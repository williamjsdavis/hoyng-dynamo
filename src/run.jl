# The run driver (dev/DECISIONS.md D6, D7, D9): a `Simulation` bundles model, stepper, noise,
# state and cache; observers see the state after every step; `run!` is the single stepping
# loop shared by record mode (aggregators and HDF5 recorders) and stream mode (aggregators
# only). The stepper and the RNG never learn which mode is running, so the same seed gives
# bitwise identical trajectories and aggregator values in both (gate T2.5).

"""
    RunSpec(; kR = 0.5, Cm = 100.0, N = 60, Δτ = 0.01, f = 6.4, schedule = ConstantForcing(f),
            τc = 0.05, Nc = 10, distribution = GaussianF(), seed = 1, stepper = :exponential,
            a0 = 1.0, T = Float64)

Every parameter of a run. `schedule` is the forcing f(τ) (`ConstantForcing`, `LinearRamp`);
`stepper` is `:exponential` or `:imex`; the initial state is `a0 · b₀`, the fundamental mode
at amplitude `a0` (a = ±1 is the nonlinear equilibrium). The defaults are HOS's standard
run (footnote 1, p. 273).
"""
struct RunSpec{T<:AbstractFloat,D<:AbstractFDistribution,S<:AbstractForcing{T}}
    kR::T
    Cm::T
    N::Int
    Δτ::T
    τc::T
    Nc::Int
    distribution::D
    schedule::S
    seed::UInt64
    stepper::Symbol
    a0::T
end

function RunSpec(; kR = 0.5, Cm = 100.0, N::Integer = 60, Δτ = 0.01, f = 6.4, schedule = nothing, τc = 0.05,
                 Nc::Integer = 10, distribution::AbstractFDistribution = GaussianF(), seed::Integer = 1,
                 stepper::Symbol = :exponential, a0 = 1.0, T::Type{<:AbstractFloat} = Float64)
    stepper in (:exponential, :imex) || throw(ArgumentError("stepper must be :exponential or :imex, got $stepper"))
    sched = schedule === nothing ? ConstantForcing(T(f)) : schedule
    sched isa AbstractForcing{T} || throw(ArgumentError("the forcing schedule must have element type $T"))
    RunSpec{T,typeof(distribution),typeof(sched)}(T(kR), T(Cm), Int(N), T(Δτ), T(τc), Int(Nc), distribution, sched,
                                                   UInt64(seed), stepper, T(a0))
end

Base.show(io::IO, s::RunSpec{T}) where {T} =
    print(io, "RunSpec{$T}(kR = $(s.kR), C_m = $(s.Cm), N = $(s.N), Δτ = $(s.Δτ), $(s.schedule), τ_c = $(s.τc), N_c = $(s.Nc), ",
          "$(nameof(typeof(s.distribution))), seed = $(s.seed), $(s.stepper), a0 = $(s.a0))")

"""
    Diagnostics(model; K = 2)

Precomputed projections for the observers: the adjoint projection vectors of the first `K`
overtones (`ã[:, k]` so that a_k = (b̃_k, b) = `dot(ã[:, k], b)`, Table I order, one member
per conjugate pair), their eigenvalues at C₀, the weights `d1w = w P̃₀²` of the D₁ integrand
(HOS Eq. 37) and the toroidal fundamental mode `T0`.
"""
struct Diagnostics{T<:AbstractFloat}
    ã::Matrix{Complex{T}}
    λ::Vector{Complex{T}}
    d1w::Vector{T}
    T0::Vector{T}
end

function Diagnostics(model::Model{T}; K::Integer = 2) where {T}
    g, m = model.grid, model.modes
    idx = [i for i in eachindex(m.λ) if imag(m.λ[i]) ≥ 0]           # one member per conjugate pair
    length(idx) ≥ K + 1 || throw(ArgumentError("only $(length(idx) - 1) overtones available, asked for $K"))
    ov = idx[2:K+1]
    ã = Matrix{Complex{T}}(undef, 2g.M, K)
    for (k, i) in enumerate(ov)
        ã[:, k] = g.W .* m.b̃[:, i]
    end
    P̃0 = view(model.b̃0, 1:g.M)
    Diagnostics{T}(ã, m.λ[ov], g.w .* P̃0 .^ 2, model.b0[g.M+1:2g.M])
end

"""
    Simulation(spec::RunSpec)

Model, stepper, binned noise, state (a0 · b₀), cache and diagnostics for one run.
"""
struct Simulation{T<:AbstractFloat,R<:RunSpec{T},S<:AbstractStepper{T},Nz<:BinnedNoise{T}}
    spec::R
    model::Model{T}
    stepper::S
    noise::Nz
    state::State{T}
    cache::Cache{T}
    diag::Diagnostics{T}
end

function Simulation(spec::RunSpec{T}) where {T}
    model = Model(; kR = spec.kR, Cm = spec.Cm, N = spec.N, T)
    noise = BinnedNoise(model, spec.Δτ; schedule = spec.schedule, τc = spec.τc, Nc = spec.Nc,
                        distribution = spec.distribution, seed = spec.seed)
    state = State(model, spec.a0 .* model.b0)
    if spec.stepper === :imex
        return Simulation(spec, model, IMEXStepper(model, spec.Δτ), noise, state, Cache(model), Diagnostics(model))
    else
        return Simulation(spec, model, ExponentialStepper(model, spec.Δτ), noise, state, Cache(model), Diagnostics(model))
    end
end

Base.show(io::IO, sim::Simulation) = print(io, "Simulation(", sim.spec, "; n = ", sim.state.n, ", τ = ", current_time(sim), ", a = ", sim.state.a, ")")

timestep(sim::Simulation) = timestep(sim.stepper)
"τ = n Δτ of the current state."
current_time(sim::Simulation) = sim.state.n * timestep(sim.stepper)
"Number of steps in `τ` diffusion times (rounded)."
steps(sim::Simulation, τ::Real) = round(Int, τ / timestep(sim.stepper))

# ---------------------------------------------------------------------------------------
# Extractors: scalar (or vector) functions of the simulation, allocation-free
# ---------------------------------------------------------------------------------------

"The dipole amplitude a = (b̃₀, b)."
@inline amplitude(sim::Simulation) = sim.state.a
"|a|, the rectified dipole amplitude (HSO Eq. 12)."
@inline absamplitude(sim::Simulation) = abs(sim.state.a)
"The complex overtone amplitude a_k = (b̃_k, b), k = 1, 2, … (HOS Eq. 15)."
@inline overtone(sim::Simulation, k::Integer) = dot(view(sim.diag.ã, :, k), sim.state.b)
@inline overtone1(sim::Simulation) = overtone(sim, 1)
@inline overtone2(sim::Simulation) = overtone(sim, 2)
"Parity 𝒫 of the toroidal field (HOS Eq. 28): −1 antisymmetric (dipolar), +1 symmetric."
@inline parity(sim::Simulation) = toroidal_parity(sim.model.grid, sim.state.b)
"The toroidal field T at the interior nodes (a view)."
@inline toroidal(sim::Simulation) = view(sim.state.b, sim.model.grid.M+1:2sim.model.grid.M)
"The poloidal potential P at the interior nodes (a view)."
@inline poloidal(sim::Simulation) = view(sim.state.b, 1:sim.model.grid.M)

"∫ T² sin θ dθ, the 'energy' of HOS Fig. 2."
@inline function toroidal_energy(sim::Simulation{T}) where {T}
    w, b, M = sim.model.grid.w, sim.state.b, sim.model.grid.M
    s = zero(T)
    @inbounds @simd for j in 1:M
        s += w[j] * b[M+j]^2
    end
    return s
end

"""
    d1_sample(sim)

The instantaneous D₁/ℱ = ∫ sin θ P̃₀² T₊² dθ with T₊ = T − a T₀ (HOS Eq. 37 before the
conditional average; Fig. 8 plots it against |a| once per τ_d).
"""
@inline function d1_sample(sim::Simulation{T}) where {T}
    d, b, a, M = sim.diag, sim.state.b, sim.state.a, sim.model.grid.M
    s = zero(T)
    @inbounds @simd for j in 1:M
        s += d.d1w[j] * (b[M+j] - a * d.T0[j])^2
    end
    return s
end

"(|a|, D₁/ℱ) for a `BinnedMean` (the D₁ measurement of HOS Fig. 8)."
@inline d1_pair(sim::Simulation) = (absamplitude(sim), d1_sample(sim))

# ---------------------------------------------------------------------------------------
# Observers
# ---------------------------------------------------------------------------------------

"""
    Observer

Something that sees the simulation after each step through `observe!(o, sim)`. Observers are
composed in Tuples or NamedTuples; `Probe` feeds an aggregator, `Recorder` and
`SnapshotRecorder` (io.jl) write HDF5 streams.
"""
abstract type Observer end

"""
    Probe(f, aggregator, every = 1)

Every `every` steps, `fit!(aggregator, f(sim))`. Steps are counted from the start of the
run (n = every, 2·every, …), so a replay of a recorded stream at the same cadence
reproduces the aggregator's values.
"""
struct Probe{F,A<:Aggregator} <: Observer
    f::F
    agg::A
    every::Int
    function Probe{F,A}(f, agg, every::Integer) where {F,A<:Aggregator}
        every ≥ 1 || throw(ArgumentError("need every ≥ 1, got $every"))
        new{F,A}(f, agg, Int(every))
    end
end

Probe(f, agg::Aggregator, every::Integer = 1) = Probe{typeof(f),typeof(agg)}(f, agg, every)

@inline function observe!(p::Probe, sim::Simulation)
    sim.state.n % p.every == 0 && fit!(p.agg, p.f(sim))
    return p
end

observe!(t::Tuple, sim::Simulation) = (foreach(o -> observe!(o, sim), t); t)
observe!(nt::NamedTuple, sim::Simulation) = (observe!(values(nt), sim); nt)
observe!(::Nothing, ::Simulation) = nothing

value(p::Probe) = value(p.agg)
value(::Observer) = nothing
value(nt::NamedTuple) = map(value, nt)

"Flush buffered output (recorders); a no-op for probes."
function flush! end
flush!(::Probe) = nothing
flush!(t::Tuple) = (foreach(flush!, t); nothing)
flush!(nt::NamedTuple) = flush!(values(nt))
flush!(::Nothing) = nothing

"Write a checkpoint of `sim` and the observers' aggregators to a run file (io.jl)."
function checkpoint! end
checkpoint!(::Nothing, ::Simulation, observers) = throw(ArgumentError("checkpoints need a run file"))

# ---------------------------------------------------------------------------------------
# The loop
# ---------------------------------------------------------------------------------------

"""
    run!(sim, nsteps, observers = nothing; file = nothing, checkpoint_every = 0)

Advance `sim` by `nsteps` steps, calling `observe!(observers, sim)` after each. With
`checkpoint_every > 0` (steps, a multiple of the redraw interval n_c) a checkpoint is written
to `file` whenever the step count is a multiple of it. Allocation-free in stream mode.
"""
function run!(sim::Simulation, nsteps::Integer, observers = nothing; file = nothing, checkpoint_every::Integer = 0)
    nsteps ≥ 0 || throw(ArgumentError("nsteps must be ≥ 0, got $nsteps"))
    if checkpoint_every > 0
        file === nothing && throw(ArgumentError("checkpoints need a run file"))
        nc = redraw_interval(sim.noise)
        checkpoint_every % nc == 0 || throw(ArgumentError("checkpoint_every = $checkpoint_every must be a multiple of the redraw interval n_c = $nc"))
    end
    st = sim.state
    for _ in 1:nsteps
        step!(st, sim.stepper, sim.noise, sim.cache)
        observe!(observers, sim)
        if checkpoint_every > 0 && st.n % checkpoint_every == 0
            checkpoint!(file, sim, observers)
        end
    end
    return sim
end

# ---------------------------------------------------------------------------------------
# The standard observer set and its summary
# ---------------------------------------------------------------------------------------

"""
    standard_observers(sim; file = nothing, snapshots = true, stream_eltype = Float64, attach = false)

The statistics of the papers as a NamedTuple of probes: the Schmitt-trigger reversal
detector, the papers' crossing-grouping rule and the exclusion-window variance on a at
every step; the moments of |a| (dv, HSO Eq. 12), the histogram of a and the D₁ measurement
(HOS Fig. 8) at 1 τ_d; the parity histogram, the toroidal-energy moments and the
autocorrelation of |a| (100 lags of n_c steps). With a `RunFile`, recorders for the
streams a (every step), a₁, a₂, parity and energy (every n_c steps), D₁/ℱ (every τ_d) and,
if `snapshots`, T and P (every τ_d) are added: record mode. `attach = true` reopens
existing streams (resume).
"""
function standard_observers(sim::Simulation{T}; file = nothing, snapshots::Bool = true,
                            stream_eltype::Type{<:AbstractFloat} = T, attach::Bool = false) where {T}
    Δτ = timestep(sim)
    nd = round(Int, 1 / Δτ)
    abs(nd * Δτ - 1) ≤ 1e-9 || throw(ArgumentError("1/Δτ = $(1 / Δτ) must be an integer for the τ_d cadence"))
    nc = redraw_interval(sim.noise)
    probes = (
        trigger = Probe(amplitude, SchmittTrigger(; Δt = Δτ)),
        groups = Probe(amplitude, CrossingGroups(; Δt = Δτ)),
        exclusion = Probe(amplitude, ExclusionVariance(; Δt = Δτ)),
        moments = Probe(amplitude, Moments()),
        abs_moments = Probe(absamplitude, Moments(), nd),
        histogram = Probe(amplitude, Histogram(T(-2.5), T(2.5), 500), nd),
        d1 = Probe(d1_pair, BinnedMean(zero(T), T(1.5), 30), nd),
        parity = Probe(parity, Histogram(-one(T), one(T), 200), nc),
        energy = Probe(toroidal_energy, Moments(), nd),
        autocorrelation = Probe(absamplitude, AutoCorrelation(100; Δt = nc * Δτ, T), nc),
    )
    file === nothing && return probes
    C = Complex{stream_eltype}
    recorders = (
        a = Recorder(file, "a", amplitude, 1; sim, eltype = stream_eltype, attach),
        a1 = Recorder(file, "a1", overtone1, nc; sim, eltype = C, attach),
        a2 = Recorder(file, "a2", overtone2, nc; sim, eltype = C, attach),
        parity_series = Recorder(file, "parity", parity, nc; sim, eltype = stream_eltype, attach),
        energy_series = Recorder(file, "energy", toroidal_energy, nc; sim, eltype = stream_eltype, attach),
        d1_series = Recorder(file, "d1", d1_sample, nd; sim, eltype = stream_eltype, attach),
    )
    snapshots || return merge(probes, recorders)
    snaps = (
        T_snapshots = SnapshotRecorder(file, "T", toroidal, nd; sim, eltype = stream_eltype, attach),
        P_snapshots = SnapshotRecorder(file, "P", poloidal, nd; sim, eltype = stream_eltype, attach),
    )
    return merge(probes, recorders, snaps)
end

"""
    summarize(observers, sim) -> NamedTuple

The headline statistics of a run with the standard observers: reversal counts and ⟨T_r⟩
with its Poisson error (Schmitt trigger, and the papers' grouping rule), the exclusion
variance of §2, dv (HSO Eq. 12), the D₁/ℱ quadratic fit of the bin averages for |a| ≤ 1.1
(HOS Eq. 61; unweighted as in the paper, with the count-weighted fit and the mean over
|a| ≤ 0.2 alongside) and the (c², γ) it implies with D₀ from Eq. 36 (Eqs. 62–63), and the
Eq. 56/58 theory for those shape parameters. The D₁ fit uses f at the last redraw
(constant forcing: f).
"""
function summarize(obs::NamedTuple, sim::Simulation{T}) where {T}
    tr, cg, ex, am, d1 = value(obs.trigger), value(obs.groups), value(obs.exclusion), value(obs.abs_moments), value(obs.d1)
    c = sim.model.constants
    ℱ = diffusion_prefactor(sim.noise.f, sim.noise.τc, sim.noise.Nc)
    sel = (d1.centers .<= 1.1) .& (d1.n .> 0)
    # the papers fit the histogram of bin averages (HOS Fig. 8, p. 291): unweighted least squares
    # on the bin means; the count-weighted fit and the direct small-|a| mean are kept alongside
    # because D₁(a) is only approximately quadratic and the intercept depends on the weighting
    c0, c2 = even_quadratic_fit(d1.centers[sel], d1.mean[sel], ones(T, count(sel)))
    w0, w2 = even_quadratic_fit(d1.centers[sel], d1.mean[sel], float.(d1.n[sel]))
    near = (d1.centers .<= 0.2) .& (d1.n .> 0)
    d1_zero = sum(d1.n[near] .* d1.mean[near]) / max(sum(d1.n[near]), 1)
    Δ = diffusion_from_fit(c0, c2; D0int = c.D0int, prefactor = ℱ)
    sp = shape_parameters(c.Λ, Δ.Δ0, Δ.Δ1)
    nan = T(NaN)
    theory = if isfinite(sp.γ) && sp.γ > 1 && sp.c2 > 0
        th = BistableTheory(sp.γ, sp.c2; Λ = c.Λ)
        (variance = relative_variance(th), mean_reversal_time = mean_reversal_time(th),
         mean_reversal_time_exact = mean_reversal_time_quadrature(th), dv = dipole_variation(th))
    else
        (variance = nan, mean_reversal_time = nan, mean_reversal_time_exact = nan, dv = nan)
    end
    return (time = tr.time, n_reversals = tr.n_reversals, n_aborted = tr.n_aborted, n_excursions = tr.n_excursions,
            mean_reversal_time = tr.mean_interval, mean_reversal_time_error = tr.mean_interval / sqrt(max(tr.n_reversals, 1)),
            n_reversals_grouped = cg.n_reversals, n_aborted_grouped = cg.n_aborted, mean_reversal_time_grouped = cg.mean_interval,
            variance = ex.variance, variance_intervals = ex.intervals, dv = am.var / am.mean^2, mean_abs = am.mean,
            d1_fit = (c0 = c0, c2 = c2), d1_fit_weighted = (c0 = w0, c2 = w2), d1_zero = d1_zero,
            prefactor = ℱ, c2 = sp.c2, γ = sp.γ, Δ0 = Δ.Δ0, Δ1 = Δ.Δ1, theory = theory)
end
