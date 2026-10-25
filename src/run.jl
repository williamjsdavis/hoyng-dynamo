# The run driver (dev/DECISIONS.md D6, D7, D9): a `Simulation` bundles model, stepper, noise,
# state and cache; observers see the state after every step; `run!` is the single stepping
# loop shared by record mode (aggregators and HDF5 recorders) and stream mode (aggregators
# only). The stepper and the RNG never learn which mode is running, so the same seed gives
# bitwise identical trajectories and aggregator values in both (gate T2.5).

"The steppers `RunSpec` accepts; `make_stepper` builds them."
const STEPPERS = (:exponential, :imex, :backward_euler, :crank_nicolson, :backward_euler_implicit_noise,
                  :crank_nicolson_implicit_noise)
"The quenching laws `RunSpec` accepts (`AbstractQuenching`); `make_quenching` builds them."
const QUENCHINGS = (:amplitude, :energy, :inverse, :local, :amplitude_floored, :energy_floored, :local_floored)

"""
    RunSpec(; kR = 0.5, Cm = 100.0, N = 60, Δτ = 0.01, f = 6.4, schedule = ConstantForcing(f),
            τc = 0.05, Nc = 10, distribution = GaussianF(), seed = 1, stepper = :exponential,
            quenching = :amplitude, a0 = 1.0, T = Float64)

Every parameter of a run. `schedule` is the forcing f(τ) (`ConstantForcing`, `LinearRamp`);
`stepper` is one of `STEPPERS`: `:exponential` (production), `:imex`, and the SG-1 sensitivity
steppers `:backward_euler`, `:crank_nicolson`, `:backward_euler_implicit_noise`,
`:crank_nicolson_implicit_noise` (`ThetaStepper`); `quenching` is one of `QUENCHINGS`: `:amplitude`
(the papers' law, q = q₀a²), the SG-2 alternatives `:energy`, `:inverse`, `:local`, or the
`_floored` variants with the α-effect suppressed but never reversed (`AbstractQuenching`; the
local law's constant is calibrated by `Simulation`); the initial
state is `a0 · b₀`, the fundamental mode at amplitude `a0` (a = ±1 is the nonlinear equilibrium). The defaults are HOS's standard
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
    quenching::Symbol
    a0::T
end

function RunSpec(; kR = 0.5, Cm = 100.0, N::Integer = 60, Δτ = 0.01, f = 6.4, schedule = nothing, τc = 0.05,
                 Nc::Integer = 10, distribution::AbstractFDistribution = GaussianF(), seed::Integer = 1,
                 stepper::Symbol = :exponential, quenching::Symbol = :amplitude, a0 = 1.0, T::Type{<:AbstractFloat} = Float64)
    stepper in STEPPERS || throw(ArgumentError("stepper must be one of $STEPPERS, got $stepper"))
    quenching in QUENCHINGS || throw(ArgumentError("quenching must be one of $QUENCHINGS, got $quenching"))
    sched = schedule === nothing ? ConstantForcing(T(f)) : schedule
    sched isa AbstractForcing{T} || throw(ArgumentError("the forcing schedule must have element type $T"))
    RunSpec{T,typeof(distribution),typeof(sched)}(T(kR), T(Cm), Int(N), T(Δτ), T(τc), Int(Nc), distribution, sched,
                                                   UInt64(seed), stepper, quenching, T(a0))
end

"""
    RunSpec(spec::RunSpec; kwargs...)

A copy of `spec` with the given fields replaced: `seed`, `f` (a new `ConstantForcing`),
`schedule`, `τc`, `Nc`, `N`, `Δτ`, `kR`, `Cm`, `distribution`, `stepper`, `quenching`, `a0`.
"""
function RunSpec(spec::RunSpec{T}; f = nothing, schedule = nothing, kR = spec.kR, Cm = spec.Cm, N = spec.N,
                 Δτ = spec.Δτ, τc = spec.τc, Nc = spec.Nc, distribution = spec.distribution, seed = spec.seed,
                 stepper = spec.stepper, quenching = spec.quenching, a0 = spec.a0) where {T}
    sched = schedule !== nothing ? schedule : (f !== nothing ? ConstantForcing(T(f)) : spec.schedule)
    return RunSpec(; kR, Cm, N, Δτ, schedule = sched, τc, Nc, distribution, seed, stepper, quenching, a0, T)
end

Base.show(io::IO, s::RunSpec{T}) where {T} =
    print(io, "RunSpec{$T}(kR = $(s.kR), C_m = $(s.Cm), N = $(s.N), Δτ = $(s.Δτ), $(s.schedule), τ_c = $(s.τc), N_c = $(s.Nc), ",
          "$(nameof(typeof(s.distribution))), seed = $(s.seed), $(s.stepper), $(s.quenching), a0 = $(s.a0))")

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
    return Simulation(spec, model, make_stepper(spec, model), noise, state, Cache(model), Diagnostics(model))
end

"The quenching law `spec.quenching` names (see `QUENCHINGS`); the local law is calibrated for `model` at `spec.Δτ`."
function make_quenching(spec::RunSpec{T}, model::Model{T}) where {T}
    q = spec.quenching
    q === :amplitude && return AmplitudeQuenching()
    q === :energy && return EnergyQuenching()
    q === :inverse && return InverseQuenching()
    q === :local && return calibrate(LocalQuenching, model, spec.Δτ)
    q === :amplitude_floored && return AmplitudeQuenching(; floored = true)
    q === :energy_floored && return EnergyQuenching(; floored = true)
    q === :local_floored && return calibrate(LocalQuenching, model, spec.Δτ; floored = true)
    throw(ArgumentError("unknown quenching law $q"))
end

"The stepper `spec.stepper` names, built for `model` with the law `spec.quenching` (see `STEPPERS`)."
function make_stepper(spec::RunSpec{T}, model::Model{T}) where {T}
    s, Δτ = spec.stepper, spec.Δτ
    quenching = make_quenching(spec, model)
    s === :exponential && return ExponentialStepper(model, Δτ; quenching)
    s === :imex && return IMEXStepper(model, Δτ; quenching)
    s === :backward_euler && return ThetaStepper(model, Δτ; θ = 1, quenching)
    s === :crank_nicolson && return ThetaStepper(model, Δτ; θ = 1 // 2, quenching)
    s === :backward_euler_implicit_noise && return ThetaStepper(model, Δτ; θ = 1, implicit_noise = true, quenching)
    s === :crank_nicolson_implicit_noise && return ThetaStepper(model, Δτ; θ = 1 // 2, implicit_noise = true, quenching)
    throw(ArgumentError("unknown stepper $s"))
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

"""
    noise_diffusion_sample(sim)

The exact diffusion coefficient of a for the current field under the binned noise,
D_b = τ_c f² Σ_k (Σ_{j∈k} ã_j g_j T_j)² (the formula of `binned_D0` with the current T in place
of T₀): the variance per unit time of the noise-driven increments of a, whose conditional
mean over |a| is what HOS's D(a) = D₀a² + D₁(a) (Eqs. 36–37) represents in the δ-correlated
limit N_c → ∞. Allocation-free (the bins are contiguous in θ).
"""
@inline function noise_diffusion_sample(sim::Simulation{T}) where {T}
    nz, ã, b, M = sim.noise, sim.model.ã, sim.state.b, sim.model.grid.M
    total = zero(T)
    acc = zero(T)
    current = @inbounds nz.bin[1]
    @inbounds for j in 1:M
        k = nz.bin[j]
        if k != current
            total += acc * acc
            acc = zero(T)
            current = k
        end
        acc += ã[j] * nz.g[j] * b[M+j]
    end
    total += acc * acc
    return nz.τc * nz.f^2 * total
end

"""
    continuum_diffusion_sample(sim)

The δ-correlated-limit counterpart of `noise_diffusion_sample`: ℱ ∫ sin θ P̃₀² T² dθ with the
full toroidal field, i.e. HOS's D₀a² + D₁ (Eqs. 36–37, cross term included) evaluated on the
current field, ℱ = πf²τ_c/(2N_c²) at the current f.
"""
@inline function continuum_diffusion_sample(sim::Simulation{T}) where {T}
    d, b, M, nz = sim.diag, sim.state.b, sim.model.grid.M, sim.noise
    s = zero(T)
    @inbounds @simd for j in 1:M
        s += d.d1w[j] * b[M+j]^2
    end
    return diffusion_prefactor(nz.f, nz.τc, nz.Nc) * s
end

"(|a|, D_b) and (|a|, D_c) for `BinnedMean`s: the exact binned and the continuum diffusion of a against |a|."
@inline noise_diffusion_pair(sim::Simulation) = (absamplitude(sim), noise_diffusion_sample(sim))
@inline continuum_diffusion_pair(sim::Simulation) = (absamplitude(sim), continuum_diffusion_sample(sim))

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
        histogram = Probe(amplitude, Histogram(T(-2.5), T(2.5), 1000), nd),     # 0.005 bins: ~100 equal-population bins for the χ² fit
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
(HOS Eq. 61, weighted by the bin counts as in the Fig. 8 caption; the unweighted fit and
the mean over |a| ≤ 0.2 alongside) and the (c², γ) it implies with D₀ from Eq. 36
(Eqs. 62–63; `c2_unweighted`, `γ_unweighted` from the unweighted fit), the χ² fits of
HSO Eq. 10 and Eq. A.1 to the |a| histogram (`histogram_fit.eq10`, `.a1`; HSO §4.3), and
the Eq. 56/58 theory for the primary shape parameters. ℱ uses f at the last redraw
(constant forcing: f).
"""
function summarize(obs::NamedTuple, sim::Simulation{T}) where {T}
    tr, cg, ex, am, d1 = value(obs.trigger), value(obs.groups), value(obs.exclusion), value(obs.abs_moments), value(obs.d1)
    c = sim.model.constants
    ℱ = diffusion_prefactor(sim.noise.f, sim.noise.τc, sim.noise.Nc)
    sel = (d1.centers .<= 1.1) .& (d1.n .> 0)
    # HOS Fig. 8 caption: a least-squares quadratic through the histogram of bin averages,
    # "weighted by the number of points in the bin", for |a| ≤ 1.1. The unweighted fit and the
    # direct small-|a| mean are kept alongside because D₁(a) is only approximately quadratic
    # and the intercept depends on the weighting (DISCREPANCIES.md DR-5)
    c0, c2 = even_quadratic_fit(d1.centers[sel], d1.mean[sel], float.(d1.n[sel]))
    u0, u2 = even_quadratic_fit(d1.centers[sel], d1.mean[sel], ones(T, count(sel)))
    near = (d1.centers .<= 0.2) .& (d1.n .> 0)
    d1_zero = sum(d1.n[near] .* d1.mean[near]) / max(sum(d1.n[near]), 1)
    Δ = diffusion_from_fit(c0, c2; D0int = c.D0int, prefactor = ℱ)
    sp = shape_parameters(c.Λ, Δ.Δ0, Δ.Δ1)
    Δu = diffusion_from_fit(u0, u2; D0int = c.D0int, prefactor = ℱ)
    spu = shape_parameters(c.Λ, Δu.Δ0, Δu.Δ1)
    # HSO §4.3: χ² fits of Eq. 10 (κ = 0) and of Eq. A.1 (κ = μ/Λ − 1 of this model) to the
    # histogram of |a| at 1 τ_d, about 100 equal-population bins
    h = value(obs.histogram)
    fe, fc = fold_histogram(h.edges, h.counts)
    hfit = if sum(fc) ≥ 1000
        (eq10 = density_fit(fe, fc, c.Λ), a1 = density_fit(fe, fc, c.Λ, c.μ))
    else
        none = (γ = T(NaN), c2 = T(NaN), chi2 = T(NaN), ν = 0, sigma = T(NaN), nbins = 0, n = sum(fc))
        (eq10 = none, a1 = none)
    end
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
            d1_fit = (c0 = c0, c2 = c2), d1_fit_unweighted = (c0 = u0, c2 = u2), d1_zero = d1_zero,
            prefactor = ℱ, Λ = c.Λ, μ = c.μ, D0int = c.D0int, c2 = sp.c2, γ = sp.γ, Δ0 = Δ.Δ0, Δ1 = Δ.Δ1,
            c2_unweighted = spu.c2, γ_unweighted = spu.γ, histogram_fit = hfit, theory = theory)
end

# ---------------------------------------------------------------------------------------
# Ensembles of seeds
# ---------------------------------------------------------------------------------------

"""
    pool!(target::NamedTuple, others) -> target

Merge the aggregators of the observer sets in `others` into those of `target` (stream-mode
sets of probes with the same keys, e.g. several seeds of `standard_observers`). The
`intervals` of a pooled `SchmittTrigger` straddle the seams between runs; take them per run
(`ensemble` does).
"""
function pool!(target::NamedTuple, others)
    for obs in others, k in keys(target)
        merge!(target[k].agg, obs[k].agg)
    end
    return target
end

"""
    ensemble(spec::RunSpec, τ; seeds = 1:3, observers = standard_observers) -> (runs, pooled, summary, intervals)

`spec` with each of `seeds`, run for `τ` diffusion times in stream mode with the observer
set `observers(sim)` (the standard observers by default; a set that extends them with
further probes, `sim -> merge(standard_observers(sim), more(sim))`, must keep their keys
for `summarize`), in parallel over the available threads (the result does not depend
on the thread count: every seed is its own simulation and the aggregators are pooled in
seed order). `runs[i]` is `(seed, sim, observers, summary, seconds)`; `pooled` holds the
pooled aggregators (a copy of the first set), `summary = summarize(pooled, sim)` and
`intervals` the polarity intervals of all runs, each run's own.
"""
function ensemble(spec::RunSpec, τ::Real; seeds = 1:3, observers = standard_observers)
    seeds = collect(seeds)
    sims = [Simulation(RunSpec(spec; seed)) for seed in seeds]
    obs = [observers(sim) for sim in sims]
    seconds = zeros(length(sims))
    Threads.@threads for i in eachindex(sims)
        seconds[i] = @elapsed run!(sims[i], steps(sims[i], τ), obs[i])
    end
    runs = [(seed = seeds[i], sim = sims[i], observers = obs[i], summary = summarize(obs[i], sims[i]),
             seconds = seconds[i]) for i in eachindex(sims)]
    pooled = pool!(deepcopy(obs[1]), obs[2:end])
    intervals = reduce(vcat, (value(r.observers.trigger).intervals for r in runs))
    return (runs = runs, pooled = pooled, summary = summarize(pooled, sims[1]), intervals = intervals)
end

"""
    segmented(spec::RunSpec, τ_segment, nsegments; seeds = 1:1) -> (runs, segments, first_reversal)

Run `spec` with each seed for `nsegments` consecutive segments of `τ_segment` diffusion
times, a fresh set of standard observers per segment (windowed statistics, e.g. of a ramped
forcing, SOH Fig. 3) and one Schmitt trigger over the whole run; seeds in parallel over the
threads. `runs[k]` is `(seed, sim, whole, segments)` with `whole` the whole-run trigger value
and `segments` that seed's per-segment summaries (ℱ at the segment's end); `segments[i]` is
the summary of segment i pooled over the seeds (its ℱ, hence its D₁-derived shape parameters,
use f at the end of the whole run); `first_reversal[k]` the time of seed k's first reversal
(NaN if none).
"""
function segmented(spec::RunSpec, τ_segment::Real, nsegments::Integer; seeds = 1:1)
    nsegments ≥ 1 || throw(ArgumentError("need at least one segment"))
    seeds = collect(seeds)
    sims = [Simulation(RunSpec(spec; seed)) for seed in seeds]
    obs = [[standard_observers(sim) for _ in 1:nsegments] for sim in sims]
    whole = [Probe(amplitude, SchmittTrigger(; Δt = timestep(sim))) for sim in sims]
    summaries = [Vector{Any}(undef, nsegments) for _ in sims]
    Threads.@threads for k in eachindex(sims)
        n = steps(sims[k], τ_segment)
        for i in 1:nsegments
            run!(sims[k], n, (whole[k], obs[k][i]))
            summaries[k][i] = summarize(obs[k][i], sims[k])
        end
    end
    runs = [(seed = seeds[k], sim = sims[k], whole = value(whole[k]), segments = [summaries[k][i] for i in 1:nsegments])
            for k in eachindex(sims)]
    segments = [summarize(pool!(deepcopy(obs[1][i]), (obs[k][i] for k in 2:length(sims))), sims[1]) for i in 1:nsegments]
    first_reversal = [isempty(r.whole.reversal_times) ? NaN : first(r.whole.reversal_times) for r in runs]
    return (runs = runs, segments = segments, first_reversal = first_reversal)
end


"""
    window_run(spec, τ0, τ1; every = 1) -> NamedTuple

Re-run `spec` from τ = 0 and record the window τ0 < τ ≤ τ1 at full resolution (every
`every` steps): the toroidal field `T` and the poloidal potential `P` at the nodes (M × n),
the dipole amplitude `a`, the overtones `a1`, `a2`, the `parity`, the toroidal `energy` and
the sample `time`s, with the grid's `θ` and the finished `sim`. The seed fixes the
trajectory, so a window chosen from the event list of a run file reproduces that run's
streams bitwise (HOS Figs. 3 and 4 are 6 τ_d windows at full resolution; `scripts/run_events.jl`).
"""
function window_run(spec::RunSpec{T}, τ0::Real, τ1::Real; every::Integer = 1) where {T}
    every ≥ 1 || throw(ArgumentError("need every ≥ 1, got $every"))
    sim = Simulation(spec)
    n0, n1 = steps(sim, τ0), steps(sim, τ1)
    0 ≤ n0 < n1 || throw(ArgumentError("need 0 ≤ τ0 < τ1, got $τ0, $τ1"))
    run!(sim, n0)
    nsamp = (n1 - n0) ÷ every
    M = sim.model.grid.M
    Tf, Pf = Matrix{T}(undef, M, nsamp), Matrix{T}(undef, M, nsamp)
    a, 𝒫, energy, time = (Vector{T}(undef, nsamp) for _ in 1:4)
    a1, a2 = Vector{Complex{T}}(undef, nsamp), Vector{Complex{T}}(undef, nsamp)
    for j in 1:nsamp
        run!(sim, every)
        Tf[:, j] .= toroidal(sim)
        Pf[:, j] .= poloidal(sim)
        a[j], a1[j], a2[j] = amplitude(sim), overtone1(sim), overtone2(sim)
        𝒫[j], energy[j], time[j] = parity(sim), toroidal_energy(sim), current_time(sim)
    end
    return (time = time, θ = sim.model.grid.θ, T = Tf, P = Pf, a = a, a1 = a1, a2 = a2, parity = 𝒫, energy = energy, sim = sim)
end
