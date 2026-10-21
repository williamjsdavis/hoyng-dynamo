# dev/DECISIONS.md D14: the hot path must be type stable and allocation-free.
#
# Register each hot function here as (f, example_args) once it exists: step! (with and
# without noise), the noise redraw and update, every aggregator's fit! and merge!, the
# closed forms of the theory module. JET checks that the call is fully inferred and
# optimised inside HoyngDynamo; the allocation test compiles with one warm-up call and
# then measures a second call.

using JET

const HOT = Vector{Tuple{Any,Tuple}}()

let model = Model(; kR = 0.5, Cm = 100.0, N = 60), cache = Cache(model)
    noise = BinnedNoise(model, 0.01; f = 6.4, seed = 1)
    ramped = BinnedNoise(model, 0.01; schedule = LinearRamp(4.0, 6.5, 0.0, 100.0), seed = 1)
    for sp in (ExponentialStepper(model, 0.01), IMEXStepper(model, 0.01)), nz in (NoNoise(), noise)
        push!(HOT, (step!, (State(model, 0.1 .* model.b0), sp, nz, cache)))
    end
    for nz in (noise, ramped), n in (0, 1)          # redraw step and quiet step
        push!(HOT, (update_noise!, (zeros(model.grid.M), nz, n)))
    end
    for dist in (GaussianF(), UniformF(), BinaryF())
        push!(HOT, (redraw!, (BinnedNoise(model, 0.01; distribution = dist, seed = 2),)))
    end
    push!(HOT, (noise_projection, (model, zeros(model.grid.M))))
end

let aggs = (Moments(), Histogram(-3.0, 3.0, 600), AutoCorrelation(500; Δt = 0.01), Sampled(Moments(), 100),
            SchmittTrigger(; Δt = 0.01), CrossingGroups(; Δt = 0.01), ExclusionVariance(; Δt = 0.01))
    for agg in aggs
        push!(HOT, (fit!, (agg, 0.9)))
    end
    push!(HOT, (fit!, (aggs, 0.9)))                                  # a Tuple of aggregators
    push!(HOT, (fit!, (BinnedMean(0.0, 1.5, 30), (0.9, 2.5))))
    push!(HOT, (smooth!, (MovingAverage(50), 0.9)))
    # merge! of the O(1)/fixed-size aggregators; the event detectors append their event lists
    for make in (Moments, () -> Histogram(-3.0, 3.0, 600), () -> AutoCorrelation(500), () -> Sampled(Moments(), 100),
                 ExclusionVariance, () -> BinnedMean(0.0, 1.5, 30))
        push!(HOT, (merge!, (make(), make())))
    end
end

# Theory module (D12): the closed forms are allocation-free; the quadratures allocate
# (QuadGK's segment heap) and are only checked for inference. This includes the nested
# escape integral, whose inner quadgk goes through a BatchIntegrand precisely so that JET
# stays clean (see the comment in theory.jl).
const INFERRED = Vector{Tuple{Any,Tuple}}()
let th = BistableTheory(2.96, 0.246), thA = BistableTheory(4.8, 0.22; μ = 1.45 * 0.527)
    for t in (th, thA)
        for f in (diffusion_constants, peak_amplitude, peak_log_density, peak_curvature, relative_variance,
                  log_mean_reversal_time, mean_reversal_time, kramers_exponent, phi0_curvature, psi_max)
            push!(HOT, (f, (t,)))
        end
        for f in (drift, diffusion, potential, psi, dpsi, d2psi, phi)
            push!(HOT, (f, (t, 0.7)))
        end
        for f in (log_normalisation, half_moments, dipole_variation, stationary_density, log_escape_time,
                  log_mean_reversal_time_quadrature)
            push!(INFERRED, (f, (t,)))
        end
    end
    push!(HOT, (log_mean_reversal_time_asymptotic, (th,)))
    push!(HOT, (logpdf, (stationary_density(th), 0.7)))
    push!(HOT, (stationary_density(th), (0.7,)))
    push!(HOT, (shape_parameters, (0.527, 0.178, 0.0438)))
end

# Run driver: the stepping loop with the standard observers in stream mode, and the extractors.
let sim = Simulation(RunSpec(; seed = 1)), obs = standard_observers(sim)
    push!(HOT, (run!, (sim, 10, obs)))
    push!(HOT, (observe!, (obs, sim)))
    for f in (amplitude, absamplitude, overtone1, overtone2, parity, toroidal_energy, d1_sample, d1_pair)
        push!(HOT, (f, (sim,)))
    end
end

# Specialised on f and the argument types so the call is inferred: a dynamically
# dispatched call would box a scalar return value and count 16 bytes against the callee.
function allocations(f::F, args::Vararg{Any,N}) where {F,N}
    f(args...)                 # warm-up: compile
    @allocated f(args...)
end

@testset "D14: type stability and allocations" begin
    for (f, args) in HOT
        JET.test_opt(f, Tuple{map(typeof, args)...}; target_modules = (HoyngDynamo,))
        @test allocations(f, args...) == 0
    end
    for (f, args) in INFERRED
        JET.test_opt(f, Tuple{map(typeof, args)...}; target_modules = (HoyngDynamo,))
        @test (@inferred f(args...)) !== nothing
    end
    @test allocations((c0, c2) -> diffusion_from_fit(c0, c2; D0int = 3.529, prefactor = 3.22e-2), 1.36, 2.00) == 0
end
