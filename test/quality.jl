# dev/DECISIONS.md D14: the hot path must be type stable and allocation-free.
#
# Register each hot function here as (f, example_args) once it exists: step! (with and
# without noise), the noise redraw and update, every aggregator's fit! and merge!. JET checks that the call is fully inferred and
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
end
