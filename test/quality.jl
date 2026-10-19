# dev/DECISIONS.md D14: the hot path must be type stable and allocation-free.
#
# Register each hot function here as (f, example_args) once it exists: step! (with and
# without noise), the noise redraw and update, every aggregator's fit!. JET checks that the call is fully inferred and
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
