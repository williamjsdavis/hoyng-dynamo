# dev/DECISIONS.md D14: the hot path must be type stable and allocation-free.
#
# Register each hot function here as (f, example_args) once it exists: step!, every
# aggregator's fit!, the noise redraw. JET checks that the call is fully inferred and
# optimised inside HoyngDynamo; the allocation test compiles with one warm-up call and
# then measures a second call.

using JET

const HOT = Vector{Tuple{Any,Tuple}}()

let model = Model(; kR = 0.5, Cm = 100.0, N = 60), cache = Cache(model)
    for sp in (ExponentialStepper(model, 0.01), IMEXStepper(model, 0.01))
        push!(HOT, (step!, (State(model, 0.1 .* model.b0), sp, NoNoise(), cache)))
    end
end

function allocations(f, args...)
    f(args...)                 # warm-up: compile
    @allocated f(args...)
end

@testset "D14: type stability and allocations" begin
    for (f, args) in HOT
        JET.test_opt(f, Tuple{map(typeof, args)...}; target_modules = (HoyngDynamo,))
        @test allocations(f, args...) == 0
    end
end
