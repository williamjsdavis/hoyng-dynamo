# Benchmarks for the hot path (dev/DECISIONS.md D14): ns per step and per aggregator
# update. Usage, from the repository root:
#     julia bench/run.jl
# Compare two commits by saving the results (BenchmarkTools.save) and using judge().
using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(path = dirname(@__DIR__))
Pkg.instantiate()

using BenchmarkTools
using HoyngDynamo

const SUITE = BenchmarkGroup()

let model = Model(; kR = 0.5, Cm = 100.0, N = 60)
    SUITE["step!"] = BenchmarkGroup()
    for (name, sp) in (("exponential", ExponentialStepper(model, 0.01)), ("imex", IMEXStepper(model, 0.01)))
        st = State(model, model.b0)
        cache = Cache(model)
        SUITE["step!"][name] = @benchmarkable step!($st, $sp, NoNoise(), $cache)
    end
end
# SUITE["aggregators"] = BenchmarkGroup()    (build step 4)

results = run(SUITE; verbose = true)
show(stdout, MIME"text/plain"(), results)
println()
