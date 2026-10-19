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
    SUITE["noise"] = BenchmarkGroup()
    for (name, sp) in (("exponential", ExponentialStepper(model, 0.01)), ("imex", IMEXStepper(model, 0.01)))
        st = State(model, model.b0)
        cache = Cache(model)
        noise = BinnedNoise(model, 0.01; f = 6.4, seed = 1)
        SUITE["step!"][name] = @benchmarkable step!($st, $sp, NoNoise(), $cache)
        SUITE["step!"][name * "+noise"] = @benchmarkable step!($st, $sp, $noise, $cache)
    end
    noise = BinnedNoise(model, 0.01; f = 6.4, seed = 1)
    v = zeros(model.grid.M)
    SUITE["noise"]["redraw step"] = @benchmarkable update_noise!($v, $noise, 0)
    SUITE["noise"]["quiet step"] = @benchmarkable update_noise!($v, $noise, 1)
    SUITE["noise"]["projection V00"] = @benchmarkable noise_projection($model, $v)
end
# SUITE["aggregators"] = BenchmarkGroup()    (build step 4)

tune!(SUITE)                 # pick evals per sample; without it sub-µs entries read as 0 ns
results = run(SUITE; verbose = true)
show(stdout, MIME"text/plain"(), results)
println()
