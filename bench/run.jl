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
    for (name, sp) in (("exponential", ExponentialStepper(model, 0.01)), ("imex", IMEXStepper(model, 0.01)),
                       ("backward_euler", BackwardEulerStepper(model, 0.01)), ("crank_nicolson", CrankNicolsonStepper(model, 0.01)))
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
let aggs = (("Moments", Moments()), ("Histogram", Histogram(-3.0, 3.0, 600)), ("AutoCorrelation L=500", AutoCorrelation(500)),
            ("Sampled", Sampled(Moments(), 100)), ("SchmittTrigger", SchmittTrigger(; Δt = 0.01)),
            ("CrossingGroups", CrossingGroups(; Δt = 0.01)), ("ExclusionVariance W=5", ExclusionVariance(; Δt = 0.01)))
    SUITE["aggregators"] = BenchmarkGroup()
    for (name, agg) in aggs
        SUITE["aggregators"][name] = @benchmarkable fit!($agg, 0.9)
    end
    tuple = map(last, aggs)
    SUITE["aggregators"]["all seven (Tuple)"] = @benchmarkable fit!($tuple, 0.9)
    SUITE["aggregators"]["BinnedMean"] = @benchmarkable fit!($(BinnedMean(0.0, 1.5, 30)), (0.9, 2.5))
end

let sim = Simulation(RunSpec(; seed = 1)), obs = standard_observers(sim)
    SUITE["run!"] = BenchmarkGroup()
    SUITE["run!"]["stream mode, 100 steps, standard observers"] = @benchmarkable run!($sim, 100, $obs)
    SUITE["run!"]["100 steps, no observers"] = @benchmarkable run!($sim, 100)
end

tune!(SUITE)                 # pick evals per sample; without it sub-µs entries read as 0 ns
results = run(SUITE; verbose = true)
show(stdout, MIME"text/plain"(), results)
println()
