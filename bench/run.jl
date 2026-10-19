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
# SUITE["step!"] = @benchmarkable ...        (build step 2)
# SUITE["aggregators"] = BenchmarkGroup()    (build step 4)

results = run(SUITE; verbose = true)
show(stdout, MIME"text/plain"(), results)
println()
