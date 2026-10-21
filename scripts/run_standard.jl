# The standard run of HOS (footnote 1, p. 273) in record mode, for the Tier 5 figures:
#
#     julia --project scripts/run_standard.jl [seed] [length in τ_d] [f]
#
# writes results/standard_f<f>_seed<seed>.h5 with every stream (a at every step, a₁, a₂,
# parity and energy every n_c steps, D₁/ℱ and the T, P snapshots every τ_d), the aggregator
# results, the summary and a checkpoint every 1000 τ_d. About 130 MB and 15 s for 5×10⁴ τ_d.

using HoyngDynamo

seed = length(ARGS) ≥ 1 ? parse(Int, ARGS[1]) : 1
τ_end = length(ARGS) ≥ 2 ? parse(Float64, ARGS[2]) : 5e4
f = length(ARGS) ≥ 3 ? parse(Float64, ARGS[3]) : 6.4

spec = RunSpec(; seed, f)
path = normpath(joinpath(@__DIR__, "..", "results", "standard_f$(f)_seed$(seed).h5"))
sim = Simulation(spec)
file = RunFile(path, spec; overwrite = true)
obs = standard_observers(sim; file)
println(spec)
t = @elapsed run!(sim, steps(sim, τ_end), obs; file, checkpoint_every = steps(sim, 1000.0))
flush!(obs)
write_results!(file, obs)
s = summarize(obs, sim)
write_summary!(file, s)
close(file)

println("wrote $path ($(round(filesize(path) / 2^20, digits = 1)) MB) in $(round(t, digits = 1)) s")
println("⟨T_r⟩ = $(round(s.mean_reversal_time, digits = 1)) ± $(round(s.mean_reversal_time_error, digits = 1)) τ_d from $(s.n_reversals) reversals " *
        "($(s.n_aborted) aborted, $(s.n_excursions) excursions); variance $(round(s.variance, digits = 3)); dv $(round(s.dv, digits = 3))")
println("D₁/ℱ = $(round(s.d1_fit.c0, digits = 3)) + $(round(s.d1_fit.c2, digits = 3)) a²  →  c² = $(round(s.c2, digits = 3)), γ = $(round(s.γ, digits = 3))")
println("theory for these shape parameters: variance $(round(s.theory.variance, digits = 3)), ⟨T_r⟩ $(round(s.theory.mean_reversal_time, digits = 1)) (Eq. 58), $(round(s.theory.mean_reversal_time_exact, digits = 1)) (Eq. 57)")
