# A run in record mode, for the Tier 5 figures; the defaults are HOS's standard run
# (footnote 1, p. 273):
#
#     julia --project scripts/run_standard.jl [seed] [length in τ_d] [f] [kR] [C_m]
#
# writes results/standard_f<f>_seed<seed>.h5 (or results/run_kR<kR>_Cm<Cm>_f<f>_seed<seed>.h5
# for another model) with every stream (a at every step, a₁, a₂, parity and energy every n_c
# steps, D₁/ℱ and the T, P snapshots every τ_d), the aggregator results, the summary and a
# checkpoint every 1000 τ_d. About 130 MB and 15 s for 5×10⁴ τ_d. HSO's two runs:
# `run_standard.jl 1 5e4 5.4` (the reference run) and `run_standard.jl 1 1e4 7 0.3 120`.

using HoyngDynamo

seed = length(ARGS) ≥ 1 ? parse(Int, ARGS[1]) : 1
τ_end = length(ARGS) ≥ 2 ? parse(Float64, ARGS[2]) : 5e4
f = length(ARGS) ≥ 3 ? parse(Float64, ARGS[3]) : 6.4
kR = length(ARGS) ≥ 4 ? parse(Float64, ARGS[4]) : 0.5
Cm = length(ARGS) ≥ 5 ? parse(Float64, ARGS[5]) : 100.0

spec = RunSpec(; seed, f, kR, Cm)
name = kR == 0.5 && Cm == 100.0 ? "standard_f$(f)_seed$(seed).h5" : "run_kR$(kR)_Cm$(round(Int, Cm))_f$(f)_seed$(seed).h5"
path = normpath(joinpath(@__DIR__, "..", "results", name))
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
hf = s.histogram_fit
println("χ² fits to the |a| histogram (HSO §4.3): Eq. 10 (γ, c²) = ($(round(hf.eq10.γ, digits = 2)), $(round(hf.eq10.c2, digits = 3))), $(round(hf.eq10.sigma, digits = 1))σ; " *
        "Eq. A.1 (γ, c²) = ($(round(hf.a1.γ, digits = 2)), $(round(hf.a1.c2, digits = 3))), $(round(hf.a1.sigma, digits = 1))σ over $(hf.eq10.nbins) bins")
