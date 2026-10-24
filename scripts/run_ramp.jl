# SOH's ramped run (Fig. 3, p. 123) in record mode: f rising linearly from 4 to 6.5 over
# 5×10⁴ τ_d, otherwise the standard model (kR 0.5, C_m 100, N_c 10, τ_c 0.05):
#
#     julia --project scripts/run_ramp.jl [seed] [length in τ_d] [f0] [f1]
#
# writes results/ramp_f<f0>-<f1>_seed<seed>.h5 with the a, a₁, a₂, parity, energy and D₁
# streams (no T, P snapshots), the results and the summary. viz/soh_fig03.py draws the
# polarity bar code and the windowed variability from it. About 80 MB and 15 s.

using HoyngDynamo

seed = length(ARGS) ≥ 1 ? parse(Int, ARGS[1]) : 1
τ_end = length(ARGS) ≥ 2 ? parse(Float64, ARGS[2]) : 5e4
f0 = length(ARGS) ≥ 3 ? parse(Float64, ARGS[3]) : 4.0
f1 = length(ARGS) ≥ 4 ? parse(Float64, ARGS[4]) : 6.5

spec = RunSpec(; seed, schedule = LinearRamp(f0, f1, 0.0, τ_end))
path = normpath(joinpath(@__DIR__, "..", "results", "ramp_f$(f0)-$(f1)_seed$(seed).h5"))
sim = Simulation(spec)
file = RunFile(path, spec; overwrite = true)
obs = standard_observers(sim; file, snapshots = false)
println(spec)
t = @elapsed run!(sim, steps(sim, τ_end), obs; file, checkpoint_every = steps(sim, 1000.0))
flush!(obs)
write_results!(file, obs)
s = summarize(obs, sim)
write_summary!(file, s)
close(file)
tr = value(obs.trigger)
println("wrote $path ($(round(filesize(path) / 2^20, digits = 1)) MB) in $(round(t, digits = 1)) s")
println("$(s.n_reversals) reversals ($(s.n_aborted) aborted, $(s.n_excursions) excursions); first reversal at τ = $(round(first(tr.reversal_times), digits = 1)) (SOH's figure: ≈ 1.9×10⁴); whole-run dv $(round(s.dv, digits = 3))")
