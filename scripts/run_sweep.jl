# The forcing sweep of HOS Table II (p. 292) in stream mode, pooled over seeds:
#
#     julia --project -t 4 scripts/run_sweep.jl [seeds=3] [scale=1]
#
# f = 4.4, 5.4, 6.4, 7.4, 8.4 at N_c = 10, τ_c = 0.05; 5×10⁴ τ_d per seed as in the paper,
# except 2×10⁵ τ_d per seed at f = 4.4, where the paper had only ~5 reversals (`scale` multiplies
# both lengths). Prints and writes results/table2.md: our Table II with the KS test of the
# polarity intervals (gate T4.11) and the theory at the measured (c², γ). Gates T4.5–T4.7,
# T4.11, T4.18 are asserted by test/tier4_stochastic.jl from the same runs.

using HoyngDynamo
using Printf
import Dates
Dates_now() = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM")

nseeds = length(ARGS) ≥ 1 ? parse(Int, ARGS[1]) : 3
scale = length(ARGS) ≥ 2 ? parse(Float64, ARGS[2]) : 1.0

# Table II, with the D₁/ℱ fit c₀ + c₂a² that the printed (c², γ) imply through Eqs. 62–63
const PAPER = Dict(4.4 => (c2 = 0.318, γ = 9.04, var56 = 0.04, var = 0.03, Tr58 = 12000.0, Tr = 10700.0, err = 5300.0, c0 = 1.22, c2fit = 0.31),
                   5.4 => (c2 = 0.299, γ = 5.15, var56 = 0.08, var = 0.07, Tr58 = 480.0, Tr = 519.0, err = 53.0, c0 = 1.34, c2fit = 0.94),
                   6.4 => (c2 = 0.246, γ = 2.96, var56 = 0.17, var = 0.16, Tr58 = 94.0, Tr = 107.0, err = 5.0, c0 = 1.36, c2fit = 2.00),
                   7.4 => (c2 = 0.240, γ = 1.82, var56 = 0.47, var = 0.30, Tr58 = 35.0, Tr = 37.6, err = 1.0, c0 = 1.62, c2fit = 3.20),
                   8.4 => (c2 = 0.228, γ = 1.09, var56 = 1.0, var = 0.50, Tr58 = 35.0, Tr = 18.7, err = 0.4, c0 = 1.99, c2fit = 5.20))

rows = String[]
push!(rows, "| f | τ_d per seed | reversals | ⟨T_r⟩ (paper) | §2 variance (paper) | dv | D₁/ℱ weighted fit (paper) | c² (paper) | γ (paper) | Eq. 56 | Eq. 58 / Eq. 57 | KS p (n > 5 τ_d) | s/seed |")
push!(rows, "|---|---|---|---|---|---|---|---|---|---|---|---|---|")
println("threads: $(Threads.nthreads()), seeds: $nseeds")
for f in (4.4, 5.4, 6.4, 7.4, 8.4)
    τ = (f == 4.4 ? 2e5 : 5e4) * scale
    e = ensemble(RunSpec(; f), τ; seeds = 1:nseeds)
    s, pp = e.summary, PAPER[f]
    ks = length(e.intervals) ≥ 2 ? ks_exponential(e.intervals; tmin = 5.0) : (p = NaN, n = 0)
    row = @sprintf("| %.1f | %.0e | %d | %.1f ± %.1f (%g ± %g) | %.3f (%.2f) | %.3f | %.2f + %.2f a² (%.2f + %.2f a²) | %.3f (%.3f) | %.2f (%.2f) | %.3f | %.1f / %.1f | %.3f (%d) | %.0f |",
                   f, τ, s.n_reversals, s.mean_reversal_time, s.mean_reversal_time_error, pp.Tr, pp.err,
                   s.variance, pp.var, s.dv, s.d1_fit.c0, s.d1_fit.c2, pp.c0, pp.c2fit, s.c2, pp.c2, s.γ, pp.γ,
                   s.theory.variance, s.theory.mean_reversal_time, s.theory.mean_reversal_time_exact, ks.p, ks.n,
                   sum(r.seconds for r in e.runs) / nseeds)
    push!(rows, row)
    println(row)
    println(@sprintf("      per seed ⟨T_r⟩: %s; unweighted fit %.2f + %.2f a² → c² %.3f, γ %.2f; D₁(0) direct %.2f",
                     join([@sprintf("%.1f (%d)", r.summary.mean_reversal_time, r.summary.n_reversals) for r in e.runs], ", "),
                     s.d1_fit_unweighted.c0, s.d1_fit_unweighted.c2, s.c2_unweighted, s.γ_unweighted, s.d1_zero))
end
path = normpath(joinpath(@__DIR__, "..", "results", "table2.md"))
mkpath(dirname(path))
open(path, "w") do io
    println(io, "# HOS Table II from this implementation ($(nseeds) seeds, $(Dates_now()))\n")
    foreach(r -> println(io, r), rows)
    println(io, "\nD₁/ℱ fitted for |a| ≤ 1.1, weighted by the bin counts (HOS Fig. 8 caption), the paper's fit being the one its (c², γ) imply; ",
                "the §2 variance is the between-reversal estimator with the ±5 τ_d exclusion, dv is HSO Eq. 12 on |a| at 1 τ_d; Eq. 56–58 are evaluated at our (c², γ), Eq. 57 is the exact escape integral. ",
                "The KS test is on the pooled polarity intervals above 5 τ_d against an exponential with fitted mean (parametric bootstrap p-value).")
end
println("wrote $path")
