# DR-5 sensitivity (dev/DISCREPANCIES.md): what moves the intercept of the D₁/ℱ fit in the
# standard run? The noise distribution, the time integrator and the resolution (Δτ halved,
# N doubled; the last two are gate T4.17), each as 3 seeds × 5×10⁴ τ_d pooled:
#
#     julia --project -t 4 scripts/dr5_sensitivity.jl
#
# About 2.5 min on 4 threads. Result 2026-10-20: converged in Δτ and N, independent of the
# distribution of F, and the integrator at Δτ = 0.01 moves the intercept by 0.26, the size of
# the gap to the paper's 1.36.

using HoyngDynamo, Printf

variants = (
    ("baseline: GaussianF, exponential, Δτ 0.01, N 60", RunSpec()),
    ("UniformF", RunSpec(; distribution = UniformF())),
    ("BinaryF", RunSpec(; distribution = BinaryF())),
    ("IMEX stepper", RunSpec(; stepper = :imex)),
    ("Δτ = 0.005 (same noise realisation)", RunSpec(; Δτ = 0.005)),
    ("N = 120", RunSpec(; N = 120)),
)
println("threads: ", Threads.nthreads())
@printf("%-48s | %-14s | %-6s | %-6s | %-16s | %-16s | %-6s | %-6s | %-13s | %s\n",
        "variant", "⟨T_r⟩", "§2 var", "dv", "weighted fit", "unweighted fit", "D₁(0)", "D₁(1)", "(c², γ)", "s/seed")
for (name, spec) in variants
    e = ensemble(spec, 5e4; seeds = 1:3)
    s = e.summary
    d = value(e.pooled.d1)
    i = findfirst(c -> abs(c - 0.975) < 1e-9, d.centers)
    d1_at_1 = d.mean[i]          # bin average at |a| = 0.95–1.0
    @printf("%-48s | %6.1f ± %4.1f | %.3f | %.3f | %.2f + %.2f a² | %.2f + %.2f a² | %.2f | %.2f | (%.3f, %.2f) | %.0f\n",
            name, s.mean_reversal_time, s.mean_reversal_time_error, s.variance, s.dv,
            s.d1_fit.c0, s.d1_fit.c2, s.d1_fit_unweighted.c0, s.d1_fit_unweighted.c2, s.d1_zero, d1_at_1, s.c2, s.γ,
            sum(r.seconds for r in e.runs) / 3)
    flush(stdout)
end
println("paper: ⟨T_r⟩ 107 ± 5; variance 0.16; fit 1.36 + 2.00 a² (D₁(1) on the curve: 3.36, bin 0.95–1.0: 3.26); (c², γ) = (0.246, 2.96)")
