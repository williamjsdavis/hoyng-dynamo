# SG-1 (dev/STRETCH_GOALS.md): how far, and in which direction, do plausible implicit schemes
# at the papers' Δτ = 0.01 move the D₁ fit and the reversal statistics (DISCREPANCIES.md DR-5)?
# Every stepper of `STEPPERS` on the standard run at the given forcings, three seeds × 5×10⁴ τ_d
# pooled (`ensemble`), preceded by the analytic expectation: each scheme's per-step factor on
# the leading modes, at C₀ (a = ±1) and at C_m (a = 0), against the exact exp(Δτλ).
#
#     julia --project -t 4 scripts/sg1_integrators.jl [f ...]        (default 6.4 5.4 7.4)
#
# Prints the tables and appends them to results/sg1_integrators.md (SG1_OUT overrides; SG1_TAU
# shortens the runs for a smoke test). About 10 min per forcing on 4 threads: the θ-method
# steppers cost 18 µs per step against 2 µs for the exponential one.

using HoyngDynamo, Printf, LinearAlgebra, Dates
BLAS.set_num_threads(1)          # the seeds run on Julia threads; keep LAPACK single-threaded

# HOS Table II (the fits are the ones the paper's (c², γ) imply, as in results/table2.md)
const PAPER = Dict(5.4 => (Tr = "519 ± 53", var = 0.07, fit = "1.34 + 0.94 a²", c2 = 0.299, γ = 5.15),
                   6.4 => (Tr = "107 ± 5", var = 0.16, fit = "1.36 + 2.00 a²", c2 = 0.246, γ = 2.96),
                   7.4 => (Tr = "37.6 ± 1", var = 0.30, fit = "1.62 + 3.20 a²", c2 = 0.240, γ = 1.82))

forcings = isempty(ARGS) ? [6.4, 5.4, 7.4] : parse.(Float64, ARGS)
seeds = 1:3
τ_end = parse(Float64, get(ENV, "SG1_TAU", "5e4"))
Δτ = 0.01
out = get(ENV, "SG1_OUT", normpath(joinpath(@__DIR__, "..", "results", "sg1_integrators.md")))
mkpath(dirname(out))
const OUT = open(out, "a")
emit(s::AbstractString = "") = (println(s); println(OUT, s); flush(OUT))

emit()
emit("## SG-1 integrator sensitivity ($(Dates.format(now(), "yyyy-mm-dd HH:MM")), $(Threads.nthreads()) threads, $(length(seeds)) seeds × $(τ_end) τ_d)")
emit()

# ---- analytic expectation ------------------------------------------------------------------
# The one-step linear maps for a frozen quenching coefficient c (c = 0 at a = ±1, c = s at
# a = 0) without noise, and their eigenvalue nearest exp(Δτλ_k) for the leading modes λ_k of
# A = G + cE: the rate ln|g|/Δτ against ℜλ_k and the frequency |arg g|/Δτ against ℑλ_k.
model = Model()
g, cst = model.grid, model.constants
M = g.M
G, E = model.G, quenching_operator(g)
L = diffusion_operator(g, cst.kR)
Lbd = zeros(2M, 2M)
Lbd[1:M, 1:M] .= L
Lbd[M+1:2M, M+1:2M] .= L
X = exp(Δτ .* G)
maps(c) = (exponential = X * (I + Δτ * c * E),
           imex = (I - Δτ * Lbd) \ (I + Δτ * (G + c * E - Lbd)),
           backward_euler = inv(I - Δτ * (G + c * E)),
           crank_nicolson = (I - Δτ / 2 * (G + c * E)) \ (I + Δτ / 2 * (G + c * E)))

emit("### Per-step factor on the leading modes, as (rate, frequency) per τ_d: ln|g|/Δτ against ℜλ, |arg g|/Δτ against ℑλ")
emit()
emit("| operator | mode | exact | exponential | IMEX | backward Euler | Crank–Nicolson |")
emit("|---|---|---|---|---|---|---|")
for (label, c) in (("C₀ (a = ±1)", 0.0), ("C_m (a = 0)", cst.s))
    λ = sort(eigvals(G + c * E); by = real, rev = true)
    λ = λ[imag.(λ) .≥ -1e-9]                  # one member per conjugate pair
    evs = map(eigvals, maps(c))
    for (k, name) in ((1, "fundamental"), (2, "first overtone"), (3, "second overtone"))
        λk = λ[k]
        cells = map(evs) do ev
            gk = ev[argmin(abs.(ev .- exp(Δτ * λk)))]
            @sprintf("%.4f, %.3f", log(abs(gk)) / Δτ, abs(angle(gk)) / Δτ)
        end
        emit(@sprintf("| %s | %s | %.4f, %.3f | %s | %s | %s | %s |", label, name, real(λk), imag(λk), cells...))
    end
end
emit()

# ---- runs ----------------------------------------------------------------------------------
for f in forcings
    emit("### f = $f")
    emit()
    emit("| stepper | ⟨T_r⟩ | reversals | §2 var | dv | weighted D₁ fit | unweighted fit | D₁(0) direct | D₁ bin 0.95–1.0 | (c², γ) | s/seed |")
    emit("|---|---|---|---|---|---|---|---|---|---|---|")
    rows = Dict{Symbol,Any}()
    for stepper in HoyngDynamo.STEPPERS
        e = ensemble(RunSpec(; f, stepper), τ_end; seeds)
        s = e.summary
        d = value(e.pooled.d1)
        i = findfirst(c -> abs(c - 0.975) < 1e-9, d.centers)
        d1_at_1 = d.mean[i]
        rows[stepper] = (s = s, d1_at_1 = d1_at_1)
        emit(@sprintf("| %s | %.1f ± %.1f | %d | %.3f | %.3f | %.2f + %.2f a² | %.2f + %.2f a² | %.2f | %.2f | (%.3f, %.2f) | %.0f |",
                      stepper, s.mean_reversal_time, s.mean_reversal_time_error, s.n_reversals, s.variance, s.dv,
                      s.d1_fit.c0, s.d1_fit.c2, s.d1_fit_unweighted.c0, s.d1_fit_unweighted.c2, s.d1_zero, d1_at_1, s.c2, s.γ,
                      sum(r.seconds for r in e.runs) / length(seeds)))
    end
    p = get(PAPER, f, nothing)
    p === nothing || emit(@sprintf("| paper | %s | | %.2f | | %s | | | | (%.3f, %.2f) | |", p.Tr, p.var, p.fit, p.c2, p.γ))
    emit()
    ref = rows[:exponential]
    emit("Shift from the exponential stepper at f = $f:")
    emit()
    emit("| stepper | Δ weighted intercept | Δ c² | Δ⟨T_r⟩ (%) | Δ dv | Δ D₁ bin 0.95–1.0 |")
    emit("|---|---|---|---|---|---|")
    for stepper in HoyngDynamo.STEPPERS
        stepper === :exponential && continue
        r = rows[stepper]
        emit(@sprintf("| %s | %+.2f | %+.3f | %+.1f | %+.3f | %+.2f |", stepper, r.s.d1_fit.c0 - ref.s.d1_fit.c0, r.s.c2 - ref.s.c2,
                      100 * (r.s.mean_reversal_time / ref.s.mean_reversal_time - 1), r.s.dv - ref.s.dv, r.d1_at_1 - ref.d1_at_1))
    end
    emit()
end
close(OUT)
println("appended to $out")
