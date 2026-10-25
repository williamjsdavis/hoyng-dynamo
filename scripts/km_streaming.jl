# Drift and diffusion of the dipole amplitude from streamed runs by online kernel-based
# regression (Davis 2023, arXiv:2307.00445; `ConditionalMoments`, `kramers_moyal`), against
# the Fokker–Planck reduction of HOS §4: S(a) = Λ(1 − a²)a[1 + κa²] (Eq. 33, HSO App. A) and
# D(a) = Δ₀a² + Δ₁ = 2D⁽²⁾(a) (Eq. 39) with (Δ₀, Δ₁) from the run's own D₁ measurement
# (Eqs. 61–63) and from HOS's printed (c², γ).
#
#     julia --project -t 4 scripts/km_streaming.jl [f = 6.4] [τ per seed = 5e4] [seeds = 3]
#
# Appends to results/km_streaming.md (env KM_OUT) and writes the conditional moments and the
# estimates to results/km_f<f>.h5 (env KM_H5) for viz/km_drift_diffusion.py. Two probes: a
# every n_c steps (Δt = τ_c) with lags to 2 τ_d, and a every step with lags to 10 steps for the
# sub-τ_c regime; 61 points on [−1.5, 1.5], bandwidth 0.1, Epanechnikov kernel.

using HoyngDynamo, HDF5, LinearAlgebra
BLAS.set_num_threads(1)

f = length(ARGS) ≥ 1 ? parse(Float64, ARGS[1]) : 6.4
τ = length(ARGS) ≥ 2 ? parse(Float64, ARGS[2]) : 5e4
nseeds = length(ARGS) ≥ 3 ? parse(Int, ARGS[3]) : 3
out = get(ENV, "KM_OUT", normpath(joinpath(@__DIR__, "..", "results", "km_streaming.md")))
h5path = get(ENV, "KM_H5", normpath(joinpath(@__DIR__, "..", "results", "km_f$(f).h5")))
mkpath(dirname(out))

spec = RunSpec(; f)
x = range(-1.5, 1.5; length = 61)
observers(sim) = merge(standard_observers(sim), drift_diffusion_observers(sim; x, h = 0.1, lags = 40),
                       (km_fine = Probe(amplitude, ConditionalMoments(x; h = 0.1, lags = 10, Δt = timestep(sim)), 1),))
t = @elapsed e = ensemble(spec, τ; seeds = 1:nseeds, observers)
s = e.summary
c = e.runs[1].sim.model.constants
Λ, μ = c.Λ, c.μ
km, kmf = value(e.pooled.km), value(e.pooled.km_fine)

# the reduced description with the run's own (c², γ) (weighted D₁ fit) and with HOS's printed values
th_run = BistableTheory(c, s.γ, s.c2)
th_run0 = BistableTheory(s.γ, s.c2; Λ)
paper = Dict(6.4 => (2.96, 0.246), 5.4 => (5.15, 0.299), 7.4 => (1.82, 0.240), 4.4 => (9.04, 0.318), 8.4 => (1.09, 0.228))
th_paper = haskey(paper, f) ? BistableTheory(c, paper[f]...) : nothing

io = open(out, "a")
emit(s) = (println(io, s); println(s))
fmt(v; d = 3) = isnan(v) ? "—" : string(round(v; digits = d))

emit("")
emit("## Drift and diffusion by online kernel-based regression: f = $f, $nseeds seeds × $τ τ_d ($(round(t, digits = 0)) s on $(Threads.nthreads()) threads)")
emit("")
emit("Run: ⟨T_r⟩ = $(fmt(s.mean_reversal_time; d = 1)) ± $(fmt(s.mean_reversal_time_error; d = 1)), dv = $(fmt(s.dv)); D₁/ℱ fit $(fmt(s.d1_fit.c0)) + $(fmt(s.d1_fit.c2)) a² → (c², γ) = ($(fmt(s.c2)), $(fmt(s.γ; d = 2))), Δ₀ = $(fmt(s.Δ0; d = 4)), Δ₁ = $(fmt(s.Δ1; d = 4)); Λ = $(fmt(Λ)), μ/Λ = $(fmt(μ / Λ; d = 2)).")
emit("Conditional moments of a: $(length(x)) points on [−1.5, 1.5], h = 0.1, Epanechnikov; coarse probe every n_c steps (Δt = $(km.lags[1])) with $(length(km.lags)) lags, fine probe every step with $(length(kmf.lags)) lags; $(km.n) and $(kmf.n) samples pooled.")
emit("")

pts = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25]
nearest(xs, p) = argmin(abs.(xs .- p))
jdx = [nearest(km.x, p) for p in pts]
lagsel = [(0.05, :fine, 5), (0.1, :fine, 10), (0.25, :coarse, 5), (0.5, :coarse, 10), (1.0, :coarse, 20), (2.0, :coarse, 40)]
S(th, a) = drift(th, a)
D(th, a) = diffusion(th, a)          # HOS's D = Δ₀(a² + c²), the coefficient of ½∂²(Dp)/∂a², i.e. 2D⁽²⁾

emit("### D⁽¹⁾(a; τ) = M⁽¹⁾/τ (direct estimation, Eq. 16) against the drift S(a)")
emit("")
emit("| a | " * join(["τ = $(l[1])" for l in lagsel], " | ") * " | S(a), run (c², γ), κ = μ/Λ − 1 | S(a), κ = 0 |")
emit("|---|" * repeat("---|", length(lagsel) + 2))
for (p, j) in zip(pts, jdx)
    cells = String[]
    for (τl, which, i) in lagsel
        d = kramers_moyal(which === :fine ? kmf : km; lags = i:i, method = :direct)
        push!(cells, "$(fmt(d.D1[j])) ± $(fmt(d.D1_err[j]))")
    end
    emit("| $p | " * join(cells, " | ") * " | $(fmt(S(th_run, p))) | $(fmt(S(th_run0, p))) |")
end
emit("")
emit("### 2D⁽²⁾(a; τ) = M⁽²*⁾/τ (conditional variance, Eq. 16 with the drift term removed) against D(a) = Δ₀a² + Δ₁")
emit("")
emit("| a | " * join(["τ = $(l[1])" for l in lagsel], " | ") * " | D(a), run (c², γ) | D(a), HOS Table II |")
emit("|---|" * repeat("---|", length(lagsel) + 2))
for (p, j) in zip(pts, jdx)
    cells = String[]
    for (τl, which, i) in lagsel
        d = kramers_moyal(which === :fine ? kmf : km; lags = i:i, method = :direct, variance = true)
        push!(cells, "$(fmt(2d.D2[j]; d = 4)) ± $(fmt(2d.D2_err[j]; d = 4))")
    end
    emit("| $p | " * join(cells, " | ") * " | $(fmt(D(th_run, p); d = 4)) | $(th_paper === nothing ? "—" : fmt(D(th_paper, p); d = 4)) |")
end
emit("")

emit("### Regression through the origin over lag windows (Eq. 17), and the shape parameters they imply")
emit("")
emit("| window (τ_d) | Λ_KM = dD⁽¹⁾/da at 0 | zero of D⁽¹⁾ (a > 0) | Δ₁ = 2D⁽²⁾(0) | Δ₀ + Δ₁ = 2D⁽²⁾(1) | even quadratic fit of 2D⁽²⁾ for \\|a\\| ≤ 1.1 | (c², γ)_KM | run D₁ fit (c², γ) |")
emit("|---|---|---|---|---|---|---|---|")
windows = [("0.05", kmf, 5:5), ("0.1–0.5", km, 2:10), ("0.25–1.0", km, 5:20), ("0.5–2.0", km, 10:40)]
curves = Dict{String,Any}()
for (label, v, lags) in windows
    d = kramers_moyal(v; lags, method = :regression, variance = true)
    curves[label] = d
    j0, jp, jm, j1 = nearest(d.x, 0.0), nearest(d.x, 0.25), nearest(d.x, -0.25), nearest(d.x, 1.0)
    Λkm = (d.D1[jp] - d.D1[jm]) / 0.5
    # first sign change of D⁽¹⁾ beyond a = 0.5
    k = findfirst(i -> d.x[i] > 0.5 && isfinite(d.D1[i]) && isfinite(d.D1[i-1]) && d.D1[i-1] > 0 ≥ d.D1[i], 2:length(d.x))
    azero = k === nothing ? NaN : (i = k + 1; d.x[i-1] + (d.x[i] - d.x[i-1]) * d.D1[i-1] / (d.D1[i-1] - d.D1[i]))
    Δ1 = 2d.D2[j0]; D1val = 2d.D2[j1]
    sel = [i for i in eachindex(d.x) if abs(d.x[i]) ≤ 1.1 && isfinite(d.D2[i])]
    q0, q2 = even_quadratic_fit(d.x[sel], 2 .* d.D2[sel], 1 ./ max.(2 .* d.D2_err[sel], 1e-6) .^ 2)
    shape = q0 > 0 && q2 > 0 ? "($(fmt(q0 / q2)), $(fmt(Λ / q2; d = 2)))" : "—"
    emit("| $label | $(fmt(Λkm)) | $(fmt(azero)) | $(fmt(Δ1; d = 4)) | $(fmt(D1val; d = 4)) | $(fmt(q0; d = 4)) + $(fmt(q2; d = 4)) a² | $shape | ($(fmt(s.c2)), $(fmt(s.γ; d = 2))) |")
end
emit("")
emit("Λ = $(fmt(Λ)), 2μ = $(fmt(2μ)); the reduced drift has S'(0) = Λ and S'(±1) = −2μ. In the Fokker–Planck convention D = 2D⁽²⁾; HOS's (c², γ) = (Δ₁/Δ₀, Λ/Δ₀).")
emit("")

# the noise diffusion of a evaluated on the fields of the run: exact for the binned noise (D_b) and in
# HOS's δ-correlated limit (D_c = ℱ∫ sin θ P̃₀² T² dθ = D₀a² + D₁ + cross term), binned in |a| once per τ_d,
# against the Kramers–Moyal 2D⁽²⁾ corrected for the piecewise-constant noise, Var(∫₀^τ ξ) = Dτ(1 − τ_c/3τ)
db, dc = value(e.pooled.d_binned), value(e.pooled.d_continuum)
τc = e.runs[1].sim.noise.τc
emit("### The noise diffusion of a on the run's fields: exact binned D_b against HOS's continuum D_c = D₀a² + D₁ (Eqs. 36–37), and the Kramers–Moyal 2D⁽²⁾ corrected for the noise correlation")
emit("")
emit("| \\|a\\| bin | samples | D_c (continuum, Eqs. 36–37) | D_b (binned, exact) | D_b/D_c | 2D⁽²⁾/(1 − τ_c/3τ) at τ = 0.1 | at τ = 0.25 | at τ = 0.5 | D(a) from the D₁ fit |")
emit("|---|---|---|---|---|---|---|---|---|")
for p in (0.025, 0.275, 0.525, 0.775, 1.025, 1.275)
    i = nearest(db.centers, p)
    cells = String[]
    for (which, l) in ((kmf, 10), (km, 5), (km, 10))
        d = kramers_moyal(which; lags = l:l, method = :direct, variance = true)
        j = nearest(d.x, p)
        corr = 1 - τc / (3 * d.lags[1])
        push!(cells, "$(fmt(2d.D2[j] / corr; d = 4)) ± $(fmt(2d.D2_err[j] / corr; d = 4))")
    end
    emit("| $(fmt(db.centers[i] - 0.025; d = 2))–$(fmt(db.centers[i] + 0.025; d = 2)) | $(db.n[i]) | $(fmt(dc.mean[i]; d = 4)) | $(fmt(db.mean[i]; d = 4)) | $(fmt(db.mean[i] / dc.mean[i]; d = 3)) | " * join(cells, " | ") * " | $(fmt(D(th_run, db.centers[i]); d = 4)) |")
end
emit("")
emit("D_c binned in |a| is the papers' D measurement with the full T (D₁ of Fig. 8 uses T − aT₀); D_b is the same for the N_c = $(e.runs[1].sim.noise.Nc) latitude bins actually used (DR-3: 5% below D_c on the pure fundamental, more where T has finer latitudinal structure).")
close(io)

h5open(h5path, "w") do h
    attrs(h)["f"] = f; attrs(h)["tau"] = τ; attrs(h)["seeds"] = nseeds; attrs(h)["Lambda"] = Λ; attrs(h)["mu"] = μ
    attrs(h)["gamma_run"] = s.γ; attrs(h)["c2_run"] = s.c2; attrs(h)["Delta0_run"] = s.Δ0; attrs(h)["Delta1_run"] = s.Δ1
    if th_paper !== nothing
        attrs(h)["gamma_paper"] = paper[f][1]; attrs(h)["c2_paper"] = paper[f][2]
    end
    for (name, v) in (("coarse", km), ("fine", kmf))
        g = create_group(h, name)
        g["x"] = collect(v.x); g["lags"] = collect(v.lags); g["W"] = v.W; g["M1"] = v.M1; g["M2"] = v.M2; g["M2star"] = v.M2star; g["neff"] = v.neff
        attrs(g)["n"] = v.n; attrs(g)["h"] = v.h
    end
    g = create_group(h, "regression")
    for (label, d) in curves
        gg = create_group(g, label)
        gg["x"] = d.x; gg["D1"] = d.D1; gg["D2"] = d.D2; gg["D1_err"] = d.D1_err; gg["D2_err"] = d.D2_err; gg["lags"] = collect(d.lags)
    end
    for (name, v) in (("d_binned", db), ("d_continuum", dc))
        g = create_group(h, name)
        g["centers"] = collect(v.centers); g["n"] = v.n; g["mean"] = v.mean; g["var"] = v.var
    end
    attrs(h)["tau_c"] = τc
end
println("wrote $out and $h5path")
