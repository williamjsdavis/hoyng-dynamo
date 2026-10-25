# Tier 6: HSO's analysis of the Sint-800 record (data/sint800.tsv) reproduced, gates T6.1–T6.5
# of dev/VALIDATION.md: dv with Monte-Carlo errors (HSO Eq. 12), the χ² fits of Eq. 10 and A.1
# to the VADM histogram (Table 1), the autocorrelation (Fig. 7) and the τ_d measurement by
# overlaying the autocorrelation of artificial Sint-800 series drawn from the two HSO runs.
#
#     julia --project scripts/sint800_comparison.jl [reference run file] [second run file]
#
# The run files default to results/standard_f5.4_seed1.h5 (run_standard.jl 1 5e4 5.4) and
# results/run_kR0.3_Cm120_f7.0_seed1.h5 (run_standard.jl 1 1e4 7 0.3 120); a missing file is
# regenerated in stream mode. Appends to results/sint800.md (env SINT_OUT) and writes the
# curves for viz/hso_fig01.py, viz/sint800_fits.py, viz/hso_fig07.py and viz/sint800_tau_d.py
# to results/sint800.h5 (env SINT_H5).

using HoyngDynamo, HDF5, Random

ref_file = length(ARGS) ≥ 1 ? ARGS[1] : normpath(joinpath(@__DIR__, "..", "results", "standard_f5.4_seed1.h5"))
second_file = length(ARGS) ≥ 2 ? ARGS[2] : normpath(joinpath(@__DIR__, "..", "results", "run_kR0.3_Cm120_f7.0_seed1.h5"))
out = get(ENV, "SINT_OUT", normpath(joinpath(@__DIR__, "..", "results", "sint800.md")))
h5path = get(ENV, "SINT_H5", normpath(joinpath(@__DIR__, "..", "results", "sint800.h5")))
mkpath(dirname(out))

const Λ1, μ1 = 0.527, 0.763               # the standard model (HSO case 1)
const Λ2, μ2 = 0.815, 1.82 * 0.815         # the second model, kR 0.3, C_m 120 (HSO Table 2)

d = read_sint800()
io = open(out, "a")
emit(s) = (println(io, s); println(s))
fmt(v; d = 3) = isnan(v) ? "—" : string(round(v; digits = d))

emit("")
emit("## Sint-800 comparison (Tier 6), $(d.n) points")
emit("")

# ---- T6.1 dv
dv = dipole_variation(d.vadm)
mc = montecarlo(dipole_variation, d.vadm, d.vadm_err; n = 5000, rng = Xoshiro(1))
emit("**T6.1** dv = $(fmt(dv; d = 4)) (HSO 0.071 ± 0.003); Monte-Carlo spread over 5000 replicas s_i + e_i N(0,1): $(fmt(mc.sigma; d = 4)) (replica mean $(fmt(mc.mean; d = 4)): the added errors add variance, so the spread, not the mean, is the error). Mean VADM $(fmt(sum(d.vadm) / d.n; d = 3))×10²² A m², median relative error $(fmt(sum(d.vadm_err ./ d.vadm) / d.n; d = 3)).")
emit("")

# ---- T6.2, T6.5 the fits
emit("**T6.2, T6.5** χ² fits of p(a) to the VADM histogram (HSO §4.2: a = (A/A_m) a_m, N fixed by the total count, χ² = Σ(N P_j − n_j)²/n_j, σ = |χ² − ν|/√(2ν)).")
emit("")
emit("| fit | bins | γ | c² | A_m (10²² A m²) | χ² | ν | σ | theory dv | theory ⟨T_r⟩ (τ_d) | paper |")
emit("|---|---|---|---|---|---|---|---|---|---|---|")
rows = Any[]
function row!(label, fit, binlabel, paper)
    push!(rows, (label, fit, binlabel))
    emit("| $label | $binlabel | $(fmt(fit.γ; d = 2)) | $(fmt(fit.c2)) | $(fmt(fit.Am; d = 2)) | $(fmt(fit.chi2; d = 1)) | $(fit.ν) | $(fmt(fit.sigma; d = 2)) | $(fmt(fit.theory.dv)) | $(fmt(fit.theory.mean_reversal_time; d = 0)) | $paper |")
end
for nb in (100, 50)
    bl = "$nb equal-population"
    row!("Eq. 10, (γ, c², A_m) free (case 1)", vadm_fit(d.vadm, Λ1; nbins = nb), bl, "6.3, 0.32, 5.98 at 1σ; dv 0.076, ⟨T_r⟩ 1090")
    row!("Eq. 10 at HSO's (6.3, 0.32), A_m free", vadm_fit(d.vadm, Λ1; nbins = nb, fixed = (6.3, 0.32)), bl, "A_m 5.98")
    row!("Eq. 10 at the run's (8, 0.50) (case 2)", vadm_fit(d.vadm, Λ1; nbins = nb, fixed = (8.0, 0.50)), bl, "1.3σ; dv 0.080, ⟨T_r⟩ 928")
    row!("Eq. 10 at (7, 0.45) (case 3)", vadm_fit(d.vadm, Λ1; nbins = nb, fixed = (7.0, 0.45)), bl, "1.1σ; dv 0.079, ⟨T_r⟩ 665")
    row!("Eq. A.1 (μ/Λ 1.45), free", vadm_fit(d.vadm, Λ1, μ1; nbins = nb), bl, "—")
    row!("Eq. A.1 at (4.8, 0.22) (case 4)", vadm_fit(d.vadm, Λ1, μ1; nbins = nb, fixed = (4.8, 0.22)), bl, "1.4σ; dv 0.071, ⟨T_r⟩ 1097")
    row!("Eq. A.1 at (5, 0.50), Λ 0.815, μ/Λ 1.82 (case 5)", vadm_fit(d.vadm, Λ2, μ2; nbins = nb, fixed = (5.0, 0.50)), bl, "5σ; dv 0.091, ⟨T_r⟩ 138")
end
row!("Eq. 43 free (HOS Fig. 10)", vadm_fit(d.vadm, Λ1; nbins = 35, equal_population = false), "35 equal-width from 0", "—")
row!("Eq. 43 at HOS's (5.15, 0.25), A_m free", vadm_fit(d.vadm, Λ1; nbins = 35, equal_population = false, fixed = (5.15, 0.25)), "35 equal-width from 0", "1σ")
row!("Eq. 43 at HOS's (5.15, 0.25), A_m free", vadm_fit(d.vadm, Λ1; nbins = 35, equal_population = false, lo = minimum(d.vadm), fixed = (5.15, 0.25)), "35 equal-width over the data range", "1σ")
row!("Eq. 43 at HOS's (5.15, 0.25), A_m free", vadm_fit(d.vadm, Λ1; nbins = 100, fixed = (5.15, 0.25)), "100 equal-population", "1σ")
row!("Eq. 4 free (SOH Fig. 4)", vadm_fit(d.vadm, Λ1; nbins = 30, equal_population = false), "30 equal-width from 0", "—")
row!("Eq. 4 at SOH's (5.2, 0.26), A_m free", vadm_fit(d.vadm, Λ1; nbins = 30, equal_population = false, fixed = (5.2, 0.26)), "30 equal-width from 0", "within 1σ, a_max 5.9")
emit("")

# ---- T6.3 autocorrelation
signed = ifelse.(d.age .> 780, -d.vadm, d.vadm)
ac_u = hso_autocorrelation(reverse(d.vadm), 50)
ac_s = hso_autocorrelation(reverse(signed), 50)
mca = montecarlo(x -> hso_autocorrelation(reverse(x), 50), d.vadm, d.vadm_err; n = 2000, rng = Xoshiro(2))
zc(ac) = (k = findfirst(<(0), ac) - 1; "$(k - 1)–$k kyr")
emit("**T6.3** autocorrelation (HSO §5.1, Σ_{i=1}^{749} x_i x_{i+k}, time reversed, mean subtracted): unsigned record: zero crossing at $(zc(ac_u)), minimum $(fmt(minimum(ac_u))) at $(argmin(ac_u) - 1) kyr; with the sign flipped beyond 780 ka (HSO Fig. 1): zero crossing at $(zc(ac_s)), minimum $(fmt(minimum(ac_s))) at $(argmin(ac_s) - 1) kyr (HSO: ≈ 27 kyr and ≈ −0.2 near 40 kyr). Monte-Carlo error bars $(fmt(mca.sigma[2]))–$(fmt(mca.sigma[end])) over lags 1–50 (HSO: of order 799^(−1/2) ≈ 0.035). Small-lag decay: e-folding over lags 1–8 of $(fmt(-1 / (sum((1:8) .* log.(ac_u[2:9])) / sum((1:8) .^ 2)); d = 1)) kyr.")
emit("")

# ---- T6.4 τ_d by overlap
function amplitude_stream(file, spec, τ)
    if isfile(file)
        f = open_run(file)
        a, every, dt = read_stream(f, "a")
        close(f)
        return a, dt, "the a stream of $(basename(file)) ($(length(a) * dt) τ_d)"
    end
    sim = Simulation(spec)
    nc = redraw_interval(sim.noise)
    probe = Probe(amplitude, Series(; sizehint = steps(sim, τ) ÷ nc + 1), nc)
    run!(sim, steps(sim, τ), (probe,))
    return value(probe).values, nc * timestep(sim), "a fresh stream-mode run ($τ τ_d, a every n_c steps)"
end
a1, dt1, src1 = amplitude_stream(ref_file, RunSpec(; f = 5.4, seed = 1), 5e4)
a2, dt2, src2 = amplitude_stream(second_file, RunSpec(; kR = 0.3, Cm = 120.0, f = 7.0, seed = 1), 1e4)
emit("**T6.4** τ_d by overlaying the autocorrelation of artificial Sint-800 series (rectified, running average τ_av, 799 points 1 kyr apart, 200 windows per τ_d; misfit = rms difference from the unsigned data curve over lags 1–25 kyr).")
emit("")
emit("| run | source | τ_av (kyr) | best τ_d per window: median (quartiles) | τ_d of least mean misfit | HSO |")
emit("|---|---|---|---|---|---|")
scans = Dict{String,Any}()
for (label, a, dt, src, paper, grid) in (("reference run f = 5.4", a1, dt1, src1, "11 ± 2", 4:0.5:24),
                                          ("second run kR 0.3, C_m 120, f 7", a2, dt2, src2, "16 ± 4", 4:0.5:30))
    for τav in (5.0, 0.0, 10.0)
        sc = tau_d_scan(a, dt, ac_u; τd = grid, τav, nwindows = 200)
        τav == 5.0 && (scans[label] = sc)
        emit("| $label | $src | $τav | $(sc.median) ($(sc.q25)–$(sc.q75)) | $(sc.best_mean) | $paper |")
    end
end
close(io)

# ---- curves for viz
h5open(h5path, "w") do h
    g = create_group(h, "data")
    g["age"] = d.age; g["vadm"] = d.vadm; g["vadm_err"] = d.vadm_err; g["signed"] = signed; g["ncores"] = d.ncores
    attrs(g)["dv"] = dv; attrs(g)["dv_sigma"] = mc.sigma
    g = create_group(h, "autocorrelation")
    g["lags"] = collect(0:50); g["unsigned"] = ac_u; g["signed"] = ac_s; g["sigma"] = mca.sigma
    g = create_group(h, "fits")
    for (i, (label, fit, binlabel)) in enumerate(rows)
        gg = create_group(g, "fit$i")
        attrs(gg)["label"] = label; attrs(gg)["bins"] = binlabel
        attrs(gg)["gamma"] = fit.γ; attrs(gg)["c2"] = fit.c2; attrs(gg)["Am"] = fit.Am; attrs(gg)["sigma"] = fit.sigma; attrs(gg)["chi2"] = fit.chi2; attrs(gg)["nu"] = fit.ν
        attrs(gg)["kappa"] = fit.κ; attrs(gg)["N"] = fit.N; attrs(gg)["am"] = fit.am
        gg["lo"] = fit.bins.lo; gg["hi"] = fit.bins.hi; gg["n"] = fit.bins.n; gg["expected"] = fit.expected
    end
    g = create_group(h, "tau_d")
    for (label, sc) in scans
        gg = create_group(g, label)
        gg["tau_d"] = sc.τd; gg["mean_misfit"] = sc.mean_misfit; gg["best"] = sc.best; gg["acor"] = sc.acor; gg["lags"] = collect(sc.lags)
        attrs(gg)["median"] = sc.median; attrs(gg)["q25"] = sc.q25; attrs(gg)["q75"] = sc.q75; attrs(gg)["best_mean"] = sc.best_mean
    end
end
println("wrote $out and $h5path")
