# Tier 4: stochastic runs, quantitative (slow). Gates T4.1–T4.18 in dev/VALIDATION.md.
# Build step 6 covers the standard run (T4.1–T4.4), build step 7 the forcing sweep (T4.5–T4.7,
# T4.11, T4.12, T4.16, T4.17, T4.18), the HSO runs (T4.8, T4.9), the SOH runs (T4.10, T4.15) and the
# event-level diagnostics (T4.13, T4.14). Definitions of
# VALIDATION.md §2: Schmitt trigger a_h = 0.5, a_e = 0.3; exclusion variance W = 5 τ_d; dv on
# |a| at 1 τ_d; D₁/ℱ bin averages of width 0.05 fitted for |a| ≤ 1.1 weighted by the bin
# counts (HOS Fig. 8 caption). Seeds run in parallel over the Julia threads
# (JULIA_NUM_THREADS=4 julia --project -e 'using Pkg; Pkg.test()').

const T4_LENGTH = 5e4                 # τ_d per seed, as in the papers
const T4_SEEDS = 1:3

function report(label, e; paper = nothing)
    s = e.summary
    per = join(["$(round(r.summary.mean_reversal_time, digits = 1)) ($(r.summary.n_reversals))" for r in e.runs], ", ")
    ks = length(e.intervals) ≥ 2 ? ks_exponential(e.intervals; tmin = 5.0) : (D = NaN, p = NaN, n = 0, mean = NaN)
    @info "$label: ⟨T_r⟩ = $(round(s.mean_reversal_time, digits = 1)) ± $(round(s.mean_reversal_time_error, digits = 1)) from $(s.n_reversals) reversals" *
          (paper === nothing ? "" : " (paper $(paper.Tr) ± $(paper.err))") * "; per seed $per; $(s.n_aborted) aborted, $(s.n_excursions) excursions; " *
          "grouping rule $(s.n_reversals_grouped); variance $(round(s.variance, digits = 3))" * (paper === nothing ? "" : " (paper $(paper.var))") *
          " over $(s.variance_intervals) intervals; dv $(round(s.dv, digits = 3)); " *
          "D₁/ℱ = $(round(s.d1_fit.c0, digits = 3)) + $(round(s.d1_fit.c2, digits = 3)) a²" * (paper === nothing ? "" : " (paper $(paper.fit.c0) + $(paper.fit.c2) a²)") *
          " (unweighted $(round(s.d1_fit_unweighted.c0, digits = 3)) + $(round(s.d1_fit_unweighted.c2, digits = 3)) a², " *
          "D₁(0) direct $(round(s.d1_zero, digits = 3))); c² = $(round(s.c2, digits = 3)), γ = $(round(s.γ, digits = 3))" *
          (paper === nothing ? "" : " (paper $(paper.c2), $(paper.γ))") * "; unweighted c² = $(round(s.c2_unweighted, digits = 3)), γ = $(round(s.γ_unweighted, digits = 3)); " *
          "theory at the measured (c², γ): variance $(round(s.theory.variance, digits = 3)), ⟨T_r⟩ Eq. 58 $(round(s.theory.mean_reversal_time, digits = 1)), Eq. 57 $(round(s.theory.mean_reversal_time_exact, digits = 1)); " *
          "KS (intervals > 5 τ_d): D = $(round(ks.D, digits = 4)), p = $(round(ks.p, digits = 3)), n = $(ks.n); " *
          "$(round(sum(r.seconds for r in e.runs), digits = 1)) s on $(Threads.nthreads()) threads"
    return ks
end

# HOS Table II (p. 292): f, c², γ, variance from run, ⟨T_r⟩ from run ± error, and the D₁/ℱ fit
# c₀ + c₂a² the printed (c², γ) imply through Eqs. 62–63 (c₀ = c²Λ/(γℱ), c₂ = Λ/(γℱ) − 3.529,
# Λ = 0.527); for f = 6.4 that is the printed 1.36 + 2.00a² to the rounding.
const TABLE_II = Dict(4.4 => (c2 = 0.318, γ = 9.04, var = 0.03, Tr = 10700.0, err = 5300.0, fit = (c0 = 1.22, c2 = 0.31)),
                      5.4 => (c2 = 0.299, γ = 5.15, var = 0.07, Tr = 519.0, err = 53.0, fit = (c0 = 1.34, c2 = 0.94)),
                      6.4 => (c2 = 0.246, γ = 2.96, var = 0.16, Tr = 107.0, err = 5.0, fit = (c0 = 1.36, c2 = 2.00)),
                      7.4 => (c2 = 0.240, γ = 1.82, var = 0.30, Tr = 37.6, err = 1.0, fit = (c0 = 1.62, c2 = 3.20)),
                      8.4 => (c2 = 0.228, γ = 1.09, var = 0.50, Tr = 18.7, err = 0.4, fit = (c0 = 1.99, c2 = 5.20)))

@testset "Tier 4: stochastic runs, quantitative (slow)" begin
    standard = ensemble(RunSpec(), T4_LENGTH; seeds = T4_SEEDS)      # shared by the first and third testsets

    @testset "T4.1–T4.4 standard run (kR 0.5, C_m 100, f 6.4): $(length(T4_SEEDS)) seeds × $(T4_LENGTH) τ_d" begin
        e = standard
        s, pp = e.summary, TABLE_II[6.4]
        ks = report("T4.1–T4.4 standard run", e; paper = pp)
        # T4.1: mean time between reversals, 2σ combined with the paper's ±5
        @test abs(s.mean_reversal_time - pp.Tr) < 2 * sqrt(pp.err^2 + s.mean_reversal_time_error^2)
        # T4.2: between-reversal relative variance (§2 exclusion estimator) and dv (HSO Eq. 12)
        @test abs(s.variance - pp.var) < 0.02
        @test abs(s.dv - pp.var) < 0.02
        # T4.3: the D₁/ℱ fit for |a| ≤ 1.1 with the paper's weighting. The intercept comes out
        # 16% above the paper's 1.36 (the unweighted fit gives 1.46, the direct small-|a| mean
        # 1.24): a known failure of the 10% gate, DISCREPANCIES.md DR-5.
        @test_broken relerr(s.d1_fit.c0, 1.36) < 0.1
        @test relerr(s.d1_fit.c2, 2.00) < 0.1
        @test relerr(s.d1_fit_unweighted.c0, 1.36) < 0.1               # the unweighted fit, for the record
        # T4.4: derived shape parameters; c² inherits the intercept's excess (0.29 against 0.246)
        @test_broken relerr(s.c2, pp.c2) < 0.1
        @test relerr(s.γ, pp.γ) < 0.1
        @test relerr(s.c2_unweighted, pp.c2) < 0.1
        # T4.11: polarity intervals above 5 τ_d are exponential (KS, bootstrap p-value)
        @test ks.p > 0.05
        @test all(r.summary.n_reversals ≥ 300 for r in e.runs)        # each seed alone has decent statistics
        # T4.13 (diagnostic, kind D): the parity histogram over the whole run, sampled every n_c
        # steps, events included; the event-window exclusion is a later piece of step 7
        pv = value(e.pooled.parity)
        cum = cumsum(pv.counts)
        median_P = pv.centers[findfirst(≥(cum[end] / 2), cum)]
        @info "T4.13 parity over the whole run (events included): median 𝒫 ≈ $(round(median_P, digits = 3)), fraction below −0.9: $(round(sum(pv.counts[pv.centers .< -0.9]) / cum[end], digits = 3))"
        @test -1 ≤ median_P < 0
    end

    sweep = Dict{Float64,Any}()                                     # the f = 5.4 ensemble is HSO's reference run

    @testset "T4.5–T4.7, T4.11, T4.18 forcing sweep (HOS Table II)" begin
        for f in (4.4, 5.4, 7.4, 8.4)
            τ = f == 4.4 ? 4 * T4_LENGTH : T4_LENGTH           # the paper had ~5 reversals at f = 4.4
            e = ensemble(RunSpec(; f), τ; seeds = T4_SEEDS)
            sweep[f] = e
            s, pp = e.summary, TABLE_II[f]
            ks = report("T4.5–T4.7 f = $f ($(length(T4_SEEDS)) × $τ τ_d)", e; paper = pp)
            # T4.5: ⟨T_r⟩ within 2σ combined
            @test abs(s.mean_reversal_time - pp.Tr) < 2 * sqrt(pp.err^2 + s.mean_reversal_time_error^2)
            # T4.6: the paper's "variance from run" within ±0.02 or 15%. dv (HSO Eq. 12 on |a|, no
            # exclusion) reproduces the column to 1% at every forcing; the §2 exclusion estimator
            # (pooled moments) is 10–13% below it at f ≥ 7.4, where the intervals are shorter than
            # 2W and the kept samples are the quiet minority (DISCREPANCIES.md DR-6)
            @test abs(s.dv - pp.var) < max(0.02, 0.15 * pp.var)
            @test abs(s.variance - pp.var) < max(0.02, 0.15 * pp.var)
            # T4.7: (c², γ) from the D₁ fit within 10%. γ passes at every forcing; c² inherits the
            # fit's intercept, which the paper's weighting reproduces at f ≥ 7.4 (broad |a|
            # distributions) and overestimates at f ≤ 6.4 (DR-5); at f = 4.4 fewer than 0.5% of the
            # samples have |a| < 0.5 and the intercept is unconstrained by either weighting
            @test relerr(s.γ, pp.γ) < 0.1
            if f ≥ 7.4
                @test relerr(s.c2, pp.c2) < 0.1
            else
                @test_broken relerr(s.c2, pp.c2) < 0.1
            end
            if 5.4 ≤ f ≤ 6.4
                @test relerr(s.c2_unweighted, pp.c2) < 0.1
            end
            # T4.11: Poisson intervals wherever there are enough of them
            ks.n ≥ 50 && @test ks.p > 0.05
            if f == 5.4
                # T4.18: HOS run 2 (519 ± 53) and the nominally identical HSO reference run (570 ± 61)
                Trs = [r.summary.mean_reversal_time for r in e.runs]
                @info "T4.18 seed scatter at f = 5.4: ⟨T_r⟩ = $(round.(Trs, digits = 1)); HOS 519 ± 53, HSO 570 ± 61"
                @test minimum(Trs) - 60 < 519 && maximum(Trs) + 60 > 570
            end
        end
    end

    @testset "T4.8, T4.9 the HSO runs; T5.13 autocorrelation rate (diagnostic)" begin
        # T4.8: HSO's reference run is HOS run 2 (f = 5.4): ⟨T_r⟩ 570 ± 61 from 88 reversals,
        # dv 0.070, and (γ, c²) = (8, 0.50) from the χ² fit of Eq. 10 to the |a| histogram
        # (HSO Table 1 case 2, §4.3); the A.1 fit is reported (HSO's case 4 fits both records)
        e = sweep[5.4]
        s = e.summary
        hf = s.histogram_fit
        @info "T4.8 HSO reference run (f = 5.4, pooled): ⟨T_r⟩ = $(round(s.mean_reversal_time, digits = 1)) ± $(round(s.mean_reversal_time_error, digits = 1)) (HSO 570 ± 61); dv $(round(s.dv, digits = 4)) (0.070); " *
              "Eq. 10 fit (γ, c²) = ($(round(hf.eq10.γ, digits = 2)), $(round(hf.eq10.c2, digits = 3))) at $(round(hf.eq10.sigma, digits = 1))σ over $(hf.eq10.nbins) bins (HSO (8, 0.50)); " *
              "Eq. A.1 fit ($(round(hf.a1.γ, digits = 2)), $(round(hf.a1.c2, digits = 3))) at $(round(hf.a1.sigma, digits = 1))σ (HSO case 4, fitted to both records: (4.8, 0.22)); " *
              "per seed Eq. 10: $(join(["($(round(r.summary.histogram_fit.eq10.γ, digits = 2)), $(round(r.summary.histogram_fit.eq10.c2, digits = 3)))" for r in e.runs], ", "))"
        @test abs(s.mean_reversal_time - 570) < 2 * sqrt(61^2 + s.mean_reversal_time_error^2)
        @test abs(s.dv - 0.070) < 0.005
        @test relerr(hf.eq10.γ, 8.0) < 0.1
        @test relerr(hf.eq10.c2, 0.50) < 0.1
        # T5.13 (diagnostic): the small-lag decay rate of the autocorrelation of |a| against
        # 2Λ (HSO Eq. 15) and 2μ (the deterministic relaxation, T1.4), from the whole-run
        # AutoCorrelation probe (lags of n_c steps to 5 τ_d), fitted over 0.1–0.6 τ_d
        ac = value(e.pooled.autocorrelation)
        c = e.runs[1].sim.model.constants
        sel = (ac.lags .≥ 0.1) .& (ac.lags .≤ 0.6) .& (ac.acor .> 0)
        rate = -(sum(ac.lags[sel] .* log.(ac.acor[sel])) - sum(ac.lags[sel]) * sum(log.(ac.acor[sel])) / count(sel)) /
                (sum(ac.lags[sel] .^ 2) - sum(ac.lags[sel])^2 / count(sel))
        first_neg = findfirst(<(0), ac.acor)
        @info "T5.13 autocorrelation of |a| at f = 5.4: decay rate over 0.1–0.6 τ_d = $(round(rate, digits = 3)) per τ_d against 2Λ = $(round(2c.Λ, digits = 3)) and 2μ = $(round(2c.μ, digits = 3)); " *
              "acor(0.5 τ_d) = $(round(ac.acor[findfirst(≥(0.5), ac.lags)], digits = 3)), acor(5 τ_d) = $(round(ac.acor[end], digits = 3))" *
              (first_neg === nothing ? "; no zero crossing within 5 τ_d" : "; first negative lag $(ac.lags[first_neg]) τ_d")
        @test rate > 0

        # T4.9: HSO's second parameter set kR 0.3, C_m 120, f 7: ⟨T_r⟩ 113 ± 5, dv 0.081
        e2 = ensemble(RunSpec(; kR = 0.3, Cm = 120.0, f = 7.0), T4_LENGTH; seeds = T4_SEEDS)
        s2 = e2.summary
        hf2 = s2.histogram_fit
        ks2 = report("T4.9 HSO second run (kR 0.3, C_m 120, f 7)", e2)
        c2 = e2.runs[1].sim.model.constants
        @info "T4.9 model constants: Λ = $(round(c2.Λ, digits = 4)) (HSO 0.815), μ/Λ = $(round(c2.μ / c2.Λ, digits = 3)) (1.82); " *
              "Eq. 10 fit (γ, c²) = ($(round(hf2.eq10.γ, digits = 2)), $(round(hf2.eq10.c2, digits = 3))) at $(round(hf2.eq10.sigma, digits = 1))σ; " *
              "Eq. A.1 fit ($(round(hf2.a1.γ, digits = 2)), $(round(hf2.a1.c2, digits = 3))) at $(round(hf2.a1.sigma, digits = 1))σ (HSO case 5, fitted to both records: (5, 0.50)); " *
              "A.2 at the A.1 fit: $(round(mean_reversal_time(BistableTheory(c2, hf2.a1.γ, hf2.a1.c2)), digits = 1)) (HSO 138 at case 5)"
        @test abs(s2.mean_reversal_time - 113) < 2 * sqrt(5^2 + s2.mean_reversal_time_error^2)
        @test abs(s2.dv - 0.081) < 0.005
        @test ks2.p > 0.05                                           # T4.11 for this model too
    end

    @testset "T4.10, T4.15 the SOH variants and the ramped run" begin
        # T4.10: SOH p. 121, f = 6 and 7 with the standard model: ⟨T_r⟩ ≈ 200 and ≈ 50, relative
        # variance between reversals 0.11 and 0.24, 20%; both variance estimators are checked
        for (f, Tr, var) in ((6.0, 200.0, 0.11), (7.0, 50.0, 0.24))
            e = ensemble(RunSpec(; f), T4_LENGTH; seeds = T4_SEEDS)
            s = e.summary
            report("T4.10 SOH variant f = $f", e)
            @test relerr(s.mean_reversal_time, Tr) < 0.2
            @test relerr(s.dv, var) < 0.2
            @test relerr(s.variance, var) < 0.2
        end
        # T4.15: SOH Fig. 3, f linear from 4 to 6.5 over 5×10⁴ τ_d: the rms variability rises
        # from 17% to 41% (20% on the rms values), the chrons shorten drastically, the first
        # chron ≈ 1.9×10⁴ τ_d in the paper's single realisation (kind J)
        r = segmented(RunSpec(; schedule = LinearRamp(4.0, 6.5, 0.0, 5e4)), 5000.0, 10; seeds = T4_SEEDS)
        rms_dv = [sqrt(sg.dv) for sg in r.segments]
        rms_ex = [sqrt(sg.variance) for sg in r.segments]
        nrev = [sg.n_reversals for sg in r.segments]
        @info "T4.15 ramp 4 → 6.5 over 5×10⁴ τ_d, $(length(T4_SEEDS)) seeds, 10 segments: rms variability √dv = $(join(string.(round.(Int, 100 .* rms_dv)), ", "))%; " *
              "√(§2 variance) = $(join(string.(round.(Int, 100 .* rms_ex)), ", "))% (SOH 17% → 41%); reversals per segment (pooled) $(join(string.(nrev), ", ")); " *
              "first reversal per seed at τ = $(join(string.(round.(Int, r.first_reversal)), ", ")) (SOH's figure ≈ 19 000); reversals per seed $(join(string.(k.whole.n_reversals for k in r.runs), ", "))"
        @test relerr(rms_dv[1], 0.17) < 0.2
        @test relerr(rms_dv[end], 0.41) < 0.2
        @test relerr(rms_ex[1], 0.17) < 0.2
        @test relerr(rms_ex[end], 0.41) < 0.2
        @test issorted(rms_dv)                                       # the rise is monotonic segment by segment
        @test nrev[1] == 0 && nrev[end] > 100                        # chrons shorten drastically
        @test all(5e3 .< r.first_reversal .< 3e4)                    # a long first chron, as in the paper
    end

    @testset "T4.13, T4.14 event-level diagnostics (standard run, $(length(T4_SEEDS)) seeds); Q9 across the sweep" begin
        # The parity series at the n_c cadence and two Schmitt triggers per seed: the §2 band
        # a_h = 0.5 and a narrower a_h = 0.3 (a_e scaled with it). The trajectories are those of
        # `standard` (same seeds), so the §2 trigger reproduces its event list exactly.
        sims = [Simulation(RunSpec(; seed)) for seed in T4_SEEDS]
        obs = map(sims) do sim
            Δτ, nc = timestep(sim), redraw_interval(sim.noise)
            (trigger = Probe(amplitude, SchmittTrigger(; Δt = Δτ)),
             narrow = Probe(amplitude, SchmittTrigger(; ah = 0.3, ae = 0.15, Δt = Δτ)),
             parity = Probe(parity, Series(; sizehint = steps(sim, T4_LENGTH) ÷ nc), nc))
        end
        Threads.@threads for k in eachindex(sims)
            run!(sims[k], steps(sims[k], T4_LENGTH), obs[k])
        end
        quantile(x, p) = (y = sort(x); y[clamp(round(Int, p * (length(y) - 1)) + 1, 1, length(y))])
        median(x) = quantile(x, 0.5)
        frac(x, c) = count(<(c), x) / length(x)
        r2(x) = round(x, digits = 2)
        r3(x) = round(x, digits = 3)
        # T4.13: the parity outside the event windows of every reversal, aborted reversal and excursion
        Pall, P0, P1 = Float64[], Float64[], Float64[]
        for (sim, o) in zip(sims, obs)
            ev, P = value(o.trigger).events, value(o.parity).values
            nc, Δτ = redraw_interval(sim.noise), timestep(sim)
            t = [(nc * i) * Δτ for i in 1:length(P)]            # step index × Δτ, as the trigger times its events
            append!(Pall, P)
            append!(P0, P[outside_events(t, ev)])
            append!(P1, P[outside_events(t, ev; pad = 1.0)])
        end
        @info "T4.13 parity outside event windows ($(length(T4_SEEDS)) seeds × $T4_LENGTH τ_d, every n_c steps): median 𝒫 = $(r3(median(P0))) over $(r3(length(P0) / length(Pall))) of the samples " *
              "(quartiles $(r3(quantile(P0, 0.25))), $(r3(quantile(P0, 0.75))); 10th percentile $(r3(quantile(P0, 0.1))); fraction below −0.9: $(r3(frac(P0, -0.9)))); " *
              "with 1 τ_d padding: $(r3(median(P1))) over $(r3(length(P1) / length(Pall))); whole run: $(r3(median(Pall)))"
        @test_broken median(P0) < -0.9                 # HOS's "𝒫 ≃ −1" taken as a median below −0.9: not met, DISCREPANCIES.md DR-7
        @test -1 < median(P0) < -0.7                   # dipole-dominated between events: E_S/E_A < 0.18
        @test abs(median(P0) - median(Pall)) < 0.05    # the events occupy too little time to move the median
        @test quantile(P0, 0.25) < -0.9                # a quarter of the time the symmetric part holds < 5% of the energy
        # T4.14: reversal durations with the §2 definition (the |a| < a_h stretch) and the narrower band
        d5 = reduce(vcat, [value(o.trigger).durations for o in obs])
        d3 = reduce(vcat, [value(o.narrow).durations for o in obs])
        @info "T4.14 reversal durations: a_h = 0.5 ($(length(d5)) reversals): median $(r2(median(d5))) τ_d, quartiles $(r2(quantile(d5, 0.25)))–$(r2(quantile(d5, 0.75))), mean $(r2(sum(d5) / length(d5))), " *
              "fraction < 0.5 τ_d $(r3(frac(d5, 0.5))), < 2 τ_d $(r3(frac(d5, 2.0))), ≥ 5 τ_d $(r3(1 - frac(d5, 5.0))) (HOS: ~3 τ_d typical, fast events ~0.5 τ_d); " *
              "a_h = 0.3 ($(length(d3)) reversals): median $(r2(median(d3))), quartiles $(r2(quantile(d3, 0.25)))–$(r2(quantile(d3, 0.75))), fraction < 0.5 τ_d $(r3(frac(d3, 0.5))), < 2 τ_d $(r3(frac(d3, 2.0))) (HSO: 0.5–2 τ_d)"
        @test length(d5) == standard.summary.n_reversals          # the same trajectories as `standard`
        @test 1.5 < median(d5) < 4.5
        @test 0.1 < frac(d5, 0.5) < 0.35
        @test 0.5 ≤ median(d3) ≤ 2.0
        # Q9: the duration distribution across the forcing sweep with the §2 definition
        for f in (4.4, 5.4, 6.4, 7.4, 8.4)
            d = f == 6.4 ? d5 : value(sweep[f].pooled.trigger).durations
            isempty(d) && continue
            @info "Q9 f = $f: $(length(d)) reversals, duration median $(r2(median(d))) τ_d, quartiles $(r2(quantile(d, 0.25)))–$(r2(quantile(d, 0.75))), " *
                  "mean $(r2(sum(d) / length(d))), fraction < 0.5 τ_d $(r3(frac(d, 0.5))), ≥ 5 τ_d $(r3(1 - frac(d, 5.0)))"
            @test 0.3 < median(d) < 10
        end
    end

    @testset "T4.17 convergence of the standard run at Δτ = 0.005, N = 120 and with a second-order scheme of another family (diagnostic)" begin
        # the Δτ = 0.005 run sees the same noise realisation (F is drawn per τ_c); the Crank–Nicolson
        # run with trapezoidal noise checks that the production stepper's results do not depend on
        # the integrator family at second order (STRETCH_GOALS.md SG-1)
        ref = standard.summary
        for (label, spec) in (("Δτ = 0.005", RunSpec(; Δτ = 0.005)), ("N = 120", RunSpec(; N = 120)),
                              ("Crank–Nicolson, trapezoidal noise", RunSpec(; stepper = :crank_nicolson_implicit_noise)))
            e = ensemble(spec, T4_LENGTH; seeds = T4_SEEDS)
            s = e.summary
            report("T4.17 $label", e)
            z = (s.mean_reversal_time - ref.mean_reversal_time) / sqrt(s.mean_reversal_time_error^2 + ref.mean_reversal_time_error^2)
            @info "T4.17 $label against Δτ = 0.01, N = 60: ⟨T_r⟩ ratio $(round(s.mean_reversal_time / ref.mean_reversal_time, digits = 3)) ($(round(z, digits = 1))σ); " *
                  "variance $(round(s.variance, digits = 3)) vs $(round(ref.variance, digits = 3)); dv $(round(s.dv, digits = 3)) vs $(round(ref.dv, digits = 3)); " *
                  "D₁ fit $(round(s.d1_fit.c0, digits = 3)) + $(round(s.d1_fit.c2, digits = 3)) a² vs $(round(ref.d1_fit.c0, digits = 3)) + $(round(ref.d1_fit.c2, digits = 3)) a²"
            @test abs(z) < 3
            @test abs(s.variance - ref.variance) < 0.02
            @test abs(s.dv - ref.dv) < 0.02
            @test relerr(s.d1_fit.c0, ref.d1_fit.c0) < 0.05
        end
    end

    @testset "T4.16 alternative quenching laws (SG-2): dipole statistics \"not materially changed\"" begin
        # HOS p. 277: q from the toroidal energy, α/(1 + q), local quenching. With h = 1 − q the
        # energy and local laws reverse the sign of α under an overtone burst and run away
        # (test/quenching.jl, DISCREPANCIES.md DR-8), so the comparison uses their floored
        # variants (α suppressed, never reversed) and the inverse law, which needs no floor.
        # Result 2026-10-24 (VALIDATION.md §4): the gate's "⟨T_r⟩, variance within ~20%" holds for
        # the floored papers' law (the floor is immaterial for it), for the ⟨T_r⟩ of the inverse law
        # only (its variability is 40–50% higher: the well is shallower, the barrier the same), and
        # for neither statistic of the floored energy and local laws, which reverse 2–7 times more
        # often. The 20% criteria stay as @test_broken where they fail.
        ref = standard.summary
        for (label, q, rate_holds) in (("papers' law, floored", :amplitude_floored, true), ("α/(1 + q)", :inverse, true),
                                       ("q from the toroidal energy, floored", :energy_floored, false), ("local quenching, floored", :local_floored, false))
            e = ensemble(RunSpec(; quenching = q), T4_LENGTH; seeds = T4_SEEDS)
            s = e.summary
            report("T4.16 $label", e)
            @info "T4.16 $label against the papers' law: ⟨T_r⟩ ratio $(round(s.mean_reversal_time / ref.mean_reversal_time, digits = 3)); " *
                  "§2 variance $(round(s.variance, digits = 3)) vs $(round(ref.variance, digits = 3)); dv $(round(s.dv, digits = 3)) vs $(round(ref.dv, digits = 3)); " *
                  "D₁ fit $(round(s.d1_fit.c0, digits = 2)) + $(round(s.d1_fit.c2, digits = 2)) a² vs $(round(ref.d1_fit.c0, digits = 2)) + $(round(ref.d1_fit.c2, digits = 2)) a²"
            @test s.n_reversals > 100
            if q === :amplitude_floored
                @test relerr(s.mean_reversal_time, ref.mean_reversal_time) < 0.05
                @test relerr(s.dv, ref.dv) < 0.05
                @test relerr(s.variance, ref.variance) < 0.05
                continue
            end
            if rate_holds
                @test relerr(s.mean_reversal_time, ref.mean_reversal_time) < 0.2
            else
                @test_broken relerr(s.mean_reversal_time, ref.mean_reversal_time) < 0.2
                @test 0.1 < s.mean_reversal_time / ref.mean_reversal_time < 0.6     # two to seven times the reversal rate, stably
            end
            @test_broken relerr(s.dv, ref.dv) < 0.2
            @test_broken relerr(s.variance, ref.variance) < 0.2
            @test s.dv > 1.2 * ref.dv                                                # every alternative raises the variability
        end
    end

    @testset "T4.12 forcing-parameter invariance at f²τ_c/N_c² = 2.05×10⁻² (diagnostic)" begin
        # HOS p. 272, Appendix E: all statistics depend on (f, τ_c, N_c) only through f²τ_c/N_c² in
        # the δ-correlated, continuum-bin limit; finite τ_c and bin width break the invariance
        # (DISCREPANCIES.md DR-3 for the bins)
        configs = ((f = 6.4, τc = 0.05, Nc = 10), (f = 4.525, τc = 0.1, Nc = 10), (f = 12.8, τc = 0.05, Nc = 20))
        es = (standard, map(c -> ensemble(RunSpec(; c.f, c.τc, c.Nc), T4_LENGTH; seeds = T4_SEEDS), configs[2:end])...)
        for (c, e) in zip(configs, es)
            report("T4.12 (f, τ_c, N_c) = ($(c.f), $(c.τc), $(c.Nc))", e)
            @test e.summary.n_reversals > 100
        end
        ref = es[1].summary
        for (c, e) in zip(configs[2:end], es[2:end])
            s = e.summary
            z = (s.mean_reversal_time - ref.mean_reversal_time) / sqrt(s.mean_reversal_time_error^2 + ref.mean_reversal_time_error^2)
            @info "T4.12 ($(c.f), $(c.τc), $(c.Nc)) against the standard run: ⟨T_r⟩ ratio $(round(s.mean_reversal_time / ref.mean_reversal_time, digits = 3)) ($(round(z, digits = 1))σ); " *
                  "variance ratio $(round(s.variance / ref.variance, digits = 3)); c² ratio $(round(s.c2 / ref.c2, digits = 3)); γ ratio $(round(s.γ / ref.γ, digits = 3))"
        end
    end
end
