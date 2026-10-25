# Tier 6: the Sint-800 comparison (gates T6.1–T6.5 in dev/VALIDATION.md), HSO §4–5 applied
# to the record in data/sint800.tsv (Guyodo & Valet 1999). T6.4 regenerates HSO's two runs in
# stream mode (the reference run f = 5.4 for 5×10⁴ τ_d, the second run kR 0.3, C_m 120, f 7 for
# 10⁴ τ_d), about 20 s on one thread. Λ = 0.527 and μ/Λ = 1.45 are the standard model's rates.

using Random

const SINT = read_sint800()
const Λ_SINT, μ_SINT = 0.527, 0.763

@testset "Tier 6: Sint-800 comparison (HSO §4–5)" begin
    d = SINT
    @test d.n == 799 && d.age[1] == 2 && d.age[end] == 800 && all(diff(d.age) .== 1)
    @test all(d.vadm .> 0) && all(d.vadm_err .> 0)

    @testset "T6.1 dv of the record with Monte-Carlo errors" begin
        dv = dipole_variation(d.vadm)
        mc = montecarlo(dipole_variation, d.vadm, d.vadm_err; n = 2000, rng = Xoshiro(1))
        @info "T6.1 dv = $(round(dv, digits = 4)) (paper 0.071 ± 0.003); Monte-Carlo spread $(round(mc.sigma, digits = 4)), replica mean $(round(mc.mean, digits = 4)); median relative error $(round(sum(d.vadm_err ./ d.vadm) / d.n, digits = 3))"
        @test abs(dv - 0.071) < 0.002
        @test 0.0015 < mc.sigma < 0.005
        @test mc.mean > dv                                   # added noise adds variance: the spread, not the mean, is the error
    end

    @testset "T6.2 χ² fit of Eq. 10 to the VADM record, ~100 equal-population bins (Table 1)" begin
        fit = vadm_fit(d.vadm, Λ_SINT; nbins = 100)
        hso = vadm_chi2(d.vadm, BistableTheory(6.3, 0.32; Λ = Λ_SINT), 5.98; bins = vadm_bins(d.vadm; nbins = 100))
        @info "T6.2 Eq. 10 optimum (γ, c², A_m) = ($(round(fit.γ, digits = 2)), $(round(fit.c2, digits = 3)), $(round(fit.Am, digits = 2))) at $(round(fit.sigma, digits = 2))σ over $(fit.nbins) bins (paper 6.3, 0.32, 5.98 at 1σ); the paper's values give $(round(hso.sigma, digits = 2))σ here; theory dv $(round(fit.theory.dv, digits = 3)), ⟨T_r⟩ $(round(fit.theory.mean_reversal_time, digits = 0)) τ_d"
        @test fit.sigma < 1.5
        @test hso.sigma < 1.5                                 # HSO's optimum lies in the same χ² valley
        @test abs(fit.Am - 5.98) < 0.2
        @test 5 < fit.γ < 7.5 && 0.2 < fit.c2 < 0.45          # the valley is degenerate in (γ, c²) (HSO p. 150)
        # cases 2–5 of Table 1: A_m alone optimised at the printed (γ, c²); the ordering of the
        # goodness of fit is the paper's (case 5, the second model, is the worst)
        c2 = vadm_fit(d.vadm, Λ_SINT; nbins = 100, fixed = (8.0, 0.50))
        c3 = vadm_fit(d.vadm, Λ_SINT; nbins = 100, fixed = (7.0, 0.45))
        c4 = vadm_fit(d.vadm, Λ_SINT, μ_SINT; nbins = 100, fixed = (4.8, 0.22))
        c5 = vadm_fit(d.vadm, 0.815, 1.82 * 0.815; nbins = 100, fixed = (5.0, 0.50))
        @info "T6.2 Table 1 cases 2–5 at the printed (γ, c²): σ = $(round(c2.sigma, digits = 2)), $(round(c3.sigma, digits = 2)), $(round(c4.sigma, digits = 2)), $(round(c5.sigma, digits = 2)) (paper 1.3, 1.1, 1.4, 5); A_m = $(round(c2.Am, digits = 2)), $(round(c3.Am, digits = 2)), $(round(c4.Am, digits = 2)), $(round(c5.Am, digits = 2)); case 4 theory dv $(round(c4.theory.dv, digits = 3)) (paper 0.071)"
        @test c2.sigma < 2 && c3.sigma < 2 && c4.sigma < 2
        @test c5.sigma > max(c2.sigma, c3.sigma, c4.sigma)
        @test abs(c4.theory.dv - 0.071) < 0.002
        # the improved Eq. A.1 with the model's μ fits at least as well as Eq. 10
        fA = vadm_fit(d.vadm, Λ_SINT, μ_SINT; nbins = 100)
        @info "T6.2 Eq. A.1 optimum (γ, c², A_m) = ($(round(fA.γ, digits = 2)), $(round(fA.c2, digits = 3)), $(round(fA.Am, digits = 2))) at $(round(fA.sigma, digits = 2))σ"
        @test fA.sigma < 1.5
    end

    @testset "T6.3 autocorrelation of the record (Fig. 7)" begin
        for (label, series) in (("unsigned", d.vadm), ("signed beyond 780 ka", ifelse.(d.age .> 780, -d.vadm, d.vadm)))
            ac = hso_autocorrelation(reverse(series), 50)
            k0 = findfirst(<(0), ac) - 1                      # first lag with a negative value
            kmin = argmin(ac) - 1
            @info "T6.3 $label: zero crossing between $(k0 - 1) and $k0 kyr (paper ≈ 27), minimum $(round(minimum(ac), digits = 3)) at $kmin kyr (paper ≈ −0.2 near 40)"
            @test ac[1] == 1 && 24 ≤ k0 ≤ 32
            @test -0.26 < minimum(ac) < -0.12 && 34 ≤ kmin ≤ 48
        end
        mc = montecarlo(x -> hso_autocorrelation(reverse(x), 50), d.vadm, d.vadm_err; n = 500, rng = Xoshiro(2))
        @info "T6.3 Monte-Carlo error bars at lags 1, 10, 25, 50: $(round.(mc.sigma[[2, 11, 26, 51]], digits = 3)) (paper: of order 799^(−1/2) ≈ 0.035)"
        @test all(0.005 .< mc.sigma[2:end] .< 0.05)
    end

    @testset "T6.5 the earlier fits: HOS Fig. 10 (35 bins) and SOH Fig. 4 (30 bins)" begin
        hos = vadm_fit(d.vadm, Λ_SINT; nbins = 35, equal_population = false)
        hos_at = vadm_fit(d.vadm, Λ_SINT; nbins = 35, equal_population = false, fixed = (5.15, 0.25))
        soh = vadm_fit(d.vadm, Λ_SINT; nbins = 30, equal_population = false)
        soh_at = vadm_fit(d.vadm, Λ_SINT; nbins = 30, equal_population = false, fixed = (5.2, 0.26))
        @info "T6.5 HOS Fig. 10, 35 equal-width bins ($(hos.nbins) non-empty): optimum ($(round(hos.γ, digits = 2)), $(round(hos.c2, digits = 3)), $(round(hos.Am, digits = 2))) at $(round(hos.sigma, digits = 2))σ; the printed (5.15, 0.25) give $(round(hos_at.sigma, digits = 2))σ with A_m = $(round(hos_at.Am, digits = 2)) (paper: 1σ)"
        @info "T6.5 SOH Fig. 4, 30 equal-width bins ($(soh.nbins) non-empty): optimum ($(round(soh.γ, digits = 2)), $(round(soh.c2, digits = 3)), $(round(soh.Am, digits = 2))) at $(round(soh.sigma, digits = 2))σ; the printed (5.2, 0.26) give $(round(soh_at.sigma, digits = 2))σ with A_m = $(round(soh_at.Am, digits = 2)) (paper: within 1σ, a_max 5.9)"
        @test soh_at.sigma < 2 && abs(soh_at.Am - 5.9) < 0.2
        @test hos.sigma < 2.5 && abs(hos.Am - 5.9) < 0.3
        @test_broken hos_at.sigma < 1.5                      # 3.2σ with our binning: DISCREPANCIES.md DR-9
        @test hos_at.sigma < 4
    end

    @testset "T6.4 τ_d by the autocorrelation overlap (τ_av = 5 kyr)" begin
        ac = hso_autocorrelation(reverse(d.vadm), 50)
        function amplitude_series(spec, τ)
            sim = Simulation(spec)
            nc = redraw_interval(sim.noise)
            probe = Probe(amplitude, Series(; sizehint = steps(sim, τ) ÷ nc + 1), nc)
            run!(sim, steps(sim, τ), (probe,))
            return value(probe).values, nc * timestep(sim)
        end
        a, dt = amplitude_series(RunSpec(; f = 5.4, seed = 1), 5e4)
        sc = tau_d_scan(a, dt, ac; τd = 4:0.5:24, nwindows = 200)
        @info "T6.4 reference run (f = 5.4, 5×10⁴ τ_d): best τ_d per 800-kyr window has median $(sc.median) kyr, quartiles $(sc.q25)–$(sc.q75); the window-averaged misfit is least at $(sc.best_mean) kyr (paper 11 ± 2)"
        @test 7 ≤ sc.median ≤ 15 && 7 ≤ sc.best_mean ≤ 15
        @test sc.q25 ≥ 5 && sc.q75 ≤ 20
        a2, dt2 = amplitude_series(RunSpec(; kR = 0.3, Cm = 120.0, f = 7.0, seed = 1), 1e4)
        sc2 = tau_d_scan(a2, dt2, ac; τd = 4:0.5:30, nwindows = 100)
        @info "T6.4 second run (kR 0.3, C_m 120, f 7, 10⁴ τ_d): median $(sc2.median) kyr, quartiles $(sc2.q25)–$(sc2.q75), least mean misfit at $(sc2.best_mean) kyr (paper 16 ± 4)"
        @test 8 ≤ sc2.median ≤ 24 && 8 ≤ sc2.best_mean ≤ 24
        @test sc2.median > sc.median                          # the larger Λ of the second model needs the longer τ_d (HSO p. 153)
    end
end
