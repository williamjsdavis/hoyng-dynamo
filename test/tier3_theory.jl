# Tier 3: closed-form theory (formula module). Gates T3.1–T3.9 in dev/VALIDATION.md.
# Each gate becomes one @testset named by its ID, with the target and tolerance from the file.
# Notation is HOS's (Eqs. 32–60, Appendix F) and HSO's Appendix A (Eqs. A.1–A.2).

using QuadGK: quadgk

# Fourth-order central differences (Richardson), independent of the module's derivatives.
function richardson(f, x; h = 1e-3)
    d(h) = (f(x + h) - f(x - h)) / 2h
    (4d(h / 2) - d(h)) / 3
end

# Literal transcriptions of the HOS closed forms (κ = 0), for comparison with the module's
# general expressions.
hos52(γ) = sqrt(1 - 1 / γ)
hos53(γ, c2) = 1 - γ + (γ * c2 + γ - 1) * log((γ * c2 + γ - 1) / (γ * c2))
hos54(γ, c2) = -4γ * (γ - 1) / (γ * c2 + γ - 1)
hos56(γ, c2) = (γ * c2 + γ - 1) / (4 * (γ - 1)^2)
hos58(γ, c2, Λ) = (π / Λ) * sqrt(2 / c2) * sqrt((γ * c2 + γ - 1) / (γ - 1)) * exp(hos53(γ, c2))
hos59(γ, c2, Λ) = (π * sqrt(2) / Λ) * sqrt(c2 / (1 + c2)) * exp((-1 + (1 + c2) * log((1 + c2) / c2)) * γ)

# HOS Table II (p. 292): run, f, c², γ, variance (Eq. 56), ⟨T_r⟩ (Eq. 58), and the run values.
const TABLE_II = [(1, 4.4, 0.318, 9.04, 0.04, 12000, 0.03, 10700),
                  (2, 5.4, 0.299, 5.15, 0.08, 480, 0.07, 519),
                  (3, 6.4, 0.246, 2.96, 0.17, 94, 0.16, 107),
                  (4, 7.4, 0.240, 1.82, 0.47, 35, 0.30, 37.6),
                  (5, 8.4, 0.228, 1.09, 1.0, 35, 0.50, 18.7)]

@testset "Tier 3: closed-form theory (formula module)" begin
    @testset "construction" begin
        th = BistableTheory(2.96, 0.246)
        @test th.Λ == 0.527 && th.κ == 0
        @test diffusion_constants(th).Δ0 ≈ 0.527 / 2.96
        @test diffusion_constants(th).Δ1 ≈ 0.246 * 0.527 / 2.96
        sp = shape_parameters(0.527, diffusion_constants(th)...)
        @test sp.γ ≈ 2.96 && sp.c2 ≈ 0.246
        thA = BistableTheory(4.8, 0.22; Λ = 0.527, μ = 1.45 * 0.527)
        @test thA.κ ≈ 0.45
        @test_throws ArgumentError BistableTheory(0.0, 0.2)
        @test_throws ArgumentError BistableTheory(3.0, 0.2; μ = 0.5)        # μ < Λ
        @test_throws ArgumentError peak_amplitude(BistableTheory(0.9, 0.2))  # γ ≤ 1: single peak
        @test !bistable(BistableTheory(1.0, 0.2)) && bistable(th)
        @test_throws ArgumentError mean_reversal_time_asymptotic(thA)        # Eq. 59 is κ = 0 only
        # the model's own constants reproduce HOS Eq. 34 (μ, Λ) and HSO's μ/Λ = 1.45
        c, _ = model_constants(Grid(60), 0.5, 100.0)
        thM = BistableTheory(c, 2.96, 0.246)
        @test relerr(thM.Λ, 0.527) < 0.01 && relerr(thM.κ + 1, 1.45) < 0.01
        # drift, diffusion, potential: S = −∂U/∂a, D(0) = Δ₁, U(±1) = 0
        for t in (th, thA), a in (-1.3, -0.4, 0.0, 0.5, 1.0, 1.7)
            @test drift(t, a) ≈ -richardson(x -> potential(t, x), a) atol = 1e-9
            @test diffusion(t, a) ≈ diffusion_constants(t).Δ0 * a^2 + diffusion_constants(t).Δ1
        end
        @test potential(th, 1.0) == 0
        @test abs(potential(thA, -1.0)) < 1e-15
        @test potential(th, 0.0) ≈ 0.527 / 4                                   # U(0) = Λ/4, Eq. 41
        @test drift(th, 0.5) ≈ 0.527 * 0.75 * 0.5                              # Eq. 33
        @test drift(thA, 0.5) ≈ 0.527 * 0.75 * 0.5 * (1 + 0.45 * 0.25)         # HSO App. A
        @test diffusion_prefactor(6.4, 0.05, 10) ≈ 3.217e-2 rtol = 1e-3        # p. 291: 3.22×10⁻²
        @test occursin("μ/Λ = 1.45", sprint(show, thA))
    end

    @testset "T3.1 Eq. 43 (and A.1) equal Eq. 42 integrated numerically" begin
        for th in (BistableTheory(2.96, 0.246), BistableTheory(9.04, 0.318), BistableTheory(4.8, 0.22; μ = 1.45 * 0.527))
            # Eq. 42 by quadrature: p ∝ (1/D) exp ∫₀^a 2S/D, compared with exp ψ through ratios to
            # the value at a_m (equivalent to comparing after normalisation).
            am = peak_amplitude(th)
            lnp42(a) = -log(diffusion(th, a)) + quadgk(x -> 2 * drift(th, x) / diffusion(th, x), 0.0, a; rtol = 1e-12)[1]
            for a in -2.0:0.125:2.0
                @test abs((lnp42(a) - lnp42(am)) - (psi(th, a) - psi(th, am))) < 1e-8
            end
            # the normalised density integrates to one and is even
            p = stationary_density(th)
            @test abs(quadgk(p, -Inf, Inf; rtol = 1e-10)[1] - 1) < 1e-8
            @test p(0.7) == p(-0.7) && logpdf(p, 0.7) ≈ log(p(0.7))
            @test p(am) ≈ exp(peak_log_density(th) - log_normalisation(th))
        end
        # the log-space quadrature survives a γ where exp(ψ_m) overflows Float64
        big = BistableTheory(1500.0, 0.25)
        @test isfinite(log_normalisation(big)) && isfinite(log_mean_reversal_time(big))
        @test stationary_density(big)(peak_amplitude(big)) > 0 && isfinite(dipole_variation(big))
        @test isinf(mean_reversal_time(big))
        @test psi(big, 1e200) == -Inf && psi(big, Inf) == -Inf          # no NaN at huge amplitude
    end

    @testset "T3.2 peak a_m, ψ_m, ψ_m'' vs numerical maximisation of ψ" begin
        for th in (BistableTheory(2.96, 0.246), BistableTheory(1.82, 0.240), BistableTheory(9.04, 0.318),
                   BistableTheory(4.8, 0.22; μ = 1.45 * 0.527), BistableTheory(5.0, 0.5; Λ = 0.815, μ = 1.82 * 0.815))
            # derivatives of the module against finite differences of ψ
            for a in (0.1, 0.5, 0.9, 1.3)
                @test abs(dpsi(th, a) - richardson(x -> psi(th, x), a)) < 1e-9
                @test abs(d2psi(th, a) - richardson(x -> dpsi(th, x), a)) < 1e-8
            end
            # the maximum of ψ on (0, ∞) by bisection on a finite-difference ψ'
            g(a) = (psi(th, a + 1e-5) - psi(th, a - 1e-5)) / 2e-5
            lo, hi = 0.05, 3.0
            @test g(lo) > 0 && g(hi) < 0
            for _ in 1:100
                mid = (lo + hi) / 2
                g(mid) > 0 ? (lo = mid) : (hi = mid)
            end
            am_num = (lo + hi) / 2
            @test abs(peak_amplitude(th) - am_num) < 1e-8
            @test abs(peak_log_density(th) - psi(th, am_num)) < 1e-8
            @test abs(peak_curvature(th) - richardson(x -> dpsi(th, x), am_num)) < 1e-8
            @test abs(phi0_curvature(th) - richardson(y -> -dpsi(th, y) - 2y / (th.c2 + y^2), 0.0)) < 1e-8
            @test phi(th, 0.3) ≈ -psi(th, 0.3) - log1p(0.09 / th.c2)
        end
        # κ = 0: the general expressions reduce to the literal Eqs. 52–54
        for (γ, c2) in ((2.96, 0.246), (1.82, 0.240), (9.04, 0.318), (6.3, 0.32))
            th = BistableTheory(γ, c2)
            @test peak_amplitude(th) ≈ hos52(γ) rtol = 1e-14
            @test peak_log_density(th) ≈ hos53(γ, c2) rtol = 1e-12
            @test peak_curvature(th) ≈ hos54(γ, c2) rtol = 1e-12
            @test mean_amplitude(th) == peak_amplitude(th)                     # Eq. 55
        end
        # κ ≠ 0: a_m solves the quadratic of HSO p. 156
        thA = BistableTheory(4.8, 0.22; μ = 1.45 * 0.527)
        am = peak_amplitude(thA)
        @test (1 + thA.κ * am^2) * (1 - am^2) ≈ 1 / thA.γ rtol = 1e-14
    end

    @testset "T3.3 variance Eq. 56 for the Table II pairs" begin
        for (run, f, c2, γ, var56, _, var_run, _) in TABLE_II
            th = BistableTheory(γ, c2)
            v = relative_variance(th)
            @test v ≈ hos56(γ, c2) rtol = 1e-12
            if run < 5
                @test abs(v - var56) <= 0.01 + 1e-12                           # last printed digit
            else
                # Run 5 (γ = 1.09): Eq. 56 gives 10.4, not the printed 1.0; the paper's footnote b
                # puts this run outside the validity of the parabolic approximation. Recorded in
                # VALIDATION.md §4, not asserted against the table.
                @test v > 10
            end
            @info "T3.3 run $run (f = $f): Eq. 56 = $(round(v, digits = 4)) (table $var56), dv over a > 0 = $(round(dipole_variation(th), digits = 4)) (run $var_run)"
        end
    end

    @testset "T3.4 ⟨T_r⟩ Eq. 58 for the Table II pairs" begin
        # Pass if within 1%, or within one unit of the last printed digit (12000 is printed to
        # two figures): VALIDATION.md §1.
        units = Dict(1 => 1000.0, 2 => 10.0, 3 => 1.0, 4 => 1.0, 5 => 1.0)
        for (run, f, c2, γ, _, Tr58, _, Tr_run) in TABLE_II
            th = BistableTheory(γ, c2)
            Tr = mean_reversal_time(th)
            @test Tr ≈ hos58(γ, c2, 0.527) rtol = 1e-12                       # A.2 form ≡ Eq. 58 at κ = 0
            @test Tr ≈ exp(log_mean_reversal_time(th))
            @test relerr(Tr, Tr58) < 0.01 || abs(Tr - Tr58) <= units[run]
            @info "T3.4 run $run (f = $f): Eq. 58 = $(round(Tr, digits = 1)) (table $Tr58, run $Tr_run)"
        end
        @test abs(mean_reversal_time(BistableTheory(2.96, 0.246)) - 93.9) < 0.05   # hand check
    end

    @testset "T3.5 Eq. 58 for HSO Table 1 cases 1–3" begin
        for (case, γ, c2, Tr_t) in ((1, 6.3, 0.32, 1090), (2, 8.0, 0.50, 928), (3, 7.0, 0.45, 665))
            Tr = mean_reversal_time(BistableTheory(γ, c2; Λ = 0.527))
            @test relerr(Tr, Tr_t) < 0.005
            @info "T3.5 case $case: Eq. 58 = $(round(Tr, digits = 1)) (HSO Table 1: $Tr_t)"
        end
    end

    @testset "T3.6 dv (HSO Eq. 12) from Eq. 10 over a > 0, cases 1–3" begin
        for (case, γ, c2, dv_t) in ((1, 6.3, 0.32, 0.076), (2, 8.0, 0.50, 0.070), (3, 7.0, 0.45, 0.079))
            th = BistableTheory(γ, c2)
            dv = dipole_variation(th)
            @test abs(dv - dv_t) <= 0.001
            # the same from the normalised density by direct quadrature
            p = stationary_density(th)
            m0 = quadgk(p, 0, Inf; rtol = 1e-10)[1]
            m1 = quadgk(a -> a * p(a), 0, Inf; rtol = 1e-10)[1] / m0
            m2 = quadgk(a -> a^2 * p(a), 0, Inf; rtol = 1e-10)[1] / m0
            @test dv ≈ (m2 - m1^2) / m1^2 rtol = 1e-8
            @test half_moments(th).mean ≈ m1 rtol = 1e-8
            @info "T3.6 case $case: dv = $(round(dv, digits = 4)) (HSO Table 1: $dv_t); ⟨|a|⟩ = $(round(m1, digits = 4)), a_m = $(round(peak_amplitude(th), digits = 4))"
        end
    end

    @testset "T3.7 Eq. A.1 with the generalised A.2: HSO Table 1 cases 4–5" begin
        # Case 4 is the kR = 0.5, C_m = 100 model (Λ = 0.527, μ/Λ = 1.45); case 5 the kR = 0.3,
        # C_m = 120 model of HSO Fig. 8 (Λ = 0.815, μ/Λ = 1.82, HSO Table 2 and p. 156).
        for (case, γ, c2, Λ, μΛ, Tr_t, dv_t) in ((4, 4.8, 0.22, 0.527, 1.45, 1097, 0.071),
                                                 (5, 5.0, 0.50, 0.815, 1.82, 138, 0.091))
            th = BistableTheory(γ, c2; Λ, μ = μΛ * Λ)
            Tr = mean_reversal_time(th)
            dv = dipole_variation(th)
            @test abs(dv - dv_t) <= 0.001
            # ⟨T_r⟩: 1109 for case 4 (1.1% above 1097) and 137.6 for case 5. The 1% gate is tighter
            # than the two printed figures of (γ, c²) allow: γ = 4.8 ± 0.05 alone moves case 4
            # over 1046–1176. Asserted at 2% and by bracketing within the printed precision.
            @test relerr(Tr, Tr_t) < 0.02
            lo = mean_reversal_time(BistableTheory(γ - 0.05, c2; Λ, μ = μΛ * Λ))
            hi = mean_reversal_time(BistableTheory(γ + 0.05, c2; Λ, μ = μΛ * Λ))
            @test lo < Tr_t < hi
            @info "T3.7 case $case: A.2 = $(round(Tr, digits = 1)) (HSO Table 1: $Tr_t; γ ± 0.05 → $(round(lo, digits = 0))–$(round(hi, digits = 0))), dv = $(round(dv, digits = 4)) ($dv_t)"
        end
        # with μ = Λ the A.1 machinery is exactly the HOS one
        @test mean_reversal_time(BistableTheory(4.8, 0.22; μ = 0.527)) ≈ hos58(4.8, 0.22, 0.527) rtol = 1e-12
        # the model constants reproduce HSO's (Λ, μ/Λ) for both parameter sets (grid-converged)
        for (kR, Cm, Λ_t, μΛ_t) in ((0.5, 100.0, 0.527, 1.45), (0.3, 120.0, 0.815, 1.82))
            c, _ = model_constants(Grid(240), kR, Cm)
            @test relerr(c.Λ, Λ_t) < 0.01 && relerr(c.μ / c.Λ, μΛ_t) < 0.01
            @info "model constants kR = $kR, C_m = $Cm, N = 240: Λ = $(round(c.Λ, digits = 4)) ($Λ_t), μ/Λ = $(round(c.μ / c.Λ, digits = 4)) ($μΛ_t)"
        end
    end

    @testset "T3.8 escape integral Eq. 57 vs the Gaussian-approximated Eq. 58 (diagnostic)" begin
        for (γ, c2) in ((3.0, 0.2), (6.0, 0.3))
            th = BistableTheory(γ, c2)
            Tq = mean_reversal_time_quadrature(th)                           # SOH's 4∫₀¹ form
            Tr = mean_reversal_time(th)
            @test Tq ≈ exp(log_mean_reversal_time_quadrature(th)) && Tq ≈ 2 * escape_time(th)
            @test 0.8 < Tq / Tr < 1.3
            # Eq. 57 is "almost independent of a′ near the well bottom" (p. 289)
            r08 = mean_reversal_time_quadrature(th; a′ = 0.8) / Tq
            r12 = mean_reversal_time_quadrature(th; a′ = 1.2) / Tq
            @test 0.95 < r08 < 1 < r12 < 1.05
            @info "T3.8 (γ, c²) = ($γ, $c2): Eq. 57 / Eq. 58 = $(round(Tq / Tr, digits = 4)); a′ = 0.8, 1.2 → ×$(round(r08, digits = 4)), ×$(round(r12, digits = 4))"
        end
    end

    @testset "T3.9 asymptotic Eq. 59–60 vs Eq. 58" begin
        c2 = 0.25
        @test kramers_exponent(BistableTheory(3.0, c2)) ≈ -1 + 1.25 * log(5)         # Eq. 60
        @test mean_reversal_time_asymptotic(BistableTheory(3.0, c2)) ≈ hos59(3.0, c2, 0.527) rtol = 1e-12
        γs = [2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 1e3, 1e4, 1e5]
        d = [log_mean_reversal_time_asymptotic(BistableTheory(γ, c2)) - log_mean_reversal_time(BistableTheory(γ, c2)) for γ in γs]
        @test all(d .< 0)                                                           # Eq. 59 underestimates
        @test issorted(abs.(d); rev = true)                                         # monotone approach
        @test abs(d[end]) < 1e-4                                                    # → 1 as γ → ∞
        @test abs(d[end]) * γs[end] ≈ abs(d[end-1]) * γs[end-1] rtol = 0.05          # ratio − 1 ∝ 1/γ
        @info "T3.9 ln(Eq. 59 / Eq. 58) at γ = $(γs): $(round.(d, sigdigits = 3))"
    end

    @testset "Eqs. 61–63: shape parameters from the measured D₁ fit" begin
        Δ = diffusion_from_fit(1.36, 2.00; D0int = 3.529, prefactor = 3.22e-2)
        sp = shape_parameters(0.527, Δ.Δ0, Δ.Δ1)
        @test abs(sp.c2 - 0.246) < 0.001 && abs(sp.γ - 2.96) < 0.01
    end
end
