# src/stats.jl: the KS test for exponential intervals (T4.11) and the χ² fit of p(a) to a
# histogram (HSO §4.3, T4.8), on synthetic data with known answers.

@testset "Statistics for the gates" begin
    rng = Xoshiro(7)

    @testset "ks_exponential" begin
        x = 5 .+ 100 .* randexp(rng, 1_400)
        k = ks_exponential(x; tmin = 5.0, nboot = 1_000, rng = Xoshiro(1))
        @test k.n == 1_400
        @test k.mean ≈ sum(x .- 5) / 1_400
        @test k.p > 0.05                            # an exponential sample passes
        @test 0 < k.D < 0.05
        y = 5 .+ 100 .* rand(rng, 1_400)            # a uniform sample fails
        @test ks_exponential(y; tmin = 5.0, nboot = 500, rng = Xoshiro(1)).p < 0.01
        @test ks_exponential([1.0, 2.0, 3.0, 9.0]; nboot = 10).n == 4
        @test_throws ArgumentError ks_exponential([1.0, 2.0]; tmin = 5.0)
        @test_throws ArgumentError ks_exponential([1.0, 2.0, 3.0]; nboot = 0)
        # deterministic for a given rng
        @test ks_exponential(x; tmin = 5.0, nboot = 200, rng = Xoshiro(3)).p == ks_exponential(x; tmin = 5.0, nboot = 200, rng = Xoshiro(3)).p
    end

    @testset "nelder_mead" begin
        rosen(x) = (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
        x, fx = nelder_mead(rosen, [-1.2, 1.0]; step = 0.5, xtol = 1e-9, ftol = 1e-14, maxiter = 5_000)
        @test fx < 1e-10
        @test x ≈ [1.0, 1.0] atol = 1e-4
        x2, f2 = nelder_mead(v -> (v[1] - 3)^2, [0.0])
        @test x2[1] ≈ 3 atol = 1e-5
    end

    @testset "fold_histogram, equal_population_bins" begin
        edges = collect(-2.0:0.5:2.0)                 # 8 bins
        counts = [1, 2, 3, 4, 10, 20, 30, 40]
        fe, fc = fold_histogram(edges, counts)
        @test fe == [0.0, 0.5, 1.0, 1.5, 2.0]
        @test fc == [14, 23, 32, 41]
        @test_throws ArgumentError fold_histogram(collect(-2.0:0.5:1.5), counts[1:7])
        @test_throws ArgumentError fold_histogram(collect(-1.5:0.5:2.5), counts)
        e = collect(0.0:0.1:1.0); c = [0, 0, 5, 5, 5, 5, 5, 5, 5, 0]
        b = equal_population_bins(e, c; nbins = 3)
        @test sum(b.n) == 35 && all(b.n .> 0)
        @test b.lo[1] == 0.2 && b.hi[end] == 0.9             # empty ends dropped
        @test all(e[b.first] .== b.lo) && all(e[b.last .+ 1] .== b.hi)
        @test_throws ArgumentError equal_population_bins(e, zeros(Int, 10))
    end

    @testset "density_fit recovers (γ, c²) from a noise-free histogram" begin
        # a histogram of |a| at 0.005 resolution (the standard observers' resolution) filled
        # with N × ∫p over each bin
        edges = collect(0.0:0.005:2.5)
        for (γ, c2, μfac) in ((8.0, 0.5, 1.0), (2.96, 0.246, 1.0), (4.8, 0.22, 1.45))
            th = BistableTheory(γ, c2; Λ = 0.527, μ = μfac * 0.527)
            ψm = psi_max(th)
            dens = [exp(psi(th, a) - ψm) for a in 0.0005:0.001:2.5]        # fine quadrature inside the bins
            cum = cumsum(dens) .* 0.001
            counts = round.(Int, 5e4 .* diff([0.0; cum[5:5:end]]) ./ cum[end])
            fit = density_fit(edges, counts, 0.527, μfac * 0.527)
            @test fit.γ ≈ γ rtol = 0.02
            @test fit.c2 ≈ c2 rtol = 0.03
            @test fit.chi2 < 0.2 * fit.ν                 # rounding noise only (σ is two-sided, so it is large here)
            @test 85 ≤ fit.nbins ≤ 115                 # "about 100" equal-population bins
        end
        # a dead dynamo: every sample in the first bin, one equal-population bin, no degrees of
        # freedom; the fit is meaningless but must not throw (σ is NaN)
        counts = zeros(Int, 500); counts[1] = 5000
        fit = density_fit(edges, counts, 0.527)
        @test fit.nbins == 1 && fit.ν < 0 && isnan(fit.sigma)
    end
end
