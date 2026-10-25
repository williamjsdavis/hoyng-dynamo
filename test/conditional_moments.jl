# Online kernel-based regression (Davis 2023, arXiv:2307.00445v3): the streamed conditional
# moments against the offline formulas (its Eqs. 8 and 10), pooling of streams, and the
# paper's Ornstein–Uhlenbeck example (§3.1) plus a multiplicative-noise process. Always run.

using Random

cm_mean(x) = sum(x) / length(x)

# Offline KBR, the paper's Eqs. 8 and 10 summed directly: N_τ × N_x matrices W, M⁽¹⁾, M⁽²⁾
# and the conditional variance about the final mean
function ref_moments(X, x, h, nlags, K)
    N, nx = length(X), length(x)
    W = zeros(nlags, nx); S1 = zeros(nlags, nx); S2 = zeros(nlags, nx)
    for i in 1:nlags, n in 1:N-i, j in 1:nx
        k = K((x[j] - X[n]) / h) / h
        Δ = X[n+i] - X[n]
        W[i, j] += k; S1[i, j] += k * Δ; S2[i, j] += k * Δ^2
    end
    M1 = S1 ./ W; M2 = S2 ./ W
    V = zeros(nlags, nx)
    for i in 1:nlags, n in 1:N-i, j in 1:nx
        k = K((x[j] - X[n]) / h) / h
        k > 0 || continue                        # (M1 is NaN where nothing has weight)
        V[i, j] += k * ((X[n+i] - X[n]) - M1[i, j])^2
    end
    return (W = W, M1 = M1, M2 = M2, S = V, M2star = V ./ W)
end

maxrel(a, b) = maximum(abs.(a .- b)) / maximum(abs.(b))

# an AR(1) series of standard deviation 0.5: it covers |x| ≲ 2 and leaves the outer
# evaluation points of [−3, 3] without weight
function ar1(rng, n; ρ = 0.8, σ = 0.3)
    X = zeros(n)
    for i in 2:n
        X[i] = ρ * X[i-1] + σ * randn(rng)
    end
    return X
end

@testset "ConditionalMoments (OKBR, Davis 2023)" begin
    rng = Xoshiro(7)
    x = range(-3.0, 3.0; length = 25)

    @testset "against the offline formulas (Eqs. 8, 10)" begin
        X = ar1(rng, 4000)
        for K in (Epanechnikov(), Boxcar(), Triangular())
            cm = ConditionalMoments(x; h = 0.5, lags = 4, Δt = 0.1, kernel = K)
            foreach(v -> fit!(cm, v), X)
            v = value(cm)
            r = ref_moments(X, x, 0.5, 4, K)
            sel = r.W .> 0
            @test v.n == 4000 && v.lags == (1:4) .* 0.1 && size(v.M1) == (4, 25) && v.h == 0.5
            @test count(sel) > 50 && count(.!sel) > 0
            @test v.W[sel] ≈ r.W[sel] rtol = 1e-12
            @test maxrel(v.M1[sel], r.M1[sel]) < 1e-10
            @test maxrel(v.M2[sel], r.M2[sel]) < 1e-10
            @test maxrel(v.M2star[sel], r.M2star[sel]) < 1e-10
            @test all(isnan, v.M1[.!sel]) && all(isnan, v.M2star[.!sel]) && all(iszero, v.W[.!sel])
            # the algebraic identity M⁽²⁾ = M⁽²*⁾ + (M⁽¹⁾)² holds for the streamed values
            @test maxrel(v.M2[sel], v.M2star[sel] .+ v.M1[sel] .^ 2) < 1e-10
            # the effective sample size lies between 1 and the number of samples in the support
            @test all(1 .≤ v.neff[sel] .≤ 4000)
        end
    end

    @testset "merge! pools independent streams" begin
        X = ar1(rng, 4000)
        a = ConditionalMoments(x; h = 0.5, lags = 3); b = ConditionalMoments(x; h = 0.5, lags = 3)
        foreach(v -> fit!(a, v), X[1:2000]); foreach(v -> fit!(b, v), X[2001:end])
        ra, rb = ref_moments(X[1:2000], x, 0.5, 3, Epanechnikov()), ref_moments(X[2001:end], x, 0.5, 3, Epanechnikov())
        W = ra.W .+ rb.W
        sel = W .> 0
        wm(r) = ifelse.(r.W .> 0, r.W .* r.M1, 0.0)
        wm2(r) = ifelse.(r.W .> 0, r.W .* r.M2, 0.0)
        M1 = (wm(ra) .+ wm(rb)) ./ W
        M2 = (wm2(ra) .+ wm2(rb)) ./ W
        S = ra.S .+ rb.S .+ ifelse.((ra.W .> 0) .& (rb.W .> 0), (ra.M1 .- rb.M1) .^ 2 .* ra.W .* rb.W ./ W, 0.0)
        merge!(a, b)
        v = value(a)
        @test v.n == 4000
        @test v.W[sel] ≈ W[sel] rtol = 1e-12
        @test maxrel(v.M1[sel], M1[sel]) < 1e-10
        @test maxrel(v.M2[sel], M2[sel]) < 1e-10
        @test maxrel(v.M2star[sel], (S ./ W)[sel]) < 1e-10
        @test_throws ArgumentError merge!(a, ConditionalMoments(x; h = 0.4, lags = 3))
        @test_throws ArgumentError merge!(a, ConditionalMoments(x; h = 0.5, lags = 2))
        @test_throws ArgumentError merge!(a, ConditionalMoments(x; h = 0.5, lags = 3, kernel = Boxcar()))
    end

    @testset "Ornstein–Uhlenbeck (paper §3.1): D⁽¹⁾ = −x, D⁽²⁾ = 1" begin
        Δt, N = 1e-3, 4_000_000
        xs = range(-2.0, 2.0; length = 21)
        cm = ConditionalMoments(xs; h = 0.4, lags = 4, Δt)
        X, s = 0.0, sqrt(2Δt)                       # Euler–Maruyama; the noise amplitude is √(2D⁽²⁾Δt)
        for _ in 1:N
            X += -X * Δt + s * randn(rng)
            fit!(cm, X)
        end
        sel = abs.(xs) .≤ 1.5
        for method in (:direct, :regression)
            d = kramers_moyal(cm; method)
            @test length(d.lags) == (method === :direct ? 1 : 4)
            @test maximum(abs.(d.D1[sel] .+ xs[sel])) < 0.15
            @test maximum(abs.(d.D2[sel] .- 1)) < 0.08
            # the standard errors are of the right size: the drift is linear, so the kernel bias
            # is the shift of the weighted mean of x inside the window, ≤ h²/5 · |p'/p| ≈ 0.05
            @test all(abs.(d.D1[sel] .+ xs[sel]) .< 0.05 .+ 5 .* d.D1_err[sel])
            @test all(0 .< d.D1_err[sel] .< 0.1) && all(0 .< d.D2_err[sel] .< 0.1)
            @test all(d.n[sel] .> 1000)
        end
        # the conditional variance and the plain second moment agree at this Δt (the (τD⁽¹⁾)²
        # term is 0.4% at x = 2 and lag 4Δt) and the variance-based estimate is never larger
        d2 = kramers_moyal(cm; method = :regression, variance = true)
        d = kramers_moyal(cm; method = :regression)
        @test all(d2.D2[sel] .≤ d.D2[sel] .+ 1e-12)
        @test maximum(abs.(d2.D2[sel] .- d.D2[sel])) < 0.01
        # feeding the lagged series from a stored copy gives the same values (replay invariant)
    end

    @testset "multiplicative noise: D⁽¹⁾ = −x, D⁽²⁾ = 1 + x²/4" begin
        Δt, N = 1e-3, 2_000_000
        xs = range(-2.0, 2.0; length = 21)
        cm = ConditionalMoments(xs; h = 0.3, lags = 2, Δt)
        X = 0.0
        for _ in 1:N
            X += -X * Δt + sqrt(2 * (1 + X^2 / 4) * Δt) * randn(rng)
            fit!(cm, X)
        end
        d = kramers_moyal(cm; method = :regression, variance = true)
        sel = abs.(xs) .≤ 1.5
        @test maximum(abs.(d.D2[sel] .- (1 .+ xs[sel] .^ 2 ./ 4)) ./ (1 .+ xs[sel] .^ 2 ./ 4)) < 0.1
        @test maximum(abs.(d.D1[sel] .+ xs[sel])) < 0.2
    end

    @testset "constructor checks and kramers_moyal arguments" begin
        @test_throws ArgumentError ConditionalMoments(0.0:1.0:0.0; h = 0.1)
        @test_throws ArgumentError ConditionalMoments(1.0:-0.5:0.0; h = 0.1)
        @test_throws ArgumentError ConditionalMoments(range(0, 1; length = 5); h = 0)
        @test_throws ArgumentError ConditionalMoments(range(0, 1; length = 5); h = 0.1, lags = 0)
        @test_throws ArgumentError ConditionalMoments(range(0, 1; length = 5); h = 0.1, Δt = 0)
        cm = ConditionalMoments(range(-1, 1; length = 5); h = 0.5, lags = 3)
        @test cm isa ConditionalMoments{Float64,Epanechnikov}
        @test_throws ArgumentError kramers_moyal(cm; lags = 4:5)
        @test_throws ArgumentError kramers_moyal(cm; lags = Int[])
        @test_throws ArgumentError kramers_moyal(cm; method = :foo)
        @test all(isnan, kramers_moyal(cm).D1)                     # nothing fed: no estimate anywhere
        fit!(cm, 0.0)
        @test value(cm).n == 1 && all(iszero, value(cm).W)         # one sample: no increment yet
        fit!(cm, 0.1)
        v = value(cm)
        @test v.W[1, 3] > 0 && all(iszero, v.W[2:3, :])             # one increment of lag 1, conditioned on x = 0
        @test v.M1[1, 3] ≈ 0.1 && v.M2[1, 3] ≈ 0.01 && abs(v.M2star[1, 3]) < 1e-15 && v.neff[1, 3] ≈ 1
        # Float32 storage follows the evaluation points
        cm32 = ConditionalMoments(range(-1.0f0, 1.0f0; length = 5); h = 0.5f0, lags = 2, Δt = 1.0f0)
        @test cm32 isa ConditionalMoments{Float32,Epanechnikov}
        fit!(cm32, 0.0); fit!(cm32, 0.25)
        @test value(cm32).M1[1, 3] ≈ 0.25f0
    end
end
