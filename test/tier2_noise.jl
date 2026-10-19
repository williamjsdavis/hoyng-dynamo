# Tier 2: noise construction. Gates T2.1–T2.5 in dev/VALIDATION.md.
# References (HOS notation): Eq. 21 (p. 272), Appendix B Eqs. B.6–B.8 (pp. 308–309), Eq. 36 (D₀),
# Appendix C Eq. C.3 (p. 310). ℱ ≡ πf²τ_c/(2N_c²) = 3.22×10⁻² for the standard run (p. 291).
# The record/stream-mode and aggregator-replay parts of T2.5 belong to build step 6 (run.jl).

using LinearAlgebra, Random

mean(x) = sum(x) / length(x)
variance(x) = (m = mean(x); sum(abs2, x .- m) / (length(x) - 1))
function correlation(x, y)
    dx, dy = x .- mean(x), y .- mean(y)
    return dot(dx, dy) / sqrt(sum(abs2, dx) * sum(abs2, dy))
end

const T2_N, T2_Δτ, T2_f, T2_τc, T2_Nc = 60, 0.01, 6.4, 0.05, 10

@testset "Tier 2: noise construction" begin
    model = Model(; kR = 0.5, Cm = 100.0, N = T2_N)
    g = model.grid
    M = g.M

    @testset "construction" begin
        nz = BinnedNoise(model, T2_Δτ; f = T2_f, seed = 1)
        @test redraw_interval(nz) == 5
        @test nz.bin == [floor(Int, θ * T2_Nc / π) + 1 for θ in g.θ]      # gap G3: bin ⌊θ_j N_c/π⌋
        @test extrema(nz.bin) == (1, T2_Nc)
        @test nz.bin == reverse(T2_Nc + 1 .- nz.bin)                       # mirror-symmetric bins
        @test nz.g ≈ 1 ./ sqrt.(2 * T2_Nc .* g.sinθ)
        @test nz.f == T2_f
        @test_throws ArgumentError BinnedNoise(model, 0.03)                # τ_c/Δτ not an integer
        @test_throws ArgumentError LinearRamp(4.0, 6.5, 10.0, 5.0)
    end

    @testset "T2.1 F sample statistics" begin
        ndraw = 10^5                                                       # redraws: 10⁶ values in 10 bins
        for dist in (GaussianF(), UniformF(), BinaryF())
            nz = BinnedNoise(model, T2_Δτ; f = T2_f, distribution = dist, seed = 7)
            Fs = Matrix{Float64}(undef, T2_Nc, ndraw)
            for m in 1:ndraw
                redraw!(nz)
                Fs[:, m] .= nz.F
            end
            x = vec(Fs)
            # 4σ: three distributions at one fixed seed; at seed 7 the binary F has a 3.4σ mean
            # over 10⁶ draws while 10⁹ draws over 20 seeds are consistent with zero (VALIDATION.md §4)
            @test abs(mean(x)) < 4 / sqrt(length(x))
            @test abs(variance(x) - 1) < 0.01                             # 1 %
            # independence across bins (45 pairs) and successive intervals (10 bins): σ = 1/√ndraw
            tol = 4.5 / sqrt(ndraw)
            @test maximum(abs(correlation(Fs[k, :], Fs[l, :])) for k in 1:T2_Nc for l in k+1:T2_Nc) < tol
            @test maximum(abs(correlation(Fs[k, 1:end-1], Fs[k, 2:end])) for k in 1:T2_Nc) < tol
            if dist isa UniformF
                @test maximum(abs, x) ≤ sqrt(3)
            elseif dist isa BinaryF
                @test all(==(1), abs.(x))
            end
        end
        # constant for n_c = 5 steps, redrawn at n mod 5 = 0, nonzero everywhere
        nz = BinnedNoise(model, T2_Δτ; f = T2_f, seed = 3)
        v = zeros(M)
        vs = [copy(update_noise!(v, nz, n)) for n in 0:19]
        @test all(!=(0), vs[1])
        for n in 1:19
            if n % 5 == 0
                @test vs[n+1] != vs[n]
            else
                @test vs[n+1] == vs[n]
            end
        end
        # nodes in one bin share F: v_j/(f g_j) is constant within each bin
        Fnode = vs[1] ./ (T2_f .* nz.g)
        for k in 1:T2_Nc
            idx = findall(==(k), nz.bin)
            @test all(Fnode[idx] .≈ Fnode[idx[1]])
        end
    end

    @testset "T2.2 correlation function of δα/α₀ (Eq. B.6)" begin
        nz = BinnedNoise(model, T2_Δτ; f = T2_f, seed = 11)
        nc = redraw_interval(nz)
        nsteps = 10^6
        nint = nsteps ÷ nc
        v = zeros(M)
        C = zeros(M, M)                                  # node covariance, one sample per interval
        jA, jB = 5, 29                                   # near the pole and near the equator
        sA, sB = Vector{Float64}(undef, nsteps), Vector{Float64}(undef, nsteps)
        for n in 0:nsteps-1
            update_noise!(v, nz, n)
            n % nc == 0 && BLAS.syr!('U', 1.0, v, C)    # upper triangle of C += v vᵀ
            sA[n+1] = v[jA]
            sB[n+1] = v[jB]
        end
        C ./= nint
        # zero lag: ⟨v_j²⟩ = f²/(2 N_c sin θ_j), within 2 %
        var_theory = T2_f^2 ./ (2 * T2_Nc .* g.sinθ)
        @test maximum(abs.(diag(C) ./ var_theory .- 1)) < 0.02
        # angle: the normalised covariance of nodes d apart, averaged over j, equals the fraction of
        # node pairs in the same bin (exact for the binned model) and approximates the triangle
        # 1 − |ζ|/θ_c of Eq. B.6 (the bins hold 5 or 6 nodes)
        θc = π / T2_Nc
        for d in 0:12
            js = 1:M-d
            measured = mean([C[j, j+d] / (T2_f^2 * nz.g[j] * nz.g[j+d]) for j in js])
            exact = mean([nz.bin[j] == nz.bin[j+d] for j in js])
            triangle = max(1 - d * g.Δθ / θc, 0.0)
            @test abs(measured - exact) < 0.01
            @test abs(exact - triangle) < 0.03
        end
        # time: autocorrelation at lag l steps is 1 − l/n_c for l ≤ n_c and 0 beyond
        for (s, j) in ((sA, jA), (sB, jB))
            c0 = sum(abs2, s) / nsteps
            @test abs(c0 / var_theory[j] - 1) < 0.02
            for l in 0:2nc
                cl = dot(view(s, 1:nsteps-l), view(s, 1+l:nsteps)) / (nsteps - l)
                @test abs(cl / c0 - max(1 - l / nc, 0)) < 0.01
            end
        end
    end

    @testset "T2.3–T2.4 D₀ measured directly" begin
        nz = BinnedNoise(model, T2_Δτ; f = T2_f, seed = 5)
        Tseg, nseg = 10.0, 10_000                        # segments of 10 τ_d = 200 τ_c
        nper = round(Int, Tseg / T2_Δτ)
        v = zeros(M)
        I = Vector{Float64}(undef, nseg)                 # ∫₀ᵀ V₀₀ dτ per segment
        n = 0
        for m in 1:nseg
            acc = 0.0
            for _ in 1:nper
                update_noise!(v, nz, n)
                acc += noise_projection(model, v)
                n += 1
            end
            I[m] = acc * T2_Δτ
        end
        D0_measured = variance(I) / Tseg
        D0_binned = binned_D0(model, nz)
        ℱ = π * T2_f^2 * T2_τc / (2 * T2_Nc^2)
        D0_formula = model.constants.D0int * ℱ           # Eq. 36 with the δ-correlated B.7
        model240 = Model(; kR = 0.5, Cm = 100.0, N = 240)
        D0_binned240 = binned_D0(model240, BinnedNoise(model240, T2_Δτ; f = T2_f))
        @test abs(mean(I)) < 4 * sqrt(variance(I) / nseg)                 # zero-mean increments
        @test abs(ℱ / 3.22e-2 - 1) < 0.01                                 # p. 291
        @test abs(D0_formula / 0.114 - 1) < 0.01                          # Eq. 36, p. 284
        @test abs(D0_measured / D0_binned - 1) < 0.05                     # 3.5σ at 10⁴ segments
        @test abs(D0_measured / D0_formula - 1) < 0.10                    # T2.3: finite bins vs B.7
        @test D0_binned < D0_formula                                      # Jensen: binning lowers D₀
        # T2.4: (D₀τ_c)^{1/2} = 0.075 (C.3): 1 % from the formula, 10 % from the measurement
        @test abs(sqrt(D0_formula * T2_τc) / 0.075 - 1) < 0.01
        @test abs(sqrt(D0_measured * T2_τc) / 0.075 - 1) < 0.10
        @info "T2.3–T2.4 (f = 6.4, N = 60, 10⁴ segments of 10 τ_d)" ℱ D0_formula D0_binned D0_binned240 D0_measured sqrt_D0τc_formula = sqrt(D0_formula * T2_τc) sqrt_D0τc_measured = sqrt(D0_measured * T2_τc)
    end

    @testset "T2.5 reproducibility" begin
        sp = ExponentialStepper(model, T2_Δτ)
        function run_a(seed, f; nsteps = 10_000)
            nz = BinnedNoise(model, T2_Δτ; f, seed)
            st, cache = State(model, model.b0), Cache(model)
            a = Vector{Float64}(undef, nsteps)
            for k in 1:nsteps
                step!(st, sp, nz, cache)
                a[k] = st.a
            end
            return a
        end
        a1, a2, a3 = run_a(42, T2_f), run_a(42, T2_f), run_a(43, T2_f)
        @test a1 == a2                                                     # bitwise
        @test a1 != a3
        @test all(isfinite, a1)
        # one noise series across f values (HOS Fig. 5): the F stream does not depend on f
        nzA = BinnedNoise(model, T2_Δτ; f = 5.4, seed = 9)
        nzB = BinnedNoise(model, T2_Δτ; f = 7.4, seed = 9)
        vA, vB = zeros(M), zeros(M)
        for n in 0:99
            update_noise!(vA, nzA, n)
            update_noise!(vB, nzB, n)
            @test nzA.F == nzB.F
            @test vB ≈ (7.4 / 5.4) .* vA rtol = 1e-14
        end
        # per-run seeds: a fixed hash of (master, index), distinct across runs and masters
        seeds = [run_seed(2026, i) for i in 0:999]
        @test seeds == [run_seed(2026, i) for i in 0:999]
        @test length(unique(seeds)) == 1000
        @test run_seed(2026, 0) != run_seed(2027, 0)
        # ramped forcing is evaluated at each redraw
        ramp = LinearRamp(4.0, 6.5, 0.0, 1.0)
        @test forcing(ramp, -1.0) == 4.0
        @test forcing(ramp, 0.5) ≈ 5.25
        @test forcing(ramp, 2.0) == 6.5
        nzR = BinnedNoise(model, T2_Δτ; schedule = ramp, seed = 1)
        vR = zeros(M)
        for n in 0:200
            update_noise!(vR, nzR, n)
            n % 5 == 0 && @test nzR.f ≈ forcing(ramp, n * T2_Δτ)
        end
        @test nzR.f == 6.5
    end
end
