# Tier 0: linear operator and eigenproblem. Gates T0.1–T0.11 in dev/VALIDATION.md.
# Reference values: HOS Table I (p. 270), pp. 269, 272, 282, 284; HSO Fig. 8 caption (p. 156).

using LinearAlgebra

const TABLE_I_C0 = [0.0 + 0.0im, -2.76 + 4.41im, -8.67 + 4.97im, -15.9 + 6.13im]
const TABLE_I_CM = [0.527 + 0.0im, -2.16 + 5.88im, -8.21 + 6.04im, -15.7 + 7.70im]

relerr(x, ref) = abs(x - ref) / abs(ref)

# ℜ and ℑ each within 1% of the printed value; a printed zero must be |·| < 0.01 (real) or
# exactly zero (imaginary part of a steady mode).
function within_1pct(λ, ref)
    okre = iszero(real(ref)) ? abs(real(λ)) < 0.01 : relerr(real(λ), real(ref)) < 0.01
    okim = iszero(imag(ref)) ? imag(λ) == 0 : relerr(imag(λ), imag(ref)) < 0.01
    return okre && okim
end

@testset "Tier 0: linear operator and eigenproblem" begin
    grids = Dict(N => Grid(N) for N in (60, 120, 240))
    consts = Dict(N => model_constants(g, 0.5, 100.0) for (N, g) in grids)

    @testset "T0.1 critical dynamo number C₀ = 57.9" begin
        c240, _ = consts[240]
        @test abs(c240.C0 - 57.9) < 0.1
        @info "T0.1" C0_60 = consts[60][1].C0 C0_120 = consts[120][1].C0 C0_240 = c240.C0
    end

    for N in (60, 240)
        g = grids[N]
        c, m0 = consts[N]
        mm = eigenmodes(g, 0.5, 100.0)

        @testset "T0.2 Table I at C₀, N = $N" begin
            λ = leading_eigenvalues(m0, 4)
            @test abs(λ[1]) < 0.01
            for i in 2:4
                @test within_1pct(λ[i], TABLE_I_C0[i])
            end
            @info "T0.2 N = $N" λ
        end

        @testset "T0.3 Table I at C_m = 100, N = $N" begin
            λ = leading_eigenvalues(mm, 4)
            @test imag(λ[1]) == 0
            @test relerr(real(λ[1]), 0.527) < 0.01
            for i in 2:4
                @test within_1pct(λ[i], TABLE_I_CM[i])
            end
            @info "T0.3 N = $N" λ
        end

        @testset "T0.4 parity pattern, N = $N" begin
            λ = leading_eigenvalues(m0, 4)
            idx = [findfirst(==(z), m0.λ) for z in λ]
            @test toroidal_parity(g, view(m0.b, :, idx[1])) ≈ -1 atol = 1e-8
            @test poloidal_parity(g, view(m0.b, :, idx[1])) ≈ +1 atol = 1e-8
            @test toroidal_parity(g, view(m0.b, :, idx[2])) ≈ +1 atol = 1e-8
            @test toroidal_parity(g, view(m0.b, :, idx[3])) ≈ -1 atol = 1e-8
            @test toroidal_parity(g, view(m0.b, :, idx[4])) ≈ +1 atol = 1e-8
            @test imag(λ[1]) == 0
            @test all(imag(λ[i]) > 0 for i in 2:4)
        end

        @testset "T0.5 bi-orthonormality, N = $N" begin
            @test biorthogonality_residual(g, m0, 10) < 1e-10
            @test biorthogonality_residual(g, mm, 10) < 1e-10
        end

        @testset "T0.6 E₀₀ = 1.048, N = $N" begin
            @test c.E00 ≈ 1.048 atol = 0.005
        end
        @testset "T0.7 ∫ sin θ (P̃₀T₀)² dθ = 3.529, N = $N" begin
            @test c.D0int ≈ 3.529 atol = 0.01
        end
        @testset "T0.8 s, q₀, μ = 0.728, 0.421, 0.763, N = $N" begin
            @test c.s ≈ 0.728 atol = 0.002
            @test c.q0 ≈ 0.421 atol = 0.002
            @test c.μ ≈ 0.763 atol = 0.002
            @test c.Λ ≈ 0.527 atol = 0.002
        end

        @testset "T0.9 spectrum of G(C₀) + sE equals spectrum of G(C_m), N = $N" begin
            G0 = dynamo_operator(g, 0.5, c.C0) + c.s * quenching_operator(g)
            Gm = dynamo_operator(g, 0.5, 100.0)
            by = z -> (-real(z), -imag(z))
            λ0 = sort(eigvals(G0); by)
            λm = sort(eigvals(Gm); by)
            @test maximum(abs, λ0[1:10] - λm[1:10]) < 1e-8
        end
    end

    @testset "T0.10 second parameter set kR = 0.3, C_m = 120: Λ = 0.815, μ/Λ = 1.82" begin
        for N in (60, 240)
            c, _ = model_constants(grids[N], 0.3, 120.0)
            @test relerr(c.Λ, 0.815) < 0.01
            @test relerr(c.μ / c.Λ, 1.82) < 0.01
            @info "T0.10 N = $N" C0 = c.C0 s = c.s Λ = c.Λ E00 = c.E00 μ_over_Λ = c.μ / c.Λ
        end
    end

    @testset "T0.11 grid convergence (diagnostic)" begin
        rows = map((60, 120, 240)) do N
            c, m0 = consts[N]
            λ0 = leading_eigenvalues(m0, 4)
            λm = leading_eigenvalues(eigenmodes(grids[N], 0.5, 100.0), 4)
            (; N, c.C0, c.Λ, c.E00, c.D0int, λ1_C0 = λ0[2], λ3_C0 = λ0[4], λ1_Cm = λm[2], λ3_Cm = λm[4])
        end
        for r in rows
            @info "T0.11" r...
        end
        @test length(rows) == 3
    end
end
