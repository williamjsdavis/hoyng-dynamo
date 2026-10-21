# Eigenmodes of G, their adjoints, the critical dynamo number and the derived constants
# (dev/DECISIONS.md D4; HOS §2.2–2.3, Table I).

"""
    Modes{T}

Eigenmodes of the discrete operator G(C) at `kR`, `C`. `λ` holds all 2M eigenvalues
ordered by decreasing real part, the ℑλ > 0 member of a conjugate pair first; `b[:, i]`
is the right eigenvector (mode) and `b̃[:, i]` the adjoint mode, normalised so that
`inner(g, b̃[:, i], b[:, j]) = δ_ij` (HOS Eqs. 11–13). The fundamental mode `i = 1` is
real, with unit toroidal energy Σ w T₀² = 1 and ∫ sin θ P₀ dθ > 0, which defines a = +1
as "normal" polarity; this scale only sets the unit of the field (HOS §11).
"""
struct Modes{T<:AbstractFloat}
    kR::T
    C::T
    λ::Vector{Complex{T}}
    b::Matrix{Complex{T}}
    b̃::Matrix{Complex{T}}
end

Base.show(io::IO, m::Modes{T}) where {T} =
    print(io, "Modes{$T}(kR = $(m.kR), C = $(m.C), λ⁰ = $(m.λ[1]), $(length(m.λ)) modes)")

"""
    eigenmodes(g, kR, C) -> Modes

Dense eigen-decomposition of G(C) (LAPACK), adjoint modes from the left eigenvectors
rescaled by the quadrature weights: if Y = V⁻¹ then b̃_i = conj(Y[i, :]) ./ W, so that
(b̃_i, b_j) = Y V = I to round-off for any discretisation.
"""
function eigenmodes(g::Grid{T}, kR, C) where {T}
    G = dynamo_operator(g, kR, C)
    F = eigen(G)
    λ = Vector{Complex{T}}(F.values)
    V = Matrix{Complex{T}}(F.vectors)
    p = sortperm(λ; by = z -> (-real(z), -imag(z)))
    λ = λ[p]
    V = V[:, p]
    Y = inv(V)
    B̃ = similar(V)
    @inbounds for i in axes(V, 2), k in axes(V, 1)
        B̃[k, i] = conj(Y[i, k]) / g.W[k]
    end
    m = Modes{T}(T(kR), T(C), λ, V, B̃)
    normalise_modes!(m, g)
    return m
end

# Every mode is scaled to unit toroidal energy ∫|T_k|² sin θ dθ = 1 (b̃_k inversely), so that
# the overtone amplitudes a_k = (b̃_k, b) are in units of the equilibrium dipole's toroidal
# r.m.s. field, and the field unit is the one in which the nonlinear equilibrium b₀ has
# unit toroidal energy. HOS fix only (b̃_i, b_j) = δ_ij and leave the unit of B to the
# constant of Eq. 18, so their |a_k| and energy scales (Figs. 2, 6) are not ours. A real
# steady fundamental mode is also made real with positive mean P ("normal" polarity).
function normalise_modes!(m::Modes{T}, g::Grid{T}) where {T}
    M = g.M
    for i in axes(m.b, 2)
        v = @view m.b[:, i]
        energy = sum(g.w[j] * abs2(v[M+j]) for j in 1:M)
        energy > 0 || continue
        c = Complex{T}(inv(sqrt(energy)))
        if i == 1 && imag(m.λ[1]) == 0 && maximum(abs ∘ imag, v) ≤ 1e-12 * maximum(abs, v)
            c *= sign(sum(g.w[j] * real(v[j]) for j in 1:M))
        end
        m.b[:, i] .*= c
        m.b̃[:, i] ./= c
    end
    return m
end

"Real fundamental mode and its adjoint, `(b₀, b̃₀)`, as `Vector{T}`."
function fundamental(m::Modes{T}) where {T}
    imag(m.λ[1]) == 0 || error("the leading eigenvalue $(m.λ[1]) is not real")
    return real.(m.b[:, 1]), real.(m.b̃[:, 1])
end

"""
    leading_eigenvalues(m, n) -> Vector{Complex}

The first `n` distinct modes in Table I's convention: one member per conjugate pair, ℑλ ≥ 0.
"""
function leading_eigenvalues(m::Modes{T}, n::Integer) where {T}
    out = Complex{T}[]
    for λ in m.λ
        imag(λ) < 0 && continue
        push!(out, λ)
        length(out) == n && break
    end
    return out
end

"""
    biorthogonality_residual(g, m, n) -> real

max |(b̃_i, b_j) − δ_ij| over the first `n` modes (HOS Eq. 12).
"""
function biorthogonality_residual(g::Grid, m::Modes, n::Integer)
    r = zero(real(eltype(m.b)))
    for i in 1:n, j in 1:n
        r = max(r, abs(inner(g, view(m.b̃, :, i), view(m.b, :, j)) - (i == j)))
    end
    return r
end

# Parity of a state vector's toroidal or poloidal part about the equator (HOS Eq. 28):
# 𝒫 = (E_S − E_A)/(E_S + E_A), −1 antisymmetric, +1 symmetric. Interior node j mirrors to M+1−j.
function _parity(g::Grid, f::AbstractVector)
    M = g.M
    ES = zero(real(eltype(f)))
    EA = zero(ES)
    @inbounds for j in 1:M
        fj, fm = f[j], f[M+1-j]
        ES += g.w[j] * abs2(fj + fm) / 4
        EA += g.w[j] * abs2(fj - fm) / 4
    end
    return (ES - EA) / (ES + EA)
end

"Parity 𝒫 of the toroidal field T of a state vector b = [P; T]: −1 dipole-like, +1 quadrupole-like."
toroidal_parity(g::Grid, b::AbstractVector) = _parity(g, view(b, g.M+1:2g.M))
"Parity of the poloidal potential P of a state vector b = [P; T]."
poloidal_parity(g::Grid, b::AbstractVector) = _parity(g, view(b, 1:g.M))

"Largest real part of the spectrum of G(C): the growth rate of the leading mode."
growth_rate(g::Grid, kR, C) = maximum(real, eigvals(dynamo_operator(g, kR, C)))

"""
    critical_dynamo_number(g, kR; bracket = (20, 100), tol = 1e-10) -> C₀

The root of ℜλ⁰(C) = 0, λ⁰ the eigenvalue of G(C) with the largest real part (HOS p. 269:
57.9 at kR = 0.5). Illinois-type regula falsi inside `bracket`; each evaluation is one
dense eigenvalue problem. Use the returned discrete value, not 57.9, so that b₀ is an
exact fixed point of the discrete model (DECISIONS.md D4).
"""
function critical_dynamo_number(g::Grid{T}, kR; bracket = (20, 100), tol = 1e-10, maxiter = 200) where {T}
    a, b = T(bracket[1]), T(bracket[2])
    fa, fb = growth_rate(g, kR, a), growth_rate(g, kR, b)
    fa < 0 < fb || throw(ArgumentError("bracket $(bracket) does not enclose C₀: ℜλ⁰ = $fa, $fb"))
    side = 0
    for _ in 1:maxiter
        c = (a * fb - b * fa) / (fb - fa)
        fc = growth_rate(g, kR, c)
        if fc == 0 || abs(b - a) < tol
            return c
        elseif fc > 0
            b, fb = c, fc
            side == 1 && (fa /= 2)
            side = 1
        else
            a, fa = c, fc
            side == -1 && (fb /= 2)
            side = -1
        end
    end
    error("critical_dynamo_number did not converge in $maxiter iterations")
end

"E₀₀ = ∫ sin θ cos θ P̃₀ T₀ dθ (HOS p. 282; 1.048 at kR = 0.5)."
function overlap_E00(g::Grid, m::Modes)
    b0, b̃0 = fundamental(m)
    M = g.M
    return sum(g.w[j] * g.cosθ[j] * b̃0[j] * b0[M+j] for j in 1:M)
end

"∫ sin θ (P̃₀ T₀)² dθ, the integral in the direct diffusion coefficient D₀ (HOS p. 284; 3.529 at kR = 0.5)."
function diffusion_integral(g::Grid, m::Modes)
    b0, b̃0 = fundamental(m)
    M = g.M
    return sum(g.w[j] * (b̃0[j] * b0[M+j])^2 for j in 1:M)
end

"""
    ModelConstants

The constants of the model equation for one (kR, C_m) on one grid (HOS §2.4, Eq. 34; DEVELOPMENT.md §1.1):
C₀ the discrete critical dynamo number, s = C_m/C₀ − 1, q₀ = s/(1 + s), Λ = λ⁰(C_m) the
small-amplitude growth rate, E₀₀ the overlap integral, μ = sE₀₀, and D₀int = ∫ sin θ (P̃₀T₀)² dθ.
"""
struct ModelConstants{T<:AbstractFloat}
    kR::T
    Cm::T
    C0::T
    s::T
    q0::T
    Λ::T
    E00::T
    μ::T
    D0int::T
end

"Compute `ModelConstants` and the modes at C₀ for grid `g`: returns `(constants, modes)`."
function model_constants(g::Grid{T}, kR, Cm) where {T}
    C0 = critical_dynamo_number(g, kR)
    m = eigenmodes(g, kR, C0)
    s = T(Cm) / C0 - 1
    E00 = overlap_E00(g, m)
    c = ModelConstants{T}(T(kR), T(Cm), C0, s, s / (1 + s), growth_rate(g, kR, Cm), E00, s * E00, diffusion_integral(g, m))
    return c, m
end
