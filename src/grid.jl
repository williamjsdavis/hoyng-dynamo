# Colatitude grid, quadrature weights and the spatial operators (dev/DECISIONS.md D2).

"""
    Grid(N; T = Float64)

Equidistant colatitude grid with `N` nodes on [0, π] including both poles,
θ_j = jπ/(N − 1), j = 0 … N − 1. The poles carry the boundary condition P = T = 0 and are
not unknowns; the `M = N − 2` interior nodes are. HOS's "60 equidistant θ points" is
`N = 60`; a grid of 60 interior points is `N = 62`.

All vector fields refer to the interior nodes. `w` holds the trapezoid weights of the
sin θ-weighted inner product of HOS Eq. 13, w_j = Δθ sin θ_j, and `W = [w; w]` the same
for a state vector b = [P; T].
"""
struct Grid{T<:AbstractFloat}
    N::Int
    M::Int
    Δθ::T
    θ::Vector{T}
    sinθ::Vector{T}
    cosθ::Vector{T}
    w::Vector{T}
    W::Vector{T}
end

function Grid(N::Integer; T::Type{<:AbstractFloat} = Float64)
    N ≥ 4 || throw(ArgumentError("need at least 4 nodes, got N = $N"))
    M = N - 2
    Δθ = T(π) / (N - 1)
    θ = [j * Δθ for j in 1:M]
    s = sin.(θ)
    w = Δθ .* s
    Grid{T}(N, M, Δθ, θ, s, cos.(θ), w, vcat(w, w))
end

Base.length(g::Grid) = g.M
Base.show(io::IO, g::Grid{T}) where {T} = print(io, "Grid{$T}(N = $(g.N), M = $(g.M) interior nodes)")

"""
    inner(g, x, y)

The weighted inner product of HOS Eq. 13 for two state vectors, (x, y) = Σ_k W_k conj(x_k) y_k,
i.e. ∫₀^π sin θ (P_x* P_y + T_x* T_y) dθ by the trapezoid rule.
"""
inner(g::Grid, x::AbstractVector, y::AbstractVector) = dot(x, Diagonal(g.W), y)

"""
    diffusion_operator(g, kR) -> Tridiagonal

Second-order conservative discretisation of L = (1/sin θ) ∂θ(sin θ ∂θ) − 1/sin²θ − (kR)²
(HOS Eq. 8) on the interior nodes with f = 0 at the poles:

    (L f)_j = [sin θ_{j+½}(f_{j+1} − f_j) − sin θ_{j−½}(f_j − f_{j−1})] / (sin θ_j Δθ²)
              − f_j / sin²θ_j − (kR)² f_j.

Self-adjoint under the sin θ-weighted inner product, like the continuous operator.
"""
function diffusion_operator(g::Grid{T}, kR) where {T}
    M, Δθ = g.M, g.Δθ
    sp = [sin(θ + Δθ / 2) for θ in g.θ]
    sm = [sin(θ - Δθ / 2) for θ in g.θ]
    k2 = T(kR)^2
    d = [-(sp[j] + sm[j]) / (g.sinθ[j] * Δθ^2) - 1 / g.sinθ[j]^2 - k2 for j in 1:M]
    du = [sp[j] / (g.sinθ[j] * Δθ^2) for j in 1:M-1]
    dl = [sm[j] / (g.sinθ[j] * Δθ^2) for j in 2:M]
    Tridiagonal(dl, d, du)
end

"""
    shear_operator(g) -> Tridiagonal

Central-difference discretisation of f ↦ ∂θ(sin θ f) with f = 0 at the poles,
(S f)_j = [sin θ_{j+1} f_{j+1} − sin θ_{j−1} f_{j−1}] / (2Δθ). The toroidal equation is
∂T/∂τ = C (S P) + L T (HOS Eq. 7).
"""
function shear_operator(g::Grid{T}) where {T}
    M, Δθ = g.M, g.Δθ
    du = [g.sinθ[j+1] / (2Δθ) for j in 1:M-1]
    dl = [-g.sinθ[j-1] / (2Δθ) for j in 2:M]
    Tridiagonal(dl, zeros(T, M), du)
end

"The α-effect operator at the critical state, f ↦ cos θ f (HOS Eq. 6 with α/α₀ = cos θ)."
alpha_operator(g::Grid) = Diagonal(g.cosθ)

"""
    dynamo_operator(g, kR, C) -> Matrix

Dense 2M×2M matrix of the linear dynamo operator G of HOS Eq. 6 acting on b = [P; T]:

    ∂P/∂τ = L P + cos θ T
    ∂T/∂τ = C ∂θ(sin θ P) + L T
"""
function dynamo_operator(g::Grid{T}, kR, C) where {T}
    M = g.M
    L = diffusion_operator(g, kR)
    G = zeros(T, 2M, 2M)
    G[1:M, 1:M] .= L
    G[1:M, M+1:2M] .= alpha_operator(g)
    G[M+1:2M, 1:M] .= T(C) .* shear_operator(g)
    G[M+1:2M, M+1:2M] .= L
    return G
end

"""
    quenching_operator(g) -> Matrix

Dense 2M×2M matrix of E, the operator multiplying s(1 − a²) in HOS Eq. 26: E b = [cos θ T; 0].
E reads only the toroidal rows and writes only the poloidal rows, so E² = 0.
"""
function quenching_operator(g::Grid{T}) where {T}
    M = g.M
    E = zeros(T, 2M, 2M)
    E[1:M, M+1:2M] .= alpha_operator(g)
    return E
end
