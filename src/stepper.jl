# Constant-memory time steppers for the model equation (dev/DECISIONS.md D3, D6, D14).
#
#     ∂b/∂τ = [ G(C₀) + s(1 − a²) E + V(τ) ] b,     a = (b̃₀, b)          (HOS Eq. 26)
#
# G is constant. The explicit part B(τ) = s(1 − a²)E + V reads only the toroidal rows and
# writes only the poloidal rows, so B² = 0 and exp(Δτ B) = I + Δτ B exactly for frozen a
# and V. Freezing a over the step is, however, explicit Euler on the amplitude equation
# and shifts the relaxation rate about a = ±1 by a relative μΔτ (found by gate T1.4), so
# the quenching coefficient is evaluated at a midpoint prediction of a instead: since
# (b̃₀, G b) = λ⁰ a = 0 at the discrete C₀, da/dτ = (ã, B b) is one dot product away.

"""
    Model(; kR = 0.5, Cm = 100.0, N = 60)

Everything about the model equation that is fixed for a run: grid, `ModelConstants`
(discrete C₀, s, Λ, E₀₀, μ, …), the modes at C₀, the real fundamental mode `b0` and its
adjoint `b̃0`, the projection weights `ã = W .* b̃0` so that a = ã ⋅ b, and the dense
matrix `G` of G(C₀).
"""
struct Model{T<:AbstractFloat}
    grid::Grid{T}
    constants::ModelConstants{T}
    modes::Modes{T}
    b0::Vector{T}
    b̃0::Vector{T}
    ã::Vector{T}
    G::Matrix{T}
end

function Model(; kR = 0.5, Cm = 100.0, N::Integer = 60, T::Type{<:AbstractFloat} = Float64)
    g = Grid(N; T)
    c, m = model_constants(g, T(kR), T(Cm))
    b0, b̃0 = fundamental(m)
    Model{T}(g, c, m, b0, b̃0, g.W .* b̃0, dynamo_operator(g, c.kR, c.C0))
end

Base.show(io::IO, m::Model{T}) where {T} =
    print(io, "Model{$T}(kR = $(m.constants.kR), C_m = $(m.constants.Cm), N = $(m.grid.N); C₀ = $(m.constants.C0), s = $(m.constants.s), Λ = $(m.constants.Λ))")

"""
    State(model, b)

The time-stepped state: the field vector `b = [P; T]` (copied), the dipole amplitude
`a = (b̃₀, b)` kept current by `step!`, and the step counter `n` (τ = n Δτ).
"""
mutable struct State{T<:AbstractFloat}
    const b::Vector{T}
    a::T
    n::Int
end

State(model::Model{T}, b::AbstractVector) where {T} = State{T}(Vector{T}(b), dot(model.ã, b), 0)

"Workspace for `step!`: `tmp` (2M) and the current noise field `v` at the interior nodes (M), zero without noise."
struct Cache{T<:AbstractFloat}
    tmp::Vector{T}
    v::Vector{T}
end

Cache(model::Model{T}) where {T} = Cache{T}(zeros(T, 2 * model.grid.M), zeros(T, model.grid.M))

abstract type AbstractStepper{T<:AbstractFloat} end

"""
    ExponentialStepper(model, Δτ)

Production stepper (DECISIONS.md D3 option b):

    b_{n+1} = X (I + Δτ B_{n+½}) b_n,    X = exp(Δτ G) precomputed,  B = s(1 − a²)E + V,

with the quenching coefficient at a midpoint prediction of a (`quenching_coefficient`).
Exact for the linear dynamo at any Δτ, so b₀ is an exact fixed point and the spectrum
is the eigensolver's by construction. With constant G this Lie form is Strang splitting
in the half-step-staggered variable, and a is invariant under X (b̃₀ is a left null
vector of G), so rates measured from the trajectory are second order in Δτ (gate T1.6).
Cost per step: one dense 2M×2M matrix–vector product plus two dot products; 2.0 µs at N = 60.
"""
struct ExponentialStepper{T} <: AbstractStepper{T}
    Δτ::T
    s::T
    M::Int
    X::Matrix{T}
    cosθ::Vector{T}
    ã::Vector{T}
end

ExponentialStepper(model::Model{T}, Δτ) where {T} =
    ExponentialStepper{T}(T(Δτ), model.constants.s, model.grid.M, exp(T(Δτ) .* model.G), model.grid.cosθ, model.ã)

"""
    IMEXStepper(model, Δτ)

Cross-check stepper (DECISIONS.md D3 option a): backward Euler on the diffusion L of both
fields, forward Euler on the α-effect, shear, quenching and noise terms:

    (I − Δτ L) P_{n+1} = P_n + Δτ [ (1 + s(1 − a_n²)) cos θ + V_n ] T_n
    (I − Δτ L) T_{n+1} = T_n + Δτ C₀ ∂θ(sin θ P_n)

Two tridiagonal solves per step; first order in Δτ: the rates of gates T1.2 and T1.4 carry
errors of order ΛΔτ and a few μΔτ (b̃₀ is not a left eigenvector of the IMEX amplification
matrix, so overtones leak into a at O(Δτ)). 1.4 µs per step at N = 60.
"""
struct IMEXStepper{T,F} <: AbstractStepper{T}
    Δτ::T
    s::T
    M::Int
    F::F
    S::Tridiagonal{T,Vector{T}}
    cosθ::Vector{T}
    ã::Vector{T}
end

function IMEXStepper(model::Model{T}, Δτ) where {T}
    g, c = model.grid, model.constants
    L = diffusion_operator(g, c.kR)
    F = lu(I - T(Δτ) * L)
    S = c.C0 * shear_operator(g)
    IMEXStepper{T,typeof(F)}(T(Δτ), c.s, g.M, F, S, g.cosθ, model.ã)
end

"""
    ThetaStepper(model, Δτ; θ = 1, implicit_noise = false)
    BackwardEulerStepper(model, Δτ; implicit_noise = false)      θ = 1
    CrankNicolsonStepper(model, Δτ; implicit_noise = false)      θ = 1/2

Sensitivity steppers (dev/STRETCH_GOALS.md SG-1): the θ-method on the whole linear operator
A = G(C₀) + cE, c = s(1 − a²) at the midpoint prediction (`quenching_coefficient`), with the
noise term V explicit or, with `implicit_noise`, treated like A:

    [I − θΔτ A] b_{n+1} = [I + (1 − θ)Δτ A] b_n + Δτ V_n b_n                (explicit noise)
    [I − θΔτ (A + V_n)] b_{n+1} = [I + (1 − θ)Δτ (A + V_n)] b_n            (implicit noise)

θ = 1 is backward Euler, first order and over-damping (an eigenmode of A with eigenvalue λ
is multiplied by 1/(1 − Δτλ) per step instead of exp(Δτλ)); θ = 1/2 is Crank–Nicolson,
second order, with the stiffest diffusion modes mapped to ≈ −1 instead of 0. Both are
A-stable. They stand in for what an unstated implicit scheme of the papers' era does to the
overtone content at Δτ = 0.01; the exponential stepper remains the production one.

Solved by block elimination. With D = I − θΔτL (tridiagonal, fixed) and the α-coupling
diagonal g = (1 + c) cos θ [+ v], the toroidal unknown satisfies the dense M×M system

    [D − (θΔτ)² C₀ S D⁻¹ diag(g)] T_{n+1} = r_T + θΔτ C₀ S D⁻¹ r_P,

factorised in place every step (LAPACK getrf/getrs; g changes with a), after which
D P_{n+1} = r_P + θΔτ g ∘ T_{n+1} is one tridiagonal solve. About 25 µs per step at N = 60.
The stepper carries its own per-step workspace, so one instance serves one simulation.
"""
struct ThetaStepper{T,F} <: AbstractStepper{T}
    Δτ::T
    s::T
    M::Int
    θ::T
    implicit_noise::Bool
    L::Tridiagonal{T,Vector{T}}              # the diffusion operator L
    S::Tridiagonal{T,Vector{T}}              # C₀ ∂θ(sin θ ·), the shear term
    F::F                                     # lu(I − θΔτ L)
    Dm::Matrix{T}                            # I − θΔτ L, dense
    Q::Matrix{T}                             # (θΔτ)² C₀ S D⁻¹, dense
    K::Matrix{T}                             # workspace: Dm − Q diag(g), factorised in place
    ipiv::Vector{LinearAlgebra.BLAS.BlasInt} # workspace: its pivots
    g::Vector{T}                             # workspace: the α-coupling diagonal
    y::Vector{T}                             # workspace: D⁻¹ r_P
    cosθ::Vector{T}
    ã::Vector{T}
end

function ThetaStepper(model::Model{T}, Δτ; θ::Real = 1, implicit_noise::Bool = false) where {T}
    1 // 2 ≤ θ ≤ 1 || throw(ArgumentError("θ must lie in [1/2, 1], got $θ"))
    g, c = model.grid, model.constants
    M = g.M
    L = diffusion_operator(g, c.kR)
    S = c.C0 * shear_operator(g)
    α = T(θ) * T(Δτ)
    D = I - α * L
    Dm = Matrix(D)
    Q = α^2 .* (Matrix(S) * inv(Dm))
    F = lu(D)
    ThetaStepper{T,typeof(F)}(T(Δτ), c.s, M, T(θ), implicit_noise, L, S, F, Dm, Q, similar(Dm),
                              Vector{LinearAlgebra.BLAS.BlasInt}(undef, M), zeros(T, M), zeros(T, M), g.cosθ, model.ã)
end

BackwardEulerStepper(model::Model, Δτ; implicit_noise::Bool = false) = ThetaStepper(model, Δτ; θ = 1, implicit_noise)
CrankNicolsonStepper(model::Model, Δτ; implicit_noise::Bool = false) = ThetaStepper(model, Δτ; θ = 1 // 2, implicit_noise)

timestep(sp::AbstractStepper) = sp.Δτ

# s(1 − a²) at the midpoint of the step, from a_{n+½} = a_n + (Δτ/2) da/dτ with
# da/dτ = (ã, B b) = Σ_j ã_j [s(1 − a_n²) cos θ_j + v_j] T_j (the P rows of B b).
@inline function quenching_coefficient(sp::AbstractStepper{T}, b::AbstractVector{T}, v::AbstractVector{T}, a::T) where {T}
    M = sp.M
    c = sp.s * (1 - a^2)
    da = zero(T)
    @inbounds @simd for j in 1:M
        da += sp.ã[j] * (c * sp.cosθ[j] + v[j]) * b[M+j]
    end
    ah = a + (sp.Δτ / 2) * da
    return sp.s * (1 - ah^2)
end

"""
    step!(state, stepper, noise, cache)

Advance `state` by one step Δτ. Allocation-free; `state.a` is current on return.
"""
function step!(st::State{T}, sp::ExponentialStepper{T}, noise, cache::Cache{T}) where {T}
    M, Δτ = sp.M, sp.Δτ
    b, tmp, v = st.b, cache.tmp, cache.v
    update_noise!(v, noise, st.n)
    c = quenching_coefficient(sp, b, v, st.a)
    @inbounds @simd for j in 1:M
        tmp[j] = b[j] + Δτ * (c * sp.cosθ[j] + v[j]) * b[M+j]
    end
    @inbounds @simd for j in M+1:2M
        tmp[j] = b[j]
    end
    mul!(b, sp.X, tmp)
    st.a = dot(sp.ã, b)
    st.n += 1
    return st
end

function step!(st::State{T}, sp::IMEXStepper{T}, noise, cache::Cache{T}) where {T}
    M, Δτ = sp.M, sp.Δτ
    b, tmp, v = st.b, cache.tmp, cache.v
    update_noise!(v, noise, st.n)
    c = quenching_coefficient(sp, b, v, st.a)
    dl, du = sp.S.dl, sp.S.du
    @inbounds begin
        @simd for j in 1:M
            tmp[j] = b[j] + Δτ * ((1 + c) * sp.cosθ[j] + v[j]) * b[M+j]
        end
        tmp[M+1] = b[M+1] + Δτ * du[1] * b[2]
        for j in 2:M-1
            tmp[M+j] = b[M+j] + Δτ * (dl[j-1] * b[j-1] + du[j] * b[j+1])
        end
        tmp[2M] = b[2M] + Δτ * dl[M-1] * b[M-1]
    end
    ldiv!(sp.F, view(tmp, 1:M))
    ldiv!(sp.F, view(tmp, M+1:2M))
    copyto!(b, tmp)
    st.a = dot(sp.ã, b)
    st.n += 1
    return st
end

# y[oy+1:oy+M] += β (A x[ox+1:ox+M]) for a tridiagonal A of size M ≥ 2.
@inline function add_tridiagonal!(y::AbstractVector{T}, oy::Int, β::T, A::Tridiagonal{T}, x::AbstractVector{T}, ox::Int, M::Int) where {T}
    dl, d, du = A.dl, A.d, A.du
    @inbounds begin
        y[oy+1] += β * (d[1] * x[ox+1] + du[1] * x[ox+2])
        for j in 2:M-1
            y[oy+j] += β * (dl[j-1] * x[ox+j-1] + d[j] * x[ox+j] + du[j] * x[ox+j+1])
        end
        y[oy+M] += β * (dl[M-1] * x[ox+M-1] + d[M] * x[ox+M])
    end
    return y
end

function step!(st::State{T}, sp::ThetaStepper{T}, noise, cache::Cache{T}) where {T}
    M, Δτ, θ = sp.M, sp.Δτ, sp.θ
    b, tmp, v = st.b, cache.tmp, cache.v
    g, y, K = sp.g, sp.y, sp.K
    update_noise!(v, noise, st.n)
    c = quenching_coefficient(sp, b, v, st.a)
    α = θ * Δτ                                   # implicit weight
    β = (1 - θ) * Δτ                             # explicit weight of A
    ν = sp.implicit_noise ? β : Δτ               # explicit weight of V
    ι = sp.implicit_noise ? one(T) : zero(T)     # V inside the implicit operator
    # r = b + β A b + ν V b into tmp: the α-coupling of the P rows here, the tridiagonal
    # parts (L on both fields, the shear) below; g is the implicit α-coupling diagonal
    @inbounds @simd for j in 1:M
        tmp[j] = b[j] + (β * (1 + c) * sp.cosθ[j] + ν * v[j]) * b[M+j]
        tmp[M+j] = b[M+j]
        g[j] = (1 + c) * sp.cosθ[j] + ι * v[j]
    end
    if θ != 1
        add_tridiagonal!(tmp, 0, β, sp.L, b, 0, M)
        add_tridiagonal!(tmp, M, β, sp.S, b, 0, M)
        add_tridiagonal!(tmp, M, β, sp.L, b, M, M)
    end
    # toroidal system K T⁺ = r_T + α C₀ S D⁻¹ r_P with K = D − Q diag(g)
    copyto!(y, 1, tmp, 1, M)
    ldiv!(sp.F, y)
    add_tridiagonal!(tmp, M, α, sp.S, y, 0, M)
    @inbounds for j in 1:M
        gj = g[j]
        @simd for i in 1:M
            K[i, j] = sp.Dm[i, j] - gj * sp.Q[i, j]
        end
    end
    LinearAlgebra.LAPACK.getrf!(K, sp.ipiv)
    LinearAlgebra.LAPACK.getrs!('N', K, sp.ipiv, view(tmp, M+1:2M))
    # poloidal system D P⁺ = r_P + α g ∘ T⁺
    @inbounds @simd for j in 1:M
        tmp[j] += α * g[j] * tmp[M+j]
    end
    ldiv!(sp.F, view(tmp, 1:M))
    copyto!(b, tmp)
    st.a = dot(sp.ã, b)
    st.n += 1
    return st
end

"The current dipole amplitude a = (b̃₀, b)."
dipole_amplitude(st::State) = st.a

"‖b − a b₀‖₂: the part of the field not in the fundamental mode."
overtone_residual(model::Model, st::State) = sqrt(sum(abs2, st.b .- st.a .* model.b0))

"""
    amplification_matrix(model, stepper) -> Matrix

The one-step map b ↦ b_{n+1} of `stepper` for the linear dynamo at C₀ (a = ±1 frozen, so
c = 0, and no noise): exp(ΔτG) for the exponential stepper, (I − ΔτL̂)⁻¹(I + Δτ(G − L̂)) with
L̂ = diag(L, L) for the IMEX, (I − θΔτG)⁻¹(I + (1 − θ)ΔτG) for the θ-method. Its eigenvalues
against exp(Δτλ_k) give each scheme's per-step bias on the modes (dev/STRETCH_GOALS.md SG-1).
"""
amplification_matrix(::Model, sp::ExponentialStepper) = copy(sp.X)

function amplification_matrix(model::Model{T}, sp::IMEXStepper{T}) where {T}
    M = model.grid.M
    L = diffusion_operator(model.grid, model.constants.kR)
    Lbd = zeros(T, 2M, 2M)
    Lbd[1:M, 1:M] .= L
    Lbd[M+1:2M, M+1:2M] .= L
    return (I - sp.Δτ * Lbd) \ (I + sp.Δτ * (model.G - Lbd))
end

amplification_matrix(model::Model{T}, sp::ThetaStepper{T}) where {T} =
    (I - sp.θ * sp.Δτ * model.G) \ (I + (1 - sp.θ) * sp.Δτ * model.G)

"""
    mode_factor(stepper, λ)

The factor the scheme applies per step to an eigenmode of the linear operator with
eigenvalue λ, for the schemes where it is a function of λ alone: exp(Δτλ) for the
exponential stepper, (1 + (1 − θ)Δτλ)/(1 − θΔτλ) for the θ-method.
"""
mode_factor(sp::ExponentialStepper, λ::Number) = exp(sp.Δτ * λ)
mode_factor(sp::ThetaStepper, λ::Number) = (1 + (1 - sp.θ) * sp.Δτ * λ) / (1 - sp.θ * sp.Δτ * λ)
