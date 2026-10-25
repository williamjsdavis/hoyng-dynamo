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
# The alternative quenching laws of STRETCH_GOALS.md SG-2 replace s(1 − a²) by node-wise
# coefficients c_j (`AbstractQuenching`, `quenching!`); B² = 0 holds for all of them.

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

"""
Workspace for `step!`: `tmp` (2M), the current noise field `v` at the interior nodes (M,
zero without noise), the node-wise quenching coefficients `c` (M) and the midpoint
toroidal field `th` (M) of the field-based quenching laws.
"""
struct Cache{T<:AbstractFloat}
    tmp::Vector{T}
    v::Vector{T}
    c::Vector{T}
    th::Vector{T}
end

Cache(model::Model{T}) where {T} = (M = model.grid.M; Cache{T}(zeros(T, 2M), zeros(T, M), zeros(T, M), zeros(T, M)))

abstract type AbstractStepper{T<:AbstractFloat} end

# ---------------------------------------------------------------------------------------
# Quenching laws (DECISIONS.md S2, STRETCH_GOALS.md SG-2)
# ---------------------------------------------------------------------------------------

"""
    AbstractQuenching

The nonlinearity of HOS Eq. 26. The α-coefficient of the poloidal equation is (1 + c) cos θ
with c = (C_m/C₀) h(q) − 1, so that c = 0 at the nonlinear equilibrium and c = s for a
vanishing field; `quenching!(cache, law, stepper, b, a)` fills `cache.c` with the node-wise
coefficients c_j at the midpoint of the step. HOS p. 277 report that the alternatives leave
the dipole statistics "not materially changed" (gate T4.16). The laws:

- `AmplitudeQuenching()`, the papers' law (HOS Eqs. 17, 19): h = 1 − q with q = q₀a², hence
  c = s(1 − a²). The midpoint is a_{n+½} = a_n + (Δτ/2)(ã, B_n b_n) (D3).
- `EnergyQuenching()`: h = 1 − q with q from the toroidal energy (HOS Eq. 18) instead of the
  fundamental mode alone, q = q₀ E_T, E_T = Σ_j w_j T_j² (= a² for a pure fundamental mode
  in the unit of `eigenmodes`), hence c = s(1 − E_T); the equilibrium is still b₀.
- `InverseQuenching()`: h = 1/(1 + q) in place of 1 − q, with q = s a² so that the
  equilibrium is again a = ±1, hence c = s(1 − a²)/(1 + s a²). The relaxation rate about
  a = ±1 is 2μ/(1 + s) instead of 2μ.
- `LocalQuenching(κ)`: Rüdiger's local law (HOS p. 270), α(θ)(1 − κ T(θ)²), hence
  c_j = s − (1 + s) κ T_j². κ only sets the field unit (HOS p. 271); `calibrate` chooses it
  so that the nonlinear equilibrium has a = 1, as for the global laws, although its shape
  is no longer the fundamental mode.

The field-based laws take the midpoint field T_{n+½} = T_n + (Δτ/2) ∂T/∂τ with
∂T/∂τ = C₀ ∂θ(sin θ P) + L T, exact because the nonlinearity and the noise act on P only.

With h = 1 − q the α-effect reverses sign once q > 1, i.e. c < −1, and the α-reversed αΩ
dynamo is supercritical for |C_eff| > C₀ (mirror symmetry of G under C → −C). The papers'
law reaches that point only at |a| > √(1 + 1/s) = 1.54, briefly, and ignores the overtones;
the energy and local laws see every overtone burst, and once their quenching measure has
reversed α the growing field drives it further: the runs diverge within a few hundred τ_d at
f = 6.4 (`VALIDATION.md` T4.16, `DISCREPANCIES.md` DR-8). `floored = true` clamps the
coefficient at c ≥ −1 (the α-effect is suppressed to zero, never reversed), which leaves the
equilibria unchanged and makes those laws stable; the inverse law needs no floor (h > 0).
"""
abstract type AbstractQuenching end
struct AmplitudeQuenching <: AbstractQuenching
    floored::Bool
end
struct EnergyQuenching <: AbstractQuenching
    floored::Bool
end
struct InverseQuenching <: AbstractQuenching end
struct LocalQuenching{T<:AbstractFloat} <: AbstractQuenching
    κ::T
    floored::Bool
end

AmplitudeQuenching(; floored::Bool = false) = AmplitudeQuenching(floored)
EnergyQuenching(; floored::Bool = false) = EnergyQuenching(floored)
LocalQuenching(κ::T; floored::Bool = false) where {T<:AbstractFloat} = LocalQuenching{T}(κ, floored)

"Whether the law clamps its coefficient at c ≥ −1 (α suppressed, never reversed)."
floored(law::AbstractQuenching) = law.floored
floored(::InverseQuenching) = false
@inline clamp_coefficient(c::T, law::AbstractQuenching) where {T} = floored(law) ? max(c, -one(T)) : c

Base.show(io::IO, q::LocalQuenching) = print(io, "LocalQuenching(κ = ", q.κ, q.floored ? ", floored)" : ")")
Base.show(io::IO, q::Union{AmplitudeQuenching,EnergyQuenching}) = print(io, nameof(typeof(q)), "(", q.floored ? "floored" : "", ")")

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

# the midpoint toroidal field T + (Δτ/2)(C₀ S P + L T) into th
@inline function midpoint_toroidal!(th::AbstractVector{T}, sp::AbstractStepper{T}, b::AbstractVector{T}) where {T}
    M, h = sp.M, sp.Δτ / 2
    @inbounds @simd for j in 1:M
        th[j] = b[M+j]
    end
    add_tridiagonal!(th, 0, h, sp.S, b, 0, M)
    add_tridiagonal!(th, 0, h, sp.L, b, M, M)
    return th
end

@inline global_coefficient(::AmplitudeQuenching, s::T, a::T) where {T} = s * (1 - a^2)
@inline global_coefficient(::InverseQuenching, s::T, a::T) where {T} = s * (1 - a^2) / (1 + s * a^2)

"""
    quenching!(cache, law, stepper, b, a) -> cache.c

Fill `cache.c` with the node-wise quenching coefficients of `law` at the midpoint of the
step from state `b`, `a` with the current noise `cache.v` (see `AbstractQuenching`).
"""
@inline function quenching!(cache::Cache{T}, law::Union{AmplitudeQuenching,InverseQuenching}, sp::AbstractStepper{T},
                            b::AbstractVector{T}, a::T) where {T}
    M, v = sp.M, cache.v
    cn = clamp_coefficient(global_coefficient(law, sp.s, a), law)
    da = zero(T)
    @inbounds @simd for j in 1:M
        da += sp.ã[j] * (cn * sp.cosθ[j] + v[j]) * b[M+j]
    end
    ah = a + (sp.Δτ / 2) * da
    return fill!(cache.c, clamp_coefficient(global_coefficient(law, sp.s, ah), law))
end

@inline function quenching!(cache::Cache{T}, law::EnergyQuenching, sp::AbstractStepper{T}, b::AbstractVector{T}, a::T) where {T}
    th = midpoint_toroidal!(cache.th, sp, b)
    E = zero(T)
    @inbounds @simd for j in 1:sp.M
        E += sp.w[j] * th[j]^2
    end
    return fill!(cache.c, clamp_coefficient(sp.s * (1 - E), law))
end

@inline function quenching!(cache::Cache{T}, law::LocalQuenching{T}, sp::AbstractStepper{T}, b::AbstractVector{T}, a::T) where {T}
    th = midpoint_toroidal!(cache.th, sp, b)
    c, s, κ = cache.c, sp.s, law.κ
    @inbounds @simd for j in 1:sp.M
        c[j] = clamp_coefficient(s - (1 + s) * κ * th[j]^2, law)
    end
    return c
end

# ---------------------------------------------------------------------------------------
# Steppers
# ---------------------------------------------------------------------------------------

"""
    ExponentialStepper(model, Δτ; quenching = AmplitudeQuenching())

Production stepper (DECISIONS.md D3 option b):

    b_{n+1} = X (I + Δτ B_{n+½}) b_n,    X = exp(Δτ G) precomputed,  B = diag(c) E + V,

with the quenching coefficients at the midpoint of the step (`quenching!`). Exact for the
linear dynamo at any Δτ, so b₀ is an exact fixed point and the spectrum is the
eigensolver's by construction. With constant G this Lie form is Strang splitting in the
half-step-staggered variable, and a is invariant under X (b̃₀ is a left null vector of G),
so rates measured from the trajectory are second order in Δτ (gate T1.6). Cost per step:
one dense 2M×2M matrix–vector product plus two dot products; 2.0 µs at N = 60.
"""
struct ExponentialStepper{T,Q<:AbstractQuenching} <: AbstractStepper{T}
    Δτ::T
    s::T
    M::Int
    q::Q
    X::Matrix{T}
    L::Tridiagonal{T,Vector{T}}      # the diffusion operator L (for the midpoint toroidal field)
    S::Tridiagonal{T,Vector{T}}      # C₀ ∂θ(sin θ ·), the shear term (likewise)
    w::Vector{T}
    cosθ::Vector{T}
    ã::Vector{T}
end

function ExponentialStepper(model::Model{T}, Δτ; quenching::AbstractQuenching = AmplitudeQuenching()) where {T}
    g, c = model.grid, model.constants
    ExponentialStepper{T,typeof(quenching)}(T(Δτ), c.s, g.M, quenching, exp(T(Δτ) .* model.G), diffusion_operator(g, c.kR),
                                            c.C0 * shear_operator(g), g.w, g.cosθ, model.ã)
end

"""
    IMEXStepper(model, Δτ; quenching = AmplitudeQuenching())

Cross-check stepper (DECISIONS.md D3 option a): backward Euler on the diffusion L of both
fields, forward Euler on the α-effect, shear, quenching and noise terms:

    (I − Δτ L) P_{n+1} = P_n + Δτ [ (1 + c_j) cos θ + V_n ] T_n
    (I − Δτ L) T_{n+1} = T_n + Δτ C₀ ∂θ(sin θ P_n)

Two tridiagonal solves per step; first order in Δτ: the rates of gates T1.2 and T1.4 carry
errors of order ΛΔτ and a few μΔτ (b̃₀ is not a left eigenvector of the IMEX amplification
matrix, so overtones leak into a at O(Δτ)). 1.4 µs per step at N = 60.
"""
struct IMEXStepper{T,F,Q<:AbstractQuenching} <: AbstractStepper{T}
    Δτ::T
    s::T
    M::Int
    q::Q
    F::F
    L::Tridiagonal{T,Vector{T}}
    S::Tridiagonal{T,Vector{T}}
    w::Vector{T}
    cosθ::Vector{T}
    ã::Vector{T}
end

function IMEXStepper(model::Model{T}, Δτ; quenching::AbstractQuenching = AmplitudeQuenching()) where {T}
    g, c = model.grid, model.constants
    L = diffusion_operator(g, c.kR)
    F = lu(I - T(Δτ) * L)
    S = c.C0 * shear_operator(g)
    IMEXStepper{T,typeof(F),typeof(quenching)}(T(Δτ), c.s, g.M, quenching, F, L, S, g.w, g.cosθ, model.ã)
end

"""
    ThetaStepper(model, Δτ; θ = 1, implicit_noise = false, quenching = AmplitudeQuenching())
    BackwardEulerStepper(model, Δτ; implicit_noise = false, quenching)      θ = 1
    CrankNicolsonStepper(model, Δτ; implicit_noise = false, quenching)      θ = 1/2

Sensitivity steppers (dev/STRETCH_GOALS.md SG-1): the θ-method on the whole linear operator
A = G(C₀) + diag(c)E with the quenching coefficients at the midpoint (`quenching!`), and
the noise term V explicit or, with `implicit_noise`, treated like A:

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
D P_{n+1} = r_P + θΔτ g ∘ T_{n+1} is one tridiagonal solve. About 18 µs per step at N = 60.
The stepper carries its own per-step workspace, so one instance serves one simulation.
"""
struct ThetaStepper{T,F,Q<:AbstractQuenching} <: AbstractStepper{T}
    Δτ::T
    s::T
    M::Int
    θ::T
    implicit_noise::Bool
    q::Q
    L::Tridiagonal{T,Vector{T}}              # the diffusion operator L
    S::Tridiagonal{T,Vector{T}}              # C₀ ∂θ(sin θ ·), the shear term
    F::F                                     # lu(I − θΔτ L)
    Dm::Matrix{T}                            # I − θΔτ L, dense
    Q::Matrix{T}                             # (θΔτ)² C₀ S D⁻¹, dense
    K::Matrix{T}                             # workspace: Dm − Q diag(g), factorised in place
    ipiv::Vector{LinearAlgebra.BLAS.BlasInt} # workspace: its pivots
    g::Vector{T}                             # workspace: the α-coupling diagonal
    y::Vector{T}                             # workspace: D⁻¹ r_P
    w::Vector{T}
    cosθ::Vector{T}
    ã::Vector{T}
end

function ThetaStepper(model::Model{T}, Δτ; θ::Real = 1, implicit_noise::Bool = false,
                      quenching::AbstractQuenching = AmplitudeQuenching()) where {T}
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
    ThetaStepper{T,typeof(F),typeof(quenching)}(T(Δτ), c.s, M, T(θ), implicit_noise, quenching, L, S, F, Dm, Q, similar(Dm),
                                                Vector{LinearAlgebra.BLAS.BlasInt}(undef, M), zeros(T, M), zeros(T, M),
                                                g.w, g.cosθ, model.ã)
end

BackwardEulerStepper(model::Model, Δτ; implicit_noise::Bool = false, quenching::AbstractQuenching = AmplitudeQuenching()) =
    ThetaStepper(model, Δτ; θ = 1, implicit_noise, quenching)
CrankNicolsonStepper(model::Model, Δτ; implicit_noise::Bool = false, quenching::AbstractQuenching = AmplitudeQuenching()) =
    ThetaStepper(model, Δτ; θ = 1 // 2, implicit_noise, quenching)

timestep(sp::AbstractStepper) = sp.Δτ
"The quenching law of a stepper."
quenching(sp::AbstractStepper) = sp.q

"""
    step!(state, stepper, noise, cache)

Advance `state` by one step Δτ. Allocation-free; `state.a` is current on return.
"""
function step!(st::State{T}, sp::ExponentialStepper{T}, noise, cache::Cache{T}) where {T}
    M, Δτ = sp.M, sp.Δτ
    b, tmp, v = st.b, cache.tmp, cache.v
    update_noise!(v, noise, st.n)
    c = quenching!(cache, sp.q, sp, b, st.a)
    @inbounds @simd for j in 1:M
        tmp[j] = b[j] + Δτ * (c[j] * sp.cosθ[j] + v[j]) * b[M+j]
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
    c = quenching!(cache, sp.q, sp, b, st.a)
    dl, du = sp.S.dl, sp.S.du
    @inbounds begin
        @simd for j in 1:M
            tmp[j] = b[j] + Δτ * ((1 + c[j]) * sp.cosθ[j] + v[j]) * b[M+j]
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

function step!(st::State{T}, sp::ThetaStepper{T}, noise, cache::Cache{T}) where {T}
    M, Δτ, θ = sp.M, sp.Δτ, sp.θ
    b, tmp, v = st.b, cache.tmp, cache.v
    g, y, K = sp.g, sp.y, sp.K
    update_noise!(v, noise, st.n)
    c = quenching!(cache, sp.q, sp, b, st.a)
    α = θ * Δτ                                   # implicit weight
    β = (1 - θ) * Δτ                             # explicit weight of A
    ν = sp.implicit_noise ? β : Δτ               # explicit weight of V
    ι = sp.implicit_noise ? one(T) : zero(T)     # V inside the implicit operator
    # r = b + β A b + ν V b into tmp: the α-coupling of the P rows here, the tridiagonal
    # parts (L on both fields, the shear) below; g is the implicit α-coupling diagonal
    @inbounds @simd for j in 1:M
        tmp[j] = b[j] + (β * (1 + c[j]) * sp.cosθ[j] + ν * v[j]) * b[M+j]
        tmp[M+j] = b[M+j]
        g[j] = (1 + c[j]) * sp.cosθ[j] + ι * v[j]
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
    calibrate(LocalQuenching, model, Δτ; κ0 = 1, τ = 100, tol = 1e-9, floored = false) -> LocalQuenching

The local law's constant such that the nonlinear equilibrium has a = 1. κ only sets the
field unit (HOS p. 271): if b* is the steady state for κ then b*/√λ is the steady state for
κλ, so κ = κ0 a*(κ0)². The steady state is reached by integrating the deterministic model
from b₀ with the exponential stepper for `τ` diffusion times; an error is thrown if a has
not settled to `tol` over the following 10 τ_d or is not positive (the equilibrium would be
oscillatory, decayed or reversed).
"""
function calibrate(::Type{LocalQuenching}, model::Model{T}, Δτ; κ0::Real = 1, τ::Real = 100, tol::Real = 1e-9, floored::Bool = false) where {T}
    sp = ExponentialStepper(model, Δτ; quenching = LocalQuenching(T(κ0)))
    st, cache = State(model, model.b0), Cache(model)
    for _ in 1:round(Int, τ / Δτ)
        step!(st, sp, NoNoise(), cache)
    end
    a1 = st.a
    for _ in 1:round(Int, 10 / Δτ)
        step!(st, sp, NoNoise(), cache)
    end
    a2 = st.a
    (a2 > 0 && abs(a2 - a1) ≤ tol * a2) || throw(ArgumentError("local quenching with κ0 = $κ0: a has not settled, $a1 → $a2 over 10 τ_d"))
    return LocalQuenching(T(κ0) * a2^2; floored)
end

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
