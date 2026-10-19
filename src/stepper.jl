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

"The current dipole amplitude a = (b̃₀, b)."
dipole_amplitude(st::State) = st.a

"‖b − a b₀‖₂: the part of the field not in the fundamental mode."
overtone_residual(model::Model, st::State) = sqrt(sum(abs2, st.b .- st.a .* model.b0))
