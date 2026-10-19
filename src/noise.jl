# The stochastic α (HOS §2.4, Eq. 21; dev/DECISIONS.md D5, D6; DEVELOPMENT.md gaps G2, G3).
#
#     δα/α₀ = f F(θ, τ) / √(2 N_c sin θ)
#
# F has zero mean and unit variance, is constant within each of the N_c latitude bins of
# width θ_c = π/N_c and within each interval τ_c, and is redrawn independently in every
# bin at the end of each interval. At the interior nodes the noise field is
#
#     v_j = f g_j F_{k_j},     g_j = 1/√(2 N_c sin θ_j),     k_j = ⌊θ_j N_c/π⌋,
#
# redrawn at the steps n with n mod n_c = 0, n_c = τ_c/Δτ (an integer). F is drawn
# independently of f and multiplied by f at application, so one seed gives one noise
# series that can be reused across forcing values (HOS p. 280, gate T2.5).
#
# Resuming from a checkpoint must start at a redraw step (n mod n_c = 0), because the
# current F is not part of the checkpoint.

"The deterministic model: V ≡ 0."
struct NoNoise end

"""
    update_noise!(v, noise, n)

Bring the noise field `v[j] = δα/α₀ at interior node j` up to date for step `n` (the step
from τ = nΔτ to (n+1)Δτ). Allocation-free. Nothing to do for `NoNoise`.
"""
@inline update_noise!(v::AbstractVector, ::NoNoise, n::Integer) = v

# ---------------------------------------------------------------- distributions of F (G2, Q7)

"Distribution of the unit-variance random function F; `GaussianF` is the default (`DECISIONS.md` D5)."
abstract type AbstractFDistribution end
"F ~ N(0, 1)."
struct GaussianF <: AbstractFDistribution end
"F uniform on [−√3, √3]."
struct UniformF <: AbstractFDistribution end
"F = ±1 with equal probability."
struct BinaryF <: AbstractFDistribution end

@inline draw(rng::AbstractRNG, ::GaussianF, ::Type{T}) where {T} = randn(rng, T)
@inline draw(rng::AbstractRNG, ::UniformF, ::Type{T}) where {T} = (2 * rand(rng, T) - 1) * sqrt(T(3))
@inline draw(rng::AbstractRNG, ::BinaryF, ::Type{T}) where {T} = rand(rng, Bool) ? one(T) : -one(T)

# ---------------------------------------------------------------- forcing schedules (D6)

"A forcing schedule f(τ), evaluated at each redraw and constant within a τ_c interval."
abstract type AbstractForcing{T<:AbstractFloat} end

"Constant forcing parameter f."
struct ConstantForcing{T} <: AbstractForcing{T}
    f::T
end

"""
    LinearRamp(f0, f1, τ0, τ1)

f(τ) = f0 for τ ≤ τ0, linear in τ to f1 at τ1, then f1 (SOH's ramped run, Fig. 3).
"""
struct LinearRamp{T} <: AbstractForcing{T}
    f0::T
    f1::T
    τ0::T
    τ1::T
    function LinearRamp{T}(f0, f1, τ0, τ1) where {T}
        τ1 > τ0 || throw(ArgumentError("need τ1 > τ0, got τ0 = $τ0, τ1 = $τ1"))
        new{T}(f0, f1, τ0, τ1)
    end
end
LinearRamp(f0, f1, τ0, τ1) = LinearRamp{float(promote_type(typeof(f0), typeof(f1), typeof(τ0), typeof(τ1)))}(f0, f1, τ0, τ1)

"The forcing parameter f at time τ."
@inline forcing(s::ConstantForcing, τ) = s.f
@inline function forcing(s::LinearRamp{T}, τ) where {T}
    x = clamp((T(τ) - s.τ0) / (s.τ1 - s.τ0), zero(T), one(T))
    return s.f0 + (s.f1 - s.f0) * x
end

# ---------------------------------------------------------------- seeds (D5)

# SplitMix64 finaliser: a fixed, version-independent 64-bit mixing function.
@inline function splitmix64(x::UInt64)
    z = x + 0x9e3779b97f4a7c15
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return z ⊻ (z >> 31)
end

"""
    run_seed(master, index)

The seed of run `index` of an experiment with the given master seed: a fixed hash of the
pair, so each run's noise series is independent of thread scheduling and run order.
"""
run_seed(master::Integer, index::Integer) = splitmix64(splitmix64(UInt64(master)) ⊻ UInt64(index))

# ---------------------------------------------------------------- the binned noise

"""
    BinnedNoise(model_or_grid, Δτ; f = 6.4, schedule = ConstantForcing(f), τc = 0.05, Nc = 10,
                distribution = GaussianF(), seed = 0)

The stochastic α of HOS Eq. 21 on the grid: F constant within each of `Nc` latitude bins and
each interval `τc`, redrawn (bins 0 … Nc−1 in order) at the steps with n mod n_c = 0,
n_c = τc/Δτ, which must be an integer. `seed` seeds a Xoshiro256** generator; the F stream
depends only on (seed, distribution, Nc, n_c), not on f. `schedule` gives f(τ), evaluated at
each redraw.
"""
mutable struct BinnedNoise{T<:AbstractFloat,D<:AbstractFDistribution,S<:AbstractForcing{T}}
    const Nc::Int
    const nc::Int
    const Δτ::T
    const τc::T
    const bin::Vector{Int}      # k_j + 1 ∈ 1:Nc
    const g::Vector{T}          # 1/√(2 N_c sin θ_j)
    const F::Vector{T}          # current F per bin
    const rng::Xoshiro
    const distribution::D
    const schedule::S
    const seed::UInt64
    f::T                        # f at the last redraw
end

function BinnedNoise(g::Grid{T}, Δτ; f = 6.4, schedule::AbstractForcing = ConstantForcing(T(f)),
                     τc = 0.05, Nc::Integer = 10, distribution::AbstractFDistribution = GaussianF(),
                     seed::Integer = 0) where {T}
    Nc ≥ 1 || throw(ArgumentError("need N_c ≥ 1, got $Nc"))
    nc = round(Int, τc / Δτ)
    (nc ≥ 1 && abs(nc * Δτ - τc) ≤ 1e-9 * τc) || throw(ArgumentError("τ_c/Δτ = $(τc / Δτ) must be a positive integer"))
    bin = [min(floor(Int, θ * Nc / π), Nc - 1) + 1 for θ in g.θ]
    gfac = [1 / sqrt(2 * Nc * s) for s in g.sinθ]
    S = typeof(schedule)
    S <: AbstractForcing{T} || throw(ArgumentError("forcing schedule must have element type $T"))
    BinnedNoise{T,typeof(distribution),S}(Nc, nc, T(Δτ), T(τc), bin, gfac, zeros(T, Nc), Xoshiro(seed),
                                          distribution, schedule, UInt64(seed), forcing(schedule, zero(T)))
end
BinnedNoise(model::Model, Δτ; kwargs...) = BinnedNoise(model.grid, Δτ; kwargs...)

Base.show(io::IO, nz::BinnedNoise{T,D,S}) where {T,D,S} =
    print(io, "BinnedNoise{$T}(N_c = $(nz.Nc), τ_c = $(nz.τc), n_c = $(nz.nc), $(D), $(nz.schedule), seed = $(nz.seed))")

"Redraw F in every bin, bins 0 … N_c−1 in order. Allocation-free."
@inline function redraw!(nz::BinnedNoise{T}) where {T}
    F, rng, d = nz.F, nz.rng, nz.distribution
    @inbounds for k in 1:nz.Nc
        F[k] = draw(rng, d, T)
    end
    return nz
end

"Steps between redraws, n_c = τ_c/Δτ."
redraw_interval(nz::BinnedNoise) = nz.nc

@inline function update_noise!(v::AbstractVector{T}, nz::BinnedNoise{T}, n::Integer) where {T}
    if n % nz.nc == 0
        redraw!(nz)
        f = forcing(nz.schedule, T(n) * nz.Δτ)
        nz.f = f
        F, g, bin = nz.F, nz.g, nz.bin
        @inbounds for j in eachindex(v)
            v[j] = f * (g[j] * F[bin[j]])
        end
    end
    return v
end

# ---------------------------------------------------------------- diagnostics

"""
    noise_projection(model, v)

V₀₀ = (b̃₀, V b₀) = ∫ sin θ (δα/α₀) P̃₀ T₀ dθ (HOS Eq. 31 with i = k = 0) for the noise field `v`.
"""
function noise_projection(model::Model{T}, v::AbstractVector{T}) where {T}
    M, ã, b0 = model.grid.M, model.ã, model.b0
    s = zero(T)
    @inbounds @simd for j in 1:M
        s += ã[j] * v[j] * b0[M+j]
    end
    return s
end

"""
    binned_D0(model, noise; f = noise.f)

The exact diffusion coefficient D₀ = 2∫₀^∞⟨V₀₀(τ)V₀₀(τ−s)⟩ds of the binned noise on the frozen
mode b₀: V₀₀ is constant over each interval τ_c and independent between intervals, so
D₀ = τ_c Var[V₀₀] = τ_c f² Σ_k (Σ_{j∈k} w_j P̃₀_j g_j T₀_j)². Compare with HOS Eq. 36, which uses the
δ-correlated limit B.7 of the correlation function (gate T2.3).
"""
function binned_D0(model::Model{T}, nz::BinnedNoise{T}; f = nz.f) where {T}
    M, ã, b0 = model.grid.M, model.ã, model.b0
    acc = zeros(T, nz.Nc)
    for j in 1:M
        acc[nz.bin[j]] += ã[j] * nz.g[j] * b0[M+j]
    end
    return nz.τc * T(f)^2 * sum(abs2, acc)
end
