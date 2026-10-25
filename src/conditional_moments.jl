# Online kernel-based regression (OKBR) for drift and diffusion functions: W. Davis (2023),
# "Reconstruction of stochastic dynamics from large streamed datasets", arXiv:2307.00445v3,
# Eqs. 11–15. The conditional moments of a scalar series X_n sampled every Δt,
#
#     M⁽ᵏ⁾(τ_i, x_j) = Σ_n K_h(x_j − X_n) [X_{n+i} − X_n]ᵏ / Σ_n K_h(x_j − X_n),    k = 1, 2,
#
# at lags τ_i = iΔt, i = 1 … N_τ, and evaluation points x_j, with the conditional variance
# M⁽²*⁾ = Σ K_h ([X_{n+i} − X_n] − M⁽¹⁾)² / Σ K_h, are accumulated in O(N_τ N_x) memory
# from a stream, one sample at a time, by the weighted Welford/West updates of the paper
# (the increment of lag i is conditioned on its starting point X_{N−i}, the Itô convention).
# The Kramers–Moyal coefficients follow from the moments offline (`kramers_moyal`).
#
# Beyond the paper's formulas (its definitions are kept): only the evaluation points inside
# the kernel's support are visited for each increment, so the cost per sample is
# O(N_τ h/Δx) instead of O(N_τ N_x); Σ K_h² is accumulated beside Σ K_h so that an effective
# sample size and standard errors are available; `merge!` pools independent streams (the
# parallel evaluation the paper's discussion mentions, after Chan et al. 1982); and the
# diffusion can be built from the conditional variance, which removes the (τD⁽¹⁾)² term of
# the finite-τ second moment. Behind the fit!/merge!/value interface of DECISIONS.md D8,
# allocation-free (D14).

"""
    AbstractKernel

A kernel K(u) of unit support (K(u) = 0 for |u| ≥ 1), callable; `K_h(u) = K(u/h)/h`.
`Epanechnikov` (the paper's default, Eq. 9), `Boxcar` (bin counting: histogram-based
regression) and `Triangular`.
"""
abstract type AbstractKernel end
struct Epanechnikov <: AbstractKernel end
struct Boxcar <: AbstractKernel end
struct Triangular <: AbstractKernel end

@inline (::Epanechnikov)(u::T) where {T<:AbstractFloat} = abs(u) < 1 ? T(0.75) * (1 - u * u) : zero(T)
@inline (::Boxcar)(u::T) where {T<:AbstractFloat} = abs(u) < 1 ? T(0.5) : zero(T)
@inline (::Triangular)(u::T) where {T<:AbstractFloat} = abs(u) < 1 ? 1 - abs(u) : zero(T)

"""
    log_lags(nmax; per_decade = 10) -> Vector{Int}

Distinct sample lags from 1 to `nmax`, about `per_decade` per decade of lag: dense where the
correlation structure of a series lives (the noise correlation time, the relaxation times)
and sparse where the moments have saturated. `log_lags(200)` is 22 lags, 1 … 200.
"""
function log_lags(nmax::Integer; per_decade::Integer = 10)
    nmax ≥ 1 || throw(ArgumentError("need nmax ≥ 1"))
    per_decade ≥ 1 || throw(ArgumentError("need per_decade ≥ 1"))
    n = ceil(Int, log10(nmax) * per_decade) + 1
    return unique!(max.(1, round.(Int, exp10.(range(0, log10(nmax); length = n)))))
end

"""
    ConditionalMoments(x::AbstractRange; h, lags = 1, Δt = 1.0, kernel = Epanechnikov())

Online kernel-based regression (Davis 2023): the conditional moments M⁽¹⁾, M⁽²⁾ and the
conditional variance M⁽²*⁾ of the increments of a streamed series at the equally spaced
evaluation points `x`, at the sample lags `lags` (an integer n for the lags 1 … n, or an
increasing vector of distinct lags in samples, e.g. `log_lags(200)`; `Δt` per sample), with
the kernel `K_h(u) = K(u/h)/h` of bandwidth `h`. The lags must straddle the correlation
times of the series: too sparse and the correlation structure is missed, too fine and the
moments do not change over the set. `fit!(cm, X)` feeds one sample; `merge!` pools
independent streams; `value` gives `(x, lags, W, W2, M1, M2, M2star, neff, n, h)` with the
lags in time units and the matrices N_τ × N_x as in the paper (`W` the cumulative kernel
weights, `W2` the sum of their squares, `neff = W²/W2` the effective number of samples; NaN
where a point has no weight). Memory O(N_τ N_x) plus a ring buffer of the longest lag; cost
per sample O(N_τ h/Δx).
"""
mutable struct ConditionalMoments{T<:AbstractFloat,K<:AbstractKernel} <: Aggregator
    const x::Vector{T}
    const x0::T
    const invdx::T
    const h::T
    const invh::T
    const kernel::K
    const lags::Vector{Int}     # the sample lags, increasing
    const L::Int                # the longest lag: length of the ring buffer
    const Δt::T
    const buf::Vector{T}        # ring of the last L samples
    const W::Matrix{T}          # N_x × N_τ internally (contiguous over x for a given lag)
    const W2::Matrix{T}
    const M1::Matrix{T}
    const M2::Matrix{T}
    const S::Matrix{T}          # Σ K_h (Δ − M⁽¹⁾)², the paper's S
    n::Int
    pos::Int                    # ring index where the next sample is stored
end

function ConditionalMoments(x::AbstractRange; h::Real, lags = 1, Δt::Real = 1.0, kernel::AbstractKernel = Epanechnikov(),
                            T::Type{<:AbstractFloat} = float(promote_type(eltype(x), typeof(h), typeof(Δt))))
    length(x) ≥ 2 || throw(ArgumentError("need at least two evaluation points"))
    dx = step(x)
    dx > 0 || throw(ArgumentError("the evaluation points must increase"))
    h > 0 || throw(ArgumentError("need a positive bandwidth h"))
    Δt > 0 || throw(ArgumentError("need Δt > 0"))
    lagidx = lags isa Integer ? collect(1:Int(lags)) : collect(Int, lags)
    (!isempty(lagidx) && lagidx[1] ≥ 1 && issorted(lagidx) && allunique(lagidx)) ||
        throw(ArgumentError("lags must be a positive integer n (the lags 1:n) or an increasing vector of distinct positive sample lags"))
    nx, nl, L = length(x), length(lagidx), lagidx[end]
    z() = zeros(T, nx, nl)
    ConditionalMoments{T,typeof(kernel)}(collect(T, x), T(first(x)), 1 / T(dx), T(h), 1 / T(h), kernel, lagidx, L, T(Δt),
                                         zeros(T, L), z(), z(), z(), z(), z(), 0, 1)
end

# One increment Δ = X_N − X_{N−i} of lag i, conditioned on its starting point Xi = X_{N−i}:
# the paper's Eqs. 11, 12 and 15 at every evaluation point inside the kernel's support
# (one point of slack on either side; the kernel is zero outside, so nothing is lost).
@inline function update_row!(cm::ConditionalMoments{T}, i::Int, Xi::T, Δ::T) where {T}
    x, W, W2, M1, M2, S = cm.x, cm.W, cm.W2, cm.M1, cm.M2, cm.S
    jlo = max(1, floor(Int, (Xi - cm.h - cm.x0) * cm.invdx))
    jhi = min(length(x), ceil(Int, (Xi + cm.h - cm.x0) * cm.invdx) + 2)
    Δ2 = Δ * Δ
    @inbounds for j in jlo:jhi
        k = cm.kernel((x[j] - Xi) * cm.invh) * cm.invh
        k > 0 || continue
        w = W[j, i] + k
        W[j, i] = w
        W2[j, i] += k * k
        m_old = M1[j, i]
        m_new = m_old + k * (Δ - m_old) / w
        M1[j, i] = m_new
        S[j, i] += k * (Δ - m_old) * (Δ - m_new)
        M2[j, i] += k * (Δ2 - M2[j, i]) / w
    end
    return nothing
end

@inline function fit!(cm::ConditionalMoments{T}, xnew::Real) where {T}
    X = T(xnew)
    cm.n += 1
    L, buf, lags, pos = cm.L, cm.buf, cm.lags, cm.pos
    @inbounds for k in eachindex(lags)
        ℓ = lags[k]
        ℓ ≤ cm.n - 1 || break            # increasing lags: the later ones have no sample yet either
        j = pos - ℓ
        j ≤ 0 && (j += L)
        Xi = buf[j]                      # X_{N−ℓ}
        update_row!(cm, k, Xi, X - Xi)
    end
    @inbounds buf[pos] = X
    cm.pos = pos == L ? 1 : pos + 1
    return cm
end

function Base.merge!(a::ConditionalMoments{T}, b::ConditionalMoments{T}) where {T}
    (a.x == b.x && a.h == b.h && a.lags == b.lags && a.Δt == b.Δt && a.kernel == b.kernel) ||
        throw(ArgumentError("conditional moments have different evaluation points, bandwidth, lags, Δt or kernel"))
    W, W2, M1, M2, S = a.W, a.W2, a.M1, a.M2, a.S
    @inbounds for idx in eachindex(W)
        Wa, Wb = W[idx], b.W[idx]
        w = Wa + Wb
        w > 0 || continue
        δ = b.M1[idx] - M1[idx]
        M1[idx] += δ * (Wb / w)
        S[idx] += b.S[idx] + δ * δ * (Wa * Wb / w)
        M2[idx] += (b.M2[idx] - M2[idx]) * (Wb / w)
        W[idx] = w
        W2[idx] += b.W2[idx]
    end
    a.n += b.n
    return a
end

function value(cm::ConditionalMoments{T}) where {T}
    W = permutedims(cm.W)
    W2 = permutedims(cm.W2)
    S = permutedims(cm.S)
    nan = T(NaN)
    has = W .> 0
    M2star = ifelse.(has, S ./ W, nan)
    neff = ifelse.(has, W .* W ./ W2, nan)
    M1 = ifelse.(has, permutedims(cm.M1), nan)
    M2 = ifelse.(has, permutedims(cm.M2), nan)
    return (x = copy(cm.x), lags = cm.lags .* cm.Δt, W = W, W2 = W2, M1 = M1, M2 = M2, M2star = M2star, neff = neff,
            n = cm.n, h = cm.h)
end

"""
    kramers_moyal(cm; lags = :, window = nothing, method = :regression, variance = false) -> NamedTuple

The drift D⁽¹⁾ and diffusion D⁽²⁾ at the evaluation points from the conditional moments of a
`ConditionalMoments` (or its `value`), Davis 2023 §3, over the lags selected by index
(`lags`, into the lag set) or by time (`window = (τ_lo, τ_hi)`, the lags with τ_lo ≤ τ ≤ τ_hi).
`method = :direct` is Eq. 16 at the first selected lag, D⁽ᵏ⁾ = M⁽ᵏ⁾/(k! τ); `:regression` is
the least-squares fit of M⁽ᵏ⁾(τ) ≃ k! τ D⁽ᵏ⁾ through the origin over the selected lags (Eq. 17). `variance = true`
builds D⁽²⁾ from the conditional variance M⁽²*⁾ instead of M⁽²⁾, which removes the (τD⁽¹⁾)²
term of the finite-τ second moment. Returns `(x, lags, D1, D2, D1_err, D2_err, n)` with the
standard errors from the conditional variance and the effective sample size at each lag
(exact for independent increments at a single lag; Gaussian increments assumed for the
error of D⁽²⁾). The moments of different lags share their increments, so the regression's
error is the coherent sum Σ τᵢσᵢ/Σ τᵢ², the fully correlated case, an upper bound that
reduces to σ/τ for one lag. In HOS's Fokker–Planck convention (Eq. 32) the drift is
S = D⁽¹⁾ and the diffusion is D = 2D⁽²⁾.
"""
kramers_moyal(cm::ConditionalMoments; kwargs...) = kramers_moyal(value(cm); kwargs...)

function kramers_moyal(v::NamedTuple; lags = Colon(), window = nothing, method::Symbol = :regression, variance::Bool = false)
    method in (:direct, :regression) || throw(ArgumentError("method must be :direct or :regression, got $method"))
    idx = if window !== nothing
        lags isa Colon || throw(ArgumentError("select the lags by index (lags) or by time (window), not both"))
        lo, hi = window
        lo ≤ hi || throw(ArgumentError("need window = (τ_lo, τ_hi) with τ_lo ≤ τ_hi"))
        findall(τ -> lo ≤ τ ≤ hi, v.lags)
    else
        lags isa Colon ? collect(eachindex(v.lags)) : collect(lags)
    end
    isempty(idx) && throw(ArgumentError("no lags selected"))
    all(1 .≤ idx .≤ length(v.lags)) || throw(ArgumentError("lags out of range 1:$(length(v.lags))"))
    method === :direct && (idx = idx[1:1])
    T = eltype(v.x)
    Mk2 = variance ? v.M2star : v.M2
    nx = length(v.x)
    D1, D2, D1_err, D2_err = (fill(T(NaN), nx) for _ in 1:4)
    n = zeros(T, nx)
    for j in 1:nx
        num1 = num2 = den = err1 = err2 = zero(T)
        used = false
        for i in idx
            w = v.W[i, j]
            (w > 0 && isfinite(v.M2star[i, j]) && v.neff[i, j] > 1) || continue
            used = true
            τ = T(v.lags[i])
            σ1 = sqrt(v.M2star[i, j] / v.neff[i, j])            # σ(M⁽¹⁾)
            σ2 = sqrt(2 / v.neff[i, j]) * v.M2star[i, j]        # σ(M⁽²⁾) for Gaussian increments
            den += τ * τ
            num1 += τ * v.M1[i, j]
            num2 += τ * Mk2[i, j]
            err1 += τ * σ1                                      # coherent sum: the lags share their increments
            err2 += τ * σ2
            n[j] = max(n[j], v.neff[i, j])
        end
        used || continue
        D1[j] = num1 / den
        D2[j] = num2 / (2den)
        D1_err[j] = err1 / den
        D2_err[j] = err2 / (2den)
    end
    return (x = copy(v.x), lags = v.lags[idx], D1 = D1, D2 = D2, D1_err = D1_err, D2_err = D2_err, n = n)
end
