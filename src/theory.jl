# The Fokker–Planck formula module (HOS §4.2–4.5, Appendices A and F; HSO Appendix A).
# Pure functions of the shape parameters (γ, c²), the growth rate Λ and the ratio μ/Λ;
# nothing here touches the simulation (dev/DECISIONS.md D12). Notation is HOS's.
#
#     ∂p/∂τ = −∂(Sp)/∂a + ½ ∂²(Dp)/∂a²                                   (HOS Eq. 32)
#     S = Λ(1 − a²)a [1 + κa²],   κ = μ/Λ − 1                           (Eq. 33; HSO App. A)
#     D = Δ₀a² + Δ₁ = Δ₀(a² + c²),   γ = Λ/Δ₀,   c² = Δ₁/Δ₀              (Eqs. 39, 44)
#     p(a) ∝ (1/D) exp ∫₀^a 2S/D da' ≡ exp ψ(a),   ψ(0) = 0              (Eq. 42)
#
# κ = 0 gives the closed forms of HOS Eqs. 43 and 51–60; κ ≠ 0 gives HSO Eq. A.1, whose
# ψ follows from the same integral with the substitution x = a'²:
#
#     ψ(a) = γ[(κ(1+c²) − 1)a² − ½κa⁴] + [γ(1+c²)(1−κc²) − 1] ln(1 + a²/c²),
#
# and the formulae for the peak, the variance and ⟨T_r⟩ (A.2) are written in terms of
# a_m, ψ_m, ψ_m'' and φ₀'' so that they reduce to Eqs. 52–58 at κ = 0 (checked by gate
# T3.2–T3.5) and remain valid for κ ≠ 0 (gate T3.7).
#
# The quadratures (normalisation, the moments over a > 0 for dv, the escape integral
# Eq. 57) evaluate every exponential relative to its maximum ψ_max = max ψ, so nothing
# overflows for large γ; the results that are themselves ∝ exp(ψ_m) are offered in
# logarithmic form as well.

"""
    BistableTheory(γ, c²; Λ = 0.527, μ = Λ)
    BistableTheory(constants::ModelConstants, γ, c²)

The reduced description of the dipole amplitude a (HOS §4.2): drift S = Λ(1 − a²)a[1 + κa²]
with κ = μ/Λ − 1, diffusion D = Δ₀a² + Δ₁ with Δ₀ = Λ/γ and Δ₁ = c²Λ/γ. With μ = Λ
(κ = 0) this is HOS Eq. 33 and the amplitude distribution is Eq. 43; with the model's μ
(HOS Eq. 34, μ/Λ = 1.45 for the standard parameters) the drift carries the next term of
HSO Appendix A and the distribution is HSO Eq. A.1.

Stored: `γ`, `c2` (= c²), `Λ`, `κ`.
"""
struct BistableTheory{T<:AbstractFloat}
    γ::T
    c2::T
    Λ::T
    κ::T
    function BistableTheory{T}(γ, c2, Λ, κ) where {T}
        γ > 0 || throw(ArgumentError("γ must be positive, got $γ"))
        c2 > 0 || throw(ArgumentError("c² must be positive, got $c2"))
        Λ > 0 || throw(ArgumentError("Λ must be positive, got $Λ"))
        κ >= 0 || throw(ArgumentError("μ/Λ must be ≥ 1 (the overtones only weaken the restoring drift), got $(κ + 1)"))
        new{T}(γ, c2, Λ, κ)
    end
end

function BistableTheory(γ::Real, c2::Real; Λ::Real = 0.527, μ::Real = Λ)
    T = float(promote_type(typeof(γ), typeof(c2), typeof(Λ), typeof(μ)))
    BistableTheory{T}(T(γ), T(c2), T(Λ), T(μ) / T(Λ) - 1)
end

BistableTheory(c::ModelConstants{T}, γ::Real, c2::Real) where {T} = BistableTheory{T}(T(γ), T(c2), c.Λ, c.μ / c.Λ - 1)

Base.show(io::IO, th::BistableTheory{T}) where {T} =
    print(io, "BistableTheory{$T}(γ = $(th.γ), c² = $(th.c2), Λ = $(th.Λ), μ/Λ = $(th.κ + 1))")

"Δ₀ = Λ/γ and Δ₁ = c²Λ/γ (HOS Eq. 44 inverted), as a NamedTuple."
diffusion_constants(th::BistableTheory) = (Δ0 = th.Λ / th.γ, Δ1 = th.c2 * th.Λ / th.γ)

"γ = Λ/Δ₀ and c² = Δ₁/Δ₀ (HOS Eq. 44), as a NamedTuple."
shape_parameters(Λ::Real, Δ0::Real, Δ1::Real) = (γ = Λ / Δ0, c2 = Δ1 / Δ0)

"""
    diffusion_prefactor(f, τc, Nc)

ℱ ≡ πf²τ_c/(2N_c²), the factor in front of the integrals of HOS Eqs. 36–37 (3.22×10⁻² for
the standard run, p. 291).
"""
diffusion_prefactor(f::Real, τc::Real, Nc::Real) = π * f^2 * τc / (2 * Nc^2)

"""
    diffusion_from_fit(c0, c2; D0int, prefactor)

Δ₀ and Δ₁ (HOS Eq. 40) from a measured D₁/ℱ ≃ c0 + c2·a² (the fit of HOS Fig. 8, Eq. 61),
the integral `D0int` of Eq. 36 (`ModelConstants.D0int`, 3.529) and the prefactor ℱ:
Δ₀ = ℱ(D0int + c2), Δ₁ = ℱ c0 (Eqs. 61–63).
"""
diffusion_from_fit(c0::Real, c2::Real; D0int::Real, prefactor::Real) =
    (Δ0 = prefactor * (D0int + c2), Δ1 = prefactor * c0)

# ---------------------------------------------------------------------------------------
# Drift, diffusion, potential, log-density and its derivatives
# ---------------------------------------------------------------------------------------

"S(a) = Λ(1 − a²)a[1 + κa²] (HOS Eq. 33; HSO Appendix A for κ ≠ 0)."
@inline drift(th::BistableTheory, a::Real) = th.Λ * (1 - a^2) * a * (1 + th.κ * a^2)

"D(a) = Δ₀(a² + c²) (HOS Eq. 39)."
@inline diffusion(th::BistableTheory, a::Real) = (th.Λ / th.γ) * (a^2 + th.c2)

"""
    potential(th, a)

The bistable potential U with −∂U/∂a = S and U(±1) = 0: ¼Λ(1 − a²)² for κ = 0 (HOS Eq. 41),
Λ[¼ + κ/12 − a²/2 + (1 − κ)a⁴/4 + κa⁶/6] in general.
"""
@inline function potential(th::BistableTheory{T}, a::Real) where {T}
    κ, x = th.κ, a^2
    th.Λ * (T(1) / 4 + κ / 12 - x / 2 + (1 - κ) * x^2 / 4 + κ * x^3 / 6)
end

"""
    psi(th, a)

ψ(a) = ln p(a) − ln p(0) (HOS Eq. 51 for κ = 0; the exponent of HSO Eq. A.1 in general).
"""
@inline function psi(th::BistableTheory, a::Real)
    γ, c2, κ = th.γ, th.c2, th.κ
    x = a^2
    u = x / c2
    isfinite(u) || return oftype(u, -Inf)        # a² overflowed: the −γa² (or −γκa⁴/2) term wins
    γ * ((κ * (1 + c2) - 1) * x - κ * x^2 / 2) + (γ * (1 + c2) * (1 - κ * c2) - 1) * log1p(u)
end

"ψ'(a) = 2a[γ(1 − a²)(1 + κa²) − 1]/(a² + c²)."
@inline function dpsi(th::BistableTheory, a::Real)
    x = a^2
    2a * (th.γ * (1 - x) * (1 + th.κ * x) - 1) / (x + th.c2)
end

"ψ''(a)."
@inline function d2psi(th::BistableTheory, a::Real)
    γ, c2, κ = th.γ, th.c2, th.κ
    x = a^2
    h = γ * (1 - x) * (1 + κ * x) - 1            # ψ' = 2a h/(a² + c²)
    dh = 2γ * a * (κ - 1 - 2κ * x)
    2 * (h + a * dh) / (x + c2) - 4x * h / (x + c2)^2
end

"φ(y) with Δ₁[D(y)p(y)]⁻¹ = exp φ(y) (HOS Appendix F): φ = −ψ − ln(1 + y²/c²)."
@inline phi(th::BistableTheory, y::Real) = -psi(th, y) - log1p(y^2 / th.c2)

"φ''(0) = −2γ/c² (HOS Appendix F; unchanged by κ)."
@inline phi0_curvature(th::BistableTheory) = -2 * th.γ / th.c2

# ---------------------------------------------------------------------------------------
# The peaks of p(a): closed forms (HOS Eqs. 52–56, generalised to κ ≠ 0)
# ---------------------------------------------------------------------------------------

"Whether p(a) is bistable (two peaks at ±a_m): γ > 1 (HOS p. 288)."
bistable(th::BistableTheory) = th.γ > 1

@inline function check_bistable(th::BistableTheory)
    bistable(th) || throw(ArgumentError("p(a) has a single peak at a = 0 for γ = $(th.γ) ≤ 1"))
    nothing
end

"""
    peak_amplitude(th)

a_m > 0, the maximum of p(a): (1 − 1/γ)^½ (HOS Eq. 52) for κ = 0; in general the positive
root of (1 + κa_m²)(1 − a_m²) = 1/γ (HSO p. 156).
"""
function peak_amplitude(th::BistableTheory{T}) where {T}
    check_bistable(th)
    κ, r = th.κ, 1 - 1 / th.γ
    sqrt(2r / ((1 - κ) + sqrt((1 - κ)^2 + 4κ * r)))        # exact for κ = 0, stable for small κ
end

"ψ_m = ψ(a_m) (HOS Eq. 53 for κ = 0)."
peak_log_density(th::BistableTheory) = psi(th, peak_amplitude(th))

"ψ_m'' = ψ''(a_m) = 4γa_m²(κ − 1 − 2κa_m²)/(a_m² + c²) (HOS Eq. 54 for κ = 0)."
function peak_curvature(th::BistableTheory)
    x = peak_amplitude(th)^2
    4 * th.γ * x * (th.κ - 1 - 2 * th.κ * x) / (x + th.c2)
end

"⟨a⟩ ≃ a_m over one peak (HOS Eq. 55)."
mean_amplitude(th::BistableTheory) = peak_amplitude(th)

"""
    relative_variance(th)

⟨(a − ⟨a⟩)²⟩/⟨a⟩² over one peak in the parabolic approximation, −1/(ψ_m'' a_m²): equals
(γc² + γ − 1)/[4(γ − 1)²] (HOS Eq. 56) for κ = 0.
"""
relative_variance(th::BistableTheory) = -1 / (peak_curvature(th) * peak_amplitude(th)^2)

# ---------------------------------------------------------------------------------------
# Mean time between reversals: closed forms (HOS Eqs. 58–60, HSO Eq. A.2)
# ---------------------------------------------------------------------------------------

"""
    log_mean_reversal_time(th)

ln⟨T_r⟩ with ⟨T_r⟩ ≈ 2⟨T_e⟩ = (4π/Δ₁)(φ₀''ψ_m'')^(−½) exp(ψ_m) (HOS Eq. F.3, HSO Eq. A.2),
which is Eq. 58 for κ = 0. Logarithmic so that it never overflows.
"""
function log_mean_reversal_time(th::BistableTheory)
    Δ1 = diffusion_constants(th).Δ1
    log(4π / Δ1) - log(phi0_curvature(th) * peak_curvature(th)) / 2 + peak_log_density(th)
end

"⟨T_r⟩ in diffusion times (HOS Eq. 58 for κ = 0, HSO Eq. A.2 in general); `Inf` if it overflows."
mean_reversal_time(th::BistableTheory) = exp(log_mean_reversal_time(th))

"K = −1 + (1 + c²) ln[(1 + c²)/c²] (HOS Eq. 60), the Kramers exponent of Eq. 59."
kramers_exponent(th::BistableTheory) = -1 + (1 + th.c2) * log((1 + th.c2) / th.c2)

@inline function check_hos_form(th::BistableTheory)
    th.κ == 0 || throw(ArgumentError("HOS Eq. 59 is the large-γ form of Eq. 58 and needs μ = Λ (κ = 0)"))
    nothing
end

"ln of HOS Eq. 59, the large-γ asymptote of Eq. 58: (π√2/Λ)[c²/(1 + c²)]^½ exp(Kγ); κ = 0 only."
function log_mean_reversal_time_asymptotic(th::BistableTheory)
    check_hos_form(th)
    log(π * sqrt(2 * one(th.Λ)) / th.Λ) + log(th.c2 / (1 + th.c2)) / 2 + kramers_exponent(th) * th.γ
end

"HOS Eq. 59, the large-γ asymptote of Eq. 58; κ = 0 only. `Inf` if it overflows."
mean_reversal_time_asymptotic(th::BistableTheory) = exp(log_mean_reversal_time_asymptotic(th))

# ---------------------------------------------------------------------------------------
# Quadratures, in log space relative to max ψ (D12)
# ---------------------------------------------------------------------------------------

"max_a ψ(a): ψ_m if bistable, else ψ(0) = 0."
psi_max(th::BistableTheory{T}) where {T} = bistable(th) ? peak_log_density(th) : zero(T)

const THEORY_RTOL = 1e-10

# ∫₀^∞ aⁿ exp[ψ(a) − ψ_max] da
function scaled_half_moment(th::BistableTheory{T}, n::Int, ψmax::T) where {T}
    I, _ = quadgk(a -> a^n * exp(psi(th, a) - ψmax), zero(T), T(Inf); rtol = THEORY_RTOL)
    return I
end

"""
    log_normalisation(th)

ln Z with Z = ∫_{−∞}^{∞} exp ψ(a) da, so that p(a) = exp[ψ(a) − ln Z] has unit area.
"""
function log_normalisation(th::BistableTheory{T}) where {T}
    ψmax = psi_max(th)
    ψmax + log(2 * scaled_half_moment(th, 0, ψmax))
end

"""
    StationaryDensity(th)

The normalised stationary amplitude distribution p(a) (HOS Eq. 43, or HSO Eq. A.1 for
κ ≠ 0), callable as `p(a)`; `logpdf(p, a)` gives ln p(a).
"""
struct StationaryDensity{T<:AbstractFloat}
    theory::BistableTheory{T}
    logZ::T
end

StationaryDensity(th::BistableTheory) = StationaryDensity(th, log_normalisation(th))
stationary_density(th::BistableTheory) = StationaryDensity(th)

@inline logpdf(p::StationaryDensity, a::Real) = psi(p.theory, a) - p.logZ
@inline (p::StationaryDensity)(a::Real) = exp(logpdf(p, a))

"""
    half_moments(th)

⟨1⟩, ⟨|a|⟩ and ⟨a²⟩ of p(a) restricted to a > 0 (HSO Eq. 12's averages), as a NamedTuple
`(mean, meansq)`, computed by quadrature.
"""
function half_moments(th::BistableTheory)
    ψmax = psi_max(th)
    m0 = scaled_half_moment(th, 0, ψmax)
    m1 = scaled_half_moment(th, 1, ψmax)
    m2 = scaled_half_moment(th, 2, ψmax)
    return (mean = m1 / m0, meansq = m2 / m0)
end

"""
    dipole_variation(th)

dv = ⟨(|a| − ⟨|a|⟩)²⟩/⟨|a|⟩² (HSO Eq. 12) evaluated on the theoretical p(a) over a > 0, by
quadrature. Unlike `relative_variance` this uses the exact distribution and includes the
reversal region near a = 0.
"""
function dipole_variation(th::BistableTheory)
    m = half_moments(th)
    (m.meansq - m.mean^2) / m.mean^2
end

"""
    log_escape_time(th; a′ = 1)

ln⟨T_e⟩ of HOS Eq. 57 (Appendix F.1), the mean time to reach a = 0 from a′ > 0 (the paper
takes a′ ≈ 1, the well bottom; a′ = Inf is allowed), evaluated by nested quadrature:

    ⟨T_e⟩ = (2/Δ₁) ∫₀^a′ dy exp φ(y) ∫_y^∞ dz exp ψ(z).

Both exponentials are taken relative to ψ_max, so the value is finite for any γ.
"""
# The integrands of the nested escape integral are callable structs. Julia's inference
# still gives up on the inner quadgk call (its recursion limit fires when QuadGK's
# internals are re-entered from inside an outer quadgk), so the inner result is asserted
# to be a T; the D14 tests check these two functions with @inferred only, not with JET.

# exp[ψ(z) − ψref]
struct ScaledDensity{T}
    th::BistableTheory{T}
    ψref::T
end

@inline (g::ScaledDensity{T})(z::T) where {T} = exp(psi(g.th, z) - g.ψref)

# ∫_y^∞ exp[ψ(z) − ψref] dz
function scaled_tail(th::BistableTheory{T}, y::T, ψref::T) where {T}
    I, _ = quadgk(ScaledDensity(th, ψref), y, T(Inf); rtol = THEORY_RTOL)
    return I::T
end

# The outer integrand of Eq. 57 divided by exp(ψ_max). The inner integral is taken relative
# to its own maximum over [y, ∞), which is ψ_max for y ≤ a_m and ψ(y) beyond the peak, so
# that the product exp φ(y) ∫_y^∞ exp ψ never meets Inf × 0 (it decays like y⁻³ for large y).
struct EscapeIntegrand{T}
    th::BistableTheory{T}
    ψmax::T
    am::T
end

function (h::EscapeIntegrand{T})(y::T) where {T}
    ψy = psi(h.th, y)
    ψref = y <= h.am ? h.ψmax : ψy
    exp(ψref - ψy - h.ψmax) * scaled_tail(h.th, y, ψref) / (1 + y^2 / h.th.c2)
end

function log_escape_time(th::BistableTheory{T}; a′::Real = 1) where {T}
    a′ > 0 || throw(ArgumentError("the starting point a′ must be positive, got $a′"))
    ψmax = psi_max(th)
    am = bistable(th) ? peak_amplitude(th) : zero(T)
    Δ1 = diffusion_constants(th).Δ1
    I, _ = quadgk(EscapeIntegrand(th, ψmax, am), zero(T), T(a′); rtol = 1e-8)
    return log(2 / Δ1) + ψmax + log(I)
end

"⟨T_e⟩ (HOS Eq. 57) by quadrature; `Inf` if it overflows."
escape_time(th::BistableTheory; a′::Real = 1) = exp(log_escape_time(th; a′))

"""
    log_mean_reversal_time_quadrature(th; a′ = 1)

ln of ⟨T_r⟩ = 2⟨T_e⟩ with the escape integral evaluated numerically (HOS Eq. 57; with
a′ = 1 this is the 4∫₀¹ form of SOH p. 122), i.e. without the Gaussian approximation
F.2 that leads to Eq. 58.
"""
log_mean_reversal_time_quadrature(th::BistableTheory; a′::Real = 1) = log(2 * one(th.γ)) + log_escape_time(th; a′)

"⟨T_r⟩ = 2⟨T_e⟩ by quadrature (HOS Eq. 57, SOH's 4∫ form at a′ = 1); `Inf` if it overflows."
mean_reversal_time_quadrature(th::BistableTheory; a′::Real = 1) = exp(log_mean_reversal_time_quadrature(th; a′))
