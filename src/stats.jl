# Small statistical tests used by the validation gates (dev/VALIDATION.md T4.11).

"""
    ks_exponential(x; tmin = 0, nboot = 2000, rng = Xoshiro(0)) -> (D, p, n, mean)

Kolmogorov–Smirnov test of the hypothesis that the values of `x` above `tmin` are `tmin`
plus an exponentially distributed variable, i.e. that the intervals form a Poisson process
once the short ones are excluded (HOS §4.5.3: polarity intervals are Poisson for
T_r ≫ Λ⁻¹; gate T4.11 uses `tmin = 5` τ_d). The exponential's mean is estimated from the
same data, so the p-value comes from a parametric bootstrap of `nboot` exponential samples
of the same size with the mean re-estimated in each (the Lilliefors correction; the
Kolmogorov distribution would overstate p). `D` is the supremum distance, `n` the number of
values above `tmin`, `mean` the fitted mean of x − tmin.
"""
function ks_exponential(x::AbstractVector{<:Real}; tmin::Real = 0, nboot::Integer = 2000,
                        rng::AbstractRNG = Xoshiro(0))
    y = sort!([Float64(v) - Float64(tmin) for v in x if v > tmin])
    n = length(y)
    n ≥ 2 || throw(ArgumentError("need at least two values above tmin = $tmin, got $n"))
    nboot ≥ 1 || throw(ArgumentError("need nboot ≥ 1, got $nboot"))
    m = sum(y) / n
    D = ks_distance_exponential(y, m)
    buf = similar(y)
    exceed = 0
    for _ in 1:nboot
        for i in eachindex(buf)
            buf[i] = m * randexp(rng)
        end
        sort!(buf)
        exceed += ks_distance_exponential(buf, sum(buf) / n) ≥ D
    end
    return (D = D, p = (exceed + 1) / (nboot + 1), n = n, mean = m)
end

"Supremum distance between the empirical distribution of the sorted `y` and Exp(mean `m`)."
function ks_distance_exponential(y::AbstractVector{Float64}, m::Float64)
    n = length(y)
    D = 0.0
    @inbounds for i in 1:n
        F = -expm1(-y[i] / m)
        D = max(D, F - (i - 1) / n, i / n - F)
    end
    return D
end
