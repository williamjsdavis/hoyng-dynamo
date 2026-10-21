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

# ---------------------------------------------------------------------------------------
# χ² fit of the stationary amplitude distribution to a histogram (HSO §4.2–4.3)
# ---------------------------------------------------------------------------------------

"""
    nelder_mead(f, x0; step = 0.5, xtol = 1e-6, ftol = 1e-10, maxiter = 2000) -> (x, f(x))

Downhill-simplex minimisation (Nelder & Mead 1965) of `f` from `x0` with an initial simplex
of edge `step`; reflection 1, expansion 2, contraction ½, shrink ½. Stops when the simplex
is smaller than `xtol` in every coordinate and the function spread is below `ftol`, or
after `maxiter` iterations.
"""
function nelder_mead(f, x0::AbstractVector{T}; step = T(0.5), xtol = T(1e-6), ftol = T(1e-10),
                     maxiter::Integer = 2000) where {T<:AbstractFloat}
    n = length(x0)
    X = [copy(x0) for _ in 1:n+1]
    for i in 1:n
        X[i+1][i] += step
    end
    F = T[f(x) for x in X]
    c = zeros(T, n)
    for _ in 1:maxiter
        order = sortperm(F)
        permute!(X, order); permute!(F, order)     # in place: X and F are never rebound (keeps inference)
        size = zero(T)
        for i in 2:n+1, k in 1:n
            size = max(size, abs(X[i][k] - X[1][k]))
        end
        size < xtol && F[end] - F[1] < ftol && break
        fill!(c, zero(T))                          # centroid of the best n vertices
        for i in 1:n
            c .+= X[i]
        end
        c ./= n
        xr = c .+ (c .- X[end]); fr = f(xr)
        if fr < F[1]
            xe = c .+ 2 .* (c .- X[end]); fe = f(xe)
            if fe < fr
                X[end] = xe; F[end] = fe
            else
                X[end] = xr; F[end] = fr
            end
        elseif fr < F[n]
            X[end] = xr; F[end] = fr
        else
            xc = fr < F[end] ? c .+ (xr .- c) ./ 2 : c .+ (X[end] .- c) ./ 2
            fc = f(xc)
            if fc < min(fr, F[end])
                X[end] = xc; F[end] = fc
            else
                for i in 2:n+1
                    X[i] = X[1] .+ (X[i] .- X[1]) ./ 2
                    F[i] = f(X[i])
                end
            end
        end
    end
    i = argmin(F)
    return (X[i], F[i])
end

"""
    fold_histogram(edges, counts) -> (edges, counts)

Fold a histogram of a whose bins are symmetric about zero onto |a| ≥ 0.
"""
function fold_histogram(edges::AbstractVector, counts::AbstractVector)
    n = length(counts)
    length(edges) == n + 1 || throw(ArgumentError("need one more edge than bins"))
    iseven(n) || throw(ArgumentError("need an even number of bins symmetric about zero"))
    h = n ÷ 2
    abs(edges[h+1]) ≤ 1e-12 * (edges[end] - edges[1]) || throw(ArgumentError("the bins must be symmetric about zero"))
    return (edges[h+1:end], counts[h+1:end] .+ counts[h:-1:1])
end

"""
    equal_population_bins(edges, counts; nbins = 100) -> (lo, hi, n, first, last)

Merge consecutive bins of a histogram into about `nbins` bins of nearly equal population
(HSO p. 149: variable widths so that no bin is empty and every bin has the same Poisson
error). Leading and trailing empty bins are dropped and a short last bin is merged into its
neighbour. `first`/`last` are the index ranges of the fine bins each merged bin covers.
"""
function equal_population_bins(edges::AbstractVector, counts::AbstractVector{<:Integer}; nbins::Integer = 100)
    nbins ≥ 1 || throw(ArgumentError("need nbins ≥ 1, got $nbins"))
    i1 = findfirst(>(0), counts)
    i1 === nothing && throw(ArgumentError("the histogram is empty"))
    i2 = something(findlast(>(0), counts))         # cannot be nothing once i1 exists; `something` tells inference so
    target = sum(counts) / nbins
    ifirst = Int[]; ilast = Int[]; n = Int[]
    acc = 0; start = i1
    for i in i1:i2
        acc += counts[i]
        nxt = i < i2 ? counts[i+1] : 0
        # close the bin at whichever side of the target is nearer: now, or after the next fine bin
        if i == i2 || acc ≥ target || (acc + nxt > target && acc + nxt - target > target - acc)
            push!(ifirst, start); push!(ilast, i); push!(n, acc)
            acc = 0; start = i + 1
        end
    end
    if length(n) ≥ 2 && n[end] < target / 2
        ilast[end-1] = ilast[end]; n[end-1] += n[end]
        pop!(ifirst); pop!(ilast); pop!(n)
    end
    return (lo = edges[ifirst], hi = edges[ilast .+ 1], n = n, first = ifirst, last = ilast)
end

"""
    density_fit(edges, counts, Λ, μ = Λ; nbins = 100, start = (γ = 5.0, c2 = 0.3)) -> NamedTuple

χ² fit of the stationary distribution p(a) ∝ exp ψ(a) of `BistableTheory(γ, c²; Λ, μ)`
(HOS Eq. 43 = HSO Eq. 10 for μ = Λ; HSO Eq. A.1 for μ > Λ) to a histogram of |a|, following
HSO §4.3: the fine bins are merged into about `nbins` bins of equal population; the model
count in each is N Σᵢ p(cᵢ)Δᵢ over the fine bins it covers (HSO take p at the centroid
times the width, the same to the order of the bin width); N is fixed by the total count;
χ² = Σⱼ (N Pⱼ − nⱼ)²/nⱼ is minimised over (γ, c²) by Nelder–Mead in (ln(γ − 1), ln c²);
the goodness of fit is quoted as σ = |χ² − ν|/√(2ν) with ν = nbins − 3 (Table 1,
footnote a). Only the shape is fitted: the a-axis of a simulation is already normalised,
its peak sitting at a_m. Returns `(γ, c2, chi2, ν, sigma, nbins, n)`.
"""
function density_fit(edges::AbstractVector{T}, counts::AbstractVector{<:Integer}, Λ::Real, μ::Real = Λ;
                     nbins::Integer = 100, start = (γ = 5.0, c2 = 0.3)) where {T<:AbstractFloat}
    b = equal_population_bins(edges, counts; nbins)
    centers = (edges[1:end-1] .+ edges[2:end]) ./ 2
    widths = diff(edges)
    total = sum(b.n)
    P = zeros(T, length(b.n))
    logp = zeros(T, length(centers))
    function chi2(x::AbstractVector{T})
        th = BistableTheory(1 + exp(x[1]), exp(x[2]); Λ, μ)
        ψmax = psi_max(th)
        @inbounds for i in eachindex(centers)
            logp[i] = psi(th, centers[i]) - ψmax
        end
        @inbounds for j in eachindex(P)
            s = zero(T)
            for i in b.first[j]:b.last[j]
                s += exp(logp[i]) * widths[i]
            end
            P[j] = s
        end
        N = total / sum(P)
        q = zero(T)                                # not named like any outer local: a closure would share it
        @inbounds for j in eachindex(P)
            q += (N * P[j] - b.n[j])^2 / b.n[j]
        end
        return q
    end
    x, χ² = nelder_mead(chi2, T[log(start.γ - 1), log(start.c2)]; step = T(0.3))
    ν = length(b.n) - 3
    return (γ = 1 + exp(x[1]), c2 = exp(x[2]), chi2 = χ², ν = ν, sigma = abs(χ² - ν) / sqrt(2ν),
            nbins = length(b.n), n = total)
end
