# Tier 6 (dev/VALIDATION.md): the Sint-800 record of Guyodo & Valet (1999) and the analysis
# HSO applied to it: dv (Eq. 12) with Monte-Carlo errors, the χ² fit of p(a) to the VADM
# histogram (§4.2), the autocorrelation (§5.1) and the τ_d measurement by overlaying the
# autocorrelation of an artificial Sint-800 series drawn from a run. Notation is HSO's where
# it differs from HOS's (dv, A_m, τ_av). The record is `data/sint800.tsv` (see data/README.md).

"""
    read_sint800(path = joinpath(pkgdir(HoyngDynamo), "data", "sint800.tsv")) -> NamedTuple

The Sint-800 record: `age` (ka, 2 … 800, one point per kyr), `ncores`, the relative stack
`sint800` with its `stderr2` and `stderr1`, and the calibrated `vadm` (10²² A m²) with its
standard error `vadm_err`; `n` = 799 points, HSO's "799 points" (p. 148). The record is
unsigned: HSO's Fig. 1 flips the sign of the points older than the reversal at 780 ka.
"""
function read_sint800(path::AbstractString = joinpath(pkgdir(HoyngDynamo), "data", "sint800.tsv"))
    age = Float64[]; ncores = Int[]; rpi = Float64[]; e2 = Float64[]; e1 = Float64[]; vadm = Float64[]; err = Float64[]
    for line in eachline(path)
        (isempty(strip(line)) || startswith(line, '#') || startswith(line, "age_ka")) && continue
        f = split(line, '\t')
        length(f) == 7 || error("malformed Sint-800 line: $line")
        push!(age, parse(Float64, f[1])); push!(ncores, round(Int, parse(Float64, f[2]))); push!(rpi, parse(Float64, f[3]))
        push!(e2, parse(Float64, f[4])); push!(e1, parse(Float64, f[5])); push!(vadm, parse(Float64, f[6])); push!(err, parse(Float64, f[7]))
    end
    return (age = age, ncores = ncores, sint800 = rpi, stderr2 = e2, stderr1 = e1, vadm = vadm, vadm_err = err, n = length(age))
end

"""
    dipole_variation(x::AbstractVector) -> dv

dv = ⟨(|x| − ⟨|x|⟩)²⟩/⟨|x|⟩² of a series (HSO Eq. 12; population averages).
"""
function dipole_variation(x::AbstractVector{<:Real})
    n = length(x)
    m = sum(abs, x) / n
    return sum(v -> (abs(v) - m)^2, x) / n / m^2
end

"""
    montecarlo(f, s, e; n = 1000, rng = Xoshiro(1)) -> (value, mean, sigma, samples)

HSO's error estimate (pp. 151–152): `f` evaluated on `n` artificial records s + e·N(0, 1)
with independent standard-normal draws; `value = f(s)`, and the mean and standard deviation
of the replicas (elementwise if `f` returns a vector).
"""
function montecarlo(f, s::AbstractVector{<:Real}, e::AbstractVector{<:Real}; n::Integer = 1000, rng::AbstractRNG = Xoshiro(1))
    length(s) == length(e) || throw(ArgumentError("the record and its errors differ in length"))
    n ≥ 2 || throw(ArgumentError("need at least two replicas"))
    samples = [f(s .+ e .* randn(rng, length(s))) for _ in 1:n]
    mean = sum(samples) / n
    sigma = sqrt.(sum(x -> (x .- mean) .^ 2, samples) ./ (n - 1))
    return (value = f(s), mean = mean, sigma = sigma, samples = samples)
end

"""
    hso_autocorrelation(x, L = 50) -> Vector

HSO's autocorrelation estimator (§5.1, p. 152): subtract the mean, c_k = Σ_{i=1}^{n−L} x_i x_{i+k}
for k = 0 … L with the same number of terms at every lag, normalised to c_0 = 1. For Sint-800
the record is fed in time order (age reversed); the model series is rectified beforehand.
"""
function hso_autocorrelation(x::AbstractVector{<:Real}, L::Integer = 50)
    n = length(x)
    1 ≤ L < n || throw(ArgumentError("need 1 ≤ L < length(x)"))
    y = x .- sum(x) / n
    c = [sum(i -> y[i] * y[i+k], 1:n-L) for k in 0:L]
    return c ./ c[1]
end

"""
    running_mean(x, w) -> Vector

Centred moving average of `x` over `w` samples (partial windows at the ends), by prefix sums.
"""
function running_mean(x::AbstractVector{<:Real}, w::Integer)
    w ≥ 1 || throw(ArgumentError("need w ≥ 1"))
    n = length(x)
    w == 1 && return float.(x)
    c = cumsum(float.(x))
    half = (w - 1) ÷ 2
    out = similar(c)
    for i in 1:n
        lo, hi = max(1, i - half), min(n, i - half + w - 1)
        out[i] = (c[hi] - (lo > 1 ? c[lo-1] : zero(eltype(c)))) / (hi - lo + 1)
    end
    return out
end

"""
    artificial_sint800(a, Δτ, τd; τav = 5.0, start = 0.0, npoints = 799, spacing = 1.0) -> Vector

HSO's recipe (§5.1): a simulated dipole amplitude `a` sampled every `Δτ` τ_d becomes an
artificial Sint-800 series once `τd` (kyr) fixes the time axis: rectify, take the running
average over `τav` kyr, pick `npoints` points `spacing` kyr apart from `start` kyr (linear
interpolation between samples) and subtract the mean. With `τav = 0` no averaging is done.
"""
function artificial_sint800(a::AbstractVector{<:Real}, Δτ::Real, τd::Real; τav::Real = 5.0, start::Real = 0.0,
                            npoints::Integer = 799, spacing::Real = 1.0)
    τd > 0 && Δτ > 0 || throw(ArgumentError("need positive τd and Δτ"))
    w = max(1, round(Int, τav / (τd * Δτ)))
    smooth = running_mean(abs.(a), w)
    return sample_series(smooth, Δτ, τd; start, npoints, spacing)
end

# `npoints` points of a smoothed series, `spacing` kyr apart from `start` kyr, mean subtracted
function sample_series(smooth::AbstractVector{<:Real}, Δτ::Real, τd::Real; start::Real, npoints::Integer, spacing::Real)
    n = length(smooth)
    out = Vector{Float64}(undef, npoints)
    for k in 1:npoints
        t = start + (k - 1) * spacing                       # kyr
        pos = t / (τd * Δτ) + 1                             # fractional sample index
        pos ≤ n + 1e-6 || throw(ArgumentError("the window reaches beyond the run ($(round(t)) kyr at τ_d = $τd kyr)"))
        i = clamp(floor(Int, pos), 1, n - 1)
        f = pos - i
        out[k] = (1 - f) * smooth[i] + f * smooth[i+1]
    end
    return out .- sum(out) / npoints
end

"""
    tau_d_scan(a, Δτ, data_acor; τd = 4:0.5:20, τav = 5.0, nwindows = 200, L = 50, lmax = 25) -> NamedTuple

HSO's τ_d measurement (§5.1, Fig. 7) made systematic: for each trial `τd` (kyr) the run `a`
(sampled every `Δτ` τ_d) is turned into `nwindows` artificial Sint-800 series starting at
evenly spaced times, their autocorrelations to lag `L` are compared with the data curve
`data_acor` by the rms difference over lags 1 … `lmax` kyr (the range HSO trust, p. 152),
and the best τd of every window is recorded. Returns the grid `τd`, the `misfit` matrix
(windows × τd), the per-window `best` τd, its `median`, `q25`, `q75`, the τd minimising
the window-averaged misfit (`best_mean`) and the mean model curves `acor` (τd × lags).
"""
function tau_d_scan(a::AbstractVector{<:Real}, Δτ::Real, data_acor::AbstractVector{<:Real}; τd = 4:0.5:20, τav::Real = 5.0,
                    nwindows::Integer = 200, L::Integer = 50, lmax::Integer = 25, npoints::Integer = 799)
    length(data_acor) ≥ lmax + 1 || throw(ArgumentError("the data autocorrelation must reach lag $lmax"))
    grid = collect(float.(τd))
    misfit = zeros(nwindows, length(grid))
    acor = zeros(length(grid), L + 1)
    span = (length(a) - 1) * Δτ                              # run length in τ_d
    for (j, τ) in enumerate(grid)
        w = max(1, round(Int, τav / (τ * Δτ)))
        smooth = running_mean(abs.(a), w)
        window = (npoints - 1) * 1.0                          # kyr
        total = span * τ * (1 - 1e-9)
        total > window || throw(ArgumentError("the run is shorter than one Sint-800 window at τ_d = $τ kyr"))
        starts = range(0, total - window; length = nwindows)
        for (i, s) in enumerate(starts)
            series = sample_series(smooth, Δτ, τ; start = s, npoints, spacing = 1.0)
            c = hso_autocorrelation(series, L)
            acor[j, :] .+= c
            misfit[i, j] = sqrt(sum(k -> (c[k+1] - data_acor[k+1])^2, 1:lmax) / lmax)
        end
    end
    acor ./= nwindows
    best = [grid[argmin(view(misfit, i, :))] for i in 1:nwindows]
    sorted = sort(best)
    q(p) = sorted[clamp(round(Int, p * (nwindows - 1)) + 1, 1, nwindows)]
    meanmisfit = vec(sum(misfit; dims = 1)) ./ nwindows
    return (τd = grid, misfit = misfit, best = best, median = q(0.5), q25 = q(0.25), q75 = q(0.75),
            best_mean = grid[argmin(meanmisfit)], mean_misfit = meanmisfit, acor = acor, lags = 0:L)
end

# ---------------------------------------------------------------------------------------
# χ² fit of p(a) to the VADM record (HSO §4.2)
# ---------------------------------------------------------------------------------------

"""
    vadm_bins(vadm; nbins = 100, equal_population = true, lo = 0) -> (lo, hi, n)

The bins of HSO §4.2: `nbins` bins of equal population (the quantiles of the record, so no
bin is empty and every bin has the same Poisson error) or of equal width over [`lo`, max]
(HOS Fig. 10 and SOH Fig. 4 bin the VADM axis from 0; `lo = minimum(vadm)` bins the data's
range); empty equal-width bins are dropped.
"""
function vadm_bins(vadm::AbstractVector{<:Real}; nbins::Integer = 100, equal_population::Bool = true, lo::Real = 0)
    nbins ≥ 2 || throw(ArgumentError("need at least two bins"))
    v = sort(float.(vadm))
    n = length(v)
    if equal_population
        edges = [v[1] * (1 - 1e-9); [(v[k] + v[k+1]) / 2 for k in round.(Int, (1:nbins-1) .* n ./ nbins)]; v[end] * (1 + 1e-9)]
    else
        lo ≤ v[1] || throw(ArgumentError("the lower edge $lo lies above the smallest value $(v[1])"))
        edges = collect(range(lo == v[1] ? lo * (1 - 1e-9) : lo, v[end] * (1 + 1e-9); length = nbins + 1))
    end
    counts = [count(x -> edges[j] ≤ x < edges[j+1], v) for j in 1:nbins]
    keep = counts .> 0
    return (lo = edges[1:end-1][keep], hi = edges[2:end][keep], n = counts[keep])
end

# ∫ exp(ψ(a) − ψmax) da over [a1, a2] by Simpson's rule on 20 subintervals
function bin_integral(th::BistableTheory{T}, ψmax, a1, a2; m::Integer = 20) where {T}
    h = (a2 - a1) / m
    s = exp(psi(th, a1) - ψmax) + exp(psi(th, a2) - ψmax)
    for k in 1:m-1
        s += (isodd(k) ? 4 : 2) * exp(psi(th, a1 + k * h) - ψmax)
    end
    return s * h / 3
end

"""
    vadm_chi2(vadm, th::BistableTheory, Am; bins = vadm_bins(vadm)) -> (chi2, ν, sigma, N, expected)

HSO §4.2 for given parameters: a = (A/A_m) a_m so that the peak of p(a) ∝ exp ψ(a) lies at
the VADM `Am`; the model count in bin j is N ∫ p da over the bin, N fixed by Σ_j N P_j = n;
χ² = Σ_j (N P_j − n_j)²/n_j; ν = nbins − 4 (three parameters and the normalisation);
σ = |χ² − ν|/√(2ν) (Table 1, footnote a).
"""
function vadm_chi2(vadm::AbstractVector{<:Real}, th::BistableTheory, Am::Real; bins = vadm_bins(vadm))
    Am > 0 || throw(ArgumentError("A_m must be positive"))
    scale = peak_amplitude(th) / Am
    ψmax = psi_max(th)
    P = [bin_integral(th, ψmax, lo * scale, hi * scale) for (lo, hi) in zip(bins.lo, bins.hi)]
    N = sum(bins.n) / sum(P)
    χ² = sum(((N * p - n)^2 / n for (p, n) in zip(P, bins.n)))
    ν = length(bins.n) - 4
    return (chi2 = χ², ν = ν, sigma = abs(χ² - ν) / sqrt(2ν), N = N, expected = N .* P)
end

"""
    vadm_fit(vadm, Λ, μ = Λ; nbins = 100, equal_population = true, lo = 0, fixed = nothing,
             start = (γ = 6.0, c2 = 0.3, Am = 6.0)) -> NamedTuple

HSO §4.2: the χ² of `vadm_chi2` minimised over (γ, c², A_m) by Nelder–Mead in
(ln(γ − 1), ln c², ln A_m), for p(a) of `BistableTheory(γ, c²; Λ, μ)` (Eq. 10 with μ = Λ,
Eq. A.1 with the model's μ) on the bins of `vadm_bins`. `fixed = (γ, c2)` optimises A_m alone (HSO's case 2). Returns
`(γ, c2, Am, chi2, ν, sigma, nbins, n, theory)` with the theoretical dv and ⟨T_r⟩ of the
fitted parameters as in Table 1.
"""
function vadm_fit(vadm::AbstractVector{<:Real}, Λ::Real, μ::Real = Λ; nbins::Integer = 100, equal_population::Bool = true,
                  lo::Real = 0, fixed = nothing, start = (γ = 6.0, c2 = 0.3, Am = 6.0))
    bins = vadm_bins(vadm; nbins, equal_population, lo)
    f = if fixed === nothing
        x -> vadm_chi2(vadm, BistableTheory(1 + exp(x[1]), exp(x[2]); Λ, μ), exp(x[3]); bins).chi2
    else
        x -> vadm_chi2(vadm, BistableTheory(fixed[1], fixed[2]; Λ, μ), exp(x[1]); bins).chi2
    end
    x0 = fixed === nothing ? [log(start.γ - 1), log(start.c2), log(start.Am)] : [log(start.Am)]
    x, _ = nelder_mead(f, x0; step = 0.3)
    γ, c2, Am = fixed === nothing ? (1 + exp(x[1]), exp(x[2]), exp(x[3])) : (float(fixed[1]), float(fixed[2]), exp(x[1]))
    th = BistableTheory(γ, c2; Λ, μ)
    r = vadm_chi2(vadm, th, Am; bins)
    ν = fixed === nothing ? r.ν : length(bins.n) - 2
    theory = (dv = dipole_variation(th), mean_reversal_time = mean_reversal_time(th))
    return (γ = γ, c2 = c2, Am = Am, chi2 = r.chi2, ν = ν, sigma = abs(r.chi2 - ν) / sqrt(2ν), nbins = length(bins.n),
            n = sum(bins.n), theory = theory, bins = bins, expected = r.expected, N = r.N, am = peak_amplitude(th), κ = th.κ)
end
