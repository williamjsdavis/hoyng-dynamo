# Running aggregators (dev/DECISIONS.md D8): unit tests on synthetic series with known
# answers, against straightforward offline reference implementations. Always run.

using Random, LinearAlgebra

ref_mean(x) = sum(x) / length(x)
ref_var(x) = (m = ref_mean(x); sum(abs2, x .- m) / (length(x) - 1))

function ref_hist(x, lo, hi, nb)
    counts = zeros(Int, nb); below = above = 0
    for v in x
        k = floor(Int, (v - lo) / (hi - lo) * nb)
        if 0 ≤ k < nb; counts[k+1] += 1 elseif k < 0; below += 1 else above += 1 end
    end
    return counts, below, above
end

# the estimator of AutoCorrelation: Σ_{i>l} x_i x_{i−l}/(n − l) − x̄²
function ref_acov(x, L)
    n = length(x); m = ref_mean(x)
    return [dot(view(x, l+1:n), view(x, 1:n-l)) / (n - l) - m^2 for l in 0:L]
end

# Schmitt trigger polarity per sample: +1/−1, 0 before the first sample outside the band
function ref_polarity(a, ah)
    p = zeros(Int, length(a)); pol = 0
    for (i, v) in enumerate(a)
        if pol == 0
            abs(v) ≥ ah && (pol = sign(v))
        elseif (pol == 1 && v < -ah) || (pol == -1 && v > ah)
            pol = -pol
        end
        p[i] = pol
    end
    return p
end
ref_flips(a, ah) = (p = ref_polarity(a, ah); [i for i in 2:length(a) if p[i] != p[i-1] && p[i-1] != 0])

# between-reversal variance with the ±W exclusion: the variance about each interval's own
# mean and |mean| pooled over the intervals by kept-sample count, then the ratio
function ref_exclusion(a, ah, W, Δt)
    w = max(1, round(Int, W / Δt)); n = length(a)
    flips = ref_flips(a, ah)
    starts = vcat(1, flips .+ w)                       # first kept sample of each interval
    stops = vcat(flips .- 1 .- w, n - w)               # last kept sample
    M2 = 0.0; absmean = 0.0; pn = 0; k = 0
    for (s, e) in zip(starts, stops)
        e - s + 1 ≥ 2 || continue
        seg = a[s:e]
        M2 += sum(abs2, seg .- ref_mean(seg)); absmean += length(seg) * abs(ref_mean(seg))
        pn += length(seg); k += 1
    end
    return (variance = (M2 / (pn - k)) / (absmean / pn)^2, n = pn, intervals = k, M2 = M2, absmean = absmean)
end
ref_exclusion_pool(r1, r2) = ((r1.M2 + r2.M2) / (r1.n + r2.n - r1.intervals - r2.intervals)) / ((r1.absmean + r2.absmean) / (r1.n + r2.n))^2

@testset "Aggregators (D8)" begin
    rng = Xoshiro(2026)

    @testset "Moments" begin
        x = randn(rng, 10_000) .+ 0.7
        m = Moments(); foreach(v -> fit!(m, v), x)
        v = value(m)
        @test v.n == 10_000
        @test v.mean ≈ ref_mean(x) rtol = 1e-12
        @test v.var ≈ ref_var(x) rtol = 1e-12
        @test v.std ≈ sqrt(ref_var(x)) rtol = 1e-12
        # merge of two halves is the whole
        m1 = Moments(); m2 = Moments()
        foreach(v -> fit!(m1, v), x[1:3000]); foreach(v -> fit!(m2, v), x[3001:end])
        merge!(m1, m2)
        @test value(m1).mean ≈ v.mean rtol = 1e-12
        @test value(m1).var ≈ v.var rtol = 1e-12
        @test value(merge!(Moments(), Moments())).n == 0
        # large offset: Welford keeps the variance (a naive Σx² would lose it)
        big = Moments(); foreach(v -> fit!(big, 1e8 + v), x)
        @test value(big).var ≈ ref_var(x) rtol = 1e-6
    end

    @testset "Histogram" begin
        x = 1.5 .* randn(rng, 20_000)
        h = Histogram(-3.0, 3.0, 60); foreach(v -> fit!(h, v), x)
        counts, below, above = ref_hist(x, -3.0, 3.0, 60)
        v = value(h)
        @test v.counts == counts && v.below == below && v.above == above
        @test v.n == 20_000 && length(v.edges) == 61 && length(v.centers) == 60
        @test v.centers[1] ≈ -2.95 && v.centers[end] ≈ 2.95
        fit!(h, NaN); @test value(h).n == 20_000                      # NaN ignored
        h2 = Histogram(-3.0, 3.0, 60); foreach(v -> fit!(h2, v), x)
        @test value(merge!(h, h2)).counts == 2 .* counts
        @test_throws ArgumentError merge!(h, Histogram(-3.0, 3.0, 30))
    end

    @testset "AutoCorrelation" begin
        # AR(1): x_n = ρ x_{n−1} + ε, autocorrelation ρ^l
        ρ, n, L = 0.9, 200_000, 20
        x = zeros(n); for i in 2:n; x[i] = ρ * x[i-1] + randn(rng); end
        ac = AutoCorrelation(L; Δt = 0.5); foreach(v -> fit!(ac, v), x)
        v = value(ac)
        @test v.n == n && v.lags == (0:L) .* 0.5
        @test v.acov ≈ ref_acov(x, L) rtol = 1e-10
        @test maximum(abs.(v.acor .- ρ .^ (0:L))) < 0.03
        # chunked accumulation: more than one chunk of 2²⁰ samples, same estimator
        y = randn(rng, (1 << 20) + 12_345) .+ 0.3
        ac2 = AutoCorrelation(5); foreach(v -> fit!(ac2, v), y)
        @test value(ac2).acov ≈ ref_acov(y, 5) rtol = 1e-8
        @test value(ac2).mean ≈ ref_mean(y) rtol = 1e-12
        # merge pools the lagged sums of independent runs
        a1 = AutoCorrelation(L); a2 = AutoCorrelation(L)
        foreach(v -> fit!(a1, v), x[1:100_000]); foreach(v -> fit!(a2, v), x[100_001:end])
        S = [dot(view(x, l+1:100_000), view(x, 1:100_000-l)) + dot(view(x, 100_000+l+1:n), view(x, 100_001:n-l)) for l in 0:L]
        merge!(a1, a2)
        @test value(a1).acov ≈ S ./ (n .- (0:L)) .- ref_mean(x)^2 rtol = 1e-10
        @test_throws ArgumentError merge!(AutoCorrelation(3), AutoCorrelation(4))
    end

    @testset "Sampled" begin
        x = randn(rng, 10_001)
        s = Sampled(Moments(), 100); foreach(v -> fit!(s, v), x)
        @test value(s).n == 101
        @test value(s).mean ≈ ref_mean(x[1:100:end]) rtol = 1e-12
        @test value(s).var ≈ ref_var(x[1:100:end]) rtol = 1e-12
        s2 = Sampled(Moments(), 100); foreach(v -> fit!(s2, v), x)
        @test value(merge!(s, s2)).n == 202 && s.count == 20_002
    end

    @testset "MovingAverage" begin
        x = randn(rng, 5_000)
        f = MovingAverage(50)
        out = [smooth!(f, v) for v in x]
        ref = [ref_mean(x[max(1, i - 49):i]) for i in 1:5_000]
        @test out ≈ ref rtol = 1e-12
    end

    @testset "SchmittTrigger" begin
        # hand-built series, Δt = 0.1: reversal (one crossing), reversal, excursion, aborted, dip
        a = [1, 1, 0.2, -0.2, -1, -1, -0.1, 0.4, 1, 1, 0.1, 0.2, 1, 0.3, -0.3, 0.3, 1, 0.4, 1]
        d = SchmittTrigger(; Δt = 0.1); foreach(v -> fit!(d, v), a)
        v = value(d)
        @test (v.n_reversals, v.n_aborted, v.n_excursions) == (2, 1, 1)
        @test v.time ≈ 1.9 && v.mean_interval ≈ 0.95
        @test v.reversal_times ≈ [0.4, 0.8] && v.intervals ≈ [0.4]
        @test v.durations ≈ [0.2, 0.2]
        @test [e.kind for e in v.events] == [:reversal, :reversal, :excursion, :aborted]
        @test v.events[3].amin ≈ 0.1 && v.events[4].crossings == 2 && isnan(v.events[3].tcross)
        # a series that starts in the negative state keeps that polarity
        d2 = SchmittTrigger(; Δt = 1.0); foreach(v -> fit!(d2, v), [-1, -1, -0.2, 0.9, 1])
        @test value(d2).n_reversals == 1 && value(d2).reversal_times == [4.0]
        # random series: the reversal count equals the number of offline polarity flips
        x = cumsum(0.3 .* randn(rng, 50_000)); x = clamp.(x, -3, 3)
        d3 = SchmittTrigger(; Δt = 0.01); foreach(v -> fit!(d3, v), x)
        @test value(d3).n_reversals == length(ref_flips(x, 0.5))
        @test all(e.kind == :reversal ? e.crossings ≥ 1 : true for e in value(d3).events)
        # a jump across the band within one sample is a reversal of zero duration
        d5 = SchmittTrigger(; Δt = 1.0); foreach(v -> fit!(d5, v), [1, 1, -1, -1, 0.2, 1])
        @test value(d5).n_reversals == 2 && value(d5).durations == [0.0, 1.0] && value(d5).reversal_times == [3.0, 5.0]
        # entering the band with the opposite sign counts as a crossing
        d6 = SchmittTrigger(; Δt = 1.0); foreach(v -> fit!(d6, v), [1, -0.2, 0.2, 1, -0.2, -1])
        @test value(d6).n_aborted == 1 && value(d6).events[1].crossings == 2
        @test value(d6).n_reversals == 1 && value(d6).events[2].crossings == 1 && value(d6).reversal_times == [5.0]
        @test all(e.kind == :aborted ? iseven(e.crossings) && e.crossings ≥ 2 : true for e in value(d3).events)
        @test all(e.kind == :excursion ? e.crossings == 0 && e.amin < 0.3 : true for e in value(d3).events)
        d4 = SchmittTrigger(; Δt = 0.01); foreach(v -> fit!(d4, v), x)
        @test value(merge!(d3, d4)).n_reversals == 2 * length(ref_flips(x, 0.5))
        @test_throws ArgumentError SchmittTrigger(; ah = 0.3, ae = 0.5)
    end

    @testset "CrossingGroups" begin
        # Δt = 0.1, W = 0.5: one crossing (reversal), two within W (aborted), three (reversal)
        s = ones(80)
        s[10:29] .= -1; s[32:49] .= -1; s[52:53] .= -1
        c = CrossingGroups(; W = 0.5, Δt = 0.1); foreach(v -> fit!(c, v), s)
        v = value(c)
        @test (v.n_reversals, v.n_aborted, v.pending) == (2, 1, 0)
        @test v.reversal_times ≈ [1.0, 5.0]
        @test v.time ≈ 8.0 && v.mean_interval ≈ 4.0
        # an open group is reported as pending
        c2 = CrossingGroups(; W = 0.5, Δt = 0.1); foreach(v -> fit!(c2, v), s[1:12])
        @test value(c2).pending == 1 && value(c2).n_reversals == 0
        @test value(merge!(c, c2)).n_reversals == 2
    end

    @testset "ExclusionVariance" begin
        # oscillation about ±1 with imposed reversals through the band
        n, Δt, W = 20_000, 0.1, 5.0
        a = 1 .+ 0.2 .* sin.(0.3 .* (1:n)) .+ 0.05 .* randn(rng, n)
        for r in (4_000, 9_000, 9_600, 15_000)
            a[r:end] .*= -1
            a[r-3:r+3] .*= [0.6, 0.3, 0.1, 0.0, 0.1, 0.3, 0.6]       # pass through the band
        end
        ev = ExclusionVariance(; ah = 0.5, W, Δt); foreach(v -> fit!(ev, v), a)
        ref = ref_exclusion(a, 0.5, W, Δt)
        v = value(ev)
        @test v.intervals == ref.intervals == 5
        @test v.n == ref.n
        @test v.variance ≈ ref.variance rtol = 1e-10
        # no reversal at all: the whole series minus the trailing W
        b = 1 .+ 0.1 .* randn(rng, 2_000)
        ev2 = ExclusionVariance(; W, Δt); foreach(v -> fit!(ev2, v), b)
        @test value(ev2).n == 1_950
        @test value(ev2).variance ≈ ref_var(b[1:1950]) / ref_mean(b[1:1950])^2 rtol = 1e-10
        # merge pools the moments by kept-sample count (not the per-interval ratios)
        ev3 = ExclusionVariance(; W, Δt); foreach(v -> fit!(ev3, v), b)
        merge!(ev, ev3)
        @test value(ev).variance ≈ ref_exclusion_pool(ref, ref_exclusion(b, 0.5, W, Δt)) rtol = 1e-10
        @test value(ev).intervals == 6
        # an interval hovering near zero does not dominate: pooled moments, not pooled ratios
        w = round(Int, W / Δt)
        c = copy(a); c[4_000+w+1:9_000-w-1] .= 0.02 .* randn(rng, 4_999 - 2w)
        ev4 = ExclusionVariance(; W, Δt); foreach(v -> fit!(ev4, v), c)
        @test value(ev4).variance ≈ ref_exclusion(c, 0.5, W, Δt).variance rtol = 1e-10
        @test value(ev4).variance < 2 * ref.variance
    end

    @testset "BinnedMean and even_quadratic_fit" begin
        x = 1.5 .* rand(rng, 50_000)
        y = 1.36 .+ 2.0 .* x .^ 2
        b = BinnedMean(0.0, 1.5, 30); foreach(xy -> fit!(b, xy), zip(x, y))
        v = value(b)
        for i in 1:30
            sel = (v.centers[i] - 0.025 .≤ x .< v.centers[i] + 0.025)
            @test v.n[i] == count(sel)
            @test v.mean[i] ≈ ref_mean(y[sel]) rtol = 1e-12
            @test v.var[i] ≈ ref_var(y[sel]) rtol = 1e-8
        end
        keep = v.centers .≤ 1.1
        c0, c2 = even_quadratic_fit(v.centers[keep], v.mean[keep], v.n[keep])
        @test abs(c0 - 1.36) < 0.01 && abs(c2 - 2.0) < 0.01
        fit!(b, (-0.1, 0.0)); fit!(b, (1.5, 0.0))
        @test value(b).below == 1 && value(b).above == 1
        b2 = BinnedMean(0.0, 1.5, 30); foreach(xy -> fit!(b2, xy), zip(x, y))
        @test value(merge!(b, b2)).n == 2 .* v.n
        @test value(b).mean ≈ v.mean rtol = 1e-12
    end

    @testset "Tuple composition" begin
        aggs = (Moments(), Histogram(-3.0, 3.0, 60), Sampled(Moments(), 10))
        x = randn(rng, 1_000)
        foreach(v -> fit!(aggs, v), x)
        v = value(aggs)
        @test v[1].n == 1_000 && v[2].n == 1_000 && v[3].n == 100
        aggs2 = (Moments(), Histogram(-3.0, 3.0, 60), Sampled(Moments(), 10))
        foreach(v -> fit!(aggs2, v), x)
        merge!(aggs, aggs2)
        @test value(aggs)[1].n == 2_000 && value(aggs)[3].n == 200
    end
end
