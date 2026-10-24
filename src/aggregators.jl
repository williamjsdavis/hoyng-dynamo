# Running aggregators (dev/DECISIONS.md D8): constant-memory statistics for arbitrarily long
# runs, behind three verbs:
#
#     fit!(agg, x)        feed one observation (a Real; a Tuple (x, y) for BinnedMean)
#     merge!(a, b)        pool b into a (Base.merge!): independent seeds or ensembles
#     value(agg)          a NamedTuple of results, computed on demand
#
# Bespoke implementation, no OnlineStats.jl (D8 option b): the generic statistics are a few
# lines each, the model-specific ones (event detector, exclusion variance, D₁ bins) are
# bespoke in any case, and every method is kept allocation-free and type-stable (D14,
# test/quality.jl). Aggregators are composed in Tuples, never in a Vector of an abstract
# type. Cadence (every k-th sample) is the `Sampled` wrapper, by an internal counter, so
# replaying a recorded series reproduces the online values exactly (D7).
#
# Precision over 10¹⁰ samples: Welford updates for means and variances, Int64 counters,
# chunked accumulation for the lagged sums of the autocorrelation.
#
# Times are in units of τ_d: each time-aware aggregator takes the sampling interval Δt of
# the series it is fed (Δτ for every-step series, 1 for once-per-τ_d series).
# Merging pools independent runs; restart segments of one run are resumed from the
# checkpointed aggregator state (D6), not merged.

abstract type Aggregator end

"""
    fit!(agg, x)

Feed one observation to `agg` (or to every aggregator of a `Tuple`). Returns `agg`.
"""
function fit! end

"""
    value(agg)

The results of `agg` as a NamedTuple (or a Tuple of them for a `Tuple` of aggregators).
"""
function value end

fit!(t::Tuple, x) = (foreach(a -> fit!(a, x), t); t)
Base.merge!(a::Tuple, b::Tuple) = (foreach(merge!, a, b); a)
value(t::Tuple) = map(value, t)

# ---------------------------------------------------------------- Moments (Welford)

"""
    Moments()

Running mean and variance by Welford's algorithm (O(1) memory). `value` gives
`(n, mean, var, std)` with the sample variance (n − 1).
"""
mutable struct Moments{T<:AbstractFloat} <: Aggregator
    n::Int
    mean::T
    M2::T
end
Moments(::Type{T} = Float64) where {T} = Moments{T}(0, zero(T), zero(T))

@inline function fit!(m::Moments{T}, x::Real) where {T}
    xT = T(x)
    m.n += 1
    δ = xT - m.mean
    m.mean += δ / m.n
    m.M2 += δ * (xT - m.mean)
    return m
end

function Base.merge!(a::Moments{T}, b::Moments{T}) where {T}
    na, nb = a.n, b.n
    n = na + nb
    n == 0 && return a
    δ = b.mean - a.mean
    a.mean += δ * (T(nb) / n)
    a.M2 += b.M2 + δ^2 * (T(na) * T(nb) / n)
    a.n = n
    return a
end

function value(m::Moments{T}) where {T}
    var = m.n > 1 ? m.M2 / (m.n - 1) : T(NaN)
    return (n = m.n, mean = m.mean, var = var, std = sqrt(var))
end

# ---------------------------------------------------------------- Histogram (fixed bins)

"""
    Histogram(lo, hi, nbins)

Counts on `nbins` equal bins of [lo, hi); samples outside are counted in `below`/`above`,
NaN is ignored. `value` gives `(edges, centers, counts, below, above, n)`.
"""
mutable struct Histogram{T<:AbstractFloat} <: Aggregator
    const lo::T
    const hi::T
    const nbins::Int
    const invwidth::T
    const counts::Vector{Int}
    below::Int
    above::Int
end
function Histogram(lo, hi, nbins::Integer)
    T = float(promote_type(typeof(lo), typeof(hi)))
    hi > lo && nbins ≥ 1 || throw(ArgumentError("need hi > lo and nbins ≥ 1"))
    Histogram{T}(T(lo), T(hi), nbins, nbins / (T(hi) - T(lo)), zeros(Int, nbins), 0, 0)
end

@inline function fit!(h::Histogram{T}, x::Real) where {T}
    xT = T(x)
    isnan(xT) && return h
    k = floor(Int, (xT - h.lo) * h.invwidth)
    if 0 ≤ k < h.nbins
        @inbounds h.counts[k+1] += 1
    elseif k < 0
        h.below += 1
    else
        h.above += 1
    end
    return h
end

function Base.merge!(a::Histogram{T}, b::Histogram{T}) where {T}
    (a.lo == b.lo && a.hi == b.hi && a.nbins == b.nbins) || throw(ArgumentError("histograms have different bins"))
    a.counts .+= b.counts
    a.below += b.below
    a.above += b.above
    return a
end

edges(h::Histogram) = range(h.lo, h.hi; length = h.nbins + 1)
centers(h::Histogram) = (e = edges(h); (e[1:end-1] .+ e[2:end]) ./ 2)
value(h::Histogram) = (edges = edges(h), centers = centers(h), counts = copy(h.counts), below = h.below,
                       above = h.above, n = sum(h.counts) + h.below + h.above)

# ---------------------------------------------------------------- AutoCorrelation (fixed max lag)

"""
    AutoCorrelation(L; Δt = 1.0)

Autocovariance of a series at lags 0 … `L` samples (a ring buffer of L samples and L + 1
lagged sums; O(L) memory). `value` gives `(lags, acov, acor, mean, n)`, with lags in time
units (`Δt` per sample) and the mean subtracted globally: acov(l) = Σ x_i x_{i−l}/(n − l) − x̄²,
which differs from the exact mean-subtracted estimate by O(L/n). Lagged sums are
accumulated in chunks of 2²⁰ samples before being added to the totals, so round-off stays
O(√n ε) rather than O(n ε) over very long series. `merge!` adds the sums (exact for
independent runs).
"""
mutable struct AutoCorrelation{T<:AbstractFloat} <: Aggregator
    const L::Int
    const Δt::T
    const buf::Vector{T}
    const S::Vector{T}        # chunk lagged sums, S[l+1] for lag l
    const Stot::Vector{T}     # flushed totals
    sum::T
    sumtot::T
    n::Int
    pos::Int                  # ring index where the next sample is stored
end
const AUTOCORR_CHUNK = 1 << 20
function AutoCorrelation(L::Integer; Δt = 1.0, T::Type{<:AbstractFloat} = Float64)
    L ≥ 1 || throw(ArgumentError("need L ≥ 1"))
    AutoCorrelation{T}(L, T(Δt), zeros(T, L), zeros(T, L + 1), zeros(T, L + 1), zero(T), zero(T), 0, 1)
end

@inline function fit!(ac::AutoCorrelation{T}, x::Real) where {T}
    xT = T(x)
    ac.n += 1
    S, buf, L = ac.S, ac.buf, ac.L
    @inbounds S[1] += xT * xT
    j = ac.pos
    @inbounds for l in 1:min(L, ac.n - 1)
        j -= 1
        j == 0 && (j = L)
        S[l+1] += xT * buf[j]
    end
    ac.sum += xT
    @inbounds buf[ac.pos] = xT
    ac.pos = ac.pos == L ? 1 : ac.pos + 1
    if ac.n % AUTOCORR_CHUNK == 0
        ac.Stot .+= S
        fill!(S, zero(T))
        ac.sumtot += ac.sum
        ac.sum = zero(T)
    end
    return ac
end

function Base.merge!(a::AutoCorrelation{T}, b::AutoCorrelation{T}) where {T}
    a.L == b.L || throw(ArgumentError("autocorrelations have different maximum lags"))
    a.Stot .+= b.Stot .+ b.S
    a.sumtot += b.sumtot + b.sum
    a.n += b.n
    return a
end

function value(ac::AutoCorrelation{T}) where {T}
    n = ac.n
    S = ac.Stot .+ ac.S
    mean = (ac.sumtot + ac.sum) / n
    acov = [l < n ? S[l+1] / (n - l) - mean^2 : T(NaN) for l in 0:ac.L]
    return (lags = (0:ac.L) .* ac.Δt, acov = acov, acor = acov ./ acov[1], mean = mean, n = n)
end

# ---------------------------------------------------------------- Sampled (cadence)

"""
    Sampled(agg, every)

Feed `agg` the first sample and every `every`-th sample after it (an internal counter).
"""
mutable struct Sampled{A<:Aggregator} <: Aggregator
    const agg::A
    const every::Int
    count::Int
end
Sampled(agg::Aggregator, every::Integer) = (every ≥ 1 || throw(ArgumentError("need every ≥ 1")); Sampled(agg, Int(every), 0))

@inline function fit!(s::Sampled, x)
    s.count += 1
    (s.count - 1) % s.every == 0 && fit!(s.agg, x)
    return s
end
Base.merge!(a::Sampled, b::Sampled) = (merge!(a.agg, b.agg); a.count += b.count; a)
value(s::Sampled) = value(s.agg)

# ---------------------------------------------------------------- MovingAverage (a stream filter)

"""
    MovingAverage(w)

Running mean of the last `w` samples (partial mean while filling); `smooth!(f, x)` feeds a
sample and returns the current average. A stream transform for feeding an aggregator
(HSO Fig. 7 smooths the rectified amplitude before its autocorrelation), not an aggregator.
"""
mutable struct MovingAverage{T<:AbstractFloat}
    const w::Int
    const buf::Vector{T}
    sum::T
    n::Int
    pos::Int
end
MovingAverage(w::Integer; T::Type{<:AbstractFloat} = Float64) = (w ≥ 1 || throw(ArgumentError("need w ≥ 1")); MovingAverage{T}(w, zeros(T, w), zero(T), 0, 1))

@inline function smooth!(f::MovingAverage{T}, x::Real) where {T}
    xT = T(x)
    f.n += 1
    @inbounds begin
        f.sum += xT - f.buf[f.pos]      # buf holds zeros while filling
        f.buf[f.pos] = xT
    end
    f.pos = f.pos == f.w ? 1 : f.pos + 1
    f.pos == 1 && (f.sum = sum(f.buf))  # exact resum once per window: no drift over long series
    return f.sum / min(f.n, f.w)
end

# ---------------------------------------------------------------- SchmittTrigger (event detector)

"""
    Event(kind, start, stop, amin, crossings)

One event of the `SchmittTrigger`: `kind` is `:reversal`, `:aborted` or `:excursion`;
`start`/`stop` the times |a| dropped below and rose above a_h; `amin` the minimum of |a|
during the event; `crossings` the number of zero crossings; `tcross` the time of the last
zero crossing (NaN if none).
"""
struct Event{T}
    kind::Symbol
    start::T
    stop::T
    amin::T
    crossings::Int
    tcross::T
end

"""
    SchmittTrigger(; ah = 0.5, ae = 0.3, Δt = 1.0)

Reversal detector of `VALIDATION.md` §2. Polarity is that of the first sample outside the
band |a| < a_h (+1 if the first sample is inside it); it flips when a crosses ∓a_h. An
*event* is a stretch with |a| < a_h: if the polarity flipped it is a reversal (its time
is the last zero crossing), otherwise an aborted reversal if a crossed zero, or an
excursion if |a| dipped below a_e. A jump from one side of the band to the other within
one sample (possible for a series sampled once per τ_d) is a reversal of zero duration.
Events are kept in a list (O(number of events));
everything else is O(1). `value` gives counts, the total time, ⟨T_r⟩ = time/reversals,
the reversal times and the polarity intervals between them.
"""
mutable struct SchmittTrigger{T<:AbstractFloat} <: Aggregator
    const ah::T
    const ae::T
    const Δt::T
    n::Int
    polarity::Int           # 0 until the first sample outside the band
    inevent::Bool
    start::Int              # sample index at which the event began
    amin::T
    crossings::Int
    lastcross::Int
    lastsign::Int
    nreversals::Int
    naborted::Int
    nexcursions::Int
    const events::Vector{Event{T}}
end
function SchmittTrigger(; ah = 0.5, ae = 0.3, Δt = 1.0, T::Type{<:AbstractFloat} = Float64)
    0 < ae < ah || throw(ArgumentError("need 0 < a_e < a_h"))
    SchmittTrigger{T}(T(ah), T(ae), T(Δt), 0, 0, false, 0, zero(T), 0, 0, 0, 0, 0, 0, Event{T}[])
end

@inline function fit!(d::SchmittTrigger{T}, x::Real) where {T}
    a = T(x)
    d.n += 1
    s = a > 0 ? 1 : (a < 0 ? -1 : 0)
    outside = abs(a) ≥ d.ah
    if d.polarity == 0                   # before the first sample outside the band
        outside && (d.polarity = s)
        d.lastsign = s == 0 ? d.lastsign : s
        return d
    end
    if d.inevent
        # zero crossings and the minimum of |a| while inside the band
        if s != 0 && d.lastsign != 0 && s != d.lastsign
            d.crossings += 1
            d.lastcross = d.n
        end
        abs(a) < d.amin && (d.amin = abs(a))
        if outside                        # the event ends: classify it
            d.inevent = false
            if s != d.polarity
                kind = :reversal
                d.polarity = s
                d.nreversals += 1
            elseif d.crossings > 0
                kind = :aborted
                d.naborted += 1
            elseif d.amin < d.ae
                kind = :excursion
                d.nexcursions += 1
            else
                kind = :dip
            end
            if kind != :dip
                tcross = d.crossings > 0 ? d.lastcross * d.Δt : T(NaN)
                push!(d.events, Event{T}(kind, d.start * d.Δt, d.n * d.Δt, d.amin, d.crossings, tcross))
            end
        end
    elseif !outside                       # an event begins, possibly with a crossing
        d.inevent = true
        d.start = d.n
        d.amin = abs(a)
        crossed = s != 0 && d.lastsign != 0 && s != d.lastsign
        d.crossings = crossed ? 1 : 0
        d.lastcross = crossed ? d.n : 0
    elseif s != d.polarity                # a jumped across the band in one sample
        d.polarity = s
        d.nreversals += 1
        push!(d.events, Event{T}(:reversal, d.n * d.Δt, d.n * d.Δt, abs(a), 1, d.n * d.Δt))
    end
    s == 0 || (d.lastsign = s)
    return d
end

function Base.merge!(a::SchmittTrigger{T}, b::SchmittTrigger{T}) where {T}
    a.n += b.n
    a.nreversals += b.nreversals
    a.naborted += b.naborted
    a.nexcursions += b.nexcursions
    append!(a.events, b.events)
    return a
end

function value(d::SchmittTrigger{T}) where {T}
    rt = [e.tcross for e in d.events if e.kind == :reversal]
    total = d.n * d.Δt
    return (n_reversals = d.nreversals, n_aborted = d.naborted, n_excursions = d.nexcursions, time = total,
            mean_interval = total / d.nreversals, reversal_times = rt, intervals = diff(rt),
            durations = [e.stop - e.start for e in d.events if e.kind == :reversal], events = copy(d.events))
end

# ---------------------------------------------------------------- Series (stored samples)

"""
    Series(; T = Float64, sizehint = 0)

Stores every sample fed to it, for diagnostics that need the whole series after the fact,
such as the parity outside the event windows of a `SchmittTrigger` (gate T4.13). Memory
grows with the run, so it serves finite diagnostics, not the open-ended statistics of a long
run; `fit!` allocates only when the vector grows past its capacity (reserve it with
`sizehint`). `merge!` appends; `value` gives `(n, values)`.
"""
struct Series{T} <: Aggregator
    values::Vector{T}
end
function Series(; T::Type = Float64, sizehint::Integer = 0)
    v = T[]
    sizehint!(v, sizehint)
    return Series{T}(v)
end
@inline fit!(s::Series{T}, x) where {T} = (push!(s.values, convert(T, x)); s)
Base.merge!(a::Series{T}, b::Series{T}) where {T} = (append!(a.values, b.values); a)
value(s::Series) = (n = length(s.values), values = copy(s.values))

# ---------------------------------------------------------------- event windows

"""
    outside_events(times, events; pad = 0) -> BitVector

`true` at every sample time that lies outside every window `[start − pad, stop + pad]` of
`events` (the event list of a `SchmittTrigger`'s value, all kinds; the window is the stretch
with |a| < a_h plus the sample that ended it). `times` must be increasing. HOS p. 276:
"except during (aborted) reversals and excursions 𝒫 ≃ −1" (gate T4.13).
"""
function outside_events(times::AbstractVector{<:Real}, events::AbstractVector{<:Event}; pad::Real = 0)
    pad ≥ 0 || throw(ArgumentError("pad must be ≥ 0, got $pad"))
    issorted(times) || throw(ArgumentError("the sample times must be increasing"))
    mask = trues(length(times))
    for e in events
        i1 = searchsortedfirst(times, e.start - pad)
        i2 = searchsortedlast(times, e.stop + pad)
        i1 ≤ i2 && (mask[i1:i2] .= false)
    end
    return mask
end

# ---------------------------------------------------------------- CrossingGroups (the papers' rule)

"""
    CrossingGroups(; W = 5.0, Δt = 1.0)

The papers' reversal count (HOS footnote 2, p. 276; `DEVELOPMENT.md` G4): zero crossings
of a that follow one another within `W` are grouped; a group with an odd number of
crossings is a reversal (timed at its first crossing), an even number an aborted
reversal. A group closes `W` after its last crossing, so `value` also reports the
crossings of a still-open group as `pending`.
"""
mutable struct CrossingGroups{T<:AbstractFloat} <: Aggregator
    const W::T
    const Δt::T
    const wn::Int             # samples in W
    n::Int
    lastsign::Int
    count::Int                # crossings in the open group (0: none open)
    first::Int                # sample index of the group's first crossing
    last::Int                 # sample index of the group's last crossing
    nreversals::Int
    naborted::Int
    const reversal_times::Vector{T}
end
CrossingGroups(; W = 5.0, Δt = 1.0, T::Type{<:AbstractFloat} = Float64) =
    CrossingGroups{T}(T(W), T(Δt), max(1, round(Int, W / Δt)), 0, 0, 0, 0, 0, 0, 0, T[])

@inline function close_group!(c::CrossingGroups{T}) where {T}
    if isodd(c.count)
        c.nreversals += 1
        push!(c.reversal_times, c.first * c.Δt)
    else
        c.naborted += 1
    end
    c.count = 0
    return c
end

@inline function fit!(c::CrossingGroups{T}, x::Real) where {T}
    a = T(x)
    c.n += 1
    s = a > 0 ? 1 : (a < 0 ? -1 : 0)
    if c.count > 0 && c.n - c.last > c.wn
        close_group!(c)
    end
    if s != 0 && c.lastsign != 0 && s != c.lastsign
        if c.count == 0
            c.first = c.n
        end
        c.count += 1
        c.last = c.n
    end
    s == 0 || (c.lastsign = s)
    return c
end

function Base.merge!(a::CrossingGroups{T}, b::CrossingGroups{T}) where {T}
    a.n += b.n
    a.nreversals += b.nreversals
    a.naborted += b.naborted
    append!(a.reversal_times, b.reversal_times)
    return a
end

value(c::CrossingGroups) = (n_reversals = c.nreversals, n_aborted = c.naborted, pending = c.count, time = c.n * c.Δt,
                            mean_interval = c.n * c.Δt / c.nreversals, reversal_times = copy(c.reversal_times))

# ---------------------------------------------------------------- ExclusionVariance (HOS/SOH style)

"""
    ExclusionVariance(; ah = 0.5, W = 5.0, Δt = 1.0)

Between-reversal relative variance of `VALIDATION.md` §2: within each polarity interval
(Schmitt trigger with threshold `ah`), the samples closer than `W` to either reversal are
dropped; the variance about each interval's own mean and the magnitude of that mean are
pooled over the intervals weighted by their kept samples, and `variance` is the ratio
Σ M2ᵢ/Σ(nᵢ−1) ÷ (Σ nᵢ|⟨a⟩ᵢ|/Σ nᵢ)². (Pooling the per-interval ratios instead lets a few
intervals whose kept part hovers near a = 0, long events the trigger does not resolve,
dominate the result at high forcing; see `DISCREPANCIES.md` DR-6.) A delay ring buffer of
W/Δt samples holds the samples that may still be discarded, so memory is O(W/Δt). `value`
gives `(variance, n, intervals)`, including the open interval.
"""
mutable struct ExclusionVariance{T<:AbstractFloat} <: Aggregator
    const ah::T
    const W::T
    const Δt::T
    const w::Int
    const buf::Vector{T}
    n::Int
    polarity::Int
    skip::Int
    filled::Int
    pos::Int
    cn::Int                   # current interval: Welford state of the kept samples
    cmean::T
    cM2::T
    pooled_M2::T              # Σ M2_i over closed intervals
    pooled_absmean::T         # Σ n_i |mean_i|
    pooled_n::Int             # Σ n_i
    intervals::Int
end
function ExclusionVariance(; ah = 0.5, W = 5.0, Δt = 1.0, T::Type{<:AbstractFloat} = Float64)
    w = max(1, round(Int, W / Δt))
    ExclusionVariance{T}(T(ah), T(W), T(Δt), w, zeros(T, w), 0, 0, 0, 0, 1, 0, zero(T), zero(T), zero(T), zero(T), 0, 0)
end

@inline function keep!(ev::ExclusionVariance{T}, a::T) where {T}
    ev.cn += 1
    δ = a - ev.cmean
    ev.cmean += δ / ev.cn
    ev.cM2 += δ * (a - ev.cmean)
    return ev
end

@inline function close_interval!(ev::ExclusionVariance{T}) where {T}
    if ev.cn ≥ 2
        ev.pooled_M2 += ev.cM2
        ev.pooled_absmean += ev.cn * abs(ev.cmean)
        ev.pooled_n += ev.cn
        ev.intervals += 1
    end
    ev.cn = 0
    ev.cmean = zero(T)
    ev.cM2 = zero(T)
    ev.filled = 0
    ev.pos = 1
    return ev
end

@inline function fit!(ev::ExclusionVariance{T}, x::Real) where {T}
    a = T(x)
    ev.n += 1
    if ev.polarity == 0
        abs(a) ≥ ev.ah && (ev.polarity = a > 0 ? 1 : -1)
    elseif (ev.polarity == 1 && a < -ev.ah) || (ev.polarity == -1 && a > ev.ah)
        close_interval!(ev)             # drops the buffered last W of the interval
        ev.polarity = -ev.polarity
        ev.skip = ev.w                  # and the first W of the next one
    end
    if ev.skip > 0
        ev.skip -= 1
        return ev
    end
    @inbounds begin
        if ev.filled == ev.w
            keep!(ev, ev.buf[ev.pos])   # the sample W ago is now safe to keep
        else
            ev.filled += 1
        end
        ev.buf[ev.pos] = a
    end
    ev.pos = ev.pos == ev.w ? 1 : ev.pos + 1
    return ev
end

function Base.merge!(a::ExclusionVariance{T}, b::ExclusionVariance{T}) where {T}
    # pool b's closed and open intervals into a's closed ones
    a.pooled_M2 += b.pooled_M2
    a.pooled_absmean += b.pooled_absmean
    a.pooled_n += b.pooled_n
    a.intervals += b.intervals
    if b.cn ≥ 2
        a.pooled_M2 += b.cM2
        a.pooled_absmean += b.cn * abs(b.cmean)
        a.pooled_n += b.cn
        a.intervals += 1
    end
    a.n += b.n
    return a
end

function value(ev::ExclusionVariance{T}) where {T}
    M2, absmean, pooled_n, intervals = ev.pooled_M2, ev.pooled_absmean, ev.pooled_n, ev.intervals
    if ev.cn ≥ 2
        M2 += ev.cM2
        absmean += ev.cn * abs(ev.cmean)
        pooled_n += ev.cn
        intervals += 1
    end
    variance = pooled_n > 0 ? (M2 / (pooled_n - intervals)) / (absmean / pooled_n)^2 : T(NaN)
    return (variance = variance, n = pooled_n, intervals = intervals)
end

# ---------------------------------------------------------------- BinnedMean (D₁ versus |a|)

"""
    BinnedMean(lo, hi, nbins)

Running mean and variance of y in fixed bins of x (Welford per bin): `fit!(b, (x, y))`.
For the D₁ measurement of HOS Fig. 8 (bin averages of D₁/ℱ against |a|). `value` gives
`(centers, n, mean, var)`.
"""
mutable struct BinnedMean{T<:AbstractFloat} <: Aggregator
    const lo::T
    const hi::T
    const nbins::Int
    const invwidth::T
    const n::Vector{Int}
    const mean::Vector{T}
    const M2::Vector{T}
    below::Int
    above::Int
end
function BinnedMean(lo, hi, nbins::Integer)
    T = float(promote_type(typeof(lo), typeof(hi)))
    hi > lo && nbins ≥ 1 || throw(ArgumentError("need hi > lo and nbins ≥ 1"))
    BinnedMean{T}(T(lo), T(hi), nbins, nbins / (T(hi) - T(lo)), zeros(Int, nbins), zeros(T, nbins), zeros(T, nbins), 0, 0)
end

@inline function fit!(b::BinnedMean{T}, xy::Tuple{Real,Real}) where {T}
    x, y = T(xy[1]), T(xy[2])
    k = floor(Int, (x - b.lo) * b.invwidth)
    if 0 ≤ k < b.nbins
        @inbounds begin
            i = k + 1
            b.n[i] += 1
            δ = y - b.mean[i]
            b.mean[i] += δ / b.n[i]
            b.M2[i] += δ * (y - b.mean[i])
        end
    elseif k < 0
        b.below += 1
    else
        b.above += 1
    end
    return b
end

function Base.merge!(a::BinnedMean{T}, b::BinnedMean{T}) where {T}
    (a.lo == b.lo && a.hi == b.hi && a.nbins == b.nbins) || throw(ArgumentError("binned means have different bins"))
    @inbounds for i in 1:a.nbins
        na, nb = a.n[i], b.n[i]
        n = na + nb
        n == 0 && continue
        δ = b.mean[i] - a.mean[i]
        a.mean[i] += δ * (T(nb) / n)
        a.M2[i] += b.M2[i] + δ^2 * (T(na) * T(nb) / n)
        a.n[i] = n
    end
    a.below += b.below
    a.above += b.above
    return a
end

function value(b::BinnedMean{T}) where {T}
    e = range(b.lo, b.hi; length = b.nbins + 1)
    var = [n > 1 ? M2 / (n - 1) : T(NaN) for (n, M2) in zip(b.n, b.M2)]
    return (centers = (e[1:end-1] .+ e[2:end]) ./ 2, n = copy(b.n), mean = copy(b.mean), var = var,
            below = b.below, above = b.above)
end

"""
    even_quadratic_fit(x, y, w) -> (c0, c2)

Weighted least-squares fit of y ≈ c0 + c2 x² (HOS Fig. 8: bin averages of D₁/ℱ fitted by a
quadratic in a, weighted by the number of points per bin).
"""
function even_quadratic_fit(x::AbstractVector, y::AbstractVector, w::AbstractVector)
    s0 = s2 = s4 = sy = sy2 = 0.0
    for (xi, yi, wi) in zip(x, y, w)
        wi > 0 || continue
        x2 = xi^2
        s0 += wi; s2 += wi * x2; s4 += wi * x2^2; sy += wi * yi; sy2 += wi * yi * x2
    end
    det = s0 * s4 - s2^2
    c0 = (sy * s4 - sy2 * s2) / det
    c2 = (s0 * sy2 - s2 * sy) / det
    return (c0, c2)
end
