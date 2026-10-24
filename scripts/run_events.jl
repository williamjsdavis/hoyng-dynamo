# HOS Figs. 3 and 4 (pp. 278–279, Colour Plates II–III): four events of the standard run at
# full resolution over 6 τ_d windows, re-run from the seed. From the event list and the
# streams of a run file written by run_standard.jl the script selects
#
#   fast        the single-crossing reversal with the shortest |a| < a_h stretch among those
#               with no other event inside the window           (Fig. 3 top: a very fast event)
#   quadrupole  among the reversals with a contiguous stretch of ≥ 1 τ_d below |a| = 0.3, the
#               one with the largest mean |a₁| over that stretch (Fig. 3 bottom: quasi-periodic,
#               the quadrupole's period 2π/ℑλ¹ ≈ 1.1 τ_d)
#   aborted     the two-crossing aborted reversal in which a reaches farthest into the opposite
#               polarity                                        (Fig. 4 top)
#   excursion   the excursion with the smallest min |a|         (Fig. 4 bottom)
#
# and re-runs the simulation from τ = 0 to each window with `window_run` (the seed fixes the
# trajectory; the window's a is checked bitwise against the recorded stream). Output:
# results/events_<run>.h5 with /windows/<name>/{time, T, P, a, a1, a2, parity, energy} and the
# event's attributes; viz/hos_fig03.py and viz/hos_fig04.py plot it. About 40 s.
#
#     julia --project scripts/run_events.jl [results/standard_f6.4_seed1.h5]

using HoyngDynamo
using HDF5
import Dates

path = length(ARGS) ≥ 1 ? ARGS[1] : normpath(joinpath(@__DIR__, "..", "results", "standard_f6.4_seed1.h5"))
file = open_run(path)
spec = RunSpec(file)
md = metadata(file)
ev = read_results(file)["trigger"]["events"]
a = read_stream(file, "a")
a1 = read_stream(file, "a1")
close(file)

kind, start, stop, amin, crossings, tcross = ev["kind"], ev["start"], ev["stop"], ev["amin"], ev["crossings"], ev["tcross"]
nev = length(kind)

"The samples of a stream with t0 ≤ τ ≤ t1 (sample i was taken at τ = i·dt)."
function segment(stream, t0, t1)
    i0 = max(ceil(Int, t0 / stream.dt - 1e-9), 1)
    i1 = min(floor(Int, t1 / stream.dt + 1e-9), length(stream.values))
    return view(stream.values, i0:i1)
end

"A 6 τ_d window on integer τ_d boundaries around an event, as HOS's τ = 219–225."
function window(start, stop)
    τ0 = floor((start + stop) / 2) - 3
    start - τ0 < 0.5 && (τ0 = floor(start - 0.5))
    return τ0, τ0 + 6
end

function clean(i)
    τ0, τ1 = window(start[i], stop[i])
    return !any(j != i && start[j] ≤ τ1 && stop[j] ≥ τ0 for j in 1:nev)
end

rev = [i for i in 1:nev if kind[i] == "reversal"]

# fast: single-crossing reversals, shortest first, the first with a clean window
fast_candidates = sort([i for i in rev if crossings[i] == 1]; by = i -> stop[i] - start[i])
fast = fast_candidates[findfirst(clean, fast_candidates)]

# quadrupole: the longest contiguous |a| < 0.3 stretch inside each reversal
function low_stretch(i; level = 0.3)
    s = segment(a, start[i], stop[i])
    best = cur = bend = 0
    for (k, v) in enumerate(s)
        cur = abs(v) < level ? cur + 1 : 0
        cur > best && (best = cur; bend = k)
    end
    i0 = max(ceil(Int, start[i] / a.dt - 1e-9), 1)
    t1 = (i0 + bend - 1) * a.dt
    return best * a.dt, t1 - (best - 1) * a.dt, t1
end
quad_candidates = Tuple{Int,Float64}[]
for i in rev
    L, t0, t1 = low_stretch(i)
    L ≥ 1 || continue
    s = segment(a1, t0, t1)
    push!(quad_candidates, (i, sum(abs, s) / length(s)))
end
quadrupole = quad_candidates[argmax(last.(quad_candidates))][1]

# aborted: two crossings, the farthest reach of a into the opposite polarity
aborted_candidates = [i for i in 1:nev if kind[i] == "aborted" && crossings[i] == 2]
function reach(i)
    pol = sign(a.values[max(round(Int, start[i] / a.dt) - 1, 1)])     # the sample before the event
    return maximum(v -> -pol * v, segment(a, start[i], stop[i]))
end
aborted = aborted_candidates[argmax(reach.(aborted_candidates))]

# excursion: the closest approach to zero without a crossing
excursion_candidates = [i for i in 1:nev if kind[i] == "excursion"]
excursion = excursion_candidates[argmin(amin[excursion_candidates])]

selected = (fast = fast, quadrupole = quadrupole, aborted = aborted, excursion = excursion)
rules = (fast = "single-crossing reversal with the shortest |a| < a_h stretch and no other event in the window",
         quadrupole = "reversal with ≥ 1 τ_d below |a| = 0.3 and the largest mean |a₁| over that stretch",
         aborted = "two-crossing aborted reversal reaching farthest into the opposite polarity",
         excursion = "excursion with the smallest min |a|")

out = normpath(joinpath(@__DIR__, "..", "results", replace(basename(path), r"^[A-Za-z]+_" => "events_")))
sim0 = Simulation(spec)
println(spec)
h5open(out, "w") do h5
    at = attributes(h5)
    for (k, v) in md
        at[k == "created" ? "source_created" : k] = v        # the run's metadata; our own creation date below
    end
    at["source"] = path
    at["Lambda"] = sim0.model.constants.Λ
    at["mu"] = sim0.model.constants.μ
    at["created"] = string(Dates.now())
    h5["grid/theta"] = collect(sim0.model.grid.θ)
    for (name, i) in pairs(selected)
        τ0, τ1 = window(start[i], stop[i])
        t = @elapsed w = window_run(spec, τ0, τ1)
        n0 = round(Int, τ0 / a.dt)
        w.a == a.values[n0+1:n0+length(w.a)] || error("the re-run window $name does not reproduce the recorded a stream")
        g = create_group(h5, "windows/$name")
        g["time"] = w.time
        g["T"] = w.T
        g["P"] = w.P
        g["a"] = w.a
        g["a1"] = w.a1
        g["a2"] = w.a2
        g["parity"] = w.parity
        g["energy"] = w.energy
        ga = attributes(g)
        ga["kind"] = kind[i]
        ga["start"] = start[i]
        ga["stop"] = stop[i]
        ga["tcross"] = tcross[i]
        ga["crossings"] = crossings[i]
        ga["amin"] = amin[i]
        ga["tau0"] = τ0
        ga["tau1"] = τ1
        ga["every"] = 1
        ga["dt"] = a.dt
        ga["rule"] = rules[name]
        println("$name: $(kind[i]) at τ = $(round(start[i], digits = 2))–$(round(stop[i], digits = 2)) ($(crossings[i]) crossings, min |a| = $(round(amin[i], digits = 4))), " *
                "window $(τ0)–$(τ1), re-run in $(round(t, digits = 1)) s, bitwise with the recorded stream")
    end
end
println("wrote $out ($(round(filesize(out) / 2^20, digits = 1)) MB)")
