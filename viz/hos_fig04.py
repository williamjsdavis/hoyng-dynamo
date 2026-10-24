"""HOS Fig. 4 (p. 279, Colour Plate III): an aborted reversal and an excursion at full resolution.

Two 6 τ_d windows with T(θ, τ) at every step, the zero contour in black, and a(τ), |a₁|, 𝒫
below: a brief sign change against none, "similar appearance" (HOS p. 277). The windows come
from results/events_f6.4_seed1.h5 (scripts/run_events.jl). Gate T5.4.
"""

from eventwindows import default_events, two_window_figure
from style import parse_args

args = parse_args(__doc__, "hos_fig04.png", default_events())
two_window_figure(args.run, ["aborted", "excursion"],
                  "HOS Fig. 4 — standard run (f = {schedule_f}, seed {seed}): an aborted reversal and an excursion, full resolution",
                  args.out)
