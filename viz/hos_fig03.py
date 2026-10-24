"""HOS Fig. 3 (p. 278, Colour Plate II): two reversals at full resolution over 6 τ_d windows.

Top: a very fast event; bottom: a quasi-periodic, quadrupolar one. Each shows T(θ, τ) at every
step with the zero contour in black, and below it a(τ), |a₁| and the parity 𝒫 with the band
±a_h dotted and an exp(Λτ) guide after the last zero crossing (the recovery of the paper's
caption). The windows come from results/events_f6.4_seed1.h5 (scripts/run_events.jl, which
documents the selection rules). Gate T5.3.
"""

from eventwindows import default_events, two_window_figure
from style import parse_args

args = parse_args(__doc__, "hos_fig03.png", default_events())
two_window_figure(args.run, ["fast", "quadrupole"],
                  "HOS Fig. 3 — standard run (f = {schedule_f}, seed {seed}): a very fast reversal and a quadrupolar one, full resolution",
                  args.out, guide=True)
