"""HOS Fig. 6 (p. 283): ℜa₁ against a = ℜa₀ along a trajectory through a reversal, f = 5.4.

The paper's window is "near τ = 680"; here the reversal nearest `--tau` (default 680) in the
run file is taken, with `--half` τ_d (default 50) on either side, at the a₁ stream's 5-step
spacing. Expected: two clouds centred at a ≈ ±1 with |ℜa₁| ≲ 0.5 and one jump between them.
Gate T5.6.
"""

import argparse

import numpy as np
import matplotlib.pyplot as plt

from runfile import RunFile
from style import FIGURES, RESULTS, savefig

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("run", nargs="?", default=str(RESULTS / "standard_f5.4_seed1.h5"), help="run file at f = 5.4")
p.add_argument("--tau", type=float, default=680.0, help="take the reversal nearest this time")
p.add_argument("--half", type=float, default=50.0, help="half-width of the window in τ_d")
p.add_argument("-o", "--out", default=str(FIGURES / "hos_fig06.png"))
args = p.parse_args()

with RunFile(args.run) as run:
    rev = run.reversal_times()
    tr = rev[np.argmin(np.abs(rev - args.tau))]
    t1, a1 = run.window("a1", tr - args.half, tr + args.half)
    a_all, every_a, dt_a = run.stream("a")
    every_1 = int(run.h5["streams/a1"].attrs["every"])
    a = a_all[(np.rint(t1 / dt_a)).astype(int) - 1]       # a at the a₁ sample steps
    f = float(run.meta["schedule_f"])
    seed = int(run.meta["seed"])

fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(a, a1.real, color="k", lw=0.4, alpha=0.8)
ax.plot(a, a1.real, ".", color="C0", ms=1.5)
ax.axhline(0, color="0.6", lw=0.4)
ax.axvline(0, color="0.6", lw=0.4)
ax.set_xlim(-2, 2)
ax.set_ylim(-1, 1)
ax.set_xlabel("a = ℜa₀")
ax.set_ylabel("ℜa₁")
ax.set_title(f"HOS Fig. 6 — f = {f:g}, seed {seed}, τ = {tr - args.half:.0f}–{tr + args.half:.0f} (reversal at {tr:.1f}), {every_1}-step spacing")
savefig(fig, args.out)
