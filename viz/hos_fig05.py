"""HOS Fig. 5 (p. 280): a(τ) over 0–2000 τ_d for f = 5.4 (top) and f = 7.4 (bottom).

The same random series as Fig. 1 (seed 1: the F stream does not depend on f), with a bar at
every reversal. The paper's panels show 6 reversal bars at f = 5.4 and several tens at
f = 7.4. Gate T5.5.

    uv run python hos_fig05.py [results/standard_f5.4_seed1.h5] [results/standard_f7.4_seed1.h5]
"""

import argparse

import matplotlib.pyplot as plt

from runfile import RunFile
from style import FIGURES, RESULTS, reversal_bars, savefig

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("low", nargs="?", default=str(RESULTS / "standard_f5.4_seed1.h5"), help="run file at f = 5.4")
p.add_argument("high", nargs="?", default=str(RESULTS / "standard_f7.4_seed1.h5"), help="run file at f = 7.4")
p.add_argument("-o", "--out", default=str(FIGURES / "hos_fig05.png"))
args = p.parse_args()
T0, T1 = 0.0, 2000.0

fig, axes = plt.subplots(2, 1, figsize=(8, 4.2), sharex=True, gridspec_kw={"hspace": 0.15})
for ax, path in zip(axes, (args.low, args.high)):
    with RunFile(path) as run:
        t, a = run.window("a", T0, T1)
        rev = run.reversal_times()
        f = float(run.meta["schedule_f"])
        seed = int(run.meta["seed"])
    nrev = int(((rev >= T0) & (rev <= T1)).sum())
    ax.plot(t, a, color="k", lw=0.4)
    reversal_bars(ax, rev, T0, T1, color="C3", lw=0.7)
    ax.axhline(0, color="0.6", lw=0.4)
    ax.set_ylim(-2, 2)
    ax.set_ylabel("a")
    ax.text(0.01, 0.95, f"f = {f:g}: {nrev} reversals in {T1 - T0:.0f} τ_d (seed {seed})", transform=ax.transAxes, va="top", fontsize=8,
            bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none", "pad": 1.5})
axes[-1].set_xlabel("τ (diffusion times)")
axes[-1].set_xlim(T0, T1)
fig.suptitle("HOS Fig. 5 — the same noise series at f = 5.4 (top) and f = 7.4 (bottom)", y=0.95)
savefig(fig, args.out)
