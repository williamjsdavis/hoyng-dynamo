"""HOS Fig. 2 (p. 275): |a₁|, |a₂|, parity 𝒫 and the toroidal energy over τ = 500–700.

Four panels at the paper's 5-step spacing (the n_c cadence of the a₁, a₂, parity and energy
streams), with a bar at every reversal. Gate T5.2.
"""

import numpy as np
import matplotlib.pyplot as plt

from runfile import RunFile
from style import parse_args, reversal_bars, savefig

args = parse_args(__doc__, "hos_fig02.png")
T0, T1 = 500.0, 700.0

with RunFile(args.run) as run:
    t1, a1 = run.window("a1", T0, T1)
    t2, a2 = run.window("a2", T0, T1)
    tp, P = run.window("parity", T0, T1)
    te, E = run.window("energy", T0, T1)
    rev = run.reversal_times()
    meta = run.meta

fig, axes = plt.subplots(4, 1, figsize=(8, 6), sharex=True, gridspec_kw={"hspace": 0.1})
panels = [(t1, np.abs(a1), r"$|a_1|$", (0, 1.5)), (t2, np.abs(a2), r"$|a_2|$", (0, 1.5)),
          (tp, P, r"parity $\mathcal{P}$", (-1.05, 1.05)), (te, E, r"energy $\int T^2 \sin\theta\, d\theta$", (0, None))]
for ax, (t, y, label, ylim) in zip(axes, panels):
    ax.plot(t, y, color="k", lw=0.6)
    reversal_bars(ax, rev, T0, T1, color="C3", lw=0.7)
    ax.set_ylabel(label)
    ax.set_ylim(*ylim)
axes[2].axhline(-1, color="0.6", lw=0.4)
axes[-1].set_xlabel("τ (diffusion times)")
axes[-1].set_xlim(T0, T1)
fig.suptitle(f"HOS Fig. 2 — standard run (seed {meta['seed']}), τ = {T0:.0f}–{T1:.0f}, 5-step spacing", y=0.93)
savefig(fig, args.out)
