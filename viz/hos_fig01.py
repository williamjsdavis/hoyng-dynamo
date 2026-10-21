"""HOS Fig. 1 (p. 274, Plate I): T(θ, τ) and the dipole amplitude a(τ) over 0–2000 τ_d.

Top: the toroidal field at the nodes, colatitude 0…π against time, from the T snapshots
(1 τ_d spacing, as the paper's 100-step spacing). Bottom: a(τ) at the same spacing, with a
bar at every reversal detected by the Schmitt trigger. Gate T5.1.
"""

import numpy as np
import matplotlib.pyplot as plt

from runfile import RunFile
from style import FIELD_CMAP, parse_args, reversal_bars, savefig

args = parse_args(__doc__, "hos_fig01.png")
T0, T1 = 0.0, 2000.0

with RunFile(args.run) as run:
    theta = run.theta
    tT, T = run.window("T", T0, T1)
    ta, a_full = run.window("a", T0, T1)
    every_a = run.stream("a")[1]
    stride = int(round(1.0 / (every_a * run.dtau)))        # 1 τ_d spacing for the plotted a
    rev = run.reversal_times()
    meta = run.meta

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 4.2), sharex=True,
                               gridspec_kw={"height_ratios": [2.2, 1], "hspace": 0.08})
vmax = np.percentile(np.abs(T), 99.5)
im = ax1.pcolormesh(tT, theta, T, cmap=FIELD_CMAP, vmin=-vmax, vmax=vmax, shading="nearest", rasterized=True)
ax1.set_ylim(np.pi, 0)
ax1.set_yticks([0, np.pi / 2, np.pi])
ax1.set_yticklabels(["0", "π/2", "π"])
ax1.set_ylabel("colatitude θ")
cb = fig.colorbar(im, ax=ax1, pad=0.01, fraction=0.03)
cb.set_label("T")

ax2.plot(ta[stride - 1::stride], a_full[stride - 1::stride], color="k", lw=0.5)
reversal_bars(ax2, rev, T0, T1, color="C3", lw=0.7)
ax2.axhline(0, color="0.6", lw=0.4)
ax2.set_ylim(-2, 2)
ax2.set_ylabel("a")
ax2.set_xlabel("τ (diffusion times)")
ax2.set_xlim(T0, T1)
fig.suptitle(f"HOS Fig. 1 — standard run (f = {meta['schedule_f']}, seed {meta['seed']}): "
             f"{np.sum((rev >= T0) & (rev <= T1))} reversals in {T1 - T0:.0f} τ_d", y=0.98)
savefig(fig, args.out)
