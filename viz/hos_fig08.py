"""HOS Fig. 8 (p. 292): the measured D₁/ℱ against |a|.

One dot per τ_d (the d1 stream and a at the same steps), the bin averages of the run's
BinnedMean (histogram), the paper's fit 1.36 + 2.00 a² and this run's weighted even-quadratic
fit for |a| ≤ 1.1 (results/summary/d1_fit). Gate T5.7.
"""

import numpy as np
import matplotlib.pyplot as plt

from runfile import RunFile
from style import parse_args, savefig

args = parse_args(__doc__, "hos_fig08.png")

with RunFile(args.run) as run:
    d1, every_d1, _ = run.stream("d1")
    a, every_a, _ = run.stream("a")
    idx = (np.arange(1, len(d1) + 1) * every_d1) // every_a - 1   # the a sample at each d1 step
    aa = np.abs(a[idx])
    bins = run.result("d1")
    fit = run.result("summary")["d1_fit"]
    meta = run.meta

fig, ax = plt.subplots(figsize=(5.5, 4))
ax.plot(aa, d1, ",", color="0.55", alpha=0.6, rasterized=True, label=f"{len(d1)} samples, 1 per τ_d")
centers, mean, n = bins["centers"], bins["mean"], bins["n"]
ok = n > 0
width = centers[1] - centers[0]
ax.bar(centers[ok], mean[ok], width=width, color="none", edgecolor="C0", lw=1.0, label="bin averages")
x = np.linspace(0, 1.1, 200)
ax.plot(x, 1.36 + 2.00 * x**2, color="C3", lw=1.4, label="HOS Eq. 61: 1.36 + 2.00 a²")
ax.plot(x, fit["c0"] + fit["c2"] * x**2, color="k", lw=1.2, ls="--",
        label=f"this run: {fit['c0']:.2f} + {fit['c2']:.2f} a²")
ax.set_xlim(0, 1.5)
ax.set_ylim(0, 15)
ax.set_xlabel("|a|")
ax.set_ylabel(r"$D_1 / \mathcal{F}$")
ax.legend(loc="upper left")
ax.set_title(f"HOS Fig. 8 — standard run (seed {meta['seed']})")
savefig(fig, args.out)
