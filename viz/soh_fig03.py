"""SOH Fig. 3 (p. 123): the polarity record of the ramped run, f from 4 to 6.5 over 5×10⁴ τ_d.

Top: polarity bar code (black normal, white reversed) from the Schmitt-trigger reversal
times, as in the paper. Below (our addition): f(τ), and in windows of 2000 τ_d the rms
variability of |a| (√dv, HSO Eq. 12) and the number of reversals. SOH: the chrons shorten
drastically and the rms variability rises roughly linearly from 17% to 41%; the first
normal chron spans roughly the first 19 000 τ_d in their realisation. Gate T5.12.
"""

import argparse

import numpy as np
import matplotlib.pyplot as plt

from runfile import RunFile
from style import RESULTS, FIGURES, savefig

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("run", nargs="?", default=str(RESULTS / "ramp_f4.0-6.5_seed1.h5"))
p.add_argument("--window", type=float, default=2000.0, help="window for the variability and the reversal count, in τ_d")
p.add_argument("-o", "--out", default=str(FIGURES / "soh_fig03.png"))
args = p.parse_args()

with RunFile(args.run) as run:
    a, _, dt = run.stream("a")
    rev = run.reversal_times()
    meta = run.meta
f0, f1, t0, t1 = (float(meta[k]) for k in ("schedule_f0", "schedule_f1", "schedule_τ0", "schedule_τ1"))
T = len(a) * dt
W = args.window
nwin = int(T // W)
tw = (np.arange(nwin) + 0.5) * W
absa = np.abs(a[::int(round(1 / dt))])              # |a| at 1 τ_d, as dv is defined
per = int(round(W))
rms = np.array([np.sqrt(absa[i * per:(i + 1) * per].var() / absa[i * per:(i + 1) * per].mean() ** 2) for i in range(nwin)])
nrev = np.array([((rev >= i * W) & (rev < (i + 1) * W)).sum() for i in range(nwin)])
f_of = lambda t: np.clip(f0 + (f1 - f0) * (t - t0) / (t1 - t0), min(f0, f1), max(f0, f1))

fig, axes = plt.subplots(3, 1, figsize=(9, 5.2), sharex=True, gridspec_kw={"height_ratios": [0.4, 1, 1], "hspace": 0.15})
ax = axes[0]
edges = np.concatenate([[0.0], rev, [T]])
pol = 1 if a[0] > 0 else -1
for s, e in zip(edges[:-1], edges[1:]):
    if pol > 0:
        ax.axvspan(s, e, color="k")
    pol = -pol
ax.set_yticks([])
ax.set_title(f"SOH Fig. 3 — f ramped {f0:g} → {f1:g} over {T:.0f} τ_d (seed {meta['seed']}): polarity (black = +), "
             f"{len(rev)} reversals, first at τ = {rev[0]:.0f} (SOH ≈ 19 000)", fontsize=8.5)
ax = axes[1]
ax.plot(tw, 100 * rms, "o-", color="k", ms=3, lw=0.8, label=f"rms variability √dv of |a| per {W:.0f} τ_d window (%)")
ax.axhline(17, color="C3", lw=0.6, ls="--"); ax.axhline(41, color="C3", lw=0.6, ls="--")
ax.text(T * 0.99, 18, "SOH: 17% at the start", ha="right", fontsize=7, color="C3")
ax.text(T * 0.01, 42, "SOH: 41% at the end", ha="left", fontsize=7, color="C3")
ax.set_ylabel("rms variability (%)")
ax.set_ylim(0, 55)
ax2 = ax.twinx()
ax2.plot(tw, f_of(tw), color="C0", lw=0.8)
ax2.set_ylabel("f", color="C0")
ax2.set_ylim(3.5, 7)
ax.legend(loc="upper left", fontsize=7)
ax = axes[2]
ax.bar(tw, nrev, width=W * 0.9, color="0.6")
ax.set_ylabel(f"reversals per {W:.0f} τ_d")
ax.set_xlabel("τ (diffusion times)")
ax.set_xlim(0, T)
savefig(fig, args.out)
