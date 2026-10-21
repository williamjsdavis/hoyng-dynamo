"""HSO Fig. 8 (p. 153): polarity bar code to 10⁴ τ_d and a(τ) to 3000 τ_d of the second run
(kR 0.3, C_m 120, f 7).

The paper's run: ⟨T_r⟩ = 113 ± 5 τ_d, dv = 0.081; of order 20–25 sign changes in the
3000 τ_d shown. Gate T5.11. Polarity from the Schmitt-trigger reversal times.
"""

import argparse

import numpy as np
import matplotlib.pyplot as plt

from runfile import RunFile
from style import RESULTS, FIGURES, reversal_bars, savefig

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("run", nargs="?", default=str(RESULTS / "run_kR0.3_Cm120_f7.0_seed1.h5"))
p.add_argument("-o", "--out", default=str(FIGURES / "hso_fig08.png"))
args = p.parse_args()
TB, TA = 10_000.0, 3000.0

with RunFile(args.run) as run:
    ta, a = run.window("a", 0.0, TA)
    rev = run.reversal_times()
    a_all, _, dt = run.stream("a")
    s = run.result("summary")
    meta = run.meta

# polarity: +1 from the start (a0 = +1), flipping at each reversal
pol0 = 1 if a_all[0] > 0 else -1
edges = np.concatenate([[0.0], rev[rev <= TB], [min(TB, len(a_all) * dt)]])
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 4), gridspec_kw={"height_ratios": [0.35, 1], "hspace": 0.45})
pol = pol0
for t0, t1 in zip(edges[:-1], edges[1:]):
    if pol > 0:
        ax1.axvspan(t0, t1, color="k")
    pol = -pol
ax1.set_xlim(0, TB)
ax1.set_yticks([])
ax1.set_xlabel("τ (diffusion times)", fontsize=8)
ax1.set_title(f"HSO Fig. 8 — kR {float(meta['kR']):g}, C_m {float(meta['Cm']):g}, f = {float(meta['schedule_f']):g}, seed {meta['seed']}: "
              f"polarity (black = +) to {TB:.0f} τ_d; {int((rev <= TB).sum())} reversals; whole run ⟨T_r⟩ = {float(s['mean_reversal_time']):.1f} ± "
              f"{float(s['mean_reversal_time_error']):.1f} (HSO 113 ± 5), dv = {float(s['dv']):.3f} (0.081)", fontsize=8)
ax2.plot(ta, a, color="k", lw=0.35)
reversal_bars(ax2, rev, 0.0, TA, color="C3", lw=0.6)
ax2.axhline(0, color="0.6", lw=0.4)
ax2.set_xlim(0, TA)
ax2.set_ylim(-2, 2)
ax2.set_ylabel("a")
ax2.set_xlabel("τ (diffusion times)")
ax2.text(0.01, 0.95, f"{int((rev <= TA).sum())} reversals in {TA:.0f} τ_d", transform=ax2.transAxes, va="top", fontsize=8,
         bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none", "pad": 1.5})
savefig(fig, args.out)
