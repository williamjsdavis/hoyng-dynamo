"""HSO Fig. 3 (p. 146): T(θ, τ) and a(τ) over 0–3000 τ_d of the reference run (f = 5.4).

HOS run 2 rerun: kR 0.5, C_m 100, f 5.4. The paper's panel shows about 7 sign changes in
3000 τ_d with the longest interval ≈ 1100 τ_d, and the run gave 88 reversals in 5×10⁴ τ_d,
⟨T_r⟩ = 570 ± 61, dv = 0.070. Gate T5.9.
"""

import numpy as np
import matplotlib.pyplot as plt

from runfile import RunFile
from style import FIELD_CMAP, RESULTS, FIGURES, reversal_bars, savefig
import argparse

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("run", nargs="?", default=str(RESULTS / "standard_f5.4_seed1.h5"))
p.add_argument("-o", "--out", default=str(FIGURES / "hso_fig03.png"))
args = p.parse_args()
T0, T1 = 0.0, 3000.0

with RunFile(args.run) as run:
    tT, Tf = run.window("T", T0, T1)
    ta, a = run.window("a", T0, T1)
    theta = run.theta
    rev = run.reversal_times()
    s = run.result("summary")
    meta = run.meta

nrev = int(((rev >= T0) & (rev <= T1)).sum())
vmax = np.abs(Tf).max()
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 5), sharex=True, gridspec_kw={"height_ratios": [1.3, 1], "hspace": 0.12})
ax1.pcolormesh(tT, theta, Tf, cmap=FIELD_CMAP, vmin=-vmax, vmax=vmax, shading="nearest", rasterized=True)
ax1.set_yticks([0, np.pi / 2, np.pi]); ax1.set_yticklabels(["0", "π/2", "π"])
ax1.set_ylabel("colatitude θ")
ax1.set_ylim(np.pi, 0)
ax1.set_title(f"HSO Fig. 3 — reference run (kR 0.5, C_m 100, f = {float(meta['schedule_f']):g}, seed {meta['seed']}): "
              f"{nrev} reversals in {T1 - T0:.0f} τ_d; whole run ⟨T_r⟩ = {float(s['mean_reversal_time']):.0f} ± "
              f"{float(s['mean_reversal_time_error']):.0f}, dv = {float(s['dv']):.3f}", fontsize=8.5)
ax2.plot(ta, a, color="k", lw=0.35)
reversal_bars(ax2, rev, T0, T1, color="C3", lw=0.7)
ax2.axhline(0, color="0.6", lw=0.4)
ax2.set_ylim(-2, 2)
ax2.set_ylabel("a")
ax2.set_xlabel("τ (diffusion times)")
ax2.set_xlim(T0, T1)
savefig(fig, args.out)
