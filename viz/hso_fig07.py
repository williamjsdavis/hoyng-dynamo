"""HSO Fig. 7 (p. 152): the model autocorrelation by HSO's recipe, for their six settings.

From the recorded a(τ) (every 0.01 τ_d): choose τ_d in kyr, rectify, running-average over
τ_av kyr, resample every 1 kyr, take 799 points, subtract the mean, C_k = Σ_{i=1}^{749} x_i x_{i+k}
for k = 0…50 kyr, normalised to C_0. Panels as in the paper: (τ_d, τ_av) = (5.9, 5), (10, 5),
(14.3, 5), (10, 0), (10, 10) and (12.5, 5) kyr; the last from a different stretch of the run
(HSO used a different artificial series). Sint-800 is not in this repository (Tier 6
deferred): its features from the resource file are marked instead (zero crossing ≈ 27 kyr,
minimum ≈ −0.2 near 40 kyr, errors ≈ 0.035). The dotted curves are exp(−2Λ s) and
exp(−2μ s) (HSO Eq. 15 and the deterministic relaxation, dev/QUESTIONS.md Q5). Gate T5.13.
"""

import argparse

import numpy as np
import matplotlib.pyplot as plt

from runfile import RunFile
from style import RESULTS, FIGURES, savefig

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("run", nargs="?", default=str(RESULTS / "standard_f5.4_seed1.h5"))
p.add_argument("--start", type=float, default=1000.0, help="start of the artificial series in τ_d")
p.add_argument("-o", "--out", default=str(FIGURES / "hso_fig07.png"))
args = p.parse_args()



def artificial_series(a, dt, tau_d, tau_av, start, npts=799):
    """HSO §5.1: rectify, running-average over tau_av kyr, resample at 1 kyr, npts points, mean removed."""
    x = np.abs(a)
    w = max(int(round(tau_av / tau_d / dt)), 1)
    if w > 1:
        x = np.convolve(x, np.ones(w) / w, mode="same")
    step = (1.0 / tau_d) / dt                       # samples per kyr
    i0 = int(round(start / dt))
    idx = i0 + np.rint(np.arange(npts) * step).astype(int)
    y = x[idx]
    return y - y.mean()


def autocorr(y, kmax=50):
    n = len(y) - kmax
    c = np.array([np.sum(y[:n] * y[k:k + n]) for k in range(kmax + 1)])
    return c / c[0]


with RunFile(args.run) as run:
    a, _, dt = run.stream("a")
    meta = run.meta
    summary = run.result("summary")
LAMBDA, MU = float(summary["Λ"]), float(summary["μ"])     # this model's rates (HSO case 1: 0.527, 0.763)

settings = [(5.9, 5.0, args.start), (10.0, 5.0, args.start), (14.3, 5.0, args.start),
            (10.0, 0.0, args.start), (10.0, 10.0, args.start), (12.5, 5.0, args.start + 5000.0)]
fig, axes = plt.subplots(2, 3, figsize=(10, 5.5), sharex=True, sharey=True)
lag = np.arange(51)
for ax, (tau_d, tau_av, start), panel in zip(axes.flat, settings, "abcdef"):
    y = artificial_series(a, dt, tau_d, tau_av, start)
    ax.plot(lag, autocorr(y), color="k", lw=1.0, label="model")
    ax.plot(lag, np.exp(-2 * LAMBDA * lag / tau_d), color="C0", lw=0.8, ls=":", label="exp(−2Λs)")
    ax.plot(lag, np.exp(-2 * MU * lag / tau_d), color="C2", lw=0.8, ls=":", label="exp(−2μs)")
    ax.axhline(0, color="0.6", lw=0.4)
    ax.axvline(27, color="C3", lw=0.6, ls="--")
    ax.plot([40], [-0.2], "x", color="C3", ms=6)
    ax.set_title(f"({panel}) τ_d = {tau_d:g} kyr, τ_av = {tau_av:g} kyr, from τ = {start:.0f}", fontsize=8)
    ax.set_ylim(-0.5, 1)
    ax.set_xlim(0, 50)
for ax in axes[1]:
    ax.set_xlabel("lag (kyr)")
for ax in axes[:, 0]:
    ax.set_ylabel("autocorrelation")
axes[0, 0].legend(loc="upper right", fontsize=7)
fig.suptitle(f"HSO Fig. 7 — model curves for the reference run (seed {meta['seed']}); red: Sint-800's zero crossing ≈ 27 kyr and minimum ≈ −0.2 near 40 kyr (resource read-off, no data)",
             fontsize=8.5, y=0.98)
savefig(fig, args.out)
