"""HSO's τ_d measurement made systematic (Tier 6, gate T6.4; HSO §5.1 and Fig. 7).

Left: the rms misfit over lags 1–25 kyr between the Sint-800 autocorrelation and the
autocorrelation of artificial Sint-800 series (rectified a, running average 5 kyr, 799 points
1 kyr apart) drawn from 200 windows of each run, averaged over the windows, against the trial
τ_d; the histograms below are the best τ_d per window. Right: the data curve with its
Monte-Carlo error bars and the window-averaged model curves at the τ_d of least misfit for
both runs, and at HSO's 11 and 16 kyr. Reads results/sint800.h5 (scripts/sint800_comparison.jl).
"""

import argparse

import h5py
import numpy as np
import matplotlib.pyplot as plt

from style import RESULTS, FIGURES, savefig

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("curves", nargs="?", default=str(RESULTS / "sint800.h5"))
p.add_argument("-o", "--out", default=str(FIGURES / "sint800_tau_d.png"))
args = p.parse_args()

with h5py.File(args.curves, "r") as h:
    lag = h["autocorrelation/lags"][...]
    data = h["autocorrelation/unsigned"][...]
    sig = h["autocorrelation/sigma"][...]
    runs = {}
    for name in h["tau_d"]:
        g = h["tau_d"][name]
        runs[name] = dict(tau=g["tau_d"][...], misfit=g["mean_misfit"][...], best=g["best"][...], acor=g["acor"][...].T,   # HDF5.jl writes τ_d × lags; h5py reads it transposed
                          lags=g["lags"][...], median=float(g.attrs["median"]), q25=float(g.attrs["q25"]),
                          q75=float(g.attrs["q75"]), best_mean=float(g.attrs["best_mean"]))

fig, axes = plt.subplots(2, 2, figsize=(10, 6.5), gridspec_kw=dict(height_ratios=[2, 1]))
paper = {"reference run f = 5.4": (11, 2, "C0"), "second run kR 0.3, C_m 120, f 7": (16, 4, "C1")}
ax = axes[0, 0]
for name, r in runs.items():
    t0, dt, color = paper.get(name, (None, None, "k"))
    ax.plot(r["tau"], r["misfit"], color=color, lw=1.2, label=f"{name}: least at {r['best_mean']:g} kyr, per-window median {r['median']:g} ({r['q25']:g}–{r['q75']:g})")
    if t0 is not None:
        ax.axvspan(t0 - dt, t0 + dt, color=color, alpha=0.12, lw=0)
        ax.axvline(t0, color=color, lw=0.6, ls="--")
ax.set_xlabel("trial τ_d (kyr)")
ax.set_ylabel("rms misfit of the autocorrelation, lags 1–25 kyr")
ax.set_title("window-averaged misfit (shaded: HSO's 11 ± 2 and 16 ± 4 kyr)", fontsize=8.5)
ax.legend(fontsize=7, loc="upper right")

ax = axes[1, 0]
for name, r in runs.items():
    color = paper.get(name, (None, None, "k"))[2]
    edges = np.concatenate([r["tau"] - 0.25, [r["tau"][-1] + 0.25]])
    ax.hist(r["best"], bins=edges, color=color, alpha=0.5, label=name)
ax.set_xlabel("best τ_d of each 800-kyr window (kyr)")
ax.set_ylabel("windows")
ax.legend(fontsize=7)

ax = axes[0, 1]
ax.errorbar(lag, data, yerr=sig, fmt="o", ms=2.5, color="k", ecolor="0.5", elinewidth=0.5, capsize=1.5, label="Sint-800 (unsigned), Monte-Carlo errors")
for name, r in runs.items():
    t0, _, color = paper.get(name, (None, None, "k"))
    j = int(np.argmin(np.abs(r["tau"] - r["best_mean"])))
    ax.plot(r["lags"], r["acor"][j], color=color, lw=1.2, label=f"{name}: τ_d = {r['tau'][j]:g} kyr (least misfit)")
    if t0 is not None and np.any(np.isclose(r["tau"], t0)):
        k = int(np.argmin(np.abs(r["tau"] - t0)))
        ax.plot(r["lags"], r["acor"][k], color=color, lw=0.8, ls=":", label=f"{name}: τ_d = {t0} kyr (HSO)")
ax.axhline(0, color="0.6", lw=0.4)
ax.set_xlim(0, 50)
ax.set_ylim(-0.5, 1)
ax.set_xlabel("lag (kyr)")
ax.set_ylabel("autocorrelation")
ax.set_title("data against the window-averaged model curves (τ_av = 5 kyr)", fontsize=8.5)
ax.legend(fontsize=6.5, loc="upper right")

ax = axes[1, 1]
for name, r in runs.items():
    color = paper.get(name, (None, None, "k"))[2]
    j = int(np.argmin(np.abs(r["tau"] - r["best_mean"])))
    ax.plot(r["lags"], r["acor"][j] - data, color=color, lw=1.0, label=name)
ax.fill_between(lag, -sig, sig, color="0.85", lw=0, label="±1σ of the data")
ax.axhline(0, color="0.6", lw=0.4)
ax.set_xlim(0, 50)
ax.set_xlabel("lag (kyr)")
ax.set_ylabel("model − data")
ax.legend(fontsize=7)
fig.suptitle("τ_d from the autocorrelation overlap (HSO §5.1): 11 ± 2 kyr (reference run) and 16 ± 4 kyr (second run) in the paper", fontsize=9, y=0.995)
savefig(fig, args.out)
