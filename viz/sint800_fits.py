"""The χ² fits of p(a) to the Sint-800 VADM histogram (Tier 6, gates T6.2 and T6.5).

Top row, the upper panels of HSO Fig. 6 (p. 150): the 50-equal-population-bin histogram
(counts per 10²² A m², 1σ Poisson error bars) with (a) Eq. 10 at our optimum and at HSO's
case 1 (6.3, 0.32, 5.98), and (b) Eq. A.1 at our optimum and at HSO's case 4 (4.8, 0.22).
Bottom row: (c) HOS Fig. 10 (p. 295), 35 equal-width bins, with HOS's (5.15, 0.25) and our
optimum; (d) SOH Fig. 4 (p. 123), 30 equal-width bins, with SOH's (5.2, 0.26) and our optimum.
The curves are N p(A a_m/A_m) a_m/A_m with the N of each fit. Reads results/sint800.h5
(scripts/sint800_comparison.jl).
"""

import argparse

import h5py
import numpy as np
import matplotlib.pyplot as plt

from style import RESULTS, FIGURES, savefig

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("curves", nargs="?", default=str(RESULTS / "sint800.h5"))
p.add_argument("-o", "--out", default=str(FIGURES / "sint800_fits.png"))
args = p.parse_args()


def psi(a, gamma, c2, kappa):
    """ln p(a) up to a constant: HOS Eq. 43 for κ = 0, HSO Eq. A.1 for κ = μ/Λ − 1."""
    return gamma * ((kappa * (1 + c2) - 1) * a**2 - 0.5 * kappa * a**4) + (gamma * (1 + c2) * (1 - kappa * c2) - 1) * np.log1p(a**2 / c2)


def curve(A, fit):
    """Expected counts per unit VADM of a fit record (attributes gamma, c2, kappa, Am, am, N)."""
    g, c2, k, Am, am, N = (float(fit.attrs[n]) for n in ("gamma", "c2", "kappa", "Am", "am", "N"))
    a = A * am / Am
    return N * np.exp(psi(a, g, c2, k) - psi(am, g, c2, k)) * am / Am


def text(v):
    return v.decode() if isinstance(v, bytes) else str(v)


with h5py.File(args.curves, "r") as h:
    fits = {text(h["fits"][k].attrs["label"]) + " | " + text(h["fits"][k].attrs["bins"]): h["fits"][k] for k in h["fits"]}
    vadm = h["data/vadm"][...]

    def pick(label_start, bins):
        for key, f in fits.items():
            lab, b = key.split(" | ")
            if lab.startswith(label_start) and b == bins:
                return f
        raise KeyError(f"{label_start} / {bins}")

    def hist(fit, ax, color="0.82"):
        lo, hi, n = fit["lo"][...], fit["hi"][...], fit["n"][...]
        w = hi - lo
        ax.bar(lo, n / w, width=w, align="edge", color=color, edgecolor="white", lw=0.3)
        ax.errorbar(0.5 * (lo + hi), n / w, yerr=np.sqrt(n) / w, fmt="none", ecolor="0.4", elinewidth=0.5, capsize=1)

    A = np.linspace(0.05, 10, 500)
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7), sharex=True)
    panels = [
        (axes[0, 0], "(a) HSO Fig. 6a top: Eq. 10, 50 equal-population bins",
         [("Eq. 10, (γ, c², A_m) free (case 1)", "50 equal-population", "k", "-", "this fit"),
          ("Eq. 10 at HSO's (6.3, 0.32), A_m free", "50 equal-population", "C3", "--", "HSO case 1 (6.3, 0.32)")]),
        (axes[0, 1], "(b) HSO Fig. 6b top: Eq. A.1 (κ = μ/Λ − 1), 50 equal-population bins",
         [("Eq. A.1 (μ/Λ 1.45), free", "50 equal-population", "k", "-", "this fit"),
          ("Eq. A.1 at (4.8, 0.22) (case 4)", "50 equal-population", "C3", "--", "HSO case 4 (4.8, 0.22)")]),
        (axes[1, 0], "(c) HOS Fig. 10: Eq. 43, 35 equal-width bins",
         [("Eq. 43 free (HOS Fig. 10)", "35 equal-width from 0", "k", "-", "this fit"),
          ("Eq. 43 at HOS's (5.15, 0.25), A_m free", "35 equal-width from 0", "C3", "--", "HOS (5.15, 0.25)")]),
        (axes[1, 1], "(d) SOH Fig. 4: Eq. 4, 30 equal-width bins",
         [("Eq. 4 free (SOH Fig. 4)", "30 equal-width from 0", "k", "-", "this fit"),
          ("Eq. 4 at SOH's (5.2, 0.26), A_m free", "30 equal-width from 0", "C3", "--", "SOH (5.2, 0.26)")]),
    ]
    for ax, title, entries in panels:
        first = pick(entries[0][0], entries[0][1])
        hist(first, ax)
        for label_start, bins, color, ls, name in entries:
            f = pick(label_start, bins)
            at = f.attrs
            ax.plot(A, curve(A, f), color=color, ls=ls, lw=1.2,
                    label=f"{name}: (γ, c²) = ({float(at['gamma']):.2f}, {float(at['c2']):.3f}), A_m = {float(at['Am']):.2f}, {float(at['sigma']):.1f}σ")
        ax.set_title(title, fontsize=8.5)
        ax.legend(loc="upper left", fontsize=7)
        ax.set_xlim(0, 10)
    for ax in axes[1]:
        ax.set_xlabel("VADM (10²² A m²)")
    for ax in axes[:, 0]:
        ax.set_ylabel("counts per 10²² A m²")
    fig.suptitle(f"Sint-800 VADM histogram ({len(vadm)} points) and the fitted p(a) (Tier 6, T6.2 and T6.5)", fontsize=9, y=0.995)
    savefig(fig, args.out)
