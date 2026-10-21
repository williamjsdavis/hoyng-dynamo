"""HSO Fig. 6, bottom panels (p. 150): the |a| histogram of the reference run with Eq. 10
and Eq. A.1 fitted (χ², ~100 equal-population bins, HSO §4.3).

Left: Eq. 10 (HOS Eq. 43, κ = 0) at this run's fit (results/summary/histogram_fit/eq10)
and at HSO's case 2 (γ, c²) = (8, 0.50). Right: Eq. A.1 (κ = μ/Λ − 1) at this run's fit
and at HSO's case 4 (4.8, 0.22), which was fitted to Sint-800 and the run together.
Counts per unit a as in the paper (peak ≈ 9×10⁴ at a ≈ 0.9 for 5×10⁴ samples). Gate T5.10.
"""

import argparse

import numpy as np
import matplotlib.pyplot as plt

from runfile import RunFile
from style import RESULTS, FIGURES, savefig

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("run", nargs="?", default=str(RESULTS / "standard_f5.4_seed1.h5"))
p.add_argument("-o", "--out", default=str(FIGURES / "hso_fig06.png"))
args = p.parse_args()


def psi(a, gamma, c2, kappa):
    """ln p(a) up to a constant: HOS Eq. 43 for κ = 0, HSO Eq. A.1 for κ = μ/Λ − 1."""
    return gamma * ((kappa * (1 + c2) - 1) * a**2 - 0.5 * kappa * a**4) + (gamma * (1 + c2) * (1 - kappa * c2) - 1) * np.log1p(a**2 / c2)


def density(a, gamma, c2, kappa):
    x = np.linspace(0, 4, 40001)
    ps = psi(x, gamma, c2, kappa); m = ps.max()
    Z = np.trapezoid(np.exp(ps - m), x)
    return np.exp(psi(a, gamma, c2, kappa) - m) / Z


with RunFile(args.run) as run:
    h = run.result("histogram")
    s = run.result("summary")
    meta = run.meta

edges, counts = h["edges"], h["counts"]
centers = 0.5 * (edges[:-1] + edges[1:]); width = edges[1] - edges[0]
pos = centers > 0
fc, fn = centers[pos], counts[pos] + counts[::-1][pos]
n = fn.sum()
fit10, fitA1 = s["histogram_fit"]["eq10"], s["histogram_fit"]["a1"]
mu_over_Lambda = float(s["μ"]) / float(s["Λ"])            # κ = μ/Λ − 1 of this model (HSO case 1: 1.45)

fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), sharey=True)
x = np.linspace(0.005, 1.5, 600)
for ax, (label, fit, kappa, paper, plabel) in zip(axes, [
        ("Eq. 10 (κ = 0)", fit10, 0.0, (8.0, 0.50), "HSO case 2: (8, 0.50), fitted to this run"),
        ("Eq. A.1 (κ = μ/Λ − 1)", fitA1, mu_over_Lambda - 1, (4.8, 0.22), "HSO case 4: (4.8, 0.22), fitted to Sint-800 and the run")]):
    ax.bar(fc, fn / width, width=width, color="0.82", edgecolor="none", label=f"{n} samples at 1 τ_d, bins of {width:g}")
    ax.plot(x, n * density(x, fit["γ"], fit["c2"], kappa), color="k", lw=1.3,
            label=f"this run: (γ, c²) = ({fit['γ']:.2f}, {fit['c2']:.3f}), {fit['sigma']:.0f}σ over {int(fit['nbins'])} bins")
    ax.plot(x, n * density(x, *paper, kappa), color="C3", lw=1.1, ls="--", label=plabel)
    ax.set_xlim(0, 1.5)
    ax.set_xlabel("|a|")
    ax.set_title(label, fontsize=9)
    ax.legend(loc="upper left", fontsize=7)
axes[0].set_ylabel("counts per unit a")
fig.suptitle(f"HSO Fig. 6 (bottom panels) — reference run f = {float(meta['schedule_f']):g}, seed {meta['seed']}", y=1.0, fontsize=9)
savefig(fig, args.out)
