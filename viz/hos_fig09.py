"""HOS Fig. 9 (p. 293): the simulated amplitude distribution p(|a|) against Eq. 43.

The histogram of a at 1 τ_d spacing (results/histogram), folded onto |a| and normalised to
unit area, against the stationary distribution of Eq. 43 with the paper's shape parameters
(γ, c²) = (2.96, 0.246) and with this run's measured ones (results/summary). The paper's own
discrepancy (histogram peak ≈ 1.35 at |a| ≈ 0.85 against a curve peak ≈ 1.1 at 0.8) should
reappear. Gate T5.8.
"""

import numpy as np
import matplotlib.pyplot as plt

from runfile import RunFile
from style import parse_args, savefig


def p43(a, gamma, c2):
    """HOS Eq. 43 on |a| ≥ 0, unit area on [0, ∞) (numerically)."""
    x = np.linspace(0, 4, 40001)
    psi = lambda u: -gamma * u**2 + (gamma * c2 + gamma - 1) * np.log1p(u**2 / c2)
    Z = np.trapezoid(np.exp(psi(x)), x)
    return np.exp(psi(np.asarray(a))) / Z


args = parse_args(__doc__, "hos_fig09.png")

with RunFile(args.run) as run:
    h = run.result("histogram")
    s = run.result("summary")
    meta = run.meta

edges, counts = h["edges"], h["counts"]
REBIN = 2                                   # the run file's bins are 0.005; draw 0.01 as the paper does
edges = edges[::REBIN]
counts = counts.reshape(-1, REBIN).sum(axis=1)
centers = 0.5 * (edges[:-1] + edges[1:])
width = edges[1] - edges[0]
# fold onto |a|: bins are symmetric about zero
pos = centers > 0
folded_c = centers[pos]
folded = counts[pos] + counts[::-1][pos]
density = folded / (folded.sum() * width)

fig, ax = plt.subplots(figsize=(5.5, 4))
ax.bar(folded_c, density, width=width, color="0.8", edgecolor="0.5", lw=0.3, label=f"simulation, {int(h['n'])} samples at 1 τ_d")
x = np.linspace(0, 1.5, 600)
ax.plot(x, p43(x, 2.96, 0.246), color="C3", lw=1.4, label="Eq. 43, (γ, c²) = (2.96, 0.246) [paper]")
ax.plot(x, p43(x, float(s["γ"]), float(s["c2"])), color="k", lw=1.0, ls="--",
        label=f"Eq. 43, (γ, c²) = ({float(s['γ']):.2f}, {float(s['c2']):.3f}) [this run]")
ipk = np.argmax(density)
ax.annotate(f"histogram peak {density[ipk]:.2f} at |a| = {folded_c[ipk]:.2f}", (folded_c[ipk], density[ipk]),
            xytext=(1.0, 1.45), fontsize=8, arrowprops={"arrowstyle": "-", "lw": 0.5})
ax.set_xlim(0, 1.5)
ax.set_ylim(0, 1.6)
ax.set_xlabel("|a|")
ax.set_ylabel("p(|a|)")
ax.legend(loc="upper left")
ax.set_title(f"HOS Fig. 9 — standard run (seed {meta['seed']})")
savefig(fig, args.out)
