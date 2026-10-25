"""The streamed conditional moments of the dipole amplitude against the lag, one curve per a
(scripts/km_streaming.jl, results/km_f<f>.h5), in the style of a conditional-moments figure:
left M⁽¹⁾/τ, right M⁽²⁾/τ, markers at the sampled lags, coloured by a, lines joining them.

The a values span both signs (default −1.4 … 1.4 in steps of 0.1), so the bistable drift is
visible as the fold of the family: a monostable drift orders the M⁽¹⁾ curves monotonically
in a, the bistable drift S(a) = Λ(1 − a²)a[1 + κa²] sends them down through zero at a = −1,
back up through zero at a = 0 and down again through zero at a = +1, and the two signs of a
are mirror images. Row 1 shows the curves on a linear lag axis, row 2 on a log axis where the
short-lag end (the noise correlation time τ_c = 0.05, dotted) is resolved, and row 3 shows the
same moments against a at a few lags with the reduced description S(a) and D(a) = Δ₀a² + Δ₁
dashed: the finite-lag drift steepens at the unstable point a = 0 and flattens at the wells
as the lag grows. The lag set is log-spaced from 0.01 to 2 τ_d, dense where the correlation
structure lives: τ_c, the overtone relaxation ≈ 0.5 and the drift time ≈ 0.7 τ_d.
Options: `--amax`, `--step` for the a values, `--taus` for the lags of row 3, `--variance` to
plot the conditional variance M⁽²*⁾ instead of M⁽²⁾.
"""

import argparse

import h5py
import numpy as np
import matplotlib.pyplot as plt

from style import RESULTS, FIGURES, savefig

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("curves", nargs="?", default=str(RESULTS / "km_f6.4.h5"))
p.add_argument("--amax", type=float, default=1.4, help="largest |a| plotted (the a values run from −amax to amax)")
p.add_argument("--step", type=float, default=0.1, help="spacing of the a values plotted")
p.add_argument("--taus", default="0.05,0.1,0.25,0.5,1,2", help="lags (τ_d) of the against-a panels; the nearest sampled lag is used")
p.add_argument("--variance", action="store_true", help="plot the conditional variance M⁽²*⁾ instead of M⁽²⁾")
p.add_argument("-o", "--out", default=None)
args = p.parse_args()

with h5py.File(args.curves, "r") as h:
    at = dict(h.attrs)
    f, τc, Λ, μ = (float(at[k]) for k in ("f", "tau_c", "Lambda", "mu"))
    Δ0, Δ1 = float(at["Delta0_run"]), float(at["Delta1_run"])
    g = h["km"]
    x, lags = g["x"][...], g["lags"][...]
    M1 = g["M1"][...].T                                   # HDF5.jl writes N_τ × N_x; h5py reads it transposed
    M2 = (g["M2star"] if args.variance else g["M2"])[...].T
    n = int(g.attrs["n"])
out = args.out or str(FIGURES / f"km_moments_f{f:g}{'_variance' if args.variance else ''}.png")

κ = μ / Λ - 1
S = lambda a: Λ * (1 - a**2) * a * (1 + κ * a**2)           # HOS Eq. 33 with HSO App. A's κ
D = lambda a: Δ0 * a**2 + Δ1                                 # HOS Eq. 61 with the run's D₁ fit; D = 2D⁽²⁾ = M⁽²⁾/τ
nearest = lambda xs, v: int(np.argmin(np.abs(xs - v)))

avals = np.arange(-args.amax, args.amax + 1e-9, args.step)
cmap = plt.get_cmap("RdYlBu_r")
norm = plt.Normalize(-args.amax, args.amax)
second = "M⁽²*⁾" if args.variance else "M⁽²⁾"
labels = ("M⁽¹⁾/τ  (τ_d⁻¹)", f"{second}/τ  (τ_d⁻¹)")

fig, axes = plt.subplots(3, 2, figsize=(10.5, 12))
fig.subplots_adjust(top=0.955, bottom=0.045, hspace=0.32, wspace=0.25)

# rows 1 and 2: the moments against the lag, one curve per a; linear and log lag axes
for row, log in ((0, False), (1, True)):
    for ax, M, label in zip(axes[row], (M1, M2), labels):
        for a0 in avals:
            y = M[:, nearest(x, a0)] / lags
            ax.plot(lags, y, "-", color="k", lw=0.5, zorder=1)
            ax.plot(lags, y, "o", ms=3.3, mfc=cmap(norm(a0)), mec="k", mew=0.4, zorder=2)
        ax.axvline(τc, color="0.6", lw=0.6, ls=":")
        ax.set_ylabel(label)
        if log:
            ax.set_xscale("log")
            ax.set_xlim(lags.min() / 1.3, lags.max() * 1.3)
            ax.set_xlabel("time-shift, τ (τ_d, log axis)")
        else:
            ax.set_xlim(0, lags.max() * 1.02)
            ax.set_xlabel("time-shift, τ (τ_d)")
    axes[row, 0].axhline(0, color="0.7", lw=0.4)
    axes[row, 1].set_ylim(bottom=0)
axes[0, 0].set_title("against the lag; colour: a, both signs", fontsize=9)
axes[0, 1].set_title("against the lag", fontsize=9)
axes[1, 0].set_title("the same, log lag axis (dotted: τ_c)", fontsize=9)
axes[1, 1].set_title("the same, log lag axis", fontsize=9)

# row 3: the moments against a at a few lags, with S(a) and D(a) of the reduced description
taus = [float(s) for s in args.taus.split(",")]
ks = sorted(set(nearest(lags, τ) for τ in taus))
lcmap = plt.get_cmap("viridis")
aa = np.linspace(x.min(), x.max(), 301)
for ax, M, label, theory, tlabel in zip(axes[2], (M1, M2), labels, (S, D),
                                        (f"S(a), κ = μ/Λ − 1 = {κ:.2f}", "D(a) = Δ₀a² + Δ₁ (run's D₁ fit)")):
    for i, k in enumerate(ks):
        c = lcmap(i / max(1, len(ks) - 1))
        ax.plot(x, M[k] / lags[k], "o-", ms=3, lw=0.8, color=c, mec="k", mew=0.3, label=f"τ = {lags[k]:g}")
    ax.plot(aa, theory(aa), "k--", lw=1.0, label=tlabel)
    ax.axhline(0, color="0.7", lw=0.4)
    for a0 in (-1, 0, 1):
        ax.axvline(a0, color="0.85", lw=0.5)
    ax.set_xlim(x.min(), x.max())
    ax.set_xlabel("a")
    ax.set_ylabel(label)
    ax.legend(fontsize=6.5, ncol=2, loc="upper right" if ax is axes[2, 1] else "lower right")
axes[2, 0].set_title("against a at fixed lags (dashed: the reduced description)", fontsize=9)
axes[2, 1].set_title("against a at fixed lags", fontsize=9)
ymax = np.nanmax([np.nanmax(M2[k] / lags[k]) for k in ks])
axes[2, 1].set_ylim(0, 1.6 * ymax)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
cb = fig.colorbar(sm, ax=axes.ravel().tolist(), pad=0.015, aspect=50)
cb.set_label("a (rows 1 and 2)")
fig.suptitle(f"Conditional moments of a, f = {f:g} ({int(at['seeds'])} seeds × {float(at['tau']):g} τ_d, {n:.1e} samples); "
             f"{len(lags)} lags from {lags.min():g} to {lags.max():g} τ_d", fontsize=9, y=0.985)
savefig(fig, out)
