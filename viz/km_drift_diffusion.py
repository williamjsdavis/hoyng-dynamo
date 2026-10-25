"""Drift and diffusion of the dipole amplitude by online kernel-based regression (Davis 2023),
against the Fokker–Planck reduction of HOS §4 (scripts/km_streaming.jl, results/km_f<f>.h5).

Top left: D⁽¹⁾(a) from the regression over two lag windows with the coherent-sum error bars,
and the drift S(a) = Λ(1 − a²)a[1 + κa²] with the run's (c², γ) (κ = μ/Λ − 1, HSO App. A) and
with κ = 0 (HOS Eq. 33). Top right: 2D⁽²⁾(a) from the conditional variance over the same
windows, D(a) = Δ₀a² + Δ₁ from the run's D₁ fit (HOS Eqs. 61–63), and the noise diffusion of
a evaluated on the run's fields: exact for the binned noise (D_b) and in HOS's δ-correlated
limit (D_c, Eqs. 36–37). Bottom: the lag dependence of M⁽¹⁾/τ and M⁽²*⁾/τ at a = 0, 0.5, 1.0
over the lag set (log axis); the dotted line marks τ_c. viz/km_moments.py shows the moments
of every a in the style of a conditional-moments figure.
"""

import argparse

import h5py
import numpy as np
import matplotlib.pyplot as plt

from style import RESULTS, FIGURES, savefig

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("curves", nargs="?", default=str(RESULTS / "km_f6.4.h5"))
p.add_argument("-o", "--out", default=None)
args = p.parse_args()

with h5py.File(args.curves, "r") as h:
    at = dict(h.attrs)
    f, Λ, μ = float(at["f"]), float(at["Lambda"]), float(at["mu"])
    γ, c2, Δ0, Δ1, τc = (float(at[k]) for k in ("gamma_run", "c2_run", "Delta0_run", "Delta1_run", "tau_c"))
    reg = {k: {n: h["regression"][k][n][...] for n in ("x", "D1", "D2", "D1_err", "D2_err")} for k in h["regression"]}
    # HDF5.jl writes the N_τ × N_x matrices column-major; h5py reads them transposed
    km = {n: (h["km"][n][...].T if h["km"][n].ndim == 2 else h["km"][n][...]) for n in ("x", "lags", "M1", "M2star", "W", "neff")}
    db = {n: h["d_binned"][n][...] for n in ("centers", "n", "mean")}
    dc = {n: h["d_continuum"][n][...] for n in ("centers", "n", "mean")}
out = args.out or str(FIGURES / f"km_drift_diffusion_f{f:g}.png")

κ = μ / Λ - 1
a = np.linspace(-1.5, 1.5, 301)
S = lambda a, kap: Λ * (1 - a**2) * a * (1 + kap * a**2)
D = Δ0 * a**2 + Δ1

fig, axes = plt.subplots(2, 2, figsize=(10, 7.5))
fig.subplots_adjust(hspace=0.35)
ax = axes[0, 0]
for key, color in (("0.1–0.5", "C0"), ("0.25–1.0", "C3")):
    r = reg[key]
    ax.errorbar(r["x"], r["D1"], yerr=r["D1_err"], fmt="o", ms=2.5, color=color, ecolor=color, elinewidth=0.6, capsize=1.5, label=f"D⁽¹⁾, lags {key} τ_d")
ax.plot(a, S(a, κ), color="k", lw=1.2, label=f"S(a), κ = μ/Λ − 1 = {κ:.2f}")
ax.plot(a, S(a, 0), color="k", lw=0.8, ls="--", label="S(a), κ = 0 (HOS Eq. 33)")
ax.axhline(0, color="0.6", lw=0.4)
ax.set_xlabel("a"); ax.set_ylabel("drift (τ_d⁻¹)")
ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.5, 1.5)
ax.legend(fontsize=7, loc="upper right")
ax.set_title("drift", fontsize=9)

ax = axes[0, 1]
for key, color in (("0.1–0.5", "C0"), ("0.25–1.0", "C3")):
    r = reg[key]
    ax.errorbar(r["x"], 2 * r["D2"], yerr=2 * r["D2_err"], fmt="o", ms=2.5, color=color, ecolor=color, elinewidth=0.6, capsize=1.5, label=f"2D⁽²⁾ (conditional variance), lags {key} τ_d")
ax.plot(a, D, color="k", lw=1.2, label=f"Δ₀a² + Δ₁ from the run's D₁ fit: (c², γ) = ({c2:.3f}, {γ:.2f})")
sel = db["n"] > 0
ax.plot(db["centers"][sel], db["mean"][sel], "s", ms=3, color="C2", label="D_b: exact binned-noise diffusion on the run's fields")
ax.plot(dc["centers"][sel], dc["mean"][sel], "^", ms=3, color="C4", label="D_c: HOS Eqs. 36–37 (δ-correlated limit) on the same fields")
ax.set_xlabel("a (|a| for the binned points)"); ax.set_ylabel("diffusion (τ_d⁻¹)")
ax.set_xlim(-1.5, 1.5); ax.set_ylim(0, 1.1 * max(np.nanmax(D), np.nanmax(dc["mean"][sel])))
ax.legend(fontsize=6.5, loc="upper left")
ax.set_title("diffusion", fontsize=9)

for ax, (name, ylabel) in zip(axes[1], (("M1", "M⁽¹⁾/τ"), ("M2star", "M⁽²*⁾/τ"))):
    for a0, color in ((0.0, "C0"), (0.5, "C2"), (1.0, "C3")):
        j = int(np.argmin(np.abs(km["x"] - a0)))
        ax.plot(km["lags"], km[name][:, j] / km["lags"], "o-", ms=3, color=color, lw=0.6, label=f"a = {a0:g}")
        if name == "M1":
            ax.axhline(S(a0, κ), color=color, lw=0.6, ls="--")
        else:
            ax.axhline(Δ0 * a0**2 + Δ1, color=color, lw=0.6, ls="--")
    ax.axvline(τc, color="0.5", lw=0.6, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel("lag τ (τ_d)"); ax.set_ylabel(ylabel)
    ax.legend(fontsize=7)
    ax.set_title(f"{ylabel} against the lag (dashed: S(a) or D(a) of the reduced description)", fontsize=8)
fig.suptitle(f"Drift and diffusion of a by online kernel-based regression, f = {f:g} ({int(at['seeds'])} seeds × {float(at['tau']):g} τ_d); Λ = {Λ:.3f}, 2μ = {2 * μ:.3f}", fontsize=9, y=0.995)
savefig(fig, out)
