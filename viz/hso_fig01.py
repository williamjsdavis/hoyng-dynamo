"""HSO Fig. 1 (p. 144): the Sint-800 VADM record (Guyodo & Valet 1999) against age, time
running right to left, with the sign flipped for the points older than the Brunhes–Matuyama
reversal at 780 ka as in the paper; the shaded band is the ±1σ error of the calibrated VADM.

Reads results/sint800.h5 written by scripts/sint800_comparison.jl (or data/sint800.tsv
directly). Gate T6.1 (dv in the title).
"""

import argparse
from pathlib import Path

import h5py
import numpy as np
import matplotlib.pyplot as plt

from style import RESULTS, FIGURES, savefig

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("curves", nargs="?", default=str(RESULTS / "sint800.h5"))
p.add_argument("-o", "--out", default=str(FIGURES / "hso_fig01.png"))
args = p.parse_args()

if Path(args.curves).exists():
    with h5py.File(args.curves, "r") as h:
        age, vadm, err, signed = (h["data"][k][...] for k in ("age", "vadm", "vadm_err", "signed"))
        dv, dv_sigma = float(h["data"].attrs["dv"]), float(h["data"].attrs["dv_sigma"])
else:
    t = np.loadtxt(Path(__file__).resolve().parent.parent / "data" / "sint800.tsv", comments="#", skiprows=3)
    age, vadm, err = t[:, 0], t[:, 5], t[:, 6]
    signed = np.where(age > 780, -vadm, vadm)
    m = vadm.mean(); dv, dv_sigma = ((vadm - m) ** 2).mean() / m**2, float("nan")

fig, ax = plt.subplots(figsize=(8, 3))
ax.fill_between(age, signed - err, signed + err, color="0.8", lw=0, label="±1σ")
ax.plot(age, signed, color="k", lw=0.8, label="Sint-800 VADM, sign flipped beyond 780 ka")
ax.axhline(0, color="0.5", lw=0.4)
ax.set_xlim(800, 0)
ax.set_ylim(-5, 10)
ax.set_xlabel("age (kyr)")
ax.set_ylabel("VADM (10²² A m²)")
ax.legend(loc="lower left", fontsize=7)
ax.set_title(f"HSO Fig. 1 — Sint-800 (Guyodo & Valet 1999), {len(age)} points, 1 per kyr; dv = {dv:.4f} ± {dv_sigma:.4f} (paper 0.071 ± 0.003)", fontsize=8.5)
savefig(fig, args.out)
