"""Shared plotting style and output helper for the paper figures (dev/DECISIONS.md D13)."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

FIGURES = Path(__file__).resolve().parent.parent / "figures"
RESULTS = Path(__file__).resolve().parent.parent / "results"

plt.rcParams.update({
    "figure.dpi": 120,
    "savefig.dpi": 150,
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "lines.linewidth": 0.8,
    "axes.linewidth": 0.6,
    "legend.frameon": False,
})

# T(θ, τ) colour plots: diverging, white at zero, as in HOS Plate I.
FIELD_CMAP = "RdBu_r"


def default_run() -> Path:
    """The standard run file written by scripts/run_standard.jl, if present."""
    return RESULTS / "standard_f6.4_seed1.h5"


def parse_args(description: str, default_name: str) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=description)
    p.add_argument("run", nargs="?", default=str(default_run()), help="HDF5 run file (default: the standard run)")
    p.add_argument("-o", "--out", default=str(FIGURES / default_name), help="output image path")
    return p.parse_args()


def savefig(fig, out: str | Path) -> Path:
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")
    return out


def reversal_bars(ax, times, t0, t1, **kw):
    """Vertical bars at the reversal times inside [t0, t1], as in HOS Figs. 1, 2, 5."""
    for t in times:
        if t0 <= t <= t1:
            ax.axvline(t, color=kw.pop("color", "k"), lw=kw.pop("lw", 0.6), alpha=kw.pop("alpha", 0.8), **kw)
