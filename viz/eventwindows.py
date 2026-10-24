"""Reader and the two-event figure for the full-resolution windows of scripts/run_events.jl.

The file holds one group per window under ``windows/<name>`` with the datasets time, T, P
(M × n), a, a1, a2, parity and energy at every step, and the event's attributes (kind, start,
stop, tcross, crossings, amin, tau0, tau1, rule); the root attributes are the run's
metadata plus Λ and μ. HOS Figs. 3 and 4 (pp. 278–279) are two such windows each.
"""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import matplotlib.pyplot as plt

from style import FIELD_CMAP, RESULTS, savefig


def default_events() -> Path:
    return RESULTS / "events_f6.4_seed1.h5"


def _decode(v):
    return v.decode() if isinstance(v, bytes) else v


class EventWindows:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.h5 = h5py.File(self.path, "r")

    def close(self):
        self.h5.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    @property
    def meta(self) -> dict:
        return {k: _decode(v) for k, v in self.h5.attrs.items()}

    @property
    def theta(self) -> np.ndarray:
        return self.h5["grid/theta"][...]

    def names(self) -> list[str]:
        return list(self.h5["windows"].keys())

    def window(self, name: str):
        """The arrays (fields as M × n) and the event attributes of one window."""
        g = self.h5["windows"][name]
        w = {k: g[k][...] for k in g}
        for k in ("T", "P"):
            if w[k].ndim == 2:
                w[k] = w[k].T                 # HDF5.jl writes M × n; h5py reads it row-major
        attrs = {k: _decode(v) for k, v in g.attrs.items()}
        return w, attrs


def event_panel(ax_field, ax_amp, w, at, theta, Lambda=None, guide=False):
    """T(θ, τ) with the zero contour, and a(τ), |a₁|, 𝒫 below it, for one window."""
    t, T = w["time"], w["T"]
    vmax = np.percentile(np.abs(T), 99.5)
    im = ax_field.pcolormesh(t, theta, T, cmap=FIELD_CMAP, vmin=-vmax, vmax=vmax, shading="nearest", rasterized=True)
    ax_field.contour(t, theta, T, levels=[0.0], colors="k", linewidths=0.5)
    ax_field.set_ylim(np.pi, 0)
    ax_field.set_yticks([0, np.pi / 2, np.pi])
    ax_field.set_yticklabels(["0", "π/2", "π"])
    ax_field.set_ylabel("colatitude θ")
    ax_field.tick_params(labelbottom=False)
    cb = ax_field.figure.colorbar(im, ax=ax_field, pad=0.01, fraction=0.03)
    cb.set_label("T")
    label = (f"{at['kind']}: τ = {at['start']:.2f}–{at['stop']:.2f}, {int(at['crossings'])} zero crossing"
             f"{'s' if at['crossings'] != 1 else ''}, min |a| = {at['amin']:.3f}")
    ax_field.text(0.01, 0.97, label, transform=ax_field.transAxes, va="top", ha="left", fontsize=7.5,
                  bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5))

    ax_amp.plot(t, w["a"], color="k", lw=0.8, label="a")
    ax_amp.plot(t, np.abs(w["a1"]), color="C0", lw=0.7, label=r"$|a_1|$")
    ax_amp.plot(t, w["parity"], color="C2", lw=0.7, ls="--", label=r"$\mathcal{P}$")
    ax_amp.axhline(0, color="0.6", lw=0.4)
    for y in (0.5, -0.5):
        ax_amp.axhline(y, color="0.6", lw=0.4, ls=":")
    for x in (at["start"], at["stop"]):
        ax_amp.axvline(x, color="0.5", lw=0.5)
    if guide and Lambda is not None and at["kind"] == "reversal" and np.isfinite(at["tcross"]):
        # HOS Fig. 3 caption: after the crossing the whole field is small and grows roughly as exp(Λτ)
        tm = at["tcross"] + 0.3
        j = int(np.searchsorted(t, tm))
        if j < len(t):
            tt = t[j:]
            g = w["a"][j] * np.exp(Lambda * (tt - tt[0]))
            keep = np.abs(g) <= 1.2
            keep.any() and ax_amp.plot(tt[keep], g[keep], color="C3", ls=":", lw=1.0, label=r"$\propto e^{\Lambda\tau}$")
    ax_amp.set_ylim(-1.6, 1.6)
    ax_amp.set_xlim(t[0], t[-1])
    ax_field.set_xlim(t[0], t[-1])
    ax_amp.set_ylabel(r"a, $|a_1|$, $\mathcal{P}$")
    ax_amp.ticklabel_format(axis="x", useOffset=False)
    ax_amp.legend(ncol=4, loc="lower right", fontsize=7, handlelength=1.5)


def two_window_figure(path, names, title, out, guide=False):
    with EventWindows(path) as ew:
        theta, meta = ew.theta, ew.meta
        fig, axes = plt.subplots(4, 1, figsize=(8, 8.5), gridspec_kw={"height_ratios": [2, 1, 2, 1], "hspace": 0.14})
        for k, name in enumerate(names):
            w, at = ew.window(name)
            event_panel(axes[2 * k], axes[2 * k + 1], w, at, theta, meta.get("Lambda"), guide)
            axes[2 * k + 1].set_xlabel("τ (diffusion times)")
        fig.suptitle(title.format(**meta), y=0.94)
        return savefig(fig, out)
