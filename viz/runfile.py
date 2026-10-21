"""Reader for HoyngDynamo HDF5 run files (dev/DECISIONS.md D7).

Layout: root attributes hold the run specification and provenance; ``streams/<name>`` are
the recorded series (attributes ``every`` in steps and ``dt`` in diffusion times; sample i,
counted from 1, was taken after step i·every, i.e. at τ = i·dt); ``results/<name>`` hold the
aggregator values; ``results/summary`` the headline statistics; ``grid/theta`` the nodes.
"""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np


def _read_group(g):
    out = {}
    for k, v in g.items():
        if isinstance(v, h5py.Group):
            out[k] = _read_group(v)
        else:
            x = v[()]
            if isinstance(x, bytes):
                x = x.decode()
            elif isinstance(x, np.ndarray) and x.dtype.kind in "OS":
                x = np.array([s.decode() if isinstance(s, bytes) else s for s in x])
            out[k] = x
    return out


class RunFile:
    """A read-only view of one run file."""

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
        return {k: (v.decode() if isinstance(v, bytes) else v) for k, v in self.h5.attrs.items()}

    @property
    def dtau(self) -> float:
        return float(self.h5.attrs["dtau"])

    @property
    def theta(self) -> np.ndarray:
        return self.h5["grid/theta"][...]

    def streams(self) -> list[str]:
        return list(self.h5["streams"].keys())

    def stream(self, name: str):
        """The recorded values, the cadence in steps and in diffusion times."""
        d = self.h5["streams"][name]
        v = d[...]
        if v.ndim == 2:                       # HDF5.jl writes M × nsamples; h5py reads it row-major
            v = v.T
        return v, int(d.attrs["every"]), float(d.attrs["dt"])

    def times(self, name: str) -> np.ndarray:
        """τ of each sample of a stream (sample i was taken after step i·every)."""
        d = self.h5["streams"][name]
        n = d.shape[-1]
        return np.arange(1, n + 1) * float(d.attrs["dt"])

    def window(self, name: str, t0: float, t1: float):
        """A stream restricted to t0 ≤ τ ≤ t1: (τ, values) with the sample axis last."""
        v, _, dt = self.stream(name)
        n = v.shape[-1]
        i0 = max(int(np.ceil(t0 / dt)) - 1, 0)
        i1 = min(int(np.floor(t1 / dt)), n)
        t = np.arange(i0 + 1, i1 + 1) * dt
        return t, v[..., i0:i1]

    def result(self, name: str) -> dict:
        return _read_group(self.h5["results"][name])

    def results(self) -> dict:
        return _read_group(self.h5["results"])

    def reversal_times(self) -> np.ndarray:
        return np.atleast_1d(self.h5["results/trigger/reversal_times"][...])

    def events(self) -> dict:
        return self.result("trigger")["events"]
