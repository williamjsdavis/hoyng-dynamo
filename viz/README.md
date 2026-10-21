# viz: plotting and visualisation for HoyngDynamo

Python side of the project (dev/DECISIONS.md D13). Reads the HDF5 output of the Julia
package and makes the figures; it is never needed to run a simulation.

Managed by [uv](https://docs.astral.sh/uv/); Python 3.12, matplotlib, numpy, h5py.

```sh
cd viz
uv sync                      # create .venv from uv.lock
uv run python <script>.py    # run a figure script in the environment
```

Layout: one script per paper figure, named after it, reading a run file and writing a PNG
to `../figures/` (git-ignored):

| Script | Figure | Needs |
|---|---|---|
| `hos_fig01.py` | HOS Fig. 1: T(θ, τ) colour plot and a(τ), 0–2000 τ_d | the `T` snapshots and the `a` stream |
| `hos_fig02.py` | HOS Fig. 2: \|a₁\|, \|a₂\|, parity, energy over τ = 500–700 | the `a1`, `a2`, `parity`, `energy` streams |
| `hos_fig08.py` | HOS Fig. 8: D₁/ℱ against \|a\| with the bin averages and fits | the `d1` and `a` streams, `results/d1`, `results/summary` |
| `hos_fig09.py` | HOS Fig. 9: simulated p(\|a\|) against Eq. 43 | `results/histogram`, `results/summary` |

Shared modules: `runfile.py` (reader for the HDF5 layout of `dev/DECISIONS.md` D7) and
`style.py` (rcParams, output paths, the reversal-bar helper). Every script takes the run file
as its first argument and defaults to `../results/standard_f6.4_seed1.h5`, written by

```sh
julia --project scripts/run_standard.jl          # from the repository root; ~15 s, ~130 MB
uv run python hos_fig01.py                       # then each figure
```
