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
| `hos_fig05.py` | HOS Fig. 5: a(τ), 0–2000 τ_d, at f = 5.4 and 7.4 with the same noise | the `a` stream of `standard_f5.4_seed1.h5` and `standard_f7.4_seed1.h5` (`run_standard.jl 1 1e4 5.4`, `… 1 2000 7.4`) |
| `hos_fig06.py` | HOS Fig. 6: ℜa₁ against a through a reversal at f = 5.4 | the `a` and `a1` streams, `results/trigger` of `standard_f5.4_seed1.h5` |
| `hso_fig03.py` | HSO Fig. 3: T(θ, τ) and a(τ), 0–3000 τ_d, reference run f = 5.4 | `standard_f5.4_seed1.h5` at 5×10⁴ τ_d (`run_standard.jl 1 5e4 5.4`) |
| `hso_fig06.py` | HSO Fig. 6 bottom: the \|a\| histogram with the Eq. 10 and A.1 χ² fits | `results/histogram`, `results/summary/histogram_fit` of the reference run |
| `hso_fig07.py` | HSO Fig. 7: model autocorrelation by HSO's recipe for six (τ_d, τ_av) settings | the `a` stream of the reference run |
| `hso_fig08.py` | HSO Fig. 8: polarity bar code and a(τ) of the second run kR 0.3, C_m 120, f 7 | `run_kR0.3_Cm120_f7.0_seed1.h5` (`run_standard.jl 1 1e4 7 0.3 120`) |
| `hos_fig08.py` | HOS Fig. 8: D₁/ℱ against \|a\| with the bin averages and fits | the `d1` and `a` streams, `results/d1`, `results/summary` |
| `hos_fig09.py` | HOS Fig. 9: simulated p(\|a\|) against Eq. 43 | `results/histogram`, `results/summary` |

Shared modules: `runfile.py` (reader for the HDF5 layout of `dev/DECISIONS.md` D7) and
`style.py` (rcParams, output paths, the reversal-bar helper). Every script takes the run file
as its first argument and defaults to `../results/standard_f6.4_seed1.h5`, written by

```sh
julia --project scripts/run_standard.jl          # from the repository root; ~15 s, ~130 MB
uv run python hos_fig01.py                       # then each figure
```
