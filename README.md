# HoyngDynamo.jl

A fast, modern implementation of the bistable-oscillator geodynamo: the one-dimensional (colatitude) stochastic αΩ mean-field dynamo of Hoyng, Ossendrijver & Schmitt (2001), Schmitt, Ossendrijver & Hoyng (2001) and Hoyng, Schmitt & Ossendrijver (2002). Random fluctuations of the α-effect drive damped overtones whose back-reaction turns the dipole amplitude into an overdamped particle in a bistable potential; secular variation and reversals follow from one mechanism.

**Status:** eigenproblem (Tier 0), deterministic time stepping (Tier 1) and the stochastic α (Tier 2) implemented and validated; 2.0 µs per step with noise. No aggregators, run driver or output yet. Progress is tracked in the gate table below and in `dev/`.

## What is here

| Path | Contents |
|---|---|
| `src/`, `test/`, `bench/` | The Julia package `HoyngDynamo`, its tier-tagged tests and the hot-path benchmarks. |
| `viz/` | Python plotting and visualisation (a uv project: matplotlib, numpy, h5py). Reads the package's HDF5 output. Not needed to simulate. |
| `dev/` | Development notes: model summary (`DEVELOPMENT.md`), validation gates (`VALIDATION.md`), open questions (`QUESTIONS.md`), implementation decisions (`DECISIONS.md`). |
| `resources/` | Verified, agent-readable summaries of the three papers: every equation, parameter, table and result with page references. Start at `resources/README.md`. |

## Using the package

```julia
] add https://github.com/<owner>/hoyng-dynamo    # once pushed
```

Run the default test suite (validation Tiers 0–3) from the repository root:

```sh
julia --project -e 'using Pkg; Pkg.test()'
HOYNG_TIERS=0,1,2,3,4 julia --project -e 'using Pkg; Pkg.test()'   # include the slow Tier 4
julia bench/run.jl                                                 # hot-path benchmarks
```

Plotting:

```sh
cd viz && uv sync && uv run python -c "import matplotlib, h5py; print('ok')"
```

## Validation status

Gates are defined in `dev/VALIDATION.md`. A = automated, D = diagnostic, J = judged by eye against the paper's figure.

| Tier | Content | Status |
|---|---|---|
| 0 | linear operator and eigenproblem (Table I, C₀ = 57.9, E₀₀, …) | passed 2026-10-18, 57 tests (`dev/VALIDATION.md` §4) |
| 1 | deterministic dynamics, noise off (fixed point, growth rate Λ, saturation, 2μ relaxation) | passed 2026-10-18, 42 tests |
| 2 | noise construction and reproducibility | passed 2026-10-19, 362 tests (record/stream equivalence tested at build step 6) |
| 3 | closed-form Fokker–Planck theory (Table II theory columns) | not run |
| 4 | stochastic runs (⟨T_r⟩, variance, D₁ fit, sweeps, Poisson intervals) | not run |
| 5 | figures of the three papers, judged by eye | not run |
| 6 | Sint-800 comparison | deferred |

## The papers

- Hoyng, P., Ossendrijver, M. A. J. H. & Schmitt, D. (2001). The geodynamo as a bistable oscillator. *Geophys. Astrophys. Fluid Dyn.* 94, 263–314. doi:10.1080/03091920108203410
- Schmitt, D., Ossendrijver, M. A. J. H. & Hoyng, P. (2001). Magnetic field reversals and secular variation in a bistable geodynamo model. *Phys. Earth Planet. Inter.* 125, 119–124. doi:10.1016/S0031-9201(01)00237-0
- Hoyng, P., Schmitt, D. & Ossendrijver, M. A. J. H. (2002). A theoretical analysis of the observed variability of the geomagnetic dipole field. *Phys. Earth Planet. Inter.* 130, 143–157. doi:10.1016/S0031-9201(02)00004-3

## Licence

MIT, see `LICENSE`.
