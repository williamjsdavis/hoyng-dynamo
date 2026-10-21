# HoyngDynamo.jl

A fast, modern implementation of the bistable-oscillator geodynamo: the one-dimensional (colatitude) stochastic αΩ mean-field dynamo of Hoyng, Ossendrijver & Schmitt (2001), Schmitt, Ossendrijver & Hoyng (2001) and Hoyng, Schmitt & Ossendrijver (2002). Random fluctuations of the α-effect drive damped overtones whose back-reaction turns the dipole amplitude into an overdamped particle in a bistable potential; secular variation and reversals follow from one mechanism.

**Status:** eigenproblem (Tier 0), deterministic time stepping (Tier 1), the stochastic α (Tier 2), the running aggregators, the Fokker–Planck theory module (Tier 3), the run driver with record and stream modes, HDF5 output and bitwise checkpoint/resume are implemented and validated. The forcing sweep reproduces HOS Table II: the mean time between reversals at all five forcings, from 19 to 11 000 τ_d, within the paper's errors, the variance column, and the Poisson statistics of the polarity intervals (Tier 4 gates T4.1–T4.7, T4.11, T4.18), at 2.3 µs per step with all statistics running; the shape parameters of the D₁ fit reproduce where the fit is well conditioned (f ≥ 7.4) and the open points are logged in `dev/DISCREPANCIES.md`. The HSO runs reproduce HSO Table 1 and 2 (the χ² fit of the stationary distribution to the amplitude histogram gives (γ, c²) = (7.9, 0.51) against their (8, 0.50); the second parameter set gives ⟨T_r⟩ 115 ± 3 against 113 ± 5). HOS Figs. 1, 2, 5, 6, 8, 9 and HSO Figs. 3, 6, 7, 8 are rendered from `viz/`. Next: the SOH variants and ramp, then the event-level diagnostics. Progress is tracked in the gate table below and in `dev/`.

## What is here

| Path | Contents |
|---|---|
| `src/`, `test/`, `bench/` | The Julia package `HoyngDynamo`, its tier-tagged tests and the hot-path benchmarks. |
| `scripts/` | Run scripts (`run_standard.jl` writes a record-mode run to `results/`, git-ignored; `run_sweep.jl` writes our Table II to `results/table2.md`; `dr5_sensitivity.jl` the convergence and sensitivity runs). |
| `viz/` | Python plotting and visualisation (a uv project: matplotlib, numpy, h5py). Reads the package's HDF5 output. Not needed to simulate. |
| `dev/` | Development notes: model summary (`DEVELOPMENT.md`), validation gates (`VALIDATION.md`), open questions (`QUESTIONS.md`), implementation decisions (`DECISIONS.md`), discrepancies between the papers and this implementation (`DISCREPANCIES.md`), stretch goals outside the build order (`STRETCH_GOALS.md`). |
| `resources/` | Verified, agent-readable summaries of the three papers: every equation, parameter, table and result with page references. Start at `resources/README.md`. |

## Using the package

```julia
] add https://github.com/<owner>/hoyng-dynamo    # once pushed
```

Run the default test suite (validation Tiers 0–3) from the repository root:

```sh
julia --project -e 'using Pkg; Pkg.test()'
HOYNG_TIERS=0,1,2,3,4 JULIA_NUM_THREADS=4 julia --project -e 'using Pkg; Pkg.test()'   # include the slow Tier 4 (seeds run in parallel)
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
| 3 | closed-form Fokker–Planck theory (Table II theory columns) | passed 2026-10-20, 301 tests |
| 4 | stochastic runs (⟨T_r⟩, variance, D₁ fit, sweeps, Poisson intervals) | standard run and forcing sweep passed 2026-10-20 (T4.1, T4.2, T4.5, T4.6, T4.11, T4.18; `HOYNG_TIERS=4`, ~5 min on 4 threads); the D₁-fit intercept and c² miss the 10% gate at f ≤ 6.4 (`@test_broken`, DR-5); convergence T4.17 and the HSO runs T4.8, T4.9 passed (the χ² histogram fit reproduces HSO Table 1); SOH runs and event diagnostics pending |
| 5 | figures of the three papers, judged by eye | HOS Figs. 1, 2, 8, 9 judged comparable by the user 2026-10-20; HOS Figs. 5, 6 and HSO Figs. 3, 6, 7, 8 rendered, judgement pending |
| 6 | Sint-800 comparison | deferred |

## The papers

- Hoyng, P., Ossendrijver, M. A. J. H. & Schmitt, D. (2001). The geodynamo as a bistable oscillator. *Geophys. Astrophys. Fluid Dyn.* 94, 263–314. doi:10.1080/03091920108203410
- Schmitt, D., Ossendrijver, M. A. J. H. & Hoyng, P. (2001). Magnetic field reversals and secular variation in a bistable geodynamo model. *Phys. Earth Planet. Inter.* 125, 119–124. doi:10.1016/S0031-9201(01)00237-0
- Hoyng, P., Schmitt, D. & Ossendrijver, M. A. J. H. (2002). A theoretical analysis of the observed variability of the geomagnetic dipole field. *Phys. Earth Planet. Inter.* 130, 143–157. doi:10.1016/S0031-9201(02)00004-3

## Licence

MIT, see `LICENSE`.
