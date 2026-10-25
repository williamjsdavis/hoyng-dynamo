# HoyngDynamo.jl

A fast, modern implementation of the bistable-oscillator geodynamo: the one-dimensional (colatitude) stochastic αΩ mean-field dynamo of Hoyng, Ossendrijver & Schmitt (2001), Schmitt, Ossendrijver & Hoyng (2001) and Hoyng, Schmitt & Ossendrijver (2002). Random fluctuations of the α-effect drive damped overtones whose back-reaction turns the dipole amplitude into an overdamped particle in a bistable potential; secular variation and reversals follow from one mechanism.

**Status:** eigenproblem (Tier 0), deterministic time stepping (Tier 1), the stochastic α (Tier 2), the running aggregators, the Fokker–Planck theory module (Tier 3), the run driver with record and stream modes, HDF5 output and bitwise checkpoint/resume are implemented and validated. The forcing sweep reproduces HOS Table II: the mean time between reversals at all five forcings, from 19 to 11 000 τ_d, within the paper's errors, the variance column, and the Poisson statistics of the polarity intervals (Tier 4 gates T4.1–T4.7, T4.11, T4.18), at 2.3 µs per step with all statistics running; the shape parameters of the D₁ fit reproduce where the fit is well conditioned (f ≥ 7.4) and the open points are logged in `dev/DISCREPANCIES.md`. The HSO runs reproduce HSO Table 1 and 2 (the χ² fit of the stationary distribution to the amplitude histogram gives (γ, c²) = (7.9, 0.51) against their (8, 0.50); the second parameter set gives ⟨T_r⟩ 115 ± 3 against 113 ± 5). SOH's ramped run reproduces the rise of the rms variability from 17% to 41% and the drastic shortening of the chrons. The event-level diagnostics complete the planned validation: the reversal duration reproduces HOS's ~3 τ_d and HSO's 0.5–2 τ_d under their respective definitions, and the parity between events sits at −0.82 where the paper says "≃ −1" (`dev/DISCREPANCIES.md` DR-7). HOS Figs. 1–6, 8, 9, HSO Figs. 3, 6, 7, 8 and SOH Fig. 3 are rendered from `viz/` and judged comparable with the originals. The first stretch goal tested the time integrator as the cause of the D₁-intercept gap (DR-5): backward-Euler and Crank–Nicolson steppers at the papers' Δτ show that a second-order scheme of another family reproduces our results and that no plausible first-order scheme reaches the paper's intercept while keeping its mean reversal time, so the integrator is excluded and the gap stays open. The second stretch goal implemented HOS's three alternative quenching laws behind a switch in the run specification: quenching by the total toroidal energy and local quenching run away as written, because an overtone burst reverses the sign of α and the α-reversed dynamo is supercritical, and once floored at α = 0 they reverse two to seven times more often; α/(1 + q) keeps the reversal rate but raises the variability by 40–50%. The paper's "no material change" is therefore recorded as not reproduced (`dev/DISCREPANCIES.md` DR-8). The Sint-800 comparison (Tier 6) is done with the record itself, found in the supplementary information of Guyodo & Valet (1999) and kept in `data/`: the dipole variation dv = 0.0706 ± 0.0027 reproduces HSO's 0.071 ± 0.003, the χ² fits of the stationary distribution to the VADM histogram land in HSO's valley with the same peak VADM, the autocorrelation has HSO's zero crossing and minimum, and a systematic version of their autocorrelation overlap gives τ_d = 12 kyr for the reference run and 20 kyr for the second against their 11 ± 2 and 16 ± 4 kyr; only HOS's 35-bin fit "at 1σ" is not reproduced (3.2σ, `dev/DISCREPANCIES.md` DR-9). A streamed estimator of drift and diffusion functions (the online kernel-based regression of Davis 2023) is available for research beyond the papers: applied to the dipole amplitude it recovers the growth rate Λ and the noise diffusion on the instantaneous field, finds the well's centre at a ≈ 0.9 rather than 1, and shows that a(τ) is not a Markov diffusion on any lag (`dev/QUESTIONS.md` Q10). The planned validation, both stretch goals and Tier 6 are complete. Progress is tracked in the gate table below and in `dev/`.

## What is here

| Path | Contents |
|---|---|
| `src/`, `test/`, `bench/` | The Julia package `HoyngDynamo`, its tier-tagged tests and the hot-path benchmarks. |
| `scripts/` | Run scripts (`run_standard.jl` writes a record-mode run to `results/`, git-ignored; `run_ramp.jl` the ramped run; `sg1_integrators.jl` the integrator-sensitivity ensembles of stretch goal SG-1; `sg2_quenching.jl` the alternative-quenching ensembles of SG-2; `run_sweep.jl` writes our Table II to `results/table2.md`; `dr5_sensitivity.jl` the convergence and sensitivity runs; `sint800_comparison.jl` the Tier 6 analysis of the Sint-800 record; `km_streaming.jl` the streamed drift and diffusion of the dipole amplitude). |
| `data/` | The Sint-800 record of Guyodo & Valet (1999) as published in their supplementary information, with its provenance (`data/README.md`). |
| `viz/` | Python plotting and visualisation (a uv project: matplotlib, numpy, h5py). Reads the package's HDF5 output. Not needed to simulate. |
| `dev/` | Development notes: model summary (`DEVELOPMENT.md`), validation gates (`VALIDATION.md`), open questions (`QUESTIONS.md`), implementation decisions (`DECISIONS.md`), discrepancies between the papers and this implementation (`DISCREPANCIES.md`), stretch goals outside the build order (`STRETCH_GOALS.md`). |
| `resources/` | Verified, agent-readable summaries of the three papers: every equation, parameter, table and result with page references. Start at `resources/README.md`. |

## Using the package

```julia
] add https://github.com/<owner>/hoyng-dynamo    # once pushed
```

Run the default test suite (validation Tiers 0–3 and 6, about 4.5 minutes) from the repository root:

```sh
julia --project -e 'using Pkg; Pkg.test()'
HOYNG_TIERS=0,1,2,3,4,6 JULIA_NUM_THREADS=4 julia --project -e 'using Pkg; Pkg.test()'   # include the slow Tier 4 (seeds run in parallel)
julia bench/run.jl                                                 # hot-path benchmarks
```

Streamed drift and diffusion of any scalar series (the online kernel-based regression of Davis 2023, `ConditionalMoments`), here of the dipole amplitude of a run:

```julia
using HoyngDynamo
spec = RunSpec(; f = 6.4)
obs = sim -> merge(standard_observers(sim), drift_diffusion_observers(sim))   # a every n_c steps, lags to 2 τ_d
e = ensemble(spec, 5e4; seeds = 1:3, observers = obs)                        # the moments are pooled over the seeds
d = kramers_moyal(e.pooled.km; lags = 5:20)                                   # D⁽¹⁾(a), D⁽²⁾(a) with error bars over lags 0.25–1 τ_d
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
| 4 | stochastic runs (⟨T_r⟩, variance, D₁ fit, sweeps, Poisson intervals) | standard run and forcing sweep passed 2026-10-20 (T4.1, T4.2, T4.5, T4.6, T4.11, T4.18; `HOYNG_TIERS=4`, ~5 min on 4 threads); the D₁-fit intercept and c² miss the 10% gate at f ≤ 6.4 (`@test_broken`, DR-5); convergence T4.17, the HSO runs T4.8, T4.9 (the χ² histogram fit reproduces HSO Table 1) and the SOH runs T4.10, T4.15 (variants and ramp) passed; the event diagnostics passed on durations (T4.14) and record the parity (T4.13, median −0.82 outside events against the paper's "≃ −1", DR-7); ~9 min on 4 threads |
| 5 | figures of the three papers, judged by eye | all thirteen figures rendered and judged comparable by the user (HOS Figs. 1, 2, 8, 9 on 2026-10-20; HOS Figs. 3–6, HSO Figs. 3, 6, 7, 8 and SOH Fig. 3 on 2026-10-24) |
| 6 | Sint-800 comparison (dv, the VADM χ² fits, the autocorrelation, τ_d by overlap) | passed 2026-10-24 (T6.1–T6.4, and T6.5 for SOH's fit; HOS's 35-bin fit recorded, DR-9); 25 tests in the default suite; `data/sint800.tsv`, `scripts/sint800_comparison.jl`, `viz/hso_fig01.py`, `sint800_fits.py`, `hso_fig07.py`, `sint800_tau_d.py` |

## The papers

- Hoyng, P., Ossendrijver, M. A. J. H. & Schmitt, D. (2001). The geodynamo as a bistable oscillator. *Geophys. Astrophys. Fluid Dyn.* 94, 263–314. doi:10.1080/03091920108203410
- Schmitt, D., Ossendrijver, M. A. J. H. & Hoyng, P. (2001). Magnetic field reversals and secular variation in a bistable geodynamo model. *Phys. Earth Planet. Inter.* 125, 119–124. doi:10.1016/S0031-9201(01)00237-0
- Hoyng, P., Schmitt, D. & Ossendrijver, M. A. J. H. (2002). A theoretical analysis of the observed variability of the geomagnetic dipole field. *Phys. Earth Planet. Inter.* 130, 143–157. doi:10.1016/S0031-9201(02)00004-3

## Licence

MIT, see `LICENSE`.
