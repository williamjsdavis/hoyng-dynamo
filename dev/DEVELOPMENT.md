# hoyng-dynamo: development notes

Working document for a fast, modern reimplementation of the bistable-oscillator geodynamo of Hoyng, Ossendrijver & Schmitt (2001, "HOS"), Schmitt, Ossendrijver & Hoyng (2001, "SOH") and Hoyng, Schmitt & Ossendrijver (2002, "HSO"). The paper content lives in `resources/`; this file holds the model summary an implementer needs in one place. The validation gates are in `VALIDATION.md`, open research questions in `QUESTIONS.md` and implementation decisions in `DECISIONS.md`. Notation is HOS's unless stated. The SOH/HSO translation is in `resources/00-series-overview-and-shared-notation.md` §3.

**Status (2026-10-19).** Reading pass over `resources/` complete. Validation gates and operational definitions accepted by the user (`VALIDATION.md`). Build step 0 (scaffold) complete the same day: package `HoyngDynamo` (Julia 1.11) with tier-tagged tests, D14 quality tests, benchmarks, CI, MIT licence and the `viz/` uv project; the smoke test passes. Build step 1 (grid, operators, eigensolver, C₀; `src/grid.jl`, `src/eigen.jl`) done the same day: Tier 0 passed, 57 tests, values in `VALIDATION.md` §4. Build step 2 (`src/stepper.jl`: `Model`, `State`, `Cache`, `ExponentialStepper`, `IMEXStepper`, allocation-free `step!`) done the same day: Tier 1 passed, 42 tests; gate T1.4 caught a first-order error in the frozen quenching coefficient, fixed with a midpoint predictor (`DECISIONS.md` D3). Build step 3 (`src/noise.jl`: `BinnedNoise` with the three F distributions, `ConstantForcing`/`LinearRamp`, SplitMix64 `run_seed`, the V₀₀ diagnostic) done 2026-10-19: Tier 2 passed, 362 tests; the binned noise's D₀ is 5% below the δ-correlated Eq. 36 (`VALIDATION.md` §4). Build step 4 (`src/aggregators.jl`: `Moments`, `Histogram`, `AutoCorrelation`, `Sampled`, `MovingAverage`, `SchmittTrigger`, `CrossingGroups`, `ExclusionVariance`, `BinnedMean`, bespoke, behind `fit!`/`merge!`/`value`) done 2026-10-19: 150 unit tests on synthetic series against offline references, every `fit!` allocation-free; the tests caught two detector defects (a jump across the band in one sample, and a crossing on the sample that enters the band). Build step 5 (`src/theory.jl`: `BistableTheory` with κ = μ/Λ − 1 selecting the HOS forms or HSO's A.1–A.2, the closed forms Eqs. 52–60, log-space quadratures for the normalisation, dv and the escape integral Eq. 57; QuadGK added as a dependency) done 2026-10-20: Tier 3 passed, 301 tests. Two findings recorded in `VALIDATION.md` §4: Table II's run-5 variance (printed 1.0) is not reproducible from Eq. 56 (10.4, the run being outside the approximation's validity by the paper's own footnote), and HSO case 4's ⟨T_r⟩ is 1.1% off within the rounding of the printed (γ, c²). The exact escape integral exceeds Eq. 58 by 8–18% for the Table II runs. Build step 6 (`src/run.jl`: `RunSpec`, `Simulation`, `Probe` observers, `standard_observers`, `run!`, `summarize`; `src/io.jl`: HDF5 run files with streams, results and bitwise checkpoint/resume; `scripts/run_standard.jl`; `viz/` figure scripts) done 2026-10-20: record and stream modes agree bitwise and a recorded a(τ) replayed through the aggregators reproduces the online values (the rest of T2.5); the standard run passes T4.1–T4.4 (pooled three seeds: ⟨T_r⟩ 104.7 ± 2.8 against 107 ± 5, variance 0.164 against 0.16, D₁/ℱ = 1.46 + 2.08a² against 1.36 + 2.00a², (c², γ) = (0.261, 2.93) against (0.246, 2.96)); HOS Figs. 1, 2, 8, 9 rendered for the user's judgement (T5.1, T5.2, T5.7, T5.8), with Fig. 9's discrepancy reproduced. Two findings: the D₁ intercept depends on the fit weighting (`DISCREPANCIES.md` DR-5), and the overtone amplitude scale is a normalisation convention, now unit toroidal energy for every mode (`DECISIONS.md` D4). Build step 7, piece 1 (the forcing sweep; `ensemble` and `pool!` in `src/run.jl`, the KS test in `src/stats.jl`, `scripts/run_sweep.jl`, HOS Figs. 5 and 6) done 2026-10-20: ⟨T_r⟩ reproduces HOS Table II at all five forcings, from 19 to 11 000 τ_d, within the paper's errors; the variance column is reproduced by dv at every forcing and by the §2 estimator within 15%; the polarity intervals are exponential at every forcing (T4.11); the (c², γ) of the D₁ fit reproduce at f ≥ 7.4 and depend on the fit's weighting at f ≤ 6.4, where the intercept is an extrapolation (T4.7 on c² recorded as `@test_broken` there). Three corrections on the way: the HOS Fig. 8 caption says the D₁ fit is count-weighted, which the step-6 gate had missed (`DISCREPANCIES.md` DR-5, the standard run's intercept and c² now fail the 10% gate); the exclusion-variance aggregator pooled per-interval ratios and was dominated by near-zero-mean intervals at high forcing (DR-6, now pooled moments); and the forcing-parameter invariance (T4.12) is broken by the bin width as DR-3 predicts, the N_c = 20 run having a 27% shorter ⟨T_r⟩ and an 18% larger D₁ intercept. The convergence gate T4.17 passed the same day (Δτ = 0.005 and N = 120 within 1σ), and `scripts/dr5_sensitivity.jl` traced the DR-5 gap to the time integrator: the noise distribution and the resolution leave the D₁ intercept unchanged, a first-order IMEX scheme at the papers' Δτ moves it by 0.26, the size of the gap, so the papers' unstated scheme is the probable cause. Step 7 piece 2 (the HSO runs) done 2026-10-20: the χ² fit of Eq. 10/A.1 to the |a| histogram (`density_fit`, `src/stats.jl`, in `summarize`) reproduces HSO Table 1 case 2 to 2% ((7.88, 0.509) against (8, 0.50)) and lands on their cases 4 and 5; T4.8 and T4.9 pass (reference run 526 ± 31 against 570 ± 61, dv 0.073 against 0.070; second run 115 ± 3 against 113 ± 5, dv 0.079 against 0.081); the autocorrelation of |a| decays at 0.87 per τ_d, near 2Λ and far from 2μ (QUESTIONS.md Q5 answered); HSO Figs. 3, 6, 7, 8 rendered. Step 7 piece 3 (the SOH runs) done 2026-10-20: `segmented` runs consecutive segments with fresh observers for windowed statistics; the f = 6 and 7 variants give ⟨T_r⟩ 185 ± 7 and 53 ± 1 with variances 0.11–0.12 and 0.21–0.24 against SOH's ≈ 200, ≈ 50, 0.11, 0.24; the ramp 4 → 6.5 over 5×10⁴ τ_d raises the rms variability from 16% to 39% (SOH 17% → 41%) and the reversals per 5000 τ_d from 0 to 149; T4.10 and T4.15 pass; SOH Fig. 3 rendered (`scripts/run_ramp.jl`). Step 7 piece 4 (the event-level diagnostics) done 2026-10-20: `Series` and `outside_events` (`src/aggregators.jl`) apply the trigger's event list to the stored parity series, `window_run` (`src/run.jl`) re-runs a window from the seed at full resolution and `scripts/run_events.jl` selects HOS's four event types by documented rules. T4.14 passes: the reversal duration has a median of 2.6 τ_d with the §2 band (HOS's "~3 τ_d typical", 19% of the reversals faster than 0.5 τ_d) and 0.8 τ_d with a_h = 0.3 (HSO's "0.5–2 τ_d"), so the two papers quote different definitions of one distribution; Q9 answered across the sweep (the median falls from 4.1 to 0.8 τ_d as f rises while the mean stays near 3–4). T4.13's median parity outside the event windows is −0.82, not below −0.9: the exclusion barely moves the whole-run value because events are rare, the level is what the inter-event |a₁| ≈ 0.2 implies, and HOS's "𝒫 ≃ −1" is a loose reading of their figure (`DISCREPANCIES.md` DR-7; recorded as `@test_broken`). HOS Figs. 3 and 4 rendered at full resolution (T5.3, T5.4). The planned validation is complete: every Tier 0–4 gate has passed or is recorded with its discrepancy, and every Tier 5 figure is rendered. Next: the stretch goals (`STRETCH_GOALS.md`, SG-1 integrator sensitivity, SG-2 alternative quenching), if the user wants them. All implementation decisions closed (`DECISIONS.md`): Julia package, type-stable and allocation-free hot path, two run modes (record and stream) over one stepping loop, running aggregators behind a fit!/merge!/value interface (OnlineStats.jl optional), HDF5 output, Python `viz/` for plots, MIT licence.

Files in `dev/`:

| File | Holds |
|---|---|
| `DEVELOPMENT.md` | This file: model summary, gaps and how they are closed, performance notes, build order, decisions log. |
| `VALIDATION.md` | Accepted validation gates, Tiers 0–6, with the tolerance rules and operational definitions they depend on. |
| `QUESTIONS.md` | Open research questions the implementation can settle. |
| `DECISIONS.md` | Implementation decisions (how, not what): options, recommendation, status. |
| `DISCREPANCIES.md` | Places where the papers' numbers cannot be reproduced from their own formulae, or where our implementation knowingly differs; one entry per finding. |
| `STRETCH_GOALS.md` | Work outside the build order, to be taken up after steps 0–7: the integrator-sensitivity test of the DR-5 hypothesis (SG-1) and the alternative quenching laws (SG-2, gate T4.16). |

Marking convention: anything labelled *(derived here)* is our own derivation or inference, not a statement in the papers. Everything else is traceable to the resource files, with journal page or equation numbers.

---

## 1. What the papers did

### 1.1 The model

A mean-field αΩ dynamo reduced to one spatial coordinate, colatitude θ ∈ [0, π]. The radial dependence is removed by a local spherical-wave ansatz cos(kr)/(kr) with kR = 0.5, which leaves a constant extra damping −(kR)² in the diffusion operator. The state is two real scalar fields, the poloidal potential P(θ, τ) and the toroidal field T(θ, τ), with P = T = 0 at both poles. Time τ is in units of the turbulent diffusion time τ_d = R²/β.

Equations actually integrated (HOS Eqs. 6–8, 21, 26–27; pp. 267–273):

```
∂P/∂τ = L P + [ (1 + s(1 − a²)) cos θ + δα/α₀ ] T
∂T/∂τ = C₀ ∂θ( sin θ P ) + L T

L      = (1/sin θ) ∂θ( sin θ ∂θ ) − 1/sin²θ − (kR)²
a(τ)   = ℜ (b̃₀, b) = ∫₀^π sin θ ( P̃₀ P + T̃₀ T ) dθ        (b̃₀ = adjoint fundamental mode at C₀)
δα/α₀  = f F(θ, τ) / √(2 N_c sin θ)
```

Three ingredients:

- **Linear operator G at the critical dynamo number C₀.** Diffusion L on both fields, α-effect cos θ feeding P from T, constant radial shear feeding T from P via C ∂θ(sin θ P). At kR = 0.5, C₀ = 57.9. The fundamental mode is a steady dipole (λ⁰ = 0 at C₀, T antisymmetric about the equator). Every overtone is oscillatory and damped, alternating in parity, quadrupole first (Table I, p. 270).
- **Global α-quenching on the dipole energy.** α_m → α_m(1 − q₀a²). Expanded about the critical state this is the term s(1 − a²) multiplying cos θ, with s = C_m/C₀ − 1 = 0.728 for the operating value C_m = 100 and q₀ = s/(1+s) = 0.421. Because the expansion basis is the critical operator's, the equilibria are exactly a = ±1 with all overtones zero (p. 273).
- **Stochastic α.** The azimuthal average runs over 2N_c sin θ convective cells per latitude circle, so α has a random part ∝ (sin θ)^(−1/2). F(θ, τ) has zero mean and unit variance, is constant within each of N_c = 10 latitude bins of width π/N_c and within each interval τ_c = 0.05, and is redrawn independently per bin at the end of each interval. HOS Appendix E proves all statistics depend on (f, τ_c, N_c) only through f²τ_c/N_c². Standard run f = 6.4, giving δα_rms/(α₀ cos θ) ≳ 2.5 (p. 296).

Parameters of the runs that appear in the series (τ_c = 0.05, N_c = 10 throughout):

| Run | kR | C_m | f | f²τ_c/N_c² | Λ | μ/Λ | Where |
|---|---|---|---|---|---|---|---|
| HOS standard (#3) | 0.5 | 100 | 6.4 | 2.05×10⁻² | 0.527 | 1.45 | HOS Table II |
| HOS #1, #2, #4, #5 | 0.5 | 100 | 4.4, 5.4, 7.4, 8.4 | 0.97, 1.46, 2.74, 3.53 ×10⁻² | 0.527 | 1.45 | HOS Table II |
| HSO reference (= HOS #2) | 0.5 | 100 | 5.4 | 1.46×10⁻² | 0.527 | 1.45 | HSO Fig. 3 |
| HSO second run | 0.3 | 120 | 7.0 | 2.45×10⁻² | 0.815 | 1.82 | HSO Fig. 8, p. 156 |
| SOH variants | 0.5 | 100 | 6, 7 | 1.80, 2.45 ×10⁻² | 0.527 | 1.45 | SOH p. 121 |
| SOH ramp | 0.5 | 100 | 4 → 6.5 linear | 0.80 → 2.11 ×10⁻² | 0.527 | 1.45 | SOH Fig. 3 |

Derived constants for kR = 0.5, C_m = 100: E₀₀ = 1.048, μ = sE₀₀ = 0.763, Λ = λ⁰(C_m) = 0.527, ∫ sin θ (P̃₀T₀)² dθ = 3.529, F ≡ πf²τ_c/2N_c² = 3.22×10⁻² at f = 6.4, D₀ = 3.529 F = 0.114, (D₀τ_c)^(1/2) = 0.075.

Eigenvalues (HOS Table I, p. 270):

| i | parity of T | λⁱ at C₀ = 57.9 | λⁱ at C_m = 100 |
|---|---|---|---|
| 0 | antisymmetric (dipole) | 0 | 0.527 |
| 1 | symmetric (quadrupole) | −2.76 + 4.41i | −2.16 + 5.88i |
| 2 | antisymmetric | −8.67 + 4.97i | −8.21 + 6.04i |
| 3 | symmetric | −15.9 + 6.13i | −15.7 + 7.70i |

*(derived here)* The spectrum of G(C₀) + sE equals the spectrum of G(C_m), because for a 2×2 block operator with equal diagonal blocks the eigenvalues depend on the off-diagonal entries only through their product, and (1 + s)C₀ = C_m. This is why Λ is both "λ⁰ at C_m" and the small-amplitude growth rate of the model equation.

### 1.2 Numerics as reported

60 equidistant θ points on [0, π], Δτ = 10⁻², runs of 5×10⁴ τ_d (5×10⁶ steps), each noise draw persisting 5 steps (HOS footnote 1 p. 273, Table II p. 292). The time integrator and the discretisation of L are not stated. The same random sequence was reused across f values when comparing runs (Fig. 5, p. 280). Mode amplitudes a_k are extracted by projecting on the adjoint modes at C₀ (Eq. 15).

### 1.3 The theory

Projecting the PDE on the eigenmodes gives coupled stochastic mode equations (HOS Eq. 29). Van Kampen's short-correlation-time expansion, followed by integrating out the overtones, yields a Fokker–Planck equation for a alone (Eqs. 32–40):

```
∂p/∂τ = −∂(Sp)/∂a + ½ ∂²(Dp)/∂a²
S ≃ Λ(1 − a²)a                      (valid for |a| ≲ 0.4; next term (μ−Λ)(1−a²)a³, HSO App. A)
D = D₀a² + D₁(a) ≃ Δ₀a² + Δ₁,       Δ₀ = D₀ + ½D₁''(0),  Δ₁ = D₁(0)
U ≃ ¼Λ(1 − a²)²
p(a) ∝ (1 + a²/c²)^(γc²+γ−1) exp(−γa²),   γ = Λ/Δ₀,  c² = Δ₁/Δ₀,  bistable iff γ > 1
```

D₀a² is the direct kick of the noise on the dipole and can never flip it (multiplicative noise vanishes at a = 0). Δ₁ = D₁(0) comes from the mean-square overtone toroidal field at a = 0 and is the sole cause of reversals (Eq. 49). Reversals are therefore Kramers escapes with f²τ_c/N_c² as effective temperature. Closed forms (Eqs. 52–60):

```
a_m   = ±(1 − 1/γ)^(1/2)
var   ≡ ⟨(a−⟨a⟩)²⟩/⟨a⟩² ≃ (γc² + γ − 1) / [4(γ − 1)²]
ψ_m   = 1 − γ + (γc² + γ − 1) ln[(γc² + γ − 1)/(γc²)]
⟨T_r⟩ ≃ (π/Λ) (2/c²)^(1/2) [(γc² + γ − 1)/(γ − 1)]^(1/2) exp(ψ_m)
      ~ (π√2/Λ) [c²/(1+c²)]^(1/2) exp(Kγ),   K = −1 + (1 + c²) ln[(1 + c²)/c²]
```

Δ₁ cannot be computed from the model parameters (HOS App. E). The papers measure D₁(a) from the run (every 1 τ_d evaluate F ∫ sin θ P̃₀² T₊² dθ with T₊ = ℜ(T − a₀T₀), bin in |a|, fit a quadratic for |a| ≤ 1.1) and then check consistency. Standard run: D₁/F = 1.36 + 2.00a², hence c² = 0.246, γ = 2.96 (Eqs. 61–63).

### 1.4 Findings

- **Simulations** (HOS §3). Dipolar field whose amplitude wobbles around ±1, sudden reversals during which the field is small and nearly quadrupolar, aborted reversals and excursions, Poisson-distributed polarity intervals, mean interval exponentially sensitive to the forcing. Fast events ~0.5 τ_d, a typical reversal ~3 τ_d (p. 277).
- **Theory vs simulation** (HOS Table II, p. 292):

| run | f | c² | γ | var (Eq. 56) | var (run) | ⟨T_r⟩ (Eq. 58) | ⟨T_r⟩ (run) |
|---|---|---|---|---|---|---|---|
| 1 | 4.4 | 0.318 | 9.04 | 0.04 | 0.03 | 12000 | 10700 ± 5300 |
| 2 | 5.4 | 0.299 | 5.15 | 0.08 | 0.07 | 480 | 519 ± 53 |
| 3 | 6.4 | 0.246 | 2.96 | 0.17 | 0.16 | 94 | 107 ± 5 |
| 4 | 7.4 | 0.240 | 1.82 | 0.47 | 0.30 | 35 | 37.6 ± 1 |
| 5 | 8.4 | 0.228 | 1.09 | 1.0 | 0.50 | 35 | 18.7 ± 0.4 |

- **SOH** adds a ramped run (f from 4 to 6.5 over 5×10⁴ τ_d) reproducing the long-term rise of the reversal rate, and the first Sint-800 fit (γ = 5.2, c² = 0.26).
- **HSO** does the quantitative data work: χ² fits of p(a) to the Sint-800 VADM histogram (optimum γ = 6.3, c² = 0.32, 1σ), the secular-variation measure dv = ⟨(|a|−⟨|a|⟩)²⟩/⟨|a|⟩² (data 0.071 ± 0.003, reference run 0.070), and an autocorrelation method giving τ_d = 11 ± 2 kyr (superseding HOS's 3 kyr). With that τ_d the reference run reverses too rarely (~6 Myr vs. 95% interval 0.22–3.3 Myr) and the second run has reversals a factor ~2 too long (8–30 kyr). The theoretical ⟨T_r⟩ overshoots the measured 570 ± 61 τ_d by almost a factor 2 for the reference run (HSO p. 151).
- **Time scales differ across papers.** SOH assumes 2 kyr, HOS infers 3 kyr, HSO measures 11 ± 2 kyr. Never carry years between papers; work in τ_d.

---

## 2. Gaps the papers leave and how each is closed

| # | Gap | Where flagged | Resolution (accepted 2026-10-18 unless noted) |
|---|---|---|---|
| G1 | Time integrator and finite-difference form of L and the shear term never stated. On 60 points the diffusion stiffness is ~4/Δθ² ≈ 1.4×10³ per τ_d, above the explicit-Euler limit 2/Δτ = 200, so some implicit or sub-stepped treatment is implied. | HOS §11 | Implementation choice, not a convention: see `DECISIONS.md` (discretisation of L, time integrator). Verified by Tiers 0–1 and the Δτ sweeps T1.6 and T4.17. |
| G2 | Distribution of F unspecified beyond zero mean, unit variance. | HOS §9, SOH §4.4, HSO §4.4 | Default Gaussian; uniform and ±1 as sensitivity options (`QUESTIONS.md` Q7). |
| G3 | Alignment of the 10 latitude bins with the 60 grid points. | HOS §11 | Assign each grid point to bin ⌊θ_j N_c/π⌋; document. |
| G4 | Reversal-counting rule: "odd number of zero crossings within a few τ_d" with no window stated. | HOS footnote 2 p. 276 | `VALIDATION.md` §2: Schmitt-trigger with threshold a_h (default 0.5) plus the paper's grouping rule for comparison. |
| G5 | How reversal transitions are excised from the between-reversal variance; HSO's dv instead uses the rectified whole record. | HOS §11, SOH §9 | `VALIDATION.md` §2: implement both: HOS-style with a ±W τ_d exclusion window (default W = 5) and HSO dv on |a|. |
| G6 | Real or complex fields. | HOS §11 | Keep P, T real (reality is preserved by the model equation). Overtone amplitudes a_k are complex because the adjoint modes are. |
| G7 | Normalisation of b₀ (sets the field unit only). | HOS §11 | Take b₀ real, bi-orthonormal; a is measured relative to it. |
| G8 | 1/√sin θ at the poles. | HOS §9 | Evaluate only at interior points; enforce P = T = 0 at the endpoints. |
| G9 | Sampling cadence of a(τ) for histograms and for the D₁ measurement. | HSO footnote 4 | Store a at every step in binary; histograms at 1 τ_d spacing for like-for-like with the papers. |


---

## 3. Performance notes

- Unknowns: 2(N−2) = 116 for N = 60 nodes including the poles. Steps per paper-length run: 5×10⁶. Per step: a few 60-vector operations, one projection for a, one noise redraw every 5 steps, and either two tridiagonal solves (IMEX) or one dense 120×120 matrix–vector product (exponential integrator); see `DECISIONS.md`. Order 10³–10⁴ flops per step, so at most ~5×10¹⁰ flops per run: seconds in compiled code, minutes in naive interpreted code.
- Everything in Tier 4 is embarrassingly parallel over seeds and forcing values.
- The eigenproblem is a one-off 120×120 (or larger for convergence) dense non-symmetric solve, plus its adjoint (left eigenvectors with the sin θ weight).
- Store a(τ) at full cadence: 5×10⁶ doubles = 40 MB per run; fine.

---

## 4. Suggested build order

0. Scaffold: Julia package (Project.toml, `src/`, `test/`, `bench/`), CI workflow, MIT licence, `viz/` uv project, `.gitignore` for `results/`. **Done 2026-10-18.**
1. Grid, L, shear operator, G(C). Eigen/adjoint solver with bi-orthonormalisation. Gates T0.1–T0.11. **Done 2026-10-18, Tier 0 passed.**
2. Deterministic integrator (f = 0). Gates T1.1–T1.6. Fix the integrator here. **Done 2026-10-18, Tier 1 passed; exponential stepper with midpoint quenching is the production integrator.**
3. Noise generator and V₀₀ diagnostic. Gates T2.1–T2.5. **Done 2026-10-19, Tier 2 passed; the record/stream and aggregator-replay parts of T2.5 are tested at step 6.**
4. Running aggregators (`DECISIONS.md` D8): interface, Welford/histogram/event detector/exclusion variance/autocorrelation, unit-tested on synthetic series with known answers. Independent of 1–3. **Done 2026-10-19, 150 tests, bespoke implementation (no OnlineStats.jl).**
5. Formula module (Tier 3), independent of 1–4; can be done first or in parallel. **Done 2026-10-20, Tier 3 passed, 301 tests.**
6. Standard stochastic run in both modes with HDF5 output (a, a_k, 𝒫, energy, D₁ measurement, reversal detection). Gates T4.1–T4.4; figures T5.1–T5.8 from `viz/`. **Done 2026-10-20: T4.1–T4.4 passed; Figs. 1, 2, 8, 9 rendered (T5.1, T5.2, T5.7, T5.8), the other Tier 5 figures need the step 7 runs.**
7. Sweeps and second parameter set. Remaining Tier 4 and Tier 5. Stretch goals (`STRETCH_GOALS.md`: the integrator-sensitivity test SG-1, alternative quenching SG-2 = T4.16) after everything else. **Piece 1 (forcing sweep: T4.5–T4.7, T4.11, T4.12, T4.18; Figs. 5, 6), piece 5 (convergence, T4.17), piece 2 (HSO runs: T4.8, T4.9, T5.13; HSO Figs. 3, 6, 7, 8) and piece 3 (SOH runs: T4.10, T4.15; SOH Fig. 3) done 2026-10-20.** **Piece 4 (event-level diagnostics: T4.13, T4.14, Q9; HOS Figs. 3, 4) done 2026-10-20. Step 7 complete.**
8. Sint-800 (Tier 6): deferred, not planned.

---

## 5. Decisions log

| Date | Decision | Status |
|---|---|---|
| 2026-10-18 | Validation gates, Tiers 0–6 (`VALIDATION.md` §3) | accepted |
| 2026-10-18 | Tolerance rules: 1% deterministic, 2σ combined for run statistics (`VALIDATION.md` §1) | accepted |
| 2026-10-18 | Operational definitions: Schmitt-trigger reversals a_h = 0.5 plus the papers' grouping rule; between-reversal variance with W = 5 τ_d exclusion; HSO dv on the rectified record; Gaussian F; bin ⌊θ_j N_c/π⌋ (`VALIDATION.md` §2) | accepted |
| 2026-10-18 | Open questions Q1–Q9 moved to `QUESTIONS.md` | recorded |
| 2026-10-18 | Tier 6 (Sint-800) deferred, not planned for now (`DECISIONS.md` S1) | accepted |
| 2026-10-18 | Alternative quenching laws (T4.16) are a stretch goal after Tiers 0–5 (`DECISIONS.md` S2; now `STRETCH_GOALS.md` SG-2) | accepted |
| 2026-10-20 | Stretch goals live in `STRETCH_GOALS.md`; the integrator-sensitivity test of the DR-5 hypothesis (SG-1) is added to them | accepted |
| 2026-10-18 | Deliverable: personal research project, built as a package others can install later (`DECISIONS.md` S3) | accepted |
| 2026-10-18 | Language: Julia, as a package; code must be type stable (`DECISIONS.md` D1, D14) | accepted |
| 2026-10-18 | Constant-memory integrator and running aggregators for arbitrarily long runs (`DECISIONS.md` D3, D8) | accepted as requirement |
| 2026-10-18 | Storage: HDF5 per run with streamed datasets (`DECISIONS.md` D7) | accepted |
| 2026-10-18 | Testing by tier and CI on the default suite (`DECISIONS.md` D11) | accepted |
| 2026-10-18 | Two run modes, record (time series to HDF5) and stream (aggregators only), identical results for the same seed (`DECISIONS.md` D7) | accepted |
| 2026-10-18 | Plotting in Python: `viz/` with a uv environment, not a simulation dependency (`DECISIONS.md` D13) | accepted |
| 2026-10-18 | Aggregators behind a fit!/merge!/value interface; OnlineStats.jl optional, allocation-checked if used (`DECISIONS.md` D8) | accepted |
| 2026-10-18 | Licence MIT (`DECISIONS.md` D10) | accepted |
| 2026-10-18 | D2–D6, D9, D12 (discretisation, integrator, eigensolver, RNG, state, parallelism, theory numerics) | adopted as recommended, no objection raised |
