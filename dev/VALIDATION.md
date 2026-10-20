# Validation gates

Accepted by the user on 2026-10-18: the gates, the tolerance rules and the operational definitions below. Targets come from the three papers via `resources/`; notation is HOS's unless stated; *(derived here)* marks our own derivation rather than a statement in a paper. Gate IDs (T0.1, …) are the ones referenced from `QUESTIONS.md` and `DECISIONS.md`.

Status: Tiers 0–3 passed (Tier 0 and Tier 1 on 2026-10-18, Tier 2 on 2026-10-19, Tier 3 on 2026-10-20); results in §4. Tiers 4–5 not yet run.

## 1. Conventions

- **Tiers** are ordered: pass a tier before spending effort on the next.
- **Kind.** A = automated numeric check with a pass/fail criterion. D = diagnostic: computed and recorded, no pass/fail, feeds `QUESTIONS.md`. J = qualitative judgement by the user, against the paper's figure.
- **Tolerances.** For deterministic quantities printed to n significant figures, pass if within ±1 in the last printed digit, relaxed to 1% because the papers' own discretisation is unknown; record both the 60-point and a grid-converged value. For run statistics, pass if |ours − theirs| < 2·√(σ_theirs² + σ_ours²), where σ_ours is from repeated seeds. For convention-dependent statistics (reversal counts, between-reversal variance) the convention used must be stated with the result.
- **Hard vs soft.** Tiers 0, 1 and 3 are hard gates: a failure blocks progress. Tiers 2, 4 are hard for A items and informative for D items. Tier 5 passes on the user's say-so. Tier 6 is deferred until the Sint-800 data are obtained.

## 2. Operational definitions

These close gaps G2–G5 and G9 of `DEVELOPMENT.md` §2. Every reported statistic quotes the definition it used.

- **Reversal detection.** Schmitt trigger on a(τ): polarity is +1 until a < −a_h, then −1 until a > +a_h. Default a_h = 0.5. A polarity flip is a reversal; a zero crossing without a flip is an aborted reversal; a dip of |a| below a_e = 0.3 without a zero crossing is an excursion. Report sensitivity to a_h ∈ {0.3, 0.5, 0.7}. Also implement the papers' rule (group zero crossings within W = 5 τ_d; odd count = reversal, even = aborted) for comparison.
- **Between-reversal variance (HOS/SOH style).** Within each polarity interval, drop W τ_d after the preceding reversal and before the next (default W = 5), compute ⟨(a−⟨a⟩)²⟩/⟨a⟩² over the remainder, pool across intervals time-weighted. Sensitivity W ∈ {2, 5, 10}.
- **dv (HSO style).** Eq. 12 on the rectified whole record |a|, no exclusion.
- **⟨T_r⟩** = run length / number of reversals, error ⟨T_r⟩/√N (Poisson), as in the papers.
- **Noise.** F ~ N(0, 1) by default. Bins: grid point j belongs to bin ⌊θ_j N_c/π⌋. Redraw at steps that are multiples of τ_c/Δτ = 5.
- **Sampling.** a(τ) is stored at every step. Histograms of a and the D₁ measurement use 1 τ_d spacing, as in the papers (HSO footnote 4; HOS Fig. 8). Overtone amplitudes a₁, a₂ use 5-step spacing when plotted (HOS Fig. 2).

## 3. Gates

### Tier 0: linear operator and eigenproblem

| ID | Check | Target | Pass criterion | Source | Kind |
|---|---|---|---|---|---|
| T0.1 | Critical dynamo number at kR = 0.5: root of ℜλ⁰(C) = 0 | C₀ = 57.9 | within 0.1 on a converged grid; record 60-point value | HOS p. 269, Table I | A |
| T0.2 | Eigenvalues of G at C₀, modes 0–3 | 0; −2.76+4.41i; −8.67+4.97i; −15.9+6.13i | ℜ and ℑ each within 1%; \|λ⁰\| < 0.01 | Table I | A |
| T0.3 | Eigenvalues of G at C_m = 100, modes 0–3 | 0.527; −2.16+5.88i; −8.21+6.04i; −15.7+7.70i | within 1% | Table I | A |
| T0.4 | Parity pattern and periodicity | mode 0: T antisym., P sym., ℑλ⁰ = 0; mode 1: T sym.; parities alternate; ℑλⁱ ≠ 0 for i ≥ 1 | exact | Table I footnotes | A |
| T0.5 | Adjoint modes and bi-orthonormality (b̃ᵢ, bⱼ) = δᵢⱼ under the sin θ inner product, first ~10 modes | — | max residual < 10⁻¹⁰ | Eqs. 11–13 | A |
| T0.6 | Overlap integral E₀₀ = ∫ sin θ cos θ P̃₀T₀ dθ | 1.048 | ±0.005 | p. 282 | A |
| T0.7 | Diffusion integral ∫ sin θ (P̃₀T₀)² dθ | 3.529 | ±0.01 | p. 284 | A |
| T0.8 | Derived constants s = C_m/C₀ − 1, q₀ = s/(1+s), μ = sE₀₀ | 0.728, 0.421, 0.763 | ±0.002 | p. 272, Eq. 34 | A |
| T0.9 | Spectrum of G(C₀) + sE equals spectrum of G(C_m) | — | eigenvalues agree to 10⁻⁸ | *(derived here)*, see `DEVELOPMENT.md` §1.1 | A |
| T0.10 | Second parameter set kR = 0.3: locate C₀(0.3), then λ⁰ at C = 120 and μ/Λ = sE₀₀/Λ | Λ = 0.815, μ/Λ = 1.82 | 1% | HSO Fig. 8 caption, p. 156 | A |
| T0.11 | Grid convergence of T0.1–T0.7 at 60, 120, 240 points | — | report; converged values within the tolerances above | — | D |

### Tier 1: deterministic nonlinear dynamics (f = 0)

| ID | Check | Target | Pass criterion | Source | Kind |
|---|---|---|---|---|---|
| T1.1 | b = ±b₀ is stationary: start at b₀, integrate 100 τ_d | a ≡ 1 | \|a − 1\| < 10⁻⁸ (or integrator tolerance) | p. 273 | A |
| T1.2 | Small-amplitude growth from a random small field | rate Λ = 0.527 (kR 0.5, C_m 100); 0.815 (kR 0.3, C_m 120) | 1% after transients | Eq. 34; HSO | A |
| T1.3 | Saturation from arbitrary initial field | \|a\| → 1, overtones → 0 | \|a\| within 10⁻⁶ of 1 by 50 τ_d | p. 273 | A |
| T1.4 | Relaxation rate of a small perturbation about a = 1 | 2μ = 1.526 (kR 0.5); 2μ = 2.97 (kR 0.3, μ = 1.483) | 1% | *(derived here)*: linearising Eq. 26 at b₀ gives dδa/dτ = −2sE₀₀ δa; overtones are driven but do not feed back at linear order. Note HSO Eq. 14–15 quote 2Λ for the noisy autocorrelation; see `QUESTIONS.md` Q5. | A |
| T1.5 | Sign symmetry b → −b | solutions map to solutions | exact to round-off | p. 273 | A |
| T1.6 | Time-step convergence of T1.2 and T1.4: Δτ = 0.01 vs 0.001 | — | report the error at Δτ = 0.01 | — | D |

### Tier 2: noise construction

| ID | Check | Target | Pass criterion | Source | Kind |
|---|---|---|---|---|---|
| T2.1 | F sample statistics: zero mean, unit variance, constant for 5 steps, independent across bins and intervals | — | mean within 4σ of 0 and variance within 1% at 10⁶ draws; successive-interval and cross-bin correlations consistent with 0 (4σ rather than 3σ because three distributions are tested at one fixed seed; see §4) | Eq. 21, p. 272 | A |
| T2.2 | Measured correlation function of δα/α₀ matches Eq. B.6: triangular in lag and in angle, amplitude f²/(2N_c √(sin θ sin θ')) | — | within 2% at zero lag | B.6, p. 308 | A |
| T2.3 | D₀ measured directly: Var[∫₀ᵀ V₀₀ dτ]/T for T ≫ τ_c, with V₀₀ from Eq. 31 on the frozen mode b₀ | 0.114 = 3.529 F at f = 6.4 | report; expect agreement within ~10% (δ-correlation of B.7 vs 10 finite bins) | Eqs. 36, B.7–B.8 | D |
| T2.4 | Short-correlation-time parameter (D₀τ_c)^(1/2) | 0.075 | from T2.3 | C.3, p. 310 | A |
| T2.5 | Reproducibility: seeded RNG reproduces a(τ) bitwise; one noise series can be reused across f values; record and stream modes (`DECISIONS.md` D7) give identical states and aggregator values, and replaying a recorded a(τ) through the aggregators reproduces the online values | — | bitwise | p. 280; D7 | A |

### Tier 3: closed-form theory (formula module, independent of the simulation)

| ID | Check | Target | Pass criterion | Source | Kind |
|---|---|---|---|---|---|
| T3.1 | Eq. 43 equals Eq. 42 evaluated numerically with S = Λ(1−a²)a, D = Δ₀a² + Δ₁ | — | 10⁻⁸ after normalisation | Eqs. 42–44 | A |
| T3.2 | Peak a_m, ψ_m, ψ_m'' (Eqs. 52–54) vs numerical maximisation of ψ | — | 10⁻⁸ | Eqs. 50–54 | A |
| T3.3 | Variance Eq. 56 for the five Table II (c², γ) pairs | 0.04, 0.08, 0.17, 0.47, 1.0 | last printed digit. Run 5 (γ = 1.09) recorded only: Eq. 56 gives 10.4 and the table's footnote b puts the run outside the validity of the approximation (§4, 2026-10-20) | Table II | A |
| T3.4 | ⟨T_r⟩ Eq. 58 for the five Table II pairs | 12000, 480, 94, 35, 35 | 1% (hand check: 93.9 for run 3) | Table II | A |
| T3.5 | Eq. 58 with Λ = 0.527 for HSO Table 1 cases 1–3 (γ, c²) = (6.3, 0.32), (8, 0.50), (7, 0.45) | 1090, 928, 665 | 0.5% (hand check: 1089, 927, 665) | HSO Table 1 | A |
| T3.6 | dv (HSO Eq. 12) from Eq. 10 over a > 0, numerically, cases 1–3 | 0.076, 0.070, 0.079 | ±0.001 | HSO Table 1 | A |
| T3.7 | Eq. A.1 and the generalised A.2 with μ/Λ = 1.45 (case 4, Λ = 0.527) and 1.82 (case 5, Λ = 0.815) | ⟨T_r⟩ 1097, 138; dv 0.071, 0.091 | 2% and within the bracket of γ ± 0.05 (the printed γ, c² have two figures and ⟨T_r⟩ is exponentially sensitive; §4, 2026-10-20); ±0.001 | HSO Table 1, App. A | A |
| T3.8 | Numerical escape integral (Eq. 57, SOH's 4∫ form) vs the Gaussian-approximated Eq. 58 at (γ, c²) = (3, 0.2) and (6, 0.3) | — | report the ratio (quantifies the approximation in F.2) | Eq. 57, F.1–F.3, Fig. 11 | D |
| T3.9 | Asymptotic Eq. 59–60 vs Eq. 58 | ratio → 1 as γ → ∞ | monotone approach | Eqs. 59–60 | A |

### Tier 4: stochastic runs, quantitative

All runs: 60 points, Δτ = 0.01, 5×10⁴ τ_d unless stated; definitions of §2 stated with each result; at least three seeds for σ_ours.

| ID | Check | Target | Pass criterion | Source | Kind |
|---|---|---|---|---|---|
| T4.1 | Standard run (kR 0.5, C_m 100, f 6.4): mean time between reversals | 107 ± 5 τ_d (468 reversals) | 2σ combined | HOS p. 291, Table II | A |
| T4.2 | Standard run: between-reversal relative variance | 0.16 | ±0.02 | Table II | A |
| T4.3 | Standard run: binned D₁/F for \|a\| ≤ 1.1, quadratic fit | 1.36 + 2.00a² | coefficients within 10% | Eq. 61, Fig. 8 | A |
| T4.4 | Derived c², γ from T4.3 and D₀ | 0.246, 2.96 | 10% | Eq. 63 | A |
| T4.5 | Sweep f = 4.4, 5.4, 7.4, 8.4: ⟨T_r⟩ | 10700 ± 5300, 519 ± 53, 37.6 ± 1, 18.7 ± 0.4 | 2σ combined; for f = 4.4 run ≥ 5×10⁵ τ_d to beat the paper's error | Table II | A |
| T4.6 | Sweep: between-reversal variance | 0.03, 0.07, 0.30, 0.50 | ±0.02 or 15% | Table II | A |
| T4.7 | Sweep: (c², γ) from the D₁ fit | (0.318, 9.04), (0.299, 5.15), (0.240, 1.82), (0.228, 1.09) | 10% | Table II | A |
| T4.8 | HSO reference run f = 5.4: ⟨T_r⟩, dv, (γ, c²) from fitting Eq. 10 to the a-histogram | 570 ± 61 (88 reversals); 0.070; (8, 0.50) | 2σ; ±0.005; 10% | HSO p. 146, Table 1 case 2 | A |
| T4.9 | HSO second run kR 0.3, C_m 120, f 7: ⟨T_r⟩, dv | 113 ± 5 τ_d; 0.081 | 2σ; ±0.005 | HSO Table 2 | A |
| T4.10 | SOH variants f = 6 and 7: ⟨T_r⟩, variance | ≈200, 0.11; ≈50, 0.24 | 20% | SOH p. 121 | A |
| T4.11 | Polarity intervals are Poisson for T ≫ Λ⁻¹ | exponential | KS test p > 0.05 on intervals > 5 τ_d | HOS p. 290 | A |
| T4.12 | Forcing-parameter invariance at f²τ_c/N_c² = 2.05×10⁻²: (f, τ_c, N_c) = (6.4, 0.05, 10), (4.525, 0.1, 10), (12.8, 0.05, 20) | same ⟨T_r⟩ and variance | within statistical error (soft: finite τ_c and bin width break exact invariance) | p. 272, App. E | D |
| T4.13 | Parity 𝒫 ≃ −1 between reversals | — | median 𝒫 outside event windows < −0.9 | Eq. 28, p. 276 | D |
| T4.14 | Reversal durations | ~3 τ_d typical, fast events ~0.5 τ_d (HOS); 0.5–2 τ_d (HSO, different definition); quadrupole period ≈ 1.1 τ_d visible | report distribution with the stated definition | HOS p. 277; HSO p. 151 | D |
| T4.15 | SOH ramped run, f linear 4 → 6.5 over 5×10⁴ τ_d | rms variability 17% → 41%; chrons shorten drastically; first chron ≈ 1.9×10⁴ τ_d [approximate read-off] | 20% on the rms values; shape by user | SOH Fig. 3, p. 123 | A/J |
| T4.16 | Robustness regression: alternative quenching (q from total toroidal energy; α/(1+q); local quenching) | dipole statistics "not materially changed" | ⟨T_r⟩, variance within ~20% | HOS pp. 277, 281 | D, stretch goal (after Tiers 0–5) |
| T4.17 | Convergence of T4.1–T4.2 at 120 points and at Δτ = 0.005 | — | within statistical error | — | D |
| T4.18 | HOS run 2 vs HSO reference run, nominally identical: 519 ± 53 vs 570 ± 61 | — | our seed-to-seed scatter at f = 5.4 should cover both | Table II; HSO p. 146 | D |

### Tier 5: figures, qualitative (user is the decision gate)

Each item: reproduce the plot with the paper's axes and run parameters, place beside the paper's figure, user judges.

| ID | Figure | What to reproduce | Features to look for | Source |
|---|---|---|---|---|
| T5.1 | HOS Fig. 1 | ℜT(θ, τ) colour plot plus a(τ), 0–2000 τ_d, standard run, 1 τ_d spacing | vertical stripes (overtone fine structure); reversal bars every ~100 τ_d; a around ±1 with occasional \|a\| ≈ 2; aborted reversals and an excursion | p. 274, Plate I |
| T5.2 | HOS Fig. 2 | \|a₁\|, \|a₂\|, parity 𝒫, energy ∫(ℜT)² sin θ dθ over 200 τ_d, 5-step spacing | \|a₁\| spikes to ~1.2, \|a₂\| to ~0.9; 𝒫 near −1 with spikes to +1 at events; reversals in clusters | p. 275 |
| T5.3 | HOS Fig. 3 | two reversals at full resolution, 6 τ_d windows | one very fast event; one quasi-periodic quadrupolar event with period ≈ 1.1 τ_d; recovery ∝ exp(Λτ) | p. 278, Plate II |
| T5.4 | HOS Fig. 4 | an aborted reversal and an excursion, 6 τ_d windows | similar appearance, sign change vs none | p. 279, Plate III |
| T5.5 | HOS Fig. 5 | a(τ), 0–2000 τ_d, f = 5.4 and 7.4, same noise series as T5.1 | ~6 reversals vs several tens | p. 280 |
| T5.6 | HOS Fig. 6 | ℜa₁ vs a scatter, f = 5.4 | two clouds at a ≈ ±1 with \|ℜa₁\| ≲ 0.5 and a jump | p. 283 |
| T5.7 | HOS Fig. 8 | D₁/F scatter (one point per τ_d) vs \|a\|, bin-average histogram, quadratic fit | rises from ≈ 1.4 at a = 0; a⁴ term appearing beyond 1.1 | p. 292 |
| T5.8 | HOS Fig. 9 | simulated p(a) vs Eq. 43 at (γ, c²) = (2.96, 0.246), unit area | the paper's discrepancy must reappear: histogram peak ≈ 1.35 at \|a\| ≈ 0.85 vs curve ≈ 1.1 at 0.8 | p. 293 |
| T5.9 | HSO Fig. 3 | T(θ, τ) and a(τ), 0–3000 τ_d, f = 5.4 | ~7 sign changes; longest interval ≈ 1100 τ_d | HSO p. 146 |
| T5.10 | HSO Fig. 6 (bottom panels) | simulation a-histogram with Eq. 10 and Eq. A.1 overlaid | A.1 visibly better than Eq. 10 | HSO p. 150 |
| T5.11 | HSO Fig. 8 | polarity bar code to 10⁴ τ_d and a(τ) to 3000 τ_d, second run | ~22 sign changes in 3000 τ_d | HSO p. 153 |
| T5.12 | SOH Fig. 3 | polarity bar code of the ramped run, 0–5×10⁴ τ_d | long first chron, then rapidly shortening | SOH p. 123 |
| T5.13 | HSO Fig. 7 (model curves only) | autocorrelation of rectified, mean-subtracted a with running average, lags to ~5 τ_d | decay rate at small lag compared with 2Λ (HSO Eq. 15) and with 2μ (T1.4); anticorrelation at long lags | HSO p. 152 |

### Tier 6: Sint-800 data comparison (deferred by the user on 2026-10-18, not planned for now; needs the Guyodo & Valet 1999 record, not in the repository)

| ID | Check | Target | Source |
|---|---|---|---|
| T6.1 | dv of the 799-point VADM record with Monte-Carlo errors | 0.071 ± 0.003 | HSO p. 148 |
| T6.2 | χ² fit of Eq. 10 with ~100 equal-population bins | γ = 6.3, c² = 0.32, A_m = 5.98×10²² A m², 1σ | HSO Table 1 case 1 |
| T6.3 | Sint-800 autocorrelation shape | zero crossing ≈ 27 kyr, minimum ≈ −0.2 near 40 kyr | HSO Fig. 7 |
| T6.4 | τ_d by autocorrelation overlap (τ_av ≈ 5 kyr) | 11 ± 2 kyr (reference run); 16 ± 4 kyr (second run) | HSO pp. 152–154 |
| T6.5 | HOS Fig. 10 / SOH Fig. 4 fits (35 and 30 bins) | (5.15, 0.25) and (5.2, 0.26) within 1σ | HOS p. 295; SOH p. 123 |

## 4. Results

### Tier 0 (2026-10-18, commit after a5d0b8c): passed, 57/57 tests

Conservative second-order finite differences, poles as Dirichlet nodes (`DECISIONS.md` D2), dense LAPACK eigensolver with weighted left eigenvectors (D4). "Paper" is the printed value.

| Quantity | Paper | N = 60 | N = 120 | N = 240 | Gate |
|---|---|---|---|---|---|
| C₀ (kR = 0.5) | 57.9 | 57.897 | 57.884 | 57.881 | T0.1 pass |
| λ¹ at C₀ | −2.76 + 4.41i | −2.757 + 4.409i | −2.757 + 4.412i | −2.757 + 4.413i | T0.2 pass |
| λ² at C₀ | −8.67 + 4.97i | −8.663 + 4.964i | — | −8.670 + 4.974i | T0.2 pass |
| λ³ at C₀ | −15.9 + 6.13i | −15.845 + 6.108i | −15.873 + 6.126i | −15.880 + 6.130i | T0.2 pass |
| Λ = λ⁰ at C_m = 100 | 0.527 | 0.5283 | 0.5277 | 0.5275 | T0.3, T0.8 pass |
| λ¹ at C_m | −2.16 + 5.88i | −2.159 + 5.875i | −2.158 + 5.880i | −2.158 + 5.882i | T0.3 pass |
| λ² at C_m | −8.21 + 6.04i | −8.201 + 6.024i | — | −8.207 + 6.039i | T0.3 pass |
| λ³ at C_m | −15.7 + 7.70i | −15.644 + 7.672i | −15.670 + 7.696i | −15.676 + 7.702i | T0.3 pass |
| E₀₀ | 1.048 | 1.0491 | 1.0483 | 1.0481 | T0.6 pass |
| ∫ sin θ (P̃₀T₀)² dθ | 3.529 | 3.5323 | 3.5294 | 3.5288 | T0.7 pass |
| C₀ (kR = 0.3) | — | 49.722 | — | 49.702 | T0.10 (recorded) |
| Λ (kR = 0.3, C_m = 120) | 0.815 | 0.8163 | — | 0.8148 | T0.10 pass |
| μ/Λ (kR = 0.3, C_m = 120) | 1.82 | 1.818 | — | 1.821 | T0.10 pass |

Also passed: parity pattern and periodicity (T0.4), bi-orthonormality residual < 10⁻¹⁰ over ten modes at both C₀ and C_m (T0.5), s, q₀, μ within ±0.002 (T0.8), spectral equivalence of G(C₀) + sE and G(C_m) to 10⁻⁸ (T0.9).

Observations. The N = 240 values agree with the printed ones to the last printed digit in every case, so the papers' eigenvalues are grid-converged; our 60-point values differ from them by 0.1–0.4% (second-order error), which suggests the papers' 60-point discretisation was either higher-order or differently placed. The grid convergence is cleanly second order (C₀: 57.897 → 57.884 → 57.881). The discrete C₀(60) = 57.897 is what the model equation will use (D4), not 57.9.

### Tier 1 (2026-10-18): passed, 42/42 tests

N = 60, Δτ = 0.01. Production stepper: exponential integrator with the quenching coefficient at a midpoint-predicted a (`DECISIONS.md` D3, correction found by T1.4). Cross-check: IMEX Euler, held to 5% instead of 1% because it is first order. "Discrete" is the value from the 60-point operator (`ModelConstants`), the exact target for the integrator; "paper" is the printed value.

| Gate | Quantity | Paper | Discrete | Exponential | IMEX | Result |
|---|---|---|---|---|---|---|
| T1.1 | \|a − 1\| after 100 τ_d from b₀ | 0 | 0 | 1×10⁻¹⁴ | 1×10⁻¹⁴ | pass (field drift 2×10⁻¹²) |
| T1.2 | Λ, kR = 0.5, C_m = 100 | 0.527 | 0.52834 | 0.52856 | 0.52125 | pass; IMEX 1.1% off the paper |
| T1.2 | Λ, kR = 0.3, C_m = 120 | 0.815 | 0.81635 | 0.81690 | 0.80997 | pass |
| T1.3 | saturation from a random field, 50 τ_d | \|a\| = 1 | — | < 10⁻⁶ | < 10⁻⁶ | pass, overtones < 10⁻⁶ |
| T1.4 | 2μ, kR = 0.5 | 1.526 *(derived)* | 1.52586 | 1.52580 | 1.4727 | pass; not 2Λ = 1.057 |
| T1.4 | 2μ, kR = 0.3 | 2.97 *(derived)* | 2.96766 | 2.96717 | 2.8282 | pass; not 2Λ = 1.633 |
| T1.5 | b → −b | exact | — | bitwise | bitwise | pass |

T1.6 (diagnostic). Relative rate errors against the discrete values, Δτ = 0.01 → 0.001:

| Stepper | Λ, kR 0.5 | 2μ, kR 0.5 | Λ, kR 0.3 | 2μ, kR 0.3 |
|---|---|---|---|---|
| exponential | 4.3×10⁻⁴ → 4.3×10⁻⁶ (×100) | 3.7×10⁻⁵ → 2.3×10⁻⁶ | 6.8×10⁻⁴ → 6.8×10⁻⁶ (×100) | 1.6×10⁻⁴ → 5.2×10⁻⁵ |
| IMEX | 1.3% → 0.13% (×10) | 3.5% → 0.36% (×10) | 0.78% → 0.08% (×10) | 4.7% → 0.49% (×10) |

The exponential stepper is second order in the growth rate; its relaxation error at Δτ = 0.001 sits on the measurement floor (round-off relative to the decayed 10⁻⁶ perturbation), so only its monotone decrease is asserted. IMEX is first order throughout; its 2μ error is large because b̃₀ is not a left eigenvector of the IMEX amplification matrix, so overtones leak into a at O(Δτ). The IMEX amplification matrix at C₀ has spectral radius 1 − 10⁻¹⁵ at Δτ = 0.01 for both kR: stable. Cost (BenchmarkTools, Apple silicon): exponential 2.0 µs per step, IMEX 1.4 µs; both `step!` methods pass the D14 type-stability and zero-allocation tests.

Observation for `QUESTIONS.md` Q5: the deterministic relaxation about a = 1 is 2μ to 10⁻⁴, and 2Λ is 31–45% lower, so the two rates are cleanly distinguishable. Whether the noisy autocorrelation follows 2Λ (HSO Eq. 15) or 2μ is for T5.13.

### Tier 2 (2026-10-19): passed, 362/362 tests

N = 60, Δτ = 0.01, N_c = 10, τ_c = 0.05 (n_c = 5 steps), f = 6.4, Gaussian F unless stated. Bins by gap G3: node j in bin ⌊θ_j N_c/π⌋, which puts 5 nodes in each polar bin and 6 in the others, mirror-symmetric about the equator. `src/noise.jl`: `BinnedNoise`, `GaussianF`/`UniformF`/`BinaryF`, `ConstantForcing`/`LinearRamp`, `run_seed`, `noise_projection` (V₀₀), `binned_D0`.

| Gate | Check | Target | Measured | Result |
|---|---|---|---|---|
| T2.1 | mean and variance of F, 10⁶ draws, each of 3 distributions | 0 ± 4σ, 1 ± 1% | Gaussian −0.4σ, 1.0005; uniform +0.9σ, 0.9996; binary −3.4σ, 1 | pass |
| T2.1 | cross-bin (45 pairs) and successive-interval (10 bins) correlations | 0 ± 4.5σ | all within 4.5σ, σ = 3.2×10⁻³ | pass |
| T2.1 | F constant for n_c = 5 steps, redrawn at n mod 5 = 0, shared by all nodes of a bin | exact | bitwise | pass |
| T2.2 | zero-lag variance ⟨v_j²⟩ = f²/(2N_c sin θ_j), all 58 nodes, 2×10⁵ independent samples | within 2% | max deviation < 2% | pass |
| T2.2 | angular structure: normalised covariance of nodes d apart (d = 0 … 12), averaged over j | same-bin fraction (exact) ± 0.01; triangle 1 − \|ζ\|/θ_c ± 0.03 | within both | pass |
| T2.2 | time structure: autocorrelation at lag l steps, nodes 5 and 29 | 1 − l/n_c (l ≤ n_c), 0 beyond, ± 0.01 | within 0.01 for l = 0 … 10 | pass |
| T2.3 | ℱ = πf²τ_c/(2N_c²) | 3.22×10⁻² | 3.217×10⁻² | pass |
| T2.3 | D₀ from Eq. 36 (δ-correlated B.7) with the 60-point integral 3.5324 | 0.114 | 0.11363 | pass |
| T2.3 | D₀ measured: Var[∫₀ᵀ V₀₀ dτ]/T, 10⁴ segments of T = 10 τ_d, V₀₀ on the frozen b₀ | 0.114 within ~10% | 0.10875 (−4.3%; statistical error 1.4%) | pass (D) |
| T2.3 | exact D₀ of the binned noise, τ_c Var[V₀₀] | — | 0.10778 at N = 60, 0.10632 at N = 240 | recorded |
| T2.4 | (D₀τ_c)^{1/2} from Eq. 36 | 0.075 | 0.0754 | pass (1%) |
| T2.4 | (D₀τ_c)^{1/2} from the measured D₀ | 0.075 within 10% | 0.0737 | pass |
| T2.5 | same seed → same a(τ), 10⁴ steps with noise; different seed → different | bitwise | bitwise | pass |
| T2.5 | one noise series across f: F stream independent of f, v(7.4)/v(5.4) = 7.4/5.4 | bitwise / 10⁻¹⁴ | bitwise / 10⁻¹⁴ | pass |
| T2.5 | `run_seed(master, i)`: deterministic, 1000 distinct seeds, distinct across masters | — | as stated | pass |
| T2.5 | `LinearRamp` evaluated at each redraw, constant within τ_c | exact | exact | pass |
| T2.5 | record and stream modes agree; replaying a(τ) through the aggregators reproduces the online values | bitwise | — | deferred to build step 6 (`run.jl`) |

Observations.

- The δ-correlated formula (Eq. 36) overestimates the D₀ of the binned noise by 5% at N = 60 and 6% at the grid-converged N = 240. For piecewise-constant F the exact value is τ_c f² Σ_k (∫_k w P̃₀T₀/√(2N_c sin θ))², and by Jensen's inequality each squared bin integral is at most θ_c times the integral of the square, with equality only if P̃₀T₀/√sin θ were constant across the bin. This is the "finite bins vs B.7" effect the gate anticipated; the measured value agrees with the exact binned value to 1% (within its 1.4% statistical error). The theory module (Tier 3) will use Eq. 36 as the papers do; the 5% is a known systematic when comparing simulation with theory.
- The 60-point bins hold 5 or 6 nodes, so the angular correlation is a staircase approximation to the triangle of Eq. B.6 (largest deviation 0.03 at d = 1–3 nodes); the time correlation is the exact triangle because every step samples the same phase pattern.
- T2.1's mean criterion was widened from 3σ to 4σ: at the fixed test seed the binary F has a −3.4σ mean over 10⁶ draws. A direct check of Julia's `rand(Xoshiro, Bool)` with 20 seeds × 5×10⁷ draws gives a pooled z of 0.13 and an rms z of 0.81, so the generator is unbiased and the event is chance (probability 7×10⁻⁴ per test).
- D14: `step!` with `BinnedNoise`, `update_noise!` (redraw and quiet steps, constant and ramped forcing), `redraw!` for all three distributions and `noise_projection` are type-stable and allocation-free. The test harness's allocation helper had to be specialised on the function: a dynamically dispatched call boxes a `Float64` return value and charged 16 bytes to `noise_projection`. Cost: the step with noise is 2.04 µs (exponential) and 1.42 µs (IMEX), unchanged from the noise-free step; a redraw step of `update_noise!` costs 56 ns, a quiet step 3 ns, the V₀₀ projection 12 ns (BenchmarkTools, Apple silicon).

### Aggregators (2026-10-19): 150 unit tests, not a tier

`src/aggregators.jl` (`DECISIONS.md` D8) is tested on synthetic series against offline reference implementations written in the test file: Welford moments against direct sums (10⁻¹²; variance preserved under a 10⁸ offset), fixed-bin histograms, the autocorrelation estimator at lags 0 … L against the direct lagged sums (10⁻¹⁰, including a series longer than one 2²⁰-sample accumulation chunk) and against ρ^l for an AR(1) series, cadence sampling, the moving average, the Schmitt-trigger detector on hand-built series (reversals, aborted reversals, excursions, dips, zero-duration band jumps, crossings on entering the band) and on a random walk (reversal count equals the offline polarity flips), the papers' crossing-grouping rule, the exclusion-window variance against an offline implementation of the §2 definition (10⁻¹⁰), binned means with the weighted even-quadratic fit (recovers 1.36 + 2.00a² to 0.01), `merge!` of every aggregator against pooling, and Tuple composition. The random-walk test found two detector defects before any simulation data existed: a jump from +a_h to −a_h within one sample (possible at τ_d sampling) was not a reversal, and a crossing on the sample that enters the band was not counted. Every `fit!` and the fixed-size `merge!` methods are type-stable and allocation-free (D14 tests). Cost per sample (BenchmarkTools, Apple silicon): `Moments` 12 ns, `Histogram` 3.5 ns, `SchmittTrigger` 3.7 ns, `CrossingGroups` 3.5 ns, `ExclusionVariance` 12 ns, `BinnedMean` 7.6 ns, `Sampled` 2.3 ns, `AutoCorrelation` with 500 lags 322 ns (fed at τ_d cadence it is negligible); all seven scalar aggregators in a Tuple 324 ns, i.e. the every-step statistics of a run add about 1% to the 2.0 µs step.

### Tier 3 (2026-10-20): passed, 301/301 tests

`src/theory.jl` (`DECISIONS.md` D12): `BistableTheory(γ, c²; Λ, μ)` with κ = μ/Λ − 1; the drift S = Λ(1 − a²)a(1 + κa²), D = Δ₀(a² + c²), the potential, ψ = ln p and its derivatives; the closed forms (a_m, ψ_m, ψ_m'', Eq. 56, ⟨T_r⟩ as HSO Eq. A.2, Eqs. 59–60) written so that one code path reduces to HOS Eqs. 52–58 at κ = 0 and remains valid for HSO's A.1; the quadratures (normalisation, the moments over a > 0 for dv, the escape integral Eq. 57) with every exponential relative to max ψ. Λ = 0.527 and κ = 0 unless stated.

| Gate | Check | Target | Computed | Result |
|---|---|---|---|---|
| T3.1 | Eq. 43 (and A.1 at κ = 0.45) vs Eq. 42 integrated numerically, a ∈ [−2, 2] | 10⁻⁸ | < 10⁻⁸ after normalisation; ∫p = 1 to 10⁻⁸ | pass |
| T3.1 | log-space quadrature at γ = 1500, where exp ψ_m overflows | finite | ln Z, ln⟨T_r⟩, dv finite; ψ(∞) = −∞, no NaN | pass |
| T3.2 | a_m, ψ_m, ψ_m'' vs bisection on a finite-difference ψ' and Richardson derivatives, five (γ, c², κ) sets | 10⁻⁸ | < 10⁻⁸; the general forms equal Eqs. 52–54 to 10⁻¹² at κ = 0 | pass |
| T3.3 | Eq. 56, Table II runs 1–4 | 0.04, 0.08, 0.17, 0.47 | 0.042, 0.083, 0.175, 0.467 | pass |
| T3.3 | Eq. 56, Table II run 5 (γ = 1.09) | 1.0 | 10.45 | recorded, see below |
| T3.4 | Eq. 58, Table II runs 1–5 | 12000, 480, 94, 35, 35 | 11843, 483.0, 93.9, 35.4, 34.8 | pass (run 1 is 1.3% off, within its last printed digit) |
| T3.5 | Eq. 58, HSO Table 1 cases 1–3 | 1090, 928, 665 | 1088.7, 927.1, 664.6 | pass (0.1%) |
| T3.6 | dv (HSO Eq. 12) from Eq. 10 over a > 0, cases 1–3 | 0.076, 0.070, 0.079 | 0.0758, 0.0696, 0.0792 | pass |
| T3.7 | A.1 with the generalised A.2, case 4 (Λ = 0.527, μ/Λ = 1.45) | 1097; dv 0.071 | 1109 (1.1%); 0.0705 | pass at 2%, see below |
| T3.7 | case 5 (Λ = 0.815, μ/Λ = 1.82) | 138; dv 0.091 | 137.6; 0.0910 | pass |
| T3.7 | model constants at N = 240: (Λ, μ/Λ) for kR = 0.5, C_m = 100 and kR = 0.3, C_m = 120 | (0.527, 1.45); (0.815, 1.82) | (0.5275, 1.446); (0.8148, 1.821) | pass (1%) |
| T3.8 | Eq. 57 at a′ = 1 (SOH's 4∫ form) over Eq. 58 at (γ, c²) = (3, 0.2), (6, 0.3) | ratio | 1.180, 1.123; a′ = 0.8 → ×0.990, ×0.998; a′ = 1.2 → ×1.005, ×1.001 | recorded (D) |
| T3.9 | ln(Eq. 59 / Eq. 58) at c² = 0.25, γ = 2, 5, 10, 20, 50, 100, 10³, 10⁴, 10⁵ | → 0 monotonically | −0.325, −0.109, −0.052, −0.026, −0.010, −0.0050, −5.0×10⁻⁴, −5.0×10⁻⁵, −5.0×10⁻⁶ | pass (Eq. 59 = Eq. 58 × (1 − 0.5/γ) at this c²) |
| Eq. 63 | (c², γ) from the D₁ fit 1.36 + 2.00a², D0int = 3.529, ℱ = 3.22×10⁻² | 0.246, 2.96 | 0.2460, 2.960 | pass |

Observations.

- Table II run 5 (f = 8.4, γ = 1.09): Eq. 56 gives (γc² + γ − 1)/[4(γ − 1)²] = 10.4, not the printed 1.0; the tracked text layer confirms the print and no reading of the formula gives 1.0. Footnote b of the table says the run is outside the validity of the parabolic approximation because γ is close to 1. For these parameters dv on the exact p(a) over a > 0 is 0.48, close to the run's 0.50. The gate asserts runs 1–4 and records run 5.
- HSO case 4: the generalised A.2 gives 1109 against the printed 1097 (1.1%). ⟨T_r⟩ is exponentially sensitive to the shape parameters: with the printed two figures, γ = 4.8 ± 0.05 spans 1046–1176 and c² = 0.22 ± 0.005 spans 1041–1183, so a 1% gate is below the resolution of the inputs. Cases 1–3 agree to 0.1%, so the formula is the same and HSO evidently evaluated those from the printed values; case 4 presumably from unrounded fit output. The test asserts 2% and that the target lies inside the γ ± 0.05 bracket. Case 5 needs Λ = 0.815, the growth rate of the kR = 0.3, C_m = 120 run of HSO Fig. 8 (with Λ = 0.527 it would be 213); our eigensolver gives Λ = 0.8148 and μ/Λ = 1.821 for that model, matching HSO's 0.815 and 1.82, an extra check of the eigenproblem at a second parameter set.
- The exact escape integral (Eq. 57) exceeds the Gaussian-approximated Eq. 58 by 18% at (3, 0.2) and 12% at (6, 0.3), and by 8%, 14%, 17% and 2% for Table II runs 1–4 (run 5, with barely a well, gives 46% of Eq. 58). This quantifies the F.2 approximation and moves the theory toward the measured values of runs 2 and 3 (552 and 110 against the runs' 519 ± 53 and 107 ± 5); it is also the sense in which HSO (p. 151) call A.2 "only approximate". The dependence on a′ is weak, as HOS p. 289 states (1% between a′ = 0.8 and 1.2).
- dv on the exact p(a) over a > 0 (HSO's definition) exceeds Eq. 56 for the first three Table II runs (0.046, 0.099, 0.196 against 0.042, 0.083, 0.175) and falls below it for run 4 (0.330 against 0.467) because it includes the reversal region near a = 0 and uses the full distribution rather than the parabola; it is not the between-reversal variance of §2. Both will be compared with the run statistics at Tier 4.
- D14: the closed forms are type-stable and allocation-free. The quadratures allocate (QuadGK's segment heap) and are checked for inference only. The nested escape integral (Eq. 57, a `quadgk` inside a `quadgk`) cannot be made JET-clean: Julia's inference recursion limit fires when QuadGK's internals are re-entered, whether the integrands are closures or callable structs and whether the inner interval is finite or not; its return type is asserted in the source, the two functions are checked with `@inferred` only, and the cost (milliseconds, a few calls per analysis) is irrelevant.

