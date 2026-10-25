# Validation gates

Accepted by the user on 2026-10-18: the gates, the tolerance rules and the operational definitions below. Targets come from the three papers via `resources/`; notation is HOS's unless stated; *(derived here)* marks our own derivation rather than a statement in a paper. Gate IDs (T0.1, …) are the ones referenced from `QUESTIONS.md` and `DECISIONS.md`.

Status: Tiers 0–3 passed (Tier 0 and Tier 1 on 2026-10-18, Tier 2 on 2026-10-19, Tier 3 on 2026-10-20); Tier 4's standard run and forcing sweep passed on 2026-10-20 (T4.1, T4.2, T4.5, T4.6, T4.11, T4.18; T4.12 recorded), with the D₁-fit intercept and c² of T4.3, T4.4, T4.7 outside the 10% gate at f ≤ 6.4 and recorded as `@test_broken` (`DISCREPANCIES.md` DR-5, attributed to the papers' unstated integrator); T4.17 (convergence), the HSO runs T4.8, T4.9 and the SOH runs T4.10, T4.15 passed the same day; T4.14 passed and T4.13 was recorded on 2026-10-20 (the median parity outside events is −0.82, not below −0.9: `DISCREPANCIES.md` DR-7); the stretch goals were done on 2026-10-24: SG-1 excludes the time integrator as the cause of DR-5, and SG-2 recorded T4.16 as not met (the alternative quenching laws change the statistics or run away, `DISCREPANCIES.md` DR-8). Every gate of Tiers 0–5 is now passed or recorded. Tier 5's T5.1, T5.2, T5.7, T5.8 were judged comparable by the user on 2026-10-20; T5.3, T5.4 (HOS Figs. 3, 4) were rendered on 2026-10-20 and judged satisfactory by the user on 2026-10-24, as were T5.5, T5.6, T5.9, T5.10, T5.11, T5.12, T5.13 (HOS Figs. 5, 6; HSO Figs. 3, 6, 7, 8; SOH Fig. 3) the same day: Tier 5 complete. Results in §4.

## 1. Conventions

- **Tiers** are ordered: pass a tier before spending effort on the next.
- **Kind.** A = automated numeric check with a pass/fail criterion. D = diagnostic: computed and recorded, no pass/fail, feeds `QUESTIONS.md`. J = qualitative judgement by the user, against the paper's figure.
- **Tolerances.** For deterministic quantities printed to n significant figures, pass if within ±1 in the last printed digit, relaxed to 1% because the papers' own discretisation is unknown; record both the 60-point and a grid-converged value. For run statistics, pass if |ours − theirs| < 2·√(σ_theirs² + σ_ours²), where σ_ours is from repeated seeds. For convention-dependent statistics (reversal counts, between-reversal variance) the convention used must be stated with the result.
- **Hard vs soft.** Tiers 0, 1 and 3 are hard gates: a failure blocks progress. Tiers 2, 4 are hard for A items and informative for D items. Tier 5 passes on the user's say-so. Tier 6 is deferred until the Sint-800 data are obtained.

## 2. Operational definitions

These close gaps G2–G5 and G9 of `DEVELOPMENT.md` §2. Every reported statistic quotes the definition it used.

- **Reversal detection.** Schmitt trigger on a(τ): polarity is +1 until a < −a_h, then −1 until a > +a_h. Default a_h = 0.5. A polarity flip is a reversal; a zero crossing without a flip is an aborted reversal; a dip of |a| below a_e = 0.3 without a zero crossing is an excursion. Report sensitivity to a_h ∈ {0.3, 0.5, 0.7}. Also implement the papers' rule (group zero crossings within W = 5 τ_d; odd count = reversal, even = aborted) for comparison.
- **Between-reversal variance (HOS/SOH style).** Within each polarity interval, drop W τ_d after the preceding reversal and before the next (default W = 5), compute ⟨(a−⟨a⟩)²⟩/⟨a⟩² over the remainder, pool across intervals time-weighted. Sensitivity W ∈ {2, 5, 10}. *Note (step 7):* the pooling is of the moments (the variance about each interval's mean and |⟨a⟩| pooled by kept samples, then the ratio), not of the per-interval ratios, which a few near-zero-mean intervals dominate at high forcing; and the paper's "variance from run" column is reproduced at every forcing by dv below, not by this estimator at f ≥ 7.4 (`DISCREPANCIES.md` DR-6).
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
| T4.16 | Robustness regression: alternative quenching (q from total toroidal energy; α/(1+q); local quenching) | dipole statistics "not materially changed" | ⟨T_r⟩, variance within ~20% | HOS pp. 277, 281 | D, stretch goal (`STRETCH_GOALS.md` SG-2); recorded 2026-10-24, not met (§4) |
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
| T2.5 | record and stream modes agree; replaying a(τ) through the aggregators reproduces the online values | bitwise | bitwise (2026-10-20, build step 6: 600 τ_d at seed 3 with 4 reversals; state, every standard aggregator, and the replay of the recorded a and D₁/ℱ streams through fresh detectors, histogram and binned means; a checkpointed run stopped at 450 τ_d and resumed from its 400 τ_d checkpoint is also bitwise identical, streams included) | pass |

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

- (Logged in `DISCREPANCIES.md`.) Table II run 5 (f = 8.4, γ = 1.09): Eq. 56 gives (γc² + γ − 1)/[4(γ − 1)²] = 10.4, not the printed 1.0; the tracked text layer confirms the print and no reading of the formula gives 1.0. Footnote b of the table says the run is outside the validity of the parabolic approximation because γ is close to 1. For these parameters dv on the exact p(a) over a > 0 is 0.48, close to the run's 0.50. The gate asserts runs 1–4 and records run 5.
- (Logged in `DISCREPANCIES.md`.) HSO case 4: the generalised A.2 gives 1109 against the printed 1097 (1.1%). ⟨T_r⟩ is exponentially sensitive to the shape parameters: with the printed two figures, γ = 4.8 ± 0.05 spans 1046–1176 and c² = 0.22 ± 0.005 spans 1041–1183, so a 1% gate is below the resolution of the inputs. Cases 1–3 agree to 0.1%, so the formula is the same and HSO evidently evaluated those from the printed values; case 4 presumably from unrounded fit output. The test asserts 2% and that the target lies inside the γ ± 0.05 bracket. Case 5 needs Λ = 0.815, the growth rate of the kR = 0.3, C_m = 120 run of HSO Fig. 8 (with Λ = 0.527 it would be 213); our eigensolver gives Λ = 0.8148 and μ/Λ = 1.821 for that model, matching HSO's 0.815 and 1.82, an extra check of the eigenproblem at a second parameter set.
- The exact escape integral (Eq. 57) exceeds the Gaussian-approximated Eq. 58 by 18% at (3, 0.2) and 12% at (6, 0.3), and by 8%, 14%, 17% and 2% for Table II runs 1–4 (run 5, with barely a well, gives 46% of Eq. 58). This quantifies the F.2 approximation and moves the theory toward the measured values of runs 2 and 3 (552 and 110 against the runs' 519 ± 53 and 107 ± 5); it is also the sense in which HSO (p. 151) call A.2 "only approximate". The dependence on a′ is weak, as HOS p. 289 states (1% between a′ = 0.8 and 1.2).
- dv on the exact p(a) over a > 0 (HSO's definition) exceeds Eq. 56 for the first three Table II runs (0.046, 0.099, 0.196 against 0.042, 0.083, 0.175) and falls below it for run 4 (0.330 against 0.467) because it includes the reversal region near a = 0 and uses the full distribution rather than the parabola; it is not the between-reversal variance of §2. Both will be compared with the run statistics at Tier 4.
- D14: the closed forms are type-stable and allocation-free. The quadratures allocate (QuadGK's segment heap) and are checked for inference only. The nested escape integral (Eq. 57, a `quadgk` inside a `quadgk`) was at first not JET-clean: Julia's inference recursion limit fires when QuadGK's adaptive routines are re-entered, whether the integrands are closures or callable structs and whether the inner interval is finite or not, and the inner call degraded to dynamic dispatch. Evaluating the inner integrand through QuadGK's `BatchIntegrand` sends it down a different method chain, so nothing is re-entered: fully inferred, identical result to the last bit, 0.3 ms instead of 1.8 ms per ⟨T_r⟩. The nesting itself is intrinsic (a two-dimensional integral over a triangle); no change of variables removes it.

### Tier 4, standard run (2026-10-20): T4.1, T4.2 passed; T4.3 and T4.4 pass on γ and c₂ but fail on the intercept and c² (`@test_broken`, DR-5)

Three seeds × 5×10⁴ τ_d in stream mode (14 s each), aggregators pooled with `merge!`; the definitions of §2 (Schmitt trigger a_h = 0.5, a_e = 0.3; exclusion variance W = 5 τ_d; dv on |a| at 1 τ_d; D₁/ℱ bin averages of width 0.05 fitted for |a| ≤ 1.1 weighted by the bin counts, HOS Fig. 8 caption). `test/tier4_stochastic.jl`, `src/run.jl` (`standard_observers`, `summarize`, `ensemble`). *Corrected at step 7:* the step-6 version of this section took the unweighted fit as the paper's procedure and reported T4.3–T4.4 as passed; the caption says the fit is count-weighted, and with it the intercept and c² miss the 10% gate (`DISCREPANCIES.md` DR-5).

| Gate | Check | Target | Measured (pooled; per seed) | Result |
|---|---|---|---|---|
| T4.1 | ⟨T_r⟩, Schmitt trigger | 107 ± 5 (468 reversals) | 104.7 ± 2.8 from 1433 reversals; 109.6 ± 5.1, 104.6 ± 4.8, 100.2 ± 4.5 | pass (2σ combined = 11.5) |
| T4.1 | the papers' grouping rule (W = 5 τ_d) | — | 432, 456, 479 reversals against the trigger's 456, 478, 499; trigger: 253, 229, 235 aborted, 2295, 2176, 2197 excursions | recorded |
| T4.2 | between-reversal relative variance, §2 estimator (pooled moments, DR-6) | 0.16 ± 0.02 | 0.142 (step 6 reported 0.164 from the pooled per-interval ratios) | pass |
| T4.2 | dv, HSO Eq. 12 on \|a\| | 0.16 ± 0.02 | 0.160 | pass |
| T4.3 | D₁/ℱ = c₀ + c₂a², \|a\| ≤ 1.1, fit of the bin averages weighted by the counts (the paper's procedure) | 1.36 + 2.00a², 10% | 1.584 + 1.91a² (16%, 4.5%) | c₀ **fails** (`@test_broken`), c₂ passes |
| T4.3 | the same fit unweighted; mean of the bins \|a\| ≤ 0.2 | — | 1.464 + 2.08a² (7.6%, 4%); 1.24 | recorded, see below |
| T4.4 | (c², γ) from the weighted fit with D0int = 3.532 (N = 60), ℱ = 3.217×10⁻², Λ = 0.5283 | (0.246, 2.96), 10% | (0.291, 3.02): 18%, 2% | c² **fails** (`@test_broken`), γ passes; unweighted (0.261, 2.93) |
| T4.4 | Eqs. 56 and 58 at the measured (c², γ) | paper: 0.17, 94 at (0.246, 2.96) | 0.181, 83.9 (Eq. 57: 97.7) against the run's 0.164 and 104.7 | recorded |
| T4.13 | parity 𝒫 every n_c steps over the whole run, events included | median < −0.9 outside event windows (D) | median −0.815; 31% of samples below −0.9 | recorded; outside the event windows −0.821, see the event-level section (DR-7) |

Observations.

- Reversal rate and variance: 104.7 ± 2.8 against 107 ± 5, and 0.142 (§2 estimator) or 0.160 (dv) against 0.16, so the standard run reproduces HOS Table II row 3 within the errors; the seed-to-seed scatter of ⟨T_r⟩ (100–110) is what the Poisson error (≈ 5) predicts. The grouping rule finds 4–5% fewer reversals than the trigger, because a pair of crossings within 5 τ_d counts as one aborted reversal.
- (Logged in `DISCREPANCIES.md` DR-5.) The D₁ intercept depends on the weighting of the fit: the bin averages are only approximately quadratic (D₁ grows faster than a² beyond |a| ≈ 0.8, the "a⁴ term" HOS see beyond 1.1), so the paper's count-weighted fit, which favours the populous bins near |a| ≈ 0.9, gives 1.58, an unweighted fit through the bin means 1.46, and the direct mean of the bins |a| ≤ 0.2 gives 1.24 ± 0.1 (1300 samples of a heavy-tailed quantity). The paper's 1.36 lies between the direct estimate and the unweighted fit; its histogram "rises from ≈ 1.4 at a = 0". c₂ is insensitive (1.9–2.1). `summarize` reports all three; the gate uses the paper's weighted fit, and the 16% excess of its intercept, equivalently a 4–6% excess of D₁ where the samples are, is the open item of DR-5.
- Theory at the measured (c², γ): Eq. 58 gives 84 against the measured 105 (the exact Eq. 57: 98); Eq. 56 gives 0.18 against 0.164. The same pattern as HOS's own comparison (94 against 107, 0.17 against 0.16), shifted by our slightly larger c².
- Parity: between reversals 𝒫 sits at −0.8 to −1 with brief spikes toward 0, and reaches +1 at the reversals (Fig. 2); the whole-run median is −0.82, and masking the event windows leaves it at −0.821 (the event-level section below, DR-7).
- Cost: 14 s per 5×10⁴ τ_d in stream mode (2.3 µs per step with all probes) and 14 s in record mode with every stream and the T, P snapshots (130 MB per run).

### Tier 4, forcing sweep (2026-10-20): T4.5, T4.6, T4.11, T4.18 passed; T4.7 passes on γ everywhere and on c² at f ≥ 7.4 (`@test_broken` at f ≤ 6.4, DR-5); T4.12 recorded

`ensemble` runs the seeds in parallel over the Julia threads and pools the aggregators (`src/run.jl`); `scripts/run_sweep.jl` writes `results/table2.md`; `test/tier4_stochastic.jl`. Three seeds per forcing, 5×10⁴ τ_d each as in the paper, 2×10⁵ τ_d each at f = 4.4 (the paper had ~5 reversals there). About 2.5 min on 4 threads. The D₁/ℱ fit is the paper's: bin averages for |a| ≤ 1.1 weighted by the bin counts; the paper's fit is the one its printed (c², γ) imply through Eqs. 62–63, which for f = 6.4 is the printed 1.36 + 2.00a².

| f | reversals | ⟨T_r⟩ (paper) | §2 variance (paper) | dv | D₁/ℱ weighted fit (paper) | c² (paper) | γ (paper) | KS p (n) |
|---|---|---|---|---|---|---|---|---|
| 4.4 | 53 | 11321 ± 1555 (10700 ± 5300) | 0.030 (0.03) | 0.031 | 1.57 − 0.01a² (1.22 + 0.31a²) | 0.447 (0.318) | 9.86 (9.04) | 0.73 (50) |
| 5.4 | 285 | 526 ± 31 (519 ± 53) | 0.068 (0.07) | 0.073 | 1.64 + 0.75a² (1.34 + 0.94a²) | 0.382 (0.299) | 5.39 (5.15) | 0.87 (282) |
| 6.4 | 1433 | 104.7 ± 2.8 (107 ± 5) | 0.142 (0.16) | 0.160 | 1.58 + 1.91a² (1.36 + 2.00a²) | 0.291 (0.246) | 3.02 (2.96) | 0.77 (1398) |
| 7.4 | 4051 | 37.0 ± 0.6 (37.6 ± 1) | 0.263 (0.30) | 0.301 | 1.62 + 3.37a² (1.62 + 3.20a²) | 0.234 (0.240) | 1.78 (1.82) | 0.37 (3707) |
| 8.4 | 7789 | 19.3 ± 0.2 (18.7 ± 0.4) | 0.457 (0.50) | 0.496 | 1.92 + 5.25a² (1.99 + 5.20a²) | 0.218 (0.228) | 1.09 (1.09) | 0.11 (6185) |

Unweighted fits: 1.41 + 0.27a², 1.39 + 1.10a², 1.46 + 2.08a², 1.58 + 3.41a², 2.04 + 4.99a², giving c² = 0.372, 0.300, 0.261, 0.228, 0.239 and γ = 9.14, 4.98, 2.93, 1.77, 1.12. Direct D₁(0) from the bins |a| ≤ 0.2: 1.04, 1.06, 1.24, 1.40, 1.50. Theory at our (c², γ): Eq. 56 0.042, 0.084, 0.178, 0.492, 10.9; Eq. 58 4909, 321, 78.6, 34.9, 35.4; Eq. 57 (exact escape integral) 5329, 369, 91.2, 35.3, 16.3.

| Gate | Check | Result |
|---|---|---|
| T4.5 | ⟨T_r⟩ at the four other forcings, 2σ combined | pass at all four (0.1σ, 0.1σ, 0.5σ, 1.3σ); the rate reproduces over three decades |
| T4.6 | variance, ±0.02 or 15% | pass: dv within 0.005 of the column everywhere; the §2 estimator within 0.02 at f ≤ 6.4 and 11–13% low at f ≥ 7.4 (DR-6) |
| T4.7 | (c², γ) from the D₁ fit, 10% | γ passes everywhere (9%, 5%, 2%, 2%, 0%); c² passes at f = 7.4, 8.4 (2.5%, 4%) and **fails** at f = 4.4, 5.4, 6.4 (41%, 28%, 18%, `@test_broken`); the unweighted c² passes at 5.4–8.4 (0.3%, 6%, 5%, 5%) and fails at 4.4 (17%) |
| T4.11 | polarity intervals above 5 τ_d exponential, KS with bootstrap p | pass at every forcing, p = 0.73, 0.87, 0.77, 0.37, 0.11 (the last at ⟨T_r⟩ = 19 τ_d, where T_r ≫ Λ⁻¹ ≈ 1.9 is weakest) |
| T4.18 | seed scatter at f = 5.4 covers HOS's 519 ± 53 and HSO's 570 ± 61 | pass: 516, 610, 472 |

Observations.

- The reversal rate is the strongest result of the sweep: ⟨T_r⟩ reproduces from 19 to 11 000 τ_d, a factor 600, with the pooled values inside the paper's errors everywhere. The exact escape integral (Eq. 57) at our (c², γ) tracks the measured rate better than Eq. 58 at f = 8.4 (16.3 against 19.3, Eq. 58 35.4) and worse at f = 4.4 (5329 against 11 321), where our c² is off.
- The D₁ fit (T4.7) is well conditioned only where the |a| distribution is broad. At f = 7.4 and 8.4 the weighted and unweighted fits agree with each other and with the paper's to 2–5%. At f ≤ 6.4 the samples concentrate near |a| ≈ 0.8–1.1 and the intercept is an extrapolation: at f = 4.4 fewer than 0.5% of the samples have |a| < 0.5, the bins below 0.5 hold 12–150 samples with means scattering from 0.05 to 3, and the weighted fit returns c₂ ≈ 0 (D₁ is flat near |a| = 1 there). The paper's fits at low forcing have the same conditioning problem, so its printed (c², γ) at f = 4.4 and 5.4 carry an uncertainty the 10% gate ignores; the excess of our weighted intercept at f = 5.4 and 6.4 (23%, 16%) over the paper's remains DR-5. γ, which Δ₀ = ℱ(c₂ + 3.529) dominates, is insensitive to all of this.
- With the sweep, the fit that reproduces the paper's (c², γ) best is the unweighted one at f = 5.4–8.4 (c² within 6%), and the count-weighted one at f ≥ 7.4.
- DR-6: the variance column of Table II is the whole-record relative variance of |a| (dv), which we reproduce to 0.005 at every forcing; the §2 exclusion estimator coincides with it at f ≤ 6.4 and falls 11–13% below it at f ≥ 7.4, where ⟨T_r⟩ is only 2–4 times the 2W = 10 τ_d exclusion. The estimator's first implementation pooled the per-interval ratios and was dominated by a few near-zero-mean intervals (0.54 and 4.3 at f = 7.4 and 8.4); it now pools the moments.
- The KS test of T4.11 is a parametric bootstrap (the mean is fitted), `src/stats.jl`; the Kolmogorov distribution would give larger p-values.
- Cost: 15 s per 5×10⁴ τ_d per seed, 54 s at 2×10⁵; the whole Tier 4 suite (standard run, sweep, invariance) runs in about 5 min on 4 threads.

### Tier 4, the HSO runs (2026-10-20): T4.8, T4.9 passed; T5.13's decay rate recorded

`test/tier4_stochastic.jl` (the f = 5.4 ensemble of the sweep is HSO's reference run, HOS run 2; the second run is a new ensemble), three seeds × 5×10⁴ τ_d each. The χ² fit of the stationary distribution to the |a| histogram follows HSO §4.3 (`density_fit`, `src/stats.jl`): about 100 equal-population bins of the 0.005-wide histogram at 1 τ_d, model counts integrated over each bin, normalisation fixed by the total, Nelder–Mead in (ln(γ − 1), ln c²), goodness of fit as σ = |χ² − ν|/√(2ν). Eq. 10 is HOS Eq. 43 (κ = 0); Eq. A.1 uses this model's κ = μ/Λ − 1.

| Gate | Check | Target (HSO) | Measured (pooled; per seed) | Result |
|---|---|---|---|---|
| T4.8 | reference run ⟨T_r⟩ | 570 ± 61 (88 reversals) | 526 ± 31 from 285 (516, 610, 472) | pass (0.6σ combined); HOS's own 519 ± 53 for the same run is the other side of ours (T4.18) |
| T4.8 | dv | 0.070 ± 0.005 | 0.0728 | pass |
| T4.8 | (γ, c²), Eq. 10 fitted to the \|a\| histogram | (8, 0.50), 10% | (7.88, 0.509) at 342σ over 109 bins; per seed (7.64, 0.47), (7.94, 0.515), (8.0, 0.532) | pass (1.5%, 2%) |
| T4.8 | (γ, c²), Eq. A.1 fitted to the histogram | HSO case 4 (4.8, 0.22) was fitted to Sint-800 and the run together | (4.81, 0.231) at 78σ | recorded: on top of case 4 |
| T4.9 | second run kR 0.3, C_m 120, f 7: ⟨T_r⟩ | 113 ± 5 | 114.9 ± 3.2 from 1306 (119, 116, 110) | pass |
| T4.9 | dv | 0.081 ± 0.005 | 0.079 | pass |
| T4.9 | model constants | Λ = 0.815, μ/Λ = 1.82 | 0.8163, 1.818 | pass (T0.10) |
| T4.9 | Eq. A.1 fit; A.2 at that fit | HSO case 5 (5, 0.50), fitted to both records; A.2 gives 138 | (5.77, 0.54) at 92σ; A.2 gives 201 against the measured 115 | recorded |
| T4.11 | Poisson intervals of the second run | p > 0.05 | p = 0.81 (1264 intervals > 5 τ_d) | pass |
| T5.13 | small-lag decay rate of the autocorrelation of \|a\|, f = 5.4, fitted over 0.1–0.6 τ_d | 2Λ = 1.057 (HSO Eq. 15) or 2μ = 1.526 (T1.4) | 0.87 per τ_d; acor(0.5 τ_d) = 0.64, acor(5 τ_d) = 0.04, no zero crossing within 5 τ_d | recorded: follows 2Λ rather than 2μ, 18% below 2Λ |

Observations.

- The histogram fits reproduce HSO's Table 1 to a few per cent wherever HSO fitted the simulation alone (case 2), and even where they fitted the simulation together with Sint-800 (cases 4 and 5) our simulation-only fits land within 0.2 in γ and 0.04 in c²: the shape of p(|a|) is reproduced in detail, including the preference for Eq. A.1 over Eq. 10 (78σ against 342σ, HSO's "many σ" of footnote 4; the σ values grow with the sample size, ours being three runs pooled).
- The histogram fit and the D₁ fit give different (γ, c²) for the same run (7.9 against 5.4 at f = 5.4; HSO's 8 against HOS's 5.15), as the papers found; the theory of Eq. 43 is only approximate (HOS Fig. 9), so each fit absorbs the error differently. Both are reported by `summarize`.
- Q5 of `QUESTIONS.md` is settled: the stochastic autocorrelation of |a| decays at 0.87 per τ_d at small lags, close to the 2Λ of HSO Eq. 15 and far from the deterministic 2μ, so HSO's method for τ_d, which measures Λ/τ_d, uses the right rate. The rate is 18% below 2Λ, so a τ_d inferred with exactly 2Λ would be biased low by that much; the paper's error bar on 11 ± 2 kyr covers it.
- A.2 overshoots the measured ⟨T_r⟩ of the second run by 1.75 (201 against 115), more than HSO's 138/113, because the fit to the simulation alone gives a larger γ than their compromise fit; at the reference run the overshoot at our A.1 fit is the papers' factor of ~2 (`DISCREPANCIES.md` DR-4).

### Tier 4, the SOH runs (2026-10-20): T4.10, T4.15 passed

`test/tier4_stochastic.jl`; `segmented` (`src/run.jl`) runs consecutive segments with fresh observers for the windowed statistics of the ramp; `scripts/run_ramp.jl` writes the record-mode ramp for SOH Fig. 3. Three seeds each; SOH's "relative variance of a between reversals" has no stated exclusion rule, so both estimators are checked.

| Gate | Check | SOH | Ours (pooled; per seed) | Result |
|---|---|---|---|---|
| T4.10 | f = 6: ⟨T_r⟩ | ≈ 200 | 185 ± 7 from 812 (198, 187, 171) | pass (20%) |
| T4.10 | f = 6: relative variance | 0.11 | dv 0.117; §2 estimator 0.106 | pass |
| T4.10 | f = 7: ⟨T_r⟩ | ≈ 50 | 52.9 ± 1.0 from 2837 (54, 53, 51) | pass |
| T4.10 | f = 7: relative variance | 0.24 | dv 0.239; §2 estimator 0.209 | pass |
| T4.15 | ramp 4 → 6.5 over 5×10⁴ τ_d: rms variability, first and last 5000 τ_d | 17% → 41% | √dv 16% → 39%; √(§2) 16% → 37%; the ten segments: 16, 17, 19, 21, 24, 26, 29, 32, 36, 39% | pass (20% on the rms values); the rise is monotonic and roughly linear in f, as the paper says |
| T4.15 | reversals per 5000 τ_d segment (pooled over three seeds) | chrons shorten drastically | 0, 2, 2, 5, 19, 27, 34, 74, 90, 149 | pass |
| T4.15 | first chron | ≈ 1.9×10⁴ τ_d in the paper's realisation [approximate read-off] | first reversal at τ = 9300, 16 000, 7200 for seeds 1–3 | recorded (kind J): the same order; a single realisation of a process whose mean time is ~10⁴ τ_d at f ≈ 4–4.5 |

Observations.

- The SOH variants are the HOS model at two more forcings and sit on the same curve as the Table II sweep: ⟨T_r⟩ falls from 185 to 53 between f = 6 and 7 while the variance doubles, and both lie within 20% of the paper's round numbers. The two variance estimators bracket the paper's value at f = 7 (0.209 and 0.239 against 0.24), the §2 estimator being the lower one for the reason of DR-6.
- The ramp reproduces SOH's picture quantitatively: the rms variability rises from 16% to 39% (17% → 41%), roughly linearly with f, and the reversal count per 5000 τ_d climbs from none to ~150, i.e. ⟨T_r⟩ from > 10⁴ to ≈ 34 τ_d. The forcing parameter f²τ_c/N_c² rises by (6.5/4)² = 2.64 over the run, SOH's "factor ~2.5".

### Tier 4, event-level diagnostics (2026-10-20): T4.14 passed; T4.13 recorded, the −0.9 criterion not met (DR-7); Q9 answered

Per seed, the §2 Schmitt trigger, a second trigger with a_h = 0.3 (a_e = 0.15) and the parity series at the n_c cadence (`Series` and `outside_events` in `src/aggregators.jl`), three seeds × 5×10⁴ τ_d on the trajectories of the standard ensemble (the §2 trigger reproduces its 1433 reversals exactly); `test/tier4_stochastic.jl`. An event window is the stretch with |a| < a_h plus the sample that ends it; every reversal, aborted reversal and excursion is masked.

| Gate | Check | Target | Measured (three seeds pooled) | Result |
|---|---|---|---|---|
| T4.13 | median 𝒫 outside the event windows (14% of the run masked) | < −0.9 (D) | −0.821; quartiles −0.924, −0.617; 10th percentile −0.967; 32% of the samples below −0.9. Windows padded by 1 τ_d (24% masked): −0.823. Whole run: −0.811 | **not met** (`@test_broken`); recorded, DR-7 |
| T4.14 | reversal duration, the \|a\| < a_h = 0.5 stretch (§2) | ~3 τ_d typical, fast events ~0.5 τ_d (HOS p. 277) | 1433 reversals: median 2.62 τ_d, quartiles 0.75–4.93, mean 3.44; 19% shorter than 0.5 τ_d, 43% shorter than 2 τ_d, 24% longer than 5 τ_d | pass (reported) |
| T4.14 | the same with a_h = 0.3 | 0.5–2 τ_d (HSO p. 151, definition unstated) | 1505 reversals: median 0.80 τ_d, quartiles 0.20–3.27; 42% below 0.5 τ_d, 65% below 2 τ_d | pass (reported) |

Observations.

- Parity: the exclusion of the events barely moves the median (−0.811 → −0.821) because the events occupy 14% of the run; between them the symmetric part of the toroidal field carries a steady E_S/E_A = (1 + 𝒫)/(1 − 𝒫) ≈ 0.10 of the energy, which is what the inter-event overtone level implies (median |a₁| = 0.21 outside events in seed 1, |a₂| = 0.14; with the unit-energy modes of `DECISIONS.md` D4, a₁ alone gives 𝒫 ≈ (2|a₁|² − 1)/(2|a₁|² + 1) ≈ −0.83). HOS Fig. 2 shows the same |a₁| ≈ 0.1–0.3 between events, so the paper's "𝒫 ≃ −1" is a loose reading of its −1…1 axis, and the −0.9 criterion was our own operationalisation. The field is within 5% of pure dipole energy (𝒫 < −0.9) a third of the time and within 2% (𝒫 < −0.967) a tenth. `DISCREPANCIES.md` DR-7.
- Durations: HOS's "~3 τ_d typical" is the median and mean of the §2 definition, and their "fast events ~0.5 τ_d" are the shortest fifth (19% below 0.5 τ_d; the fastest, τ = 9304.63–9304.69 in seed 1, takes a from −a_h to +a_h in six steps). HSO's "(0.5–2) τ_d" is what the narrower band gives (median 0.8 τ_d). The two quoted ranges differ by definition, not by model. The distribution is broad, a factor 6.5 between the quartiles, because 31% of the reversals have three or more zero crossings and the longest last 10–18 τ_d (seed 1: 316 single-crossing reversals, 96 with three, 29 with five, 15 with seven or more).
- Q9 (`QUESTIONS.md`), the duration across the forcing sweep with the §2 definition, three seeds each: the median falls from 4.1 τ_d at f = 4.4 to 0.76 τ_d at f = 8.4 while the mean stays at 3–4 τ_d, because fast single-crossing reversals become the majority as the forcing rises while the long multi-crossing ones persist.

| f | reversals | median (τ_d) | quartiles | mean | fraction < 0.5 τ_d | fraction ≥ 5 τ_d |
|---|---|---|---|---|---|---|
| 4.4 | 53 | 4.14 | 2.28–5.73 | 4.16 | 0% | 36% |
| 5.4 | 285 | 3.44 | 1.60–5.81 | 4.00 | 10% | 32% |
| 6.4 | 1433 | 2.62 | 0.75–4.93 | 3.44 | 19% | 24% |
| 7.4 | 4051 | 1.71 | 0.37–4.50 | 3.16 | 30% | 22% |
| 8.4 | 7789 | 0.76 | 0.19–3.73 | 2.84 | 42% | 18% |

- Cost: 15 s per seed in parallel; the stored parity series is 8 MB per seed.

### Tier 4, convergence T4.17 and the DR-5 sensitivity (2026-10-20): passed

`scripts/dr5_sensitivity.jl`; the two resolution variants are also a testset in `test/tier4_stochastic.jl`. Standard run, three seeds × 5×10⁴ τ_d each; the Δτ = 0.005 run sees the same noise realisation (the F stream is drawn per τ_c).

| Variant | ⟨T_r⟩ | §2 variance | dv | weighted D₁ fit | D₁ bin 0.95–1.0 | (c², γ) |
|---|---|---|---|---|---|---|
| baseline (GaussianF, exponential, Δτ 0.01, N 60) | 104.7 ± 2.8 | 0.142 | 0.160 | 1.58 + 1.91a² | 3.33 | (0.291, 3.02) |
| Δτ = 0.005 | 106.9 ± 2.9 | 0.142 | 0.160 | 1.57 + 1.89a² | 3.32 | (0.289, 3.03) |
| N = 120 | 108.5 ± 2.9 | 0.139 | 0.156 | 1.58 + 1.89a² | 3.32 | (0.292, 3.03) |
| UniformF | 106.5 ± 2.8 | 0.145 | 0.162 | 1.55 + 1.98a² | 3.37 | (0.281, 2.98) |
| BinaryF | 111.9 ± 3.1 | 0.145 | 0.162 | 1.56 + 1.98a² | 3.39 | (0.283, 2.98) |
| IMEX stepper, Δτ 0.01 | 98.3 ± 2.5 | 0.135 | 0.153 | 1.84 + 1.78a² | 3.43 | (0.345, 3.09) |

- T4.17 passes: halving Δτ and doubling N move ⟨T_r⟩ by 2% and 4% (under 1σ), the variances by ≤ 0.004 and the D₁ fit by 0.01. The production integrator is converged at the papers' resolution.
- The distribution of F does not matter (the three agree within errors), as HOS state for the "other details" of the cell model.
- The time integrator does: a first-order IMEX scheme at the same Δτ raises the D₁ intercept by 0.26 and c² by 0.05 while changing ⟨T_r⟩ by 6% and dv by 0.007. This is the size of the DR-5 gap and the most probable explanation of it: the paper's unstated scheme carries an O(|λ₁|Δτ) bias on the overtone content that the extrapolated intercept amplifies, while ⟨T_r⟩ and the variance, which the gates reproduce, are an order of magnitude less sensitive. See `DISCREPANCIES.md` DR-5.

### Stretch goal SG-1, integrator sensitivity (2026-10-24): done; the DR-5 integrator hypothesis is not supported

`scripts/sg1_integrators.jl` (`-t 4`, 35 min; the tables are also in `results/sg1_integrators.md`). Two further steppers in `src/stepper.jl` (`DECISIONS.md` D3 note): `ThetaStepper`, the θ-method on the whole linear operator G(C₀) + s(1 − a²)E, as backward Euler (θ = 1, first order) and Crank–Nicolson (θ = 1/2, second order), each with the noise term explicit or treated like the operator (`implicit_noise`; trapezoidal in the Crank–Nicolson case). Block elimination to one dense (N−2)×(N−2) LU per step, 18 µs against 2 µs for the exponential stepper; Tier 1 passes on both (backward Euler at the IMEX's 5%, Crank–Nicolson at 1%), the block solve is checked against the dense θ-method formula to 10⁻¹² (`test/stepper.jl`). Standard run at f = 6.4, 5.4 and 7.4, six steppers, three seeds × 5×10⁴ τ_d pooled; the Crank–Nicolson run with trapezoidal noise is now part of the T4.17 testset as a second-order control from another integrator family.

**Analytic expectation.** The per-step factor each scheme applies to the leading modes, as (damping rate, frequency) per τ_d from ln|g|/Δτ and |arg g|/Δτ, against the exact (ℜλ, ℑλ):

| operator | mode | exact | exponential | IMEX | backward Euler | Crank–Nicolson |
|---|---|---|---|---|---|---|
| C₀ (a = ±1) | first overtone | (−2.757, 4.409) | (−2.757, 4.409) | (−2.567, 4.332) | (−2.811, 4.288) | (−2.756, 4.409) |
| C_m (a = 0) | fundamental | (0.528, 0) | (0.529, 0) | (0.521, 0) | (0.530, 0) | (0.528, 0) |
| C_m (a = 0) | first overtone | (−2.159, 5.875) | (−2.158, 5.875) | (−1.899, 5.724) | (−2.301, 5.745) | (−2.157, 5.874) |
| C_m (a = 0) | second overtone | (−8.201, 6.024) | (−8.201, 6.024) | (−7.644, 5.716) | (−8.037, 5.562) | (−8.198, 6.033) |

The IMEX under-damps the first overtone by 12% and backward Euler over-damps it by 7%, both lowering its frequency by 2–3%; Crank–Nicolson is exact to 0.1%, the exponential stepper to 0.02% (the splitting of the sE term).

**Runs.** Weighted and unweighted quadratic fits of D₁/ℱ for |a| ≤ 1.1, the direct mean for |a| ≤ 0.2, the bin average at |a| = 0.95–1.0, and the shape parameters:

| f = 6.4 | ⟨T_r⟩ | §2 var | dv | weighted D₁ fit | unweighted fit | D₁(0) direct | D₁ bin 0.95–1.0 | (c², γ) |
|---|---|---|---|---|---|---|---|---|
| exponential | 104.7 ± 2.8 | 0.142 | 0.160 | 1.58 + 1.91a² | 1.46 + 2.08a² | 1.24 | 3.33 | (0.291, 3.02) |
| IMEX | 98.3 ± 2.5 | 0.135 | 0.153 | 1.84 + 1.78a² | 1.71 + 1.97a² | 1.41 | 3.43 | (0.345, 3.09) |
| backward Euler, noise explicit | 110.9 ± 3.0 | 0.135 | 0.152 | 1.58 + 1.72a² | 1.46 + 1.88a² | 1.24 | 3.13 | (0.300, 3.13) |
| backward Euler, noise implicit | 129.5 ± 3.8 | 0.136 | 0.151 | 1.45 + 1.67a² | 1.31 + 1.86a² | 1.11 | 3.08 | (0.278, 3.15) |
| Crank–Nicolson, noise explicit | 100.5 ± 2.6 | 0.143 | 0.161 | 1.64 + 1.88a² | 1.52 + 2.05a² | 1.27 | 3.38 | (0.303, 3.03) |
| Crank–Nicolson, noise trapezoidal | 107.9 ± 2.9 | 0.143 | 0.160 | 1.57 + 1.87a² | 1.45 + 2.04a² | 1.24 | 3.30 | (0.290, 3.04) |
| paper | 107 ± 5 | 0.16 | | 1.36 + 2.00a² | | | 3.26 (on the curve) | (0.246, 2.96) |

| f = 5.4 | ⟨T_r⟩ | §2 var | dv | weighted D₁ fit | unweighted fit | D₁(0) direct | D₁ bin 0.95–1.0 | (c², γ) |
|---|---|---|---|---|---|---|---|---|
| exponential | 526 ± 31 | 0.068 | 0.073 | 1.64 + 0.75a² | 1.39 + 1.10a² | 1.06 | 2.32 | (0.382, 5.39) |
| IMEX | 497 ± 29 | 0.066 | 0.071 | 1.79 + 0.67a² | 1.53 + 1.06a² | 1.25 | 2.33 | (0.427, 5.49) |
| backward Euler, noise explicit | 588 ± 37 | 0.066 | 0.070 | 1.62 + 0.64a² | 1.38 + 0.97a² | 1.10 | 2.19 | (0.388, 5.53) |
| backward Euler, noise implicit | 704 ± 48 | 0.065 | 0.069 | 1.54 + 0.62a² | 1.29 + 0.98a² | 1.00 | 2.11 | (0.370, 5.55) |
| Crank–Nicolson, noise explicit | 500 ± 29 | 0.068 | 0.073 | 1.67 + 0.74a² | 1.43 + 1.09a² | 1.08 | 2.33 | (0.391, 5.41) |
| Crank–Nicolson, noise trapezoidal | 542 ± 33 | 0.068 | 0.073 | 1.63 + 0.72a² | 1.38 + 1.08a² | 1.05 | 2.31 | (0.383, 5.42) |
| paper | 519 ± 53 | 0.07 | | 1.34 + 0.94a² (implied) | | | | (0.299, 5.15) |

| f = 7.4 | ⟨T_r⟩ | §2 var | dv | weighted D₁ fit | unweighted fit | D₁(0) direct | D₁ bin 0.95–1.0 | (c², γ) |
|---|---|---|---|---|---|---|---|---|
| exponential | 37.0 ± 0.6 | 0.263 | 0.301 | 1.62 + 3.37a² | 1.58 + 3.41a² | 1.40 | 4.66 | (0.234, 1.78) |
| IMEX | 34.7 ± 0.5 | 0.241 | 0.280 | 2.03 + 3.15a² | 1.98 + 3.21a² | 1.79 | 4.79 | (0.304, 1.84) |
| backward Euler, noise explicit | 39.3 ± 0.6 | 0.249 | 0.285 | 1.61 + 3.05a² | 1.57 + 3.10a² | 1.41 | 4.35 | (0.245, 1.87) |
| backward Euler, noise implicit | 45.0 ± 0.8 | 0.256 | 0.290 | 1.38 + 2.99a² | 1.34 + 3.02a² | 1.13 | 4.10 | (0.211, 1.88) |
| Crank–Nicolson, noise explicit | 35.5 ± 0.5 | 0.261 | 0.298 | 1.72 + 3.30a² | 1.69 + 3.33a² | 1.51 | 4.51 | (0.253, 1.80) |
| Crank–Nicolson, noise trapezoidal | 37.9 ± 0.6 | 0.265 | 0.302 | 1.58 + 3.29a² | 1.55 + 3.33a² | 1.38 | 4.58 | (0.232, 1.80) |
| paper | 37.6 ± 1 | 0.30 | | 1.62 + 3.20a² (implied) | | | | (0.240, 1.82) |

- The second-order control agrees with the production stepper everywhere: Crank–Nicolson with trapezoidal noise moves the weighted intercept by −0.01 to −0.04, c² by ≤ 0.002, ⟨T_r⟩ by +2.5 to +3.1% (under 1σ) and dv by ≤ 0.001. The standard-run results are converged across integrator families, not only in Δτ and N (T4.17 now asserts this).
- The first-order schemes move the intercept in both directions, by amounts set by their treatment of the noise term as much as by their overtone damping. The IMEX, which under-damps the first overtone by 12%, raises it by +0.25, +0.16 and +0.41 at f = 6.4, 5.4, 7.4. Backward Euler with the noise explicit over-damps that mode by 7% and lowers D₁ at |a| ≈ 1 by 6–7%, yet leaves the intercept unchanged (−0.01, −0.02, 0.00) while lengthening ⟨T_r⟩ by 6–12%. Fully implicit backward Euler lowers the intercept by −0.14, −0.10 and −0.24 but lengthens ⟨T_r⟩ by 24%, 34% and 21%: the implicit solve attenuates the noise injection itself, which acts like a smaller effective forcing (every D₁ value falls, 7–12% at |a| ≈ 1, with the intercept). Crank–Nicolson with the noise explicit, first order in that term alone, shifts the intercept by +0.06, +0.03, +0.11 and ⟨T_r⟩ by −4 to −5%.
- The hypothesis fails its own test: no scheme reaches the paper's 1.36 ± 0.1 at f = 6.4 (1.34 at 5.4) while keeping ⟨T_r⟩ within 107 ± 5 (519 ± 53). The nearest, fully implicit backward Euler, gives 1.45 and 1.54 at ⟨T_r⟩ = 129.5 ± 3.8 and 704 ± 48, 3.6σ and 2.6σ from the paper's values; at f = 7.4, where the exponential stepper already matches the paper's intercept (1.62), the same scheme would spoil both numbers (1.38; 45.0 against 37.6 ± 1). The variances are insensitive to the scheme (dv within 0.01, §2 variance within 0.02), as DR-5 already noted of the IMEX.
- The shape of D₁(a) is also not an integrator effect: c₂/c₀ of the weighted fit is 1.15–1.25 for every scheme but the IMEX (0.97) at f = 6.4, against the paper's 1.47. The unweighted fit through the bin means is within 0.05–0.10 of the paper's intercept at all three forcings (1.39, 1.46, 1.58 against 1.34, 1.36, 1.62), the weighted one only at f = 7.4.
- Conclusion for DR-5: the papers' unstated time-integration scheme does not explain the intercept gap, within the family of implicit schemes a 2000-era code at Δτ = 0.01 could have used (the diffusion stiffness rules out explicit ones). The remaining candidates are the fit procedure (weighting, the range |a| ≤ 1.1, the paper's histogram binning), the alignment of the noise bins with the grid, and the finite-difference forms of L and of the shear term; DR-5 stays open as a reproduction gap.

### Stretch goal SG-2, alternative quenching laws, gate T4.16 (2026-10-24): recorded, not met; the alternatives are not immaterial

`scripts/sg2_quenching.jl` (`-t 4`, 8 min; the tables are also in `results/sg2_quenching.md`). The nonlinearity is now a pluggable law (`AbstractQuenching` in `src/stepper.jl`, `DECISIONS.md` D3 note, selected by `RunSpec(quenching = …)`): the papers' `AmplitudeQuenching` (c = s(1 − a²); bitwise the previous code, `test/quenching.jl`) and HOS's three alternatives of p. 277, `EnergyQuenching` (q from the toroidal energy Σ w T² instead of a² alone, c = s(1 − E_T)), `InverseQuenching` (α/(1 + q), c = s(1 − a²)/(1 + s a²)) and `LocalQuenching` (α(θ)(1 − κ T(θ)²), c_j = s − (1 + s) κ T_j², with κ = 0.551 calibrated so that the equilibrium has a = 1; the sin θ-mean of q over b₀ would need 2q₀ = 0.842). `floored = true` clamps c ≥ −1: the α-effect is suppressed to zero, never reversed. Standard run at f = 6.4, 5.4 and 7.4, three seeds × 5×10⁴ τ_d pooled. The deterministic equilibria first:

| law | coefficient c | relaxation rate about a = 1 | equilibrium ‖b − a b₀‖ | equilibrium toroidal energy | parity |
|---|---|---|---|---|---|
| amplitude (papers) | s(1 − a²) | 1.526 (= 2μ) | 10⁻¹² | 1 | −1 |
| energy | s(1 − Σ w T²) | 1.077 | 10⁻¹² | 1 | −1 |
| inverse | s(1 − a²)/(1 + s a²) | 0.883 (= 2μ/(1 + s)) | 10⁻¹² | 1 | −1 |
| local | s − (1 + s) κ T(θ)² | 2.403 | 0.68 | 1.21 | −1 |

| f = 6.4 | ⟨T_r⟩ | reversals | §2 var | dv | weighted D₁ fit | D₁(0) direct | D₁ bin 0.95–1.0 | (c², γ) |
|---|---|---|---|---|---|---|---|---|
| amplitude | 104.7 ± 2.8 | 1433 | 0.142 | 0.160 | 1.58 + 1.91a² | 1.24 | 3.33 | (0.291, 3.02) |
| amplitude, floored | 104.7 ± 2.8 | 1432 | 0.144 | 0.161 | 1.59 + 1.90a² | 1.23 | 3.35 | (0.292, 3.02) |
| inverse | 98.2 ± 2.5 | 1528 | 0.202 | 0.223 | 1.31 + 2.36a² | 1.13 | 3.61 | (0.223, 2.79) |
| energy | diverged at τ = 475, ≈ 143, 192 (seeds 1–3) | | | | | | | |
| energy, floored | 36.6 ± 0.6 | 4102 | 0.204 | 0.233 | 2.41 + 0.48a² | 2.52 | 2.95 | (0.600, 4.10) |
| local | diverged at τ = 402, 116 (seeds 1, 2) | | | | | | | |
| local, floored | 32.5 ± 0.5 | 4618 | 0.252 | 0.283 | 2.36 + 1.35a² | 2.25 | 3.65 | (0.483, 3.37) |
| paper | 107 ± 5 | | | 0.16 | 1.36 + 2.00a² | | | (0.246, 2.96) |

| f = 5.4 | ⟨T_r⟩ | §2 var | dv | weighted D₁ fit | (c², γ) | | f = 7.4 | ⟨T_r⟩ | §2 var | dv | weighted D₁ fit | (c², γ) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| amplitude | 526 ± 31 | 0.068 | 0.073 | 1.64 + 0.75a² | (0.382, 5.39) | | amplitude | 37.0 ± 0.6 | 0.263 | 0.301 | 1.62 + 3.37a² | (0.234, 1.78) |
| amplitude, floored | 526 ± 31 | 0.068 | 0.073 | 1.63 + 0.76a² | (0.381, 5.38) | | amplitude, floored | 37.4 ± 0.6 | 0.274 | 0.311 | 1.65 + 3.33a² | (0.240, 1.79) |
| inverse | 497 ± 29 | 0.103 | 0.108 | 1.23 + 1.29a² | (0.255, 4.79) | | inverse | 35.9 ± 0.6 | 0.365 | 0.419 | 1.47 + 3.87a² | (0.198, 1.66) |
| energy | diverged (τ = 1556, seed 1) | | | | | | energy | diverged | | | | |
| energy, floored | 104.8 ± 2.8 | 0.110 | 0.123 | 2.30 − 0.35a² | (0.723, 7.25) | | energy, floored | 18.4 ± 0.2 | 0.355 | 0.387 | 2.54 + 2.04a² | (0.456, 2.20) |
| local | diverged (τ = 558, seed 1) | | | | | | local | diverged | | | | |
| local, floored | 75.5 ± 1.7 | 0.150 | 0.170 | 2.32 + 0.38a² | (0.593, 5.90) | | local, floored | 18.3 ± 0.2 | 0.417 | 0.450 | 2.39 + 3.13a² | (0.359, 1.84) |

- The energy and local laws with h = 1 − q run away at every forcing: an overtone burst lifts their quenching measure past the point where the α-effect reverses sign (E_T > 1 + 1/s = 2.37, or κ T² > 1 locally), and the α-reversed αΩ dynamo is as supercritical as the normal one (the operator is mirror-symmetric under C → −C: Λ = 0.528 at C = −100), so the growing field drives the reversal further and the run overflows within a diffusion time (`DISCREPANCIES.md` DR-8). The papers' law escapes because q depends on a alone: its brief reversed-α episodes (|a| up to 2.5, E_T up to 17 in the standard run) are immaterial, the floored version agreeing with it to 0.1% in ⟨T_r⟩ and 1–4% in the variances.
- The inverse law keeps the reversal rate but not the variability: ⟨T_r⟩ −6%, −6%, −3% at f = 6.4, 5.4, 7.4 (within the gate), the §2 variance and dv +40–50% (0.202 and 0.223 against 0.142 and 0.160 at f = 6.4). The restoring rate about a = ±1 is 2μ/(1 + s) = 0.88 instead of 1.53, so the well is shallower, while the barrier near a = 0, where c → s under either law, is the same. Its D₁ fit at f = 6.4, 1.31 + 2.36a² with (c², γ) = (0.223, 2.79), happens to lie near the paper's 1.36 + 2.00a², (0.246, 2.96), but its variance does not match the paper's 0.16, which the papers' law reproduces.
- The floored energy and local laws change the statistics materially: ⟨T_r⟩ falls to 36.6 and 32.5 τ_d at f = 6.4 (−65%, −69%), 105 and 76 at f = 5.4 (−80%, −86%) and 18.4 and 18.3 at f = 7.4 (−50%, −51%); the variances rise by 35–135%; D₁(0) doubles (2.5 against 1.2) while D₁ at |a| ≈ 1 changes little, so the quadratic fits flatten. A quenching that counts the overtone energy, or the local field, switches the α-effect off exactly during the bursts that drive reversals, removing the dipole's supercritical restoring force when it is needed, so the dipole reverses two to seven times more often.
- Gate T4.16 (kind D): recorded, not met. The inverse law meets the ⟨T_r⟩ half of "⟨T_r⟩, variance within ~20%" and fails the variance half; the energy and local laws fail both, or diverge. HOS's "no material change in the dipole statistics" (p. 277) is reproduced only for the reversal rate under law (2); the paper does not say how its alternatives were bounded, which modes entered q, or how long they ran (DR-8; the resource file carries a Note). The testset keeps the 20% criteria as `@test_broken` beside assertions on what holds: the inverse law's ⟨T_r⟩ within 20%, the floored laws stable with two to seven times the reversal rate, every alternative raising dv, and the floored papers' law within 5% of the unfloored.

### Tier 4, forcing-parameter invariance T4.12 (2026-10-20): recorded, the invariance is broken by the bin width

Three configurations at f²τ_c/N_c² = 2.05×10⁻², three seeds × 5×10⁴ τ_d each (`test/tier4_stochastic.jl`). HOS p. 272 and Appendix E: all statistics depend on (f, τ_c, N_c) only through f²τ_c/N_c², in the δ-correlated, continuum-bin limit.

| (f, τ_c, N_c) | ⟨T_r⟩ | ratio to standard | §2 variance | dv | D₁/ℱ weighted fit | (c², γ) | exact binned D₀ / Eq. 36 |
|---|---|---|---|---|---|---|---|
| (6.4, 0.05, 10), standard | 104.7 ± 2.8 (1433) | 1 | 0.142 | 0.160 | 1.58 + 1.91a² | (0.291, 3.02) | 0.9485 |
| (4.525, 0.1, 10) | 114.8 ± 3.2 (1307) | 1.096 (2.4σ) | 0.134 | 0.153 | 1.59 + 1.72a² | (0.302, 3.13) | 0.9485 |
| (12.8, 0.05, 20) | 75.9 ± 1.7 (1975) | 0.726 (−8.8σ) | 0.153 | 0.175 | 1.87 + 2.12a² | (0.330, 2.91) | 1.0008 |

- Doubling τ_c at fixed f²τ_c/N_c² lengthens ⟨T_r⟩ by 10% (2.4σ) and lowers the variance and the D₁ slope slightly: a finite-correlation-time effect, the noise no longer being δ-correlated on the overtones' damping times (|ℜλ₁| ≈ 6, so τ_c = 0.1 is 0.6 of a damping time).
- Doubling N_c at fixed f²τ_c/N_c² shortens ⟨T_r⟩ by 27% (8.8σ) and raises the D₁ intercept by 18% and dv by 9%. The mechanism is DR-3's: with N_c = 10 the piecewise-constant F delivers 5.2% less D₀ than Eq. 36 (Jensen's inequality over the bin), and with N_c = 20 the deficit is gone (ratio 1.0008, `binned_D0`); the overtones, which have finer latitudinal structure than the dipole, lose more at N_c = 10 and so gain more at N_c = 20, which is the 18% in D₁. The reversal rate's sensitivity to the forcing (d ln⟨T_r⟩/d ln f² ≈ −3.6 from Table II) turns the 5% in D₀ and the 18% in D₁ into the 27%.
- So the invariance holds only in the limit the paper states it for; at the paper's own N_c = 10 the results are specific to that bin width, and comparisons with the paper must use N_c = 10, τ_c = 0.05, as every gate here does. `DISCREPANCIES.md` DR-3 carries the addendum.

### Tier 5, standard-run figures (2026-10-20): T5.1, T5.2, T5.7, T5.8 rendered, the user's judgement pending

`scripts/run_standard.jl` (seed 1) and `viz/hos_fig01.py`, `hos_fig02.py`, `hos_fig08.py`, `hos_fig09.py` write `figures/hos_fig0*.png`. The summariser's reading against the features listed in §3:

- T5.1 (HOS Fig. 1): vertical stripes of overtone fine structure in T(θ, τ); 27 reversal bars in 0–2000 τ_d (paper: of order 20 at a mean spacing of order 100); a around ±1 with excursions to |a| ≈ 1.8; aborted reversals and excursions visible as dips. The colour scale is ±1.7 in our field unit.
- T5.2 (HOS Fig. 2): |a₁| spikes to 1.3–1.4 and |a₂| to 1.0 at the events (paper: ~1.2 and ~0.9), 0.1–0.5 between them; 𝒫 near −1 with spikes to +1 at the reversals and brief spikes toward 0 in between; the reversals of the window come in clusters (506; 538; 652 and 670). The overtone scale is in the unit-toroidal-energy convention of `eigenmodes` (`DECISIONS.md` D4); with it the spike heights match the paper's. The energy axis reaches 7 in our field unit while the paper's runs to 30, so their field unit is 2–3 times ours (the constant of HOS Eq. 18 is not given).
- T5.7 (HOS Fig. 8): the scatter, the bin averages rising from ≈ 1.4 at a = 0 and the a⁴ departure beyond 1.1 all as in the paper; the fitted curve is 1.45 + 2.10a² for this seed against 1.36 + 2.00a².
- T5.8 (HOS Fig. 9): the paper's discrepancy reappears exactly: histogram peak 1.42 at |a| = 0.84 (paper ≈ 1.35 at 0.85) against the Eq. 43 curve peak 1.12 at 0.8 (paper ≈ 1.1 at 0.8), the curve above the histogram for |a| < 0.5 and below it near the peak.
- The user judged Figs. 1, 2, 8, 9 "very comparable" with the originals (2026-10-20). Fig. 8 now also draws the unweighted fit (DR-5).

### Tier 5, event figures (2026-10-20): T5.3, T5.4 passed, judged satisfactory by the user 2026-10-24

`scripts/run_events.jl` selects four events from the event list of `results/standard_f6.4_seed1.h5` by documented rules (the single-crossing reversal with the shortest |a| < a_h stretch and no other event in its window; among reversals with ≥ 1 τ_d below |a| = 0.3, the one with the largest mean |a₁| over that stretch; the two-crossing aborted reversal reaching farthest into the opposite polarity; the excursion with the smallest min |a|), re-runs the simulation from the seed to each 6 τ_d window with `window_run` (`src/run.jl`) at full resolution, checks the window's a against the recorded stream bitwise, and writes `results/events_f6.4_seed1.h5`; `viz/hos_fig03.py` and `hos_fig04.py` draw T(θ, τ) with the zero contour and, below it, a, |a₁| and 𝒫. The paper's windows (τ = 219–225, 1133–1139, 823–829, 1280–1286) belong to its own realisation; ours are chosen by the rules above.

- T5.3 (HOS Fig. 3): top, the reversal at τ = 9304.63–9304.69: a goes from −1.3 to +1.5 in about 0.3 τ_d while |a₁| peaks at 1.3 and 𝒫 reaches +0.8, the zero contour of T sweeping across the whole colatitude range in a fraction of a diffusion time, the paper's "very fast event". Bottom, the reversal at τ = 38258.23–38260.63 (5 crossings): after the first crossing the field is small for 2 τ_d (|a| ≲ 0.3, |a₁| falling from 1.1 to 0.2) with 𝒫 swinging between +1 and −1 and the zero contour of T bulging at intervals of about 1 τ_d, the quadrupole's period 2π/5.88 = 1.07 τ_d; a then recovers along the dotted exp(Λτ) guide over the last 2 τ_d, the roughly exponential recovery of the caption.
- T5.4 (HOS Fig. 4): the aborted reversal at τ = 276.5–280.21 (a to −0.44 and back, 2 crossings) and the excursion at τ = 13711.72–13714.21 (min |a| = 0.0006, no crossing) look alike, as the paper says: a sudden rise of |a₁| to 0.5–1, 𝒫 to +1, the toroidal field weak and symmetric for 1–3 τ_d, then the dipole regrows, with or without a sign change.
- The user judged Figs. 3 and 4 satisfactory (2026-10-24): T5.3 and T5.4 pass.

### Tier 5, sweep figures (2026-10-20): T5.5, T5.6 passed, judged good by the user 2026-10-24

`scripts/run_standard.jl 1 1e4 5.4` and `… 1 2000 7.4` write `results/standard_f5.4_seed1.h5` and `standard_f7.4_seed1.h5`; `viz/hos_fig05.py`, `hos_fig06.py`.

- T5.5 (HOS Fig. 5): with the same noise series as Fig. 1 (seed 1), 6 reversals in 0–2000 τ_d at f = 5.4, the paper's count exactly, including a short reversed interval near 1760–1840 (paper: 1760–1790); 58 at f = 7.4 (paper: several tens).
- T5.6 (HOS Fig. 6): ℜa₁ against a over ±50 τ_d around the reversal at τ = 833 (the one nearest the paper's 680): two clouds centred at a ≈ ±0.9 with |ℜa₁| ≲ 0.5 and a single jump through a ≈ 0, as described.

### Tier 5, HSO and SOH figures (2026-10-20): T5.9, T5.10, T5.11, T5.12, T5.13 passed, judged good by the user 2026-10-24

`scripts/run_standard.jl 1 5e4 5.4` (the reference run, 129 MB) and `… 1 1e4 7 0.3 120` (the second run, 26 MB); `viz/hso_fig03.py`, `hso_fig06.py`, `hso_fig07.py`, `hso_fig08.py`.

- T5.9 (HSO Fig. 3): 7 reversals in 0–3000 τ_d, the paper's count; the longest interval shown is ≈ 880 τ_d (paper ≈ 1100); seed 1 of the whole run gives ⟨T_r⟩ = 515 ± 52 from 97 reversals and dv = 0.072 (paper 570 ± 61 from 88, 0.070). The T(θ, τ) panel has the same vertical striping and polarity blocks as HOS Fig. 1.
- T5.10 (HSO Fig. 6, bottom panels): the |a| histogram (5×10⁴ samples, peak 8.7×10⁴ counts per unit a at |a| ≈ 0.93; paper ≈ 9×10⁴ at ≈ 0.9) with the χ² fits. Eq. 10 gives (γ, c²) = (7.64, 0.470) at 104σ over 110 equal-population bins, and the curve is indistinguishable from HSO's case 2 (8, 0.50); Eq. A.1 gives (4.66, 0.204) at 21σ, on top of HSO's case 4 (4.8, 0.22), which they fitted to Sint-800 and the run together. A.1 is visibly the better description, as HSO say ("many σ" for either, their footnote 4).
- T5.11 (HSO Fig. 8): the polarity bar code to 10⁴ τ_d and a(τ) to 3000 τ_d of the kR 0.3, C_m 120, f 7 run: 26 reversals in 3000 τ_d (paper: of order 20–25 read off the figure, 27 from the mean rate); seed 1 over 10⁴ τ_d gives ⟨T_r⟩ = 118 ± 13 and dv = 0.078 (paper 113 ± 5, 0.081).
- T5.12 (SOH Fig. 3; `scripts/run_ramp.jl`, `viz/soh_fig03.py`): the polarity bar code of the ramped run over 0–5×10⁴ τ_d shows a first normal chron of 9300 τ_d (the paper's realisation: ≈ 19 000), a reversed chron to 14 000, then chrons shortening drastically to a few tens of τ_d by the end, 124 reversals in all; the added panels show the rms variability per 2000 τ_d window rising near-linearly from 16% to 40% between SOH's 17% and 41% lines, and the reversals per window climbing from 0–1 to 22.
- T5.13 (HSO Fig. 7, model curves only; Sint-800 is not in the repository): HSO's recipe applied to the reference run (rectify, running average over τ_av, resample at 1 kyr, 799 points, lags to 50 kyr) for their six (τ_d, τ_av) settings. The behaviour they describe reappears: at τ_av = 5 kyr the curve decays faster with τ_d = 5.9 kyr than with 14.3; with τ_av = 0 the small-lag decay runs between exp(−2Λs) and exp(−2μs) before the smoothing flattens it; τ_av = 10 kyr over-smooths; the zero crossing with τ_d = 10 kyr is near 37 kyr against Sint-800's ≈ 27 kyr, and a different stretch of the run (panel f, from τ = 6000) shows the anticorrelation to −0.3 that HSO say the model "generally" reproduces. The quantitative decay rate is settled from the whole-run autocorrelation in the Tier 4 HSO section (QUESTIONS.md Q5).

