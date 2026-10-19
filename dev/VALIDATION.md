# Validation gates

Accepted by the user on 2026-10-18: the gates, the tolerance rules and the operational definitions below. Targets come from the three papers via `resources/`; notation is HOS's unless stated; *(derived here)* marks our own derivation rather than a statement in a paper. Gate IDs (T0.1, …) are the ones referenced from `QUESTIONS.md` and `DECISIONS.md`.

Status: no gate has been run yet. Results will be recorded against each gate once code exists.

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
| T2.1 | F sample statistics: zero mean, unit variance, constant for 5 steps, independent across bins and intervals | — | mean within 3σ of 0 and variance within 1% at 10⁶ draws; successive-interval and cross-bin correlations consistent with 0 | Eq. 21, p. 272 | A |
| T2.2 | Measured correlation function of δα/α₀ matches Eq. B.6: triangular in lag and in angle, amplitude f²/(2N_c √(sin θ sin θ')) | — | within 2% at zero lag | B.6, p. 308 | A |
| T2.3 | D₀ measured directly: Var[∫₀ᵀ V₀₀ dτ]/T for T ≫ τ_c, with V₀₀ from Eq. 31 on the frozen mode b₀ | 0.114 = 3.529 F at f = 6.4 | report; expect agreement within ~10% (δ-correlation of B.7 vs 10 finite bins) | Eqs. 36, B.7–B.8 | D |
| T2.4 | Short-correlation-time parameter (D₀τ_c)^(1/2) | 0.075 | from T2.3 | C.3, p. 310 | A |
| T2.5 | Reproducibility: seeded RNG reproduces a(τ) bitwise; one noise series can be reused across f values; record and stream modes (`DECISIONS.md` D7) give identical states and aggregator values, and replaying a recorded a(τ) through the aggregators reproduces the online values | — | bitwise | p. 280; D7 | A |

### Tier 3: closed-form theory (formula module, independent of the simulation)

| ID | Check | Target | Pass criterion | Source | Kind |
|---|---|---|---|---|---|
| T3.1 | Eq. 43 equals Eq. 42 evaluated numerically with S = Λ(1−a²)a, D = Δ₀a² + Δ₁ | — | 10⁻⁸ after normalisation | Eqs. 42–44 | A |
| T3.2 | Peak a_m, ψ_m, ψ_m'' (Eqs. 52–54) vs numerical maximisation of ψ | — | 10⁻⁸ | Eqs. 50–54 | A |
| T3.3 | Variance Eq. 56 for the five Table II (c², γ) pairs | 0.04, 0.08, 0.17, 0.47, 1.0 | last printed digit | Table II | A |
| T3.4 | ⟨T_r⟩ Eq. 58 for the five Table II pairs | 12000, 480, 94, 35, 35 | 1% (hand check: 93.9 for run 3) | Table II | A |
| T3.5 | Eq. 58 with Λ = 0.527 for HSO Table 1 cases 1–3 (γ, c²) = (6.3, 0.32), (8, 0.50), (7, 0.45) | 1090, 928, 665 | 0.5% (hand check: 1089, 927, 665) | HSO Table 1 | A |
| T3.6 | dv (HSO Eq. 12) from Eq. 10 over a > 0, numerically, cases 1–3 | 0.076, 0.070, 0.079 | ±0.001 | HSO Table 1 | A |
| T3.7 | Eq. A.1 and the generalised A.2 with μ/Λ = 1.45 (case 4) and 1.82 (case 5) | ⟨T_r⟩ 1097, 138; dv 0.071, 0.091 | 1%; ±0.001 | HSO Table 1, App. A | A |
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
