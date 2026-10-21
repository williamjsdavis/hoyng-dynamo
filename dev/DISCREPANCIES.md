# Discrepancies between the papers and this implementation

A running log of places where a number in the papers cannot be reproduced from the papers' own formulae and parameters, or where our implementation knowingly differs from the papers' formulae. Each entry says what the paper prints, what we get, how we interpret the difference and how the validation gate treats it. Notation is HOS's. Gate IDs refer to `VALIDATION.md`; the full numbers are in its §4. Add new entries at the end; do not rewrite old ones when a later finding changes the interpretation, append to them.

## Paper values we cannot reproduce

### DR-1. HOS Table II, run 5: the Eq. 56 variance is printed as 1.0 but evaluates to 10.4

- **Where.** HOS p. 292, Table II, row 5 (f = 8.4, c² = 0.228, γ = 1.09), column "variance (56)".
- **Paper.** 1.0. The run's measured relative variance is 0.50. Footnote b: "Run 5 is outside the validity range of the asymptotic analysis of Section 4.5 because γ is so close to 1."
- **Ours.** Eq. 56, (γc² + γ − 1)/[4(γ − 1)²] = 0.33852/0.0324 = 10.45. Rows 1–4 reproduce from the same formula to the printed precision (0.042, 0.083, 0.175, 0.467 against 0.04, 0.08, 0.17, 0.47), and all five ⟨T_r⟩ values reproduce from Eq. 58 (11843, 483, 93.9, 35.4, 34.8 against 12000, 480, 94, 35, 35), so the formula and the parameters are right.
- **Checked.** The tracked text layer `resources/raw-text/hoyng-2001-gafd.txt` reads "1.0" for this cell, so the resource transcription is faithful. No PDF was opened. No reading of Eq. 56 gives 1.0; the exact-distribution dv over a > 0 for these parameters is 0.48, close to the run's 0.50.
- **Interpretation.** The parabolic approximation behind Eq. 56 has no meaning for γ this close to 1 (the well is barely there: a_m = 0.29), which the footnote says in other words. The printed 1.0 may be a placeholder or a slip; it is not a value of Eq. 56.
- **Gate.** T3.3 asserts rows 1–4 and records row 5. Status 2026-10-20: open, nothing further to do unless a later paper or erratum explains the cell.

### DR-2. HSO Table 1, case 4: ⟨T_r⟩ from Eq. A.2 is printed as 1097 but evaluates to 1109

- **Where.** HSO p. 149, Table 1, row 4 (Eq. A.1 fitted to both Sint-800 and the Fig. 3 run; γ = 4.8, c² = 0.22; μ/Λ = 1.45, Λ = 0.527), column ⟨T_r⟩, footnote b "from relation (A.2)". HSO say the explicit expressions used for the table are "straightforward but suppressed" (p. 156).
- **Paper.** 1097 τ_d.
- **Ours.** 1109 τ_d (1.1% high) from our generalisation of A.2 to κ = μ/Λ − 1 ≠ 0 (`src/theory.jl`: a_m from the quadratic of HSO p. 156, ψ_m and ψ_m'' from the A.1 exponent, φ₀'' = −2γ/c² unchanged). Cases 1–3 (Eq. 10, κ = 0) reproduce to 0.1% (1088.7, 927.1, 664.6 against 1090, 928, 665) and case 5 to 0.3% (137.6 against 138, with Λ = 0.815 of the Fig. 8 run), so the machinery is the same as HSO's. The dv column reproduces for all five cases to ±0.001.
- **Interpretation.** ⟨T_r⟩ is exponentially sensitive to the shape parameters, and the table prints them to two figures: γ = 4.8 ± 0.05 alone moves case 4 over 1046–1176 and c² = 0.22 ± 0.005 over 1041–1183. Cases 1–3 were evidently evaluated from the printed values (they are round numbers chosen by hand); case 4 is a fit optimum, presumably evaluated from the unrounded fit output. A 1% target is below the resolution of the printed inputs.
- **Gate.** T3.7 asserts 2% and that 1097 lies inside the γ ± 0.05 bracket, instead of the 1% originally written. Status 2026-10-20: closed as a rounding effect.
- **Also noted.** Case 5 needs Λ = 0.815 (HSO Table 2, the kR = 0.3, C_m = 120 model); with the reference run's Λ = 0.527 it would be 213. Our eigensolver gives Λ = 0.8148 and μ/Λ = 1.821 for that model (HSO: 0.815, 1.82).

## Known systematic differences between our simulation and the papers' formulae

### DR-3. D₀ of the binned noise is 5–6% below the δ-correlated Eq. 36

- **Where.** HOS Eq. 36 with the correlation function B.7, p. 284 and p. 309; our Tier 2 results, `VALIDATION.md` §4.
- **Paper.** D₀ = ℱ ∫ sin θ (P̃₀T₀)² dθ = 0.1136 for the standard run (ℱ = 3.22×10⁻², integral 3.529), i.e. (D₀τ_c)^½ = 0.075.
- **Ours.** The exact D₀ of F piecewise constant on N_c = 10 bins is τ_c f² Σ_k (∫_k w P̃₀T₀/√(2N_c sin θ))² = 0.1078 at N = 60 and 0.1063 at N = 240; measured directly from V₀₀ on the frozen b₀: 0.1087 ± 1.4%.
- **Interpretation.** Jensen's inequality: each squared bin integral is at most θ_c times the integral of the square, with equality only if P̃₀T₀/√sin θ were constant across the bin. Eq. 36 is the continuum limit N_c → ∞ at fixed f²τ_c/N_c². The papers' own simulations carry the same 5% (they used the same bins), so their measured D₁ and the fitted (c², γ) absorb it; it only matters when a theoretical D₀ is combined with a measured D₁.
- **Gate.** T2.3 passes on the measured value against the exact binned formula; the theory module uses Eq. 36 as the papers do. Status 2026-10-19: closed, a known systematic to quote when comparing simulation with theory.

### DR-4. The Gaussian-approximated ⟨T_r⟩ (Eq. 58, A.2) is 8–18% below the exact escape integral (Eq. 57)

- **Where.** HOS Appendix F, Eqs. F.1–F.3; SOH p. 122 (the 4∫₀¹ form); HSO p. 151 ("only approximate").
- **Paper.** Eq. 58 is Eq. 57 with the Gaussian approximation F.2 of both integrands; the papers quote Eq. 58/A.2 as "the" theoretical ⟨T_r⟩ and note the factor-2 gap to the f = 5.4 run.
- **Ours.** Eq. 57 at a′ = 1 exceeds Eq. 58 by 8%, 14%, 17%, 2% for Table II runs 1–4 (12761, 552, 110, 36.1 against 11843, 483, 93.9, 35.4) and by 18% and 12% at (γ, c²) = (3, 0.2) and (6, 0.3); for run 5 it is 46% of Eq. 58. The dependence on a′ is 1% between 0.8 and 1.2.
- **Interpretation.** Not an error, a quantification of the F.2 approximation; it moves the theory toward the measured values of runs 2 and 3 (519 ± 53, 107 ± 5) and is part of the gap discussed in `QUESTIONS.md` Q2. When Tier 4 compares measured reversal times with theory, both Eq. 58 and Eq. 57 are reported.
- **Gate.** T3.8 (diagnostic). Status 2026-10-20: recorded.

### DR-5. The intercept of the D₁ fit depends on how the bin averages are weighted

- **Where.** HOS Eq. 61 and Fig. 8 (p. 291–292): "histogram of bin averages, fitted curve 1.36 + 2.00a² for a ≤ 1.1"; our Tier 4 results, `VALIDATION.md` §4.
- **Paper.** D₁/ℱ ≃ 1.36 + 2.00a², hence c² = 0.246, γ = 2.96 (Eq. 63); the histogram "rises from ≈ 1.4 at a = 0" (figure read-off).
- **Ours.** Three seeds × 5×10⁴ τ_d pooled: an unweighted least-squares fit through the bin means for |a| ≤ 1.1 gives 1.46 + 2.08a²; weighting the bins by their counts gives 1.58 + 1.91a²; the direct mean of the samples with |a| ≤ 0.2 gives 1.24 ± 0.1. Hence c² = 0.261 (unweighted) or 0.291 (weighted) against 0.246; γ is 2.93 or 3.02 against 2.96 either way.
- **Interpretation.** D₁(a) is only approximately quadratic: the bin averages bend upward beyond |a| ≈ 0.8 (HOS note the a⁴ term beyond 1.1), so the intercept is an extrapolation whose value depends on the weighting. The paper does not say how it fitted; a curve through the histogram is an unweighted fit, which is what the gate uses, and the paper's 1.36 lies between our unweighted fit and our direct small-|a| estimate. The physically relevant number, D₁(0) = Δ₁, is best estimated directly from the small-|a| bins once enough samples exist; the forcing sweep (step 7) will pool more.
- **Gate.** T4.3 and T4.4 pass at 10% with the unweighted fit (7.6% and 6%); the weighted fit would fail both. `summarize` reports all three estimates. Status 2026-10-20: open; revisit with the sweep statistics (T4.7).

