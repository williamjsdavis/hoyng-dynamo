# Verification log

How the resource files were produced and checked. Dates are 2026-10-15.

## Pipeline

1. **Extraction.** Text layers extracted with PyMuPDF; every page rendered to a 2× PNG; equation-dense regions re-rendered at 3–6× as crops. Papers 2 and 3 are vector PDFs with clean text. Paper 1 is a scan whose OCR layer omits every displayed equation and garbles Greek letters, so it was treated as image-only for anything mathematical.
2. **Writing.** One agent per paper wrote the resource from the page images, following a fixed 15-section template, keeping a scratch list of every transcribed equation, table and caption with its journal page.
3. **Independent verification.** A separate agent per paper (two for Paper 1, split at journal page 290/291 with disjoint edit rights) re-read every page, built a checklist of every equation, number, page reference, figure description and glossary entry in the resource, checked each against its own zoomed crops, fixed errors in place, and reported. Verifiers were told to assume errors were present and to treat omissions as errors.
4. **Coordinator spot checks.** The coordinating session independently compared selected pages with the final text (listed below).
5. **DOIs** confirmed against the Crossref API.

## Paper 2 — Schmitt, Ossendrijver & Hoyng (2001), PEPI 125, 119–124

- Verifier checked ~140 items: the 5 equations (Eqs. 1–4, the escape-time formula, U(a)), ~45 numeric values, ~60 page references, 4 figure descriptions, 25 glossary rows, all 26 reference entries, frontmatter and page mapping.
- Errors: 0 equation, 0 numeric, 0 page-reference.
- Changes: 7 additions of content present in the paper but missing from the draft (omitted citations on p. 119; the criterion for choosing C on p. 120; a qualifier that the length unit R is implied rather than stated; the axisymmetry argument for VADM ∝ a₀ on p. 123; the Discussion's conclusions on p. 124 that reversals are ultimately nonlinear in origin and have no unique cause; the stated reason the method generalises to 2D). The operator-ordering question in Eq. (1) was resolved: the print `C(∂/∂θ) sin θ` can only mean C ∂θ(sin θ P) for a radial-shear Ω-effect, and the resource keeps the flag that the paper does not spell this out.
- The writer's own derived numbers (forcing-parameter values, the Eq. 4 peak position, the one-well variance, the ln T_r slope) were recomputed and reproduced.
- Coordinator spot check: p. 122 (Eqs. 3–4, U(a), T_r formula, Fig. 2 caption) matches the resource exactly.
- Verifier confidence: high.

## Paper 3 — Hoyng, Schmitt & Ossendrijver (2002), PEPI 130, 143–157

- Verifier checked ~210 items: 17 equations (1–15, A.1, A.2) plus the unnumbered appendix relations, every value in Tables 1–2 and their footnotes, both run captions (Figs. 3, 8), all Fig. 7 panel labels, ~90 page references (all outline ranges tested against page headers), 9 figure descriptions, 50 glossary rows, the 33-entry reference list, frontmatter and page mapping.
- Errors: 0 equation, 0 numeric, 1 page-reference (a sentence spanning pp. 145–146 was cited as p. 146 only).
- Changes: that fix plus 2 short additions (the p. 144 argument that Poissonian chron lengths make rapid-variability-induced reversals "almost unavoidable"; the p. 156 remark that the exponent of Eq. A.1 contains an a⁴ term).
- Independent re-derivations: Eq. (10) from Eq. (9); C_α/C_Ω ≈ 8.5×10⁻⁴; β/η ≈ 27 and 19 for the two cases; Fig. 3 has 7 sign changes with longest interval ≈ 1100 τ_d; Fig. 8 ≈ 22 sign changes in 0–3000 τ_d (the resource's "20–25 [uncertain]" is kept as a visual estimate).
- Typesetting slips in the original are transcribed as printed and flagged: Poisson law without a minus sign (p. 143), Eq. (1) with an italic x for the cross product, Eqs. (13)–(14) written with D and Δ₀ rather than their square roots, Fig. 2 caption attributing poloidal-from-toroidal generation to differential rotation, McFadden & Merrill (1995) cited on p. 148 but absent from the reference list.
- Coordinator spot check: Table 1 (p. 149), all five rows and three footnotes, matches the resource exactly.
- Verifier confidence: high.

## Paper 1 — Hoyng, Ossendrijver & Schmitt (2001), GAFD 94, 263–314

Verified in two halves because of its length (55 PDF pages, 63 numbered equations plus 43 appendix equations).

**Journal pages 291–314 and the figure inventory (agent verifier).** ~260 items: all 47 equations in scope ((61)–(63), (A.1)–(A.19), (B.1)–(B.13), (C.1)–(C.3), (D.1)–(D.3), (E.1)–(E.2), (F.1)–(F.3), the unnumbered identity on p. 314) from independent 4× crops; Table I (16 values, 3 footnotes) and Table II (35 values, 2 footnotes); all 11 figure captions and the 3 colour plates; every number in Sections 5–6; all 60 reference entries; all outline page ranges; frontmatter and page mapping (tested on 8 PDF pages).
- Errors: 0 equation, 2 numeric (both figure read-offs), 0 page-reference.
- Changes (9): Fig. 5 shows six reversal bars for f = 5.4, not four; Fig. 2 bar positions refined and kept `[uncertain]`; Table I footnotes now quoted exactly as printed ("~ i(i+1)" with no minus sign; "i ≥ 0"); Fig. 10 panel title and error bars added; reference count corrected to 60; one journal name corrected (Van Kampen 1976, Phys. Reports 24C); the paper's actual wording on additive-noise models restored; Appendix D bra-ket notation changed to the paper's round brackets; one inference in section 15 marked as such.
- Confirmed as printed: Eq. (A.4) lacks the factor s present in (A.7)/(A.12); the p. 296 line has no relation symbol before "2.5"; the sum 5.52 ≈ 2.00 + 3.529.

**Journal pages 263–290 (coordinator, in-session, after the first agent assigned to this half was cut off by a usage limit before recording any finding).** Checked against the 47 zoomed crops that agent had rendered plus one new crop: Eqs. (1)–(60) character by character (operator matrix, boundary conditions, inner product, adjoint relations, Eq. (16) definitions, quenching Eqs. (17)–(19), fluctuation Eqs. (20)–(21), model Eqs. (22)–(27), parity Eq. (28), mode Eqs. (29)–(31), Fokker–Planck Eqs. (32)–(40), potential and stationary solution Eqs. (41)–(45), reversal mechanism Eqs. (46)–(49), parabolic approximation and statistics Eqs. (50)–(60)); Table I in full; footnotes 1–3; Fig. 1–7 captions; the numbers kR = 0.5, C₀ = 57.9, C_m = 100, s = 0.728, q₀ = 0.421, τ_c = 0.05, N_c = 10, 60 grid points, Δτ = 10⁻², E₀₀ = 1.048, μ = 0.763, Λ = 0.527, the integral 3.529, 2π/5.88 ≃ 1.1, ~0.5 τ_d fast events, ~3 τ_d reversal duration, τ_d ~ 3 kyr, 100 τ_d ~ 3×10⁵ yr, the aborted-reversal and excursion times in the Fig. 1 caption, and the Fig. 6–7 parameters (f = 5.4, τ ≈ 680; γ = 3, c² = 0.2); the implementation guide and glossary against the above.
- Errors: 0 equation, 0 numeric, 0 page-reference.
- Changes (2): added footnote 3 of p. 296 (maximum possible |α| is the r.m.s. convective speed, for maximally helical convection) to section 4.4 and the α_m row of section 4.6; "4N_c²/π ≈ 130" changed to the printed "~ 130".
- Two issues the first agent had flagged as out of scope (Fig. 5 reversal count; the Table I footnote sign) were found already corrected in the file.

Overall confidence for Paper 1: high. Residual approximations are figure read-offs, all marked.

## Overview file

`00-series-overview-and-shared-notation.md` was written by the coordinator from the three verified resources only; every number in it is traceable to one of them. The forcing-parameter conversion table was computed from f, τ_c = 0.05 and N_c = 10 (the SOH column is π² times the HOS column); the year conversions in section 5 are products of the quoted dimensionless values with the quoted τ_d.

## Implementation-time checks (2026-10-20)

Not a transcription pass; recorded because two Notes were added to the resources. While implementing the theory module (`src/theory.jl`, validation gates T3.3 and T3.7), Eq. (56) for Table II run 5 of Paper 1 gave 10.4 against the printed 1.0. The tracked text layer `raw-text/hoyng-2001-gafd.txt` was checked and reads "1.0", so the resource is faithful to the print; no PDF was opened. A **Note:** to that effect was added after the Table II discussion in the Paper 1 file (section 6.8), and a **Note:** on the reproduction of Table 1's ⟨T_r⟩ column (cases 1–3 and 5 to 0.3%, case 4 to 1.1%) after the case 4 bullet in the Paper 3 file (section 6.4). Nothing transcribed from the papers was changed. Details in `dev/DISCREPANCIES.md`.

A third **Note:** was added on 2026-10-20 after Eq. (28) in the Paper 1 file (section 4, the parity): the paper's "except during (aborted) reversals and excursions 𝒫 ≃ −1" is qualitative, and our reproduction gives a median of −0.82 outside the events, consistent with the overtone amplitudes the paper's own Fig. 2 shows. Nothing transcribed was changed; no PDF was opened. Details in `dev/DISCREPANCIES.md` DR-7.

A fourth **Note:** was added on 2026-10-24 after the alternative-nonlinearity paragraph of the Paper 1 file (section 4.5): alternatives (1) and (3) of p. 277, quenching by the total toroidal energy and local quenching with h = 1 − q, run away in our reproduction because the α-effect reverses sign under an overtone burst, and that the paper's "no material change" is reproduced only for the reversal rate under alternative (2): the floored forms of (1) and (3) reverse two to seven times more often, and (2) raises the variability by 40–50%. Nothing transcribed was changed; no PDF was opened. Details in `dev/DISCREPANCIES.md` DR-8 and `dev/VALIDATION.md` T4.16.
