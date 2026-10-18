# Paper resources: the bistable-oscillator geodynamo series (Hoyng, Schmitt & Ossendrijver, 2001–2002)

Structured, verified summaries of the three PDFs in `papers/`, written so that a coding agent can load one file and have the paper's model, equations, parameters, results and caveats as context without opening the PDF. All three papers are by the same group and describe one model from three angles; read this file first, then `00-series-overview-and-shared-notation.md`, then the paper you need.

## What is here

| File | Contents |
|---|---|
| `README.md` | This index: what each paper contributes, reading order, conventions, provenance. |
| `VERIFICATION.md` | How each file was produced and checked, item counts, every correction made, and the residual caveats. |
| `00-series-overview-and-shared-notation.md` | The common model in one place, a unified symbol table reconciling notation differences between the papers, the headline numbers side by side, and which paper to consult for what. |
| `hoyng-2001-gafd-geodynamo-bistable-oscillator.md` | Paper 1 (chronological). Hoyng, Ossendrijver & Schmitt, *Geophys. Astrophys. Fluid Dyn.* 94, 263–314 (2001). The long paper: full model, simulations, Fokker–Planck derivation, reversal statistics, six appendices. |
| `schmitt-2001-pepi-bistable-geodynamo-reversals.md` | Paper 2. Schmitt, Ossendrijver & Hoyng, *Phys. Earth Planet. Inter.* 125, 119–124 (2001). Six-page summary of the method with the first Sint-800 comparison and a ramped-forcing run. |
| `hoyng-2002-pepi-variability-geomagnetic-dipole.md` | Paper 3. Hoyng, Schmitt & Ossendrijver, *Phys. Earth Planet. Inter.* 130, 143–157 (2002). Quantitative confrontation with the Sint-800 VADM record: χ² fits of the amplitude distribution, autocorrelation-based estimate of the turbulent diffusion time, constraint accounting. |
| `raw-text/*.txt` | The raw text layer of each PDF, for keyword search only. Headers state provenance and quality. The GAFD file is OCR of a scan and has no equations. |

Source PDFs (unchanged):

| Resource | PDF | Pages | DOI |
|---|---|---|---|
| Paper 1 | `papers/The_Geodynamo_as_a_Bistable_Oscillator.pdf` | 55 (52 numbered + 3 colour plates) | 10.1080/03091920108203410 |
| Paper 2 | `papers/Magnetic_Field_Reversals_and_Secular_Var.pdf` | 6 | 10.1016/S0031-9201(01)00237-0 |
| Paper 3 | `papers/geophysical_1_Breno.pdf` | 15 | 10.1016/S0031-9201(02)00004-3 |

DOIs were confirmed against Crossref on 2026-10-15.

## The series in three sentences

The geodynamo is modelled as a one-dimensional (colatitude only) αΩ mean-field dynamo whose α-effect carries a rapidly fluctuating random part, because the azimuthal average that defines the mean field runs over a finite number of convective cells. The fluctuations excite damped, oscillatory overtones, and their back-reaction on the fundamental dipole mode turns the dipole amplitude into an overdamped particle in a bistable potential ∝ (1 − a²)², so the field wobbles in one well (secular variation) and occasionally hops to the other (a reversal). This single mechanism links the reversal rate to the secular variation, reproduces the Poisson statistics of polarity intervals, and fits the Sint-800 virtual-axial-dipole-moment distribution.

## Which file to read for what

| You need | Go to |
|---|---|
| A fast, correct mental model of the mechanism | Paper 2, sections 1–2 and 6; then overview |
| The governing equations, boundary conditions, eigenmode problem, Table of eigenvalues | Paper 1, section 4 and 12 (glossary) |
| How the stochastic α is constructed and what the "forcing parameter" is | Paper 1, section 4.4; overview section on notation (the three papers write it differently) |
| The Fokker–Planck equation, drift/diffusion coefficients, the potential, p(a), variance and reversal-rate formulas | Paper 1, section 6.3–6.7 (derivations in 6.9); Paper 2, section 6 for the compact form |
| Numbers to validate a reimplementation against | Paper 1, sections 4.6, 6.1–6.2 and 11; Paper 3, section 4.6 (the two runs used there) |
| The Sint-800 data set and how it was binned, fitted and resampled | Paper 3, sections 5 and 6.2–6.4 |
| The turbulent diffusion time τ_d in years and why the three papers assume different values | Paper 3, section 6.6; overview "time scales" table |
| Known discrepancies and open problems | Paper 3, sections 6.7 and 15; Paper 1, section 15 |

## Conventions used in the resource files

- Each file follows the same 15-section template, so the same information is in the same place in every file (model in §4, methods in §5, results in §6, figures in §7, implementation guide in §11, glossary in §12, outline in §13).
- Math is LaTeX. Equation numbers are the paper's own, e.g. "(Eq. 43)". Page references are journal page numbers as printed; the frontmatter of each file gives the PDF-page to journal-page mapping.
- Every quantitative statement carries a page reference. Values read off figure axes rather than stated in text are marked approximate or `[uncertain]`. Nothing is marked `[unreadable]` in the final versions.
- Sentences starting with **Note:** are commentary by the summariser (reconstructions, consistency checks, implementation remarks). Everything else is the paper's own content.
- Section 11 of each file states explicitly where the paper does not specify something an implementer would need (time integrator, random-number distribution, reversal-counting window, etc.) rather than filling the gap silently.

## Provenance and verification

- Text was extracted with PyMuPDF. Papers 2 and 3 are vector PDFs with clean text; Paper 1 is a scanned PDF whose OCR layer drops every displayed equation, so for all three papers equations, tables and figure captions were transcribed from rendered page images (2× pages, 4–5× zoomed crops), not from the text layer.
- Each resource was written by one agent from the full page images and then checked independently: a separate verifier re-read every page and compared every equation, number and page reference against its own zoomed crops, fixing errors in place. Papers 2 and 3 and journal pages 291–314 of Paper 1 were verified by agents; journal pages 263–290 of Paper 1 (the model and the statistical analysis) were verified by the coordinating session itself after the agent assigned to them was cut off by a usage limit. Across all three papers the verifiers found no equation transcription errors; the corrections were figure read-offs, one page range, wording fidelity and a few omissions. Details and counts are in `VERIFICATION.md`.
- The overview file was written last, from the verified resources only.
- Residual risk: a handful of figure read-offs remain approximate (marked), and typesetting slips in the originals (e.g. a missing relation symbol on p. 296 of Paper 1, a missing minus sign in the Poisson law on p. 143 of Paper 3) are transcribed as printed and flagged with a Note rather than silently corrected.
