---
title: "Series overview and shared notation: the bistable-oscillator geodynamo papers"
type: overview
covers:
  - hoyng-2001-gafd-geodynamo-bistable-oscillator.md
  - schmitt-2001-pepi-bistable-geodynamo-reversals.md
  - hoyng-2002-pepi-variability-geomagnetic-dipole.md
purpose: "One place for the common model, a concordance of the three papers' notations, the headline numbers side by side, and a map of which paper holds which equation."
canonical_notation: "Hoyng, Ossendrijver & Schmitt (2001), GAFD 94 (the long paper). Where the other two papers differ, the table in section 3 gives the translation."
---

# Series overview and shared notation

Three papers by the same group (P. Hoyng, D. Schmitt, M.A.J.H. Ossendrijver) describe one model from three angles. In chronological order of submission:

| # | Short name | Citation | Role | Resource file |
|---|---|---|---|---|
| 1 | **HOS** | Hoyng, Ossendrijver & Schmitt, *Geophys. Astrophys. Fluid Dyn.* 94, 263–314 (2001). DOI 10.1080/03091920108203410 | The long paper: model definition, simulations, Fokker–Planck derivation, reversal statistics, six appendices, first comparison with data | `hoyng-2001-gafd-geodynamo-bistable-oscillator.md` |
| 2 | **SOH** | Schmitt, Ossendrijver & Hoyng, *Phys. Earth Planet. Inter.* 125, 119–124 (2001). DOI 10.1016/S0031-9201(01)00237-0 | Six-page summary of the method, the Sint-800 amplitude-distribution fit, and a ramped-forcing run explaining the long-term rise of the reversal rate | `schmitt-2001-pepi-bistable-geodynamo-reversals.md` |
| 3 | **HSO** | Hoyng, Schmitt & Ossendrijver, *Phys. Earth Planet. Inter.* 130, 143–157 (2002). DOI 10.1016/S0031-9201(02)00004-3 | Quantitative confrontation with Sint-800: χ² fits, improved p(a), autocorrelation-based turbulent diffusion time, constraint accounting | `hoyng-2002-pepi-variability-geomagnetic-dipole.md` |

These short names (HOS, SOH, HSO) are used below and in the HSO resource; the other two resources call the papers "GAFD 2001", "PEPI 2001" and "PEPI 2002". Page references below are journal pages of the paper named.

## 1. The model in one page (HOS notation)

**Mean-field induction equation** with azimuthal averaging, so ⟨B⟩ is the axisymmetric field and the dynamo coefficients inherit a random part from the finite number of convective cells on a latitude circle (HOS Eq. 1, p. 265; HSO Eq. 1, p. 145):

$$\frac{\partial\langle\mathbf B\rangle}{\partial t}=\nabla\times\left[\mathbf v\times\;+\;\alpha\;-\;(\beta+\eta)\nabla\times\right]\langle\mathbf B\rangle$$

**Reduction to one dimension.** Poloidal/toroidal split, local radial ansatz cos(kr)/(kr) with kR = 0.5, differential rotation with constant radial shear Ω′, αΩ approximation. The field is b = (P, T)ᵀ on colatitude θ ∈ [0, π] with P = T = 0 at the poles (HOS Eqs. 3–9, pp. 267–268; SOH Eq. 1, p. 120):

$$\frac{\partial\mathbf b}{\partial\tau}=G\mathbf b,\qquad G=\begin{pmatrix}L & \alpha/\alpha_0\\ C\,\partial_\theta\sin\theta & L\end{pmatrix},\qquad L=\frac1{\sin\theta}\partial_\theta\sin\theta\,\partial_\theta-\frac1{\sin^2\theta}-(kR)^2$$

$$C=\frac{\alpha_0R^4\Omega'}{\beta^2},\qquad \tau=\frac{t}{\tau_d},\qquad \tau_d=\frac{R^2}{\beta},\qquad \frac{\alpha}{\alpha_0}=\cos\theta$$

The shear entry acts on P as C ∂θ(sin θ P). For kR = 0.5 the critical dynamo number is C₀ = 57.9; the model operates at C_m = 100 (HOS p. 269; SOH writes C_crit = 58, C = 100).

**Eigenmodes.** G b_i = λ^i b_i with adjoint modes b̃_i and bi-orthonormality (b̃_i, b_j) = δ_ij under the inner product ∫₀^π sin θ (P₁*P₂ + T₁*T₂) dθ. Mode 0 at C₀ is a non-oscillatory dipole (λ⁰ = 0); every overtone is oscillatory and damped (HOS Table I, p. 270):

| i | parity of T | λ^i at C₀ = 57.9 | λ^i at C_m = 100 |
|---|---|---|---|
| 0 | antisymmetric (dipole) | 0 | 0.527 |
| 1 | symmetric (quadrupole) | −2.76 + 4.41i | −2.16 + 5.88i |
| 2 | antisymmetric | −8.67 + 4.97i | −8.21 + 6.04i |
| 3 | symmetric | −15.9 + 6.13i | −15.7 + 7.70i |

**Nonlinearity and fluctuations** (HOS Eqs. 17–27, pp. 269–273; HSO Eqs. 2–3, p. 145): global α-quenching on the fundamental-mode energy, plus a random part,

$$\alpha=\alpha_m(1-q)\cos\theta+\delta\alpha,\qquad q=q_0a^2,\qquad \delta\alpha=\alpha_0\frac{f\,F(\theta,\tau)}{\sqrt{2N_c\sin\theta}}$$

with a ≡ ℜa₀ the dipole amplitude in units of its equilibrium value (so the equilibria are a = ±1), s ≡ C_m/C₀ − 1 = 0.728, q₀ = 0.421, F a zero-mean unit-variance random function redrawn every τ_c in each of N_c latitude bins. The model equation actually integrated is

$$\frac{\partial\mathbf b}{\partial\tau}=\left[G+s(1-a^2)E+V\right]\mathbf b,\qquad E=\begin{pmatrix}0&\cos\theta\\0&0\end{pmatrix},\quad V=\begin{pmatrix}0&\delta\alpha/\alpha_0\\0&0\end{pmatrix}$$

with G at C₀. Standard run: f = 6.4, τ_c = 0.05, N_c = 10, 60-point θ grid, Δτ = 10⁻², 5×10⁴ τ_d (HOS footnote 1 p. 273, Table II p. 292).

**Mode equations and Fokker–Planck reduction** (HOS Eqs. 29–44, pp. 282–286; SOH Eqs. 2–4, pp. 121–122; HSO Eqs. 5–11, pp. 146–147):

$$\frac{da_i}{d\tau}=\lambda^ia_i+s(1-a^2)E_{ik}a_k+V_{ik}a_k$$

$$\frac{\partial p}{\partial\tau}=-\frac{\partial}{\partial a}Sp+\frac12\frac{\partial^2}{\partial a^2}Dp,\qquad S\simeq\Lambda(1-a^2)a,\qquad D\simeq\Delta_0a^2+\Delta_1$$

$$U\simeq\tfrac14\Lambda(1-a^2)^2,\qquad p(a)\propto\left(1+\frac{a^2}{c^2}\right)^{\gamma c^2+\gamma-1}e^{-\gamma a^2},\qquad \gamma=\frac{\Lambda}{\Delta_0},\quad c^2=\frac{\Delta_1}{\Delta_0}$$

Λ = 0.527 is the growth rate of the dipole mode at C_m. The term Δ₀a² is the direct effect of δα on the dipole and can never reverse it; Δ₁ = D₁(0) comes from the overtones' mean-square toroidal field during reversals and is the sole cause of reversals (HOS Eq. 49, p. 287). Bistability needs γ > 1.

**Closed-form statistics** (HOS Eqs. 52–60, pp. 288–290):

$$\langle a\rangle\simeq\left(1-\frac1\gamma\right)^{1/2},\qquad \frac{\langle(a-\langle a\rangle)^2\rangle}{\langle a\rangle^2}\simeq\frac{\gamma c^2+\gamma-1}{4(\gamma-1)^2}$$

$$\langle T_r\rangle\simeq\frac\pi\Lambda\left(\frac2{c^2}\right)^{1/2}\left(\frac{\gamma c^2+\gamma-1}{\gamma-1}\right)^{1/2}e^{\psi_m}\;\sim\;\frac{\pi\sqrt2}{\Lambda}\left(\frac{c^2}{1+c^2}\right)^{1/2}e^{K\gamma},\qquad K=-1+(1+c^2)\ln\frac{1+c^2}{c^2}$$

The reversal rate is a Kramers escape rate with f²τ_c/N_c² as effective temperature; polarity intervals are Poisson distributed. Larger forcing means more secular variation and more frequent reversals, which is the observed correlation.

## 2. What each paper adds

- **HOS** derives everything above (Appendices A–F), validates the theory against the simulations (standard run: predicted variance 0.17 vs measured 0.16; predicted ⟨T_r⟩ 94 vs measured 107 ± 5 τ_d), argues from reversal duration that τ_d ~ 3 kyr, shows that a reversal rate of once per ~100 τ_d implies a between-reversal relative variance ~0.16, matching VADM records (0.1–0.2), and fits p(a) to Sint-800 with γ = 5.15, c² = 0.25 (pp. 291–296).
- **SOH** restates the method in six pages, gives the escape-time integral in the form T_r ≃ 4∫₀¹ dx [D(x)p(x)]⁻¹ ∫ₓ^∞ p(y) dy, fits Sint-800 within 1σ with γ = 5.2, c² = 0.26, and adds the ramped run (f from 4 to 6.5 over 5×10⁴ diffusion times) that reproduces the fall of chron lengths from 10⁷ to 10⁵ yr with a factor ~2.5 increase of the forcing parameter (p. 123).
- **HSO** fits p(a) and the improved form (A.1) to Sint-800 and to simulated a(τ) by χ² (Table 1: γ = 6.3, c² = 0.32 optimal for Sint-800; degenerate direction with γ and c² rising together), defines the secular dipole variation dv (Eq. 12; Sint-800 gives 0.071 ± 0.003, reference run 0.070), introduces the autocorrelation technique and obtains τ_d = 11 ± 2 kyr (reference run) or 16 ± 4 kyr (second run), and finds that no single run satisfies the reversal-rate, secular-variation and reversal-duration constraints at once (pp. 149–154).

## 3. Notation concordance

Canonical column is HOS. Blank means the paper does not use the quantity.

| Quantity | HOS (GAFD 2001) | SOH (PEPI 2001) | HSO (PEPI 2002) | Remarks |
|---|---|---|---|---|
| Turbulent magnetic diffusivity | β | η | β | **SOH's η is the turbulent diffusivity**, not the molecular one; HOS/HSO use η for the molecular value (~1.3×10⁴ cm² s⁻¹). |
| Turbulent diffusion time | τ_d = R²/β | "diffusion time" (= R²/η in its notation, formula not written) | τ_d = R²/β | Unit of all dimensionless times in all three papers. |
| Dynamo number, critical / operating | C = α₀R⁴Ω′/β²; C₀ = 57.9, C_m = 100 | C = R_αR_Ω, R_α = α₀R/η, R_Ω = (dΩ/dr)R³/η; C_crit = 58, C = 100 | C_m = α_mR⁴Ω′/β² (Eq. 4) | Same product; SOH rounds C₀. |
| Reference α | α₀ (critical state), α_m (operating) | α₀ | α_m | HOS: α_m(1 − q₀) = α₀. |
| Supercriticality | s = C_m/C₀ − 1 = 0.728 | absorbed into E_ik | absorbed into E_ik (footnote 3) | HOS Eq. 29 has s(1−a²)E_ik; SOH Eq. 2 and HSO Eq. 5 have (1−a₀²)E_ik. |
| Eigenvalues | λ^i (superscript) | λ_i (subscript) | λ^i | Same modes, ordered by decreasing growth rate; all at C₀. |
| Dipole amplitude | a ≡ ℜa₀ | a = a₀ | a ≡ a₀ | Normalised to the nonlinear equilibrium (±1). |
| Fluctuation amplitude | δα = α₀ fF/√(2N_c sin θ) | δα/α₀ = fF/√N, N = 2N_c sin θ | δα ∝ fF/√(2N_c sin θ) | Identical. |
| Cell size in latitude | N_c cells on [0, π] (N_c = 10) | θ_c = π/N_c (θ_c = π/10) | N_c | — |
| Correlation time | τ_c (= 0.05 τ_d) | τ_c (= 5×10⁻²) | τ_c (= 0.05) | — |
| Forcing parameter | f²τ_c/N_c² ("effective forcing parameter"; f alone is called "the forcing parameter") | f²θ_c²τ_c | f²τ_c/N_c² | **See section 4: SOH's value is π² times the HOS/HSO value.** |
| Diffusion prefactor | F ≡ πf²τ_c/2N_c² (3.22×10⁻² standard) | — | — | HOS Eqs. 36–37; not the random function F(θ,τ). |
| Drift, diffusion | S, D = D₀a² + D₁(a) ≃ Δ₀a² + Δ₁ | S, D(a) | S, D ≃ Δ₀a² + Δ₁ | — |
| Growth rate of dipole at C_m | Λ = 0.527; μ = sE₀₀ = 0.763 (bare) | Λ | Λ (0.527 case 1, 0.815 case 2); μ/Λ = 1.45, 1.82 | — |
| Potential | U ≃ ¼Λ(1−a²)² | U(a) = ¼Λ(1−a²)² | U ≃ ¼Λ(1−a²)² | Identical. |
| Shape parameters of p(a) | γ = Λ/Δ₀, c² = Δ₁/Δ₀ | γ, c² ("measures of U(0), D(0)") | γ, c² (Eq. 11) | Identical definitions. |
| Secular variation measure | ⟨(a−⟨a⟩)²⟩/⟨a⟩² between reversals | same, "relative variance" | dv = ⟨(\|a\|−⟨\|a\|⟩)²⟩/⟨\|a\|⟩² (Eq. 12) | HSO uses the rectified record \|a\|; for a well-separated bimodal p(a) the two agree. SOH's Fig. 3 quotes r.m.s. values (square roots). |
| Mean time between reversals | ⟨T_r⟩, Eq. 58 | T_r, escape-time integral | ⟨T_r⟩, Eq. A.2 | Eq. A.2 of HSO is the first form of HOS Appendix F. |
| Parity | 𝒫 = (ℰ_S − ℰ_A)/(ℰ_S + ℰ_A) | — | — | HOS Eq. 28. |
| VADM scale | Sint-800 histogram, 35 bins | a_max = 5.9×10²² A m², 30 bins | A_m = 5.98–6.14×10²² A m², ~100 variable-width bins | VADM at the peak of p. |

## 4. The forcing parameter: three spellings, one quantity

HOS and HSO use f²τ_c/N_c²; SOH uses f²θ_c²τ_c with θ_c = π/N_c, so

$$f^2\theta_c^2\tau_c=\pi^2\,\frac{f^2\tau_c}{N_c^2}.$$

HOS proves (Appendix E) that all statistical properties depend on f, τ_c, N_c only through this combination. Values for the runs that appear in the series (τ_c = 0.05, N_c = 10 throughout):

| f | f²τ_c/N_c² (HOS, HSO) | f²θ_c²τ_c (SOH) | Where used |
|---|---|---|---|
| 4.0 | 8.0×10⁻³ | 0.079 | start of SOH ramped run |
| 4.4 | 9.7×10⁻³ | 0.096 | HOS Table II |
| 5.4 | 1.46×10⁻² | 0.144 | HOS run #2; HSO reference run (Fig. 3); SOH Sint-800 match |
| 6.0 | 1.80×10⁻² | 0.178 | SOH (T_r ≈ 200) |
| 6.4 | 2.05×10⁻² | 0.202 | HOS standard run; SOH Fig. 1 |
| 6.5 | 2.11×10⁻² | 0.209 | end of SOH ramped run |
| 7.0 | 2.45×10⁻² | 0.242 | SOH (T_r ≈ 50); HSO second run (with kR = 0.3, C_m = 120) |
| 7.4 | 2.74×10⁻² | 0.270 | HOS Fig. 5, Table II |
| 8.4 | 3.53×10⁻² | 0.348 | HOS Table II |

(The SOH column is computed from the HOS column by the π² factor; SOH itself prints only f values.)

## 5. Time scales: why the three papers give different years

All model results are in diffusion times τ_d. The papers convert to years with different τ_d:

| Paper | τ_d assumed or inferred | Basis |
|---|---|---|
| SOH | 2000 yr | "estimated as 2000 years" (Fig. 1 caption, p. 121), no derivation |
| HOS | ~3 kyr | model reversals last ~3 τ_d; observed ≲10 kyr (p. 277, 293) |
| HSO | 11 ± 2 kyr (reference run), 16 ± 4 kyr (second run) | autocorrelation of Sint-800 matched to the simulated amplitude (pp. 152–153); supersedes the earlier values |

Consequences for the same dimensionless results:

| Dimensionless result | SOH (2 kyr) | HOS (3 kyr) | HSO (11 kyr) |
|---|---|---|---|
| Standard run ⟨T_r⟩ ≈ 100–107 τ_d | ≈ 2×10⁵ yr | ≈ 3×10⁵ yr | ≈ 1.1 Myr |
| Run f = 5.4, ⟨T_r⟩ ≈ 500–570 τ_d | ≈ 1 Myr | ≈ 1.5 Myr | ≈ 6 Myr (HSO Table 2: 5–7.5 Myr) |
| τ_c = 0.05 τ_d | 100 yr | 150 yr | ~550 yr |
| Reversal duration ~3 τ_d | ~6 kyr | ~9 kyr | 8–30 kyr (HSO, second run) |

This is the origin of HSO's central difficulty: with τ_d from the autocorrelation, the run that fits the Sint-800 amplitude distribution and secular variation reverses too rarely (6 Myr vs. the 95% interval 0.22–3.3 Myr inferred from one reversal in 800 kyr), while the run that reverses often enough has reversals a factor ~2 too long (HSO pp. 153–154). HSO also notes that the physical α_m and β/η shift accordingly: HOS quotes β ~ 10⁶ cm² s⁻¹, β/η ~ 80, α_m ~ −3×10⁻³ cm s⁻¹; HSO quotes β/η ~ 20–30 and α_m ~ −2.5×10⁻⁴ cm s⁻¹ (HOS p. 295–296; HSO p. 155).

## 6. Headline numbers side by side

| Observable | Data | HOS | SOH | HSO |
|---|---|---|---|---|
| Relative variance between reversals (standard run, f = 6.4) | 0.1 (last 130 kyr), 0.15 (Matuyama), 0.2 (Gauss) [HOS p. 294]; 0.15–0.20 [SOH p. 121] | 0.16 measured, 0.17 predicted | 0.16 | — |
| Same, run f = 5.4 | Sint-800: 0.07 [SOH]; dv = 0.071 ± 0.003 [HSO] | — | 0.07 → f = 5.4 | dv = 0.070 (reference run) |
| ⟨T_r⟩, f = 6.4 | (2–3)×10⁵ yr at present | 107 ± 5 τ_d (468 reversals); theory 94 | ≈ 100 diffusion times | — |
| ⟨T_r⟩, f = 5.4 | Brunhes chron 780 kyr | ~500 τ_d | 500 diffusion times | 570 ± 61 τ_d (88 reversals); theory (A.2) almost 2× larger |
| p(a) fit to Sint-800 | 800 points, 1 kyr spacing | γ = 5.15, c² = 0.25 (1σ) | γ = 5.2, c² = 0.26 (1σ, "not optimal") | γ = 6.3, c² = 0.32 (1σ optimum); γ = 4.8, c² = 0.22 with Eq. A.1 (1.4σ) |
| Standard-run shape parameters | — | c² ≃ 0.246, γ ≃ 2.96 (measured D₁) | — | — |
| Reversal duration | ≲10 kyr | ~3 τ_d | fast, few diffusion times | 8–30 kyr (second run), factor ~2 too long |
| Polarity-interval distribution | Poisson | Poisson (verified, no details) | Poisson | Poisson assumed in the 95% interval |

## 7. Where each result lives (equation map)

| Concept | HOS | SOH | HSO |
|---|---|---|---|
| Induction equation | Eq. 1 | — | Eq. 1 |
| 1-D operator G, L, C | Eqs. 5–9 | Eq. 1 (matrix form, with Δ for L) | — |
| Eigenproblem, adjoint, bi-orthonormality | Eqs. 10–15, Table I | — | — |
| α-quenching, fluctuation model | Eqs. 17–21 | unnumbered (p. 120) | Eqs. 2–3 |
| Model equation with s, E, V | Eqs. 22–27 | — | — |
| Parity diagnostic | Eq. 28 | — | — |
| Mode equations | Eq. 29 | Eq. 2 | Eq. 5 |
| Fokker–Planck equation; S, D | Eqs. 32–40 | Eq. 3 | Eqs. 6–7 |
| Potential U | Eq. 41, Fig. 7 | Fig. 2 caption | Eq. 8 |
| Stationary p(a); γ, c² | Eqs. 42–45 | Eq. 4 | Eqs. 9–11; improved form Eq. A.1 |
| Reversal mechanism (D₁(0) only) | Eqs. 46–49 | prose, p. 122 | prose, Section 3 |
| Parabolic approximation, ⟨a⟩, variance | Eqs. 50–56 | — | a_m (p. 149), dv (Eq. 12) |
| Mean time between reversals | Eqs. 57–60, App. F | unnumbered escape-time integral | Eq. A.2 |
| D₁ measured from a run | Eq. 61, Fig. 8 | — | — |
| Standard-run c², γ | Eq. 63 | — | — |
| Autocorrelation / Langevin form | — | — | Eqs. 13–15 |
| χ² fitting procedure | — | — | Section 4.2, Table 1 |
| Mean-field self-consistency estimates | Section 5.2.3–5.2.4 | — | Section 7 |
| Derivation details | Appendices A–F | — | Appendix A |

## 8. Cross-paper caveats worth knowing

- **SOH's "η" is HOS/HSO's "β".** Do not read SOH's turbulent diffusivity as the molecular value.
- **Factor s.** HOS's mode equation carries an explicit s = 0.728 in front of (1−a²)E_ik; SOH and HSO absorb it into E_ik. Numerical E_ik values are therefore not interchangeable between HOS and the other two without the factor.
- **Forcing parameter** differs by π² between SOH and the others (section 4).
- **Secular variation** is "between reversals" in HOS and SOH, and the rectified whole-record dv in HSO; HSO's own simulation gives 0.070 by its definition for the run that SOH assigns 0.07 by the other definition.
- **τ_d** changes by a factor 3–5 across the series (section 5); never carry years across papers, only diffusion times.
- **Sint-800 fits** are not all optima: HOS and SOH quote parameters chosen to lie within 1σ of both the data and a particular run; HSO's Table 1 case 1 is the actual χ² optimum.
- **Typesetting slips transcribed as printed and flagged in the resources:** HOS Table I footnote b prints "~ i(i+1)" without a minus sign and footnote c "i ≥ 0" where i ≥ 1 is meant; HOS p. 296 omits the relation symbol before "2.5" (p. 299 prints ≳ 2.5); HOS Eq. (A.4) lacks the factor s that (A.7) and (A.12) carry; HSO p. 143 prints the Poisson law without a minus sign; HSO Eqs. (13)–(14) use D and Δ₀ where the Itô form has their square roots; HSO Fig. 2 caption attributes poloidal-from-toroidal generation to differential rotation; HSO cites McFadden & Merrill (1995) but omits it from the reference list.

## 9. Reading paths

- **Understand the mechanism in 20 minutes:** SOH resource sections 1, 2, 6, then section 1 of this file.
- **Reimplement the model:** HOS resource sections 4, 5, 11 and 12, with Table I and the standard-run numbers in section 4.6 as validation targets; then HSO section 4.6 for the second run (kR = 0.3, C_m = 120, Λ = 0.815, f = 7).
- **Reproduce the data comparison:** HSO resource sections 5, 6.2–6.7 and 11 (binning, χ² definition, resampling, autocorrelation recipe); the Sint-800 record itself (Guyodo & Valet 1999) is not included in this repository.
- **Audit a derivation:** HOS resource section 6.9 (appendix-by-appendix) and the verification notes in `VERIFICATION.md`.
