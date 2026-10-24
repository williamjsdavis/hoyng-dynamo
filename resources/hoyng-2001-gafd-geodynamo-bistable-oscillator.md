---
title: "The geodynamo as a bistable oscillator"
authors: ["P. Hoyng", "M. A. J. H. Ossendrijver", "D. Schmitt"]
year: 2001
journal: "Geophysical & Astrophysical Fluid Dynamics"
volume: 94
pages: "263–314"
doi: "10.1080/03091920108203410"
source_pdf: papers/The_Geodynamo_as_a_Bistable_Oscillator.pdf
pdf_pages: 55
page_mapping: "PDF page 1 = journal page 263; journal page = PDF page + 262 for PDF pages 1–52 (verified on PDF pages 1, 5, 6, 13, 15, 29, 36, 37, 39, 40). PDF pages 53–55 are unnumbered Colour Plates I–III (colour versions of Figs. 1, 3 and 4, referring back to pp. 274, 278, 279)."
series: "Hoyng–Schmitt–Ossendrijver bistable geodynamo series, paper 1 of 3 (chronological)"
keywords: ["Geodynamo", "Reversals", "Secular variation", "Dynamo theory", "Helicity fluctuations", "Stochastic processes", "Bistable oscillator"]
extraction_notes: "Scanned PDF with an OCR text layer; prose taken from the text layer and checked against page images, all equations (1)–(63), (A.1)–(A.19), (B.1)–(B.13), (C.1)–(C.3), (D.1)–(D.3), (E.1)–(E.2), (F.1)–(F.3), Tables I–II and all figure captions transcribed from 4x zoomed page images. Received 1 June 2000; in final form 22 December 2000."
---

# The geodynamo as a bistable oscillator

## 1. One-paragraph summary

The paper builds a deliberately simple nonlinear mean-field αΩ dynamo for the Earth's outer core (one spatial co-ordinate, colatitude θ; poloidal and toroidal scalars P, T) and adds to the α-effect a random component δα that fluctuates on the eddy turnover time τ_c. The fluctuations are argued to be unavoidable once the mean field is defined by an azimuthal average over a finite number of convective cells. Long simulations (5×10^4 turbulent diffusion times per run, i.e. ~150 Myr with τ_d ≈ 3 kyr; Table II p. 292, p. 298) produce a dominantly dipolar field whose amplitude varies rapidly, with sudden polarity reversals (during which the field becomes almost quadrupolar), aborted reversals, excursions and Poisson-distributed polarity intervals (pp. 273–281). A Fokker–Planck equation for the amplitude a of the fundamental dipole mode is derived (Appendices A–B); it shows that a behaves as a heavily damped, randomly forced particle in a bistable potential U ≃ ¼Λ(1−a²)², and that reversals are caused solely by diffusion contributed by the overtones (pp. 282–287). Closed-form expressions follow for the stationary amplitude distribution p(a), the relative variance of a between reversals and the mean time between reversals ⟨T_r⟩ (Kramers escape, exponential in the forcing parameter f²τ_c/N_c²) (pp. 288–290). Validation against simulations (Table II) and against geomagnetic data shows that the present reversal rate (once per (2–3)×10^5 yr) implies a ~15 % secular variation of the axial dipole moment, as observed, and that p(a) fits the Sint-800 VADM distribution (pp. 291–296).

## 2. Key takeaways

- The model is the αΩ reduction (Eqs. 5–8) with kR = 0.5, critical dynamo number C_0 = 57.9 and operating value C_m = 100, giving supercriticality s = C_m/C_0 − 1 = 0.728 and equilibrium quenching q_0 = 0.421 (p. 269, p. 272).
- The fundamental mode at C_0 is a non-periodic dipole (λ^0 = 0); all overtones are periodic and damped, e.g. λ^1 = −2.76 + 4.41i at C_0 and −2.16 + 5.88i at C_m (Table I, p. 270). At C_m the fundamental mode grows at Λ = 0.527 per diffusion time (Eq. 34, p. 283).
- Fluctuations: δα = α_0 f F(θ,τ)/√(2N_c sin θ), F random with zero mean and unit variance, redrawn once every τ_c independently in each latitude interval θ_c = π/N_c; standard run τ_c = 0.05, N_c = 10, f = 6.4; all statistics depend on f, τ_c, N_c only through f²τ_c/N_c² (Eq. 21, pp. 271–272; Appendix E, p. 311).
- Standard run (60-point θ grid, Δτ = 10^−2, 5×10^4 diffusion times): 468 reversals, ⟨T_r⟩ = 107 ± 5 τ_d, relative variance of a between reversals 0.16 (p. 291, Table II p. 292).
- The dipole amplitude obeys a Fokker–Planck equation ∂p/∂τ = −∂(Sp)/∂a + ½∂²(Dp)/∂a² with drift S ≃ Λ(1−a²)a and diffusion D = D_0a² + D_1(a); D_0a² is the direct effect of δα on a and never causes reversals; D_1(a), from the overtones' mean-square toroidal field ⟨T_+²⟩|_a, is what makes reversals possible (Eqs. 32–37, 49; pp. 282–287).
- Stationary amplitude distribution p(a) ∝ (1 + a²/c²)^{γc²+γ−1} exp(−γa²), with γ = Λ/Δ_0 (dimensionless barrier height) and c² = Δ_1/Δ_0 (dimensionless diffusion at a = 0); bistable only for γ > 1 (Eqs. 43–44, pp. 285–286, 288).
- Relative variance between reversals ⟨(a−⟨a⟩)²⟩/⟨a⟩² ≃ (γc²+γ−1)/[4(γ−1)²] (Eq. 56, p. 289); mean time between reversals ⟨T_r⟩ ≃ (π/Λ)(2/c²)^{1/2}[(γc²+γ−1)/(γ−1)]^{1/2} exp(ψ_m) ~ exp(Kγ) (Eqs. 58–60, p. 290): a Kramers escape rate with f²τ_c/N_c² acting as effective temperature.
- Standard case: c² ≃ 0.246, γ ≃ 2.96, predicted variance 0.17 (measured 0.16), predicted ⟨T_r⟩ ≃ 94 (measured 107 ± 5) (Eq. 63, p. 291). Halving f²τ_c/N_c² (f from 6.4 to 4.4) raises ⟨T_r⟩ by a factor ~100, to ~20 Myr (p. 299).
- A typical model reversal lasts ~3τ_d; with the observed ≲10 kyr this gives τ_d ~ 3 kyr and ⟨T_r⟩ ≈ 100τ_d ≈ 3×10^5 yr (p. 277, p. 293).
- Geomagnetic test: a reversal rate of once per ~100τ_d implies a relative variance ~0.16 of the dipole amplitude between reversals; VADM records give 0.1 (last 130 kyr), 0.15 (Matuyama, 2–2.5 Myr) and 0.2 (Gauss, 2.5–3 Myr) (p. 294).
- Sint-800 VADM histogram (800 points, 1 kyr spacing, 35 bins) is fitted by Eq. (43) with γ = 5.15, c² = 0.25 (χ² at 1σ), close to run #2 (f = 5.4; γ = 5.15, c² = 0.3, χ² at 2σ), whose ⟨T_r⟩ ~ 500τ_d ~ 1.5 Myr (Fig. 10, p. 295).
- Self-consistency estimates: β ~ 10^6 cm² s^−1, η ~ 1.3×10^4 cm² s^−1, β/η ~ 80; C_Ω ~ −100, C_α ~ −1; α_m ~ −3×10^−3 cm s^−1; τ_c = 0.05 ↔ 150 yr; λ_c ~ 500 km; the r.m.s. fluctuation exceeds the mean α by a factor ≳ 2.5 (pp. 295–296, 299).

## 3. Context and motivation

The geomagnetic field has existed for at least 3.5×10^9 yr, is largely dipolar with the axis within ≲10° of the rotation axis, and reverses polarity suddenly: during a reversal the dipole intensity becomes small, the axis may swing rapidly, and the whole process is completed in ≲10 kyr; the present mean interval between reversals is (2–3)×10^5 yr. Reversals can also abort (aborted reversal / excursion). The authors cite Jacobs (1994) and Merrill et al. (1996) as the standard documentation (p. 264).

State of the art described (p. 264): three-dimensional hydromagnetic simulations (Glatzmaier and Roberts 1995a,b, 1996; Kageyama and Sato 1997; Kuang and Bloxham 1997; Olson et al. 1999) show self-sustained dynamo action; the solid inner core stabilises the field against too-frequent reversals (Hollerbach and Jones 1993); in the 2.5-D model of Sarson and Jones (1999) stabilisation comes from a convection-driven meridional circulation whose fluctuations cause reversals; Glatzmaier et al. (1999) found that reversal frequency and polar paths depend strongly on the CMB heat-flux pattern. The idea that rapid variability induces reversals goes back to Parker (1969) and Levy (1972), and later Olson (1983), Olson and Hagee (1990), Le Mouël et al. (1997), Narteau et al. (2000) (see also Gibbons 1998), all within the mean-field concept; the Poissonian distribution of polarity intervals points the same way (pp. 264–265).

What is new (p. 265): (1) a theoretical framework for computing the statistical properties of mean-field dynamos "of the geo-type" (an αΩ model with a non-periodic dipolar fundamental mode), illustrated on a specific model; (2) a detailed comparison of the predictions with the observed geomagnetic field. The αΩ choice is motivated by evidence that the inner core rotates ≲1° yr^−1 faster than the outer core (Song and Richards 1996; Su et al. 1996; Vidale et al. 2000; Richards 2000; dissenting: Souriau et al. 1997; Laske and Masters 1999) and by the simulations of Glatzmaier and Roberts (1996) (p. 266).

The origin of the fluctuations (p. 266): the average ⟨·⟩ is interpreted as an azimuthal average (Braginskii 1965a,b), i.e. an average over an incomplete ensemble containing only the system and its azimuthally rotated states. Because only a finite number of convective elements lie on a circle, and these are renewed every correlation time τ_c, the averaged parameters α, β and the mean flow v must have a fluctuating component on the time scale τ_c. Fluctuations in α are expected to dominate because u·∇×u can take either sign while |u|² cannot (Otmianowska-Mazur et al. 1997). The same idea explained the phase–amplitude correlation of the solar cycle (Hoyng 1993, 1996; Ossendrijver et al. 1996). Schmitt et al. (2000; = the PEPI 2001 paper of this series) had already applied it to the geodynamo in short form; the present paper is the detailed numerical and theoretical study (p. 266).

## 4. The model

### 4.1 Physical setup and geometry

An axisymmetric mean-field dynamo in a spherical shell representing the outer core, "very similar" to Schmitt and Schüssler (1989) (p. 267). Spherical polar co-ordinates (r, θ, φ). The mean field is split into poloidal and toroidal parts (Eq. 3), and the radial dependence is removed by a local approximation to a radial spherical wave (Eq. 4) with wavenumber k; the radial co-ordinate then becomes a free parameter set equal to the reference radius R (p. 268). The mean flow is a differential rotation v = Ω × r = Ω r sin θ e_φ with constant radial shear dΩ/dr ≡ Ω′ (pp. 267–268). The αΩ approximation is made (the α-term is dropped from the toroidal equation). The remaining problem is one-dimensional in colatitude θ ∈ [0, π] with boundary conditions P = T = 0 at θ = 0, π (p. 268). The mean α-effect is α ∝ cos θ, antisymmetric about the equator (Eq. 9).

**Note:** the model has no inner core and no explicit radial structure; the paper lists a 2-D (r, θ) generalisation with an inner core as the most desirable extension (p. 298).

### 4.2 Governing equations

Mean-field induction equation (Moffatt 1978; Krause and Rädler 1980), p. 265:

$$\frac{\partial\langle\mathbf B\rangle}{\partial t}=\nabla\times\left[\mathbf v\times\;+\;\alpha\;-\;(\beta+\eta)\nabla\times\right]\langle\mathbf B\rangle \qquad \text{(Eq. 1)}$$

⟨B⟩ is the mean (azimuthally averaged) field; v the mean flow; u the turbulent convection (total flow v + u); η the molecular resistivity; β the turbulent resistivity; α the helicity-driven α-effect. With correlation time τ_c (eddy turnover time) and correlation length λ_c (eddy size):

$$\alpha\simeq-\tfrac13\langle\mathbf u\cdot\nabla\times\mathbf u\rangle\tau_c;\qquad \beta\simeq\tfrac13\langle|\mathbf u|^2\rangle\tau_c \qquad \text{(Eq. 2)}$$

Mean-field validity requires the eddy magnetic Reynolds number uλ_c/η ≫ 1, i.e. η ≪ β; η is dropped relative to β in (1) (p. 265).

Field decomposition and radial ansatz (p. 267):

$$\langle\mathbf B\rangle=\nabla\times A\,\mathbf e_\varphi+B\,\mathbf e_\varphi \qquad \text{(Eq. 3)}$$

$$A=\frac{\alpha_0R^2}{\beta}\,P(\theta,t)\,\frac{\cos kr}{kr};\qquad B=T(\theta,t)\,\frac{\cos kr}{kr} \qquad \text{(Eq. 4)}$$

A e_φ is the poloidal vector potential (field ⊥ e_φ), B the toroidal field (∥ e_φ). The factor α_0R²/β makes P and T have the same dimension; β is constant; α_0 and R are reference values of α and of the radial position; cos(kr)/(kr) crudely models radial turbulent diffusion (p. 267).

The 1-D αΩ system (Ossendrijver et al. 1996), pp. 267–268:

$$\frac{\partial\mathbf b}{\partial\tau}=G\mathbf b\quad\text{with}\quad \mathbf b(\theta,\tau)=\begin{pmatrix}P\\T\end{pmatrix} \qquad \text{(Eq. 5)}$$

$$G=\begin{pmatrix}L & \alpha/\alpha_0\\ C\,(\partial/\partial\theta)\sin\theta & L\end{pmatrix} \qquad \text{(Eq. 6)}$$

$$C=\alpha_0R^4\Omega'/\beta^2,\qquad \tau=t/\tau_d,\qquad \tau_d=R^2/\beta \qquad \text{(Eq. 7)}$$

$$L=\frac1{\sin\theta}\frac{\partial}{\partial\theta}\sin\theta\frac{\partial}{\partial\theta}-\frac1{\sin^2\theta}-(kR)^2 \qquad \text{(Eq. 8)}$$

C is the dynamo number, τ the dimensionless time in units of the turbulent diffusion time τ_d, L the diffusion operator (the −(kR)² term is the radial diffusion from the ansatz (4)), Ω′ = dΩ/dr the (constant) differential rotation. Boundary conditions P = T = 0 at θ = 0, π (p. 268). Mean α-effect:

$$\alpha/\alpha_0=\cos\theta \qquad \text{(Eq. 9)}$$

α_0 being the value of α at the North pole (p. 268).

**Note:** in the toroidal equation the operator C(∂/∂θ) sin θ acts on P, i.e. the shear term is C ∂(sin θ P)/∂θ; the upper-right entry α/α_0 multiplies T (α-effect generating poloidal from toroidal field).

Eigenvalue problem (Section 2.1.1, pp. 268–269):

$$G\mathbf b_i=\lambda^i\mathbf b_i;\qquad \mathbf b_i=\begin{pmatrix}P_i\\T_i\end{pmatrix} \qquad \text{(Eq. 10)}$$

$$\tilde G\tilde{\mathbf b}_i=\lambda^{i*}\tilde{\mathbf b}_i;\qquad \tilde{\mathbf b}_i=\begin{pmatrix}\tilde P_i\\ \tilde T_i\end{pmatrix} \qquad \text{(Eq. 11)}$$

$$(\tilde{\mathbf b}_i,\mathbf b_j)=\delta_{ij} \qquad \text{(Eq. 12)}$$

$$(\mathbf b_1,\mathbf b_2)=\int_0^\pi d\theta\,\sin\theta\,(P_1^*P_2+T_1^*T_2) \qquad \text{(Eq. 13)}$$

$$\mathbf b(\theta,\tau)=\sum_k a_k(\tau)\mathbf b_k(\theta) \qquad \text{(Eq. 14)}$$

$$a_i=(\tilde{\mathbf b}_i,\mathbf b) \qquad \text{(Eq. 15)}$$

G is not self-adjoint; eigenvalues λ^i and eigenmodes b_i are in general complex; G̃ is the adjoint operator, b̃_i its eigenfunctions, * complex conjugation; {b_i}, {b̃_i} form a bi-orthonormal set, which makes (15) follow from (14). Details of the adjoint problem: Kleeorin and Ruzmaikin (1984), Ossendrijver (1996), Ossendrijver et al. (1996) (p. 269).

The requirement of a non-periodic dipole fundamental mode restricts 0 < kR < 1 and C > 0. For given kR there is a critical C_0 below which all modes are damped. Two dynamo numbers are distinguished (p. 269):

$$C_0\equiv\alpha_0R^4\Omega'/\beta^2,\qquad C_m\equiv\alpha_mR^4\Omega'/\beta^2 \qquad \text{(Eq. 16)}$$

with α_0 the reference α associated with the critical state and α_m the (larger) unquenched value at which the model operates. The model has three free parameters: kR, C_m and a fluctuation strength. Chosen: kR = 0.5 (so C_0 = 57.9 and 1/k = 2R is the radial length scale), C_m = 100 (p. 269).

Nonlinearity (Section 2.2.1, pp. 269–270): local α-quenching α_m → α_m(1 − const·⟨B⟩²) (Rüdiger 1973) is replaced, for tractability, by global quenching α_m → α_m(1 − const·∫⟨B⟩² d³r) (Brandenburg et al. 1989):

$$\alpha=\alpha_m(1-q)\cos\theta \qquad \text{(Eq. 17)}$$

$$q\simeq\text{const}'\int_0^\pi(\Re T)^2\sin\theta\,d\theta \qquad \text{(Eq. 18)}$$

$$\;\;\simeq\text{const}''(\Re a_0)^2 \qquad \text{(Eq. 19)}$$

(18) uses that the toroidal field dominates the poloidal one in the αΩ approximation and that quenching is done by the real field; (19) keeps only the fundamental mode's contribution (coefficient a_0), usually dominant (p. 270).

Nonlinear equilibrium (Section 2.2.2, pp. 270–271): with quenching function h(q) (= 1 − q here) G has h(q) cos θ as upper-right element and C = C_m lower-left. Since global quenching leaves the angular shape of α intact, the dynamo has a unique nonlinear equilibrium, independent of the shape of h: the field evolves under an effective dynamo number C_m h(q) and the equilibrium coincides with the linear critical solution, C_m h(q_0) = C_0, i.e. α_m h(q_0) = α_0; it is stable if dh/dq < 0. The constants in (18)–(19) only fix the unit of ⟨B⟩; the choice of q matters little as long as it is quadratic in the field (p. 271).

Fluctuations (Section 2.2.3, pp. 271–272):

$$\alpha=\alpha_m(1-q)\cos\theta+\delta\alpha(\theta,\tau) \qquad \text{(Eq. 20)}$$

$$\delta\alpha=\alpha_0\,\frac{f\,F(\theta,\tau)}{\sqrt{2N_c\sin\theta}} \qquad \text{(Eq. 21)}$$

See Section 4.4 below for the definitions of f, F, N_c.

Model equations (Section 2.3, pp. 272–273). Supercriticality parameter:

$$s\equiv\frac{C_m}{C_0}-1=\frac{q_0}{1-q_0} \qquad \text{(Eq. 22)}$$

(since C_m(1 − q_0) = C_0); for kR = 0.5, C_m = 100: s = 0.728, q_0 = 0.421 (p. 272).

$$q=q_0(\Re a_0)^2 \qquad \text{(Eq. 23)}$$

with a_0 now measured in units of its nonlinear-equilibrium value, and

$$a\equiv\Re a_0 \qquad \text{(Eq. 24)}$$

so that (20) becomes

$$\frac{\alpha}{\alpha_0}=\cos\theta+s(1-a^2)\cos\theta+\frac{\delta\alpha}{\alpha_0} \qquad \text{(Eq. 25)}$$

and (5) is replaced by the model equation

$$\frac{\partial\mathbf b}{\partial\tau}=\left[G+s(1-a^2)E+V\right]\mathbf b \qquad \text{(Eq. 26)}$$

$$E=\begin{pmatrix}0&\cos\theta\\0&0\end{pmatrix};\qquad V=\begin{pmatrix}0&\delta\alpha/\alpha_0\\0&0\end{pmatrix} \qquad \text{(Eq. 27)}$$

Here G is the operator (6) at the critical dynamo number C_0 (with α/α_0 = cos θ), and a = ℜ(b̃_0, b) with b̃_0 the adjoint fundamental mode at C_0. The first term is the linear evolution, the second the nonlinearity (and, for small a, the supercriticality), the third the fluctuations. Because the nonlinearity is even in the field there are two equilibria, a = +1 and a = −1 (B → −B invariance) (p. 273). The expansion in eigenfunctions at C_0 is an expansion around the nonlinear equilibrium and therefore efficient (p. 273).

Parity diagnostic (Brandenburg et al. 1989), p. 276:

$$\mathcal P=\frac{\mathcal E_S-\mathcal E_A}{\mathcal E_S+\mathcal E_A} \qquad \text{(Eq. 28)}$$

with ℰ_{S,A} = ∫_0^π (ℜT_{S,A})² sin θ dθ, T_S(θ) = T(θ) + T(π−θ), T_A(θ) = T(θ) − T(π−θ). Except during (aborted) reversals and excursions 𝒫 ≃ −1 (antisymmetric, dipolar). **Note:** in our reproduction of the standard run the median 𝒫 outside the (aborted) reversals and excursions is −0.82, with 𝒫 < −0.9 a third of the time; the symmetric part then carries about 10% of the toroidal energy, which is what the |a₁| ≈ 0.1–0.3 between events in Fig. 2 implies, so "≃ −1" is a loose description of the figure rather than a value. See `dev/DISCREPANCIES.md` DR-7.

The coupled mode equations (29)–(31), Fokker–Planck equation (32)–(41) and reversal statistics (42)–(60) are derived results and are given in Section 6.

### 4.3 Non-dimensionalisation, units and time scales

- Length: R, the reference radius (set equal to the radial position at which the local radial approximation is made); radial scale of the field 1/k = 2R for kR = 0.5 (p. 269). In Section 5.2.3 R ≈ 3500 km (p. 295).
- Time: τ = t/τ_d with τ_d = R²/β, the turbulent diffusion time (Eq. 7). The paper fixes τ_d empirically: a typical model reversal lasts ~3τ_d; the observed duration is ≲10 kyr, hence τ_d ~ 3 kyr (p. 277, p. 293). With R ~ 3500 km this gives β ~ 10^6 cm² s^−1 (p. 295). All simulation times (figures' "t [diffusion times]") and ⟨T_r⟩ are in units of τ_d.
- Field: P and T have the same dimension (factor α_0R²/β in Eq. 4). In the nonlinear model the field unit is such that the fundamental-mode amplitude a_0 = 1 in nonlinear equilibrium (Eq. 23, p. 272); all a_k in Eq. (29) are in this unit (p. 282). The constants in (18)–(19) merely set the unit of ⟨B⟩ (p. 271).
- α: in units of α_0 (Eq. 25). Rates (λ^i, Λ, μ, S, D) are per diffusion time.
- Conversion used by the paper: τ_c = 0.05 τ_d ↔ 150 yr; the oscillation time about equilibrium (2Λ)^−1 ≈ 1 τ_d ≈ 3 kyr; overtone fine structure (ℜλ)^−1 ≲ 1 kyr; ⟨T_r⟩ ≈ 100τ_d ≈ 3×10^5 yr; ⟨T_r⟩ ≈ 500τ_d ≈ 1.5 Myr (run #2); 5×10^4 τ_d of simulation (one Table II run) ≈ 150 Myr (pp. 277, 295–296, 298).

### 4.4 Fluctuation model

Construction (pp. 271–272): helicity fluctuations with correlation length ("cell size") λ_c and correlation time τ_c are homogeneously distributed in the outer core. δα, being an azimuthal average, is proportional to N^{−1/2}, where N = 2N_c sin θ is the number of cells on a circle of constant latitude; the total number of cells on a sphere is 4N_c²/π. F(θ,τ) is a random function of θ for given τ and of τ for given θ, with zero mean and unit variance. In the numerical experiments F is updated once every correlation time τ_c, independently in every latitude interval θ_c = π/N_c. The parameter f measures the unknown relative strength of the fluctuations: in nonlinear equilibrium (α_m(1 − q_0) = α_0) each cell contributes a mean α_0 cos θ and an r.m.s. fluctuation f α_0 per cell (p. 272). The fluctuations are not magnetically quenched in the reported runs (p. 272; reason in Section 3.3, p. 281: quenching merely reduces the effective f and can be offset by raising f).

Scaling: the mean part of α scales as cos θ, the fluctuating part as (sin θ)^{−1/2}, i.e. increasing towards the poles (no problem there because the field vanishes at the poles). Most other details (cell size, distribution) are irrelevant (p. 272). All statistical properties depend only on f²τ_c/N_c², the third model parameter (p. 272; proof in Appendix E, p. 311). Values: τ_c = 0.05, N_c = 10 (p. 272); f = 6.4 standard, also 4.4, 5.4, 7.4, 8.4 (Table II, p. 292).

The "forcing parameter": the paper calls f the forcing parameter (Fig. 1 caption, p. 274; p. 277) and f²τ_c/N_c² the "effective forcing parameter" (p. 294); the latter acts like a temperature in the Kramers analogy (p. 290). For the standard run πf²τ_c/2N_c² = 3.22×10^−2 (p. 291) and πf²τ_c²/2N_c² ≃ 1.6×10^−3 (p. 310).

Correlation function implied by the model (Appendix B.2, Eq. B.6–B.7, p. 308):

$$\left\langle\frac{\delta\alpha(\theta,\tau)}{\alpha_0}\frac{\delta\alpha(\theta-\zeta,\tau-\sigma)}{\alpha_0}\right\rangle=\frac{f^2\langle F(\theta,\tau)F(\theta-\zeta,\tau-\sigma)\rangle}{2N_c\sqrt{\sin\theta\sin(\theta-\zeta)}}=\frac{f^2(1-|\zeta|/\theta_c)(1-|\sigma|/\tau_c)}{2N_c\sqrt{\sin\theta\sin(\theta-\zeta)}} \qquad \text{(Eq. B.6)}$$

(zero for |ζ| > θ_c or |σ| > τ_c), approximated when τ_c and θ_c are small compared with the relevant scales by

$$\simeq\frac{f^2\tau_c\theta_c}{2N_c\sin\theta}\,\delta(\zeta)\,\delta(\sigma) \qquad \text{(Eq. B.7)}$$

Magnitude: in equilibrium the ratio of fluctuating to systematic term in (20) is δα_rms/(α_0 cos θ) = (α_0 f/√(2N_c sin θ))/(α_0 cos θ) ≳ 2.5 for the standard run (p. 296, where the relation symbol between the expression and "2.5" is missing in print; p. 299 prints "≳ 2.5"). Roughly the same number explains the solar phase–amplitude correlation (Ossendrijver et al. 1996), "probably a coincidence" (p. 299).

Physical numbers (Section 5.2.4, p. 296): τ_c = 0.05 ↔ 150 yr, eddy size λ_c ~ uτ_c ~ 500 km with convective speed u ~ 10^−2 cm s^−1; N_c = 10 ↔ 20 cells on a great circle, 4N_c²/π ~ 130 on a sphere. The value α_m ~ −3×10^−3 cm s^−1 is called reasonable because |α| is then everywhere smaller than the convective speed u ~ 10^−2 cm s^−1; footnote 3 (p. 296) adds that the maximum possible value of |α| is the r.m.s. convective flow speed, attained for maximally helical convection.

### 4.5 Nonlinearity / saturation

Global α-quenching, Eqs. (17), (19), (23): α_m → α_m(1 − q) with q = q_0 a², a = ℜa_0 in equilibrium units. The effective dynamo number is C_m(1 − q); equilibrium at C_m(1 − q_0) = C_0 (p. 271). In the model equation (26) the nonlinearity appears as s(1 − a²)E: for |a| < 1 the dynamo is supercritical (effective growth), for |a| > 1 subcritical. Alternatives tested with no material change in the dipole statistics (Section 3.3, p. 277): (1) α as in (17) with q from (18) (full toroidal energy instead of only a_0); (2) α = α_m cos θ/(1 + q) with q from (19); (3) local instead of global quenching; also changes of C_m and kR, and magnetic quenching of δα. Magnetic quenching of the fluctuations reduces the reversal rate and the r.m.s. variability, because it lowers the effective f; this can be offset by raising f (p. 281). Hence the "standard case defined by (26)" without quenching of δα.

The paper stresses that during a reversal 1 − a² ≃ 1 in (26)/(29), so the dynamo then evolves as a linear system; it argues this may hold for the geodynamo too (field down by a factor ≥ 4, Lorentz force 10–20 times smaller) (p. 300).

### 4.6 Parameter values used

| Symbol | Meaning | Value(s) | Where |
|---|---|---|---|
| kR | radial wavenumber × reference radius | 0.5 (1/k = 2R) | p. 269; footnote 1 p. 273 |
| C_0 | critical dynamo number (for kR = 0.5) | 57.9 | p. 269; Table I p. 270 |
| C_m | operating dynamo number | 100 | p. 269; footnote 1 p. 273 |
| s | supercriticality C_m/C_0 − 1 | 0.728 | p. 272 |
| q_0 | equilibrium quenching | 0.421 | p. 272 |
| α/α_0 (mean) | angular profile | cos θ | Eq. 9 p. 268 |
| τ_c | correlation time of δα (units τ_d) | 0.05 (= 5×10^−2) | p. 272; footnote 1 p. 273 |
| N_c | number of latitude cells on [0, π] (θ_c = π/N_c) | 10 | p. 272 |
| f | fluctuation strength (r.m.s. fα_0 per cell) | 6.4 (standard); 4.4, 5.4, 7.4, 8.4 (Table II) | Fig. 1 p. 274; Table II p. 292 |
| πf²τ_c/2N_c² | diffusion prefactor, standard run | 3.22×10^−2 | p. 291 |
| grid | equidistant θ-grid on [0, π] | 60 points | footnote 1 p. 273 |
| Δτ | time step | 10^−2 | footnote 1 p. 273 |
| run length | per run in Table II | 5×10^4 τ_d | Table II p. 292 |
| E_00 | overlap integral ∫ sin θ cos θ P̃_0 T_0 dθ | 1.048 | p. 282 |
| ∫ sin θ (P̃_0T_0)² dθ | integral in Eq. 36 | 3.529 | p. 284; p. 310 |
| μ | sE_00 | 0.763 | Eq. 34 p. 283 |
| Λ | λ^0 at C = C_m (growth rate) | 0.527 | Eq. 34 p. 283; Table I |
| λ^i (C_0) | eigenvalues i = 0..3 | 0; −2.76+4.41i; −8.67+4.97i; −15.9+6.13i | Table I p. 270 |
| λ^i (C_m) | eigenvalues i = 0..3 | 0.527; −2.16+5.88i; −8.21+6.04i; −15.7+7.70i | Table I p. 270 |
| D_1 fit | standard run | (πf²τ_c/2N_c²)(1.36 + 2.00a²), valid |a| ≤ 1.1 | Eq. 61 p. 291; Fig. 8 p. 292 |
| c², γ | standard run | 0.246, 2.96 | Eq. 63 p. 291 |
| τ_d | diffusion time (physical) | ~3 kyr | p. 277, 293 |
| R | reference radius | ~3500 km | p. 295 |
| β | turbulent resistivity | ~10^6 cm² s^−1 | p. 295 |
| η | molecular resistivity (σ = 6×10^5 S m^−1) | ~1.3×10^4 cm² s^−1 | p. 296 |
| ΔR | outer minus inner core radius | 2300 km | p. 296 |
| Ω′ | shear from inner-core super-rotation ≲1° yr^−1 | ~ −1° yr^−1 / ΔR | p. 296 |
| C_Ω, C_α | R³Ω′/β, α_mR/β | ~ −100, ~ −1 | p. 296 |
| α_m | unquenched reference α (|α| < u everywhere; its maximum possible value is the r.m.s. convective speed, footnote 3) | ~ −3×10^−3 cm s^−1 | p. 296 |
| u | convective speed | ~10^−2 cm s^−1 | p. 296 |
| λ_c | eddy size ~ uτ_c | ~500 km | p. 296 |

## 5. Methods

**Numerical scheme (Section 3, footnote 1, p. 273).** Equation (26) is solved on a 60-point equidistant θ-grid on [0, π] with time step Δτ = 10^−2, kR = 0.5, C_m = 100, τ_c = 5×10^−2, N_c = 10. The time-stepping algorithm and the spatial discretisation of L are not stated. Runs extend over 5×10^4 diffusion times each (Table II, p. 292); simulations "extending over ~10^5 turbulent diffusion times or more are quite feasible" (p. 273), and the discussion mentions runs "up to ~150 Myr" (p. 298). Plotted data in Fig. 1 have a spacing of 100 time steps (= 1τ_d); Fig. 2 uses 5 time steps; Figs. 3–4 are at full resolution (pp. 274–279). The same random series is used for the three values of f in Figs. 1 and 5 (p. 280). Mode amplitudes a_k are extracted from the simulated b via the projection (15) with the adjoint eigenfunctions at C_0 (Fig. 2 caption, p. 275).

**Eigenmode expansion.** Eigenfunctions and eigenvalues of G at C_0 (and C_m) are computed, including the adjoint problem (Kleeorin and Ruzmaikin 1984; Ossendrijver 1996; Ossendrijver et al. 1996), bi-orthonormalised according to (12)–(13). Modes are ordered by decreasing growth rate ℜλ^i, which is "a decreasing series, asymptotically ~ i(i+1)" (printed without a minus sign, i.e. −ℜλ^i grows like i(i+1)) because the eigenvalue problem is parabolic (Table I footnote b, p. 270). Mode 0 is the non-periodic dipole (T antisymmetric about the equator); all overtones i ≥ 1 are periodic (ℑλ^i ≠ 0; footnote c misprints this as "for i ≥ 0") and alternate in parity (Table I).

**Statistical analysis (Section 4, Appendices A–F).** (i) Projection of (26) onto the modes gives the coupled stochastic mode equations (29). (ii) The complex system is rewritten in real "supervectors" (Appendix A.1). (iii) Van Kampen's (1976, 1992) cumulant/short-correlation-time expansion for multiplicative noise yields a Fokker–Planck equation for the joint distribution Π of all mode amplitudes (A.9–A.15), under the conditions |𝒱|τ_c ≪ 1 (A.10, verified in Appendix C) and "no appreciable evolution during a correlation time" (exp(±s𝒜) ≈ 1). (iv) Integration over all overtone amplitudes gives the one-variable Fokker–Planck equation (32) for p(a) with drift and diffusion coefficients expressed as conditional averages (A.17–A.19). (v) These are evaluated in Appendix B using the δ-correlated model (B.7) and the identity M_μννσ = 0 (Appendix D, completeness of eigenfunctions). (vi) Stationary solution, parabolic approximation around the peaks, and the mean-first-passage-time formula (Gardiner 1990; Van Kampen 1992) give p(a), the relative variance, and ⟨T_r⟩ (Appendix F).

**Fitting procedure for D_1 (Section 5.1, p. 291).** Since c², γ cannot yet be expressed in the model parameters (Appendix E), D_1 is measured from the standard run: at every diffusion time (to ensure independent draws) the instantaneous value of (37) is evaluated with T_+ known and the averaging omitted, giving 5×10^4 points of D_1/F (F = πf²τ_c/2N_c²) versus |a| (Fig. 8). A histogram of bin averages is least-squares fitted by a quadratic in a, weighted by the number of points per bin, for a ≤ 1.1 (for larger a an a⁴ term becomes appreciable) (Fig. 8 caption, p. 292). g_1(a) was measured the same way to test (33) (p. 292).

**Observational data.** (1) Three VADM records for the relative variance of the dipole between reversals: the last 130 kyr (Fig. 4.9 of Merrill et al. 1996), the Matuyama chron 2–2.5 Myr ago and the Gauss chron 2.5–3 Myr ago (both from Valet and Meynadier 1993); a is taken proportional to the VADM ordinate since the poloidal dipole field is ∝ a (p. 294). (2) The Sint-800 record (Guyodo and Valet 1999): 800 kyr of VADM, 800 points with 1000-yr spacing; its amplitude distribution is binned in 35 bins; data older than 780 kyr (the last reversal) are folded to positive values using p(−a) = p(a) (Fig. 10 caption, p. 295). The comparison is explicitly "preliminary"; no retuning of kR, C_m was attempted (p. 295).

## 6. Results

### 6.1 Simulated behaviour of the standard case (Section 3.1–3.2, pp. 273–277)

Figures 1–2 (f = 6.4) share with the geomagnetic record (Jacobs 1994; Merrill et al. 1996): (1) a dominantly dipolar field of rapidly varying amplitude with fine structure (vertical stripes in the T(θ,τ) plot) from transient excitation of overtones; a fluctuates around its equilibrium value; (2) occasional sudden reversals, with f tuned so that the mean time between reversals is about 100τ_d; (3) aborted reversals (a changes sign briefly, then returns) and excursions (pp. 275–276). Definitions (footnote 2, p. 276): a reversal is an odd number of zero crossings of a within a few τ_d; an aborted reversal an even number in rapid succession; an excursion has a "approaching zero" without sign change; the authors call this subjective and note ambiguities (three closely spaced crossings could be an aborted event plus a reversal).

Phenomenology of (aborted) reversals from Figs. 3–4 (pp. 276–277):
1. During the event the field tends to become symmetric (𝒫 ↑ +1) because overtones dominate; mode #1 (quadrupole) frequently stands out; quasi-periodic behaviour with period of order τ_d can appear because all overtones are periodic (mode #1 alone would give 2π/5.88 ≃ 1.1τ_d, Table I).
2. Onset: a decrease of |a| combined with a sudden enhancement (≲0.5τ_d) of modes #1 and/or #2, or rather of the overtone toroidal field T_+ (Eq. 38); but the realisation of δα must be "favorable" — not all jumps in T_+ produce a reversal.
3. Just after the zero crossing the whole field is often small (all a_k small) and enters a recovery phase: with |a| ≪ 1 the nonlinearity is off and (26) gives roughly exponential growth, time scale ~Λ^−1 (Fig. 3 caption, p. 278).
4. Fast events last ~0.5τ_d; typical reversal duration ~3τ_d. With the observed ≲10 kyr: τ_d ~ 3 kyr, so 100τ_d ~ 3×10^5 yr, about the observed recent mean interval (p. 277).

### 6.2 Scaling with the forcing (Section 3.3, pp. 277–281)

The reversal rate is extremely sensitive to f (Fig. 5, p. 280: for f = 5.4 the 2000τ_d window shows 6 reversal bars, at τ ≈ 540, 690, 1180, 1660 and a short reversed interval ≈ 1760–1790; for f = 7.4 several tens — counts read off the figure, not stated in the text). The dominant scaling of mean properties (mean time between reversals, r.m.s. variability between reversals) is through f²τ_c/N_c², and the reversal rate depends exponentially on it (p. 277). Other nonlinearities, other C_m and kR, and quenching of δα do not materially change the dipole statistics (p. 277, 281).

### 6.3 Coupled mode equations (Section 4.1, p. 282)

Taking d/dτ of (15) with (26) and (14):

$$\frac{da_i}{d\tau}=\lambda^ia_i+s(1-a^2)E_{ik}a_k+V_{ik}a_k \qquad \text{(Eq. 29)}$$

with a = ℜa_0; summation convention on repeated lower indices (in E_ik a_k and V_ik a_k, not in λ^i a_i); all a_k in units of the equilibrium a_0.

$$E_{ik}=(\tilde{\mathbf b}_i,E\mathbf b_k)=\int_0^\pi d\theta\,\sin\theta\cos\theta\,\tilde P_i^*T_k \qquad \text{(Eq. 30)}$$

$$V_{ik}(\tau)=(\tilde{\mathbf b}_i,V\mathbf b_k)=\int_0^\pi d\theta\,\sin\theta\,(\delta\alpha/\alpha_0)\,\tilde P_i^*T_k \qquad \text{(Eq. 31)}$$

E_00 = ∫ sin θ cos θ P̃_0T_0 dθ = 1.048 plays a special role. The λ^i are those at C_0 (Table I). The second term drives the system to the equilibrium a_0 = ±1, a_i = 0 (i ≥ 1); the third is multiplicative noise (p. 282). Figure 6 (p. 283) shows ℜa_1 versus a = ℜa_0 for f = 5.4: the overtone fluctuates rapidly while a stays tied to ±1 with occasional jumps.

### 6.4 Fokker–Planck equation and its coefficients (Section 4.2, pp. 282–284)

$$\frac{\partial p}{\partial\tau}=-\frac{\partial}{\partial a}Sp+\frac12\frac{\partial^2}{\partial a^2}Dp \qquad \text{(Eq. 32)}$$

p(a, τ) is the probability density of a = ℜa_0; S the drift, D the diffusion coefficient (Risken 1984; Gardiner 1990; Van Kampen 1992).

$$S=\mu(1-a^2)a+g_1(a)\simeq\Lambda(1-a^2)a \qquad \text{(Eq. 33)}$$

$$\mu\equiv sE_{00}=0.763,\qquad \Lambda\equiv\lambda^0\big|_{C=C_m}=0.527 \qquad \text{(Eq. 34)}$$

g_1(a) is the overtone contribution to the drift (Eq. B.3); the second form of (33) is a small-a approximation, accurate for |a| ≲ 0.4 (p. 292; Appendix B.1).

$$D=D_0a^2+D_1(a) \qquad \text{(Eq. 35)}$$

$$D_0=\frac{\pi f^2\tau_c}{2N_c^2}\int_0^\pi d\theta\,\sin\theta\,(\tilde P_0T_0)^2 \qquad \text{(Eq. 36)}$$

$$D_1=\frac{\pi f^2\tau_c}{2N_c^2}\int_0^\pi d\theta\,\sin\theta\,\tilde P_0^2\,\langle T_+^2\rangle\big|_a \qquad \text{(Eq. 37)}$$

The integral in (36) equals 3.529 (p. 284). T_+ is the toroidal field of the overtones,

$$T_+=\Re\sum_{k\ge1}a_kT_k=\Re(T-a_0T_0) \qquad \text{(Eq. 38)}$$

and ⟨·⟩|_a is a conditional average at fixed a. D_0a² is the direct effect of δα on a (term V_00a_0 in (29)); D_1(a) is the indirect effect through stochastic forcing of the overtones (terms V_0k a_k, k ≥ 1). Expansion to order a²:

$$D\simeq\Delta_0a^2+\Delta_1 \qquad \text{(Eq. 39)}$$

$$\Delta_0=D_0+\tfrac12D_1''(0);\qquad \Delta_1=D_1(0) \qquad \text{(Eq. 40)}$$

using D_1′(0) = 0 because D_1(−a) = D_1(a) (p. 284). For the standard run D_1 > D_0a² for |a| ≲ 1 (Fig. 7 caption, p. 285).

### 6.5 Particle-in-well analogy and amplitude distribution (Section 4.3, pp. 284–286)

The only remarkable feature is that D depends on the total non-dipole toroidal field, so all dependence on individual overtones has disappeared (p. 284; Appendix E). Equation (32) also governs the position of a heavily damped, randomly forced particle in the bistable potential defined by −∂U/∂a = S:

$$U\simeq\tfrac14\Lambda(1-a^2)^2 \qquad \text{(Eq. 41)}$$

The particle sits near a = ±1 and performs a random motion (secular variation) with occasional jumps (reversals): "the geodynamo is thus pictured as a bistable oscillator". The central hill at a = 0 is due to the supercriticality, the outer walls to the nonlinearity (p. 285). Stationary solution:

$$p(a)\propto\frac1D\exp\left(\int_0^a\frac{2S}{D}\,da'\right) \qquad \text{(Eq. 42)}$$

$$\;\;\propto\left(1+\frac{a^2}{c^2}\right)^{\gamma c^2+\gamma-1}\exp(-\gamma a^2) \qquad \text{(Eq. 43)}$$

$$\gamma=\Lambda/\Delta_0;\qquad c^2=\Delta_1/\Delta_0 \qquad \text{(Eq. 44)}$$

γ and c² are dimensionless measures of U(0) and D(0). A steeply rising power law times a sharply decreasing exponential gives two peaks near a = ±1 (Fig. 7). Since c² ≈ D_1(0)/D_0 (p. 286):

$$c^2\propto\left(\frac{\text{r.m.s. toroidal field at a reversal}}{\text{r.m.s. toroidal field between reversals}}\right)^2 \qquad \text{(Eq. 45)}$$

**Note:** (43) follows from (42) with S = Λ(1−a²)a and D = Δ_0a² + Δ_1: 2S/D = 2Λ(1−a²)a/(Δ_0(a²+c²)) integrates to −γa² + (γc²+γ)ln(1+a²/c²), and the 1/D prefactor supplies the extra −1 in the exponent.

### 6.6 The reversal mechanism (Section 4.4, pp. 286–287)

Three steps: (1) Without overtones (a_k = 0 for k ≥ 1; λ^0 = 0), (29) reduces to

$$\frac{da_0}{d\tau}=s(1-a^2)E_{00}a_0+V_{00}(\tau)a_0 \qquad \text{(Eq. 46)}$$

The noise is multiplicative and vanishes for small a_0; since E_00 and V_00 are real,

$$\frac{d\ln a}{d\tau}=s(1-a^2)E_{00}+V_{00}(\tau) \qquad \text{(Eq. 47)}$$

so ln a performs a restricted random walk and a cannot change sign: no reversals. (2) If a goes rapidly to zero through favourable fluctuations, the damped overtones do not follow immediately, and

$$\frac{da}{d\tau}=\Re\sum_{k\ge1}(sE_{0k}+V_{0k})a_k \qquad \text{(Eq. 48)}$$

a is then temporarily driven by additive noise and may reverse. (3) A reversal cannot be slow: if the approach is slow the overtones decay, the system falls back to (47), and the result is an excursion ("the flip becomes a flop").

In Fokker–Planck terms: starting from p(a,0) = δ(a−1), a quasi-equilibrium forms in the well on a time scale Λ^−1, then probability creeps diffusively over the hill, filling the other well, towards (43). Writing (32) as ∂p/∂τ = −∂F/∂a, the probability current is F = Sp − ½∂(Dp)/∂a; since S and dD/da vanish at a = 0,

$$F=\tfrac12D(0)\frac{\partial p}{\partial a}\Big|_{a=0}=\tfrac12D_1(0)\frac{\partial p}{\partial a}\Big|_{a=0} \qquad \text{(Eq. 49)}$$

D_0a² contributes nothing: reversals are caused by the small but finite D_1(0), set by ⟨T_+²⟩|_{a=0}, the mean-square non-dipole toroidal field during reversals. Since this favours no overtone, it will in the long run be impossible to attribute reversals to any particular overtone (p. 287).

### 6.7 Statistics of field reversals (Section 4.5, pp. 288–290)

Parabolic approximation to p(a) (Van Kampen 1992), valid for sufficiently large γ:

$$p(a)\propto\exp[\psi(a)]\simeq\exp\left[\psi_m+\tfrac12\psi_m''(a-a_m)^2\right] \qquad \text{(Eq. 50)}$$

$$\psi(a)=-\gamma a^2+(\gamma c^2+\gamma-1)\ln(1+a^2/c^2) \qquad \text{(Eq. 51)}$$

$$a_m=\pm\left(1-\frac1\gamma\right)^{1/2} \qquad \text{(Eq. 52)}$$

$$\psi_m=1-\gamma+(\gamma c^2+\gamma-1)\ln\left(\frac{\gamma c^2+\gamma-1}{\gamma c^2}\right) \qquad \text{(Eq. 53)}$$

$$\psi_m''=-\frac{4\gamma(\gamma-1)}{\gamma c^2+\gamma-1} \qquad \text{(Eq. 54)}$$

γ > 1 is required, otherwise p(a) has a single peak at a = 0 (p. 288).

Variability between reversals (4.5.1): integrating over a single peak (particle in the well at a = 1),

$$\langle a\rangle=\frac{\int_{-\infty}^\infty a\,p(a)\,da}{\int_{-\infty}^\infty p(a)\,da}\simeq a_m=\left(1-\frac1\gamma\right)^{1/2} \qquad \text{(Eq. 55)}$$

The mean of |a| is smaller than 1 (asymmetry of the well and D increasing with a: the random motion is slower for |a| < 1) (pp. 288–289).

$$\frac{\langle(a-\langle a\rangle)^2\rangle}{\langle a\rangle^2}\simeq\frac{\gamma c^2+\gamma-1}{4(\gamma-1)^2} \qquad \text{(Eq. 56)}$$

independent of the units of a, hence directly comparable to the measured relative variability of the dipole between reversals (p. 289).

Mean reversal rate (4.5.2): the mean escape time ⟨T_e⟩ from 0 ≤ a < ∞ starting at a′ > 0 (Gardiner 1990; Van Kampen 1992),

$$\langle T_e\rangle=2\int_0^{a'}\frac{dy}{D(y)p(y)}\int_y^\infty p(z)\,dz \qquad \text{(Eq. 57)}$$

almost independent of a′ near the well bottom (a′ ≈ 1). At escape the particle is at a = 0; naively it falls with probability ½ into either well, so the mean waiting time for a reversal would be 2⟨T_e⟩ = ⟨T_r⟩ (Van Kampen 1992, exercise p. 47); the probabilities are in fact unequal (Appendix F), so

$$\langle T_r\rangle\approx2\langle T_e\rangle\simeq\frac\pi\Lambda\left(\frac2{c^2}\right)^{1/2}\left(\frac{\gamma c^2+\gamma-1}{\gamma-1}\right)^{1/2}\exp(\psi_m) \qquad \text{(Eq. 58)}$$

$$\;\;\sim\frac{\pi\sqrt2}{\Lambda}\left(\frac{c^2}{1+c^2}\right)^{1/2}\exp(K\gamma) \qquad \text{(Eq. 59)}$$

$$K=-1+(1+c^2)\ln\left(\frac{1+c^2}{c^2}\right) \qquad \text{(Eq. 60)}$$

(59) is the large-γ asymptotic expansion of (58); it shows ⟨T_r⟩ ∝ exp[const·U(0)/(f²τ_c/N_c²)] (with Appendix E), explaining the very strong f-dependence. Kramers (1940) analogy: a thermal particle with energy ~kT must acquire U(0) − U(1) = U(0) to escape, with rate ∝ exp[−U(0)/kT]; f²τ_c/N_c² acts as effective temperature (p. 290).

Distribution of polarity intervals (4.5.3): by standard theory the lengths T_r of polarity intervals are Poisson distributed, ∝ exp(−T_r/⟨T_r⟩) for T_r ≫ Λ^−1; verified from a number of long runs, no details given (p. 290).

### 6.8 Internal consistency (Section 5.1, pp. 291–293)

Measured D_1 for the standard case (Fig. 8):

$$D_1\simeq\frac{\pi f^2\tau_c}{2N_c^2}(1.36+2.00a^2) \qquad \text{(Eq. 61)}$$

$$D\simeq\frac{\pi f^2\tau_c}{2N_c^2}(1.36+5.52a^2) \qquad \text{(Eq. 62)}$$

(adding D_0a², i.e. 5.52 ≃ 2.00 + 3.529). With πf²τ_c/2N_c² = 3.22×10^−2:

$$c^2\simeq0.246;\qquad \gamma\simeq2.96 \qquad \text{(Eq. 63)}$$

Then (56) and (58) give ⟨(a−⟨a⟩)²⟩/⟨a⟩² ≃ 0.17 and ⟨T_r⟩ ≃ 94; the simulation had 468 reversals in 5×10^4 τ_d, i.e. ⟨T_r⟩ = 107 ± 5, and a measured relative variance 0.16 (p. 291). Table II (p. 292):

| run | f | c² | γ | variance (56) | variance from run | ⟨T_r⟩ (58) | ⟨T_r⟩ from run |
|---|---|---|---|---|---|---|---|
| 1 | 4.4 | 0.318 | 9.04 | 0.04 | 0.03 | 12000 | 10700 ± 5300 |
| 2 | 5.4 | 0.299 | 5.15 | 0.08 | 0.07 | 480 | 519 ± 53 |
| 3 (standard) | 6.4 | 0.246 | 2.96 | 0.17 | 0.16 | 94 | 107 ± 5 |
| 4 | 7.4 | 0.240 | 1.82 | 0.47 | 0.30 | 35 | 37.6 ± 1 |
| 5 | 8.4 | 0.228 | 1.09 | 1.0 | 0.50 | 35 | 18.7 ± 0.4 |

All runs: N_c = 10, τ_c = 5×10^−2, 5×10^4 diffusion times. Run 5 is outside the validity range of the asymptotic analysis because γ is close to 1. The agreement is called satisfactory given that (58) is only approximate (p. 291). **Note:** Eq. (56) evaluated with the printed c² = 0.228, γ = 1.09 of run 5 gives 10.4, not the 1.0 printed in the table (the text layer confirms the print); rows 1–4 and all five ⟨T_r⟩ values reproduce from (56) and (58) to the printed precision. Footnote b already excludes run 5 from the parabolic approximation. See `dev/DISCREPANCIES.md` DR-1. Appendix E's crude scaling γ ~ (f²τ_c/N_c²)^−1, c² ~ const (Eq. E.1) is "not too bad as a first approximation" (p. 312); indeed γ·f² ≈ 175, 150, 121, 100, 77 across runs 1–5 (**Note:** computed from Table II; the drift of this product shows ⟨T_+²⟩|_{a=0} does depend on the forcing).

Figure 9 (p. 293): the simulated p(a) versus (43) agrees less well; the discrepancy is traced to the approximate drift (33), which is correct only for |a| ≲ 0.4; including the next term in g_1(a) (Appendix B.1; Hoyng et al. 2001, the PEPI 2002 paper) virtually removes it and has little effect on ⟨T_r⟩ (p. 293).

### 6.9 Derivation details (Appendices A–F, pp. 303–314)

**Appendix A (Fokker–Planck derivation).** A.1 Notation: a complex system ȧ_i = B_ik a_k (A.1) is written with a real supervector a = (ℜa_0, ℜa_1, …, ℑa_0, ℑa_1, …)^T and supermatrix

$$\mathbf B=\begin{pmatrix}\Re B & -\Im B\\ \Im B & \Re B\end{pmatrix} \qquad \text{(Eq. A.2)}$$

so that ȧ_μ = B_μν a_ν (A.3); Greek indices, boldface for 2∞-dimensional objects. The mode equations become

$$\dot{\mathbf a}_\mu=\left[\mathbf L_{\mu\nu}+(1-\mathbf a_0^2)\mathbf E_{\mu\nu}+\mathbf V_{\mu\nu}\right]\mathbf a_\nu \qquad \text{(Eq. A.4)}$$

(printed without the factor s in front of (1 − a_0²); (A.7) and (A.12) include s), with L the block matrix of diag(ℜλ^i) and ∓diag(ℑλ^i) (A.5), ℜλ^0 = ℑλ^0 = 0. A.2: the density f(a,τ) of copy systems obeys the continuity equation ∂f/∂τ = −(∂/∂a_μ)ȧ_μ f = [𝒜 + 𝒱(t)]f (A.6) with 𝒜 = −(∂/∂a_μ)[L_μν + s(1−a_0²)E_μν]a_ν (A.7), 𝒱 = −(∂/∂a_μ)V_μν a_ν (A.8). For Π(a,τ) ≡ ⟨f⟩ (Van Kampen 1976, 1992):

$$\frac{\partial\Pi}{\partial\tau}=\left[\mathcal A+\int_0^\infty ds\,\langle\mathcal V(\tau)\exp(s\mathcal A)\mathcal V(\tau-s)\rangle\exp(-s\mathcal A)\right]\Pi \qquad \text{(Eq. A.9)}$$

valid if |𝒱|τ_c ≪ 1 (A.10). Approximating exp(±s𝒜) ≈ 1 (no appreciable evolution during τ_c; justified for low-order modes, effectively assuming high-order modes contribute little, p. 305) gives

$$\frac{\partial\Pi}{\partial\tau}=-\frac{\partial}{\partial\mathbf a_\mu}\left(\mathbf N_{\mu\nu}\mathbf a_\nu-\mathbf M_{\mu\nu\rho\sigma}\mathbf a_\nu\frac{\partial}{\partial\mathbf a_\rho}\mathbf a_\sigma\right)\Pi \qquad \text{(Eq. A.11)}$$

$$\mathbf N=\mathbf L+s(1-\mathbf a_0^2)\mathbf E \qquad \text{(Eq. A.12)}$$

$$\mathbf M_{\mu\nu\rho\sigma}=\int_0^\infty ds\,\langle\mathbf V_{\mu\nu}(\tau)\mathbf V_{\rho\sigma}(\tau-s)\rangle \qquad \text{(Eq. A.13)}$$

(stationary fluctuations, M independent of τ). With the identity a_ν ∂/∂a_ρ = (∂/∂a_ρ)a_ν − δ_νρ (A.14) and M_μννσ = 0 (Appendix D):

$$\frac{\partial\Pi}{\partial\tau}=-\frac{\partial}{\partial\mathbf a_\mu}\left(\mathbf N_{\mu\nu}\mathbf a_\nu-\mathbf M_{\mu\nu\rho\sigma}\frac{\partial}{\partial\mathbf a_\rho}\mathbf a_\nu\mathbf a_\sigma\right)\Pi \qquad \text{(Eq. A.15)}$$

A.3: integrating over all a_κ̄ (κ ≠ 0), p(a_0,τ) = ∏_κ̄ ∫ da_κ̄ Π (A.16), yields (32) with

$$S=\mathbf N_{00}a+\mathbf N_{0\bar\nu}\langle\mathbf a_{\bar\nu}\rangle|_a \qquad \text{(Eq. A.17)}$$

$$D=2\left(\mathbf M_{0000}a^2+\mathbf M_{0\bar\nu0\bar\sigma}\langle\mathbf a_{\bar\nu}\mathbf a_{\bar\sigma}\rangle|_a+\mathbf M_{000\bar\sigma}a\langle\mathbf a_{\bar\sigma}\rangle|_a+\mathbf M_{0\bar\nu00}a\langle\mathbf a_{\bar\nu}\rangle|_a\right) \qquad \text{(Eq. A.18)}$$

with conditional averages ⟨a_ν̄ a_σ̄⟩|_a = ∏∫da_κ̄ a_ν̄ a_σ̄ Π / ∏∫da_κ̄ Π (A.19). Computing these from the stationary Π gave no useful results (p. 307).

**Appendix B (drift and diffusion).** B.1: N_00 = L_00 + s(1−a_0²)E_00 = μ(1−a²) (λ^0 = 0, μ = sE_00). N_0ν̄ a_ν̄ = s(1−a_0²)ℜΣ_{k≥1}E_0k a_k (B.1), using L_0ν̄ = 0 and ℑE_00 = 0. With T_+ ≡ ℜΣ_{k≥1}a_kT_k = ℜ(T − a_0T_0) (B.2),

$$g_1(a)=s(1-a^2)\int_0^\pi d\theta\,\sin\theta\cos\theta\,\tilde P_0\langle T_+\rangle|_a \qquad \text{(Eq. B.3)}$$

⟨T_+⟩|_a vanishes at a = 0 by symmetry and is ∝ a for small a, so g_1(a)/(1−a²) = Aa(1 + Ba² + …); lowest order g_1 = A(1−a²)a; the total S = const·(1−a²)a must have const = λ^0 at C = C_m, which proves (33) for small a (p. 307). **Note:** hence A = Λ − μ = −0.236 for the standard parameters; the overtones reduce the effective restoring rate from μ = 0.763 to Λ = 0.527.

B.2: D_0 = 2M_0000 with M_0000 = ∫_0^∞ ds⟨V_00(τ)V_00(τ−s)⟩ (B.4) = ∫ds ∫∫dθdψ sin θ sin ψ P̃_0(θ)T_0(θ)P̃_0(ψ)T_0(ψ)⟨(δα(θ,τ)/α_0)(δα(ψ,τ−s)/α_0)⟩ (B.5); with (B.7),

$$\mathbf M_{0000}=\frac{\pi f^2\tau_c}{2N_c^2}\int_0^\infty ds\,\delta(s)\int_0^\pi\!\!\int_0^\pi d\theta\,d\psi\,\sin\psi\,\tilde P_0(\theta)T_0(\theta)\tilde P_0(\psi)T_0(\psi)\,\delta(\theta-\psi) \qquad \text{(Eq. B.8)}$$

giving (36) since ∫_0^∞δ(s)ds = ½ (p. 309). Next, M_0ν̄0σ̄⟨a_ν̄a_σ̄⟩|_a = ∫ds⟨V_0ν̄(τ)V_0σ̄(τ−s)⟩⟨a_ν̄a_σ̄⟩|_a (B.9); V_0ν̄a_ν̄ = ℜΣ_{k≥1}V_0k a_k (B.10, ℑV_00 = 0) = ∫dθ sin θ P̃_0 T_+ δα/α_0 (B.11); hence

$$\mathbf M_{0\bar\nu0\bar\sigma}\langle\mathbf a_{\bar\nu}\mathbf a_{\bar\sigma}\rangle|_a=\frac{\pi f^2\tau_c}{2N_c^2}\int_0^\infty ds\,\delta(s)\int\!\!\int d\theta\,d\psi\,\sin\psi\,\tilde P_0(\theta)\tilde P_0(\psi)\langle T_+(\theta)T_+(\psi)\rangle|_a\,\delta(\theta-\psi) \qquad \text{(Eq. B.12)}$$

which leads to (37) (factor 2), assuming low-order modes dominate T_+ so that a δ-peaked θ-correlation is justified. The remaining term g_2(a) = 2a(M_000σ̄⟨a_σ̄⟩|_a + M_0ν̄00⟨a_ν̄⟩|_a), both parts equal,

$$g_2(a)=2a\frac{\pi f^2\tau_c}{2N_c^2}\int_0^\pi d\theta\,\sin\theta\,\tilde P_0^2T_0\langle T_+\rangle|_a \qquad \text{(Eq. B.13)}$$

is numerically always small relative to D and is omitted from the main text (p. 309).

**Appendix C (short correlation time).** |𝒱|² ~ ⟨|V_ik|²⟩ = ∫∫dθdψ sin θ sin ψ (P̃_i^*T_k)_θ(P̃_i^*T_k)^*_ψ⟨(δα(θ,τ)/α_0)(δα(ψ,τ)/α_0)⟩ (C.1), with correlation function (f²/2N_c sin θ)θ_cδ(θ−ψ), gives |𝒱|² ~ (πf²/2N_c²)∫dθ sin θ|P̃_i^*T_k|² (C.2), largest for i = k = 0 where it equals D_0/τ_c. Condition:

$$|\mathcal V|\tau_c\sim(D_0\tau_c)^{1/2}\ll1 \qquad \text{(Eq. C.3)}$$

For the standard run πf²τ_c²/2N_c² ≃ 1.6×10^−3 and the integral 3.529 give (D_0τ_c)^{1/2} ≃ 7.5×10^−2 (p. 310).

**Appendix D (proof of M_μννσ = 0).** M_μννσ = ∫ds⟨V(τ)·V(τ−s)⟩_μσ; in the real representation V^τ·V^{τ−s} has the block form (D.1) with ℜ{V^τ·V^{τ−s}} and ∓ℑ{…}. Using bra–ket notation (Merzbacher 1970) b = |b), b̃ = (b| (round brackets in the paper, to distinguish from the averaging brackets ⟨·⟩), (12) reads (i|j) = δ_ij and completeness |k)(k| = 1, explicitly {|k)(k|}_nm = δ_nm δ(θ−ψ) (D.2) (e.g. Σ_k T̃_k^*(ψ)P_k(θ) = 0 for the n,m = 1,2 element). Since V from (27) satisfies V(1)V(2) = 0 (nilpotent 2×2 structure), 0 = (i|V(1)V(2)|j) = (i|V(1)|k)(k|V(2)|j) = V_ik(τ_1)V_kj(τ_2) (D.3), so V^τ·V^{τ−s} = 0 in any equivalent representation (p. 311).

**Appendix E (scaling of c² and γ).** c², γ are combinations of Λ, Δ_0, Δ_1, which derive from (A.15); that equation depends on kR, C_m and on the overall factor f²τ_c/N_c² in M (e.g. B.8), proving the dependence on f, τ_c, N_c only through f²τ_c/N_c², but the functional form is unknown because ⟨T_+²⟩|_{a=0} is not known as a function of f²τ_c/N_c². Taking it constant:

$$\gamma\sim\left(\frac{f^2\tau_c}{N_c^2}\right)^{-1},\qquad c^2\sim\text{constant} \qquad \text{(Eq. E.1)}$$

The dependence of c², γ on ⟨T_+²⟩ only suggests a direct route via the Ansatz b = a b_0 + b_+ (E.2), b_+ the overtone field, inserted in ḃ = Ob: taking the inner product with b̃_0 gives ȧ = (b̃_0, Ob_0) + (b̃_0, Ob_+), leading to two coupled equations for a and b_+; no credible derivation of (32) was extracted from them, but their existence proves the dipole statistics depend only on the collective effect of all overtones (p. 312).

**Appendix F (evaluation of ⟨T_r⟩).** Write p(a) = exp[ψ(a)], [D(y)p(y)]^−1 = Δ_1^−1 exp[φ(y)] with φ(y) = γy² − γ(c²+1)ln(1 + y²/c²), φ(0) = φ′(0) = 0, φ″(0) = −2γ/c² ≡ φ_0″. Then

$$\langle T_e\rangle=\frac2{\Delta_1}\int_0^{a'}dy\,\exp[\phi(y)]\int_y^\infty dz\,\exp[\psi(z)] \qquad \text{(Eq. F.1)}$$

$$\;\;\simeq\frac2{\Delta_1}\int_0^\infty dy\,\exp\left(\tfrac12\phi_0''y^2\right)\int_{-\infty}^\infty dz\,\exp\left[\psi_m+\tfrac12\psi_m''(z-a_m)^2\right] \qquad \text{(Eq. F.2)}$$

$$\;\;=\frac{2\pi}{\Delta_1}\left(\frac1{\phi_0''\psi_m''}\right)^{1/2}\exp(\psi_m) \qquad \text{(Eq. F.3)}$$

The Gaussian approximation of exp[φ(y)] is adequate when a′ ≈ 1 (Fig. 11); a′ → ∞ is allowed since exp(−γy²/c²) → 0; the lower limit y of the inner integral is replaced by 0 because the outer integrand is concentrated in [0, c/√γ], and after expansion of ψ the lower limit becomes −∞. Substituting φ_0″, (54) and (44) gives (58); substituting (53) and expanding asymptotically gives (59); the only non-trivial step is

$$\left(\frac{\gamma c^2+\gamma-1}{\gamma c^2}\right)^{\gamma c^2+\gamma}=\left(\frac{1+c^2}{c^2}\right)^{\gamma(1+c^2)}\left(1-\frac1{\gamma(1+c^2)}\right)^{\gamma(1+c^2)}\sim\left(\frac{1+c^2}{c^2}\right)^{\gamma(1+c^2)}e^{-1}$$

for large γ (p. 314, unnumbered). Why ⟨T_r⟩ = 2⟨T_e⟩ is only approximate: a = a_0 = 0 would be the separatrix only if a_0 = 0 remained zero without fluctuations, but for a_0 = 0 and V_ik = 0 Eq. (29) gives da_0/dτ = sΣ_{k≠0}E_0k a_k ≠ 0, a residual pulling on the dipole by the overtones that makes escape to left and right unequally probable (p. 313).

## 7. Figures and tables

| Item | Page | What it shows | Parameters / run | Key readable content |
|---|---|---|---|---|
| Fig. 1 | 274 (Colour Plate I, PDF 53) | Top: colour/grey plot of toroidal field ℜT(θ, t), colatitude 0…π vs t = 0…2000 diffusion times. Bottom: normalised dipole amplitude a = ℜa_0 (−2…2). Vertical bars mark reversals. Data spacing 100 time steps. | standard case f = 6.4; footnote-1 parameters | reversal bars at a mean spacing of order 100τ_d; a hovers around ±1 with excursions to |a| ≈ 2; aborted reversals at τ = 825, 863, 1048, also 1475; excursion at τ = 648; vertical stripes = overtone fine structure |
| Fig. 2 | 275 | Four panels for τ = 500–700 (spacing 5 steps): |a_1| and |a_2| (0–1.5), parity 𝒫 (−1…1), "energy" ∫(ℜT)² sin θ dθ (0–30) | standard run | |a_1| spikes to ~1.2, |a_2| to ~0.9; 𝒫 mostly near −1 with spikes to +1; reversal bars at τ ≈ 535, 557 and 660, 678, 685 (two clusters) [uncertain] |
| Fig. 3 | 278 (Colour Plate II, PDF 54) | Two reversals at full resolution; T(θ,t) with black zero line; top τ = 219–225, bottom τ = 1133–1139 | standard run | top: very fast event near τ ≈ 222; bottom: quasi-periodic (quadrupole, period 1.1), recovery after τ = 1135 ≈ exponential with time scale ~Λ^−1 |
| Fig. 4 | 279 (Colour Plate III, PDF 55) | Aborted reversal (τ = 823–829) and excursion at τ ≃ 1282 (τ = 1280–1286) | standard run | dipole sign change vs none; similar appearance |
| Fig. 5 | 280 | a(t), 0–2000 τ_d, for f = 5.4 (top) and f = 7.4 (bottom); same random series as Fig. 1 | footnote-1 parameters | f = 5.4: 6 reversal bars (sign changes near τ ≈ 540, 690, 1180, 1660, plus a short reversed interval ≈ 1760–1790 marked by two bars) [uncertain]; f = 7.4: several tens |
| Fig. 6 | 283 | ℜa_1 (−1…1) vs a = ℜa_0 (−2…2), trajectory | f = 5.4, near τ = 680 | two clouds centred at a ≈ ±1 with |ℜa_1| ≲ 0.5 and one jump |
| Fig. 7 | 285 | Sketch: U(a) = ¼Λ(1−a²)², dashed p(a) from (43), particle in the right well, arrow over the hill; inset diagram δα → a (D_0a²) and δα → overtones → a (D_1(a)) | γ = 3, c² = 0.2 | p(a) peaks near a ≈ ±0.8 |
| Fig. 8 | 292 | 5×10^4 dots of D_1/F (F = πf²τ_c/2N_c²) vs |a| (0–1.5), histogram of bin averages, fitted curve 1.36 + 2.00a² for a ≤ 1.1 | standard run | D_1/F 0–15; histogram rises from ≈1.4 at a = 0 |
| Fig. 9 | 293 | Simulated p(a) (histogram) vs theoretical (43) (curve), both unit area, |a| 0–1.5 | standard run, γ = 2.96, c² = 0.246 | histogram peak ≈1.35 at |a| ≈ 0.85; curve peak ≈1.1 at |a| ≈ 0.8 [uncertain] |
| Fig. 10 | 295 | "Amplitude Distribution": Sint-800 VADM histogram (35 bins, each with an error bar; x: VADM 0–10 ×10^22 A m², y: number of data points 0–80) with fitted p(a) | fit γ = 5.15, c² = 0.25 | peak ≈70 points near 5.5–6×10^22 A m² [uncertain] |
| Fig. 11 | 313 | exp[φ(y)] = Δ_1[D(y)p(y)]^−1 and ∫_y^∞ p(z)dz vs y (0–2) | γ = 3, c² = 0.2 | exp[φ] falls steeply from 1 to a minimum near y ≈ 1 and rises again; ∫p decays by y ≈ 1.5 |
| Table I | 270 | Eigenvalues λ^i, parity of T_i, at C_0 = 57.9 and C_m = 100, kR = 0.5, α/α_0 = cos θ | — | see §4.6; footnotes: (a) s/a = symmetric/antisymmetric about θ = π/2, P and T have opposite symmetry, P and P̃ (and T and T̃) the same; (b) ℜλ^i a decreasing series, printed "asymptotically ~ i(i+1)" (no minus sign in print); (c) all overtones periodic, printed "ℑλ^i ≠ 0 for i ≥ 0" (evidently i ≥ 1) |
| Table II | 292 | Comparison numerical vs analytic (c², γ, variance, ⟨T_r⟩) for f = 4.4–8.4 | N_c = 10, τ_c = 0.05, 5×10^4 τ_d | see §6.8 |

## 8. Comparison with observations

Observed features summarised by the authors (Jacobs 1994; Merrill et al. 1996), p. 293: ⟨T_r⟩ has decreased over the last 160 Myr and is now (2–3)×10^5 yr; a reversal takes ≲10 kyr; no correlation between lengths of successive polarity intervals nor between polarity and length; interval lengths follow a non-stationary Poisson process.

| Observable | Data | Model | Stated agreement |
|---|---|---|---|
| Reversal duration | ≲10 kyr | ~3τ_d (fast events ~0.5τ_d) | used to calibrate τ_d ~ 3 kyr (p. 277, 293) |
| Mean time between reversals | (2–3)×10^5 yr now | ⟨T_r⟩ ≈ 100τ_d for f = 6.4 → ~3×10^5 yr | matches by tuning f (p. 277, 293); the rate itself is not predicted since f²τ_c/N_c² is free (p. 294) |
| Polarity-interval statistics | Poisson, independent of polarity | Poisson for T_r ≫ Λ^−1; independent of sign by construction | reproduced (pp. 290, 293) |
| Rapid variability during reversals | reported (Coe et al. 1995) | quasi-periodic overtone behaviour, period ~τ_d (~3 kyr) | "may be related" (p. 294, 297) |
| Field structure during reversal | — | nearly symmetric/quadrupolar | qualitative (p. 297) |
| Secular variation vs reversal rate | trend suggested in the record (McFadden and Merrill 1995) and in simulations (Glatzmaier et al. 1999) | larger forcing → larger variance and more reversals (Table II, Fig. 5) | simple theoretical explanation (p. 294) |
| Relative variance of dipole between reversals | 0.1 (last 130 kyr, Merrill et al. 1996 Fig. 4.9), 0.15 (Matuyama 2–2.5 Myr), 0.2 (Gauss 2.5–3 Myr) (Valet and Meynadier 1993) | ~0.16 for once-per-100τ_d reversals (Table II) | "reasonable agreement" (p. 294); abstract: present rate implies ~15 % secular variation |
| Dipole amplitude distribution | Sint-800 histogram (800 points, 35 bins) | p(a) from (43): γ = 5.15, c² = 0.25 fits at χ² ↔ 1σ; run-#2 values (5.15, 0.3) at 2σ | "surprisingly good" with no retuning of kR, C_m; the geodynamo over the past 800 kyr behaved like run #2 (p. 295) |
| Time since last reversal | 780 kyr | run #2: ⟨T_r⟩ ~ 500τ_d ~ 1.5 Myr | "would be consistent … but this requires more careful statistics" (p. 295) |
| Cretaceous superchron / long-term rate increase | long polarity bias ~100 Myr ago; slow increase of reversal rate and SV since | factor-2 decrease of f²τ_c/N_c² → ⟨T_r⟩ ×100 (~20 Myr); slow factor-2 increase of f²τ_c/N_c² since | offered as explanation; elaborated in Schmitt et al. (2000) (pp. 299–300) |
| Mean-field / αΩ validity | σ = 6×10^5 S m^−1; inner core super-rotation ≲1° yr^−1 | β/η ~ 80; |C_α| ~ 1 ≪ |C_Ω| ~ 100 | "self-consistent", not proven (p. 296) |

## 9. Assumptions, caveats and limitations

Stated by the authors:
- Mean-field theory is applied indirectly to the geodynamo; its validity is only argued to be self-consistent (β/η ~ 80, |C_α| ≪ |C_Ω|), "estimates of this kind may always be stretched one way or the other" (pp. 295–296).
- The model is 1-D in θ with a local radial ansatz; r is a free parameter set to R; no inner core (pp. 267–268, 298). The results are claimed not to depend on the model being one-dimensional, only on eigenvalues and overlap integrals (p. 298).
- The αΩ approximation and α ∝ cos θ are "customary" choices (p. 268); η is neglected relative to β (p. 265).
- Nonlinearity: global α-quenching with q from the fundamental mode only, chosen for tractability; the authors tested alternatives and found no material change (pp. 270, 277, 281). The constants in (18)–(19) only set the field unit (p. 271).
- Fluctuation model: only δα fluctuates (not β or v, "never been studied yet", p. 299); no magnetic quenching of δα (p. 272, 281); homogeneously distributed renewing cells; details of cell size/distribution irrelevant except the cos θ mean and (sin θ)^{−1/2} fluctuation scaling (p. 272). The required fluctuation level (δα_rms/α_0 cos θ ≳ 2.5) is large and it is unknown whether it can be attained; no prediction of f exists from theory or simulations (p. 299).
- Reversal/excursion definitions are "dictated by theoretical convenience rather than by geophysical practice" and subjective (p. 276, footnote 2).
- Fokker–Planck derivation assumptions: short correlation time |𝒱|τ_c ~ (D_0τ_c)^{1/2} ≃ 0.075 ≪ 1 (Appendix C); exp(±s𝒜) ≈ 1, i.e. high-order modes (with time scales comparable to τ_c) contribute little (p. 305); δ-peaked θ-correlation, justified if low-order modes dominate T_+ (p. 309); g_2(a) neglected as numerically small (p. 309).
- Drift approximation S ≃ Λ(1−a²)a holds only for |a| ≲ 0.4; it degrades p(a) (Fig. 9) but barely affects ⟨T_r⟩ (pp. 292–293).
- Parabolic approximation and (58)–(59) require sufficiently large γ (γ > 1 strictly); run 5 (γ = 1.09) is outside the validity range (p. 288, Table II).
- ⟨T_r⟩ = 2⟨T_e⟩ is only approximate because a = 0 is not a separatrix (residual pulling by overtones) (p. 313).
- c² and γ cannot be expressed in the model parameters; ⟨T_+²⟩|_{a=0} as a function of f²τ_c/N_c² is unknown, so a numerical fit of D_1 is needed (p. 291, Appendix E).
- The Sint-800 comparison is preliminary; a parameter study is deferred (p. 295).
- It remains unclear why a comparable α² dynamo study (Hoyng and Van Geffen 1993) showed no reversals (p. 298).

**Note:** the time unit τ_d ≈ 3 kyr is calibrated from the reversal duration, not derived; the later paper (PEPI 2002) obtains τ_d ≈ 10–15 kyr from the Sint-800 autocorrelation, so all physical-time conversions here (150-yr τ_c, 3×10^5-yr ⟨T_r⟩ for 100τ_d, 20-Myr superchron) scale with that choice.
**Note:** with Δτ = 10^−2 and τ_c = 0.05, each random draw of F persists for only 5 time steps; the noise is therefore nearly white relative to the mode time scales (overtone decay 0.1–0.5τ_d) but not relative to Δτ, which matters for reproducing D_0 numerically.
**Note:** the fluctuation amplitude ∝ (sin θ)^{−1/2} diverges at the poles; the paper relies on P = T = 0 there. A discrete implementation must avoid evaluating 1/√sin θ exactly at θ = 0, π.
**Note:** the paper never states whether F is Gaussian, uniform or two-valued; only zero mean and unit variance are specified.

## 10. Relation to the other two papers in this set

- This paper (GAFD 2001, received June 2000, final form December 2000) is the long paper with the full model definition, simulations, the Fokker–Planck derivation (Appendices A–F) and the first observational comparison. It explicitly refers to the short companion Schmitt, Ossendrijver and Hoyng (submitted 2000; = PEPI 125, 119–124, 2001) as having "recently applied these ideas in the context of the geodynamo" and delivered the amplitude distribution, the reversal picture and the SV–reversal-rate relation (p. 266), and as elaborating the slow increase of reversal rate since the Cretaceous via a factor-2 change of f²τ_c/N_c² (p. 300).
- Paper 2 (Schmitt, Ossendrijver, Hoyng, PEPI 2001) is the compact 6-page summary of the same method: Fokker–Planck equation for the dipole moment, bistable potential, normal/reversed minima, SV–reversal-rate relation and the Sint-800 amplitude distribution. Everything it states is derived here; it adds no new model ingredients.
- Paper 3 (Hoyng, Schmitt, Ossendrijver, PEPI 130, 143–157, 2002; cited here as "submitted (2001)") is the detailed Sint-800 analysis announced on p. 295 ("A parameter study is now in progress, and the results will be reported elsewhere (Hoyng et al., 2001)") and on p. 293 (the improved g_1(a) correction). It adds the autocorrelation function (mentioned here as computable but not computed, p. 298), determines τ_d = R²/β ≈ 10–15 kyr from Sint-800, reproduces the VADM distribution, SV and mean reversal interval, and finds the reversal duration a factor ~2 too long.
- Division of labour: full derivations — this paper; data analysis (autocorrelation, τ_d, duration discrepancy) — PEPI 2002; compact summary — PEPI 2001.

## 11. Implementation guide for coding agents

**System to solve.** Two real or complex scalar fields P(θ,τ), T(θ,τ) on θ ∈ [0, π]:

∂P/∂τ = L P + [cos θ (1 + s(1−a²)) + δα(θ,τ)/α_0] T
∂T/∂τ = C_0 ∂(sin θ P)/∂θ + L T

with L from Eq. (8) (including −(kR)², kR = 0.5), C_0 = 57.9, s = 0.728, boundary conditions P = T = 0 at θ = 0, π, and a = ℜ(b̃_0, b) = ℜ∫_0^π sin θ (P̃_0^* P + T̃_0^* T) dθ evaluated at every step from the current field. The fields can be kept real (**Note:** inference, not a statement of the paper — with real G, E, V and real initial data (26) preserves reality; the paper only says the fields "are in general complex" and that quenching is done by the real field ℜT, p. 270; the overtone amplitudes a_k = (b̃_k, b) are complex regardless, because b̃_k is). The paper's own simulations use the equivalent operator form (26)–(27).

**Eigenvalue problem.** Discretise G (Eq. 6) at C = C_0 with α/α_0 = cos θ on the same grid; solve for eigenpairs of G and of its adjoint G̃ (or use left eigenvectors of the discrete matrix, with the sin θ weight of inner product (13)). Order modes by decreasing ℜλ^i; check against Table I: at C_0 the leading eigenvalue is 0 (real, dipole: T antisymmetric about θ = π/2, P symmetric), then −2.76 ± 4.41i (quadrupole), −8.67 ± 4.97i, −15.9 ± 6.13i; at C_m = 100: 0.527, −2.16 ± 5.88i, −8.21 ± 6.04i, −15.7 ± 7.70i (Table I lists only ℑλ^i > 0; the conjugate partners follow from G being real). First locate C_0 by scanning C until ℜλ^0 = 0 (target 57.9 for kR = 0.5). Normalise with (b̃_i, b_j) = δ_ij; take b_0 real. Validation integrals: E_00 = ∫ sin θ cos θ P̃_0 T_0 dθ = 1.048 and ∫ sin θ (P̃_0T_0)² dθ = 3.529 (both are invariant under the rescaling b_0 → cb_0, b̃_0 → b̃_0/c permitted by (12), so they test the adjoint solution and the bi-orthonormalisation rather than the field unit).

**Field unit.** With b_0 real and bi-orthonormal, b = ±b_0 is an exact stationary solution of (26) without noise (a = ±1, G b_0 = 0). The scale of b_0 is a free choice that only sets the unit of the physical field; a is always measured relative to it.

**Stochastic α.** δα/α_0 = f F(θ,τ)/√(2N_c sin θ) with N_c = 10, f = 6.4 (standard). F: zero mean, unit variance; constant within each of N_c latitude bins of width π/N_c and within each interval of length τ_c = 0.05 (= 5 time steps at Δτ = 10^−2); redrawn independently per bin at the end of each τ_c. The distribution of F is not specified (gap). Whether the latitude-bin boundaries are aligned with grid points is not specified (gap). The paper used the same random sequence across f values when comparing runs (Fig. 5).

**Numerics.** 60 equidistant θ points on [0, π], Δτ = 10^−2 (footnote 1, p. 273). The time integrator and the finite-difference form of L and of the shear term are not stated (gap). The largest diffusive eigenvalue on a 60-point grid is ~4/Δθ² ≈ 1.4×10³ per τ_d, above the explicit-Euler stability limit 2/Δτ = 200, so the diffusion operator presumably needs implicit or sub-stepped treatment — the paper does not say how this was handled (**Note:** mine, not the paper's). Run length for statistics: 5×10^4 τ_d (≈ 470 reversals at f = 6.4).

**Diagnostics as defined in the paper.**
- Dipole amplitude a(τ) = ℜ(b̃_0, b) (Eq. 15/24); overtone amplitudes a_k = (b̃_k, b); |a_1|, |a_2| plotted in Fig. 2.
- Parity 𝒫 from Eq. (28) with ℰ_{S,A} = ∫(ℜT_{S,A})² sin θ dθ; "energy" ∫(ℜT)² sin θ dθ.
- Reversal counting (footnote 2, p. 276): group zero crossings of a that occur within "a few τ_d"; an odd number in a group = one reversal, an even number = aborted reversal; approach to zero without crossing = excursion. The grouping window is not specified (gap). ⟨T_r⟩ = run length / number of reversals, with the quoted uncertainties (e.g. 107 ± 5) consistent with ⟨T_r⟩/√N for Poisson counts.
- Secular-variation measure: ⟨(a−⟨a⟩)²⟩/⟨a⟩² computed over intervals between reversals (Eq. 56 is the theory; "from run" in Table II is the measurement). How the reversal transitions themselves were excluded is not specified (gap).
- Amplitude distribution: histogram of consecutive values of a (Fig. 9); compare with (43) after normalising both to unit area.
- D_1 measurement: every 1τ_d evaluate F·∫ sin θ P̃_0² T_+² dθ with T_+ = ℜ(T − a_0T_0), bin in |a|, fit a quadratic for a ≤ 1.1 (Fig. 8). Expected for the standard run: D_1/F ≈ 1.36 + 2.00a².
- Derived parameters: Δ_1 = D_1(0), Δ_0 = D_0 + ½D_1″(0) with D_0 = F·3.529; c² = Δ_1/Δ_0, γ = Λ/Δ_0 with Λ = 0.527; then (52)–(56), (58)–(60).

**Expected outputs (standard run).** c² ≃ 0.246, γ ≃ 2.96; predicted variance 0.17 vs measured 0.16; predicted ⟨T_r⟩ ≃ 94 vs measured 107 ± 5 (468 reversals in 5×10^4 τ_d); typical reversal duration ~3τ_d; 𝒫 ≃ −1 between reversals; quadrupole period ≈ 1.1τ_d during some reversals; (D_0τ_c)^{1/2} ≃ 0.075. Other f values: Table II (§6.8).

**Pitfalls.**
- Time unit: all rates are per τ_d; physical conversion uses τ_d ~ 3 kyr in this paper (10–15 kyr in PEPI 2002).
- The λ^i used in (29) are those at C_0, not C_m; but the drift constant in (33) is Λ = λ^0 at C_m (0.527), not μ = sE_00 (0.763) — the difference is the overtone contribution g_1.
- The expansion a_k is in modes of the critical (C_0) operator, so the equilibrium is exactly a_0 = ±1 with all a_k = 0.
- p(a) in (43) is unnormalised; (50)–(56) use the parabolic approximation valid only for γ well above 1.
- (58) counts reversals, not escapes: ⟨T_r⟩ ≈ 2⟨T_e⟩ with the caveat of Appendix F. Aborted reversals (even crossings) are not reversals.
- The (sin θ)^{−1/2} factor must never be evaluated at θ = 0, π.
- The reported statistical results were insensitive to the type of nonlinearity and to C_m, kR; this can be used as a regression test (p. 277).

## 12. Symbol glossary

| Symbol | Meaning | Dimension / unit | First defined |
|---|---|---|---|
| $\langle\mathbf B\rangle$ | mean (azimuthally averaged) magnetic field | field | Eq. 1, p. 265 |
| $\mathbf v,\ \mathbf u$ | mean flow, turbulent convection | velocity | p. 265 |
| $\alpha,\ \beta,\ \eta$ | α-effect, turbulent resistivity, molecular resistivity | cm s^−1; cm² s^−1 | Eqs. 1–2, p. 265 |
| $\tau_c,\ \lambda_c$ | correlation time (eddy turnover), correlation length (eddy size) | τ_c in units of τ_d (0.05); λ_c ~ 500 km | p. 265, 272, 296 |
| $A,\ B$ | poloidal vector potential (A e_φ), toroidal field | field | Eq. 3, p. 267 |
| $P(\theta,t),\ T(\theta,t)$ | poloidal, toroidal scalar amplitudes (same dimension) | field units | Eq. 4, p. 267 |
| $k,\ R$ | radial wavenumber, reference radius (kR = 0.5) | 1/length, length | Eq. 4, p. 267 |
| $\alpha_0,\ \alpha_m$ | reference α at North pole for critical / operating state | cm s^−1 | p. 268, 269 |
| $\Omega,\ \Omega'$ | angular velocity, dΩ/dr (constant) | s^−1, s^−1 cm^−1 | p. 267–268 |
| $\mathbf b=(P,T)^T$ | field vector | — | Eq. 5 |
| $G,\ L$ | dynamo operator, diffusion operator | per τ_d | Eqs. 6, 8 |
| $C,\ C_0,\ C_m$ | dynamo number; critical (57.9); operating (100) | dimensionless | Eqs. 7, 16 |
| $\tau,\ \tau_d$ | dimensionless time; turbulent diffusion time R²/β | —; ~3 kyr | Eq. 7 |
| $\mathbf b_i,\ \lambda^i$ | eigenmodes, eigenvalues of G | —; per τ_d | Eq. 10 |
| $\tilde G,\ \tilde{\mathbf b}_i$ | adjoint operator and eigenmodes (P̃_i, T̃_i) | — | Eq. 11 |
| $(\cdot,\cdot)$ | inner product with weight sin θ | — | Eq. 13 |
| $a_k$ | mode amplitudes (b̃_k, b), in units of equilibrium a_0 | dimensionless | Eq. 15; p. 272, 282 |
| $q,\ q_0$ | quenching measure; its equilibrium value (0.421) | dimensionless | Eqs. 18–19, 22 |
| $h(q)$ | quenching function (= 1 − q) | — | p. 270–271 |
| $\delta\alpha,\ F,\ f,\ N_c,\ N,\ \theta_c$ | α fluctuation; unit-variance random function; strength; cells per π in latitude (10); cells on a latitude circle 2N_c sin θ; cell width π/N_c | δα in units of α_0 (cm s^−1); the rest dimensionless | Eq. 21, p. 271–272 |
| $s$ | supercriticality C_m/C_0 − 1 (0.728) | — | Eq. 22 |
| $a$ | ℜa_0, dipole amplitude | — | Eq. 24 |
| $E,\ V$ | nonlinearity and noise operators (2×2) | — | Eq. 27 |
| $\mathcal P,\ \mathcal E_{S,A},\ T_{S,A}$ | parity; symmetric/antisymmetric toroidal energies and fields | — | Eq. 28, p. 276 |
| $E_{ik},\ V_{ik}$ | overlap matrices (E_00 = 1.048) | — | Eqs. 30–31 |
| $p(a,\tau)$ | probability density of a | — | Eq. 32 |
| $S,\ D$ | drift, diffusion coefficients | per τ_d | Eqs. 32–35 |
| $\mu,\ \Lambda$ | sE_00 = 0.763; λ^0 at C_m = 0.527 | per τ_d | Eq. 34 |
| $g_1(a),\ g_2(a)$ | overtone corrections to S and D | per τ_d | Eq. 33, B.3, B.13 |
| $D_0,\ D_1(a)$ | direct and overtone diffusion coefficients | per τ_d | Eqs. 36–37 |
| $T_+$ | toroidal field of overtones ℜ(T − a_0T_0) | field units | Eq. 38 |
| $\langle\cdot\rangle\vert_a$ | conditional average at fixed a | — | Eq. 37, A.19 |
| $\Delta_0,\ \Delta_1$ | D_0 + ½D_1″(0); D_1(0) | per τ_d | Eq. 40 |
| $U(a)$ | bistable potential ¼Λ(1−a²)² | per τ_d | Eq. 41 |
| $\gamma,\ c^2$ | Λ/Δ_0 (barrier measure); Δ_1/Δ_0 (diffusion at a = 0) | — | Eq. 44 |
| $F$ (p. 287) | probability current Sp − ½∂(Dp)/∂a in ∂p/∂τ = −∂F/∂a (distinct from the random function F(θ,τ) and from the prefactor F of Fig. 8) | per τ_d | Eq. 49 |
| $\psi(a),\ a_m,\ \psi_m,\ \psi_m''$ | log-density, peak position, value and curvature | — | Eqs. 50–54 |
| $\langle T_e\rangle,\ \langle T_r\rangle$ | mean escape time, mean time between reversals | τ_d | Eqs. 57–58 |
| $K$ | −1 + (1 + c²) ln((1 + c²)/c²) | — | Eq. 60 |
| $F$ (p. 292) | πf²τ_c/2N_c² (= 3.22×10^−2 standard) | per τ_d | Fig. 8 caption |
| $\mathbf a,\ \mathbf B,\ \mathbf L,\ \mathbf E,\ \mathbf V,\ \mathbf N,\ \mathbf M$ | real supervector / supermatrices | — | Eqs. A.2–A.5, A.12–A.13 |
| $\mathcal A,\ \mathcal V$ | Liouville operators (systematic, noise) | per τ_d | Eqs. A.7–A.8 |
| $\Pi(\mathbf a,\tau)$ | joint probability of all mode amplitudes | — | p. 304 |
| $\phi(y),\ \phi_0''$ | γy² − γ(c²+1)ln(1 + y²/c²); −2γ/c² | — | p. 312 |
| $C_\alpha,\ C_\Omega$ | α_mR/β, R³Ω′/β (~ −1, ~ −100) | — | p. 296 |

## 13. Section-by-section outline

| Section | Pages | Content |
|---|---|---|
| Abstract, keywords | 263 | Summary; bistable potential ∝(1−a²)²; present reversal rate implies ~15 % SV; Sint-800 agreement. |
| 1. Introduction | 264–266 | Observed properties; 3-D simulations; helicity-fluctuation reversal ideas (Parker 1969 …); mean-field Eq. (1)–(2); azimuthal average as origin of fluctuations; relation to Schmitt et al. (2000); outline. |
| 2. Dynamo model — 2.1 The linear system | 267–269 | Eqs. (3)–(9); 2.1.1 Eigenfunctions and parameters: Eqs. (10)–(16), kR = 0.5, C_0 = 57.9, C_m = 100. |
| 2.2 Nonlinear feedback and fluctuations | 269–272 | 2.2.1 Nonlinearity (Eqs. 17–19, Table I on p. 270); 2.2.2 Nonlinear equilibrium (C_m h(q_0) = C_0); 2.2.3 Fluctuations (Eqs. 20–21, cell model, f²τ_c/N_c²). |
| 2.3 Model equations | 272–273 | Eqs. (22)–(27); two equilibria a = ±1. |
| 3. Simulations — 3.1 Standard case | 273–276 | Footnote 1 (numerics); Figs. 1–2; three geomagnetic-like features; parity (28); footnote 2 (reversal definitions). |
| 3.2 Reversals and aborted reversals | 276–277 | Four-point phenomenology; Figs. 3–4 (pp. 278–279); τ_d ~ 3 kyr. |
| 3.3 Other cases and scaling | 277–281 | Fig. 5 (p. 280); scaling with f²τ_c/N_c²; alternative nonlinearities; quenching of δα. |
| 4. Statistical analysis (intro) | 281 | Recap; Kramers analogy announced; method independent of model details. |
| 4.1 Coupled mode equations | 282 | Eqs. (29)–(31); E_00 = 1.048. |
| 4.2 Fokker–Planck equation | 282–284 | Fig. 6 (p. 283); Eqs. (32)–(40); integral 3.529. |
| 4.3 Particle-in-well analogy | 284–286 | Eqs. (41)–(45); Fig. 7 (p. 285). |
| 4.4 The reversal mechanism | 286–287 | Eqs. (46)–(49); three-step argument; probability current. |
| 4.5 The statistics of field reversals | 288–290 | Eqs. (50)–(54); 4.5.1 variability (55)–(56); 4.5.2 mean reversal rate (57)–(60), Kramers; 4.5.3 Poisson intervals. |
| 5. Model validation — 5.1 Internal consistency | 291–293 | D_1 fit (61)–(63); Fig. 8, Table II (p. 292); Fig. 9 (p. 293) and its explanation. |
| 5.2 Comparison with the geodynamo | 293–296 | Observed reversal properties; 5.2.1 SV vs reversal rate (VADM variances 0.1/0.15/0.2); 5.2.2 Sint-800 (Fig. 10, p. 295); 5.2.3 mean-field approximation (β, η, C_α, C_Ω); 5.2.4 characteristic numbers (why C_m = 100; τ_c, N_c, fluctuation ratio). |
| 6. Discussion and summary | 297–300 | Key elements; three main results; relation to Parker, Mazaud–Laj, Narteau; generality; unresolved issues; additive vs multiplicative noise; superchron explanation; role of nonlinearities. Acknowledgements p. 300. |
| References | 301–303 | 60 entries. |
| Appendix A. Derivation of the Fokker–Planck equation | 303–307 | A.1 Notation (A.1–A.5); A.2 Probability distribution of {a_ν} (A.6–A.15); A.3 Probability distribution of a_0 (A.16–A.19). |
| Appendix B. Drift and diffusion coefficients | 307–309 | B.1 Drift (B.1–B.3); B.2 Diffusion (B.4–B.13). |
| Appendix C. Short correlation time | 310 | (C.1)–(C.3); (D_0τ_c)^{1/2} ≃ 0.075. |
| Appendix D. Proof of M_μννσ = 0 | 310–311 | (D.1)–(D.3), completeness relation. |
| Appendix E. The scaling of c² and γ | 311–312 | (E.1)–(E.2). |
| Appendix F. Evaluation of ⟨T_r⟩ | 312–314 | (F.1)–(F.3), Fig. 11 (p. 313), unnumbered asymptotic identity (p. 314), why ⟨T_r⟩ ≠ 2⟨T_e⟩ exactly. |
| Colour Plates I–III | PDF 53–55 | Colour versions of Figs. 1, 3, 4. |

## 14. Key references cited

- Parker, E. N. (1969), Astrophys. J. 158, 815–827 — the original proposal that helicity fluctuations can reverse the geodynamo; this paper is "a generalisation" of it (stationary fluctuation level instead of a single jump in α) (pp. 264, 297–298).
- Levy, E. H. (1972), Astrophys. J. 171, 635–642; Olson (1983); Olson and Hagee (1990); Le Mouël, Allègre and Narteau (1997); Narteau et al. (2000), PEPI 120, 271–287; Gibbons (1998) — mean-field helicity-fluctuation reversal models that this work extends; Narteau's hierarchical turbulence model contrasted with the simpler cell model (pp. 265, 272, 297).
- Moffatt (1978); Krause and Rädler (1980) — mean-field induction equation and α, β expressions (p. 265).
- Braginskii (1965a,b) — azimuthal average interpretation of the mean field (p. 266).
- Schmitt and Schüssler (1989), Astron. Astrophys. 223, 343–351 — the αΩ model geometry adopted and the complicated C-dependence of eigenvalues (pp. 267, 269).
- Ossendrijver, Hoyng and Schmitt (1996), Astron. Astrophys. 313, 938–948; Hoyng (1993), A&A 272, 321–339; Hoyng (1996) — solar-dynamo α-fluctuation work (phase–amplitude correlation) from which the method and the 1-D equations derive (pp. 266, 267, 271, 299).
- Kleeorin and Ruzmaikin (1984); Ossendrijver (1996, PhD thesis) — normal and adjoint eigenvalue problems (p. 269).
- Brandenburg et al. (1989), A&A 213, 411–422 — global α-quenching and the parity definition (pp. 270, 276).
- Rüdiger (1973) — local α-quenching (p. 270).
- Van Kampen (1976) Phys. Reports 24C, 171–228; Van Kampen (1992) Stochastic Processes in Physics and Chemistry; Gardiner (1990); Risken (1984) — stochastic-process theory: multiplicative-noise expansion (A.9), parabolic approximation, escape times, Kramers (pp. 283, 288–290, 305).
- Hoyng and Van Geffen (1993), GAFD 68, 203–236 — comparable α² study without reversals, flagged as unexplained (p. 298).
- Jacobs (1994); Merrill, McElhinny and McFadden (1996) — observational facts (reversal rate, duration, Poisson statistics, conductivity 6×10^5 S m^−1, VADM Fig. 4.9) (pp. 264, 293–296).
- Valet and Meynadier (1993), Nature 366, 234–238 — VADM records for Matuyama and Gauss chrons (p. 294).
- Guyodo and Valet (1999), Nature 399, 249–252 — Sint-800 record (p. 294–295).
- McFadden and Merrill (1995), PEPI 91, 253–260 — SV–reversal-rate trend and "fundamental transitions" (pp. 294, 297, 299).
- Glatzmaier et al. (1999), Nature 401, 885–890; Glatzmaier and Roberts (1995a,b, 1996) — 3-D simulations, CMB heat-flux control of reversals, inner-core rotation (pp. 264, 266, 294, 300).
- Crossley, Jensen and Jacobs (1986), PEPI 42, 143–153; Meinel and Brandenburg (1990), A&A 238, 369–376; Farrell and Ioannou (1999a,b), ApJ 522 — additive-noise dynamo studies; "interesting in itself", but the paper cautions that a physical basis for the additive term must be ascertained, since it makes Eq. (1) inhomogeneous and fixes the field magnitude even without nonlinearities (p. 299).
- Gubbins (1999), Geophys. J. Int. 137, F1–F3 — reversal duration set by inner-core decay time; motivates a 2-D extension (p. 298).
- Merzbacher (1970) — bra–ket notation and completeness relation used in Appendix D (p. 311).

## 15. Open questions and unresolved points

- Two "main unresolved issues" (p. 299): whether the required large fluctuation level (δα_rms/α_0 cos θ ≳ 2.5) can be attained in the core, and the scaling of ⟨T_+²⟩|_{a=0} (hence c², γ) with f²τ_c/N_c², which would remove the need for the numerical fit of D_1 (pp. 291, 298, 311–312).
- No credible derivation of (32) from the two-field Ansatz b = a b_0 + b_+ (Appendix E, p. 312).
- Why an α² dynamo in a comparable study (Hoyng and Van Geffen 1993) showed no reversals (p. 298).
- Fluctuations in β and v have never been studied (p. 299).
- The discrepancy between simulated and theoretical p(a) (Fig. 9) is attributed to the drift approximation; the corrected g_1 is deferred to Hoyng et al. (2001/2002) (p. 293).
- ⟨T_r⟩ ≈ 2⟨T_e⟩ is only approximate; the predicted ⟨T_r⟩ (94) is 12 % below the measured 107 ± 5 for the standard run, and run 4 (predicted variance 0.47 vs 0.30) shows the parabolic approximation deteriorating as γ → 1 (Table II; **Note:** the run-4 remark is an inference from the table — the paper itself flags only run 5, footnote b of Table II).
- The Sint-800 fit prefers c² = 0.25 to run #2's 0.3; a parameter study in kR and C_m is announced (p. 295). The consistency of run #2's ⟨T_r⟩ ~ 1.5 Myr with the 780-kyr age of the last reversal "requires more careful statistics" (p. 295).
- Autocorrelation of the dipole amplitude and the frequency of excursions of given magnitude are listed as computable but not computed (p. 298).
- A 2-D (r, θ) model with an inner core is identified as the most desirable generalisation, to test Gubbins' (1999) suggestion that reversal duration is gauged by the inner-core resistive decay time (p. 298).
- Whether the Cretaceous superchron and the slow post-Cretaceous rise in reversal rate and SV reflect a factor-2 change of f²τ_c/N_c² — possibly tied to inner-core growth or CMB heat-flux changes — is raised as speculation (pp. 299–300).
