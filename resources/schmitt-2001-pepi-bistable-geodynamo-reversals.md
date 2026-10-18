---
title: "Magnetic field reversals and secular variation in a bistable geodynamo model"
authors: ["D. Schmitt", "M.A.J.H. Ossendrijver", "P. Hoyng"]
year: 2001
journal: "Physics of the Earth and Planetary Interiors"
volume: 125
pages: "119–124"
doi: "10.1016/S0031-9201(01)00237-0"
source_pdf: "papers/Magnetic_Field_Reversals_and_Secular_Var.pdf"
pdf_pages: 6
page_mapping: "PDF page 1 = journal page 119; journal page = PDF page + 118"
series: "Hoyng–Schmitt–Ossendrijver bistable geodynamo series, paper 2 of 3 (chronological)"
keywords: ["Dynamo model", "Helicity fluctuations", "Reversals", "Secular variation"]
extraction_notes: "Vector PDF with a good text layer; equations (1)–(4), the unnumbered T_r formula, the potential U(a) and all figure captions were re-transcribed from 2x page images and 4–5x zoomed crops. Values read off figure axes rather than captions are marked as approximate."
---

# Magnetic field reversals and secular variation in a bistable geodynamo model

## 1. One-paragraph summary

A compact account of a stochastic mean-field geodynamo and a method for computing its statistics. The model is a one-dimensional $\alpha\Omega$ dynamo on a sphere of radius $R$ whose axisymmetric mean field depends on colatitude $\theta$ only (Eq. 1, p. 120). The $\alpha$-effect carries rapid random fluctuations $\delta\alpha$ representing the finite number of convective eddies per latitude circle; growth is limited by $\alpha$-quenching on the energy of the fundamental dipole mode. Runs at dynamo number $C=100$ ($C_{\rm crit}=58$) show a predominantly dipolar field of constant polarity with sudden reversals, aborted reversals and Poisson-distributed reversal intervals (p. 121). Expanding in eigenmodes of the marginal linear operator gives coupled mode equations (Eq. 2); averaging over the fluctuations gives a Fokker–Planck equation (Eq. 3) for the dipole amplitude $a$, equivalent to an overdamped particle in the bistable potential $U(a)=\tfrac14\Lambda(1-a^2)^2$ (p. 122). Its stationary solution (Eq. 4) and an escape-time formula for the mean reversal time $T_r$ explain why secular variation grows while $T_r$ shrinks exponentially with the forcing parameter $f^2\theta_c^2\tau_c$. With a 2000-yr diffusion time, $f=6.4$ gives $T_r\approx100$ diffusion times and relative variance 0.16 (p. 121); the Sint-800 VADM histogram is fitted within $1\sigma$ by Eq. (4) with $\gamma=5.2$, $c^2=0.26$ (p. 123).

## 2. Key takeaways

- 1D $\alpha\Omega$ mean-field dynamo, field a function of $\theta$ only, $\alpha=\alpha_0\cos\theta$, constant radial shear; for $C>0$ the first excited mode is a non-oscillatory dipole, all overtones (dipolar and quadrupolar) are oscillatory and need considerably larger amplification; $C$ is chosen so that the fundamental dipole is supercritical while the overtones are damped (p. 120).
- Standard run: $C=100$, $C_{\rm crit}=58$, $f=6.4$, $\theta_c=\pi/10$, $\tau_c=5\times10^{-2}$; run length 25 × the 2000 diffusion times displayed, i.e. 50,000 (Fig. 1 caption, p. 121).
- Time unit = turbulent diffusion time, "estimated as 2000 years" (p. 121); so $T_r\approx100$ diffusion times $\approx2\times10^5$ yr, 500 diffusion times $\approx10^6$ yr, 50,000 diffusion times $\approx10^8$ yr (pp. 121, 123).
- Relative variance of $a$ between reversals: 0.11 ($f=6$, $T_r\approx200$), 0.16 ($f=6.4$, $T_r\approx100$), 0.24 ($f=7$, $T_r\approx50$); observed secular variation 0.15–0.20 for records dated 2–3 Myr (p. 121).
- All statistics depend only on the forcing parameter $f^2\theta_c^2\tau_c$ (p. 120); $T_r$ depends exponentially on it (pp. 121–122).
- Fokker–Planck Eq. (3) with drift $S\approx\Lambda(1-a^2)a$ is the equation of a heavily damped particle in $U(a)=\tfrac14\Lambda(1-a^2)^2$; symmetry from $\pm B$ invariance, central hill from supercriticality, side walls from the nonlinearity (p. 122, Fig. 2).
- Stationary $p(a)\simeq(1+a^2/c^2)^{\gamma c^2+\gamma-1}\exp(-\gamma a^2)$ (Eq. 4), $\gamma$ and $c^2$ measuring $U(0)$ and $D(0)$ (p. 122).
- $T_r\simeq4\int_0^1dx\,[D(x)p(x)]^{-1}\int_x^\infty dy\,p(y)$ (unnumbered, p. 122).
- Ramping $f$ from 4 to 6.5 over 50,000 diffusion times shortens chrons drastically while the rms variability rises from 17% to 41% (relative variance 0.03 → 0.17); the observed fall of chron lengths from $10^7$ to $10^5$ yr over the last $10^8$ yr is attributed to a factor $\sim2.5$ rise of the forcing parameter (Fig. 3, p. 123).
- Sint-800 VADM relative variance 0.07 ⇒ $f=5.4$ ⇒ $T_r=500$ diffusion times ($\sim1$ Myr), consistent with the 780,000-yr Brunhes chron; Eq. (4) fits the 30-bin VADM histogram within $1\sigma$ with $\gamma=5.2$, $c^2=0.26$, $a_{\max}=5.9\times10^{22}$ A m$^2$ (p. 123).
- $a_0\leftrightarrow$ VADM is argued accurate to $<10\%$: $l=1$ dominates $a_0$, $l=3$ contributes ~20%, the next antisymmetric mode $a_2$ is ~10% of $a_0$, and both fall off faster with height (p. 123).
- Reversals are always fast, caused by interaction of the dipole mode with the overtones; during a reversal $a_0\simeq0$, the system is effectively linear and quadrupolar components appear (pp. 121, 122, 124); the authors conclude that no unique physical cause for a reversal exists in the model, and that this may hold for the geodynamo too (p. 124).

## 3. Context and motivation

The paper starts from the paleomagnetic facts to be explained: reversals at irregular intervals of $10^5$–$10^7$ yr (Jacobs, 1994; Merrill et al., 1996), lasting only a few thousand years (Bogue and Merrill, 1992); considerable fluctuation of field strength between reversals (Valet and Meynadier, 1993; Guyodo and Valet, 1999); and excursions, i.e. aborted reversals (Hoffman, 1981; Gubbins, 1999) (p. 119).

State of the art as described: dynamo action in the liquid outer core (Gubbins and Roberts, 1987; Fearn, 1998); the 3D simulation of Glatzmaier and Roberts (1995), where the outer-core dynamo generates a dipole that "wants to oscillate" on about 100 yr while the conducting inner core, diffusing on a much longer time scale, stabilises it (Hollerbach and Jones, 1993), so a reversal needs the outer-core field to stay reversed long enough to diffuse into the inner core (p. 119). Earlier kinematic schemes (Parker, 1969; Levy, 1972) invoked changes in the latitudinal distribution of cyclonic eddies; Sarson (2000) reviews reversal models from dynamo calculations (p. 120).

The objective is the long-term statistics of the field and a physical explanation of the relation between secular variation and reversal records. This needs runs with many reversals, which 3D codes cannot provide (few reversals, unrealistic parameters), so the authors use mean-field theory with three ingredients from simulations and observations: differential rotation in the outer core, helicity, and large fluctuations, allowing runs of typically $10^8$ yr (p. 120). The random $\alpha$ fluctuations were earlier shown to explain the phase–amplitude correlation of the solar cycle (Hoyng, 1993; Ossendrijver et al., 1996) (p. 120). The paper describes the method, defers all derivations to Hoyng et al. (2001) (pp. 120–121), and applies it to the Sint-800 record (Guyodo and Valet, 1999) and to the long-term reversal-rate trend (McFadden and Merrill, 1995). **Note:** The paper does not state what is new relative to GAFD 2001; from the two abstracts, the Sint-800 distribution comparison appears in both, while the ramped-forcing run (Fig. 3) and the mapping of the Sint-800 variance to $f=5.4$, $T_r=500$ are specific to this paper.

## 4. The model

### 4.1 Physical setup and geometry

A spherical dynamo of fixed radius $R$ with the mean field depending on colatitude $\theta$ only, after the one-dimensional thin-layer model of Schmitt and Schüssler (1989) (p. 120). The field is split into the toroidal part $T$ and the vector potential $P$ of the poloidal field. $\alpha=\alpha_0\cos\theta$ parametrises helical convection (poloidal from toroidal); differential rotation with constant radial shear $d\Omega/dr=\mathrm{const}$, assumed "for reasons of simplicity", makes toroidal from poloidal field; the dynamo is of $\alpha\Omega$ type (Soward, 1991) (p. 120). Mean quantities are longitudinal averages, so the mean field is the axisymmetric ($m=0$) component of the real field (p. 120). **Note:** Fig. 1(a) refers to "a representative radial position, $R$, in the outer core"; the numerical value of $R$ is not given.

### 4.2 Governing equations

Non-dimensional mean-field induction equation (Eq. 1, p. 120):

$$
\frac{\partial}{\partial t}\begin{pmatrix} P \\ T \end{pmatrix}
=\begin{pmatrix} \Delta & (\alpha+\delta\alpha)/\alpha_0 \\[2pt] C\,(\partial/\partial\theta)\sin\theta & \Delta \end{pmatrix}
\begin{pmatrix} P \\ T \end{pmatrix}
\qquad (\text{Eq. }1)
$$

$P$: poloidal vector potential; $T$: toroidal field; $\Delta$: diffusion operator; $\alpha=\alpha_0\cos\theta$; $\delta\alpha$: random fluctuation of $\alpha$; $C$: dynamo number; $t$: time in diffusion times. The dynamo number (p. 120):

$$
C=R_\alpha R_\Omega,\qquad R_\alpha=\frac{\alpha_0R}{\eta},\qquad R_\Omega=\frac{d\Omega}{dr}\frac{R^3}{\eta},
$$

with $\eta$ the turbulent magnetic diffusivity. **Note:** The explicit $\theta$-form of $\Delta$ and the pole boundary conditions are not given here (see Schmitt and Schüssler, 1989; Hoyng et al., 2001). The $\Omega$-entry is printed as the operator $C(\partial/\partial\theta)\sin\theta$ acting on $P$, i.e. $C\,\partial_\theta(\sin\theta P)$ read left to right, consistent with $B_r\propto(\sin\theta)^{-1}\partial_\theta(\sin\theta P)$ sheared by $d\Omega/dr$ and multiplied by $\sin\theta$; the paper does not spell the ordering out. [Verified against the page image: the print is exactly $C(\partial/\partial\theta)\sin\theta$, and for $\Omega=\Omega(r)$ the $\Omega$-effect in $\partial_tB_\phi$ is $r\sin\theta\,B_r\,d\Omega/dr=(d\Omega/dr)\,\partial_\theta(\sin\theta P)$, so $\partial_\theta(\sin\theta P)$ is the only ordering consistent with the stated model; the alternative $\sin\theta\,\partial_\theta P$ would have been printed $C\sin\theta(\partial/\partial\theta)$.]

Mode-coefficient equation (Eq. 2, p. 121), from decomposing the field into eigenfunctions of the unperturbed ($\delta\alpha=0$) linear operator at marginal excitation $C_{\rm crit}=58$:

$$
\frac{\partial a_i}{\partial t}=\lambda_ia_i+(1-a_0^2)\sum_kE_{ik}a_k+\sum_kV_{ik}(t)\,a_k
\qquad (\text{Eq. }2)
$$

$a_i$: amplitude of eigenmode $i$; $a_0$: amplitude of the fundamental dipole mode; $\lambda_i$: eigenvalues at $C_{\rm crit}$; $E_{ik}$, $V_{ik}(t)$: overlap integrals of eigenmodes $i$ and $k$. The first term is the linear evolution, the second the supercritical excitation plus global $\alpha$-quenching by the fundamental mode, which drive the system to its nonlinear equilibrium, the third the fluctuations in $\alpha$; $E_{ik}$ and $V_{ik}$ couple all eigenmodes (p. 121). Derivations are in Hoyng et al. (2001).

Fokker–Planck equation for the distribution $p(a)$ of $a=a_0$, obtained by averaging Eq. (2) over the fluctuations (Eq. 3, p. 122):

$$
\frac{\partial p}{\partial t}=-\frac{\partial}{\partial a}Sp+\frac12\frac{\partial^2}{\partial a^2}Dp
\qquad (\text{Eq. }3)
$$

Drift $S\approx\Lambda(1-a^2)a$ from supercriticality and nonlinear quenching, $\Lambda$ being the growth rate of the fundamental mode at the model's $C$; diffusion coefficient $D(a)$ describes the average effect of the fluctuations and "scales with the forcing parameter"; the overtones' effect is included in both coefficients (p. 122). Eq. (3) is that of a heavily damped particle in a bistable potential $U$ with $-\partial U/\partial a=S$ (p. 122); the Fig. 2 caption gives

$$
U(a)=\tfrac14\,\Lambda\,(1-a^2)^2 .
$$

Stationary solution (Eq. 4, p. 122):

$$
p(a)\propto\frac1D\exp\!\left(\int_0^a\frac{2S}{D}\,dx\right)
\simeq\left(1+\frac{a^2}{c^2}\right)^{\gamma c^2+\gamma-1}\exp(-\gamma a^2)
\qquad (\text{Eq. }4)
$$

where the last expression is "a useful approximation" and $\gamma$, $c^2$ are dimensionless measures of $U(0)$ and $D(0)$ (p. 122). Mean time between reversals, from the mean escape rate from a well (unnumbered, p. 122):

$$
T_r\simeq4\int_0^1dx\,[D(x)p(x)]^{-1}\int_x^\infty dy\,p(y).
$$

**Note:** The right-hand form of Eq. (4) follows exactly from the left-hand form if $D(a)=D(0)(1+a^2/c^2)$ and $\gamma\equiv\Lambda c^2/D(0)$: then $\int_0^a2S/D\,dx=-\gamma a^2+\gamma(1+c^2)\ln(1+a^2/c^2)$, and the $1/D$ prefactor supplies the $-1$ in the exponent. This reconstruction is mine, not the paper's; it agrees with $\gamma$ measuring $U(0)=\Lambda/4$ against $D(0)$ and $c^2$ measuring $D(0)$ ($c$ being the amplitude where the $a$-dependent part of $D$ equals its constant part).

### 4.3 Non-dimensionalisation, units and time scales

- Length unit: $R$ (p. 120; implied by the non-dimensional form, the correlation length $\theta_cR$ and the Reynolds numbers, not stated explicitly).
- Time unit: the turbulent diffusion time, "estimated as 2000 years" (Fig. 1 caption, p. 121). **Note:** The paper does not write the formula; with $R$ as length unit and $\eta$ the turbulent diffusivity it is $R^2/\eta$. The 2000-yr estimate is not justified; the 2002 PEPI paper (Section 10) later infers 10–15 kyr from the Sint-800 autocorrelation.
- Field unit: $a$ is "normalised" (p. 121); the normalisation is not stated, but Fig. 1(b) and the wells of $U$ at $\pm1$ (Fig. 2) imply the nonlinear equilibrium is $a=\pm1$. The absolute field ($\sim0.5$ G) cannot be predicted without a fully nonlinear model (p. 124).
- Conversions used: $T_r\approx100$ diffusion times is "typical" for 2–3 Myr ago (p. 121); 500 diffusion times $\approx1$ Myr; 50,000 diffusion times $\approx10^8$ yr (p. 123). **Note:** $\tau_c=5\times10^{-2}$ then corresponds to 100 yr.

### 4.4 Fluctuation model

Because the number of convective eddies on a latitude circle is finite, the longitudinal average leaves rapid residual fluctuations $\delta\alpha$ in latitude and time (p. 120):

$$
\frac{\delta\alpha}{\alpha_0}=\frac{fF}{\sqrt N},\qquad N=2N_c\sin\theta,\qquad N_c=\frac{\pi}{\theta_c},
$$

$f$: mean relative amplitude; $\theta_cR$: correlation length; $\tau_c$: correlation time of individual eddies; $F$: random number with zero mean and unit variance, updated once every $\tau_c$ independently for every eddy; $N$: eddies on a latitude circle; $N_c$: cells from pole to pole. All statistical properties depend only on $f^2\theta_c^2\tau_c$, the forcing parameter (p. 120). For the standard run the longitudinally averaged fluctuation strength is $\delta\alpha/\alpha_0\approx2$ (p. 121). **Note:** $\theta_c=\pi/10$ gives $N_c=10$ and $N=20\sin\theta$; at the equator $f/\sqrt N=6.4/\sqrt{20}\approx1.4$, larger toward the poles. The distribution of $F$ and the placement of cells on the $\theta$ grid are not specified.

### 4.5 Nonlinearity / saturation

Lorentz-force back-reaction is modelled by letting $\alpha$ decrease with the total energy of the fundamental mode ($\alpha$-quenching; Stix, 1992) (p. 120). In Eq. (2) this is the global factor $(1-a_0^2)$ on the supercritical term; in Eq. (3) it gives the side walls of $U$ (p. 122). The noise is "linear multiplicative" (p. 124). **Note:** With $\delta\alpha=0$ and only the fundamental mode, Eq. (2) gives $\dot a_0=\lambda_0a_0+(1-a_0^2)E_{00}a_0$; with $\lambda_0=0$ at marginal excitation the fixed points are $a_0=\pm1$ and the small-amplitude growth rate is $E_{00}=\Lambda$, consistent with the normalisation and with $S$. The exact quenching function is not given here.

### 4.6 Parameter values used

| Symbol | Meaning | Value(s) | Where |
|---|---|---|---|
| $C$ | dynamo number $R_\alpha R_\Omega$ | 100 | Fig. 1 caption, p. 121 |
| $C_{\rm crit}$ | marginal dynamo number | 58 | p. 121 |
| $f$ | mean relative fluctuation amplitude | 6.4 (Fig. 1); 6 and 7 (variants); 5.4 (Sint-800 match); 4 → 6.5 linear ramp (Fig. 3) | pp. 121, 123 |
| $\theta_c$ | correlation angle (length $\theta_cR$) | $\pi/10$ | p. 121 |
| $\tau_c$ | correlation time | $5\times10^{-2}$ diffusion times | p. 121 |
| $f^2\theta_c^2\tau_c$ | forcing parameter | not tabulated; **Note:** 0.202 ($f=6.4$), 0.144 ($f=5.4$), 0.079 → 0.209 (ramp), my arithmetic | — |
| time unit | turbulent diffusion time | 2000 yr (estimate) | p. 121 |
| run length | | 50,000 diffusion times (Fig. 1 run and Fig. 3 run) | pp. 121, 123 |
| $\delta\alpha/\alpha_0$ | averaged fluctuation strength | $\approx2$ | p. 121 |
| $\gamma$, $c^2$ | Eq. (4) parameters, Sint-800 fit | 5.2, 0.26 | Fig. 4 caption, p. 123 |
| $a_{\max}$ | VADM at the maximum of $p$ | $5.9\times10^{22}$ A m$^2$ | p. 123 |
| $\Lambda$ | growth rate of fundamental at $C=100$ | not given | p. 122 |

**Note:** The Fig. 3 caption gives only $f$ and the run length; $C$, $\theta_c$, $\tau_c$ for that run are presumably those of Fig. 1 but are not stated.

## 5. Methods

Numerics: no grid, time step or scheme is given. The Fig. 1 run lasted 50,000 diffusion times, of which the first 2000 are shown (p. 121); the Fig. 3 run covers 50,000 diffusion times with $f$ rising linearly from 4 to 6.5 (p. 123). $F$ is redrawn every $\tau_c$ independently per eddy (p. 120).

Analysis: (i) eigenmode decomposition at $C_{\rm crit}=58$ giving Eq. (2) "by standard techniques" (p. 121); (ii) averaging over the fluctuations to obtain Eq. (3), with the overtones folded into $S$ and $D$ (p. 122); (iii) stationary solution Eq. (4) and the escape-time formula for $T_r$ from standard Fokker–Planck theory (Van Kampen, 1992) (p. 122); (iv) run diagnostics: normalised $a$, its relative variance between reversals, the mean reversal time and the interval distribution (Poissonian, p. 121).

Observations: the Sint-800 record of the virtual axial dipole moment (VADM) over the past 800,000 yr (Guyodo and Valet, 1999), the VADM being the surface value of the $l=1$, $m=0$ harmonic of the poloidal field (p. 123). The paper bins it into 30 bins with error bars (Fig. 4) and quotes a VADM relative variance of 0.07 (p. 123); no other pre-processing is mentioned and the record's time resolution is not stated. Also used: secular variation 0.15–0.20 from records dated 2–3 Myr (Valet and Meynadier, 1993) and the chron-length trend of McFadden and Merrill (1995) (pp. 121, 123).

## 6. Results

### 6.1 Reference run (Fig. 1, p. 121)

With $C=100$, $C_{\rm crit}=58$, $f=6.4$, $\theta_c=\pi/10$, $\tau_c=5\times10^{-2}$ the field fluctuates rapidly but keeps its polarity for long times, with occasional sudden reversals and aborted reversals (e.g. $t=648$, 825, 863, 1048). Reversal intervals are Poissonian, as observed (Merrill et al., 1996). Between reversals the parity is predominantly antisymmetric (dipolar); at a reversal quadrupolar components come in and the field strength becomes small, as for the Earth (Hoffman, 1991; McFadden et al., 1991); the field tends to become quasi-periodic during a reversal because the oscillatory overtones dominate temporarily, possibly related to the rapid variability observed during reversals (Coe et al., 1995) (p. 121).

| Run | $f$ | $T_r$ (diffusion times) | relative variance of $a$ between reversals | Page |
|---|---|---|---|---|
| variant | 6 | $\approx200$ | 0.11 | 121 |
| reference | 6.4 | $\approx100$ | 0.16 | 121 |
| variant | 7 | $\approx50$ | 0.24 | 121 |
| Sint-800 match | 5.4 | 500 | 0.07 (inferred from VADM) | 123 |

The value 0.16 is "close to" the observed 0.15–0.20 at 2–3 Myr, and $T_r\approx100$ diffusion times "is typical for that period" (p. 121).

### 6.2 Mechanism and the bistable picture (pp. 121–122, Fig. 2)

The fluctuations perturb the fundamental mode and excite otherwise damped overtones; "reversals are always fast events, caused by the interaction of the dipole mode with the overtones" (p. 122); during a reversal $a_0\simeq0$ and the model evolves linearly; since the fluctuations are a manifestation of nonlinear turbulence, the amplitude variations and reversals are ultimately caused by nonlinear effects, and no unique physical cause for a reversal exists in the model (p. 124). In the Fokker–Planck picture the dipole amplitude is a heavily damped particle in $U(a)=\tfrac14\Lambda(1-a^2)^2$, symmetric about $a=0$ by $\pm B$ invariance, with a central hill from supercriticality and side walls from the nonlinearity; it sits most of the time in one well (normal or reversed polarity), random-walks around the minimum, and occasionally jumps over the hill: "the geodynamo is thus pictured as a bistable oscillator" (p. 122). $p(a)$ from Eq. (4) is bimodal with peaks near the wells and "turns out to describe the simulated and observed dipole distribution rather well" (p. 122). **Note:** For $\gamma=5.2$, $c^2=0.26$ the exponent in Eq. (4) is 5.55 and the peak is at $a\approx0.90$; integrating that form over one well gives a relative variance of about 0.09 (my evaluation), somewhat above the quoted 0.07, in line with the caption's remark that the fit is "not the optimal possibility".

### 6.3 Secular variation versus reversal rate (pp. 122–123)

From Eq. (4) one obtains e.g. the relative variance $\langle(a-\langle a\rangle)^2\rangle/\langle a\rangle^2$, "the relative secular variation in geophysical parlance"; $T_r$ from the escape-time integral depends exponentially on $f^2\theta_c^2\tau_c$ (p. 122). Numerically, the variability between reversals increases with forcing while $T_r$ decreases sharply, as in the geomagnetic record: a more strongly kicked particle oscillates more widely in its well and jumps more often. Tuning the forcing so that $T_r$ has the observed value, the model reproduces the VADM secular variation (cf. Fig. 1), and this explains the lower secular variation during long superchrons such as the late Cretaceous one (McFadden and Merrill, 1995) (pp. 122–123). **Note:** Using the four $(f,T_r)$ pairs above with $\theta_c=\pi/10$, $\tau_c=0.05$, $\ln T_r$ falls from 6.2 to 3.9 as the forcing parameter rises from 0.144 to 0.242, roughly linearly (slope about $-25$); $T_r$ halves per $\approx0.03$ of forcing parameter in this range. This is my fit to four points.

### 6.4 Ramped-forcing run (Fig. 3, p. 123)

A 50,000-diffusion-time ($\sim10^8$ yr) run with $f$ rising linearly from 4 to 6.5 shows a drastic decrease of the mean time between reversals and an approximately linear rise of the rms variability of $a$ from 17% to 41% (relative variance 0.03 to 0.17), in "qualitative agreement" with the Earth over the last $10^8$ yr. The observed fall of chron lengths from $10^7$ to $10^5$ yr is explained by a factor $\sim2.5$ rise of $f^2\theta_c^2\tau_c$, "presumably associated with small gradual changes in the mean properties of the convection in the outer core"; the sensitivity of the reversal rate to the outer boundary condition found by Glatzmaier et al. (1999) is read the same way (p. 123). **Note:** $(6.5/4)^2=2.64$. In the figure the first normal chron spans roughly the first 19,000 diffusion times [approximate reading].

### 6.5 Sint-800 comparison (p. 123, Fig. 4)

Argument for VADM $\propto a_0$: by axisymmetry $m=0$, and $a_0$ gives the time dependence of both the toroidal and the poloidal field of the fundamental mode; $l=1$ dominates $a_0$; $l=3$ contributes about 20% and higher $l$ are negligible; the next antisymmetric mode $a_2$, only temporarily excited and otherwise decaying, has mean amplitude about 10% of $a_0$; both fall off with height faster than the $l=1$ part, so the error is estimated below 10%, about the VADM observational error (p. 123). The VADM relative variance 0.07 corresponds to $f=5.4$, giving $T_r=500$ diffusion times ($\sim1$ Myr), consistent with the Brunhes chron, now 780,000 yr long. The binned Sint-800 distribution is matched by Eq. (4) within $1\sigma$ with $\gamma=5.2$, $c^2=0.26$, the $a$–VADM scale chosen so that the maxima coincide at $a_{\max}=5.9\times10^{22}$ A m$^2$; the fit "is not the optimal possibility but one that lies within $1\sigma$ of the Sint data and the numerical run with $f=5.4$". The match "supports, by implication, the underlying physics of the model" (p. 123).

## 7. Figures and tables

Four figures, no tables.

| Figure | Page | What is plotted | Run / parameters | Key content |
|---|---|---|---|---|
| Fig. 1(a) | 121 | Colour plot of field strength vs. time (0–2000 diffusion times) and colatitude (0–$\pi$) at a representative radius; red positive, blue negative | $C=100$, $C_{\rm crit}=58$, $f=6.4$, $\theta_c=\pi/10$, $\tau_c=5\times10^{-2}$ | Dipolar field between reversals; polarity blocks with sharp transitions |
| Fig. 1(b) | 121 | Normalised dipole amplitude $a$ ($-2$ to 2) vs. time | same | $a$ around $\pm1$ with excursions toward 0; aborted reversals at $t=648$, 825, 863, 1048; of order 15 polarity changes in 2000 diffusion times [approximate count, consistent with $T_r\approx100$]; relative variance 0.16 |
| Fig. 2 | 122 | $U(a)=\tfrac14\Lambda(1-a^2)^2$ (solid) and $p(a)$ from Eq. (4) (dashed) vs. $a\in[-2,2]$, with particle cartoon and jump arrow | $\gamma$, $c^2$ "corresponding to Fig. 1" (values not given) | Wells at $\pm1$, hill at 0; $p(a)$ peaks near $a\approx\pm0.85$ [approximate reading] |
| Fig. 3 | 123 | Polarity record (black normal, white reversed) vs. time, 0–50,000 diffusion times | $f$ ramped 4 → 6.5 | Chrons shorten drastically; rms variability 17% → 41% (relative variance 0.03 → 0.17) |
| Fig. 4 | 123 | Sint-800 VADM histogram, 30 bins, $0$–$10\times10^{22}$ A m$^2$; y-axis "No. of points / $10^{22}$ A m$^2$", 0–250; $1\sigma$ error bars; theoretical $p(a)$ overlaid | $\gamma=5.2$, $c^2=0.26$, $a_{\max}=5.9\times10^{22}$ A m$^2$; run with $f=5.4$ | Histogram peak about 230 at VADM $\approx5.5$–6, fitted curve peak about 190 near 5.9 [approximate readings]; fit within $1\sigma$ |

## 8. Comparison with observations

| Observable | Model | Data (as quoted) | Stated agreement | Page |
|---|---|---|---|---|
| Relative secular variation | 0.16 ($f=6.4$) | 0.15–0.20 at 2–3 Myr (Valet and Meynadier, 1993) | "close" | 121 |
| Mean time between reversals | $\approx100$ diffusion times ($\approx2\times10^5$ yr) | "typical for that period" | agreement | 121 |
| Reversal-interval distribution | Poissonian | Poissonian (Merrill et al., 1996) | agreement | 121 |
| Field during reversals | small, quadrupolar, quasi-periodic | quadrupolar transitional fields; rapid variability | qualitative | 121 |
| Reversal duration | "always fast" (no number) | a few thousand years | qualitative | 119, 122 |
| SV–reversal-rate relation | SV up, $T_r$ down with forcing | lower SV in superchrons | explained | 122–123 |
| Chron-length trend | factor $\sim2.5$ forcing rise (Fig. 3) | $10^7\to10^5$ yr over last $10^8$ yr | "qualitative agreement" | 123 |
| Sint-800 VADM variance and chron | $f=5.4$: 0.07, $T_r=500$ diffusion times $\sim1$ Myr | 0.07; Brunhes 780,000 yr | consistent | 123 |
| VADM distribution | Eq. (4), $\gamma=5.2$, $c^2=0.26$ | 30-bin Sint-800 histogram | within $1\sigma$ | 123 |
| Absolute field | not predicted | $\sim0.5$ G | needs fully nonlinear model | 124 |

## 9. Assumptions, caveats and limitations

Stated by the authors:
- The model is "clearly highly simplified": 1D in $\theta$, constant radial shear, $\alpha\propto\cos\theta$, mean-field approximation (pp. 120, 124).
- The method requires a non-oscillatory, predominantly dipolar field, a dominant stationary dipole mode and a source of fluctuations (pp. 119, 124).
- Treating nonlinear turbulence as noise is valid on the eddy-turnover time; on longer time scales nonlinearities remain essential (p. 124).
- Not all properties follow; the absolute field magnitude needs a fully nonlinear model (p. 124).
- VADM $\propto a_0$ only "within the error bounds" ($<10\%$) (p. 123); the Fig. 4 fit is "not the optimal possibility" (p. 123).
- $S$ and the right-hand side of Eq. (4) are approximations (p. 122); the origin of the factor $\sim2.5$ forcing increase is only "presumably" slow convective change (p. 123).

**Note:** Further points for implementers: the 2000-yr diffusion time underlies every year conversion, and the 2002 PEPI paper's 10–15 kyr would rescale $T_r$, $\tau_c$ and the $10^8$-yr run length by a factor 5–7.5; "relative variance between reversals" is not defined operationally (how transitions are excluded); the 0.15–0.20 (2–3 Myr) and 0.07 (Sint-800) values are matched by changing $f$ (6.4 vs 5.4), a tuning rather than a prediction; Fig. 2's $\gamma$, $c^2$ are unspecified.

## 10. Relation to the other two papers in this set

This is paper 2 of 3 chronologically (received 22 December 2000) and is the compact summary.

- Hoyng, Ossendrijver and Schmitt (2001), GAFD 94, 263–314, "The geodynamo as a bistable oscillator": the long paper with full derivations. This paper defers to it twice: "Full details of the model are given in Hoyng et al. (2001)" (p. 120) and "Details of the derivation of this [Eq. 2] as well as the other equations of this section can be found in Hoyng et al. (2001)" (p. 121). Its abstract states the same picture (bistable potential $\propto(1-a^2)^2$, reversals induced solely by the overtones, theoretical expressions for the amplitude distribution, between-reversal variance and mean reversal rate, Poisson intervals robust to the type of nonlinearity) and adds a claim absent here: a present reversal rate of once per $(2$–$3)\times10^5$ yr implies a secular variation of the axial dipole moment of $\sim15\%$. The GAFD abstract also reports agreement of the theoretical distribution with Sint-800, so the Fig. 4 comparison is shared; the Fig. 3 ramped run and the chron-length argument appear, from the abstracts, to be specific to this paper.
- Hoyng, Schmitt and Ossendrijver (2002), PEPI 130, 143–157, "A theoretical analysis of the observed variability of the geomagnetic dipole field": the detailed Sint-800 analysis. Per its abstract it predicts the VADM distribution and autocorrelation function, uses the autocorrelation to determine the turbulent diffusion time $\tau_d=R^2/\beta\approx10$–15 kyr, reproduces the secular variation and the mean time between reversals, finds the reversal duration a factor $\sim2$ too long, and shows mean-field theory to be self-consistent. It extends Section 5 of the present paper, and its diffusion-time result supersedes the 2000-yr estimate assumed here.

GAFD 2001 = derivations; this paper = compact summary with the Sint-800 comparison; PEPI 2002 = detailed data analysis with autocorrelation and diffusion time.

## 11. Implementation guide for coding agents

System to solve: the 1D PDEs of Eq. (1) for $P(\theta,t)$, $T(\theta,t)$ on $0\le\theta\le\pi$, i.e. $\partial_tP=\Delta P+[(\alpha+\delta\alpha)/\alpha_0]T$ and $\partial_tT=C\,\partial_\theta(\sin\theta P)+\Delta T$ [operator ordering as read from the printed matrix; Note in 4.2], with $\alpha/\alpha_0=\cos\theta$, time in diffusion times, length in $R$. Not specified in this paper: the explicit $\theta$-form of $\Delta$, the pole boundary conditions, the grid, the time stepper, and the exact quenching function (only "decreasing function of the total energy of the fundamental mode", appearing as $(1-a_0^2)$ in Eq. 2). Take these from Schmitt and Schüssler (1989) and Hoyng et al. (2001).

Eigenvalue problem: compute the eigenmodes of the $\delta\alpha=0$ operator at $C_{\rm crit}=58$ (p. 121). For $C>0$ the first excited mode is the non-oscillatory, equatorially antisymmetric dipole (the fundamental, index 0); all others are oscillatory and damped at this $C$ (p. 120). $a_2$ is "the next antisymmetric dynamo mode" (p. 123), which suggests the index runs over dipolar and quadrupolar modes in order of excitation, $a_1$ then being quadrupolar [uncertain; not stated]. Normalise so that $|a_0|=1$ at equilibrium (Fig. 1b, Fig. 2). $E_{ik}$, $V_{ik}$ are not given here.

Stochastic $\alpha$: $N_c=\pi/\theta_c$ cells from pole to pole (10 for $\theta_c=\pi/10$), $N=2N_c\sin\theta$ eddies per latitude circle, $\delta\alpha/\alpha_0=fF/\sqrt N$ with $F$ zero-mean unit-variance, redrawn every $\tau_c$ independently per eddy (p. 120). The distribution of $F$ and the mapping of cells to grid points are not specified. Checks: $\delta\alpha/\alpha_0\approx2$ overall for $f=6.4$ (p. 121); statistics invariant under changes of $f$, $\theta_c$, $\tau_c$ at fixed $f^2\theta_c^2\tau_c$ (p. 120).

Standard run: $C=100$, $f=6.4$, $\theta_c=\pi/10$, $\tau_c=0.05$, 50,000 diffusion times (p. 121); variants $f=6$, 7, 5.4 and a linear ramp 4 → 6.5 over 50,000 diffusion times (pp. 121, 123).

Diagnostics as defined by the paper:
- $a(t)=a_0(t)$, normalised fundamental-mode amplitude, identified with the VADM (pp. 121, 123).
- Relative variance between reversals $\langle(a-\langle a\rangle)^2\rangle/\langle a\rangle^2$ = relative secular variation (pp. 121–122); its square root is the "rms variability" (17–41% ↔ 0.03–0.17, p. 123). The rule for excluding transitions is not given.
- Reversal: $a$ crosses zero and settles in the other well; aborted reversal: $a$ approaches zero and returns ($t=648$, 825, 863, 1048 in Fig. 1). No threshold or minimum dwell time is specified; choose and document one.
- $T_r$: mean time between reversals; intervals should be Poissonian (p. 121).
- $p(a)$: histogram of $a$ (30 bins for Sint-800, p. 123); fit Eq. (4); compare $T_r$ with the escape-time integral if $D(a)$ is modelled (not given here).

Validation targets ($C=100$, $\theta_c=\pi/10$, $\tau_c=0.05$):

| $f$ | $T_r$ (diffusion times) | relative variance | Page |
|---|---|---|---|
| 5.4 | 500 | 0.07 | 123 |
| 6 | $\approx200$ | 0.11 | 121 |
| 6.4 | $\approx100$ | 0.16 | 121 |
| 7 | $\approx50$ | 0.24 | 121 |
| 4 → 6.5 ramp | drastic decrease over 50,000 | 0.03 → 0.17 | 123 |

Sint-800 target: Eq. (4) with $\gamma=5.2$, $c^2=0.26$, peak at VADM $5.9\times10^{22}$ A m$^2$ (p. 123).

Pitfalls: 1 diffusion time = 2000 yr in this paper only (later 10–15 kyr), so keep internal times in diffusion times; $C_{\rm crit}=58$ is the eigenbasis reference and $C=100$ the run value, with $\Lambda$ not given; $U(0)=\Lambda/4$ matters only relative to $D(0)$ through $\gamma$; Eq. (4) is a stationary one-variable approximation with the overtones averaged into $S$ and $D$ and should not be expected to describe reversal transients; relative variances are "between reversals" and would be inflated by including transitions; the Fig. 3 percentages are square roots of relative variances.

## 12. Symbol glossary

| Symbol | Meaning | Dimensionless / units | First defined |
|---|---|---|---|
| $R$ | radius of the dynamo; length unit | length | p. 120 |
| $\theta$ | colatitude | rad | p. 120 |
| $T$, $P$ | toroidal field; poloidal vector potential | normalised | Eq. (1), p. 120 |
| $\Delta$ | diffusion operator | dimensionless | Eq. (1), p. 120 |
| $\alpha$, $\alpha_0$, $\delta\alpha$ | $\alpha$-effect $\alpha_0\cos\theta$, its amplitude, its fluctuation | $\alpha/\alpha_0$ dimensionless | p. 120 |
| $d\Omega/dr$ | constant radial shear | 1/(length·time) | p. 120 |
| $\eta$ | turbulent magnetic diffusivity | length$^2$/time | p. 120 |
| $C$, $R_\alpha$, $R_\Omega$ | dynamo number; $\alpha$ and $\Omega$ Reynolds numbers | dimensionless | p. 120 |
| $C_{\rm crit}$ | marginal dynamo number (58) | dimensionless | p. 121 |
| $f$ | mean relative amplitude of $\delta\alpha$ | dimensionless | p. 120 |
| $\theta_c$, $\tau_c$ | eddy correlation angle and time | rad; diffusion times | p. 120 |
| $F$ | random number, zero mean, unit variance | dimensionless | p. 120 |
| $N$, $N_c$ | eddies per latitude circle; cells pole to pole | integers | p. 120 |
| $f^2\theta_c^2\tau_c$ | forcing parameter | dimensionless | p. 120 |
| $a_i$, $a_0$ ($=a$), $a_2$ | mode amplitudes; fundamental dipole; next antisymmetric mode | normalised ($a_0=\pm1$ at equilibrium) | Eq. (2), p. 121; p. 123 |
| $\lambda_i$ | eigenvalues at $C_{\rm crit}$ | 1/diffusion time | Eq. (2), p. 121 |
| $E_{ik}$, $V_{ik}(t)$ | overlap integrals of eigenmodes | 1/diffusion time | Eq. (2), p. 121 |
| $p(a)$ | probability distribution of $a$ | 1/(unit of $a$) | p. 121 |
| $S$, $D$ | drift, diffusion coefficients | 1/diffusion time ($a$ normalised) | Eq. (3), p. 122 |
| $\Lambda$ | growth rate of the fundamental mode at $C$ | 1/diffusion time | p. 122 |
| $U(a)$ | bistable potential $\tfrac14\Lambda(1-a^2)^2$ | 1/diffusion time ($a$ normalised) | Fig. 2, p. 122 |
| $\gamma$, $c^2$ | Eq. (4) constants measuring $U(0)$, $D(0)$ | dimensionless | Eq. (4), p. 122 |
| $T_r$ | mean time between reversals | diffusion times | pp. 121–122 |
| VADM, $a_{\max}$ | virtual axial dipole moment; VADM at the peak of $p$ | A m$^2$ ($10^{22}$ A m$^2$ in Fig. 4) | p. 123 |
| $l$, $m$ | spherical-harmonic degree, order | integers | p. 123 |

## 13. Section-by-section outline

| Section | Pages | Content |
|---|---|---|
| Abstract | 119 | Fokker–Planck method for the dipole moment; needs a non-oscillatory dipolar field and helicity fluctuations; bistable potential; reproduces the SV–reversal-rate relation and the Sint-800 distribution. |
| 1. Introduction | 119–120 | Observational facts; 3D simulations and inner-core stabilisation; kinematic reversal schemes; objective of long-term statistics via mean-field theory with differential rotation, helicity and large fluctuations; $\sim10^8$-yr runs. |
| 2. Dynamo model | 120 | 1D $\alpha\Omega$ model, Eq. (1), dynamo number, mode structure, $\alpha$-quenching, longitudinal averaging, construction of $\delta\alpha$, forcing parameter. |
| 3. Numerical results | 120–121 | Fig. 1 run; Poisson intervals; dipolar parity between reversals, quadrupolar and quasi-periodic behaviour during reversals. |
| 4. Statistical analysis | 121–122 | Eigenmode expansion at $C_{\rm crit}=58$, Eq. (2); Fokker–Planck Eq. (3); drift and diffusion; potential and Fig. 2; stationary Eq. (4). |
| 5. Comparison with geomagnetic record | 122–123 | Relative variance; $T_r$ formula; SV vs. $T_r$ trend; superchrons; Fig. 3 ramp and chron-length trend; Glatzmaier et al. (1999); Sint-800 and the $a_0\leftrightarrow$VADM argument; $f=5.4$; Fig. 4 fit. |
| 6. Discussion | 124 | Linear multiplicative noise; linear evolution during reversals; noise idealisation valid on the eddy-turnover time; no unique cause of a reversal; absolute field not predicted; method applicable to 2D models as a prototype. |
| References | 124 | 26 entries. |

## 14. Key references cited

- Hoyng, Ossendrijver, Schmitt (2001), GAFD 94, 263–314: full model details and derivations of Eqs. (2)–(4) (pp. 120, 121).
- Schmitt and Schüssler (1989), A&A 223, 343–351: the 1D thin-layer dynamo and its mode structure (p. 120).
- Hoyng (1993), A&A 272, 321–339; Ossendrijver, Hoyng, Schmitt (1996), A&A 313, 938–948: random $\alpha$ fluctuations in the solar dynamo (p. 120).
- Guyodo and Valet (1999), Nature 399, 249–252: the Sint-800 VADM record (pp. 119, 123).
- Valet and Meynadier (1993), Nature 366, 234–238: secular variation 0.15–0.20 at 2–3 Myr (p. 121).
- McFadden and Merrill (1995), PEPI 91, 253–260: lower SV in superchrons; chron-length trend (pp. 122–123).
- Merrill, McElhinny, McFadden (1996), The Magnetic Field of the Earth: Poisson reversal statistics (pp. 119, 121).
- Glatzmaier and Roberts (1995), Nature 377, 203–209: 3D reversal simulation (p. 119).
- Glatzmaier, Coe, Hongre, Roberts (1999), Nature 401, 885–890: reversal frequency vs. boundary condition (p. 123).
- Hollerbach and Jones (1993), Nature 365, 541–543: stabilising inner core (p. 119).
- Parker (1969), ApJ 158, 815–827; Levy (1972), ApJ 171, 635–642: kinematic reversal schemes (p. 120).
- Van Kampen (1992), Stochastic Processes in Physics and Chemistry: Fokker–Planck theory and escape rates (p. 122).
- Hoffman (1991), Nature 354, 273–277; McFadden et al. (1991), JGR 96, 3923–3933: quadrupolar transitional fields (p. 121).
- Coe, Prévot, Camps (1995), Nature 374, 687–692: rapid change during a reversal (p. 121).
- Stix (1992), A&A 20, 9–12 (as printed): $\alpha$-quenching (p. 120); Soward (1991), GAFD 62, 191–209: $\alpha\Omega$ classification (p. 120).

## 15. Open questions and unresolved points

- The cause of the inferred factor $\sim2.5$ rise of the forcing parameter over $10^8$ yr is only presumed (slow change of mean convection) (p. 123).
- Application of the method to more realistic (2D) models is proposed, not done; the authors argue it carries over because Eq. (2) involves only the eigenvalues and certain overlap integrals of the eigenfunctions, the only requirements being a source of fluctuations and a dominant stationary dipole mode (p. 124).
- The absolute field magnitude ($\sim0.5$ G) is out of reach (p. 124).
- The Fig. 4 fit is non-optimal; no best-fit parameters or goodness-of-fit statistic are reported (p. 123).
- Reversals are "always fast" but no model duration is quoted against the observed few kyr (pp. 119, 122). **Note:** The 2002 PEPI paper later finds the model duration a factor $\sim2$ too long.
- The 2000-yr diffusion time is an unexamined estimate on which all year conversions rest (p. 121). **Note:** The 2002 PEPI paper's autocorrelation analysis gives 10–15 kyr.
- Reconciling the 0.16 (2–3 Myr) and 0.07 (Sint-800) secular-variation levels other than by changing $f$ is not discussed.
