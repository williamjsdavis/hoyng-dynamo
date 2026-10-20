---
title: "A theoretical analysis of the observed variability of the geomagnetic dipole field"
authors: ["P. Hoyng", "D. Schmitt", "M.A.J.H. Ossendrijver"]
year: 2002
journal: "Physics of the Earth and Planetary Interiors"
volume: 130
pages: "143–157"
doi: "10.1016/S0031-9201(02)00004-3"
source_pdf: papers/geophysical_1_Breno.pdf
pdf_pages: 15
page_mapping: "PDF page 1 = journal page 143; journal page = PDF page + 142"
series: "Hoyng–Schmitt–Ossendrijver bistable geodynamo series, paper 3 of 3 (chronological)"
keywords: ["Geodynamo", "Reversals", "Secular variation", "Sint-800 record", "Helicity fluctuations", "Bistable oscillator"]
extraction_notes: "Vector PDF with a good text layer; all 17 numbered equations (1–15, A.1, A.2), both tables and all figure captions were transcribed from 2x page images and 4x zoomed crops, not from the text layer (which flattens fractions, sub/superscripts and integral limits). No equations were unreadable."
---

# A theoretical analysis of the observed variability of the geomagnetic dipole field

## 1. One-paragraph summary

The paper confronts the authors' stochastic αΩ mean-field geodynamo (developed in Hoyng, Ossendrijver & Schmitt 2001, "HOS", and summarised in Schmitt, Ossendrijver & Hoyng 2001, "SOH") with the Sint-800 record of the virtual axial dipole moment (VADM) of the past 800 kyr (Guyodo & Valet 1999). In the model the α-effect has a rapidly fluctuating part, which excites overtones and produces a variable dipole amplitude $a$ with occasional reversals; the amplitude obeys a Fokker–Planck equation equivalent to a heavily damped particle in a bistable potential $U\simeq\tfrac14\Lambda(1-a^2)^2$. The stationary amplitude distribution $p(a)$ (Eq. 10, improved in Eq. A.1) is fitted by $\chi^2$ to the binned Sint-800 VADM and to simulated $a(\tau)$; fits at the 1–1.4σ level are obtained (Table 1, p. 149), and the secular dipole variation $dv$ (Eq. 12) of the data, $0.071\pm0.003$, is matched (model 0.070). A new autocorrelation technique—overlaying the autocorrelation of the Sint-800 series with that of the simulated amplitude—yields the turbulent diffusion time $\tau_d=R^2/\beta=11\pm2$ kyr (p. 153), i.e. $\tau_d$ is about 10–15 kyr. With this $\tau_d$ the mean time between reversals ($\sim6$ Myr) is too long for the reference run; a second run fixes that at the cost of a reversal duration of 8–30 kyr, a factor ~2 too long. Mean-field theory is argued to be self-consistent ($\beta/\eta\sim20$–30, $C_\alpha/C_\Omega\sim10^{-3}$, p. 155).

## 2. Key takeaways

- The Sint-800 VADM record (799 points, one per 1000 yr, individual errors ~10%, one reversal at 780 kyr) has secular dipole variation $dv=\langle(|a|-\langle|a|\rangle)^2\rangle/\langle|a|\rangle^2=0.071\pm0.003$; the uncertainty was obtained by Monte-Carlo resampling $s_i\to s_i+e_iN(0,1)$ (pp. 144, 148, 151).
- The dipole amplitude $a$ behaves like a strongly damped particle in the bistable potential $U\simeq\tfrac14\Lambda(1-a^2)^2$ (Eq. 8), with drift $S\simeq\Lambda(1-a^2)a$ and diffusion $D\simeq\Delta_0a^2+\Delta_1$ (Eq. 7); the stationary distribution is $p(a)\propto(1+a^2/c^2)^{\gamma c^2+\gamma-1}\exp(-\gamma a^2)$ with $\gamma=\Lambda/\Delta_0$, $c^2=\Delta_1/\Delta_0$ (Eqs. 10–11, p. 147).
- The optimal fit of Eq. (10) to Sint-800 has $\gamma=6.3$, $c^2=0.32$, peak VADM $A_m=5.98\times10^{22}$ A m², $\chi^2$ at 1σ (Table 1 case 1, p. 149); the fit to the simulated $a(\tau)$ of the reference run gives $\gamma=8$, $c^2=0.50$ (case 2). A compromise $\gamma=7$, $c^2=0.45$ fits Sint-800 at 1.1σ (case 3); with the improved Eq. (A.1), $\gamma=4.8$, $c^2=0.22$ fits at 1.4σ and gives $dv=0.071$ (case 4).
- Fits are individually sensitive to $\gamma$ and $c^2$ but degenerate: increasing $\gamma$ can be offset by increasing $c^2$, so a wide range of ($\gamma$, $c^2$) is compatible with Sint-800 (p. 150).
- Reference run (Fig. 3; run #2 of HOS): $kR=0.5$, $C_m=100$, $\Lambda=0.527$, $f^2\tau_c/N_c^2=1.46\times10^{-2}$ ($f=5.4$, $\tau_c=0.05$, $N_c=10$), $\Delta\tau=0.01$, length $5\times10^4\tau_d$, 88 reversals, $\langle T_r\rangle=570\pm61\,\tau_d$, $dv=0.070$ (p. 146).
- Linearising the Langevin equation about $a=1$ gives an autocorrelation $\langle\tilde a(\tau)\tilde a(\tau+s)\rangle\propto\exp(-2\Lambda s)$ (Eq. 15, p. 152): the observed variability time scale is the roll-down time $\Lambda^{-1}$ of the potential well, independent of $\tau_c$ as long as $\Lambda\tau_c\ll1$.
- Matching the autocorrelation functions of Sint-800 and of the rectified simulated amplitude (running average $\tau_{av}\sim5$ kyr) gives $\tau_d=R^2/\beta=11\pm2$ kyr for the reference run and $16\pm4$ kyr for the second run; this is a model-dependent measurement of $\Lambda/\tau_d$ (pp. 152–153).
- With $\tau_d\sim11$ kyr the reference run implies $\langle T_r\rangle\approx6$ Myr (Table 2: 5–7.5 Myr), outside the 95% confidence interval 220 kyr $<\langle T_r\rangle<$ 3.3 Myr inferred from the single reversal in 800 kyr (pp. 151, 153).
- Second run (Fig. 8): $kR=0.3$, $C_m=120$, $\Lambda=0.815$, $f^2\tau_c/N_c^2=2.45\times10^{-2}$ ($f=7$): $\langle T_r\rangle=(113\pm5)\tau_d\approx1.8$ Myr (inside the interval), $dv=0.081$, but the (A.1) fit to Sint-800 is 5σ and the reversal duration is 8–30 kyr (pp. 153–154).
- Either $\langle T_r\rangle$ (case 1) or the reversal duration (case 2) is too long; a $\tau_d$ smaller by a factor 2–3 would remove the discrepancy (p. 154).
- Theoretical $\langle T_r\rangle$ from Eq. (A.2) is almost a factor 2 larger than the measured $570\pm61\,\tau_d$ (except case 3); the authors attribute this to (A.2) being approximate (p. 151).
- Mean-field self-consistency: with $R=3500$ km and $\eta\sim1.3\times10^4$ cm² s⁻¹, $\beta/\eta\sim30$ (case 1) and $\sim20$ (case 2); $\alpha_m\sim-2.5\times10^{-4}$ cm s⁻¹ and $C_\alpha/C_\Omega\sim10^{-3}$ justify the αΩ approximation (p. 155).

## 3. Context and motivation

The paper starts from the paleomagnetic facts (p. 143): the geomagnetic field has always been approximately dipolar and aligned with the rotation axis; reversals are fast events lasting $\lesssim10$ kyr; the intervals $T_r$ between reversals are Poisson distributed (printed as $\propto\exp(T_r/\langle T_r\rangle)$—**Note:** the minus sign is evidently missing in print); and $\langle T_r\rangle$ has decreased from $\sim10^7$ yr in the Cretaceous to $(2$–$3)\times10^5$ yr during the past 10 Myr (Merrill et al. 1996; Valet & Meynadier 1993). On shorter time scales the dipole amplitude fluctuates considerably (Fig. 1, the Sint-800 record of Guyodo & Valet 1999).

State of the art as described by the authors (p. 144): three-dimensional simulations (Glatzmaier & Roberts 1995a,b; Kageyama & Sato 1997; Kuang & Bloxham 1997; Olson et al. 1999; Glatzmaier et al. 1999; Christensen et al. 1999) reproduce a largely dipolar field with $B\simeq0.5$ G at the surface, sudden reversals and rapid variability, but run at Ekman numbers $E\sim10^{-6}$ whereas the core is believed to have $E\sim10^{-15}$ (see Zhang & Schubert 2000), and they are extremely expensive. Many simpler models reproduce the observations qualitatively (Mazaud & Laj 1989; Narteau et al. 2000), but a *quantitative* comparison of model predictions with observed geodynamo statistics "has received little attention"; the paper aims to bridge that gap.

Prior stochastic ideas the model builds on (p. 144): Parker (1969) and Levy (1972) considered a burst in helicity followed by a jump in the mean flow (the Parker–Levy mechanism; see also Gibbons 1998); Olson (1983) found reversals due to jumps in the α-effect; Narteau et al. (2000) due to a steady level of helicity fluctuations. Because the observed distribution of intervals between reversals is Poissonian, the authors regard it as "almost unavoidable that reversals are induced by some kind of rapid variability" (p. 144). Ossendrijver et al. (1996) and Hoyng (1993, 1996) used the same kind of α-fluctuations to explain the phase–amplitude correlation of the solar cycle (periodic dynamos react differently, pp. 145–146). The two key elements of the present model are (1) a stationary level of random fluctuations in α over very long times, and (2) numerical *and* theoretical evaluation of the statistics of the dipole amplitude via a Fokker–Planck equation (p. 144).

What is new in this paper relative to HOS and SOH: a systematic, quantitative confrontation with Sint-800 ("a systematic comparison of the model with the observations had not yet been undertaken", p. 148): the binning and error analysis of the VADM distribution, the $\chi^2$ fitting procedure and its degeneracies, an improved $p(a)$ (Appendix A), the autocorrelation method for measuring $\tau_d$, an attempt to satisfy all constraints simultaneously (Section 6 of the paper), and the mean-field self-consistency estimate.

## 4. The model

### 4.1 Physical setup and geometry

- Mean-field (azimuthally averaged) axisymmetric dynamo in the liquid outer core, of αΩ type (Fig. 2, p. 144; "on general grounds", see also Richards 2000).
- The mean flow is a differential rotation $\boldsymbol v=\boldsymbol\Omega\times\boldsymbol r=\Omega r\sin\theta\,\boldsymbol e_\varphi$ with $d\Omega/dr$ constant (p. 145). Differential rotation makes toroidal from poloidal field; the α-effect of the helical convection makes poloidal from toroidal field (Fig. 2 caption). **Note:** the Fig. 2 caption as printed says differential rotation "generates new poloidal field from toroidal field", an evident slip for toroidal from poloidal.
- The radial dependence of the field is modelled as a spherical wave $\propto\exp(ikr)/kr$, so that $1/k$ measures the radial length scale and $r$ becomes a parameter set to a representative radial position $R$. The model is one-dimensional in colatitude $\theta$ (plus time $t$) and "catches the dynamo action at that position $R$ in the outer core, with crude allowance for radial turbulent diffusion" (p. 145). Details are in HOS.
- Turbulent diffusion coefficient $\beta$ is constant and $\gg$ the molecular resistivity $\eta$ (p. 145).
- The fundamental mode must be a non-periodic (stationary) dipole, which restricts $0<kR<1$ and $C_m>0$ (p. 145).

### 4.2 Governing equations

Mean-field induction equation (Moffatt 1978; Krause & Rädler 1980), p. 145:

$$\frac{\partial\langle\boldsymbol B\rangle}{\partial t}=\nabla\times\left[\boldsymbol v\times{}+\alpha-(\beta+\eta)\nabla\times\right]\langle\boldsymbol B\rangle \qquad\text{(Eq. 1)}$$

Here $\langle\cdot\rangle$ is the azimuthal average, so $\langle\boldsymbol B\rangle$ is the axisymmetric part of $\boldsymbol B$; the bracket is an operator acting on $\langle\boldsymbol B\rangle$ (i.e. $\nabla\times(\boldsymbol v\times\langle\boldsymbol B\rangle)+\nabla\times(\alpha\langle\boldsymbol B\rangle)-\nabla\times[(\beta+\eta)\nabla\times\langle\boldsymbol B\rangle]$); $\alpha$ is the dynamo coefficient. **Note:** the print shows "$\boldsymbol v x$" with an italic $x$; it stands for the cross product.

Dynamo coefficient with mean and fluctuating parts, p. 145:

$$\alpha=\alpha_m(1-q)\cos\theta+\delta\alpha(\theta,\tau) \qquad\text{(Eq. 2)}$$

$\alpha_m$ is the amplitude of the steady part, $q\propto a_0^2$ with $a_0$ the amplitude of the fundamental mode (global α-quenching, Brandenburg et al. 1989: the steady part of α decreases with the energy of the fundamental mode, which is approximately the total magnetic energy), $\delta\alpha$ the rapidly fluctuating part, $\tau$ dimensionless time.

Fluctuating part, p. 145:

$$\delta\alpha\propto\frac{fF(\theta,\tau)}{\sqrt{2N_c\sin\theta}} \qquad\text{(Eq. 3)}$$

$f$ gauges the strength of the fluctuations, $F$ is a random function of position and time with zero mean and unit variance, $N_c$ sets the cell size (see 4.4).

Dynamo number, p. 145:

$$C_m=\frac{\alpha_mR^4}{\beta^2}\frac{d\Omega}{dr} \qquad\text{(Eq. 4)}$$

Mode equations (p. 146). The field is expanded in eigenmodes of the unperturbed dynamo equation at marginal excitation (eigenvalues $\lambda^i$; $\lambda^0=0$, $\mathrm{Re}\,\lambda^i<0$ for $i\ge1$):

$$\frac{da_i}{d\tau}=\lambda^ia_i+(1-a_0^2)E_{ik}a_k+V_{ik}(\tau)a_k \qquad\text{(Eq. 5)}$$

with summation over repeated lower indices; $a_i$ is the amplitude of mode $i$, normalised to the amplitude $a_0$ of the dipole mode in non-linear equilibrium; the first term is the linear evolution, the second the supercritical excitation and α-quenching, the third the α-fluctuations that cause all variability; $E_{ik}$ and $V_{ik}$ are overlap integrals of eigenmodes $i$ and $k$ (the factor $s$ of HOS has been absorbed in $E_{ik}$, footnote 3).

Fokker–Planck equation for the probability distribution $p(a,\tau)$ of the dipole amplitude (p. 146), obtained by averaging Eq. (5) over the fluctuations:

$$\frac{\partial p}{\partial\tau}=-\frac{\partial}{\partial a}Sp+\frac12\frac{\partial^2}{\partial a^2}Dp,\qquad a\equiv a_0 \qquad\text{(Eq. 6)}$$

with drift and diffusion coefficients, p. 146:

$$S\simeq\Lambda(1-a^2)a+\text{correction term},\qquad D\simeq\Delta_0a^2+\Delta_1 \qquad\text{(Eq. 7)}$$

$\Lambda$ is the linear growth rate of the fundamental mode at small amplitude ($a\ll1$); the correction term is treated in Appendix A; $\Delta_0$, $\Delta_1$ are constants that depend on $kR$, $C_m$ and $f^2\tau_c/N_c^2$ "but the functional dependence is unfortunately unknown" (p. 147).

Bistable potential (p. 147), with $-dU/da=S$ and the correction term ignored:

$$U\simeq\tfrac14\Lambda(1-a^2)^2 \qquad\text{(Eq. 8)}$$

Stationary solution of Eq. (6), p. 147:

$$p(a)\propto\frac1D\exp\left(\int_0^a\frac{2S}{D}\,da'\right) \qquad\text{(Eq. 9)}$$

which, ignoring the correction term, reduces to

$$p(a)\propto\left(1+\frac{a^2}{c^2}\right)^{\gamma c^2+\gamma-1}\exp(-\gamma a^2) \qquad\text{(Eq. 10)}$$

with the two dimensionless parameters

$$\gamma=\frac{\Lambda}{\Delta_0},\qquad c^2=\frac{\Delta_1}{\Delta_0} \qquad\text{(Eq. 11)}$$

described as "dimensionless measures of $U(0)$ and $D(0)$" (p. 147). $p(a)$ has two peaks near $a=\pm1$ (Fig. 4); its maximum is at $a_m=(1-1/\gamma)^{1/2}$ (p. 149).

Secular dipole variation (p. 148):

$$dv=\frac{\langle(|a|-\langle|a|\rangle)^2\rangle}{\langle|a|\rangle^2} \qquad\text{(Eq. 12)}$$

i.e. the variance of the dipole amplitude relative to the square of its mean, with $\langle|a|\rangle=\int_0^\infty ap(a)\,da\big/\int_0^\infty p(a)\,da$ etc.

Langevin (Itô) form of Eq. (6) (p. 151; Gardiner 1990): $\dot a=S+DL(\tau)$ with $L$ a rapidly varying Langevin force, i.e.

$$\frac{da}{d\tau}\simeq\Lambda(1-a^2)a+DL(\tau) \qquad\text{(Eq. 13)}$$

and for small oscillations about $a=1$, $a=1+\delta a$:

$$\frac{d\delta a}{d\tau}\simeq-2\Lambda\,\delta a+\Delta_0L(\tau) \qquad\text{(Eq. 14)}$$

whence, after rectifying the data to eliminate reversals and subtracting the mean (series $\tilde a$), p. 152:

$$\langle\tilde a(\tau)\tilde a(\tau+s)\rangle\propto\exp(-2\Lambda s) \qquad\text{(Eq. 15)}$$

$s$ is the time lag in units of $\tau_d$. **Note:** Eqs. (13)–(14) are printed with noise amplitude $D$ and $\Delta_0$ respectively; in the standard Itô correspondence to Eq. (6) the noise amplitude is $\sqrt D$, and $D(1)=\Delta_0+\Delta_1$. Neither point affects Eq. (15), whose decay rate $2\Lambda$ comes only from the drift.

Improved distribution (Appendix A, p. 156), obtained from Eq. (9) with $S=\Lambda(1-a^2)a+(\mu-\Lambda)(1-a^2)a^3$:

$$p(a)\propto\frac{1}{a^2+c^2}\exp\left\{\gamma\int_0^{a^2}\frac{1+(\mu/\Lambda-1)x}{x+c^2}(1-x)\,dx\right\} \qquad\text{(Eq. A.1)}$$

$\mu$ is the "bare" supercritical growth-rate parameter of HOS (HOS Eq. 33: $S=\mu(1-a^2)a+A(1-a^2)a=\Lambda(1-a^2)a$), and $\mu/\Lambda$ is known for a given model: 1.45 for the run of Fig. 3 and 1.82 for the run of Fig. 8. Setting $\mu/\Lambda-1=0$ regains Eq. (10); since the highest-order term in the integrand is $\propto x$, the exponent of (A.1) contains an $a^4$ term (explicit evaluation is "straightforward but suppressed", p. 156). The maximum $a_m$ of (A.1) follows from the quadratic $\{1+(\mu/\Lambda-1)a_m^2\}(1-a_m^2)=1/\gamma$ (p. 156).

Mean time between reversals (p. 156):

$$\langle T_r\rangle\approx\frac{4\pi}{\Delta_1}\left(\frac{1}{\phi_0''\psi_m''}\right)^{1/2}\exp(\psi_m)=\frac{4\pi}{\Delta_1}\left(\frac{1}{\phi_0''\psi_m''}\right)^{1/2}\frac{p(a_m)}{p(0)} \qquad\text{(Eq. A.2)}$$

where $\psi(a)$ and $\phi(a)$ are defined in Appendix A of HOS, $\phi_0''=\phi''(0)$, $\psi_m$ and $\psi_m''$ are $\psi$ and $\psi''$ at $a=a_m$, and the second form follows from $p(a)\propto\exp[\psi(a)]$ with $\psi(0)=0$. The first form was given in Appendix F of HOS. (A.2) remains valid when the correction term is included; explicit expressions for $\langle T_r\rangle$ and $dv$ derived from it were used for Table 1 but are omitted in the paper.

### 4.3 Non-dimensionalisation, units and time scales

- Dimensionless time $\tau=t/\tau_d$ with $\tau_d=R^2/\beta$ the turbulent diffusion time (p. 145). All rates in the statistical theory ($\Lambda$, $\Delta_0$, $\Delta_1$, $\lambda^i$) are therefore in units of $\tau_d^{-1}$; the amplitude $a$ is dimensionless (normalised to the non-linear equilibrium amplitude of the dipole mode), so $\gamma$ and $c^2$ are dimensionless ratios (p. 147). "Basically we measure $\Lambda/\tau_d$ and $\Lambda$ is model dependent" (p. 153).
- Eddy turnover time $\tau_c$ (fast); the simulations use $\tau_c=0.05\,\tau_d$ (pp. 146, 153).
- Conversion to years is *not assumed* but measured: $\tau_d=11\pm2$ kyr (reference run, p. 153) and $\simeq16\pm4$ kyr (second run, p. 153). HOS had concluded $\tau_d\sim3$ kyr, which this paper supersedes (p. 151). The abstract summarises "$\tau_d$ is about 10–15 kyr" (p. 143).
- Reversal duration is $(0.5$–$2)\tau_d$ in the model, i.e. 5–20 kyr for $\tau_d\approx10$ kyr (p. 151) and 8–30 kyr for the second run (p. 154).
- Field unit: the VADM data are in $10^{22}$ A m²; the model amplitude is mapped to VADM through the fitted peak position $A_m$ (Section 4.2 of the paper, p. 149): $a=(A/A_m)(1-1/\gamma)^{1/2}$.
- For the physical estimates of p. 155 the paper uses cgs: $R=3500$ km, $\eta\sim1.3\times10^4$ cm² s⁻¹, $\Omega'=d\Omega/dr\simeq-2.4\times10^{-18}$ (cm s)⁻¹.

### 4.4 Fluctuation model

- The helical convection at radius $R$ consists of cells of angular size $\pi/N_c$, homogeneously distributed over the sphere, $4N_c^2/\pi$ in total. Each cell contributes a mean part $\alpha_m(1-q)\cos\theta$ and a fluctuating part $\propto fF$ (p. 145).
- Azimuthal averaging reduces the random contribution of the $N=2N_c\sin\theta$ cells on a latitude circle by $N^{-1/2}$, giving Eq. (3) with its $(\sin\theta)^{-1/2}$ scaling (fluctuations increase towards the poles).
- Update rule: $F$ is updated every turnover time $\tau_c$, independently in each latitude interval $\Delta\theta=\pi/N_c$ (p. 145). $F$ has zero mean and unit variance; its distribution is not specified in this paper.
- $\beta$ is kept constant: fluctuations in $\beta$ are believed to be less important (p. 145).
- The statistical properties of the model depend on $kR$, $C_m$ and only one further combination, the forcing parameter $f^2\tau_c/N_c^2$ (p. 145). Reference run: $f^2\tau_c/N_c^2=1.46\times10^{-2}$ with $f=5.4$, $\tau_c=0.05$, $N_c=10$ (p. 146); second run: $2.45\times10^{-2}$ with $f=7$, $\tau_c=0.05$, $N_c=10$ (p. 153).
- The required fluctuation level is large, $\delta\alpha_{rms}/\langle\alpha\rangle\sim2.5$ (p. 148, quoting HOS/SOH).
- Robustness claims (p. 145): independence of $\delta\alpha$ between latitude intervals, the $(\sin\theta)^{-1/2}$ scaling and the precise form of the quenching are not essential; "any kind of 'rapid rattling' inside the core is likely to induce a similar behaviour of variability and reversals"; the convection model is much simpler than Narteau et al. (2000) because the statistics of $\langle\boldsymbol B\rangle$ depend only weakly on the turbulence model.

### 4.5 Nonlinearity / saturation

Global α-quenching: $\alpha_m(1-q)$ with $q\propto a_0^2$ (Eq. 2, p. 145). In the mode equations this appears as the factor $(1-a_0^2)$ multiplying $E_{ik}a_k$ (Eq. 5); in the Fokker–Planck description it yields the cubic drift $S\simeq\Lambda(1-a^2)a$ and the double-well potential (Eq. 8). The next-order correction (Appendix A) adds $(\mu-\Lambda)(1-a^2)a^3$ to $S$, derived from the requirement that the mean overtone toroidal field $\langle T_+\rangle$ vanishes at $a=1$ (numerically verified by storing $T_+(\theta)$ each time $a$ passes through 1 in the Fig. 3 run; this gives $B=-1$ in $g_1(a)/(1-a^2)=Aa(1+Ba^2+\cdots)$, p. 156). The constant of proportionality in $q\propto a_0^2$ is not given in this paper.

### 4.6 Parameter values used

| Symbol | Meaning | Value(s) | Where |
|---|---|---|---|
| $kR$ | radial wavenumber × radius (radial scale of field) | 0.5 (case 1); 0.3 (case 2); allowed $0<kR<1$ | p. 145; Fig. 3 p. 146; Fig. 8 p. 153; Table 2 p. 154 |
| $C_m$ | dynamo number (Eq. 4) | 100 (case 1); 120 (case 2); $C_m>0$ | same |
| $\Lambda$ | linear growth rate of dipole mode, units $\tau_d^{-1}$ | 0.527 (case 1); 0.815 (case 2) | Fig. 3, Fig. 8 captions |
| $f^2\tau_c/N_c^2$ | forcing parameter | $1.46\times10^{-2}$ (case 1); $2.45\times10^{-2}$ (case 2) | same; Table 2 |
| $f$ | fluctuation strength | 5.4 (case 1); 7 (case 2) | same |
| $\tau_c$ | eddy turnover time (units of $\tau_d$) | 0.05 | same |
| $N_c$ | number of cells across $\pi$ in colatitude | 10 | same |
| $\Delta\tau$ | timestep | 0.01 | same |
| run length | | $5\times10^4\tau_d$ (case 1); case 2 polarity shown to $10^4\tau_d$ | Fig. 3 caption; Fig. 8 |
| $\mu/\Lambda$ | ratio of bare to effective growth rate | 1.45 (case 1); 1.82 (case 2) | p. 156 |
| $\tau_d$ | turbulent diffusion time | $11\pm2$ kyr (case 1); $16\pm4$ kyr (case 2) | pp. 153–154, Table 2 |
| $\tau_{av}$ | running-average window for artificial Sint-800 series | ~5 kyr best; 0 too small, 10 too large | p. 152 |
| $R$ | representative core radius | 3500 km (used for $\beta=R^2/\tau_d$) | p. 155 |
| $\Delta R$ | depth of outer core | 2300 km | p. 155 |
| $\sigma$ (conductivity) | core conductivity | $6\times10^5$ S m⁻¹ (Merrill et al. 1996, p. 274) | p. 155 |
| $\eta$ | molecular magnetic diffusivity | $\sim1.3\times10^4$ cm² s⁻¹ | p. 155 |
| $\Omega'=d\Omega/dr$ | radial shear | $\simeq-2.4\times10^{-18}$ (cm s)⁻¹ ($-1°$ per year over 2300 km, Richards 2000) | p. 155 |
| $\alpha_m$ | mean α amplitude | $\sim-2.5\times10^{-4}$ cm s⁻¹ (from $\alpha_m\Omega'=C_m/\tau_d^2$) | p. 155 |
| $\beta/\eta$ | eddy magnetic Reynolds number | ~30 (case 1); ~20 (case 2) | p. 155, Table 2 |
| $C_\alpha/C_\Omega$ | $(\alpha_mR/\beta)/(\Omega'R^3/\beta)$ | $\sim10^{-3}$ | p. 155 |

## 5. Methods

**Numerical model.** One-dimensional (colatitude) αΩ mean-field dynamo at radius $R$, integrated with timestep $\Delta\tau=0.01$ with $F$ in Eq. (3) redrawn every $\tau_c=0.05$ (i.e. every 5 steps) independently in each of the $N_c=10$ latitude intervals (pp. 145–146). Long runs ($5\times10^4\tau_d$, "e.g. $10^5\tau_d$ can be performed in a short time", p. 146). Output: the toroidal field $T(\theta,\tau)$ and the normalised amplitude $a(\tau)$ of the fundamental (dipole) mode (Fig. 3). The reference run is run #2 of HOS. Grid resolution, boundary conditions and the projection used to extract $a(\tau)$ are not restated here (see HOS).

**Statistical theory.** Expansion in eigenmodes at marginal excitation → coupled mode equations (Eq. 5) → Fokker–Planck equation (Eq. 6) for $p(a)$ → stationary solution (Eqs. 9–10, A.1) → $\langle T_r\rangle$ from (A.2) and $dv$ from (12). The theory applies to dynamos of any dimensionality with a stationary fundamental mode (pp. 146, 155).

**Observational data: Sint-800** (Guyodo & Valet 1999; Fig. 1, p. 144). The record gives the VADM in units of $10^{22}$ A m² over the past 800 kyr with one data point per 1000 yr (799 points used, p. 148); the paper does not describe how the record was constructed. Time runs from right to left in Fig. 1; the known reversal at $t=780$ kyr is "allowed for" (the plotted VADM goes negative beyond ~780 kyr). The 1σ errors $e_i$ of each point $s_i$ are known and as a rule $e_i/s_i\sim0.1$ (p. 148). The data of Valet & Meynadier (1993), covering 4 Myr, were not used because they are less reliable (p. 155).

**Binning and error model** (p. 148, Fig. 5). Histograms depend strongly on the binning method; apparent over-abundances near 4, 5 and $6\times10^{22}$ A m² are of unclear origin (Guyodo, private communication) but not significant: since $e_i$ exceeds the bin width, the count $n_j$ in bin $j$ has Poisson uncertainty $n_j^{1/2}$, i.e. a constant relative error $n_j^{-1/2}=(799/50)^{-1/2}\cong25\%$ per bin for 50 equal-population bins—as large as the peaks.

**$\chi^2$ fitting of $p(a)$ to Sint-800** (Section 4.2, p. 149). With bin centroids $A_j$, widths $w_j$ and counts $n_j$: choose $\gamma$, $c^2$ and $A_m$ (VADM at which the distribution peaks); rescale $a_j=(A_j/A_m)(1-1/\gamma)^{1/2}$ (since the maximum of Eq. 10 is at $a_m=(1-1/\gamma)^{1/2}$); fix the normalisation $N$ by $N\sum_jp_jw_j=799$ where $p_j=p(a_j)$; evaluate $\chi^2=\sum_j(Np_jw_j-n_j)^2/n_j$; iterate over $(\gamma,c^2,A_m)$ to the minimum. Variable bin widths (equal population, as in Fig. 5 left) and ~100 bins are used to avoid empty bins. Goodness of fit is quoted as the number of "σ", defined as $|\chi^2-\nu|/(2\nu)^{1/2}$ with $\nu$ the effective number of degrees of freedom (Table 1, footnote a). No maximum-likelihood fitting was attempted (p. 155).

**Fitting $p(a)$ to simulation data** (Section 4.3, pp. 149–150). Histogram $a(\tau)$; because the $a$-axis is already normalised (peak at $(1-1/\gamma)^{1/2}$), only $\gamma$ and $c^2$ (plus normalisation) are free. The minimum $\chi^2$ is "many σ" because the fit uses $5\times10^4$ data points and small systematic theory–simulation differences become significant (footnote 4, p. 150). **Note:** $5\times10^4$ points over a $5\times10^4\tau_d$ run suggests one sample per $\tau_d$, but the sampling is not stated.

**Monte-Carlo error on $dv$ and on the autocorrelation** (pp. 151–152): artificial Sint-800 sets are generated by replacing each $s_i$ by $s_i+e_iN(0,1)$ with independent standard-normal draws; the spread gives $\pm0.003$ on $dv$ and the error bars in Fig. 7.

**Autocorrelation method for $\tau_d$** (Section 5.1, p. 152). Simulated $a_j$ (spacing $0.01\tau_d$) are turned into an artificial Sint-800 series: choose $\tau_d$ (fixes the absolute time axis), rectify the data, take a running average over $\tau_{av}$ years (mimicking the physical time resolution of Sint-800), pick 799 consecutive points spaced 1000 yr, subtract the mean. For the 799 Sint-800 points, reverse the direction of time and subtract the mean. For either series $\{x_i\}$ compute $\sum_{i=1}^{799-50}x_ix_{i+k}$ for $k=0$–50 and rescale to 1 at $k=0$. Sint-800 error bars are of order $799^{-1/2}\simeq0.035$ and are trustworthy up to lags of about 25 kyr. $\tau_d$ is adjusted until the model curve overlaps the data curve (Fig. 7). Different starting points in the long run give different artificial series (Fig. 7f), so the quality of overlap varies.

## 6. Results

### 6.1 The bistable-oscillator picture and its status (Section 3, pp. 146–148)

The dipole amplitude behaves as a strongly damped, randomly forced particle in the well $U\simeq\tfrac14\Lambda(1-a^2)^2$ (Fig. 4); the wells are the normal and reversed polarity states; velocity $da/d\tau$ is proportional to the force $-dU/da$. Overtones are excited because $\tau_c\ll\tau_d$; most fluctuations do not produce a reversal; reversals are rare because chance must produce the right sequence of α fluctuations—"a random process without a clearly identifiable cause, like radioactive decay", so that analysis of single events (Parker–Levy) is meaningless (p. 147). From HOS: reversals *must* occur as soon as overtones are excited (general proof, independent of overtone geometry, eigenvalues and dimensionality); the mean reversal rate depends on an overlap integral of the poloidal field of the fundamental mode and the rms toroidal overtone field; reversals are necessarily fast, lasting about one $\tau_d$ (p. 147).

Status before this paper (Section 3.1, pp. 147–148): (1) Sint-800 is well described by Eq. (10), itself a reasonable approximation to simulated $a(\tau)$; required $\delta\alpha_{rms}/\langle\alpha\rangle\sim2.5$; chron lengths are Poissonian. (2) The observed negative trend between $dv$ and $\langle T_r\rangle$ (McFadden & Merrill 1995) is explained: stronger forcing increases the oscillation amplitude in the well ($dv\uparrow$) and the jump frequency ($\langle T_r\rangle\downarrow$); $\langle T_r\rangle$ depends exponentially on $f^2\tau_c/N_c^2$, which also sets $dv$; requiring $\langle T_r\rangle\sim3\times10^5$ yr (typical of the past few Myr) gives $dv\simeq0.15$, about the observed value for that period. Properties such as the fraction of time below a threshold or the excursion frequency could be computed but have not been (p. 148).

### 6.2 Sint-800 amplitude distribution and errors (Section 4, p. 148)

Fig. 5 shows the 799-point histogram in 50 equal-population bins (left) and on a 30-point equidistant grid (right). Secondary peaks near 4, 5 and $6\times10^{22}$ A m² are within the ~25% per-bin Poisson error and are not significant. $dv$ computed from the data is $0.071\pm0.003$.

### 6.3 Fitting like to like (Section 4.1, p. 149)

VADM is the $\ell=1$, $m=0$ poloidal surface-field component. All eigenmodes have $m=0$; the fundamental (dipole) mode is equatorially antisymmetric with dominant $\ell=1$ and decreasing $\ell=3,5,\dots$; the first overtone is symmetric (quadrupolar, even $\ell$); the second overtone is antisymmetric (octupolar) with a sizable $\ell=1$ part but $a_2\ll a_0$. Hence VADM $\propto$ a linear combination of $a_0,a_2,a_4,\dots$ dominated by $a_0$; the error of omitting overtones is ~10% (SOH), comparable to the Sint-800 errors, so Sint-800 is fitted directly to $a_0\equiv a$. Two comparisons are made simultaneously: theory vs. simulation (is the theory trustworthy?) and theory vs. Sint-800 (does the model describe the data?). Direct fitting of simulation to data is impractical.

### 6.4 The fits (Sections 4.2–4.4; Table 1, p. 149; Fig. 6, p. 150)

| Case | Type of fit | $\gamma$ | $c^2$ | $A_m$ ($10^{22}$ A m²) | $\chi^2$ (σ)ᵃ | $\langle T_r\rangle$ ($\tau_d$)ᵇ | $dv$ᶜ | Figure |
|---|---|---|---|---|---|---|---|---|
| 1 | (10) to Sint-800 | 6.3 | 0.32 | 5.98 | 1 | 1090 | 0.076 | – |
| 2 | (10) to $a(\tau)$ of Fig. 3 | 8 | 0.50 | 6.07 | 1.3 | 928 | 0.070 | – |
| 3 | (10) to both | 7 | 0.45 | 6.04 | 1.1 | 665 | 0.079 | 6a |
| 4 | (A.1) to both | 4.8 | 0.22 | 6.14 | 1.4 | 1097 | 0.071 | 6b |
| 5 | As 4 but $a(\tau)$ of Fig. 8 | 5 | 0.50 | 6.10 | 5 | 138 | 0.091 | 9 |

ᵃ of the fit of (10) or (A.1) to Sint-800, quoted as $|\chi^2-\nu|/(2\nu)^{1/2}$; ᵇ from relation (A.2); ᶜ computed with (12) and (10) or (A.1). (Table 1, p. 149.)

- Case 1 is the three-parameter optimum to Sint-800 (1σ). Case 2 fixes $\gamma$, $c^2$ from the simulation histogram and optimises only $A_m$ and the normalisation against Sint-800: still 1.3σ, "quite acceptable" (p. 150). Hence a wide range of parameters is consistent with Sint-800, including those of the Fig. 3 run. The fits are sensitive to $\gamma$ or $c^2$ separately, but an increase in $\gamma$ can be largely offset by increasing $c^2$ (cases 1–4 show $c^2$ increasing with $\gamma$), p. 150.
- Going from the Sint-800 fit to a simulation is not straightforward: two numbers ($\gamma,c^2$) cannot fix three model parameters ($kR$, $C_m$, $f^2\tau_c/N_c^2$)—that could be cured with extra information such as $\langle T_r\rangle$—but "the real snag" is that no useful expressions exist for $\gamma$ and $c^2$ in terms of the model parameters (p. 149).
- Case 3 (Fig. 6a): compromise values between cases 1 and 2; $\chi^2\cong1.1\sigma$ for Sint-800.
- Case 4 (Fig. 6b): with the improved (A.1), compromise values; the representation of the simulation improved visibly while Sint-800 is not significantly worse ($\chi^2\cong1.4\sigma$). Conclusion: Sint-800, the Fig. 3 simulation and the theory "are in good agreement with one another" (p. 151); theoretical $dv$ values in Table 1 match the simulation's 0.070 and the data's $0.071\pm0.003$. **Note:** Evaluating (A.2) with the printed parameters reproduces the ⟨T_r⟩ column for cases 1–3 to 0.1% (1089, 927, 665 with $\Lambda=0.527$) and case 5 to 0.3% (137.6, with $\Lambda=0.815$ of the Fig. 8 run), but case 4 only to 1.1% (1109 against 1097); ⟨T_r⟩ changes by 6% across the rounding interval of the printed $\gamma=4.8$, so the printed value was presumably computed from unrounded fit parameters. The $dv$ column reproduces to ±0.001 for all five cases. See `dev/DISCREPANCIES.md` DR-2.

### 6.5 Time scales: mean reversal time and reversal duration (Section 4.5, p. 151)

- Measured $\langle T_r\rangle=570\pm61\,\tau_d$ in the Fig. 3 run vs. theoretical values 928–1097 $\tau_d$ from (A.2) (except case 3, 665): almost a factor 2; attributed to (A.2) being only approximate (HOS), "not a serious problem as long as we have no accurate theoretical expression".
- Current $\langle T_r\rangle$: one reversal in 800 kyr gives only the 95% confidence interval 220 kyr $<\langle T_r\rangle<$ 3.3 Myr. With the HOS value $\tau_d\sim3$ kyr, $(570\pm61)\tau_d$ gave 1.5–1.9 Myr, inside the interval. Extending to the past 10 Myr would give $(2$–$3)\times10^5$ yr, but for consistency the dipole variability and the reversal rate must refer to the same period (the past 800 kyr).
- Sint-800 actually permits measuring $\tau_d$ (Section 5): $\tau_d\sim11$ kyr rather than 3 kyr, so $\langle T_r\rangle\approx6$ Myr, outside the 95% interval.
- Reversal duration is $(0.5$–$2)\tau_d$, i.e. 5–20 kyr: "a little long but still acceptable".

### 6.6 Turbulent diffusion time from autocorrelation (Section 5, pp. 151–153; Fig. 7, p. 152)

- Why $\tau_d$ is poorly known: for the periodic solar dynamo there are two relations, the period $P\simeq(\alpha\,d\Omega/dr)^{-1/2}$ and the dynamo number (Eq. 4), so $\alpha$, $\beta$, $\tau_d$ can be estimated; the geodynamo is non-periodic, so only one combination (Eq. 4) is known (p. 151).
- Idea: after scaling both axes, $a(\tau)$ of Fig. 3 and the VADM of Fig. 1 should be "the same" process with a different realisation of the forcing; the autocorrelation function is the right tool, and choosing $\tau_d$ fixes the absolute time axis of the model.
- Prediction (Eq. 15): exponential decay with rate $2\Lambda$ (roll-down time $\Lambda^{-1}$), independent of $\tau_c$ for $\Lambda\tau_c\ll1$ (well satisfied; exact result in Gardiner 1990, p. 131). The typical time scale seen in Fig. 1 is thus the time to restore the non-linear equilibrium of the geodynamo (p. 152).
- Fig. 7 (all panels show Sint-800 dashed, same in each; model drawn): (a) $\tau_d=5.9$, $\tau_{av}=5.0$ kyr — model decays too fast; (b) $\tau_d=10.0$, $\tau_{av}=5.0$ — close overlap; (c) $\tau_d=14.3$, $\tau_{av}=5.0$ — model decays too slowly; (d) $\tau_d=10.0$, $\tau_{av}=0$ — too small an average, follows Eq. (15) well at small $s$; (e) $\tau_d=10.0$, $\tau_{av}=10.0$ — too much averaging; (f) $\tau_d=12.5$, $\tau_{av}=5.0$ — best fit for a *different* artificial series, showing that the fit is not always as good as in (b). The anticorrelation of Sint-800 at lags $t\gtrsim25$ kyr is generally reproduced by the model (p. 152).
- Result: $\tau_d=R^2/\beta=11\pm2$ kyr, with physical time resolution of Sint-800 of several kyr ($\tau_{av}\sim5$ kyr) (p. 153). Model dependent: what is measured is $\Lambda/\tau_d$.

### 6.7 Can all constraints be accommodated? (Section 6, pp. 153–154; Table 2, p. 154; Figs. 8–9)

Other models (different $kR$, $C_m$, $f^2\tau_c/N_c^2$) were explored to reduce $\langle T_r\rangle$ while keeping the rest; "we have not been completely successful. Many compromises are possible and we found no clear route to an optimal model." One of the better ones is case 2 (Fig. 8/9, Table 1 row 5):

| Case | Figure | $kR$ | $C_m$ | $f^2\tau_c/N_c^2$ | $\tau_d$ (kyr) | $\tau_r$ (kyr)ᵃ | $\langle T_r\rangle$ (Myr)ᵇ | $dv$ᶜ | $\beta/\eta$ |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 3, 6b | 0.5 | 100 | $1.46\times10^{-2}$ | $11\pm2$ | 5–20 | 5–7.5 | 0.070 | ~30 |
| 2 | 8, 9 | 0.3 | 120 | $2.45\times10^{-2}$ | $16\pm4$ | 8–30 | 1.3–2.2 | 0.081 | ~20 |

ᵃ duration of a reversal; ᵇ from run and the value of $\tau_d$; ᶜ from run; the Sint-800 value is $0.071\pm0.003$. (Table 2, p. 154.)

For case 2: $\tau_d\simeq16\pm4$ kyr from the autocorrelation method, larger mainly because $\Lambda$ is larger (0.815); measured $\langle T_r\rangle=(113\pm5)\tau_d$ now agrees rather well with the theoretical 138 $\tau_d$; corresponds to ~1.8 Myr, inside (0.22, 3.3) Myr. The smaller $\langle T_r\rangle$ arises because $p(0)$ is larger (the amplitude is more often near $a=0$) and $D\simeq\Delta_1$ is larger there (cf. Eq. A.2). The price: the (A.1) fit to Sint-800 is 5σ (Table 1 case 5, Fig. 9), though "on the whole the fits are not unreasonable"; $dv=0.081$ is somewhat larger than $0.071\pm0.003$; the reversal duration is 8–30 kyr, "rather long" (p. 154).

### 6.8 Discussion and summary (Section 7, pp. 154–155)

- Philosophy: full 3-D MHD is "strictly speaking" the only complete approach, but risks losing the perspective on basic mechanisms; simulations and simple models should be developed jointly (the Brownian-motion/Langevin analogy).
- Either $\langle T_r\rangle$ (case 1) or the reversal duration (case 2) is too long; both trace back to $\tau_d$: "If that had been a factor 2–3 times smaller the model would describe the data perfectly." This may indicate model inadequacy, or overlooked systematics in constructing the artificial Sint-800 series (only the finite time resolution was allowed for, via the running average).
- Variability and reversals are "inseparably linked—you cannot have one without the other"; however small the fluctuations, there is always a next reversal, with waiting time growing exponentially with $(f^2\tau_c/N_c^2)^{-1}$—possibly what happened in the Cretaceous superchron.
- Transfer to 2-D models with an independent radial coordinate is affirmed provided the fundamental dipole mode is stationary, because the mode equations (5) do not depend on dimensionality. Mechanisms favouring stationary αΩ action: spatial separation of induction effects (Deinzer et al. 1974; Gibbons 1998); α² dynamos are usually non-oscillatory (Roberts 1972); meridional circulation (Roberts 1972; Sarson & Jones 1999); a conducting solid inner core (Hollerbach & Jones 1993; Glatzmaier & Roberts 1995a,b). A model with a passive conducting inner core has been constructed with promising initial tests; it could test whether reversal duration is set by the inner-core Ohmic time (Gubbins 1999) or by $\tau_d$ of the outer core (p. 155).
- Fitting data determines only some characteristics of the bistable potential (slopes, height of the central hill); trade-offs between parameters limit the information gained on the "substrate" (fluctuation-free) model.

### 6.9 Mean-field self-consistency (Section 7, p. 155)

The eddy magnetic Reynolds number is $R_m=u\lambda_c/\eta\sim\beta/\eta$ (since $u\tau_c\sim\lambda_c$). With $\beta=R^2/\tau_d$, $R=3500$ km, core conductivity $6\times10^5$ S m⁻¹ (Merrill et al. 1996, p. 274) hence $\eta\sim1.3\times10^4$ cm² s⁻¹, the resulting $\beta/\eta$ (~30 and ~20 for the two cases, Table 2) "are perhaps just large enough to support the notion of field lines frozen into the eddy convection." From Eq. (4), $\alpha_m\Omega'=C_m/\tau_d^2$; with $\Omega'=d\Omega/dr\simeq-2.4\times10^{-18}$ (cm s)⁻¹ (a differential rotation of $-1°$ per year over $\Delta R=2300$ km, Richards 2000) this gives $\alpha_m\sim-2.5\times10^{-4}$ cm s⁻¹. Then $C_\alpha\equiv\alpha_mR/\beta$ and $C_\Omega\equiv\Omega'R^3/\beta$ give $C_\alpha/C_\Omega\sim10^{-3}$, justifying the αΩ approximation. **Note:** $C_\alpha/C_\Omega=\alpha_m/(\Omega'R^2)$ is independent of $\beta$; with the quoted numbers it evaluates to $\approx9\times10^{-4}$.

## 7. Figures and tables

| Item | Page | What it shows (axes, quantities) | Run / data | Key readings |
|---|---|---|---|---|
| Fig. 1 | 144 | Sint-800 VADM ($10^{22}$ A m², axis −5 to 10) vs. age (kyr, 0–800); time runs right to left | Guyodo & Valet 1999 | 1 point/1000 yr; values mostly 4–8; sign flipped for the reversal at 780 kyr (dips to about −4); errors ~10%; $dv=0.071\pm0.003$ |
| Fig. 2 | 144 | Schematic of the αΩ model: stationary dipolar fundamental mode (heavy lines), toroidal loops, helical cells; dashed loop with wrong orientation = random irregularity in α | — | Conceptual only |
| Fig. 3 | 146 | Top: toroidal field at radius $R$ vs. colatitude (0, π/2, π) and time (0–3000 $\tau_d$), colour; bottom: dipole amplitude $a(\tau)$ (−2 to 2) | Case 1 (run #2 of HOS): $kR=0.5$, $C_m=100$, $\Lambda=0.527$, $f^2\tau_c/N_c^2=1.46\times10^{-2}$ ($f=5.4$, $\tau_c=0.05$, $N_c=10$), $\Delta\tau=0.01$ | Run $5\times10^4\tau_d$, 88 reversals, $\langle T_r\rangle=570\pm61$, $dv=0.070$; $a$ oscillates around ±1; about 7 sign changes (8 polarity intervals) in the 3000 $\tau_d$ shown, the longest interval ≈1100 $\tau_d$ |
| Fig. 4 | 147 | Bistable potential $U(a)=\tfrac14\Lambda(1-a^2)^2$ (solid) and $p(a)$ (dashed) vs. $a$ (−2 to 2); particle in the right well with arrows to the left well | — | Wells at $a=\pm1$; $p(a)$ peaked near ±1 |
| Fig. 5 | 148 | Sint-800 VADM histograms: left, 50 equal-population bins (counts per $10^{22}$ A m², 0–400); right, 30-point equidistant grid (0–250) | 799 Sint-800 points | Left: spiky, peak ~350 near 6; right: smooth, peak ~225 near 5.5–6; $dv=0.071\pm0.003$ |
| Fig. 6 | 150 | (a) left column: Eq. (10) fitted to Sint-800 with 1σ error bars (top; counts per $10^{22}$ A m², 0–250, VADM 0–10) and to the amplitude histogram of Fig. 3 (bottom; counts per unit $a$, 0–$10^5$, $a$ 0–1.5). (b) right column: same with Eq. (A.1). Horizontal axes scaled so the two $p(a)$ curves coincide | Cases 3 (a) and 4 (b) | VADM peak ~200–225 at ~6; amplitude peak ~$9\times10^4$ at $a\approx0.9$; (b) visibly better for the simulation |
| Fig. 7 | 152 | Autocorrelation (−0.5 to 1) vs. lag (0–50 kyr): Sint-800 dashed with error bars (same in all panels), model drawn; panels (a) $\tau_d=5.9$, $\tau_{av}=5.0$; (b) 10.0, 5.0; (c) 14.3, 5.0; (d) 10.0, 0.0; (e) 10.0, 10.0; (f) 12.5, 5.0 (kyr) | Case 1 run | Sint-800 crosses zero near 27 kyr and reaches ≈−0.2 near 40 kyr; (b) overlaps best; (f) is a different artificial series |
| Fig. 8 | 153 | Top: polarity bar code (0–10 000 $\tau_d$); bottom: $a(\tau)$ (0–3000 $\tau_d$) | Case 2: $kR=0.3$, $C_m=120$, $\Lambda=0.815$, $f^2\tau_c/N_c^2=2.45\times10^{-2}$ ($f=7$, $\tau_c=0.05$, $N_c=10$), $\Delta\tau=0.01$ | $\langle T_r\rangle=113\pm5$, $dv=0.081$; many more reversals than Fig. 3 (of order 20–25 sign changes in the 3000 $\tau_d$ shown [uncertain]; the average rate is 3000/113 ≈ 27) |
| Fig. 9 | 153 | Same layout as Fig. 6b (VADM top, amplitude bottom) using Eq. (A.1) and the amplitudes of Fig. 8 | Case 5 of Table 1 | VADM peak ~210 at ~6; amplitude peak ~$8.5\times10^4$ at $a\approx0.95$; fit to Sint-800 is 5σ |
| Table 1 | 149 | Fitting parameters: $\gamma$, $c^2$, $A_m$, $\chi^2$(σ), theoretical $\langle T_r\rangle$, $dv$ for five fit cases | — | Reproduced in §6.4 |
| Table 2 | 154 | Model parameters and derived physical quantities for the two cases | — | Reproduced in §6.7 |

## 8. Comparison with observations

| Observable | Data (source, page) | Model, case 1 (Fig. 3) | Model, case 2 (Fig. 8) | Stated agreement |
|---|---|---|---|---|
| VADM distribution shape | Sint-800 histogram (Fig. 5) | Eq. (10): 1σ (case 1), 1.3σ (case 2), 1.1σ (case 3); Eq. (A.1): 1.4σ (case 4) | Eq. (A.1): 5σ (case 5) | Good for case 1; "not unreasonable" for case 2 (pp. 150–151, 154) |
| Secular dipole variation $dv$ | $0.071\pm0.003$ (p. 148) | 0.070 (run), 0.071–0.079 (theory) | 0.081 (run), 0.091 (theory) | Matches for case 1; somewhat large for case 2 (pp. 151, 154) |
| Autocorrelation function | Sint-800, lags to 50 kyr, errors ~0.035 (Fig. 7) | Overlaps for $\tau_d\approx10$–11 kyr, $\tau_{av}\sim5$ kyr; anticorrelation beyond 25 kyr reproduced | overlaps for $\tau_d\simeq16$ kyr | Good, with realisation dependence (pp. 152–153) |
| Mean time between reversals | 95% CI 0.22–3.3 Myr from one reversal in 800 kyr (p. 151); (2–3)×10⁵ yr over past 10 Myr | 5–7.5 Myr (570 $\tau_d$ × 11 kyr) — outside | 1.3–2.2 Myr (113 $\tau_d$ × 16 kyr) — inside | Case 1 too long; case 2 acceptable (pp. 151, 153) |
| Reversal duration | $\lesssim10$ kyr (p. 143) | 5–20 kyr | 8–30 kyr | Factor ~2 too long (abstract; p. 154) |
| Chron-length distribution | Poissonian | Poissonian (HOS) | — | Agreement (p. 148) |
| $dv$–$\langle T_r\rangle$ anticorrelation | McFadden & Merrill 1995 | explained qualitatively; $\langle T_r\rangle\sim3\times10^5$ yr ↔ $dv\simeq0.15$ | — | Agreement (p. 148) |
| Eddy magnetic Reynolds number | — | $\beta/\eta\sim30$ | $\sim20$ | "perhaps just large enough" for mean-field theory (p. 155) |

## 9. Assumptions, caveats and limitations

Stated by the authors:
- One-dimensional model in colatitude at a representative radius $R$ with radial structure $\exp(ikr)/kr$; applicability to 2-D models requires a stationary fundamental dipole mode, deferred to Section 7 (pp. 145, 155).
- Constant $\beta$ (no fluctuations in $\beta$); $\delta\alpha$ independent per latitude interval and $\propto(\sin\theta)^{-1/2}$ for analytical convenience, claimed non-essential (p. 145).
- The functional dependence of $\Delta_0$, $\Delta_1$ (hence $\gamma$, $c^2$) on $kR$, $C_m$, $f^2\tau_c/N_c^2$ is unknown (pp. 147, 149).
- Eq. (10)/(A.1) are approximations; the simulation histogram deviates "formally significantly" from (10) when $5\times10^4$ points are used (footnote 4, p. 150). Relation (A.2) for $\langle T_r\rangle$ is only approximate (factor ~2, p. 151).
- VADM is identified with $a_0$ alone, neglecting overtone contributions (~10% error, p. 149).
- Sint-800 has ~10% errors per point, possible systematic effects (peaks near 4, 5, 6×10²² A m², origin unclear), finite physical time resolution of several kyr, and only one reversal, so $\langle T_r\rangle$ is weakly constrained (pp. 148, 151, 153–155). Possible overlooked systematics in constructing the artificial series (p. 154).
- $\tau_d$ from the autocorrelation method is model dependent ($\Lambda/\tau_d$ is what is measured, p. 153), and the quality of the overlap depends on the realisation (Fig. 7f).
- Binning-based $\chi^2$ fits rather than maximum likelihood (p. 155); the "σ" measure is $|\chi^2-\nu|/\sqrt{2\nu}$.
- Mean-field theory at $\beta/\eta\sim20$–30 is only marginally justified (p. 155).
- Whether "any rapid rattling" produces reversals is likely but "has not been proven" (p. 147).

**Note:** additional points for implementers. (i) Eqs. (13)–(14) use $D$ and $\Delta_0$ rather than $\sqrt D$ as noise amplitude; use the Fokker–Planck form (6)–(7) as the authoritative definition. (ii) $dv$ is a variance ratio; the rms relative variation is $\sqrt{dv}\approx0.27$ for Sint-800, while the "~15%" quoted in the GAFD abstract corresponds to $dv\simeq0.15$ (p. 148). (iii) The paper does not say whether the Sint-800 autocorrelation uses the signed series of Fig. 1 or the unsigned VADM; the model series is rectified, so like-to-like suggests the unsigned record. (iv) The theoretical $\langle T_r\rangle$ column of Table 1 and the measured 570 $\tau_d$ differ by ~2; the physical $\langle T_r\rangle$ in Table 2 uses the *measured* run value.

## 10. Relation to the other two papers in this set

- **HOS = Hoyng, Ossendrijver & Schmitt (2001), GAFD 94, 263–314** (received June 2000) is the long paper with the full derivations: the 1-D model equations, the eigenmode expansion and the mode equations (its factor $s$ is absorbed here into $E_{ik}$), the Fokker–Planck derivation, the proof that overtone excitation implies reversals, the expressions for $\langle T_r\rangle$ (its Appendix F; first form of Eq. A.2 here), the functions $\phi$, $\psi$ (its Appendix A), the correction function $g_1(a)$ (its Appendix B.1, Eqs. 33 and B3 referred to on p. 156), the elaborate fitting method of its Section 5.1, run #2 (reproduced here as Fig. 3), the conclusion $\tau_d\sim3$ kyr (superseded here by 11 ± 2 kyr), the $dv$–reversal-rate relation and the result that $\langle T_r\rangle\sim(2$–$3)\times10^5$ yr implies $dv\simeq0.15$. This paper cites HOS for "details on all matters" of Section 3 and takes from it the model, the theory and the reference run; it gives back a quantitative test of (10)/(A.1) against Sint-800, the numerical finding $\langle T_+\rangle|_{a=1}=0$ that fixes $B=-1$ in $g_1$, and the measured $\tau_d$.
- **SOH = Schmitt, Ossendrijver & Hoyng (2001), PEPI 125, 119–124** (received Dec 2000) is the compact 6-page summary of the method (Fokker–Planck equation for the dipole moment, bistable potential from supercritical excitation and nonlinear limitation) with a first Sint-800 comparison of the amplitude distribution and the secular-variation/reversal-rate relation. This paper cites SOH for the ~10% estimate of the error made by identifying VADM with $a_0$ (p. 149) and for the earlier finding that the Fig. 3 run corresponds closely to Sint-800.
- **This paper (received July 2001)** is the detailed data-analysis paper: the Sint-800 binning and error model, the $\chi^2$ fitting procedure and its degeneracy, the improved $p(a)$ (A.1), the autocorrelation technique and the resulting $\tau_d\approx10$–15 kyr, the attempt to satisfy all constraints with a second run, and the mean-field self-consistency estimate. It supersedes HOS's $\tau_d\sim3$ kyr and consequently revises the physical $\langle T_r\rangle$ of the reference run from 1.5–1.9 Myr to ~6 Myr.

## 11. Implementation guide for coding agents

**What this paper specifies (and what it does not).** The paper gives the model structure (Eqs. 1–4), the mode-equation form (5), the Fokker–Planck/potential description (6–11), the fitting and autocorrelation procedures, and two complete parameter sets. It does *not* restate: the explicit 1-D equations in $\theta$ for the poloidal and toroidal field components after inserting $\exp(ikr)/kr$; boundary conditions at the poles; the grid in $\theta$ or the time integrator; the constant in $q\propto a_0^2$; the proportionality constant in Eq. (3); the probability distribution of $F$; the definitions and normalisation of the eigenmodes, $E_{ik}$, $V_{ik}$; how $a_0(\tau)$ is projected out of the field; the reversal-counting criterion; the functions $\phi$, $\psi$ in (A.2). All of these must be taken from HOS (GAFD 94, 263) or decided and documented by the implementer.

**System to solve.** Eq. (1) with $\boldsymbol v=\Omega r\sin\theta\,\boldsymbol e_\varphi$, constant $d\Omega/dr$, radial dependence $\exp(ikr)/kr$ at $r=R$, constant $\beta\gg\eta$, and $\alpha$ from Eqs. (2)–(3), in dimensionless time $\tau=t/\tau_d$ ($\tau_d=R^2/\beta$). Free parameters: $kR\in(0,1)$, $C_m>0$ (Eq. 4), and the fluctuation triple $(f,\tau_c,N_c)$ entering the statistics only via $f^2\tau_c/N_c^2$.

**Standard runs.** Case 1 (Fig. 3; run #2 of HOS): $kR=0.5$, $C_m=100$, $f=5.4$, $\tau_c=0.05$, $N_c=10$ ($f^2\tau_c/N_c^2=1.46\times10^{-2}$), $\Delta\tau=0.01$, length $5\times10^4\tau_d$; expected $\Lambda=0.527$, $\mu/\Lambda=1.45$. Case 2 (Fig. 8): $kR=0.3$, $C_m=120$, $f=7$, $\tau_c=0.05$, $N_c=10$ ($2.45\times10^{-2}$), $\Delta\tau=0.01$; expected $\Lambda=0.815$, $\mu/\Lambda=1.82$.

**Constructing the stochastic α.** Divide colatitude into $N_c=10$ intervals of width $\pi/N_c$. Every $\tau_c=0.05$ (every 5 steps of $\Delta\tau=0.01$) draw a new zero-mean, unit-variance $F$ independently in each interval; set $\delta\alpha\propto fF/\sqrt{2N_c\sin\theta}$ (Eq. 3). The constant of proportionality is not given here; the physically motivated size is $\delta\alpha_{rms}/\langle\alpha\rangle\sim2.5$ (p. 148). Keep $\beta$ constant.

**Eigenvalue problem.** Linearise at marginal excitation (the $C_m$ at which the fundamental dipole mode has $\lambda^0=0$); order modes by $\mathrm{Re}\,\lambda^i$ (all $<0$ for $i\ge1$). The fundamental mode is equatorially antisymmetric (dipole, dominant $\ell=1$), the first overtone symmetric (quadrupole), the second antisymmetric (octupole) (p. 149). Normalise amplitudes to the non-linear equilibrium amplitude of the dipole mode so that $a=\pm1$ are the well bottoms (p. 146).

**Diagnostics as defined in the paper.**
- Dipole amplitude $a(\tau)=a_0(\tau)$; record at least once per $\tau_d$ (footnote 4 implies $5\times10^4$ histogram points for the $5\times10^4\tau_d$ run).
- Reversals: count sign changes of the polarity (Fig. 8 top shows a polarity bar code); $\langle T_r\rangle=$ run length / number of reversals ($5\times10^4/88\approx570\,\tau_d$). The paper does not state how aborted reversals/excursions are excluded; HOS should be consulted. Reversal duration: expect $(0.5$–$2)\tau_d$.
- $dv$ from Eq. (12) using $|a|$ (p. 148). Expected: 0.070 (case 1), 0.081 (case 2).
- Amplitude histogram and fit of Eq. (10) with the peak fixed at $(1-1/\gamma)^{1/2}$: expected $\gamma=8$, $c^2=0.50$ for case 1 (Table 1 case 2).
- Autocorrelation: rectify, running-average over $\tau_{av}$ (try 5 kyr for a chosen $\tau_d$), resample at 1 kyr, take 799 points, subtract mean, compute $\sum_{i=1}^{749}x_ix_{i+k}$, $k=0..50$, normalise to 1 at $k=0$. Expected e-folding rate $2\Lambda$ per $\tau_d$ at small lags (Eq. 15), i.e. an e-folding lag of $\approx0.95\tau_d\approx10$ kyr for case 1 with $\tau_d=11$ kyr (**Note:** derived, not stated). Compare with the Sint-800 curve (zero crossing ≈27 kyr, minimum ≈−0.2 near 40 kyr, Fig. 7).

**Fitting Eq. (10)/(A.1) to VADM data.** Follow Section 4.2 exactly: equal-population bins (~100), parameters $(\gamma,c^2,A_m)$, $a_j=(A_j/A_m)(1-1/\gamma)^{1/2}$ for (10) (for (A.1) use the quadratic for $a_m$), normalisation $N\sum_jp_jw_j=799$, $\chi^2=\sum_j(Np_jw_j-n_j)^2/n_j$, report $|\chi^2-\nu|/\sqrt{2\nu}$. Expected optimum for Sint-800: $\gamma=6.3$, $c^2=0.32$, $A_m=5.98\times10^{22}$ A m² (1σ). **Note:** Eq. (10) follows from (9) with $S=\Lambda(1-a^2)a$, $D=\Delta_0(a^2+c^2)$ via the substitution $x=a'^2$: $\int_0^a2S/D\,da'=\gamma\int_0^{a^2}(1-x)/(x+c^2)\,dx=-\gamma a^2+\gamma(1+c^2)\ln(1+a^2/c^2)$; the same substitution with $S=\Lambda(1-a^2)a[1+(\mu/\Lambda-1)a^2]$ gives (A.1), and $d\ln p/da=0$ reproduces $a_m^2=1-1/\gamma$ and the quadratic of p. 156. These checks confirm the transcription.

**Validation targets.** Case 1: $\langle T_r\rangle=570\pm61\,\tau_d$, $dv=0.070$, $(\gamma,c^2)\approx(8,0.50)$, $\mu/\Lambda=1.45$, autocorrelation overlap with Sint-800 at $\tau_d\approx10$–11 kyr with $\tau_{av}\approx5$ kyr. Case 2: $\langle T_r\rangle=113\pm5\,\tau_d$, $dv=0.081$, $(\gamma,c^2)\approx(5,0.50)$ with (A.1), $\mu/\Lambda=1.82$, $\tau_d\approx16$ kyr. Data: $dv=0.071\pm0.003$; Table 1 case 1 fit.

**Pitfalls.**
- Time unit: $\tau_d$ is a fitted quantity (11 ± 2 or 16 ± 4 kyr), not an input; HOS's 3 kyr is obsolete. Physical $\langle T_r\rangle$ scales linearly with $\tau_d$.
- The theoretical $\langle T_r\rangle$ from (A.2) is ~2× the measured value (except case 3); do not use (A.2) as a validation target for the simulation.
- $dv$ is a variance ratio ($\sqrt{dv}$ is the relative rms).
- Histogram shape depends strongly on binning (Fig. 5); use equal-population bins and the Poisson per-bin error $n_j^{1/2}$.
- $\chi^2$ against a $5\times10^4$-point simulation histogram will be "many σ" due to small systematic errors; this is expected (footnote 4).
- Degeneracy: $\gamma\uparrow$ compensated by $c^2\uparrow$; report both and the $\chi^2$ surface.
- Autocorrelation: $\tau_{av}=0$ is too little smoothing and 10 kyr too much; error bars are reliable only up to ~25 kyr lag; results vary with the starting point in the run (Fig. 7f).
- The Sint-800 time axis in Fig. 1 runs right to left; the paper reverses the direction of time before computing the autocorrelation.

## 12. Symbol glossary

| Symbol | Meaning | Units | First defined |
|---|---|---|---|
| $\langle\boldsymbol B\rangle$ | azimuthally averaged (axisymmetric) magnetic field | field | Eq. 1, p. 145 |
| $\boldsymbol v=\Omega r\sin\theta\,\boldsymbol e_\varphi$ | mean (differential-rotation) flow | velocity | p. 145 |
| $\alpha$, $\alpha_m$ | dynamo coefficient; amplitude of its steady part | velocity (cm s⁻¹) | Eq. 2; p. 155 |
| $q\propto a_0^2$ | α-quenching factor | dimensionless | p. 145 |
| $\delta\alpha$ | fluctuating part of α | velocity | Eqs. 2–3 |
| $\beta$, $\eta$ | turbulent, molecular magnetic diffusivity | cm² s⁻¹ | p. 145 |
| $k$, $kR$ | radial wavenumber; dimensionless radial scale | cm⁻¹; dimensionless | p. 145 |
| $R$ | representative radial position (3500 km) | cm | pp. 145, 155 |
| $\theta$ | colatitude | rad | p. 145 |
| $t$, $\tau=t/\tau_d$ | physical, dimensionless time | yr; dimensionless | p. 145 |
| $\tau_d=R^2/\beta$ | turbulent diffusion time | kyr | p. 145 |
| $\tau_c$ | eddy turnover time (update interval of $F$) | units of $\tau_d$ | pp. 144–146 |
| $f$, $F(\theta,\tau)$ | fluctuation strength; random function, zero mean, unit variance | dimensionless | p. 145 |
| $N_c$; $N=2N_c\sin\theta$ | cells across $\pi$; cells on a latitude circle | integer | p. 145 |
| $C_m$ | dynamo number | dimensionless | Eq. 4 |
| $f^2\tau_c/N_c^2$ | forcing parameter | dimensionless | p. 145 |
| $\lambda^i$ | eigenvalues at marginal excitation ($\lambda^0=0$) | $\tau_d^{-1}$ | p. 146 |
| $a_i$; $a\equiv a_0$ | mode amplitudes; dipole amplitude normalised to equilibrium | dimensionless | Eq. 5–6 |
| $E_{ik}$, $V_{ik}(\tau)$ | overlap integrals (quenching; fluctuations) | $\tau_d^{-1}$ | Eq. 5 |
| $p(a,\tau)$, $p(a)$ | probability distribution of $a$; its stationary form | — | Eqs. 6, 9–10 |
| $S$, $D$ | drift, diffusion coefficients | $\tau_d^{-1}$ | Eqs. 6–7 |
| $\Lambda$ | linear growth rate of dipole mode at $a\ll1$ | $\tau_d^{-1}$ | Eq. 7 |
| $\Delta_0$, $\Delta_1$ | coefficients of $D=\Delta_0a^2+\Delta_1$ | $\tau_d^{-1}$ | Eq. 7 |
| $U(a)$ | bistable potential, $-dU/da=S$ | $\tau_d^{-1}$ | Eq. 8 |
| $\gamma=\Lambda/\Delta_0$, $c^2=\Delta_1/\Delta_0$ | shape parameters of $p(a)$ | dimensionless | Eq. 11 |
| $a_m=(1-1/\gamma)^{1/2}$ | location of maximum of (10) | dimensionless | p. 149 |
| $A$, $A_m$, $A_j$, $w_j$, $n_j$ | VADM; VADM at peak; bin centroid, width, count | $10^{22}$ A m² | p. 149 |
| $N$ (Section 4.2) | normalisation of $p$ | — | p. 149 |
| $\nu$ | effective degrees of freedom | — | Table 1 |
| $dv$ | secular dipole variation | dimensionless | Eq. 12 |
| $T_r$, $\langle T_r\rangle$ | interval between reversals; its mean | $\tau_d$ or yr | p. 143 |
| $\tau_r$ | duration of a reversal | kyr | Table 2 |
| $s_i$, $e_i$ | Sint-800 value and 1σ error of point $i$ | $10^{22}$ A m² | p. 148 |
| $L(\tau)$ | Langevin force | — | Eq. 13 |
| $\delta a$, $\tilde a$ | deviation from $a=1$; rectified, mean-subtracted series | dimensionless | Eqs. 14–15 |
| $s$ | time lag | $\tau_d$ (or kyr) | Eq. 15 |
| $\tau_{av}$ | running-average window | kyr | p. 152 |
| $x_i$ | mean-subtracted series for autocorrelation | — | p. 152 |
| $\mu$ | bare supercritical growth parameter (HOS) | $\tau_d^{-1}$ | p. 156 |
| $g_1(a)$, $A$, $B$ | correction function and its expansion coefficients (HOS) | — | p. 156 |
| $\tilde P_0$, $T_+(\theta)$ | poloidal fundamental-mode function; toroidal overtone field (HOS notation) | — | p. 156 |
| $\phi(a)$, $\psi(a)$ | functions in the $\langle T_r\rangle$ formula (HOS App. A); $p\propto e^{\psi}$ | — | Eq. A.2 |
| $R_m=u\lambda_c/\eta$ | eddy magnetic Reynolds number ($\sim\beta/\eta$) | dimensionless | p. 155 |
| $\Omega'=d\Omega/dr$ | radial shear | (cm s)⁻¹ | p. 155 |
| $C_\alpha=\alpha_mR/\beta$, $C_\Omega=\Omega'R^3/\beta$ | α and Ω dynamo numbers | dimensionless | p. 155 |
| $E$ | Ekman number | dimensionless | p. 144 |

## 13. Section-by-section outline

| Section | Pages | Content |
|---|---|---|
| Abstract, keywords | 143 | Summary; $\tau_d\approx10$–15 kyr; reversal duration factor ~2 too long; mean-field theory self-consistent. |
| 1. Introduction | 143–145 | Paleomagnetic facts; 3-D simulations and their parameter gap; need for quantitative simple models; prior stochastic reversal ideas; the two key elements of the model; footnote defining HOS and SOH. Figs. 1–2. |
| 2. Introducing the model | 145–146 | Eqs. (1)–(4); geometry; α-quenching; convection-cell fluctuation model; three control parameters; requirement of a stationary dipole mode. |
| 3. The bistable oscillator picture | 146–147 | Fig. 3 reference run; mode equations (5); Fokker–Planck (6)–(7); potential (8), Fig. 4; Parker–Levy vs. continuous fluctuations; HOS results on reversals; stationary solution (9)–(11). |
| 3.1 Status | 147–148 | $dv$ (12); main results of HOS/SOH; $dv$–$\langle T_r\rangle$ relation. |
| 4. The distribution of the dipole amplitude | 148 | Binning of Sint-800 (Fig. 5); Poisson error per bin ~25%. |
| 4.1 Fitting like to like | 149 | VADM vs. $a_0$; overtone content; two simultaneous comparisons. |
| 4.2 Optimal fit to the Sint-800 data | 149 | $\chi^2$ procedure; Table 1 case 1. |
| 4.3 Best fit to the simulation data | 149–150 | Why one cannot simply simulate the fitted parameters; fit to $a(\tau)$ (case 2); degeneracy of $\gamma$, $c^2$; footnote 4. |
| 4.4 Compromises and improvements | 150–151 | Cases 3 and 4; Fig. 6; agreement of $dv$; Monte-Carlo error. |
| 4.5 Time scales | 151 | Factor-2 discrepancy of theoretical $\langle T_r\rangle$; 95% confidence interval; $\tau_d$ 3 vs. 11 kyr; reversal duration 5–20 kyr. |
| 5. The turbulent diffusion time of the geodynamo | 151–152 | Why $\tau_d$ is unknown; idea of the autocorrelation method; Eqs. (13)–(15). |
| 5.1 Autocorrelation functions | 152–153 | Procedure; Fig. 7; $\tau_d=11\pm2$ kyr; model dependence. |
| 6. Can we accommodate all constraints? | 153–154 | Second run (Figs. 8–9, Table 1 row 5); $\tau_d\simeq16\pm4$ kyr; $\langle T_r\rangle\approx1.8$ Myr; costs (5σ fit, $dv=0.081$, 8–30 kyr reversals). |
| 7. Discussion and summary | 154–155 | Philosophy; summary of model; Table 2; sources of discrepancy; Valet–Meynadier data; maximum likelihood; carry-over to 2-D; inner core; mean-field self-consistency. |
| Acknowledgements | 155 | Guyodo, Gill, Langereis, Haverkamp. |
| Appendix A. Improved expression for $p(a)$ | 156 | $g_1(a)$, $B=-1$, Eq. (A.1), $\mu/\Lambda$ values, quadratic for $a_m$, Eq. (A.2) and its sensitivity to $p(0)$ and $\Delta_1$. |
| References | 156–157 | 33 entries (McFadden & Merrill 1995, cited on p. 148, is missing from the list). |

## 14. Key references cited

- Hoyng, Ossendrijver & Schmitt (2001), GAFD 94, 263–314 (HOS): source of the model, the statistical theory, run #2, and all derivations referred to (Eqs. 33, B3, Appendices A, B.1, F).
- Schmitt, Ossendrijver & Hoyng (2001), PEPI 125, 119–124 (SOH): compact summary; the ~10% estimate of overtone contamination of the VADM.
- Guyodo & Valet (1999), Nature 399, 249–252: the Sint-800 VADM record analysed throughout (799 points, errors $e_i$).
- Valet & Meynadier (1993), Nature 366, 234–238: 4 Myr intensity/reversal record; source of $\langle T_r\rangle\sim(2$–$3)\times10^5$ yr; not used because less reliable.
- Merrill, McElhinny & McFadden (1996), *The Magnetic Field of the Earth*: reversal-rate history; core conductivity $6\times10^5$ S m⁻¹ (p. 274).
- Moffatt (1978); Krause & Rädler (1980): the mean-field induction equation (Eq. 1).
- Brandenburg et al. (1989), A&A 213, 411: global α-quenching used in Eq. (2).
- Parker (1969), ApJ 158, 815; Levy (1972), ApJ 171, 635; Gibbons (1998), PEPI 106, 129: the Parker–Levy single-event reversal mechanism contrasted with the continuous-fluctuation picture.
- Olson (1983), PEPI 33, 260: reversals from jumps in the α-effect.
- Narteau et al. (2000), PEPI 120, 271: reversals from a steady level of helicity fluctuations in a multiple-scale model (more elaborate turbulence model than used here).
- Van Kampen (1992); Gardiner (1990): Fokker–Planck/Langevin theory; Gardiner p. 131 for the exact autocorrelation of the bistable process.
- Ossendrijver, Hoyng & Schmitt (1996), A&A 313, 938; Hoyng (1993, 1996): the same α-fluctuation mechanism applied to the (periodic) solar dynamo.
- Glatzmaier & Roberts (1995a,b); Kuang & Bloxham (1997); Olson et al. (1999); Christensen et al. (1999); Glatzmaier et al. (1999); Kageyama & Sato (1997): 3-D geodynamo simulations; Zhang & Schubert (2000) for the parameter-regime caveat.
- Roberts (1972); Deinzer et al. (1974); Sarson & Jones (1999); Hollerbach & Jones (1993); Gubbins (1999); Richards (2000): mechanisms for stationary dynamo action, the inner core, and the differential-rotation estimate.
- McFadden & Merrill (1995): cited for the observed anticorrelation between secular variation and $\langle T_r\rangle$ (**Note:** this item is cited in the text, p. 148, but is absent from the printed reference list).

## 15. Open questions and unresolved points

- Which of $\langle T_r\rangle$ (case 1) or reversal duration (case 2) is wrong, and why: a $\tau_d$ smaller by a factor 2–3 would reconcile everything; the authors leave open whether the model is inadequate or whether systematics in Sint-800 or in the construction of the artificial series are responsible (p. 154).
- The functional dependence of $\Delta_0$, $\Delta_1$ (and hence $\gamma$, $c^2$) on $kR$, $C_m$, $f^2\tau_c/N_c^2$ is unknown; without it one cannot go from fitted $(\gamma,c^2)$ to a simulation (pp. 147, 149).
- Relation (A.2) overestimates $\langle T_r\rangle$ by ~2 for case 1 but not for case 2; no accurate theoretical expression for $\langle T_r\rangle$ exists (p. 151).
- No optimal model was found: "many compromises are possible and we found no clear route to an optimal model" (p. 153).
- Origin of the apparent peaks near 4, 5, 6×10²² A m² in the Sint-800 histogram is unclear (p. 148).
- Whether "any kind of rapid rattling" produces reversals is likely but unproven (p. 147).
- Statistics not yet computed from $p(a)$: fraction of time below a threshold, frequency of excursions of given magnitude (p. 148).
- Future work: apply the analysis to the 4 Myr Valet & Meynadier (1993) data to constrain $\langle T_r\rangle$ better; use maximum-likelihood fitting to avoid binning; develop 2-D models with an independent radial coordinate and a conducting inner core (initial tests "look promising") to test whether reversal duration is set by the inner-core Ohmic time (Gubbins 1999) or by $\tau_d$ (p. 155).
- Mean-field applicability at $\beta/\eta\sim20$–30 "has never really been settled" (p. 155).
