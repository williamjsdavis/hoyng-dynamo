# Open questions

Things the papers leave open that a fast implementation can address with sweeps the authors could not afford. These are research questions, not scope: none of them blocks the build. Gate IDs refer to `VALIDATION.md`; *(derived here)* marks our own derivation. Status 2026-10-18: all open.

When a question is settled, write the answer under it with the runs that settled it and move it to the Settled section.

## Open

- **Q1. Scaling of Δ₁ (i.e., ⟨T₊²⟩ at a = 0) with f²τ_c/N_c².** HOS App. E assumes it constant, giving γ ∝ (f²τ_c/N_c²)⁻¹; Table II shows γf² drifting from 175 to 77 across runs 1–5. A dense sweep in forcing would give the functional form and remove the need to measure D₁ from each run.
- **Q2. The factor-2 gap between Eq. 58/A.2 and measured ⟨T_r⟩ at f = 5.4** (HSO p. 151), and the 12% gap at f = 6.4. Longer runs separate statistics from theory error; T3.8 separates the Gaussian approximation from the drift approximation.
- **Q3. Exactness of the f²τ_c/N_c² invariance** at finite τ_c and bin width (T4.12).
- **Q4. Resolution and time-step convergence** of the reported statistics (never reported; T0.11, T1.6, T4.17).
- **Q5. Relaxation rate near a = 1.** Deterministically it is 2μ = 1.526 *(derived here)*; HSO Eq. 15 uses 2Λ = 1.054 for the noisy autocorrelation and reports the simulated curve "follows Eq. 15 well at small lag" (Fig. 7d). T5.13 measures which rate the stochastic autocorrelation actually follows. This bears on HSO's τ_d = 11 ± 2 kyr, since what they measure is Λ/τ_d.
- **Q6. The Fig. 9 discrepancy** (T5.8) and whether Eq. A.1 removes it (T5.10).
- **Q7. Dependence on the distribution of F** (Gaussian vs uniform vs ±1) at fixed variance. The Fokker–Planck theory only sees the second moment; rare-event tails might not.
- **Q8. Seed-to-seed scatter** at fixed parameters (T4.18), to calibrate how much of the 519 vs 570 difference is realisation noise.
- **Q9. Reversal duration** with a fixed, documented definition across the forcing sweep, since HOS (~3 τ_d) and HSO (0.5–2 τ_d) differ.

## Settled

None yet.
