# Implementation decisions

How the accepted model and gates get implemented. Scope (*what*) is fixed by `DEVELOPMENT.md` §1 (equations, parameter sets) and `VALIDATION.md` (gates, tolerance rules, operational definitions); nothing here changes a definition or a target. Each decision lists the options, a recommendation with its reason, which gates it touches, and a status. Status 2026-10-18: all decisions closed. Scope items S1–S3 closed; D1 (Julia), D7 (HDF5, two run modes), D8 (running aggregators; OnlineStats.jl optional, allocation-checked if used), D10 (MIT), D11 (testing, CI), D13 (Python `viz/`) and D14 (type stability) accepted by the user; D2–D6, D9 and D12 adopted as recommended without objection, to be revisited only if a gate fails. *(derived here)* marks our own analysis.

## 0. Scope check

Fixed: the model equation (HOS Eq. 26 form), the run table, the seven validation tiers, the tolerance rules, the operational definitions (reversal detection, variance window, dv, F distribution, bin assignment, sampling cadence).

Three items were scope rather than implementation. The user closed them on 2026-10-18:

- **S1. Sint-800 record.** Tier 6 is deferred and not planned for now. Obtaining the Guyodo & Valet (1999) VADM series becomes a task only if Tier 6 is reopened.
- **S2. Alternative quenching laws** (T4.16: q from total toroidal energy, α/(1+q), local quenching) are a stretch goal. Keep a quenching-variant switch in the parameters object (D10) so they can be added without restructuring, and implement them only after Tiers 0–5 pass.
- **S3. Form of the deliverable.** A personal research project, built as a package that others can install and use later. Consequences: package structure from the first commit, a public API, docstrings and a README, and a licence (D10). S3 constrains the language only to "installable by others", which every option in D1 satisfies.

## 1. Decisions

### D1. Language and numerical stack

Need: ~10⁹ time steps over the whole project (sweeps, seeds, the 5×10⁵ τ_d run for f = 4.4), so ≲ 1 µs per step; a dense non-symmetric eigensolver; quadrature, χ² fitting, KS test, autocorrelation; plots; reproducible RNG; thread-parallel sweeps.

| Option | For | Against |
|---|---|---|
| (a) Julia | One language for all tiers; compiled loops without ceremony; LAPACK, QuadGK, Optim, Distributions, HDF5, CairoMakie all mature; native threads; Xoshiro RNG | JIT latency per session; smaller user base if S3 is a public package |
| (b) Rust core + Python analysis | Fastest, no JIT; a Python wheel (PyO3/maturin) is the most reusable form for geophysics users; scipy/matplotlib for analysis | Two languages, FFI layer; eigenproblem and expm need faer/nalgebra or LAPACK bindings, or live on the Python side |
| (c) Python + NumPy with Numba for the stepper | Lowest friction if the user is Python-first; scipy/matplotlib; one language | Numba compile time and fragility; everything hot must be Numba-compatible; "modern/fast" only by courtesy |
| (d) Python + JAX | `lax.scan` stepper at ~µs/step, `vmap` over seeds/forcings for free, counter-based RNG built in | float64 needs a flag; CPU-only here; debugging inside `scan` is awkward |

Recommendation: (a) Julia, structured as a Julia package (installable from the repository URL with `Pkg.add`, registrable later). S3 only requires that others can install it, which (a) satisfies; choose (b) instead if reaching Python users matters more than a single language. Cost estimates below assume compiled code; (c) and (d) are within 2× of them.

Affects: everything. Status: accepted 2026-10-18, Julia for simulation, analysis and theory. The user's condition: the code must be type stable so the hot loop compiles to efficient machine code; see D14. Plotting and visualisation are in Python under `viz/` (D13), not in the Julia package.

### D2. Spatial discretisation

Grid: θ_j = jπ/(N−1), j = 0 … N−1, poles are nodes with P = T = 0; unknowns are the interior j = 1 … N−2, so 2(N−2) = 116 for N = 60. "60 interior points" is the same code with N = 62 (T0.11 tries both). Bin index ⌊θ_j N_c/π⌋ is evaluated at interior nodes only (closes G8).

Operator L, options:

- (a) Conservative second-order differences: (LP)_j = [sin θ_{j+½}(P_{j+1} − P_j) − sin θ_{j−½}(P_j − P_{j−1})] / (sin θ_j Δθ²) − P_j/sin²θ_j − (kR)² P_j. Self-adjoint under the discrete sin θ inner product, mirroring the continuous operator; tridiagonal.
- (b) Non-conservative second order (expand L = ∂θ² + cot θ ∂θ − 1/sin²θ − (kR)²). Same order; loses discrete self-adjointness.
- (c) Fourth-order stencils. Faster grid convergence; wider band; not needed for 1% gates.
- (d) Spectral in P¹_l(cos θ), where L is diagonal with eigenvalue −l(l+1) − (kR)². Exact diffusion, tridiagonal couplings in l for cos θ and ∂θ(sin θ ·); but the noise is piecewise constant in θ with a 1/√sin θ factor, so its spectral representation converges only algebraically, and comparability with the papers' 60-point grid is lost.

Shear term: central differences on the product, C [sin θ_{j+1}P_{j+1} − sin θ_{j−1}P_{j−1}] / (2Δθ). α-effect and noise: pointwise.

Recommendation: (a) with N a parameter; (d) later as an independent cross-check of Table I only (T0.2–T0.3), where the noise does not enter.

Affects: Tier 0 values at N = 60, T0.11. Status: adopted 2026-10-18 as recommended.

### D3. Time integrator

The noise is a piecewise-constant coefficient (redrawn every 5 steps), so the model is an ODE with time-dependent coefficients, not an SDE: ordinary ODE methods apply, there is no Itô/Stratonovich question, and the Fokker–Planck theory is the τ_c → 0 limit. Write Eq. 26 as ∂b/∂τ = G b + B(τ) b with G = G(C₀) constant and B = s(1 − a²)E + V(τ).

*(derived here)* E and V both read the T rows and write the P rows, so E² = V² = EV = VE = 0. Hence exp(Δτ B) = I + Δτ B exactly for frozen a and F: an explicit Euler update of the quenching and noise terms is not an approximation *for frozen coefficients*. **Correction (gate T1.4, 2026-10-18):** freezing a over the step is explicit Euler on the amplitude equation da/dτ = (ã, B b) and shifts the relaxation rate about a = ±1 by a relative μΔτ (0.8% at kR = 0.5, 1.5% at kR = 0.3, both measured). The production stepper therefore evaluates s(1 − a²) at a midpoint prediction a_{n+½} = a_n + (Δτ/2)(ã, B_n b_n), which is one extra dot product because (b̃₀, G b) = λ⁰a = 0 at the discrete C₀. With it the rate errors are second order: 4×10⁻⁴ (Λ) and 1×10⁻⁴ (2μ) at Δτ = 0.01, falling 100× per decade.

Options:

- (a) IMEX Euler: backward Euler on L (two tridiagonal solves per step), forward Euler on α-effect, shear, quenching, noise. Cost ~0.3 µs/step. First-order error in every eigenvalue of G, about ΛΔτ/2 ≈ 0.3% on Λ at Δτ = 0.01. Stability relies on the implicit diffusion damping the explicitly treated P→T coupling at high wavenumber; plausible at Δτ = 0.01 but must be shown (T1.6).
- (b) Exponential integrator on G plus exact update of B: b_{n+1} = X (I + Δτ B_n) b_n with X = exp(Δτ G) precomputed once by scaling-and-squaring (`expm`), dense 116×116. Cost ~3 µs/step, ~15 s per 5×10⁴ τ_d run. Exact for the linear dynamo at any Δτ, so b₀ is an exact fixed point (T1.1) and the integrator's spectrum is the eigensolver's by construction. *(derived here)* The leading-eigenvalue error of the split map vanishes at O(Δτ²) because ⟨b̃₀ᵐ, [G, sE] b₀ᵐ⟩ = 0 for the eigenvectors of G + sE, so T1.2 should pass with margin. Strang splitting costs nothing extra: with constant G the half-exponentials of consecutive steps merge, leaving one X per step with a evaluated at the staggered midpoint.
- (c) Crank–Nicolson on L with a second-order explicit scheme on the rest. Second order, but CN damps the stiffest modes weakly (amplification → −1), which with impulsive noise redraws leaves high-wavenumber ringing. Avoid.
- (d) Exponential of L only (two 58×58 dense products) with explicit coupling. Cheaper than (b), same stability question as (a).
- (e) Adaptive library ODE solvers. Awkward with the step-tied redraw schedule; no benefit.

**Constant memory.** The user requires arbitrarily long runs, so the stepper must run in O(1) memory with no history in RAM, and statistics must come from running aggregators (D8). This is orthogonal to the choice of stepper: every option above mutates a fixed state and workspace, and the only growing objects are stored time series, which D7 streams to disk or omits. The stepper API is `step!(state, cache, params, rng)`, mutating, allocation-free after warm-up (D14). Trade-offs:

| Stepper | Cost per step | Working memory | Error in Λ at Δτ = 0.01 | Remarks |
|---|---|---|---|---|
| (a) IMEX Euler | 1.4 µs measured (N = 60) | a few KB (tridiagonal factors, four vectors) | 1.3% measured (ΛΔτ, first order); 3.5–4.7% on 2μ because b̃₀ is not a left eigenvector of the IMEX amplification matrix | amplification matrix at C₀ has spectral radius 1 − 10⁻¹⁵ at Δτ = 0.01: stable (T1.6) |
| (b) exp(ΔτG) with midpoint-predicted B | 2.0 µs measured (N = 60; one 116×116 gemv, two dot products) | ~110 KB for X plus vectors; fits in L2 | 4×10⁻⁴ measured, second order | exact linear dynamics; b₀ a fixed point to 10⁻¹⁴ over 10⁴ steps; sign symmetry bitwise |
| (d) exp(ΔτL) with explicit coupling | ~0.6 µs (two 58×58 gemv) | ~55 KB | first order | same stability question as (a) |
| Galerkin in M eigenmodes (HOS Eq. 29) | O(M²) per step, O(N_c M²) per redraw | O(N_c M²) | truncation in M | not the papers' computation; a later research tool for Q1 |
| (e) DifferentialEquations.jl | adaptive | allocation-free only with care | adaptive | coefficients jump every 5 steps, which defeats adaptivity; no benefit |

Recommendation: (b) as the production stepper, (a) alongside as the cheap cross-check (same operator assembly). T1.6 and T4.17 compare them. Implemented 2026-10-18 as `ExponentialStepper` and `IMEXStepper` in `src/stepper.jl`; both type stable and allocation-free (D14 tests). Keep Eq. 26's form, G(C₀) with coefficient 1 + s(1 − a²), rather than G(C_m) − s a² E, so that b₀ is an eigenvector of the operator that is exponentiated. Do not build X from the eigendecomposition (non-normal, possibly ill-conditioned); use `expm`.

Affects: Tier 1, T4.17, run cost. Status: adopted 2026-10-18 as recommended.

### D4. Eigenproblem, adjoint modes, normalisation, C₀

- Solver: dense non-symmetric (LAPACK `geev`) on the 2(N−2) interior matrix; N up to 240 for T0.11 is still trivial.
- Adjoint modes: left eigenvectors of the discrete matrix, rescaled by the inverse quadrature weights so that (b̃_i, b_j) = Σ_j w_j sin θ_j (P̃_i* P_j + T̃_i* T_j) = δ_ij with trapezoid weights w_j = Δθ on the interior. Bi-orthogonality is then exact to round-off for any discretisation (T0.5); the discretisation only affects how close P̃₀, T̃₀ are to the continuous adjoint (T0.6, T0.7).
- Mode ordering: by decreasing ℜλ; complex pairs keep the ℑλ > 0 member (Table I convention); parity from the symmetric/antisymmetric energy fractions of T.
- Normalisation (closes G7): b₀ real; scale so that Σ_j w_j sin θ_j T₀_j² = 1; sign so that P₀ > 0 at the equator, which defines a = +1 as "normal" polarity; b̃₀ scaled by (b̃₀, b₀) = 1. E₀₀ and ∫ sin θ (P̃₀T₀)² dθ are invariant under this scaling.
- C₀: Brent's method on ℜλ⁰(C) with λ⁰ the leading eigenvalue, bracket [40, 80], tolerance 10⁻¹⁰. Use the discrete C₀(N) in the model, not 57.9, so that T1.1 holds exactly; report C₀(N) against 57.9 (T0.1). Define runs by (kR, C_m) and derive s = C_m/C₀(N) − 1, so that G(C₀) + sE ≅ G(C_m) holds in the discrete system (T0.9).

Recommendation: as listed. Affects: Tier 0, T1.1. Status: adopted 2026-10-18 as recommended.

### D5. Random numbers

- Generator: Xoshiro256** (Julia, Rust `rand_xoshiro`), PCG64 (NumPy) or Philox (JAX), following D1. All are fine; record the generator name and version with the output.
- Streams: one master seed per experiment; per-run seed derived by hashing (master seed, run index) so results are independent of thread scheduling and run order.
- Gaussian draws via the library's ziggurat; distribution switch gaussian | uniform(±√3) | binary(±1) for Q7.
- Draw order fixed: at each redraw, bins 0 … N_c−1 in order. F is drawn independently of f and multiplied by f at application, so one stream serves several f values (T2.5, HOS Fig. 5).

Recommendation: as listed. Affects: Tier 2, T4.12, T4.18. Status: adopted 2026-10-18 as recommended.

Implementation (2026-10-19, `src/noise.jl`): `Xoshiro(seed)` per noise object; `run_seed(master, index)` is a SplitMix64 hash of the pair written out in the package, so it is fixed across Julia versions; the distributions are the dispatch types `GaussianF`, `UniformF`, `BinaryF`. The F stream depends only on (seed, distribution, N_c, n_c). Julia documents that the sequence of a seeded generator may change between minor Julia versions, so the Julia version is part of the run metadata (D7).

### D6. State layout and step bookkeeping

- State: one vector b of length 2(N−2), P then T, Float64. Precomputed per node: sin θ, cos θ, noise factor g_j = 1/√(2N_c sin θ_j), bin index k_j.
- Projections: a = ã · b with ã = W b̃₀ precomputed (one dot product per step); likewise complex ã₁, ã₂ for the overtones.
- Schedules by integer step count n, never by floating τ: redraw when n mod n_c = 0 with n_c = τ_c/Δτ (require it to be an integer), snapshots when n mod n_d = 0 with n_d = 1/Δτ.
- Ramped forcing (SOH): f(τ) as a piecewise-linear table evaluated at each redraw, constant within a τ_c interval; schedules are small structs (`ConstantForcing`, `LinearRamp`) dispatched on type (D14).
- Step counter `Int64` (10¹⁸ steps), τ = nΔτ computed from it, never accumulated.
- Checkpoint = (b, n, RNG state, aggregator states), written to the HDF5 file at a fixed cadence so an arbitrarily long run can be stopped and resumed bitwise. The current F is not stored, so the cadence is a multiple of n_c and a resume starts at a redraw step (noted in `noise.jl`, 2026-10-19).

Recommendation: as listed. Affects: all runs. Status: adopted 2026-10-18 as recommended.

### D7. Run modes, output format, cadence, precision

**Run modes (user's decision).** Two modes, one stepping loop:

- **Record mode**: stepper, aggregators and stream writers. The time series (a every step, a₁ and a₂ every n_c steps, P and T snapshots every τ_d, any subset) are written to HDF5. For paper-length runs, figures and offline re-analysis.
- **Stream mode**: stepper and aggregators only; the output is the aggregator values and checkpoints. For arbitrarily long runs, sweeps and ensembles.

The mode is nothing more than the tuple of observers handed to `run!`: aggregators in both modes, writers added in record mode. The stepper and the RNG sequence do not know which mode is running, so the same seed and parameters give bitwise identical state trajectories and aggregator values in both, and replaying a recorded a(τ) through the aggregators reproduces the online values. Both invariants are tests (D11, T2.5).

Contents per run: metadata (every parameter, seed, generator, N, Δτ, integrator, F distribution, code git hash, date); the aggregator results (D8) and checkpoints (D6); optionally the streams a at every step, a₁, a₂ every n_c steps, and P and T snapshots every τ_d. Streams are written from fixed buffers into chunked, extendible datasets (e.g. 10⁴ steps per chunk), so RAM stays constant whatever the run length; a run with streams off keeps only the aggregators.

| Option | For | Against |
|---|---|---|
| HDF5, one file per run | Cross-language, chunked, attributes for metadata, readable from Python whatever D1 is | Dependency |
| Raw binary + JSON sidecar | Zero dependencies | Ad hoc layout, no partial reads |
| JLD2 / NPZ | Native, simple | Single-language |
| Arrow/Parquet | Good for event tables | Poor fit for the field snapshots |

Size: a every step in Float64 is 40 MB per 5×10⁴ τ_d, snapshots 46 MB; a 5×10⁵ τ_d run is ten times that. Proposal: Float32 for the per-step a stream (10⁻⁷ relative precision, far below anything measured from it; T1.1 is checked in memory), Float64 for everything else; a "light" mode without snapshots for long sweeps; `results/` git-ignored.

Recommendation: HDF5 with the layout above (HDF5.jl); `viz/` reads the same files with h5py. Affects: Tier 4–5 analysis, disk. Status: accepted 2026-10-18, including the two modes.

### D8. Running aggregators (constant-memory statistics)

The user's requirement: simulate arbitrarily long series and compute every statistic with running aggregators, so memory does not grow with run length.

- **Interface.** Each aggregator is a concrete struct with `fit!(agg, obs)` called once per step (slower cadences by an internal counter), `merge!(a, b)` to combine seeds or restart segments, and `value(agg)`. The observation is a view of the simulation state (b, a, step index); each aggregator derives what it needs. A run holds its aggregators in a `Tuple`, not a `Vector` of an abstract type (D14).
- **Aggregators the gates need, with their memory:**
  - mean and variance of a and of |a|, hence dv (Welford): O(1);
  - histogram of a on fixed bins (e.g. [−3, 3] in 600 bins) at τ_d cadence; HSO's equal-population binning is derived from it afterwards: O(bins);
  - event detector: Schmitt-trigger state, time of last crossing, counts of reversals, aborted reversals and excursions, fixed histograms of polarity-interval lengths and event durations, optionally the event list itself, which grows only with the number of reversals: O(1) + O(N_rev). The papers' grouping rule and the (a_h, W) sensitivity variants run as extra instances at negligible cost;
  - between-reversal variance with the W τ_d exclusion: samples pass through a delay ring buffer of W/Δτ entries before entering the current interval's Welford accumulator; at a reversal the buffer is discarded and the next W τ_d skipped; intervals are pooled time-weighted: O(W/Δτ);
  - D₁ accumulator: at each τ_d, F ∫ sin θ P̃₀² T₊² dθ added to the running mean of its |a| bin: O(bins);
  - autocorrelation of |a| − ⟨|a|⟩ to a fixed maximum lag (HSO Fig. 7): ring buffer and lagged sums: O(L_max);
  - parity, toroidal energy and |a₁|, |a₂| running statistics: O(1).
- **Precision over 10¹⁰ and more samples:** Welford or Kahan accumulation, Int64 counters, Float64 throughout.
- **Offline path.** A stored a(τ) stream (D7) is analysed by replaying it through the same aggregators, so online and offline results cannot disagree.
- **Ensembles.** `merge!` across seeds gives ensemble statistics and the seed-to-seed σ the tolerance rule needs.

Interface options:

- (a) Adopt OnlineStats.jl's `OnlineStat` interface (`fit!`, `merge!`, `value`, `Series`): Mean, Variance, Hist, AutoCov and Lag exist; the model-specific aggregators (event detector, exclusion variance, D₁ bins) are written as subtypes. One well-maintained dependency, mergeable by construction.
- (b) A bespoke interface with the same three verbs and hand-written Welford and histograms. Full control, no dependency, a few hundred lines to write and test.

Recommendation: the three-verb interface (`fit!`, `merge!`, `value`) is fixed. Whether it is OnlineStats.jl's or our own is decided when the aggregators are written, on whichever makes more sense then (dependency weight, allocation behaviour, fit with the state view); mixing is fine as long as everything composes and merges uniformly. Affects: Tier 4, Tier 5, memory of long runs. Status: accepted 2026-10-18. User's conditions: OnlineStats.jl is optional, not required; if any OnlineStats type is used, its methods are not trusted to be allocation-free, and the default test suite asserts `@allocated fit!(o, x) == 0` and `@allocated merge!(o, o′) == 0` after warm-up for each one; an allocating type is replaced by a bespoke implementation behind the same interface.

Implementation (2026-10-19, `src/aggregators.jl`): option (b), bespoke, no OnlineStats.jl. Reasons: the generic statistics are a few lines each; the model-specific aggregators are bespoke in any case; no dependency to allocation-audit; `merge!` extends `Base.merge!`, `fit!` and `value` are the package's own. Types: `Moments` (Welford), `Histogram` (fixed bins, NaN ignored), `AutoCorrelation` (lags 0 … L, ring buffer, lagged sums accumulated in 2²⁰-sample chunks, global-mean subtraction), `Sampled` (cadence by an internal counter), `MovingAverage` (a stream filter, `smooth!`), `SchmittTrigger` (events with kind, start, stop, min |a|, crossings, time of the last crossing; the reversal time is the last zero crossing; a jump across the band in one sample is a zero-duration reversal), `CrossingGroups` (the papers' rule), `ExclusionVariance` (delay ring buffer, own trigger, pooled by kept-sample count, open interval included), `BinnedMean` (Welford per bin) with `even_quadratic_fit` for HOS Fig. 8. Tuples of aggregators compose through `fit!`/`merge!`/`value`. Time-aware aggregators take the sampling interval Δt of the series they are fed. Merging pools independent runs; restart segments resume from the checkpointed state (D6). The event detectors' `merge!` appends event lists and so allocates; every other `fit!` and `merge!` is allocation-free (`test/quality.jl`).

### D9. Parallelism

Each run is single-threaded (the 116×116 product is below the size where BLAS threading helps; set BLAS threads to 1). Sweeps are a list of run specifications executed by a thread or process pool; a CLI taking one specification file also allows shell-level parallelism. Determinism independent of worker count follows from D5.

Recommendation: as stated. Affects: Tier 4 wall-clock. Status: adopted 2026-10-18 as recommended.

### D10. Architecture and repository layout

Julia package modules: `grid` (nodes, weights, L, shear), `eigen` (G(C), modes, adjoints, C₀), `noise` (F streams, schedules), `stepper` (D3 integrators), `aggregators` (D8), `io` (HDF5 writers, checkpoints), `run` (parameters and observers → output file), `theory` (Fokker–Planck formula module, Tier 3), `cli`. No plotting in the package (D13).

Parameters object: kR, C_m, f or f-schedule, τ_c, N_c, N, Δτ, run length, seed, integrator, F distribution, quenching variant (S2). C₀ and s are derived (D4).

Layout: `src/`, `test/`, `bench/`, `runs/` (committed run specifications), `results/` (ignored), `viz/` (Python plotting, D13), `figures/` (committed PNGs of validation figures, small), `data/` (Sint-800, if Tier 6 is ever reopened), `dev/`, `resources/`, `papers/` (ignored, unchanged). Each output file carries the git hash of the code and a hash of its specification.

Recommendation: as stated. Per S3, build it as a package from the first commit: package metadata, a public API (parameters → run → aggregator values and HDF5 output), docstrings, a README that carries the gate status, and a licence. Status: accepted 2026-10-18; licence MIT.

### D11. Testing, gate automation, CI

- Tests tagged by tier. Tiers 0–3 and the short Tier 1 integrations (10⁴ steps) form the default suite, target under two minutes. Tier 4 runs are a separate "slow" script that writes a results table (gate ID, ours, theirs, σ, pass/fail) using the same tolerance functions as the tests; its summary goes into `VALIDATION.md` as a status column.
- Tier 5: Python scripts in `viz/` (D13) render our figures from the HDF5 outputs to `figures/`, named after the paper figure (`HOS_fig01.png`, …); the paper's figures are not copied into the repository (copyright, and the `papers/` rule); the user compares against the PDF.
- CI (GitHub Actions) runs the default suite on push; the slow suite is manual.

Recommendation: as stated. Status: accepted 2026-10-18. The default suite also carries the type-stability and allocation tests of D14, the aggregator allocation checks of D8, and the mode-equivalence tests of D7 (record and stream modes bitwise identical; replay of a recorded stream reproduces the online aggregator values).

### D12. Theory-module numerics

Closed forms (Eqs. 52–60, HSO A.1–A.2) as pure functions of (γ, c², Λ, μ/Λ). Numerical counterparts (normalisation of Eq. 43, the escape integral Eq. 57, the dv integral) by adaptive quadrature with the integrand evaluated in log space relative to its maximum, since exp(2Φ) overflows for large γ. T3.1 and T3.8 cross-check closed form against quadrature. Tier 6 fitting: χ² over equal-population bins minimised by Nelder–Mead, errors by the same Monte-Carlo resampling HSO used.

Recommendation: as stated; QuadGK and Optim in Julia (the fits stay on the Julia side so that `viz/` is plotting only). Status: adopted 2026-10-18 as recommended.

Implementation (2026-10-20, `src/theory.jl`): `BistableTheory{T}(γ, c², Λ, κ)` with κ = μ/Λ − 1 selects between the HOS forms (κ = 0) and HSO's Eq. A.1 through one ψ(a) = γ[(κ(1 + c²) − 1)a² − ½κa⁴] + [γ(1 + c²)(1 − κc²) − 1] ln(1 + a²/c²); the peak, the variance and ⟨T_r⟩ are written in terms of a_m, ψ_m, ψ_m'' and φ₀'' (HSO Eq. A.2, with a_m from the quadratic of HSO p. 156) so that one code path reduces to Eqs. 52–58 at κ = 0 (checked to 10⁻¹²) and stays valid for κ ≠ 0. Every quadrature integrand is exp[ψ − max ψ]; ⟨T_r⟩, ⟨T_e⟩ and Z have `log_` forms that never overflow. The nested escape integral (Eq. 57) evaluates its inner integrand through QuadGK's `BatchIntegrand`: a plain `quadgk` inside a `quadgk` trips Julia's inference recursion limit (the adaptive routines are re-entered) and degrades to dynamic dispatch, while the batched inner call goes through different methods and stays fully inferred (D14), with the same result to the last bit. Constructors from `ModelConstants` connect the module to the eigensolver's Λ and μ, and `diffusion_from_fit`/`shape_parameters` implement Eqs. 61–63 for the Tier 4 D₁ measurement. Dependency: QuadGK 2. Optim is not added until Tier 6 fitting is planned.

### D13. Figures and visualisation (Python, `viz/`)

User's decision: plotting and visualisation live in Python, not Julia. `viz/` is a small, self-contained project managed by uv (`pyproject.toml`, `uv.lock`) with matplotlib, numpy and h5py; it reads the HDF5 outputs of D7 and is never a dependency of the simulation. One script or function per paper figure, paper axes and run parameters, a shared style module. T(θ, τ) colour plots use a diverging colour map with a > 0 mapped to the same colour sense as HOS Plate I. Interactive exploration (notebooks) also goes here.

Status: accepted 2026-10-18.

### D14. Performance discipline (type stability, allocations)

The user's requirement: type-stable code so that Julia compiles the hot loop to efficient machine code.

- Parametric concrete structs (`Params{T}`, `State{T}`, `Cache{T}`), no abstract or untyped fields; `T = Float64` by default.
- Hot path (`step!`, every aggregator's `fit!`, the noise redraw) allocation-free: preallocated workspace, in-place `mul!`, `randn!(rng, buf)`, `@views`; no closures over mutable locals; forcing schedules and the F-distribution switch as small structs dispatched on type.
- Heterogeneous collections (aggregators) as `Tuple`s, recursed with a type-stable `foreach`.
- Tests enforce it in the default suite: `Test.@inferred` on the hot functions, `@allocated step!(…) == 0` after warm-up, `JET.@test_opt` on `step!` and `fit!`. A `bench/` directory with BenchmarkTools records ns per step and per aggregator; a regression beyond ~20% is a CI warning.
- Linear algebra: `mul!` with BLAS threads set to 1, benchmarked against a plain `@inbounds @simd` loop for the 116×116 product; keep the faster.
- Latency: a PrecompileTools workload for the stepper and the eigensolver once the API settles.

Status: accepted 2026-10-18 as a requirement; the measures listed are the proposal.

## 2. Summary

| ID | Decision | Recommendation | Status |
|---|---|---|---|
| D1 | Language and stack | Julia, as a Julia package | accepted 2026-10-18 |
| D2 | Spatial discretisation | Conservative 2nd-order FD, poles as Dirichlet nodes, N parameter; spectral as Tier 0 cross-check | adopted 2026-10-18 |
| D3 | Time integrator | Exponential integrator exp(ΔτG) with exact (nilpotent) update of quenching and noise; IMEX Euler as cross-check; O(1)-memory `step!` | adopted 2026-10-18 |
| D4 | Eigenproblem and adjoints | LAPACK geev, weighted left eigenvectors, discrete C₀(N), s derived from C_m | adopted 2026-10-18 |
| D5 | RNG | Library generator per D1, hashed per-run seeds, Gaussian default with switch | adopted 2026-10-18; implemented 2026-10-19 |
| D6 | State and bookkeeping | Single Float64 vector, precomputed projection weights, Int64 step schedules, checkpointable | adopted 2026-10-18 |
| D7 | Run modes and output | Record mode (time series to HDF5) and stream mode (aggregators only), same results for the same seed; HDF5 per run with streamed extendible datasets and checkpoints | accepted 2026-10-18 |
| D8 | Running aggregators | Constant-memory online statistics behind a fit!/merge!/value interface; OnlineStats.jl optional, allocation-checked if used; offline = replay | accepted 2026-10-18; implemented bespoke 2026-10-19 |
| D9 | Parallelism | Single-threaded runs, pool over runs | adopted 2026-10-18 |
| D10 | Architecture | Module list above; package from the first commit; `runs/` committed, `results/` ignored; licence MIT | accepted 2026-10-18 |
| D11 | Testing and CI | Tier-tagged tests; slow script writes gate table; CI for default suite | accepted 2026-10-18 |
| D12 | Theory numerics | Log-space quadrature; closed forms cross-checked | adopted 2026-10-18; implemented 2026-10-20 |
| D13 | Figures | Python `viz/` with a uv environment (matplotlib, h5py); one script per paper figure; not a simulation dependency | accepted 2026-10-18 |
| D14 | Performance discipline | Type-stable, allocation-free hot path, enforced by tests and benchmarks | accepted 2026-10-18 |
