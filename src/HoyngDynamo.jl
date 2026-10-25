"""
    HoyngDynamo

A fast, modern implementation of the one-dimensional (colatitude) stochastic αΩ
mean-field geodynamo of Hoyng, Ossendrijver & Schmitt (2001), Schmitt, Ossendrijver &
Hoyng (2001) and Hoyng, Schmitt & Ossendrijver (2002): the "bistable oscillator"
model of the geodynamo.

The model, the validation gates and the implementation decisions are documented in
`dev/` (`DEVELOPMENT.md`, `VALIDATION.md`, `DECISIONS.md`). Notation follows the 2001
GAFD paper ("HOS") unless a docstring says otherwise.

Package layout (`dev/DECISIONS.md` D10), in build order:

- `grid.jl`         colatitude nodes, quadrature weights, the operators L and shear
- `eigen.jl`        G(C), eigenmodes and weighted adjoints, the critical dynamo number C₀
- `noise.jl`        the stochastic α: binned, piecewise-constant F(θ, τ), forcing schedules
- `stepper.jl`      constant-memory time steppers (exponential integrator, IMEX Euler, θ-method) and the quenching laws
- `aggregators.jl`  running statistics behind `fit!`, `merge!`, `value`
- `theory.jl`       the Fokker–Planck formula module (closed forms and quadratures)
- `stats.jl`        statistics for the gates: KS test for exponential intervals, χ² fit of p(a) to a histogram
- `sint800.jl`      the Sint-800 record (`data/`) and HSO's analysis of it: dv, the VADM χ² fit, autocorrelation, τ_d
- `run.jl`          `RunSpec`, `Simulation`, observers and the `run!` loop (record or stream mode)
- `io.jl`           HDF5 run files: streams, results, checkpoints

Plotting lives outside the package, in `viz/` (Python). Investigations built on the package live
in sibling repositories (`hoyng-dynamo-sde`: drift and diffusion functions from streamed runs).
"""
module HoyngDynamo

using LinearAlgebra
using Random
using QuadGK: quadgk, BatchIntegrand
using HDF5
using Serialization: serialize, deserialize
import Dates

export Grid, inner, diffusion_operator, shear_operator, alpha_operator, dynamo_operator,
    quenching_operator
export Modes, eigenmodes, fundamental, leading_eigenvalues, biorthogonality_residual,
    toroidal_parity, poloidal_parity, growth_rate, critical_dynamo_number, overlap_E00,
    diffusion_integral, ModelConstants, model_constants
export NoNoise, update_noise!, AbstractFDistribution, GaussianF, UniformF, BinaryF,
    AbstractForcing, ConstantForcing, LinearRamp, forcing, run_seed, BinnedNoise, redraw!,
    redraw_interval, noise_projection, binned_D0
export Model, State, Cache, AbstractStepper, ExponentialStepper, IMEXStepper, ThetaStepper, BackwardEulerStepper,
    CrankNicolsonStepper, step!, timestep, dipole_amplitude, overtone_residual, amplification_matrix, mode_factor,
    AbstractQuenching, AmplitudeQuenching, EnergyQuenching, InverseQuenching, LocalQuenching, quenching!, quenching, calibrate, floored
export Aggregator, fit!, value, Moments, Histogram, AutoCorrelation, Sampled, MovingAverage, smooth!,
    Series, Event, SchmittTrigger, outside_events, CrossingGroups, ExclusionVariance, BinnedMean, even_quadratic_fit
export BistableTheory, diffusion_constants, shape_parameters, diffusion_prefactor, diffusion_from_fit,
    drift, diffusion, potential, psi, dpsi, d2psi, phi, phi0_curvature, bistable, peak_amplitude,
    peak_log_density, peak_curvature, mean_amplitude, relative_variance, log_mean_reversal_time,
    mean_reversal_time, kramers_exponent, log_mean_reversal_time_asymptotic, mean_reversal_time_asymptotic,
    psi_max, log_normalisation, StationaryDensity, stationary_density, logpdf, half_moments,
    dipole_variation, log_escape_time, escape_time, log_mean_reversal_time_quadrature,
    mean_reversal_time_quadrature
export RunSpec, Diagnostics, Simulation, current_time, steps, amplitude, absamplitude, overtone, overtone1, overtone2,
    parity, toroidal, poloidal, toroidal_energy, d1_sample, d1_pair, noise_diffusion_sample, continuum_diffusion_sample,
    noise_diffusion_pair, continuum_diffusion_pair, Observer, Probe, observe!, flush!, run!,
    standard_observers, summarize, pool!, ensemble, segmented, window_run
export ks_exponential, nelder_mead, fold_histogram, equal_population_bins, density_fit
export read_sint800, montecarlo, hso_autocorrelation, running_mean, artificial_sint800, tau_d_scan, vadm_bins, vadm_chi2, vadm_fit
export RunFile, open_run, metadata, Recorder, SnapshotRecorder, read_stream, streams, write_results!, write_summary!,
    read_results, checkpoint!, restore!, checkpoint_step

include("grid.jl")
include("eigen.jl")
include("stepper.jl")
include("noise.jl")
include("aggregators.jl")
include("theory.jl")
include("stats.jl")
include("sint800.jl")
include("run.jl")
include("io.jl")

end # module
