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
- `stepper.jl`      constant-memory time steppers (exponential integrator, IMEX Euler)
- `aggregators.jl`  running statistics behind `fit!`, `merge!`, `value`
- `io.jl`           HDF5 streams and checkpoints
- `run.jl`          parameters and observers → a run in record or stream mode
- `theory.jl`       the Fokker–Planck formula module (closed forms and quadratures)

Plotting lives outside the package, in `viz/` (Python).
"""
module HoyngDynamo

using LinearAlgebra
using Random

export Grid, inner, diffusion_operator, shear_operator, alpha_operator, dynamo_operator,
    quenching_operator
export Modes, eigenmodes, fundamental, leading_eigenvalues, biorthogonality_residual,
    toroidal_parity, poloidal_parity, growth_rate, critical_dynamo_number, overlap_E00,
    diffusion_integral, ModelConstants, model_constants

include("grid.jl")
include("eigen.jl")

end # module
