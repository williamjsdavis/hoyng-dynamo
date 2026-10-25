using HoyngDynamo
using Test

# Tier selection, see dev/VALIDATION.md. Default: the fast tiers 0–3 and the Sint-800
# comparison (Tier 6, under a minute). Tier 4 is slow and runs on request:
#     HOYNG_TIERS=0,1,2,3,4,6 julia --project -e 'using Pkg; Pkg.test()'
const TIERS = Set(parse.(Int, split(get(ENV, "HOYNG_TIERS", "0,1,2,3,6"), ',')))
tier(t) = t in TIERS
relerr(x, ref) = abs(x - ref) / abs(ref)        # shared by the tier files

@testset "HoyngDynamo" begin
    @testset "package loads" begin
        @test HoyngDynamo isa Module
    end
    tier(0) && include("tier0_eigenproblem.jl")
    tier(1) && include("tier1_deterministic.jl")
    tier(2) && include("tier2_noise.jl")
    tier(3) && include("tier3_theory.jl")
    tier(4) && include("tier4_stochastic.jl")
    tier(6) && include("tier6_sint800.jl")
    include("stepper.jl")
    include("quenching.jl")
    include("aggregators.jl")
    include("conditional_moments.jl")
    include("stats.jl")
    include("run_io.jl")
    include("quality.jl")
end
