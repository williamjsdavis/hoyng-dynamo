# Tier 1: deterministic nonlinear dynamics, f = 0. Gates T1.1–T1.6 in dev/VALIDATION.md.
# Reference values: HOS p. 273 (fixed point), Eq. 34 (Λ = 0.527), HSO p. 156 (Λ = 0.815);
# 2μ = 1.526 (kR 0.5) and 2.97 (kR 0.3) are our derivation (DEVELOPMENT.md §1.1, T1.4).

using LinearAlgebra, Random


function integrate!(st, sp, cache, nsteps)
    for _ in 1:nsteps
        step!(st, sp, NoNoise(), cache)
    end
    return st
end

# a(τ) at steps 0 … nsteps
function trace_a(st, sp, cache, nsteps)
    a = Vector{Float64}(undef, nsteps + 1)
    a[1] = st.a
    for k in 1:nsteps
        step!(st, sp, NoNoise(), cache)
        a[k+1] = st.a
    end
    return a
end

# exponential rate of |x| between sample indices i1 and i2 (steps i1−1 and i2−1)
rate_between(x, Δτ, i1, i2) = log(abs(x[i2]) / abs(x[i1])) / ((i2 - i1) * Δτ)

# T1.2 protocol: tiny random field, rate of a between τ = 10 and τ = 15 (overtones dead, a ≪ 1)
function measured_growth_rate(model, sp; rng = Xoshiro(1))
    Δτ = timestep(sp)
    st = State(model, 1e-10 .* randn(rng, 2 * model.grid.M))
    a = trace_a(st, sp, Cache(model), round(Int, 15 / Δτ))
    return rate_between(a, Δτ, round(Int, 10 / Δτ) + 1, round(Int, 15 / Δτ) + 1)
end

# T1.4 protocol: b = (1 + δ) b₀, decay rate (positive) of a − 1 between τ = 1 and τ = 3.
# δ = 10⁻⁶ keeps the nonlinear correction (relative O(δ)) below the Δτ² splitting error.
function measured_relaxation_rate(model, sp; δ = 1e-6)
    Δτ = timestep(sp)
    st = State(model, (1 + δ) .* model.b0)
    a = trace_a(st, sp, Cache(model), round(Int, 3 / Δτ))
    return -rate_between(a .- 1, Δτ, round(Int, 1 / Δτ) + 1, round(Int, 3 / Δτ) + 1)
end

# The 1% gates apply to the production stepper and the second-order Crank–Nicolson; the
# first-order IMEX and backward Euler cross-checks are held to 5%.
second_order(name) = Symbol(name) in (:exponential, :crank_nicolson)
gate_tol(name) = second_order(name) ? 0.01 : 0.05

@testset "Tier 1: deterministic nonlinear dynamics, f = 0" begin
    cases = ((kR = 0.5, Cm = 100.0, Λ = 0.527, twoμ = 1.526), (kR = 0.3, Cm = 120.0, Λ = 0.815, twoμ = 2.97))
    for case in cases
        model = Model(; case.kR, case.Cm, N = 60)
        c = model.constants
        M = model.grid.M
        steppers = (exponential = ExponentialStepper(model, 0.01), imex = IMEXStepper(model, 0.01),
                    backward_euler = BackwardEulerStepper(model, 0.01), crank_nicolson = CrankNicolsonStepper(model, 0.01))

        for (name, sp) in pairs(steppers)
            label = "$name, kR = $(case.kR), C_m = $(case.Cm)"

            @testset "T1.1 b₀ is a fixed point over 100 τ_d ($label)" begin
                st = integrate!(State(model, model.b0), sp, Cache(model), 10_000)
                @test abs(st.a - 1) < 1e-8
                @test norm(st.b - model.b0) < 1e-8
                @info "T1.1 $label" a_minus_1 = st.a - 1 field_drift = norm(st.b - model.b0)
            end

            @testset "T1.2 small-amplitude growth rate Λ ($label)" begin
                rate = measured_growth_rate(model, sp)
                @test relerr(rate, case.Λ) < gate_tol(name)
                @info "T1.2 $label" rate paper = case.Λ discrete_Λ = c.Λ err_vs_discrete = relerr(rate, c.Λ)
            end

            @testset "T1.3 saturation at |a| = 1 from an arbitrary field ($label)" begin
                st = integrate!(State(model, 0.1 .* randn(Xoshiro(2), 2M)), sp, Cache(model), 5_000)
                @test abs(abs(st.a) - 1) < 1e-6
                @test overtone_residual(model, st) < 1e-6
            end

            @testset "T1.4 relaxation about a = 1 at rate 2μ ($label)" begin
                rate = measured_relaxation_rate(model, sp)
                @test relerr(rate, case.twoμ) < gate_tol(name)
                @test relerr(rate, 2c.Λ) > 0.2           # it is not 2Λ (QUESTIONS.md Q5)
                @info "T1.4 $label" rate twoμ_paper = case.twoμ twoμ_discrete = 2c.μ twoΛ = 2c.Λ
            end

            @testset "T1.5 sign symmetry b → −b ($label)" begin
                b = 0.1 .* randn(Xoshiro(3), 2M)
                a_plus = trace_a(State(model, b), sp, Cache(model), 2_000)
                a_minus = trace_a(State(model, -b), sp, Cache(model), 2_000)
                @test maximum(abs, a_plus .+ a_minus) ≤ 1e-12
                @info "T1.5 $label" bitwise = (a_plus == -a_minus)
            end
        end

        @testset "T1.6 time-step convergence (kR = $(case.kR))" begin
            for (name, make) in (("exponential", Δτ -> ExponentialStepper(model, Δτ)), ("imex", Δτ -> IMEXStepper(model, Δτ)),
                                 ("backward_euler", Δτ -> BackwardEulerStepper(model, Δτ)), ("crank_nicolson", Δτ -> CrankNicolsonStepper(model, Δτ)))
                errs = map((0.01, 0.001)) do Δτ
                    sp = make(Δτ)
                    (; Δτ, growth = relerr(measured_growth_rate(model, sp), c.Λ),
                        relaxation = relerr(measured_relaxation_rate(model, sp), 2c.μ))
                end
                ratios = (growth = errs[1].growth / errs[2].growth, relaxation = errs[1].relaxation / errs[2].relaxation)
                @info "T1.6 $name kR = $(case.kR)" errs ratios
                if second_order(name)
                    # second order: the growth-rate error falls ×100 per decade. The relaxation
                    # error cannot show that: at Δτ = 0.001 it sits on the measurement floor
                    # (round-off ~10⁻¹⁴ relative to the decayed perturbation δ e^{−2μ·3} ~ 10⁻¹⁰),
                    # so only its size at Δτ = 0.01 and its monotone decrease are asserted.
                    @test errs[1].growth < 1e-3 && errs[1].relaxation < 1e-3
                    @test ratios.growth > 30
                    @test errs[2].relaxation < errs[1].relaxation
                else                          # first order: ×~10 per decade in the growth rate
                    @test errs[1].growth < 5e-2 && errs[1].relaxation < 5e-2
                    @test 5 < ratios.growth < 20
                    if name == "imex"             # overtones leak into a at O(Δτ): b̃₀ is not a left eigenvector of the IMEX map
                        @test 5 < ratios.relaxation < 20
                    else                          # backward Euler keeps b̃₀ as a left eigenvector (b̃₀ᵀG = 0), and the
                        # midpoint-predicted quenching makes its relaxation rate second order like the exponential's
                        @test errs[1].relaxation < 1e-3 && errs[2].relaxation < errs[1].relaxation
                    end
                end
            end
            # stability diagnostic for IMEX Euler at Δτ = 0.01: amplification matrix at C₀
            g = model.grid
            L = diffusion_operator(g, c.kR)
            Lbd = [Matrix(L) zeros(M, M); zeros(M, M) Matrix(L)]
            amp = (I - 0.01 * Lbd) \ (I + 0.01 * (model.G - Lbd))
            μs = abs.(eigvals(amp))
            @info "T1.6 IMEX amplification at C₀, kR = $(case.kR)" max_modulus = maximum(μs) n_above_one = count(>(1 + 1e-9), μs)
            @test maximum(μs) ≤ 1 + 1e-9
        end
    end
end
