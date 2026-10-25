# The θ-method steppers (dev/STRETCH_GOALS.md SG-1): the block-elimination solve against the
# dense θ-method formula, the analytic per-mode factors, and the agreement of every stepper
# with the exponential one to the order of its truncation error over two diffusion times
# with noise. Always run; the Tier 1 gates also run on backward Euler and Crank–Nicolson.

using LinearAlgebra, Random

@testset "θ-method steppers (SG-1)" begin
    model = Model(; kR = 0.5, Cm = 100.0, N = 60)
    g = model.grid
    M = g.M
    Δτ = 0.01
    G, E = model.G, quenching_operator(g)
    rng = Xoshiro(11)

    @testset "one step equals the dense θ-method (θ = $θ, implicit noise $ι)" for θ in (1, 1 // 2), ι in (false, true)
        sp = ThetaStepper(model, Δτ; θ, implicit_noise = ι)
        @test timestep(sp) == Δτ
        b = 0.3 .* randn(rng, 2M)
        v = 0.5 .* randn(rng, M)
        st = State(model, b)
        cache = Cache(model)
        cache.v .= v                                     # NoNoise leaves cache.v as set
        c = HoyngDynamo.quenching_coefficient(sp, b, v, st.a)
        V = zeros(2M, 2M)
        V[1:M, M+1:2M] .= Diagonal(v)
        A = G .+ c .* E
        α, β = θ * Δτ, (1 - θ) * Δτ
        expected = ι ? (I - α * (A + V)) \ ((I + β * (A + V)) * b) : (I - α * A) \ ((I + β * A) * b + Δτ * (V * b))
        step!(st, sp, NoNoise(), cache)
        @test st.b ≈ expected rtol = 1e-12
        @test st.a ≈ dot(model.ã, expected) rtol = 1e-12
        @test st.n == 1
    end

    @testset "constructor" begin
        @test_throws ArgumentError ThetaStepper(model, Δτ; θ = 0.4)
        @test_throws ArgumentError ThetaStepper(model, Δτ; θ = 1.1)
        @test BackwardEulerStepper(model, Δτ).θ == 1
        @test CrankNicolsonStepper(model, Δτ).θ == 0.5
        @test CrankNicolsonStepper(model, Δτ; implicit_noise = true).implicit_noise
        @test_throws ArgumentError RunSpec(; stepper = :forward_euler)
        for s in HoyngDynamo.STEPPERS
            @test Simulation(RunSpec(; stepper = s)).stepper isa AbstractStepper
        end
        @test Simulation(RunSpec(; stepper = :backward_euler_implicit_noise)).stepper.implicit_noise
    end

    @testset "amplification matrix and per-mode factors" begin
        λ = model.modes.λ[1:6]                          # fundamental and the first overtone pairs
        for (sp, order) in ((BackwardEulerStepper(model, Δτ), 1), (CrankNicolsonStepper(model, Δτ), 2), (ExponentialStepper(model, Δτ), 0))
            ev = eigvals(amplification_matrix(model, sp))
            for λk in λ
                gk = mode_factor(sp, λk)
                @test minimum(abs.(ev .- gk)) < 1e-9     # the predicted factor is an eigenvalue of the one-step map
                order == 0 || @test abs(gk - exp(Δτ * λk)) < 2 * abs(Δτ * λk)^(order + 1) + 1e-14   # the floor is for λ₀ ≈ 0
            end
        end
        λ1 = λ[2]
        exact = abs(exp(Δτ * λ1))
        be, cn = abs(mode_factor(BackwardEulerStepper(model, Δτ), λ1)), abs(mode_factor(CrankNicolsonStepper(model, Δτ), λ1))
        @info "per-step damping of the first overtone at C₀" exact backward_euler = be crank_nicolson = cn
        @test be < exact                                 # backward Euler over-damps
        @test abs(cn - exact) < 1e-4                     # Crank–Nicolson is second order: |Δτλ₁|³/12 ≈ 2×10⁻⁵
        @test maximum(abs, eigvals(amplification_matrix(model, IMEXStepper(model, Δτ)))) ≤ 1 + 1e-9
    end

    @testset "all steppers agree with the exponential one over 2 τ_d with noise" begin
        steppers = (imex = IMEXStepper(model, Δτ),
                    backward_euler = BackwardEulerStepper(model, Δτ),
                    backward_euler_implicit_noise = BackwardEulerStepper(model, Δτ; implicit_noise = true),
                    crank_nicolson = CrankNicolsonStepper(model, Δτ),
                    crank_nicolson_implicit_noise = CrankNicolsonStepper(model, Δτ; implicit_noise = true),
                    exponential = ExponentialStepper(model, Δτ))
        traces = map(steppers) do sp
            noise = BinnedNoise(model, Δτ; f = 6.4, seed = 5)
            st = State(model, 0.9 .* model.b0)
            cache = Cache(model)
            a = Float64[]
            for _ in 1:200
                step!(st, sp, noise, cache)
                push!(a, st.a)
            end
            (a = a, b = copy(st.b))
        end
        ref = traces.exponential
        devs = map(tr -> maximum(abs, tr.a .- ref.a), traces)
        @info "maximum deviation in a from the exponential stepper over 2 τ_d" devs...
        for (name, tr) in pairs(traces)
            @test all(isfinite, tr.b)
            # Only Crank–Nicolson with trapezoidal noise is second order throughout; with the
            # noise term explicit it is forward Euler on that term, first order like the rest.
            @test devs[name] < (name === :crank_nicolson_implicit_noise ? 1e-2 : 1e-1)
        end
        @test devs.crank_nicolson_implicit_noise < minimum(d for (n, d) in pairs(devs) if n ∉ (:exponential, :crank_nicolson_implicit_noise))
        @test traces.crank_nicolson.a != traces.crank_nicolson_implicit_noise.a      # the noise treatment matters
        @test traces.backward_euler.a != traces.backward_euler_implicit_noise.a
    end
end
