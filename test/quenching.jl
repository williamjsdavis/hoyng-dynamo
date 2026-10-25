# The quenching laws (dev/STRETCH_GOALS.md SG-2, gate T4.16): the papers' amplitude law is
# unchanged by the refactor (regression against values recorded before it), each law's
# coefficients match the dense formula, every law has the equilibrium a = ±1, and the
# inverse law relaxes at the predicted 2μ/(1 + s). Always run.

using LinearAlgebra, Random

@testset "quenching laws (SG-2)" begin
    model = Model(; kR = 0.5, Cm = 100.0, N = 60)
    g, cst = model.grid, model.constants
    M, s, Δτ = g.M, cst.s, 0.01
    G, E = model.G, quenching_operator(g)

    @testset "the amplitude law reproduces the trajectories recorded before the quenching switch" begin
        # a and Σb after n steps of the standard run (seed 1) and of the kR 0.3 run, recorded on
        # 2026-10-24 with the previous scalar quenching code; bitwise on the recording machine
        refs = ((:exponential, 10_000, -0.5137507890414534, 2.4549402704351837),
                (:imex, 2_000, 0.09218890861895997, 1.4228560711135876),
                (:backward_euler, 2_000, 0.278287644844738, 4.547502574318473),
                (:crank_nicolson_implicit_noise, 2_000, -0.7581345621479251, -10.628596361325453))
        for (stepper, n, a, Σb) in refs
            sim = Simulation(RunSpec(; stepper))
            run!(sim, n)
            @test isapprox(sim.state.a, a; rtol = 1e-9)
            @test isapprox(sum(sim.state.b), Σb; rtol = 1e-9)
            @info "amplitude-law regression, $stepper, $n steps" bitwise = (sim.state.a == a && sum(sim.state.b) == Σb)
        end
        sim = Simulation(RunSpec(; kR = 0.3, Cm = 120.0, f = 7.0, seed = 2))
        run!(sim, 5_000)
        @test isapprox(sim.state.a, -0.8963057338159734; rtol = 1e-9)
    end

    @testset "coefficients against the dense formulas" begin
        rng = Xoshiro(7)
        sp = ExponentialStepper(model, Δτ)
        cache = Cache(model)
        # a mixed state with overtones, and a frozen noise field
        b = 0.8 .* model.b0 .+ 0.3 .* real.(model.modes.b[:, 2]) .+ 0.05 .* randn(rng, 2M)
        v = 0.4 .* randn(rng, M)
        cache.v .= v
        a = dot(model.ã, b)
        V = zeros(2M, 2M)
        V[1:M, M+1:2M] .= Diagonal(v)
        Th = b[M+1:2M] .+ (Δτ / 2) .* (G * b)[M+1:2M]              # midpoint toroidal field
        # amplitude law: c = s(1 − a_h²) with a_h = a + (Δτ/2) ã ⋅ (s(1 − a²)E + V) b
        ah = a + (Δτ / 2) * dot(model.ã, (s * (1 - a^2) .* E .+ V) * b)
        c = quenching!(cache, AmplitudeQuenching(), sp, b, a)
        @test all(==(c[1]), c)
        @test c[1] ≈ s * (1 - ah^2) rtol = 1e-12
        # inverse law: c = s(1 − a_h²)/(1 + s a_h²), its own a_h
        cn = s * (1 - a^2) / (1 + s * a^2)
        ah = a + (Δτ / 2) * dot(model.ã, (cn .* E .+ V) * b)
        c = quenching!(cache, InverseQuenching(), sp, b, a)
        @test c[1] ≈ s * (1 - ah^2) / (1 + s * ah^2) rtol = 1e-12
        # energy law: c = s(1 − Σ w T_h²)
        c = quenching!(cache, EnergyQuenching(), sp, b, a)
        @test all(==(c[1]), c)
        @test c[1] ≈ s * (1 - sum(g.w .* Th .^ 2)) rtol = 1e-12
        @test !isapprox(c[1], s * (1 - ah^2); rtol = 1e-3)             # the overtones make it differ from the amplitude law
        # local law: c_j = s − (1 + s) κ T_h,j²
        c = quenching!(cache, LocalQuenching(0.55), sp, b, a)
        @test c ≈ s .- (1 + s) * 0.55 .* Th .^ 2 rtol = 1e-12
        # for a pure fundamental mode E_T = a² and G b₀ = 0, so the energy law's midpoint field is
        # the field itself and c = s(1 − a²) exactly; the amplitude law's midpoint a_h differs from
        # a at O(Δτ) because the nonlinearity acts on P, which a sees at once and T only via the shear
        fill!(cache.v, 0)
        b = 0.6 .* model.b0
        ce = quenching!(cache, EnergyQuenching(), sp, b, 0.6)[1]
        ca = quenching!(cache, AmplitudeQuenching(), sp, b, 0.6)[1]
        @test ce ≈ s * (1 - 0.36) rtol = 1e-10
        @test 0 < ce - ca < 5e-3
        # the same coefficients through the other steppers
        for other in (IMEXStepper(model, Δτ), BackwardEulerStepper(model, Δτ))
            @test quenching!(cache, LocalQuenching(0.55), other, b, 0.6) ≈ quenching!(cache, LocalQuenching(0.55), sp, b, 0.6) rtol = 1e-14
        end
    end

    function settle(law; b = model.b0, τ = 200.0, stepper = ExponentialStepper(model, Δτ; quenching = law))
        st, cache = State(model, b), Cache(model)
        for _ in 1:round(Int, τ / Δτ)
            step!(st, stepper, NoNoise(), cache)
        end
        return st
    end

    loc = calibrate(LocalQuenching, model, Δτ)
    @info "local quenching constant" κ = loc.κ twice_q0 = 2cst.q0

    @testset "equilibria: a = ±1 for every law ($name)" for (name, law) in
            (("amplitude", AmplitudeQuenching()), ("energy", EnergyQuenching()), ("inverse", InverseQuenching()), ("local", loc))
        st = settle(law)                                                 # from b₀
        @test abs(st.a - 1) < 1e-8
        st_r = settle(law; b = 0.1 .* randn(Xoshiro(2), 2M))            # from an arbitrary field
        @test abs(abs(st_r.a) - 1) < 1e-6
        st_m = settle(law; b = -model.b0)                                # the reversed polarity
        @test abs(st_m.a + 1) < 1e-8
        res = overtone_residual(model, st)
        P = HoyngDynamo.toroidal_parity(g, st.b)
        @info "equilibrium of the $name law" overtone_residual = res parity = P
        if law isa LocalQuenching
            @test res > 0.1                                              # the local law's equilibrium is not the fundamental mode
            @test P ≈ -1 atol = 1e-6                                     # but it is antisymmetric
            @test law.κ < 2cst.q0                                        # the sin θ-mean of q over b₀ would need κ = 2q₀
        else
            @test res < 1e-8                                             # the global laws keep b₀
        end
        # IMEX and backward Euler reach the same equilibrium amplitude to O(Δτ) accuracy
        for other in (IMEXStepper(model, Δτ; quenching = law), BackwardEulerStepper(model, Δτ; quenching = law))
            @test abs(settle(law; stepper = other).a - 1) < 1e-3
        end
    end

    @testset "relaxation rates about a = 1" begin
        function relaxation_rate(law; δ = 1e-6)
            sp = ExponentialStepper(model, Δτ; quenching = law)
            b_eq = settle(law).b
            st, cache = State(model, (1 + δ) .* b_eq), Cache(model)
            x = Float64[]
            for _ in 1:300
                step!(st, sp, NoNoise(), cache)
                push!(x, st.a - 1)
            end
            return -log(abs(x[300]) / abs(x[100])) / 2.0
        end
        rates = (amplitude = relaxation_rate(AmplitudeQuenching()), energy = relaxation_rate(EnergyQuenching()),
                 inverse = relaxation_rate(InverseQuenching()), local_law = relaxation_rate(loc))
        @info "relaxation rates about a = 1" rates... twoμ = 2cst.μ
        @test relerr(rates.amplitude, 2cst.μ) < 0.01
        @test relerr(rates.inverse, 2cst.μ / (1 + s)) < 0.01             # dc/da at a = 1 is −2s/(1 + s) instead of −2s
        @test rates.energy < rates.amplitude                             # T lags the quenching, which acts on P
        @test rates.local_law > rates.amplitude
    end

    @testset "the floor: α suppressed, never reversed" begin
        sp = ExponentialStepper(model, Δτ)
        cache = Cache(model)
        big = 3.0 .* model.b0 .+ 2.0 .* real.(model.modes.b[:, 2])             # a strong burst: E_T ≫ 1 + 1/s
        a = dot(model.ã, big)
        for (law, flaw) in ((AmplitudeQuenching(), AmplitudeQuenching(; floored = true)), (EnergyQuenching(), EnergyQuenching(; floored = true)),
                            (LocalQuenching(0.55), LocalQuenching(0.55; floored = true)))
            @test !floored(law) && floored(flaw)
            @test minimum(quenching!(cache, law, sp, big, a)) < -1                    # the α-effect has reversed
            @test minimum(quenching!(cache, flaw, sp, big, a)) == -1                  # clamped at α = 0
            @test maximum(quenching!(cache, flaw, sp, big, a)) ≤ s
        end
        @test !floored(InverseQuenching())
        @test minimum(quenching!(cache, InverseQuenching(), sp, big, a)) > -1          # h = 1/(1 + q) never reverses α
        # the floor is inactive at the equilibria: same coefficients and same local constant
        for (law, flaw) in ((AmplitudeQuenching(), AmplitudeQuenching(; floored = true)), (EnergyQuenching(), EnergyQuenching(; floored = true)), (loc, LocalQuenching(loc.κ; floored = true)))
            b_eq = settle(law).b
            @test quenching!(cache, law, sp, b_eq, dot(model.ã, b_eq)) == quenching!(cache, flaw, sp, b_eq, dot(model.ã, b_eq))
        end
        @test calibrate(LocalQuenching, model, Δτ; floored = true).κ == loc.κ
        # the unfloored energy and local laws run away at f = 6.4 (DR-8): seed 1 of the energy law
        # by τ ≈ 475, seed 2 of the local law by τ ≈ 116; the floored ones and the papers' law do not
        function diverges_before(spec, τ)
            sim = Simulation(spec)
            for _ in 1:round(Int, τ / 1.0)
                run!(sim, steps(sim, 1.0))
                (isfinite(amplitude(sim)) && abs(amplitude(sim)) < 20) || return true
            end
            return false
        end
        @test diverges_before(RunSpec(; quenching = :energy, seed = 1), 600)
        @test diverges_before(RunSpec(; quenching = :local, seed = 2), 200)
        @test !diverges_before(RunSpec(; quenching = :energy_floored, seed = 1), 600)
        @test !diverges_before(RunSpec(; quenching = :local_floored, seed = 2), 200)
        @test !diverges_before(RunSpec(; seed = 1), 600)
    end

    @testset "RunSpec and run files" begin
        @test_throws ArgumentError RunSpec(; quenching = :global)
        @test RunSpec(; quenching = :energy_floored).quenching === :energy_floored
        @test floored(quenching(Simulation(RunSpec(; quenching = :local_floored)).stepper))
        @test !floored(quenching(Simulation(RunSpec(; quenching = :local)).stepper))
        @test RunSpec().quenching === :amplitude
        for q in HoyngDynamo.QUENCHINGS
            sim = Simulation(RunSpec(; quenching = q))
            @test quenching(sim.stepper) isa AbstractQuenching
            run!(sim, 10)
            @test isfinite(sim.state.a)
        end
        @test quenching(Simulation(RunSpec(; quenching = :local)).stepper).κ ≈ loc.κ
        @test quenching(Simulation(RunSpec(; quenching = :inverse, stepper = :imex)).stepper) isa InverseQuenching
        spec = RunSpec(; quenching = :energy, seed = 3)
        @test RunSpec(spec; seed = 4).quenching === :energy
        @test occursin("energy", string(spec))
        path = joinpath(mktempdir(), "quenching.h5")
        file = RunFile(path, spec; overwrite = true)
        close(file)
        file = open_run(path)
        @test RunSpec(file) === spec
        @test metadata(file)["quenching"] == "energy"
        close(file)
    end
end
