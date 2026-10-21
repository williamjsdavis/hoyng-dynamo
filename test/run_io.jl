# Run driver and HDF5 output (dev/DECISIONS.md D6, D7): the record/stream-mode and
# aggregator-replay parts of gate T2.5, checkpoint/resume, metadata round trip. Fast, always run.

const RUN_STEPS = 60_000          # 600 τ_d at Δτ = 0.01: seed 3 has 4 reversals and 2 aborted ones

@testset "Run driver and HDF5 output (D6, D7)" begin
    spec = RunSpec(; seed = 3)

    @testset "construction and extractors" begin
        sim = Simulation(spec)
        @test sim.state.n == 0 && sim.state.a ≈ 1
        @test amplitude(sim) == sim.state.a && absamplitude(sim) == abs(sim.state.a)
        @test abs(overtone1(sim)) < 1e-10 && abs(overtone2(sim)) < 1e-10          # pure fundamental: a_k = 0
        @test parity(sim) ≈ -1 atol = 1e-12                                         # dipole: T antisymmetric
        @test abs(d1_sample(sim)) < 1e-20                                           # T₊ = 0
        @test toroidal_energy(sim) ≈ sum(sim.model.grid.w .* toroidal(sim) .^ 2)
        @test length(toroidal(sim)) == length(poloidal(sim)) == sim.model.grid.M
        @test steps(sim, 1.0) == 100 && current_time(sim) == 0 && timestep(sim) == 0.01
        # a known overtone content: b = b₀ + c b₁ + c̄ b̄₁ (real) has a₁ = (b̃₁, b) = c and a = 1
        m = sim.model.modes
        idx = findall(λ -> imag(λ) ≥ 0, m.λ)
        c = 0.3 + 0.1im
        b = sim.model.b0 .+ 2 .* real.(c .* m.b[:, idx[2]])
        sim.state.b .= b
        sim.state.a = dot(sim.model.ã, b)
        @test overtone1(sim) ≈ c atol = 1e-10
        @test abs(overtone2(sim)) < 1e-10
        @test abs(sim.state.a - 1) < 1e-10
        @test sim.diag.λ[1] ≈ m.λ[idx[2]] && imag(sim.diag.λ[1]) > 0
        @test d1_sample(sim) ≈ sum(sim.model.grid.w .* view(sim.model.b̃0, 1:sim.model.grid.M) .^ 2 .* (toroidal(sim) .- sim.state.a .* sim.diag.T0) .^ 2)
        @test d1_pair(sim) == (absamplitude(sim), d1_sample(sim))
        @test occursin("seed = 3", sprint(show, spec)) && occursin("n = 0", sprint(show, sim))
        # the IMEX path and the ramp schedule build too
        simI = Simulation(RunSpec(; stepper = :imex, schedule = LinearRamp(4.0, 6.5, 0.0, 100.0)))
        run!(simI, 20)
        @test simI.state.n == 20 && simI.stepper isa IMEXStepper
        @test_throws ArgumentError RunSpec(; stepper = :rk4)
        @test_throws ArgumentError Probe(amplitude, Moments(), 0)
        @test_throws ArgumentError run!(sim, -1)
    end

    @testset "T2.5 record and stream modes agree bitwise; replay reproduces the aggregators" begin
        sim_s = Simulation(spec)
        obs_s = standard_observers(sim_s)
        run!(sim_s, RUN_STEPS, obs_s)
        path = tempname() * ".h5"
        sim_r = Simulation(spec)
        file = RunFile(path, spec)
        obs_r = standard_observers(sim_r; file)
        run!(sim_r, RUN_STEPS, obs_r)
        flush!(obs_r)
        write_results!(file, obs_r)
        summary = summarize(obs_r, sim_r)
        write_summary!(file, summary)
        @test sim_s.state.b == sim_r.state.b && sim_s.state.a == sim_r.state.a
        @test sim_s.state.n == sim_r.state.n == RUN_STEPS
        for k in keys(obs_s)
            @test isequal(value(obs_s[k]), value(obs_r[k]))
        end
        @test value(obs_s.trigger).n_reversals ≥ 1                      # the test is only meaningful with events
        # replay the recorded streams through fresh aggregators at the same cadence
        a = read_stream(file, "a")
        @test a.every == 1 && a.dt == 0.01 && length(a.values) == RUN_STEPS && a.values[end] == sim_s.state.a
        tr = SchmittTrigger(; Δt = 0.01)
        ex = ExclusionVariance(; Δt = 0.01)
        cg = CrossingGroups(; Δt = 0.01)
        h = Histogram(-2.5, 2.5, 1000)
        for (i, x) in enumerate(a.values)
            fit!(tr, x); fit!(ex, x); fit!(cg, x)
            i % 100 == 0 && fit!(h, x)
        end
        @test isequal(value(tr), value(obs_s.trigger))
        @test isequal(value(ex), value(obs_s.exclusion))
        @test isequal(value(cg), value(obs_s.groups))
        @test isequal(value(h), value(obs_s.histogram))
        d1 = read_stream(file, "d1")
        @test d1.every == 100 && length(d1.values) == RUN_STEPS ÷ 100
        bm = BinnedMean(0.0, 1.5, 30)
        for (i, y) in enumerate(d1.values)
            fit!(bm, (abs(a.values[100i]), y))
        end
        @test isequal(value(bm), value(obs_s.d1))
        a1 = read_stream(file, "a1")
        @test eltype(a1.values) == ComplexF64 && length(a1.values) == RUN_STEPS ÷ 5 && a1.values[end] == overtone1(sim_s)
        p = read_stream(file, "parity")
        @test p.every == 5 && p.values[end] == parity(sim_s)
        Tsnap = read_stream(file, "T")
        @test size(Tsnap.values) == (sim_s.model.grid.M, RUN_STEPS ÷ 100) && Tsnap.values[:, end] == toroidal(sim_s)
        @test Set(streams(file)) == Set(["a", "a1", "a2", "parity", "energy", "d1", "T", "P"])
        # results and summary are written and readable
        res = read_results(file)
        tv = value(obs_s.trigger)
        @test res["trigger"]["n_reversals"] == tv.n_reversals && res["trigger"]["reversal_times"] == tv.reversal_times
        @test res["trigger"]["events"]["kind"] == [String(e.kind) for e in tv.events]
        @test res["histogram"]["counts"] == value(obs_s.histogram).counts
        @test res["d1"]["mean"] == value(obs_s.d1).mean
        @test res["summary"]["n_reversals"] == tv.n_reversals && res["summary"]["theory"]["variance"] == summary.theory.variance
        @test summary.mean_reversal_time == tv.mean_interval
        @test summary.prefactor ≈ 3.217e-2 rtol = 1e-3
        # metadata round trip
        @test RunSpec(file) === spec
        md = metadata(file)
        @test md["seed"] == 3 && md["N"] == 60 && md["schedule_f"] == 6.4 && md["distribution"] == "GaussianF"
        @test read(file.h5["grid/theta"]) == sim_s.model.grid.θ
        @test occursin(path, sprint(show, file))
        close(file)
        @test !isopen(file)
        rm(path)
    end

    @testset "checkpoint and resume are bitwise" begin
        # reference: one straight run in record mode
        sim_a = Simulation(spec)
        path_a = tempname() * ".h5"
        file_a = RunFile(path_a, spec)
        obs_a = standard_observers(sim_a; file = file_a)
        run!(sim_a, RUN_STEPS, obs_a)
        flush!(obs_a)
        # interrupted run: checkpoints every 20000 steps, stopped at 45000, after the checkpoint at 40000
        sim_b = Simulation(spec)
        path_b = tempname() * ".h5"
        file_b = RunFile(path_b, spec)
        obs_b = standard_observers(sim_b; file = file_b)
        run!(sim_b, 45_000, obs_b; file = file_b, checkpoint_every = 20_000)
        flush!(obs_b)                                       # the streams hold 45000 samples, the checkpoint 40000
        @test checkpoint_step(file_b) == 40_000
        @test length(read_stream(file_b, "a").values) == 45_000
        close(file_b)
        # resume into the reopened file from the stored specification
        file_b = open_run(path_b; mode = "r+")
        sim_c = Simulation(RunSpec(file_b))
        obs_c = standard_observers(sim_c; file = file_b, attach = true)
        restore!(sim_c, obs_c, file_b)
        @test sim_c.state.n == 40_000 && current_time(sim_c) == 400
        @test length(read_stream(file_b, "a").values) == 40_000 && size(read_stream(file_b, "T").values, 2) == 400
        run!(sim_c, RUN_STEPS - 40_000, obs_c)
        flush!(obs_c)
        @test sim_c.state.b == sim_a.state.b && sim_c.state.a == sim_a.state.a && sim_c.state.n == RUN_STEPS
        for k in keys(obs_a)
            @test isequal(value(obs_c[k]), value(obs_a[k]))
        end
        for name in streams(file_a)
            @test read_stream(file_b, name).values == read_stream(file_a, name).values
        end
        close(file_a); close(file_b)
        rm(path_a); rm(path_b)
        # refusals: cadence off the redraw grid, checkpoint without a file, checkpoint off a redraw step
        sim_d = Simulation(spec)
        path_d = tempname() * ".h5"
        file_d = RunFile(path_d, spec)
        @test_throws ArgumentError run!(sim_d, 10, nothing; file = file_d, checkpoint_every = 7)
        @test_throws ArgumentError run!(sim_d, 10, nothing; checkpoint_every = 5)
        run!(sim_d, 3)
        @test_throws ArgumentError checkpoint!(file_d, sim_d, nothing)
        @test checkpoint_step(file_d) === nothing
        @test_throws ArgumentError restore!(sim_d, nothing, file_d)
        close(file_d); rm(path_d)
    end

    @testset "short runs, empty results, refusals" begin
        sim = Simulation(spec)
        path = tempname() * ".h5"
        file = RunFile(path, spec)
        obs = standard_observers(sim; file, snapshots = false, stream_eltype = Float32)
        run!(sim, 50, obs)
        flush!(obs)
        write_results!(file, obs)
        write_summary!(file, summarize(obs, sim))               # nothing to fit yet: NaNs, no error
        res = read_results(file)
        @test res["trigger"]["n_reversals"] == 0 && isempty(res["trigger"]["reversal_times"]) && isempty(res["trigger"]["events"]["kind"])
        @test eltype(read_stream(file, "a").values) == Float32 && length(read_stream(file, "a").values) == 50
        @test eltype(read_stream(file, "a1").values) == ComplexF32
        @test length(read_stream(file, "d1").values) == 0 && !("T" in streams(file))
        @test_throws ArgumentError Recorder(file, "a", amplitude, 1; sim)                  # exists
        @test_throws ArgumentError Recorder(file, "zzz", amplitude, 1; sim, attach = true)  # missing
        @test_throws ArgumentError Recorder(file, "a", amplitude, 7; sim, attach = true)    # wrong cadence
        close(file)
        @test_throws ArgumentError RunFile(path, spec)                                      # no overwrite
        f2 = RunFile(path, spec; overwrite = true)
        close(f2)
        rm(path)
    end
end
