# Tier 4: stochastic runs, quantitative (slow). Gates T4.1–T4.18 in dev/VALIDATION.md.
# Build step 6 covers the standard run (T4.1–T4.4, T4.13 as a diagnostic); the sweeps and
# the other parameter sets (T4.5–T4.18) are build step 7. Definitions of VALIDATION.md §2:
# Schmitt trigger a_h = 0.5, a_e = 0.3; exclusion variance W = 5 τ_d; dv on |a| at 1 τ_d.

const T4_LENGTH = 5e4                 # τ_d per seed, as in the papers
const T4_SEEDS = (1, 2, 3)

# Pool the aggregators of several runs with the standard observers into the first.
function pool!(target::NamedTuple, others)
    for obs in others, k in keys(target)
        merge!(target[k].agg, obs[k].agg)
    end
    return target
end

@testset "Tier 4: stochastic runs, quantitative (slow)" begin
    @testset "T4.1–T4.4 standard run (kR 0.5, C_m 100, f 6.4): $(length(T4_SEEDS)) seeds × $(T4_LENGTH) τ_d" begin
        runs = map(T4_SEEDS) do seed
            sim = Simulation(RunSpec(; seed))
            obs = standard_observers(sim)
            t = @elapsed run!(sim, steps(sim, T4_LENGTH), obs)
            s = summarize(obs, sim)
            @info "T4 seed $seed ($(round(t, digits = 1)) s): ⟨T_r⟩ = $(round(s.mean_reversal_time, digits = 1)) ± $(round(s.mean_reversal_time_error, digits = 1)) " *
                  "($(s.n_reversals) reversals, $(s.n_aborted) aborted, $(s.n_excursions) excursions; grouping rule: $(s.n_reversals_grouped) reversals, " *
                  "$(s.n_aborted_grouped) aborted); variance $(round(s.variance, digits = 3)) over $(s.variance_intervals) intervals; dv $(round(s.dv, digits = 3)); " *
                  "D₁/ℱ = $(round(s.d1_fit.c0, digits = 2)) + $(round(s.d1_fit.c2, digits = 2)) a²; c² = $(round(s.c2, digits = 3)), γ = $(round(s.γ, digits = 2))"
            (sim = sim, obs = obs, summary = s)
        end
        pooled = pool!(runs[1].obs, (r.obs for r in runs[2:end]))
        s = summarize(pooled, runs[1].sim)
        @info "T4 pooled ($(length(T4_SEEDS)) × $(T4_LENGTH) τ_d): ⟨T_r⟩ = $(round(s.mean_reversal_time, digits = 1)) ± $(round(s.mean_reversal_time_error, digits = 1)) " *
              "($(s.n_reversals) reversals; paper 107 ± 5 from 468); variance $(round(s.variance, digits = 3)) (paper 0.16); dv $(round(s.dv, digits = 3)); " *
              "D₁/ℱ = $(round(s.d1_fit.c0, digits = 3)) + $(round(s.d1_fit.c2, digits = 3)) a² (paper 1.36 + 2.00 a²; count-weighted fit $(round(s.d1_fit_weighted.c0, digits = 3)) + $(round(s.d1_fit_weighted.c2, digits = 3)) a²; mean over |a| ≤ 0.2: $(round(s.d1_zero, digits = 3))); c² = $(round(s.c2, digits = 3)) (0.246), γ = $(round(s.γ, digits = 3)) (2.96); " *
              "theory for these (c², γ): variance $(round(s.theory.variance, digits = 3)), ⟨T_r⟩ Eq. 58 $(round(s.theory.mean_reversal_time, digits = 1)), Eq. 57 $(round(s.theory.mean_reversal_time_exact, digits = 1)), dv $(round(s.theory.dv, digits = 3))"
        # T4.1: mean time between reversals, 2σ combined with the paper's ±5
        @test abs(s.mean_reversal_time - 107) < 2 * sqrt(5^2 + s.mean_reversal_time_error^2)
        # T4.2: between-reversal relative variance
        @test abs(s.variance - 0.16) < 0.02
        # T4.3: the D₁/ℱ fit for |a| ≤ 1.1
        @test relerr(s.d1_fit.c0, 1.36) < 0.1
        @test relerr(s.d1_fit.c2, 2.00) < 0.1
        # T4.4: derived shape parameters
        @test relerr(s.c2, 0.246) < 0.1
        @test relerr(s.γ, 2.96) < 0.1
        # seed-to-seed scatter of ⟨T_r⟩ (feeds T4.18) and the grouping rule against the trigger
        Trs = [r.summary.mean_reversal_time for r in runs]
        @info "T4 seed scatter: ⟨T_r⟩ = $(round.(Trs, digits = 1)); grouping-rule counts $(join(string.(r.summary.n_reversals_grouped for r in runs), ", ")) vs trigger $(join(string.(r.summary.n_reversals for r in runs), ", "))"
        @test all(r.summary.n_reversals ≥ 300 for r in runs)        # each seed alone has decent statistics
        # T4.13 (diagnostic, kind D): the parity histogram over the whole run, sampled every n_c
        # steps, events included; the event-window exclusion is build step 7
        pv = value(pooled.parity)
        cum = cumsum(pv.counts)
        median_P = pv.centers[findfirst(≥(cum[end] / 2), cum)]
        @info "T4.13 parity over the whole run (events included): median 𝒫 ≈ $(round(median_P, digits = 3)), fraction below −0.9: $(round(sum(pv.counts[pv.centers .< -0.9]) / cum[end], digits = 3))"
        @test -1 ≤ median_P < 0
    end
end
