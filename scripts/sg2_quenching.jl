# SG-2 (dev/STRETCH_GOALS.md), gate T4.16: do the alternative quenching laws leave the dipole
# statistics "not materially changed" (HOS p. 277)? The four laws of `QUENCHINGS` on the
# standard run at the given forcings, three seeds × 5×10⁴ τ_d pooled (`ensemble`), preceded by
# each law's deterministic equilibrium (relaxation rate about a = 1, overtone content, parity):
#
#     julia --project -t 4 scripts/sg2_quenching.jl [f ...]        (default 6.4 5.4 7.4)
#
# Prints the tables and appends them to results/sg2_quenching.md (SG2_OUT overrides; SG2_TAU
# shortens the runs for a smoke test). About 2 min per forcing on 4 threads.

using HoyngDynamo, Printf, LinearAlgebra, Dates

forcings = isempty(ARGS) ? [6.4, 5.4, 7.4] : parse.(Float64, ARGS)
seeds = 1:3
τ_end = parse(Float64, get(ENV, "SG2_TAU", "5e4"))
Δτ = 0.01
out = get(ENV, "SG2_OUT", normpath(joinpath(@__DIR__, "..", "results", "sg2_quenching.md")))
mkpath(dirname(out))
const OUT = open(out, "a")
emit(s::AbstractString = "") = (println(s); println(OUT, s); flush(OUT))
unwrap(e) = e isa TaskFailedException ? unwrap(e.task.result) : (e isa CompositeException && !isempty(e) ? unwrap(first(e)) : e)

emit()
emit("## SG-2 alternative quenching laws, gate T4.16 ($(Dates.format(now(), "yyyy-mm-dd HH:MM")), $(Threads.nthreads()) threads, $(length(seeds)) seeds × $(τ_end) τ_d)")
emit()

# ---- deterministic equilibria ---------------------------------------------------------------
model = Model()
g, cst = model.grid, model.constants
M = g.M
loc = calibrate(LocalQuenching, model, Δτ)
laws = ["amplitude" => AmplitudeQuenching(), "energy" => EnergyQuenching(), "inverse" => InverseQuenching(), "local" => loc]

function settle(law; b = model.b0, τ = 200.0)
    sp = ExponentialStepper(model, Δτ; quenching = law)
    st, cache = State(model, b), Cache(model)
    for _ in 1:round(Int, τ / Δτ)
        step!(st, sp, NoNoise(), cache)
    end
    return st
end

function relaxation_rate(law; δ = 1e-6)
    sp = ExponentialStepper(model, Δτ; quenching = law)
    st, cache = State(model, (1 + δ) .* settle(law).b), Cache(model)
    x = Float64[]
    for _ in 1:300
        step!(st, sp, NoNoise(), cache)
        push!(x, st.a - 1)
    end
    return -log(abs(x[300]) / abs(x[100])) / 2.0
end

emit("Local law constant κ = $(round(loc.κ, digits = 4)) for an equilibrium at a = 1 (the sin θ-mean of q over b₀ would need 2q₀ = $(round(2cst.q0, digits = 4))). 2μ = $(round(2cst.μ, digits = 4)), 2μ/(1 + s) = $(round(2cst.μ / (1 + cst.s), digits = 4)).")
emit()
emit("| law | coefficient c | relaxation rate about a = 1 | equilibrium overtone residual ‖b − a b₀‖ | equilibrium toroidal energy | parity |")
emit("|---|---|---|---|---|---|")
formulas = Dict("amplitude" => "s(1 − a²)", "energy" => "s(1 − Σ w T²)", "inverse" => "s(1 − a²)/(1 + s a²)", "local" => "s − (1 + s) κ T(θ)²")
for (name, law) in laws
    st = settle(law)
    ET = sum(g.w .* st.b[M+1:2M] .^ 2)
    emit(@sprintf("| %s | %s | %.4f | %.3e | %.4f | %.4f |", name, formulas[name], relaxation_rate(law), overtone_residual(model, st), ET,
                  HoyngDynamo.toroidal_parity(g, st.b)))
end
emit()

# ---- runs ----------------------------------------------------------------------------------
for f in forcings
    emit("### f = $f")
    emit()
    emit("| law | ⟨T_r⟩ | reversals | aborted | §2 var | dv | weighted D₁ fit | unweighted fit | D₁(0) direct | D₁ bin 0.95–1.0 | (c², γ) | s/seed |")
    emit("|---|---|---|---|---|---|---|---|---|---|---|---|")
    rows = Dict{Symbol,Any}()
    for q in HoyngDynamo.QUENCHINGS
        e = try
            ensemble(RunSpec(; f, quenching = q), τ_end; seeds)
        catch err
            # a law may run away (the energy law reverses the sign of α under a strong overtone
            # burst, see VALIDATION.md T4.16); record it and go on
            emit("| $q | diverged ($(first(split(sprint(showerror, unwrap(err)), '\n')))) | | | | | | | | | | |")
            continue
        end
        s = e.summary
        d = value(e.pooled.d1)
        i = findfirst(c -> abs(c - 0.975) < 1e-9, d.centers)
        rows[q] = (s = s, d1_at_1 = d.mean[i])
        emit(@sprintf("| %s | %.1f ± %.1f | %d | %d | %.3f | %.3f | %.2f + %.2f a² | %.2f + %.2f a² | %.2f | %.2f | (%.3f, %.2f) | %.0f |",
                      q, s.mean_reversal_time, s.mean_reversal_time_error, s.n_reversals, s.n_aborted, s.variance, s.dv,
                      s.d1_fit.c0, s.d1_fit.c2, s.d1_fit_unweighted.c0, s.d1_fit_unweighted.c2, s.d1_zero, d.mean[i], s.c2, s.γ,
                      sum(r.seconds for r in e.runs) / length(seeds)))
    end
    emit()
    ref = rows[:amplitude]
    emit("Change from the amplitude law at f = $f (gate T4.16: ⟨T_r⟩ and variance within ~20%):")
    emit()
    emit("| law | Δ⟨T_r⟩ (%) | Δ §2 var (%) | Δ dv (%) | Δ weighted intercept | Δ c² | Δ γ |")
    emit("|---|---|---|---|---|---|---|")
    for q in HoyngDynamo.QUENCHINGS
        (q === :amplitude || !haskey(rows, q)) && continue
        r = rows[q]
        emit(@sprintf("| %s | %+.1f | %+.1f | %+.1f | %+.2f | %+.3f | %+.2f |", q, 100 * (r.s.mean_reversal_time / ref.s.mean_reversal_time - 1),
                      100 * (r.s.variance / ref.s.variance - 1), 100 * (r.s.dv / ref.s.dv - 1), r.s.d1_fit.c0 - ref.s.d1_fit.c0,
                      r.s.c2 - ref.s.c2, r.s.γ - ref.s.γ))
    end
    emit()
end
close(OUT)
println("appended to $out")
