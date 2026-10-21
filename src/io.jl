# HDF5 run files (dev/DECISIONS.md D7) and checkpoints (D6).
#
# Layout of a run file:
#   /                 attributes: every RunSpec field, package and Julia versions, git commit, date
#   /grid/theta       the interior colatitudes
#   /streams/<name>   extendible, chunked datasets written from fixed buffers (attributes
#                     `every` in steps and `dt` in τ_d); `T`, `P` are M × nsnapshots
#   /results/<name>   the `value` of each probe's aggregator, one group per probe
#   /results/summary  the NamedTuple of `summarize`
#   /checkpoint       `blob` (serialised state, RNG and aggregators), attributes `n`, `time`

"""
    RunFile(path, spec; overwrite = false)
    open_run(path; mode = "r")

An HDF5 run file. `RunFile` creates it and writes the metadata of `spec`; `open_run` reopens
an existing one (`mode = "r+"` to resume into it). Close with `close`.
"""
mutable struct RunFile
    const path::String
    h5::HDF5.File
end

function RunFile(path::AbstractString, spec::RunSpec; overwrite::Bool = false)
    isfile(path) && !overwrite && throw(ArgumentError("$path exists; pass overwrite = true to replace it"))
    mkpath(dirname(abspath(path)))
    h5 = h5open(path, "w")
    file = RunFile(String(path), h5)
    write_metadata!(file, spec)
    create_group(h5, "streams")
    create_group(h5, "results")
    return file
end

open_run(path::AbstractString; mode::AbstractString = "r") = RunFile(String(path), h5open(path, mode))
Base.close(f::RunFile) = close(f.h5)
Base.isopen(f::RunFile) = isopen(f.h5)
Base.flush(f::RunFile) = flush(f.h5)
Base.show(io::IO, f::RunFile) = print(io, "RunFile(\"", f.path, "\"", isopen(f) ? "" : ", closed", ")")

const RUN_FILE_LAYOUT = 1

function git_commit()
    try
        return readchomp(`git -C $(pkgdir(HoyngDynamo)) rev-parse HEAD`)
    catch
        return "unknown"
    end
end

function write_metadata!(file::RunFile, spec::RunSpec{T}) where {T}
    a = attributes(file.h5)
    a["layout_version"] = RUN_FILE_LAYOUT
    a["package_version"] = string(pkgversion(HoyngDynamo))
    a["julia_version"] = string(VERSION)
    a["git_commit"] = git_commit()
    a["created"] = string(Dates.now())
    a["float_type"] = string(T)
    a["kR"] = spec.kR; a["Cm"] = spec.Cm; a["N"] = spec.N; a["dtau"] = spec.Δτ
    a["tau_c"] = spec.τc; a["Nc"] = spec.Nc; a["seed"] = spec.seed
    a["stepper"] = String(spec.stepper); a["a0"] = spec.a0
    a["distribution"] = string(nameof(typeof(spec.distribution)))
    a["schedule"] = string(nameof(typeof(spec.schedule)))
    for name in fieldnames(typeof(spec.schedule))
        a["schedule_" * String(name)] = getfield(spec.schedule, name)
    end
    g = Grid(spec.N; T)
    create_group(file.h5, "grid")["theta"] = g.θ
    return file
end

"All root attributes of a run file as a Dict."
metadata(file::RunFile) = Dict{String,Any}(k => read(attributes(file.h5)[k]) for k in keys(attributes(file.h5)))

const DISTRIBUTIONS = Dict("GaussianF" => GaussianF(), "UniformF" => UniformF(), "BinaryF" => BinaryF())

"""
    RunSpec(file::RunFile)

The run specification stored in a run file's metadata (for resuming or re-running).
"""
function RunSpec(file::RunFile)
    m = metadata(file)
    T = m["float_type"] == "Float32" ? Float32 : Float64
    sched = if m["schedule"] == "ConstantForcing"
        ConstantForcing(T(m["schedule_f"]))
    elseif m["schedule"] == "LinearRamp"
        LinearRamp(T(m["schedule_f0"]), T(m["schedule_f1"]), T(m["schedule_τ0"]), T(m["schedule_τ1"]))
    else
        error("unknown forcing schedule $(m["schedule"])")
    end
    RunSpec(; kR = m["kR"], Cm = m["Cm"], N = m["N"], Δτ = m["dtau"], schedule = sched, τc = m["tau_c"], Nc = m["Nc"],
            distribution = DISTRIBUTIONS[m["distribution"]], seed = m["seed"], stepper = Symbol(m["stepper"]), a0 = m["a0"], T)
end

# ---------------------------------------------------------------------------------------
# Stream recorders
# ---------------------------------------------------------------------------------------

function stream_dataset(file::RunFile, name::AbstractString, eltype::Type, dims::Tuple, maxdims::Tuple, chunk::Tuple,
                        every::Integer, dt::Real, attach::Bool)
    g = file.h5["streams"]
    if haskey(g, name)
        attach || throw(ArgumentError("stream \"$name\" already exists in $(file.path)"))
        d = g[name]
        read(attributes(d)["every"]) == every || throw(ArgumentError("stream \"$name\" was recorded every $(read(attributes(d)["every"])) steps, not $every"))
        return d
    end
    attach && throw(ArgumentError("stream \"$name\" not found in $(file.path)"))
    d = create_dataset(g, name, eltype, (dims, maxdims); chunk = chunk)
    attributes(d)["every"] = Int(every)
    attributes(d)["dt"] = Float64(dt)
    return d
end

"""
    Recorder(file, name, f, every; sim, eltype = Float64, chunk = 10_000, attach = false)

An observer that appends `f(sim)` every `every` steps to the stream `name` of `file`,
through a fixed buffer of `chunk` samples. `flush!` writes a partial buffer.
"""
mutable struct Recorder{T,F} <: Observer
    const f::F
    const every::Int
    const buf::Vector{T}
    const dset::HDF5.Dataset
    fill::Int
    written::Int
end

function Recorder(file::RunFile, name::AbstractString, f, every::Integer; sim::Simulation, eltype::Type = Float64,
                  chunk::Integer = 10_000, attach::Bool = false)
    every ≥ 1 || throw(ArgumentError("need every ≥ 1"))
    d = stream_dataset(file, name, eltype, (0,), (-1,), (Int(chunk),), every, every * timestep(sim), attach)
    Recorder{eltype,typeof(f)}(f, Int(every), Vector{eltype}(undef, chunk), d, 0, 0)
end

@inline function observe!(r::Recorder, sim::Simulation)
    if sim.state.n % r.every == 0
        r.fill += 1
        @inbounds r.buf[r.fill] = r.f(sim)
        r.fill == length(r.buf) && flush!(r)
    end
    return r
end

function flush!(r::Recorder)
    r.fill == 0 && return nothing
    n = r.written + r.fill
    HDF5.set_extent_dims(r.dset, (n,))
    r.dset[r.written+1:n] = r.buf[1:r.fill]
    r.written = n
    r.fill = 0
    return nothing
end

"""
    SnapshotRecorder(file, name, f, every; sim, eltype = Float64, chunk = 100, attach = false)

Like `Recorder` for a vector-valued `f(sim)` (a field at the nodes): the stream is an
M × nsnapshots dataset.
"""
mutable struct SnapshotRecorder{T,F} <: Observer
    const f::F
    const every::Int
    const buf::Matrix{T}
    const dset::HDF5.Dataset
    fill::Int
    written::Int
end

function SnapshotRecorder(file::RunFile, name::AbstractString, f, every::Integer; sim::Simulation, eltype::Type = Float64,
                          chunk::Integer = 100, attach::Bool = false)
    every ≥ 1 || throw(ArgumentError("need every ≥ 1"))
    M = length(f(sim))
    d = stream_dataset(file, name, eltype, (M, 0), (M, -1), (M, Int(chunk)), every, every * timestep(sim), attach)
    SnapshotRecorder{eltype,typeof(f)}(f, Int(every), Matrix{eltype}(undef, M, chunk), d, 0, 0)
end

@inline function observe!(r::SnapshotRecorder, sim::Simulation)
    if sim.state.n % r.every == 0
        r.fill += 1
        v = r.f(sim)
        @inbounds for j in eachindex(v)
            r.buf[j, r.fill] = v[j]
        end
        r.fill == size(r.buf, 2) && flush!(r)
    end
    return r
end

function flush!(r::SnapshotRecorder)
    r.fill == 0 && return nothing
    n = r.written + r.fill
    HDF5.set_extent_dims(r.dset, (size(r.buf, 1), n))
    r.dset[:, r.written+1:n] = r.buf[:, 1:r.fill]
    r.written = n
    r.fill = 0
    return nothing
end

"""
    read_stream(file, name) -> (values, every, dt)

A stream and its cadence (steps and τ_d).
"""
function read_stream(file::RunFile, name::AbstractString)
    d = file.h5["streams"][name]
    return (values = read(d), every = read(attributes(d)["every"]), dt = read(attributes(d)["dt"]))
end

"The names of the streams in a run file."
streams(file::RunFile) = keys(file.h5["streams"])

# ---------------------------------------------------------------------------------------
# Results: aggregator values as groups of datasets
# ---------------------------------------------------------------------------------------

write_value!(g, k::AbstractString, v::Real) = (g[k] = v; nothing)
write_value!(g, k::AbstractString, v::Complex) = (g[k] = v; nothing)
function write_value!(g, k::AbstractString, v::AbstractVector{<:Union{Real,Complex}})
    isempty(v) ? create_dataset(g, k, eltype(v), (0,)) : (g[k] = collect(v))
    return nothing
end
write_string_vector!(g, k::AbstractString, v::Vector{String}) = (isempty(v) ? create_dataset(g, k, String, (0,)) : (g[k] = v); nothing)
write_value!(g, k::AbstractString, v::Symbol) = (g[k] = String(v); nothing)
write_value!(g, k::AbstractString, v::AbstractString) = (g[k] = String(v); nothing)
write_value!(g, k::AbstractString, v::NamedTuple) = (write_value!(create_group(g, k), v); nothing)

function write_value!(g, k::AbstractString, events::AbstractVector{<:Event})
    e = create_group(g, k)
    write_string_vector!(e, "kind", [String(ev.kind) for ev in events])
    write_value!(e, "start", [ev.start for ev in events])
    write_value!(e, "stop", [ev.stop for ev in events])
    write_value!(e, "amin", [ev.amin for ev in events])
    write_value!(e, "crossings", [ev.crossings for ev in events])
    write_value!(e, "tcross", [ev.tcross for ev in events])
    return nothing
end

function write_value!(g, nt::NamedTuple)
    for (k, v) in pairs(nt)
        write_value!(g, String(k), v)
    end
    return nothing
end

"""
    write_results!(file, observers)

Write `value` of every probe in the NamedTuple `observers` to `/results/<name>`, replacing
earlier results of the same name.
"""
function write_results!(file::RunFile, observers::NamedTuple)
    g = file.h5["results"]
    for (k, o) in pairs(observers)
        o isa Probe || continue
        name = String(k)
        haskey(g, name) && delete_object(g, name)
        sub = create_group(g, name)
        attributes(sub)["aggregator"] = string(nameof(typeof(o.agg)))
        attributes(sub)["every"] = o.every
        write_value!(sub, value(o))
    end
    flush(file.h5)
    return file
end

"Write a `summarize` NamedTuple to `/results/summary`."
function write_summary!(file::RunFile, summary::NamedTuple)
    g = file.h5["results"]
    haskey(g, "summary") && delete_object(g, "summary")
    write_value!(create_group(g, "summary"), summary)
    flush(file.h5)
    return file
end

read_results(file::RunFile) = read(file.h5["results"])

# ---------------------------------------------------------------------------------------
# Checkpoints
# ---------------------------------------------------------------------------------------

aggregator_states(p::Probe) = p.agg
aggregator_states(::Observer) = nothing
aggregator_states(t::Tuple) = map(aggregator_states, t)
aggregator_states(nt::NamedTuple) = map(aggregator_states, nt)
aggregator_states(::Nothing) = nothing

"""
    checkpoint!(file, sim, observers)

Flush the recorders and write the state (b, a, n), the noise RNG and every probe's
aggregator to `/checkpoint` so that the run can be resumed bitwise with `restore!`. Must
fall on a redraw step (n a multiple of n_c): the current F is not stored (D6).
"""
function checkpoint!(file::RunFile, sim::Simulation, observers)
    nc = redraw_interval(sim.noise)
    sim.state.n % nc == 0 || throw(ArgumentError("checkpoints must fall on a redraw step (n = $(sim.state.n) is not a multiple of n_c = $nc)"))
    flush!(observers)
    payload = (n = sim.state.n, a = sim.state.a, b = copy(sim.state.b), rng = copy(sim.noise.rng),
               aggregators = aggregator_states(observers))
    io = IOBuffer()
    serialize(io, payload)
    bytes = take!(io)
    h5 = file.h5
    g = haskey(h5, "checkpoint") ? h5["checkpoint"] : create_group(h5, "checkpoint")
    if !haskey(g, "blob")
        create_dataset(g, "blob", UInt8, ((0,), (-1,)); chunk = (1 << 16,))
    end
    d = g["blob"]
    HDF5.set_extent_dims(d, (length(bytes),))
    d[1:length(bytes)] = bytes
    attrs(g)["n"] = sim.state.n                 # attrs overwrites; attributes would refuse the second time
    attrs(g)["time"] = current_time(sim)
    flush(h5)
    return file
end

# Copy the state of `src` into `dst` field by field: const fields are arrays (copied into)
# or nested aggregators (recursed into); the rest are assigned.
function restore!(dst::A, src::A) where {A}
    for i in 1:fieldcount(A)
        v = getfield(src, i)
        if Base.isconst(A, i)
            restore_field!(getfield(dst, i), v, fieldname(A, i))
        else
            setfield!(dst, i, v)
        end
    end
    return dst
end
restore!(dst, src) = throw(ArgumentError("cannot restore a $(typeof(dst)) from a checkpointed $(typeof(src))"))

restore_field!(d::AbstractVector, v::AbstractVector, _) = (length(d) == length(v) || resize!(d, length(v)); copyto!(d, v); nothing)
restore_field!(d::Aggregator, v::Aggregator, _) = (restore!(d, v); nothing)
restore_field!(d, v, name) = (isequal(d, v) || throw(ArgumentError("field $name differs between the run ($d) and the checkpoint ($v)")); nothing)

restore_observer!(p::Probe, agg, n) = (restore!(p.agg, agg); nothing)
function restore_observer!(r::Recorder, ::Nothing, n)
    r.fill = 0
    r.written = n ÷ r.every
    HDF5.set_extent_dims(r.dset, (r.written,))
    return nothing
end
function restore_observer!(r::SnapshotRecorder, ::Nothing, n)
    r.fill = 0
    r.written = n ÷ r.every
    HDF5.set_extent_dims(r.dset, (size(r.buf, 1), r.written))
    return nothing
end
restore_observers!(obs::Tuple, saved::Tuple, n) = (foreach((o, s) -> restore_observer!(o, s, n), obs, saved); nothing)
restore_observers!(obs::NamedTuple, saved::NamedTuple, n) = restore_observers!(values(obs), values(saved), n)
restore_observers!(::Nothing, ::Nothing, n) = nothing

"""
    restore!(sim, observers, file)

Resume from the checkpoint in `file`: the state, the noise RNG and the probes' aggregators
are restored, and the recorders (created with `attach = true`) are truncated to the
checkpointed length. `sim` and `observers` must have been built from the same `RunSpec`
and observer set as the checkpointed run; the file must be open for writing (`mode = "r+"`).
"""
function restore!(sim::Simulation, observers, file::RunFile)
    haskey(file.h5, "checkpoint") || throw(ArgumentError("$(file.path) has no checkpoint"))
    g = file.h5["checkpoint"]
    payload = deserialize(IOBuffer(read(g["blob"])))
    nc = redraw_interval(sim.noise)
    payload.n % nc == 0 || error("checkpoint at n = $(payload.n) is not a redraw step")
    length(payload.b) == length(sim.state.b) || throw(ArgumentError("checkpoint has $(length(payload.b)) state entries, the simulation $(length(sim.state.b))"))
    copyto!(sim.state.b, payload.b)
    sim.state.a = payload.a
    sim.state.n = payload.n
    copy!(sim.noise.rng, payload.rng)
    restore_observers!(observers, payload.aggregators, payload.n)
    return sim
end

"The step count of the checkpoint in a run file, or `nothing`."
checkpoint_step(file::RunFile) = haskey(file.h5, "checkpoint") ? read(attributes(file.h5["checkpoint"])["n"]) : nothing
