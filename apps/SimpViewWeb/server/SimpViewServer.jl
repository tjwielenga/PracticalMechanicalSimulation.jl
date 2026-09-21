"""Local model-preview bridge for the browser-based SimpView interface."""
module SimpViewServer

using HTTP
using JSON
using TOML
using PracticalMechanicalSimulation

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(REPOSITORY_ROOT, "src", "viewer", "ViewerData.jl"))
include(joinpath(REPOSITORY_ROOT, "src", "viewer", "StoredResultViewer.jl"))
include(joinpath(REPOSITORY_ROOT, "src", "viewer",
    "PortableViewerDocument.jl"))

using .StoredResultViewer: model_mechanism_result
using .StoredResultViewer: simulation_mechanism_result
using .PortableViewerDocument: viewer_document, write_graphics

export preview_model_document, start_model_run, model_run_snapshot,
       serve_simpview

const MAXIMUM_MODEL_BYTES = 20 * 1024 * 1024
const RUNS = Dict{String,Any}()
const RUNS_LOCK = ReentrantLock()
const RUN_NUMBER = Ref(0)
const MODEL_CACHE = Dict{Tuple{String,UInt,Int},Any}()
const MODEL_CACHE_LIMIT = 8
allowed_origin(origin) = occursin(
    r"^http://(?:127\.0\.0\.1|localhost):\d+$", origin)

function model_dimension(source)
    document = TOML.parse(source)
    model = get(document, "model", nothing)
    model isa AbstractDict || throw(ArgumentError(
        "model preview requires a [model] table"))
    dimension = Symbol(get(model, "dimension", ""))
    dimension in (:planar, :spatial) || throw(ArgumentError(
        "model.dimension must be 'planar' or 'spatial'"))
    dimension
end

function load_uploaded_model(source, name)
    extension = lowercase(splitext(name)[2])
    if extension == ".lua"
        return load_spatial_model(IOBuffer(source); format = :lua,
            source_directory = REPOSITORY_ROOT, source_label = name)
    elseif extension == ".toml"
        dimension = model_dimension(source)
        return dimension == :spatial ?
            load_spatial_model(IOBuffer(source); format = :toml,
                source_directory = REPOSITORY_ROOT, source_label = name) :
            load_planar_model(IOBuffer(source);
                source_directory = REPOSITORY_ROOT)
    end
    throw(ArgumentError("model preview requires a .toml or .lua file"))
end

model_cache_key(source, name) =
    (String(name), hash(source), ncodeunits(source))

function cache_model!(source, name, loaded)
    lock(RUNS_LOCK) do
        MODEL_CACHE[model_cache_key(source, name)] = loaded
        while length(MODEL_CACHE) > MODEL_CACHE_LIMIT
            delete!(MODEL_CACHE, first(keys(MODEL_CACHE)))
        end
    end
    loaded
end

function cached_model(source, name)
    lock(RUNS_LOCK) do
        get(MODEL_CACHE, model_cache_key(source, name), nothing)
    end
end

"""Build entered and consistent-IC viewer choices from uploaded model text."""
function preview_model_document(source::AbstractString,
        name::AbstractString = "model.toml")
    ncodeunits(source) <= MAXIMUM_MODEL_BYTES || throw(ArgumentError(
        "model source is larger than the 20 MB preview limit"))
    model_name = basename(name)
    loaded = cached_model(source, model_name)
    if isnothing(loaded)
        loaded = cache_model!(source, model_name,
            load_uploaded_model(source, model_name))
    end
    results = [
        model_mechanism_result(loaded; configuration = :entered),
        model_mechanism_result(loaded; configuration = :consistent),
    ]
    document = viewer_document(results;
        labels = ["Model input", "Consistent initial conditions"],
        choice_name = "Configuration", include_signals = false)
    document["title"] = loaded.title
    start = loaded.simulation.start_time
    finish = loaded.simulation.end_time
    document["runnable"] = true
    document["run_settings"] = Dict(
        "start_time" => start,
        "end_time" => finish,
        "frames_per_second" => 60.0,
        "analysis_mode" => String(requested_analysis_mode(loaded)),
    )
    document
end

mutable struct ModelRun
    id::String
    result_name::String
    result_path::String
    loaded::Any
    viewer_loaded::Any
    analysis_mode::Symbol
    starting_result::Any
    status::Symbol
    message::String
    revision::Int
    times::Vector{Float64}
    states::Vector{Vector{Float64}}
    static_events::Vector{Any}
    dynamic_started::Bool
    result::Any
    writer::Any
    task::Union{Nothing,Task}
    lock::ReentrantLock
end

function next_run_id()
    lock(RUNS_LOCK) do
        RUN_NUMBER[] += 1
        string(time_ns(), '-', RUN_NUMBER[])
    end
end

function requested_analysis_mode(loaded)
    mode = loaded.analysis.mode
    if mode == :automatic
        return loaded isa LoadedPlanarModel &&
            loaded.analysis.degrees_of_freedom == 0 ? :kinematic : :dynamic
    end
    mode
end

function loaded_for_analysis(loaded, mode::Symbol)
    mode == :initial_conditions && return loaded
    actual_mode = mode == :dynamic && loaded isa LoadedPlanarModel &&
        loaded.analysis.degrees_of_freedom == 0 ? :kinematic : mode
    initialization = mode in (:dynamic, :static) ? :none :
        loaded.analysis.initialization
    settings = merge(loaded.analysis,
        (; mode = actual_mode, initialization))
    fields = fieldnames(typeof(loaded))
    values = ntuple(length(fields)) do index
        fields[index] == :analysis ? settings : getfield(loaded, fields[index])
    end
    typeof(loaded)(values...)
end

function completed_starting_result(id, loaded)
    isnothing(id) && return nothing
    prior = find_run(String(id))
    isnothing(prior) && throw(ArgumentError(
        "the selected starting result is no longer available"))
    lock(prior.lock) do
        prior.status == :complete || throw(ArgumentError(
            "the selected starting analysis is not complete"))
        !isnothing(prior.result) &&
            prior.result.analysis_mode in (:static, :dynamic, :kinematic) ||
            throw(ArgumentError(
                "the selected result cannot start the requested analysis"))
        prior.loaded.model_source == loaded.model_source || throw(ArgumentError(
            "the starting result belongs to a different model"))
        continuing = prior.result.analysis_mode in (:dynamic, :kinematic)
        (; state = copy(last(prior.result.states)),
           time = Float64(last(prior.result.times)),
           restore_declared_velocities = !continuing,
           times = continuing ? Float64.(prior.result.times) : Float64[],
           states = continuing ? copy.(prior.result.states) :
               Vector{Vector{Float64}}())
    end
end

function initial_consistency_count(loaded)
    conditions = loaded.initial_conditions
    sum(Int(getproperty(conditions, name)) for name in
        (:position_corrections, :velocity_corrections,
         :acceleration_corrections) if hasproperty(conditions, name))
end

function initial_conditions_result(run::ModelRun)
    time = run.loaded.simulation.start_time
    record_run_event!(run, (;
        kind = :begin, loaded = run.loaded,
        analysis_mode = :initial_conditions, times = [time]))
    state = copy(run.loaded.initial_values)
    record_run_event!(run, (; kind = :sample, time, state))
    (; times = [time], states = [state], loaded = run.loaded,
       analysis_mode = :initial_conditions,
       initial_consistency_iterations = initial_consistency_count(run.loaded))
end

function dynamic_result_from_previous(run::ModelRun, end_time, samples)
    loaded = run.loaded
    start = run.starting_result.time
    times = samples == 1 ? [start] :
        collect(range(start, end_time; length = samples))
    state = copy(run.starting_result.state)
    record_run_event!(run, (;
        kind = :begin, loaded, analysis_mode = run.analysis_mode, times,
        initial_state = state))
    continued = if loaded isa LoadedSpatialModel
        runner = PracticalMechanicalSimulation.SpatialSimulationRunner
        runner.set_spatial_analysis_stage!(loaded, :dynamic;
            state, time = start)
        runner.prepare_spatial_dynamic_state!(state, loaded, start;
            restore_declared_velocities =
                run.starting_result.restore_declared_velocities)
        result = runner.run_spatial_implicit_model(loaded, times;
            initial_state = state,
            sample_progress = event -> record_run_event!(run, event))
        merge(result, (; static_progress_events = Any[]))
    else
        if run.starting_result.restore_declared_velocities
            for variable in loaded.layout.catalog.variables
                variable.kind in
                    (:velocity, :angular_velocity, :relative_velocity) ||
                    continue
                state[variable.index] = loaded.initial_values[variable.index]
            end
        end
        PracticalMechanicalSimulation.PlanarModelIO.correct_initial_velocities!(
            state, loaded.model, loaded.initial_variable_weights,
            Set(loaded.imposed_initial_variable_indices), start)
        PracticalMechanicalSimulation.SimulationRunner.run_implicit_model(
            loaded, times, run.analysis_mode; initial_state = state,
            sample_progress = event -> record_run_event!(run, event))
    end
    isempty(run.starting_result.times) && return continued
    merge(continued, (;
        times = [run.starting_result.times; continued.times[2:end]],
        states = [run.starting_result.states; continued.states[2:end]]))
end

function modal_result_from_previous(run::ModelRun)
    loaded = run.loaded
    time = run.starting_result.time
    state = copy(run.starting_result.state)
    record_run_event!(run, (;
        kind = :begin, loaded, analysis_mode = :modal, times = [time],
        initial_state = state))
    if loaded isa LoadedSpatialModel
        runner = PracticalMechanicalSimulation.SpatialSimulationRunner
        return runner.run_spatial_modal_model(
            loaded, time; initial_state = state)
    end
    PracticalMechanicalSimulation.SimulationRunner.run_modal_model(
        loaded, time; initial_state = state)
end

function result_file_name(model_name)
    stem = first(splitext(basename(model_name)))
    isempty(stem) && (stem = "simpview-result")
    stem * ".simp"
end

function record_run_event!(run::ModelRun, event)
    lock(run.lock) do
        if event.kind == :begin
            run.analysis_mode = event.analysis_mode
            seeded = !isempty(run.times)
            initial_time = seeded ? first(run.times) : first(event.times)
            initial_state = seeded ? first(run.states) :
                (hasproperty(event, :initial_state) ? event.initial_state :
                    event.loaded.initial_values)
            run.writer = begin_incremental_result(run.result_path,
                event.loaded, event.analysis_mode; overwrite = true,
                initial_time, initial_state)
            if seeded
                for index in 2:length(run.times)
                    record_incremental_event!(run.writer, (;
                        kind = :sample, time = run.times[index],
                        state = run.states[index]))
                end
            end
        else
            isnothing(run.writer) || record_incremental_event!(run.writer, event)
            if event.kind in (:sample, :failure)
                time = Float64(event.time)
                state = Float64.(event.state)
                if !isempty(run.times) && run.times[end] == time
                    run.times[end] = time
                    run.states[end] = state
                else
                    push!(run.times, time)
                    push!(run.states, state)
                end
                run.dynamic_started = true
            elseif event.kind == :static
                push!(run.static_events, event)
            end
        end
        run.revision += 1
    end
    yield()
    nothing
end

function run_mechanisms(run::ModelRun)
    snapshot = lock(run.lock) do
        (; result = run.result, loaded = run.viewer_loaded,
           analysis_mode = run.analysis_mode,
           times = copy(run.times), states = copy.(run.states),
           static_events = copy(run.static_events),
           dynamic_started = run.dynamic_started,
           status = run.status, message = run.message,
           revision = run.revision, result_name = run.result_name)
    end
    results, labels, choice_name = if !isnothing(snapshot.result)
        result = snapshot.result
        if result.analysis_mode == :modal
            count = length(result.eigenvalues)
            mechanisms = [simulation_mechanism_result(result; mode = index)
                for index in 1:count]
            names = ["Mode $index — " *
                string(round(result.natural_frequencies_hz[index];
                    sigdigits = 6)) * " Hz" for index in 1:count]
            mechanisms, names, "Mode"
        elseif result.analysis_mode in (:dynamic, :kinematic) &&
                hasproperty(result, :static_progress_events) &&
                !isempty(result.static_progress_events)
            [simulation_mechanism_result(result;
                 analysis_mode = :static),
             simulation_mechanism_result(result)],
                ["Static initialization", "Dynamic"], "Analysis"
        else
            [simulation_mechanism_result(result)], ["Result"], "Result"
        end
    elseif !snapshot.dynamic_started && !isempty(snapshot.static_events)
        partial = (; loaded = snapshot.loaded, analysis_mode = :static,
            times = [snapshot.loaded.simulation.start_time],
            states = [copy(last(snapshot.static_events).state)],
            static_progress_events = snapshot.static_events)
        [simulation_mechanism_result(partial)], ["Static progress"], "Result"
    elseif !isempty(snapshot.states)
        partial = (; loaded = snapshot.loaded,
            analysis_mode = snapshot.analysis_mode,
            times = snapshot.times, states = snapshot.states,
            static_progress_events = snapshot.static_events)
        if snapshot.analysis_mode in (:dynamic, :kinematic) &&
                !isempty(snapshot.static_events)
            [simulation_mechanism_result(partial;
                 analysis_mode = :static),
             simulation_mechanism_result(partial)],
                ["Static initialization", "Running dynamics"], "Analysis"
        else
            [simulation_mechanism_result(partial)],
                ["Running result"], "Result"
        end
    else
        [model_mechanism_result(snapshot.loaded; configuration = :entered)],
            ["Model input"], "Configuration"
    end
    (; results, labels, choice_name, snapshot)
end

function model_run_document(run::ModelRun)
    prepared = run_mechanisms(run)
    document = viewer_document(prepared.results;
        labels = prepared.labels, choice_name = prepared.choice_name,
        include_signals = true)
    for (choice, label) in zip(document["choices"], prepared.labels)
        startswith(label, "Static") &&
            (choice["time_label"] = "static history sample")
    end
    document["run_status"] = String(prepared.snapshot.status)
    document["run_message"] = prepared.snapshot.message
    document["run_revision"] = prepared.snapshot.revision
    document["result_name"] = prepared.snapshot.result_name
    document
end

function finish_run_file!(run::ModelRun, result)
    finish_incremental_result!(run.writer, :complete)
    write_result(run.result_path, result; overwrite = true)
    prepared = run_mechanisms(run)
    write_graphics(run.result_path, prepared.results;
        labels = prepared.labels, choice_name = prepared.choice_name)
end

function fail_run_file!(run::ModelRun, exception)
    isnothing(run.writer) && return
    failure_time = isempty(run.times) ? run.loaded.simulation.start_time :
        last(run.times)
    finish_incremental_result!(run.writer, :failed;
        message = sprint(showerror, exception), failure_time)
    prepared = run_mechanisms(run)
    write_graphics(run.result_path, prepared.results;
        labels = prepared.labels, choice_name = prepared.choice_name)
end

function execute_model_run!(run::ModelRun, end_time, samples)
    lock(run.lock) do
        run.status = :running
        run.revision += 1
    end
    try
        result = if run.analysis_mode == :initial_conditions
            initial_conditions_result(run)
        elseif run.analysis_mode == :modal && !isnothing(run.starting_result)
            modal_result_from_previous(run)
        elseif !isnothing(run.starting_result)
            dynamic_result_from_previous(run, end_time, samples)
        elseif run.loaded isa LoadedSpatialModel
            run_spatial_model(run.loaded; end_time, samples,
                result_progress = event -> record_run_event!(run, event))
        else
            run_planar_model(run.loaded; end_time, samples,
                result_progress = event -> record_run_event!(run, event))
        end
        lock(run.lock) do
            run.result = result
            finish_run_file!(run, result)
            run.status = :complete
            run.message = "Analysis complete"
            run.revision += 1
        end
    catch exception
        lock(run.lock) do
            run.status = exception isa InterruptException ? :interrupted : :failed
            run.message = sprint(showerror, exception)
            try
                fail_run_file!(run, exception)
            catch file_error
                run.message *= "; result finalization failed: " *
                    sprint(showerror, file_error)
            end
            run.revision += 1
        end
    end
    nothing
end

"""Start a model analysis in the background and return its run identifier."""
function start_model_run(source::AbstractString, name::AbstractString;
        end_time = nothing, frames_per_second = 60.0,
        analysis = nothing, initial_run_id = nothing)
    ncodeunits(source) <= MAXIMUM_MODEL_BYTES || throw(ArgumentError(
        "model source is larger than the 20 MB run limit"))
    viewer_loaded = cached_model(source, name)
    if isnothing(viewer_loaded)
        viewer_loaded = cache_model!(source, name,
            load_uploaded_model(source, name))
    end
    base_loaded = viewer_loaded isa LoadedSpatialModel ?
        load_spatial_model(IOBuffer(viewer_loaded.model_source);
            source_directory = REPOSITORY_ROOT) :
        load_planar_model(IOBuffer(viewer_loaded.model_source);
            source_directory = REPOSITORY_ROOT)
    requested = isnothing(analysis) ? requested_analysis_mode(base_loaded) :
        Symbol(analysis)
    requested in (:initial_conditions, :static, :dynamic, :modal) ||
        throw(ArgumentError(
            "analysis must be initial_conditions, static, dynamic, or modal"))
    loaded = loaded_for_analysis(base_loaded, requested)
    actual_mode = requested == :dynamic && loaded.analysis.mode == :kinematic ?
        :kinematic : requested
    starting_result = requested in (:dynamic, :modal) ?
        completed_starting_result(initial_run_id, loaded) : nothing
    start = isnothing(starting_result) ? loaded.simulation.start_time :
        starting_result.time
    finish = isnothing(end_time) ? loaded.simulation.end_time : Float64(end_time)
    finish >= start || throw(ArgumentError(
        "end time must not precede the model start time"))
    fps = Float64(frames_per_second)
    isfinite(fps) && fps > 0 || throw(ArgumentError(
        "frames per second must be positive and finite"))
    samples = requested == :dynamic ?
        max(1, round(Int, (finish - start) * fps) + 1) : 1
    requested == :dynamic || (finish = start)
    id = next_run_id()
    directory = mktempdir(; prefix = "simpview-run-")
    result_name = result_file_name(name)
    run = ModelRun(id, result_name, joinpath(directory, result_name), loaded,
        viewer_loaded,
        actual_mode, starting_result, :queued, "Waiting to start", 0,
        isnothing(starting_result) ? Float64[] : copy(starting_result.times),
        isnothing(starting_result) ? Vector{Vector{Float64}}() :
            copy.(starting_result.states),
        Any[], false, nothing, nothing, nothing,
        ReentrantLock())
    lock(RUNS_LOCK) do
        RUNS[id] = run
    end
    run.task = Threads.@spawn execute_model_run!(run, finish, samples)
    (; id, result_name, analysis_mode = actual_mode,
       end_time = finish, frames_per_second = fps, samples)
end

function find_run(id)
    lock(RUNS_LOCK) do
        get(RUNS, String(id), nothing)
    end
end

"""Return a complete browser snapshot, or `nothing` when unchanged."""
function model_run_snapshot(id::AbstractString, after_revision = -1)
    run = find_run(id)
    isnothing(run) && throw(KeyError(id))
    current_revision = lock(run.lock) do
        run.revision
    end
    current_revision <= after_revision && return nothing
    document = model_run_document(run)
    revision = Int(document["run_revision"])
    status = Symbol(document["run_status"])
    message = String(document["run_message"])
    (; id = run.id, revision, status, message,
       analysis_mode = run.analysis_mode,
       result_ready = status in (:complete, :failed, :interrupted) &&
            isfile(run.result_path),
       result_name = run.result_name, document)
end

function response_headers(request)
    origin = HTTP.header(request, "Origin", "")
    headers = ["Content-Type" => "application/json; charset=utf-8",
        "Cache-Control" => "no-store"]
    allowed_origin(origin) && append!(headers, [
        "Access-Control-Allow-Origin" => origin,
        "Access-Control-Allow-Headers" => "Content-Type",
        "Access-Control-Allow-Methods" => "GET, POST, OPTIONS",
        "Vary" => "Origin",
    ])
    headers
end

json_response(request, status, value) =
    HTTP.Response(status, response_headers(request), JSON.json(value))

function request_after_revision(target)
    query = length(split(target, '?'; limit = 2)) == 2 ?
        split(target, '?'; limit = 2)[2] : ""
    for field in split(query, '&')
        pair = split(field, '='; limit = 2)
        length(pair) == 2 && pair[1] == "after" &&
            return parse(Int, pair[2])
    end
    -1
end

function result_response(request, run)
    isfile(run.result_path) || return json_response(request, 404,
        Dict("error" => "the run does not yet have a result file"))
    headers = response_headers(request)
    filter!(pair -> first(pair) != "Content-Type", headers)
    append!(headers, [
        "Content-Type" => "application/x-hdf5",
        "Content-Disposition" =>
            "attachment; filename=\"$(run.result_name)\"",
        "Cache-Control" => "no-store",
    ])
    HTTP.Response(200, headers, read(run.result_path))
end

function request_handler(request)
    raw_target = String(request.target)
    target = first(split(raw_target, '?'; limit = 2))
    request.method == "OPTIONS" &&
        return HTTP.Response(204, response_headers(request))
    if request.method == "GET" && target == "/api/health"
        return json_response(request, 200,
            Dict("service" => "SimpView", "status" => "ready"))
    end
    if request.method == "POST" && target == "/api/preview"
        length(request.body) <= MAXIMUM_MODEL_BYTES ||
            return json_response(request, 413,
                Dict("error" => "model upload is too large"))
        try
            payload = JSON.parse(String(request.body))
            source = get(payload, "source", nothing)
            name = get(payload, "name", nothing)
            source isa AbstractString && name isa AbstractString ||
                throw(ArgumentError(
                    "preview request requires string name and source fields"))
            return json_response(request, 200,
                preview_model_document(source, name))
        catch error
            return json_response(request, 400,
                Dict("error" => sprint(showerror, error)))
        end
    end
    if request.method == "POST" && target == "/api/run"
        length(request.body) <= MAXIMUM_MODEL_BYTES ||
            return json_response(request, 413,
                Dict("error" => "model run request is too large"))
        try
            payload = JSON.parse(String(request.body))
            source = get(payload, "source", nothing)
            name = get(payload, "name", nothing)
            source isa AbstractString && name isa AbstractString ||
                throw(ArgumentError(
                    "run request requires string name and source fields"))
            started = start_model_run(source, name;
                end_time = get(payload, "end_time", nothing),
                frames_per_second = get(payload,
                    "frames_per_second", 60.0),
                analysis = get(payload, "analysis", nothing),
                initial_run_id = get(payload, "initial_run_id", nothing))
            return json_response(request, 202, started)
        catch error
            return json_response(request, 400,
                Dict("error" => sprint(showerror, error)))
        end
    end
    run_match = match(r"^/api/run/([^/]+)$", target)
    if request.method == "GET" && !isnothing(run_match)
        run = find_run(run_match.captures[1])
        isnothing(run) && return json_response(request, 404,
            Dict("error" => "unknown SimpView run"))
        snapshot = try
            model_run_snapshot(run.id, request_after_revision(raw_target))
        catch error
            return json_response(request, 500,
                Dict("error" => sprint(showerror, error)))
        end
        isnothing(snapshot) && return HTTP.Response(204,
            response_headers(request))
        return json_response(request, 200, snapshot)
    end
    result_match = match(r"^/api/run/([^/]+)/result$", target)
    if request.method == "GET" && !isnothing(result_match)
        run = find_run(result_match.captures[1])
        isnothing(run) && return json_response(request, 404,
            Dict("error" => "unknown SimpView run"))
        return result_response(request, run)
    end
    json_response(request, 404, Dict("error" => "unknown SimpView endpoint"))
end

"""Serve the local preview API on the loopback interface."""
function serve_simpview(; host = "127.0.0.1", port = 8123)
    println("SimpView model service ready at http://$host:$port")
    HTTP.serve(request_handler, host, port; verbose = false)
end

end
