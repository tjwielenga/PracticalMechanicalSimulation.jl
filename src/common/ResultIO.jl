"""
    ResultIO

Versioned HDF5 storage for `.simp` results. Files contain portable numeric
arrays and metadata rather than serialized Julia objects: canonical variable
and equation catalogs, active status, the original TOML model, analysis and
solver diagnostics, health snapshots, state-selection changes, and optional
modal data.

Keeping the source model in each result supports inspection, extraction, and
name-and-type checked initialization of a later analysis.
"""
module ResultIO

using HDF5

export StoredBodyReference, StoredHealthSnapshot, StoredStateSelectionChange,
       StoredStaticSnapshot, StoredSimulationResult, IncrementalResultWriter,
       begin_incremental_result, record_incremental_event!,
       finish_incremental_result!, finalize_result!, write_result, read_result,
       HealthPeakCollector, record_health_step!, finish_health_episode!

const RESULT_FORMAT = "PracticalMechanicalSimulation"
const RESULT_FORMAT_VERSION = 2
const HEALTH_WARNING_ERROR = 10.0
const HEALTH_CHECK_ERROR = 25.0
const HEALTH_AMPLIFICATION = 10.0

"""One accepted active state retained at the peak of an unhealthy episode."""
struct RawHealthPeak
    time::Float64
    active_values::Vector{Float64}
    physical_error::Float64
    controlled_error::Float64
    dominant_local_index::Int
    order::Int
    step_size::Float64
end

"""Largest monitored predictor error, retained without a state copy."""
struct PhysicalErrorPeak
    time::Float64
    error::Float64
    dominant_local_index::Int
end

"""
Bounded collector for physical-error diagnostics from accepted BDF steps.

Only the largest state in each contiguous unhealthy episode is copied. Scalar
information for the overall largest monitored error is retained separately.
"""
mutable struct HealthPeakCollector
    peaks::Vector{RawHealthPeak}
    current::Union{Nothing,RawHealthPeak}
    maximum::Union{Nothing,PhysicalErrorPeak}
    previous_order::Int
    previous_step::Float64
end

HealthPeakCollector() = HealthPeakCollector(
    RawHealthPeak[], nothing, nothing, 0, 0.0)

function finish_health_episode!(collector::HealthPeakCollector)
    isnothing(collector.current) || push!(collector.peaks, collector.current)
    collector.current = nothing
    collector
end

"""Inspect one accepted step and retain only bounded health information."""
function record_health_step!(collector::HealthPeakCollector, time, values,
        order, step, controlled_error, physical_error,
        dominant_local_index)
    valid = isfinite(physical_error) && isfinite(controlled_error) &&
        dominant_local_index > 0
    if valid && (isnothing(collector.maximum) ||
            physical_error > collector.maximum.error)
        collector.maximum = PhysicalErrorPeak(Float64(time),
            Float64(physical_error), Int(dominant_local_index))
    end
    stable_formula = collector.previous_order != 0 &&
        order == collector.previous_order && begin
            ratio = iszero(collector.previous_step) ? 1.0 :
                abs(step / collector.previous_step)
            2 / 3 <= ratio <= 3 / 2
        end
    amplification = valid ? physical_error / max(controlled_error, 0.1) : 0.0
    unhealthy = valid && stable_formula &&
        physical_error > HEALTH_WARNING_ERROR &&
        amplification > HEALTH_AMPLIFICATION
    if unhealthy
        if isnothing(collector.current) ||
                physical_error > collector.current.physical_error
            collector.current = RawHealthPeak(Float64(time), collect(values),
                Float64(physical_error), Float64(controlled_error),
                Int(dominant_local_index), Int(order), Float64(step))
        end
    else
        finish_health_episode!(collector)
    end
    collector.previous_order = Int(order)
    collector.previous_step = Float64(step)
    collector
end

"""One accepted state retained at the peak of an unhealthy episode."""
struct StoredHealthSnapshot
    time::Float64
    values::Vector{Float64}
    physical_error::Float64
    controlled_error::Float64
    amplification::Float64
    dominant_variable::String
    dominant_variable_index::Int
    order::Int
    step_size::Float64
    severity::Symbol
    selected_variables::Vector{String}
end

"""One initial or runtime choice of the independent velocity coordinates."""
struct StoredStateSelectionChange
    time::Float64
    reason::Symbol
    previous::Vector{String}
    selected::Vector{String}
end

"""One retained state from static relaxation or Newton correction."""
struct StoredStaticSnapshot
    phase::Symbol
    status::Symbol
    model_time::Float64
    pseudo_time::Union{Nothing,Float64}
    iteration::Int
    relaxation_cycle::Int
    values::Vector{Float64}
    force_imbalance::Float64
    torque_imbalance::Float64
    constraint_error::Float64
    equivalent_acceleration::Float64
    reciprocal_condition::Float64
    mass_regularized::Bool
    message::String
end

"""Constant transform from a modeled body reference frame to its CM marker."""
struct StoredBodyReference
    name::String
    center_of_mass_marker::Union{Nothing,String}
    center_of_mass_position::Vector{Float64}
    center_of_mass_orientation::Matrix{Float64}
end

"""
    StoredSimulationResult

Portable data read from a versioned `.simp` HDF5 file.

Rows of `values` correspond to `times`; columns correspond to the parallel
variable metadata arrays. The structure also retains equation metadata,
active/inactive flags, the verbatim TOML model, analysis information, and
solver statistics without serializing executable Julia objects.
"""
struct StoredSimulationResult
    title::String
    times::Vector{Float64}
    values::Matrix{Float64}
    variable_names::Vector{String}
    variable_components::Vector{String}
    variable_kinds::Vector{String}
    variable_levels::Vector{Int}
    variable_scales::Vector{Float64}
    variable_active::BitVector
    equation_names::Vector{String}
    equation_components::Vector{String}
    equation_kinds::Vector{String}
    equation_levels::Vector{Int}
    equation_active::BitVector
    model_source::String
    body_references::Vector{StoredBodyReference}
    initial_consistency_iterations::Int
    static_initialization_iterations::Union{Nothing,Int}
    static_relaxation_cycles::Union{Nothing,Int}
    static_mass_regularized::Union{Nothing,Bool,BitVector}
    analysis_mode::Symbol
    degrees_of_freedom::Int
    solver_statistics::NamedTuple
    health_snapshots::Vector{StoredHealthSnapshot}
    state_selection_changes::Vector{StoredStateSelectionChange}
    static_snapshots::Vector{StoredStaticSnapshot}
    modal_eigenvalues::Vector{ComplexF64}
    modal_natural_frequencies_hz::Vector{Float64}
    modal_damped_frequencies_hz::Vector{Float64}
    modal_damping_ratios::Vector{Float64}
    modal_mode_shapes::Matrix{ComplexF64}
    modal_equation_errors::Vector{Float64}
    modal_shift::ComplexF64
    output_precision::Symbol
    status::Symbol
    status_message::String
    failure_time::Union{Nothing,Float64}
end

"""Buffered writer that keeps an in-progress `.simp` file valid on disk."""
mutable struct IncrementalResultWriter
    path::String
    variable_count::Int
    storage_type::DataType
    compression::Int
    flush_interval::Int
    sample_times::Vector{Float64}
    sample_states::Vector{Vector{Float64}}
    static_events::Vector{Any}
    sample_count::Int
    static_count::Int
    last_sample_time::Union{Nothing,Float64}
    closed::Bool
end

function selected_variable_names(result, time)
    if hasproperty(result, :state_selection_changes) &&
            !isempty(result.state_selection_changes)
        changes = result.state_selection_changes
        selected = copy(first(changes).selected)
        for change in Iterators.drop(changes, 1)
            if change.time == time
                selected = change.previous
                break
            elseif change.time < time
                selected = change.selected
            else
                break
            end
        end
        names = String[]
        for velocity in selected
            text = String(velocity)
            push!(names, text)
            suffix = endswith(text, ".V_x") ? ".R_x" :
                endswith(text, ".V_y") ? ".R_y" : ".theta"
            push!(names, replace(text, r"\.(V_x|V_y|omega)$" => suffix))
        end
        return names
    end
    solution = result.solution
    isnothing(solution.differential_vars) && return String[]
    variables = result.loaded.layout.catalog.variables
    [begin
         canonical = result.loaded.active_variable_indices[local_index]
         variable = variables[canonical]
         string(variable.component, ".", variable.name)
     end for local_index in findall(solution.differential_vars)]
end

function health_episode_peaks(result)
    hasproperty(result, :solution) || return Int[]
    solution = result.solution
    monitor = solution.error_monitor
    isnothing(monitor) && return Int[]
    peaks = Int[]
    peak = 0
    for index in 2:length(solution.t)
        physical = monitor.maximum_errors[index]
        controlled = monitor.controlled_errors[index]
        valid = isfinite(physical) && isfinite(controlled)
        amplification = valid ? physical / max(controlled, 0.1) : 0.0
        stable_formula = index > 2 && begin
            previous_step = abs(solution.steps[index - 1])
            ratio = iszero(previous_step) ? 1.0 :
                abs(solution.steps[index]) / previous_step
            solution.orders[index] == solution.orders[index - 1] &&
                2 / 3 <= ratio <= 3 / 2
        end
        unhealthy = valid && stable_formula &&
            physical > HEALTH_WARNING_ERROR &&
            amplification > HEALTH_AMPLIFICATION
        if unhealthy
            if peak == 0 || physical > monitor.maximum_errors[peak]
                peak = index
            end
        elseif peak != 0
            push!(peaks, peak)
            peak = 0
        end
    end
    peak != 0 && push!(peaks, peak)
    peaks
end

function health_snapshots(result)
    if hasproperty(result, :health_step_peaks)
        peaks = result.health_step_peaks
        isempty(peaks) && return StoredHealthSnapshot[]
        active = result.loaded.active_variable_indices
        variables = result.loaded.layout.catalog.variables
        return [begin
            canonical_variable = active[peak.dominant_local_index]
            variable = variables[canonical_variable]
            values = copy(first(result.states))
            values[active] .= peak.active_values
            StoredHealthSnapshot(peak.time, values, peak.physical_error,
                peak.controlled_error,
                peak.physical_error / max(peak.controlled_error, 0.1),
                string(variable.component, ".", variable.name),
                canonical_variable, peak.order, peak.step_size,
                peak.physical_error > HEALTH_CHECK_ERROR ? :check : :warning,
                selected_variable_names(result, peak.time))
        end for peak in peaks]
    end
    indices = health_episode_peaks(result)
    isempty(indices) && return StoredHealthSnapshot[]
    solution = result.solution
    monitor = solution.error_monitor
    active = result.loaded.active_variable_indices
    variables = result.loaded.layout.catalog.variables
    [begin
         local_variable = monitor.maximum_error_indices[index]
         canonical_variable = active[local_variable]
         variable = variables[canonical_variable]
         physical = monitor.maximum_errors[index]
         controlled = monitor.controlled_errors[index]
         values = copy(first(result.states))
         values[active] .= solution.y[index]
         StoredHealthSnapshot(solution.t[index], values, physical,
             controlled, physical / max(controlled, 0.1),
             string(variable.component, ".", variable.name),
             canonical_variable, solution.orders[index],
             solution.steps[index],
             physical > HEALTH_CHECK_ERROR ? :check : :warning,
             selected_variable_names(result, solution.t[index]))
     end for index in indices]
end

function output_storage_type(precision)
    value = Symbol(precision)
    value == :single && return Float32
    value == :double && return Float64
    throw(ArgumentError("output precision must be single or double"))
end

function result_output_precision(result, override = nothing)
    isnothing(override) ? result.loaded.simulation.output_precision :
        Symbol(override)
end

function history_matrix(states, ::Type{T} = Float64) where {T<:AbstractFloat}
    isempty(states) && return Matrix{T}(undef, 0, 0)
    values = Matrix{T}(undef, length(states), length(first(states)))
    for (row, state) in enumerate(states)
        values[row, :] .= state
    end
    values
end

function replace_attribute!(object, name, value)
    object_attributes = attributes(object)
    haskey(object_attributes, name) && delete_attribute(object, name)
    object_attributes[name] = value
end

function write_body_references!(model, loaded)
    hasproperty(loaded, :body_reference_frames) || return
    references = loaded.body_reference_frames
    names = sort!(collect(keys(references)))
    body_references = create_group(model, "body_references")
    body_references["name"] = string.(names)
    body_references["center_of_mass_marker"] = [
        isnothing(references[name].center_of_mass_marker) ? "" :
            string(references[name].center_of_mass_marker)
        for name in names]
    positions = zeros(Float64, length(names), 3)
    orientations = zeros(Float64, length(names), 3, 3)
    for (index, name) in enumerate(names)
        positions[index, :] .= references[name].center_of_mass_position
        orientations[index, :, :] .=
            references[name].center_of_mass_orientation
    end
    body_references["center_of_mass_position"] = positions
    body_references["center_of_mass_orientation"] = orientations
end

function write_catalog!(file, loaded)
    variables = loaded.layout.catalog.variables
    equations = loaded.layout.catalog.equations
    metadata = create_group(file, "variables")
    metadata["name"] = string.(getproperty.(variables, :name))
    metadata["component"] = string.(getproperty.(variables, :component))
    metadata["kind"] = string.(getproperty.(variables, :kind))
    metadata["level"] = getproperty.(variables, :level)
    metadata["scale"] = getproperty.(variables, :scale)
    active_variables = Set(loaded.active_variable_indices)
    metadata["active"] = UInt8[
        variable.index in active_variables for variable in variables]
    equation_metadata = create_group(file, "equations")
    equation_metadata["name"] = string.(getproperty.(equations, :name))
    equation_metadata["component"] =
        string.(getproperty.(equations, :component))
    equation_metadata["kind"] = string.(getproperty.(equations, :kind))
    equation_metadata["level"] = getproperty.(equations, :level)
    active_equations = Set(loaded.active_equation_indices)
    equation_metadata["active"] = UInt8[
        equation.index in active_equations for equation in equations]
end

function create_extendible_dataset(parent, name, ::Type{T}, dimensions,
        maximum_dimensions; chunk, compression = nothing) where {T}
    keywords = isnothing(compression) ? (; chunk) :
        (; chunk, shuffle = (), deflate = compression)
    create_dataset(parent, name, datatype(T),
        (dimensions, maximum_dimensions); keywords...)
end

function create_static_progress!(diagnostics, variable_count, storage_type,
        compression)
    progress = create_group(diagnostics, "static_progress")
    for name in ("phase", "status", "message")
        create_extendible_dataset(progress, name, String, (0,), (-1,);
            chunk = (16,))
    end
    for name in ("model_time", "pseudo_time", "force_imbalance",
            "torque_imbalance", "constraint_error",
            "equivalent_acceleration", "reciprocal_condition")
        create_extendible_dataset(progress, name, Float64, (0,), (-1,);
            chunk = (16,))
    end
    for name in ("iteration", "relaxation_cycle")
        create_extendible_dataset(progress, name, Int, (0,), (-1,);
            chunk = (16,))
    end
    create_extendible_dataset(progress, "mass_regularized", UInt8,
        (0,), (-1,); chunk = (16,))
    create_extendible_dataset(progress, "values", storage_type,
        (0, variable_count), (-1, variable_count);
        chunk = (8, max(variable_count, 1)), compression)
end

"""
    begin_incremental_result(path, loaded, analysis_mode; overwrite=false,
                             initial_time=loaded.simulation.start_time,
                             initial_state=loaded.initial_values)

Create a valid, initially `running` result file before analysis starts. The
initial canonical state is stored immediately so even an early failure leaves
a viewable model configuration.
"""
function begin_incremental_result(path::AbstractString, loaded,
        analysis_mode::Symbol; overwrite = false, compression = 3,
        flush_interval = 8, initial_time = loaded.simulation.start_time,
        initial_state = loaded.initial_values,
        output_precision = loaded.simulation.output_precision)
    isfile(path) && !overwrite &&
        throw(ArgumentError("result file already exists: $path"))
    variable_count = length(loaded.layout.catalog.variables)
    precision = Symbol(output_precision)
    storage_type = output_storage_type(precision)
    h5open(path, "w") do file
        attributes(file)["format"] = RESULT_FORMAT
        attributes(file)["format_version"] = RESULT_FORMAT_VERSION
        attributes(file)["title"] = loaded.title
        attributes(file)["status"] = "running"
        attributes(file)["status_message"] = ""
        attributes(file)["output_precision"] = String(precision)
        model = create_group(file, "model")
        model["toml"] = loaded.model_source
        write_body_references!(model, loaded)
        simulation = create_group(file, "simulation")
        create_extendible_dataset(simulation, "time", Float64,
            (0,), (-1,); chunk = (256,))
        results = create_group(file, "results")
        create_extendible_dataset(results, "values", storage_type,
            (0, variable_count), (-1, variable_count);
            chunk = (8, max(variable_count, 1)), compression)
        write_catalog!(file, loaded)
        diagnostics = create_group(file, "diagnostics")
        attributes(diagnostics)["analysis_mode"] = String(analysis_mode)
        attributes(diagnostics)["degrees_of_freedom"] =
            loaded.analysis.degrees_of_freedom
        diagnostics["initial_consistency_iterations"] = 0
        create_static_progress!(diagnostics, variable_count, storage_type,
            compression)
    end
    writer = IncrementalResultWriter(abspath(path), variable_count,
        storage_type,
        compression, max(Int(flush_interval), 1), Float64[],
        Vector{Float64}[], Any[], 0, 0, nothing, false)
    record_incremental_event!(writer, (; kind = :sample,
        time = Float64(initial_time),
        state = copy(initial_state)))
    flush_incremental_result!(writer)
    writer
end

function append_rows!(dataset, values::AbstractVector)
    isempty(values) && return
    old = size(dataset, 1)
    HDF5.set_extent_dims(dataset, (old + length(values),))
    dataset[(old + 1):(old + length(values))] = values
end

function append_rows!(dataset, values::AbstractMatrix)
    size(values, 1) == 0 && return
    old = size(dataset, 1)
    HDF5.set_extent_dims(dataset, (old + size(values, 1), size(dataset, 2)))
    dataset[(old + 1):(old + size(values, 1)), :] = values
end

function event_value(event, name, default)
    hasproperty(event, name) ? getproperty(event, name) : default
end

function flush_incremental_result!(writer::IncrementalResultWriter)
    writer.closed && return writer
    isempty(writer.sample_times) && isempty(writer.static_events) && return writer
    h5open(writer.path, "r+") do file
        if !isempty(writer.sample_times)
            append_rows!(file["simulation/time"], writer.sample_times)
            append_rows!(file["results/values"],
                history_matrix(writer.sample_states, writer.storage_type))
            writer.sample_count += length(writer.sample_times)
            empty!(writer.sample_times)
            empty!(writer.sample_states)
        end
        if !isempty(writer.static_events)
            progress = file["diagnostics/static_progress"]
            events = writer.static_events
            append_rows!(progress["phase"],
                string.(getproperty.(events, :phase)))
            append_rows!(progress["status"],
                string.(getproperty.(events, :status)))
            append_rows!(progress["model_time"], Float64.(
                getproperty.(events, :model_time)))
            append_rows!(progress["pseudo_time"], Float64[
                isnothing(event_value(event, :pseudo_time, nothing)) ? NaN :
                    event.pseudo_time for event in events])
            append_rows!(progress["iteration"], Int.(
                getproperty.(events, :iteration)))
            append_rows!(progress["relaxation_cycle"], Int.(
                getproperty.(events, :relaxation_cycle)))
            for name in (:force_imbalance, :torque_imbalance,
                    :constraint_error, :equivalent_acceleration,
                    :reciprocal_condition)
                append_rows!(progress[String(name)], Float64[
                    event_value(event, name, NaN) for event in events])
            end
            append_rows!(progress["mass_regularized"], UInt8[
                event_value(event, :mass_regularized, false)
                for event in events])
            append_rows!(progress["message"], String[
                String(event_value(event, :message, "")) for event in events])
            append_rows!(progress["values"], history_matrix(
                [event.state for event in events], writer.storage_type))
            writer.static_count += length(events)
            empty!(writer.static_events)
        end
    end
    writer
end

function replace_last_sample!(writer, time, state)
    if !isempty(writer.sample_times)
        writer.sample_times[end] = time
        writer.sample_states[end] = copy(state)
    elseif writer.sample_count > 0
        h5open(writer.path, "r+") do file
            file["simulation/time"][writer.sample_count] = time
            file["results/values"][writer.sample_count, :] =
                writer.storage_type.(state)
        end
    end
end

"""
    record_incremental_event!(writer, event)

Buffer one accepted dynamic sample, failure state, or static-progress event.
Samples at the same time replace one another so an initial placeholder can be
updated by the corrected state. Reaching the flush interval, or recording a
failure, writes the buffer while leaving the HDF5 file valid and readable.
"""
function record_incremental_event!(writer::IncrementalResultWriter, event)
    writer.closed && return writer
    if event.kind in (:sample, :failure)
        time = Float64(event.time)
        state = Float64.(event.state)
        length(state) == writer.variable_count || throw(DimensionMismatch(
            "incremental result state has the wrong length"))
        if !isnothing(writer.last_sample_time) &&
                time == writer.last_sample_time
            replace_last_sample!(writer, time, state)
        else
            push!(writer.sample_times, time)
            push!(writer.sample_states, state)
            writer.last_sample_time = time
        end
    elseif event.kind == :static
        push!(writer.static_events, event)
    end
    if length(writer.sample_times) >= writer.flush_interval ||
            length(writer.static_events) >= writer.flush_interval ||
            event.kind == :failure
        flush_incremental_result!(writer)
    end
    writer
end

function write_solver_statistics!(file, statistics)
    solver = haskey(file["diagnostics"], "solver") ?
        file["diagnostics/solver"] : create_group(file["diagnostics"], "solver")
    for name in (:accepted_steps, :rejected_steps, :residual_evaluations,
            :jacobian_evaluations, :factorizations,
            :symbolic_factorizations, :newton_iterations,
            :corrector_failures, :root_evaluations, :events_found,
            :history_restarts, :state_reselections)
        hasproperty(statistics, name) || continue
        replace_attribute!(solver, String(name), Int(getproperty(statistics, name)))
    end
end

"""
    finish_incremental_result!(writer, status; ...)

Flush pending data, store terminal status and solver statistics, and close the
logical writer. The file remains the same result whether `status` is
`:complete`, `:failed`, or `:interrupted`; accepted history is never renamed
to a separate partial-result format.
"""
function finish_incremental_result!(writer::IncrementalResultWriter,
        status::Symbol; message = "", failure_time = nothing,
        statistics = nothing)
    writer.closed && return writer.path
    flush_incremental_result!(writer)
    h5open(writer.path, "r+") do file
        replace_attribute!(file, "status", String(status))
        replace_attribute!(file, "status_message", String(message))
        isnothing(failure_time) || replace_attribute!(file, "failure_time",
            Float64(failure_time))
        isnothing(statistics) || write_solver_statistics!(file, statistics)
    end
    writer.closed = true
    writer.path
end

"""
    finalize_result!(path; message="run finalized by the user")

Mark an abandoned `running` result as `interrupted`, retaining every sample
and static correction already flushed to the file. This repairs the metadata
after a hard process termination; it cannot recover an accepted state that had
not yet reached the file. Already interrupted or failed results are left
unchanged, while a completed result is rejected.
"""
function finalize_result!(path::AbstractString;
        message = "run finalized by the user")
    h5open(path, "r+") do file
        file_attributes = attributes(file)
        haskey(file_attributes, "format") &&
            read(file_attributes["format"]) == RESULT_FORMAT ||
            throw(ArgumentError(
                "not a PracticalMechanicalSimulation result file"))
        status = haskey(file_attributes, "status") ?
            Symbol(read(file_attributes["status"])) : :complete
        status in (:failed, :interrupted) && return path
        status == :running || throw(ArgumentError(
            "only a running result can be finalized; current status is $status"))
        times = read(file["simulation/time"])
        failure_time = isempty(times) ? nothing : Float64(last(times))
        replace_attribute!(file, "status", "interrupted")
        replace_attribute!(file, "status_message", String(message))
        isnothing(failure_time) || replace_attribute!(file, "failure_time",
            failure_time)
    end
    path
end

"""
    write_result(path, result; overwrite=false, compression=3,
                 output_precision=nothing)

Write a completed simulation result as compact, versioned HDF5.

Only requested output samples are stored. Existing files require
`overwrite=true`; `compression` is the HDF5 deflate level used for the canonical
history matrix. Returns `path`.
"""
function write_result(path::AbstractString, result; overwrite = false,
        compression = 3, output_precision = nothing)
    isfile(path) && !overwrite &&
        throw(ArgumentError("result file already exists: $path"))
    variables = result.loaded.layout.catalog.variables
    equations = result.loaded.layout.catalog.equations
    precision = result_output_precision(result, output_precision)
    storage_type = output_storage_type(precision)
    values = history_matrix(result.states, storage_type)
    chunk_rows = max(1, min(size(values, 1), 256))
    chunk_columns = max(1, size(values, 2))
    h5open(path, "w") do file
        attributes(file)["format"] = RESULT_FORMAT
        attributes(file)["format_version"] = RESULT_FORMAT_VERSION
        attributes(file)["title"] = result.loaded.title
        attributes(file)["status"] = "complete"
        attributes(file)["status_message"] = ""
        attributes(file)["output_precision"] = String(precision)
        model = create_group(file, "model")
        model["toml"] = result.loaded.model_source
        if hasproperty(result.loaded, :body_reference_frames)
            references = result.loaded.body_reference_frames
            names = sort!(collect(keys(references)))
            body_references = create_group(model, "body_references")
            body_references["name"] = string.(names)
            body_references["center_of_mass_marker"] = [
                isnothing(references[name].center_of_mass_marker) ? "" :
                    string(references[name].center_of_mass_marker)
                for name in names]
            positions = zeros(Float64, length(names), 3)
            orientations = zeros(Float64, length(names), 3, 3)
            for (index, name) in enumerate(names)
                positions[index, :] .=
                    references[name].center_of_mass_position
                orientations[index, :, :] .=
                    references[name].center_of_mass_orientation
            end
            body_references["center_of_mass_position"] = positions
            body_references["center_of_mass_orientation"] = orientations
        end
        simulation = create_group(file, "simulation")
        simulation["time"] = Float64.(result.times)
        results = create_group(file, "results")
        results["values", chunk = (chunk_rows, chunk_columns),
                shuffle = (), deflate = compression] = values
        metadata = create_group(file, "variables")
        metadata["name"] = string.(getproperty.(variables, :name))
        metadata["component"] = string.(getproperty.(variables, :component))
        metadata["kind"] = string.(getproperty.(variables, :kind))
        metadata["level"] = getproperty.(variables, :level)
        metadata["scale"] = getproperty.(variables, :scale)
        active_variables = Set(hasproperty(result, :static_variable_indices) ?
            result.static_variable_indices : result.loaded.active_variable_indices)
        metadata["active"] = UInt8[
            variable.index in active_variables for variable in variables]
        equation_metadata = create_group(file, "equations")
        equation_metadata["name"] = string.(getproperty.(equations, :name))
        equation_metadata["component"] =
            string.(getproperty.(equations, :component))
        equation_metadata["kind"] = string.(getproperty.(equations, :kind))
        equation_metadata["level"] = getproperty.(equations, :level)
        active_equations = Set(hasproperty(result, :static_equation_indices) ?
            result.static_equation_indices : result.loaded.active_equation_indices)
        equation_metadata["active"] = UInt8[
            equation.index in active_equations for equation in equations]
        diagnostics = create_group(file, "diagnostics")
        attributes(diagnostics)["analysis_mode"] = String(result.analysis_mode)
        attributes(diagnostics)["degrees_of_freedom"] =
            result.loaded.analysis.degrees_of_freedom
        diagnostics["initial_consistency_iterations"] =
            result.initial_consistency_iterations
        if hasproperty(result, :static_initialization_iterations) &&
                !isnothing(result.static_initialization_iterations)
            diagnostics["static_initialization_iterations"] =
                result.static_initialization_iterations
        end
        if hasproperty(result, :static_relaxation_cycles) &&
                !isnothing(result.static_relaxation_cycles) &&
                result.static_relaxation_cycles isa Integer
            diagnostics["static_relaxation_cycles"] =
                result.static_relaxation_cycles
        end
        if hasproperty(result, :static_mass_regularized) &&
                result.static_mass_regularized isa Bool
            diagnostics["static_mass_regularized"] =
                UInt8(result.static_mass_regularized)
        end
        if hasproperty(result, :static_progress_events) &&
                !isempty(result.static_progress_events)
            events = result.static_progress_events
            progress = create_group(diagnostics, "static_progress")
            progress["phase"] = string.(getproperty.(events, :phase))
            progress["status"] = string.(getproperty.(events, :status))
            progress["model_time"] = Float64.(
                getproperty.(events, :model_time))
            progress["pseudo_time"] = Float64[
                isnothing(event.pseudo_time) ? NaN : event.pseudo_time
                for event in events]
            progress["iteration"] = Int.(getproperty.(events, :iteration))
            progress["relaxation_cycle"] = Int.(
                getproperty.(events, :relaxation_cycle))
            progress["force_imbalance"] = Float64.(
                getproperty.(events, :force_imbalance))
            progress["torque_imbalance"] = Float64.(
                getproperty.(events, :torque_imbalance))
            progress["constraint_error"] = Float64.(
                getproperty.(events, :constraint_error))
            progress["equivalent_acceleration"] = Float64.(
                getproperty.(events, :equivalent_acceleration))
            progress["reciprocal_condition"] = Float64.(
                getproperty.(events, :reciprocal_condition))
            progress["mass_regularized"] = UInt8.(
                getproperty.(events, :mass_regularized))
            progress["message"] = String.(getproperty.(events, :message))
            progress_values = history_matrix(getproperty.(events, :state),
                storage_type)
            progress["values", chunk = (min(length(events), 64),
                    size(progress_values, 2)), shuffle = (),
                    deflate = compression] = progress_values
        end
        if hasproperty(result, :solution)
            statistics = result.solution.stats
            solver = create_group(diagnostics, "solver")
            attributes(solver)["accepted_steps"] = statistics.accepted_steps
            attributes(solver)["rejected_steps"] = statistics.rejected_steps
            attributes(solver)["residual_evaluations"] =
                statistics.residual_evaluations
            attributes(solver)["jacobian_evaluations"] =
                statistics.jacobian_evaluations
            attributes(solver)["factorizations"] = statistics.factorizations
            attributes(solver)["symbolic_factorizations"] =
                statistics.symbolic_factorizations
            attributes(solver)["newton_iterations"] = statistics.newton_iterations
            attributes(solver)["corrector_failures"] =
                statistics.corrector_failures
            attributes(solver)["root_evaluations"] = statistics.root_evaluations
            attributes(solver)["events_found"] = statistics.events_found
            attributes(solver)["history_restarts"] = statistics.history_restarts
            attributes(solver)["state_reselections"] =
                statistics.state_reselections
            if hasproperty(result, :state_selection_changes) &&
                    !isempty(result.state_selection_changes)
                changes = create_group(diagnostics, "state_selection")
                changes["time"] = getproperty.(
                    result.state_selection_changes, :time)
                changes["reason"] = string.(getproperty.(
                    result.state_selection_changes, :reason))
                changes["previous"] = [join(string.(change.previous), '\n')
                    for change in result.state_selection_changes]
                changes["selected"] = [join(string.(change.selected), '\n')
                    for change in result.state_selection_changes]
            end
            snapshots = health_snapshots(result)
            if !isempty(snapshots)
                health = create_group(diagnostics, "health")
                health["time"] = getproperty.(snapshots, :time)
                health["physical_error"] =
                    getproperty.(snapshots, :physical_error)
                health["controlled_error"] =
                    getproperty.(snapshots, :controlled_error)
                health["amplification"] =
                    getproperty.(snapshots, :amplification)
                health["dominant_variable"] =
                    getproperty.(snapshots, :dominant_variable)
                health["dominant_variable_index"] =
                    getproperty.(snapshots, :dominant_variable_index)
                health["order"] = getproperty.(snapshots, :order)
                health["step_size"] = getproperty.(snapshots, :step_size)
                health["severity"] = string.(getproperty.(snapshots, :severity))
                health["selected_variables"] =
                    join.(getproperty.(snapshots, :selected_variables), '\n')
                snapshot_values = history_matrix(
                    getproperty.(snapshots, :values), storage_type)
                health["values", chunk = (min(length(snapshots), 64),
                        size(snapshot_values, 2)), shuffle = (),
                        deflate = compression] = snapshot_values
            end
        elseif hasproperty(result, :static_iterations)
            static = create_group(diagnostics, "static")
            static["newton_iterations"] = result.static_iterations
            hasproperty(result, :static_relaxation_cycles) &&
                (static["relaxation_cycles"] = result.static_relaxation_cycles)
            hasproperty(result, :static_mass_regularized) &&
                (static["mass_regularized"] =
                    UInt8.(result.static_mass_regularized))
            attributes(static)["continuation_solutions"] =
                result.continuation_solutions
        end
        if result.analysis_mode == :modal
            modal = create_group(file, "modal")
            modal["eigenvalue_real"] = real.(result.eigenvalues)
            modal["eigenvalue_imaginary"] = imag.(result.eigenvalues)
            modal["natural_frequency_hz"] = result.natural_frequencies_hz
            modal["damped_frequency_hz"] = result.damped_frequencies_hz
            modal["damping_ratio"] = result.damping_ratios
            modal["mode_shape_real"] = storage_type.(real.(result.mode_shapes))
            modal["mode_shape_imaginary"] =
                storage_type.(imag.(result.mode_shapes))
            modal["equation_error"] = result.equation_errors
            attributes(modal)["shift_real"] = real(result.shift)
            attributes(modal)["shift_imaginary"] = imag(result.shift)
            attributes(modal)["normalization"] =
                "maximum scaled physical displacement"
            attributes(modal)["sparse_factorizations"] =
                result.sparse_factorizations
        end
    end
    path
end

function read_string_vector(dataset)
    String.(read(dataset))
end

function read_health_snapshots(file, variable_count)
    haskey(file["diagnostics"], "health") ||
        return StoredHealthSnapshot[]
    health = file["diagnostics/health"]
    times = Float64.(read(health["time"]))
    values = Matrix{Float64}(read(health["values"]))
    size(values) == (length(times), variable_count) ||
        throw(DimensionMismatch(
            "stored health snapshots have the wrong dimensions"))
    physical_errors = Float64.(read(health["physical_error"]))
    controlled_errors = Float64.(read(health["controlled_error"]))
    amplifications = Float64.(read(health["amplification"]))
    dominant_variables = read_string_vector(health["dominant_variable"])
    dominant_indices = Int.(read(health["dominant_variable_index"]))
    orders = Int.(read(health["order"]))
    step_sizes = Float64.(read(health["step_size"]))
    severities = Symbol.(read_string_vector(health["severity"]))
    selected = [isempty(value) ? String[] : String.(split(value, '\n'))
        for value in read_string_vector(health["selected_variables"])]
    metadata = (physical_errors, controlled_errors, amplifications,
        dominant_variables, dominant_indices, orders, step_sizes,
        severities, selected)
    all(length(values) == length(times) for values in metadata) ||
        throw(DimensionMismatch(
            "stored health snapshot metadata have inconsistent lengths"))
    [StoredHealthSnapshot(times[index], collect(@view(values[index, :])),
         physical_errors[index], controlled_errors[index],
         amplifications[index], dominant_variables[index],
         dominant_indices[index], orders[index], step_sizes[index],
         severities[index], selected[index])
     for index in eachindex(times)]
end

function read_static_snapshots(file, variable_count)
    haskey(file["diagnostics"], "static_progress") ||
        return StoredStaticSnapshot[]
    progress = file["diagnostics/static_progress"]
    phases = Symbol.(read_string_vector(progress["phase"]))
    statuses = Symbol.(read_string_vector(progress["status"]))
    model_times = Float64.(read(progress["model_time"]))
    stored_pseudo_times = Float64.(read(progress["pseudo_time"]))
    pseudo_times = Union{Nothing,Float64}[
        isnan(value) ? nothing : value for value in stored_pseudo_times]
    iterations = Int.(read(progress["iteration"]))
    cycles = Int.(read(progress["relaxation_cycle"]))
    force_imbalances = Float64.(read(progress["force_imbalance"]))
    torque_imbalances = Float64.(read(progress["torque_imbalance"]))
    constraint_errors = Float64.(read(progress["constraint_error"]))
    accelerations = Float64.(read(progress["equivalent_acceleration"]))
    conditions = Float64.(read(progress["reciprocal_condition"]))
    regularized = BitVector(!iszero(value) for value in
        read(progress["mass_regularized"]))
    messages = read_string_vector(progress["message"])
    values = Matrix{Float64}(read(progress["values"]))
    count = length(phases)
    metadata = (statuses, model_times, pseudo_times, iterations, cycles,
        force_imbalances, torque_imbalances, constraint_errors,
        accelerations, conditions, regularized, messages)
    all(length(field) == count for field in metadata) ||
        throw(DimensionMismatch(
            "stored static-progress metadata have inconsistent lengths"))
    size(values) == (count, variable_count) || throw(DimensionMismatch(
        "stored static-progress states have the wrong dimensions"))
    [StoredStaticSnapshot(phases[index], statuses[index], model_times[index],
         pseudo_times[index], iterations[index], cycles[index],
         collect(@view(values[index, :])), force_imbalances[index],
         torque_imbalances[index], constraint_errors[index],
         accelerations[index], conditions[index], regularized[index],
         messages[index]) for index in 1:count]
end

function read_state_selection_changes(file)
    haskey(file["diagnostics"], "state_selection") ||
        return StoredStateSelectionChange[]
    changes = file["diagnostics/state_selection"]
    times = Float64.(read(changes["time"]))
    reasons = Symbol.(read_string_vector(changes["reason"]))
    previous = [isempty(value) ? String[] : String.(split(value, '\n'))
        for value in read_string_vector(changes["previous"])]
    selected = [isempty(value) ? String[] : String.(split(value, '\n'))
        for value in read_string_vector(changes["selected"])]
    all(length(values) == length(times)
        for values in (reasons, previous, selected)) ||
        throw(DimensionMismatch(
            "stored state-selection metadata have inconsistent lengths"))
    [StoredStateSelectionChange(times[index], reasons[index],
        previous[index], selected[index]) for index in eachindex(times)]
end

function read_body_references(file)
    haskey(file["model"], "body_references") ||
        return StoredBodyReference[]
    references = file["model/body_references"]
    names = read_string_vector(references["name"])
    marker_names = read_string_vector(references["center_of_mass_marker"])
    positions = Matrix{Float64}(read(references["center_of_mass_position"]))
    orientations = Array{Float64,3}(
        read(references["center_of_mass_orientation"]))
    size(positions) == (length(names), 3) || throw(DimensionMismatch(
        "stored body-reference positions have the wrong dimensions"))
    size(orientations) == (length(names), 3, 3) ||
        throw(DimensionMismatch(
            "stored body-reference orientations have the wrong dimensions"))
    length(marker_names) == length(names) || throw(DimensionMismatch(
        "stored body-reference marker names have the wrong length"))
    [StoredBodyReference(names[index],
         isempty(marker_names[index]) ? nothing : marker_names[index],
         collect(@view(positions[index, :])),
         Matrix(@view(orientations[index, :, :])))
     for index in eachindex(names)]
end

"""
    read_result(path) -> StoredSimulationResult

Read and validate a PracticalMechanicalSimulation HDF5 result. The format name
and version are checked before portable arrays and metadata are returned.
"""
function read_result(path::AbstractString)
    h5open(path, "r") do file
        file_attributes = attributes(file)
        haskey(file_attributes, "format") &&
            read(file_attributes["format"]) == RESULT_FORMAT ||
            throw(ArgumentError("not a PracticalMechanicalSimulation result file"))
        version = Int(read(file_attributes["format_version"]))
        version in (1, RESULT_FORMAT_VERSION) ||
            throw(ArgumentError("unsupported result format version $version"))
        status = haskey(file_attributes, "status") ?
            Symbol(read(file_attributes["status"])) : :complete
        status_message = haskey(file_attributes, "status_message") ?
            String(read(file_attributes["status_message"])) : ""
        failure_time = haskey(file_attributes, "failure_time") ?
            Float64(read(file_attributes["failure_time"])) : nothing
        output_precision = haskey(file_attributes, "output_precision") ?
            Symbol(read(file_attributes["output_precision"])) :
            eltype(file["results/values"]) == Float32 ? :single : :double
        diagnostic_attributes = attributes(file["diagnostics"])
        solver_statistics = if haskey(file["diagnostics"], "solver")
            solver_attributes = attributes(file["diagnostics/solver"])
            optional_statistic(name) = haskey(solver_attributes, name) ?
                Int(read(solver_attributes[name])) : 0
            (; accepted_steps = Int(read(solver_attributes["accepted_steps"])),
               rejected_steps = Int(read(solver_attributes["rejected_steps"])),
               residual_evaluations = Int(read(solver_attributes["residual_evaluations"])),
               jacobian_evaluations = Int(read(solver_attributes["jacobian_evaluations"])),
               factorizations = Int(read(solver_attributes["factorizations"])),
               symbolic_factorizations =
                   optional_statistic("symbolic_factorizations"),
               newton_iterations = Int(read(solver_attributes["newton_iterations"])),
               corrector_failures = optional_statistic("corrector_failures"),
               root_evaluations = optional_statistic("root_evaluations"),
               events_found = optional_statistic("events_found"),
               history_restarts = optional_statistic("history_restarts"),
               state_reselections = optional_statistic("state_reselections"))
        elseif haskey(file["diagnostics"], "static")
            static = file["diagnostics/static"]
            (; newton_iterations = Int.(read(static["newton_iterations"])),
               relaxation_cycles = haskey(static, "relaxation_cycles") ?
                   Int.(read(static["relaxation_cycles"])) : Int[],
               continuation_solutions = Int(read(
                   attributes(static)["continuation_solutions"])))
        else
            (;)
        end
        modal_values = if haskey(file, "modal")
            modal = file["modal"]
            modal_attributes = attributes(modal)
            (; eigenvalues = ComplexF64.(read(modal["eigenvalue_real"]),
                    read(modal["eigenvalue_imaginary"])),
               natural_frequencies_hz =
                    Float64.(read(modal["natural_frequency_hz"])),
               damped_frequencies_hz =
                    Float64.(read(modal["damped_frequency_hz"])),
               damping_ratios = Float64.(read(modal["damping_ratio"])),
               mode_shapes = ComplexF64.(read(modal["mode_shape_real"]),
                    read(modal["mode_shape_imaginary"])),
               equation_errors = Float64.(read(modal["equation_error"])),
               shift = complex(
                    Float64(read(modal_attributes["shift_real"])),
                    Float64(read(modal_attributes["shift_imaginary"]))))
        else
            (; eigenvalues = ComplexF64[], natural_frequencies_hz = Float64[],
               damped_frequencies_hz = Float64[], damping_ratios = Float64[],
               mode_shapes = Matrix{ComplexF64}(undef, 0, 0),
               equation_errors = Float64[], shift = 0.0 + 0.0im)
        end
        variable_names = read_string_vector(file["variables/name"])
        health = read_health_snapshots(file, length(variable_names))
        state_selection_changes = read_state_selection_changes(file)
        static_snapshots = read_static_snapshots(file, length(variable_names))
        body_references = read_body_references(file)
        StoredSimulationResult(
            String(read(file_attributes["title"])),
            Float64.(read(file["simulation/time"])),
            Matrix{Float64}(read(file["results/values"])),
            variable_names,
            read_string_vector(file["variables/component"]),
            read_string_vector(file["variables/kind"]),
            Int.(read(file["variables/level"])),
            Float64.(read(file["variables/scale"])),
            BitVector(!iszero(value) for value in read(file["variables/active"])),
            read_string_vector(file["equations/name"]),
            read_string_vector(file["equations/component"]),
            read_string_vector(file["equations/kind"]),
            Int.(read(file["equations/level"])),
            BitVector(!iszero(value) for value in read(file["equations/active"])),
            String(read(file["model/toml"])),
            body_references,
            Int(read(file["diagnostics/initial_consistency_iterations"])),
            haskey(file["diagnostics"], "static_initialization_iterations") ?
                Int(read(file["diagnostics/static_initialization_iterations"])) :
                nothing,
            haskey(file["diagnostics"], "static_relaxation_cycles") ?
                Int(read(file["diagnostics/static_relaxation_cycles"])) : nothing,
            haskey(file["diagnostics"], "static_mass_regularized") ?
                !iszero(read(file["diagnostics/static_mass_regularized"])) :
                haskey(file, "diagnostics/static/mass_regularized") ?
                    BitVector(!iszero(value) for value in
                        read(file["diagnostics/static/mass_regularized"])) :
                    nothing,
            haskey(diagnostic_attributes, "analysis_mode") ?
                Symbol(read(diagnostic_attributes["analysis_mode"])) : :kinematic,
            haskey(diagnostic_attributes, "degrees_of_freedom") ?
                Int(read(diagnostic_attributes["degrees_of_freedom"])) : 0,
            solver_statistics,
            health,
            state_selection_changes,
            static_snapshots,
            modal_values.eigenvalues,
            modal_values.natural_frequencies_hz,
            modal_values.damped_frequencies_hz,
            modal_values.damping_ratios,
            modal_values.mode_shapes,
            modal_values.equation_errors,
            modal_values.shift,
            output_precision,
            status,
            status_message,
            failure_time)
    end
end

end
