"""
    Sim3D

Julia model-building interface for the spatial mechanism simulator. Builders
create the same hierarchical model document used by TOML and Lua input, then
pass it through the ordinary spatial loader and validation path.
"""
module Sim3D

using TOML
import Base: run
import ..SpatialModelIO: load_spatial_model
import ..SpatialSimulationRunner: run_spatial_model
import ..ResultIO: write_result, read_result

export Model, ElementRef, VariableRef, document, element!, marker!, graphic!,
       set_properties!, analysis!, simulation!, state_selection!,
       initial_conditions!, parameters!, graphics!, variable, load_model,
       simulate, run, write_model, save_result,
       ground!, rigid_body!, gravity!, applied_force!, applied_torque!,
       spanning_force!, spherical!, perp!, inplane!, inline!, hinge!, orient!,
       revolute!, fixed!, rotational_motion!, translational_motion!,
       spanning_motion!, span!, directed_distance!, bushing!, plane_contact!,
       surface_friction!, revolute_friction!, translational_friction!,
       inplane_friction!, rolling_tire!, coupler!, gear_pair!,
       rack_and_pinion!, pulley!, belt!, belt_span!, equation_component!

const RESERVED_TABLES = Set(("model", "parameters", "analysis", "simulation",
    "state_selection", "initial_conditions", "graphics"))
const VIEWER_RUNTIME = Ref{Any}(nothing)

"""Mutable Julia description of one spatial model."""
mutable struct Model
    document::Dict{String,Any}
    source_directory::String
end

"""Typed reference to a named element in a [`Model`](@ref)."""
struct ElementRef
    model::Model
    name::String
    kind::String
end

"""Reference to a named variable belonging to a model element."""
struct VariableRef
    model::Model
    name::String
end

Base.string(reference::ElementRef) = reference.name
Base.string(reference::VariableRef) = reference.name

function Base.show(io::IO, reference::ElementRef)
    print(io, "Sim3D.", reference.kind, "(\"", reference.name, "\")")
end

function Base.show(io::IO, reference::VariableRef)
    print(io, "Sim3D.variable(\"", reference.name, "\")")
end

function model_name(name)
    text = string(name)
    isempty(text) && throw(ArgumentError("model and element names cannot be empty"))
    any(isempty, split(text, '.')) && throw(ArgumentError(
        "model and element names cannot contain empty path segments"))
    text
end

"""
    Model(name; title=string(name), source_directory=pwd(), kwargs...)

Create an empty spatial model document. Additional keywords are stored in the
top-level `model` table.
"""
function Model(name; title = string(name), source_directory = pwd(), kwargs...)
    name_text = model_name(name)
    model_table = Dict{String,Any}(
        "name" => name_text,
        "title" => String(title),
        "dimension" => "spatial")
    for (key, value) in kwargs
        isnothing(value) || (model_table[string(key)] = model_value(nothing, value))
    end
    Model(Dict{String,Any}("model" => model_table),
        abspath(String(source_directory)))
end

same_model(model::Model, reference::Union{ElementRef,VariableRef}) =
    reference.model === model || throw(ArgumentError(
        "a Sim3D reference belongs to a different model"))

function model_value(model::Union{Nothing,Model}, reference::ElementRef)
    isnothing(model) || same_model(model, reference)
    reference.name
end

function model_value(model::Union{Nothing,Model}, reference::VariableRef)
    isnothing(model) || same_model(model, reference)
    reference.name
end

model_value(model, value::Symbol) = String(value)
model_value(model, value::Tuple) = model_value(model, collect(value))
model_value(model, value::AbstractVector) =
    Any[model_value(model, item) for item in value]
model_value(model, value::AbstractMatrix) =
    [Any[model_value(model, value[row, column]) for column in axes(value, 2)]
     for row in axes(value, 1)]
model_value(model, value::NamedTuple) =
    Dict{String,Any}(string(key) => model_value(model, item)
        for (key, item) in pairs(value))
model_value(model, value::AbstractDict) =
    Dict{String,Any}(string(key) => model_value(model, item)
        for (key, item) in value)
model_value(model, value) = value

function table_at!(document, path::AbstractVector{<:AbstractString}; create = true)
    table = document
    for segment in path
        if create
            child = get!(table, segment) do
                Dict{String,Any}()
            end
        else
            child = get(table, segment, nothing)
            isnothing(child) && return nothing
        end
        child isa AbstractDict || throw(ArgumentError(
            "'$segment' is not a model table"))
        table = child
    end
    table
end

function insert_table!(model::Model, qualified_name::String, table)
    segments = split(qualified_name, '.')
    first(segments) in RESERVED_TABLES && throw(ArgumentError(
        "'$qualified_name' conflicts with a reserved model table"))
    parent = table_at!(model.document, segments[1:end - 1])
    haskey(parent, last(segments)) && throw(ArgumentError(
        "model element '$qualified_name' is already defined"))
    parent[last(segments)] = table
    table
end

function update_table!(model::Model,
        path::AbstractVector{<:AbstractString}; kwargs...)
    table = table_at!(model.document, path)
    for (key, value) in kwargs
        isnothing(value) || (table[string(key)] = model_value(model, value))
    end
    model
end

"""
    element!(model, kind, name; kwargs...)

Add any supported spatial element. Typed convenience builders call this same
function; it is also the escape hatch for newly added element types.
"""
function element!(model::Model, kind, name; kwargs...)
    qualified_name = model_name(name)
    table = Dict{String,Any}("type" => string(kind))
    for (key, value) in kwargs
        isnothing(value) || (table[string(key)] = model_value(model, value))
    end
    insert_table!(model, qualified_name, table)
    ElementRef(model, qualified_name, string(kind))
end

"""Add a marker beneath a rigid body or ground element."""
function marker!(owner::ElementRef, name; kwargs...)
    owner.kind in ("ground", "rigid_body") || throw(ArgumentError(
        "markers must be owned by ground or a rigid body"))
    local_name = model_name(name)
    occursin('.', local_name) && throw(ArgumentError(
        "a marker's local name cannot contain a period"))
    element!(owner.model, "marker", "$(owner.name).$local_name"; kwargs...)
end

"""Return a qualified variable reference belonging to an element."""
function variable(owner::ElementRef, name)
    local_name = model_name(name)
    occursin('.', local_name) && throw(ArgumentError(
        "a variable's local name cannot contain a period"))
    VariableRef(owner.model, "$(owner.name).$local_name")
end

"""Merge fields into an existing element table."""
function set_properties!(reference::ElementRef; kwargs...)
    table = table_at!(reference.model.document, split(reference.name, '.');
        create = false)
    isnothing(table) && throw(ArgumentError(
        "model element '$(reference.name)' is not defined"))
    for (key, value) in kwargs
        isnothing(value) ||
            (table[string(key)] = model_value(reference.model, value))
    end
    reference
end

analysis!(model::Model; kwargs...) = update_table!(model, ["analysis"]; kwargs...)
state_selection!(model::Model; kwargs...) =
    update_table!(model, ["state_selection"]; kwargs...)
initial_conditions!(model::Model; kwargs...) =
    update_table!(model, ["initial_conditions"]; kwargs...)
parameters!(model::Model; kwargs...) = update_table!(model, ["parameters"]; kwargs...)
graphics!(model::Model; kwargs...) = update_table!(model, ["graphics"]; kwargs...)

"""
    simulation!(model; frames_per_second=nothing, kwargs...)

Set simulation fields. `frames_per_second` is a Julia-side convenience that
sets `output_samples` from the effective start and end times.
"""
function simulation!(model::Model; frames_per_second = nothing, kwargs...)
    update_table!(model, ["simulation"]; kwargs...)
    if !isnothing(frames_per_second)
        fps = Float64(frames_per_second)
        fps > 0 || throw(ArgumentError("frames_per_second must be positive"))
        table = model.document["simulation"]
        start_time = Float64(get(table, "start_time", 0.0))
        end_time = Float64(get(table, "end_time", 1.0))
        end_time >= start_time || throw(ArgumentError(
            "simulation end_time must not precede start_time"))
        table["output_samples"] = round(Int,
            (end_time - start_time) * fps) + 1
    end
    model
end

"""Attach one simple graphics table directly to an element."""
function graphics!(reference::ElementRef; kwargs...)
    update_table!(reference.model,
        [split(reference.name, '.'); "graphics"]; kwargs...)
    reference
end

"""Attach a named graphics child beneath an element's graphics table."""
function graphic!(reference::ElementRef, name; kwargs...)
    local_name = model_name(name)
    occursin('.', local_name) && throw(ArgumentError(
        "a graphic's local name cannot contain a period"))
    update_table!(reference.model,
        [split(reference.name, '.'); "graphics"; local_name]; kwargs...)
    reference
end

for (function_name, element_kind) in (
        (:ground!, "ground"),
        (:rigid_body!, "rigid_body"),
        (:gravity!, "gravity"),
        (:applied_force!, "applied_force"),
        (:applied_torque!, "applied_torque"),
        (:spanning_force!, "spanning_force"),
        (:spherical!, "spherical"),
        (:perp!, "perp"),
        (:inplane!, "inplane"),
        (:inline!, "inline"),
        (:hinge!, "hinge"),
        (:orient!, "orient"),
        (:revolute!, "revolute"),
        (:fixed!, "fixed"),
        (:rotational_motion!, "rotational_motion"),
        (:translational_motion!, "translational_motion"),
        (:spanning_motion!, "spanning_motion"),
        (:span!, "span"),
        (:directed_distance!, "directed_distance"),
        (:bushing!, "bushing"),
        (:plane_contact!, "plane_contact"),
        (:surface_friction!, "surface_friction"),
        (:revolute_friction!, "revolute_friction"),
        (:translational_friction!, "translational_friction"),
        (:inplane_friction!, "inplane_friction"),
        (:rolling_tire!, "rolling_tire"),
        (:coupler!, "coupler"),
        (:gear_pair!, "gear_pair"),
        (:rack_and_pinion!, "rack_and_pinion"),
        (:pulley!, "pulley"),
        (:belt!, "belt"),
        (:belt_span!, "belt_span"),
        (:equation_component!, "equation_component"))
    @eval begin
        function $(function_name)(model::Model, name; kwargs...)
            element!(model, $element_kind, name; kwargs...)
        end
    end
end

"""Return a detached copy of the hierarchical model document."""
document(model::Model) = deepcopy(model.document)

"""Validate and allocate a Julia-built model through the ordinary loader."""
function load_model(model::Model)
    load_spatial_model(document(model);
        source_directory = model.source_directory,
        source_label = "Julia Sim3D model '$(model.document["model"]["name"])'")
end

function load_spatial_model(model::Model;
        source_directory = model.source_directory,
        source_label = "Julia Sim3D model '$(model.document["model"]["name"])'")
    load_spatial_model(document(model);
        source_directory, source_label)
end

"""Run a Julia-built model and return the ordinary Sim3D analysis result."""
function run_spatial_model(model::Model; kwargs...)
    run_spatial_model(load_model(model); kwargs...)
end

simulate(model::Model; kwargs...) = run_spatial_model(model; kwargs...)
run(model::Model; kwargs...) = simulate(model; kwargs...)

"""Write the generated portable model document as TOML."""
function write_model(path::AbstractString, model::Model; overwrite = false)
    isfile(path) && !overwrite && throw(ArgumentError(
        "model file already exists: $path"))
    open(path, "w") do io
        TOML.print(io, document(model))
    end
    String(path)
end

function viewer_runtime()
    !isnothing(VIEWER_RUNTIME[]) && return VIEWER_RUNTIME[]
    runtime = Module(:Sim3DViewerRuntime, true, true)
    Core.eval(runtime, :(using PracticalMechanicalSimulation))
    viewer_directory = normpath(joinpath(@__DIR__, "..", "viewer"))
    Base.include(runtime, joinpath(viewer_directory, "ViewerData.jl"))
    Base.include(runtime, joinpath(viewer_directory, "StoredResultViewer.jl"))
    Base.include(runtime, joinpath(viewer_directory,
        "PortableViewerDocument.jl"))
    VIEWER_RUNTIME[] = runtime
    runtime
end

function attach_viewer_graphics!(path)
    runtime = viewer_runtime()
    stored_viewer = Base.invokelatest(getproperty,
        runtime, :StoredResultViewer)
    portable_viewer = Base.invokelatest(getproperty,
        runtime, :PortableViewerDocument)
    mechanism_result = Base.invokelatest(getproperty,
        stored_viewer, :stored_mechanism_result)
    write_graphics = Base.invokelatest(getproperty,
        portable_viewer, :write_graphics)
    stored = read_result(path)
    if stored.analysis_mode == :modal
        count = length(stored.modal_eigenvalues)
        results = [Base.invokelatest(mechanism_result, path;
            mode = number) for number in 1:count]
        labels = ["Mode $number — " * string(round(
            stored.modal_natural_frequencies_hz[number]; sigdigits = 6)) *
            " Hz" for number in 1:count]
        Base.invokelatest(write_graphics, path, results;
            labels, choice_name = "Mode")
    else
        result = Base.invokelatest(mechanism_result, path)
        Base.invokelatest(write_graphics, path, result;
            label = "Result")
    end
    path
end

"""
    save_result(path, result; viewer_graphics=true, kwargs...)

Write an ordinary `.simp` result and, by default, add the portable graphics
section read directly by SimpView. Remaining keywords are passed to
[`write_result`](@ref).
"""
function save_result(path::AbstractString, result;
        viewer_graphics = true, kwargs...)
    write_result(path, result; kwargs...)
    if viewer_graphics
        try
            attach_viewer_graphics!(path)
        catch error
            @warn("the result was saved, but SimpView graphics could not be prepared",
                exception = (error, catch_backtrace()))
        end
    end
    String(path)
end

end
