"""Shared hierarchical-document builder used by the Sim2D and Sim3D APIs."""
module JuliaModelBuilder

using TOML
import ..ResultIO: write_result, read_result

export Model, ElementRef, VariableRef, document, element!, marker!, graphic!,
       set_properties!, analysis!, simulation!, state_selection!,
       initial_conditions!, parameters!, graphics!, variable, write_model,
       save_result, normalized_document_value

const RESERVED_TABLES = Set(("model", "parameters", "analysis", "simulation",
    "state_selection", "initial_conditions", "graphics"))
const VIEWER_RUNTIME = Ref{Any}(nothing)

"""Mutable Julia description of a planar or spatial model."""
mutable struct Model{D}
    document::Dict{String,Any}
    source_directory::String
end

"""Typed reference to a named element in a [`Model`](@ref)."""
struct ElementRef{D}
    model::Model{D}
    name::String
    kind::String
end

"""Reference to a named variable belonging to a model element."""
struct VariableRef{D}
    model::Model{D}
    name::String
end

Base.string(reference::ElementRef) = reference.name
Base.string(reference::VariableRef) = reference.name

api_name(::Model{:planar}) = "Sim2D"
api_name(::Model{:spatial}) = "Sim3D"

function Base.show(io::IO, reference::ElementRef)
    print(io, api_name(reference.model), ".", reference.kind,
        "(\"", reference.name, "\")")
end

function Base.show(io::IO, reference::VariableRef)
    print(io, api_name(reference.model), ".variable(\"",
        reference.name, "\")")
end

function model_name(name)
    text = string(name)
    isempty(text) && throw(ArgumentError(
        "model and element names cannot be empty"))
    any(isempty, split(text, '.')) && throw(ArgumentError(
        "model and element names cannot contain empty path segments"))
    text
end

"""Create an empty model document for dimension `D`."""
function Model{D}(name; title = string(name), source_directory = pwd(),
        kwargs...) where {D}
    D in (:planar, :spatial) || throw(ArgumentError(
        "Julia model dimension must be :planar or :spatial"))
    name_text = model_name(name)
    model_table = Dict{String,Any}(
        "name" => name_text,
        "title" => String(title),
        "dimension" => String(D))
    for (key, value) in kwargs
        isnothing(value) ||
            (model_table[string(key)] = model_value(nothing, value))
    end
    Model{D}(Dict{String,Any}("model" => model_table),
        abspath(String(source_directory)))
end

same_model(model::Model, reference::Union{ElementRef,VariableRef}) =
    reference.model === model || throw(ArgumentError(
        "a $(api_name(model)) reference belongs to a different model"))

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

normalized_document_value(value) = model_value(nothing, value)

function table_at!(document, path::AbstractVector{<:AbstractString};
        create = true)
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

"""Add any element supported by the dimension-specific model loader."""
function element!(model::Model{D}, kind, name; kwargs...) where {D}
    qualified_name = model_name(name)
    table = Dict{String,Any}("type" => string(kind))
    for (key, value) in kwargs
        isnothing(value) || (table[string(key)] = model_value(model, value))
    end
    insert_table!(model, qualified_name, table)
    ElementRef{D}(model, qualified_name, string(kind))
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
function variable(owner::ElementRef{D}, name) where {D}
    local_name = model_name(name)
    occursin('.', local_name) && throw(ArgumentError(
        "a variable's local name cannot contain a period"))
    VariableRef{D}(owner.model, "$(owner.name).$local_name")
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

analysis!(model::Model; kwargs...) = update_table!(model, ["analysis"];
    kwargs...)
state_selection!(model::Model; kwargs...) =
    update_table!(model, ["state_selection"]; kwargs...)
initial_conditions!(model::Model; kwargs...) =
    update_table!(model, ["initial_conditions"]; kwargs...)
parameters!(model::Model; kwargs...) = update_table!(model, ["parameters"];
    kwargs...)
graphics!(model::Model; kwargs...) = update_table!(model, ["graphics"];
    kwargs...)

"""Set simulation fields, optionally deriving sample count from frame rate."""
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

"""Return a detached copy of the hierarchical model document."""
document(model::Model) = deepcopy(model.document)

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
    runtime = Module(:JuliaBuilderViewerRuntime, true, true)
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

"""Write a `.simp` result and optionally attach native SimpView graphics."""
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
