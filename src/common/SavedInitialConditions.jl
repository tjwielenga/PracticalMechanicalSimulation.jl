"""Shared name-and-type checked initialization from stored `.simp` results."""
module SavedInitialConditions

using TOML
using ..ResultIO

export disabled_saved_initial_conditions, apply_saved_initial_conditions!

disabled_saved_initial_conditions() =
    (; enabled = false, result_path = nothing, sample = nothing,
       requested_sample = nothing, source_time = nothing,
       include_velocities = false, transferred_variables = Symbol[])

function saved_result_sample(sample, stored)
    count = length(stored.times)
    count > 0 || throw(ArgumentError(
        "saved initial-condition result contains no samples"))
    if sample == "last"
        return count
    elseif sample == "first"
        return 1
    elseif sample == "static"
        isnothing(stored.static_initialization_iterations) &&
            throw(ArgumentError("initial_conditions.sample = 'static' " *
                "requires a dynamic result that used static initialization"))
        return 1
    elseif sample isa Integer && !(sample isa Bool)
        1 <= sample <= count || throw(ArgumentError(
            "initial_conditions.sample must be between 1 and $count"))
        return Int(sample)
    end
    throw(ArgumentError(
        "initial_conditions.sample must be 'first', 'last', 'static', " *
        "or a sample number"))
end

function collect_typed_tables!(result, table, path = String[])
    for (key, value) in table
        value isa AbstractDict || continue
        element_path = [path; key]
        haskey(value, "type") &&
            (result[Symbol(join(element_path, "."))] = value)
        collect_typed_tables!(result, value, element_path)
    end
    result
end

function saved_element_kinds(stored, dimension)
    isempty(stored.model_source) && throw(ArgumentError(
        "saved initial-condition result does not contain its source model"))
    document = try
        TOML.parse(stored.model_source)
    catch error
        throw(ArgumentError(
            "saved initial-condition result contains invalid model TOML: " *
            sprint(showerror, error)))
    end
    model = get(document, "model", Dict{String,Any}())
    get(model, "dimension", "") == dimension || throw(ArgumentError(
        "saved initial-condition result is not a $dimension model"))
    elements = collect_typed_tables!(Dict{Symbol,Any}(), document)
    Dict(name => string(get(table, "type", ""))
        for (name, table) in elements)
end

"""
    apply_saved_initial_conditions!(initial, layout, target_kinds, document,
        source_directory; dimension, configuration_kinds, velocity_kinds,
        element_types)

Transfer compatible configuration variables and optional velocities from a
stored result. Components, element types, canonical names, and variable kinds
must agree. The returned metadata records the resolved source and transferred
variables; accelerations, reactions, and work variables are never copied.
"""
function apply_saved_initial_conditions!(initial, layout, target_kinds,
        document, source_directory; dimension, configuration_kinds,
        velocity_kinds, element_types)
    haskey(document, "initial_conditions") ||
        return disabled_saved_initial_conditions()
    table = document["initial_conditions"]
    table isa AbstractDict || throw(ArgumentError(
        "initial_conditions must be a TOML table"))
    haskey(table, "result") || throw(ArgumentError(
        "initial_conditions.result is required"))
    result_name = table["result"]
    result_name isa AbstractString || throw(ArgumentError(
        "initial_conditions.result must be a filename"))
    include_velocities = get(table, "include_velocities", false)
    include_velocities isa Bool || throw(ArgumentError(
        "initial_conditions.include_velocities must be Boolean"))
    result_path = isabspath(result_name) ? normpath(result_name) :
        normpath(joinpath(source_directory, result_name))
    stored = read_result(result_path)
    sample = get(table, "sample", "last")
    sample_index = saved_result_sample(sample, stored)
    column_count = size(stored.values, 2)
    metadata_lengths = (length(stored.variable_names),
        length(stored.variable_components), length(stored.variable_kinds))
    all(==(column_count), metadata_lengths) || throw(ArgumentError(
        "saved initial-condition result has inconsistent variable metadata"))
    size(stored.values, 1) == length(stored.times) || throw(ArgumentError(
        "saved initial-condition result has inconsistent sample data"))
    source_kinds = saved_element_kinds(stored, dimension)
    source_columns = Dict{Tuple{Symbol,Symbol},Tuple{Symbol,Int}}()
    for index in eachindex(stored.variable_names)
        key = (Symbol(stored.variable_components[index]),
               Symbol(stored.variable_names[index]))
        haskey(source_columns, key) && throw(ArgumentError(
            "saved initial-condition result has duplicate variable " *
            "'$(key[1]).$(key[2])'"))
        source_columns[key] = (Symbol(stored.variable_kinds[index]), index)
    end
    transferred = Symbol[]
    for variable in layout.catalog.variables
        transferable = variable.kind in configuration_kinds ||
            (include_velocities && variable.kind in velocity_kinds)
        transferable || continue
        component = variable.component
        haskey(target_kinds, component) || continue
        string(target_kinds[component]) in element_types || continue
        haskey(source_kinds, component) || continue
        string(target_kinds[component]) == source_kinds[component] ||
            throw(ArgumentError("initial-condition component '$component' " *
                "has type '$(source_kinds[component])' in the result but " *
                "'$(target_kinds[component])' in the model"))
        key = (component, variable.name)
        haskey(source_columns, key) || continue
        source_kind, column = source_columns[key]
        same_user_state = source_kind in (:user_state_hold,
            :user_state_steady) && variable.kind in (:user_state_hold,
            :user_state_steady)
        (source_kind == variable.kind || same_user_state) || throw(ArgumentError(
            "initial-condition variable '$component.$(variable.name)' has " *
            "kind '$source_kind' in the result but '$(variable.kind)' " *
            "in the model"))
        value = stored.values[sample_index, column]
        isfinite(value) || throw(ArgumentError(
            "initial-condition variable '$component.$(variable.name)' " *
            "is not finite"))
        initial[variable.index] = value
        push!(transferred, Symbol(component, :., variable.name))
    end
    isempty(transferred) && throw(ArgumentError(
        "saved initial-condition result has no compatible configuration variables"))
    (; enabled = true, result_path, sample = sample_index,
       requested_sample = sample, source_time = stored.times[sample_index],
       include_velocities, transferred_variables = transferred)
end

end
