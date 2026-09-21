"""
    PlanarModelIO

Front end for hierarchical planar TOML models. The loader validates names and
references, declares and allocates the complete canonical system, constructs
markers and elements, corrects the initial configuration and velocity, detects
redundant ideal constraints, and selects independent physical velocities by
pivoted QR.

Automatic selection requires model derivatives, so loading is intentionally a
two-pass operation: a provisional assembly supplies the constraint partials;
the final assembly activates the chosen state equations while retaining the
dormant candidates needed for recovery during integration.
"""
module PlanarModelIO

using LinearAlgebra
using TOML
using ..InputUnits
using ..ScalarExpressions: compile_time_expression, compile_motion_expression,
    expression_uses_model_variables, constant_scalar_law
import ..ScalarExpressions
using ..AutomaticAnalysis
using ..PlanarAppliedForces
using ..PlanarDirectedDistances
using ..PlanarComponentAssembly
using ..PlanarEquationComponents
using ..PlanarFrictionForces
using ..PlanarModeling
using ..SavedInitialConditions
using ..JuliaModelBuilder: normalized_document_value
using ..AssemblyExpansion: expand_model_assemblies, parse_lua_model

export LoadedPlanarModel, PlanarModelMarker, load_planar_model,
       compile_time_expression

const RESERVED_TABLES = Set(("model", "parameters", "analysis", "simulation",
                             "state_selection", "initial_conditions",
                             "graphics"))

const CONFIGURATION_KINDS = Set((:position, :orientation, :relative_position,
                                 :user_state_hold, :user_state_steady))
const VELOCITY_KINDS = Set((:velocity, :angular_velocity, :relative_velocity))
const SAVED_STATE_ELEMENT_TYPES = Set(("rigid_body", "revolute",
                                       "distance_coordinate",
                                       "equation_component",
                                       "surface_friction",
                                       "revolute_friction",
                                       "translational_friction",
                                       "inplane_friction"))

const PLANAR_ELEMENT_TYPES = Set(("ground", "rigid_body", "marker",
    "floating_marker", "revolute", "inplane", "perp", "translational",
    "fixed", "gear_pair", "rack_and_pinion", "distance_coordinate",
    "coupler", "span", "pulley", "belt", "belt_span",
    "rotational_motion", "translational_motion", "gravity",
    "applied_force", "applied_torque", "bushing", "plane_contact",
    "spanning_force", "torsional_spring_damper", "surface_friction",
    "revolute_friction", "translational_friction", "inplane_friction",
    "equation_component"))

const EXPRESSION_KINEMATIC_KINDS = Set((
    :position, :orientation, :relative_position,
    :velocity, :angular_velocity, :relative_velocity))

expression_variable_supported(variable) =
    variable.kind in EXPRESSION_KINEMATIC_KINDS ||
    variable.kind in (:user_algebraic, :user_state_hold,
        :user_state_steady) ||
    (variable.kind == :applied_geometry && variable.name == :length) ||
    (variable.kind == :applied_rate && variable.name == :length_rate)

"""
Compile a restricted scalar expression that may use qualified kinematic model
variables. The returned law evaluates from `(t, z)` and obtains its local
gradient with forward-mode dual numbers.
"""
function compile_model_expression(source::AbstractString, parameters, layout)
    ScalarExpressions.compile_model_expression(source, parameters, layout,
        expression_variable_supported)
end

function linear_spanning_law(layout, name, stiffness, damping, free_length)
    variables = component_variable_indices(layout, name)
    dependencies = [variables[3], variables[6]]
    value = (t, z) -> -stiffness * (z[variables[3]] - free_length) -
        damping * z[variables[6]]
    gradient = (t, z) -> [-stiffness, -damping]
    PlanarScalarLaw(value, gradient, dependencies)
end

"""One model marker, represented internally by point and orientation views."""
struct PlanarModelMarker{B,P,O}
    name::Symbol
    owner::B
    point::P
    orientation::O
end

local_marker_angle(marker::PlanarBodyOrientationMarker) = marker.angle_offset
local_marker_angle(marker::PlanarGroundOrientationMarker) = marker.angle

function initialize_spanning_force!(initial, force, time)
    element = force.element
    marker_1 = PlanarAppliedForces.point_marker_kinematics(
        element.marker_1, initial)
    marker_2 = PlanarAppliedForces.point_marker_kinematics(
        element.marker_2, initial)
    spanning = marker_2.position - marker_1.position
    length = norm(spanning)
    length > 0 || throw(ArgumentError(
        "spanning force '$(force.name)' has zero initial length"))
    unit = spanning ./ length
    length_rate = dot(unit, marker_2.velocity - marker_1.velocity)
    initial[element.spanning_variables] .= spanning
    initial[element.length_variable] = length
    initial[element.unit_variables] .= unit
    initial[element.length_rate_variable] = length_rate
    initial[element.force_variable] = force.active[] ?
        force.law(time, initial) : 0.0
    isfinite(initial[element.force_variable]) || throw(ArgumentError(
        "spanning force '$(force.name)' expression is not finite initially"))
    initial[element.global_force_variables] .=
        -unit .* initial[element.force_variable]
    return initial
end

function initialize_span_measure!(initial, measure)
    marker_1 = PlanarAppliedForces.point_marker_kinematics(
        measure.marker_1, initial)
    marker_2 = PlanarAppliedForces.point_marker_kinematics(
        measure.marker_2, initial)
    spanning = marker_2.position - marker_1.position
    distance = norm(spanning)
    distance > 0 || throw(ArgumentError(
        "span '$(measure.name)' requires initially separated markers"))
    unit = spanning ./ distance
    velocity = dot(unit, marker_2.velocity - marker_1.velocity)
    initial[measure.spanning_variables] .= spanning
    initial[measure.length_variable] = distance
    initial[measure.unit_variables] .= unit
    initial[measure.length_rate_variable] = velocity
    values = PlanarComponentAssembly.planar_span_measure_values(
        measure, initial)
    initial[measure.length_acceleration_variable] =
        dot(values.unit, values.relative_acceleration) +
        values.transverse_speed_squared / values.length
    initial
end

function initialize_torsional_force!(initial, spring)
    element = spring.element
    relative = PlanarAppliedForces.relative_rotation(element, initial)
    initial[element.torque_variable] =
        -element.stiffness * (relative.angle - element.free_angle) -
        element.damping * relative.angular_velocity
    return initial
end

function select_belt_tangent(pulley_1, pulley_2, hint_1, hint_2, initial)
    center_1 = PlanarAppliedForces.point_marker_kinematics(
        pulley_1.center_marker, initial).position
    center_2 = PlanarAppliedForces.point_marker_kinematics(
        pulley_2.center_marker, initial).position
    hint_position_1 = PlanarAppliedForces.point_marker_kinematics(
        hint_1.point, initial).position
    hint_position_2 = PlanarAppliedForces.point_marker_kinematics(
        hint_2.point, initial).position
    candidates = NamedTuple[]
    for tangent_kind in (1, -1), tangent_side in (1, -1)
        geometry = try
            belt_tangent_geometry(center_1, pulley_1.pitch_radius,
                center_2, pulley_2.pitch_radius, tangent_kind, tangent_side)
        catch error
            error isa DomainError || rethrow()
            continue
        end
        score = sum(abs2, geometry.point_1 - hint_position_1) +
            sum(abs2, geometry.point_2 - hint_position_2)
        push!(candidates, (; tangent_kind, tangent_side, geometry, score))
    end
    isempty(candidates) && throw(ArgumentError(
        "no common tangent exists between pulleys '$(pulley_1.name)' and '$(pulley_2.name)'"))
    selected = candidates[argmin(getproperty.(candidates, :score))]
    length(candidates) > 1 && begin
        scores = sort(getproperty.(candidates, :score))
        isapprox(scores[1], scores[2]; rtol = 1.0e-10,
            atol = 1.0e-12) && throw(ArgumentError(
            "belt tangent hints are ambiguous between pulleys '$(pulley_1.name)' and '$(pulley_2.name)'"))
    end
    selected
end

function initialize_belt_span!(initial, span)
    values = belt_span_values(span, initial)
    initial[span.point_1_variables] .= values.actual_point_1
    initial[span.point_2_variables] .= values.actual_point_2
    initial[span.tangent_variables] .= values.actual_tangent
    initial[span.length_variable] = values.actual_length
    initial[span.extension_variable] = values.expected_extension
    initial[span.extension_rate_variable] = values.expected_extension_rate
    initial[span.tension_variable] = span.stiffness * values.expected_extension +
        span.damping * values.expected_extension_rate
    initial[span.force_variables] .=
        initial[span.tension_variable] .* values.actual_tangent
    initial
end

function belt_reference_path_length(specifications)
    straight_length = sum(spec.geometry.length for spec in specifications)
    wrap_length = 0.0
    for index in eachindex(specifications)
        incoming = specifications[index]
        outgoing = specifications[mod1(index + 1, length(specifications))]
        incoming.pulley_2 === outgoing.pulley_1 || throw(ArgumentError(
            "ordered belt spans must form a closed pulley loop"))
        radius = incoming.pulley_2.pitch_radius
        beta_in = incoming.geometry.beta_2
        beta_out = outgoing.geometry.beta_1
        sign(beta_in) == sign(beta_out) || throw(ArgumentError(
            "belt tangent hints produce incompatible wrap directions at pulley '$(incoming.pulley_2.name)'"))
        normal_in = incoming.geometry.normal_2
        normal_out = outgoing.geometry.normal_1
        angle_in = atan(normal_in[2], normal_in[1])
        angle_out = atan(normal_out[2], normal_out[1])
        circulation = sign(beta_in)
        wrap_angle = mod(circulation * (angle_out - angle_in), 2pi)
        wrap_length += radius * wrap_angle
    end
    straight_length + wrap_length
end

"""
    LoadedPlanarModel

Validated and allocated planar model returned by [`load_planar_model`](@ref).

The object retains the canonical layout and executable implicit model together
with named bodies, markers, connections, forces, and drivers.
`entered_initial_values` preserves the values supplied by the model before
consistency correction, while `initial_values` contains the corrected values;
both use canonical variable order. `active_variable_indices` and
`active_equation_indices` omit redundant ideal-constraint families and dormant
candidate state equations. `model_source` is the verbatim TOML used for later
result-file storage and model extraction.
"""
struct LoadedPlanarModel
    title::String
    layout::ModelLayout
    model::ExecutableAnalysisModel
    bodies::Dict{Symbol,Any}
    markers::Dict{Symbol,Any}
    connections::Dict{Symbol,Any}
    measures::Dict{Symbol,Any}
    forces::Dict{Symbol,Any}
    equation_components::Dict{Symbol,PlanarEquationComponent}
    drivers::Dict{Symbol,Any}
    analysis::NamedTuple
    simulation::NamedTuple
    state_selection::NamedTuple
    initial_conditions::NamedTuple
    initial_condition_weights::Dict{Symbol,NamedTuple}
    initial_variable_weights::Vector{Float64}
    imposed_initial_variable_indices::Vector{Int}
    entered_initial_values::Vector{Float64}
    initial_values::Vector{Float64}
    active_variable_indices::Vector{Int}
    active_equation_indices::Vector{Int}
    model_source::String
end

function parse_simulation(document)
    table = get(document, "simulation", Dict{String,Any}())
    start_time = Float64(get(table, "start_time", 0.0))
    end_time = Float64(get(table, "end_time", 5.0))
    output_samples = Int(get(table, "output_samples", 81))
    output_precision = Symbol(get(table, "output_precision", "single"))
    end_time >= start_time ||
        throw(ArgumentError("simulation.end_time must not precede start_time"))
    output_samples >= 1 ||
        throw(ArgumentError("simulation.output_samples must be positive"))
    output_precision in (:single, :double) || throw(ArgumentError(
        "simulation.output_precision must be single or double"))
    relative_tolerance = Float64(get(table, "relative_tolerance", 1.0e-7))
    absolute_tolerance = Float64(get(table, "absolute_tolerance", 1.0e-9))
    initial_step = Float64(get(table, "initial_step", 1.0e-6))
    maximum_step = Float64(get(table, "maximum_step", 0.005))
    all(>(0), (relative_tolerance, absolute_tolerance,
               initial_step, maximum_step)) ||
        throw(ArgumentError("simulation tolerances and step sizes must be positive"))
    (; start_time, end_time, output_samples, output_precision,
       relative_tolerance,
       absolute_tolerance, initial_step, maximum_step)
end

function parse_analysis(document)
    table = get(document, "analysis", Dict{String,Any}())
    haskey(table, "deficit") && throw(ArgumentError(
        "analysis.deficit is not a model option; the current StateSelected " *
        "formulation has deficit zero"))
    mode = Symbol(get(table, "mode", "automatic"))
    mode in (:automatic, :dynamic, :kinematic, :static, :modal) ||
        throw(ArgumentError(
            "analysis.mode must be automatic, dynamic, kinematic, static, or modal"))
    initialization = Symbol(get(table, "initialization", "none"))
    initialization in (:none, :static_equilibrium) || throw(ArgumentError(
        "analysis.initialization must be none or static_equilibrium"))
    initialization == :static_equilibrium && mode in (:kinematic, :static) &&
        throw(ArgumentError("static_equilibrium initialization requires " *
            "dynamic or automatic analysis"))
    static_method = Symbol(get(table, "static_method", "newton"))
    static_method in (:newton, :dynamic_relaxation) || throw(ArgumentError(
        "analysis.static_method must be newton or dynamic_relaxation"))
    relaxation_duration = Float64(get(table, "relaxation_duration", 0.5))
    relaxation_duration > 0 || throw(ArgumentError(
        "analysis.relaxation_duration must be positive"))
    relaxation_reduction_factor =
        Float64(get(table, "relaxation_reduction_factor", 0.25))
    0 <= relaxation_reduction_factor < 1 || throw(ArgumentError(
        "analysis.relaxation_reduction_factor must be at least zero and less than one"))
    relaxation_cycles = Int(get(table, "relaxation_cycles", 8))
    relaxation_cycles >= 1 || throw(ArgumentError(
        "analysis.relaxation_cycles must be positive"))
    number_of_modes = Int(get(table, "modes", 10))
    number_of_modes >= 1 || throw(ArgumentError(
        "analysis.modes must be positive"))
    frequency_shift_hz = Float64(get(table, "frequency_shift_hz", 0.0))
    frequency_shift_hz >= 0 || throw(ArgumentError(
        "analysis.frequency_shift_hz must be nonnegative"))
    modal_tolerance = Float64(get(table, "modal_tolerance", 1.0e-9))
    modal_tolerance > 0 || throw(ArgumentError(
        "analysis.modal_tolerance must be positive"))
    (; mode, initialization, static_method, relaxation_duration,
       relaxation_reduction_factor, relaxation_cycles,
       number_of_modes, frequency_shift_hz, modal_tolerance,
       formulation = :StateSelected, deficit = 0)
end

function collect_elements!(result, table, path = String[])
    for (key, value) in table
        key in RESERVED_TABLES && isempty(path) && continue
        value isa AbstractDict || continue
        element_path = [path; key]
        if haskey(value, "type")
            name = Symbol(join(element_path, "."))
            haskey(result, name) && throw(ArgumentError("duplicate element '$name'"))
            result[name] = value
        end
        collect_elements!(result, value, element_path)
    end
    result
end

function require_vector(table, field, n; default = nothing, label = "element")
    value = get(table, field, default)
    value isa Vector && length(value) == n ||
        throw(ArgumentError("$label.$field must contain $n values"))
    all(x -> x isa Number, value) || throw(ArgumentError("$label.$field must be numeric"))
    Float64.(value)
end

function planar_orientation(table, label)
    has_orientation = haskey(table, "orientation")
    has_angle = haskey(table, "angle")
    has_orientation && has_angle && throw(ArgumentError(
        "$label must not specify both orientation and angle"))
    value = has_orientation ? table["orientation"] : get(table, "angle", 0.0)
    angle_value(value, "$label.orientation")
end

function required(collection, name, kind)
    key = Symbol(name)
    haskey(collection, key) || throw(ArgumentError("unknown $kind '$name'"))
    collection[key]
end

function parent_names(name::Symbol)
    parts = split(String(name), ".")
    [Symbol(join(parts[1:k], ".")) for k in (length(parts) - 1):-1:1]
end

function marker_owner(name, elements, bodies)
    for parent in parent_names(name)
        haskey(bodies, parent) && return (:body, parent)
        haskey(elements, parent) && get(elements[parent], "type", "") == "ground" &&
            return (:ground, parent)
    end
    throw(ArgumentError("marker '$name' must be nested under a rigid_body or ground"))
end

function parse_state_selection(document, valid_velocity_names)
    table = get(document, "state_selection", Dict{String,Any}())
    method = Symbol(get(table, "method", "automatic"))
    method in (:automatic, :preferred) ||
        throw(ArgumentError("state_selection.method must be 'automatic' or 'preferred'"))
    preferred = Symbol.(get(table, "preferred_velocities", String[]))
    valid = Set(valid_velocity_names)
    all(name -> name in valid, preferred) ||
        throw(ArgumentError("state_selection contains an unknown preferred velocity"))
    method == :preferred && isempty(preferred) &&
        throw(ArgumentError("preferred state selection requires preferred_velocities"))
    maximum_condition_number = Float64(get(table,
        "maximum_condition_number", inv(sqrt(eps(Float64)))))
    maximum_condition_number > 1 || throw(ArgumentError(
        "state_selection.maximum_condition_number must be greater than one"))
    requested = Symbol.(get(table, "_requested_preferred_velocities",
                            String.(preferred)))
    (; method, preferred_velocities = requested,
       selected_velocities = preferred,
       allow_fallback = Bool(get(table, "allow_fallback", true)),
       maximum_condition_number,
       fallback_used = Bool(get(table, "_fallback_used", false)),
       dependent_condition = haskey(table, "_dependent_condition") ?
            Float64(table["_dependent_condition"]) : nothing,
       qr_rank = haskey(table, "_qr_rank") ? Int(table["_qr_rank"]) : nothing,
       qr_tolerance = haskey(table, "_qr_tolerance") ?
            Float64(table["_qr_tolerance"]) : nothing,
       qr_pivots = Symbol.(get(table, "_qr_pivots", String[])),
       body_characteristic_lengths = Dict(Symbol(name) => Float64(value)
            for (name, value) in get(table, "_body_characteristic_lengths",
                                     Dict{String,Any}())),
       qr_row_scales = Float64.(get(table, "_qr_row_scales", Float64[])),
       qr_column_scales = Float64.(get(table, "_qr_column_scales", Float64[])),
       position_corrections = Int(get(table, "_position_corrections", 0)),
       velocity_corrections = Int(get(table, "_velocity_corrections", 0)),
       imposed_initial_variables = Symbol.(get(table,
            "_imposed_initial_variables", String[])),
       redundant_velocity_equations =
            Symbol.(get(table, "_redundant_velocity_equations", String[])))
end

function body_characteristic_lengths(bodies, body_names, markers, overrides)
    marker_radii = Dict(name => Float64[] for name in body_names)
    for marker in values(markers)
        isnothing(marker.owner) && continue
        marker.point isa PlanarBodyPointMarker || continue
        push!(marker_radii[marker.owner.name], norm(marker.point.r_body))
    end
    natural = Dict{Symbol,Union{Nothing,Float64}}()
    for name in body_names
        radii = marker_radii[name]
        marker_length = isempty(radii) ? 0.0 : maximum(radii)
        body = bodies[name]
        gyration_length = body.mass > 0 && body.inertia > 0 ?
            sqrt(body.inertia / body.mass) : nothing
        natural[name] = marker_length > 0 ? marker_length : gyration_length
    end
    positive = Float64[value for value in values(natural)
                       if !isnothing(value) && value > 0]
    model_length = isempty(positive) ? 1.0 : maximum(positive)
    minimum_length = sqrt(eps(Float64)) * model_length
    Dict(name => haskey(overrides, name) ? overrides[name] :
        (isnothing(natural[name]) ? model_length :
         max(natural[name], minimum_length)) for name in body_names)
end

"""
Select independent planar velocities from the scaled constraint partial.

Body translational columns are scaled by characteristic length so they can be
compared with angular rates. Column-pivoted QR chooses state velocities. A
second QR of the transpose finds redundant equation rows independently of that
choice. Both results are retained for preferred-state checks and diagnostics.
"""
function automatic_velocity_selection(model, bodies, body_names, markers,
        connections, connection_names, characteristic_length_overrides,
        initial, time)
    catalog = model.catalog
    rows = [equation.index for equation in catalog.equations
            if equation.level == 1 &&
               equation.kind in (:constraint, :coordinate_relation, :motion)]
    columns = Int[]
    names = Symbol[]
    lengths = body_characteristic_lengths(bodies, body_names, markers,
        characteristic_length_overrides)
    column_scales = Float64[]
    for body_name in body_names
        body = bodies[body_name]
        append!(columns, body.velocity_variables)
        push!(columns, body.angular_velocity_variable)
        append!(names, (Symbol(body_name, :., :V_x),
                        Symbol(body_name, :., :V_y),
                        Symbol(body_name, :., :omega)))
        append!(column_scales, (lengths[body_name], lengths[body_name], 1.0))
    end
    for connection_name in connection_names
        connection = connections[connection_name]
        if connection isa PlanarRevoluteJointComponent &&
                !isempty(connection.rotation_variables)
            push!(columns, connection.rotation_variables[2])
            push!(names, Symbol(connection_name, :., :omega))
            push!(column_scales, 1.0)
        elseif connection isa PlanarDistanceCoordinateComponent
            push!(columns, connection.velocity_variable)
            push!(names, Symbol(connection_name, :., :velocity))
            push!(column_scales, maximum(values(lengths)))
        end
    end
    selection = AnalysisSelection(VelocityIC(), columns, rows)
    D = evaluate_analysis_jacobian(model, selection, time, initial,
                                   zeros(length(initial)), 0.0)
    column_scaled = D * Diagonal(column_scales)
    row_norms = [norm(row) for row in eachrow(column_scaled)]
    row_scales = [value > 0 ? inv(value) : 1.0 for value in row_norms]
    scaled_D = Diagonal(row_scales) * column_scaled
    if isempty(rows)
        return copy(names), (; rank = 0, tolerance = 0.0,
            pivots = copy(names), independent_rows = Int[],
            redundant_rows = Int[], D, scaled_D, column_names = names,
            body_lengths = lengths, row_scales, column_scales)
    end
    factorization = qr(scaled_D, ColumnNorm())
    pivots = collect(factorization.p)
    diagonal = abs.(diag(factorization.R))
    tolerance = isempty(diagonal) ? 0.0 :
        max(size(scaled_D)...) * eps(Float64) * maximum(diagonal)
    numerical_rank = count(>(tolerance), diagonal)
    independent = pivots[(numerical_rank + 1):end]
    dependent = pivots[1:numerical_rank]
    row_factorization = qr(transpose(scaled_D[:, dependent]), ColumnNorm())
    row_pivots = collect(row_factorization.p)
    independent_rows = row_pivots[1:numerical_rank]
    redundant_rows = row_pivots[(numerical_rank + 1):end]
    return names[independent], (; rank = numerical_rank, tolerance,
        pivots = names[pivots], independent_rows = rows[independent_rows],
        redundant_rows = rows[redundant_rows], D,
        scaled_D, column_names = names, body_lengths = lengths,
        row_scales, column_scales)
end

function check_preferred_velocities(preferred, diagnostics,
        maximum_condition_number)
    length(unique(preferred)) == length(preferred) ||
        throw(ArgumentError("preferred velocities must be unique"))
    independent = Int[]
    for name in preferred
        index = findfirst(==(name), diagnostics.column_names)
        isnothing(index) && throw(ArgumentError(
            "unknown preferred velocity '$name'"))
        push!(independent, index)
    end
    expected = size(diagnostics.scaled_D, 2) - diagnostics.rank
    length(independent) == expected || return (acceptable = false,
        condition = Inf, reason = "expected $expected preferred velocities, got $(length(independent))")
    dependent = setdiff(collect(axes(diagnostics.scaled_D, 2)), independent)
    block = diagnostics.scaled_D[:, dependent]
    values = svdvals(block)
    block_rank = count(>(diagnostics.tolerance), values)
    condition = isempty(values) ? 1.0 :
        (last(values) == 0 ? Inf : first(values) / last(values))
    block_rank == diagnostics.rank || return (acceptable = false,
        condition, reason = "the complementary dependent block is rank deficient")
    condition <= maximum_condition_number || return (acceptable = false,
        condition, reason = "the complementary dependent block has condition number $condition")
    (; acceptable = true, condition, reason = "")
end

function constraint_family_indices(layout, velocity_equation_indices)
    inactive_variables = Int[]
    inactive_equations = Int[]
    for velocity_index in velocity_equation_indices
        velocity = layout.catalog.equations[velocity_index]
        velocity.kind == :constraint || continue
        matches = [family for family in values(layout.constraint_families)
                   if family.velocity_equation == velocity_index]
        length(matches) == 1 || error(
            "Constraint velocity equation $(velocity.component).$(velocity.name) " *
            "must belong to exactly one scalar constraint family")
        family = only(matches)
        append!(inactive_equations, (family.position_equation,
            family.velocity_equation, family.acceleration_equation))
        push!(inactive_variables, family.reaction_variable)
    end
    sort!(unique!(inactive_variables)), sort!(unique!(inactive_equations))
end

function weighted_initial_projection!(initial, model, variable_kinds,
        equation_level, equation_kinds, variable_weights, imposed_variables,
        time, label; tolerance = 1.0e-11, maximum_iterations = 30)
    catalog = model.catalog
    variables = [variable.index for variable in catalog.variables
        if variable.kind in variable_kinds &&
           variable.index ∉ imposed_variables]
    equations = [equation.index for equation in catalog.equations
        if equation.level == equation_level &&
           equation.kind in equation_kinds]
    isempty(equations) && return 0
    inverse_sqrt_weights = 1 ./ sqrt.(variable_weights[variables])
    selection = AnalysisSelection(Dynamics(), variables, equations)
    derivative = zeros(length(initial))
    residual = zeros(length(equations))
    imposed_names = [string(variable.component, ".", variable.name)
        for variable in catalog.variables
        if variable.index in imposed_variables &&
           variable.kind in variable_kinds]
    imposed_detail = isempty(imposed_names) ? "" :
        "; imposed values: " * join(imposed_names, ", ")
    for iteration in 0:maximum_iterations
        evaluate_analysis_equations!(residual, model, selection,
            time, initial, derivative)
        norm(residual, Inf) <= tolerance && return iteration
        isempty(variables) && throw(ArgumentError(
            "$label cannot satisfy its equations because all available " *
            "variables are imposed$imposed_detail"))
        iteration == maximum_iterations && break
        jacobian = evaluate_analysis_jacobian(model, selection,
            time, initial, derivative, 0.0)
        row_qr = qr(transpose(jacobian), ColumnNorm())
        diagonal = abs.(diag(row_qr.R))
        rank_tolerance = isempty(diagonal) ? 0.0 :
            max(size(jacobian)...) * eps(Float64) * maximum(diagonal)
        row_rank = count(>(rank_tolerance), diagonal)
        active_rows = collect(row_qr.p)[1:row_rank]
        active_jacobian = jacobian[active_rows, :]
        active_residual = residual[active_rows]
        isempty(active_rows) && break
        scaled = active_jacobian * Diagonal(inverse_sqrt_weights)
        correction = inverse_sqrt_weights .* (scaled' *
            ((scaled * scaled') \ active_residual))
        initial[variables] .-= correction
    end
    throw(ArgumentError("$label failed to converge$imposed_detail"))
end

"""Project the supplied configuration onto all level-zero model equations."""
function correct_initial_configuration!(initial, model, variable_weights,
        imposed_variables, time)
    weighted_initial_projection!(initial, model, CONFIGURATION_KINDS, 0,
        Set((:constraint, :normalization, :coordinate_relation, :motion)),
        variable_weights, imposed_variables, time,
        "weighted initial-configuration correction")
end

"""Project the supplied velocities onto all level-one model equations."""
function correct_initial_velocities!(initial, model, variable_weights,
        imposed_variables, time)
    weighted_initial_projection!(initial, model, VELOCITY_KINDS, 1,
        Set((:constraint, :coordinate_relation, :motion)),
        variable_weights, imposed_variables, time,
        "weighted initial-velocity correction")
end

function body_ic_weights(name, table)
    old_fields = ("position_ic_weights", "orientation_ic_weight",
        "velocity_ic_weights", "angular_velocity_ic_weight")
    old = findfirst(field -> haskey(table, field), old_fields)
    isnothing(old) || throw(ArgumentError(
        "body '$name'.$(old_fields[old]) has been replaced by " *
        "the scalar ic_weight_scale"))
    scale = finite_number(get(table, "ic_weight_scale", 1.0),
        "body '$name'.ic_weight_scale")
    scale > 0 || throw(ArgumentError(
        "body '$name'.ic_weight_scale must be positive"))
    (; scale)
end

function initial_table(table, label)
    value = get(table, "initial", Dict{String,Any}())
    value isa AbstractDict || throw(ArgumentError(
        "$label.initial must be a table"))
    value
end

function impose_table(table, label)
    initial = initial_table(table, label)
    value = get(initial, "impose", Dict{String,Any}())
    value isa AbstractDict || throw(ArgumentError(
        "$label.initial.impose must be a table"))
    value
end

function finite_number(value, label)
    value isa Number || throw(ArgumentError("$label must be numeric"))
    result = Float64(value)
    isfinite(result) || throw(ArgumentError("$label must be finite"))
    result
end

const PLANAR_ANALYSIS_STAGES = (:static, :dynamic, :modal)

function active_during_stages(table, label)
    has_active = haskey(table, "active_during")
    has_inactive = haskey(table, "inactive_during")
    has_active && has_inactive && throw(ArgumentError(
        "$label may specify active_during or inactive_during, but not both"))
    !has_active && !has_inactive && return PLANAR_ANALYSIS_STAGES

    field = has_active ? "active_during" : "inactive_during"
    specification = table[field]
    names = if specification isa AbstractString
        [String(specification)]
    elseif specification isa Vector &&
            all(value -> value isa AbstractString, specification)
        String.(specification)
    else
        throw(ArgumentError("$label.$field must be a stage name or " *
            "an array of stage names"))
    end
    isempty(names) && throw(ArgumentError(
        "$label.$field must contain at least one stage"))
    stages = Symbol.(lowercase.(strip.(names)))
    if :always in stages
        has_active || throw(ArgumentError(
            "$label.inactive_during cannot use always"))
        length(stages) == 1 || throw(ArgumentError(
            "$label.active_during cannot combine always with other stages"))
        return PLANAR_ANALYSIS_STAGES
    end
    all(stage -> stage in PLANAR_ANALYSIS_STAGES, stages) ||
        throw(ArgumentError("$label.$field stages must be static, " *
            "dynamic, or modal"))
    length(unique(stages)) == length(stages) || throw(ArgumentError(
        "$label.$field must not repeat a stage"))
    has_active ? Tuple(stages) :
        Tuple(stage for stage in PLANAR_ANALYSIS_STAGES if stage ∉ stages)
end

function positive_weight(value, label)
    result = finite_number(value, label)
    result > 0 || throw(ArgumentError("$label must be positive"))
    result
end

function friction_parameters(table, label; angular = false,
        effective_radius = false)
    stiffness = finite_number(get(table, "stiffness", NaN),
        "$label.stiffness")
    damping = finite_number(get(table, "damping", 0.0),
        "$label.damping")
    preload = finite_number(get(table, "preload", 0.0),
        "$label.preload")
    static_coefficient = finite_number(
        get(table, "static_coefficient", NaN),
        "$label.static_coefficient")
    dynamic_coefficient = finite_number(
        get(table, "dynamic_coefficient", NaN),
        "$label.dynamic_coefficient")
    transition_speed = finite_number(get(table, "transition_speed",
        angular ? 0.1 : 0.01), "$label.transition_speed")
    release_time = finite_number(get(table, "release_time", 0.01),
        "$label.release_time")
    stiffness > 0 && damping >= 0 && preload >= 0 &&
        static_coefficient >= dynamic_coefficient >= 0 &&
        transition_speed > 0 && release_time > 0 || throw(ArgumentError(
            "$label requires positive stiffness, nonnegative damping and " *
            "preload, static_coefficient >= dynamic_coefficient >= 0, " *
            "and positive transition_speed and release_time"))
    radius = if effective_radius
        value = finite_number(get(table, "effective_radius", NaN),
            "$label.effective_radius")
        value > 0 || throw(ArgumentError(
            "$label.effective_radius must be positive"))
        value
    else
        nothing
    end
    (; stiffness, damping, preload, static_coefficient, dynamic_coefficient,
       transition_speed, release_time, effective_radius = radius)
end

function apply_body_initial_impose!(initial, imposed_variables, body, table,
        name)
    initial_spec = initial_table(table, "body '$name'")
    unknown_initial = setdiff(Set(keys(initial_spec)), Set(("impose",)))
    isempty(unknown_initial) || throw(ArgumentError(
        "body '$name'.initial has unknown field '$(first(unknown_initial))'"))
    imposed = impose_table(table, "body '$name'")
    indices = Dict(
        "R_x" => body.position_variables[1],
        "R_y" => body.position_variables[2],
        "theta" => body.orientation_variable,
        "V_x" => body.velocity_variables[1],
        "V_y" => body.velocity_variables[2],
        "omega" => body.angular_velocity_variable)
    for (field, value) in imposed
        haskey(indices, field) || throw(ArgumentError(
            "body '$name'.initial.impose has unknown state '$field'"))
        index = indices[field]
        initial[index] = field == "theta" ?
            angle_value(value, "body '$name'.initial.impose.theta") :
            finite_number(value, "body '$name'.initial.impose.$field")
        push!(imposed_variables, index)
    end
    nothing
end

function relative_initial_specification(table, label, kind)
    initial = initial_table(table, label)
    coordinate, velocity = kind == :rotation ?
        ("angle", "omega") : ("distance", "velocity")
    allowed = Set((coordinate, velocity,
        coordinate * "_weight", velocity * "_weight", "impose"))
    unknown = setdiff(Set(keys(initial)), allowed)
    isempty(unknown) || throw(ArgumentError(
        "$label.initial has unknown field '$(first(unknown))'"))
    imposed = impose_table(table, label)
    imposed_allowed = Set((coordinate, velocity))
    unknown_imposed = setdiff(Set(keys(imposed)), imposed_allowed)
    isempty(unknown_imposed) || throw(ArgumentError(
        "$label.initial.impose has unknown state '$(first(unknown_imposed))'"))
    for field in keys(imposed)
        haskey(initial, field) && throw(ArgumentError(
            "$label.initial.$field and initial.impose.$field must not both be specified"))
    end
    coordinate_value(value, field_label) = kind == :rotation ?
        angle_value(value, field_label) : finite_number(value, field_label)
    coordinate_guess = haskey(initial, coordinate) ?
        coordinate_value(initial[coordinate], "$label.initial.$coordinate") :
        nothing
    velocity_guess = haskey(initial, velocity) ?
        finite_number(initial[velocity], "$label.initial.$velocity") : nothing
    coordinate_imposed = haskey(imposed, coordinate) ?
        coordinate_value(imposed[coordinate],
            "$label.initial.impose.$coordinate") : nothing
    velocity_imposed = haskey(imposed, velocity) ?
        finite_number(imposed[velocity],
            "$label.initial.impose.$velocity") : nothing
    coordinate_weight = positive_weight(get(initial,
        coordinate * "_weight", 1.0),
        "$label.initial.$(coordinate)_weight")
    velocity_weight = positive_weight(get(initial,
        velocity * "_weight", 1.0),
        "$label.initial.$(velocity)_weight")
    (; coordinate, velocity, coordinate_guess, velocity_guess,
       coordinate_imposed, velocity_imposed, coordinate_weight,
       velocity_weight)
end

function apply_relative_initial_specification!(initial, variable_weights,
        imposed_variables, coordinate_index, velocity_index, specification)
    if !isnothing(specification.coordinate_guess)
        initial[coordinate_index] = specification.coordinate_guess
    end
    if !isnothing(specification.velocity_guess)
        initial[velocity_index] = specification.velocity_guess
    end
    variable_weights[coordinate_index] = specification.coordinate_weight
    variable_weights[velocity_index] = specification.velocity_weight
    if !isnothing(specification.coordinate_imposed)
        initial[coordinate_index] = specification.coordinate_imposed
        push!(imposed_variables, coordinate_index)
    end
    if !isnothing(specification.velocity_imposed)
        initial[velocity_index] = specification.velocity_imposed
        push!(imposed_variables, velocity_index)
    end
    nothing
end

function driver_laws(table, parameters, name)
    function_name = get(table, "function", "expression")
    if function_name == "constant_speed"
        angle_0 = angle_value(get(table, "initial_angle", 0.0),
            "driver '$name'.initial_angle")
        omega_specification = table["angular_velocity"]
        omega = omega_specification isa Number ? Float64(omega_specification) :
            compile_time_expression(string(omega_specification), parameters)(0.0)
        return (t -> atan(sin(angle_0 + omega*t), cos(angle_0 + omega*t)),
                t -> omega, t -> 0.0)
    elseif function_name == "expression"
        haskey(table, "angle") || throw(ArgumentError(
            "driver '$name' requires angle"))
        for field in ("angular_velocity", "angular_acceleration")
            haskey(table, field) && throw(ArgumentError(
                "driver '$name'.$field is differentiated automatically " *
                "from angle and must be omitted"))
        end
        return compile_motion_expression(string(table["angle"]), parameters)
    end
    throw(ArgumentError("unsupported rotational-motion function '$function_name'"))
end

function translational_driver_laws(table, parameters, name)
    function_name = get(table, "function", "expression")
    if function_name == "constant_speed"
        distance_0 = Float64(get(table, "initial_distance", 0.0))
        haskey(table, "velocity") || throw(ArgumentError(
            "translational driver '$name' requires velocity"))
        velocity_specification = table["velocity"]
        velocity = velocity_specification isa Number ?
            Float64(velocity_specification) :
            compile_time_expression(
                string(velocity_specification), parameters)(0.0)
        return (t -> distance_0 + velocity * t,
                t -> velocity, t -> 0.0)
    elseif function_name == "expression"
        haskey(table, "distance") || throw(ArgumentError(
            "translational driver '$name' requires distance"))
        for field in ("velocity", "acceleration")
            haskey(table, field) && throw(ArgumentError(
                "translational driver '$name'.$field is differentiated " *
                "automatically from distance and must be omitted"))
        end
        return compile_motion_expression(
            string(table["distance"]), parameters)
    end
    throw(ArgumentError(
        "unsupported translational-motion function '$function_name'"))
end

"""
    load_planar_model(path::AbstractString) -> LoadedPlanarModel
    load_planar_model(input::IO; format=:toml, source_directory=pwd()) -> LoadedPlanarModel
    load_planar_model(document::AbstractDict; source_directory=pwd()) -> LoadedPlanarModel

Parse, validate, allocate, and assemble a planar TOML or Lua model.

Before selecting states, the loader corrects the supplied configuration and
velocity onto the position- and velocity-level equations while respecting IC
weights and imposed values. It then uses the scaled velocity-constraint partial
matrix to select independent physical velocities and suppress any redundant
scalar ideal-constraint families. The input stream is consumed but not closed;
the path form opens and closes its own file.
"""
function load_planar_text(source::AbstractString, format::Symbol,
        source_directory::AbstractString, source_label::AbstractString)
    format in (:toml, :lua) || throw(ArgumentError(
        "planar model format must be :toml or :lua"))
    document = if format == :lua
        parse_lua_model(source;
            label = "Lua model '$source_label'",
            source_directory,
            element_types = PLANAR_ELEMENT_TYPES)
    else
        TOML.parse(source)
    end
    model_table = get(document, "model", Dict{String,Any}())
    get(model_table, "dimension", "") == "planar" ||
        throw(ArgumentError("model.dimension must be 'planar'"))
    has_assembly_imports = !isempty(get(model_table, "assemblies", String[]))
    document = expand_model_assemblies(document, source_directory;
        dimension = "planar")
    stored_source = format == :lua || has_assembly_imports ?
        sprint(io -> TOML.print(io, document)) : String(source)
    load_planar_document(document, stored_source;
        source_directory)
end

function load_planar_model(path::AbstractString; format = nothing,
        source_directory = nothing, source_label = nothing)
    full_path = abspath(path)
    source = read(full_path, String)
    detected_format = endswith(lowercase(full_path), ".lua") ? :lua : :toml
    selected_format = isnothing(format) ? detected_format : Symbol(format)
    directory = isnothing(source_directory) ? dirname(full_path) :
        abspath(String(source_directory))
    label = isnothing(source_label) ? path : String(source_label)
    load_planar_text(source, selected_format, directory, label)
end

function load_planar_model(input::IO; format = :toml,
        source_directory = pwd(), source_label = "input")
    source = read(input, String)
    load_planar_text(source, Symbol(format),
        abspath(String(source_directory)), String(source_label))
end

function load_planar_model(document::AbstractDict; source_directory = pwd())
    normalized = normalized_document_value(document)
    source = sprint(io -> TOML.print(io, normalized))
    load_planar_document(normalized, source;
        source_directory = abspath(String(source_directory)))
end

function load_planar_document(document, source = ""; initial_override = nothing,
        entered_initial_override = nothing,
        initial_conditions_override = nothing, source_directory = pwd())
    model_table = get(document, "model", Dict{String,Any}())
    get(model_table, "dimension", "") == "planar" ||
        throw(ArgumentError("model.dimension must be 'planar'"))
    title = string(get(model_table, "title", get(model_table, "name", "Planar model")))
    parameters = get(document, "parameters", Dict{String,Any}())
    analysis = parse_analysis(document)
    simulation = parse_simulation(document)
    elements = collect_elements!(Dict{Symbol,Any}(), document)
    kinds = Dict(name => get(table, "type", "") for (name, table) in elements)
    names_of(types) = sort!([name for (name, kind) in kinds if kind in types])
    body_names = names_of(("rigid_body",))
    fixed_marker_names = names_of(("marker",))
    floating_marker_names = names_of(("floating_marker",))
    marker_names = sort!([fixed_marker_names; floating_marker_names])
    connection_names = names_of((
        "revolute", "inplane", "perp", "translational", "fixed",
        "gear_pair", "rack_and_pinion", "distance_coordinate", "coupler"))
    distance_coordinate_names = names_of(("distance_coordinate",))
    coupler_names = names_of(("coupler",))
    pulley_names = names_of(("pulley",))
    belt_names = names_of(("belt",))
    belt_span_names = names_of(("belt_span",))
    driver_names = names_of(("rotational_motion", "translational_motion"))
    bushing_names = names_of(("bushing",))
    contact_names = names_of(("plane_contact",))
    surface_friction_names = names_of(("surface_friction",))
    spanning_names = names_of(("spanning_force",))
    span_measure_names = names_of(("span",))
    torsional_names = names_of(("torsional_spring_damper",))
    revolute_friction_names = names_of(("revolute_friction",))
    translational_friction_names = names_of(("translational_friction",))
    inplane_friction_names = names_of(("inplane_friction",))
    friction_names = sort!([surface_friction_names; revolute_friction_names;
        translational_friction_names; inplane_friction_names])
    equation_component_names = names_of(("equation_component",))
    stateful_applied_force_names = [name for name in names_of(("applied_force",))
        if (haskey(elements[name], "expression") && expression_uses_model_variables(
            string(elements[name]["expression"]))) ||
           haskey(elements[name], "active_during") ||
           haskey(elements[name], "inactive_during")]
    stateful_applied_torque_names = [name for name in names_of(("applied_torque",))
        if (haskey(elements[name], "expression") && expression_uses_model_variables(
            string(elements[name]["expression"]))) ||
           haskey(elements[name], "active_during") ||
           haskey(elements[name], "inactive_during")]
    force_names = [names_of(("gravity", "applied_force",
                             "applied_torque",
                             "bushing", "plane_contact",
                             "spanning_force",
                             "torsional_spring_damper",
                             "revolute_friction", "translational_friction",
                             "inplane_friction"));
                   surface_friction_names]
    all(kind -> kind in PLANAR_ELEMENT_TYPES, values(kinds)) ||
        throw(ArgumentError("unsupported element type"))
    isempty(body_names) && throw(ArgumentError("model requires at least one rigid_body"))
    rotation_coordinate_names = Symbol[]
    for name in connection_names
        kinds[name] == "revolute" || continue
        setting = get(elements[name], "rotation_coordinates", false)
        setting isa Bool || throw(ArgumentError(
            "revolute joint '$name' rotation_coordinates must be Boolean"))
        !setting && !isempty(initial_table(elements[name],
            "revolute joint '$name'")) && throw(ArgumentError(
            "revolute joint '$name'.initial requires " *
            "rotation_coordinates = true"))
        setting && push!(rotation_coordinate_names, name)
    end
    valid_velocity_names = Symbol[
        Symbol(name, :., component) for name in body_names
        for component in (:V_x, :V_y, :omega)]
    append!(valid_velocity_names,
        [Symbol(name, :., :omega) for name in rotation_coordinate_names])
    append!(valid_velocity_names,
        [Symbol(name, :., :velocity) for name in distance_coordinate_names])
    state_table = get(document, "state_selection", Dict{String,Any}())
    state_method = get(state_table, "method", "automatic")
    selection_complete = Bool(get(state_table, "_selection_complete", false))
    declared_preferred_names =
        Symbol.(get(state_table, "preferred_velocities", String[]))
    preferred_names = selection_complete ? declared_preferred_names : Symbol[]
    selected_by_body = Dict(name => Symbol[] for name in body_names)
    coordinate_velocity_names = Dict(
        Symbol(name, :., :omega) => name for name in rotation_coordinate_names)
    merge!(coordinate_velocity_names, Dict(
        Symbol(name, :., :velocity) => name
        for name in distance_coordinate_names))
    selected_by_connection = Dict(name => false
        for name in [rotation_coordinate_names; distance_coordinate_names])
    for preferred in preferred_names
        if haskey(coordinate_velocity_names, preferred)
            coordinate_name = coordinate_velocity_names[preferred]
            selected_by_connection[coordinate_name] && throw(ArgumentError(
                "duplicate preferred velocity '$preferred'"))
            selected_by_connection[coordinate_name] = true
            continue
        end
        parts = split(String(preferred), ".")
        component = Symbol(last(parts))
        component in (:V_x, :V_y, :omega) || throw(ArgumentError(
            "preferred velocity must end in V_x, V_y, or omega"))
        body_name = Symbol(join(parts[1:end-1], "."))
        body_name in body_names ||
            throw(ArgumentError("unknown preferred-state body '$body_name'"))
        component in selected_by_body[body_name] &&
            throw(ArgumentError("duplicate preferred velocity '$preferred'"))
        push!(selected_by_body[body_name], component)
    end

    registrations = Dict{Symbol,ComponentRegistration}()
    for name in body_names
        registrations[name] = planar_body_registration(name)
    end
    connection_component_names = Symbol[]
    for name in connection_names
        if kinds[name] == "translational"
            inplane_name = Symbol(name, ".inplane")
            perp_name = Symbol(name, ".perp")
            registrations[inplane_name] =
                inplane_constraint_registration(inplane_name)
            registrations[perp_name] = perp_constraint_registration(perp_name)
            append!(connection_component_names, (inplane_name, perp_name))
            continue
        elseif kinds[name] == "fixed"
            revolute_name = Symbol(name, ".revolute")
            perp_name = Symbol(name, ".perp")
            registrations[revolute_name] =
                revolute_joint_registration(revolute_name)
            registrations[perp_name] = perp_constraint_registration(perp_name)
            append!(connection_component_names, (revolute_name, perp_name))
            continue
        end
        registrations[name] = if kinds[name] == "revolute"
            revolute_joint_registration(name;
                rotation_coordinates = name in rotation_coordinate_names)
        elseif kinds[name] == "distance_coordinate"
            distance_coordinate_registration(name)
        elseif kinds[name] == "coupler"
            coordinate_coupler_registration(name)
        elseif kinds[name] == "inplane"
            inplane_constraint_registration(name)
        elseif kinds[name] == "perp"
            perp_constraint_registration(name)
        elseif kinds[name] == "rack_and_pinion"
            rack_and_pinion_registration(name)
        else
            gear_pair_registration(name)
        end
        push!(connection_component_names, name)
    end
    for name in driver_names
        registrations[name] = kinds[name] == "rotational_motion" ?
            rotational_motion_registration(name) :
            translational_motion_registration(name)
    end
    for name in bushing_names
        registrations[name] = bushing_registration(name)
    end
    for name in contact_names
        registrations[name] = plane_contact_registration(name)
    end
    for name in surface_friction_names
        registrations[name] = planar_surface_friction_registration(name)
    end
    for name in spanning_names
        registrations[name] = spanning_force_registration(name)
    end
    for name in span_measure_names
        registrations[name] = span_measure_registration(name)
    end
    for name in torsional_names
        registrations[name] = torsional_spring_registration(name)
    end
    for name in revolute_friction_names
        registrations[name] = planar_revolute_friction_registration(name)
    end
    for name in translational_friction_names
        registrations[name] = planar_translational_friction_registration(name)
    end
    for name in inplane_friction_names
        registrations[name] = planar_inplane_friction_registration(name)
    end
    for name in stateful_applied_force_names
        registrations[name] = applied_force_registration(name)
    end
    for name in stateful_applied_torque_names
        registrations[name] = applied_torque_registration(name)
    end
    for name in belt_span_names
        registrations[name] = belt_span_registration(name)
    end
    for name in equation_component_names
        registrations[name] = planar_equation_registration(name,
            elements[name])
    end
    builder = ModelLayoutBuilder()
    for name in (body_names..., connection_component_names..., driver_names...,
                 bushing_names..., contact_names..., spanning_names...,
                 span_measure_names...,
                 torsional_names..., stateful_applied_force_names...,
                 stateful_applied_torque_names..., belt_span_names...,
                 friction_names...,
                 equation_component_names...)
        allocate_component_variables!(builder, registrations[name])
    end
    for name in body_names
        allocate_component_equation_block!(builder, registrations[name], :balance)
        for component in (:V_x, :V_y, :omega)
            block = component == :omega ? :selected_state :
                Symbol(:selected_, component)
            allocate_component_equation_block!(builder, registrations[name], block)
        end
    end
    for name in (connection_component_names..., driver_names...)
        for block in (:acceleration, :velocity, :position)
            allocate_component_equation_block!(builder, registrations[name], block)
        end
    end
    for name in rotation_coordinate_names
        for block in (:rotation_acceleration, :rotation_velocity,
                      :rotation_position)
            allocate_component_equation_block!(builder, registrations[name], block)
        end
        allocate_component_equation_block!(builder, registrations[name],
            :selected_rotation_state)
    end
    for name in distance_coordinate_names
        allocate_component_equation_block!(builder, registrations[name],
            :selected_distance_state)
    end
    for name in bushing_names
        allocate_component_equation_block!(builder, registrations[name], :load)
    end
    for name in contact_names
        allocate_component_equation_block!(builder, registrations[name], :contact)
    end
    for name in spanning_names
        for block in (:geometry, :rate, :load)
            allocate_component_equation_block!(builder, registrations[name], block)
        end
    end
    for name in span_measure_names
        for block in (:geometry, :velocity, :acceleration)
            allocate_component_equation_block!(builder, registrations[name],
                block)
        end
    end
    for name in torsional_names
        allocate_component_equation_block!(builder, registrations[name], :load)
    end
    for name in friction_names
        allocate_component_equation_block!(builder, registrations[name],
            :friction)
    end
    for name in (stateful_applied_force_names...,
                 stateful_applied_torque_names...)
        allocate_component_equation_block!(builder, registrations[name], :load)
    end
    for name in belt_span_names
        for block in (:geometry, :rate, :load)
            allocate_component_equation_block!(builder, registrations[name], block)
        end
    end
    for name in equation_component_names
        for block in (:differential, :algebraic)
            haskey(registrations[name].equation_blocks, block) || continue
            allocate_component_equation_block!(builder, registrations[name],
                block)
        end
    end
    layout = finish_layout(builder, collect(values(registrations)))

    bodies = Dict{Symbol,Any}()
    weights = Dict{Symbol,NamedTuple}()
    variable_weights = ones(Float64, length(layout.catalog.variables))
    imposed_variables = Set{Int}()
    characteristic_length_overrides = Dict{Symbol,Float64}()
    initial = zeros(Float64, length(layout.catalog.variables))
    equation_components = Dict{Symbol,PlanarEquationComponent}()
    for name in equation_component_names
        component = allocated_planar_equation_component(layout, name,
            elements[name], parameters)
        initialize_planar_equation_component!(initial, component,
            elements[name])
        equation_components[name] = component
    end
    for name in body_names
        table = elements[name]
        haskey(table, "mass") ||
            throw(ArgumentError("body '$name' requires mass"))
        haskey(table, "inertia") ||
            throw(ArgumentError("body '$name' requires inertia"))
        mass, inertia = Float64(table["mass"]), Float64(table["inertia"])
        mass >= 0 ||
            throw(ArgumentError("body '$name' mass must be nonnegative"))
        inertia >= 0 ||
            throw(ArgumentError("body '$name' inertia must be nonnegative"))
        if haskey(table, "characteristic_length")
            characteristic_length = Float64(table["characteristic_length"])
            characteristic_length > 0 || throw(ArgumentError(
                "body '$name' characteristic_length must be positive"))
            characteristic_length_overrides[name] = characteristic_length
        end
        body = allocated_planar_body(layout, name, mass, inertia;
            selected_states = selected_by_body[name],
            retain_state_candidates = true)
        bodies[name] = body
        initial[body.position_variables] .= require_vector(table, "position", 2;
            default = [0.0, 0.0], label = String(name))
        initial[body.orientation_variable] =
            planar_orientation(table, "body '$name'")
        initial[body.velocity_variables] .= require_vector(table, "velocity", 2;
            default = [0.0, 0.0], label = String(name))
        initial[body.angular_velocity_variable] =
            Float64(get(table, "angular_velocity", 0.0))
        body_weights = body_ic_weights(name, table)
        weights[name] = body_weights
        apply_body_initial_impose!(initial, imposed_variables, body, table,
            name)
    end
    initial_conditions = if isnothing(initial_override)
        apply_saved_initial_conditions!(initial, layout, kinds, document,
            source_directory; dimension = "planar",
            configuration_kinds = CONFIGURATION_KINDS,
            velocity_kinds = VELOCITY_KINDS,
            element_types = SAVED_STATE_ELEMENT_TYPES)
    else
        length(initial_override) == length(initial) ||
            throw(DimensionMismatch("initial override does not match model variables"))
        initial .= initial_override
        isnothing(initial_conditions_override) ?
            disabled_saved_initial_conditions() : initial_conditions_override
    end
    saved_initial_values = initial_conditions.enabled ? copy(initial) : nothing

    markers = Dict{Symbol,Any}()
    for name in fixed_marker_names
        table = elements[name]
        owner_kind, owner_name = marker_owner(name, elements, bodies)
        position = require_vector(table, "position", 2;
            default = [0.0, 0.0], label = String(name))
        angle = planar_orientation(table, "marker '$name'")
        if owner_kind == :body
            body = bodies[owner_name]
            markers[name] = PlanarModelMarker(name, body, planar_point_marker(body, position),
                PlanarBodyOrientationMarker(body.orientation_variable,
                    body.angular_velocity_variable, body.balance_equations[3], angle))
        else
            markers[name] = PlanarModelMarker(name, nothing,
                PlanarGroundPointMarker(position), PlanarGroundOrientationMarker(angle))
        end
    end

    floating_followers = Dict{Symbol,PlanarModelMarker}()
    for name in floating_marker_names
        table = elements[name]
        owner_kind, owner_name = marker_owner(name, elements, bodies)
        owner_kind == :body || throw(ArgumentError(
            "floating marker '$name' must be nested under a rigid_body"))
        haskey(table, "follows") || throw(ArgumentError(
            "floating marker '$name' requires follows"))
        follower_name = Symbol(table["follows"])
        follower = required(markers, follower_name, "marker")
        follower.point isa PlanarFloatingPointMarker && throw(ArgumentError(
            "floating marker '$name' must follow an ordinary marker"))
        body = bodies[owner_name]
        angle = planar_orientation(table, "floating marker '$name'")
        markers[name] = PlanarModelMarker(name, body,
            planar_floating_point_marker(body, follower.point),
            PlanarBodyOrientationMarker(body.orientation_variable,
                body.angular_velocity_variable, body.balance_equations[3], angle))
        floating_followers[name] = follower
    end

    body_lengths = body_characteristic_lengths(bodies, body_names, markers,
        characteristic_length_overrides)
    positive_masses = [Float64(body.mass) for body in values(bodies)
        if body.mass > 0]
    reference_mass = isempty(positive_masses) ? 1.0 : minimum(positive_masses)
    mass_floor = 1.0e-6 * reference_mass
    for name in body_names
        body = bodies[name]
        scale = weights[name].scale
        effective_mass = max(Float64(body.mass), mass_floor)
        translational_weight = scale * effective_mass
        rotational_weight = translational_weight * body_lengths[name]^2
        weights[name] = (; scale, effective_mass,
            characteristic_length = body_lengths[name],
            translational_weight, rotational_weight)
        variable_weights[body.position_variables] .= translational_weight
        variable_weights[body.orientation_variable] = rotational_weight
        variable_weights[body.velocity_variables] .= translational_weight
        variable_weights[body.angular_velocity_variable] = rotational_weight
    end

    measures = Dict{Symbol,Any}()
    for name in span_measure_names
        table = elements[name]
        endpoints = get(table, "markers", Any[])
        length(endpoints) == 2 || throw(ArgumentError(
            "span '$name' requires two markers"))
        marker_1 = required(markers, endpoints[1], "marker")
        marker_2 = required(markers, endpoints[2], "marker")
        (marker_1.point isa PlanarFloatingPointMarker ||
         marker_2.point isa PlanarFloatingPointMarker) &&
            throw(ArgumentError(
                "span '$name' cannot use a floating marker"))
        measure = allocated_span_measure(layout, name, marker_1, marker_2)
        initialize_span_measure!(initial, measure)
        measures[name] = measure
    end

    connections = Dict{Symbol,Any}()
    relative_initial_specifications = Dict{Symbol,NamedTuple}()
    for name in connection_names
        kinds[name] in ("gear_pair", "rack_and_pinion", "coupler") && continue
        table = elements[name]
        endpoints = get(table, "markers", Any[])
        length(endpoints) == 2 || throw(ArgumentError("connection '$name' requires two markers"))
        a, b = required(markers, endpoints[1], "marker"), required(markers, endpoints[2], "marker")
        if kinds[name] == "revolute"
            (a.point isa PlanarFloatingPointMarker ||
             b.point isa PlanarFloatingPointMarker) && throw(ArgumentError(
                "revolute joint '$name' cannot use a floating marker"))
            rotation_coordinates = name in rotation_coordinate_names
            connection = allocated_revolute_joint(layout, name, a.owner,
                a.point, b.owner, b.point;
                rotation_coordinates,
                selected_state = get(selected_by_connection, name, false),
                rotation_marker_a = a.orientation,
                rotation_marker_b = b.orientation,
                retain_state_candidate = true)
            connections[name] = connection
            if rotation_coordinates
                alpha, omega, theta = connection.rotation_variables
                initial[alpha] = (isnothing(a.owner) ? 0.0 :
                    initial[a.owner.angular_acceleration_variable]) -
                    (isnothing(b.owner) ? 0.0 :
                    initial[b.owner.angular_acceleration_variable])
                initial[omega] =
                    PlanarAppliedForces.marker_angular_velocity(a.orientation,
                        initial) -
                    PlanarAppliedForces.marker_angular_velocity(b.orientation,
                        initial)
                initial[theta] = PlanarAppliedForces.marker_angle(a.orientation,
                    initial) - PlanarAppliedForces.marker_angle(b.orientation,
                    initial)
                specification = relative_initial_specification(
                    table, "revolute joint '$name'", :rotation)
                relative_initial_specifications[name] = specification
                apply_relative_initial_specification!(initial,
                    variable_weights, imposed_variables, theta, omega,
                    specification)
                weights[name] = (; angle = specification.coordinate_weight,
                    omega = specification.velocity_weight)
            end
        elseif kinds[name] == "inplane"
            (a.point isa PlanarFloatingPointMarker ||
             b.point isa PlanarFloatingPointMarker) && throw(ArgumentError(
                "inplane constraint '$name' cannot use a floating marker"))
            connections[name] = allocated_inplane_constraint(layout, name,
                a.owner, a.point, b.owner, b.point, b.orientation)
        elseif kinds[name] == "perp"
            connections[name] = allocated_perp_constraint(layout, name,
                a.owner, a.orientation, b.owner, b.orientation)
        elseif kinds[name] == "translational"
            (a.point isa PlanarFloatingPointMarker ||
             b.point isa PlanarFloatingPointMarker) && throw(ArgumentError(
                "translational joint '$name' cannot use a floating marker"))
            connections[name] = allocated_translational_joint(layout, name,
                a.owner, a.point, a.orientation,
                b.owner, b.point, b.orientation)
        elseif kinds[name] == "fixed"
            (a.point isa PlanarFloatingPointMarker ||
             b.point isa PlanarFloatingPointMarker) && throw(ArgumentError(
                "fixed joint '$name' cannot use a floating marker"))
            connections[name] = allocated_fixed_joint(layout, name,
                a.owner, a.point, a.orientation,
                b.owner, b.point, b.orientation)
        elseif kinds[name] == "distance_coordinate"
            (a.point isa PlanarFloatingPointMarker ||
             b.point isa PlanarFloatingPointMarker) && throw(ArgumentError(
                "distance coordinate '$name' cannot use a floating marker"))
            coordinate = allocated_distance_coordinate(layout, name,
                a.owner, a.point, b.owner, b.point, b.orientation;
                selected_state = get(selected_by_connection, name, false),
                retain_state_candidate = true)
            connections[name] = coordinate
            initial[coordinate.distance_variable] =
                PlanarComponentAssembly.inplane_position(coordinate, initial)
            initial[coordinate.velocity_variable] =
                PlanarComponentAssembly.inplane_velocity(coordinate, initial)
            initial[coordinate.acceleration_variable] =
                PlanarComponentAssembly.inplane_acceleration(coordinate, initial)
            specification = relative_initial_specification(
                table, "distance coordinate '$name'", :translation)
            relative_initial_specifications[name] = specification
            apply_relative_initial_specification!(initial, variable_weights,
                imposed_variables, coordinate.distance_variable,
                coordinate.velocity_variable, specification)
            weights[name] = (; distance = specification.coordinate_weight,
                velocity = specification.velocity_weight)
        end
    end

    for name in coupler_names
        table = elements[name]
        port_specs = get(table, "coordinates", Any[])
        length(port_specs) >= 2 || throw(ArgumentError(
            "coupler '$name' requires at least two coordinates"))
        all(spec -> spec isa AbstractString, port_specs) ||
            throw(ArgumentError(
                "coupler '$name' coordinates must be strings"))
        length(unique(port_specs)) == length(port_specs) ||
            throw(ArgumentError(
                "coupler '$name' coordinates must be unique"))
        coordinates = Any[]
        coordinate_kinds = Symbol[]
        for specification in port_specs
            parts = split(String(specification), ".")
            length(parts) >= 2 || throw(ArgumentError(
                "coupler '$name' coordinate '$specification' must name an element port"))
            port = Symbol(last(parts))
            component_name = Symbol(join(parts[1:end-1], "."))
            component = required(connections, component_name,
                "coordinate element")
            if port == :rotation &&
                    component isa PlanarRevoluteJointComponent &&
                    !isempty(component.rotation_variables)
                push!(coordinate_kinds, :rotation)
            elseif port == :distance &&
                    component isa PlanarDistanceCoordinateComponent
                push!(coordinate_kinds, :distance)
            else
                throw(ArgumentError(
                    "coupler '$name' coordinate '$specification' is not an available coordinate port"))
            end
            push!(coordinates, component)
        end
        coordinate_kind = length(unique(coordinate_kinds)) == 1 ?
            first(coordinate_kinds) : :mixed
        raw_coefficients = get(table, "coefficients", Any[])
        length(raw_coefficients) == length(coordinates) ||
            throw(ArgumentError(
                "coupler '$name' requires one coefficient per coordinate"))
        all(value -> value isa Number, raw_coefficients) ||
            throw(ArgumentError(
                "coupler '$name' coefficients must be numeric"))
        coefficients = Float64.(raw_coefficients)
        all(isfinite, coefficients) || throw(ArgumentError(
            "coupler '$name' coefficients must be finite"))
        count(value -> !iszero(value), coefficients) >= 2 || throw(ArgumentError(
            "coupler '$name' requires at least two nonzero coefficients"))
        initial_sum = sum(factor * initial[
            PlanarComponentAssembly.coordinate_variables(coordinate)[1]]
            for (coordinate, factor) in zip(coordinates, coefficients))
        offset_specification = get(table, "offset", "initial")
        offset = if offset_specification isa Number
            Float64(offset_specification)
        elseif offset_specification == "initial"
            initial_sum
        else
            throw(ArgumentError(
                "coupler '$name' offset must be numeric or 'initial'"))
        end
        isfinite(offset) || throw(ArgumentError(
            "coupler '$name' offset must be finite"))
        connections[name] = allocated_coordinate_coupler(layout, name,
            coordinates, coefficients, offset, coordinate_kind)
    end

    for name in connection_names
        kinds[name] == "gear_pair" || continue
        table = elements[name]
        joint_names = get(table, "joints", Any[])
        length(joint_names) == 2 || throw(ArgumentError(
            "gear pair '$name' requires two joints"))
        joint_names[1] != joint_names[2] || throw(ArgumentError(
            "gear pair '$name' requires two different joints"))
        joints = [required(connections, joint_name, "joint")
                  for joint_name in joint_names]
        all(joint -> joint isa PlanarRevoluteJointComponent, joints) ||
            throw(ArgumentError(
                "gear pair '$name' joints must be revolute joints"))

        haskey(table, "contact_marker") || throw(ArgumentError(
            "gear pair '$name' requires contact_marker"))
        contact = required(markers, table["contact_marker"], "contact marker")
        contact.point isa PlanarFloatingPointMarker && throw(ArgumentError(
            "gear pair '$name' contact_marker must be fixed to its carrier"))

        gear_bodies = Any[]
        gear_centers = Any[]
        floating = PlanarModelMarker[]
        for (index, joint) in enumerate(joints)
            gear_body = joint.body_a
            isnothing(gear_body) && throw(ArgumentError(
                "gear pair '$name' joint '$(joint.name)' must list its " *
                "gear marker first"))
            push!(gear_bodies, gear_body)
            push!(gear_centers, joint.marker_b)
            floating_name = Symbol(name, ".contact_", index)
            haskey(markers, floating_name) && throw(ArgumentError(
                "gear pair '$name' generated marker '$floating_name' " *
                "conflicts with an existing marker"))
            marker = PlanarModelMarker(floating_name, gear_body,
                planar_floating_point_marker(gear_body, contact.point),
                joint.rotation_marker_a)
            markers[floating_name] = marker
            push!(floating, marker)
        end
        gear_bodies[1] !== gear_bodies[2] || throw(ArgumentError(
            "gear pair '$name' joints must identify two different gear bodies"))

        contact_position = PlanarAppliedForces.point_marker_kinematics(
            contact.point, initial).position
        center_positions = [PlanarAppliedForces.point_marker_kinematics(
            center, initial).position for center in gear_centers]
        radial_vectors = [contact_position - center
                          for center in center_positions]
        radii = norm.(radial_vectors)
        all(>(0), radii) || throw(ArgumentError(
            "gear pair '$name' contact point must differ from both joint centers"))
        normalized_dot = dot(radial_vectors[1], radial_vectors[2]) /
            (radii[1] * radii[2])
        isapprox(abs(normalized_dot), 1.0; atol = 1.0e-8) ||
            throw(ArgumentError(
                "gear pair '$name' joint centers and contact point must be collinear"))
        eligible_references = findall(
            index -> joints[index].body_b === contact.owner,
            eachindex(joints))
        isempty(eligible_references) && throw(ArgumentError(
            "gear pair '$name' contact_marker must be fixed to the base body " *
            "of at least one referenced revolute joint"))
        reference_index = first(eligible_references)
        reference_joint = joints[reference_index]
        reference_radial = radial_vectors[reference_index]
        reference_direction = reference_radial ./ radii[reference_index]
        signed_radii = [
            dot(radial_vectors[1], reference_direction),
            dot(radial_vectors[2], reference_direction),
        ]

        carrier_angle = isnothing(contact.owner) ? 0.0 :
            initial[contact.owner.orientation_variable]
        initial_phase = signed_radii[1] *
            (PlanarAppliedForces.marker_angle(floating[1].orientation,
                initial) - carrier_angle) -
            signed_radii[2] *
            (PlanarAppliedForces.marker_angle(floating[2].orientation,
                initial) - carrier_angle)
        phase_specification = get(table, "phase", "initial")
        phase = if phase_specification isa Number
            Float64(phase_specification)
        elseif phase_specification == "initial"
            initial_phase
        else
            throw(ArgumentError(
                "gear pair '$name' phase must be numeric or 'initial'"))
        end
        isfinite(phase) || throw(ArgumentError(
            "gear pair '$name' phase must be finite"))
        connections[name] = allocated_gear_pair(layout, name,
            gear_bodies[1], floating[1].point, floating[1].orientation,
            gear_bodies[2], floating[2].point, floating[2].orientation,
            joints[1], joints[2], contact.owner, contact.point,
            reference_joint.marker_b, reference_joint.rotation_marker_b,
            signed_radii[1], signed_radii[2], phase)
    end

    for name in connection_names
        kinds[name] == "rack_and_pinion" || continue
        table = elements[name]
        joint_names = get(table, "joints", Any[])
        length(joint_names) == 2 || throw(ArgumentError(
            "rack and pinion '$name' requires two joints"))
        translational = required(connections, joint_names[1], "joint")
        revolute = required(connections, joint_names[2], "joint")
        translational isa PlanarTranslationalJoint || throw(ArgumentError(
            "rack and pinion '$name' first joint must be translational"))
        revolute isa PlanarRevoluteJointComponent || throw(ArgumentError(
            "rack and pinion '$name' second joint must be revolute"))
        haskey(table, "pitch_radius") || throw(ArgumentError(
            "rack and pinion '$name' requires pitch_radius"))
        pitch_radius = Float64(table["pitch_radius"])
        pitch_radius > 0 || throw(ArgumentError(
            "rack and pinion '$name' pitch_radius must be positive"))

        rack_body = translational.inplane.geometry.body_i
        carrier_body = translational.inplane.geometry.body_j
        rack_body isa PlanarRigidBodyComponent || throw(ArgumentError(
            "rack and pinion '$name' translational joint must put the rack first"))
        carrier_center, pinion_body, pinion_orientation,
        pinion_carrier_orientation =
            if revolute.body_a === carrier_body &&
                    revolute.body_b isa PlanarRigidBodyComponent
                (revolute.marker_a, revolute.body_b,
                 revolute.rotation_marker_b, revolute.rotation_marker_a)
            elseif revolute.body_b === carrier_body &&
                    revolute.body_a isa PlanarRigidBodyComponent
                (revolute.marker_b, revolute.body_a,
                 revolute.rotation_marker_a, revolute.rotation_marker_b)
            else
                throw(ArgumentError(
                    "rack and pinion '$name' joints must share a carrier"))
            end
        pinion_body !== rack_body || throw(ArgumentError(
            "rack and pinion '$name' requires different rack and pinion bodies"))
        carrier_orientation = translational.perp.orientation_j
        local_axis_angle = local_marker_angle(carrier_orientation)
        local_normal = [-sin(local_axis_angle), cos(local_axis_angle)]
        contact_point = if isnothing(carrier_body)
            PlanarGroundPointMarker(carrier_center.position +
                pitch_radius .* local_normal)
        else
            planar_point_marker(carrier_body, carrier_center.r_body +
                pitch_radius .* local_normal)
        end

        carrier_contact_name = Symbol(name, ".carrier_contact")
        rack_contact_name = Symbol(name, ".rack_contact")
        pinion_contact_name = Symbol(name, ".pinion_contact")
        any(haskey(markers, generated) for generated in
            (carrier_contact_name, rack_contact_name,
             pinion_contact_name)) && throw(ArgumentError(
            "rack and pinion '$name' generated contact marker name is already used"))
        rack_contact = planar_floating_point_marker(rack_body, contact_point)
        pinion_contact = planar_floating_point_marker(
            pinion_body, contact_point)
        markers[carrier_contact_name] = PlanarModelMarker(
            carrier_contact_name, carrier_body, contact_point,
            carrier_orientation)
        markers[rack_contact_name] = PlanarModelMarker(
            rack_contact_name, rack_body, rack_contact,
            translational.perp.orientation_i)
        markers[pinion_contact_name] = PlanarModelMarker(
            pinion_contact_name, pinion_body, pinion_contact,
            pinion_orientation)

        initial_inplane = PlanarComponentAssembly.inplane_values(
            translational.inplane, initial)
        initial_axis = -initial_inplane.direction.normal
        initial_rack_position = dot(initial_inplane.separation, initial_axis)
        initial_relative_angle = PlanarAppliedForces.marker_angle(
            pinion_orientation, initial) -
            PlanarAppliedForces.marker_angle(
                pinion_carrier_orientation, initial)
        initial_phase = initial_rack_position +
            pitch_radius * initial_relative_angle
        phase_specification = get(table, "phase", "initial")
        phase = if phase_specification isa Number
            Float64(phase_specification)
        elseif phase_specification == "initial"
            initial_phase
        else
            throw(ArgumentError(
                "rack and pinion '$name' phase must be numeric or 'initial'"))
        end
        isfinite(phase) || throw(ArgumentError(
            "rack and pinion '$name' phase must be finite"))
        connections[name] = allocated_rack_and_pinion(layout, name,
            rack_body, rack_contact,
            pinion_body, pinion_contact, pinion_orientation,
            translational, revolute, carrier_body, contact_point,
            pinion_carrier_orientation, pitch_radius, phase)
    end

    for name in pulley_names
        table = elements[name]
        haskey(table, "joint") || throw(ArgumentError(
            "pulley '$name' requires joint"))
        joint = required(connections, table["joint"], "revolute joint")
        joint isa PlanarRevoluteJointComponent || throw(ArgumentError(
            "pulley '$name' joint must be revolute"))
        body_name = if haskey(table, "body")
            Symbol(table["body"])
        else
            parents = [parent for parent in parent_names(name)
                if haskey(bodies, parent)]
            isempty(parents) && throw(ArgumentError(
                "pulley '$name' requires body unless nested beneath its rigid body"))
            first(parents)
        end
        body = required(bodies, body_name, "pulley body")
        center_marker = if joint.body_a === body
            joint.marker_a
        elseif joint.body_b === body
            joint.marker_b
        else
            throw(ArgumentError(
                "pulley '$name' body '$body_name' is not attached to joint '$(table["joint"])'"))
        end
        haskey(table, "pitch_radius") || throw(ArgumentError(
            "pulley '$name' requires pitch_radius"))
        pitch_radius = Float64(table["pitch_radius"])
        pitch_radius > 0 || throw(ArgumentError(
            "pulley '$name' pitch_radius must be positive"))
        connections[name] = PlanarPulleyComponent(
            name, body, joint, center_marker, pitch_radius)
    end

    belt_forces = Dict{Symbol,Vector{Any}}()
    claimed_spans = Symbol[]
    for belt_name in belt_names
        table = elements[belt_name]
        ordered_names = Symbol.(get(table, "spans", String[]))
        length(ordered_names) >= 2 || throw(ArgumentError(
            "belt '$belt_name' requires at least two ordered spans"))
        length(unique(ordered_names)) == length(ordered_names) ||
            throw(ArgumentError("belt '$belt_name' spans must be unique"))
        specifications = NamedTuple[]
        for span_name in ordered_names
            span_name in belt_span_names || throw(ArgumentError(
                "belt '$belt_name' references unknown belt span '$span_name'"))
            push!(claimed_spans, span_name)
            span_table = elements[span_name]
            endpoints = get(span_table, "pulleys", Any[])
            length(endpoints) == 2 || throw(ArgumentError(
                "belt span '$span_name' requires two pulleys"))
            pulley_1 = required(connections, endpoints[1], "pulley")
            pulley_2 = required(connections, endpoints[2], "pulley")
            (pulley_1 isa PlanarPulleyComponent &&
             pulley_2 isa PlanarPulleyComponent) || throw(ArgumentError(
                "belt span '$span_name' endpoints must be pulleys"))
            pulley_1 !== pulley_2 || throw(ArgumentError(
                "belt span '$span_name' requires two different pulleys"))
            hints = get(span_table, "near_points", Any[])
            length(hints) == 2 || throw(ArgumentError(
                "belt span '$span_name' requires two near_points"))
            hint_1 = required(markers, hints[1], "near-point marker")
            hint_2 = required(markers, hints[2], "near-point marker")
            selected = select_belt_tangent(
                pulley_1, pulley_2, hint_1, hint_2, initial)
            stiffness = Float64(get(span_table, "stiffness",
                get(table, "stiffness", 0.0)))
            stiffness > 0 || throw(ArgumentError(
                "belt span '$span_name' stiffness must be positive"))
            damping_time_scale = Float64(get(span_table,
                "damping_time_scale", get(table, "damping_time_scale", 0.0)))
            damping_time_scale >= 0 || throw(ArgumentError(
                "belt span '$span_name' damping_time_scale must be nonnegative"))
            damping = Float64(get(span_table, "damping",
                get(table, "damping", damping_time_scale * stiffness)))
            damping >= 0 || throw(ArgumentError(
                "belt span '$span_name' damping must be nonnegative"))
            push!(specifications, (; name = span_name, pulley_1, pulley_2,
                selected.tangent_kind, selected.tangent_side,
                geometry = selected.geometry, stiffness, damping,
                damping_time_scale))
        end
        path_length = belt_reference_path_length(specifications)
        has_initial_tension = haskey(table, "initial_tension")
        has_free_length = haskey(table, "free_length")
        has_initial_tension && has_free_length && throw(ArgumentError(
            "belt '$belt_name' cannot specify both initial_tension and free_length"))
        initial_tension = if has_initial_tension
            specified = Float64(table["initial_tension"])
            specified >= 0 || throw(ArgumentError(
                "belt '$belt_name' initial_tension must be nonnegative"))
            specified
        elseif has_free_length
            free = Float64(table["free_length"])
            free > 0 || throw(ArgumentError(
                "belt '$belt_name' free_length must be positive"))
            (path_length - free) /
                sum(inv(spec.stiffness) for spec in specifications)
        else
            0.0
        end
        if initial_tension < 0
            @warn "belt '$belt_name' begins with negative tension" initial_tension
        end
        free_length = has_free_length ? Float64(table["free_length"]) :
            path_length - sum(initial_tension / spec.stiffness
                              for spec in specifications)
        spans = Any[]
        for spec in specifications
            geometry = spec.geometry
            span = allocated_belt_span(layout, spec.name, belt_name,
                spec.pulley_1, spec.pulley_2, spec.tangent_kind,
                spec.tangent_side, geometry.beta_1, geometry.beta_2,
                geometry.length, geometry.tangent,
                initial[spec.pulley_1.body.orientation_variable],
                initial[spec.pulley_2.body.orientation_variable],
                initial_tension / spec.stiffness, spec.stiffness,
                spec.damping, spec.damping_time_scale)
            initialize_belt_span!(initial, span)
            push!(spans, span)
        end
        connections[belt_name] = PlanarBeltComponent(
            belt_name, spans, initial_tension, free_length)
        belt_forces[belt_name] = spans
    end
    sort!(claimed_spans) == sort!(belt_span_names) || throw(ArgumentError(
        "every belt_span must belong to exactly one belt"))

    forces = Dict{Symbol,Any}()
    merge!(forces, belt_forces)
    for name in force_names
        table, kind = elements[name], kinds[name]
        staged = haskey(table, "active_during") ||
            haskey(table, "inactive_during")
        stage_supported = kind in ("applied_force", "applied_torque",
            "spanning_force", "bushing", "plane_contact")
        staged && !stage_supported && throw(ArgumentError(
            "force '$name' of type '$kind' does not support " *
            "active_during or inactive_during"))
        if kind == "gravity"
            acceleration = require_vector(table, "acceleration", 2; label = String(name))
            targets = get(table, "bodies", string.(body_names))
            forces[name] = [PlanarGravityComponent(Symbol(name, :_, Symbol(target)),
                required(bodies, target, "body"), acceleration) for target in targets]
        elseif kind == "applied_force"
            active_during = active_during_stages(table,
                "applied force '$name'")
            endpoints = get(table, "markers", Any[])
            length(endpoints) == 2 || throw(ArgumentError(
                "applied force '$name' requires an application marker and a direction marker"))
            application = required(markers, endpoints[1], "marker")
            application.owner === nothing && throw(ArgumentError(
                "applied force '$name' application marker must belong to a body"))
            direction = required(markers, endpoints[2], "marker")
            has_force = haskey(table, "force")
            has_expression = haskey(table, "expression")
            xor(has_force, has_expression) || throw(ArgumentError(
                "applied force '$name' requires exactly one of force or expression"))
            law = has_force ? constant_scalar_law(table["force"]) :
                compile_model_expression(string(table["expression"]),
                    parameters, layout)
            reaction_point = nothing
            if haskey(table, "reaction_body")
                reaction_body_name = Symbol(table["reaction_body"])
                reaction_body = required(bodies, reaction_body_name,
                    "reaction body")
                reaction_body === application.owner && throw(ArgumentError(
                    "applied force '$name' reaction_body must differ from the application body"))
                reaction_name = Symbol(name, ".reaction")
                haskey(markers, reaction_name) && throw(ArgumentError(
                    "applied force '$name' generated marker '$reaction_name' already exists"))
                reaction_point = planar_floating_point_marker(
                    reaction_body, application.point)
                markers[reaction_name] = PlanarModelMarker(reaction_name,
                    reaction_body, reaction_point,
                    PlanarBodyOrientationMarker(
                        reaction_body.orientation_variable,
                        reaction_body.angular_velocity_variable,
                        reaction_body.balance_equations[3], 0.0))
            end
            component = if name in stateful_applied_force_names
                allocated_applied_force(layout, name, application.point,
                    direction.owner, direction.orientation, reaction_point,
                    law; active_during)
            else
                axis = PlanarDirectedAxis(
                    direction.owner, direction.orientation)
                PlanarAppliedForceComponent(name, application.point,
                    axis, reaction_point, law)
            end
            if component.magnitude_variable != 0
                initial[component.magnitude_variable] = component.active[] ?
                    law(simulation.start_time, initial) : 0.0
                isfinite(initial[component.magnitude_variable]) ||
                    throw(ArgumentError(
                        "applied force '$name' expression is not finite initially"))
            end
            forces[name] = [component]
        elseif kind == "spanning_force"
            active_during = active_during_stages(table,
                "spanning force '$name'")
            endpoints = get(table, "markers", Any[])
            length(endpoints) == 2 || throw(ArgumentError(
                "spanning force '$name' requires two markers"))
            marker_1 = required(markers, endpoints[1], "marker")
            marker_2 = required(markers, endpoints[2], "marker")
            (marker_1.point isa PlanarFloatingPointMarker ||
             marker_2.point isa PlanarFloatingPointMarker) &&
                throw(ArgumentError(
                    "spanning force '$name' cannot use a floating marker"))
            has_constant = haskey(table, "force")
            has_expression = haskey(table, "expression")
            has_linear = haskey(table, "stiffness")
            count(identity, (has_constant, has_expression, has_linear)) == 1 ||
                throw(ArgumentError("spanning force '$name' requires exactly " *
                    "one of force, expression, or spring-damper coefficients"))
            linear_extras = any(haskey(table, field)
                for field in ("damping", "free_length"))
            !has_linear && linear_extras && throw(ArgumentError(
                "spanning force '$name' damping and free_length require stiffness"))
            stiffness = damping = free_length = 0.0
            law = if has_constant
                constant_scalar_law(table["force"])
            elseif has_expression
                compile_model_expression(string(table["expression"]),
                    parameters, layout)
            else
                haskey(table, "damping") || throw(ArgumentError(
                    "spanning force '$name' spring-damper law requires damping"))
                haskey(table, "free_length") || throw(ArgumentError(
                    "spanning force '$name' spring-damper law requires free_length"))
                stiffness = Float64(table["stiffness"])
                damping = Float64(table["damping"])
                free_length = Float64(table["free_length"])
                min(stiffness, damping, free_length) >= 0 ||
                    throw(ArgumentError("spanning force '$name' spring-damper " *
                        "coefficients must be nonnegative"))
                linear_spanning_law(layout, name, stiffness, damping, free_length)
            end
            force = allocated_spanning_force(layout, name,
                marker_1.point, marker_2.point, law;
                stiffness, damping, free_length, active_during)
            initialize_spanning_force!(initial, force, simulation.start_time)
            forces[name] = [force]
        elseif kind == "torsional_spring_damper"
            endpoints = get(table, "markers", Any[])
            length(endpoints) == 2 || throw(ArgumentError(
                "torsional spring-damper '$name' requires two markers"))
            marker_1 = required(markers, endpoints[1], "marker")
            marker_2 = required(markers, endpoints[2], "marker")
            haskey(table, "stiffness") || throw(ArgumentError(
                "torsional spring-damper '$name' requires stiffness"))
            stiffness = Float64(table["stiffness"])
            damping_time_scale =
                Float64(get(table, "damping_time_scale", 0.0))
            damping = haskey(table, "damping") ?
                Float64(table["damping"]) :
                damping_time_scale * stiffness
            min(stiffness, damping, damping_time_scale) >= 0 ||
                throw(ArgumentError(
                    "torsional spring-damper '$name' coefficients must be nonnegative"))
            free_angle = angle_value(get(table, "free_angle", 0.0),
                "torsional spring-damper '$name'.free_angle")
            spring = allocated_torsional_spring_damper(layout, name,
                marker_1, marker_2, stiffness, damping,
                damping_time_scale, free_angle)
            initialize_torsional_force!(initial, spring)
            forces[name] = [spring]
        elseif kind == "bushing"
            active_during = active_during_stages(table, "bushing '$name'")
            endpoints = get(table, "markers", Any[])
            length(endpoints) == 2 ||
                throw(ArgumentError("bushing '$name' requires two markers"))
            marker_1 = required(markers, endpoints[1], "marker")
            marker_2 = required(markers, endpoints[2], "marker")
            translational_stiffness = require_vector(table,
                "translational_stiffness", 2; label = String(name))
            all(>=(0), translational_stiffness) || throw(ArgumentError(
                "bushing '$name' translational stiffness must be nonnegative"))
            damping_time_scale =
                Float64(get(table, "damping_time_scale", 0.0))
            damping_time_scale >= 0 || throw(ArgumentError(
                "bushing '$name' damping_time_scale must be nonnegative"))
            translational_damping = haskey(table, "translational_damping") ?
                require_vector(table, "translational_damping", 2;
                    label = String(name)) :
                damping_time_scale .* translational_stiffness
            all(>=(0), translational_damping) || throw(ArgumentError(
                "bushing '$name' translational damping must be nonnegative"))
            rotational_stiffness =
                Float64(get(table, "rotational_stiffness", 0.0))
            rotational_damping = haskey(table, "rotational_damping") ?
                Float64(table["rotational_damping"]) :
                damping_time_scale * rotational_stiffness
            min(rotational_stiffness, rotational_damping) >= 0 ||
                throw(ArgumentError(
                    "bushing '$name' rotational coefficients must be nonnegative"))
            free_position = require_vector(table, "free_position", 2;
                default = [0.0, 0.0], label = String(name))
            free_angle = angle_value(get(table, "free_angle", 0.0),
                "bushing '$name'.free_angle")
            forces[name] = [allocated_planar_bushing(layout, name,
                marker_1, marker_2, translational_stiffness,
                translational_damping, rotational_stiffness,
                rotational_damping, damping_time_scale,
                free_position, free_angle; active_during)]
        elseif kind == "plane_contact"
            active_during = active_during_stages(table,
                "plane contact '$name'")
            endpoints = get(table, "markers", Any[])
            length(endpoints) == 2 || throw(ArgumentError(
                "plane contact '$name' requires two markers"))
            marker_1 = required(markers, endpoints[1], "marker")
            marker_2 = required(markers, endpoints[2], "marker")
            haskey(table, "radius") || throw(ArgumentError(
                "plane contact '$name' requires radius"))
            radius = Float64(table["radius"])
            radius > 0 || throw(ArgumentError(
                "plane contact '$name' radius must be positive"))
            stiffness = Float64(get(table, "stiffness", 0.0))
            stiffness > 0 || throw(ArgumentError(
                "plane contact '$name' stiffness must be positive"))
            damping_factor = Float64(get(table, "damping_factor", 0.0))
            damping_factor >= 0 || throw(ArgumentError(
                "plane contact '$name' damping_factor must be nonnegative"))
            forces[name] = [allocated_plane_contact(layout, name,
                marker_1, marker_2, radius, stiffness, damping_factor;
                active_during)]
        elseif kind == "surface_friction"
            contact_name = get(table, "contact", nothing)
            contact_name isa AbstractString || throw(ArgumentError(
                "surface friction '$name'.contact must name a plane_contact"))
            contact_collection = get(forces, Symbol(contact_name), nothing)
            contact = contact_collection isa AbstractVector &&
                length(contact_collection) == 1 ? only(contact_collection) :
                nothing
            contact isa PlanarPlaneContactComponent || throw(ArgumentError(
                "surface friction '$name'.contact must name an existing " *
                "plane_contact"))
            haskey(table, "preload") && throw(ArgumentError(
                "surface friction '$name' cannot specify preload; its " *
                "capacity comes from the contact normal force"))
            settings = friction_parameters(table,
                "surface friction '$name'")
            friction = allocated_planar_surface_friction(layout, name,
                contact, settings.stiffness, settings.damping,
                settings.static_coefficient, settings.dynamic_coefficient,
                settings.transition_speed, settings.release_time)
            initialize_planar_surface_friction!(initial, friction;
                reset_anchor = true)
            forces[name] = [friction]
        elseif kind == "revolute_friction"
            joint_name = get(table, "joint", nothing)
            joint_name isa AbstractString || throw(ArgumentError(
                "revolute friction '$name'.joint must name a revolute joint"))
            joint = get(connections, Symbol(joint_name), nothing)
            joint isa PlanarRevoluteJointComponent || throw(ArgumentError(
                "revolute friction '$name'.joint must name an existing " *
                "revolute joint"))
            settings = friction_parameters(table,
                "revolute friction '$name'"; angular = true,
                effective_radius = true)
            friction = allocated_planar_revolute_friction(layout, name,
                joint, settings.stiffness, settings.damping,
                settings.effective_radius, settings.preload,
                settings.static_coefficient, settings.dynamic_coefficient,
                settings.transition_speed, settings.release_time)
            initialize_planar_revolute_friction!(initial, friction;
                reset_anchor = true)
            forces[name] = [friction]
        elseif kind == "translational_friction"
            joint_name = get(table, "joint", nothing)
            joint_name isa AbstractString || throw(ArgumentError(
                "translational friction '$name'.joint must name a " *
                "translational joint"))
            joint = get(connections, Symbol(joint_name), nothing)
            joint isa PlanarTranslationalJoint || throw(ArgumentError(
                "translational friction '$name'.joint must name an " *
                "existing translational joint"))
            settings = friction_parameters(table,
                "translational friction '$name'")
            friction = allocated_planar_translational_friction(layout,
                name, joint, settings.stiffness, settings.damping,
                settings.preload, settings.static_coefficient,
                settings.dynamic_coefficient, settings.transition_speed,
                settings.release_time)
            initialize_planar_translational_friction!(initial, friction;
                reset_anchor = true)
            forces[name] = [friction]
        elseif kind == "inplane_friction"
            constraint_name = get(table, "constraint", nothing)
            constraint_name isa AbstractString || throw(ArgumentError(
                "inplane friction '$name'.constraint must name an " *
                "inplane constraint"))
            constraint = get(connections, Symbol(constraint_name), nothing)
            constraint isa PlanarInplaneConstraint || throw(ArgumentError(
                "inplane friction '$name'.constraint must name an existing " *
                "standalone inplane constraint"))
            settings = friction_parameters(table,
                "inplane friction '$name'")
            friction = allocated_planar_inplane_friction(layout, name,
                constraint, settings.stiffness, settings.damping,
                settings.preload, settings.static_coefficient,
                settings.dynamic_coefficient, settings.transition_speed,
                settings.release_time)
            initialize_planar_inplane_friction!(initial, friction;
                reset_anchor = true)
            forces[name] = [friction]
        else
            active_during = active_during_stages(table,
                "applied torque '$name'")
            endpoints = get(table, "markers", Any[])
            length(endpoints) == 2 || throw(ArgumentError("force '$name' requires two markers"))
            marker_a = required(markers, endpoints[1], "marker")
            marker_b = required(markers, endpoints[2], "marker")
            a = marker_a.orientation
            b = marker_b.orientation
            has_torque = haskey(table, "torque")
            has_expression = haskey(table, "expression")
            xor(has_torque, has_expression) || throw(ArgumentError(
                "applied torque '$name' requires exactly one of torque or expression"))
            law = has_torque ? constant_scalar_law(table["torque"]) :
                compile_model_expression(string(table["expression"]),
                    parameters, layout)
            component = if name in stateful_applied_torque_names
                allocated_applied_torque(layout, name, a, b, law,
                    marker_a.point, marker_b.point; active_during)
            else
                PlanarAppliedTorqueComponent(name, a, b, law,
                    marker_a.point, marker_b.point)
            end
            if component.torque_variable != 0
                initial[component.torque_variable] = component.active[] ?
                    law(simulation.start_time, initial) : 0.0
                isfinite(initial[component.torque_variable]) ||
                    throw(ArgumentError(
                        "applied torque '$name' expression is not finite initially"))
            end
            forces[name] = [component]
        end
    end

    drivers = Dict{Symbol,Any}()
    for name in driver_names
        table = elements[name]
        endpoints = get(table, "markers", Any[])
        length(endpoints) == 2 || throw(ArgumentError("driver '$name' requires two markers"))
        a, b = required(markers, endpoints[1], "marker"), required(markers, endpoints[2], "marker")
        driver = if kinds[name] == "rotational_motion"
            motion, velocity, acceleration =
                driver_laws(table, parameters, name)
            allocated_rotational_motion_generator(layout, name,
                a.owner, a.orientation, b.owner, b.orientation,
                motion, velocity, acceleration)
        else
            (a.point isa PlanarFloatingPointMarker ||
             b.point isa PlanarFloatingPointMarker) && throw(ArgumentError(
                "translational driver '$name' cannot use a floating marker"))
            motion, velocity, acceleration =
                translational_driver_laws(table, parameters, name)
            allocated_translational_motion_generator(layout, name,
                a.owner, a.point, b.owner, b.point, b.orientation,
                motion, velocity, acceleration)
        end
        drivers[name] = driver
        if driver isa PlanarRotationalMotionGenerator
            initial[driver.angle_variable] = driver.motion(simulation.start_time)
            initial[driver.angular_velocity_variable] =
                driver.motion_derivative(simulation.start_time)
            initial[driver.angular_acceleration_variable] =
                driver.motion_second_derivative(simulation.start_time)
        else
            initial[driver.distance_variable] = driver.motion(simulation.start_time)
            initial[driver.velocity_variable] =
                driver.motion_derivative(simulation.start_time)
            initial[driver.acceleration_variable] =
                driver.motion_second_derivative(simulation.start_time)
        end
    end
    if !isnothing(initial_override)
        length(initial_override) == length(initial) ||
            throw(DimensionMismatch("initial override does not match model variables"))
        initial .= initial_override
        for measure in values(measures)
            initialize_span_measure!(initial, measure)
        end
    elseif initial_conditions.enabled
        transferred = Set(initial_conditions.transferred_variables)
        for variable in layout.catalog.variables
            qualified = Symbol(variable.component, :., variable.name)
            qualified in transferred || continue
            initial[variable.index] = saved_initial_values[variable.index]
        end
    end
    if isnothing(initial_override)
        for name in body_names
            apply_body_initial_impose!(initial, imposed_variables,
                bodies[name], elements[name], name)
        end
        for (name, specification) in relative_initial_specifications
            coordinate = connections[name]
            coordinate_index, velocity_index, _ =
                PlanarComponentAssembly.coordinate_variables(coordinate)
            apply_relative_initial_specification!(initial, variable_weights,
                imposed_variables, coordinate_index, velocity_index,
                specification)
        end
    end
    for name in friction_names
        friction = only(forces[name])
        if friction isa PlanarSurfaceFriction
            initialize_planar_surface_friction!(initial, friction;
                reset_anchor = true)
        elseif friction isa PlanarRevoluteFriction
            initialize_planar_revolute_friction!(initial, friction;
                reset_anchor = true)
        elseif friction isa PlanarTranslationalFriction
            initialize_planar_translational_friction!(initial, friction;
                reset_anchor = true)
        elseif friction isa PlanarInplaneFriction
            initialize_planar_inplane_friction!(initial, friction;
                reset_anchor = true)
        end
    end

    entered_initial_values = isnothing(entered_initial_override) ?
        copy(initial) : copy(entered_initial_override)

    bushing_components = [only(forces[name]) for name in bushing_names]
    contact_components = [only(forces[name]) for name in contact_names]
    spanning_components = [only(forces[name]) for name in spanning_names]
    torsional_components = [only(forces[name]) for name in torsional_names]
    friction_components = [only(forces[name]) for name in friction_names]
    applied_load_components = [only(forces[name]) for name in
        (stateful_applied_force_names..., stateful_applied_torque_names...)]
    belt_span_components = isempty(belt_forces) ? Any[] :
        collect(Iterators.flatten(values(belt_forces)))
    owned_components = [collect(values(bodies)); collect(values(connections));
                        collect(values(drivers)); collect(values(measures));
                        bushing_components;
                        contact_components; spanning_components;
                        torsional_components; applied_load_components;
                        belt_span_components; friction_components;
                        collect(values(equation_components))]
    contribution_components = Any[]
    for collection in values(forces)
        append!(contribution_components, collection)
    end
    append!(contribution_components, values(connections))
    append!(contribution_components, values(drivers))
    model = assemble_planar_model(layout, owned_components,
        contribution_components)
    selected_count = sum(length, values(selected_by_body)) +
        count(values(selected_by_connection))
    if !selection_complete
        position_corrections = correct_initial_configuration!(initial, model,
            variable_weights, imposed_variables, simulation.start_time)
        velocity_corrections = correct_initial_velocities!(initial, model,
            variable_weights, imposed_variables, simulation.start_time)
        automatic_choice, diagnostics = automatic_velocity_selection(
            model, bodies, body_names, markers, connections, connection_names,
            characteristic_length_overrides, initial, simulation.start_time)
        chosen = automatic_choice
        fallback_used = false
        dependent_condition = nothing
        settings = parse_state_selection(document, valid_velocity_names)
        if state_method == "preferred"
            check = check_preferred_velocities(declared_preferred_names,
                diagnostics, settings.maximum_condition_number)
            dependent_condition = check.condition
            if check.acceptable
                chosen = declared_preferred_names
            elseif settings.allow_fallback
                fallback_used = true
            else
                throw(ArgumentError("preferred state selection is unusable: " *
                    check.reason))
            end
        end
        revised = deepcopy(document)
        revised_state = Dict{String,Any}(
            "method" => state_method,
            "preferred_velocities" => String.(chosen),
            "allow_fallback" => settings.allow_fallback,
            "maximum_condition_number" => settings.maximum_condition_number,
            "_qr_rank" => diagnostics.rank,
            "_qr_tolerance" => diagnostics.tolerance,
            "_qr_pivots" => String.(diagnostics.pivots),
            "_body_characteristic_lengths" => Dict(String(name) => value
                for (name, value) in diagnostics.body_lengths),
            "_qr_row_scales" => diagnostics.row_scales,
            "_qr_column_scales" => diagnostics.column_scales,
            "_position_corrections" => position_corrections,
            "_velocity_corrections" => velocity_corrections,
            "_imposed_initial_variables" => [string(
                layout.catalog.variables[index].component, ".",
                layout.catalog.variables[index].name)
                for index in sort!(collect(imposed_variables))],
            "_fallback_used" => fallback_used,
            "_redundant_velocity_equations" => [string(
                layout.catalog.equations[index].component, ".",
                layout.catalog.equations[index].name)
                for index in diagnostics.redundant_rows
                if layout.catalog.equations[index].kind == :constraint],
            "_selection_complete" => true)
        state_method == "preferred" && begin
            revised_state["_requested_preferred_velocities"] =
                String.(declared_preferred_names)
            revised_state["_dependent_condition"] = dependent_condition
        end
        revised["state_selection"] = revised_state
        return load_planar_document(revised, source;
            initial_override = initial,
            entered_initial_override = entered_initial_values,
            initial_conditions_override = initial_conditions,
            source_directory)
    end
    inactive_variables = Int[]
    inactive_equations = Int[]
    chosen, diagnostics = automatic_velocity_selection(
        model, bodies, body_names, markers, connections, connection_names,
        characteristic_length_overrides, initial, simulation.start_time)
    mechanical_degrees_of_freedom = length(chosen)
    inactive_variables, inactive_equations = constraint_family_indices(
        layout, diagnostics.redundant_rows)
    candidate_state_equations = Int[]
    selected_state_equations = Int[]
    for body in values(bodies)
        append!(candidate_state_equations,
            Iterators.flatten(values(body.candidate_state_equations)))
        append!(selected_state_equations,
            Iterators.flatten(body.state_equations))
    end
    for connection in values(connections)
        if connection isa PlanarRevoluteJointComponent ||
                connection isa PlanarDistanceCoordinateComponent
            append!(candidate_state_equations,
                connection.candidate_state_equations)
            append!(selected_state_equations, connection.state_equations)
        end
    end
    append!(inactive_equations, setdiff(candidate_state_equations,
        selected_state_equations))
    sort!(unique!(inactive_equations))
    selected_count == mechanical_degrees_of_freedom || begin
        mechanical_degrees_of_freedom == 0 && selected_count == 0 ||
            throw(ArgumentError("model has $mechanical_degrees_of_freedom independent states; " *
                "specify that many preferred velocities"))
    end
    active_variables = setdiff(collect(eachindex(layout.catalog.variables)),
                               inactive_variables)
    active_equations = setdiff(collect(eachindex(layout.catalog.equations)),
                               inactive_equations)
    length(active_variables) == length(active_equations) ||
        throw(ArgumentError("active dynamic system is not square"))
    auxiliary_state_count = count(variable -> variable.kind in
        (:user_state_hold, :user_state_steady), layout.catalog.variables)
    analysis = merge(analysis,
        (; degrees_of_freedom = mechanical_degrees_of_freedom +
            auxiliary_state_count))
    LoadedPlanarModel(title, layout, model, bodies, markers, connections,
        measures, forces, equation_components, drivers, analysis, simulation,
        parse_state_selection(document, valid_velocity_names),
        initial_conditions, weights, variable_weights,
        sort!(collect(imposed_variables)), entered_initial_values, initial,
        active_variables,
        active_equations, source)
end

end
