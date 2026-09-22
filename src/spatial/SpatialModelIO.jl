"""
    SpatialModelIO

TOML reader for the spatial rigid-body modeler. The initial vertical slice
supports rigid bodies, ground and body-fixed markers, spherical, revolute, and
fixed joints, perpendicular-axis, inplane, inline, hinge, and orient constraints, and
gravity, applied forces, joint-based applied torques, and rotational and
translational and spanning motion generators, spatial bushings, sphere-plane
contacts, plus span and directed-distance measurements, linear coordinate
couplers, ideal gear pairs, spur rack-and-pinion sets, and steady or transient
rolling tires.

Initial assembly uses rank-revealing projection and pivoted QR to remove
redundant scalar ideal-constraint families before state selection.
"""
module SpatialModelIO

using LinearAlgebra
using TOML
using ..InputUnits
using ..AssemblyExpansion: expand_model_assemblies, parse_lua_model
using ..ScalarExpressions: ScalarLaw, constant_scalar_law,
    compile_time_expression, compile_motion_expression
import ..ScalarExpressions
using ..AutomaticAnalysis
using ..SpatialComponentAssembly
using ..SpatialModeling
using ..SpatialDirectedDistances
using ..SpatialConstraints
using ..SpatialCoordinateCouplers
using ..SpatialGearPairs
using ..SpatialRackAndPinions
using ..SpatialSpans
using ..SpatialBelts
using ..SpatialAppliedForces
using ..SpatialBushings
using ..SpatialPlaneContacts
using ..SpatialCurveContacts
using ..PlanarCurveContacts: PlanarClosedCurve
using ..SpatialFrictionForces
using ..SpatialTires
using ..SpatialEquationComponents
using ..SpatialMotionGenerators
using ..SavedInitialConditions

export LoadedSpatialModel, load_spatial_model

const SPATIAL_CONFIGURATION_KINDS = Set((
    :position, :orientation_parameter, :relative_position, :elastic_position,
    :internal_state,
    :user_state_hold, :user_state_steady))
const SPATIAL_VELOCITY_KINDS = Set((
    :velocity, :angular_velocity, :relative_velocity, :elastic_velocity))
const SPATIAL_SAVED_STATE_ELEMENT_TYPES = Set((
    "rigid_body", "flexible_beam", "hinge", "revolute", "inline", "rolling_tire",
    "surface_friction", "revolute_friction", "translational_friction",
    "inplane_friction",
    "equation_component"))
const SPATIAL_ELEMENT_TYPES = Set(("ground", "rigid_body", "flexible_beam",
    "marker", "gravity",
    "applied_force", "applied_torque", "spanning_force", "spherical",
    "perp", "inplane", "inline", "hinge", "orient", "revolute",
    "fixed", "rotational_motion", "translational_motion",
    "spanning_motion", "span", "directed_distance", "bushing",
    "curve", "curve_contact", "flat_follower_contact",
    "plane_contact", "surface_friction", "revolute_friction",
    "translational_friction", "inplane_friction",
    "rolling_tire", "coupler", "gear_pair", "rack_and_pinion",
    "pulley", "belt", "belt_span", "equation_component"))

"""
    LoadedSpatialModel

Validated and allocated spatial model returned by [`load_spatial_model`](@ref).

The object retains the canonical layout, executable components, named model
objects, analysis settings, and both the entered and corrected initial states.
`body_reference_frames` records any reference-frame-to-CM transform.
`active_variable_indices` and `active_equation_indices` omit redundant ideal
constraint families and dormant candidate state equations without removing
them from the canonical layout. The expanded TOML in `model_source` is the
self-contained description stored in a result, including models entered
through Lua assemblies.
"""
struct LoadedSpatialModel
    title::String
    layout::ModelLayout
    model::ExecutableAnalysisModel
    bodies::Dict{Symbol,Any}
    body_reference_frames::Dict{Symbol,NamedTuple}
    markers::Dict{Symbol,Any}
    grounds::Set{Symbol}
    measures::Dict{Symbol,Any}
    forces::Dict{Symbol,Any}
    equation_components::Dict{Symbol,SpatialEquationComponent}
    connections::Dict{Symbol,Any}
    drivers::Dict{Symbol,Any}
    analysis::NamedTuple
    simulation::NamedTuple
    state_selection::NamedTuple
    initial_conditions::NamedTuple
    initial_condition_weights::Dict{Symbol,NamedTuple}
    initial_variable_weights::Vector{Float64}
    imposed_initial_variable_indices::Vector{Int}
    imposed_initial_orientations::Vector{Symbol}
    entered_initial_values::Vector{Float64}
    initial_values::Vector{Float64}
    active_variable_indices::Vector{Int}
    active_equation_indices::Vector{Int}
    model_source::String
end

const SPATIAL_EXPRESSION_KINEMATIC_KINDS = Set((
    :position, :orientation, :relative_position,
    :elastic_position, :velocity, :angular_velocity, :relative_velocity,
    :elastic_velocity))

spatial_expression_variable_supported(variable) =
    variable.kind in SPATIAL_EXPRESSION_KINEMATIC_KINDS ||
    variable.kind in (:user_algebraic, :user_state_hold,
        :user_state_steady) ||
    (variable.kind == :internal_state && variable.name in
        (:longitudinal_deformation, :lateral_deformation)) ||
    (variable.kind == :applied_geometry && variable.name in
        (:length, :gap, :deflection, :camber_angle, :station, :curvature,
         :contact_x, :contact_y, :contact_z,
         :normal_x, :normal_y, :normal_z)) ||
    (variable.kind == :applied_rate && variable.name in
        (:length_rate, :gap_rate, :station_rate, :deflection_rate,
         :forward_velocity,
         :lateral_velocity, :longitudinal_slip_velocity,
         :lateral_slip_velocity, :slip_ratio, :slip_angle)) ||
    (variable.kind == :applied_load && variable.name == :normal_force)

function compile_spatial_model_expression(source, parameters, layout)
    ScalarExpressions.compile_model_expression(source, parameters, layout,
        spatial_expression_variable_supported)
end

function spatial_linear_spanning_law(layout, name, stiffness, damping,
        free_length)
    variables = component_variable_indices(layout, name)
    dependencies = [variables[4], variables[8]]
    value = (t, z) -> -stiffness * (z[variables[4]] - free_length) -
        damping * z[variables[8]]
    gradient = (t, z) -> [-stiffness, -damping]
    ScalarLaw(value, gradient, dependencies)
end

function spatial_linear_tire_normal_law(layout, name, stiffness, damping)
    variables = component_variable_indices(layout, name)
    dependencies = [variables[1], variables[2]]
    value = (t, z) -> stiffness * z[variables[1]] +
        damping * z[variables[2]]
    gradient = (t, z) -> [stiffness, damping]
    ScalarLaw(value, gradient, dependencies)
end

function tire_load_curve(table, field, name)
    rows = get(table, field, nothing)
    rows isa Vector && length(rows) >= 2 || throw(ArgumentError(
        "rolling tire '$name'.$field requires at least two [load, value] rows"))
    curve = NTuple{2,Float64}[]
    for row in rows
        row isa Vector && length(row) == 2 || throw(ArgumentError(
            "rolling tire '$name'.$field rows must be [load, value]"))
        push!(curve, (finite_number(row[1], "rolling tire '$name'.$field load"),
            finite_number(row[2], "rolling tire '$name'.$field value")))
    end
    curve[1] == (0.0, 0.0) && all(curve[index][1] > curve[index - 1][1] &&
        curve[index][2] > curve[index - 1][2] for index in 2:length(curve)) ||
        throw(ArgumentError("rolling tire '$name'.$field must begin at " *
            "[0, 0] and increase strictly in load and value"))
    curve
end

function spatial_linear_torque_law(hinge, stiffness, damping, free_angle)
    _, omega, theta = hinge.rotation_variables
    dependencies = [theta, omega]
    value = (t, z) -> -stiffness * (z[theta] - free_angle[]) -
        damping * z[omega]
    gradient = (t, z) -> [-stiffness, -damping]
    ScalarLaw(value, gradient, dependencies)
end

function spatial_rotational_driver_laws(table, parameters, name)
    function_name = get(table, "function", "expression")
    if function_name == "constant_speed"
        angle_0 = angle_value(get(table, "initial_angle", 0.0),
            "rotational motion '$name'.initial_angle")
        haskey(table, "angular_velocity") || throw(ArgumentError(
            "rotational motion '$name' requires angular_velocity"))
        specification = table["angular_velocity"]
        omega = specification isa Number ? Float64(specification) :
            compile_time_expression(string(specification), parameters)(0.0)
        isfinite(omega) || throw(ArgumentError(
            "rotational motion '$name'.angular_velocity must be finite"))
        return (t -> angle_0 + omega * t, t -> omega, t -> 0.0)
    elseif function_name == "expression"
        haskey(table, "angle") || throw(ArgumentError(
            "rotational motion '$name' requires angle"))
        for field in ("angular_velocity", "angular_acceleration")
            haskey(table, field) && throw(ArgumentError(
                "rotational motion '$name'.$field is differentiated " *
                "automatically from angle and must be omitted"))
        end
        return compile_motion_expression(string(table["angle"]), parameters)
    end
    throw(ArgumentError(
        "unsupported rotational-motion function '$function_name'"))
end

function spatial_translational_driver_laws(table, parameters, name;
        label = "translational motion")
    function_name = get(table, "function", "expression")
    if function_name == "constant_speed"
        distance_0 = Float64(get(table, "initial_distance", 0.0))
        isfinite(distance_0) || throw(ArgumentError(
            "$label '$name'.initial_distance must be finite"))
        haskey(table, "velocity") || throw(ArgumentError(
            "$label '$name' requires velocity"))
        specification = table["velocity"]
        velocity = specification isa Number ? Float64(specification) :
            compile_time_expression(string(specification), parameters)(0.0)
        isfinite(velocity) || throw(ArgumentError(
            "$label '$name'.velocity must be finite"))
        return (t -> distance_0 + velocity * t,
            t -> velocity, t -> 0.0)
    elseif function_name == "expression"
        haskey(table, "distance") || throw(ArgumentError(
            "$label '$name' requires distance"))
        for field in ("velocity", "acceleration")
            haskey(table, field) && throw(ArgumentError(
                "$label '$name'.$field is differentiated automatically " *
                "from distance and must be omitted"))
        end
        return compile_motion_expression(
            string(table["distance"]), parameters)
    end
    throw(ArgumentError(
        "unsupported $label function '$function_name'"))
end

function spatial_spanning_driver_laws(table, parameters, name)
    spatial_translational_driver_laws(table, parameters, name;
        label = "spanning motion")
end

function source_text_and_directory(source)
    if source isa IO
        return read(source, String), pwd(), :toml
    elseif source isa AbstractString
        path = abspath(source)
        format = endswith(lowercase(path), ".lua") ? :lua : :toml
        return read(path, String), dirname(path), format
    end
    throw(ArgumentError("spatial model source must be a path or IO"))
end

function normalize_transferred_orientations!(initial, bodies,
        transferred_variables)
    transferred = Set(transferred_variables)
    for body in values(bodies)
        parameter_names = [Symbol(body.name, :., name)
            for name in (:p_0, :p_1, :p_2, :p_3)]
        present = [name in transferred for name in parameter_names]
        any(present) || continue
        all(present) || throw(ArgumentError(
            "saved initial-condition result has an incomplete orientation " *
            "for body '$(body.name)'"))
        parameters = @view initial[body.euler_parameter_variables]
        magnitude = norm(parameters)
        magnitude > sqrt(eps(Float64)) || throw(ArgumentError(
            "saved initial-condition orientation for body '$(body.name)' " *
            "has zero magnitude"))
        parameters ./= magnitude
    end
    initial
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

function numeric_vector(table, field, length_required; default = nothing,
        label = "table")
    value = get(table, field, default)
    isnothing(value) && throw(ArgumentError("$label requires $field"))
    value isa Vector && length(value) == length_required &&
        all(item -> item isa Number, value) || throw(ArgumentError(
            "$label.$field must contain $length_required numbers"))
    result = Float64.(value)
    all(isfinite, result) || throw(ArgumentError(
        "$label.$field must contain finite numbers"))
    result
end

function orientation_value(value, label)
    matrix = if value isa Vector && length(value) == 4 &&
            (value[1] isa Number || value[1] isa AbstractString) &&
            all(item -> item isa Number, value[2:4])
        angle = angle_value(value[1], "$label angle")
        axis_angle_rotation(angle, value[2:4])
    elseif value isa Vector && length(value) == 3 &&
            all(row -> row isa Vector && length(row) == 3 &&
                all(item -> item isa Number, row), value)
        reduce(vcat, permutedims.([Float64.(row) for row in value]))
    else
        throw(ArgumentError("$label must be [angle, axis_x, " *
            "axis_y, axis_z] or a 3 by 3 numeric matrix; angle may be " *
            "radians or a quoted degree value"))
    end
    matrix_to_euler_parameters(matrix)
    matrix
end

function orientation_matrix(table, field, label)
    value = get(table, field, nothing)
    isnothing(value) && return Matrix{Float64}(I, 3, 3)
    orientation_value(value, "$label.$field")
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

const SPATIAL_ANALYSIS_STAGES = (:static, :dynamic, :modal)

function active_during_stages(table, label)
    has_active = haskey(table, "active_during")
    has_inactive = haskey(table, "inactive_during")
    has_active && has_inactive && throw(ArgumentError(
        "$label may specify active_during or inactive_during, but not both"))
    !has_active && !has_inactive && return SPATIAL_ANALYSIS_STAGES

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
        return SPATIAL_ANALYSIS_STAGES
    end
    all(stage -> stage in SPATIAL_ANALYSIS_STAGES, stages) ||
        throw(ArgumentError("$label.$field stages must be static, " *
            "dynamic, or modal"))
    length(unique(stages)) == length(stages) || throw(ArgumentError(
        "$label.$field must not repeat a stage"))
    has_active ? Tuple(stages) :
        Tuple(stage for stage in SPATIAL_ANALYSIS_STAGES if stage ∉ stages)
end

function body_ic_weights(name, table)
    scale = finite_number(get(table, "ic_weight_scale", 1.0),
        "body '$name'.ic_weight_scale")
    scale > 0 || throw(ArgumentError(
        "body '$name'.ic_weight_scale must be positive"))
    (; scale)
end

function apply_body_initial_impose!(initial, imposed_variables,
        imposed_orientations, body, table, name)
    label = "body '$name'"
    initial_spec = initial_table(table, label)
    unknown_initial = setdiff(Set(keys(initial_spec)), Set(("impose",)))
    isempty(unknown_initial) || throw(ArgumentError(
        "$label.initial has unknown field '$(first(unknown_initial))'"))
    imposed = impose_table(table, label)
    indices = Dict(
        "R_x" => body.position_variables[1],
        "R_y" => body.position_variables[2],
        "R_z" => body.position_variables[3],
        "V_x" => body.velocity_variables[1],
        "V_y" => body.velocity_variables[2],
        "V_z" => body.velocity_variables[3],
        "omega_x" => body.angular_velocity_variables[1],
        "omega_y" => body.angular_velocity_variables[2],
        "omega_z" => body.angular_velocity_variables[3])
    allowed = Set([collect(keys(indices)); "orientation"])
    unknown = setdiff(Set(keys(imposed)), allowed)
    isempty(unknown) || throw(ArgumentError(
        "$label.initial.impose has unknown state '$(first(unknown))'"))
    for (field, value) in imposed
        if field == "orientation"
            matrix = orientation_value(value,
                "$label.initial.impose.orientation")
            initial[body.euler_parameter_variables] .=
                matrix_to_euler_parameters(matrix)
            push!(imposed_orientations, name)
        else
            index = indices[field]
            initial[index] = finite_number(value,
                "$label.initial.impose.$field")
            push!(imposed_variables, index)
        end
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
    coordinate_weight = finite_number(get(initial,
        coordinate * "_weight", 1.0),
        "$label.initial.$(coordinate)_weight")
    velocity_weight = finite_number(get(initial,
        velocity * "_weight", 1.0),
        "$label.initial.$(velocity)_weight")
    coordinate_weight > 0 && velocity_weight > 0 || throw(ArgumentError(
        "$label initial-condition weights must be positive"))
    coordinate_specified = any(haskey(initial, field)
        for field in (coordinate, coordinate * "_weight")) ||
        haskey(imposed, coordinate)
    velocity_specified = any(haskey(initial, field)
        for field in (velocity, velocity * "_weight")) ||
        haskey(imposed, velocity)
    (; coordinate, velocity, coordinate_guess, velocity_guess,
       coordinate_imposed, velocity_imposed, coordinate_weight,
       velocity_weight, coordinate_specified, velocity_specified)
end

function apply_relative_initial_specification!(initial, variable_weights,
        imposed_variables, coordinate_index, velocity_index, specification)
    !isnothing(specification.coordinate_guess) &&
        (initial[coordinate_index] = specification.coordinate_guess)
    !isnothing(specification.velocity_guess) &&
        (initial[velocity_index] = specification.velocity_guess)
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

function inertia_matrix(table, name)
    haskey(table, "inertia") || throw(ArgumentError(
        "body '$name' requires inertia"))
    value = table["inertia"]
    matrix = if value isa Vector && length(value) == 3 &&
            all(item -> item isa Number, value)
        Matrix(Diagonal(Float64.(value)))
    elseif value isa Vector && length(value) == 3 &&
            all(row -> row isa Vector && length(row) == 3 &&
                all(item -> item isa Number, row), value)
        reduce(vcat, permutedims.([Float64.(row) for row in value]))
    else
        throw(ArgumentError("body '$name'.inertia must contain three " *
            "principal values or a 3 by 3 matrix"))
    end
    all(isfinite, matrix) || throw(ArgumentError(
        "body '$name'.inertia must be finite"))
    issymmetric(matrix) || throw(ArgumentError(
        "body '$name'.inertia must be symmetric"))
    isposdef(Symmetric(matrix)) || throw(ArgumentError(
        "a free spatial body requires positive-definite inertia"))
    matrix
end

function simulation_settings(document)
    table = get(document, "simulation", Dict{String,Any}())
    table isa AbstractDict || throw(ArgumentError(
        "simulation must be a TOML table"))
    start_time = Float64(get(table, "start_time", 0.0))
    end_time = Float64(get(table, "end_time", 1.0))
    output_samples = Int(get(table, "output_samples", 201))
    output_precision = Symbol(get(table, "output_precision", "single"))
    relative_tolerance = Float64(get(table, "relative_tolerance", 1.0e-7))
    absolute_tolerance = Float64(get(table, "absolute_tolerance", 1.0e-9))
    initial_step = Float64(get(table, "initial_step", 1.0e-6))
    maximum_step = Float64(get(table, "maximum_step", 0.01))
    end_time >= start_time || throw(ArgumentError(
        "simulation.end_time must not precede start_time"))
    output_samples >= 1 || throw(ArgumentError(
        "simulation.output_samples must be positive"))
    output_precision in (:single, :double) || throw(ArgumentError(
        "simulation.output_precision must be single or double"))
    relative_tolerance > 0 && absolute_tolerance > 0 || throw(ArgumentError(
        "simulation tolerances must be positive"))
    initial_step > 0 && maximum_step > 0 || throw(ArgumentError(
        "simulation step sizes must be positive"))
    (; start_time, end_time, output_samples, output_precision,
       relative_tolerance,
       absolute_tolerance, initial_step, maximum_step)
end

function analysis_settings(document, degrees_of_freedom)
    table = get(document, "analysis", Dict{String,Any}())
    table isa AbstractDict || throw(ArgumentError(
        "analysis must be a TOML table"))
    haskey(table, "deficit") && throw(ArgumentError(
        "analysis.deficit is not a model option; the current StateSelected " *
        "formulation has deficit zero"))
    mode = Symbol(get(table, "mode", "automatic"))
    mode in (:automatic, :dynamic, :static, :modal) || throw(ArgumentError(
        "spatial analysis.mode must be automatic, dynamic, static, or modal"))
    initialization = Symbol(get(table, "initialization", "none"))
    initialization in (:none, :static_equilibrium) || throw(ArgumentError(
        "spatial analysis.initialization must be none or static_equilibrium"))
    initialization == :static_equilibrium && mode == :static &&
        throw(ArgumentError("static_equilibrium initialization requires " *
            "dynamic, modal, or automatic analysis"))
    static_method = Symbol(get(table, "static_method", "newton"))
    static_method in (:newton, :dynamic_relaxation) || throw(ArgumentError(
        "spatial analysis.static_method must be newton or dynamic_relaxation"))
    static_tolerance = haskey(table, "static_tolerance") ?
        Float64(table["static_tolerance"]) : nothing
    (isnothing(static_tolerance) ||
        (isfinite(static_tolerance) && static_tolerance > 0)) ||
        throw(ArgumentError(
            "spatial analysis.static_tolerance must be positive"))
    relaxation_duration = Float64(get(table, "relaxation_duration", 0.5))
    relaxation_duration > 0 || throw(ArgumentError(
        "analysis.relaxation_duration must be positive"))
    relaxation_reduction_factor =
        Float64(get(table, "relaxation_reduction_factor", 0.25))
    0 <= relaxation_reduction_factor < 1 || throw(ArgumentError(
        "analysis.relaxation_reduction_factor must be at least zero and less than one"))
    haskey(table, "relaxation_cycles") && throw(ArgumentError(
        "analysis.relaxation_cycles was renamed relaxation_max_cycles"))
    relaxation_min_cycles = Int(get(table, "relaxation_min_cycles", 1))
    relaxation_min_cycles >= 1 || throw(ArgumentError(
        "analysis.relaxation_min_cycles must be positive"))
    relaxation_max_cycles = Int(get(table, "relaxation_max_cycles", 20))
    relaxation_max_cycles >= relaxation_min_cycles || throw(ArgumentError(
        "analysis.relaxation_max_cycles must be at least relaxation_min_cycles"))
    relaxation_polish = get(table, "relaxation_polish", true)
    relaxation_polish isa Bool || throw(ArgumentError(
        "analysis.relaxation_polish must be true or false"))
    handoff_acceleration = Float64(get(table, "handoff_acceleration", 0.01))
    isfinite(handoff_acceleration) && handoff_acceleration > 0 ||
        throw(ArgumentError("analysis.handoff_acceleration must be positive"))
    handoff_speed = Float64(get(table, "handoff_speed", 0.01))
    isfinite(handoff_speed) && handoff_speed > 0 ||
        throw(ArgumentError("analysis.handoff_speed must be positive"))
    handoff_correction = Float64(get(table, "handoff_correction", 0.01))
    isfinite(handoff_correction) && handoff_correction > 0 ||
        throw(ArgumentError("analysis.handoff_correction must be positive"))
    number_of_modes = Int(get(table, "modes", 10))
    number_of_modes >= 1 || throw(ArgumentError(
        "analysis.modes must be positive"))
    frequency_shift_hz = Float64(get(table, "frequency_shift_hz", 0.0))
    frequency_shift_hz >= 0 || throw(ArgumentError(
        "analysis.frequency_shift_hz must be nonnegative"))
    modal_tolerance = Float64(get(table, "modal_tolerance", 1.0e-9))
    modal_tolerance > 0 || throw(ArgumentError(
        "analysis.modal_tolerance must be positive"))
    (; mode, initialization, static_method, static_tolerance,
       relaxation_duration,
       relaxation_reduction_factor, relaxation_min_cycles,
       relaxation_max_cycles, relaxation_polish,
       handoff_acceleration, handoff_speed,
       handoff_correction,
       number_of_modes, frequency_shift_hz, modal_tolerance,
       formulation = :StateSelected, deficit = 0, degrees_of_freedom)
end

function marker_owner(name, bodies, grounds)
    text = String(name)
    occursin('.', text) || throw(ArgumentError(
        "spatial marker '$name' must be nested beneath a body or ground"))
    candidates = [candidate for candidate in [collect(keys(bodies));
            collect(grounds)] if startswith(text, String(candidate) * ".")]
    isempty(candidates) && throw(ArgumentError(
        "spatial marker '$name' has no body or ground owner"))
    owner = argmax(candidate -> length(String(candidate)), candidates)
    if haskey(bodies, owner)
        return bodies[owner]
    elseif owner in grounds
        return nothing
    end
    throw(ArgumentError("spatial marker '$name' has no body or ground owner"))
end

function spatial_body_reference_frames(body_names, typed)
    references = Dict{Symbol,NamedTuple}()
    for name in body_names
        body_table = typed[name]
        specification = get(body_table, "center_of_mass", nothing)
        if isnothing(specification)
            references[name] = (;
                center_of_mass_marker = nothing,
                center_of_mass_position = zeros(3),
                center_of_mass_orientation = Matrix{Float64}(I, 3, 3))
            continue
        end
        specification isa AbstractString || throw(ArgumentError(
            "body '$name'.center_of_mass must name one of its markers"))
        marker_name = Symbol(specification)
        haskey(typed, marker_name) && typed[marker_name]["type"] == "marker" ||
            throw(ArgumentError(
                "body '$name'.center_of_mass references unknown marker " *
                "'$specification'"))
        candidates = [candidate for candidate in body_names
            if startswith(String(marker_name), String(candidate) * ".")]
        owner_name = isempty(candidates) ? nothing :
            argmax(candidate -> length(String(candidate)), candidates)
        owner_name == name || throw(ArgumentError(
            "body '$name'.center_of_mass marker '$specification' must belong " *
            "to that body"))
        marker_table = typed[marker_name]
        position = numeric_vector(marker_table, "position", 3;
            default = zeros(3), label = "marker '$marker_name'")
        orientation = orientation_matrix(
            marker_table, "orientation", "marker '$marker_name'")
        references[name] = (;
            center_of_mass_marker = marker_name,
            center_of_mass_position = position,
            center_of_mass_orientation = orientation)
    end
    references
end

function required_marker(markers, specification, label)
    specification isa AbstractString || throw(ArgumentError(
        "$label marker name must be a string"))
    name = Symbol(specification)
    haskey(markers, name) || throw(ArgumentError(
        "$label references unknown marker '$name'"))
    markers[name]
end

function imposed_detail(layout, imposed_variables, imposed_orientations)
    names = [string(variable.component, ".", variable.name)
        for variable in layout.catalog.variables
        if variable.index in imposed_variables]
    append!(names, [string(name, ".orientation")
        for name in imposed_orientations])
    isempty(names) ? "" : "; imposed values: " * join(sort!(names), ", ")
end

function weighted_constraint_projection!(values, model, variable_indices,
        equation_indices, variable_weights, imposed_variables, layout, label;
        time = 0.0)
    isempty(equation_indices) && return 0
    variables = [index for index in variable_indices
        if index ∉ imposed_variables]
    selection = AnalysisSelection(Dynamics(), variables,
        equation_indices)
    rates = zeros(length(values))
    equations = zeros(length(equation_indices))
    evaluate_analysis_equations!(equations, model, selection, time,
        values, rates)
    norm(equations, Inf) <= 1.0e-12 && return 0
    detail = imposed_detail(layout, imposed_variables, Set{Symbol}())
    isempty(variables) && throw(ArgumentError(
        "$label cannot satisfy its equations because all correction " *
        "variables are imposed$detail"))
    jacobian = Matrix(evaluate_analysis_jacobian(model, selection, time,
        values, rates, 0.0))
    inverse_sqrt_weights = 1 ./ sqrt.(variable_weights[variables])
    scaled_jacobian = jacobian * Diagonal(inverse_sqrt_weights)
    row_factorization = qr(transpose(scaled_jacobian), ColumnNorm())
    diagonal = abs.(diag(row_factorization.R))
    tolerance = isempty(diagonal) ? 0.0 :
        max(size(scaled_jacobian)...) * eps(Float64) * maximum(diagonal)
    numerical_rank = count(>(tolerance), diagonal)
    active_rows = collect(row_factorization.p)[1:numerical_rank]
    isempty(active_rows) && throw(ArgumentError(
        "$label has no independent correction equations$detail"))
    factorization = svd(scaled_jacobian[active_rows, :])
    correction = -inverse_sqrt_weights .*
        (factorization \ equations[active_rows])
    values[variables] .+= correction
    evaluate_analysis_equations!(equations, model, selection, time,
        values, rates)
    norm(equations, Inf) <= 1.0e-10 || throw(ArgumentError(
        "$label could not be made consistent"))
    1
end

function apply_spatial_position_correction!(values, bodies, descriptors,
        correction, body_lengths)
    rotations = Dict(name => zeros(3) for name in keys(bodies))
    for (column, descriptor) in enumerate(descriptors)
        if descriptor.kind == :rotation
            rotations[descriptor.body][descriptor.component] =
                correction[column] / body_lengths[descriptor.body]
        else
            values[descriptor.variable] += correction[column]
        end
    end
    for (name, rotation_correction) in rotations
        angle = norm(rotation_correction)
        iszero(angle) && continue
        body = bodies[name]
        parameters = @view values[body.euler_parameter_variables]
        orientation = rotation_matrix(parameters)
        increment = axis_angle_rotation(angle, rotation_correction)
        values[body.euler_parameter_variables] .=
            matrix_to_euler_parameters(orientation * increment)
    end
    nothing
end

"""Correct body origins and orientations against all position constraints."""
function spatial_position_projection!(values, model, bodies, markers,
        equation_indices, relative_coordinate_indices, variable_weights,
        body_weights, imposed_variables, imposed_orientations, layout, label;
        maximum_iterations = 20, time = 0.0)
    isempty(equation_indices) && return 0
    body_names = sort!(collect(keys(bodies)))
    body_lengths = spatial_body_lengths(bodies, markers)
    variable_indices = reduce(vcat,
        [[collect(bodies[name].position_variables);
          collect(bodies[name].euler_parameter_variables);
          (bodies[name] isa SpatialFlexibleBeamComponent ?
              collect(bodies[name].elastic_position_variables) : Int[])]
         for name in body_names]; init = Int[])
    append!(variable_indices, relative_coordinate_indices)
    selection = AnalysisSelection(Dynamics(), variable_indices,
        equation_indices)
    local_rows = Dict(index => row
        for (row, index) in enumerate(variable_indices))
    descriptors = NamedTuple[]
    correction_weights = Float64[]
    for name in body_names
        body = bodies[name]
        weights = body_weights[name]
        for component in 1:3
            variable = body.position_variables[component]
            variable in imposed_variables && continue
            push!(descriptors, (; kind = :translation, body = name,
                component, variable))
            push!(correction_weights, weights.weight)
        end
        if name ∉ imposed_orientations
            for component in 1:3
                push!(descriptors, (; kind = :rotation, body = name,
                    component, variable = 0))
                push!(correction_weights, weights.weight)
            end
        end
        if body isa SpatialFlexibleBeamComponent
            for variable in body.elastic_position_variables
                variable in imposed_variables && continue
                push!(descriptors, (; kind = :elastic, body = name,
                    component = 0, variable))
                push!(correction_weights, variable_weights[variable])
            end
        end
    end
    for variable in relative_coordinate_indices
        variable in imposed_variables && continue
        push!(descriptors, (; kind = :relative, body = Symbol(),
            component = 0, variable))
        push!(correction_weights, variable_weights[variable])
    end
    detail = imposed_detail(layout, imposed_variables, imposed_orientations)
    rates = zeros(length(values))
    equations = zeros(length(equation_indices))
    for iteration in 0:maximum_iterations
        evaluate_analysis_equations!(equations, model, selection, time,
            values, rates)
        norm(equations, Inf) <= 1.0e-12 && return iteration
        iteration == maximum_iterations && break

        canonical_jacobian = Matrix(evaluate_analysis_jacobian(model,
            selection, time, values, rates, 0.0))
        isempty(descriptors) && throw(ArgumentError(
            "$label cannot satisfy its equations because all correction " *
            "variables are imposed$detail"))
        correction_map = zeros(eltype(values), length(variable_indices),
            length(descriptors))
        for (column, descriptor) in enumerate(descriptors)
            if descriptor.kind == :rotation
                body = bodies[descriptor.body]
                parameters = @view values[body.euler_parameter_variables]
                parameter_rows = local_rows[
                    body.euler_parameter_variables.start]:local_rows[
                    body.euler_parameter_variables.stop]
                correction_map[parameter_rows, column] .=
                    0.5 .* quaternion_rate_matrix(parameters)[:,
                        descriptor.component] ./ body_lengths[descriptor.body]
            else
                correction_map[local_rows[descriptor.variable], column] = 1.0
            end
        end
        inverse_sqrt_weights = 1 ./ sqrt.(correction_weights)
        jacobian = canonical_jacobian * correction_map *
            Diagonal(inverse_sqrt_weights)
        row_factorization = qr(transpose(jacobian), ColumnNorm())
        diagonal = abs.(diag(row_factorization.R))
        tolerance = isempty(diagonal) ? 0.0 :
            max(size(jacobian)...) * eps(Float64) * maximum(diagonal)
        numerical_rank = count(>(tolerance), diagonal)
        active_rows = collect(row_factorization.p)[1:numerical_rank]
        isempty(active_rows) && break
        factorization = svd(jacobian[active_rows, :])
        correction = -inverse_sqrt_weights .*
            (factorization \ equations[active_rows])
        apply_spatial_position_correction!(values, bodies, descriptors,
            correction, body_lengths)
    end
    throw(ArgumentError("$label could not be made consistent"))
end

function spatial_body_lengths(bodies, markers)
    result = Dict{Symbol,Float64}()
    for (name, body) in bodies
        offsets = [norm(marker.position_body) for marker in values(markers)
            if marker isa Union{SpatialBodyMarker,SpatialFlexibleBeamMarker} &&
                marker.body === body &&
                norm(marker.position_body) > 0]
        fallback = body.mass > 0 ?
            sqrt(maximum(diag(body.inertia)) / body.mass) : 1.0
        result[name] = isempty(offsets) ? max(fallback, 1.0) :
            maximum(offsets)
    end
    result
end

"""
Select independent spatial velocities from the scaled constraint partial.

Column-pivoted QR identifies velocity columns not required to satisfy the
constraints; these become states. A separate pivoted QR of the transpose,
restricted to the dependent columns, identifies independent and redundant
constraint rows. The returned diagnostics preserve the scaling and pivots for
preferred-state validation and result reporting.
"""
function automatic_spatial_velocity_selection(model, velocity_indices,
        velocity_rows, column_scales, canonical, time = 0.0)
    if isempty(velocity_rows)
        return copy(velocity_indices),
            (; rank = 0, tolerance = 0.0,
               scaled_partial = zeros(0, length(velocity_indices)),
               pivots = collect(eachindex(velocity_indices)),
               independent_rows = Int[], redundant_rows = Int[],
               row_scales = Float64[])
    end
    selection = AnalysisSelection(VelocityIC(), velocity_indices,
        velocity_rows)
    derivative = zeros(length(canonical))
    partial = Matrix(evaluate_analysis_jacobian(model, selection, time,
        canonical, derivative, 0.0))
    column_scaled = partial * Diagonal(column_scales)
    row_norms = [norm(row) for row in eachrow(column_scaled)]
    row_scales = [value > 0 ? inv(value) : 1.0 for value in row_norms]
    scaled_partial = Diagonal(row_scales) * column_scaled
    factorization = qr(scaled_partial, ColumnNorm())
    diagonal = abs.(diag(factorization.R))
    tolerance = isempty(diagonal) ? 0.0 :
        max(size(scaled_partial)...) * eps(Float64) * maximum(diagonal)
    rank_value = count(>(tolerance), diagonal)
    pivots = collect(factorization.p)
    selected = velocity_indices[pivots[(rank_value + 1):end]]
    independent_rows, redundant_rows = if rank_value == 0
        (Int[], copy(velocity_rows))
    else
        dependent = pivots[1:rank_value]
        row_factorization = qr(
            transpose(scaled_partial[:, dependent]), ColumnNorm())
        row_pivots = collect(row_factorization.p)
        (velocity_rows[row_pivots[1:rank_value]],
         velocity_rows[row_pivots[(rank_value + 1):end]])
    end
    selected, (; rank = rank_value, tolerance, scaled_partial,
        pivots, independent_rows, redundant_rows, row_scales)
end

function spatial_constraint_family_indices(layout, velocity_equation_indices)
    inactive_variables = Int[]
    inactive_equations = Int[]
    for velocity_index in velocity_equation_indices
        velocity = layout.catalog.equations[velocity_index]
        velocity.kind == :constraint || throw(ArgumentError(
            "spatial equation $(velocity.component).$(velocity.name) is " *
            "redundant but does not belong to an ideal constraint family"))
        matches = [family for family in values(layout.constraint_families)
            if family.velocity_equation == velocity_index]
        length(matches) == 1 || error(
            "Constraint velocity equation $(velocity.component)." *
            "$(velocity.name) must belong to exactly one scalar " *
            "constraint family")
        family = only(matches)
        append!(inactive_equations, (family.position_equation,
            family.velocity_equation, family.acceleration_equation))
        push!(inactive_variables, family.reaction_variable)
    end
    sort!(unique!(inactive_variables)), sort!(unique!(inactive_equations))
end

function state_selection(model, layout, bodies, markers, connections,
        drivers, canonical, document)
    velocity_indices = Int[]
    column_scales = Float64[]
    coordinate_indices = Dict{Int,Int}()
    state_equations = Dict{Int,Vector{Int}}()
    body_lengths = spatial_body_lengths(bodies, markers)
    for name in sort!(collect(keys(bodies)))
        body = bodies[name]
        for component in 1:3
            velocity = body.velocity_variables[component]
            push!(velocity_indices, velocity)
            push!(column_scales, body_lengths[name])
            coordinate_indices[velocity] = body.position_variables[component]
            state_equations[velocity] = [
                body.acceleration_state_equations[component],
                body.position_state_equations[component]]
        end
        for component in 1:3
            velocity = body.angular_velocity_variables[component]
            push!(velocity_indices, velocity)
            push!(column_scales, 1.0)
            coordinate_indices[velocity] = body.pseudo_angle_variables[component]
            state_equations[velocity] = [
                body.angular_acceleration_state_equations[component],
                body.pseudo_angle_state_equations[component]]
        end
        if body isa SpatialFlexibleBeamComponent
            for component in 1:6
                velocity = body.elastic_velocity_variables[component]
                push!(velocity_indices, velocity)
                push!(column_scales, component <= 3 ? body_lengths[name] : 1.0)
                coordinate_indices[velocity] =
                    body.elastic_position_variables[component]
                state_equations[velocity] = [
                    body.elastic_acceleration_state_equations[component],
                    body.elastic_position_state_equations[component]]
            end
        end
    end
    for name in sort!(collect(keys(connections)))
        hinge = connection_hinge(connections[name])
        if !isnothing(hinge) && !isempty(hinge.rotation_variables)
            alpha, velocity, coordinate = hinge.rotation_variables
            push!(velocity_indices, velocity)
            push!(column_scales, 1.0)
            coordinate_indices[velocity] = coordinate
            state_equations[velocity] = copy(hinge.state_equations)
        end
        inline = connection_inline(connections[name])
        if !isnothing(inline) && !isempty(inline.translation_variables)
            acceleration, velocity, coordinate =
                inline.translation_variables
            push!(velocity_indices, velocity)
            push!(column_scales, maximum(values(body_lengths)))
            coordinate_indices[velocity] = coordinate
            state_equations[velocity] = copy(inline.state_equations)
        end
    end
    velocity_rows = [equation.index for equation in layout.catalog.equations
        if equation.level == 1 && equation.kind == :constraint]
    for connection in values(connections)
        hinge = connection_hinge(connection)
        !isnothing(hinge) && !isempty(hinge.rotation_equations) &&
            push!(velocity_rows, hinge.rotation_equations[2])
        inline = connection_inline(connection)
        !isnothing(inline) && !isempty(inline.translation_equations) &&
            push!(velocity_rows, inline.translation_equations[2])
    end
    for driver in values(drivers)
        if driver isa SpatialRotationalMotionGenerator
            push!(velocity_rows, driver.velocity_equation)
        elseif driver isa SpatialTranslationalMotionGenerator
            push!(velocity_rows, first(driver.velocity_equations))
        elseif driver isa SpatialSpanningMotionGenerator
            push!(velocity_rows, driver.span.velocity_equation)
        end
    end
    sort!(unique!(velocity_rows))
    automatic_selected, diagnostics = automatic_spatial_velocity_selection(
        model, velocity_indices, velocity_rows, column_scales, canonical)
    rank_value = diagnostics.rank
    tolerance = diagnostics.tolerance
    scaled_partial = diagnostics.scaled_partial
    redundant_variables, redundant_equations =
        spatial_constraint_family_indices(layout,
            diagnostics.redundant_rows)
    active_velocity_rows = setdiff(
        velocity_rows, diagnostics.redundant_rows)
    candidate_names = [Symbol(layout.catalog.variables[index].component, ".",
        layout.catalog.variables[index].name) for index in velocity_indices]
    settings = get(document, "state_selection", Dict{String,Any}())
    settings isa AbstractDict || throw(ArgumentError(
        "state_selection must be a TOML table"))
    method = String(get(settings, "method", "automatic"))
    method in ("automatic", "preferred") || throw(ArgumentError(
        "spatial state_selection.method must be 'automatic' or 'preferred'"))
    requested_names = Symbol.(get(settings, "preferred_velocities", String[]))
    allow_fallback = Bool(get(settings, "allow_fallback", true))
    maximum_condition_number = Float64(get(
        settings, "maximum_condition_number", 1.0e8))
    maximum_condition_number > 0 || throw(ArgumentError(
        "state_selection.maximum_condition_number must be positive"))
    selected = automatic_selected
    fallback_used = false
    dependent_condition = nothing
    if method == "preferred"
        length(unique(requested_names)) == length(requested_names) ||
            throw(ArgumentError("preferred velocities must be unique"))
        expected = length(velocity_indices) - rank_value
        reason = nothing
        preferred_columns = Int[]
        length(requested_names) > expected &&
            (reason = "at most $expected preferred velocities may be " *
                "specified, got $(length(requested_names))")
        for name in requested_names
            column = findfirst(==(name), candidate_names)
            isnothing(column) && throw(ArgumentError(
                "unknown preferred velocity '$name'"))
            push!(preferred_columns, column)
        end
        completed_columns = Int[]
        if isnothing(reason)
            available_columns = setdiff(
                collect(eachindex(velocity_indices)), preferred_columns)
            available_factorization = qr(
                scaled_partial[:, available_columns], ColumnNorm())
            available_pivots = collect(available_factorization.p)
            available_diagonal = abs.(diag(available_factorization.R))
            available_tolerance = isempty(available_diagonal) ? 0.0 :
                max(size(scaled_partial)...) * eps(Float64) *
                    maximum(available_diagonal)
            available_rank = count(>(available_tolerance),
                available_diagonal)
            available_rank == rank_value ||
                (reason = "the requested velocities cannot all be states")
            dependent_columns = isnothing(reason) ?
                available_columns[available_pivots[1:rank_value]] : Int[]
            completed_columns = isnothing(reason) ?
                vcat(preferred_columns,
                    setdiff(available_columns, dependent_columns)) : Int[]
            block = scaled_partial[:, dependent_columns]
            singular_values = svdvals(block)
            block_tolerance = isempty(singular_values) ? 0.0 :
                max(size(block)...) * eps(Float64) * maximum(singular_values)
            block_rank = count(>(block_tolerance), singular_values)
            condition = isempty(singular_values) ? 1.0 :
                (last(singular_values) == 0 ? Inf :
                 first(singular_values) / last(singular_values))
            dependent_condition = condition
            block_rank == rank_value ||
                (reason = "the complementary dependent block is rank deficient")
            isnothing(reason) && condition > maximum_condition_number &&
                (reason = "the complementary dependent block has condition " *
                    "number $condition")
        end
        if isnothing(reason)
            selected = velocity_indices[completed_columns]
        elseif allow_fallback
            fallback_used = true
        else
            throw(ArgumentError(
                "preferred state selection is unusable: $reason"))
        end
    end
    selected_equations = reduce(vcat,
        [state_equations[index] for index in selected]; init = Int[])
    candidate_equations = reduce(vcat, collect(Base.values(state_equations));
        init = Int[])
    inactive_equations = setdiff(candidate_equations, selected_equations)
    inactive_variables = Int[]
    append!(inactive_variables, redundant_variables)
    append!(inactive_equations, redundant_equations)
    sort!(unique!(inactive_variables))
    sort!(unique!(inactive_equations))
    selected_names = [Symbol(layout.catalog.variables[index].component, ".",
        layout.catalog.variables[index].name) for index in selected]
    selected_coordinates = [coordinate_indices[index] for index in selected]
    redundant_velocity_equations = [Symbol(
        layout.catalog.equations[index].component, ".",
        layout.catalog.equations[index].name)
        for index in diagnostics.redundant_rows]
    (; velocity_indices, velocity_rows = active_velocity_rows,
       coordinate_indices, state_equations,
       candidate_names, selected_velocity_indices = selected,
       selected_coordinate_indices = selected_coordinates,
       selected_velocities = selected_names, inactive_variables,
       inactive_equations,
       rank = rank_value, tolerance, body_lengths, column_scales,
       qr_pivots = diagnostics.pivots,
       independent_velocity_equations = diagnostics.independent_rows,
       redundant_velocity_equations,
       qr_row_scales = diagnostics.row_scales,
       method = Symbol(method), requested_velocities = requested_names,
       allow_fallback, maximum_condition_number, fallback_used,
       dependent_condition)
end

function initialize_spatial_accelerations!(values, model, layout, time = 0.0;
        inactive_variables = Int[], inactive_equations = Int[])
    complete = select_analysis(layout.catalog, AccelerationIC())
    relative_accelerations = [variable.index
        for variable in layout.catalog.variables
        if variable.kind == :relative_acceleration]
    relative_acceleration_equations = [equation.index
        for equation in layout.catalog.equations
        if equation.kind == :coordinate_relation && equation.level == 2]
    motion_acceleration_equations = [equation.index
        for equation in layout.catalog.equations
        if equation.kind == :motion && equation.level == 2]
    variable_indices = setdiff(
        [complete.variable_indices; relative_accelerations],
        inactive_variables)
    equation_indices = setdiff(
        [complete.equation_indices; relative_acceleration_equations;
         motion_acceleration_equations], inactive_equations)
    selection = AnalysisSelection(AccelerationIC(), variable_indices,
        equation_indices)
    length(selection.variable_indices) == length(selection.equation_indices) ||
        throw(ArgumentError("spatial acceleration initialization is not square"))
    derivatives = zeros(length(values))
    equations = zeros(length(selection.equation_indices))
    evaluate_analysis_equations!(equations, model, selection, time,
        values, derivatives)
    initial_norm = norm(equations, Inf)
    initial_norm <= 1.0e-12 && return 0
    residual_norm = initial_norm
    for iteration in 1:6
        jacobian = evaluate_analysis_sparse_jacobian(model, selection, time,
            values, derivatives, 0.0)
        correction = try
            -(jacobian \ equations)
        catch error
            error isa LinearAlgebra.SingularException || rethrow()
            throw(ArgumentError(
                "spatial acceleration initialization is singular"))
        end
        values[selection.variable_indices] .+= correction
        evaluate_analysis_equations!(equations, model, selection, time,
            values, derivatives)
        residual_norm = norm(equations, Inf)
        residual_norm <= 1.0e-10 && return iteration
    end
    throw(ArgumentError(
        "spatial acceleration initialization did not converge; " *
        "equation norm fell from $initial_norm to $residual_norm"))
end

"""
    load_spatial_model(source; format=nothing, source_directory=nothing,
        source_label=nothing)

Load a spatial model from TOML, Lua, or an `AbstractDict` model document. A
path determines its format and directory. The optional keywords let an
uploaded input stream retain a declared format, assembly search directory,
and diagnostic filename. Dictionary keys and symbolic values are normalized
to the string representation used by TOML before ordinary model validation.
"""
function normalized_model_document_value(value::AbstractDict)
    Dict{String,Any}(string(key) => normalized_model_document_value(item)
        for (key, item) in value)
end

normalized_model_document_value(value::NamedTuple) =
    normalized_model_document_value(Dict(pairs(value)))
normalized_model_document_value(value::Tuple) =
    normalized_model_document_value(collect(value))
normalized_model_document_value(value::AbstractVector) =
    Any[normalized_model_document_value(item) for item in value]
normalized_model_document_value(value::AbstractMatrix) =
    [Any[normalized_model_document_value(value[row, column])
         for column in axes(value, 2)] for row in axes(value, 1)]
normalized_model_document_value(value::Symbol) = String(value)
normalized_model_document_value(value) = value

function load_spatial_model(source; format = nothing,
        source_directory = nothing, source_label = nothing)
    is_document = source isa AbstractDict
    is_document && !isnothing(format) && throw(ArgumentError(
        "format is not used when loading a spatial model document"))
    text, detected_directory, detected_format = if is_document
        ("", pwd(), :document)
    else
        source_text_and_directory(source)
    end
    source_format = isnothing(format) ? detected_format : Symbol(format)
    source_format in (:toml, :lua, :document) || throw(ArgumentError(
        "spatial model format must be :toml or :lua"))
    model_directory = isnothing(source_directory) ? detected_directory :
        abspath(String(source_directory))
    label = isnothing(source_label) ?
        (is_document ? "Julia model document" : string(source)) :
        String(source_label)
    document = if source_format == :document
        normalized_model_document_value(source)
    elseif source_format == :lua
        parse_lua_model(text;
            label = "Lua model '$label'", source_directory = model_directory,
            element_types = SPATIAL_ELEMENT_TYPES)
    else
        TOML.parse(text)
    end
    model_table = get(document, "model", nothing)
    model_table isa AbstractDict || throw(ArgumentError(
        "spatial model requires a [model] table"))
    get(model_table, "dimension", "") == "spatial" || throw(ArgumentError(
        "model.dimension must be 'spatial'"))
    document = expand_model_assemblies(document, model_directory;
        dimension = "spatial")
    text = sprint(io -> TOML.print(io, document))
    model_table = document["model"]
    title = String(get(model_table, "title",
        get(model_table, "name", "Spatial model")))
    parameters = get(document, "parameters", Dict{String,Any}())
    parameters isa AbstractDict || throw(ArgumentError(
        "parameters must be a TOML table"))
    for (name, value) in parameters
        value isa Number && isfinite(value) || throw(ArgumentError(
            "parameter '$name' must be a finite number"))
    end
    simulation = simulation_settings(document)
    typed = collect_typed_tables!(Dict{Symbol,Any}(), document)
    for (name, table) in typed
        kind = String(table["type"])
        kind in SPATIAL_ELEMENT_TYPES || throw(ArgumentError(
            "spatial element '$name' has unsupported type '$kind'"))
    end
    kinds = Dict(name => String(table["type"]) for (name, table) in typed)
    body_names = sort!([name for (name, table) in typed
        if table["type"] in ("rigid_body", "flexible_beam")])
    beam_names = Set(name for name in body_names
        if typed[name]["type"] == "flexible_beam")
    body_reference_frames = spatial_body_reference_frames(body_names, typed)
    spherical_names = sort!([name for (name, table) in typed
        if table["type"] == "spherical"])
    perp_names = sort!([name for (name, table) in typed
        if table["type"] == "perp"])
    inplane_names = sort!([name for (name, table) in typed
        if table["type"] == "inplane"])
    inline_names = sort!([name for (name, table) in typed
        if table["type"] == "inline"])
    hinge_names = sort!([name for (name, table) in typed
        if table["type"] == "hinge"])
    orient_names = sort!([name for (name, table) in typed
        if table["type"] == "orient"])
    revolute_names = sort!([name for (name, table) in typed
        if table["type"] == "revolute"])
    fixed_names = sort!([name for (name, table) in typed
        if table["type"] == "fixed"])
    coupler_names = sort!([name for (name, table) in typed
        if table["type"] == "coupler"])
    gear_pair_names = sort!([name for (name, table) in typed
        if table["type"] == "gear_pair"])
    rack_and_pinion_names = sort!([name for (name, table) in typed
        if table["type"] == "rack_and_pinion"])
    pulley_names = sort!([name for (name, table) in typed
        if table["type"] == "pulley"])
    belt_names = sort!([name for (name, table) in typed
        if table["type"] == "belt"])
    belt_span_names = sort!([name for (name, table) in typed
        if table["type"] == "belt_span"])
    spanning_force_names = sort!([name for (name, table) in typed
        if table["type"] == "spanning_force"])
    applied_force_names = sort!([name for (name, table) in typed
        if table["type"] == "applied_force"])
    applied_torque_names = sort!([name for (name, table) in typed
        if table["type"] == "applied_torque"])
    bushing_names = sort!([name for (name, table) in typed
        if table["type"] == "bushing"])
    curve_names = sort!([name for (name, table) in typed
        if table["type"] == "curve"])
    curve_contact_names = sort!([name for (name, table) in typed
        if table["type"] in ("curve_contact", "flat_follower_contact")])
    plane_contact_names = sort!([name for (name, table) in typed
        if table["type"] == "plane_contact"])
    surface_friction_names = sort!([name for (name, table) in typed
        if table["type"] == "surface_friction"])
    revolute_friction_names = sort!([name for (name, table) in typed
        if table["type"] == "revolute_friction"])
    translational_friction_names = sort!([name for (name, table) in typed
        if table["type"] == "translational_friction"])
    inplane_friction_names = sort!([name for (name, table) in typed
        if table["type"] == "inplane_friction"])
    tire_names = sort!([name for (name, table) in typed
        if table["type"] == "rolling_tire"])
    equation_component_names = sort!([name for (name, table) in typed
        if table["type"] == "equation_component"])
    transient_tires = Dict{Symbol,Bool}()
    for name in tire_names
        table = typed[name]
        tangential_model = lowercase(string(get(table, "tangential_model",
            "expression")))
        tangential_model in ("expression", "bristle") ||
            throw(ArgumentError("rolling tire '$name'.tangential_model must " *
                "be 'expression' or 'bristle'"))
        has_longitudinal = haskey(table, "longitudinal_relaxation_length")
        has_lateral = haskey(table, "lateral_relaxation_length")
        has_longitudinal == has_lateral || throw(ArgumentError(
            "rolling tire '$name' must specify both longitudinal and " *
            "lateral relaxation lengths, or neither"))
        tangential_model == "bristle" && has_longitudinal && throw(
            ArgumentError("rolling tire '$name' bristle mode cannot use " *
                "expression relaxation lengths"))
        transient_tires[name] = has_longitudinal || tangential_model == "bristle"
    end
    rotational_motion_names = sort!([name for (name, table) in typed
        if table["type"] == "rotational_motion"])
    translational_motion_names = sort!([name for (name, table) in typed
        if table["type"] == "translational_motion"])
    spanning_motion_names = sort!([name for (name, table) in typed
        if table["type"] == "spanning_motion"])
    span_measure_names = sort!([name for (name, table) in typed
        if table["type"] == "span"])
    directed_distance_names = sort!([name for (name, table) in typed
        if table["type"] == "directed_distance"])
    measure_names = sort!([span_measure_names; directed_distance_names])
    torque_joint_names = Dict{Symbol,Symbol}()
    connection_name_set = Set([hinge_names; revolute_names])
    for name in applied_torque_names
        joint_value = get(typed[name], "joint", nothing)
        joint_value isa AbstractString || throw(ArgumentError(
            "applied torque '$name'.joint must name a hinge or revolute"))
        joint = Symbol(joint_value)
        joint in connection_name_set || throw(ArgumentError(
            "applied torque '$name' names unknown hinge or revolute " *
            "'$joint_value'"))
        torque_joint_names[name] = joint
    end
    driver_joint_names = Dict{Symbol,Symbol}()
    for name in rotational_motion_names
        joint_value = get(typed[name], "joint", nothing)
        joint_value isa AbstractString || throw(ArgumentError(
            "rotational motion '$name'.joint must name a hinge or revolute"))
        joint = Symbol(joint_value)
        joint in connection_name_set || throw(ArgumentError(
            "rotational motion '$name' names unknown hinge or revolute " *
            "'$joint_value'"))
        joint in values(driver_joint_names) && throw(ArgumentError(
            "hinge or revolute '$joint' has more than one rotational motion"))
        haskey(typed[joint], "initial") && throw(ArgumentError(
            "driven hinge or revolute '$joint' cannot also specify initial values"))
        driver_joint_names[name] = joint
    end
    pulley_joint_names = Dict{Symbol,Symbol}()
    for name in pulley_names
        joint_value = get(typed[name], "joint", nothing)
        joint_value isa AbstractString || throw(ArgumentError(
            "pulley '$name'.joint must name a revolute joint"))
        joint = Symbol(joint_value)
        joint in Set(revolute_names) || throw(ArgumentError(
            "pulley '$name' names unknown revolute joint '$joint_value'"))
        pulley_joint_names[name] = joint
    end
    rotation_coordinate_names = Set(name for name in
            [hinge_names; revolute_names] if begin
        setting = get(typed[name], "rotation_coordinates", false)
        setting isa Bool || throw(ArgumentError(
            "'$name'.rotation_coordinates must be Boolean"))
        setting
    end)
    union!(rotation_coordinate_names, values(torque_joint_names))
    union!(rotation_coordinate_names, values(driver_joint_names))
    union!(rotation_coordinate_names, values(pulley_joint_names))
    translation_coordinate_names = Set(name for name in inline_names if begin
        setting = get(typed[name], "translation_coordinates", false)
        setting isa Bool || throw(ArgumentError(
            "'$name'.translation_coordinates must be Boolean"))
        setting
    end)
    rack_and_pinion_joint_names = Dict{Symbol,Tuple{Symbol,Symbol}}()
    for name in rack_and_pinion_names
        raw_joints = get(typed[name], "joints", nothing)
        raw_joints isa Vector && length(raw_joints) == 2 &&
            all(joint -> joint isa AbstractString, raw_joints) ||
            throw(ArgumentError("rack and pinion '$name'.joints must " *
                "contain an inline constraint and a revolute joint"))
        inline_name, revolute_name = Symbol.(raw_joints)
        inline_name in inline_names || throw(ArgumentError(
            "rack and pinion '$name' first joint must be an inline constraint"))
        revolute_name in revolute_names || throw(ArgumentError(
            "rack and pinion '$name' second joint must be a revolute joint"))
        rack_and_pinion_joint_names[name] = (inline_name, revolute_name)
        push!(translation_coordinate_names, inline_name)
        push!(rotation_coordinate_names, revolute_name)
    end
    coupler_coordinate_specs = Dict{Symbol,Vector{Tuple{Symbol,Symbol}}}()
    for name in coupler_names
        raw_specs = get(typed[name], "coordinates", Any[])
        length(raw_specs) >= 2 || throw(ArgumentError(
            "coupler '$name' requires at least two coordinates"))
        all(specification -> specification isa AbstractString, raw_specs) ||
            throw(ArgumentError(
                "coupler '$name' coordinates must be strings"))
        length(unique(raw_specs)) == length(raw_specs) || throw(ArgumentError(
            "coupler '$name' coordinates must be unique"))
        parsed = Tuple{Symbol,Symbol}[]
        for raw_specification in raw_specs
            specification = String(raw_specification)
            parts = split(specification, ".")
            length(parts) >= 2 || throw(ArgumentError(
                "coupler '$name' coordinate '$specification' must name an element port"))
            port = Symbol(last(parts))
            component = Symbol(join(parts[1:end-1], "."))
            if port == :rotation &&
                    (component in hinge_names || component in revolute_names)
                push!(rotation_coordinate_names, component)
            elseif port == :distance && component in inline_names
                push!(translation_coordinate_names, component)
            else
                throw(ArgumentError(
                    "coupler '$name' coordinate '$specification' is not an available coordinate port"))
            end
            push!(parsed, (component, port))
        end
        coupler_coordinate_specs[name] = parsed
    end
    for name in [hinge_names; revolute_names]
        haskey(typed[name], "initial") && name ∉ rotation_coordinate_names &&
            throw(ArgumentError("$(typed[name]["type"]) '$name'.initial " *
                "requires rotation_coordinates = true"))
    end
    for name in inline_names
        haskey(typed[name], "initial") && name ∉ translation_coordinate_names &&
            throw(ArgumentError("inline constraint '$name'.initial requires " *
                "translation_coordinates = true"))
    end
    base_connection_names = sort!([spherical_names; perp_names; inplane_names;
        inline_names; hinge_names; orient_names; revolute_names; fixed_names])
    connection_names = sort!([base_connection_names; gear_pair_names;
        rack_and_pinion_names; coupler_names])
    isempty(body_names) && throw(ArgumentError(
        "spatial model requires at least one rigid body or flexible beam"))
    ground_names = Set(name for (name, table) in typed
        if table["type"] == "ground")

    registrations = Dict(name => (name in beam_names ?
        spatial_flexible_beam_registration(name) : spatial_body_registration(name))
        for name in body_names)
    for name in spherical_names
        registrations[name] = spherical_joint_registration(name)
    end
    for name in perp_names
        registrations[name] = perp_constraint_registration(name)
    end
    for name in inplane_names
        registrations[name] = inplane_constraint_registration(name)
    end
    for name in inline_names
        registrations[name] = inline_constraint_registration(name;
            translation_coordinates = name in translation_coordinate_names)
    end
    for name in hinge_names
        registrations[name] = hinge_constraint_registration(name;
            rotation_coordinates = name in rotation_coordinate_names)
    end
    for name in orient_names
        registrations[name] = orient_constraint_registration(name)
    end
    for name in revolute_names
        registrations[name] = revolute_joint_registration(name;
            rotation_coordinates = name in rotation_coordinate_names)
    end
    for name in fixed_names
        registrations[name] = fixed_joint_registration(name)
    end
    for name in gear_pair_names
        registrations[name] = spatial_gear_pair_registration(name)
    end
    for name in rack_and_pinion_names
        registrations[name] = spatial_rack_and_pinion_registration(name)
    end
    for name in coupler_names
        registrations[name] = spatial_coordinate_coupler_registration(name)
    end
    for name in spanning_force_names
        registrations[name] = spatial_spanning_force_registration(name)
    end
    for name in belt_span_names
        registrations[name] = spatial_belt_span_registration(name)
    end
    for name in applied_force_names
        registrations[name] = spatial_applied_force_registration(name)
    end
    for name in applied_torque_names
        registrations[name] = spatial_applied_torque_registration(name)
    end
    for name in bushing_names
        registrations[name] = spatial_bushing_registration(name)
    end
    for name in curve_contact_names
        registrations[name] = spatial_curve_contact_registration(name)
    end
    for name in plane_contact_names
        registrations[name] = spatial_plane_contact_registration(name)
    end
    for name in surface_friction_names
        registrations[name] = spatial_surface_friction_registration(name)
    end
    for name in revolute_friction_names
        registrations[name] = spatial_revolute_friction_registration(name)
    end
    for name in translational_friction_names
        registrations[name] = spatial_translational_friction_registration(name)
    end
    for name in inplane_friction_names
        registrations[name] = spatial_inplane_friction_registration(name)
    end
    for name in tire_names
        registrations[name] = spatial_tire_registration(name;
            transient = transient_tires[name])
    end
    for name in equation_component_names
        registrations[name] = spatial_equation_registration(name, typed[name])
    end
    for name in rotational_motion_names
        registrations[name] = spatial_rotational_motion_registration(name)
    end
    for name in translational_motion_names
        registrations[name] = spatial_translational_motion_registration(name)
    end
    for name in spanning_motion_names
        registrations[name] = spatial_spanning_motion_registration(name)
    end
    for name in span_measure_names
        registrations[name] = spatial_span_measure_registration(name)
    end
    for name in directed_distance_names
        registrations[name] = spatial_directed_distance_registration(name)
    end
    builder = ModelLayoutBuilder()
    for name in (body_names..., connection_names..., belt_span_names...,
            spanning_force_names...,
            applied_force_names..., applied_torque_names..., bushing_names...,
            curve_contact_names...,
            plane_contact_names...,
            surface_friction_names...,
            revolute_friction_names...,
            translational_friction_names...,
            inplane_friction_names...,
            tire_names..., equation_component_names...,
            rotational_motion_names..., translational_motion_names...,
            spanning_motion_names..., measure_names...)
        registration = registrations[name]
        allocate_component_variables!(builder, registration)
    end
    for name in body_names
        registration = registrations[name]
        for block in (:balance, :selected_state, :orientation)
            allocate_component_equation_block!(builder, registration, block)
        end
    end
    for name in connection_names
        registration = registrations[name]
        if name in gear_pair_names
            for block in (:phase_acceleration, :phase_velocity,
                    :phase_position, :acceleration, :velocity, :position)
                allocate_component_equation_block!(builder, registration,
                    block)
            end
        elseif name in coupler_names || name in rack_and_pinion_names
            for block in (:acceleration, :velocity, :position)
                allocate_component_equation_block!(builder, registration,
                    block)
            end
        else
            allocate_component_equation_block!(builder, registration,
                :constraint)
        end
    end
    for name in rotational_motion_names
        allocate_component_equation_block!(builder, registrations[name],
            :constraint)
    end
    for name in translational_motion_names
        for block in (:position, :velocity, :acceleration)
            allocate_component_equation_block!(builder, registrations[name],
                block)
        end
    end
    for name in spanning_motion_names
        for block in (:span_geometry, :span_velocity, :span_acceleration,
                :position, :velocity, :acceleration, :load)
            allocate_component_equation_block!(builder, registrations[name],
                block)
        end
    end
    for name in span_measure_names
        for block in (:span_geometry, :span_velocity, :span_acceleration)
            allocate_component_equation_block!(builder, registrations[name],
                block)
        end
    end
    for name in directed_distance_names
        for block in (:position, :velocity, :acceleration)
            allocate_component_equation_block!(builder, registrations[name],
                block)
        end
    end
    for name in rotation_coordinate_names
        for block in (:rotation_acceleration, :rotation_velocity,
                :rotation_position, :selected_rotation_state)
            allocate_component_equation_block!(builder, registrations[name],
                block)
        end
    end
    for name in translation_coordinate_names
        for block in (:translation_acceleration, :translation_velocity,
                :translation_position, :selected_translation_state)
            allocate_component_equation_block!(builder, registrations[name],
                block)
        end
    end
    for name in spanning_force_names
        for block in (:span_geometry, :span_velocity, :span_acceleration,
                :load)
            allocate_component_equation_block!(builder, registrations[name],
                block)
        end
    end
    for name in belt_span_names
        for block in (:geometry, :rate, :load)
            allocate_component_equation_block!(builder, registrations[name],
                block)
        end
    end
    for name in applied_force_names
        allocate_component_equation_block!(builder, registrations[name],
            :load)
    end
    for name in applied_torque_names
        allocate_component_equation_block!(builder, registrations[name],
            :load)
    end
    for name in bushing_names
        for block in (:kinematics, :load)
            allocate_component_equation_block!(builder, registrations[name],
                block)
        end
    end
    for name in curve_contact_names
        allocate_component_equation_block!(builder, registrations[name],
            :contact)
    end
    for name in plane_contact_names
        allocate_component_equation_block!(builder, registrations[name],
            :contact)
    end
    for name in surface_friction_names
        allocate_component_equation_block!(builder, registrations[name],
            :friction)
    end
    for name in revolute_friction_names
        allocate_component_equation_block!(builder, registrations[name],
            :friction)
    end
    for name in translational_friction_names
        allocate_component_equation_block!(builder, registrations[name],
            :friction)
    end
    for name in inplane_friction_names
        allocate_component_equation_block!(builder, registrations[name],
            :friction)
    end
    for name in tire_names
        blocks = transient_tires[name] ?
            (:kinematics, :load, :deformation) : (:kinematics, :load)
        for block in blocks
            allocate_component_equation_block!(builder, registrations[name],
                block)
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
    initial = zeros(Float64, length(layout.catalog.variables))
    equation_components = Dict{Symbol,SpatialEquationComponent}()
    for name in equation_component_names
        component = allocated_spatial_equation_component(layout, name,
            typed[name], parameters)
        initialize_spatial_equation_component!(initial, component, typed[name])
        equation_components[name] = component
    end
    variable_weights = ones(Float64, length(initial))
    initial_condition_weights = Dict{Symbol,NamedTuple}()
    imposed_variables = Set{Int}()
    imposed_orientations = Set{Symbol}()
    for name in body_names
        table = typed[name]
        haskey(table, "mass") || throw(ArgumentError(
            "body or beam '$name' requires mass"))
        mass = Float64(table["mass"])
        isfinite(mass) && mass > 0 || throw(ArgumentError(
            "a free spatial body requires positive mass"))
        reference = body_reference_frames[name]
        body = if name in beam_names
            isnothing(reference.center_of_mass_marker) || throw(ArgumentError(
                "flexible beam '$name' uses its floating center frame and " *
                "cannot specify center_of_mass"))
            required_positive(field) = begin
                haskey(table, field) || throw(ArgumentError(
                    "flexible beam '$name' requires $field"))
                value = Float64(table[field])
                isfinite(value) && value > 0 || throw(ArgumentError(
                    "flexible beam '$name'.$field must be positive"))
                value
            end
            length = required_positive("length")
            area = required_positive("area")
            elastic_modulus = required_positive("elastic_modulus")
            shear_modulus = required_positive("shear_modulus")
            second_moment_y = required_positive("second_moment_y")
            second_moment_z = required_positive("second_moment_z")
            torsion_constant = required_positive("torsion_constant")
            shear_y = finite_number(get(table,
                "shear_coefficient_y", 5 / 6),
                "flexible beam '$name'.shear_coefficient_y")
            shear_z = finite_number(get(table,
                "shear_coefficient_z", 5 / 6),
                "flexible beam '$name'.shear_coefficient_z")
            shear_y > 0 && shear_z > 0 || throw(ArgumentError(
                "flexible beam '$name' shear coefficients must be positive"))
            damping_time_scale = finite_number(get(table,
                "damping_time_scale", 0.0),
                "flexible beam '$name'.damping_time_scale")
            damping_time_scale >= 0 || throw(ArgumentError(
                "flexible beam '$name'.damping_time_scale must be nonnegative"))
            specified_inertia = haskey(table, "inertia") ?
                inertia_matrix(table, name) : nothing
            SpatialModeling.allocated_spatial_flexible_beam(
                layout, name, mass, length,
                area, elastic_modulus, shear_modulus, second_moment_y,
                second_moment_z, torsion_constant, shear_y, shear_z,
                damping_time_scale; inertia = specified_inertia)
        else
            inertia_at_cm = inertia_matrix(table, name)
            inertia_in_body = reference.center_of_mass_orientation *
                inertia_at_cm * transpose(reference.center_of_mass_orientation)
            allocated_spatial_body(layout, name, mass, inertia_in_body)
        end
        bodies[name] = body
        reference_position = numeric_vector(table, "position", 3;
            default = zeros(3), label = "body '$name'")
        reference_velocity = numeric_vector(table, "velocity", 3;
            default = zeros(3), label = "body '$name'")
        angular_velocity = numeric_vector(
            table, "angular_velocity", 3; default = zeros(3),
            label = "body '$name'")
        initial[body.angular_velocity_variables] .= angular_velocity
        orientation = orientation_matrix(table, "orientation", "body '$name'")
        initial[body.euler_parameter_variables] .=
            matrix_to_euler_parameters(orientation)
        if body isa SpatialFlexibleBeamComponent
            initial[body.elastic_position_variables] .= numeric_vector(table,
                "elastic_position", 6; default = zeros(6),
                label = "flexible beam '$name'")
            initial[body.elastic_velocity_variables] .= numeric_vector(table,
                "elastic_velocity", 6; default = zeros(6),
                label = "flexible beam '$name'")
        end
        initial_condition_weights[name] = body_ic_weights(name, table)

        # An imposed orientation or angular velocity participates in the
        # reference-origin to CM conversion. Reapplying the impositions after
        # conversion preserves any explicitly imposed canonical R or V values.
        apply_body_initial_impose!(initial, imposed_variables,
            imposed_orientations, body, table, name)
        orientation = rotation_matrix(
            @view initial[body.euler_parameter_variables])
        angular_velocity =
            collect(@view initial[body.angular_velocity_variables])
        center_offset = reference.center_of_mass_position
        initial[body.position_variables] .=
            reference_position + orientation * center_offset
        initial[body.velocity_variables] .= reference_velocity +
            orientation * cross(angular_velocity, center_offset)
        apply_body_initial_impose!(initial, imposed_variables,
            imposed_orientations, body, table, name)
    end

    initial_conditions = apply_saved_initial_conditions!(initial, layout,
        kinds, document, model_directory; dimension = "spatial",
        configuration_kinds = SPATIAL_CONFIGURATION_KINDS,
        velocity_kinds = SPATIAL_VELOCITY_KINDS,
        element_types = SPATIAL_SAVED_STATE_ELEMENT_TYPES)
    normalize_transferred_orientations!(initial, bodies,
        initial_conditions.transferred_variables)
    if initial_conditions.enabled
        for name in body_names
            apply_body_initial_impose!(initial, imposed_variables,
                imposed_orientations, bodies[name], typed[name], name)
        end
    end
    saved_initial_values = initial_conditions.enabled ? copy(initial) : nothing
    transferred_variables = Set(initial_conditions.transferred_variables)

    markers = Dict{Symbol,Any}()
    for name in sort!(collect(beam_names))
        beam = bodies[name]
        for node in (:end_i, :cm, :end_j)
            marker = spatial_flexible_beam_marker(beam, node)
            haskey(typed, marker.name) && throw(ArgumentError(
                "flexible beam '$name' automatically owns marker '$(marker.name)'"))
            markers[marker.name] = marker
        end
    end
    marker_names = sort!([name for (name, table) in typed
        if table["type"] == "marker"])
    for name in marker_names
        table = typed[name]
        owner = marker_owner(name, bodies, ground_names)
        position = numeric_vector(table, "position", 3;
            default = zeros(3), label = "marker '$name'")
        orientation = orientation_matrix(
            table, "orientation", "marker '$name'")
        markers[name] = if isnothing(owner)
            SpatialGroundMarker(name, position, orientation)
        else
            owner isa SpatialFlexibleBeamComponent && throw(ArgumentError(
                "flexible beam '$(owner.name)' currently provides generated " *
                "end_i, cm, and end_j markers; arbitrary beam markers are " *
                "not yet supported"))
            center_offset =
                body_reference_frames[owner.name].center_of_mass_position
            SpatialBodyMarker(
                name, owner, position - center_offset, orientation)
        end
    end

    curves = Dict{Symbol,SpatialCurveDefinition}()
    for name in curve_names
        table = typed[name]
        marker_name = get(table, "marker", nothing)
        marker_name isa AbstractString || throw(ArgumentError(
            "curve '$name'.marker must name a body or ground marker"))
        marker = required_marker(markers, marker_name, "curve '$name'")
        marker isa Union{SpatialBodyMarker,SpatialGroundMarker} ||
            throw(ArgumentError("curve '$name' marker must be fixed to a " *
                "body or ground"))
        get(table, "closed", true) == true || throw(ArgumentError(
            "curve '$name' must be closed"))
        points = get(table, "points", nothing)
        points isa Vector || throw(ArgumentError(
            "curve '$name'.points must be an array of [x, y] points"))
        profile = PlanarClosedCurve(points)
        half_width = finite_number(get(table, "half_width", 0.04),
            "curve '$name'.half_width")
        half_width > 0 || throw(ArgumentError(
            "curve '$name'.half_width must be positive"))
        curves[name] = SpatialCurveDefinition(name, profile, marker, half_width)
    end

    body_lengths = spatial_body_lengths(bodies, markers)
    positive_masses = [Float64(body.mass) for body in values(bodies)
        if body.mass > 0]
    reference_mass = isempty(positive_masses) ? 1.0 : minimum(positive_masses)
    mass_floor = 1.0e-6 * reference_mass
    for name in body_names
        body = bodies[name]
        scale = initial_condition_weights[name].scale
        effective_mass = max(Float64(body.mass), mass_floor)
        weight = scale * effective_mass
        initial_condition_weights[name] =
            (; scale, effective_mass, weight)
        variable_weights[body.position_variables] .= weight
        variable_weights[body.velocity_variables] .= weight
        variable_weights[body.angular_velocity_variables] .=
            weight * body_lengths[name]^2
        if body isa SpatialFlexibleBeamComponent
            variable_weights[body.elastic_position_variables[1:3]] .= weight
            variable_weights[body.elastic_velocity_variables[1:3]] .= weight
            variable_weights[body.elastic_position_variables[4:6]] .=
                weight * body_lengths[name]^2
            variable_weights[body.elastic_velocity_variables[4:6]] .=
                weight * body_lengths[name]^2
        end
    end

    measures = Dict{Symbol,Any}()
    for name in span_measure_names
        table = typed[name]
        endpoints = get(table, "markers", nothing)
        endpoints isa Vector && length(endpoints) == 2 || throw(ArgumentError(
            "span '$name'.markers must contain two marker names"))
        marker_1 = required_marker(markers, endpoints[1], "span '$name'")
        marker_2 = required_marker(markers, endpoints[2], "span '$name'")
        initial_separation = spatial_marker_position(marker_2, initial) -
            spatial_marker_position(marker_1, initial)
        norm(initial_separation) > 0 || throw(ArgumentError(
            "span '$name' requires initially separated markers"))
        measure = allocated_spatial_span_measure(layout, name,
            marker_1, marker_2)
        initialize_spatial_span_measure!(initial, measure)
        measures[name] = measure
    end
    for name in directed_distance_names
        table = typed[name]
        endpoints = get(table, "markers", nothing)
        endpoints isa Vector && length(endpoints) == 2 || throw(ArgumentError(
            "directed distance '$name'.markers must contain two marker names"))
        marker_i = required_marker(markers, endpoints[1],
            "directed distance '$name'")
        marker_j = required_marker(markers, endpoints[2],
            "directed distance '$name'")
        measure = allocated_spatial_directed_distance_measure(layout, name,
            marker_i, marker_j)
        initialize_spatial_directed_distance_measure!(initial, measure)
        measures[name] = measure
    end

    connections = Dict{Symbol,Any}()
    for name in spherical_names
        table = typed[name]
        endpoints = get(table, "markers", nothing)
        endpoints isa Vector && length(endpoints) == 2 || throw(ArgumentError(
            "spherical joint '$name'.markers must contain two marker names"))
        marker_a = required_marker(markers, endpoints[1],
            "spherical joint '$name'")
        marker_b = required_marker(markers, endpoints[2],
            "spherical joint '$name'")
        marker_a isa SpatialGroundMarker &&
            marker_b isa SpatialGroundMarker && throw(ArgumentError(
                "spherical joint '$name' cannot connect ground to ground"))
        connections[name] = allocated_spherical_joint(layout, name,
            marker_a, marker_b)
    end
    for name in perp_names
        table = typed[name]
        endpoints = get(table, "markers", nothing)
        endpoints isa Vector && length(endpoints) == 2 || throw(ArgumentError(
            "perp constraint '$name'.markers must contain two marker names"))
        marker_i = required_marker(markers, endpoints[1],
            "perp constraint '$name'")
        marker_j = required_marker(markers, endpoints[2],
            "perp constraint '$name'")
        marker_i isa SpatialGroundMarker &&
            marker_j isa SpatialGroundMarker && throw(ArgumentError(
                "perp constraint '$name' cannot connect ground to ground"))
        connections[name] = allocated_perp_constraint(layout, name,
            marker_i, marker_j)
    end
    for name in inplane_names
        table = typed[name]
        endpoints = get(table, "markers", nothing)
        endpoints isa Vector && length(endpoints) == 2 || throw(ArgumentError(
            "inplane constraint '$name'.markers must contain two marker names"))
        marker_i = required_marker(markers, endpoints[1],
            "inplane constraint '$name'")
        marker_j = required_marker(markers, endpoints[2],
            "inplane constraint '$name'")
        marker_i isa SpatialGroundMarker &&
            marker_j isa SpatialGroundMarker && throw(ArgumentError(
                "inplane constraint '$name' cannot connect ground to ground"))
        connections[name] = allocated_inplane_constraint(layout, name,
            marker_i, marker_j)
    end
    for name in inline_names
        table = typed[name]
        endpoints = get(table, "markers", nothing)
        endpoints isa Vector && length(endpoints) == 2 || throw(ArgumentError(
            "inline constraint '$name'.markers must contain two marker names"))
        marker_i = required_marker(markers, endpoints[1],
            "inline constraint '$name'")
        marker_j = required_marker(markers, endpoints[2],
            "inline constraint '$name'")
        marker_i isa SpatialGroundMarker &&
            marker_j isa SpatialGroundMarker && throw(ArgumentError(
                "inline constraint '$name' cannot connect ground to ground"))
        connections[name] = allocated_inline_constraint(layout, name,
            marker_i, marker_j;
            translation_coordinates = name in translation_coordinate_names)
    end
    for name in hinge_names
        table = typed[name]
        endpoints = get(table, "markers", nothing)
        endpoints isa Vector && length(endpoints) == 2 || throw(ArgumentError(
            "hinge constraint '$name'.markers must contain two marker names"))
        marker_i = required_marker(markers, endpoints[1],
            "hinge constraint '$name'")
        marker_j = required_marker(markers, endpoints[2],
            "hinge constraint '$name'")
        marker_i isa SpatialGroundMarker &&
            marker_j isa SpatialGroundMarker && throw(ArgumentError(
                "hinge constraint '$name' cannot connect ground to ground"))
        connections[name] = allocated_hinge_constraint(layout, name,
            marker_i, marker_j;
            rotation_coordinates = name in rotation_coordinate_names)
    end
    for name in orient_names
        table = typed[name]
        endpoints = get(table, "markers", nothing)
        endpoints isa Vector && length(endpoints) == 2 || throw(ArgumentError(
            "orient constraint '$name'.markers must contain two marker names"))
        marker_i = required_marker(markers, endpoints[1],
            "orient constraint '$name'")
        marker_j = required_marker(markers, endpoints[2],
            "orient constraint '$name'")
        marker_i isa SpatialGroundMarker &&
            marker_j isa SpatialGroundMarker && throw(ArgumentError(
                "orient constraint '$name' cannot connect ground to ground"))
        connections[name] = allocated_orient_constraint(layout, name,
            marker_i, marker_j)
    end
    for name in revolute_names
        table = typed[name]
        endpoints = get(table, "markers", nothing)
        endpoints isa Vector && length(endpoints) == 2 || throw(ArgumentError(
            "revolute joint '$name'.markers must contain two marker names"))
        marker_a = required_marker(markers, endpoints[1],
            "revolute joint '$name'")
        marker_b = required_marker(markers, endpoints[2],
            "revolute joint '$name'")
        marker_a isa SpatialGroundMarker &&
            marker_b isa SpatialGroundMarker && throw(ArgumentError(
                "revolute joint '$name' cannot connect ground to ground"))
        connections[name] = allocated_revolute_joint(layout, name,
            marker_a, marker_b;
            rotation_coordinates = name in rotation_coordinate_names)
    end
    for name in fixed_names
        table = typed[name]
        endpoints = get(table, "markers", nothing)
        endpoints isa Vector && length(endpoints) == 2 || throw(ArgumentError(
            "fixed joint '$name'.markers must contain two marker names"))
        marker_a = required_marker(markers, endpoints[1],
            "fixed joint '$name'")
        marker_b = required_marker(markers, endpoints[2],
            "fixed joint '$name'")
        marker_a isa SpatialGroundMarker &&
            marker_b isa SpatialGroundMarker && throw(ArgumentError(
                "fixed joint '$name' cannot connect ground to ground"))
        connections[name] = allocated_fixed_joint(layout, name,
            marker_a, marker_b)
    end

    # Pulleys refer to completed revolute joints. The pulley body is required
    # to be the first side so its continuous joint coordinate has the expected
    # sign in the no-slip belt relation.
    for name in pulley_names
        table = typed[name]
        joint = connections[pulley_joint_names[name]]
        body_specification = get(table, "body", nothing)
        body_specification isa AbstractString || throw(ArgumentError(
            "pulley '$name'.body must name its rigid body"))
        body_name = Symbol(body_specification)
        haskey(bodies, body_name) || throw(ArgumentError(
            "pulley '$name' names unknown body '$body_specification'"))
        body = bodies[body_name]
        joint.marker_a isa SpatialBodyMarker &&
            joint.marker_a.body === body || throw(ArgumentError(
            "pulley '$name' body must be the first marker body of revolute " *
            "joint '$(joint.name)'"))
        radius = finite_number(get(table, "pitch_radius", NaN),
            "pulley '$name'.pitch_radius")
        radius > 0 || throw(ArgumentError(
            "pulley '$name'.pitch_radius must be positive"))
        connections[name] = SpatialPulleyComponent(name, body, joint,
            joint.marker_a, radius)
    end

    # Gear pairs refer to completed revolute joints. Their generated floating
    # markers belong to the gears but follow the carrier contact point.
    for name in gear_pair_names
        table = typed[name]
        joint_specs = get(table, "joints", nothing)
        joint_specs isa Vector && length(joint_specs) == 2 &&
            all(value -> value isa AbstractString, joint_specs) ||
            throw(ArgumentError(
                "gear pair '$name'.joints must contain two revolute joint names"))
        length(unique(joint_specs)) == 2 || throw(ArgumentError(
            "gear pair '$name'.joints must name two different revolute joints"))
        joint_names = Symbol.(joint_specs)
        all(joint -> haskey(connections, joint) &&
            connections[joint] isa SpatialRevoluteJoint, joint_names) ||
            throw(ArgumentError(
                "gear pair '$name'.joints must reference spatial revolute joints"))
        contact_spec = get(table, "contact_marker", nothing)
        contact_marker = required_marker(markers, contact_spec,
            "gear pair '$name'.contact_marker")
        contact_marker isa SpatialFloatingMarker && throw(ArgumentError(
            "gear pair '$name'.contact_marker must be fixed to its carrier"))
        joints = (connections[joint_names[1]], connections[joint_names[2]])
        generated_names = (Symbol(name, ".contact_1"),
            Symbol(name, ".contact_2"))
        any(generated -> haskey(markers, generated), generated_names) &&
            throw(ArgumentError(
                "gear pair '$name' generated contact marker name is already in use"))
        floating = (SpatialFloatingMarker(generated_names[1],
                joints[1].marker_a.body, contact_marker),
            SpatialFloatingMarker(generated_names[2],
                joints[2].marker_a.body, contact_marker))
        markers[generated_names[1]] = floating[1]
        markers[generated_names[2]] = floating[2]
        connections[name] = allocated_spatial_gear_pair(layout, name,
            joints[1], joints[2], contact_marker, floating[1], floating[2],
            initial; phase = get(table, "phase", "initial"))
    end

    # A spatial spur rack and pinion uses the inline base frame as its pitch
    # frame: x is the pinion axis, y points toward contact, and z is rack
    # travel. The carrier contact and two body-owned floating markers are
    # generated from the referenced joints and pitch radius.
    for name in rack_and_pinion_names
        table = typed[name]
        inline_name, revolute_name = rack_and_pinion_joint_names[name]
        inline = connections[inline_name]
        revolute = connections[revolute_name]
        pitch_radius = finite_number(get(table, "pitch_radius", NaN),
            "rack and pinion '$name'.pitch_radius")
        pitch_radius > 0 || throw(ArgumentError(
            "rack and pinion '$name'.pitch_radius must be positive"))
        alignment_tolerance = finite_number(get(table,
            "alignment_tolerance", 1.0e-5),
            "rack and pinion '$name'.alignment_tolerance")
        alignment_tolerance > 0 || throw(ArgumentError(
            "rack and pinion '$name'.alignment_tolerance must be positive"))
        geometry = spatial_rack_and_pinion_geometry(inline, revolute,
            pitch_radius, initial; alignment_tolerance)

        generated_names = (Symbol(name, ".carrier_contact"),
            Symbol(name, ".rack_contact"),
            Symbol(name, ".pinion_contact"))
        any(generated -> haskey(markers, generated), generated_names) &&
            throw(ArgumentError("rack and pinion '$name' generated " *
                "contact marker name is already in use"))
        carrier_contact = if isnothing(geometry.carrier)
            SpatialGroundMarker(generated_names[1],
                collect(geometry.contact_position),
                Matrix(geometry.guide_orientation))
        else
            carrier = geometry.carrier
            carrier_orientation = rotation_matrix(
                @view initial[carrier.euler_parameter_variables])
            position_body = transpose(carrier_orientation) *
                (geometry.contact_position -
                 initial[carrier.position_variables])
            orientation_body = transpose(carrier_orientation) *
                geometry.guide_orientation
            SpatialBodyMarker(generated_names[1], carrier,
                collect(position_body), Matrix(orientation_body))
        end
        rack_contact = SpatialFloatingMarker(generated_names[2],
            geometry.rack_body, carrier_contact)
        pinion_contact = SpatialFloatingMarker(generated_names[3],
            geometry.pinion_body, carrier_contact)
        markers[generated_names[1]] = carrier_contact
        markers[generated_names[2]] = rack_contact
        markers[generated_names[3]] = pinion_contact
        connections[name] = allocated_spatial_rack_and_pinion(layout, name,
            inline, revolute, carrier_contact, rack_contact, pinion_contact,
            pitch_radius, initial; alignment_tolerance,
            phase = get(table, "phase", "initial"))
    end

    relative_coordinate_indices = Int[]
    relative_velocity_indices = Int[]
    relative_position_rows = Int[]
    relative_velocity_rows = Int[]
    relative_initial_specifications = Dict{Symbol,NamedTuple}()
    for name in gear_pair_names
        gear = connections[name]
        append!(relative_coordinate_indices, gear.position_variables)
        append!(relative_velocity_indices, gear.velocity_variables)
        append!(relative_position_rows, gear.position_definition_equations)
        append!(relative_velocity_rows, gear.velocity_definition_equations)
    end
    for name in base_connection_names
        connection = connections[name]
        hinge = connection_hinge(connection)
        if !isnothing(hinge) && !isempty(hinge.rotation_variables)
            initialize_hinge_coordinates!(initial, hinge, 1)
            _, velocity, coordinate = hinge.rotation_variables
            specification = relative_initial_specification(typed[name],
                "$(typed[name]["type"]) '$name'", :rotation)
            relative_initial_specifications[name] = specification
            apply_relative_initial_specification!(initial, variable_weights,
                imposed_variables, coordinate, velocity, specification)
            coordinate_name = Symbol(name, :.,
                layout.catalog.variables[coordinate].name)
            velocity_name = Symbol(name, :.,
                layout.catalog.variables[velocity].name)
            if specification.coordinate_specified ||
                    coordinate_name in transferred_variables
                push!(relative_coordinate_indices, coordinate)
                push!(relative_position_rows, hinge.rotation_equations[3])
            end
            if specification.velocity_specified ||
                    velocity_name in transferred_variables
                push!(relative_velocity_indices, velocity)
                push!(relative_velocity_rows, hinge.rotation_equations[2])
            end
        end
        inline = connection_inline(connection)
        if !isnothing(inline) && !isempty(inline.translation_variables)
            initialize_inline_coordinates!(initial, inline, 1)
            _, velocity, coordinate = inline.translation_variables
            specification = relative_initial_specification(typed[name],
                "inline constraint '$name'", :translation)
            relative_initial_specifications[name] = specification
            apply_relative_initial_specification!(initial, variable_weights,
                imposed_variables, coordinate, velocity, specification)
            coordinate_name = Symbol(name, :.,
                layout.catalog.variables[coordinate].name)
            velocity_name = Symbol(name, :.,
                layout.catalog.variables[velocity].name)
            if specification.coordinate_specified ||
                    coordinate_name in transferred_variables
                push!(relative_coordinate_indices, coordinate)
                push!(relative_position_rows, inline.translation_equations[3])
            end
            if specification.velocity_specified ||
                    velocity_name in transferred_variables
                push!(relative_velocity_indices, velocity)
                push!(relative_velocity_rows, inline.translation_equations[2])
            end
        end
    end

    belt_forces = Dict{Symbol,SpatialBeltSpanComponent}()
    claimed_belt_spans = Symbol[]
    for belt_name in belt_names
        table = typed[belt_name]
        ordered_specs = get(table, "spans", nothing)
        ordered_specs isa Vector && length(ordered_specs) >= 2 &&
            all(item -> item isa AbstractString, ordered_specs) ||
            throw(ArgumentError(
                "belt '$belt_name'.spans must contain at least two names"))
        ordered_names = Symbol.(ordered_specs)
        length(unique(ordered_names)) == length(ordered_names) ||
            throw(ArgumentError("belt '$belt_name' spans must be unique"))
        specifications = NamedTuple[]
        for span_name in ordered_names
            span_name in belt_span_names || throw(ArgumentError(
                "belt '$belt_name' references unknown belt span '$span_name'"))
            push!(claimed_belt_spans, span_name)
            span_table = typed[span_name]
            pulley_specs = get(span_table, "pulleys", nothing)
            pulley_specs isa Vector && length(pulley_specs) == 2 &&
                all(item -> item isa AbstractString, pulley_specs) ||
                throw(ArgumentError(
                    "belt span '$span_name'.pulleys must contain two names"))
            pulley_names_for_span = Symbol.(pulley_specs)
            all(pulley -> haskey(connections, pulley) &&
                    connections[pulley] isa SpatialPulleyComponent,
                pulley_names_for_span) || throw(ArgumentError(
                "belt span '$span_name' endpoints must be spatial pulleys"))
            pulley_1 = connections[pulley_names_for_span[1]]
            pulley_2 = connections[pulley_names_for_span[2]]
            pulley_1 !== pulley_2 || throw(ArgumentError(
                "belt span '$span_name' requires two different pulleys"))
            hint_specs = get(span_table, "near_points", nothing)
            hint_specs isa Vector && length(hint_specs) == 2 &&
                all(item -> item isa AbstractString, hint_specs) ||
                throw(ArgumentError(
                    "belt span '$span_name'.near_points must contain two markers"))
            hint_1 = required_marker(markers, hint_specs[1],
                "belt span '$span_name'.near_points")
            hint_2 = required_marker(markers, hint_specs[2],
                "belt span '$span_name'.near_points")
            tolerance = finite_number(get(span_table,
                "feasibility_tolerance", get(table,
                    "feasibility_tolerance", 1.0e-5)),
                "belt span '$span_name'.feasibility_tolerance")
            tolerance > 0 || throw(ArgumentError(
                "belt span '$span_name'.feasibility_tolerance must be positive"))
            selected = select_spatial_belt_tangent(pulley_1, pulley_2,
                hint_1, hint_2, initial; tolerance)
            stiffness = finite_number(get(span_table, "stiffness",
                get(table, "stiffness", 0.0)),
                "belt span '$span_name'.stiffness")
            stiffness > 0 || throw(ArgumentError(
                "belt span '$span_name'.stiffness must be positive"))
            damping_time_scale = finite_number(get(span_table,
                "damping_time_scale", get(table,
                    "damping_time_scale", 0.0)),
                "belt span '$span_name'.damping_time_scale")
            damping_time_scale >= 0 || throw(ArgumentError(
                "belt span '$span_name'.damping_time_scale must be nonnegative"))
            damping = finite_number(get(span_table, "damping",
                get(table, "damping", damping_time_scale * stiffness)),
                "belt span '$span_name'.damping")
            damping >= 0 || throw(ArgumentError(
                "belt span '$span_name'.damping must be nonnegative"))
            push!(specifications, (; name = span_name, pulley_1, pulley_2,
                geometry = selected.geometry, selected, stiffness, damping,
                damping_time_scale, tolerance))
        end
        path_length = spatial_belt_reference_path_length(specifications)
        has_initial_tension = haskey(table, "initial_tension")
        has_free_length = haskey(table, "free_length")
        has_initial_tension && has_free_length && throw(ArgumentError(
            "belt '$belt_name' cannot specify both initial_tension and free_length"))
        initial_tension = if has_initial_tension
            value = finite_number(table["initial_tension"],
                "belt '$belt_name'.initial_tension")
            value >= 0 || throw(ArgumentError(
                "belt '$belt_name'.initial_tension must be nonnegative"))
            value
        elseif has_free_length
            free = finite_number(table["free_length"],
                "belt '$belt_name'.free_length")
            free > 0 || throw(ArgumentError(
                "belt '$belt_name'.free_length must be positive"))
            (path_length - free) /
                sum(inv(spec.stiffness) for spec in specifications)
        else
            0.0
        end
        initial_tension < 0 && @warn(
            "belt '$belt_name' begins with negative tension",
            initial_tension = initial_tension)
        free_length = has_free_length ? Float64(table["free_length"]) :
            path_length - sum(initial_tension / spec.stiffness
                for spec in specifications)
        spans = SpatialBeltSpanComponent[]
        for spec in specifications
            initial_extension = initial_tension / spec.stiffness
            span = allocated_spatial_belt_span(layout, spec.name, belt_name,
                spec.pulley_1, spec.pulley_2, spec.selected,
                initial, initial_extension,
                spec.stiffness, spec.damping, spec.damping_time_scale;
                feasibility_tolerance = spec.tolerance)
            initialize_spatial_belt_span!(initial, span)
            belt_forces[spec.name] = span
            push!(spans, span)
        end
        connections[belt_name] = SpatialBeltComponent(belt_name, spans,
            initial_tension, free_length)
    end
    sort!(claimed_belt_spans) == sort!(belt_span_names) || throw(ArgumentError(
        "every spatial belt_span must belong to exactly one belt"))

    drivers = Dict{Symbol,Any}()
    driven_connection_names = Set(values(driver_joint_names))
    for name in rotational_motion_names
        joint_name = driver_joint_names[name]
        hinge = connection_hinge(connections[joint_name])
        motion, velocity, acceleration = spatial_rotational_driver_laws(
            typed[name], parameters, name)
        driver = allocated_spatial_rotational_motion_generator(layout, name,
            hinge, motion, velocity, acceleration)
        alpha, omega, theta = hinge.rotation_variables
        initial[theta] = motion(simulation.start_time)
        initial[omega] = velocity(simulation.start_time)
        initial[alpha] = acceleration(simulation.start_time)
        all(isfinite, initial[[alpha, omega, theta]]) || throw(ArgumentError(
            "rotational motion '$name' is not finite initially"))
        push!(relative_coordinate_indices, theta)
        push!(relative_position_rows, hinge.rotation_equations[3])
        push!(relative_velocity_indices, omega)
        push!(relative_velocity_rows, hinge.rotation_equations[2])
        push!(imposed_variables, theta)
        push!(imposed_variables, omega)
        drivers[name] = driver
    end
    for name in translational_motion_names
        table = typed[name]
        endpoints = get(table, "markers", nothing)
        endpoints isa Vector && length(endpoints) == 2 || throw(ArgumentError(
            "translational motion '$name'.markers must contain two marker names"))
        marker_i = required_marker(markers, endpoints[1],
            "translational motion '$name'")
        marker_j = required_marker(markers, endpoints[2],
            "translational motion '$name'")
        marker_i isa SpatialGroundMarker &&
            marker_j isa SpatialGroundMarker && throw(ArgumentError(
                "translational motion '$name' cannot connect ground to ground"))
        motion, velocity, acceleration = spatial_translational_driver_laws(
            table, parameters, name)
        driver = allocated_spatial_translational_motion_generator(layout,
            name, marker_i, marker_j, motion, velocity, acceleration)
        initial[driver.distance_variable] = motion(simulation.start_time)
        initial[driver.velocity_variable] = velocity(simulation.start_time)
        initial[driver.acceleration_variable] =
            acceleration(simulation.start_time)
        indices = [driver.distance_variable, driver.velocity_variable,
            driver.acceleration_variable]
        all(isfinite, initial[indices]) || throw(ArgumentError(
            "translational motion '$name' is not finite initially"))
        push!(relative_coordinate_indices, driver.distance_variable)
        push!(relative_position_rows, first(driver.position_equations))
        push!(relative_velocity_indices, driver.velocity_variable)
        push!(relative_velocity_rows, first(driver.velocity_equations))
        push!(imposed_variables, driver.distance_variable)
        push!(imposed_variables, driver.velocity_variable)
        drivers[name] = driver
    end
    for name in spanning_motion_names
        table = typed[name]
        endpoints = get(table, "markers", nothing)
        endpoints isa Vector && length(endpoints) == 2 || throw(ArgumentError(
            "spanning motion '$name'.markers must contain two marker names"))
        marker_1 = required_marker(markers, endpoints[1],
            "spanning motion '$name'")
        marker_2 = required_marker(markers, endpoints[2],
            "spanning motion '$name'")
        marker_1 isa SpatialGroundMarker &&
            marker_2 isa SpatialGroundMarker && throw(ArgumentError(
                "spanning motion '$name' cannot connect ground to ground"))
        initial_separation = spatial_marker_position(marker_2, initial) -
            spatial_marker_position(marker_1, initial)
        norm(initial_separation) > 0 || throw(ArgumentError(
            "spanning motion '$name' requires initially separated markers " *
            "to establish its line of action"))
        motion, velocity, acceleration = spatial_spanning_driver_laws(
            table, parameters, name)
        driver = allocated_spatial_spanning_motion_generator(layout, name,
            marker_1, marker_2, motion, velocity, acceleration)
        initialize_spatial_span!(initial, driver.span)
        initial[driver.span.distance_variable] = motion(simulation.start_time)
        initial[driver.span.velocity_variable] = velocity(simulation.start_time)
        initial[driver.span.acceleration_variable] =
            acceleration(simulation.start_time)
        indices = [driver.span.distance_variable, driver.span.velocity_variable,
            driver.span.acceleration_variable]
        all(isfinite, initial[indices]) &&
            initial[driver.span.distance_variable] > 0 || throw(ArgumentError(
                "spanning motion '$name' must be positive and finite initially"))
        push!(relative_coordinate_indices, driver.span.distance_variable)
        push!(relative_position_rows, driver.span.distance_equation)
        push!(relative_velocity_indices, driver.span.velocity_variable)
        push!(relative_velocity_rows, driver.span.velocity_equation)
        push!(imposed_variables, driver.span.distance_variable)
        push!(imposed_variables, driver.span.velocity_variable)
        drivers[name] = driver
    end

    # Construct couplers after initializing any motion generators so an
    # `initial` offset records the actual starting values of driven ports.
    for name in coupler_names
        table = typed[name]
        coordinates = Any[]
        coordinate_kinds = Symbol[]
        for (component_name, port) in coupler_coordinate_specs[name]
            connection = connections[component_name]
            if port == :rotation
                coordinate = connection_hinge(connection)
                isnothing(coordinate) && error(
                    "coupler '$name' lost rotation coordinate '$component_name'")
                push!(coordinate_kinds, :rotation)
            else
                coordinate = connection_inline(connection)
                isnothing(coordinate) && error(
                    "coupler '$name' lost distance coordinate '$component_name'")
                push!(coordinate_kinds, :distance)
            end
            push!(coordinates, coordinate)
        end
        raw_coefficients = get(table, "coefficients", Any[])
        length(raw_coefficients) == length(coordinates) || throw(ArgumentError(
            "coupler '$name' requires one coefficient per coordinate"))
        all(value -> value isa Number, raw_coefficients) || throw(ArgumentError(
            "coupler '$name' coefficients must be numeric"))
        coefficients = Float64.(raw_coefficients)
        all(isfinite, coefficients) || throw(ArgumentError(
            "coupler '$name' coefficients must be finite"))
        count(value -> !iszero(value), coefficients) >= 2 || throw(ArgumentError(
            "coupler '$name' requires at least two nonzero coefficients"))
        initial_sum = sum(factor * initial[
            spatial_coordinate_variables(coordinate)[1]]
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
        coordinate_kind = length(unique(coordinate_kinds)) == 1 ?
            first(coordinate_kinds) : :mixed
        connections[name] = allocated_spatial_coordinate_coupler(layout,
            name, coordinates, coefficients, offset, coordinate_kind)
    end

    forces = Dict{Symbol,Any}(belt_forces)
    for name in curve_contact_names
        table = typed[name]
        curve_name = get(table, "curve", nothing)
        curve_name isa AbstractString || throw(ArgumentError(
            "curve contact '$name'.curve must name a curve"))
        curve_definition = get(curves, Symbol(curve_name), nothing)
        isnothing(curve_definition) && throw(ArgumentError(
            "curve contact '$name' references unknown curve '$curve_name'"))
        follower_kind = table["type"] == "curve_contact" ? :roller : :flat
        marker_field = follower_kind == :roller ? "roller_marker" :
            "follower_marker"
        marker_name = get(table, marker_field, nothing)
        marker_name isa AbstractString || throw(ArgumentError(
            "curve contact '$name'.$marker_field must name a marker"))
        follower_marker = required_marker(markers, marker_name,
            "curve contact '$name'")
        follower_marker isa Union{SpatialBodyMarker,SpatialGroundMarker} ||
            throw(ArgumentError("curve contact '$name' follower marker must " *
                "be fixed to a body or ground"))
        curve_marker = curve_definition.marker
        curve_marker isa SpatialGroundMarker &&
            follower_marker isa SpatialGroundMarker && throw(ArgumentError(
                "curve contact '$name' cannot connect ground to ground"))
        curve_marker isa SpatialBodyMarker &&
            follower_marker isa SpatialBodyMarker &&
            curve_marker.body === follower_marker.body && throw(ArgumentError(
                "curve contact '$name' markers must belong to different bodies"))
        radius = follower_kind == :roller ? finite_number(
            get(table, "radius", NaN), "curve contact '$name'.radius") : 0.0
        follower_kind == :roller && radius <= 0 && throw(ArgumentError(
            "curve contact '$name'.radius must be positive"))
        side_name = lowercase(string(get(table, "side", "outside")))
        side_name in ("outside", "inside") || throw(ArgumentError(
            "curve contact '$name'.side must be 'outside' or 'inside'"))
        has_stiffness = haskey(table, "stiffness")
        has_expression = haskey(table, "expression")
        xor(has_stiffness, has_expression) || throw(ArgumentError(
            "curve contact '$name' requires exactly one of stiffness or expression"))
        predefined_fields = ("damping_factor", "transition_depth")
        has_expression && any(haskey(table, field)
            for field in predefined_fields) && throw(ArgumentError(
            "curve contact '$name' damping_factor and transition_depth " *
            "cannot be combined with expression"))
        stiffness = has_stiffness ? finite_number(table["stiffness"],
            "curve contact '$name'.stiffness") : 0.0
        has_stiffness && stiffness <= 0 && throw(ArgumentError(
            "curve contact '$name'.stiffness must be positive"))
        damping_factor = finite_number(get(table, "damping_factor", 0.0),
            "curve contact '$name'.damping_factor")
        damping_factor >= 0 || throw(ArgumentError(
            "curve contact '$name'.damping_factor must be nonnegative"))
        transition_depth = finite_number(get(table, "transition_depth", 0.0),
            "curve contact '$name'.transition_depth")
        transition_depth >= 0 || throw(ArgumentError(
            "curve contact '$name'.transition_depth must be nonnegative"))
        axis_tolerance = finite_number(get(table, "axis_tolerance", 1.0e-6),
            "curve contact '$name'.axis_tolerance")
        axis_tolerance > 0 || throw(ArgumentError(
            "curve contact '$name'.axis_tolerance must be positive"))
        law = has_expression ? compile_spatial_model_expression(
            string(table["expression"]), parameters, layout) : nothing
        active_during = active_during_stages(table, "curve contact '$name'")
        component = allocated_spatial_curve_contact(layout, name,
            curve_definition.profile, curve_marker, follower_marker,
            radius, curve_definition.half_width, stiffness, damping_factor; law,
            expression = has_expression, transition_depth,
            side = Symbol(side_name), follower_kind, axis_tolerance,
            active_during)
        station = haskey(table, "initial_station") ? finite_number(
            table["initial_station"],
            "curve contact '$name'.initial_station") : nothing
        initialize_spatial_curve_contact!(initial, component,
            simulation.start_time; station)
        forces[name] = component
    end
    for name in plane_contact_names
        table = typed[name]
        endpoints = get(table, "markers", nothing)
        endpoints isa Vector && length(endpoints) == 2 || throw(ArgumentError(
            "plane contact '$name'.markers must contain a sphere-center " *
            "marker and a plane marker"))
        sphere_marker = required_marker(markers, endpoints[1],
            "plane contact '$name'")
        plane_marker = required_marker(markers, endpoints[2],
            "plane contact '$name'")
        sphere_marker isa Union{SpatialBodyMarker,SpatialGroundMarker} ||
            throw(ArgumentError("plane contact '$name' sphere marker must " *
                "be fixed to a body or ground"))
        plane_marker isa Union{SpatialBodyMarker,SpatialGroundMarker} ||
            throw(ArgumentError("plane contact '$name' plane marker must " *
                "be fixed to a body or ground"))
        sphere_marker isa SpatialGroundMarker &&
            plane_marker isa SpatialGroundMarker && throw(ArgumentError(
                "plane contact '$name' cannot connect ground to ground"))
        sphere_marker isa SpatialBodyMarker &&
            plane_marker isa SpatialBodyMarker &&
            sphere_marker.body === plane_marker.body && throw(ArgumentError(
                "plane contact '$name' markers must belong to different bodies"))
        radius = finite_number(get(table, "radius", NaN),
            "plane contact '$name'.radius")
        radius > 0 || throw(ArgumentError(
            "plane contact '$name'.radius must be positive"))
        has_stiffness = haskey(table, "stiffness")
        has_expression = haskey(table, "expression")
        xor(has_stiffness, has_expression) || throw(ArgumentError(
            "plane contact '$name' requires exactly one of stiffness or expression"))
        predefined_fields = ("damping_factor", "transition_depth")
        has_expression && any(haskey(table, field)
            for field in predefined_fields) && throw(ArgumentError(
            "plane contact '$name' damping_factor and transition_depth " *
            "cannot be combined with expression"))
        stiffness = has_stiffness ? finite_number(table["stiffness"],
            "plane contact '$name'.stiffness") : 0.0
        has_stiffness && stiffness <= 0 && throw(ArgumentError(
            "plane contact '$name'.stiffness must be positive"))
        damping_factor = finite_number(get(table, "damping_factor", 0.0),
            "plane contact '$name'.damping_factor")
        damping_factor >= 0 || throw(ArgumentError(
            "plane contact '$name'.damping_factor must be nonnegative"))
        transition_depth = finite_number(get(table, "transition_depth", 0.0),
            "plane contact '$name'.transition_depth")
        transition_depth >= 0 || throw(ArgumentError(
            "plane contact '$name'.transition_depth must be nonnegative"))
        law = has_expression ? compile_spatial_model_expression(
            string(table["expression"]), parameters, layout) : nothing
        active_during = active_during_stages(table, "plane contact '$name'")
        component = allocated_spatial_plane_contact(layout, name,
            sphere_marker, plane_marker, radius, stiffness, damping_factor;
            law, expression = has_expression, transition_depth, active_during)
        initialize_spatial_plane_contact!(initial, component,
            simulation.start_time)
        forces[name] = component
    end
    for name in surface_friction_names
        table = typed[name]
        contact_name = get(table, "contact", nothing)
        contact_name isa AbstractString || throw(ArgumentError(
            "surface friction '$name'.contact must name a plane_contact"))
        contact = get(forces, Symbol(contact_name), nothing)
        contact isa SpatialPlaneContactComponent || throw(ArgumentError(
            "surface friction '$name'.contact must name an existing plane_contact"))
        stiffness = finite_number(get(table, "stiffness", NaN),
            "surface friction '$name'.stiffness")
        damping = finite_number(get(table, "damping", 0.0),
            "surface friction '$name'.damping")
        mu_static = finite_number(get(table, "static_coefficient", NaN),
            "surface friction '$name'.static_coefficient")
        mu_dynamic = finite_number(get(table, "dynamic_coefficient", NaN),
            "surface friction '$name'.dynamic_coefficient")
        speed = finite_number(get(table, "transition_speed", 0.01),
            "surface friction '$name'.transition_speed")
        release = finite_number(get(table, "release_time", 0.01),
            "surface friction '$name'.release_time")
        stiffness > 0 && damping >= 0 &&
            mu_static >= mu_dynamic >= 0 && speed > 0 &&
            release > 0 || throw(ArgumentError(
            "surface friction '$name' requires positive stiffness, " *
            "static_coefficient >= dynamic_coefficient >= 0, " *
            "positive transition_speed and release_time"))
        friction = allocated_spatial_surface_friction(layout, name, contact,
            stiffness, damping, mu_static, mu_dynamic, speed, release)
        initialize_spatial_surface_friction!(initial, friction;
            reset_anchor = true)
        forces[name] = friction
    end
    for name in revolute_friction_names
        table = typed[name]
        joint_name = get(table, "joint", nothing)
        joint_name isa AbstractString || throw(ArgumentError(
            "revolute friction '$name'.joint must name a revolute joint"))
        joint = get(connections, Symbol(joint_name), nothing)
        joint isa SpatialRevoluteJoint || throw(ArgumentError(
            "revolute friction '$name'.joint must name an existing " *
            "revolute joint"))
        stiffness = finite_number(get(table, "stiffness", NaN),
            "revolute friction '$name'.stiffness")
        damping = finite_number(get(table, "damping", 0.0),
            "revolute friction '$name'.damping")
        radius = finite_number(get(table, "effective_radius", NaN),
            "revolute friction '$name'.effective_radius")
        preload = finite_number(get(table, "preload", 0.0),
            "revolute friction '$name'.preload")
        mu_static = finite_number(get(table, "static_coefficient", NaN),
            "revolute friction '$name'.static_coefficient")
        mu_dynamic = finite_number(get(table, "dynamic_coefficient", NaN),
            "revolute friction '$name'.dynamic_coefficient")
        speed = finite_number(get(table, "transition_speed", 0.1),
            "revolute friction '$name'.transition_speed")
        release = finite_number(get(table, "release_time", 0.01),
            "revolute friction '$name'.release_time")
        stiffness > 0 && damping >= 0 && radius > 0 && preload >= 0 &&
            mu_static >= mu_dynamic >= 0 && speed > 0 && release > 0 ||
            throw(ArgumentError("revolute friction '$name' requires " *
                "positive stiffness and effective_radius, nonnegative " *
                "damping and preload, static_coefficient >= " *
                "dynamic_coefficient >= 0, positive transition_speed " *
                "and release_time"))
        friction = allocated_spatial_revolute_friction(layout, name, joint,
            stiffness, damping, radius, preload, mu_static, mu_dynamic,
            speed, release)
        initialize_spatial_revolute_friction!(initial, friction;
            reset_anchor = true)
        forces[name] = friction
    end
    for name in translational_friction_names
        table = typed[name]
        joint_name = get(table, "joint", nothing)
        joint_name isa AbstractString || throw(ArgumentError(
            "translational friction '$name'.joint must name an inline constraint"))
        joint = get(connections, Symbol(joint_name), nothing)
        joint isa SpatialInlineConstraint || throw(ArgumentError(
            "translational friction '$name'.joint must name an existing " *
            "inline constraint"))
        stiffness = finite_number(get(table, "stiffness", NaN),
            "translational friction '$name'.stiffness")
        damping = finite_number(get(table, "damping", 0.0),
            "translational friction '$name'.damping")
        preload = finite_number(get(table, "preload", 0.0),
            "translational friction '$name'.preload")
        mu_static = finite_number(get(table, "static_coefficient", NaN),
            "translational friction '$name'.static_coefficient")
        mu_dynamic = finite_number(get(table, "dynamic_coefficient", NaN),
            "translational friction '$name'.dynamic_coefficient")
        speed = finite_number(get(table, "transition_speed", 0.01),
            "translational friction '$name'.transition_speed")
        release = finite_number(get(table, "release_time", 0.01),
            "translational friction '$name'.release_time")
        stiffness > 0 && damping >= 0 && preload >= 0 &&
            mu_static >= mu_dynamic >= 0 && speed > 0 && release > 0 ||
            throw(ArgumentError("translational friction '$name' requires " *
                "positive stiffness, nonnegative damping and preload, " *
                "static_coefficient >= dynamic_coefficient >= 0, " *
                "positive transition_speed and release_time"))
        friction = allocated_spatial_translational_friction(layout, name,
            joint, stiffness, damping, preload, mu_static, mu_dynamic,
            speed, release)
        initialize_spatial_translational_friction!(initial, friction;
            reset_anchor = true)
        forces[name] = friction
    end
    for name in inplane_friction_names
        table = typed[name]
        constraint_name = get(table, "constraint", nothing)
        constraint_name isa AbstractString || throw(ArgumentError(
            "inplane friction '$name'.constraint must name an inplane constraint"))
        constraint = get(connections, Symbol(constraint_name), nothing)
        constraint isa SpatialInplaneConstraint || throw(ArgumentError(
            "inplane friction '$name'.constraint must name an existing " *
            "inplane constraint"))
        stiffness = finite_number(get(table, "stiffness", NaN),
            "inplane friction '$name'.stiffness")
        damping = finite_number(get(table, "damping", 0.0),
            "inplane friction '$name'.damping")
        preload = finite_number(get(table, "preload", 0.0),
            "inplane friction '$name'.preload")
        mu_static = finite_number(get(table, "static_coefficient", NaN),
            "inplane friction '$name'.static_coefficient")
        mu_dynamic = finite_number(get(table, "dynamic_coefficient", NaN),
            "inplane friction '$name'.dynamic_coefficient")
        speed = finite_number(get(table, "transition_speed", 0.01),
            "inplane friction '$name'.transition_speed")
        release = finite_number(get(table, "release_time", 0.01),
            "inplane friction '$name'.release_time")
        stiffness > 0 && damping >= 0 && preload >= 0 &&
            mu_static >= mu_dynamic >= 0 && speed > 0 && release > 0 ||
            throw(ArgumentError("inplane friction '$name' requires " *
                "positive stiffness, nonnegative damping and preload, " *
                "static_coefficient >= dynamic_coefficient >= 0, " *
                "positive transition_speed and release_time"))
        friction = allocated_spatial_inplane_friction(layout, name,
            constraint, stiffness, damping, preload, mu_static, mu_dynamic,
            speed, release)
        initialize_spatial_inplane_friction!(initial, friction;
            reset_anchor = true)
        forces[name] = friction
    end
    for name in tire_names
        table = typed[name]
        endpoints = get(table, "markers", nothing)
        endpoints isa Vector && length(endpoints) == 2 || throw(ArgumentError(
            "rolling tire '$name'.markers must contain a wheel-center " *
            "marker and a road marker"))
        wheel_marker = required_marker(markers, endpoints[1],
            "rolling tire '$name'")
        road_marker = required_marker(markers, endpoints[2],
            "rolling tire '$name'")
        wheel_marker isa SpatialBodyMarker || throw(ArgumentError(
            "rolling tire '$name' wheel marker must belong to a body"))
        road_marker isa SpatialBodyMarker &&
            road_marker.body === wheel_marker.body && throw(ArgumentError(
                "rolling tire '$name' markers must belong to different bodies"))
        radius = finite_number(get(table, "radius", NaN),
            "rolling tire '$name'.radius")
        radius > 0 || throw(ArgumentError(
            "rolling tire '$name'.radius must be positive"))
        regularization_speed = finite_number(
            get(table, "regularization_speed", 0.1),
            "rolling tire '$name'.regularization_speed")
        regularization_speed > 0 || throw(ArgumentError(
            "rolling tire '$name'.regularization_speed must be positive"))
        bristle = lowercase(string(get(table, "tangential_model",
            "expression"))) == "bristle"
        if !bristle
            any(haskey(table, field) for field in
                ("patch_length_by_load", "cornering_stiffness_by_load",
                 "longitudinal_slip_stiffness_by_load",
                 "longitudinal_relaxation_fraction",
                 "lateral_relaxation_fraction", "shear_release_time")) &&
                throw(ArgumentError("rolling tire '$name' bristle settings " *
                    "require tangential_model = 'bristle'"))
        end
        longitudinal_relaxation_length = !bristle && transient_tires[name] ?
            finite_number(table["longitudinal_relaxation_length"],
                "rolling tire '$name'.longitudinal_relaxation_length") : 0.0
        lateral_relaxation_length = !bristle && transient_tires[name] ?
            finite_number(table["lateral_relaxation_length"],
                "rolling tire '$name'.lateral_relaxation_length") : 0.0
        !bristle && transient_tires[name] &&
            min(longitudinal_relaxation_length,
                lateral_relaxation_length) <= 0 && throw(ArgumentError(
                "rolling tire '$name' relaxation lengths must be positive"))

        has_normal_expression = haskey(table, "normal_expression")
        has_normal_stiffness = haskey(table, "normal_stiffness")
        xor(has_normal_expression, has_normal_stiffness) ||
            throw(ArgumentError("rolling tire '$name' requires exactly one " *
                "of normal_expression or normal_stiffness"))
        normal_stiffness = normal_damping = 0.0
        normal_law = if has_normal_expression
            any(haskey(table, field) for field in
                ("normal_damping", "normal_damping_time_scale")) &&
                throw(ArgumentError("rolling tire '$name' normal damping " *
                    "fields require normal_stiffness"))
            compile_spatial_model_expression(
                string(table["normal_expression"]), parameters, layout)
        else
            normal_stiffness = finite_number(table["normal_stiffness"],
                "rolling tire '$name'.normal_stiffness")
            normal_stiffness > 0 || throw(ArgumentError(
                "rolling tire '$name'.normal_stiffness must be positive"))
            has_damping = haskey(table, "normal_damping")
            has_time_scale = haskey(table, "normal_damping_time_scale")
            has_damping && has_time_scale && throw(ArgumentError(
                "rolling tire '$name' may specify normal_damping or " *
                "normal_damping_time_scale, but not both"))
            normal_damping = has_damping ? finite_number(
                table["normal_damping"],
                "rolling tire '$name'.normal_damping") :
                normal_stiffness * finite_number(
                    get(table, "normal_damping_time_scale", 0.0),
                    "rolling tire '$name'.normal_damping_time_scale")
            normal_damping >= 0 || throw(ArgumentError(
                "rolling tire '$name'.normal_damping must be nonnegative"))
            spatial_linear_tire_normal_law(layout, name,
                normal_stiffness, normal_damping)
        end
        if bristle
            !has_normal_expression || throw(ArgumentError(
                "rolling tire '$name' bristle mode requires normal_stiffness"))
            any(haskey(table, field) for field in
                ("longitudinal_expression", "lateral_expression")) &&
                throw(ArgumentError("rolling tire '$name' bristle mode " *
                    "cannot use tangential expressions"))
            lowercase(string(get(table, "friction_limit", "ellipse"))) ==
                "ellipse" ||
                throw(ArgumentError("rolling tire '$name' bristle mode " *
                    "requires friction_limit = 'ellipse'"))
        else
            haskey(table, "longitudinal_expression") || throw(ArgumentError(
                "rolling tire '$name' requires longitudinal_expression"))
            haskey(table, "lateral_expression") || throw(ArgumentError(
                "rolling tire '$name' requires lateral_expression"))
        end
        longitudinal_law = bristle ? constant_scalar_law(0.0) :
            compile_spatial_model_expression(
                string(table["longitudinal_expression"]), parameters, layout)
        lateral_law = bristle ? constant_scalar_law(0.0) :
            compile_spatial_model_expression(
                string(table["lateral_expression"]), parameters, layout)

        friction_specification = lowercase(string(
            get(table, "friction_limit", bristle ? "ellipse" : "none")))
        friction_specification in ("none", "ellipse") ||
            throw(ArgumentError("rolling tire '$name'.friction_limit must " *
                "be 'none' or 'ellipse'"))
        friction_limit = Symbol(friction_specification)
        has_mu_longitudinal = haskey(table, "mu_longitudinal")
        has_mu_lateral = haskey(table, "mu_lateral")
        if friction_limit == :ellipse
            has_mu_longitudinal && has_mu_lateral || throw(ArgumentError(
                "rolling tire '$name' ellipse limit requires " *
                "mu_longitudinal and mu_lateral"))
        elseif has_mu_longitudinal || has_mu_lateral
            throw(ArgumentError("rolling tire '$name' friction coefficients " *
                "require friction_limit = 'ellipse'"))
        end
        mu_longitudinal = friction_limit == :ellipse ? finite_number(
            table["mu_longitudinal"],
            "rolling tire '$name'.mu_longitudinal") : Inf
        mu_lateral = friction_limit == :ellipse ? finite_number(
            table["mu_lateral"],
            "rolling tire '$name'.mu_lateral") : Inf
        friction_limit == :none || min(mu_longitudinal, mu_lateral) > 0 ||
            throw(ArgumentError("rolling tire '$name' friction " *
            "coefficients must be positive"))

        patch_curve = bristle ? tire_load_curve(table,
            "patch_length_by_load", name) : NTuple{2,Float64}[]
        longitudinal_curve = bristle ? tire_load_curve(table,
            "longitudinal_slip_stiffness_by_load", name) :
            NTuple{2,Float64}[]
        lateral_curve = bristle ? tire_load_curve(table,
            "cornering_stiffness_by_load", name) : NTuple{2,Float64}[]
        longitudinal_fraction = finite_number(get(table,
            "longitudinal_relaxation_fraction", 1.0),
            "rolling tire '$name'.longitudinal_relaxation_fraction")
        lateral_fraction = finite_number(get(table,
            "lateral_relaxation_fraction", 1.0),
            "rolling tire '$name'.lateral_relaxation_fraction")
        shear_release_time = finite_number(get(table,
            "shear_release_time", 0.01),
            "rolling tire '$name'.shear_release_time")
        min(longitudinal_fraction, lateral_fraction, shear_release_time) > 0 ||
            throw(ArgumentError("rolling tire '$name' bristle fractions " *
                "and release time must be positive"))

        initial_axes = spatial_marker_orientation(wheel_marker, initial)[:, 3]
        initial_normal = spatial_marker_orientation(road_marker, initial)[:, 3]
        norm(cross(initial_axes, initial_normal)) > sqrt(eps(Float64)) ||
            throw(ArgumentError("rolling tire '$name' axle must not be " *
                "parallel to the road normal"))
        component = allocated_spatial_tire(layout, name, wheel_marker,
            road_marker, radius, regularization_speed, normal_law,
            longitudinal_law, lateral_law; normal_stiffness,
            normal_damping, normal_expression = has_normal_expression,
            friction_limit, mu_longitudinal, mu_lateral,
            longitudinal_relaxation_length, lateral_relaxation_length,
            bristle, patch_curve, longitudinal_curve, lateral_curve,
            longitudinal_fraction, lateral_fraction, shear_release_time)
        initialize_spatial_tire!(initial, component, simulation.start_time)
        forces[name] = component
    end
    for name in bushing_names
        table = typed[name]
        endpoints = get(table, "markers", nothing)
        endpoints isa Vector && length(endpoints) == 2 || throw(ArgumentError(
            "bushing '$name'.markers must contain two marker names"))
        marker_i = required_marker(markers, endpoints[1], "bushing '$name'")
        marker_j = required_marker(markers, endpoints[2], "bushing '$name'")
        marker_i isa SpatialGroundMarker &&
            marker_j isa SpatialGroundMarker && throw(ArgumentError(
                "bushing '$name' cannot connect ground to ground"))
        marker_i isa SpatialBodyMarker && marker_j isa SpatialBodyMarker &&
            marker_i.body === marker_j.body && throw(ArgumentError(
                "bushing '$name' markers must belong to different bodies"))

        translational_stiffness = numeric_vector(table,
            "translational_stiffness", 3; label = "bushing '$name'")
        rotational_stiffness = numeric_vector(table,
            "rotational_stiffness", 3; default = zeros(3),
            label = "bushing '$name'")
        damping_time_scale = finite_number(
            get(table, "damping_time_scale", 0.0),
            "bushing '$name'.damping_time_scale")
        translational_damping = haskey(table, "translational_damping") ?
            numeric_vector(table, "translational_damping", 3;
                label = "bushing '$name'") :
            damping_time_scale .* translational_stiffness
        rotational_damping = haskey(table, "rotational_damping") ?
            numeric_vector(table, "rotational_damping", 3;
                label = "bushing '$name'") :
            damping_time_scale .* rotational_stiffness
        all(>=(0), translational_stiffness) &&
            all(>=(0), rotational_stiffness) &&
            all(>=(0), translational_damping) &&
            all(>=(0), rotational_damping) &&
            damping_time_scale >= 0 || throw(ArgumentError(
                "bushing '$name' stiffness, damping, and " *
                "damping_time_scale must be nonnegative"))
        active_during = active_during_stages(table, "bushing '$name'")

        reaction_marker = nothing
        if marker_j isa SpatialBodyMarker
            generated_name = Symbol(name, ".reaction")
            haskey(markers, generated_name) && throw(ArgumentError(
                "bushing '$name' generated marker '$generated_name' " *
                "already exists"))
            reaction_marker = SpatialFloatingMarker(
                generated_name, marker_j.body, marker_i)
            markers[generated_name] = reaction_marker
        end
        component = allocated_spatial_bushing(layout, name, marker_i,
            marker_j, reaction_marker, translational_stiffness,
            translational_damping, rotational_stiffness,
            rotational_damping, damping_time_scale, active_during)
        initialize_spatial_bushing!(initial, component)
        forces[name] = component
    end
    for name in applied_torque_names
        table = typed[name]
        active_during = active_during_stages(
            table, "applied torque '$name'")
        joint_name = torque_joint_names[name]
        hinge = connection_hinge(connections[joint_name])
        has_constant = haskey(table, "torque")
        has_expression = haskey(table, "expression")
        has_linear = haskey(table, "stiffness")
        count(identity, (has_constant, has_expression, has_linear)) == 1 ||
            throw(ArgumentError("applied torque '$name' requires exactly " *
                "one of torque, expression, or spring-damper coefficients"))
        spring_fields = ("damping", "damping_time_scale", "free_angle")
        !has_linear && any(haskey(table, field) for field in spring_fields) &&
            throw(ArgumentError("applied torque '$name' damping, " *
                "damping_time_scale, and free_angle require stiffness"))
        stiffness = damping = 0.0
        initial_free_angle = false
        free_angle = Ref(0.0)
        law = if has_constant
            constant_scalar_law(finite_number(table["torque"],
                "applied torque '$name'.torque"))
        elseif has_expression
            compile_spatial_model_expression(string(table["expression"]),
                parameters, layout)
        else
            stiffness = finite_number(table["stiffness"],
                "applied torque '$name'.stiffness")
            has_damping = haskey(table, "damping")
            has_time_scale = haskey(table, "damping_time_scale")
            has_damping && has_time_scale && throw(ArgumentError(
                "applied torque '$name' may specify damping or " *
                "damping_time_scale, but not both"))
            damping_time_scale = has_time_scale ? finite_number(
                table["damping_time_scale"],
                "applied torque '$name'.damping_time_scale") : 0.0
            damping = has_damping ? finite_number(table["damping"],
                "applied torque '$name'.damping") :
                stiffness * damping_time_scale
            min(stiffness, damping, damping_time_scale) >= 0 ||
                throw(ArgumentError("applied torque '$name' spring-damper " *
                    "coefficients must be nonnegative"))
            free_specification = get(table, "free_angle", 0.0)
            free_angle_reference = if free_specification isa AbstractString &&
                    lowercase(strip(free_specification)) == "initial"
                initial_free_angle = true
                free_angle
            else
                Ref(angle_value(free_specification,
                    "applied torque '$name'.free_angle"))
            end
            free_angle = free_angle_reference
            spatial_linear_torque_law(hinge, stiffness, damping,
                free_angle_reference)
        end
        component = allocated_spatial_applied_torque(layout, name, hinge,
            law; initial_free_angle, stiffness, damping, free_angle,
            active_during)
        initialize_spatial_applied_torque!(initial, component,
            simulation.start_time; reset_free_angle = true)
        forces[name] = component
    end
    for name in applied_force_names
        table = typed[name]
        active_during = active_during_stages(
            table, "applied force '$name'")
        endpoints = get(table, "markers", nothing)
        endpoints isa Vector && length(endpoints) == 2 || throw(ArgumentError(
            "applied force '$name'.markers must contain an application " *
            "marker and a direction marker"))
        application = required_marker(markers, endpoints[1],
            "applied force '$name'")
        application isa SpatialBodyMarker || throw(ArgumentError(
            "applied force '$name' application marker must belong to a body"))
        direction = required_marker(markers, endpoints[2],
            "applied force '$name'")
        has_force = haskey(table, "force")
        has_expression = haskey(table, "expression")
        xor(has_force, has_expression) || throw(ArgumentError(
            "applied force '$name' requires exactly one of force or expression"))
        law = if has_force
            constant_scalar_law(finite_number(table["force"],
                "applied force '$name'.force"))
        else
            compile_spatial_model_expression(string(table["expression"]),
                parameters, layout)
        end
        reaction_marker = nothing
        if haskey(table, "reaction_body")
            reaction_name = Symbol(table["reaction_body"])
            haskey(bodies, reaction_name) || throw(ArgumentError(
                "applied force '$name' names unknown reaction body " *
                "'$(table["reaction_body"])'"))
            reaction_body = bodies[reaction_name]
            reaction_body === application.body && throw(ArgumentError(
                "applied force '$name' reaction_body must differ from the " *
                "application body"))
            generated_name = Symbol(name, ".reaction")
            haskey(markers, generated_name) && throw(ArgumentError(
                "applied force '$name' generated marker '$generated_name' " *
                "already exists"))
            reaction_marker = SpatialFloatingMarker(
                generated_name, reaction_body, application)
            markers[generated_name] = reaction_marker
        end
        component = allocated_spatial_applied_force(layout, name,
            application, direction, reaction_marker, law, active_during)
        initialize_spatial_applied_force!(initial, component,
            simulation.start_time)
        forces[name] = component
    end
    for name in spanning_force_names
        table = typed[name]
        active_during = active_during_stages(
            table, "spanning force '$name'")
        endpoints = get(table, "markers", nothing)
        endpoints isa Vector && length(endpoints) == 2 || throw(ArgumentError(
            "spanning force '$name'.markers must contain two marker names"))
        marker_1 = required_marker(markers, endpoints[1],
            "spanning force '$name'")
        marker_2 = required_marker(markers, endpoints[2],
            "spanning force '$name'")
        marker_1 isa SpatialGroundMarker &&
            marker_2 isa SpatialGroundMarker && throw(ArgumentError(
                "spanning force '$name' cannot connect ground to ground"))
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
            value = finite_number(table["force"],
                "spanning force '$name'.force")
            constant_scalar_law(value)
        elseif has_expression
            compile_spatial_model_expression(string(table["expression"]),
                parameters, layout)
        else
            haskey(table, "damping") || throw(ArgumentError(
                "spanning force '$name' spring-damper law requires damping"))
            haskey(table, "free_length") || throw(ArgumentError(
                "spanning force '$name' spring-damper law requires free_length"))
            stiffness = finite_number(table["stiffness"],
                "spanning force '$name'.stiffness")
            damping = finite_number(table["damping"],
                "spanning force '$name'.damping")
            free_length = finite_number(table["free_length"],
                "spanning force '$name'.free_length")
            min(stiffness, damping, free_length) >= 0 || throw(ArgumentError(
                "spanning force '$name' spring-damper coefficients must be nonnegative"))
            spatial_linear_spanning_law(layout, name, stiffness, damping,
                free_length)
        end
        force = allocated_spatial_spanning_force(layout, name,
            marker_1, marker_2, law; stiffness, damping, free_length,
            active_during)
        initialize_spatial_spanning_force!(initial, force,
            simulation.start_time)
        forces[name] = force
    end
    gravity_names = sort!([name for (name, table) in typed
        if table["type"] == "gravity"])
    for name in gravity_names
        table = typed[name]
        acceleration = numeric_vector(table, "acceleration", 3;
            label = "gravity '$name'")
        body_list = get(table, "bodies", nothing)
        body_list isa Vector && !isempty(body_list) &&
            all(item -> item isa AbstractString, body_list) ||
            throw(ArgumentError("gravity '$name'.bodies must name one or more bodies"))
        for (index, body_name) in enumerate(body_list)
            body_symbol = Symbol(body_name)
            haskey(bodies, body_symbol) || throw(ArgumentError(
                "gravity '$name' names unknown body '$body_name'"))
            force_name = length(body_list) == 1 ? name : Symbol(name, ".", index)
            forces[force_name] = SpatialGravityComponent(
                force_name, bodies[body_symbol], acceleration)
        end
    end

    blocks = ExecutableEquationBlock[]
    contributions = EquationContribution[]
    for body in values(bodies)
        append!(blocks, executable_blocks(body))
    end
    for connection in values(connections)
        append!(blocks, executable_blocks(connection))
        append!(contributions, equation_contributions(connection))
    end
    for driver in values(drivers)
        append!(blocks, executable_blocks(driver))
        append!(contributions, equation_contributions(driver))
    end
    for measure in values(measures)
        append!(blocks, executable_blocks(measure))
        append!(contributions, equation_contributions(measure))
    end
    for force in values(forces)
        append!(blocks, executable_blocks(force))
        append!(contributions, equation_contributions(force))
    end
    for component in values(equation_components)
        append!(blocks, executable_blocks(component))
    end
    model = ExecutableAnalysisModel(layout.catalog, blocks, contributions)

    if initial_conditions.enabled
        for variable in layout.catalog.variables
            qualified = Symbol(variable.component, :., variable.name)
            qualified in transferred_variables || continue
            initial[variable.index] = saved_initial_values[variable.index]
        end
        normalize_transferred_orientations!(initial, bodies,
            initial_conditions.transferred_variables)
        for name in body_names
            apply_body_initial_impose!(initial, imposed_variables,
                imposed_orientations, bodies[name], typed[name], name)
        end
        for (name, specification) in relative_initial_specifications
            hinge = connection_hinge(connections[name])
            if !isnothing(hinge)
                _, velocity, coordinate = hinge.rotation_variables
            else
                inline = connection_inline(connections[name])
                _, velocity, coordinate = inline.translation_variables
            end
            apply_relative_initial_specification!(initial, variable_weights,
                imposed_variables, coordinate, velocity, specification)
        end
        for (name, driver) in drivers
            if driver isa SpatialRotationalMotionGenerator
                hinge = connection_hinge(connections[driver_joint_names[name]])
                alpha, omega, theta = hinge.rotation_variables
                initial[theta] = driver.motion(simulation.start_time)
                initial[omega] = driver.motion_derivative(simulation.start_time)
                initial[alpha] =
                    driver.motion_second_derivative(simulation.start_time)
            end
        end
    end

    # Retain the configuration supplied by the model (including any selected
    # saved-result values and imposed initial values) before the consistency
    # projections move bodies to satisfy the joint equations.
    entered_initial_values = copy(initial)

    velocity_indices = reduce(vcat,
        [[collect(body.velocity_variables);
          collect(body.angular_velocity_variables);
          (body isa SpatialFlexibleBeamComponent ?
              collect(body.elastic_velocity_variables) : Int[])]
         for body in values(bodies)]; init = Int[])
    append!(velocity_indices, relative_velocity_indices)
    position_rows = [equation.index for equation in layout.catalog.equations
        if equation.kind == :constraint && equation.level == 0]
    append!(position_rows, relative_position_rows)
    sort!(position_rows)
    velocity_rows = [equation.index for equation in layout.catalog.equations
        if equation.kind == :constraint && equation.level == 1]
    append!(velocity_rows, relative_velocity_rows)
    sort!(velocity_rows)
    position_corrections = spatial_position_projection!(initial, model,
        bodies, markers, position_rows, relative_coordinate_indices,
        variable_weights, initial_condition_weights, imposed_variables,
        imposed_orientations, layout,
        "spatial initial position"; time = simulation.start_time)
    for (name, specification) in relative_initial_specifications
        name in driven_connection_names && continue
        coordinate = if !isnothing(connection_hinge(connections[name]))
            connection_hinge(connections[name]).rotation_variables[3]
        else
            connection_inline(connections[name]).translation_variables[3]
        end
        coordinate_name = Symbol(name, :.,
            layout.catalog.variables[coordinate].name)
        (specification.coordinate_specified ||
            coordinate_name in transferred_variables) && continue
        hinge = connection_hinge(connections[name])
        if !isnothing(hinge)
            initial[hinge.rotation_variables[3]] = hinge_angle(hinge, initial)
        else
            inline = connection_inline(connections[name])
            initial[inline.translation_variables[3]] =
                inline_distance(inline, initial)
        end
    end
    velocity_corrections = weighted_constraint_projection!(initial, model,
        velocity_indices, velocity_rows, variable_weights,
        imposed_variables, layout,
        "spatial initial velocity"; time = simulation.start_time)
    for (name, specification) in relative_initial_specifications
        name in driven_connection_names && continue
        velocity = if !isnothing(connection_hinge(connections[name]))
            connection_hinge(connections[name]).rotation_variables[2]
        else
            connection_inline(connections[name]).translation_variables[2]
        end
        velocity_name = Symbol(name, :.,
            layout.catalog.variables[velocity].name)
        (specification.velocity_specified ||
            velocity_name in transferred_variables) && continue
        hinge = connection_hinge(connections[name])
        if !isnothing(hinge)
            initial[hinge.rotation_variables[2]] =
                hinge_angular_velocity(hinge, initial)
        else
            inline = connection_inline(connections[name])
            initial[inline.translation_variables[2]] =
                inline_velocity(inline, initial)
        end
    end
    for measure in values(measures)
        if measure isa SpatialSpanMeasure
            initialize_spatial_span_measure!(initial, measure)
        else
            initialize_spatial_directed_distance_measure!(initial, measure)
        end
    end
    for force in values(forces)
        if force isa SpatialBeltSpanComponent
            initialize_spatial_belt_span!(initial, force)
        elseif force isa SpatialSpanningForceComponent
            initialize_spatial_spanning_force!(initial, force,
                simulation.start_time)
        elseif force isa SpatialAppliedForceComponent
            initialize_spatial_applied_force!(initial, force,
                simulation.start_time)
        elseif force isa SpatialAppliedTorqueComponent
            initialize_spatial_applied_torque!(initial, force,
                simulation.start_time; reset_free_angle = true)
        elseif force isa SpatialBushingComponent
            initialize_spatial_bushing!(initial, force)
        elseif force isa SpatialPlaneContactComponent
            initialize_spatial_plane_contact!(initial, force)
        elseif force isa SpatialCurveContactComponent
            initialize_spatial_curve_contact!(initial, force;
                station = initial[force.station_variable])
        elseif force isa SpatialSurfaceFriction
            initialize_spatial_surface_friction!(initial, force;
                reset_anchor = true)
        elseif force isa SpatialRevoluteFriction
            initialize_spatial_revolute_friction!(initial, force;
                reset_anchor = true)
        elseif force isa SpatialTranslationalFriction
            initialize_spatial_translational_friction!(initial, force;
                reset_anchor = true)
        elseif force isa SpatialInplaneFriction
            initialize_spatial_inplane_friction!(initial, force;
                reset_anchor = true)
        elseif force isa SpatialTireComponent
            initialize_spatial_tire!(initial, force,
                simulation.start_time)
        end
    end

    imposed_names = [Symbol(string(variable.component, ".", variable.name))
        for variable in layout.catalog.variables
        if variable.index in imposed_variables]
    append!(imposed_names, [Symbol(string(name, ".orientation"))
        for name in imposed_orientations])
    sort!(imposed_names)
    selection = merge(state_selection(
        model, layout, bodies, markers, connections, drivers, initial,
        document),
        (; imposed_initial_variables = imposed_names))
    active_variables = setdiff(collect(eachindex(layout.catalog.variables)),
        selection.inactive_variables)
    active_equations = setdiff(collect(eachindex(layout.catalog.equations)),
        selection.inactive_equations)
    length(active_variables) == length(active_equations) || error(
        "spatial dynamic system must be square")

    acceleration_corrections = initialize_spatial_accelerations!(initial,
        model, layout, simulation.start_time;
        inactive_variables = selection.inactive_variables,
        inactive_equations = selection.inactive_equations)
    for connection in values(connections)
        hinge = connection_hinge(connection)
        if !isnothing(hinge) && !isempty(hinge.rotation_variables)
            initial[hinge.rotation_variables[1]] =
                hinge_angular_acceleration(hinge, initial)
        end
        inline = connection_inline(connection)
        if !isnothing(inline) && !isempty(inline.translation_variables)
            initial[inline.translation_variables[1]] =
                inline_acceleration(inline, initial)
        end
    end
    degrees_of_freedom = length(selection.selected_velocity_indices)
    analysis = analysis_settings(document, degrees_of_freedom)
    if any(force isa SpatialTireComponent && force.bristle
            for force in values(forces)) &&
            (analysis.mode == :static ||
             analysis.initialization == :static_equilibrium)
        analysis.static_method == :dynamic_relaxation &&
            !analysis.relaxation_polish || throw(ArgumentError(
                "bristle tires need static_method = 'dynamic_relaxation' " *
                "and relaxation_polish = false; static Newton cannot " *
                "determine parked shear from the final pose"))
    end
    initial_conditions = merge(initial_conditions,
        (; position_corrections, velocity_corrections,
         acceleration_corrections))
    LoadedSpatialModel(title, layout, model, bodies, body_reference_frames,
        markers, ground_names, measures, forces, equation_components,
        connections, drivers,
        analysis, simulation,
        selection, initial_conditions, initial_condition_weights,
        variable_weights, sort!(collect(imposed_variables)),
        sort!(collect(imposed_orientations)), entered_initial_values, initial,
        active_variables, active_equations, text)
end

end
