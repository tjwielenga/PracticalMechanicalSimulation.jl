"""
    PlanarModeling

Bridge between allocated canonical indices and the planar component library.
The `allocated_*` constructors turn a completed `ModelLayout` into executable
components. [`assemble_planar_model`](@ref) then combines their owned blocks
and additive contributions into one `ExecutableAnalysisModel`.

This layer contains no TOML policy: it can also be used by programmatic models
and tests that construct the same components directly.
"""
module PlanarModeling

using LinearAlgebra
using ..AutomaticAnalysis
using ..PlanarAppliedForces
using ..PlanarDirectedDistances
using ..PlanarComponentAssembly

export allocated_planar_body, planar_point_marker, planar_floating_point_marker,
       allocated_revolute_joint, allocated_inplane_constraint,
       allocated_perp_constraint,
       allocated_translational_joint, allocated_fixed_joint,
       allocated_gear_pair, allocated_rack_and_pinion,
       allocated_belt_span,
       allocated_rotational_motion_generator,
       allocated_translational_motion_generator,
       allocated_distance_coordinate, allocated_span_measure,
       allocated_coordinate_coupler,
       allocated_applied_force, allocated_applied_torque,
       allocated_planar_bushing,
       allocated_plane_contact, allocated_spanning_force,
       allocated_torsional_spring_damper,
       assemble_planar_model,
       solve_planar_analysis!, predict_planar_configuration!

"""Construct an applied scalar force with a component-local load equation."""
function allocated_applied_force(layout, name, application_marker,
        direction_body, direction_marker, reaction_marker, magnitude;
        active_during = (:static, :dynamic, :modal))
    axis = PlanarDirectedAxis(direction_body, direction_marker)
    PlanarAppliedForceComponent(name, application_marker, axis,
        reaction_marker, magnitude, active_during,
        Ref(:dynamic in active_during),
        only(component_variable_indices(layout, name)),
        only(component_equation_indices(layout, name, :load)))
end

"""Construct an applied scalar torque with a component-local load equation."""
function allocated_applied_torque(layout, name, marker_a, marker_b, torque,
        point_a = nothing, point_b = nothing;
        active_during = (:static, :dynamic, :modal))
    PlanarAppliedTorqueComponent(name, marker_a, marker_b, torque,
        point_a, point_b, active_during, Ref(:dynamic in active_during),
        only(component_variable_indices(layout, name)),
        only(component_equation_indices(layout, name, :load)))
end

"""Construct a planar body from its locations in a completed model layout."""
function allocated_planar_body(layout, name, mass, inertia;
        selected_state = false, selected_states = Symbol[],
        retain_state_candidates = false)
    variables = component_variable_indices(layout, name)
    selected = selected_state ? [:omega] : collect(selected_states)
    blocks = Dict(:V_x => :selected_V_x, :V_y => :selected_V_y,
                  :omega => :selected_state)
    state_equations = [component_equation_indices(layout, name, blocks[kind])
                       for kind in selected]
    state_variables = NTuple{3,Int}[]
    for kind in selected
        if kind == :V_x
            push!(state_variables, (variables[1], variables[4], variables[7]))
        elseif kind == :V_y
            push!(state_variables, (variables[2], variables[5], variables[8]))
        elseif kind == :omega
            push!(state_variables, (variables[3], variables[6], variables[9]))
        else
            throw(ArgumentError("unknown planar selected state '$kind'"))
        end
    end
    candidate_equations = retain_state_candidates ? Dict(
        kind => component_equation_indices(layout, name, blocks[kind])
        for kind in keys(blocks)) : Dict{Symbol,UnitRange{Int}}()
    candidate_variables = retain_state_candidates ? Dict(
        :V_x => (variables[1], variables[4], variables[7]),
        :V_y => (variables[2], variables[5], variables[8]),
        :omega => (variables[3], variables[6], variables[9])) :
        Dict{Symbol,NTuple{3,Int}}()
    return PlanarRigidBodyComponent(name, mass, inertia,
        variables[1:2], variables[3], variables[4:5], variables[6],
        variables[7:8], variables[9],
        component_equation_indices(layout, name, :balance), state_equations,
        state_variables, candidate_equations, candidate_variables)
end

"""Construct an allocated elastic tangent span belonging to a planar belt."""
function allocated_belt_span(layout, name, belt_name, pulley_1, pulley_2,
        tangent_kind, tangent_side, beta_1, beta_2, reference_length,
        reference_tangent, reference_angle_1, reference_angle_2,
        initial_extension, stiffness, damping, damping_time_scale)
    variables = component_variable_indices(layout, name)
    geometry = component_equation_indices(layout, name, :geometry)
    rate = only(component_equation_indices(layout, name, :rate))
    load = component_equation_indices(layout, name, :load)
    PlanarBeltSpanComponent(name, belt_name, pulley_1, pulley_2,
        tangent_kind, tangent_side, beta_1, beta_2, reference_length,
        collect(reference_tangent), reference_angle_1, reference_angle_2,
        initial_extension, stiffness, damping, damping_time_scale,
        variables[1:2], variables[3:4], variables[5:6], variables[7],
        variables[8], variables[9], variables[10], variables[11:12],
        first(geometry):last(geometry), rate, first(load):last(load))
end


"""Construct a planar bushing from its allocated load variables and equations."""
function allocated_planar_bushing(layout, name, marker_1, marker_2,
        translational_stiffness, translational_damping,
        rotational_stiffness, rotational_damping, damping_time_scale,
        free_position, free_angle;
        active_during = (:static, :dynamic, :modal))
    variables = component_variable_indices(layout, name)
    PlanarBushingComponent(name, marker_1, marker_2,
        collect(translational_stiffness), collect(translational_damping),
        rotational_stiffness, rotational_damping, damping_time_scale,
        collect(free_position), free_angle, active_during,
        Ref(:dynamic in active_during), variables[1:2], variables[3],
        component_equation_indices(layout, name, :load))
end

"""Construct a one-sided plane contact from its allocated variables."""
function allocated_plane_contact(layout, name, marker_1, marker_2, radius,
        stiffness, damping_factor;
        law = nothing, expression = false, transition_depth = 0.0,
        active_during = (:static, :dynamic, :modal))
    variables = component_variable_indices(layout, name)
    axis = PlanarDirectedAxis(marker_2.owner, marker_2.orientation)
    geometry = PlanarDirectedDistance(marker_1.owner, marker_1.point,
        marker_2.owner, marker_2.point, axis)
    PlanarPlaneContactComponent(name, marker_1, marker_2, geometry, law,
        Bool(expression), radius, stiffness, damping_factor,
        transition_depth, active_during,
        Ref(:dynamic in active_during),
        variables[1], variables[2], variables[3],
        variables[4:5], component_equation_indices(layout, name, :contact))
end

"""Construct an allocated nine-equation spanning force."""
function allocated_spanning_force(layout, name, marker_1, marker_2,
        law; stiffness = 0.0, damping = 0.0, free_length = 0.0,
        active_during = (:static, :dynamic, :modal))
    variables = component_variable_indices(layout, name)
    geometry = component_equation_indices(layout, name, :geometry)
    rate = only(component_equation_indices(layout, name, :rate))
    load = component_equation_indices(layout, name, :load)
    element = PlanarSpanningSpringDamper(marker_1, marker_2,
        stiffness, damping, free_length,
        variables[1:2], variables[3], variables[4:5], variables[6],
        variables[7], variables[8:9],
        geometry[1:2], geometry[3], geometry[4:5], rate,
        load[1], load[2:3])
    PlanarSpanningForceComponent(name, element, law, active_during,
        Ref(:dynamic in active_during))
end

"""Construct a reaction-free marker-to-marker span measurement."""
function allocated_span_measure(layout, name, marker_1, marker_2)
    variables = component_variable_indices(layout, name)
    geometry = component_equation_indices(layout, name, :geometry)
    PlanarSpanMeasureComponent(name, marker_1.owner, marker_1.point,
        marker_2.owner, marker_2.point,
        variables[1:2], variables[3], variables[4:5], variables[6],
        variables[7], geometry[1:2], geometry[3], geometry[4:5],
        only(component_equation_indices(layout, name, :velocity)),
        only(component_equation_indices(layout, name, :acceleration)))
end

"""Construct an allocated one-equation torsional spring-damper."""
function allocated_torsional_spring_damper(layout, name, marker_1, marker_2,
        stiffness, damping, damping_time_scale, free_angle)
    torque_variable = only(component_variable_indices(layout, name))
    constitutive_equation = only(
        component_equation_indices(layout, name, :load))
    element = PlanarTorsionalSpringDamper(
        marker_1.orientation, marker_2.orientation,
        stiffness, damping, free_angle,
        torque_variable, constitutive_equation)
    PlanarTorsionalSpringComponent(name, marker_1, marker_2,
        element, damping_time_scale)
end

"""Construct a point marker using a body's allocated variables and balances."""
function planar_point_marker(body, local_position)
    return PlanarBodyPointMarker(body.position_variables,
        body.orientation_variable, body.velocity_variables,
        body.angular_velocity_variable, body.balance_equations[1:2],
        body.balance_equations[3], collect(local_position))
end

"""Construct a body-owned point whose translation follows another marker."""
function planar_floating_point_marker(body, follower)
    PlanarFloatingPointMarker(body.position_variables,
        body.balance_equations[1:2], body.balance_equations[3], follower)
end

"""Construct a revolute joint from its locations in a completed layout."""
function allocated_revolute_joint(layout, name,
        body_a, marker_a, body_b, marker_b;
        rotation_coordinates = false, selected_state = false,
        rotation_marker_a = nothing, rotation_marker_b = nothing,
        retain_state_candidate = false)
    selected_state && !rotation_coordinates && throw(ArgumentError(
        "a revolute relative state requires rotation coordinates"))
    rotation_coordinates &&
        (isnothing(rotation_marker_a) || isnothing(rotation_marker_b)) &&
        throw(ArgumentError(
            "revolute rotation coordinates require two orientation markers"))
    variables = component_variable_indices(layout, name)
    rotation_variables = rotation_coordinates ? collect(variables[3:5]) : Int[]
    rotation_equations = rotation_coordinates ? [
        only(component_equation_indices(layout, name, :rotation_acceleration)),
        only(component_equation_indices(layout, name, :rotation_velocity)),
        only(component_equation_indices(layout, name, :rotation_position)),
    ] : Int[]
    state_equations = selected_state ? collect(component_equation_indices(
        layout, name, :selected_rotation_state)) : Int[]
    candidate_state_equations = rotation_coordinates && retain_state_candidate ?
        collect(component_equation_indices(
            layout, name, :selected_rotation_state)) : Int[]
    selected_variables = selected_state ?
        [(rotation_variables[1], rotation_variables[2],
          rotation_variables[3])] : NTuple{3,Int}[]
    return PlanarRevoluteJointComponent(name, body_a, marker_a,
        body_b, marker_b, first(variables):(first(variables) + 1),
        component_equation_indices(layout, name, :acceleration),
        component_equation_indices(layout, name, :velocity),
        component_equation_indices(layout, name, :position),
        rotation_marker_a, rotation_marker_b, rotation_variables,
        rotation_equations, state_equations,
        selected_variables, candidate_state_equations)
end

"""Construct an in-plane constraint from its allocated scalar locations."""
function allocated_inplane_constraint(layout, name,
        body_i, marker_i, body_j, marker_j, orientation_j)
    axis = PlanarDirectedAxis(body_j, orientation_j)
    geometry = PlanarDirectedDistance(
        body_i, marker_i, body_j, marker_j, axis)
    return PlanarInplaneConstraint(name, geometry,
        only(component_variable_indices(layout, name)),
        only(component_equation_indices(layout, name, :acceleration)),
        only(component_equation_indices(layout, name, :velocity)),
        only(component_equation_indices(layout, name, :position)))
end

"""Construct an allocated planar perpendicular-axis constraint."""
function allocated_perp_constraint(layout, name,
        body_i, orientation_i, body_j, orientation_j)
    PlanarPerpConstraint(name, body_i, orientation_i,
        body_j, orientation_j,
        only(component_variable_indices(layout, name)),
        only(component_equation_indices(layout, name, :acceleration)),
        only(component_equation_indices(layout, name, :velocity)),
        only(component_equation_indices(layout, name, :position)))
end

"""Construct a planar translational joint from allocated inplane and perp primitives."""
function allocated_translational_joint(layout, name,
        body_i, point_i, orientation_i, body_j, point_j, orientation_j)
    inplane_name = Symbol(name, ".inplane")
    perp_name = Symbol(name, ".perp")
    inplane = allocated_inplane_constraint(layout, inplane_name,
        body_i, point_i, body_j, point_j, orientation_j)
    perp = allocated_perp_constraint(layout, perp_name,
        body_i, orientation_i, body_j, orientation_j)
    PlanarTranslationalJoint(name, inplane, perp)
end

"""Construct a planar fixed joint from allocated revolute and perp primitives."""
function allocated_fixed_joint(layout, name,
        body_i, point_i, orientation_i, body_j, point_j, orientation_j)
    revolute_name = Symbol(name, ".revolute")
    perp_name = Symbol(name, ".perp")
    revolute = allocated_revolute_joint(layout, revolute_name,
        body_i, point_i, body_j, point_j)
    perp = allocated_perp_constraint(layout, perp_name,
        body_i, orientation_i, body_j, orientation_j)
    PlanarFixedJoint(name, revolute, perp)
end

"""Construct an allocated external or internal planar gear-pair constraint."""
function allocated_gear_pair(layout, name, body_1, marker_1, orientation_1,
        body_2, marker_2, orientation_2, joint_1, joint_2,
        carrier_body, contact_marker, force_reference_marker,
        carrier_orientation, radius_1, radius_2, phase)
    PlanarGearPairComponent(name, body_1, marker_1, orientation_1,
        body_2, marker_2, orientation_2, joint_1, joint_2,
        carrier_body, contact_marker, force_reference_marker,
        carrier_orientation, radius_1, radius_2, phase,
        only(component_variable_indices(layout, name)),
        only(component_equation_indices(layout, name, :acceleration)),
        only(component_equation_indices(layout, name, :velocity)),
        only(component_equation_indices(layout, name, :position)))
end

"""Construct an allocated ideal planar rack-and-pinion constraint."""
function allocated_rack_and_pinion(layout, name,
        rack_body, rack_marker, pinion_body, pinion_marker, pinion_orientation,
        translational_joint, revolute_joint, carrier_body, contact_marker,
        pinion_carrier_orientation, pitch_radius, phase)
    PlanarRackAndPinionComponent(name, rack_body, rack_marker,
        pinion_body, pinion_marker, pinion_orientation,
        translational_joint, revolute_joint, carrier_body, contact_marker,
        pinion_carrier_orientation, pitch_radius, phase,
        only(component_variable_indices(layout, name)),
        only(component_equation_indices(layout, name, :acceleration)),
        only(component_equation_indices(layout, name, :velocity)),
        only(component_equation_indices(layout, name, :position)))
end

"""Construct an allocated rotational motion generator."""
function allocated_rotational_motion_generator(layout, name,
        body_a, marker_a, body_b, marker_b,
        motion, motion_derivative, motion_second_derivative)
    variables = component_variable_indices(layout, name)
    return PlanarRotationalMotionGenerator(name, body_a, marker_a,
        body_b, marker_b, variables[1], variables[2], variables[3],
        variables[4], component_equation_indices(layout, name, :position),
        component_equation_indices(layout, name, :velocity),
        component_equation_indices(layout, name, :acceleration),
        motion, motion_derivative, motion_second_derivative)
end

"""Construct a prescribed distance along the second marker's local y axis."""
function allocated_translational_motion_generator(layout, name,
        body_i, marker_i, body_j, marker_j, orientation_j,
        motion, motion_derivative, motion_second_derivative)
    variables = component_variable_indices(layout, name)
    axis = PlanarDirectedAxis(body_j, orientation_j)
    geometry = PlanarDirectedDistance(
        body_i, marker_i, body_j, marker_j, axis)
    PlanarTranslationalMotionGenerator(name, geometry,
        variables[1], variables[2], variables[3], variables[4],
        component_equation_indices(layout, name, :position),
        component_equation_indices(layout, name, :velocity),
        component_equation_indices(layout, name, :acceleration),
        motion, motion_derivative, motion_second_derivative)
end

"""Construct an allocated, unprescribed relative-distance coordinate."""
function allocated_distance_coordinate(layout, name,
        body_i, marker_i, body_j, marker_j, orientation_j;
        selected_state = false, retain_state_candidate = false)
    variables = component_variable_indices(layout, name)
    state_equations = selected_state ? collect(component_equation_indices(
        layout, name, :selected_distance_state)) : Int[]
    candidate_state_equations = retain_state_candidate ? collect(
        component_equation_indices(layout, name, :selected_distance_state)) :
        Int[]
    selected_variables = selected_state ?
        [(variables[3], variables[2], variables[1])] : NTuple{3,Int}[]
    axis = PlanarDirectedAxis(body_j, orientation_j)
    geometry = PlanarDirectedDistance(
        body_i, marker_i, body_j, marker_j, axis)
    PlanarDistanceCoordinateComponent(name, geometry,
        variables[1], variables[2], variables[3],
        only(component_equation_indices(layout, name, :acceleration)),
        only(component_equation_indices(layout, name, :velocity)),
        only(component_equation_indices(layout, name, :position)),
        state_equations, selected_variables, candidate_state_equations)
end

"""Construct one allocated linear constraint among scalar coordinates."""
function allocated_coordinate_coupler(layout, name, coordinates,
        coefficients, offset, coordinate_kind)
    PlanarCoordinateCoupler(name, collect(coordinates), collect(coefficients),
        offset, coordinate_kind,
        only(component_variable_indices(layout, name)),
        only(component_equation_indices(layout, name, :acceleration)),
        only(component_equation_indices(layout, name, :velocity)),
        only(component_equation_indices(layout, name, :position)))
end

"""Combine owned component blocks and load contributions into one model."""
function assemble_planar_model(layout, equation_components,
                               contribution_components)
    blocks = ExecutableEquationBlock[]
    for component in equation_components
        append!(blocks, executable_blocks(component))
    end
    contributions = EquationContribution[]
    for component in contribution_components
        append!(contributions, equation_contributions(component))
    end
    return ExecutableAnalysisModel(layout.catalog, blocks, contributions)
end

"""Newton solve one square selected analysis in canonical storage."""
function solve_planar_analysis!(canonical, model::ExecutableAnalysisModel,
        selection::AnalysisSelection, t;
        tolerance = 1.0e-11, maximum_iterations = 15)
    values = analysis_values(canonical, selection)
    derivative = zeros(eltype(canonical), length(canonical))
    corrections = 0
    for iteration in 1:maximum_iterations
        expand_analysis_values!(canonical, values, selection)
        equations = zeros(eltype(canonical), length(selection.equation_indices))
        evaluate_analysis_equations!(equations, model, selection,
            t, canonical, derivative)
        norm(equations, Inf) <= tolerance && break
        jacobian = evaluate_analysis_sparse_jacobian(model, selection,
            t, canonical, derivative, 0.0)
        values .-= jacobian \ equations
        corrections = iteration
    end
    expand_analysis_values!(canonical, values, selection)
    equations = zeros(eltype(canonical), length(selection.equation_indices))
    evaluate_analysis_equations!(equations, model, selection,
        t, canonical, derivative)
    norm(equations, Inf) <= tolerance ||
        error("$(typeof(selection.policy)) failed to converge")
    return corrections
end

function solve_planar_analysis!(canonical, assembly,
        policy::AutomaticAnalysis.AnalysisPolicy, t; kwargs...)
    selection = select_analysis(assembly.layout.catalog, policy)
    return solve_planar_analysis!(canonical, assembly.model, selection, t;
        kwargs...)
end

"""Second-order prediction of every supplied body configuration."""
function predict_planar_configuration!(canonical, bodies, step)
    for body in bodies
        canonical[body.position_variables] .+=
            step .* canonical[body.velocity_variables] .+
            (step^2 / 2) .* canonical[body.acceleration_variables]
        canonical[body.orientation_variable] +=
            step * canonical[body.angular_velocity_variable] +
            (step^2 / 2) * canonical[body.angular_acceleration_variable]
    end
    return canonical
end

end
