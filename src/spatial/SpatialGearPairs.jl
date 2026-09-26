"""Ideal spatial gear pairs with carrier-relative phase coordinates."""
module SpatialGearPairs

using LinearAlgebra
using ..AutomaticAnalysis
using ..SpatialComponentAssembly
using ..SpatialModeling
using ..SpatialDirectedDistances
using ..SpatialConstraints
using ..SpatialConstraints: marker_angular_kinematics

import ..SpatialComponentAssembly: component_registration,
    executable_blocks, equation_contributions

export SpatialGearPair, spatial_gear_pair_registration,
       allocated_spatial_gear_pair, initialize_spatial_gear_coordinates!,
       spatial_gear_contact_direction, spatial_gear_position,
       spatial_gear_velocity, spatial_gear_acceleration

struct SpatialGearSide{J,B,V,T}
    joint::J
    body::B
    axis_carrier::V
    radial_carrier::V
    radius::T
end

"""
An ideal gear pair whose pitch geometry is fixed in the contact marker's
carrier frame. The two floating markers apply opposite forces at the carrier
contact point.
"""
struct SpatialGearPair{S1,S2,C,F1,F2,V,T}
    name::Symbol
    side_1::S1
    side_2::S2
    contact_marker::C
    contact_1::F1
    contact_2::F2
    direction_carrier::V
    phase::T
    reaction_variable::Int
    acceleration_variables::UnitRange{Int}
    velocity_variables::UnitRange{Int}
    position_variables::UnitRange{Int}
    acceleration_definition_equations::UnitRange{Int}
    velocity_definition_equations::UnitRange{Int}
    position_definition_equations::UnitRange{Int}
    acceleration_equation::Int
    velocity_equation::Int
    position_equation::Int
end

function spatial_gear_pair_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:lambda, :reaction, 2),
        VariableDeclaration(:alpha_1, :relative_acceleration, 2),
        VariableDeclaration(:alpha_2, :relative_acceleration, 2),
        VariableDeclaration(:omega_1, :relative_velocity, 1),
        VariableDeclaration(:omega_2, :relative_velocity, 1),
        VariableDeclaration(:theta_1, :relative_position, 0),
        VariableDeclaration(:theta_2, :relative_position, 0),
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:phase_acceleration, [
            EquationDeclaration(:alpha_1_definition, :coordinate_relation,
                2, :carrier_relative_rotation),
            EquationDeclaration(:alpha_2_definition, :coordinate_relation,
                2, :carrier_relative_rotation)]),
        EquationBlockDeclaration(:phase_velocity, [
            EquationDeclaration(:omega_1_definition, :coordinate_relation,
                1, :carrier_relative_rotation),
            EquationDeclaration(:omega_2_definition, :coordinate_relation,
                1, :carrier_relative_rotation)]),
        EquationBlockDeclaration(:phase_position, [
            EquationDeclaration(:theta_1_definition, :coordinate_relation,
                0, :carrier_relative_rotation),
            EquationDeclaration(:theta_2_definition, :coordinate_relation,
                0, :carrier_relative_rotation)]),
        EquationBlockDeclaration(:acceleration, [
            EquationDeclaration(:Phi_ddot, :constraint, 2, :gear_pair)]),
        EquationBlockDeclaration(:velocity, [
            EquationDeclaration(:Phi_dot, :constraint, 1, :gear_pair)]),
        EquationBlockDeclaration(:position, [
            EquationDeclaration(:Phi, :constraint, 0, :gear_pair)]),
    ]
    families = ConstraintFamilyDeclaration[
        ConstraintFamilyDeclaration(:gear_pair, :lambda, :Phi, :Phi_dot,
            :Phi_ddot),
    ]
    ComponentRegistration(name, variables, blocks, families)
end

component_registration(gear::SpatialGearPair) =
    spatial_gear_pair_registration(gear.name)

marker_owner(::SpatialGroundMarker) = nothing
marker_owner(marker::SpatialBodyMarker) = marker.body

function carrier_vector(marker, local_vector, z)
    spatial_marker_orientation(marker, z) * local_vector
end

function carrier_angular_kinematics(marker, z)
    marker isa SpatialGroundMarker && return (
        omega = zeros(eltype(z), 3), alpha = zeros(eltype(z), 3))
    body = marker.body
    orientation = rotation_matrix(@view z[body.euler_parameter_variables])
    (; omega = orientation * z[body.angular_velocity_variables],
       alpha = orientation * z[body.angular_acceleration_variables])
end

function physical_gear_phase(gear, side, z)
    axis = carrier_vector(gear.contact_marker, side.axis_carrier, z)
    radial = carrier_vector(gear.contact_marker, side.radial_carrier, z)
    gear_x = marker_axis_kinematics(side.joint.marker_a, z, 1).direction
    sine = dot(axis, cross(radial, gear_x))
    cosine = dot(radial, gear_x)
    atan(sine, cosine)
end

function physical_gear_velocity(gear, side, z)
    axis = carrier_vector(gear.contact_marker, side.axis_carrier, z)
    gear_motion = marker_angular_kinematics(side.joint.marker_a, z)
    carrier_motion = carrier_angular_kinematics(gear.contact_marker, z)
    dot(gear_motion.omega - carrier_motion.omega, axis)
end

function physical_gear_acceleration(gear, side, z)
    axis = carrier_vector(gear.contact_marker, side.axis_carrier, z)
    gear_motion = marker_angular_kinematics(side.joint.marker_a, z)
    carrier_motion = carrier_angular_kinematics(gear.contact_marker, z)
    axis_velocity = cross(carrier_motion.omega, axis)
    dot(gear_motion.alpha - carrier_motion.alpha, axis) +
        dot(gear_motion.omega - carrier_motion.omega, axis_velocity)
end

spatial_gear_contact_direction(gear::SpatialGearPair, z) =
    carrier_vector(gear.contact_marker, gear.direction_carrier, z)

function spatial_gear_position(gear::SpatialGearPair, z)
    radii = (gear.side_1.radius, gear.side_2.radius)
    radii[1] * z[gear.position_variables[1]] -
        radii[2] * z[gear.position_variables[2]] - gear.phase
end

function spatial_gear_velocity(gear::SpatialGearPair, z)
    radii = (gear.side_1.radius, gear.side_2.radius)
    radii[1] * physical_gear_velocity(gear, gear.side_1, z) -
        radii[2] * physical_gear_velocity(gear, gear.side_2, z)
end

function spatial_gear_acceleration(gear::SpatialGearPair, z)
    radii = (gear.side_1.radius, gear.side_2.radius)
    radii[1] * physical_gear_acceleration(gear, gear.side_1, z) -
        radii[2] * physical_gear_acceleration(gear, gear.side_2, z)
end

function body_rotational_dependencies(body)
    [collect(body.angular_acceleration_variables);
     collect(body.angular_velocity_variables);
     collect(body.euler_parameter_variables)]
end

function phase_dependencies(gear)
    columns = [body_rotational_dependencies(gear.side_1.body);
               body_rotational_dependencies(gear.side_2.body)]
    owner = marker_owner(gear.contact_marker)
    isnothing(owner) || append!(columns, body_rotational_dependencies(owner))
    sort!(unique!(columns))
end

function wrapped_angle_difference(angle, physical)
    sine_angle, cosine_angle = sincos(angle)
    sine_physical, cosine_physical = sincos(physical)
    atan(sine_angle * cosine_physical - cosine_angle * sine_physical,
        cosine_angle * cosine_physical + sine_angle * sine_physical)
end

function gear_phase_definition_vector(gear, z)
    [wrapped_angle_difference(z[gear.position_variables[1]],
         physical_gear_phase(gear, gear.side_1, z)),
     wrapped_angle_difference(z[gear.position_variables[2]],
         physical_gear_phase(gear, gear.side_2, z)),
     z[gear.velocity_variables[1]] -
         physical_gear_velocity(gear, gear.side_1, z),
     z[gear.velocity_variables[2]] -
         physical_gear_velocity(gear, gear.side_2, z),
     z[gear.acceleration_variables[1]] -
         physical_gear_acceleration(gear, gear.side_1, z),
     z[gear.acceleration_variables[2]] -
         physical_gear_acceleration(gear, gear.side_2, z)]
end

function executable_blocks(gear::SpatialGearPair)
    definition_rows = [collect(gear.position_definition_equations);
        collect(gear.velocity_definition_equations);
        collect(gear.acceleration_definition_equations)]
    definition_variables = [collect(gear.position_variables);
        collect(gear.velocity_variables); collect(gear.acceleration_variables)]
    dependencies = sort!(unique!([phase_dependencies(gear);
        definition_variables]))
    definitions! = function (equations, t, z, zdot)
        equations[definition_rows] .= gear_phase_definition_vector(gear, z)
    end
    definitions_jacobian! = function (jacobian, t, z, zdot, coefficient)
        jacobian[definition_rows, dependencies] .+=
            spatial_local_state_jacobian(
            state -> gear_phase_definition_vector(gear, state), z,
            dependencies)
    end

    coupling_rows = [gear.position_equation, gear.velocity_equation,
        gear.acceleration_equation]
    coupling_dependencies = sort!(unique!([phase_dependencies(gear);
        collect(gear.position_variables)]))
    coupling! = function (equations, t, z, zdot)
        equations[gear.position_equation] = spatial_gear_position(gear, z)
        equations[gear.velocity_equation] = spatial_gear_velocity(gear, z)
        equations[gear.acceleration_equation] =
            spatial_gear_acceleration(gear, z)
    end
    coupling_vector = state -> [spatial_gear_position(gear, state),
        spatial_gear_velocity(gear, state),
        spatial_gear_acceleration(gear, state)]
    coupling_jacobian! = function (jacobian, t, z, zdot, coefficient)
        jacobian[coupling_rows, coupling_dependencies] .+=
            spatial_local_state_jacobian(coupling_vector, z,
                coupling_dependencies)
    end
    [ExecutableEquationBlock(gear.name, :phase_definitions,
         definition_rows, definitions!, definitions_jacobian!),
     ExecutableEquationBlock(gear.name, :gear_constraint,
         coupling_rows, coupling!, coupling_jacobian!)]
end

function initialize_spatial_gear_coordinates!(z, gear::SpatialGearPair,
        level = 2)
    if level >= 0
        z[gear.position_variables[1]] =
            physical_gear_phase(gear, gear.side_1, z)
        z[gear.position_variables[2]] =
            physical_gear_phase(gear, gear.side_2, z)
    end
    if level >= 1
        z[gear.velocity_variables[1]] =
            physical_gear_velocity(gear, gear.side_1, z)
        z[gear.velocity_variables[2]] =
            physical_gear_velocity(gear, gear.side_2, z)
    end
    if level >= 2
        z[gear.acceleration_variables[1]] =
            physical_gear_acceleration(gear, gear.side_1, z)
        z[gear.acceleration_variables[2]] =
            physical_gear_acceleration(gear, gear.side_2, z)
    end
    z
end

function body_force_contribution(body, z, point, force)
    parameters = @view z[body.euler_parameter_variables]
    orientation = rotation_matrix(parameters)
    lever_body = transpose(orientation) *
        (point - z[body.position_variables])
    force_body = transpose(orientation) * force
    [-force; -cross(lever_body, force_body)]
end

function gear_reaction_vector(gear, z)
    point = spatial_marker_position(gear.contact_marker, z)
    force = z[gear.reaction_variable] .*
        spatial_gear_contact_direction(gear, z)
    [body_force_contribution(gear.side_1.body, z, point, force);
     body_force_contribution(gear.side_2.body, z, point, -force)]
end

function reaction_dependencies(gear)
    columns = Int[gear.reaction_variable]
    for side in (gear.side_1, gear.side_2)
        append!(columns, side.body.position_variables)
        append!(columns, side.body.euler_parameter_variables)
    end
    owner = marker_owner(gear.contact_marker)
    if !isnothing(owner)
        append!(columns, owner.position_variables)
        append!(columns, owner.euler_parameter_variables)
    end
    sort!(unique!(columns))
end

function equation_contributions(gear::SpatialGearPair)
    rows = [collect(gear.side_1.body.balance_equations);
            collect(gear.side_2.body.balance_equations)]
    columns = reaction_dependencies(gear)
    residual! = function (equations, t, z, zdot)
        equations[rows] .+= gear_reaction_vector(gear, z)
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        jacobian[rows, columns] .+= spatial_local_state_jacobian(
            state -> gear_reaction_vector(gear, state), z, columns)
    end
    [EquationContribution(gear.name, :reaction_to_bodies, rows,
        residual!, jacobian!)]
end

function allocated_spatial_gear_pair(layout, name, joint_1, joint_2,
        contact_marker, contact_1, contact_2, initial;
        tangent_tolerance = 1.0e-5, phase = "initial")
    joints = (joint_1, joint_2)
    all(joint -> joint isa SpatialRevoluteJoint, joints) || throw(ArgumentError(
        "gear pair '$name' joints must be spatial revolute joints"))
    all(joint -> joint.marker_a isa SpatialBodyMarker, joints) ||
        throw(ArgumentError(
            "gear pair '$name' must list each gear marker first in its revolute joint"))
    bodies = (joint_1.marker_a.body, joint_2.marker_a.body)
    bodies[1] !== bodies[2] || throw(ArgumentError(
        "gear pair '$name' requires two different gear bodies"))
    carrier = marker_owner(contact_marker)
    any(marker_owner(joint.marker_b) === carrier for joint in joints) ||
        throw(ArgumentError(
            "gear pair '$name' contact_marker must be fixed to the base body of at least one referenced revolute joint"))

    contact_position = spatial_marker_position(contact_marker, initial)
    contact_orientation = spatial_marker_orientation(contact_marker, initial)
    axes = [marker_axis_kinematics(joint.marker_b, initial, 3).direction
        for joint in joints]
    centers = [spatial_marker_position(joint.marker_b, initial)
        for joint in joints]
    radial_vectors = [contact_position - center for center in centers]
    tangents = [cross(axis, radial) for (axis, radial) in
        zip(axes, radial_vectors)]
    tangent_norms = norm.(tangents)
    for index in 1:2
        scale = max(norm(radial_vectors[index]), 1.0)
        tangent_norms[index] > 1.0e-12 * scale || throw(ArgumentError(
            "gear pair '$name' contact point lies on the axis of joint '$(joints[index].name)'"))
    end
    unit_tangents = tangents ./ tangent_norms
    alignment = clamp(abs(dot(unit_tangents[1], unit_tangents[2])), 0.0, 1.0)
    mismatch = acos(alignment)
    if mismatch > tangent_tolerance
        degrees = round(rad2deg(mismatch); sigdigits = 6)
        throw(ArgumentError(
            "gear pair '$name' pitch tangents are not collinear " *
            "(angular mismatch $degrees degrees)"))
    end
    direction = unit_tangents[1]
    radii = [dot(axis, cross(radial, direction))
        for (axis, radial) in zip(axes, radial_vectors)]
    all(radius -> abs(radius) > 1.0e-12, radii) || throw(ArgumentError(
        "gear pair '$name' has a zero effective pitch radius"))

    sides = [SpatialGearSide(joints[index], bodies[index],
        transpose(contact_orientation) * axes[index],
        transpose(contact_orientation) *
            (radial_vectors[index] -
             dot(radial_vectors[index], axes[index]) .* axes[index]) /
            norm(radial_vectors[index] -
                 dot(radial_vectors[index], axes[index]) .* axes[index]),
        radii[index]) for index in 1:2]
    variables = component_variable_indices(layout, name)
    gear = SpatialGearPair(name, sides[1], sides[2], contact_marker,
        contact_1, contact_2,
        transpose(contact_orientation) * direction, 0.0, variables[1],
        variables[2:3], variables[4:5], variables[6:7],
        component_equation_indices(layout, name, :phase_acceleration),
        component_equation_indices(layout, name, :phase_velocity),
        component_equation_indices(layout, name, :phase_position),
        only(component_equation_indices(layout, name, :acceleration)),
        only(component_equation_indices(layout, name, :velocity)),
        only(component_equation_indices(layout, name, :position)))
    initialize_spatial_gear_coordinates!(initial, gear, 1)
    initial_phase = radii[1] * initial[gear.position_variables[1]] -
        radii[2] * initial[gear.position_variables[2]]
    phase_value = if phase isa Number
        Float64(phase)
    elseif phase == "initial"
        initial_phase
    else
        throw(ArgumentError(
            "gear pair '$name' phase must be numeric or 'initial'"))
    end
    isfinite(phase_value) || throw(ArgumentError(
        "gear pair '$name' phase must be finite"))
    SpatialGearPair(gear.name, gear.side_1, gear.side_2,
        gear.contact_marker, gear.contact_1, gear.contact_2,
        gear.direction_carrier, phase_value, gear.reaction_variable,
        gear.acceleration_variables, gear.velocity_variables,
        gear.position_variables, gear.acceleration_definition_equations,
        gear.velocity_definition_equations, gear.position_definition_equations,
        gear.acceleration_equation, gear.velocity_equation,
        gear.position_equation)
end

end
