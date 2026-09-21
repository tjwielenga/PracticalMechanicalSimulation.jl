"""Compliant scalar friction on planar joint primitives."""
module PlanarFrictionForces

using ForwardDiff
using LinearAlgebra
using ..AutomaticAnalysis
using ..FrictionLaws
using ..PlanarAppliedForces: PlanarBodyOrientationMarker,
    PlanarGroundOrientationMarker, PlanarBodyPointMarker,
    marker_angle, marker_angular_velocity, add_marker_torque!,
    add_marker_torque_jacobian!, add_body_point_force!
import ..PlanarAppliedForces
using ..PlanarDirectedDistances
using ..PlanarComponentAssembly

import ..PlanarComponentAssembly: component_registration,
    executable_blocks, equation_contributions

export PlanarRevoluteFriction, planar_revolute_friction_registration,
    allocated_planar_revolute_friction,
    initialize_planar_revolute_friction!,
    set_planar_revolute_friction_stage!, revolute_friction_rate,
    revolute_bearing_load,
    PlanarTranslationalFriction,
    planar_translational_friction_registration,
    allocated_planar_translational_friction,
    initialize_planar_translational_friction!,
    set_planar_translational_friction_stage!,
    translational_friction_rate, translational_guide_load,
    PlanarInplaneFriction, planar_inplane_friction_registration,
    allocated_planar_inplane_friction,
    initialize_planar_inplane_friction!,
    set_planar_inplane_friction_stage!, inplane_friction_rate,
    inplane_normal_load

function local_state_jacobian(function_value, z, columns)
    isempty(columns) && return zeros(eltype(z), length(function_value(z)), 0)
    inputs = collect(z[columns])
    ForwardDiff.jacobian(inputs) do local_values
        state = Vector{eltype(local_values)}(undef, length(z))
        state .= z
        state[columns] .= local_values
        function_value(state)
    end
end

smoothed_magnitude(value) = sqrt(value^2 + 1.0e-18) - 1.0e-9

mutable struct PlanarRevoluteFriction{J}
    name::Symbol
    joint::J
    stiffness::Float64
    damping::Float64
    effective_radius::Float64
    preload::Float64
    static_coefficient::Float64
    dynamic_coefficient::Float64
    transition_speed::Float64
    release_time::Float64
    stage::Symbol
    anchor::Float64
    shear_variable::Int
    slip_variable::Int
    bearing_load_variable::Int
    torque_variable::Int
    friction_equations::UnitRange{Int}
end

function planar_revolute_friction_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:shear, :user_state_steady, 0),
        VariableDeclaration(:slip, :applied_rate, 1),
        VariableDeclaration(:bearing_load, :applied_load, 2),
        VariableDeclaration(:torque, :applied_load, 2)]
    equations = EquationDeclaration[
        EquationDeclaration(:shear, :user_differential_steady, 0, :friction),
        EquationDeclaration(:slip, :applied_definition, 1, :friction),
        EquationDeclaration(:bearing_load, :applied_definition, 2, :friction),
        EquationDeclaration(:torque, :applied_definition, 2, :friction)]
    ComponentRegistration(name, variables,
        [EquationBlockDeclaration(:friction, equations)])
end

component_registration(friction::PlanarRevoluteFriction) =
    planar_revolute_friction_registration(friction.name)

function allocated_planar_revolute_friction(layout, name, joint,
        stiffness, damping, effective_radius, preload, static_coefficient,
        dynamic_coefficient, transition_speed, release_time)
    variables = component_variable_indices(layout, name)
    PlanarRevoluteFriction(name, joint, stiffness, damping, effective_radius,
        preload, static_coefficient, dynamic_coefficient, transition_speed,
        release_time, :dynamic, 0.0, variables[1], variables[2],
        variables[3], variables[4],
        component_equation_indices(layout, name, :friction))
end

function revolute_angle(friction, z)
    joint = friction.joint
    difference = marker_angle(joint.rotation_marker_a, z) -
        marker_angle(joint.rotation_marker_b, z)
    atan(sin(difference), cos(difference))
end

function revolute_slip(friction, z)
    joint = friction.joint
    marker_angular_velocity(joint.rotation_marker_a, z) -
        marker_angular_velocity(joint.rotation_marker_b, z)
end

function revolute_bearing_load(friction, z)
    reaction = z[friction.joint.reaction_variables]
    sqrt(sum(abs2, reaction) + 1.0e-18) - 1.0e-9 + friction.preload
end

function revolute_friction_capacity(friction, z, slip_squared)
    scalar_friction_coefficient(friction, slip_squared) *
        friction.effective_radius *
        max(zero(eltype(z)), z[friction.bearing_load_variable])
end

function revolute_friction_rate(friction, z)
    shear = z[friction.shear_variable]
    slip = revolute_slip(friction, z)
    capacity = revolute_friction_capacity(friction, z, slip^2)
    scalar_bristle_rate(friction, shear, slip, capacity)
end

function calculated_revolute_friction_torque(friction, z)
    shear = z[friction.shear_variable]
    slip = revolute_slip(friction, z)
    capacity = revolute_friction_capacity(friction, z, slip^2)
    scalar_bristle_force(friction, shear, slip, capacity)
end

function initialize_planar_revolute_friction!(state, friction;
        reset_anchor = false)
    reset_anchor &&
        (friction.anchor = revolute_angle(friction, state) -
            state[friction.shear_variable])
    state[friction.slip_variable] = revolute_slip(friction, state)
    state[friction.bearing_load_variable] = revolute_bearing_load(friction, state)
    state[friction.torque_variable] =
        calculated_revolute_friction_torque(friction, state)
    state
end

function set_planar_revolute_friction_stage!(friction, stage)
    friction.stage = stage
    friction
end

function revolute_friction_residual(friction, z, zdot)
    shear = z[friction.shear_variable]
    shear_equation = friction.stage == :static ?
        shear - (revolute_angle(friction, z) - friction.anchor) :
        zdot[friction.shear_variable] - revolute_friction_rate(friction, z)
    [shear_equation;
     z[friction.slip_variable] - revolute_slip(friction, z);
     z[friction.bearing_load_variable] - revolute_bearing_load(friction, z);
     z[friction.torque_variable] -
        calculated_revolute_friction_torque(friction, z)]
end

function orientation_dependencies(marker::PlanarBodyOrientationMarker)
    [marker.theta_variable, marker.omega_variable]
end
orientation_dependencies(::PlanarGroundOrientationMarker) = Int[]

function executable_blocks(friction::PlanarRevoluteFriction)
    joint = friction.joint
    columns = sort!(unique!([
        orientation_dependencies(joint.rotation_marker_a);
        orientation_dependencies(joint.rotation_marker_b);
        collect(joint.reaction_variables);
        friction.shear_variable;
        friction.slip_variable;
        friction.bearing_load_variable;
        friction.torque_variable]))
    residual! = (equations, t, z, zdot) ->
        (equations[friction.friction_equations] .=
            revolute_friction_residual(friction, z, zdot))
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        partials = local_state_jacobian(
            state -> revolute_friction_residual(friction, state, zdot),
            z, columns)
        jacobian[friction.friction_equations, columns] .+= partials
        friction.stage == :static ||
            (jacobian[first(friction.friction_equations),
                friction.shear_variable] += coefficient)
    end
    [ExecutableEquationBlock(friction.name, :friction,
        collect(friction.friction_equations), residual!, jacobian!)]
end

function equation_contributions(friction::PlanarRevoluteFriction)
    joint = friction.joint
    marker_a = joint.rotation_marker_a
    marker_b = joint.rotation_marker_b
    rows = Int[]
    marker_a isa PlanarBodyOrientationMarker &&
        push!(rows, marker_a.torque_equation)
    marker_b isa PlanarBodyOrientationMarker &&
        push!(rows, marker_b.torque_equation)
    residual! = function (equations, t, z, zdot)
        torque = z[friction.torque_variable]
        add_marker_torque!(equations, marker_a, torque)
        add_marker_torque!(equations, marker_b, -torque)
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        add_marker_torque_jacobian!(jacobian, marker_a,
            friction.torque_variable, 1)
        add_marker_torque_jacobian!(jacobian, marker_b,
            friction.torque_variable, -1)
    end
    [EquationContribution(friction.name, :torque_to_bodies, unique(rows),
        residual!, jacobian!)]
end

abstract type AbstractPlanarTangentialFriction end

mutable struct PlanarTranslationalFriction{J} <:
        AbstractPlanarTangentialFriction
    name::Symbol
    joint::J
    stiffness::Float64
    damping::Float64
    preload::Float64
    static_coefficient::Float64
    dynamic_coefficient::Float64
    transition_speed::Float64
    release_time::Float64
    stage::Symbol
    anchor::Float64
    shear_variable::Int
    slip_variable::Int
    load_variable::Int
    force_variable::Int
    global_force_variables::UnitRange{Int}
    friction_equations::UnitRange{Int}
end

mutable struct PlanarInplaneFriction{C} <: AbstractPlanarTangentialFriction
    name::Symbol
    constraint::C
    stiffness::Float64
    damping::Float64
    preload::Float64
    static_coefficient::Float64
    dynamic_coefficient::Float64
    transition_speed::Float64
    release_time::Float64
    stage::Symbol
    anchor::Float64
    shear_variable::Int
    slip_variable::Int
    load_variable::Int
    force_variable::Int
    global_force_variables::UnitRange{Int}
    friction_equations::UnitRange{Int}
end

function tangential_friction_registration(name::Symbol, load_name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:shear, :user_state_steady, 0),
        VariableDeclaration(:slip, :applied_rate, 1),
        VariableDeclaration(load_name, :applied_load, 2),
        VariableDeclaration(:force, :applied_load, 2),
        VariableDeclaration(:F_x, :applied_load, 2),
        VariableDeclaration(:F_y, :applied_load, 2)]
    equations = EquationDeclaration[
        EquationDeclaration(:shear, :user_differential_steady, 0, :friction),
        EquationDeclaration(:slip, :applied_definition, 1, :friction),
        EquationDeclaration(load_name, :applied_definition, 2, :friction),
        EquationDeclaration(:force, :applied_definition, 2, :friction),
        EquationDeclaration(:global_force_x, :applied_definition, 2, :friction),
        EquationDeclaration(:global_force_y, :applied_definition, 2, :friction)]
    ComponentRegistration(name, variables,
        [EquationBlockDeclaration(:friction, equations)])
end

planar_translational_friction_registration(name::Symbol) =
    tangential_friction_registration(name, :guide_load)
planar_inplane_friction_registration(name::Symbol) =
    tangential_friction_registration(name, :normal_load)

component_registration(friction::PlanarTranslationalFriction) =
    planar_translational_friction_registration(friction.name)
component_registration(friction::PlanarInplaneFriction) =
    planar_inplane_friction_registration(friction.name)

function allocated_planar_translational_friction(layout, name, joint,
        stiffness, damping, preload, static_coefficient, dynamic_coefficient,
        transition_speed, release_time)
    variables = component_variable_indices(layout, name)
    PlanarTranslationalFriction(name, joint, stiffness, damping, preload,
        static_coefficient, dynamic_coefficient, transition_speed,
        release_time, :dynamic, 0.0, variables[1], variables[2], variables[3],
        variables[4], variables[5:6],
        component_equation_indices(layout, name, :friction))
end

function allocated_planar_inplane_friction(layout, name, constraint,
        stiffness, damping, preload, static_coefficient, dynamic_coefficient,
        transition_speed, release_time)
    variables = component_variable_indices(layout, name)
    PlanarInplaneFriction(name, constraint, stiffness, damping, preload,
        static_coefficient, dynamic_coefficient, transition_speed,
        release_time, :dynamic, 0.0, variables[1], variables[2], variables[3],
        variables[4], variables[5:6],
        component_equation_indices(layout, name, :friction))
end

friction_constraint(friction::PlanarTranslationalFriction) =
    friction.joint.inplane
friction_constraint(friction::PlanarInplaneFriction) = friction.constraint

function tangent_kinematics(friction, z)
    constraint = friction_constraint(friction)
    values = directed_distance_values(constraint.geometry, z;
        include_acceleration = false)
    tangent = .-values.axis.transverse
    tangent_rate = values.axis.unit .* values.axis.omega
    position = dot(values.separation, tangent)
    slip = dot(values.relative_velocity, tangent) +
        dot(values.separation, tangent_rate)
    (; values, tangent, position, slip)
end

function tangential_normal_load(friction, z)
    reaction = z[friction_constraint(friction).reaction_variable]
    smoothed_magnitude(reaction) + friction.preload
end

translational_guide_load(friction::PlanarTranslationalFriction, z) =
    tangential_normal_load(friction, z)
inplane_normal_load(friction::PlanarInplaneFriction, z) =
    tangential_normal_load(friction, z)

function tangential_friction_capacity(friction, z, slip_squared)
    scalar_friction_coefficient(friction, slip_squared) *
        max(zero(eltype(z)), z[friction.load_variable])
end

function tangential_friction_rate(friction, z)
    shear = z[friction.shear_variable]
    slip = tangent_kinematics(friction, z).slip
    capacity = tangential_friction_capacity(friction, z, slip^2)
    scalar_bristle_rate(friction, shear, slip, capacity)
end

translational_friction_rate(friction::PlanarTranslationalFriction, z) =
    tangential_friction_rate(friction, z)
inplane_friction_rate(friction::PlanarInplaneFriction, z) =
    tangential_friction_rate(friction, z)

function calculated_tangential_friction_force(friction, z)
    shear = z[friction.shear_variable]
    slip = tangent_kinematics(friction, z).slip
    capacity = tangential_friction_capacity(friction, z, slip^2)
    scalar_bristle_force(friction, shear, slip, capacity)
end

function initialize_planar_tangential_friction!(state, friction;
        reset_anchor = false)
    kinematics = tangent_kinematics(friction, state)
    reset_anchor &&
        (friction.anchor = kinematics.position - state[friction.shear_variable])
    state[friction.slip_variable] = kinematics.slip
    state[friction.load_variable] = tangential_normal_load(friction, state)
    state[friction.force_variable] =
        calculated_tangential_friction_force(friction, state)
    state[friction.global_force_variables] .=
        state[friction.force_variable] .* kinematics.tangent
    state
end

initialize_planar_translational_friction!(state, friction;
        reset_anchor = false) =
    initialize_planar_tangential_friction!(state, friction; reset_anchor)
initialize_planar_inplane_friction!(state, friction;
        reset_anchor = false) =
    initialize_planar_tangential_friction!(state, friction; reset_anchor)

function set_planar_tangential_friction_stage!(friction, stage)
    friction.stage = stage
    friction
end
set_planar_translational_friction_stage!(friction, stage) =
    set_planar_tangential_friction_stage!(friction, stage)
set_planar_inplane_friction_stage!(friction, stage) =
    set_planar_tangential_friction_stage!(friction, stage)

function tangential_friction_residual(friction, z, zdot)
    kinematics = tangent_kinematics(friction, z)
    shear = z[friction.shear_variable]
    shear_equation = friction.stage == :static ?
        shear - (kinematics.position - friction.anchor) :
        zdot[friction.shear_variable] - tangential_friction_rate(friction, z)
    force = calculated_tangential_friction_force(friction, z)
    [shear_equation;
     z[friction.slip_variable] - kinematics.slip;
     z[friction.load_variable] - tangential_normal_load(friction, z);
     z[friction.force_variable] - force;
     z[friction.global_force_variables] .-
        z[friction.force_variable] .* kinematics.tangent]
end

function executable_blocks(friction::AbstractPlanarTangentialFriction)
    constraint = friction_constraint(friction)
    columns = sort!(unique!([
        directed_distance_dependencies(constraint.geometry);
        constraint.reaction_variable;
        friction.shear_variable;
        friction.slip_variable;
        friction.load_variable;
        friction.force_variable;
        collect(friction.global_force_variables)]))
    residual! = (equations, t, z, zdot) ->
        (equations[friction.friction_equations] .=
            tangential_friction_residual(friction, z, zdot))
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        partials = local_state_jacobian(
            state -> tangential_friction_residual(friction, state, zdot),
            z, columns)
        jacobian[friction.friction_equations, columns] .+= partials
        friction.stage == :static ||
            (jacobian[first(friction.friction_equations),
                friction.shear_variable] += coefficient)
    end
    [ExecutableEquationBlock(friction.name, :friction,
        collect(friction.friction_equations), residual!, jacobian!)]
end

function equation_contributions(friction::AbstractPlanarTangentialFriction)
    constraint = friction_constraint(friction)
    first_marker = constraint.geometry.marker_i
    second_marker = constraint.geometry.marker_j
    rows = Int[]
    for marker in (first_marker, second_marker)
        marker isa PlanarBodyPointMarker || continue
        append!(rows, marker.force_equations)
        push!(rows, marker.torque_equation)
    end
    residual! = function (equations, t, z, zdot)
        kinematics = tangent_kinematics(friction, z).values
        force = z[friction.global_force_variables]
        add_body_point_force!(equations, first_marker,
            kinematics.marker_i, force)
        add_body_point_force!(equations, second_marker,
            kinematics.marker_j, -force)
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        kinematics = tangent_kinematics(friction, z).values
        force = z[friction.global_force_variables]
        PlanarAppliedForces.add_point_force_partials!(jacobian,
            first_marker, 1, kinematics.marker_i,
            friction.global_force_variables, force)
        PlanarAppliedForces.add_point_force_partials!(jacobian,
            second_marker, -1, kinematics.marker_j,
            friction.global_force_variables, force)
    end
    [EquationContribution(friction.name, :force_to_bodies, unique(rows),
        residual!, jacobian!)]
end

end
