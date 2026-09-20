"""Compliant friction on spatial contacts and joint primitives."""
module SpatialFrictionForces

using ForwardDiff
using LinearAlgebra
using ..AutomaticAnalysis
using ..SpatialComponentAssembly
using ..SpatialModeling
using ..SpatialConstraints
using ..SpatialAppliedForces
using ..SpatialPlaneContacts

import ..SpatialComponentAssembly: component_registration,
    executable_blocks, equation_contributions

export SpatialSurfaceFriction, spatial_surface_friction_registration,
    allocated_spatial_surface_friction, initialize_spatial_surface_friction!,
    set_spatial_surface_friction_stage!, surface_friction_rates,
    SpatialRevoluteFriction, spatial_revolute_friction_registration,
    allocated_spatial_revolute_friction,
    initialize_spatial_revolute_friction!,
    set_spatial_revolute_friction_stage!, revolute_friction_rate,
    SpatialTranslationalFriction,
    spatial_translational_friction_registration,
    allocated_spatial_translational_friction,
    initialize_spatial_translational_friction!,
    set_spatial_translational_friction_stage!, translational_friction_rate,
    SpatialInplaneFriction, spatial_inplane_friction_registration,
    allocated_spatial_inplane_friction,
    initialize_spatial_inplane_friction!,
    set_spatial_inplane_friction_stage!, inplane_friction_rates

mutable struct SpatialSurfaceFriction{C}
    name::Symbol
    contact::C
    stiffness::Float64
    damping::Float64
    static_coefficient::Float64
    dynamic_coefficient::Float64
    transition_speed::Float64
    release_time::Float64
    stage::Symbol
    anchor::Vector{Float64}
    shear_variables::UnitRange{Int}
    velocity_variables::UnitRange{Int}
    force_variables::UnitRange{Int}
    global_force_variables::UnitRange{Int}
    friction_equations::UnitRange{Int}
end

function spatial_surface_friction_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:shear_x, :user_state_steady, 0),
        VariableDeclaration(:shear_y, :user_state_steady, 0),
        VariableDeclaration(:slip_x, :applied_rate, 1),
        VariableDeclaration(:slip_y, :applied_rate, 1),
        VariableDeclaration(:friction_x, :applied_load, 2),
        VariableDeclaration(:friction_y, :applied_load, 2),
        [VariableDeclaration(Symbol(:F_, axis), :applied_load, 2)
            for axis in (:x, :y, :z)]...]
    equations = EquationDeclaration[
        EquationDeclaration(:shear_x, :user_differential_steady, 0, :friction),
        EquationDeclaration(:shear_y, :user_differential_steady, 0, :friction),
        EquationDeclaration(:slip_x, :applied_definition, 1, :friction),
        EquationDeclaration(:slip_y, :applied_definition, 1, :friction),
        EquationDeclaration(:friction_x, :applied_definition, 2, :friction),
        EquationDeclaration(:friction_y, :applied_definition, 2, :friction),
        [EquationDeclaration(Symbol(:global_force_, axis),
            :applied_definition, 2, :friction)
            for axis in (:x, :y, :z)]...]
    ComponentRegistration(name, variables,
        [EquationBlockDeclaration(:friction, equations)])
end

component_registration(friction::SpatialSurfaceFriction) =
    spatial_surface_friction_registration(friction.name)

function allocated_spatial_surface_friction(layout, name, contact,
        stiffness, damping, static_coefficient, dynamic_coefficient,
        transition_speed, release_time)
    variables = component_variable_indices(layout, name)
    SpatialSurfaceFriction(name, contact, stiffness, damping, static_coefficient,
        dynamic_coefficient, transition_speed, release_time, :dynamic,
        zeros(2), variables[1:2], variables[3:4], variables[5:6],
        variables[7:9], component_equation_indices(layout, name, :friction))
end

function tangent_axes(friction, z)
    orientation = spatial_marker_orientation(friction.contact.plane_marker, z)
    orientation[:, 1:2]
end

function tangent_position(friction, z)
    contact = friction.contact
    axes = tangent_axes(friction, z)
    transpose(axes) * (spatial_marker_position(contact.sphere_marker, z) -
        spatial_marker_position(contact.plane_marker, z))
end

function contact_slip(friction, z)
    contact = friction.contact
    values = SpatialPlaneContacts.spatial_plane_contact_kinematics(contact, z)
    first = contact.sphere_marker
    second = contact.plane_marker
    velocity = spatial_marker_velocity(first, z) -
        spatial_marker_velocity(second, z)
    # The sphere's surface moves relative to its center when it rotates.
    if first isa SpatialBodyMarker
        body = first.body
        omega = spatial_marker_orientation(first, z) *
            (transpose(first.orientation_body) * z[body.angular_velocity_variables])
        velocity -= cross(omega, contact.radius .* values.normal)
    end
    if second isa SpatialBodyMarker
        body = second.body
        orientation = rotation_matrix(z[body.euler_parameter_variables])
        omega = orientation * z[body.angular_velocity_variables]
        velocity -= cross(omega, values.contact_point -
            spatial_marker_position(second, z))
    end
    transpose(tangent_axes(friction, z)) * velocity
end

function friction_capacity(friction, z, speed_squared)
    normal_force = max(zero(eltype(z)),
        z[friction.contact.normal_force_variable])
    scalar_friction_coefficient(friction, speed_squared) * normal_force
end

function surface_friction_rates(friction, z)
    shear = z[friction.shear_variables]
    velocity = contact_slip(friction, z)
    capacity = friction_capacity(friction, z, sum(abs2, velocity))
    vector_bristle_rate(friction, shear, velocity, capacity)
end

function calculated_friction_force(friction, z)
    shear = z[friction.shear_variables]
    normal_force = z[friction.contact.normal_force_variable]
    normal_force > 0 || return zero.(shear)
    velocity = contact_slip(friction, z)
    capacity = friction_capacity(friction, z, sum(abs2, velocity))
    vector_bristle_force(friction, shear, velocity, capacity)
end

function initialize_spatial_surface_friction!(state, friction;
        reset_anchor = false)
    reset_anchor && (friction.anchor .= tangent_position(friction, state) .-
        state[friction.shear_variables])
    state[friction.velocity_variables] .= contact_slip(friction, state)
    state[friction.force_variables] .= calculated_friction_force(friction, state)
    state[friction.global_force_variables] .= tangent_axes(friction, state) *
        state[friction.force_variables]
    state
end

function set_spatial_surface_friction_stage!(friction, stage)
    friction.stage = stage
    friction
end

function friction_residual(friction, z, zdot)
    shear = z[friction.shear_variables]
    rates = friction.stage == :static ?
        shear .- (tangent_position(friction, z) .- friction.anchor) :
        zdot[friction.shear_variables] .- surface_friction_rates(friction, z)
    [rates;
     z[friction.velocity_variables] .- contact_slip(friction, z);
     z[friction.force_variables] .- calculated_friction_force(friction, z);
     z[friction.global_force_variables] .-
        tangent_axes(friction, z) * z[friction.force_variables]]
end

function executable_blocks(friction::SpatialSurfaceFriction)
    contact = friction.contact
    columns = sort!(unique!([
        SpatialPlaneContacts.contact_kinematic_dependencies(contact.sphere_marker);
        SpatialPlaneContacts.contact_kinematic_dependencies(contact.plane_marker);
        contact.normal_force_variable;
        collect(friction.shear_variables);
        collect(friction.velocity_variables);
        collect(friction.force_variables);
        collect(friction.global_force_variables)]))
    residual! = (equations, t, z, zdot) ->
        (equations[friction.friction_equations] .=
            friction_residual(friction, z, zdot))
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        partials = SpatialPlaneContacts.local_state_jacobian(
            state -> friction_residual(friction, state, zdot), z, columns)
        jacobian[friction.friction_equations, columns] .+= partials
        if friction.stage != :static
            for (row, variable) in zip(friction.friction_equations[1:2],
                    friction.shear_variables)
                jacobian[row, variable] += coefficient
            end
        end
    end
    [ExecutableEquationBlock(friction.name, :friction,
        collect(friction.friction_equations), residual!, jacobian!)]
end

function equation_contributions(friction::SpatialSurfaceFriction)
    contact = friction.contact
    rows = Int[]
    first = contact.sphere_marker
    second = contact.plane_marker
    first isa SpatialBodyMarker && append!(rows, first.body.balance_equations)
    second isa SpatialBodyMarker && append!(rows, second.body.balance_equations)
    columns = sort!(unique!([
        SpatialPlaneContacts.contact_configuration_dependencies(first);
        SpatialPlaneContacts.contact_configuration_dependencies(second);
        collect(friction.global_force_variables)]))
    contribution = function (z)
        point = SpatialPlaneContacts.spatial_plane_contact_kinematics(
            contact, z).contact_point
        force = z[friction.global_force_variables]
        values = eltype(z)[]
        first isa SpatialBodyMarker && append!(values,
            SpatialPlaneContacts.body_force_contribution(first.body, z,
                point, force))
        second isa SpatialBodyMarker && append!(values,
            SpatialPlaneContacts.body_force_contribution(second.body, z,
                point, -force))
        values
    end
    residual! = (equations, t, z, zdot) ->
        (equations[rows] .+= contribution(z))
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        partials = SpatialPlaneContacts.local_state_jacobian(contribution,
            z, columns)
        jacobian[rows, columns] .+= partials
    end
    [EquationContribution(friction.name, :load_to_bodies, rows,
        residual!, jacobian!)]
end

"""Common two-direction bristle law for plane contact and inplane guides."""
function vector_bristle_rate(friction, shear, velocity, capacity)
    capacity <= 1.0e-12 && return -shear ./ friction.release_time
    speed_squared = sum(abs2, velocity)
    # Avoid differentiating sqrt at zero slip; the sticking branch is regular.
    speed = speed_squared > 0 ? sqrt(speed_squared) : zero(speed_squared)
    velocity .- (friction.stiffness * speed / capacity) .* shear
end

function vector_bristle_force(friction, shear, velocity, capacity)
    capacity <= 1.0e-12 && return zero.(shear)
    rate = friction.stage == :static ? zero.(shear) :
        vector_bristle_rate(friction, shear, velocity, capacity)
    trial = -friction.stiffness .* shear .- friction.damping .* rate
    magnitude = norm(trial)
    factor = magnitude > capacity ? capacity / magnitude : one(magnitude)
    factor .* trial
end

"""Common scalar bristle law; the joint supplies slip and normal capacity."""
scalar_friction_coefficient(friction, slip_squared) =
    friction.dynamic_coefficient +
    (friction.static_coefficient - friction.dynamic_coefficient) *
        exp(-slip_squared / friction.transition_speed^2)

function scalar_bristle_rate(friction, shear, slip, capacity)
    capacity <= 1.0e-12 && return -shear / friction.release_time
    slip - (friction.stiffness * abs(slip) / capacity) * shear
end

function scalar_bristle_force(friction, shear, slip, capacity)
    capacity <= 1.0e-12 && return zero(shear)
    rate = friction.stage == :static ? zero(shear) :
        scalar_bristle_rate(friction, shear, slip, capacity)
    trial = -friction.stiffness * shear - friction.damping * rate
    clamp(trial, -capacity, capacity)
end

"""
One angular bristle state opposing a revolute joint's free rotation.

The radial component of the spherical-joint reaction estimates the bearing
normal load. `preload` supplies an independent normal load when appropriate.
Axial thrust friction is not included in this first bearing model.
"""
mutable struct SpatialRevoluteFriction{J}
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
    anchor::Base.RefValue{Float64}
    shear_variable::Int
    slip_variable::Int
    bearing_load_variable::Int
    torque_variable::Int
    global_torque_variables::UnitRange{Int}
    friction_equations::UnitRange{Int}
end

function spatial_revolute_friction_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:shear, :user_state_steady, 0),
        VariableDeclaration(:slip, :applied_rate, 1),
        VariableDeclaration(:bearing_load, :applied_load, 2),
        VariableDeclaration(:torque, :applied_load, 2),
        [VariableDeclaration(Symbol(:T_, axis), :applied_load, 2)
            for axis in (:x, :y, :z)]...]
    equations = EquationDeclaration[
        EquationDeclaration(:shear, :user_differential_steady, 0, :friction),
        EquationDeclaration(:slip, :applied_definition, 1, :friction),
        EquationDeclaration(:bearing_load, :applied_definition, 2, :friction),
        EquationDeclaration(:torque, :applied_definition, 2, :friction),
        [EquationDeclaration(Symbol(:global_torque_, axis),
            :applied_definition, 2, :friction)
            for axis in (:x, :y, :z)]...]
    ComponentRegistration(name, variables,
        [EquationBlockDeclaration(:friction, equations)])
end

component_registration(friction::SpatialRevoluteFriction) =
    spatial_revolute_friction_registration(friction.name)

function allocated_spatial_revolute_friction(layout, name, joint,
        stiffness, damping, effective_radius, preload, static_coefficient,
        dynamic_coefficient, transition_speed, release_time)
    variables = component_variable_indices(layout, name)
    SpatialRevoluteFriction(name, joint, stiffness, damping, effective_radius,
        preload, static_coefficient, dynamic_coefficient, transition_speed,
        release_time, :dynamic, Ref(0.0), variables[1], variables[2],
        variables[3], variables[4], variables[5:7],
        component_equation_indices(layout, name, :friction))
end

function revolute_friction_axis(friction, z)
    SpatialConstraints.marker_axis_kinematics(
        friction.joint.hinge.marker_j, z, 3).direction
end

function revolute_bearing_load(friction, z)
    reaction = z[friction.joint.spherical.reaction_variables]
    axis = revolute_friction_axis(friction, z)
    radial = reaction .- dot(reaction, axis) .* axis
    # The small regularization makes the normal-load Jacobian defined when
    # a bearing is unloaded; it is negligible at engineering force scales.
    regularization = 1.0e-9
    sqrt(sum(abs2, radial) + regularization^2) - regularization +
        friction.preload
end

function revolute_friction_capacity(friction, z, slip_squared)
    scalar_friction_coefficient(friction, slip_squared) *
        friction.effective_radius *
        max(zero(eltype(z)), z[friction.bearing_load_variable])
end

function revolute_friction_rate(friction, z)
    shear = z[friction.shear_variable]
    slip = SpatialConstraints.hinge_angular_velocity(friction.joint.hinge, z)
    capacity = revolute_friction_capacity(friction, z, slip^2)
    scalar_bristle_rate(friction, shear, slip, capacity)
end

function calculated_revolute_friction_torque(friction, z)
    slip = SpatialConstraints.hinge_angular_velocity(friction.joint.hinge, z)
    capacity = revolute_friction_capacity(friction, z, slip^2)
    shear = z[friction.shear_variable]
    scalar_bristle_force(friction, shear, slip, capacity)
end

function initialize_spatial_revolute_friction!(state, friction;
        reset_anchor = false)
    reset_anchor && (friction.anchor[] =
        SpatialConstraints.hinge_angle(friction.joint.hinge, state) -
        state[friction.shear_variable])
    state[friction.slip_variable] =
        SpatialConstraints.hinge_angular_velocity(friction.joint.hinge, state)
    state[friction.bearing_load_variable] = revolute_bearing_load(friction, state)
    state[friction.torque_variable] =
        calculated_revolute_friction_torque(friction, state)
    state[friction.global_torque_variables] .=
        state[friction.torque_variable] .* revolute_friction_axis(friction, state)
    state
end

function set_spatial_revolute_friction_stage!(friction, stage)
    friction.stage = stage
    friction
end

function revolute_friction_residual(friction, z, zdot)
    shear = z[friction.shear_variable]
    shear_equation = friction.stage == :static ?
        shear - (SpatialConstraints.hinge_angle(friction.joint.hinge, z) -
            friction.anchor[]) :
        zdot[friction.shear_variable] - revolute_friction_rate(friction, z)
    slip = SpatialConstraints.hinge_angular_velocity(friction.joint.hinge, z)
    bearing_load = revolute_bearing_load(friction, z)
    torque = calculated_revolute_friction_torque(friction, z)
    [shear_equation;
     z[friction.slip_variable] - slip;
     z[friction.bearing_load_variable] - bearing_load;
     z[friction.torque_variable] - torque;
     z[friction.global_torque_variables] .-
        z[friction.torque_variable] .* revolute_friction_axis(friction, z)]
end

function executable_blocks(friction::SpatialRevoluteFriction)
    joint = friction.joint
    columns = sort!(unique!([
        SpatialPlaneContacts.contact_kinematic_dependencies(joint.marker_a);
        SpatialPlaneContacts.contact_kinematic_dependencies(joint.marker_b);
        collect(joint.spherical.reaction_variables);
        friction.shear_variable;
        friction.slip_variable;
        friction.bearing_load_variable;
        friction.torque_variable;
        collect(friction.global_torque_variables)]))
    residual! = (equations, t, z, zdot) ->
        (equations[friction.friction_equations] .=
            revolute_friction_residual(friction, z, zdot))
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        partials = SpatialPlaneContacts.local_state_jacobian(
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

function equation_contributions(friction::SpatialRevoluteFriction)
    marker_i = friction.joint.hinge.marker_i
    marker_j = friction.joint.hinge.marker_j
    rows = unique!([SpatialConstraints.body_torque_rows(marker_i);
                    SpatialConstraints.body_torque_rows(marker_j)])
    residual! = function (equations, t, z, zdot)
        torque = z[friction.global_torque_variables]
        SpatialAppliedForces.add_applied_torque_to_body!(
            equations, z, marker_i, torque, 1)
        SpatialAppliedForces.add_applied_torque_to_body!(
            equations, z, marker_j, torque, -1)
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        torque = z[friction.global_torque_variables]
        SpatialAppliedForces.add_applied_torque_jacobian!(jacobian, z,
            marker_i, torque, 1, friction.global_torque_variables)
        SpatialAppliedForces.add_applied_torque_jacobian!(jacobian, z,
            marker_j, torque, -1, friction.global_torque_variables)
    end
    [EquationContribution(friction.name, :torque_to_bodies, rows,
        residual!, jacobian!)]
end

"""Axial bristle friction for an inline guide, with transverse guide load."""
mutable struct SpatialTranslationalFriction{J}
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
    anchor::Base.RefValue{Float64}
    shear_variable::Int
    slip_variable::Int
    guide_load_variable::Int
    force_variable::Int
    global_force_variables::UnitRange{Int}
    friction_equations::UnitRange{Int}
end

function spatial_translational_friction_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:shear, :user_state_steady, 0),
        VariableDeclaration(:slip, :applied_rate, 1),
        VariableDeclaration(:guide_load, :applied_load, 2),
        VariableDeclaration(:force, :applied_load, 2),
        [VariableDeclaration(Symbol(:F_, axis), :applied_load, 2)
            for axis in (:x, :y, :z)]...]
    equations = EquationDeclaration[
        EquationDeclaration(:shear, :user_differential_steady, 0, :friction),
        EquationDeclaration(:slip, :applied_definition, 1, :friction),
        EquationDeclaration(:guide_load, :applied_definition, 2, :friction),
        EquationDeclaration(:force, :applied_definition, 2, :friction),
        [EquationDeclaration(Symbol(:global_force_, axis),
            :applied_definition, 2, :friction)
            for axis in (:x, :y, :z)]...]
    ComponentRegistration(name, variables,
        [EquationBlockDeclaration(:friction, equations)])
end

component_registration(friction::SpatialTranslationalFriction) =
    spatial_translational_friction_registration(friction.name)

function allocated_spatial_translational_friction(layout, name, joint,
        stiffness, damping, preload, static_coefficient, dynamic_coefficient,
        transition_speed, release_time)
    variables = component_variable_indices(layout, name)
    SpatialTranslationalFriction(name, joint, stiffness, damping, preload,
        static_coefficient, dynamic_coefficient, transition_speed,
        release_time, :dynamic, Ref(0.0), variables[1], variables[2],
        variables[3], variables[4], variables[5:7],
        component_equation_indices(layout, name, :friction))
end

function translational_guide_load(friction, z)
    joint = friction.joint
    λx = z[joint.inplane_x.reaction_variable]
    λy = z[joint.inplane_y.reaction_variable]
    regularization = 1.0e-9
    sqrt(λx^2 + λy^2 + regularization^2) - regularization +
        friction.preload
end

function translational_friction_capacity(friction, z, slip_squared)
    scalar_friction_coefficient(friction, slip_squared) *
        max(zero(eltype(z)), z[friction.guide_load_variable])
end

function translational_friction_rate(friction, z)
    shear = z[friction.shear_variable]
    slip = SpatialConstraints.inline_velocity(friction.joint, z)
    capacity = translational_friction_capacity(friction, z, slip^2)
    scalar_bristle_rate(friction, shear, slip, capacity)
end

function calculated_translational_friction_force(friction, z)
    slip = SpatialConstraints.inline_velocity(friction.joint, z)
    capacity = translational_friction_capacity(friction, z, slip^2)
    shear = z[friction.shear_variable]
    scalar_bristle_force(friction, shear, slip, capacity)
end

function initialize_spatial_translational_friction!(state, friction;
        reset_anchor = false)
    reset_anchor && (friction.anchor[] =
        SpatialConstraints.inline_distance(friction.joint, state) -
        state[friction.shear_variable])
    state[friction.slip_variable] =
        SpatialConstraints.inline_velocity(friction.joint, state)
    state[friction.guide_load_variable] =
        translational_guide_load(friction, state)
    state[friction.force_variable] =
        calculated_translational_friction_force(friction, state)
    state[friction.global_force_variables] .=
        state[friction.force_variable] .*
        SpatialConstraints.inline_axis(friction.joint, state)
    state
end

function set_spatial_translational_friction_stage!(friction, stage)
    friction.stage = stage
    friction
end

function translational_friction_residual(friction, z, zdot)
    shear = z[friction.shear_variable]
    shear_equation = friction.stage == :static ?
        shear - (SpatialConstraints.inline_distance(friction.joint, z) -
            friction.anchor[]) :
        zdot[friction.shear_variable] - translational_friction_rate(friction, z)
    slip = SpatialConstraints.inline_velocity(friction.joint, z)
    guide_load = translational_guide_load(friction, z)
    force = calculated_translational_friction_force(friction, z)
    [shear_equation;
     z[friction.slip_variable] - slip;
     z[friction.guide_load_variable] - guide_load;
     z[friction.force_variable] - force;
     z[friction.global_force_variables] .-
        z[friction.force_variable] .*
        SpatialConstraints.inline_axis(friction.joint, z)]
end

function executable_blocks(friction::SpatialTranslationalFriction)
    joint = friction.joint
    columns = sort!(unique!([
        SpatialPlaneContacts.contact_kinematic_dependencies(joint.marker_i);
        SpatialPlaneContacts.contact_kinematic_dependencies(joint.marker_j);
        joint.inplane_x.reaction_variable;
        joint.inplane_y.reaction_variable;
        friction.shear_variable;
        friction.slip_variable;
        friction.guide_load_variable;
        friction.force_variable;
        collect(friction.global_force_variables)]))
    residual! = (equations, t, z, zdot) ->
        (equations[friction.friction_equations] .=
            translational_friction_residual(friction, z, zdot))
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        partials = SpatialPlaneContacts.local_state_jacobian(
            state -> translational_friction_residual(friction, state, zdot),
            z, columns)
        jacobian[friction.friction_equations, columns] .+= partials
        friction.stage == :static ||
            (jacobian[first(friction.friction_equations),
                friction.shear_variable] += coefficient)
    end
    [ExecutableEquationBlock(friction.name, :friction,
        collect(friction.friction_equations), residual!, jacobian!)]
end

function equation_contributions(friction::SpatialTranslationalFriction)
    first_marker = friction.joint.marker_i
    second_marker = friction.joint.marker_j
    rows = Int[]
    first_marker isa SpatialBodyMarker &&
        append!(rows, first_marker.body.balance_equations)
    second_marker isa SpatialBodyMarker &&
        append!(rows, second_marker.body.balance_equations)
    columns = sort!(unique!([
        SpatialPlaneContacts.contact_configuration_dependencies(first_marker);
        SpatialPlaneContacts.contact_configuration_dependencies(second_marker);
        collect(friction.global_force_variables)]))
    contribution = function (z)
        point = spatial_marker_position(first_marker, z)
        force = z[friction.global_force_variables]
        values = eltype(z)[]
        first_marker isa SpatialBodyMarker && append!(values,
            SpatialPlaneContacts.body_force_contribution(first_marker.body,
                z, point, force))
        second_marker isa SpatialBodyMarker && append!(values,
            SpatialPlaneContacts.body_force_contribution(second_marker.body,
                z, point, -force))
        values
    end
    residual! = (equations, t, z, zdot) ->
        (equations[rows] .+= contribution(z))
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        partials = SpatialPlaneContacts.local_state_jacobian(contribution,
            z, columns)
        jacobian[rows, columns] .+= partials
    end
    [EquationContribution(friction.name, :load_to_bodies, rows,
        residual!, jacobian!)]
end

"""Two tangential bristle states acting on an ideal inplane constraint."""
mutable struct SpatialInplaneFriction{C}
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
    anchor::Vector{Float64}
    shear_variables::UnitRange{Int}
    slip_variables::UnitRange{Int}
    normal_load_variable::Int
    force_variables::UnitRange{Int}
    global_force_variables::UnitRange{Int}
    friction_equations::UnitRange{Int}
end

function spatial_inplane_friction_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:shear_x, :user_state_steady, 0),
        VariableDeclaration(:shear_y, :user_state_steady, 0),
        VariableDeclaration(:slip_x, :applied_rate, 1),
        VariableDeclaration(:slip_y, :applied_rate, 1),
        VariableDeclaration(:normal_load, :applied_load, 2),
        VariableDeclaration(:friction_x, :applied_load, 2),
        VariableDeclaration(:friction_y, :applied_load, 2),
        [VariableDeclaration(Symbol(:F_, axis), :applied_load, 2)
            for axis in (:x, :y, :z)]...]
    equations = EquationDeclaration[
        EquationDeclaration(:shear_x, :user_differential_steady, 0, :friction),
        EquationDeclaration(:shear_y, :user_differential_steady, 0, :friction),
        EquationDeclaration(:slip_x, :applied_definition, 1, :friction),
        EquationDeclaration(:slip_y, :applied_definition, 1, :friction),
        EquationDeclaration(:normal_load, :applied_definition, 2, :friction),
        EquationDeclaration(:friction_x, :applied_definition, 2, :friction),
        EquationDeclaration(:friction_y, :applied_definition, 2, :friction),
        [EquationDeclaration(Symbol(:global_force_, axis),
            :applied_definition, 2, :friction)
            for axis in (:x, :y, :z)]...]
    ComponentRegistration(name, variables,
        [EquationBlockDeclaration(:friction, equations)])
end

component_registration(friction::SpatialInplaneFriction) =
    spatial_inplane_friction_registration(friction.name)

function allocated_spatial_inplane_friction(layout, name, constraint,
        stiffness, damping, preload, static_coefficient, dynamic_coefficient,
        transition_speed, release_time)
    variables = component_variable_indices(layout, name)
    SpatialInplaneFriction(name, constraint, stiffness, damping, preload,
        static_coefficient, dynamic_coefficient, transition_speed,
        release_time, :dynamic, zeros(2), variables[1:2], variables[3:4],
        variables[5], variables[6:7], variables[8:10],
        component_equation_indices(layout, name, :friction))
end

inplane_markers(friction) = (friction.constraint.geometry.marker_i,
    friction.constraint.geometry.marker_j)

function inplane_tangent_axes(friction, z)
    _, second = inplane_markers(friction)
    spatial_marker_orientation(second, z)[:, 1:2]
end

function inplane_tangent_position(friction, z)
    first, second = inplane_markers(friction)
    transpose(inplane_tangent_axes(friction, z)) *
        (spatial_marker_position(first, z) -
            spatial_marker_position(second, z))
end

function inplane_slip(friction, z)
    first, second = inplane_markers(friction)
    displacement = spatial_marker_position(first, z) -
        spatial_marker_position(second, z)
    velocity = spatial_marker_velocity(first, z) -
        spatial_marker_velocity(second, z)
    if second isa SpatialBodyMarker
        body = second.body
        orientation = rotation_matrix(z[body.euler_parameter_variables])
        omega = orientation * z[body.angular_velocity_variables]
        velocity -= cross(omega, displacement)
    end
    transpose(inplane_tangent_axes(friction, z)) * velocity
end

function inplane_normal_load(friction, z)
    reaction = z[friction.constraint.reaction_variable]
    regularization = 1.0e-9
    sqrt(reaction^2 + regularization^2) - regularization + friction.preload
end

function inplane_friction_capacity(friction, z, speed_squared)
    scalar_friction_coefficient(friction, speed_squared) *
        max(zero(eltype(z)), z[friction.normal_load_variable])
end

function inplane_friction_rates(friction, z)
    shear = z[friction.shear_variables]
    velocity = inplane_slip(friction, z)
    capacity = inplane_friction_capacity(friction, z,
        sum(abs2, velocity))
    vector_bristle_rate(friction, shear, velocity, capacity)
end

function calculated_inplane_friction_force(friction, z)
    velocity = inplane_slip(friction, z)
    capacity = inplane_friction_capacity(friction, z,
        sum(abs2, velocity))
    shear = z[friction.shear_variables]
    vector_bristle_force(friction, shear, velocity, capacity)
end

function initialize_spatial_inplane_friction!(state, friction;
        reset_anchor = false)
    reset_anchor && (friction.anchor .=
        inplane_tangent_position(friction, state) .-
        state[friction.shear_variables])
    state[friction.slip_variables] .= inplane_slip(friction, state)
    state[friction.normal_load_variable] = inplane_normal_load(friction, state)
    state[friction.force_variables] .=
        calculated_inplane_friction_force(friction, state)
    state[friction.global_force_variables] .=
        inplane_tangent_axes(friction, state) * state[friction.force_variables]
    state
end

function set_spatial_inplane_friction_stage!(friction, stage)
    friction.stage = stage
    friction
end

function inplane_friction_residual(friction, z, zdot)
    shear = z[friction.shear_variables]
    shear_equations = friction.stage == :static ?
        shear .- (inplane_tangent_position(friction, z) .-
            friction.anchor) :
        zdot[friction.shear_variables] .- inplane_friction_rates(friction, z)
    [shear_equations;
     z[friction.slip_variables] .- inplane_slip(friction, z);
     z[friction.normal_load_variable] - inplane_normal_load(friction, z);
     z[friction.force_variables] .-
        calculated_inplane_friction_force(friction, z);
     z[friction.global_force_variables] .-
        inplane_tangent_axes(friction, z) * z[friction.force_variables]]
end

function executable_blocks(friction::SpatialInplaneFriction)
    first, second = inplane_markers(friction)
    columns = sort!(unique!([
        SpatialPlaneContacts.contact_kinematic_dependencies(first);
        SpatialPlaneContacts.contact_kinematic_dependencies(second);
        friction.constraint.reaction_variable;
        collect(friction.shear_variables);
        collect(friction.slip_variables);
        friction.normal_load_variable;
        collect(friction.force_variables);
        collect(friction.global_force_variables)]))
    residual! = (equations, t, z, zdot) ->
        (equations[friction.friction_equations] .=
            inplane_friction_residual(friction, z, zdot))
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        partials = SpatialPlaneContacts.local_state_jacobian(
            state -> inplane_friction_residual(friction, state, zdot),
            z, columns)
        jacobian[friction.friction_equations, columns] .+= partials
        if friction.stage != :static
            for (row, variable) in zip(friction.friction_equations[1:2],
                    friction.shear_variables)
                jacobian[row, variable] += coefficient
            end
        end
    end
    [ExecutableEquationBlock(friction.name, :friction,
        collect(friction.friction_equations), residual!, jacobian!)]
end

function equation_contributions(friction::SpatialInplaneFriction)
    first, second = inplane_markers(friction)
    rows = Int[]
    first isa SpatialBodyMarker && append!(rows, first.body.balance_equations)
    second isa SpatialBodyMarker && append!(rows, second.body.balance_equations)
    columns = sort!(unique!([
        SpatialPlaneContacts.contact_configuration_dependencies(first);
        SpatialPlaneContacts.contact_configuration_dependencies(second);
        collect(friction.global_force_variables)]))
    contribution = function (z)
        point = spatial_marker_position(first, z)
        force = z[friction.global_force_variables]
        values = eltype(z)[]
        first isa SpatialBodyMarker && append!(values,
            SpatialPlaneContacts.body_force_contribution(first.body,
                z, point, force))
        second isa SpatialBodyMarker && append!(values,
            SpatialPlaneContacts.body_force_contribution(second.body,
                z, point, -force))
        values
    end
    residual! = (equations, t, z, zdot) ->
        (equations[rows] .+= contribution(z))
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        partials = SpatialPlaneContacts.local_state_jacobian(contribution,
            z, columns)
        jacobian[rows, columns] .+= partials
    end
    [EquationContribution(friction.name, :load_to_bodies, rows,
        residual!, jacobian!)]
end

end
