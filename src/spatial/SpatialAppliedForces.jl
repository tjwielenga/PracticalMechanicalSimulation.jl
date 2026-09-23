"""
Applied forces and torques for spatial models.

The named application side receives the action load and an optional second
side receives the equal-and-opposite reaction. Constitutive magnitudes and
their global vectors remain explicit canonical variables, so their definitions
and reactions participate in initialization, statics, dynamics, and modes.
"""
module SpatialAppliedForces

using LinearAlgebra
using ..AutomaticAnalysis
using ..ScalarExpressions: ScalarLaw
using ..SpatialComponentAssembly
using ..SpatialModeling
using ..SpatialSpans
using ..SpatialDirectedDistances
using ..SpatialConstraints: body_force_rows, body_torque_rows,
    add_reaction_to_body!,
    add_reaction_jacobian!

import ..SpatialComponentAssembly: component_registration,
    executable_blocks, equation_contributions

export SpatialAppliedForceComponent, spatial_applied_force_registration,
       allocated_spatial_applied_force, initialize_spatial_applied_force!,
       spatial_applied_force_value, set_spatial_applied_force_stage!,
       SpatialDirectedTorqueComponent, spatial_directed_torque_registration,
       allocated_spatial_directed_torque,
       initialize_spatial_directed_torque!, spatial_directed_torque_value,
       set_spatial_directed_torque_stage!,
       SpatialAppliedTorqueComponent, spatial_applied_torque_registration,
       allocated_spatial_applied_torque, initialize_spatial_applied_torque!,
       spatial_applied_torque_value, set_spatial_applied_torque_stage!,
       SpatialSpanningForceComponent, spatial_spanning_force_registration,
       allocated_spatial_spanning_force, initialize_spatial_spanning_force!,
       spatial_spanning_values, set_spatial_spanning_force_stage!

"""
Force applied at the first marker along the second marker's local z-axis. An
optional floating marker carries the opposite force on a reaction body.
"""
struct SpatialAppliedForceComponent{A,D,R,L,S}
    name::Symbol
    application_marker::A
    direction_axis::D
    reaction_marker::R
    magnitude::L
    active_during::S
    active::Base.RefValue{Bool}
    magnitude_variable::Int
    global_force_variables::UnitRange{Int}
    magnitude_equation::Int
    global_force_equations::UnitRange{Int}
end

"""
Torque applied to a body marker along a direction marker's local z-axis. An
optional floating marker carries the opposite torque on a reaction body.
"""
struct SpatialDirectedTorqueComponent{A,D,R,L,S}
    name::Symbol
    application_marker::A
    direction_axis::D
    reaction_marker::R
    magnitude::L
    active_during::S
    active::Base.RefValue{Bool}
    magnitude_variable::Int
    global_torque_variables::UnitRange{Int}
    magnitude_equation::Int
    global_torque_equations::UnitRange{Int}
end

"""
Torque applied to the first side of a hinge or revolute joint around the
second marker's local z-axis. The second side receives the opposite torque.
"""
struct SpatialAppliedTorqueComponent{H,L,T,R,S}
    name::Symbol
    hinge::H
    magnitude::L
    initial_free_angle::Bool
    stiffness::T
    damping::T
    free_angle::R
    active_during::S
    active::Base.RefValue{Bool}
    magnitude_variable::Int
    global_torque_variables::UnitRange{Int}
    magnitude_equation::Int
    global_torque_equations::UnitRange{Int}
end

function spatial_applied_torque_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:torque, :applied_load, 2);
        [VariableDeclaration(Symbol(:T_, axis), :applied_load, 2)
            for axis in (:x, :y, :z)]
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:load, EquationDeclaration[
            EquationDeclaration(:magnitude, :applied_definition, 2,
                :constitutive_torque);
            [EquationDeclaration(Symbol(:global_torque_, axis),
                :applied_definition, 2, :global_torque)
                for axis in (:x, :y, :z)]
        ]),
    ]
    ComponentRegistration(name, variables, blocks)
end

component_registration(torque::SpatialAppliedTorqueComponent) =
    spatial_applied_torque_registration(torque.name)

function allocated_spatial_applied_torque(layout, name, hinge,
        magnitude::ScalarLaw; initial_free_angle = false,
        stiffness = 0.0, damping = 0.0, free_angle = Ref(0.0),
        active_during = (:static, :dynamic, :modal))
    variables = component_variable_indices(layout, name)
    equations = component_equation_indices(layout, name, :load)
    SpatialAppliedTorqueComponent(name, hinge, magnitude,
        initial_free_angle, stiffness, damping, free_angle,
        active_during, Ref(:dynamic in active_during),
        variables[1], variables[2:4], equations[1], equations[2:4])
end

"""Select whether this applied torque is active in the analysis stage."""
function set_spatial_applied_torque_stage!(torque::SpatialAppliedTorqueComponent,
        stage)
    torque.active[] = stage in torque.active_during
    torque
end

function initialize_spatial_applied_torque!(initial, torque, time;
        reset_free_angle = false)
    if reset_free_angle && torque.initial_free_angle
        torque.free_angle[] =
            initial[torque.hinge.rotation_variables[3]]
    end
    magnitude = torque.active[] ? torque.magnitude(time, initial) : 0.0
    isfinite(magnitude) || throw(ArgumentError(
        "applied torque '$(torque.name)' is not finite initially"))
    direction = spatial_marker_orientation(
        torque.hinge.marker_j, initial)[:, 3]
    initial[torque.magnitude_variable] = magnitude
    initial[torque.global_torque_variables] .= magnitude .* direction
    initial
end

spatial_applied_torque_value(torque::SpatialAppliedTorqueComponent, z) =
    collect(@view z[torque.global_torque_variables])

function executable_blocks(torque::SpatialAppliedTorqueComponent)
    residual! = function (equations, t, z, zdot)
        magnitude = z[torque.magnitude_variable]
        direction = marker_axis_kinematics(torque.hinge.marker_j, z, 3)
        equations[torque.magnitude_equation] =
            magnitude - (torque.active[] ? torque.magnitude(t, z) : 0.0)
        equations[torque.global_torque_equations] .=
            z[torque.global_torque_variables] .-
            magnitude .* direction.direction
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        magnitude = z[torque.magnitude_variable]
        direction = marker_axis_kinematics(torque.hinge.marker_j, z, 3)
        jacobian[torque.magnitude_equation, torque.magnitude_variable] += 1
        if torque.active[]
            for (column, partial) in zip(torque.magnitude.dependencies,
                    torque.magnitude.gradient(t, z))
                jacobian[torque.magnitude_equation, column] -= partial
            end
        end
        for component in 1:3
            row = torque.global_torque_equations[component]
            jacobian[row, torque.global_torque_variables[component]] += 1
            jacobian[row, torque.magnitude_variable] -=
                direction.direction[component]
        end
        if !isnothing(direction.body)
            jacobian[torque.global_torque_equations,
                direction.body.euler_parameter_variables] .-=
                magnitude .* direction.direction_parameters
        end
    end
    ExecutableEquationBlock[ExecutableEquationBlock(torque.name, :load,
        [torque.magnitude_equation;
         collect(torque.global_torque_equations)], residual!, jacobian!)]
end

function add_applied_torque_to_body!(equations, z, marker, global_torque,
        sign)
    marker isa SpatialGroundMarker && return nothing
    body = marker.body
    parameters = @view z[body.euler_parameter_variables]
    torque_body = transpose(rotation_matrix(parameters)) *
        (sign .* global_torque)
    equations[body.balance_equations[4:6]] .-= torque_body
    if marker isa SpatialFlexibleBeamMarker
        equations[body.balance_equations[7:12]] .-=
            transpose(marker.orientation_shape) * torque_body
    end
    nothing
end

function add_applied_torque_jacobian!(jacobian, z, marker, global_torque,
        sign, torque_variables)
    marker isa SpatialGroundMarker && return nothing
    body = marker.body
    parameters = @view z[body.euler_parameter_variables]
    orientation = rotation_matrix(parameters)
    rows = body.balance_equations[4:6]
    jacobian[rows, torque_variables] .-=
        sign .* transpose(orientation)
    jacobian[rows, body.euler_parameter_variables] .-=
        rotation_transpose_vector_jacobian(
            parameters, sign .* global_torque)
    if marker isa SpatialFlexibleBeamMarker
        elastic_rows = body.balance_equations[7:12]
        shape_transpose = transpose(marker.orientation_shape)
        jacobian[elastic_rows, torque_variables] .-=
            sign .* shape_transpose * transpose(orientation)
        jacobian[elastic_rows, body.euler_parameter_variables] .-=
            shape_transpose * rotation_transpose_vector_jacobian(
                parameters, sign .* global_torque)
    end
    nothing
end

function equation_contributions(torque::SpatialAppliedTorqueComponent)
    rows = unique!([body_torque_rows(torque.hinge.marker_i);
                    body_torque_rows(torque.hinge.marker_j)])
    residual! = function (equations, t, z, zdot)
        global_torque = @view z[torque.global_torque_variables]
        add_applied_torque_to_body!(equations, z, torque.hinge.marker_i,
            global_torque, 1)
        add_applied_torque_to_body!(equations, z, torque.hinge.marker_j,
            global_torque, -1)
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        global_torque = @view z[torque.global_torque_variables]
        add_applied_torque_jacobian!(jacobian, z, torque.hinge.marker_i,
            global_torque, 1, torque.global_torque_variables)
        add_applied_torque_jacobian!(jacobian, z, torque.hinge.marker_j,
            global_torque, -1, torque.global_torque_variables)
    end
    EquationContribution[EquationContribution(torque.name,
        :torque_to_bodies, rows, residual!, jacobian!)]
end

function spatial_applied_force_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:force, :applied_load, 2);
        [VariableDeclaration(Symbol(:F_, axis), :applied_load, 2)
            for axis in (:x, :y, :z)]
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:load, EquationDeclaration[
            EquationDeclaration(:magnitude, :applied_definition, 2,
                :constitutive_force);
            [EquationDeclaration(Symbol(:global_force_, axis),
                :applied_definition, 2, :global_force)
                for axis in (:x, :y, :z)]
        ]),
    ]
    ComponentRegistration(name, variables, blocks)
end

component_registration(force::SpatialAppliedForceComponent) =
    spatial_applied_force_registration(force.name)

function allocated_spatial_applied_force(layout, name, application_marker,
        direction_marker, reaction_marker, magnitude::ScalarLaw,
        active_during = (:static, :dynamic, :modal))
    variables = component_variable_indices(layout, name)
    equations = component_equation_indices(layout, name, :load)
    axis = SpatialDirectedAxis(direction_marker, 3)
    SpatialAppliedForceComponent(name, application_marker, axis,
        reaction_marker, magnitude, active_during,
        Ref(:dynamic in active_during), variables[1], variables[2:4],
        equations[1], equations[2:4])
end

"""Select whether this applied force is active in the analysis stage."""
function set_spatial_applied_force_stage!(force::SpatialAppliedForceComponent,
        stage)
    force.active[] = stage in force.active_during
    force
end

function initialize_spatial_applied_force!(initial, force, time)
    magnitude = force.active[] ? force.magnitude(time, initial) : 0.0
    isfinite(magnitude) || throw(ArgumentError(
        "applied force '$(force.name)' is not finite initially"))
    direction = directed_axis_values(force.direction_axis, initial).direction
    initial[force.magnitude_variable] = magnitude
    initial[force.global_force_variables] .= magnitude .* direction
    initial
end

spatial_applied_force_value(force::SpatialAppliedForceComponent, z) =
    collect(@view z[force.global_force_variables])

function executable_blocks(force::SpatialAppliedForceComponent)
    residual! = function (equations, t, z, zdot)
        magnitude = z[force.magnitude_variable]
        direction = directed_axis_values(force.direction_axis, z)
        equations[force.magnitude_equation] =
            magnitude - (force.active[] ? force.magnitude(t, z) : 0.0)
        equations[force.global_force_equations] .=
            z[force.global_force_variables] .- magnitude .* direction.direction
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        magnitude = z[force.magnitude_variable]
        direction = directed_axis_values(force.direction_axis, z)
        jacobian[force.magnitude_equation, force.magnitude_variable] += 1
        if force.active[]
            for (column, partial) in zip(force.magnitude.dependencies,
                    force.magnitude.gradient(t, z))
                jacobian[force.magnitude_equation, column] -= partial
            end
        end
        for component in 1:3
            row = force.global_force_equations[component]
            jacobian[row, force.global_force_variables[component]] += 1
            jacobian[row, force.magnitude_variable] -=
                direction.direction[component]
        end
        if !isnothing(direction.body)
            jacobian[force.global_force_equations,
                direction.body.euler_parameter_variables] .-=
                magnitude .* direction.direction_parameters
        end
    end
    ExecutableEquationBlock[ExecutableEquationBlock(force.name, :load,
        [force.magnitude_equation;
         collect(force.global_force_equations)], residual!, jacobian!)]
end

function add_floating_force!(equations, z, marker::SpatialFloatingMarker,
        global_force)
    body = marker.body
    parameters = @view z[body.euler_parameter_variables]
    orientation = rotation_matrix(parameters)
    lever_global = spatial_marker_position(marker, z) -
        z[body.position_variables]
    force_body = transpose(orientation) * global_force
    lever_body = transpose(orientation) * lever_global
    equations[body.balance_equations[1:3]] .-= global_force
    equations[body.balance_equations[4:6]] .-=
        cross(lever_body, force_body)
    nothing
end

function add_floating_force_jacobian!(jacobian, z,
        marker::SpatialFloatingMarker, global_force, force_variables, sign)
    body = marker.body
    follower = marker.follower
    parameters = @view z[body.euler_parameter_variables]
    orientation = rotation_matrix(parameters)
    transpose_orientation = transpose(orientation)
    lever_global = spatial_marker_position(marker, z) -
        z[body.position_variables]
    signed_force = sign .* global_force
    lever_body = transpose_orientation * lever_global
    force_body = transpose_orientation * signed_force
    force_rows = body.balance_equations[1:3]
    torque_rows = body.balance_equations[4:6]

    jacobian[force_rows, force_variables] .-=
        sign .* Matrix{eltype(z)}(I, 3, 3)
    jacobian[torque_rows, force_variables] .-=
        sign .* skew(lever_body) * transpose_orientation
    jacobian[torque_rows, body.position_variables] .-=
        skew(force_body) * transpose_orientation

    lever_parameters = rotation_transpose_vector_jacobian(
        parameters, lever_global)
    force_parameters = rotation_transpose_vector_jacobian(
        parameters, signed_force)
    jacobian[torque_rows, body.euler_parameter_variables] .+=
        skew(force_body) * lever_parameters -
        skew(lever_body) * force_parameters

    follower_values = marker_point_kinematics(follower, z)
    if !isnothing(follower_values.body)
        follower_body = follower_values.body
        jacobian[torque_rows, follower_body.position_variables] .+=
            skew(force_body) * transpose_orientation
        jacobian[torque_rows, follower_body.euler_parameter_variables] .+=
            skew(force_body) * transpose_orientation *
            follower_values.position_parameters
    end
    nothing
end

function equation_contributions(force::SpatialAppliedForceComponent)
    rows = body_force_rows(force.application_marker)
    if !isnothing(force.reaction_marker)
        append!(rows, collect(force.reaction_marker.body.balance_equations))
    end
    unique!(rows)
    residual! = function (equations, t, z, zdot)
        global_force = @view z[force.global_force_variables]
        add_reaction_to_body!(equations, z, force.application_marker,
            global_force, 1)
        if !isnothing(force.reaction_marker)
            add_floating_force!(equations, z, force.reaction_marker,
                -global_force)
        end
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        global_force = @view z[force.global_force_variables]
        add_reaction_jacobian!(jacobian, z, force.application_marker,
            global_force, 1, force.global_force_variables)
        if !isnothing(force.reaction_marker)
            add_floating_force_jacobian!(jacobian, z, force.reaction_marker,
                global_force, force.global_force_variables, -1)
        end
    end
    EquationContribution[EquationContribution(force.name, :force_to_bodies,
        rows, residual!, jacobian!)]
end

spatial_directed_torque_registration(name::Symbol) =
    spatial_applied_torque_registration(name)

component_registration(torque::SpatialDirectedTorqueComponent) =
    spatial_directed_torque_registration(torque.name)

function allocated_spatial_directed_torque(layout, name, application_marker,
        direction_marker, reaction_marker, magnitude::ScalarLaw,
        active_during = (:static, :dynamic, :modal))
    variables = component_variable_indices(layout, name)
    equations = component_equation_indices(layout, name, :load)
    axis = SpatialDirectedAxis(direction_marker, 3)
    SpatialDirectedTorqueComponent(name, application_marker, axis,
        reaction_marker, magnitude, active_during,
        Ref(:dynamic in active_during), variables[1], variables[2:4],
        equations[1], equations[2:4])
end

"""Select whether this directed torque is active in the analysis stage."""
function set_spatial_directed_torque_stage!(
        torque::SpatialDirectedTorqueComponent, stage)
    torque.active[] = stage in torque.active_during
    torque
end

function initialize_spatial_directed_torque!(initial, torque, time)
    magnitude = torque.active[] ? torque.magnitude(time, initial) : 0.0
    isfinite(magnitude) || throw(ArgumentError(
        "directed torque '$(torque.name)' is not finite initially"))
    direction = directed_axis_values(torque.direction_axis, initial).direction
    initial[torque.magnitude_variable] = magnitude
    initial[torque.global_torque_variables] .= magnitude .* direction
    initial
end

spatial_directed_torque_value(torque::SpatialDirectedTorqueComponent, z) =
    collect(@view z[torque.global_torque_variables])

function executable_blocks(torque::SpatialDirectedTorqueComponent)
    residual! = function (equations, t, z, zdot)
        magnitude = z[torque.magnitude_variable]
        direction = directed_axis_values(torque.direction_axis, z)
        equations[torque.magnitude_equation] = magnitude -
            (torque.active[] ? torque.magnitude(t, z) : 0.0)
        equations[torque.global_torque_equations] .=
            z[torque.global_torque_variables] .-
            magnitude .* direction.direction
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        magnitude = z[torque.magnitude_variable]
        direction = directed_axis_values(torque.direction_axis, z)
        jacobian[torque.magnitude_equation, torque.magnitude_variable] += 1
        if torque.active[]
            for (column, partial) in zip(torque.magnitude.dependencies,
                    torque.magnitude.gradient(t, z))
                jacobian[torque.magnitude_equation, column] -= partial
            end
        end
        for component in 1:3
            row = torque.global_torque_equations[component]
            jacobian[row, torque.global_torque_variables[component]] += 1
            jacobian[row, torque.magnitude_variable] -=
                direction.direction[component]
        end
        if !isnothing(direction.body)
            jacobian[torque.global_torque_equations,
                direction.body.euler_parameter_variables] .-=
                magnitude .* direction.direction_parameters
        end
    end
    ExecutableEquationBlock[ExecutableEquationBlock(torque.name, :load,
        [torque.magnitude_equation;
         collect(torque.global_torque_equations)], residual!, jacobian!)]
end

function equation_contributions(torque::SpatialDirectedTorqueComponent)
    rows = body_torque_rows(torque.application_marker)
    if !isnothing(torque.reaction_marker)
        append!(rows,
            collect(torque.reaction_marker.body.balance_equations[4:6]))
    end
    unique!(rows)
    residual! = function (equations, t, z, zdot)
        global_torque = @view z[torque.global_torque_variables]
        add_applied_torque_to_body!(equations, z,
            torque.application_marker, global_torque, 1)
        if !isnothing(torque.reaction_marker)
            add_applied_torque_to_body!(equations, z,
                torque.reaction_marker, global_torque, -1)
        end
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        global_torque = @view z[torque.global_torque_variables]
        add_applied_torque_jacobian!(jacobian, z,
            torque.application_marker, global_torque, 1,
            torque.global_torque_variables)
        if !isnothing(torque.reaction_marker)
            add_applied_torque_jacobian!(jacobian, z,
                torque.reaction_marker, global_torque, -1,
                torque.global_torque_variables)
        end
    end
    EquationContribution[EquationContribution(torque.name,
        :torque_to_bodies, rows, residual!, jacobian!)]
end

"""Axial force acting through a shared spatial span coordinate."""
struct SpatialSpanningForceComponent{S,L,T,A}
    name::Symbol
    span::S
    law::L
    stiffness::T
    damping::T
    free_length::T
    active_during::A
    active::Base.RefValue{Bool}
    force_variable::Int
    global_force_variables::UnitRange{Int}
    force_equation::Int
    global_force_equations::UnitRange{Int}
end

function spatial_spanning_force_registration(name::Symbol)
    variables = spatial_span_variable_declarations(distance_name = :length,
        velocity_name = :length_rate,
        acceleration_name = :length_acceleration)
    push!(variables, VariableDeclaration(:force, :applied_load, 2))
    append!(variables, [VariableDeclaration(Symbol(:F_, axis),
        :applied_load, 2) for axis in (:x, :y, :z)])
    blocks = spatial_span_equation_blocks()
    push!(blocks, EquationBlockDeclaration(:load, EquationDeclaration[
            EquationDeclaration(:scalar_force, :applied_definition, 2,
                :constitutive_force);
            [EquationDeclaration(Symbol(:global_force_, axis),
                :applied_definition, 2, :global_force)
                for axis in (:x, :y, :z)];
        ]))
    ComponentRegistration(name, variables, blocks)
end

component_registration(force::SpatialSpanningForceComponent) =
    spatial_spanning_force_registration(force.name)

function allocated_spatial_spanning_force(layout, name, marker_1, marker_2,
        law::ScalarLaw; stiffness = 0.0, damping = 0.0, free_length = 0.0,
        active_during = (:static, :dynamic, :modal))
    variables = component_variable_indices(layout, name)
    load = component_equation_indices(layout, name, :load)
    span = allocated_spatial_span(layout, name, marker_1, marker_2)
    SpatialSpanningForceComponent(name, span, law,
        Float64(stiffness), Float64(damping), Float64(free_length),
        active_during, Ref(:dynamic in active_during),
        variables[10], variables[11:13], load[1], load[2:4])
end

"""Select whether this spanning force is active in the analysis stage."""
function set_spatial_spanning_force_stage!(force::SpatialSpanningForceComponent,
        stage)
    force.active[] = stage in force.active_during
    force
end

function spatial_spanning_values(force::SpatialSpanningForceComponent, z)
    values = spatial_span_values(force.span, z)
    (; values..., spanning = values.separation, length = values.distance,
       unit = values.direction, length_rate = values.velocity,
       length_acceleration = values.acceleration,
       scalar_force = z[force.force_variable],
       global_force = @view(z[force.global_force_variables]))
end

function initialize_spatial_spanning_force!(initial, force, time)
    initialize_spatial_span!(initial, force.span)
    initial[force.force_variable] =
        force.active[] ? force.law(time, initial) : 0.0
    isfinite(initial[force.force_variable]) || throw(ArgumentError(
        "spanning force '$(force.name)' is not finite initially"))
    initial[force.global_force_variables] .=
        -initial[force.span.direction_variables] .* initial[force.force_variable]
    initial
end

function executable_blocks(force::SpatialSpanningForceComponent)
    load! = function (equations, t, z, zdot)
        values = spatial_spanning_values(force, z)
        equations[force.force_equation] = values.scalar_force -
            (force.active[] ? force.law(t, z) : 0.0)
        equations[force.global_force_equations] .= values.global_force .+
            values.unit .* values.scalar_force
    end
    load_jacobian! = function (jacobian, t, z, zdot, coefficient)
        values = spatial_spanning_values(force, z)
        jacobian[force.force_equation, force.force_variable] += 1
        if force.active[]
            for (column, partial) in zip(force.law.dependencies,
                    force.law.gradient(t, z))
                jacobian[force.force_equation, column] -= partial
            end
        end
        for component in 1:3
            row = force.global_force_equations[component]
            jacobian[row, force.global_force_variables[component]] += 1
            jacobian[row, force.span.direction_variables[component]] +=
                values.scalar_force
            jacobian[row, force.force_variable] += values.unit[component]
        end
    end
    blocks = spatial_span_executable_blocks(force.span)
    push!(blocks, ExecutableEquationBlock(force.name, :load,
        [force.force_equation; collect(force.global_force_equations)],
        load!, load_jacobian!))
    blocks
end

function equation_contributions(force::SpatialSpanningForceComponent)
    rows = unique!([body_force_rows(force.span.marker_1);
                    body_force_rows(force.span.marker_2)])
    residual! = function (equations, t, z, zdot)
        global_force = @view z[force.global_force_variables]
        add_reaction_to_body!(equations, z, force.span.marker_1,
            global_force, 1)
        add_reaction_to_body!(equations, z, force.span.marker_2,
            global_force, -1)
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        global_force = @view z[force.global_force_variables]
        add_reaction_jacobian!(jacobian, z, force.span.marker_1,
            global_force, 1, force.global_force_variables)
        add_reaction_jacobian!(jacobian, z, force.span.marker_2,
            global_force, -1, force.global_force_variables)
    end
    EquationContribution[EquationContribution(force.name, :load_to_bodies,
        rows, residual!, jacobian!)]
end

end
