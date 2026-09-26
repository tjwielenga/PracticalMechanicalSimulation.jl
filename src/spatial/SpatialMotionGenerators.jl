"""Prescribed rotational, translational, and spanning motion elements."""
module SpatialMotionGenerators

using LinearAlgebra
using ..AutomaticAnalysis
using ..SpatialComponentAssembly
using ..SpatialModeling: is_flexible_marker
using ..SpatialSpans
using ..SpatialDirectedDistances
using ..SpatialConstraints: body_force_rows, body_torque_rows,
    directed_distance_reaction_rows, add_directed_distance_reaction!,
    add_directed_distance_reaction_jacobian!, add_reaction_to_body!,
    add_reaction_jacobian!, add_perp_reaction_to_body!,
    add_perp_reaction_jacobian!, constraint_dependencies, add_ad_jacobian!

import ..SpatialComponentAssembly: component_registration,
    executable_blocks, equation_contributions

export SpatialRotationalMotionGenerator,
       spatial_rotational_motion_registration,
       allocated_spatial_rotational_motion_generator,
       SpatialTranslationalMotionGenerator,
       spatial_translational_motion_registration,
       allocated_spatial_translational_motion_generator,
       SpatialSpanningMotionGenerator,
       spatial_spanning_motion_registration,
       allocated_spatial_spanning_motion_generator

"""Prescribed rotation around the second marker's z-axis."""
struct SpatialRotationalMotionGenerator{H,F,F1,F2}
    name::Symbol
    hinge::H
    reaction_variable::Int
    acceleration_equation::Int
    velocity_equation::Int
    position_equation::Int
    motion::F
    motion_derivative::F1
    motion_second_derivative::F2
end

"""Prescribed signed marker distance along the second marker's z-axis."""
struct SpatialTranslationalMotionGenerator{G,F,F1,F2}
    name::Symbol
    geometry::G
    distance_variable::Int
    velocity_variable::Int
    acceleration_variable::Int
    reaction_variable::Int
    position_equations::UnitRange{Int}
    velocity_equations::UnitRange{Int}
    acceleration_equations::UnitRange{Int}
    motion::F
    motion_derivative::F1
    motion_second_derivative::F2
end

"""Prescribed positive distance between two spatial marker points."""
struct SpatialSpanningMotionGenerator{S,F,F1,F2}
    name::Symbol
    span::S
    reaction_variable::Int
    global_force_variables::UnitRange{Int}
    position_equations::UnitRange{Int}
    velocity_equations::UnitRange{Int}
    acceleration_equations::UnitRange{Int}
    global_force_equations::UnitRange{Int}
    motion::F
    motion_derivative::F1
    motion_second_derivative::F2
end

function spatial_rotational_motion_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:torque, :reaction, 2),
    ]
    equations = EquationDeclaration[
        EquationDeclaration(:prescribed_angular_acceleration, :motion, 2,
            :prescribed_rotation),
        EquationDeclaration(:prescribed_angular_velocity, :motion, 1,
            :prescribed_rotation),
        EquationDeclaration(:prescribed_angle, :motion, 0,
            :prescribed_rotation),
    ]
    ComponentRegistration(name, variables,
        [EquationBlockDeclaration(:constraint, equations)])
end

component_registration(generator::SpatialRotationalMotionGenerator) =
    spatial_rotational_motion_registration(generator.name)

function spatial_translational_motion_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:distance_m, :relative_position, 0),
        VariableDeclaration(:velocity_m, :relative_velocity, 1),
        VariableDeclaration(:acceleration_m, :relative_acceleration, 2),
        VariableDeclaration(:force_m, :reaction, 2),
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:position, EquationDeclaration[
            EquationDeclaration(:relative_distance, :coordinate_relation, 0,
                :relative_translation),
            EquationDeclaration(:prescribed_distance, :motion, 0,
                :prescribed_translation),
        ]),
        EquationBlockDeclaration(:velocity, EquationDeclaration[
            EquationDeclaration(:relative_velocity, :coordinate_relation, 1,
                :relative_translation),
            EquationDeclaration(:prescribed_velocity, :motion, 1,
                :prescribed_translation),
        ]),
        EquationBlockDeclaration(:acceleration, EquationDeclaration[
            EquationDeclaration(:relative_acceleration,
                :coordinate_relation, 2, :relative_translation),
            EquationDeclaration(:prescribed_acceleration, :motion, 2,
                :prescribed_translation),
        ]),
    ]
    ComponentRegistration(name, variables, blocks)
end

component_registration(generator::SpatialTranslationalMotionGenerator) =
    spatial_translational_motion_registration(generator.name)

function spatial_spanning_motion_registration(name::Symbol)
    variables = spatial_span_variable_declarations(
        distance_name = :distance_m, velocity_name = :velocity_m,
        acceleration_name = :acceleration_m)
    push!(variables, VariableDeclaration(:force_m, :reaction, 2))
    append!(variables, [VariableDeclaration(Symbol(:F_, axis), :reaction, 2)
        for axis in (:x, :y, :z)])
    blocks = spatial_span_equation_blocks()
    append!(blocks, EquationBlockDeclaration[
        EquationBlockDeclaration(:position, EquationDeclaration[
            EquationDeclaration(:prescribed_distance, :motion, 0,
                :prescribed_distance),
        ]),
        EquationBlockDeclaration(:velocity, EquationDeclaration[
            EquationDeclaration(:prescribed_velocity, :motion, 1,
                :prescribed_distance),
        ]),
        EquationBlockDeclaration(:acceleration, EquationDeclaration[
            EquationDeclaration(:prescribed_acceleration, :motion, 2,
                :prescribed_distance),
        ]),
        EquationBlockDeclaration(:load,
            [EquationDeclaration(Symbol(:global_force_, axis),
                :applied_definition, 2, :global_force)
                for axis in (:x, :y, :z)]),
    ])
    ComponentRegistration(name, variables, blocks)
end

component_registration(generator::SpatialSpanningMotionGenerator) =
    spatial_spanning_motion_registration(generator.name)

function allocated_spatial_rotational_motion_generator(layout, name, hinge,
        motion, motion_derivative, motion_second_derivative)
    reaction = only(component_variable_indices(layout, name))
    equations = component_equation_indices(layout, name, :constraint)
    SpatialRotationalMotionGenerator(name, hinge, reaction,
        equations[1], equations[2], equations[3], motion,
        motion_derivative, motion_second_derivative)
end

function allocated_spatial_translational_motion_generator(layout, name,
        marker_i, marker_j, motion, motion_derivative,
        motion_second_derivative)
    variables = component_variable_indices(layout, name)
    position_equations = component_equation_indices(layout, name, :position)
    velocity_equations = component_equation_indices(layout, name, :velocity)
    acceleration_equations = component_equation_indices(
        layout, name, :acceleration)
    geometry = SpatialDirectedDistance(
        marker_i, marker_j, SpatialDirectedAxis(marker_j, 3))
    SpatialTranslationalMotionGenerator(name, geometry,
        variables[1], variables[2], variables[3], variables[4],
        position_equations, velocity_equations, acceleration_equations,
        motion, motion_derivative, motion_second_derivative)
end

function allocated_spatial_spanning_motion_generator(layout, name,
        marker_1, marker_2, motion, motion_derivative,
        motion_second_derivative)
    variables = component_variable_indices(layout, name)
    span = allocated_spatial_span(layout, name, marker_1, marker_2)
    SpatialSpanningMotionGenerator(name, span, variables[10], variables[11:13],
        component_equation_indices(layout, name, :position),
        component_equation_indices(layout, name, :velocity),
        component_equation_indices(layout, name, :acceleration),
        component_equation_indices(layout, name, :load),
        motion, motion_derivative, motion_second_derivative)
end

function positive_spanning_motion(generator, time)
    value = generator.motion(time)
    isfinite(value) && value > 0 || throw(DomainError(value,
        "spanning motion '$(generator.name)' must remain positive"))
    value
end
function executable_blocks(generator::SpatialRotationalMotionGenerator)
    alpha, omega, theta = generator.hinge.rotation_variables
    residual! = function (equations, t, z, zdot)
        equations[generator.acceleration_equation] =
            z[alpha] - generator.motion_second_derivative(t)
        equations[generator.velocity_equation] =
            z[omega] - generator.motion_derivative(t)
        equations[generator.position_equation] =
            z[theta] - generator.motion(t)
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        jacobian[generator.acceleration_equation, alpha] += 1
        jacobian[generator.velocity_equation, omega] += 1
        jacobian[generator.position_equation, theta] += 1
    end
    ExecutableEquationBlock[ExecutableEquationBlock(generator.name,
        :constraint, [generator.acceleration_equation,
            generator.velocity_equation, generator.position_equation],
        residual!, jacobian!)]
end

function executable_blocks(generator::SpatialTranslationalMotionGenerator)
    position_relation, position_prescription = generator.position_equations
    velocity_relation, velocity_prescription = generator.velocity_equations
    acceleration_relation, acceleration_prescription =
        generator.acceleration_equations
    residual! = function (equations, t, z, zdot)
        values = directed_distance_values(generator.geometry, z)
        equations[position_relation] = z[generator.distance_variable] -
            values.position
        equations[position_prescription] = z[generator.distance_variable] -
            generator.motion(t)
        equations[velocity_relation] = z[generator.velocity_variable] -
            values.velocity
        equations[velocity_prescription] = z[generator.velocity_variable] -
            generator.motion_derivative(t)
        equations[acceleration_relation] =
            z[generator.acceleration_variable] -
            values.acceleration
        equations[acceleration_prescription] =
            z[generator.acceleration_variable] -
            generator.motion_second_derivative(t)
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        jacobian[position_relation, generator.distance_variable] += 1
        jacobian[position_prescription, generator.distance_variable] += 1
        jacobian[velocity_relation, generator.velocity_variable] += 1
        jacobian[velocity_prescription, generator.velocity_variable] += 1
        jacobian[acceleration_relation,
            generator.acceleration_variable] += 1
        jacobian[acceleration_prescription,
            generator.acceleration_variable] += 1
        directed_distance_jacobian!(jacobian, z, generator.geometry,
            acceleration_relation, velocity_relation, position_relation, -1)
    end
    rows = [collect(generator.position_equations);
            collect(generator.velocity_equations);
            collect(generator.acceleration_equations)]
    ExecutableEquationBlock[ExecutableEquationBlock(generator.name,
        :constraint, rows, residual!, jacobian!)]
end

function executable_blocks(generator::SpatialSpanningMotionGenerator)
    position_prescription = only(generator.position_equations)
    velocity_prescription = only(generator.velocity_equations)
    acceleration_prescription = only(generator.acceleration_equations)
    prescription! = function (equations, t, z, zdot)
        equations[position_prescription] =
            z[generator.span.distance_variable] -
            positive_spanning_motion(generator, t)
        equations[velocity_prescription] =
            z[generator.span.velocity_variable] -
            generator.motion_derivative(t)
        equations[acceleration_prescription] =
            z[generator.span.acceleration_variable] -
            generator.motion_second_derivative(t)
    end
    prescription_jacobian! = function (jacobian, t, z, zdot, coefficient)
        jacobian[position_prescription,
            generator.span.distance_variable] += 1
        jacobian[velocity_prescription,
            generator.span.velocity_variable] += 1
        jacobian[acceleration_prescription,
            generator.span.acceleration_variable] += 1
    end
    load! = function (equations, t, z, zdot)
        direction = @view z[generator.span.direction_variables]
        equations[generator.global_force_equations] .=
            z[generator.global_force_variables] .+
            z[generator.reaction_variable] .* direction
    end
    load_jacobian! = function (jacobian, t, z, zdot, coefficient)
        direction = @view z[generator.span.direction_variables]
        force_rows = generator.global_force_equations
        jacobian[force_rows, generator.global_force_variables] .+= I(3)
        jacobian[force_rows, generator.reaction_variable] .+= direction
        jacobian[force_rows, generator.span.direction_variables] .+=
            z[generator.reaction_variable] .* I(3)
    end
    blocks = spatial_span_executable_blocks(generator.span)
    push!(blocks, ExecutableEquationBlock(generator.name, :constraint,
        [position_prescription, velocity_prescription,
            acceleration_prescription], prescription!,
        prescription_jacobian!))
    push!(blocks, ExecutableEquationBlock(generator.name, :load,
        collect(generator.global_force_equations), load!, load_jacobian!))
    blocks
end

function equation_contributions(generator::SpatialRotationalMotionGenerator)
    marker_i = generator.hinge.marker_i
    marker_j = generator.hinge.marker_j
    rows = unique!([body_torque_rows(marker_i); body_torque_rows(marker_j)])
    residual! = function (equations, t, z, zdot)
        axis = marker_axis_kinematics(marker_j, z, 3)
        torque = z[generator.reaction_variable] .* axis.direction
        add_perp_reaction_to_body!(equations, z, marker_i, torque, 1)
        add_perp_reaction_to_body!(equations, z, marker_j, torque, -1)
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        if is_flexible_marker(marker_i) || is_flexible_marker(marker_j)
            columns = [generator.reaction_variable;
                       constraint_dependencies(marker_i, marker_j)]
            values = function (local_z)
                local_equations = zeros(eltype(local_z), maximum(rows))
                axis = marker_axis_kinematics(marker_j, local_z, 3)
                torque = local_z[generator.reaction_variable] .*
                    axis.direction
                add_perp_reaction_to_body!(local_equations, local_z,
                    marker_i, torque, 1)
                add_perp_reaction_to_body!(local_equations, local_z,
                    marker_j, torque, -1)
                local_equations[rows]
            end
            add_ad_jacobian!(jacobian, z, rows, columns, values)
            return nothing
        end
        axis = marker_axis_kinematics(marker_j, z, 3)
        normal_derivatives = Tuple{Any,Any}[]
        if !isnothing(axis.body)
            push!(normal_derivatives,
                (axis.body, axis.direction_parameters))
        end
        add_perp_reaction_jacobian!(jacobian, z, generator, marker_i, 1,
            axis.direction, normal_derivatives)
        add_perp_reaction_jacobian!(jacobian, z, generator, marker_j, -1,
            axis.direction, normal_derivatives)
    end
    EquationContribution[EquationContribution(generator.name,
        :torque_to_bodies, rows, residual!, jacobian!)]
end

function equation_contributions(generator::SpatialTranslationalMotionGenerator)
    rows = directed_distance_reaction_rows(generator.geometry)
    residual! = (equations, t, z, zdot) ->
        add_directed_distance_reaction!(equations, z, generator.geometry,
            generator.reaction_variable)
    jacobian! = (jacobian, t, z, zdot, coefficient) ->
        add_directed_distance_reaction_jacobian!(jacobian, z,
            generator.geometry, generator.reaction_variable)
    EquationContribution[EquationContribution(generator.name,
        :force_to_bodies, rows, residual!, jacobian!)]
end

function equation_contributions(generator::SpatialSpanningMotionGenerator)
    rows = unique!([body_force_rows(generator.span.marker_1);
                    body_force_rows(generator.span.marker_2)])
    residual! = function (equations, t, z, zdot)
        force = @view z[generator.global_force_variables]
        add_reaction_to_body!(equations, z, generator.span.marker_1, force, 1)
        add_reaction_to_body!(equations, z, generator.span.marker_2, force, -1)
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        force = @view z[generator.global_force_variables]
        add_reaction_jacobian!(jacobian, z, generator.span.marker_1, force, 1,
            generator.global_force_variables)
        add_reaction_jacobian!(jacobian, z, generator.span.marker_2, force, -1,
            generator.global_force_variables)
    end
    EquationContribution[EquationContribution(generator.name,
        :force_to_bodies, rows, residual!, jacobian!)]
end

end
