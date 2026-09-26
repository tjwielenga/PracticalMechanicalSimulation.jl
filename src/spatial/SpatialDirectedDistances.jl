"""Shared directed-axis and point-to-plane kinematics for spatial elements."""
module SpatialDirectedDistances

using ForwardDiff
using LinearAlgebra
using StaticArrays: SVector
using ..AutomaticAnalysis
using ..SpatialComponentAssembly
using ..SpatialModeling

import ..SpatialComponentAssembly: component_registration,
    executable_blocks, equation_contributions

export SpatialDirectedAxis, SpatialDirectedDistance,
       SpatialDirectedDistanceMeasure,
       spatial_directed_distance_registration,
       allocated_spatial_directed_distance_measure,
       initialize_spatial_directed_distance_measure!,
       marker_axis_kinematics, marker_point_kinematics,
       directed_axis_values, directed_distance_values,
       directed_distance_direction, directed_distance_position,
       directed_distance_velocity, directed_distance_acceleration,
       directed_distance_jacobian!

"""One selected axis of an oriented spatial marker."""
struct SpatialDirectedAxis{M}
    marker::M
    index::Int

    function SpatialDirectedAxis(marker::M, index::Integer) where {M}
        index in 1:3 || throw(ArgumentError(
            "a spatial directed axis must be 1, 2, or 3"))
        new{M}(marker, Int(index))
    end
end

"""Signed distance from one marker point to a plane through another."""
struct SpatialDirectedDistance{A,B,D}
    marker_i::A
    marker_j::B
    axis::D
end

"""Reaction-free directed distance, velocity, and acceleration output."""
struct SpatialDirectedDistanceMeasure{G}
    name::Symbol
    geometry::G
    distance_variable::Int
    velocity_variable::Int
    acceleration_variable::Int
    position_equation::Int
    velocity_equation::Int
    acceleration_equation::Int
end

function spatial_directed_distance_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:distance, :relative_position, 0),
        VariableDeclaration(:velocity, :relative_velocity, 1),
        VariableDeclaration(:acceleration, :relative_acceleration, 2),
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:position, [EquationDeclaration(
            :distance_definition, :coordinate_relation, 0,
            :directed_distance)]),
        EquationBlockDeclaration(:velocity, [EquationDeclaration(
            :velocity_definition, :coordinate_relation, 1,
            :directed_distance)]),
        EquationBlockDeclaration(:acceleration, [EquationDeclaration(
            :acceleration_definition, :coordinate_relation, 2,
            :directed_distance)]),
    ]
    ComponentRegistration(name, variables, blocks)
end

component_registration(measure::SpatialDirectedDistanceMeasure) =
    spatial_directed_distance_registration(measure.name)

function allocated_spatial_directed_distance_measure(layout, name,
        marker_i, marker_j)
    variables = component_variable_indices(layout, name)
    geometry = SpatialDirectedDistance(
        marker_i, marker_j, SpatialDirectedAxis(marker_j, 3))
    SpatialDirectedDistanceMeasure(name, geometry,
        variables[1], variables[2], variables[3],
        only(component_equation_indices(layout, name, :position)),
        only(component_equation_indices(layout, name, :velocity)),
        only(component_equation_indices(layout, name, :acceleration)))
end

function marker_axis_kinematics(marker::SpatialGroundMarker, z, axis)
    direction = SVector{3}(@view marker.orientation[:, axis])
    zeros_vector = zero(SVector{3,eltype(z)})
    (; direction, velocity = zeros_vector,
       acceleration = zeros_vector, body = nothing,
       direction_parameters = nothing, velocity_parameters = nothing,
       velocity_omega = nothing, acceleration_parameters = nothing,
       acceleration_omega = nothing, acceleration_alpha = nothing)
end

function marker_axis_kinematics(marker::SpatialBodyMarker, z, axis)
    body = marker.body
    parameters = SVector{4}(@view z[body.euler_parameter_variables])
    orientation = rotation_matrix(parameters)
    local_direction = SVector{3}(@view marker.orientation_body[:, axis])
    omega = SVector{3}(@view z[body.angular_velocity_variables])
    alpha = SVector{3}(@view z[body.angular_acceleration_variables])
    local_velocity = cross(omega, local_direction)
    local_acceleration = cross(alpha, local_direction) .+
        cross(omega, local_velocity)
    direction = orientation * local_direction
    velocity = orientation * local_velocity
    acceleration = orientation * local_acceleration
    direction_parameters = rotation_vector_jacobian(
        parameters, local_direction)
    velocity_parameters = rotation_vector_jacobian(
        parameters, local_velocity)
    velocity_omega = -orientation * skew(local_direction)
    acceleration_parameters = rotation_vector_jacobian(
        parameters, local_acceleration)
    acceleration_omega = orientation *
        (-skew(local_velocity) - skew(omega) * skew(local_direction))
    acceleration_alpha = -orientation * skew(local_direction)
    (; direction, velocity, acceleration, body,
       direction_parameters, velocity_parameters, velocity_omega,
       acceleration_parameters, acceleration_omega,
       acceleration_alpha)
end

function marker_axis_kinematics(marker::SpatialFlexibleBeamMarker, z, axis)
    body = marker.body
    parameters = @view z[body.euler_parameter_variables]
    orientation = rotation_matrix(parameters)
    base_direction = @view marker.orientation_body[:, axis]
    elastic = @view z[body.elastic_position_variables]
    elastic_rate = @view z[body.elastic_velocity_variables]
    elastic_acceleration = @view z[body.elastic_acceleration_variables]
    phi = marker.orientation_shape * elastic
    phi_rate = marker.orientation_shape * elastic_rate
    phi_acceleration = marker.orientation_shape * elastic_acceleration
    local_direction = base_direction + cross(phi, base_direction)
    local_elastic_rate = cross(phi_rate, base_direction)
    local_velocity = cross(
        @view(z[body.angular_velocity_variables]), local_direction) +
        local_elastic_rate
    omega = @view z[body.angular_velocity_variables]
    alpha = @view z[body.angular_acceleration_variables]
    local_acceleration = cross(alpha, local_direction) +
        cross(omega, cross(omega, local_direction)) +
        2 .* cross(omega, local_elastic_rate) +
        cross(phi_acceleration, base_direction)
    direction = orientation * local_direction
    velocity = orientation * local_velocity
    acceleration = orientation * local_acceleration
    (; direction, velocity, acceleration, body,
       direction_parameters = nothing, velocity_parameters = nothing,
       velocity_omega = nothing, acceleration_parameters = nothing,
       acceleration_omega = nothing, acceleration_alpha = nothing)
end

function marker_point_kinematics(marker::SpatialGroundMarker, z)
    zero_vector = zero(SVector{3,eltype(z)})
    (; position = SVector{3}(marker.position), velocity = zero_vector,
       acceleration = zero_vector, body = nothing,
       position_parameters = nothing, velocity_parameters = nothing,
       velocity_omega = nothing, acceleration_parameters = nothing,
       acceleration_omega = nothing, acceleration_alpha = nothing)
end

function marker_point_kinematics(marker::SpatialBodyMarker, z)
    body = marker.body
    parameters = SVector{4}(@view z[body.euler_parameter_variables])
    orientation = rotation_matrix(parameters)
    offset = SVector{3}(marker.position_body)
    omega = SVector{3}(@view z[body.angular_velocity_variables])
    alpha = SVector{3}(@view z[body.angular_acceleration_variables])
    local_velocity = cross(omega, offset)
    local_acceleration = cross(alpha, offset) .+
        cross(omega, local_velocity)
    position = SVector{3}(@view z[body.position_variables]) .+
        orientation * offset
    velocity = SVector{3}(@view z[body.velocity_variables]) .+
        orientation * local_velocity
    acceleration = SVector{3}(@view z[body.acceleration_variables]) .+
        orientation * local_acceleration
    position_parameters = rotation_vector_jacobian(parameters, offset)
    velocity_parameters = rotation_vector_jacobian(
        parameters, local_velocity)
    velocity_omega = -orientation * skew(offset)
    acceleration_parameters = rotation_vector_jacobian(
        parameters, local_acceleration)
    acceleration_omega = orientation *
        (-skew(local_velocity) - skew(omega) * skew(offset))
    acceleration_alpha = -orientation * skew(offset)
    (; position, velocity, acceleration, body, position_parameters,
       velocity_parameters, velocity_omega, acceleration_parameters,
       acceleration_omega, acceleration_alpha)
end

function marker_point_kinematics(marker::SpatialFlexibleBeamMarker, z)
    body = marker.body
    position = spatial_marker_position(marker, z)
    velocity = spatial_marker_velocity(marker, z)
    acceleration = spatial_marker_acceleration(marker, z)
    (; position, velocity, acceleration, body,
       position_parameters = nothing, velocity_parameters = nothing,
       velocity_omega = nothing, acceleration_parameters = nothing,
       acceleration_omega = nothing, acceleration_alpha = nothing)
end

directed_axis_values(axis::SpatialDirectedAxis, z) =
    marker_axis_kinematics(axis.marker, z, axis.index)

function directed_distance_values(geometry::SpatialDirectedDistance, z)
    first = marker_point_kinematics(geometry.marker_i, z)
    second = marker_point_kinematics(geometry.marker_j, z)
    direction = directed_axis_values(geometry.axis, z)
    separation = first.position - second.position
    relative_velocity = first.velocity - second.velocity
    relative_acceleration = first.acceleration - second.acceleration
    position = dot(separation, direction.direction)
    velocity = dot(relative_velocity, direction.direction) +
        dot(separation, direction.velocity)
    acceleration = dot(relative_acceleration, direction.direction) +
        2dot(relative_velocity, direction.velocity) +
        dot(separation, direction.acceleration)
    (; first, second, direction, separation, relative_velocity,
       relative_acceleration, position, velocity, acceleration)
end

directed_distance_direction(geometry::SpatialDirectedDistance, z) =
    directed_axis_values(geometry.axis, z).direction
directed_distance_position(geometry::SpatialDirectedDistance, z) =
    directed_distance_values(geometry, z).position
directed_distance_velocity(geometry::SpatialDirectedDistance, z) =
    directed_distance_values(geometry, z).velocity
directed_distance_acceleration(geometry::SpatialDirectedDistance, z) =
    directed_distance_values(geometry, z).acceleration

function add_point_jacobian!(jacobian, values, point, sign,
        acceleration_row, velocity_row, position_row, multiplier)
    isnothing(point.body) && return nothing
    body = point.body
    direction = values.direction
    signed = multiplier * sign
    jacobian[position_row, body.position_variables] .+=
        signed .* direction.direction
    jacobian[velocity_row, body.position_variables] .+=
        signed .* direction.velocity
    jacobian[acceleration_row, body.position_variables] .+=
        signed .* direction.acceleration
    jacobian[velocity_row, body.velocity_variables] .+=
        signed .* direction.direction
    jacobian[acceleration_row, body.velocity_variables] .+=
        2signed .* direction.velocity
    jacobian[acceleration_row, body.acceleration_variables] .+=
        signed .* direction.direction

    jacobian[position_row, body.euler_parameter_variables] .+= signed .* (
        transpose(point.position_parameters) * direction.direction)
    jacobian[velocity_row, body.euler_parameter_variables] .+= signed .* (
        transpose(point.velocity_parameters) * direction.direction .+
        transpose(point.position_parameters) * direction.velocity)
    jacobian[acceleration_row, body.euler_parameter_variables] .+= signed .* (
        transpose(point.acceleration_parameters) * direction.direction .+
        2transpose(point.velocity_parameters) * direction.velocity .+
        transpose(point.position_parameters) * direction.acceleration)
    jacobian[velocity_row, body.angular_velocity_variables] .+= signed .* (
        transpose(point.velocity_omega) * direction.direction)
    jacobian[acceleration_row, body.angular_velocity_variables] .+= signed .* (
        transpose(point.acceleration_omega) * direction.direction .+
        2transpose(point.velocity_omega) * direction.velocity)
    jacobian[acceleration_row,
        body.angular_acceleration_variables] .+= signed .* (
        transpose(point.acceleration_alpha) * direction.direction)
    nothing
end

function add_axis_jacobian!(jacobian, values, acceleration_row,
        velocity_row, position_row, multiplier)
    direction = values.direction
    isnothing(direction.body) && return nothing
    body = direction.body
    separation = values.separation
    relative_velocity = values.relative_velocity
    relative_acceleration = values.relative_acceleration
    jacobian[position_row, body.euler_parameter_variables] .+= multiplier .* (
        transpose(direction.direction_parameters) * separation)
    jacobian[velocity_row, body.euler_parameter_variables] .+= multiplier .* (
        transpose(direction.direction_parameters) * relative_velocity .+
        transpose(direction.velocity_parameters) * separation)
    jacobian[acceleration_row, body.euler_parameter_variables] .+= multiplier .* (
        transpose(direction.direction_parameters) * relative_acceleration .+
        2transpose(direction.velocity_parameters) * relative_velocity .+
        transpose(direction.acceleration_parameters) * separation)
    jacobian[velocity_row, body.angular_velocity_variables] .+= multiplier .* (
        transpose(direction.velocity_omega) * separation)
    jacobian[acceleration_row, body.angular_velocity_variables] .+= multiplier .* (
        2transpose(direction.velocity_omega) * relative_velocity .+
        transpose(direction.acceleration_omega) * separation)
    jacobian[acceleration_row,
        body.angular_acceleration_variables] .+= multiplier .* (
        transpose(direction.acceleration_alpha) * separation)
    nothing
end

"""Add the analytical Jacobian of all three directed-distance levels."""
function directed_distance_jacobian!(jacobian, z,
        geometry::SpatialDirectedDistance, acceleration_row, velocity_row,
        position_row, multiplier = 1)
    markers = (geometry.marker_i, geometry.marker_j, geometry.axis.marker)
    if any(is_flexible_marker, markers)
        columns = sort!(unique!(reduce(vcat,
            (spatial_marker_dependency_indices(marker) for marker in markers);
            init = Int[])))
        rows = [acceleration_row, velocity_row, position_row]
        initial = collect(z[columns])
        block = ForwardDiff.jacobian(initial) do local_values
            local_z = z .+ zero(eltype(local_values))
            local_z[columns] .= local_values
            distance_values = directed_distance_values(geometry, local_z)
            multiplier .* [distance_values.acceleration,
                distance_values.velocity, distance_values.position]
        end
        jacobian[rows, columns] .+= block
        return nothing
    end
    values = directed_distance_values(geometry, z)
    add_point_jacobian!(jacobian, values, values.first, 1,
        acceleration_row, velocity_row, position_row, multiplier)
    add_point_jacobian!(jacobian, values, values.second, -1,
        acceleration_row, velocity_row, position_row, multiplier)
    add_axis_jacobian!(jacobian, values, acceleration_row, velocity_row,
        position_row, multiplier)
    nothing
end

function initialize_spatial_directed_distance_measure!(initial,
        measure::SpatialDirectedDistanceMeasure)
    values = directed_distance_values(measure.geometry, initial)
    initial[measure.distance_variable] = values.position
    initial[measure.velocity_variable] = values.velocity
    initial[measure.acceleration_variable] = values.acceleration
    initial
end

function executable_blocks(measure::SpatialDirectedDistanceMeasure)
    residual! = function (equations, t, z, zdot)
        values = directed_distance_values(measure.geometry, z)
        equations[measure.position_equation] =
            z[measure.distance_variable] - values.position
        equations[measure.velocity_equation] =
            z[measure.velocity_variable] - values.velocity
        equations[measure.acceleration_equation] =
            z[measure.acceleration_variable] - values.acceleration
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        jacobian[measure.position_equation, measure.distance_variable] += 1
        jacobian[measure.velocity_equation, measure.velocity_variable] += 1
        jacobian[measure.acceleration_equation,
            measure.acceleration_variable] += 1
        directed_distance_jacobian!(jacobian, z, measure.geometry,
            measure.acceleration_equation, measure.velocity_equation,
            measure.position_equation, -1)
    end
    rows = [measure.acceleration_equation, measure.velocity_equation,
            measure.position_equation]
    ExecutableEquationBlock[ExecutableEquationBlock(measure.name,
        :measurement, rows, residual!, jacobian!)]
end

equation_contributions(::SpatialDirectedDistanceMeasure) =
    EquationContribution[]

end
