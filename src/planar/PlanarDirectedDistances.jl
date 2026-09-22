"""Shared directed-axis and point-to-plane kinematics for planar elements."""
module PlanarDirectedDistances

using LinearAlgebra
using ..PlanarAppliedForces

export PlanarDirectedAxis, PlanarDirectedDistance,
       directed_axis_values, directed_axis_dependencies,
       directed_distance_values, directed_distance_position,
       directed_distance_velocity, directed_distance_acceleration,
       directed_distance_dependencies

"""The oriented local y-axis of a planar marker and its owning body."""
struct PlanarDirectedAxis{B,O}
    body::B
    orientation::O
end

"""Signed position of one marker from the plane through another marker."""
struct PlanarDirectedDistance{BI,MI,BJ,MJ,A}
    body_i::BI
    marker_i::MI
    body_j::BJ
    marker_j::MJ
    axis::A
end

function directed_axis_values(axis::PlanarDirectedAxis, z)
    angle = PlanarAppliedForces.marker_angle(axis.orientation, z)
    omega = PlanarAppliedForces.marker_angular_velocity(axis.orientation, z)
    alpha = isnothing(axis.body) ? zero(eltype(z)) :
        axis.orientation isa PlanarFlexibleOrientationMarker ?
            PlanarAppliedForces.marker_angular_acceleration(
                axis.orientation, z) :
            z[axis.body.angular_acceleration_variable]
    unit = [-sin(angle), cos(angle)]
    transverse = [-unit[2], unit[1]]
    unit_rate = transverse .* omega
    unit_acceleration = transverse .* alpha .- unit .* omega^2
    (; unit, transverse, unit_rate, unit_acceleration, omega, alpha)
end

directed_axis_dependencies(axis::PlanarDirectedAxis) =
    isnothing(axis.body) ? Int[] :
        axis.orientation isa PlanarFlexibleOrientationMarker ?
            [axis.body.orientation_variable;
             collect(axis.body.elastic_position_variables)] :
            [axis.body.orientation_variable]

function marker_acceleration(body, marker, values, z)
    isnothing(body) && return zeros(eltype(z), 2)
    marker isa PlanarFloatingPointMarker && throw(ArgumentError(
        "a directed-distance acceleration cannot use a floating marker"))
    marker isa PlanarFlexiblePointMarker &&
        return PlanarAppliedForces.point_marker_acceleration(marker, z)
    z[body.acceleration_variables] .+
        values.d .* z[body.angular_acceleration_variable] .-
        values.r .* z[body.angular_velocity_variable]^2
end

function directed_distance_values(geometry::PlanarDirectedDistance, z;
        include_acceleration = true)
    marker_i = PlanarAppliedForces.point_marker_kinematics(
        geometry.marker_i, z)
    marker_j = PlanarAppliedForces.point_marker_kinematics(
        geometry.marker_j, z)
    axis = directed_axis_values(geometry.axis, z)
    separation = marker_i.position - marker_j.position
    relative_velocity = marker_i.velocity - marker_j.velocity
    relative_acceleration = include_acceleration ?
        marker_acceleration(geometry.body_i, geometry.marker_i, marker_i, z) -
        marker_acceleration(geometry.body_j, geometry.marker_j, marker_j, z) :
        nothing
    position = dot(separation, axis.unit)
    velocity = dot(relative_velocity, axis.unit) +
        dot(separation, axis.unit_rate)
    acceleration = include_acceleration ?
        dot(relative_acceleration, axis.unit) +
            2dot(relative_velocity, axis.unit_rate) +
            dot(separation, axis.unit_acceleration) : nothing
    (; marker_i, marker_j, axis, separation, relative_velocity,
       relative_acceleration, position, velocity, acceleration)
end

directed_distance_position(geometry::PlanarDirectedDistance, z) =
    directed_distance_values(geometry, z).position
directed_distance_velocity(geometry::PlanarDirectedDistance, z) =
    directed_distance_values(geometry, z).velocity
directed_distance_acceleration(geometry::PlanarDirectedDistance, z) =
    directed_distance_values(geometry, z).acceleration

function body_dependencies(body)
    isnothing(body) && return Int[]
    dependencies = [collect(body.acceleration_variables);
        body.angular_acceleration_variable;
        collect(body.velocity_variables);
        body.angular_velocity_variable;
        collect(body.position_variables);
        body.orientation_variable]
    hasproperty(body, :elastic_acceleration_variables) && append!(dependencies,
        [collect(body.elastic_acceleration_variables);
         collect(body.elastic_velocity_variables);
         collect(body.elastic_position_variables)])
    dependencies
end

directed_distance_dependencies(geometry::PlanarDirectedDistance) = unique([
    body_dependencies(geometry.body_i);
    body_dependencies(geometry.body_j)])

end
