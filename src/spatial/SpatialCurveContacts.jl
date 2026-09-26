"""Marker-fixed extruded cam profiles with roller and flat followers."""
module SpatialCurveContacts

using LinearAlgebra
using ..AutomaticAnalysis
using ..PlanarCurveContacts: PlanarClosedCurve, curve_point
using ..SpatialComponentAssembly
using ..SpatialModeling

import ..SpatialComponentAssembly: component_registration,
    executable_blocks, equation_contributions

export SpatialCurveDefinition, SpatialCurveContactComponent,
       spatial_curve_contact_registration, allocated_spatial_curve_contact,
       initialize_spatial_curve_contact!, set_spatial_curve_contact_stage!,
       spatial_curve_contact_values, spatial_curve_contact_gap,
       spatial_curve_contact_damping_surface

"""A closed planar profile fixed to a spatial marker and extruded along local z."""
struct SpatialCurveDefinition{C,M,T}
    name::Symbol
    profile::C
    marker::M
    half_width::T
end

"""One-sided roller or flat contact with an extruded spatial cam profile."""
struct SpatialCurveContactComponent{C,M,F,L,T,S}
    name::Symbol
    curve::C
    curve_marker::M
    follower_marker::F
    follower_kind::Symbol
    law::L
    expression::Bool
    radius::T
    curve_half_width::T
    stiffness::T
    damping_factor::T
    transition_depth::T
    normal_side::Int
    axis_tolerance::T
    active_during::S
    active::Base.RefValue{Bool}
    station_variable::Int
    station_rate_variable::Int
    curvature_variable::Int
    gap_variable::Int
    gap_rate_variable::Int
    normal_force_variable::Int
    contact_point_variables::UnitRange{Int}
    normal_variables::UnitRange{Int}
    global_force_variables::UnitRange{Int}
    contact_equations::UnitRange{Int}
end

function spatial_curve_contact_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:station, :applied_geometry, 0),
        VariableDeclaration(:station_rate, :applied_rate, 1),
        VariableDeclaration(:curvature, :applied_geometry, 0),
        VariableDeclaration(:gap, :applied_geometry, 0),
        VariableDeclaration(:gap_rate, :applied_rate, 1),
        VariableDeclaration(:normal_force, :applied_load, 2),
        [VariableDeclaration(Symbol(:contact_, axis), :applied_geometry, 0)
            for axis in (:x, :y, :z)]...,
        [VariableDeclaration(Symbol(:normal_, axis), :applied_geometry, 0)
            for axis in (:x, :y, :z)]...,
        [VariableDeclaration(Symbol(:F_, axis), :applied_load, 2)
            for axis in (:x, :y, :z)]...,
    ]
    equations = EquationDeclaration[
        EquationDeclaration(:tangent, :applied_definition, 0,
            :contact_station),
        EquationDeclaration(:tangent_rate, :applied_definition, 1,
            :contact_station_rate),
        EquationDeclaration(:curvature, :applied_definition, 0,
            :contact_curvature),
        EquationDeclaration(:gap, :applied_definition, 0, :contact_gap),
        EquationDeclaration(:gap_rate, :applied_definition, 1,
            :contact_gap_rate),
        EquationDeclaration(:normal_force, :applied_definition, 2,
            :contact_force),
        [EquationDeclaration(Symbol(:contact_, axis), :applied_definition, 0,
            :contact_point) for axis in (:x, :y, :z)]...,
        [EquationDeclaration(Symbol(:normal_, axis), :applied_definition, 0,
            :contact_normal) for axis in (:x, :y, :z)]...,
        [EquationDeclaration(Symbol(:global_force_, axis),
            :applied_definition, 2, :global_force)
            for axis in (:x, :y, :z)]...,
    ]
    ComponentRegistration(name, variables,
        [EquationBlockDeclaration(:contact, equations)])
end

component_registration(contact::SpatialCurveContactComponent) =
    spatial_curve_contact_registration(contact.name)

function allocated_spatial_curve_contact(layout, name, curve,
        curve_marker, follower_marker, radius, curve_half_width,
        stiffness, damping_factor;
        law = nothing, expression = false, transition_depth = 0.0,
        side = :outside, follower_kind = :roller,
        axis_tolerance = 1.0e-6,
        active_during = (:static, :dynamic, :modal))
    follower_kind in (:roller, :flat) || throw(ArgumentError(
        "curve contact '$name' follower_kind must be roller or flat"))
    side in (:outside, :inside) || throw(ArgumentError(
        "curve contact '$name' side must be 'outside' or 'inside'"))
    variables = component_variable_indices(layout, name)
    SpatialCurveContactComponent(name, curve, curve_marker, follower_marker,
        follower_kind, law, Bool(expression), Float64(radius),
        Float64(curve_half_width), Float64(stiffness), Float64(damping_factor),
        Float64(transition_depth), side == :outside ? 1 : -1,
        Float64(axis_tolerance), active_during,
        Ref(:dynamic in active_during), variables[1], variables[2],
        variables[3], variables[4], variables[5], variables[6],
        variables[7:9], variables[10:12], variables[13:15],
        component_equation_indices(layout, name, :contact))
end

function set_spatial_curve_contact_stage!(contact::SpatialCurveContactComponent,
        stage)
    contact.active[] = stage in contact.active_during
    contact
end

function marker_angular_velocity(::SpatialGroundMarker, z)
    zeros(eltype(z), 3)
end

function marker_angular_velocity(marker::SpatialBodyMarker, z)
    body = marker.body
    rotation_matrix(@view(z[body.euler_parameter_variables])) *
        @view(z[body.angular_velocity_variables])
end

marker_angular_velocity(marker::SpatialFloatingMarker, z) =
    marker_angular_velocity(marker.follower, z)

function effective_contact_penetration(contact, gap)
    penetration = max(-gap, zero(gap))
    depth = contact.transition_depth
    if !(depth > 0) || penetration >= depth
        return penetration - (depth > 0 ? depth / 2 : zero(depth))
    end
    x = penetration / depth
    depth * x^6 * (21 + x * (-60 + x * (67.5 + x * (-35 + 7x))))
end

function spatial_curve_contact_kinematics(
        contact::SpatialCurveContactComponent, z)
    station = z[contact.station_variable]
    station_rate = z[contact.station_rate_variable]
    local_values = curve_point(contact.curve, station)

    curve_origin = spatial_marker_position(contact.curve_marker, z)
    curve_origin_velocity = spatial_marker_velocity(contact.curve_marker, z)
    curve_orientation = spatial_marker_orientation(contact.curve_marker, z)
    curve_omega = marker_angular_velocity(contact.curve_marker, z)
    follower_position = spatial_marker_position(contact.follower_marker, z)
    follower_velocity = spatial_marker_velocity(contact.follower_marker, z)
    follower_orientation = spatial_marker_orientation(
        contact.follower_marker, z)
    follower_omega = marker_angular_velocity(contact.follower_marker, z)

    local_point = [local_values.position; zero(station)]
    local_first = [local_values.derivative; zero(station)]
    local_second = [local_values.second_derivative; zero(station)]
    point_offset = curve_orientation * local_point
    base_point = curve_origin + point_offset
    first = curve_orientation * local_first
    second = curve_orientation * local_second
    speed = norm(first)
    speed > sqrt(eps(real(float(one(speed))))) || throw(DomainError(
        speed, "spatial curve-contact tangent is undefined"))
    tangent = first ./ speed
    tangent_station = (second .- tangent .* dot(tangent, second)) ./ speed
    normal_sign = contact.normal_side * contact.curve.outward_sign
    local_tangent = local_values.derivative ./ norm(local_values.derivative)
    local_curve_normal = normal_sign .* [local_tangent[2],
        -local_tangent[1], zero(station)]
    curve_normal = curve_orientation * local_curve_normal
    local_tangent_station = (local_values.second_derivative .-
        local_tangent .* dot(local_tangent, local_values.second_derivative)) ./
        norm(local_values.derivative)
    local_normal_station = normal_sign .* [local_tangent_station[2],
        -local_tangent_station[1], zero(station)]
    curve_normal_station = curve_orientation * local_normal_station

    axis = curve_orientation[:, 3]
    axis_rate = cross(curve_omega, axis)
    base_material_velocity = curve_origin_velocity +
        cross(curve_omega, point_offset)
    moving_base_velocity = base_material_velocity + first .* station_rate
    axial = dot(follower_position - base_point, axis)
    axial_rate = dot(follower_velocity - moving_base_velocity, axis) +
        dot(follower_position - base_point, axis_rate)
    contact_point = base_point + axial .* axis
    contact_velocity = moving_base_velocity + axial_rate .* axis +
        axial .* axis_rate
    separation = follower_position - contact_point
    separation_rate = follower_velocity - contact_velocity
    tangent_rate = cross(curve_omega, tangent) +
        tangent_station .* station_rate

    normal, normal_rate, tangent_error, tangent_error_rate, gap =
        if contact.follower_kind == :roller
            curve_normal_rate = cross(curve_omega, curve_normal) +
                curve_normal_station .* station_rate
            (curve_normal, curve_normal_rate,
             dot(separation, tangent),
             dot(separation_rate, tangent) + dot(separation, tangent_rate),
             dot(separation, curve_normal) - contact.radius)
        else
            follower_normal = follower_orientation[:, 2]
            follower_normal_rate = cross(follower_omega, follower_normal)
            (follower_normal, follower_normal_rate,
             dot(tangent, follower_normal),
             dot(tangent_rate, follower_normal) +
                dot(tangent, follower_normal_rate),
             dot(separation, follower_normal))
        end
    gap_rate = dot(separation_rate, normal) + dot(separation, normal_rate)
    curvature = dot(tangent_station ./ speed, curve_normal)
    penetration = max(-gap, zero(gap))
    effective_penetration = effective_contact_penetration(contact, gap)
    damping_multiplier = max(zero(gap),
        one(gap) - contact.damping_factor * gap_rate)
    (; curve_origin, follower_position, contact_point, tangent, normal,
       tangent_error, tangent_error_rate, curvature, gap, gap_rate,
       penetration, effective_penetration, damping_multiplier,
       axis_alignment = dot(axis, follower_orientation[:, 3]))
end

function calculated_contact_force(contact, time, z, kinematics)
    contact.active[] || return zero(kinematics.gap)
    contact.expression && return contact.law(time, z)
    contact.stiffness * kinematics.effective_penetration *
        kinematics.damping_multiplier
end

function spatial_curve_contact_values(contact::SpatialCurveContactComponent,
        z, time = 0.0)
    kinematics = spatial_curve_contact_kinematics(contact, z)
    normal_force = calculated_contact_force(contact, time, z, kinematics)
    global_force = normal_force .* kinematics.normal
    (; kinematics..., normal_force, global_force)
end

spatial_curve_contact_gap(contact::SpatialCurveContactComponent, z) =
    z[contact.gap_variable]

function spatial_curve_contact_damping_surface(
        contact::SpatialCurveContactComponent, z)
    contact.expression && return one(eltype(z))
    gap = z[contact.gap_variable]
    gap < 0 ? one(gap) -
        contact.damping_factor * z[contact.gap_rate_variable] : one(gap)
end

function contact_kinematic_dependencies(marker)
    marker isa SpatialGroundMarker && return Int[]
    body = marker.body
    [collect(body.position_variables);
     collect(body.velocity_variables);
     collect(body.angular_velocity_variables);
     collect(body.euler_parameter_variables)]
end

function contact_configuration_dependencies(marker)
    marker isa SpatialGroundMarker && return Int[]
    body = marker.body
    [collect(body.position_variables); collect(body.euler_parameter_variables)]
end

function contact_residual(contact, time, z)
    values = spatial_curve_contact_kinematics(contact, z)
    force = calculated_contact_force(contact, time, z, values)
    [values.tangent_error;
     values.tangent_error_rate;
     z[contact.curvature_variable] - values.curvature;
     z[contact.gap_variable] - values.gap;
     z[contact.gap_rate_variable] - values.gap_rate;
     z[contact.normal_force_variable] - force;
     z[contact.contact_point_variables] .- values.contact_point;
     z[contact.normal_variables] .- values.normal;
     z[contact.global_force_variables] .-
        z[contact.normal_force_variable] .* values.normal]
end

function executable_blocks(contact::SpatialCurveContactComponent)
    dependencies = sort!(unique!([
        contact_kinematic_dependencies(contact.curve_marker);
        contact_kinematic_dependencies(contact.follower_marker);
        contact.station_variable; contact.station_rate_variable;
        contact.curvature_variable; contact.gap_variable;
        contact.gap_rate_variable; contact.normal_force_variable;
        collect(contact.contact_point_variables);
        collect(contact.normal_variables);
        collect(contact.global_force_variables);
        contact.expression ? contact.law.dependencies : Int[]]))
    residual! = function (equations, t, z, zdot)
        equations[contact.contact_equations] .= contact_residual(contact, t, z)
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        partials = spatial_local_state_jacobian(
            state -> contact_residual(contact, t, state), z, dependencies)
        jacobian[contact.contact_equations, dependencies] .+= partials
    end
    [ExecutableEquationBlock(contact.name, :contact,
        collect(contact.contact_equations), residual!, jacobian!)]
end

function body_force_contribution(body, z, point, global_force)
    orientation = rotation_matrix(@view(z[body.euler_parameter_variables]))
    lever_body = transpose(orientation) *
        (point - z[body.position_variables])
    force_body = transpose(orientation) * global_force
    [-global_force; -cross(lever_body, force_body)]
end

function marker_body(marker)
    marker isa SpatialGroundMarker && return nothing
    marker.body
end

function contact_body_contribution(contact, z)
    values = spatial_curve_contact_kinematics(contact, z)
    global_force = @view z[contact.global_force_variables]
    contributions = eltype(z)[]
    follower_body = marker_body(contact.follower_marker)
    isnothing(follower_body) || append!(contributions,
        body_force_contribution(follower_body, z, values.contact_point,
            global_force))
    curve_body = marker_body(contact.curve_marker)
    isnothing(curve_body) || append!(contributions,
        body_force_contribution(curve_body, z, values.contact_point,
            -global_force))
    contributions
end

function equation_contributions(contact::SpatialCurveContactComponent)
    rows = Int[]
    follower_body = marker_body(contact.follower_marker)
    curve_body = marker_body(contact.curve_marker)
    isnothing(follower_body) || append!(rows, follower_body.balance_equations)
    isnothing(curve_body) || append!(rows, curve_body.balance_equations)
    dependencies = sort!(unique!([
        contact_configuration_dependencies(contact.curve_marker);
        contact_configuration_dependencies(contact.follower_marker);
        contact.station_variable;
        collect(contact.global_force_variables)]))
    residual! = function (equations, t, z, zdot)
        equations[rows] .+= contact_body_contribution(contact, z)
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        partials = spatial_local_state_jacobian(
            state -> contact_body_contribution(contact, state), z,
            dependencies)
        jacobian[rows, dependencies] .+= partials
    end
    [EquationContribution(contact.name, :load_to_bodies, rows,
        residual!, jacobian!)]
end

function closest_station(contact, state)
    follower = spatial_marker_position(contact.follower_marker, state)
    origin = spatial_marker_position(contact.curve_marker, state)
    orientation = spatial_marker_orientation(contact.curve_marker, state)
    samples = max(64, 8 * size(contact.curve.points, 2))
    stations = range(0.0, contact.curve.length; length = samples + 1)[1:end-1]
    station = if contact.follower_kind == :roller
        distances = [begin
            local_point = curve_point(contact.curve, candidate).position
            offset = follower - (origin + orientation * [local_point; 0.0])
            in_plane = offset - dot(offset, orientation[:, 3]) .* orientation[:, 3]
            sum(abs2, in_plane)
        end for candidate in stations]
        stations[argmin(distances)]
    else
        follower_normal = spatial_marker_orientation(
            contact.follower_marker, state)[:, 2]
        projections = [begin
            local_point = curve_point(contact.curve, candidate).position
            dot(origin + orientation * [local_point; 0.0], follower_normal)
        end for candidate in stations]
        stations[argmax(projections)]
    end
    for _ in 1:12
        values = curve_point(contact.curve, station)
        point = origin + orientation * [values.position; 0.0]
        first = orientation * [values.derivative; 0.0]
        second = orientation * [values.second_derivative; 0.0]
        residual, derivative = if contact.follower_kind == :roller
            separation = follower - point
            (dot(separation, first),
             -dot(first, first) + dot(separation, second))
        else
            follower_normal = spatial_marker_orientation(
                contact.follower_marker, state)[:, 2]
            (dot(first, follower_normal), dot(second, follower_normal))
        end
        abs(derivative) > eps(Float64) || break
        increment = clamp(residual / derivative,
            -contact.curve.length / 8, contact.curve.length / 8)
        station -= increment
        abs(increment) <= 1.0e-12 * max(contact.curve.length, 1.0) && break
    end
    mod(station, contact.curve.length)
end

function initialize_spatial_curve_contact!(state,
        contact::SpatialCurveContactComponent, time = 0.0;
        station = nothing)
    curve_axis = spatial_marker_orientation(contact.curve_marker, state)[:, 3]
    follower_axis = spatial_marker_orientation(
        contact.follower_marker, state)[:, 3]
    misalignment = norm(cross(curve_axis, follower_axis))
    misalignment <= contact.axis_tolerance || throw(ArgumentError(
        "curve contact '$(contact.name)' requires parallel cam and follower " *
        "marker z axes; initial sine misalignment is $misalignment"))
    state[contact.station_variable] = isnothing(station) ?
        closest_station(contact, state) : Float64(station)
    state[contact.station_rate_variable] = 0.0
    base = spatial_curve_contact_kinematics(contact, state)
    trial = copy(state)
    trial[contact.station_rate_variable] = 1.0
    unit = spatial_curve_contact_kinematics(contact, trial)
    coefficient = unit.tangent_error_rate - base.tangent_error_rate
    abs(coefficient) > eps(Float64) &&
        (state[contact.station_rate_variable] =
            -base.tangent_error_rate / coefficient)
    values = spatial_curve_contact_kinematics(contact, state)
    state[contact.curvature_variable] = values.curvature
    state[contact.gap_variable] = values.gap
    state[contact.gap_rate_variable] = values.gap_rate
    state[contact.contact_point_variables] .= values.contact_point
    state[contact.normal_variables] .= values.normal
    force = calculated_contact_force(contact, time, state, values)
    isfinite(force) || throw(ArgumentError(
        "curve contact '$(contact.name)' force is not finite initially"))
    state[contact.normal_force_variable] = force
    state[contact.global_force_variables] .= force .* values.normal
    state
end

end
