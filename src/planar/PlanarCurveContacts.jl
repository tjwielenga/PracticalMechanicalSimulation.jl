"""Smooth closed profiles and roller contact for planar models."""
module PlanarCurveContacts

using ForwardDiff
using LinearAlgebra
using ..AutomaticAnalysis
using ..PlanarAppliedForces
using ..PlanarComponentAssembly

import ..PlanarComponentAssembly: component_registration,
    executable_blocks, equation_contributions

export PlanarClosedCurve, PlanarCurveContactComponent,
       planar_curve_contact_registration, allocated_planar_curve_contact,
       initialize_planar_curve_contact!, set_planar_curve_contact_stage!,
       curve_point, curve_contact_values, curve_contact_gap,
       curve_contact_damping_surface

"""
A periodic, twice-continuously-differentiable cubic spline in a marker frame.

`station` is chord length around the supplied control polygon. Evaluation is
periodic, so an integrated or predicted station may pass through either end
without a discontinuity.
"""
struct PlanarClosedCurve{T}
    points::Matrix{T}
    knots::Vector{T}
    second_derivatives::Matrix{T}
    length::T
    outward_sign::T
end

function PlanarClosedCurve(raw_points::AbstractVector)
    length(raw_points) >= 4 || throw(ArgumentError(
        "a closed planar curve requires at least four points"))
    points = [Float64.(point) for point in raw_points]
    all(point -> length(point) == 2 && all(isfinite, point), points) ||
        throw(ArgumentError("planar curve points must be finite [x, y] pairs"))
    if isapprox(points[1], points[end]; rtol = 0.0, atol = 1.0e-14)
        pop!(points)
    end
    length(points) >= 4 || throw(ArgumentError(
        "a closed planar curve requires at least four distinct points"))
    count = length(points)
    matrix = reduce(hcat, points)
    intervals = [norm(points[mod1(index + 1, count)] - points[index])
                 for index in 1:count]
    all(>(sqrt(eps(Float64))), intervals) || throw(ArgumentError(
        "consecutive planar curve points must be distinct"))
    knots = [0.0; cumsum(intervals)]
    total = knots[end]

    system = zeros(Float64, count, count)
    right = zeros(Float64, count, 2)
    for index in 1:count
        previous = mod1(index - 1, count)
        following = mod1(index + 1, count)
        h_previous = intervals[previous]
        h_next = intervals[index]
        system[index, previous] = h_previous
        system[index, index] = 2(h_previous + h_next)
        system[index, following] = h_next
        right[index, :] .= 6 .* ((points[following] - points[index]) ./ h_next .-
            (points[index] - points[previous]) ./ h_previous)
    end
    second = Matrix(transpose(system \ right))
    twice_area = sum(points[index][1] * points[mod1(index + 1, count)][2] -
        points[mod1(index + 1, count)][1] * points[index][2]
        for index in 1:count)
    abs(twice_area) > sqrt(eps(Float64)) || throw(ArgumentError(
        "planar curve control polygon must enclose a nonzero area"))
    PlanarClosedCurve(matrix, knots, second, total, sign(twice_area))
end

function curve_point(curve::PlanarClosedCurve, station)
    wrapped = mod(station, curve.length)
    index = min(searchsortedlast(curve.knots, wrapped),
        size(curve.points, 2))
    next = mod1(index + 1, size(curve.points, 2))
    left = curve.knots[index]
    right = curve.knots[index + 1]
    h = right - left
    a = (right - wrapped) / h
    b = (wrapped - left) / h
    first_point = @view curve.points[:, index]
    second_point = @view curve.points[:, next]
    first_second = @view curve.second_derivatives[:, index]
    second_second = @view curve.second_derivatives[:, next]
    position = a .* first_point .+ b .* second_point .+
        ((a^3 - a) .* first_second .+ (b^3 - b) .* second_second) .* (h^2 / 6)
    derivative = (second_point - first_point) ./ h .+
        ((-3a^2 + 1) .* first_second .+ (3b^2 - 1) .* second_second) .* (h / 6)
    second_derivative = a .* first_second .+ b .* second_second
    (; position, derivative, second_derivative)
end

"""Circular roller in one-sided contact with a closed marker-fixed curve."""
struct PlanarCurveContactComponent{C,M,R,L,T,S}
    name::Symbol
    curve::C
    curve_marker::M
    roller_marker::R
    law::L
    expression::Bool
    radius::T
    stiffness::T
    damping_factor::T
    transition_depth::T
    normal_side::Int
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

function planar_curve_contact_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:station, :applied_geometry, 0),
        VariableDeclaration(:station_rate, :applied_rate, 1),
        VariableDeclaration(:curvature, :applied_geometry, 0),
        VariableDeclaration(:gap, :applied_geometry, 0),
        VariableDeclaration(:gap_rate, :applied_rate, 1),
        VariableDeclaration(:normal_force, :applied_load, 2),
        VariableDeclaration(:contact_x, :applied_geometry, 0),
        VariableDeclaration(:contact_y, :applied_geometry, 0),
        VariableDeclaration(:normal_x, :applied_geometry, 0),
        VariableDeclaration(:normal_y, :applied_geometry, 0),
        VariableDeclaration(:F_x, :applied_load, 2),
        VariableDeclaration(:F_y, :applied_load, 2),
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
        EquationDeclaration(:contact_x, :applied_definition, 0,
            :contact_point),
        EquationDeclaration(:contact_y, :applied_definition, 0,
            :contact_point),
        EquationDeclaration(:normal_x, :applied_definition, 0,
            :contact_normal),
        EquationDeclaration(:normal_y, :applied_definition, 0,
            :contact_normal),
        EquationDeclaration(:global_force_x, :applied_definition, 2,
            :global_force),
        EquationDeclaration(:global_force_y, :applied_definition, 2,
            :global_force),
    ]
    ComponentRegistration(name, variables,
        [EquationBlockDeclaration(:contact, equations)])
end

component_registration(contact::PlanarCurveContactComponent) =
    planar_curve_contact_registration(contact.name)

function allocated_planar_curve_contact(layout, name, curve, curve_marker,
        roller_marker, radius, stiffness, damping_factor;
        law = nothing, expression = false, transition_depth = 0.0,
        side = :outside,
        active_during = (:static, :dynamic, :modal))
    side in (:outside, :inside) || throw(ArgumentError(
        "curve contact '$name' side must be 'outside' or 'inside'"))
    variables = component_variable_indices(layout, name)
    PlanarCurveContactComponent(name, curve, curve_marker, roller_marker,
        law, Bool(expression), Float64(radius), Float64(stiffness),
        Float64(damping_factor), Float64(transition_depth),
        side == :outside ? 1 : -1, active_during,
        Ref(:dynamic in active_during), variables[1], variables[2],
        variables[3], variables[4], variables[5], variables[6],
        variables[7:8], variables[9:10], variables[11:12],
        component_equation_indices(layout, name, :contact))
end

function set_planar_curve_contact_stage!(contact::PlanarCurveContactComponent,
        stage)
    contact.active[] = stage in contact.active_during
    contact
end

perpendicular(vector) = [-vector[2], vector[1]]

function effective_contact_penetration(contact, gap)
    penetration = max(-gap, zero(gap))
    depth = contact.transition_depth
    if !(depth > 0) || penetration >= depth
        return penetration - (depth > 0 ? depth / 2 : zero(depth))
    end
    x = penetration / depth
    depth * x^6 * (21 + x * (-60 + x * (67.5 + x * (-35 + 7x))))
end

function curve_contact_kinematics(contact::PlanarCurveContactComponent, z)
    station = z[contact.station_variable]
    station_rate = z[contact.station_rate_variable]
    local_values = curve_point(contact.curve, station)
    curve_origin = PlanarAppliedForces.point_marker_kinematics(
        contact.curve_marker.point, z)
    roller = PlanarAppliedForces.point_marker_kinematics(
        contact.roller_marker.point, z)
    angle = PlanarAppliedForces.marker_angle(contact.curve_marker.orientation, z)
    omega = PlanarAppliedForces.marker_angular_velocity(
        contact.curve_marker.orientation, z)
    cosine, sine = cos(angle), sin(angle)
    rotation = [cosine -sine; sine cosine]
    local_point = local_values.position
    point_offset = rotation * local_point
    first = rotation * local_values.derivative
    second = rotation * local_values.second_derivative
    speed = norm(first)
    speed > sqrt(eps(real(float(one(speed))))) || throw(DomainError(
        speed, "curve contact tangent is undefined"))
    tangent = first ./ speed
    tangent_station = (second .- tangent .* dot(tangent, second)) ./ speed
    normal_sign = contact.normal_side * contact.curve.outward_sign
    normal = normal_sign .* [tangent[2], -tangent[1]]
    normal_station = normal_sign .* [tangent_station[2], -tangent_station[1]]
    contact_point = curve_origin.position + point_offset
    point_material_velocity = curve_origin.velocity +
        omega .* perpendicular(point_offset)
    separation = roller.position - contact_point
    moving_point_velocity = point_material_velocity + first .* station_rate
    separation_rate = roller.velocity - moving_point_velocity
    tangent_rate = omega .* perpendicular(tangent) .+
        tangent_station .* station_rate
    normal_rate = omega .* perpendicular(normal) .+
        normal_station .* station_rate
    tangent_error = dot(separation, tangent)
    tangent_error_rate = dot(separation_rate, tangent) +
        dot(separation, tangent_rate)
    gap = dot(separation, normal) - contact.radius
    gap_rate = dot(separation_rate, normal) + dot(separation, normal_rate)
    curvature = dot(tangent_station ./ speed, normal)
    penetration = max(-gap, zero(gap))
    effective_penetration = effective_contact_penetration(contact, gap)
    damping_multiplier = max(zero(gap),
        one(gap) - contact.damping_factor * gap_rate)
    (; curve_origin, roller, point_offset, contact_point, tangent, normal,
       tangent_error, tangent_error_rate, curvature, gap, gap_rate,
       penetration, effective_penetration, damping_multiplier)
end

function calculated_contact_force(contact, time, z, kinematics)
    contact.active[] || return zero(kinematics.gap)
    contact.expression && return contact.law(time, z)
    contact.stiffness * kinematics.effective_penetration *
        kinematics.damping_multiplier
end

function curve_contact_values(contact::PlanarCurveContactComponent, z,
        time = 0.0)
    kinematics = curve_contact_kinematics(contact, z)
    normal_force = calculated_contact_force(contact, time, z, kinematics)
    global_force = normal_force .* kinematics.normal
    (; kinematics..., normal_force, global_force)
end

curve_contact_gap(contact::PlanarCurveContactComponent, z) =
    z[contact.gap_variable]

function curve_contact_damping_surface(contact::PlanarCurveContactComponent, z)
    contact.expression && return one(eltype(z))
    gap = z[contact.gap_variable]
    gap < 0 ? one(gap) -
        contact.damping_factor * z[contact.gap_rate_variable] : one(gap)
end

function body_kinematic_dependencies(marker)
    body = marker.owner
    isnothing(body) && return Int[]
    [collect(body.position_variables); body.orientation_variable;
     collect(body.velocity_variables); body.angular_velocity_variable]
end

function body_configuration_dependencies(marker)
    body = marker.owner
    isnothing(body) && return Int[]
    [collect(body.position_variables); body.orientation_variable]
end

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

function contact_residual(contact, time, z)
    values = curve_contact_kinematics(contact, z)
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

function executable_blocks(contact::PlanarCurveContactComponent)
    dependencies = sort!(unique!([
        body_kinematic_dependencies(contact.curve_marker);
        body_kinematic_dependencies(contact.roller_marker);
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
        partials = local_state_jacobian(
            state -> contact_residual(contact, t, state), z, dependencies)
        jacobian[contact.contact_equations, dependencies] .+= partials
    end
    [ExecutableEquationBlock(contact.name, :contact,
        collect(contact.contact_equations), residual!, jacobian!)]
end

function add_curve_body_force!(equations, marker, values, force)
    body = marker.owner
    isnothing(body) && return nothing
    equations[body.balance_equations[1:2]] .-= force
    lever = values.contact_point - values.curve_origin.position +
        values.curve_origin.r
    equations[body.balance_equations[3]] -= dot(perpendicular(lever), force)
    nothing
end

function contact_body_contribution(contact, z)
    values = curve_contact_kinematics(contact, z)
    force = @view z[contact.global_force_variables]
    rows = Int[]
    contributions = eltype(z)[]
    if !isnothing(contact.roller_marker.owner)
        body = contact.roller_marker.owner
        local_equations = zeros(eltype(z), maximum(body.balance_equations))
        PlanarAppliedForces.add_body_point_force!(local_equations,
            contact.roller_marker.point, values.roller, force)
        append!(rows, body.balance_equations)
        append!(contributions, local_equations[body.balance_equations])
    end
    if !isnothing(contact.curve_marker.owner)
        body = contact.curve_marker.owner
        local_equations = zeros(eltype(z), maximum(body.balance_equations))
        add_curve_body_force!(local_equations, contact.curve_marker,
            values, -force)
        append!(rows, body.balance_equations)
        append!(contributions, local_equations[body.balance_equations])
    end
    rows, contributions
end

function equation_contributions(contact::PlanarCurveContactComponent)
    rows = Int[]
    !isnothing(contact.roller_marker.owner) &&
        append!(rows, contact.roller_marker.owner.balance_equations)
    !isnothing(contact.curve_marker.owner) &&
        append!(rows, contact.curve_marker.owner.balance_equations)
    dependencies = sort!(unique!([
        body_configuration_dependencies(contact.curve_marker);
        body_configuration_dependencies(contact.roller_marker);
        contact.station_variable;
        collect(contact.global_force_variables)]))
    function contribution_values(state)
        _, values = contact_body_contribution(contact, state)
        values
    end
    residual! = function (equations, t, z, zdot)
        equations[rows] .+= contribution_values(z)
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        partials = local_state_jacobian(contribution_values, z, dependencies)
        jacobian[rows, dependencies] .+= partials
    end
    [EquationContribution(contact.name, :load_to_bodies, rows,
        residual!, jacobian!)]
end

function closest_station(contact, state)
    roller = PlanarAppliedForces.point_marker_kinematics(
        contact.roller_marker.point, state).position
    origin = PlanarAppliedForces.point_marker_kinematics(
        contact.curve_marker.point, state).position
    angle = PlanarAppliedForces.marker_angle(contact.curve_marker.orientation,
        state)
    cosine, sine = cos(angle), sin(angle)
    rotation = [cosine -sine; sine cosine]
    samples = max(64, 8 * size(contact.curve.points, 2))
    stations = range(0.0, contact.curve.length; length = samples + 1)[1:end-1]
    distances = [sum(abs2, roller - (origin + rotation *
        curve_point(contact.curve, station).position)) for station in stations]
    station = stations[argmin(distances)]
    for _ in 1:12
        values = curve_point(contact.curve, station)
        point = origin + rotation * values.position
        first = rotation * values.derivative
        second = rotation * values.second_derivative
        separation = roller - point
        residual = dot(separation, first)
        derivative = -dot(first, first) + dot(separation, second)
        abs(derivative) > eps(Float64) || break
        increment = clamp(residual / derivative,
            -contact.curve.length / 8, contact.curve.length / 8)
        station -= increment
        abs(increment) <= 1.0e-12 * max(contact.curve.length, 1.0) && break
    end
    mod(station, contact.curve.length)
end

function initialize_planar_curve_contact!(state,
        contact::PlanarCurveContactComponent, time = 0.0;
        station = nothing)
    state[contact.station_variable] = isnothing(station) ?
        closest_station(contact, state) : Float64(station)
    # Tangency differentiated at fixed body motion is linear in station rate.
    state[contact.station_rate_variable] = 0.0
    base = curve_contact_kinematics(contact, state)
    trial = copy(state)
    trial[contact.station_rate_variable] = 1.0
    unit = curve_contact_kinematics(contact, trial)
    coefficient = unit.tangent_error_rate - base.tangent_error_rate
    abs(coefficient) > eps(Float64) &&
        (state[contact.station_rate_variable] =
            -base.tangent_error_rate / coefficient)
    values = curve_contact_kinematics(contact, state)
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
