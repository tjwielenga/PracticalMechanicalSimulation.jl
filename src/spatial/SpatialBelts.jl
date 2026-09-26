"""Elastic, no-slip belts tangent to spatial pulley pitch circles."""
module SpatialBelts

using ForwardDiff
using LinearAlgebra
using ..AutomaticAnalysis
using ..SpatialComponentAssembly
using ..SpatialModeling
using ..SpatialConstraints

import ..SpatialComponentAssembly: component_registration,
    executable_blocks, equation_contributions

export SpatialPulleyComponent, SpatialBeltComponent,
       SpatialBeltSpanComponent, spatial_belt_span_registration,
       spatial_belt_tangent_geometry, select_spatial_belt_tangent,
       spatial_belt_reference_path_length, allocated_spatial_belt_span,
       initialize_spatial_belt_span!, spatial_belt_span_values

"""A circular pitch surface carried by the first side of a revolute joint."""
struct SpatialPulleyComponent{B,J,M,T}
    name::Symbol
    body::B
    joint::J
    center_marker::M
    pitch_radius::T
end

"""One ordered, closed belt assembled from elastic tangent spans."""
struct SpatialBeltComponent{S,T}
    name::Symbol
    spans::Vector{S}
    initial_tension::T
    free_length::T
end

"""One elastic, no-slip straight span between two spatial pulleys."""
struct SpatialBeltSpanComponent{P1,P2,T}
    name::Symbol
    belt_name::Symbol
    pulley_1::P1
    pulley_2::P2
    parallel_axes::Bool
    branch_1::Int
    branch_2::Int
    reference_length::T
    reference_phase_1::T
    reference_phase_2::T
    beta_1::T
    beta_2::T
    initial_extension::T
    stiffness::T
    damping::T
    damping_time_scale::T
    feasibility_tolerance::T
    point_1_variables::UnitRange{Int}
    point_2_variables::UnitRange{Int}
    tangent_variables::UnitRange{Int}
    length_variable::Int
    extension_variable::Int
    extension_rate_variable::Int
    tension_variable::Int
    force_variables::UnitRange{Int}
    geometry_equations::UnitRange{Int}
    rate_equation::Int
    load_equations::UnitRange{Int}
end

function spatial_belt_span_registration(name::Symbol)
    variables = VariableDeclaration[
        [VariableDeclaration(Symbol(:P1_, axis), :applied_geometry, 0)
            for axis in (:x, :y, :z)];
        [VariableDeclaration(Symbol(:P2_, axis), :applied_geometry, 0)
            for axis in (:x, :y, :z)];
        [VariableDeclaration(Symbol(:t_, axis), :applied_geometry, 0)
            for axis in (:x, :y, :z)];
        VariableDeclaration(:length, :applied_geometry, 0);
        VariableDeclaration(:extension, :applied_geometry, 0);
        VariableDeclaration(:extension_rate, :applied_rate, 1);
        VariableDeclaration(:tension, :applied_load, 2);
        [VariableDeclaration(Symbol(:F_, axis), :applied_load, 2)
            for axis in (:x, :y, :z)]
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:geometry, EquationDeclaration[
            [EquationDeclaration(Symbol(:point_1_, axis),
                :applied_definition, 0, :belt_tangent)
                for axis in (:x, :y, :z)];
            [EquationDeclaration(Symbol(:point_2_, axis),
                :applied_definition, 0, :belt_tangent)
                for axis in (:x, :y, :z)];
            [EquationDeclaration(Symbol(:tangent_, axis),
                :applied_definition, 0, :belt_tangent)
                for axis in (:x, :y, :z)];
            EquationDeclaration(:length, :applied_definition, 0,
                :belt_length);
            EquationDeclaration(:extension, :applied_definition, 0,
                :belt_extension)
        ]),
        EquationBlockDeclaration(:rate, EquationDeclaration[
            EquationDeclaration(:extension_rate, :applied_definition, 1,
                :belt_extension_rate)
        ]),
        EquationBlockDeclaration(:load, EquationDeclaration[
            EquationDeclaration(:tension, :applied_definition, 2,
                :belt_tension);
            [EquationDeclaration(Symbol(:force_, axis),
                :applied_definition, 2, :global_force)
                for axis in (:x, :y, :z)]
        ])
    ]
    ComponentRegistration(name, variables, blocks)
end

component_registration(span::SpatialBeltSpanComponent) =
    spatial_belt_span_registration(span.name)
component_registration(::SpatialPulleyComponent) = nothing
component_registration(::SpatialBeltComponent) = nothing
executable_blocks(::SpatialPulleyComponent) = ExecutableEquationBlock[]
equation_contributions(::SpatialPulleyComponent) = EquationContribution[]
executable_blocks(::SpatialBeltComponent) = ExecutableEquationBlock[]
equation_contributions(::SpatialBeltComponent) = EquationContribution[]

unit_vector(vector, label) = begin
    magnitude = norm(vector)
    magnitude > zero(magnitude) || throw(DomainError(magnitude,
        "$label must have positive length"))
    vector ./ magnitude
end

function pulley_circle_kinematics(pulley::SpatialPulleyComponent, z)
    marker = pulley.center_marker
    body = pulley.body
    orientation = spatial_marker_orientation(marker, z)
    body_orientation = rotation_matrix(@view z[body.euler_parameter_variables])
    (; center = spatial_marker_position(marker, z),
       velocity = spatial_marker_velocity(marker, z),
       axis = orientation[:, 3],
       omega = body_orientation * z[body.angular_velocity_variables])
end

function parallel_tangent_geometry(center_1, axis_1, radius_1,
        center_2, axis_2, radius_2, tangent_kind, tangent_side, tolerance)
    alignment = dot(axis_1, axis_2)
    abs(alignment) >= 1 - tolerance || throw(DomainError(alignment,
        "parallel-axis belt tangent axes are no longer parallel"))
    separation = center_2 - center_1
    axial_offset = dot(separation, axis_1)
    scale = max(norm(separation), radius_1, radius_2, one(axial_offset))
    abs(axial_offset) <= tolerance * scale || throw(DomainError(axial_offset,
        "parallel pulley pitch circles are not coplanar"))
    planar_separation = separation - axial_offset .* axis_1
    center_distance = norm(planar_separation)
    center_distance > zero(center_distance) || throw(DomainError(
        center_distance, "belt pulley centers must be distinct"))
    center_unit = planar_separation ./ center_distance
    center_perpendicular = cross(axis_1, center_unit)
    ratio = (radius_1 - tangent_kind * radius_2) / center_distance
    abs(ratio) <= one(ratio) || throw(DomainError(ratio,
        "selected belt tangent is geometrically infeasible"))
    transverse = tangent_side * sqrt(max(zero(ratio), one(ratio) - ratio^2))
    normal_1 = ratio .* center_unit .+ transverse .* center_perpendicular
    normal_2 = tangent_kind .* normal_1
    point_1 = center_1 .+ radius_1 .* normal_1
    point_2 = center_2 .+ radius_2 .* normal_2
    delta = point_2 - point_1
    length = norm(delta)
    tangent = unit_vector(delta, "belt tangent span")
    tangency_error = max(abs(dot(tangent, normal_1)),
        abs(dot(tangent, normal_2)), abs(dot(tangent, axis_1)),
        abs(dot(tangent, axis_2)))
    tangency_error <= 10tolerance || throw(DomainError(tangency_error,
        "parallel pulley tangent failed its feasibility check"))
    beta_1 = radius_1 * dot(cross(axis_1, normal_1), tangent)
    beta_2 = radius_2 * dot(cross(axis_2, normal_2), tangent)
    (; point_1, point_2, normal_1, normal_2, tangent, length,
       axis_1, axis_2, beta_1, beta_2, feasibility_error = tangency_error)
end

function angled_tangent_geometry(center_1, axis_1, radius_1,
        center_2, axis_2, radius_2, radial_side_1, radial_side_2, tolerance)
    axis_cross = cross(axis_1, axis_2)
    axis_cross_norm = norm(axis_cross)
    axis_cross_norm > tolerance || throw(DomainError(axis_cross_norm,
        "angled-axis belt tangent axes have become parallel"))
    candidate_tangent = axis_cross ./ axis_cross_norm
    normal_1 = radial_side_1 .* unit_vector(cross(candidate_tangent, axis_1),
        "first pulley tangent radius")
    normal_2 = radial_side_2 .* unit_vector(cross(candidate_tangent, axis_2),
        "second pulley tangent radius")
    point_1 = center_1 .+ radius_1 .* normal_1
    point_2 = center_2 .+ radius_2 .* normal_2
    delta = point_2 - point_1
    length = norm(delta)
    tangent = unit_vector(delta, "belt tangent span")
    feasibility_error = max(abs(dot(tangent, normal_1)),
        abs(dot(tangent, normal_2)), abs(dot(tangent, axis_1)),
        abs(dot(tangent, axis_2)))
    feasibility_error <= tolerance || throw(DomainError(feasibility_error,
        "no straight common tangent exists for the current pulley geometry"))
    beta_1 = radius_1 * dot(cross(axis_1, normal_1), tangent)
    beta_2 = radius_2 * dot(cross(axis_2, normal_2), tangent)
    (; point_1, point_2, normal_1, normal_2, tangent, length,
       axis_1, axis_2, beta_1, beta_2, feasibility_error)
end

"""
Return one selected common tangent between two oriented pitch circles.

For parallel axes `branch_1` is the planar tangent kind and `branch_2` is its
side. For angled axes the two branches choose the radial direction on each
pitch circle. The latter geometry is feasible only when the line between the
candidate points is tangent to both circles.
"""
function spatial_belt_tangent_geometry(pulley_1, pulley_2, z,
        parallel_axes, branch_1, branch_2; tolerance = 1.0e-5)
    first = pulley_circle_kinematics(pulley_1, z)
    second = pulley_circle_kinematics(pulley_2, z)
    geometry = if parallel_axes
        parallel_tangent_geometry(first.center, first.axis,
            pulley_1.pitch_radius, second.center, second.axis,
            pulley_2.pitch_radius, branch_1, branch_2, tolerance)
    else
        angled_tangent_geometry(first.center, first.axis,
            pulley_1.pitch_radius, second.center, second.axis,
            pulley_2.pitch_radius, branch_1, branch_2, tolerance)
    end
    (; geometry..., pulley_1_kinematics = first,
       pulley_2_kinematics = second)
end

function select_spatial_belt_tangent(pulley_1, pulley_2, hint_1, hint_2,
        initial; tolerance = 1.0e-5)
    first = pulley_circle_kinematics(pulley_1, initial)
    second = pulley_circle_kinematics(pulley_2, initial)
    axis_cross = norm(cross(first.axis, second.axis))
    parallel_axes = axis_cross <= tolerance
    hint_position_1 = spatial_marker_position(hint_1, initial)
    hint_position_2 = spatial_marker_position(hint_2, initial)
    candidates = NamedTuple[]
    for branch_1 in (1, -1), branch_2 in (1, -1)
        geometry = try
            spatial_belt_tangent_geometry(pulley_1, pulley_2, initial,
                parallel_axes, branch_1, branch_2; tolerance)
        catch error
            error isa DomainError || rethrow()
            continue
        end
        score = sum(abs2, geometry.point_1 - hint_position_1) +
            sum(abs2, geometry.point_2 - hint_position_2)
        push!(candidates, (; parallel_axes, branch_1, branch_2,
            geometry, score))
    end
    isempty(candidates) && throw(ArgumentError(
        "no feasible straight common tangent exists between pulleys " *
        "'$(pulley_1.name)' and '$(pulley_2.name)'"))
    selected = candidates[argmin(getproperty.(candidates, :score))]
    if length(candidates) > 1
        scores = sort(getproperty.(candidates, :score))
        isapprox(scores[1], scores[2]; rtol = 1.0e-10,
            atol = 1.0e-12) && throw(ArgumentError(
            "belt tangent hints are ambiguous between pulleys " *
            "'$(pulley_1.name)' and '$(pulley_2.name)'"))
    end
    selected
end

function pulley_contact_phase(pulley, normal, z)
    hinge = pulley.joint.hinge
    isempty(hinge.rotation_variables) && error(
        "pulley revolute joint must expose rotation coordinates")
    theta = z[hinge.rotation_variables[3]]
    base_orientation = spatial_marker_orientation(pulley.joint.marker_b, z)
    axis = base_orientation[:, 3]
    reference = base_orientation[:, 1]
    projected = normal - dot(normal, axis) .* axis
    projected = projected ./ norm(projected)
    contact_angle = atan(dot(axis, cross(reference, projected)),
        dot(reference, projected))
    theta - contact_angle
end

function spatial_belt_span_values(span::SpatialBeltSpanComponent, z)
    geometry = spatial_belt_tangent_geometry(span.pulley_1, span.pulley_2,
        z, span.parallel_axes, span.branch_1, span.branch_2;
        tolerance = span.feasibility_tolerance)
    phase_1 = pulley_contact_phase(span.pulley_1, geometry.normal_1, z)
    phase_2 = pulley_contact_phase(span.pulley_2, geometry.normal_2, z)
    expected_extension = geometry.length - span.reference_length +
        span.beta_2 * (phase_2 - span.reference_phase_2) -
        span.beta_1 * (phase_1 - span.reference_phase_1) +
        span.initial_extension
    radial_1 = geometry.point_1 - geometry.pulley_1_kinematics.center
    radial_2 = geometry.point_2 - geometry.pulley_2_kinematics.center
    surface_velocity_1 = geometry.pulley_1_kinematics.velocity +
        cross(geometry.pulley_1_kinematics.omega, radial_1)
    surface_velocity_2 = geometry.pulley_2_kinematics.velocity +
        cross(geometry.pulley_2_kinematics.omega, radial_2)
    expected_extension_rate = dot(surface_velocity_2, geometry.tangent) -
        dot(surface_velocity_1, geometry.tangent)
    (; geometry, actual_point_1 = geometry.point_1,
       actual_point_2 = geometry.point_2,
       actual_normal_1 = geometry.normal_1,
       actual_normal_2 = geometry.normal_2,
       actual_tangent = geometry.tangent,
       actual_length = geometry.length,
       phase_1, phase_2, expected_extension,
       expected_extension_rate, surface_velocity_1, surface_velocity_2,
       point_1 = @view(z[span.point_1_variables]),
       point_2 = @view(z[span.point_2_variables]),
       tangent_variable = @view(z[span.tangent_variables]),
       length_variable = z[span.length_variable],
       extension = z[span.extension_variable],
       extension_rate = z[span.extension_rate_variable],
       tension = z[span.tension_variable],
       force = @view(z[span.force_variables]))
end

function allocated_spatial_belt_span(layout, name, belt_name, pulley_1,
        pulley_2, selected, initial_state, initial_extension, stiffness, damping,
        damping_time_scale; feasibility_tolerance = 1.0e-5)
    variables = component_variable_indices(layout, name)
    geometry_rows = component_equation_indices(layout, name, :geometry)
    selected_geometry = selected.geometry
    phase_1 = pulley_contact_phase(pulley_1,
        selected_geometry.normal_1, initial_state)
    phase_2 = pulley_contact_phase(pulley_2,
        selected_geometry.normal_2, initial_state)
    SpatialBeltSpanComponent(name, belt_name, pulley_1, pulley_2,
        selected.parallel_axes, selected.branch_1, selected.branch_2,
        selected_geometry.length, phase_1, phase_2,
        selected_geometry.beta_1, selected_geometry.beta_2,
        initial_extension, stiffness, damping, damping_time_scale,
        feasibility_tolerance, variables[1:3], variables[4:6],
        variables[7:9], variables[10], variables[11], variables[12],
        variables[13], variables[14:16], geometry_rows,
        only(component_equation_indices(layout, name, :rate)),
        component_equation_indices(layout, name, :load))
end

function initialize_spatial_belt_span!(initial, span)
    values = spatial_belt_span_values(span, initial)
    initial[span.point_1_variables] .= values.actual_point_1
    initial[span.point_2_variables] .= values.actual_point_2
    initial[span.tangent_variables] .= values.actual_tangent
    initial[span.length_variable] = values.actual_length
    initial[span.extension_variable] = values.expected_extension
    initial[span.extension_rate_variable] = values.expected_extension_rate
    initial[span.tension_variable] = span.stiffness * values.expected_extension +
        span.damping * values.expected_extension_rate
    initial[span.force_variables] .=
        initial[span.tension_variable] .* values.actual_tangent
    initial
end

function spatial_belt_reference_path_length(specifications)
    straight_length = sum(spec.geometry.length for spec in specifications)
    wrap_length = 0.0
    for index in eachindex(specifications)
        incoming = specifications[index]
        outgoing = specifications[mod1(index + 1, length(specifications))]
        incoming.pulley_2 === outgoing.pulley_1 || throw(ArgumentError(
            "ordered spatial belt spans must form a closed pulley loop"))
        pulley = incoming.pulley_2
        beta_in = incoming.geometry.beta_2
        beta_out = outgoing.geometry.beta_1
        sign(beta_in) == sign(beta_out) || throw(ArgumentError(
            "belt tangent hints produce incompatible wrap directions at " *
            "pulley '$(pulley.name)'"))
        axis = incoming.geometry.axis_2
        normal_in = incoming.geometry.normal_2
        normal_out = outgoing.geometry.normal_1
        signed_angle = atan(dot(axis, cross(normal_in, normal_out)),
            dot(normal_in, normal_out))
        wrap_angle = mod(sign(beta_in) * signed_angle, 2pi)
        wrap_length += pulley.pitch_radius * wrap_angle
    end
    straight_length + wrap_length
end

function belt_body_dependencies(pulley; velocity = false)
    columns = Int[]
    bodies = Any[pulley.body]
    base = pulley.joint.marker_b
    base isa SpatialBodyMarker && push!(bodies, base.body)
    for body in unique(bodies)
        append!(columns, body.position_variables)
        append!(columns, body.euler_parameter_variables)
        velocity && append!(columns, body.velocity_variables)
        velocity && append!(columns, body.angular_velocity_variables)
    end
    append!(columns, pulley.joint.hinge.rotation_variables)
    sort!(unique!(columns))
end

function local_jacobian(function_value, z, columns)
    inputs = collect(z[columns])
    ForwardDiff.jacobian(inputs) do local_values
        state = Vector{eltype(local_values)}(undef, length(z))
        state .= z
        state[columns] .= local_values
        function_value(state)
    end
end

function geometry_residual(span, z)
    values = spatial_belt_span_values(span, z)
    [values.point_1 - values.actual_point_1;
     values.point_2 - values.actual_point_2;
     values.tangent_variable - values.actual_tangent;
     values.length_variable - values.actual_length;
     values.extension - values.expected_extension]
end

function rate_residual(span, z)
    values = spatial_belt_span_values(span, z)
    [values.extension_rate - values.expected_extension_rate]
end

function load_residual(span, z)
    values = spatial_belt_span_values(span, z)
    [values.tension - span.stiffness * values.extension -
         span.damping * values.extension_rate;
     values.force - values.tension .* values.tangent_variable]
end

function executable_blocks(span::SpatialBeltSpanComponent)
    geometry_rows = collect(span.geometry_equations)
    rate_rows = [span.rate_equation]
    load_rows = collect(span.load_equations)
    geometry_dependencies = sort!(unique!([
        belt_body_dependencies(span.pulley_1);
        belt_body_dependencies(span.pulley_2);
        collect(span.point_1_variables); collect(span.point_2_variables);
        collect(span.tangent_variables); span.length_variable;
        span.extension_variable]))
    rate_dependencies = sort!(unique!([
        belt_body_dependencies(span.pulley_1; velocity = true);
        belt_body_dependencies(span.pulley_2; velocity = true);
        span.extension_rate_variable]))
    load_dependencies = [collect(span.tangent_variables);
        span.extension_variable; span.extension_rate_variable;
        span.tension_variable; collect(span.force_variables)]
    geometry! = (equations, t, z, zdot) ->
        (equations[geometry_rows] .= geometry_residual(span, z))
    rate! = (equations, t, z, zdot) ->
        (equations[rate_rows] .= rate_residual(span, z))
    load! = (equations, t, z, zdot) ->
        (equations[load_rows] .= load_residual(span, z))
    geometry_jacobian! = (jacobian, t, z, zdot, coefficient) ->
        (jacobian[geometry_rows, geometry_dependencies] .+=
            local_jacobian(state -> geometry_residual(span, state), z,
                geometry_dependencies))
    rate_jacobian! = (jacobian, t, z, zdot, coefficient) ->
        (jacobian[rate_rows, rate_dependencies] .+=
            local_jacobian(state -> rate_residual(span, state), z,
                rate_dependencies))
    load_jacobian! = (jacobian, t, z, zdot, coefficient) ->
        (jacobian[load_rows, load_dependencies] .+=
            local_jacobian(state -> load_residual(span, state), z,
                load_dependencies))
    [ExecutableEquationBlock(span.name, :geometry, geometry_rows,
         geometry!, geometry_jacobian!),
     ExecutableEquationBlock(span.name, :rate, rate_rows,
         rate!, rate_jacobian!),
     ExecutableEquationBlock(span.name, :load, load_rows,
         load!, load_jacobian!)]
end

function add_belt_point_force!(equations, body, z, point, force)
    orientation = rotation_matrix(@view z[body.euler_parameter_variables])
    lever_body = transpose(orientation) *
        (point - z[body.position_variables])
    force_body = transpose(orientation) * force
    equations[body.balance_equations[1:3]] .-= force
    equations[body.balance_equations[4:6]] .-=
        cross(lever_body, force_body)
    nothing
end

function equation_contributions(span::SpatialBeltSpanComponent)
    body_1 = span.pulley_1.body
    body_2 = span.pulley_2.body
    rows = [collect(body_1.balance_equations);
            collect(body_2.balance_equations)]
    columns = sort!(unique!([
        collect(body_1.position_variables);
        collect(body_1.euler_parameter_variables);
        collect(body_2.position_variables);
        collect(body_2.euler_parameter_variables);
        collect(span.point_1_variables); collect(span.point_2_variables);
        collect(span.force_variables)]))
    residual! = function (equations, t, z, zdot)
        add_belt_point_force!(equations, body_1, z,
            @view(z[span.point_1_variables]), @view(z[span.force_variables]))
        add_belt_point_force!(equations, body_2, z,
            @view(z[span.point_2_variables]), -@view(z[span.force_variables]))
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        # Build the two six-row body contributions without relying on their
        # global equation numbers inside the differentiated local function.
        local_value = function (state)
            output = zeros(eltype(state), 12)
            for (offset, body, point, force) in
                    ((0, body_1, @view(state[span.point_1_variables]),
                        @view(state[span.force_variables])),
                     (6, body_2, @view(state[span.point_2_variables]),
                        -@view(state[span.force_variables])))
                orientation = rotation_matrix(
                    @view state[body.euler_parameter_variables])
                lever_body = transpose(orientation) *
                    (point - state[body.position_variables])
                force_body = transpose(orientation) * force
                output[offset .+ (1:3)] .-= force
                output[offset .+ (4:6)] .-=
                    cross(lever_body, force_body)
            end
            output
        end
        jacobian[rows, columns] .+= local_jacobian(local_value, z, columns)
    end
    [EquationContribution(span.name, :load_to_bodies, rows,
        residual!, jacobian!)]
end

end
