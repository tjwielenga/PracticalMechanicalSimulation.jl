"""Ideal spatial spur rack-and-pinion constraints and contact reactions."""
module SpatialRackAndPinions

using LinearAlgebra
using ..AutomaticAnalysis
using ..SpatialComponentAssembly
using ..SpatialModeling
using ..SpatialDirectedDistances: marker_axis_kinematics
using ..SpatialConstraints

import ..SpatialComponentAssembly: component_registration,
    executable_blocks, equation_contributions

export SpatialRackAndPinion, spatial_rack_and_pinion_registration,
       spatial_rack_and_pinion_geometry,
       allocated_spatial_rack_and_pinion,
       spatial_rack_and_pinion_position,
       spatial_rack_and_pinion_velocity,
       spatial_rack_and_pinion_acceleration,
       spatial_rack_contact_direction

"""
An ideal spur rack and pinion coupling one inline distance to one revolute
angle. Generated floating markers apply the equal-and-opposite tooth forces
at the carrier-fixed pitch contact point.
"""
struct SpatialRackAndPinion{I,R,B1,B2,C,F1,F2,T}
    name::Symbol
    inline::I
    revolute::R
    rack_body::B1
    pinion_body::B2
    contact_marker::C
    rack_contact::F1
    pinion_contact::F2
    pitch_radius::T
    rolling_radius::T
    phase::T
    reaction_variable::Int
    acceleration_equation::Int
    velocity_equation::Int
    position_equation::Int
end

function spatial_rack_and_pinion_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:lambda, :reaction, 2),
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:acceleration, [EquationDeclaration(
            :Phi_ddot, :constraint, 2, :rack_and_pinion)]),
        EquationBlockDeclaration(:velocity, [EquationDeclaration(
            :Phi_dot, :constraint, 1, :rack_and_pinion)]),
        EquationBlockDeclaration(:position, [EquationDeclaration(
            :Phi, :constraint, 0, :rack_and_pinion)]),
    ]
    families = ConstraintFamilyDeclaration[
        ConstraintFamilyDeclaration(:rack_and_pinion, :lambda, :Phi,
            :Phi_dot, :Phi_ddot),
    ]
    ComponentRegistration(name, variables, blocks, families)
end

component_registration(component::SpatialRackAndPinion) =
    spatial_rack_and_pinion_registration(component.name)

fixed_marker_owner(::SpatialGroundMarker) = nothing
fixed_marker_owner(marker::SpatialBodyMarker) = marker.body

"""
Validate the ordinary spur arrangement and return its carrier-fixed contact
geometry. The inline base marker uses `x` for the pinion axis, `y` for the
contact radius, and `z` for rack travel.
"""
function spatial_rack_and_pinion_geometry(inline, revolute, pitch_radius,
        initial; alignment_tolerance = 1.0e-5)
    inline isa SpatialInlineConstraint || throw(ArgumentError(
        "spatial rack and pinion first joint must be an inline constraint"))
    revolute isa SpatialRevoluteJoint || throw(ArgumentError(
        "spatial rack and pinion second joint must be a revolute joint"))
    inline.marker_i isa SpatialBodyMarker || throw(ArgumentError(
        "spatial rack and pinion inline constraint must list the rack marker first"))
    revolute.marker_a isa SpatialBodyMarker || throw(ArgumentError(
        "spatial rack and pinion revolute joint must list the pinion marker first"))
    pitch_radius > 0 || throw(ArgumentError(
        "spatial rack and pinion pitch_radius must be positive"))
    alignment_tolerance > 0 || throw(ArgumentError(
        "spatial rack and pinion alignment_tolerance must be positive"))

    carrier = fixed_marker_owner(inline.marker_j)
    fixed_marker_owner(revolute.marker_b) === carrier || throw(ArgumentError(
        "spatial rack and pinion joints must share a carrier"))
    rack_body = inline.marker_i.body
    pinion_body = revolute.marker_a.body
    rack_body !== pinion_body || throw(ArgumentError(
        "spatial rack and pinion requires different rack and pinion bodies"))

    guide_orientation = spatial_marker_orientation(inline.marker_j, initial)
    guide_x = collect(@view guide_orientation[:, 1])
    guide_y = collect(@view guide_orientation[:, 2])
    guide_z = collect(@view guide_orientation[:, 3])
    pinion_axis = marker_axis_kinematics(
        revolute.marker_b, initial, 3).direction
    alignment = clamp(dot(pinion_axis, guide_x), -1.0, 1.0)
    mismatch = acos(abs(alignment))
    mismatch <= alignment_tolerance || throw(ArgumentError(
        "spatial rack and pinion pinion axis must be parallel or " *
        "antiparallel to the inline base marker's x-axis " *
        "(angular mismatch $(round(rad2deg(mismatch); sigdigits = 6)) degrees)"))

    center = spatial_marker_position(revolute.marker_b, initial)
    contact_position = center + pitch_radius .* guide_y
    rolling_radius = pitch_radius * alignment
    (; carrier, rack_body, pinion_body, guide_orientation,
       guide_x, guide_y, guide_z, pinion_axis, contact_position,
       rolling_radius)
end

function allocated_spatial_rack_and_pinion(layout, name, inline, revolute,
        contact_marker, rack_contact, pinion_contact, pitch_radius, initial;
        alignment_tolerance = 1.0e-5, phase = "initial")
    geometry = spatial_rack_and_pinion_geometry(
        inline, revolute, pitch_radius, initial; alignment_tolerance)
    variables = component_variable_indices(layout, name)
    component = SpatialRackAndPinion(name, inline, revolute,
        geometry.rack_body, geometry.pinion_body, contact_marker,
        rack_contact, pinion_contact, pitch_radius, geometry.rolling_radius,
        0.0, only(variables),
        only(component_equation_indices(layout, name, :acceleration)),
        only(component_equation_indices(layout, name, :velocity)),
        only(component_equation_indices(layout, name, :position)))
    initial_phase = spatial_rack_and_pinion_position(component, initial)
    phase_value = if phase isa Number
        Float64(phase)
    elseif phase == "initial"
        initial_phase
    else
        throw(ArgumentError(
            "rack and pinion '$name' phase must be numeric or 'initial'"))
    end
    isfinite(phase_value) || throw(ArgumentError(
        "rack and pinion '$name' phase must be finite"))
    SpatialRackAndPinion(component.name, component.inline,
        component.revolute, component.rack_body, component.pinion_body,
        component.contact_marker, component.rack_contact,
        component.pinion_contact, component.pitch_radius,
        component.rolling_radius, phase_value, component.reaction_variable,
        component.acceleration_equation, component.velocity_equation,
        component.position_equation)
end

function rack_and_pinion_variables(component)
    inline = component.inline.translation_variables
    rotation = component.revolute.hinge.rotation_variables
    isempty(inline) && error(
        "rack and pinion inline coordinate was not allocated")
    isempty(rotation) && error(
        "rack and pinion revolute coordinate was not allocated")
    (; qddot = inline[1], qdot = inline[2], q = inline[3],
       alpha = rotation[1], omega = rotation[2], theta = rotation[3])
end

function spatial_rack_and_pinion_position(component::SpatialRackAndPinion, z)
    variables = rack_and_pinion_variables(component)
    z[variables.q] - component.rolling_radius * z[variables.theta] -
        component.phase
end

function spatial_rack_and_pinion_velocity(component::SpatialRackAndPinion, z)
    variables = rack_and_pinion_variables(component)
    z[variables.qdot] - component.rolling_radius * z[variables.omega]
end

function spatial_rack_and_pinion_acceleration(
        component::SpatialRackAndPinion, z)
    variables = rack_and_pinion_variables(component)
    z[variables.qddot] - component.rolling_radius * z[variables.alpha]
end

function executable_blocks(component::SpatialRackAndPinion)
    variables = rack_and_pinion_variables(component)
    rows = (component.position_equation, component.velocity_equation,
        component.acceleration_equation)
    columns = ((variables.q, variables.theta),
        (variables.qdot, variables.omega),
        (variables.qddot, variables.alpha))
    names = (:position, :velocity, :acceleration)
    offsets = (component.phase, 0.0, 0.0)
    blocks = ExecutableEquationBlock[]
    for (name, row, pair, offset) in zip(names, rows, columns, offsets)
        residual! = function (equations, t, z, zdot)
            equations[row] = z[pair[1]] -
                component.rolling_radius * z[pair[2]] - offset
        end
        jacobian! = function (jacobian, t, z, zdot, coefficient)
            jacobian[row, pair[1]] += 1
            jacobian[row, pair[2]] -= component.rolling_radius
        end
        push!(blocks, ExecutableEquationBlock(component.name, name, [row],
            residual!, jacobian!))
    end
    blocks
end

spatial_rack_contact_direction(component::SpatialRackAndPinion, z) =
    inline_axis(component.inline, z)

function body_force_contribution(body, z, point, force)
    parameters = @view z[body.euler_parameter_variables]
    orientation = rotation_matrix(parameters)
    lever_body = transpose(orientation) *
        (point - z[body.position_variables])
    force_body = transpose(orientation) * force
    [-force; -cross(lever_body, force_body)]
end

function rack_and_pinion_reaction_vector(component, z)
    point = spatial_marker_position(component.contact_marker, z)
    force = z[component.reaction_variable] .*
        spatial_rack_contact_direction(component, z)
    [body_force_contribution(component.rack_body, z, point, force);
     body_force_contribution(component.pinion_body, z, point, -force)]
end

function reaction_dependencies(component)
    columns = Int[component.reaction_variable]
    for body in (component.rack_body, component.pinion_body)
        append!(columns, body.position_variables)
        append!(columns, body.euler_parameter_variables)
    end
    carrier = fixed_marker_owner(component.contact_marker)
    if !isnothing(carrier)
        append!(columns, carrier.position_variables)
        append!(columns, carrier.euler_parameter_variables)
    end
    sort!(unique!(columns))
end

function equation_contributions(component::SpatialRackAndPinion)
    rows = [collect(component.rack_body.balance_equations);
            collect(component.pinion_body.balance_equations)]
    columns = reaction_dependencies(component)
    residual! = function (equations, t, z, zdot)
        equations[rows] .+= rack_and_pinion_reaction_vector(component, z)
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        jacobian[rows, columns] .+= spatial_local_state_jacobian(
            state -> rack_and_pinion_reaction_vector(component, state),
            z, columns)
    end
    [EquationContribution(component.name, :reaction_to_bodies, rows,
        residual!, jacobian!)]
end

end
