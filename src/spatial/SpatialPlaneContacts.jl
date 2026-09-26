"""One-sided compliant sphere-plane contact for spatial models."""
module SpatialPlaneContacts

using LinearAlgebra
using ..AutomaticAnalysis
using ..SpatialComponentAssembly
using ..SpatialModeling
using ..SpatialDirectedDistances
using ..ScalarExpressions: ScalarLaw

import ..SpatialComponentAssembly: component_registration,
    executable_blocks, equation_contributions

export SpatialPlaneContactComponent, spatial_plane_contact_registration,
       allocated_spatial_plane_contact, initialize_spatial_plane_contact!,
       set_spatial_plane_contact_stage!, spatial_plane_contact_values,
       spatial_plane_contact_gap, spatial_plane_contact_damping_surface

"""
One-sided compliant contact between a marker-centered sphere and an oriented
plane. The first marker locates the sphere center. The second marker's local
z-axis is the outward plane normal.
"""
struct SpatialPlaneContactComponent{M1,M2,G,L,T,S}
    name::Symbol
    sphere_marker::M1
    plane_marker::M2
    geometry::G
    law::L
    expression::Bool
    radius::T
    stiffness::T
    damping_factor::T
    transition_depth::T
    active_during::S
    active::Base.RefValue{Bool}
    gap_variable::Int
    gap_rate_variable::Int
    normal_force_variable::Int
    global_force_variables::UnitRange{Int}
    contact_equations::UnitRange{Int}
end

function spatial_plane_contact_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:gap, :applied_geometry, 0),
        VariableDeclaration(:gap_rate, :applied_rate, 1),
        VariableDeclaration(:normal_force, :applied_load, 2),
        [VariableDeclaration(Symbol(:F_, axis), :applied_load, 2)
            for axis in (:x, :y, :z)]...
    ]
    equations = EquationDeclaration[
        EquationDeclaration(:gap, :applied_definition, 0, :contact_gap),
        EquationDeclaration(:gap_rate, :applied_definition, 1,
            :contact_gap_rate),
        EquationDeclaration(:normal_force, :applied_definition, 2,
            :contact_force),
        [EquationDeclaration(Symbol(:global_force_, axis),
            :applied_definition, 2, :global_force)
            for axis in (:x, :y, :z)]...
    ]
    ComponentRegistration(name, variables,
        [EquationBlockDeclaration(:contact, equations)])
end

component_registration(contact::SpatialPlaneContactComponent) =
    spatial_plane_contact_registration(contact.name)

function allocated_spatial_plane_contact(layout, name, sphere_marker,
        plane_marker, radius, stiffness, damping_factor;
        law = nothing, expression = false,
        transition_depth = 0.0,
        active_during = (:static, :dynamic, :modal))
    variables = component_variable_indices(layout, name)
    equations = component_equation_indices(layout, name, :contact)
    geometry = SpatialDirectedDistance(sphere_marker, plane_marker,
        SpatialDirectedAxis(plane_marker, 3))
    SpatialPlaneContactComponent(name, sphere_marker, plane_marker, geometry,
        law, Bool(expression), Float64(radius), Float64(stiffness),
        Float64(damping_factor),
        Float64(transition_depth), active_during,
        Ref(:dynamic in active_during), variables[1],
        variables[2], variables[3], variables[4:6], equations)
end

"""Select whether the contact force is active in the current analysis stage."""
function set_spatial_plane_contact_stage!(contact::SpatialPlaneContactComponent,
        stage)
    contact.active[] = stage in contact.active_during
    contact
end

"""Effective penetration with a fifth-derivative-smooth contact entry."""
function effective_contact_penetration(
        contact::SpatialPlaneContactComponent, gap)
    penetration = max(-gap, zero(gap))
    depth = contact.transition_depth
    if !(depth > 0) || penetration >= depth
        return penetration - (depth > 0 ? depth / 2 : zero(depth))
    end
    x = penetration / depth
    # Integral of the ninth-degree smootherstep stiffness
    # 126x^5 - 420x^6 + 540x^7 - 315x^8 + 70x^9.
    depth * x^6 * (21 + x * (-60 + x * (67.5 + x * (-35 + 7x))))
end

function spatial_plane_contact_kinematics(
        contact::SpatialPlaneContactComponent, z)
    geometry = directed_distance_values(contact.geometry, z)
    normal = geometry.direction.direction
    gap = geometry.position - contact.radius
    gap_rate = geometry.velocity
    penetration = max(-gap, zero(gap))
    effective_penetration = effective_contact_penetration(contact, gap)
    damping_multiplier = max(zero(gap),
        one(gap) - contact.damping_factor * gap_rate)
    contact_point = geometry.first.position - geometry.position .* normal
    (; geometry, normal, gap, gap_rate, penetration, effective_penetration,
       damping_multiplier, contact_point)
end

function calculated_plane_contact_force(contact, time, z, kinematics)
    contact.active[] || return zero(kinematics.gap)
    contact.expression && return contact.law(time, z)
    contact.stiffness * kinematics.effective_penetration *
        kinematics.damping_multiplier
end

function spatial_plane_contact_values(contact::SpatialPlaneContactComponent,
        z, time = 0.0)
    kinematics = spatial_plane_contact_kinematics(contact, z)
    normal_force = calculated_plane_contact_force(
        contact, time, z, kinematics)
    normal = kinematics.normal
    global_force = normal_force .* normal
    (; kinematics..., normal_force, global_force)
end

spatial_plane_contact_gap(contact::SpatialPlaneContactComponent, z) =
    z[contact.gap_variable]

function spatial_plane_contact_damping_surface(
        contact::SpatialPlaneContactComponent, z)
    contact.expression && return one(eltype(z))
    gap = z[contact.gap_variable]
    gap < 0 ? one(gap) -
        contact.damping_factor * z[contact.gap_rate_variable] : one(gap)
end

function initialize_spatial_plane_contact!(initial,
        contact::SpatialPlaneContactComponent, time = 0.0)
    kinematics = spatial_plane_contact_kinematics(contact, initial)
    initial[contact.gap_variable] = kinematics.gap
    initial[contact.gap_rate_variable] = kinematics.gap_rate
    normal_force = calculated_plane_contact_force(
        contact, time, initial, kinematics)
    isfinite(normal_force) || throw(ArgumentError(
        "plane contact '$(contact.name)' force is not finite initially"))
    initial[contact.normal_force_variable] = normal_force
    initial[contact.global_force_variables] .=
        normal_force .* kinematics.normal
    initial
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
    [collect(body.position_variables);
     collect(body.euler_parameter_variables)]
end

function contact_residual(contact, time, z)
    values = spatial_plane_contact_kinematics(contact, z)
    gap = z[contact.gap_variable]
    gap_rate = z[contact.gap_rate_variable]
    calculated_force = calculated_plane_contact_force(
        contact, time, z, values)
    [gap - values.gap;
     gap_rate - values.gap_rate;
     z[contact.normal_force_variable] - calculated_force;
     z[contact.global_force_variables] .-
        z[contact.normal_force_variable] .* values.normal]
end

function executable_blocks(contact::SpatialPlaneContactComponent)
    dependencies = sort!(unique!([
        contact_kinematic_dependencies(contact.sphere_marker);
        contact_kinematic_dependencies(contact.plane_marker);
        contact.gap_variable;
        contact.gap_rate_variable;
        contact.normal_force_variable;
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
    parameters = @view z[body.euler_parameter_variables]
    orientation = rotation_matrix(parameters)
    lever_body = transpose(orientation) *
        (point - z[body.position_variables])
    force_body = transpose(orientation) * global_force
    [-global_force; -cross(lever_body, force_body)]
end

function contact_body_contribution(contact, z)
    values = spatial_plane_contact_kinematics(contact, z)
    global_force = @view z[contact.global_force_variables]
    contributions = eltype(z)[]
    if contact.sphere_marker isa SpatialBodyMarker
        append!(contributions, body_force_contribution(
            contact.sphere_marker.body, z, values.contact_point, global_force))
    end
    if contact.plane_marker isa SpatialBodyMarker
        append!(contributions, body_force_contribution(
            contact.plane_marker.body, z, values.contact_point, -global_force))
    end
    contributions
end

function equation_contributions(contact::SpatialPlaneContactComponent)
    rows = Int[]
    contact.sphere_marker isa SpatialBodyMarker &&
        append!(rows, contact.sphere_marker.body.balance_equations)
    contact.plane_marker isa SpatialBodyMarker &&
        append!(rows, contact.plane_marker.body.balance_equations)
    dependencies = sort!(unique!([
        contact_configuration_dependencies(contact.sphere_marker);
        contact_configuration_dependencies(contact.plane_marker);
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

end
