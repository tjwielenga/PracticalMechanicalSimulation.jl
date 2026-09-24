"""Ideal constraint components for the spatial modeler."""
module SpatialConstraints

using ForwardDiff
using LinearAlgebra
using ..AutomaticAnalysis
using ..SpatialComponentAssembly
using ..SpatialModeling
using ..SpatialDirectedDistances

import ..SpatialComponentAssembly: component_registration,
    executable_blocks, equation_contributions

export SpatialSphericalJoint, SpatialPerpConstraint, SpatialCVPhaseConstraint,
       SpatialHingeConstraint,
       SpatialOrientConstraint, SpatialRevoluteJoint, SpatialInplaneConstraint,
       SpatialFixedJoint, SpatialConstantVelocityJoint,
       SpatialInlineConstraint, SpatialCylindricalJoint,
       SpatialTranslationalJoint,
       spherical_joint_registration, perp_constraint_registration,
       cv_phase_constraint_registration,
       hinge_constraint_registration, orient_constraint_registration,
       revolute_joint_registration, fixed_joint_registration,
       constant_velocity_joint_registration,
       inplane_constraint_registration, inline_constraint_registration,
       cylindrical_joint_registration, translational_joint_registration,
       allocated_spherical_joint, allocated_perp_constraint,
       allocated_cv_phase_constraint,
       allocated_hinge_constraint, allocated_orient_constraint,
       allocated_revolute_joint, allocated_fixed_joint,
       allocated_constant_velocity_joint,
       allocated_inplane_constraint, allocated_inline_constraint,
       allocated_cylindrical_joint, allocated_translational_joint,
       perp_normal, perp_position, perp_velocity, perp_acceleration,
       cv_phase_position, cv_phase_velocity, cv_phase_acceleration,
       inplane_normal, inplane_position, inplane_velocity,
       inplane_acceleration,
       directed_distance_reaction_rows, add_directed_distance_reaction!,
       add_directed_distance_reaction_jacobian!,
       inline_axis, inline_distance, inline_velocity, inline_acceleration,
       hinge_angle, hinge_angular_velocity, hinge_angular_acceleration,
       initialize_hinge_coordinates!, initialize_inline_coordinates!,
       connection_hinge, connection_inline

"""Three coincident-point constraints and their global reaction vector."""
struct SpatialSphericalJoint{A,B}
    name::Symbol
    marker_a::A
    marker_b::B
    reaction_variables::UnitRange{Int}
    acceleration_equations::UnitRange{Int}
    velocity_equations::UnitRange{Int}
    position_equations::UnitRange{Int}
end

"""One perpendicular-axis constraint and its scalar torque reaction."""
struct SpatialPerpConstraint{A,B}
    name::Symbol
    marker_i::A
    marker_j::B
    axis_i::Int
    axis_j::Int
    reaction_variable::Int
    acceleration_equation::Int
    velocity_equation::Int
    position_equation::Int
end

"""One ideal constant-velocity phase constraint between two shaft frames."""
struct SpatialCVPhaseConstraint{A,B}
    name::Symbol
    marker_i::A
    marker_j::B
    reaction_variable::Int
    acceleration_equation::Int
    velocity_equation::Int
    position_equation::Int
end

"""Two perpendicular-axis primitives leaving relative rotation about `z_j`."""
struct SpatialHingeConstraint{A,B,P,Q}
    name::Symbol
    marker_i::A
    marker_j::B
    perp_xz::P
    perp_yz::Q
    rotation_variables::Vector{Int}
    rotation_equations::Vector{Int}
    state_equations::Vector{Int}
end

"""Three perpendicular-axis primitives fixing relative marker orientation."""
struct SpatialOrientConstraint{A,B,H,P}
    name::Symbol
    marker_i::A
    marker_j::B
    hinge::H
    perp_xy::P
end

"""A spherical joint and hinge constraint sharing the same two markers."""
struct SpatialRevoluteJoint{A,B,S,H}
    name::Symbol
    marker_a::A
    marker_b::B
    spherical::S
    hinge::H
end

"""A spherical joint and orient constraint sharing the same two markers."""
struct SpatialFixedJoint{A,B,S,O}
    name::Symbol
    marker_a::A
    marker_b::B
    spherical::S
    orient::O
end

"""A spherical joint and constant-velocity phase constraint at one point."""
struct SpatialConstantVelocityJoint{A,B,S,C}
    name::Symbol
    marker_i::A
    marker_j::B
    spherical::S
    phase::C
end

"""One point constrained to a plane normal to a selected axis of marker `j`."""
struct SpatialInplaneConstraint{G}
    name::Symbol
    geometry::G
    reaction_variable::Int
    acceleration_equation::Int
    velocity_equation::Int
    position_equation::Int
end

"""Two inplane primitives leaving translation along `z_j`."""
struct SpatialInlineConstraint{A,B,X,Y,Z}
    name::Symbol
    marker_i::A
    marker_j::B
    inplane_x::X
    inplane_y::Y
    axial_geometry::Z
    translation_variables::Vector{Int}
    translation_equations::Vector{Int}
    state_equations::Vector{Int}
end

"""An inline and hinge constraint leaving translation and rotation about `z_j`."""
struct SpatialCylindricalJoint{A,B,I,H}
    name::Symbol
    marker_i::A
    marker_j::B
    inline::I
    hinge::H
end

"""An inline and orient constraint leaving only translation along `z_j`."""
struct SpatialTranslationalJoint{A,B,I,O}
    name::Symbol
    marker_i::A
    marker_j::B
    inline::I
    orient::O
end

function spherical_joint_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:F_x, :reaction, 2),
        VariableDeclaration(:F_y, :reaction, 2),
        VariableDeclaration(:F_z, :reaction, 2),
    ]
    axes = (:x, :y, :z)
    equations = EquationDeclaration[
        [EquationDeclaration(Symbol(:a_, axis), :constraint, 2,
            :spherical_acceleration) for axis in axes];
        [EquationDeclaration(Symbol(:v_, axis), :constraint, 1,
            :spherical_velocity) for axis in axes];
        [EquationDeclaration(Symbol(:p_, axis), :constraint, 0,
            :spherical_position) for axis in axes];
    ]
    families = ConstraintFamilyDeclaration[
        ConstraintFamilyDeclaration(Symbol(:spherical_, axis),
            Symbol(:F_, axis), Symbol(:p_, axis), Symbol(:v_, axis),
            Symbol(:a_, axis)) for axis in axes]
    ComponentRegistration(name, variables,
        [EquationBlockDeclaration(:constraint, equations)], families)
end

function perp_constraint_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:lambda, :reaction, 2),
    ]
    equations = EquationDeclaration[
        EquationDeclaration(:Phi_ddot, :constraint, 2, :perp_acceleration),
        EquationDeclaration(:Phi_dot, :constraint, 1, :perp_velocity),
        EquationDeclaration(:Phi, :constraint, 0, :perp_position),
    ]
    families = ConstraintFamilyDeclaration[
        ConstraintFamilyDeclaration(:perp, :lambda, :Phi, :Phi_dot,
            :Phi_ddot)]
    ComponentRegistration(name, variables,
        [EquationBlockDeclaration(:constraint, equations)], families)
end

function cv_phase_constraint_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:lambda, :reaction, 2),
    ]
    equations = EquationDeclaration[
        EquationDeclaration(:Phi_ddot, :constraint, 2,
            :cv_phase_acceleration),
        EquationDeclaration(:Phi_dot, :constraint, 1,
            :cv_phase_velocity),
        EquationDeclaration(:Phi, :constraint, 0, :cv_phase_position),
    ]
    families = ConstraintFamilyDeclaration[
        ConstraintFamilyDeclaration(:cv_phase, :lambda, :Phi, :Phi_dot,
            :Phi_ddot)]
    ComponentRegistration(name, variables,
        [EquationBlockDeclaration(:constraint, equations)], families)
end

function inplane_constraint_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:lambda, :reaction, 2),
    ]
    equations = EquationDeclaration[
        EquationDeclaration(:Phi_ddot, :constraint, 2,
            :inplane_acceleration),
        EquationDeclaration(:Phi_dot, :constraint, 1,
            :inplane_velocity),
        EquationDeclaration(:Phi, :constraint, 0, :inplane_position),
    ]
    families = ConstraintFamilyDeclaration[
        ConstraintFamilyDeclaration(:inplane, :lambda, :Phi, :Phi_dot,
            :Phi_ddot)]
    ComponentRegistration(name, variables,
        [EquationBlockDeclaration(:constraint, equations)], families)
end

function inline_constraint_registration(name::Symbol;
        translation_coordinates = false)
    variables = VariableDeclaration[
        VariableDeclaration(:lambda_x, :reaction, 2),
        VariableDeclaration(:lambda_y, :reaction, 2),
    ]
    equations = EquationDeclaration[
        EquationDeclaration(:Phi_x_ddot, :constraint, 2,
            :inline_acceleration),
        EquationDeclaration(:Phi_y_ddot, :constraint, 2,
            :inline_acceleration),
        EquationDeclaration(:Phi_x_dot, :constraint, 1, :inline_velocity),
        EquationDeclaration(:Phi_y_dot, :constraint, 1, :inline_velocity),
        EquationDeclaration(:Phi_x, :constraint, 0, :inline_position),
        EquationDeclaration(:Phi_y, :constraint, 0, :inline_position),
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:constraint, equations),
    ]
    if translation_coordinates
        append!(variables, VariableDeclaration[
            VariableDeclaration(:acceleration, :relative_acceleration, 2),
            VariableDeclaration(:velocity, :relative_velocity, 1),
            VariableDeclaration(:distance, :relative_position, 0),
        ])
        append!(blocks, EquationBlockDeclaration[
            EquationBlockDeclaration(:translation_acceleration,
                [EquationDeclaration(:acceleration_definition,
                    :coordinate_relation, 2, :relative_translation)]),
            EquationBlockDeclaration(:translation_velocity,
                [EquationDeclaration(:velocity_definition,
                    :coordinate_relation, 1, :relative_translation)]),
            EquationBlockDeclaration(:translation_position,
                [EquationDeclaration(:distance_definition,
                    :coordinate_relation, 0, :relative_translation)]),
            EquationBlockDeclaration(:selected_translation_state,
                [EquationDeclaration(:acceleration_state, :state_equation, 2,
                     :relative_translation_state),
                 EquationDeclaration(:velocity_state, :state_equation, 1,
                     :relative_translation_state)]),
        ])
    end
    families = ConstraintFamilyDeclaration[
        ConstraintFamilyDeclaration(:inline_x, :lambda_x, :Phi_x,
            :Phi_x_dot, :Phi_x_ddot),
        ConstraintFamilyDeclaration(:inline_y, :lambda_y, :Phi_y,
            :Phi_y_dot, :Phi_y_ddot)]
    ComponentRegistration(name, variables, blocks, families)
end

function hinge_constraint_registration(name::Symbol;
        rotation_coordinates = false)
    variables = VariableDeclaration[
        VariableDeclaration(:lambda_xz, :reaction, 2),
        VariableDeclaration(:lambda_yz, :reaction, 2),
    ]
    equations = EquationDeclaration[
        EquationDeclaration(:Phi_xz_ddot, :constraint, 2,
            :hinge_acceleration),
        EquationDeclaration(:Phi_yz_ddot, :constraint, 2,
            :hinge_acceleration),
        EquationDeclaration(:Phi_xz_dot, :constraint, 1, :hinge_velocity),
        EquationDeclaration(:Phi_yz_dot, :constraint, 1, :hinge_velocity),
        EquationDeclaration(:Phi_xz, :constraint, 0, :hinge_position),
        EquationDeclaration(:Phi_yz, :constraint, 0, :hinge_position),
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:constraint, equations),
    ]
    if rotation_coordinates
        append!(variables, VariableDeclaration[
            VariableDeclaration(:alpha, :relative_acceleration, 2),
            VariableDeclaration(:omega, :relative_velocity, 1),
            VariableDeclaration(:theta, :relative_position, 0),
        ])
        append!(blocks, EquationBlockDeclaration[
            EquationBlockDeclaration(:rotation_acceleration,
                [EquationDeclaration(:alpha_definition,
                    :coordinate_relation, 2, :relative_rotation)]),
            EquationBlockDeclaration(:rotation_velocity,
                [EquationDeclaration(:omega_definition,
                    :coordinate_relation, 1, :relative_rotation)]),
            EquationBlockDeclaration(:rotation_position,
                [EquationDeclaration(:theta_definition,
                    :coordinate_relation, 0, :relative_rotation)]),
            EquationBlockDeclaration(:selected_rotation_state,
                [EquationDeclaration(:alpha_state, :state_equation, 2,
                     :relative_rotation_state),
                 EquationDeclaration(:omega_state, :state_equation, 1,
                     :relative_rotation_state)]),
        ])
    end
    families = ConstraintFamilyDeclaration[
        ConstraintFamilyDeclaration(:hinge_xz, :lambda_xz, :Phi_xz,
            :Phi_xz_dot, :Phi_xz_ddot),
        ConstraintFamilyDeclaration(:hinge_yz, :lambda_yz, :Phi_yz,
            :Phi_yz_dot, :Phi_yz_ddot)]
    ComponentRegistration(name, variables, blocks, families)
end

function orient_constraint_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:lambda_xz, :reaction, 2),
        VariableDeclaration(:lambda_yz, :reaction, 2),
        VariableDeclaration(:lambda_xy, :reaction, 2),
    ]
    equations = EquationDeclaration[
        EquationDeclaration(:Phi_xz_ddot, :constraint, 2,
            :orient_acceleration),
        EquationDeclaration(:Phi_yz_ddot, :constraint, 2,
            :orient_acceleration),
        EquationDeclaration(:Phi_xy_ddot, :constraint, 2,
            :orient_acceleration),
        EquationDeclaration(:Phi_xz_dot, :constraint, 1,
            :orient_velocity),
        EquationDeclaration(:Phi_yz_dot, :constraint, 1,
            :orient_velocity),
        EquationDeclaration(:Phi_xy_dot, :constraint, 1,
            :orient_velocity),
        EquationDeclaration(:Phi_xz, :constraint, 0,
            :orient_position),
        EquationDeclaration(:Phi_yz, :constraint, 0,
            :orient_position),
        EquationDeclaration(:Phi_xy, :constraint, 0,
            :orient_position),
    ]
    families = ConstraintFamilyDeclaration[
        ConstraintFamilyDeclaration(:orient_xz, :lambda_xz, :Phi_xz,
            :Phi_xz_dot, :Phi_xz_ddot),
        ConstraintFamilyDeclaration(:orient_yz, :lambda_yz, :Phi_yz,
            :Phi_yz_dot, :Phi_yz_ddot),
        ConstraintFamilyDeclaration(:orient_xy, :lambda_xy, :Phi_xy,
            :Phi_xy_dot, :Phi_xy_ddot)]
    ComponentRegistration(name, variables,
        [EquationBlockDeclaration(:constraint, equations)], families)
end

function revolute_joint_registration(name::Symbol;
        rotation_coordinates = false)
    spherical = spherical_joint_registration(name)
    hinge = hinge_constraint_registration(name; rotation_coordinates)
    variables = [spherical.variables; hinge.variables]
    constraint_equations = [
        spherical.equation_blocks[:constraint].equations;
        hinge.equation_blocks[:constraint].equations;
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:constraint, constraint_equations),
    ]
    if rotation_coordinates
        for block in (:rotation_acceleration, :rotation_velocity,
                :rotation_position, :selected_rotation_state)
            push!(blocks, hinge.equation_blocks[block])
        end
    end
    ComponentRegistration(name, variables, blocks,
        [spherical.constraint_families; hinge.constraint_families])
end

function fixed_joint_registration(name::Symbol)
    spherical = spherical_joint_registration(name)
    orient = orient_constraint_registration(name)
    variables = [spherical.variables; orient.variables]
    constraint_equations = [
        spherical.equation_blocks[:constraint].equations;
        orient.equation_blocks[:constraint].equations;
    ]
    ComponentRegistration(name, variables,
        [EquationBlockDeclaration(:constraint, constraint_equations)],
        [spherical.constraint_families; orient.constraint_families])
end

function constant_velocity_joint_registration(name::Symbol)
    spherical = spherical_joint_registration(name)
    phase = cv_phase_constraint_registration(name)
    constraint_equations = [
        spherical.equation_blocks[:constraint].equations;
        phase.equation_blocks[:constraint].equations;
    ]
    ComponentRegistration(name, [spherical.variables; phase.variables],
        [EquationBlockDeclaration(:constraint, constraint_equations)],
        [spherical.constraint_families; phase.constraint_families])
end

function cylindrical_joint_registration(name::Symbol;
        translation_coordinates = false, rotation_coordinates = false)
    inline = inline_constraint_registration(name; translation_coordinates)
    hinge = hinge_constraint_registration(name; rotation_coordinates)
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:constraint, [
            inline.equation_blocks[:constraint].equations;
            hinge.equation_blocks[:constraint].equations;
        ]),
    ]
    if translation_coordinates
        for block in (:translation_acceleration, :translation_velocity,
                :translation_position, :selected_translation_state)
            push!(blocks, inline.equation_blocks[block])
        end
    end
    if rotation_coordinates
        for block in (:rotation_acceleration, :rotation_velocity,
                :rotation_position, :selected_rotation_state)
            push!(blocks, hinge.equation_blocks[block])
        end
    end
    ComponentRegistration(name, [inline.variables; hinge.variables], blocks,
        [inline.constraint_families; hinge.constraint_families])
end

function translational_joint_registration(name::Symbol;
        translation_coordinates = false)
    inline = inline_constraint_registration(name; translation_coordinates)
    orient = orient_constraint_registration(name)
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:constraint, [
            inline.equation_blocks[:constraint].equations;
            orient.equation_blocks[:constraint].equations;
        ]),
    ]
    if translation_coordinates
        for block in (:translation_acceleration, :translation_velocity,
                :translation_position, :selected_translation_state)
            push!(blocks, inline.equation_blocks[block])
        end
    end
    ComponentRegistration(name, [inline.variables; orient.variables], blocks,
        [inline.constraint_families; orient.constraint_families])
end

component_registration(joint::SpatialSphericalJoint) =
    spherical_joint_registration(joint.name)
component_registration(constraint::SpatialPerpConstraint) =
    perp_constraint_registration(constraint.name)
component_registration(constraint::SpatialCVPhaseConstraint) =
    cv_phase_constraint_registration(constraint.name)
component_registration(constraint::SpatialInplaneConstraint) =
    inplane_constraint_registration(constraint.name)
component_registration(constraint::SpatialInlineConstraint) =
    inline_constraint_registration(constraint.name;
        translation_coordinates = !isempty(constraint.translation_variables))
component_registration(constraint::SpatialHingeConstraint) =
    hinge_constraint_registration(constraint.name;
        rotation_coordinates = !isempty(constraint.rotation_variables))
component_registration(constraint::SpatialOrientConstraint) =
    orient_constraint_registration(constraint.name)
component_registration(joint::SpatialRevoluteJoint) =
    revolute_joint_registration(joint.name;
        rotation_coordinates = !isempty(joint.hinge.rotation_variables))
component_registration(joint::SpatialFixedJoint) =
    fixed_joint_registration(joint.name)
component_registration(joint::SpatialConstantVelocityJoint) =
    constant_velocity_joint_registration(joint.name)
component_registration(joint::SpatialCylindricalJoint) =
    cylindrical_joint_registration(joint.name;
        translation_coordinates = !isempty(joint.inline.translation_variables),
        rotation_coordinates = !isempty(joint.hinge.rotation_variables))
component_registration(joint::SpatialTranslationalJoint) =
    translational_joint_registration(joint.name;
        translation_coordinates = !isempty(joint.inline.translation_variables))

function allocated_spherical_joint(layout, name, marker_a, marker_b)
    variables = component_variable_indices(layout, name)
    equations = component_equation_indices(layout, name, :constraint)
    SpatialSphericalJoint(name, marker_a, marker_b, variables,
        equations[1:3], equations[4:6], equations[7:9])
end

function allocated_perp_constraint(layout, name, marker_i, marker_j)
    variable = only(component_variable_indices(layout, name))
    equations = component_equation_indices(layout, name, :constraint)
    SpatialPerpConstraint(name, marker_i, marker_j, 1, 2, variable,
        equations[1], equations[2], equations[3])
end

function allocated_cv_phase_constraint(layout, name, marker_i, marker_j)
    variable = only(component_variable_indices(layout, name))
    equations = component_equation_indices(layout, name, :constraint)
    SpatialCVPhaseConstraint(name, marker_i, marker_j, variable,
        equations[1], equations[2], equations[3])
end

function allocated_inplane_constraint(layout, name, marker_i, marker_j;
        normal_axis = 3)
    normal_axis in 1:3 || throw(ArgumentError(
        "an inplane normal axis must be 1, 2, or 3"))
    variable = only(component_variable_indices(layout, name))
    equations = component_equation_indices(layout, name, :constraint)
    axis = SpatialDirectedAxis(marker_j, normal_axis)
    geometry = SpatialDirectedDistance(marker_i, marker_j, axis)
    SpatialInplaneConstraint(name, geometry, variable,
        equations[1], equations[2], equations[3])
end

function allocated_inline_constraint(layout, name, marker_i, marker_j;
        translation_coordinates = false)
    variables = component_variable_indices(layout, name)
    equations = component_equation_indices(layout, name, :constraint)
    geometry_x = SpatialDirectedDistance(
        marker_i, marker_j, SpatialDirectedAxis(marker_j, 1))
    geometry_y = SpatialDirectedDistance(
        marker_i, marker_j, SpatialDirectedAxis(marker_j, 2))
    inplane_x = SpatialInplaneConstraint(Symbol(name, ".x"), geometry_x,
        variables[1], equations[1], equations[3], equations[5])
    inplane_y = SpatialInplaneConstraint(Symbol(name, ".y"), geometry_y,
        variables[2], equations[2], equations[4], equations[6])
    translation_variables = translation_coordinates ?
        collect(variables[3:5]) : Int[]
    translation_equations = translation_coordinates ? [
        only(component_equation_indices(layout, name,
            :translation_acceleration)),
        only(component_equation_indices(layout, name,
            :translation_velocity)),
        only(component_equation_indices(layout, name,
            :translation_position)),
    ] : Int[]
    state_equations = translation_coordinates ? collect(
        component_equation_indices(layout, name,
            :selected_translation_state)) : Int[]
    axial_geometry = SpatialDirectedDistance(
        marker_i, marker_j, SpatialDirectedAxis(marker_j, 3))
    SpatialInlineConstraint(name, marker_i, marker_j, inplane_x, inplane_y,
        axial_geometry, translation_variables, translation_equations,
        state_equations)
end

function allocated_hinge_constraint(layout, name, marker_i, marker_j;
        rotation_coordinates = false)
    variables = component_variable_indices(layout, name)
    equations = component_equation_indices(layout, name, :constraint)
    perp_xz = SpatialPerpConstraint(Symbol(name, ".xz"), marker_i, marker_j,
        1, 3, variables[1], equations[1], equations[3], equations[5])
    perp_yz = SpatialPerpConstraint(Symbol(name, ".yz"), marker_i, marker_j,
        2, 3, variables[2], equations[2], equations[4], equations[6])
    rotation_variables = rotation_coordinates ? collect(variables[3:5]) : Int[]
    rotation_equations = rotation_coordinates ? [
        only(component_equation_indices(layout, name,
            :rotation_acceleration)),
        only(component_equation_indices(layout, name, :rotation_velocity)),
        only(component_equation_indices(layout, name, :rotation_position)),
    ] : Int[]
    state_equations = rotation_coordinates ? collect(
        component_equation_indices(layout, name,
            :selected_rotation_state)) : Int[]
    SpatialHingeConstraint(name, marker_i, marker_j, perp_xz, perp_yz,
        rotation_variables, rotation_equations, state_equations)
end

function allocated_orient_constraint(layout, name, marker_i, marker_j)
    variables = component_variable_indices(layout, name)
    equations = component_equation_indices(layout, name, :constraint)
    perp_xz = SpatialPerpConstraint(Symbol(name, ".xz"), marker_i, marker_j,
        1, 3, variables[1], equations[1], equations[4], equations[7])
    perp_yz = SpatialPerpConstraint(Symbol(name, ".yz"), marker_i, marker_j,
        2, 3, variables[2], equations[2], equations[5], equations[8])
    hinge = SpatialHingeConstraint(name, marker_i, marker_j,
        perp_xz, perp_yz, Int[], Int[], Int[])
    perp_xy = SpatialPerpConstraint(Symbol(name, ".xy"), marker_i, marker_j,
        1, 2, variables[3], equations[3], equations[6], equations[9])
    SpatialOrientConstraint(name, marker_i, marker_j, hinge, perp_xy)
end

function allocated_revolute_joint(layout, name, marker_a, marker_b;
        rotation_coordinates = false)
    variables = component_variable_indices(layout, name)
    equations = component_equation_indices(layout, name, :constraint)
    spherical = SpatialSphericalJoint(name, marker_a, marker_b,
        variables[1:3], equations[1:3], equations[4:6], equations[7:9])
    hinge_equations = equations[10:15]
    perp_xz = SpatialPerpConstraint(Symbol(name, ".xz"), marker_a, marker_b,
        1, 3, variables[4], hinge_equations[1], hinge_equations[3],
        hinge_equations[5])
    perp_yz = SpatialPerpConstraint(Symbol(name, ".yz"), marker_a, marker_b,
        2, 3, variables[5], hinge_equations[2], hinge_equations[4],
        hinge_equations[6])
    rotation_variables = rotation_coordinates ? collect(variables[6:8]) : Int[]
    rotation_equations = rotation_coordinates ? [
        only(component_equation_indices(layout, name,
            :rotation_acceleration)),
        only(component_equation_indices(layout, name, :rotation_velocity)),
        only(component_equation_indices(layout, name, :rotation_position)),
    ] : Int[]
    state_equations = rotation_coordinates ? collect(
        component_equation_indices(layout, name,
            :selected_rotation_state)) : Int[]
    hinge = SpatialHingeConstraint(name, marker_a, marker_b,
        perp_xz, perp_yz, rotation_variables, rotation_equations,
        state_equations)
    SpatialRevoluteJoint(name, marker_a, marker_b, spherical, hinge)
end

function allocated_fixed_joint(layout, name, marker_a, marker_b)
    variables = component_variable_indices(layout, name)
    equations = component_equation_indices(layout, name, :constraint)
    spherical = SpatialSphericalJoint(name, marker_a, marker_b,
        variables[1:3], equations[1:3], equations[4:6], equations[7:9])
    orient_equations = equations[10:18]
    perp_xz = SpatialPerpConstraint(Symbol(name, ".xz"), marker_a, marker_b,
        1, 3, variables[4], orient_equations[1], orient_equations[4],
        orient_equations[7])
    perp_yz = SpatialPerpConstraint(Symbol(name, ".yz"), marker_a, marker_b,
        2, 3, variables[5], orient_equations[2], orient_equations[5],
        orient_equations[8])
    hinge = SpatialHingeConstraint(name, marker_a, marker_b,
        perp_xz, perp_yz, Int[], Int[], Int[])
    perp_xy = SpatialPerpConstraint(Symbol(name, ".xy"), marker_a, marker_b,
        1, 2, variables[6], orient_equations[3], orient_equations[6],
        orient_equations[9])
    orient = SpatialOrientConstraint(name, marker_a, marker_b, hinge, perp_xy)
    SpatialFixedJoint(name, marker_a, marker_b, spherical, orient)
end

function allocated_constant_velocity_joint(layout, name, marker_i, marker_j)
    variables = component_variable_indices(layout, name)
    equations = component_equation_indices(layout, name, :constraint)
    spherical = SpatialSphericalJoint(name, marker_i, marker_j,
        variables[1:3], equations[1:3], equations[4:6], equations[7:9])
    phase = SpatialCVPhaseConstraint(name, marker_i, marker_j, variables[4],
        equations[10], equations[11], equations[12])
    SpatialConstantVelocityJoint(name, marker_i, marker_j, spherical, phase)
end

function allocated_cylindrical_joint(layout, name, marker_i, marker_j;
        translation_coordinates = false, rotation_coordinates = false)
    variables = component_variable_indices(layout, name)
    equations = component_equation_indices(layout, name, :constraint)

    geometry_x = SpatialDirectedDistance(
        marker_i, marker_j, SpatialDirectedAxis(marker_j, 1))
    geometry_y = SpatialDirectedDistance(
        marker_i, marker_j, SpatialDirectedAxis(marker_j, 2))
    inplane_x = SpatialInplaneConstraint(Symbol(name, ".x"), geometry_x,
        variables[1], equations[1], equations[3], equations[5])
    inplane_y = SpatialInplaneConstraint(Symbol(name, ".y"), geometry_y,
        variables[2], equations[2], equations[4], equations[6])
    translation_variables = translation_coordinates ?
        collect(variables[3:5]) : Int[]
    translation_equations = translation_coordinates ? [
        only(component_equation_indices(layout, name,
            :translation_acceleration)),
        only(component_equation_indices(layout, name, :translation_velocity)),
        only(component_equation_indices(layout, name, :translation_position)),
    ] : Int[]
    translation_state_equations = translation_coordinates ? collect(
        component_equation_indices(layout, name,
            :selected_translation_state)) : Int[]
    axial_geometry = SpatialDirectedDistance(
        marker_i, marker_j, SpatialDirectedAxis(marker_j, 3))
    inline = SpatialInlineConstraint(name, marker_i, marker_j,
        inplane_x, inplane_y, axial_geometry, translation_variables,
        translation_equations, translation_state_equations)

    hinge_variable_start = translation_coordinates ? 6 : 3
    hinge_equations = equations[7:12]
    perp_xz = SpatialPerpConstraint(Symbol(name, ".xz"), marker_i, marker_j,
        1, 3, variables[hinge_variable_start], hinge_equations[1],
        hinge_equations[3], hinge_equations[5])
    perp_yz = SpatialPerpConstraint(Symbol(name, ".yz"), marker_i, marker_j,
        2, 3, variables[hinge_variable_start + 1], hinge_equations[2],
        hinge_equations[4], hinge_equations[6])
    rotation_variables = rotation_coordinates ? collect(
        variables[(hinge_variable_start + 2):(hinge_variable_start + 4)]) : Int[]
    rotation_equations = rotation_coordinates ? [
        only(component_equation_indices(layout, name,
            :rotation_acceleration)),
        only(component_equation_indices(layout, name, :rotation_velocity)),
        only(component_equation_indices(layout, name, :rotation_position)),
    ] : Int[]
    rotation_state_equations = rotation_coordinates ? collect(
        component_equation_indices(layout, name,
            :selected_rotation_state)) : Int[]
    hinge = SpatialHingeConstraint(name, marker_i, marker_j,
        perp_xz, perp_yz, rotation_variables, rotation_equations,
        rotation_state_equations)
    SpatialCylindricalJoint(name, marker_i, marker_j, inline, hinge)
end

function allocated_translational_joint(layout, name, marker_i, marker_j;
        translation_coordinates = false)
    variables = component_variable_indices(layout, name)
    equations = component_equation_indices(layout, name, :constraint)

    geometry_x = SpatialDirectedDistance(
        marker_i, marker_j, SpatialDirectedAxis(marker_j, 1))
    geometry_y = SpatialDirectedDistance(
        marker_i, marker_j, SpatialDirectedAxis(marker_j, 2))
    inplane_x = SpatialInplaneConstraint(Symbol(name, ".x"), geometry_x,
        variables[1], equations[1], equations[3], equations[5])
    inplane_y = SpatialInplaneConstraint(Symbol(name, ".y"), geometry_y,
        variables[2], equations[2], equations[4], equations[6])
    translation_variables = translation_coordinates ?
        collect(variables[3:5]) : Int[]
    translation_equations = translation_coordinates ? [
        only(component_equation_indices(layout, name,
            :translation_acceleration)),
        only(component_equation_indices(layout, name, :translation_velocity)),
        only(component_equation_indices(layout, name, :translation_position)),
    ] : Int[]
    translation_state_equations = translation_coordinates ? collect(
        component_equation_indices(layout, name,
            :selected_translation_state)) : Int[]
    axial_geometry = SpatialDirectedDistance(
        marker_i, marker_j, SpatialDirectedAxis(marker_j, 3))
    inline = SpatialInlineConstraint(name, marker_i, marker_j,
        inplane_x, inplane_y, axial_geometry, translation_variables,
        translation_equations, translation_state_equations)

    orient_variable_start = translation_coordinates ? 6 : 3
    orient_equations = equations[7:15]
    perp_xz = SpatialPerpConstraint(Symbol(name, ".xz"), marker_i, marker_j,
        1, 3, variables[orient_variable_start], orient_equations[1],
        orient_equations[4], orient_equations[7])
    perp_yz = SpatialPerpConstraint(Symbol(name, ".yz"), marker_i, marker_j,
        2, 3, variables[orient_variable_start + 1], orient_equations[2],
        orient_equations[5], orient_equations[8])
    hinge = SpatialHingeConstraint(name, marker_i, marker_j,
        perp_xz, perp_yz, Int[], Int[], Int[])
    perp_xy = SpatialPerpConstraint(Symbol(name, ".xy"), marker_i, marker_j,
        1, 2, variables[orient_variable_start + 2], orient_equations[3],
        orient_equations[6], orient_equations[9])
    orient = SpatialOrientConstraint(name, marker_i, marker_j, hinge, perp_xy)
    SpatialTranslationalJoint(name, marker_i, marker_j, inline, orient)
end

connection_hinge(constraint::SpatialHingeConstraint) = constraint
connection_hinge(joint::SpatialRevoluteJoint) = joint.hinge
connection_hinge(joint::SpatialCylindricalJoint) = joint.hinge
connection_hinge(::Any) = nothing
connection_inline(constraint::SpatialInlineConstraint) = constraint
connection_inline(joint::SpatialCylindricalJoint) = joint.inline
connection_inline(joint::SpatialTranslationalJoint) = joint.inline
connection_inline(::Any) = nothing

inline_axis(constraint::SpatialInlineConstraint, z) =
    directed_distance_direction(constraint.axial_geometry, z)
inline_distance(constraint::SpatialInlineConstraint, z) =
    directed_distance_position(constraint.axial_geometry, z)
inline_velocity(constraint::SpatialInlineConstraint, z) =
    directed_distance_velocity(constraint.axial_geometry, z)
inline_acceleration(constraint::SpatialInlineConstraint, z) =
    directed_distance_acceleration(constraint.axial_geometry, z)

inplane_normal(constraint::SpatialInplaneConstraint, z) =
    directed_distance_direction(constraint.geometry, z)

inplane_position(constraint::SpatialInplaneConstraint, z) =
    directed_distance_position(constraint.geometry, z)
inplane_velocity(constraint::SpatialInplaneConstraint, z) =
    directed_distance_velocity(constraint.geometry, z)
inplane_acceleration(constraint::SpatialInplaneConstraint, z) =
    directed_distance_acceleration(constraint.geometry, z)

function perp_kinematics(constraint::SpatialPerpConstraint, z)
    first = marker_axis_kinematics(
        constraint.marker_i, z, constraint.axis_i)
    second = marker_axis_kinematics(
        constraint.marker_j, z, constraint.axis_j)
    (; first, second)
end

function perp_normal(constraint::SpatialPerpConstraint, z)
    kinematics = perp_kinematics(constraint, z)
    cross(kinematics.first.direction, kinematics.second.direction)
end

function perp_position(constraint::SpatialPerpConstraint, z)
    kinematics = perp_kinematics(constraint, z)
    dot(kinematics.first.direction, kinematics.second.direction)
end

function perp_velocity(constraint::SpatialPerpConstraint, z)
    kinematics = perp_kinematics(constraint, z)
    first, second = kinematics.first, kinematics.second
    dot(first.velocity, second.direction) +
        dot(first.direction, second.velocity)
end

function perp_acceleration(constraint::SpatialPerpConstraint, z)
    kinematics = perp_kinematics(constraint, z)
    first, second = kinematics.first, kinematics.second
    dot(first.acceleration, second.direction) +
        2dot(first.velocity, second.velocity) +
        dot(first.direction, second.acceleration)
end

function cv_phase_kinematics(constraint::SpatialCVPhaseConstraint, z)
    first_x = marker_axis_kinematics(constraint.marker_i, z, 1)
    first_z = marker_axis_kinematics(constraint.marker_i, z, 3)
    second_x = marker_axis_kinematics(constraint.marker_j, z, 1)
    second_z = marker_axis_kinematics(constraint.marker_j, z, 3)
    axis_sum = first_z.direction + second_z.direction
    dot(axis_sum, axis_sum) > 1.0e-12 || throw(ArgumentError(
        "cv_phase constraint '$(constraint.name)' shaft axes are nearly opposite"))
    phase_cross = cross(first_x.direction, second_x.direction)
    axis_sum_velocity = first_z.velocity + second_z.velocity
    axis_sum_magnitude = norm(axis_sum)
    bisector = axis_sum ./ axis_sum_magnitude
    bisector_velocity = (axis_sum_velocity .-
        bisector .* dot(bisector, axis_sum_velocity)) ./ axis_sum_magnitude
    (; first_x, first_z, second_x, second_z, axis_sum, phase_cross,
       bisector, bisector_velocity)
end

function cv_phase_position(constraint::SpatialCVPhaseConstraint, z)
    k = cv_phase_kinematics(constraint, z)
    dot(k.axis_sum, k.phase_cross)
end

function cv_phase_velocity(constraint::SpatialCVPhaseConstraint, z)
    k = cv_phase_kinematics(constraint, z)
    first = marker_angular_kinematics(constraint.marker_i, z)
    second = marker_angular_kinematics(constraint.marker_j, z)
    dot(first.omega - second.omega, k.bisector)
end

function cv_phase_acceleration(constraint::SpatialCVPhaseConstraint, z)
    k = cv_phase_kinematics(constraint, z)
    first = marker_angular_kinematics(constraint.marker_i, z)
    second = marker_angular_kinematics(constraint.marker_j, z)
    dot(first.alpha - second.alpha, k.bisector) +
        dot(first.omega - second.omega, k.bisector_velocity)
end

function cv_phase_reaction_directions(constraint::SpatialCVPhaseConstraint, z)
    k = cv_phase_kinematics(constraint, z)
    (; first = k.bisector, second = -k.bisector)
end

function marker_angular_kinematics(marker::SpatialGroundMarker, z)
    zero_vector = zeros(eltype(z), 3)
    (; omega = zero_vector, alpha = copy(zero_vector), body = nothing,
       orientation = nothing, omega_parameters = nothing,
       alpha_parameters = nothing)
end

function marker_angular_kinematics(marker::SpatialBodyMarker, z)
    body = marker.body
    parameters = @view z[body.euler_parameter_variables]
    orientation = rotation_matrix(parameters)
    omega_body = @view z[body.angular_velocity_variables]
    alpha_body = @view z[body.angular_acceleration_variables]
    omega = orientation * omega_body
    alpha = orientation * alpha_body
    omega_parameters = rotation_vector_jacobian(parameters, omega_body)
    alpha_parameters = rotation_vector_jacobian(parameters, alpha_body)
    (; omega, alpha, body, orientation, omega_parameters, alpha_parameters)
end

function marker_angular_kinematics(marker::SpatialFlexibleBeamMarker, z)
    body = marker.body
    parameters = @view z[body.euler_parameter_variables]
    orientation = rotation_matrix(parameters)
    omega_body = @view z[body.angular_velocity_variables]
    alpha_body = @view z[body.angular_acceleration_variables]
    elastic_rate = marker.orientation_shape *
        (@view z[body.elastic_velocity_variables])
    elastic_acceleration = marker.orientation_shape *
        (@view z[body.elastic_acceleration_variables])
    omega = orientation * (omega_body + elastic_rate)
    alpha = orientation * (alpha_body + elastic_acceleration +
        cross(omega_body, elastic_rate))
    (; omega, alpha, body, orientation,
       omega_parameters = nothing, alpha_parameters = nothing)
end

function add_ad_jacobian!(jacobian, z, rows, columns, values)
    columns = sort!(unique!(collect(columns)))
    (isempty(rows) || isempty(columns)) && return nothing
    initial = collect(z[columns])
    block = ForwardDiff.jacobian(initial) do local_values
        local_z = z .+ zero(eltype(local_values))
        local_z[columns] .= local_values
        collect(values(local_z))
    end
    jacobian[rows, columns] .+= block
    nothing
end

constraint_dependencies(markers...) = sort!(unique!(reduce(vcat,
    (spatial_marker_dependency_indices(marker) for marker in markers);
    init = Int[])))

function hinge_angle_geometry(constraint::SpatialHingeConstraint, z)
    first_x = marker_axis_kinematics(constraint.marker_i, z, 1)
    second_x = marker_axis_kinematics(constraint.marker_j, z, 1)
    second_z = marker_axis_kinematics(constraint.marker_j, z, 3)
    cosine = dot(second_x.direction, first_x.direction)
    sine = dot(second_z.direction,
        cross(second_x.direction, first_x.direction))
    (; first_x, second_x, second_z, cosine, sine)
end

"""Wrapped physical hinge angle, positive from `x_j` toward `x_i` about `z_j`."""
function hinge_angle(constraint::SpatialHingeConstraint, z)
    geometry = hinge_angle_geometry(constraint, z)
    atan(geometry.sine, geometry.cosine)
end

function hinge_angular_velocity(constraint::SpatialHingeConstraint, z)
    first = marker_angular_kinematics(constraint.marker_i, z)
    second = marker_angular_kinematics(constraint.marker_j, z)
    axis = marker_axis_kinematics(constraint.marker_j, z, 3)
    dot(first.omega - second.omega, axis.direction)
end

function hinge_angular_acceleration(constraint::SpatialHingeConstraint, z)
    first = marker_angular_kinematics(constraint.marker_i, z)
    second = marker_angular_kinematics(constraint.marker_j, z)
    axis = marker_axis_kinematics(constraint.marker_j, z, 3)
    dot(first.alpha - second.alpha, axis.direction) +
        dot(first.omega - second.omega, axis.velocity)
end

function initialize_hinge_coordinates!(z, constraint::SpatialHingeConstraint,
        level)
    isempty(constraint.rotation_variables) && return nothing
    alpha, omega, theta = constraint.rotation_variables
    level >= 0 && (z[theta] = hinge_angle(constraint, z))
    level >= 1 && (z[omega] = hinge_angular_velocity(constraint, z))
    level >= 2 && (z[alpha] = hinge_angular_acceleration(constraint, z))
    nothing
end


initialize_hinge_coordinates!(z, joint::SpatialRevoluteJoint, level) =
    initialize_hinge_coordinates!(z, joint.hinge, level)
initialize_hinge_coordinates!(z, joint::SpatialCylindricalJoint, level) =
    initialize_hinge_coordinates!(z, joint.hinge, level)

function initialize_inline_coordinates!(z,
        constraint::SpatialInlineConstraint, level)
    isempty(constraint.translation_variables) && return nothing
    acceleration, velocity, distance = constraint.translation_variables
    level >= 0 && (z[distance] = inline_distance(constraint, z))
    level >= 1 && (z[velocity] = inline_velocity(constraint, z))
    level >= 2 && (z[acceleration] = inline_acceleration(constraint, z))
    nothing
end

initialize_inline_coordinates!(z, joint::SpatialCylindricalJoint, level) =
    initialize_inline_coordinates!(z, joint.inline, level)
initialize_inline_coordinates!(z, joint::SpatialTranslationalJoint, level) =
    initialize_inline_coordinates!(z, joint.inline, level)

function perp_constraints!(equations, z, constraint)
    equations[constraint.acceleration_equation] =
        perp_acceleration(constraint, z)
    equations[constraint.velocity_equation] = perp_velocity(constraint, z)
    equations[constraint.position_equation] = perp_position(constraint, z)
    nothing
end

function cv_phase_constraints!(equations, z, constraint)
    equations[constraint.acceleration_equation] =
        cv_phase_acceleration(constraint, z)
    equations[constraint.velocity_equation] = cv_phase_velocity(constraint, z)
    equations[constraint.position_equation] = cv_phase_position(constraint, z)
    nothing
end

function cv_phase_constraints_jacobian!(jacobian, z, constraint)
    rows = [constraint.acceleration_equation,
            constraint.velocity_equation,
            constraint.position_equation]
    values = local_z -> [cv_phase_acceleration(constraint, local_z),
                         cv_phase_velocity(constraint, local_z),
                         cv_phase_position(constraint, local_z)]
    add_ad_jacobian!(jacobian, z, rows,
        constraint_dependencies(constraint.marker_i, constraint.marker_j),
        values)
    nothing
end

function add_perp_axis_jacobian!(jacobian, constraint, kinematics, side)
    item = side == :first ? kinematics.first : kinematics.second
    isnothing(item.body) && return nothing
    other = side == :first ? kinematics.second : kinematics.first
    body = item.body
    parameter_columns = body.euler_parameter_variables
    omega_columns = body.angular_velocity_variables
    alpha_columns = body.angular_acceleration_variables
    position_row = constraint.position_equation
    velocity_row = constraint.velocity_equation
    acceleration_row = constraint.acceleration_equation

    jacobian[position_row, parameter_columns] .+=
        transpose(item.direction_parameters) * other.direction
    jacobian[velocity_row, parameter_columns] .+=
        transpose(item.velocity_parameters) * other.direction .+
        transpose(item.direction_parameters) * other.velocity
    jacobian[velocity_row, omega_columns] .+=
        transpose(item.velocity_omega) * other.direction
    jacobian[acceleration_row, parameter_columns] .+=
        transpose(item.acceleration_parameters) * other.direction .+
        2transpose(item.velocity_parameters) * other.velocity .+
        transpose(item.direction_parameters) * other.acceleration
    jacobian[acceleration_row, omega_columns] .+=
        transpose(item.acceleration_omega) * other.direction .+
        2transpose(item.velocity_omega) * other.velocity
    jacobian[acceleration_row, alpha_columns] .+=
        transpose(item.acceleration_alpha) * other.direction
    nothing
end

function perp_constraints_jacobian!(jacobian, z, constraint)
    if is_flexible_marker(constraint.marker_i) ||
            is_flexible_marker(constraint.marker_j)
        rows = [constraint.acceleration_equation,
                constraint.velocity_equation,
                constraint.position_equation]
        values = local_z -> [perp_acceleration(constraint, local_z),
                             perp_velocity(constraint, local_z),
                             perp_position(constraint, local_z)]
        add_ad_jacobian!(jacobian, z, rows,
            constraint_dependencies(
                constraint.marker_i, constraint.marker_j), values)
        return nothing
    end
    kinematics = perp_kinematics(constraint, z)
    add_perp_axis_jacobian!(jacobian, constraint, kinematics, :first)
    add_perp_axis_jacobian!(jacobian, constraint, kinematics, :second)
    nothing
end

function inplane_constraints!(equations, z, constraint)
    values = directed_distance_values(constraint.geometry, z)
    equations[constraint.acceleration_equation] =
        values.acceleration
    equations[constraint.velocity_equation] =
        values.velocity
    equations[constraint.position_equation] =
        values.position
    nothing
end

function inplane_constraints_jacobian!(jacobian, z, constraint,
        multiplier = 1)
    directed_distance_jacobian!(jacobian, z, constraint.geometry,
        constraint.acceleration_equation, constraint.velocity_equation,
        constraint.position_equation, multiplier)
end

function spherical_constraints!(equations, z, joint)
    equations[joint.acceleration_equations] .=
        spatial_marker_acceleration(joint.marker_a, z) .-
        spatial_marker_acceleration(joint.marker_b, z)
    equations[joint.velocity_equations] .=
        spatial_marker_velocity(joint.marker_a, z) .-
        spatial_marker_velocity(joint.marker_b, z)
    equations[joint.position_equations] .=
        spatial_marker_position(joint.marker_a, z) .-
        spatial_marker_position(joint.marker_b, z)
    nothing
end

function add_marker_constraint_jacobian!(jacobian, z, marker, sign,
        acceleration_rows, velocity_rows, position_rows)
    marker isa SpatialGroundMarker && return nothing
    body = marker.body
    parameters = @view z[body.euler_parameter_variables]
    orientation = rotation_matrix(parameters)
    offset = marker.position_body
    omega = @view z[body.angular_velocity_variables]
    alpha = @view z[body.angular_acceleration_variables]
    local_velocity = cross(omega, offset)
    local_acceleration = cross(alpha, offset) .+
        cross(omega, local_velocity)
    identity = Matrix{eltype(z)}(I, 3, 3)

    jacobian[position_rows, body.position_variables] .+= sign .* identity
    jacobian[position_rows, body.euler_parameter_variables] .+=
        sign .* rotation_vector_jacobian(parameters, offset)

    jacobian[velocity_rows, body.velocity_variables] .+= sign .* identity
    jacobian[velocity_rows, body.angular_velocity_variables] .+=
        -sign .* orientation * skew(offset)
    jacobian[velocity_rows, body.euler_parameter_variables] .+=
        sign .* rotation_vector_jacobian(parameters, local_velocity)

    jacobian[acceleration_rows, body.acceleration_variables] .+=
        sign .* identity
    jacobian[acceleration_rows,
        body.angular_acceleration_variables] .+=
        -sign .* orientation * skew(offset)
    omega_derivative = -skew(local_velocity) .-
        skew(omega) * skew(offset)
    jacobian[acceleration_rows, body.angular_velocity_variables] .+=
        sign .* orientation * omega_derivative
    jacobian[acceleration_rows, body.euler_parameter_variables] .+=
        sign .* rotation_vector_jacobian(parameters, local_acceleration)
    nothing
end

function spherical_constraints_jacobian!(jacobian, z, joint)
    if is_flexible_marker(joint.marker_a) || is_flexible_marker(joint.marker_b)
        rows = [collect(joint.acceleration_equations);
                collect(joint.velocity_equations);
                collect(joint.position_equations)]
        values = local_z -> [
            spatial_marker_acceleration(joint.marker_a, local_z) .-
                spatial_marker_acceleration(joint.marker_b, local_z);
            spatial_marker_velocity(joint.marker_a, local_z) .-
                spatial_marker_velocity(joint.marker_b, local_z);
            spatial_marker_position(joint.marker_a, local_z) .-
                spatial_marker_position(joint.marker_b, local_z)]
        add_ad_jacobian!(jacobian, z, rows,
            constraint_dependencies(joint.marker_a, joint.marker_b), values)
        return nothing
    end
    for (marker, sign) in ((joint.marker_a, 1), (joint.marker_b, -1))
        add_marker_constraint_jacobian!(jacobian, z, marker, sign,
            joint.acceleration_equations, joint.velocity_equations,
            joint.position_equations)
    end
    nothing
end

function executable_blocks(joint::SpatialSphericalJoint)
    rows = collect(joint.acceleration_equations.start:
        joint.position_equations.stop)
    ExecutableEquationBlock[
        ExecutableEquationBlock(joint.name, :constraint, rows,
            (e, t, z, zd) -> spherical_constraints!(e, z, joint),
            (J, t, z, zd, c) ->
                spherical_constraints_jacobian!(J, z, joint)),
    ]
end

function executable_blocks(constraint::SpatialPerpConstraint)
    rows = [constraint.acceleration_equation,
            constraint.velocity_equation,
            constraint.position_equation]
    ExecutableEquationBlock[
        ExecutableEquationBlock(constraint.name, :constraint, rows,
            (e, t, z, zd) -> perp_constraints!(e, z, constraint),
            (J, t, z, zd, c) ->
                perp_constraints_jacobian!(J, z, constraint)),
    ]
end

function executable_blocks(constraint::SpatialCVPhaseConstraint)
    rows = [constraint.acceleration_equation,
            constraint.velocity_equation,
            constraint.position_equation]
    ExecutableEquationBlock[
        ExecutableEquationBlock(constraint.name, :constraint, rows,
            (e, t, z, zd) -> cv_phase_constraints!(e, z, constraint),
            (J, t, z, zd, c) ->
                cv_phase_constraints_jacobian!(J, z, constraint)),
    ]
end

function executable_blocks(constraint::SpatialInplaneConstraint)
    rows = [constraint.acceleration_equation,
            constraint.velocity_equation,
            constraint.position_equation]
    ExecutableEquationBlock[
        ExecutableEquationBlock(constraint.name, :constraint, rows,
            (e, t, z, zd) -> inplane_constraints!(e, z, constraint),
            (J, t, z, zd, c) ->
                inplane_constraints_jacobian!(J, z, constraint)),
    ]
end

function executable_blocks(constraint::SpatialInlineConstraint)
    blocks = [executable_blocks(constraint.inplane_x);
              executable_blocks(constraint.inplane_y)]
    isempty(constraint.translation_variables) && return blocks
    acceleration, velocity, distance = constraint.translation_variables
    acceleration_equation, velocity_equation, distance_equation =
        constraint.translation_equations
    geometry = constraint.axial_geometry

    coordinates! = function (equations, t, z, zdot)
        values = directed_distance_values(geometry, z)
        equations[acceleration_equation] =
            z[acceleration] - values.acceleration
        equations[velocity_equation] =
            z[velocity] - values.velocity
        equations[distance_equation] =
            z[distance] - values.position
    end
    coordinates_jacobian! = function (jacobian, t, z, zdot, coefficient)
        jacobian[acceleration_equation, acceleration] += 1
        jacobian[velocity_equation, velocity] += 1
        jacobian[distance_equation, distance] += 1
        directed_distance_jacobian!(jacobian, z, geometry,
            acceleration_equation, velocity_equation, distance_equation, -1)
    end
    state! = function (equations, t, z, zdot)
        equations[constraint.state_equations[1]] =
            z[acceleration] - zdot[velocity]
        equations[constraint.state_equations[2]] =
            z[velocity] - zdot[distance]
    end
    state_jacobian! = function (jacobian, t, z, zdot, coefficient)
        jacobian[constraint.state_equations[1], acceleration] += 1
        jacobian[constraint.state_equations[1], velocity] -= coefficient
        jacobian[constraint.state_equations[2], velocity] += 1
        jacobian[constraint.state_equations[2], distance] -= coefficient
    end
    append!(blocks, ExecutableEquationBlock[
        ExecutableEquationBlock(constraint.name, :translation_coordinates,
            copy(constraint.translation_equations), coordinates!,
            coordinates_jacobian!),
        ExecutableEquationBlock(constraint.name, :selected_translation_state,
            copy(constraint.state_equations), state!, state_jacobian!),
    ])
    blocks
end

function executable_blocks(constraint::SpatialHingeConstraint)
    blocks = [executable_blocks(constraint.perp_xz);
              executable_blocks(constraint.perp_yz)]
    isempty(constraint.rotation_variables) && return blocks
    alpha, omega, theta = constraint.rotation_variables
    alpha_equation, omega_equation, theta_equation =
        constraint.rotation_equations

    angle! = function (equations, t, z, zdot)
        geometry = hinge_angle_geometry(constraint, z)
        sine_theta, cosine_theta = sincos(z[theta])
        numerator = sine_theta * geometry.cosine -
            cosine_theta * geometry.sine
        denominator = cosine_theta * geometry.cosine +
            sine_theta * geometry.sine
        equations[theta_equation] = atan(numerator, denominator)
    end
    angle_jacobian! = function (jacobian, t, z, zdot, coefficient)
        geometry = hinge_angle_geometry(constraint, z)
        first_x, second_x, second_z = geometry.first_x,
            geometry.second_x, geometry.second_z
        sine_theta, cosine_theta = sincos(z[theta])
        numerator = sine_theta * geometry.cosine -
            cosine_theta * geometry.sine
        denominator = cosine_theta * geometry.cosine +
            sine_theta * geometry.sine
        magnitude_squared = numerator^2 + denominator^2
        magnitude_squared > eps(Float64) || throw(DomainError(
            magnitude_squared, "hinge angle is undefined"))
        cosine_factor = (denominator * sine_theta -
            numerator * cosine_theta) / magnitude_squared
        sine_factor = (-denominator * cosine_theta -
            numerator * sine_theta) / magnitude_squared
        jacobian[theta_equation, theta] += 1
        if !isnothing(first_x.body)
            gradient = cosine_factor .* (transpose(
                first_x.direction_parameters) * second_x.direction) .+
                sine_factor .* (transpose(first_x.direction_parameters) *
                    cross(second_z.direction, second_x.direction))
            jacobian[theta_equation,
                first_x.body.euler_parameter_variables] .+= gradient
        end
        if !isnothing(second_x.body)
            gradient = cosine_factor .* (transpose(
                second_x.direction_parameters) * first_x.direction) .+
                sine_factor .* (
                    transpose(second_x.direction_parameters) *
                        cross(first_x.direction, second_z.direction) .+
                    transpose(second_z.direction_parameters) *
                        cross(second_x.direction, first_x.direction))
            jacobian[theta_equation,
                second_x.body.euler_parameter_variables] .+= gradient
        end
    end

    velocity! = function (equations, t, z, zdot)
        equations[omega_equation] = z[omega] -
            hinge_angular_velocity(constraint, z)
    end
    velocity_jacobian! = function (jacobian, t, z, zdot, coefficient)
        first = marker_angular_kinematics(constraint.marker_i, z)
        second = marker_angular_kinematics(constraint.marker_j, z)
        axis = marker_axis_kinematics(constraint.marker_j, z, 3)
        difference = first.omega - second.omega
        jacobian[omega_equation, omega] += 1
        if !isnothing(first.body)
            jacobian[omega_equation,
                first.body.euler_parameter_variables] .-=
                transpose(first.omega_parameters) * axis.direction
            jacobian[omega_equation,
                first.body.angular_velocity_variables] .-=
                transpose(first.orientation) * axis.direction
        end
        if !isnothing(second.body)
            jacobian[omega_equation,
                second.body.euler_parameter_variables] .-=
                -transpose(second.omega_parameters) * axis.direction .+
                transpose(axis.direction_parameters) * difference
            jacobian[omega_equation,
                second.body.angular_velocity_variables] .+=
                transpose(second.orientation) * axis.direction
        end
    end

    acceleration! = function (equations, t, z, zdot)
        equations[alpha_equation] = z[alpha] -
            hinge_angular_acceleration(constraint, z)
    end
    acceleration_jacobian! = function (jacobian, t, z, zdot, coefficient)
        first = marker_angular_kinematics(constraint.marker_i, z)
        second = marker_angular_kinematics(constraint.marker_j, z)
        axis = marker_axis_kinematics(constraint.marker_j, z, 3)
        omega_difference = first.omega - second.omega
        alpha_difference = first.alpha - second.alpha
        jacobian[alpha_equation, alpha] += 1
        if !isnothing(first.body)
            jacobian[alpha_equation,
                first.body.euler_parameter_variables] .-=
                transpose(first.alpha_parameters) * axis.direction .+
                transpose(first.omega_parameters) * axis.velocity
            jacobian[alpha_equation,
                first.body.angular_acceleration_variables] .-=
                transpose(first.orientation) * axis.direction
            jacobian[alpha_equation,
                first.body.angular_velocity_variables] .-=
                transpose(first.orientation) * axis.velocity
        end
        if !isnothing(second.body)
            parameter_gradient =
                -transpose(second.alpha_parameters) * axis.direction .+
                transpose(axis.direction_parameters) * alpha_difference .-
                transpose(second.omega_parameters) * axis.velocity .+
                transpose(axis.velocity_parameters) * omega_difference
            omega_gradient =
                -transpose(second.orientation) * axis.velocity .+
                transpose(axis.velocity_omega) * omega_difference
            jacobian[alpha_equation,
                second.body.euler_parameter_variables] .-=
                parameter_gradient
            jacobian[alpha_equation,
                second.body.angular_acceleration_variables] .+=
                transpose(second.orientation) * axis.direction
            jacobian[alpha_equation,
                second.body.angular_velocity_variables] .-=
                omega_gradient
        end
    end
    state! = function (equations, t, z, zdot)
        equations[constraint.state_equations[1]] =
            z[alpha] - zdot[omega]
        equations[constraint.state_equations[2]] =
            z[omega] - zdot[theta]
    end
    state_jacobian! = function (jacobian, t, z, zdot, coefficient)
        jacobian[constraint.state_equations[1], alpha] += 1
        jacobian[constraint.state_equations[1], omega] -= coefficient
        jacobian[constraint.state_equations[2], omega] += 1
        jacobian[constraint.state_equations[2], theta] -= coefficient
    end
    append!(blocks, ExecutableEquationBlock[
        ExecutableEquationBlock(constraint.name, :rotation_acceleration,
            [alpha_equation], acceleration!, acceleration_jacobian!),
        ExecutableEquationBlock(constraint.name, :rotation_velocity,
            [omega_equation], velocity!, velocity_jacobian!),
        ExecutableEquationBlock(constraint.name, :rotation_position,
            [theta_equation], angle!, angle_jacobian!),
        ExecutableEquationBlock(constraint.name, :selected_rotation_state,
            copy(constraint.state_equations), state!, state_jacobian!),
    ])
    blocks
end


function executable_blocks(joint::SpatialRevoluteJoint)
    [executable_blocks(joint.spherical);
     executable_blocks(joint.hinge)]
end

function executable_blocks(constraint::SpatialOrientConstraint)
    [executable_blocks(constraint.hinge);
     executable_blocks(constraint.perp_xy)]
end

function executable_blocks(joint::SpatialFixedJoint)
    [executable_blocks(joint.spherical);
     executable_blocks(joint.orient)]
end

function executable_blocks(joint::SpatialConstantVelocityJoint)
    [executable_blocks(joint.spherical);
     executable_blocks(joint.phase)]
end


function executable_blocks(joint::SpatialCylindricalJoint)
    [executable_blocks(joint.inline);
     executable_blocks(joint.hinge)]
end

function executable_blocks(joint::SpatialTranslationalJoint)
    [executable_blocks(joint.inline);
     executable_blocks(joint.orient)]
end

function body_force_rows(marker)
    marker isa Union{SpatialBodyMarker,SpatialFlexibleBeamMarker} || return Int[]
    collect(marker.body.balance_equations)
end

function body_torque_rows(marker)
    marker isa Union{SpatialBodyMarker,SpatialFlexibleBeamMarker} || return Int[]
    rows = collect(marker.body.balance_equations[4:6])
    marker isa SpatialFlexibleBeamMarker &&
        append!(rows, marker.body.balance_equations[7:12])
    rows
end

function add_reaction_to_body!(equations, z, marker, reaction, sign)
    marker isa SpatialGroundMarker && return nothing
    body = marker.body
    force = sign .* reaction
    parameters = @view z[body.euler_parameter_variables]
    force_body = transpose(rotation_matrix(parameters)) * force
    offset = marker.position_body
    if marker isa SpatialFlexibleBeamMarker
        offset = offset + marker.translation_shape *
            (@view z[body.elastic_position_variables])
    end
    torque_body = cross(offset, force_body)
    equations[body.balance_equations[1:3]] .-= force
    equations[body.balance_equations[4:6]] .-= torque_body
    if marker isa SpatialFlexibleBeamMarker
        equations[body.balance_equations[7:12]] .-=
            transpose(marker.translation_shape) * force_body
    end
    nothing
end

function add_reaction_jacobian!(jacobian, z, marker, reaction, sign,
        reaction_variables)
    marker isa SpatialGroundMarker && return nothing
    if marker isa SpatialFlexibleBeamMarker
        rows = collect(marker.body.balance_equations)
        columns = [collect(reaction_variables);
                   spatial_marker_dependency_indices(marker)]
        values = function (local_z)
            local_equations = zeros(eltype(local_z), maximum(rows))
            add_reaction_to_body!(local_equations, local_z, marker,
                @view(local_z[reaction_variables]), sign)
            local_equations[rows]
        end
        add_ad_jacobian!(jacobian, z, rows, columns, values)
        return nothing
    end
    body = marker.body
    parameters = @view z[body.euler_parameter_variables]
    orientation = rotation_matrix(parameters)
    offset_cross = skew(marker.position_body)
    force = sign .* reaction
    force_rows = body.balance_equations[1:3]
    torque_rows = body.balance_equations[4:6]
    jacobian[force_rows, reaction_variables] .-=
        sign .* Matrix{eltype(z)}(I, 3, 3)
    jacobian[torque_rows, reaction_variables] .-=
        sign .* offset_cross * transpose(orientation)
    jacobian[torque_rows, body.euler_parameter_variables] .-=
        offset_cross * rotation_transpose_vector_jacobian(parameters, force)
    nothing
end

function equation_contributions(joint::SpatialSphericalJoint)
    rows = unique!([body_force_rows(joint.marker_a);
                    body_force_rows(joint.marker_b)])
    residual! = function (equations, t, z, zdot)
        reaction = @view z[joint.reaction_variables]
        add_reaction_to_body!(equations, z, joint.marker_a, reaction, 1)
        add_reaction_to_body!(equations, z, joint.marker_b, reaction, -1)
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        reaction = @view z[joint.reaction_variables]
        add_reaction_jacobian!(jacobian, z, joint.marker_a, reaction, 1,
            joint.reaction_variables)
        add_reaction_jacobian!(jacobian, z, joint.marker_b, reaction, -1,
            joint.reaction_variables)
    end
    EquationContribution[EquationContribution(joint.name, :reaction,
        rows, residual!, jacobian!)]
end

function add_perp_reaction_to_body!(equations, z, marker, torque, sign)
    marker isa SpatialGroundMarker && return nothing
    body = marker.body
    parameters = @view z[body.euler_parameter_variables]
    torque_body = transpose(rotation_matrix(parameters)) * (sign .* torque)
    equations[body.balance_equations[4:6]] .-= torque_body
    if marker isa SpatialFlexibleBeamMarker
        equations[body.balance_equations[7:12]] .-=
            transpose(marker.orientation_shape) * torque_body
    end
    nothing
end

function add_perp_reaction_jacobian!(jacobian, z, constraint, marker, sign,
        normal, normal_derivatives)
    marker isa SpatialGroundMarker && return nothing
    if marker isa SpatialFlexibleBeamMarker
        rows = body_torque_rows(marker)
        columns = [constraint.reaction_variable;
                   constraint_dependencies(
                       constraint.marker_i, constraint.marker_j)]
        values = function (local_z)
            local_equations = zeros(eltype(local_z), maximum(rows))
            local_torque = local_z[constraint.reaction_variable] .*
                perp_normal(constraint, local_z)
            add_perp_reaction_to_body!(local_equations, local_z, marker,
                local_torque, sign)
            local_equations[rows]
        end
        add_ad_jacobian!(jacobian, z, rows, columns, values)
        return nothing
    end
    body = marker.body
    parameters = @view z[body.euler_parameter_variables]
    orientation = rotation_matrix(parameters)
    reaction = z[constraint.reaction_variable]
    rows = body.balance_equations[4:6]
    jacobian[rows, constraint.reaction_variable] .-=
        sign .* transpose(orientation) * normal
    jacobian[rows, body.euler_parameter_variables] .-=
        rotation_transpose_vector_jacobian(
            parameters, sign .* reaction .* normal)
    for (axis_body, derivative) in normal_derivatives
        jacobian[rows, axis_body.euler_parameter_variables] .-=
            sign .* reaction .* transpose(orientation) * derivative
    end
    nothing
end

function equation_contributions(constraint::SpatialPerpConstraint)
    rows = unique!([body_torque_rows(constraint.marker_i);
                    body_torque_rows(constraint.marker_j)])
    residual! = function (equations, t, z, zdot)
        torque = z[constraint.reaction_variable] .* perp_normal(constraint, z)
        add_perp_reaction_to_body!(
            equations, z, constraint.marker_i, torque, 1)
        add_perp_reaction_to_body!(
            equations, z, constraint.marker_j, torque, -1)
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        if is_flexible_marker(constraint.marker_i) ||
                is_flexible_marker(constraint.marker_j)
            columns = [constraint.reaction_variable;
                       constraint_dependencies(
                           constraint.marker_i, constraint.marker_j)]
            values = function (local_z)
                local_equations = zeros(eltype(local_z), maximum(rows))
                torque = local_z[constraint.reaction_variable] .*
                    perp_normal(constraint, local_z)
                add_perp_reaction_to_body!(local_equations, local_z,
                    constraint.marker_i, torque, 1)
                add_perp_reaction_to_body!(local_equations, local_z,
                    constraint.marker_j, torque, -1)
                local_equations[rows]
            end
            add_ad_jacobian!(jacobian, z, rows, columns, values)
            return nothing
        end
        kinematics = perp_kinematics(constraint, z)
        first, second = kinematics.first, kinematics.second
        normal = cross(first.direction, second.direction)
        normal_derivatives = Tuple{Any,Any}[]
        if !isnothing(first.body)
            push!(normal_derivatives, (first.body,
                -skew(second.direction) * first.direction_parameters))
        end
        if !isnothing(second.body)
            push!(normal_derivatives, (second.body,
                skew(first.direction) * second.direction_parameters))
        end
        add_perp_reaction_jacobian!(jacobian, z, constraint,
            constraint.marker_i, 1, normal, normal_derivatives)
        add_perp_reaction_jacobian!(jacobian, z, constraint,
            constraint.marker_j, -1, normal, normal_derivatives)
    end
    EquationContribution[EquationContribution(constraint.name, :reaction,
        rows, residual!, jacobian!)]
end

function equation_contributions(constraint::SpatialCVPhaseConstraint)
    rows = unique!([body_torque_rows(constraint.marker_i);
                    body_torque_rows(constraint.marker_j)])
    residual! = function (equations, t, z, zdot)
        directions = cv_phase_reaction_directions(constraint, z)
        reaction = z[constraint.reaction_variable]
        add_perp_reaction_to_body!(equations, z, constraint.marker_i,
            reaction .* directions.first, 1)
        add_perp_reaction_to_body!(equations, z, constraint.marker_j,
            reaction .* directions.second, 1)
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        columns = [constraint.reaction_variable;
                   constraint_dependencies(
                       constraint.marker_i, constraint.marker_j)]
        values = function (local_z)
            local_equations = zeros(eltype(local_z), maximum(rows))
            directions = cv_phase_reaction_directions(constraint, local_z)
            reaction = local_z[constraint.reaction_variable]
            add_perp_reaction_to_body!(local_equations, local_z,
                constraint.marker_i, reaction .* directions.first, 1)
            add_perp_reaction_to_body!(local_equations, local_z,
                constraint.marker_j, reaction .* directions.second, 1)
            local_equations[rows]
        end
        add_ad_jacobian!(jacobian, z, rows, columns, values)
    end
    EquationContribution[EquationContribution(constraint.name, :reaction,
        rows, residual!, jacobian!)]
end

function equation_contributions(constraint::SpatialOrientConstraint)
    [equation_contributions(constraint.hinge);
     equation_contributions(constraint.perp_xy)]
end

function equation_contributions(joint::SpatialFixedJoint)
    [equation_contributions(joint.spherical);
     equation_contributions(joint.orient)]
end

function equation_contributions(joint::SpatialConstantVelocityJoint)
    [equation_contributions(joint.spherical);
     equation_contributions(joint.phase)]
end


function equation_contributions(joint::SpatialCylindricalJoint)
    [equation_contributions(joint.inline);
     equation_contributions(joint.hinge)]
end

function equation_contributions(joint::SpatialTranslationalJoint)
    [equation_contributions(joint.inline);
     equation_contributions(joint.orient)]
end

directed_distance_reaction_rows(geometry::SpatialDirectedDistance) =
    unique!([body_force_rows(geometry.marker_i);
             body_force_rows(geometry.marker_j)])

function add_directed_distance_reaction!(equations, z, geometry,
        reaction_variable, scale = 1)
    values = directed_distance_values(geometry, z)
    reaction = scale * z[reaction_variable]
    force = reaction .* values.direction.direction
    add_reaction_to_body!(equations, z, geometry.marker_i, force, 1)

    marker_j = geometry.marker_j
    marker_j isa SpatialGroundMarker && return nothing
    body = marker_j.body
    parameters = @view z[body.euler_parameter_variables]
    orientation = rotation_matrix(parameters)
    force_body = reaction .*
        marker_j.orientation_body[:, geometry.axis.index]
    lever_body = transpose(orientation) *
        (values.first.position - z[body.position_variables])
    equations[body.balance_equations[1:3]] .+= force
    equations[body.balance_equations[4:6]] .+=
        cross(lever_body, force_body)
    nothing
end

function add_directed_distance_reaction_jacobian!(jacobian, z, geometry,
        reaction_variable, scale = 1)
    values = directed_distance_values(geometry, z)
    marker_i, marker_j = geometry.marker_i, geometry.marker_j
    normal = values.direction.direction
    reaction = scale * z[reaction_variable]
    force = reaction .* normal

    if marker_i isa SpatialBodyMarker
        body = marker_i.body
        parameters = @view z[body.euler_parameter_variables]
        orientation = rotation_matrix(parameters)
        force_rows = body.balance_equations[1:3]
        torque_rows = body.balance_equations[4:6]
        offset_cross = skew(marker_i.position_body)
        jacobian[force_rows, reaction_variable] .-= scale .* normal
        jacobian[torque_rows, reaction_variable] .-=
            scale .* offset_cross * transpose(orientation) * normal
        jacobian[torque_rows, body.euler_parameter_variables] .-=
            offset_cross *
                rotation_transpose_vector_jacobian(parameters, force)
        if !isnothing(values.direction.body)
            normal_columns = values.direction.body.euler_parameter_variables
            normal_derivative = values.direction.direction_parameters
            jacobian[force_rows, normal_columns] .-=
                reaction .* normal_derivative
            jacobian[torque_rows, normal_columns] .-=
                reaction .* offset_cross * transpose(orientation) *
                    normal_derivative
        end
    end

    if marker_j isa SpatialBodyMarker
        body = marker_j.body
        parameters = @view z[body.euler_parameter_variables]
        orientation = rotation_matrix(parameters)
        position = @view z[body.position_variables]
        lever_global = values.first.position - position
        lever_body = transpose(orientation) * lever_global
        normal_body = marker_j.orientation_body[:, geometry.axis.index]
        force_rows = body.balance_equations[1:3]
        torque_rows = body.balance_equations[4:6]
        jacobian[force_rows, reaction_variable] .+= scale .* normal
        jacobian[force_rows, body.euler_parameter_variables] .+=
            reaction .* values.direction.direction_parameters
        jacobian[torque_rows, reaction_variable] .+=
            scale .* cross(lever_body, normal_body)

        lever_gradient = -reaction .* skew(normal_body)
        jacobian[torque_rows, body.position_variables] .+=
            lever_gradient * (-transpose(orientation))
        jacobian[torque_rows, body.euler_parameter_variables] .+=
            lever_gradient * rotation_transpose_vector_jacobian(
                parameters, lever_global)
        if marker_i isa SpatialBodyMarker
            first_body = marker_i.body
            jacobian[torque_rows, first_body.position_variables] .+=
                lever_gradient * transpose(orientation)
            jacobian[torque_rows,
                first_body.euler_parameter_variables] .+=
                lever_gradient * transpose(orientation) *
                    values.first.position_parameters
        end
    end
    nothing
end

add_inplane_reaction!(equations, z, constraint) =
    add_directed_distance_reaction!(equations, z, constraint.geometry,
        constraint.reaction_variable)

add_inplane_reaction_jacobian!(jacobian, z, constraint) =
    add_directed_distance_reaction_jacobian!(jacobian, z,
        constraint.geometry, constraint.reaction_variable)

function equation_contributions(constraint::SpatialInplaneConstraint)
    rows = directed_distance_reaction_rows(constraint.geometry)
    residual! = (equations, t, z, zdot) ->
        add_inplane_reaction!(equations, z, constraint)
    jacobian! = (jacobian, t, z, zdot, coefficient) ->
        add_inplane_reaction_jacobian!(jacobian, z, constraint)
    EquationContribution[EquationContribution(constraint.name, :reaction,
        rows, residual!, jacobian!)]
end

function equation_contributions(constraint::SpatialInlineConstraint)
    [equation_contributions(constraint.inplane_x);
     equation_contributions(constraint.inplane_y)]
end

function equation_contributions(constraint::SpatialHingeConstraint)
    [equation_contributions(constraint.perp_xz);
     equation_contributions(constraint.perp_yz)]
end


function equation_contributions(joint::SpatialRevoluteJoint)
    [equation_contributions(joint.spherical);
     equation_contributions(joint.hinge)]
end

end
