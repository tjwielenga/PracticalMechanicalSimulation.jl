"""
    PlanarComponentAssembly

Planar component types, canonical declarations, and executable equation and
Jacobian callbacks. A component keeps only its physical data and allocated
canonical indices. It owns equations for its local definitions or constraints
and contributes forces and reactions to the balance equations owned by bodies.

New element implementations normally provide a registration method, an
allocated component type, owned executable blocks, and any additive equation
contributions. See `architecture/planar/adding-planar-component.md` for the full contract.
"""
module PlanarComponentAssembly

using LinearAlgebra
using ..AutomaticAnalysis
using ..ScalarExpressions: ScalarLaw
using ..PlanarAppliedForces
using ..PlanarDirectedDistances

export PlanarRigidBodyComponent, PlanarGravityComponent,
       PlanarScalarLaw, PlanarAppliedForceComponent,
       PlanarConstantTorqueComponent, PlanarAppliedTorqueComponent,
       PlanarRevoluteJointComponent, PlanarGroundRevoluteJointComponent,
       PlanarSpanningForceComponent, PlanarSpanMeasureComponent,
       PlanarTorsionalSpringComponent,
       PlanarBushingComponent,
       PlanarPlaneContactComponent,
       PlanarRotationalMotionGenerator, PlanarTranslationalMotionGenerator,
       PlanarDistanceCoordinateComponent, PlanarCoordinateCoupler,
       PlanarInplaneConstraint, PlanarPerpConstraint,
       PlanarTranslationalJoint, PlanarFixedJoint,
       PlanarGearPairComponent, PlanarRackAndPinionComponent,
       PlanarPulleyComponent, PlanarBeltComponent, PlanarBeltSpanComponent,
       planar_body_registration, revolute_joint_registration,
       applied_force_registration, applied_torque_registration,
       spanning_force_registration, span_measure_registration,
       torsional_spring_registration,
       bushing_registration,
       plane_contact_registration, perp_constraint_registration,
       gear_pair_registration, rack_and_pinion_registration,
       belt_span_registration,
       rotational_motion_registration, translational_motion_registration,
       distance_coordinate_registration, coordinate_coupler_registration,
       inplane_constraint_registration,
       component_registration,
       executable_blocks, equation_contributions, applied_force_value,
       applied_torque_value, plane_contact_gap,
       plane_contact_damping_surface, belt_tangent_geometry,
       belt_span_values,
       set_planar_applied_force_stage!, set_planar_applied_torque_stage!,
       set_planar_spanning_force_stage!, set_planar_bushing_stage!,
       set_planar_plane_contact_stage!

"""
Allocated planar rigid body in the unreduced canonical system.

The body always retains acceleration, velocity, and configuration variables.
`state_equations` selects the present independent coordinates, while the
candidate dictionaries retain dormant alternatives for runtime reselection.
A zero mass or inertia is allowed; the assembled model must still obtain a
nonsingular square system through constraints or force balance.
"""
struct PlanarRigidBodyComponent{T}
    name::Symbol
    mass::T
    inertia::T
    acceleration_variables::UnitRange{Int}
    angular_acceleration_variable::Int
    velocity_variables::UnitRange{Int}
    angular_velocity_variable::Int
    position_variables::UnitRange{Int}
    orientation_variable::Int
    balance_equations::UnitRange{Int}
    state_equations::Vector{UnitRange{Int}}
    selected_state_variables::Vector{NTuple{3,Int}}
    candidate_state_equations::Dict{Symbol,UnitRange{Int}}
    candidate_state_variables::Dict{Symbol,NTuple{3,Int}}
end

function PlanarRigidBodyComponent(name, mass, inertia, acceleration_variables,
        angular_acceleration_variable, velocity_variables,
        angular_velocity_variable, position_variables, orientation_variable,
        balance_equations, state_equations::UnitRange{Int})
    selected = isempty(state_equations) ? NTuple{3,Int}[] :
        [(angular_acceleration_variable, angular_velocity_variable,
          orientation_variable)]
    PlanarRigidBodyComponent(name, mass, inertia, acceleration_variables,
        angular_acceleration_variable, velocity_variables,
        angular_velocity_variable, position_variables, orientation_variable,
        balance_equations, isempty(state_equations) ? UnitRange{Int}[] :
        [state_equations], selected, Dict{Symbol,UnitRange{Int}}(),
        Dict{Symbol,NTuple{3,Int}}())
end

"""Constant global gravitational acceleration applied to one planar body."""
struct PlanarGravityComponent{B,T}
    name::Symbol
    body::B
    acceleration::Vector{T}
end

const PlanarScalarLaw = ScalarLaw

"""
Force applied at one point marker along the local y axis of an orientation
marker. An optional floating marker carries the equal-and-opposite force on a
second body.
"""
struct PlanarAppliedForceComponent{P,A,R,F,S}
    name::Symbol
    application_marker::P
    direction_axis::A
    reaction_marker::R
    magnitude::F
    active_during::S
    active::Base.RefValue{Bool}
    magnitude_variable::Int
    magnitude_equation::Int
end


PlanarAppliedForceComponent(name, application_marker, direction_axis,
        reaction_marker, magnitude) =
    PlanarAppliedForceComponent(name, application_marker, direction_axis,
        reaction_marker, magnitude, (:static, :dynamic, :modal), Ref(true),
        0, 0)

function set_planar_applied_force_stage!(force::PlanarAppliedForceComponent,
        stage)
    force.active[] = stage in force.active_during
    force
end

"""Constant action and reaction torque between two orientation markers."""
struct PlanarConstantTorqueComponent{MA,MB,T}
    name::Symbol
    marker_a::MA
    marker_b::MB
    torque::T
end

"""Prescribed or state-dependent torque between two orientation markers."""
struct PlanarAppliedTorqueComponent{MA,MB,F,PA,PB,S}
    name::Symbol
    marker_a::MA
    marker_b::MB
    torque::F
    point_a::PA
    point_b::PB
    active_during::S
    active::Base.RefValue{Bool}
    torque_variable::Int
    torque_equation::Int
end

PlanarAppliedTorqueComponent(name, marker_a, marker_b, torque) =
    PlanarAppliedTorqueComponent(name, marker_a, marker_b, torque,
        nothing, nothing, (:static, :dynamic, :modal), Ref(true), 0, 0)

PlanarAppliedTorqueComponent(name, marker_a, marker_b, torque,
        point_a, point_b) =
    PlanarAppliedTorqueComponent(name, marker_a, marker_b, torque,
        point_a, point_b, (:static, :dynamic, :modal), Ref(true), 0, 0)

function set_planar_applied_torque_stage!(torque::PlanarAppliedTorqueComponent,
        stage)
    torque.active[] = stage in torque.active_during
    torque
end

"""
Two coincident-point constraints with an optional relative rotation coordinate.

When present, `rotation_variables` and `rotation_equations` explicitly define
relative angular acceleration, velocity, and angle. The coordinate becomes a
state candidate only when `candidate_state_equations` is retained.
"""
struct PlanarRevoluteJointComponent{BA,MA,BB,MB,OA,OB}
    name::Symbol
    body_a::BA
    marker_a::MA
    body_b::BB
    marker_b::MB
    reaction_variables::UnitRange{Int}
    acceleration_equations::UnitRange{Int}
    velocity_equations::UnitRange{Int}
    position_equations::UnitRange{Int}
    rotation_marker_a::OA
    rotation_marker_b::OB
    rotation_variables::Vector{Int}
    rotation_equations::Vector{Int}
    state_equations::Vector{Int}
    selected_state_variables::Vector{NTuple{3,Int}}
    candidate_state_equations::Vector{Int}
end

function PlanarRevoluteJointComponent(name, body_a, marker_a, body_b, marker_b,
        reaction_variables, acceleration_equations, velocity_equations,
        position_equations)
    PlanarRevoluteJointComponent(name, body_a, marker_a, body_b, marker_b,
        reaction_variables, acceleration_equations, velocity_equations,
        position_equations, nothing, nothing, Int[], Int[], Int[],
        NTuple{3,Int}[], Int[])
end

const PlanarGroundRevoluteJointComponent = PlanarRevoluteJointComponent

"""Axial marker-to-marker force with explicit geometry and a scalar law."""
struct PlanarSpanningForceComponent{E,L,S}
    name::Symbol
    element::E
    law::L
    active_during::S
    active::Base.RefValue{Bool}
end

function set_planar_spanning_force_stage!(force::PlanarSpanningForceComponent,
        stage)
    force.active[] = stage in force.active_during
    force
end

"""Reaction-free marker-to-marker span measurement."""
struct PlanarSpanMeasureComponent{B1,M1,B2,M2}
    name::Symbol
    body_1::B1
    marker_1::M1
    body_2::B2
    marker_2::M2
    spanning_variables::UnitRange{Int}
    length_variable::Int
    unit_variables::UnitRange{Int}
    length_rate_variable::Int
    length_acceleration_variable::Int
    spanning_equations::UnitRange{Int}
    length_equation::Int
    unit_equations::UnitRange{Int}
    length_rate_equation::Int
    length_acceleration_equation::Int
end

"""Marker-to-marker torsional spring-damper with one explicit torque."""
struct PlanarTorsionalSpringComponent{M1,M2,E,T}
    name::Symbol
    marker_1::M1
    marker_2::M2
    element::E
    damping_time_scale::T
end

"""Linear planar bushing expressed in the second marker's frame."""
struct PlanarBushingComponent{M1,M2,T,S}
    name::Symbol
    marker_1::M1
    marker_2::M2
    translational_stiffness::Vector{T}
    translational_damping::Vector{T}
    rotational_stiffness::T
    rotational_damping::T
    damping_time_scale::T
    free_position::Vector{T}
    free_angle::T
    active_during::S
    active::Base.RefValue{Bool}
    force_variables::UnitRange{Int}
    torque_variable::Int
    load_equations::UnitRange{Int}
end

function set_planar_bushing_stage!(bushing::PlanarBushingComponent, stage)
    bushing.active[] = stage in bushing.active_during
    bushing
end

"""
One-sided compliant contact between a marker-centered sphere and an oriented
plane.

The first marker locates the sphere center. The second marker's local y axis is
the outward plane normal. Positive gap is separation; negative gap is
penetration.
"""
struct PlanarPlaneContactComponent{M1,M2,G,T,S}
    name::Symbol
    marker_1::M1
    marker_2::M2
    geometry::G
    radius::T
    stiffness::T
    damping_factor::T
    active_during::S
    active::Base.RefValue{Bool}
    gap_variable::Int
    gap_rate_variable::Int
    normal_force_variable::Int
    global_force_variables::UnitRange{Int}
    contact_equations::UnitRange{Int}
end

function set_planar_plane_contact_stage!(contact::PlanarPlaneContactComponent,
        stage)
    contact.active[] = stage in contact.active_during
    contact
end

struct PlanarRotationalMotionGenerator{BA,MA,BB,MB,F,F1,F2}
    name::Symbol
    body_a::BA
    marker_a::MA
    body_b::BB
    marker_b::MB
    angle_variable::Int
    angular_velocity_variable::Int
    angular_acceleration_variable::Int
    torque_variable::Int
    position_equations::UnitRange{Int}
    velocity_equations::UnitRange{Int}
    acceleration_equations::UnitRange{Int}
    motion::F
    motion_derivative::F1
    motion_second_derivative::F2
end

"""Prescribed signed marker distance along the second marker's local y axis."""
struct PlanarTranslationalMotionGenerator{G,F,F1,F2}
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

"""Measured signed marker distance along the second marker's local y axis."""
struct PlanarDistanceCoordinateComponent{G}
    name::Symbol
    geometry::G
    distance_variable::Int
    velocity_variable::Int
    acceleration_variable::Int
    acceleration_equation::Int
    velocity_equation::Int
    position_equation::Int
    state_equations::Vector{Int}
    selected_state_variables::Vector{NTuple{3,Int}}
    candidate_state_equations::Vector{Int}
end

"""One ideal linear relation among relative rotation or distance coordinates."""
struct PlanarCoordinateCoupler{C,T}
    name::Symbol
    coordinates::Vector{C}
    coefficients::Vector{T}
    offset::T
    coordinate_kind::Symbol
    reaction_variable::Int
    acceleration_equation::Int
    velocity_equation::Int
    position_equation::Int
end

struct PlanarInplaneConstraint{G}
    name::Symbol
    geometry::G
    reaction_variable::Int
    acceleration_equation::Int
    velocity_equation::Int
    position_equation::Int
end

"""Primitive planar orientation constraint between two marker axes."""
struct PlanarPerpConstraint{BI,OI,BJ,OJ}
    name::Symbol
    body_i::BI
    orientation_i::OI
    body_j::BJ
    orientation_j::OJ
    reaction_variable::Int
    acceleration_equation::Int
    velocity_equation::Int
    position_equation::Int
end

"""Planar translational joint composed of inplane and perp primitives."""
struct PlanarTranslationalJoint{I,P}
    name::Symbol
    inplane::I
    perp::P
end

"""Planar fixed joint composed of revolute and perp primitives."""
struct PlanarFixedJoint{R,P}
    name::Symbol
    revolute::R
    perp::P
end

"""Ideal planar gear pair with signed pitch radii derived from its carrier."""
struct PlanarGearPairComponent{B1,M1,O1,B2,M2,O2,J1,J2,BC,MC,MR,OC,T}
    name::Symbol
    body_1::B1
    marker_1::M1
    orientation_1::O1
    body_2::B2
    marker_2::M2
    orientation_2::O2
    joint_1::J1
    joint_2::J2
    carrier_body::BC
    contact_marker::MC
    force_reference_marker::MR
    carrier_orientation::OC
    radius_1::T
    radius_2::T
    phase::T
    reaction_variable::Int
    acceleration_equation::Int
    velocity_equation::Int
    position_equation::Int
end

"""Ideal planar rack-and-pinion coupling located by two carrier joints."""
struct PlanarRackAndPinionComponent{BR,MR,BG,MG,OG,JT,JR,BC,MC,OC,T}
    name::Symbol
    rack_body::BR
    rack_marker::MR
    pinion_body::BG
    pinion_marker::MG
    pinion_orientation::OG
    translational_joint::JT
    revolute_joint::JR
    carrier_body::BC
    contact_marker::MC
    pinion_carrier_orientation::OC
    pitch_radius::T
    phase::T
    reaction_variable::Int
    acceleration_equation::Int
    velocity_equation::Int
    position_equation::Int
end

"""Circular pitch surface carried by one side of a revolute joint."""
struct PlanarPulleyComponent{B,J,M,T}
    name::Symbol
    body::B
    joint::J
    center_marker::M
    pitch_radius::T
end

"""One ordered closed belt assembled from local elastic tangent spans."""
struct PlanarBeltComponent{S,T}
    name::Symbol
    spans::Vector{S}
    initial_tension::T
    free_length::T
end

"""Elastic no-slip belt span tangent to two planar pulley pitch circles."""
struct PlanarBeltSpanComponent{P1,P2,T}
    name::Symbol
    belt_name::Symbol
    pulley_1::P1
    pulley_2::P2
    tangent_kind::Int
    tangent_side::Int
    beta_1::T
    beta_2::T
    reference_length::T
    reference_tangent::Vector{T}
    reference_angle_1::T
    reference_angle_2::T
    initial_extension::T
    stiffness::T
    damping::T
    damping_time_scale::T
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

function planar_body_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:a_x, :acceleration, 2),
        VariableDeclaration(:a_y, :acceleration, 2),
        VariableDeclaration(:alpha, :angular_acceleration, 2),
        VariableDeclaration(:V_x, :velocity, 1),
        VariableDeclaration(:V_y, :velocity, 1),
        VariableDeclaration(:omega, :angular_velocity, 1),
        VariableDeclaration(:R_x, :position, 0),
        VariableDeclaration(:R_y, :position, 0),
        VariableDeclaration(:theta, :orientation, 0),
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:balance, EquationDeclaration[
            EquationDeclaration(:sum_F_x, :balance, 2, :body_balance),
            EquationDeclaration(:sum_F_y, :balance, 2, :body_balance),
            EquationDeclaration(:sum_T, :balance, 2, :body_balance),
        ]),
        EquationBlockDeclaration(:selected_state, EquationDeclaration[
            EquationDeclaration(:alpha_state, :state_equation, 2,
                                :angular_state),
            EquationDeclaration(:omega_state, :state_equation, 1,
                                :angular_state),
        ]),
        EquationBlockDeclaration(:selected_V_x, EquationDeclaration[
            EquationDeclaration(:a_x_state, :state_equation, 2,
                                :translational_state),
            EquationDeclaration(:V_x_state, :state_equation, 1,
                                :translational_state),
        ]),
        EquationBlockDeclaration(:selected_V_y, EquationDeclaration[
            EquationDeclaration(:a_y_state, :state_equation, 2,
                                :translational_state),
            EquationDeclaration(:V_y_state, :state_equation, 1,
                                :translational_state),
        ]),
    ]
    return ComponentRegistration(name, variables, blocks)
end

component_registration(body::PlanarRigidBodyComponent) =
    planar_body_registration(body.name)

function applied_force_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:F, :applied_load, 2),
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:load, EquationDeclaration[
            EquationDeclaration(:magnitude, :applied_definition, 2,
                :constitutive_force),
        ]),
    ]
    ComponentRegistration(name, variables, blocks)
end

component_registration(force::PlanarAppliedForceComponent) =
    force.magnitude_variable == 0 ? nothing :
        applied_force_registration(force.name)

function applied_torque_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:T, :applied_load, 2),
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:load, EquationDeclaration[
            EquationDeclaration(:magnitude, :applied_definition, 2,
                :constitutive_torque),
        ]),
    ]
    ComponentRegistration(name, variables, blocks)
end

component_registration(torque::PlanarAppliedTorqueComponent) =
    torque.torque_variable == 0 ? nothing :
        applied_torque_registration(torque.name)

function revolute_joint_registration(name::Symbol;
        rotation_coordinates = false)
    variables = VariableDeclaration[
        VariableDeclaration(:lambda_x, :reaction, 2),
        VariableDeclaration(:lambda_y, :reaction, 2),
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:acceleration, EquationDeclaration[
            EquationDeclaration(:Phi_ddot_x, :constraint, 2, :pin),
            EquationDeclaration(:Phi_ddot_y, :constraint, 2, :pin),
        ]),
        EquationBlockDeclaration(:velocity, EquationDeclaration[
            EquationDeclaration(:Phi_dot_x, :constraint, 1, :pin),
            EquationDeclaration(:Phi_dot_y, :constraint, 1, :pin),
        ]),
        EquationBlockDeclaration(:position, EquationDeclaration[
            EquationDeclaration(:Phi_x, :constraint, 0, :pin),
            EquationDeclaration(:Phi_y, :constraint, 0, :pin),
        ]),
    ]
    families = ConstraintFamilyDeclaration[
        ConstraintFamilyDeclaration(:x, :lambda_x, :Phi_x,
            :Phi_dot_x, :Phi_ddot_x),
        ConstraintFamilyDeclaration(:y, :lambda_y, :Phi_y,
            :Phi_dot_y, :Phi_ddot_y),
    ]
    if rotation_coordinates
        append!(variables, VariableDeclaration[
            VariableDeclaration(:alpha, :relative_acceleration, 2),
            VariableDeclaration(:omega, :relative_velocity, 1),
            VariableDeclaration(:theta, :relative_position, 0),
        ])
        append!(blocks, EquationBlockDeclaration[
            EquationBlockDeclaration(:rotation_acceleration,
                EquationDeclaration[
                    EquationDeclaration(:alpha_definition,
                        :coordinate_relation, 2, :relative_rotation),
                ]),
            EquationBlockDeclaration(:rotation_velocity,
                EquationDeclaration[
                    EquationDeclaration(:omega_definition,
                        :coordinate_relation, 1, :relative_rotation),
                ]),
            EquationBlockDeclaration(:rotation_position,
                EquationDeclaration[
                    EquationDeclaration(:theta_definition,
                        :coordinate_relation, 0, :relative_rotation),
                ]),
            EquationBlockDeclaration(:selected_rotation_state,
                EquationDeclaration[
                    EquationDeclaration(:alpha_state, :state_equation, 2,
                        :relative_rotation_state),
                    EquationDeclaration(:omega_state, :state_equation, 1,
                        :relative_rotation_state),
                ]),
        ])
    end
    return ComponentRegistration(name, variables, blocks, families)
end

component_registration(joint::PlanarRevoluteJointComponent) =
    revolute_joint_registration(joint.name;
        rotation_coordinates = !isempty(joint.rotation_variables))

function spanning_force_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:s_x, :applied_geometry, 0),
        VariableDeclaration(:s_y, :applied_geometry, 0),
        VariableDeclaration(:length, :applied_geometry, 0),
        VariableDeclaration(:u_x, :applied_geometry, 0),
        VariableDeclaration(:u_y, :applied_geometry, 0),
        VariableDeclaration(:length_rate, :applied_rate, 1),
        VariableDeclaration(:force, :applied_load, 2),
        VariableDeclaration(:F_x, :applied_load, 2),
        VariableDeclaration(:F_y, :applied_load, 2),
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:geometry, EquationDeclaration[
            EquationDeclaration(:span_x, :applied_definition, 0, :span),
            EquationDeclaration(:span_y, :applied_definition, 0, :span),
            EquationDeclaration(:length, :applied_definition, 0, :length),
            EquationDeclaration(:unit_x, :applied_definition, 0, :unit),
            EquationDeclaration(:unit_y, :applied_definition, 0, :unit),
        ]),
        EquationBlockDeclaration(:rate, EquationDeclaration[
            EquationDeclaration(:length_rate, :applied_definition, 1,
                                :length_rate),
        ]),
        EquationBlockDeclaration(:load, EquationDeclaration[
            EquationDeclaration(:scalar_force, :applied_definition, 2,
                                :constitutive_force),
            EquationDeclaration(:global_force_x, :applied_definition, 2,
                                :global_force),
            EquationDeclaration(:global_force_y, :applied_definition, 2,
                                :global_force),
        ]),
    ]
    return ComponentRegistration(name, variables, blocks)
end

component_registration(force::PlanarSpanningForceComponent) =
    spanning_force_registration(force.name)

function span_measure_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:s_x, :applied_geometry, 0),
        VariableDeclaration(:s_y, :applied_geometry, 0),
        VariableDeclaration(:distance, :relative_position, 0),
        VariableDeclaration(:u_x, :applied_geometry, 0),
        VariableDeclaration(:u_y, :applied_geometry, 0),
        VariableDeclaration(:velocity, :relative_velocity, 1),
        VariableDeclaration(:acceleration, :relative_acceleration, 2),
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:geometry, EquationDeclaration[
            EquationDeclaration(:span_x, :applied_definition, 0, :span),
            EquationDeclaration(:span_y, :applied_definition, 0, :span),
            EquationDeclaration(:distance, :applied_definition, 0,
                :spanning_distance),
            EquationDeclaration(:unit_x, :applied_definition, 0, :unit),
            EquationDeclaration(:unit_y, :applied_definition, 0, :unit),
        ]),
        EquationBlockDeclaration(:velocity, EquationDeclaration[
            EquationDeclaration(:velocity, :applied_definition, 1,
                :spanning_distance),
        ]),
        EquationBlockDeclaration(:acceleration, EquationDeclaration[
            EquationDeclaration(:acceleration, :applied_definition, 2,
                :spanning_distance),
        ]),
    ]
    ComponentRegistration(name, variables, blocks)
end

component_registration(measure::PlanarSpanMeasureComponent) =
    span_measure_registration(measure.name)

function torsional_spring_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:T, :applied_load, 2),
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:load, EquationDeclaration[
            EquationDeclaration(:torque, :applied_definition, 2,
                :constitutive_torque),
        ]),
    ]
    ComponentRegistration(name, variables, blocks)
end

component_registration(spring::PlanarTorsionalSpringComponent) =
    torsional_spring_registration(spring.name)

function bushing_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:F_x, :applied_load, 2),
        VariableDeclaration(:F_y, :applied_load, 2),
        VariableDeclaration(:T, :applied_load, 2),
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:load, EquationDeclaration[
            EquationDeclaration(:force_x, :applied_definition, 2, :bushing_force),
            EquationDeclaration(:force_y, :applied_definition, 2, :bushing_force),
            EquationDeclaration(:torque, :applied_definition, 2, :bushing_torque),
        ]),
    ]
    ComponentRegistration(name, variables, blocks)
end

component_registration(bushing::PlanarBushingComponent) =
    bushing_registration(bushing.name)

function plane_contact_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:gap, :applied_geometry, 0),
        VariableDeclaration(:gap_rate, :applied_rate, 1),
        VariableDeclaration(:normal_force, :applied_load, 2),
        VariableDeclaration(:F_x, :applied_load, 2),
        VariableDeclaration(:F_y, :applied_load, 2),
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:contact, EquationDeclaration[
            EquationDeclaration(:gap, :applied_definition, 0, :contact_gap),
            EquationDeclaration(:gap_rate, :applied_definition, 1,
                                :contact_gap_rate),
            EquationDeclaration(:normal_force, :applied_definition, 2,
                                :contact_force),
            EquationDeclaration(:global_force_x, :applied_definition, 2,
                                :global_force),
            EquationDeclaration(:global_force_y, :applied_definition, 2,
                                :global_force),
        ]),
    ]
    ComponentRegistration(name, variables, blocks)
end

component_registration(contact::PlanarPlaneContactComponent) =
    plane_contact_registration(contact.name)

function rotational_motion_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:theta_m, :relative_position, 0),
        VariableDeclaration(:omega_m, :relative_velocity, 1),
        VariableDeclaration(:alpha_m, :relative_acceleration, 2),
        VariableDeclaration(:tau_m, :reaction, 2),
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:position, EquationDeclaration[
            EquationDeclaration(:relative_angle, :coordinate_relation, 0,
                                :relative_rotation),
            EquationDeclaration(:prescribed_angle, :motion, 0,
                                :prescribed_rotation),
        ]),
        EquationBlockDeclaration(:velocity, EquationDeclaration[
            EquationDeclaration(:relative_angular_velocity,
                                :coordinate_relation, 1,
                                :relative_rotation),
            EquationDeclaration(:prescribed_angular_velocity, :motion, 1,
                                :prescribed_rotation),
        ]),
        EquationBlockDeclaration(:acceleration, EquationDeclaration[
            EquationDeclaration(:relative_angular_acceleration,
                                :coordinate_relation, 2,
                                :relative_rotation),
            EquationDeclaration(:prescribed_angular_acceleration, :motion, 2,
                                :prescribed_rotation),
        ]),
    ]
    return ComponentRegistration(name, variables, blocks)
end

component_registration(generator::PlanarRotationalMotionGenerator) =
    rotational_motion_registration(generator.name)

function translational_motion_registration(name::Symbol)
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
                                :coordinate_relation, 2,
                                :relative_translation),
            EquationDeclaration(:prescribed_acceleration, :motion, 2,
                                :prescribed_translation),
        ]),
    ]
    ComponentRegistration(name, variables, blocks)
end

component_registration(generator::PlanarTranslationalMotionGenerator) =
    translational_motion_registration(generator.name)

function distance_coordinate_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:distance, :relative_position, 0),
        VariableDeclaration(:velocity, :relative_velocity, 1),
        VariableDeclaration(:acceleration, :relative_acceleration, 2),
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:position, EquationDeclaration[
            EquationDeclaration(:distance_definition, :coordinate_relation,
                                0, :relative_translation),
        ]),
        EquationBlockDeclaration(:velocity, EquationDeclaration[
            EquationDeclaration(:velocity_definition, :coordinate_relation,
                                1, :relative_translation),
        ]),
        EquationBlockDeclaration(:acceleration, EquationDeclaration[
            EquationDeclaration(:acceleration_definition,
                                :coordinate_relation, 2,
                                :relative_translation),
        ]),
        EquationBlockDeclaration(:selected_distance_state,
            EquationDeclaration[
                EquationDeclaration(:acceleration_state, :state_equation, 2,
                                    :relative_translation_state),
                EquationDeclaration(:velocity_state, :state_equation, 1,
                                    :relative_translation_state),
            ]),
    ]
    ComponentRegistration(name, variables, blocks)
end

component_registration(coordinate::PlanarDistanceCoordinateComponent) =
    distance_coordinate_registration(coordinate.name)

function coordinate_coupler_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:lambda, :reaction, 2),
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:acceleration, EquationDeclaration[
            EquationDeclaration(:Phi_ddot, :constraint, 2,
                                :coordinate_coupler),
        ]),
        EquationBlockDeclaration(:velocity, EquationDeclaration[
            EquationDeclaration(:Phi_dot, :constraint, 1,
                                :coordinate_coupler),
        ]),
        EquationBlockDeclaration(:position, EquationDeclaration[
            EquationDeclaration(:Phi, :constraint, 0,
                                :coordinate_coupler),
        ]),
    ]
    families = ConstraintFamilyDeclaration[
        ConstraintFamilyDeclaration(:coordinate_coupler, :lambda, :Phi,
            :Phi_dot, :Phi_ddot),
    ]
    ComponentRegistration(name, variables, blocks, families)
end

component_registration(coupler::PlanarCoordinateCoupler) =
    coordinate_coupler_registration(coupler.name)

function inplane_constraint_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:lambda, :reaction, 2),
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:acceleration, EquationDeclaration[
            EquationDeclaration(:Phi_ddot, :constraint, 2, :inplane),
        ]),
        EquationBlockDeclaration(:velocity, EquationDeclaration[
            EquationDeclaration(:Phi_dot, :constraint, 1, :inplane),
        ]),
        EquationBlockDeclaration(:position, EquationDeclaration[
            EquationDeclaration(:Phi, :constraint, 0, :inplane),
        ]),
    ]
    families = ConstraintFamilyDeclaration[
        ConstraintFamilyDeclaration(:normal, :lambda, :Phi,
            :Phi_dot, :Phi_ddot),
    ]
    return ComponentRegistration(name, variables, blocks, families)
end

component_registration(constraint::PlanarInplaneConstraint) =
    inplane_constraint_registration(constraint.name)

function perp_constraint_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:lambda, :reaction, 2),
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:acceleration, EquationDeclaration[
            EquationDeclaration(:Phi_ddot, :constraint, 2, :perp),
        ]),
        EquationBlockDeclaration(:velocity, EquationDeclaration[
            EquationDeclaration(:Phi_dot, :constraint, 1, :perp),
        ]),
        EquationBlockDeclaration(:position, EquationDeclaration[
            EquationDeclaration(:Phi, :constraint, 0, :perp),
        ]),
    ]
    families = ConstraintFamilyDeclaration[
        ConstraintFamilyDeclaration(:perp, :lambda, :Phi,
            :Phi_dot, :Phi_ddot),
    ]
    ComponentRegistration(name, variables, blocks, families)
end

component_registration(constraint::PlanarPerpConstraint) =
    perp_constraint_registration(constraint.name)

function gear_pair_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:lambda, :reaction, 2),
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:acceleration, EquationDeclaration[
            EquationDeclaration(:Phi_ddot, :constraint, 2, :gear_pair),
        ]),
        EquationBlockDeclaration(:velocity, EquationDeclaration[
            EquationDeclaration(:Phi_dot, :constraint, 1, :gear_pair),
        ]),
        EquationBlockDeclaration(:position, EquationDeclaration[
            EquationDeclaration(:Phi, :constraint, 0, :gear_pair),
        ]),
    ]
    families = ConstraintFamilyDeclaration[
        ConstraintFamilyDeclaration(:gear_pair, :lambda, :Phi,
            :Phi_dot, :Phi_ddot),
    ]
    ComponentRegistration(name, variables, blocks, families)
end

component_registration(gear::PlanarGearPairComponent) =
    gear_pair_registration(gear.name)

function rack_and_pinion_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:lambda, :reaction, 2),
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:acceleration, EquationDeclaration[
            EquationDeclaration(:Phi_ddot, :constraint, 2,
                                :rack_and_pinion),
        ]),
        EquationBlockDeclaration(:velocity, EquationDeclaration[
            EquationDeclaration(:Phi_dot, :constraint, 1,
                                :rack_and_pinion),
        ]),
        EquationBlockDeclaration(:position, EquationDeclaration[
            EquationDeclaration(:Phi, :constraint, 0, :rack_and_pinion),
        ]),
    ]
    families = ConstraintFamilyDeclaration[
        ConstraintFamilyDeclaration(:rack_and_pinion, :lambda, :Phi,
            :Phi_dot, :Phi_ddot),
    ]
    ComponentRegistration(name, variables, blocks, families)
end

component_registration(component::PlanarRackAndPinionComponent) =
    rack_and_pinion_registration(component.name)

function belt_span_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:P1_x, :applied_geometry, 0),
        VariableDeclaration(:P1_y, :applied_geometry, 0),
        VariableDeclaration(:P2_x, :applied_geometry, 0),
        VariableDeclaration(:P2_y, :applied_geometry, 0),
        VariableDeclaration(:t_x, :applied_geometry, 0),
        VariableDeclaration(:t_y, :applied_geometry, 0),
        VariableDeclaration(:ell, :applied_geometry, 0),
        VariableDeclaration(:extension, :applied_geometry, 0),
        VariableDeclaration(:extension_rate, :applied_rate, 1),
        VariableDeclaration(:tension, :applied_load, 2),
        VariableDeclaration(:F_x, :applied_load, 2),
        VariableDeclaration(:F_y, :applied_load, 2),
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:geometry, EquationDeclaration[
            EquationDeclaration(:point_1_x, :applied_definition, 0, :belt_tangent),
            EquationDeclaration(:point_1_y, :applied_definition, 0, :belt_tangent),
            EquationDeclaration(:point_2_x, :applied_definition, 0, :belt_tangent),
            EquationDeclaration(:point_2_y, :applied_definition, 0, :belt_tangent),
            EquationDeclaration(:tangent_x, :applied_definition, 0, :belt_tangent),
            EquationDeclaration(:tangent_y, :applied_definition, 0, :belt_tangent),
            EquationDeclaration(:length, :applied_definition, 0, :belt_length),
            EquationDeclaration(:extension, :applied_definition, 0, :belt_extension),
        ]),
        EquationBlockDeclaration(:rate, EquationDeclaration[
            EquationDeclaration(:extension_rate, :applied_definition, 1,
                :belt_extension_rate),
        ]),
        EquationBlockDeclaration(:load, EquationDeclaration[
            EquationDeclaration(:tension, :applied_definition, 2,
                :belt_tension),
            EquationDeclaration(:force_x, :applied_definition, 2,
                :global_force),
            EquationDeclaration(:force_y, :applied_definition, 2,
                :global_force),
        ]),
    ]
    ComponentRegistration(name, variables, blocks)
end

component_registration(span::PlanarBeltSpanComponent) =
    belt_span_registration(span.name)
component_registration(::PlanarPulleyComponent) = nothing
component_registration(::PlanarBeltComponent) = nothing

function body_balance!(equations, z, body)
    equations[body.balance_equations[1:2]] .=
        body.mass .* z[body.acceleration_variables]
    equations[body.balance_equations[3]] =
        body.inertia * z[body.angular_acceleration_variable]
end

function body_balance_jacobian!(jacobian, body)
    for k in 1:2
        jacobian[body.balance_equations[k], body.acceleration_variables[k]] +=
            body.mass
    end
    jacobian[body.balance_equations[3], body.angular_acceleration_variable] +=
        body.inertia
end

function body_states!(equations, z, zdot, body)
    pairs = isempty(body.candidate_state_equations) ?
        zip(body.state_equations, body.selected_state_variables) :
        ((body.candidate_state_equations[name], variables)
         for (name, variables) in body.candidate_state_variables)
    for (rows, (acceleration, velocity, position)) in pairs
        equations[rows[1]] = z[acceleration] - zdot[velocity]
        equations[rows[2]] = z[velocity] - zdot[position]
    end
end

function body_states_jacobian!(jacobian, coefficient, body)
    pairs = isempty(body.candidate_state_equations) ?
        zip(body.state_equations, body.selected_state_variables) :
        ((body.candidate_state_equations[name], variables)
         for (name, variables) in body.candidate_state_variables)
    for (rows, (acceleration, velocity, position)) in pairs
        jacobian[rows[1], acceleration] += 1
        jacobian[rows[1], velocity] -= coefficient
        jacobian[rows[2], velocity] += 1
        jacobian[rows[2], position] -= coefficient
    end
end

function executable_blocks(body::PlanarRigidBodyComponent)
    owned_state_equations = isempty(body.candidate_state_equations) ?
        collect(Iterators.flatten(body.state_equations)) :
        collect(Iterators.flatten(values(body.candidate_state_equations)))
    return ExecutableEquationBlock[
        ExecutableEquationBlock(body.name, :balance,
            collect(body.balance_equations),
            (e, t, z, zd) -> body_balance!(e, z, body),
            (J, t, z, zd, c) -> body_balance_jacobian!(J, body)),
        ExecutableEquationBlock(body.name, :selected_state,
            owned_state_equations,
            (e, t, z, zd) -> body_states!(e, z, zd, body),
            (J, t, z, zd, c) -> body_states_jacobian!(J, c, body)),
    ]
end

function joint_marker_acceleration(body::PlanarRigidBodyComponent,
                                   marker, values, z)
    return z[body.acceleration_variables] .+
        values.d .* z[body.angular_acceleration_variable] .-
        values.r .* z[body.angular_velocity_variable]^2
end

joint_marker_acceleration(::Nothing, marker, values, z) = zeros(eltype(z), 2)

function add_joint_marker_acceleration_jacobian!(jacobian, rows, sign,
        body::PlanarRigidBodyComponent, values, z)
    alpha = z[body.angular_acceleration_variable]
    omega = z[body.angular_velocity_variable]
    for k in 1:2
        jacobian[rows[k], body.acceleration_variables[k]] += sign
        jacobian[rows[k], body.angular_acceleration_variable] +=
            sign * values.d[k]
        jacobian[rows[k], body.angular_velocity_variable] +=
            sign * (-2 * values.r[k] * omega)
        jacobian[rows[k], body.orientation_variable] +=
            sign * (-values.r[k] * alpha - values.d[k] * omega^2)
    end
end

add_joint_marker_acceleration_jacobian!(jacobian, rows, sign, ::Nothing,
                                        values, z) = nothing

function add_joint_marker_velocity_jacobian!(jacobian, rows, sign,
        body::PlanarRigidBodyComponent, values, z)
    omega = z[body.angular_velocity_variable]
    for k in 1:2
        jacobian[rows[k], body.velocity_variables[k]] += sign
        jacobian[rows[k], body.angular_velocity_variable] +=
            sign * values.d[k]
        jacobian[rows[k], body.orientation_variable] +=
            sign * (-values.r[k] * omega)
    end
end

add_joint_marker_velocity_jacobian!(jacobian, rows, sign, ::Nothing,
                                    values, z) = nothing

function add_joint_marker_position_jacobian!(jacobian, rows, sign,
        body::PlanarRigidBodyComponent, values, z)
    for k in 1:2
        jacobian[rows[k], body.position_variables[k]] += sign
        jacobian[rows[k], body.orientation_variable] +=
            sign * values.d[k]
    end
end

add_joint_marker_position_jacobian!(jacobian, rows, sign, ::Nothing,
                                    values, z) = nothing

function executable_blocks(joint::PlanarRevoluteJointComponent)
    acceleration! = function (equations, t, z, zdot)
        values_a = PlanarAppliedForces.point_marker_kinematics(joint.marker_a, z)
        values_b = PlanarAppliedForces.point_marker_kinematics(joint.marker_b, z)
        equations[joint.acceleration_equations] .=
            joint_marker_acceleration(joint.body_a, joint.marker_a, values_a, z) .-
            joint_marker_acceleration(joint.body_b, joint.marker_b, values_b, z)
    end
    acceleration_jacobian! = function (jacobian, t, z, zdot, coefficient)
        values_a = PlanarAppliedForces.point_marker_kinematics(joint.marker_a, z)
        values_b = PlanarAppliedForces.point_marker_kinematics(joint.marker_b, z)
        add_joint_marker_acceleration_jacobian!(jacobian,
            joint.acceleration_equations, 1, joint.body_a, values_a, z)
        add_joint_marker_acceleration_jacobian!(jacobian,
            joint.acceleration_equations, -1, joint.body_b, values_b, z)
    end
    velocity! = function (equations, t, z, zdot)
        values_a = PlanarAppliedForces.point_marker_kinematics(joint.marker_a, z)
        values_b = PlanarAppliedForces.point_marker_kinematics(joint.marker_b, z)
        equations[joint.velocity_equations] .= values_a.velocity .- values_b.velocity
    end
    velocity_jacobian! = function (jacobian, t, z, zdot, coefficient)
        values_a = PlanarAppliedForces.point_marker_kinematics(joint.marker_a, z)
        values_b = PlanarAppliedForces.point_marker_kinematics(joint.marker_b, z)
        add_joint_marker_velocity_jacobian!(jacobian,
            joint.velocity_equations, 1, joint.body_a, values_a, z)
        add_joint_marker_velocity_jacobian!(jacobian,
            joint.velocity_equations, -1, joint.body_b, values_b, z)
    end
    position! = function (equations, t, z, zdot)
        values_a = PlanarAppliedForces.point_marker_kinematics(joint.marker_a, z)
        values_b = PlanarAppliedForces.point_marker_kinematics(joint.marker_b, z)
        equations[joint.position_equations] .= values_a.position .- values_b.position
    end
    position_jacobian! = function (jacobian, t, z, zdot, coefficient)
        values_a = PlanarAppliedForces.point_marker_kinematics(joint.marker_a, z)
        values_b = PlanarAppliedForces.point_marker_kinematics(joint.marker_b, z)
        add_joint_marker_position_jacobian!(jacobian,
            joint.position_equations, 1, joint.body_a, values_a, z)
        add_joint_marker_position_jacobian!(jacobian,
            joint.position_equations, -1, joint.body_b, values_b, z)
    end
    blocks = ExecutableEquationBlock[
        ExecutableEquationBlock(joint.name, :acceleration,
            collect(joint.acceleration_equations), acceleration!,
            acceleration_jacobian!),
        ExecutableEquationBlock(joint.name, :velocity,
            collect(joint.velocity_equations), velocity!, velocity_jacobian!),
        ExecutableEquationBlock(joint.name, :position,
            collect(joint.position_equations), position!, position_jacobian!),
    ]
    isempty(joint.rotation_variables) && return blocks

    alpha, omega, theta = joint.rotation_variables
    alpha_equation, omega_equation, theta_equation = joint.rotation_equations
    marker_alpha(body, z) = isnothing(body) ? zero(eltype(z)) :
        z[body.angular_acceleration_variable]
    rotation_coordinates! = function (equations, t, z, zdot)
        equations[alpha_equation] = z[alpha] -
            (marker_alpha(joint.body_a, z) - marker_alpha(joint.body_b, z))
        equations[omega_equation] = z[omega] -
            (PlanarAppliedForces.marker_angular_velocity(
                 joint.rotation_marker_a, z) -
             PlanarAppliedForces.marker_angular_velocity(
                 joint.rotation_marker_b, z))
        equations[theta_equation] = z[theta] -
            (PlanarAppliedForces.marker_angle(joint.rotation_marker_a, z) -
             PlanarAppliedForces.marker_angle(joint.rotation_marker_b, z))
    end
    rotation_coordinates_jacobian! =
            function (jacobian, t, z, zdot, coefficient)
        jacobian[alpha_equation, alpha] += 1
        jacobian[omega_equation, omega] += 1
        jacobian[theta_equation, theta] += 1
        if !isnothing(joint.body_a)
            jacobian[alpha_equation,
                joint.body_a.angular_acceleration_variable] -= 1
            jacobian[omega_equation,
                joint.body_a.angular_velocity_variable] -= 1
            jacobian[theta_equation,
                joint.body_a.orientation_variable] -= 1
        end
        if !isnothing(joint.body_b)
            jacobian[alpha_equation,
                joint.body_b.angular_acceleration_variable] += 1
            jacobian[omega_equation,
                joint.body_b.angular_velocity_variable] += 1
            jacobian[theta_equation,
                joint.body_b.orientation_variable] += 1
        end
    end
    push!(blocks, ExecutableEquationBlock(joint.name, :rotation_coordinates,
        copy(joint.rotation_equations), rotation_coordinates!,
        rotation_coordinates_jacobian!))

    state_equations = isempty(joint.candidate_state_equations) ?
        joint.state_equations : joint.candidate_state_equations
    if !isempty(state_equations)
        state! = function (equations, t, z, zdot)
            equations[state_equations[1]] = z[alpha] - zdot[omega]
            equations[state_equations[2]] = z[omega] - zdot[theta]
        end
        state_jacobian! = function (jacobian, t, z, zdot, coefficient)
            jacobian[state_equations[1], alpha] += 1
            jacobian[state_equations[1], omega] -= coefficient
            jacobian[state_equations[2], omega] += 1
            jacobian[state_equations[2], theta] -= coefficient
        end
        push!(blocks, ExecutableEquationBlock(joint.name,
            :selected_rotation_state, copy(state_equations), state!,
            state_jacobian!))
    end
    blocks
end

function equation_contributions(gravity::PlanarGravityComponent)
    body = gravity.body
    residual! = function (equations, t, z, zdot)
        equations[body.balance_equations[1:2]] .-=
            body.mass .* gravity.acceleration
    end
    return EquationContribution[
        EquationContribution(gravity.name, :gravity_to_body,
            collect(body.balance_equations), residual!,
            (J, t, z, zd, c) -> nothing),
    ]
end

function executable_blocks(force::PlanarAppliedForceComponent)
    force.magnitude_variable == 0 && return ExecutableEquationBlock[]
    residual! = function (equations, t, z, zdot)
        equations[force.magnitude_equation] =
            z[force.magnitude_variable] -
                (force.active[] ? force.magnitude(t, z) : 0.0)
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        jacobian[force.magnitude_equation, force.magnitude_variable] += 1
        if force.active[]
            partials = force.magnitude.gradient(t, z)
            for (column, partial) in zip(force.magnitude.dependencies, partials)
                jacobian[force.magnitude_equation, column] -= partial
            end
        end
    end
    ExecutableEquationBlock[ExecutableEquationBlock(force.name, :load,
        [force.magnitude_equation], residual!, jacobian!)]
end

function executable_blocks(torque::PlanarAppliedTorqueComponent)
    torque.torque_variable == 0 && return ExecutableEquationBlock[]
    residual! = function (equations, t, z, zdot)
        equations[torque.torque_equation] =
            z[torque.torque_variable] -
                (torque.active[] ? torque.torque(t, z) : 0.0)
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        jacobian[torque.torque_equation, torque.torque_variable] += 1
        if torque.active[]
            partials = torque.torque.gradient(t, z)
            for (column, partial) in zip(torque.torque.dependencies, partials)
                jacobian[torque.torque_equation, column] -= partial
            end
        end
    end
    ExecutableEquationBlock[ExecutableEquationBlock(torque.name, :load,
        [torque.torque_equation], residual!, jacobian!)]
end

function applied_force_value(force::PlanarAppliedForceComponent, t, z)
    direction = directed_axis_values(force.direction_axis, z).unit
    magnitude = !force.active[] ? zero(eltype(z)) :
        force.magnitude_variable == 0 ? force.magnitude(t) :
            z[force.magnitude_variable]
    magnitude .* direction
end

applied_torque_value(torque::PlanarAppliedTorqueComponent, t, z) =
    !torque.active[] ? zero(eltype(z)) :
        torque.torque_variable == 0 ? torque.torque(t) :
            z[torque.torque_variable]

function add_applied_force!(equations, force::PlanarAppliedForceComponent,
        t, z)
    global_force = applied_force_value(force, t, z)
    application = PlanarAppliedForces.point_marker_kinematics(
        force.application_marker, z)
    PlanarAppliedForces.add_body_point_force!(equations,
        force.application_marker, application, global_force)
    if !isnothing(force.reaction_marker)
        reaction = PlanarAppliedForces.point_marker_kinematics(
            force.reaction_marker, z)
        PlanarAppliedForces.add_body_point_force!(equations,
            force.reaction_marker, reaction, -global_force)
    end
    return nothing
end

function applied_force_dependencies(marker::PlanarBodyPointMarker)
    [collect(marker.position_variables); marker.theta_variable]
end

applied_force_dependencies(::PlanarGroundPointMarker) = Int[]

function applied_force_dependencies(marker::PlanarFloatingPointMarker)
    unique([collect(marker.owner_position_variables);
            applied_force_dependencies(marker.follower)])
end

function applied_force_target_rows(marker::PlanarBodyPointMarker)
    [collect(marker.force_equations); marker.torque_equation]
end

applied_force_target_rows(::PlanarGroundPointMarker) = Int[]

function applied_force_target_rows(marker::PlanarFloatingPointMarker)
    [collect(marker.force_equations); marker.torque_equation]
end

function equation_contributions(force::PlanarAppliedForceComponent)
    target_rows = unique([applied_force_target_rows(force.application_marker);
                          isnothing(force.reaction_marker) ? Int[] :
                              applied_force_target_rows(force.reaction_marker)])
    dependencies = unique([
        applied_force_dependencies(force.application_marker);
        isnothing(force.reaction_marker) ? Int[] :
            applied_force_dependencies(force.reaction_marker);
        directed_axis_dependencies(force.direction_axis);
        force.magnitude_variable == 0 ? Int[] : [force.magnitude_variable]])
    residual! = (equations, t, z, zdot) ->
        add_applied_force!(equations, force, t, z)
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        add_contribution_finite_difference!(jacobian, target_rows,
            (equations, trial) -> add_applied_force!(
                equations, force, t, trial), z, dependencies)
    end
    return EquationContribution[
        EquationContribution(force.name, :force_to_bodies,
            target_rows, residual!, jacobian!),
    ]
end

function equation_contributions(torque::PlanarConstantTorqueComponent)
    residual! = function (equations, t, z, zdot)
        PlanarAppliedForces.add_marker_torque!(
            equations, torque.marker_a, torque.torque)
        PlanarAppliedForces.add_marker_torque!(
            equations, torque.marker_b, -torque.torque)
    end
    target_rows = Int[]
    for marker in (torque.marker_a, torque.marker_b)
        marker isa PlanarBodyOrientationMarker || continue
        push!(target_rows, marker.torque_equation)
    end
    return EquationContribution[
        EquationContribution(torque.name, :torque_to_bodies,
            unique(target_rows), residual!, (J, t, z, zd, c) -> nothing),
    ]
end

function equation_contributions(torque::PlanarAppliedTorqueComponent)
    residual! = function (equations, t, z, zdot)
        value = applied_torque_value(torque, t, z)
        PlanarAppliedForces.add_marker_torque!(equations, torque.marker_a, value)
        PlanarAppliedForces.add_marker_torque!(equations, torque.marker_b, -value)
    end
    target_rows = Int[]
    for marker in (torque.marker_a, torque.marker_b)
        marker isa PlanarBodyOrientationMarker || continue
        push!(target_rows, marker.torque_equation)
    end
    jacobian! = if torque.torque_variable == 0
        (J, t, z, zd, c) -> nothing
    else
        function (jacobian, t, z, zdot, coefficient)
            PlanarAppliedForces.add_marker_torque_jacobian!(jacobian,
                torque.marker_a, torque.torque_variable, 1)
            PlanarAppliedForces.add_marker_torque_jacobian!(jacobian,
                torque.marker_b, torque.torque_variable, -1)
        end
    end
    return EquationContribution[
        EquationContribution(torque.name, :torque_to_bodies,
            unique(target_rows), residual!, jacobian!),
    ]
end

function equation_contributions(joint::PlanarRevoluteJointComponent)
    residual! = function (equations, t, z, zdot)
        reaction = z[joint.reaction_variables]
        values_a = PlanarAppliedForces.point_marker_kinematics(joint.marker_a, z)
        values_b = PlanarAppliedForces.point_marker_kinematics(joint.marker_b, z)
        PlanarAppliedForces.add_body_point_force!(equations, joint.marker_a,
            values_a, reaction)
        PlanarAppliedForces.add_body_point_force!(equations, joint.marker_b,
            values_b, -reaction)
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        reaction = z[joint.reaction_variables]
        values_a = PlanarAppliedForces.point_marker_kinematics(joint.marker_a, z)
        values_b = PlanarAppliedForces.point_marker_kinematics(joint.marker_b, z)
        PlanarAppliedForces.add_point_force_partials!(jacobian,
            joint.marker_a, 1, values_a, joint.reaction_variables, reaction)
        PlanarAppliedForces.add_point_force_partials!(jacobian,
            joint.marker_b, -1, values_b, joint.reaction_variables, reaction)
    end
    target_rows = Int[]
    for marker in (joint.marker_a, joint.marker_b)
        marker isa PlanarBodyPointMarker || continue
        append!(target_rows, marker.force_equations)
        push!(target_rows, marker.torque_equation)
    end
    return EquationContribution[
        EquationContribution(joint.name, :reaction_to_bodies,
            unique(target_rows),
            residual!, jacobian!),
    ]
end

directed_geometry(component) = component.geometry

function inplane_direction(component, z)
    axis = directed_axis_values(directed_geometry(component).axis, z)
    (; unit = axis.unit, normal = axis.transverse,
       omega = axis.omega, alpha = axis.alpha)
end

function inplane_values(component, z)
    values = directed_distance_values(directed_geometry(component), z)
    direction = (; unit = values.axis.unit,
                 normal = values.axis.transverse,
                 omega = values.axis.omega,
                 alpha = values.axis.alpha)
    (; values.marker_i, values.marker_j, direction, values.separation,
       values.relative_velocity, values.relative_acceleration,
       unit_rate = values.axis.unit_rate,
       unit_acceleration = values.axis.unit_acceleration)
end

inplane_position(component, z) =
    directed_distance_position(directed_geometry(component), z)
inplane_velocity(component, z) =
    directed_distance_velocity(directed_geometry(component), z)
inplane_acceleration(component, z) =
    directed_distance_acceleration(directed_geometry(component), z)

function body_dependency_indices(body::PlanarRigidBodyComponent)
    return [collect(body.acceleration_variables);
        body.angular_acceleration_variable;
        collect(body.velocity_variables);
        body.angular_velocity_variable;
        collect(body.position_variables);
        body.orientation_variable]
end
body_dependency_indices(::Nothing) = Int[]

function inplane_dependency_indices(constraint)
    directed_distance_dependencies(directed_geometry(constraint))
end

function add_scalar_finite_difference!(jacobian, row, scalar_function,
        z, dependencies)
    for column in dependencies
        step = sqrt(eps(real(float(one(eltype(z)))))) *
            max(abs(z[column]), one(eltype(z)))
        z_plus = copy(z)
        z_minus = copy(z)
        z_plus[column] += step
        z_minus[column] -= step
        jacobian[row, column] +=
            (scalar_function(z_plus) - scalar_function(z_minus)) / (2step)
    end
end

function add_contribution_finite_difference!(jacobian, rows,
        residual!, z, dependencies)
    for column in dependencies
        step = sqrt(eps(real(float(one(eltype(z)))))) *
            max(abs(z[column]), one(eltype(z)))
        plus, minus = copy(z), copy(z)
        plus[column] += step
        minus[column] -= step
        equations_plus = zeros(eltype(z), size(jacobian, 1))
        equations_minus = zeros(eltype(z), size(jacobian, 1))
        residual!(equations_plus, plus)
        residual!(equations_minus, minus)
        for row in rows
            jacobian[row, column] +=
                (equations_plus[row] - equations_minus[row]) / (2step)
        end
    end
end

function executable_blocks(constraint::PlanarInplaneConstraint)
    dependencies = inplane_dependency_indices(constraint)
    functions = (inplane_acceleration, inplane_velocity, inplane_position)
    names = (:acceleration, :velocity, :position)
    rows = (constraint.acceleration_equation,
            constraint.velocity_equation,
            constraint.position_equation)
    blocks = ExecutableEquationBlock[]
    for (name, row, scalar_function) in zip(names, rows, functions)
        residual! = (equations, t, z, zdot) ->
            (equations[row] = scalar_function(constraint, z))
        jacobian! = function (jacobian, t, z, zdot, coefficient)
            add_scalar_finite_difference!(jacobian, row,
                trial -> scalar_function(constraint, trial), z, dependencies)
        end
        push!(blocks, ExecutableEquationBlock(constraint.name, name,
            [row], residual!, jacobian!))
    end
    return blocks
end

function add_inplane_reaction!(equations, constraint, z)
    values = inplane_values(constraint, z)
    force = z[constraint.reaction_variable] .* values.direction.unit
    geometry = directed_geometry(constraint)
    PlanarAppliedForces.add_body_point_force!(equations,
        geometry.marker_i, values.marker_i, force)
    PlanarAppliedForces.add_body_point_force!(equations,
        geometry.marker_j, values.marker_j, -force)
end

function equation_contributions(constraint::PlanarInplaneConstraint)
    target_rows = Int[]
    geometry = directed_geometry(constraint)
    for marker in (geometry.marker_i, geometry.marker_j)
        marker isa PlanarBodyPointMarker || continue
        append!(target_rows, marker.force_equations)
        push!(target_rows, marker.torque_equation)
    end
    target_rows = unique(target_rows)
    residual! = (equations, t, z, zdot) ->
        add_inplane_reaction!(equations, constraint, z)
    dependencies = unique([inplane_dependency_indices(constraint);
                           constraint.reaction_variable])
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        add_contribution_finite_difference!(jacobian, target_rows,
            (equations, trial) ->
                add_inplane_reaction!(equations, constraint, trial),
            z, dependencies)
    end
    return EquationContribution[
        EquationContribution(constraint.name, :reaction_to_bodies,
            target_rows, residual!, jacobian!),
    ]
end

perp_marker_alpha(body::PlanarRigidBodyComponent, z) =
    z[body.angular_acceleration_variable]
perp_marker_alpha(::Nothing, z) = zero(eltype(z))

function perp_relative_angle(constraint::PlanarPerpConstraint, z)
    PlanarAppliedForces.marker_angle(constraint.orientation_i, z) -
        PlanarAppliedForces.marker_angle(constraint.orientation_j, z)
end

perp_position(constraint::PlanarPerpConstraint, z) =
    sin(perp_relative_angle(constraint, z))

perp_velocity(constraint::PlanarPerpConstraint, z) =
    PlanarAppliedForces.marker_angular_velocity(
        constraint.orientation_i, z) -
    PlanarAppliedForces.marker_angular_velocity(
        constraint.orientation_j, z)

perp_acceleration(constraint::PlanarPerpConstraint, z) =
    perp_marker_alpha(constraint.body_i, z) -
    perp_marker_alpha(constraint.body_j, z)

function executable_blocks(constraint::PlanarPerpConstraint)
    position! = (equations, t, z, zdot) ->
        (equations[constraint.position_equation] =
            perp_position(constraint, z))
    velocity! = (equations, t, z, zdot) ->
        (equations[constraint.velocity_equation] =
            perp_velocity(constraint, z))
    acceleration! = (equations, t, z, zdot) ->
        (equations[constraint.acceleration_equation] =
            perp_acceleration(constraint, z))
    position_jacobian! = function (jacobian, t, z, zdot, coefficient)
        cosine = cos(perp_relative_angle(constraint, z))
        add_motion_position_partial!(jacobian, constraint.position_equation,
            constraint.orientation_i, cosine)
        add_motion_position_partial!(jacobian, constraint.position_equation,
            constraint.orientation_j, -cosine)
    end
    velocity_jacobian! = function (jacobian, t, z, zdot, coefficient)
        add_motion_velocity_partial!(jacobian, constraint.velocity_equation,
            constraint.orientation_i, 1)
        add_motion_velocity_partial!(jacobian, constraint.velocity_equation,
            constraint.orientation_j, -1)
    end
    acceleration_jacobian! = function (jacobian, t, z, zdot, coefficient)
        add_motion_acceleration_partial!(jacobian,
            constraint.acceleration_equation, constraint.body_i, 1)
        add_motion_acceleration_partial!(jacobian,
            constraint.acceleration_equation, constraint.body_j, -1)
    end
    ExecutableEquationBlock[
        ExecutableEquationBlock(constraint.name, :acceleration,
            [constraint.acceleration_equation], acceleration!,
            acceleration_jacobian!),
        ExecutableEquationBlock(constraint.name, :velocity,
            [constraint.velocity_equation], velocity!, velocity_jacobian!),
        ExecutableEquationBlock(constraint.name, :position,
            [constraint.position_equation], position!, position_jacobian!),
    ]
end

function add_perp_reaction!(equations, constraint::PlanarPerpConstraint, z)
    torque = z[constraint.reaction_variable]
    PlanarAppliedForces.add_marker_torque!(
        equations, constraint.orientation_i, torque)
    PlanarAppliedForces.add_marker_torque!(
        equations, constraint.orientation_j, -torque)
    return nothing
end

function equation_contributions(constraint::PlanarPerpConstraint)
    residual! = (equations, t, z, zdot) ->
        add_perp_reaction!(equations, constraint, z)
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        PlanarAppliedForces.add_marker_torque_jacobian!(jacobian,
            constraint.orientation_i, constraint.reaction_variable, 1)
        PlanarAppliedForces.add_marker_torque_jacobian!(jacobian,
            constraint.orientation_j, constraint.reaction_variable, -1)
    end
    target_rows = Int[]
    for marker in (constraint.orientation_i, constraint.orientation_j)
        marker isa PlanarBodyOrientationMarker || continue
        push!(target_rows, marker.torque_equation)
    end
    EquationContribution[EquationContribution(constraint.name,
        :reaction_to_bodies, unique(target_rows), residual!, jacobian!)]
end

executable_blocks(joint::PlanarTranslationalJoint) =
    [executable_blocks(joint.inplane); executable_blocks(joint.perp)]

equation_contributions(joint::PlanarTranslationalJoint) =
    [equation_contributions(joint.inplane);
     equation_contributions(joint.perp)]

executable_blocks(joint::PlanarFixedJoint) =
    [executable_blocks(joint.revolute); executable_blocks(joint.perp)]

equation_contributions(joint::PlanarFixedJoint) =
    [equation_contributions(joint.revolute);
     equation_contributions(joint.perp)]

gear_marker_alpha(body::PlanarRigidBodyComponent, z) =
    z[body.angular_acceleration_variable]
gear_marker_alpha(::Nothing, z) = zero(eltype(z))

gear_carrier_angle(body::PlanarRigidBodyComponent, z) =
    z[body.orientation_variable]
gear_carrier_angle(::Nothing, z) = zero(eltype(z))

gear_coefficients(gear::PlanarGearPairComponent) =
    (gear.radius_1, -gear.radius_2)

function gear_position(gear::PlanarGearPairComponent, z)
    coefficient_1, coefficient_2 = gear_coefficients(gear)
    carrier_angle = gear_carrier_angle(gear.carrier_body, z)
    return coefficient_1 * (PlanarAppliedForces.marker_angle(
        gear.orientation_1, z) - carrier_angle) +
        coefficient_2 * (PlanarAppliedForces.marker_angle(
            gear.orientation_2, z) - carrier_angle) - gear.phase
end

function gear_velocity(gear::PlanarGearPairComponent, z)
    coefficient_1, coefficient_2 = gear_coefficients(gear)
    carrier_velocity = PlanarAppliedForces.marker_angular_velocity(
        gear.carrier_orientation, z)
    return coefficient_1 * (PlanarAppliedForces.marker_angular_velocity(
        gear.orientation_1, z) - carrier_velocity) +
        coefficient_2 * (PlanarAppliedForces.marker_angular_velocity(
            gear.orientation_2, z) - carrier_velocity)
end

function gear_acceleration(gear::PlanarGearPairComponent, z)
    coefficient_1, coefficient_2 = gear_coefficients(gear)
    carrier_acceleration = gear_marker_alpha(gear.carrier_body, z)
    return coefficient_1 * (gear_marker_alpha(gear.body_1, z) -
        carrier_acceleration) +
        coefficient_2 * (gear_marker_alpha(gear.body_2, z) -
        carrier_acceleration)
end

function executable_blocks(gear::PlanarGearPairComponent)
    coefficient_1, coefficient_2 = gear_coefficients(gear)
    carrier_coefficient = -(coefficient_1 + coefficient_2)
    position! = (equations, t, z, zdot) ->
        (equations[gear.position_equation] = gear_position(gear, z))
    velocity! = (equations, t, z, zdot) ->
        (equations[gear.velocity_equation] = gear_velocity(gear, z))
    acceleration! = (equations, t, z, zdot) ->
        (equations[gear.acceleration_equation] = gear_acceleration(gear, z))
    position_jacobian! = function (jacobian, t, z, zdot, coefficient)
        add_motion_position_partial!(jacobian, gear.position_equation,
            gear.orientation_1, coefficient_1)
        add_motion_position_partial!(jacobian, gear.position_equation,
            gear.orientation_2, coefficient_2)
        add_motion_position_partial!(jacobian, gear.position_equation,
            gear.carrier_orientation, carrier_coefficient)
    end
    velocity_jacobian! = function (jacobian, t, z, zdot, coefficient)
        add_motion_velocity_partial!(jacobian, gear.velocity_equation,
            gear.orientation_1, coefficient_1)
        add_motion_velocity_partial!(jacobian, gear.velocity_equation,
            gear.orientation_2, coefficient_2)
        add_motion_velocity_partial!(jacobian, gear.velocity_equation,
            gear.carrier_orientation, carrier_coefficient)
    end
    acceleration_jacobian! = function (jacobian, t, z, zdot, coefficient)
        add_motion_acceleration_partial!(jacobian, gear.acceleration_equation,
            gear.body_1, coefficient_1)
        add_motion_acceleration_partial!(jacobian, gear.acceleration_equation,
            gear.body_2, coefficient_2)
        add_motion_acceleration_partial!(jacobian, gear.acceleration_equation,
            gear.carrier_body, carrier_coefficient)
    end
    ExecutableEquationBlock[
        ExecutableEquationBlock(gear.name, :acceleration,
            [gear.acceleration_equation], acceleration!,
            acceleration_jacobian!),
        ExecutableEquationBlock(gear.name, :velocity,
            [gear.velocity_equation], velocity!, velocity_jacobian!),
        ExecutableEquationBlock(gear.name, :position,
            [gear.position_equation], position!, position_jacobian!),
    ]
end

"""Return the unit direction of the gear-pair force acting on the first gear."""
function gear_force_direction(gear::PlanarGearPairComponent, z)
    contact = PlanarAppliedForces.point_marker_kinematics(
        gear.contact_marker, z).position
    reference = PlanarAppliedForces.point_marker_kinematics(
        gear.force_reference_marker, z).position
    radial = contact - reference
    radius = norm(radial)
    radius > 0 || throw(ArgumentError(
        "gear pair '$(gear.name)' contact point coincides with its force-reference joint"))
    [-radial[2], radial[1]] ./ radius
end

function add_gear_reaction!(equations, gear, z)
    marker_1 = PlanarAppliedForces.point_marker_kinematics(gear.marker_1, z)
    marker_2 = PlanarAppliedForces.point_marker_kinematics(gear.marker_2, z)
    force = z[gear.reaction_variable] .* gear_force_direction(gear, z)
    PlanarAppliedForces.add_body_point_force!(
        equations, gear.marker_1, marker_1, force)
    PlanarAppliedForces.add_body_point_force!(
        equations, gear.marker_2, marker_2, -force)
    return nothing
end

function equation_contributions(gear::PlanarGearPairComponent)
    target_rows = unique([
        collect(gear.marker_1.force_equations);
        gear.marker_1.torque_equation;
        collect(gear.marker_2.force_equations);
        gear.marker_2.torque_equation;
    ])
    residual! = (equations, t, z, zdot) ->
        add_gear_reaction!(equations, gear, z)
    dependencies = unique([
        body_dependency_indices(gear.body_1);
        body_dependency_indices(gear.body_2);
        body_dependency_indices(gear.carrier_body);
        gear.reaction_variable;
    ])
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        add_contribution_finite_difference!(jacobian, target_rows,
            (equations, trial) -> add_gear_reaction!(equations, gear, trial),
            z, dependencies)
    end
    EquationContribution[EquationContribution(gear.name,
        :reaction_to_bodies, target_rows, residual!, jacobian!)]
end

function rack_translation_values(component::PlanarRackAndPinionComponent, z)
    values = inplane_values(component.translational_joint.inplane, z)
    axis = -values.direction.normal
    axis_rate = values.direction.unit .* values.direction.omega
    axis_acceleration = values.direction.unit .* values.direction.alpha .-
        axis .* values.direction.omega^2
    position = dot(values.separation, axis)
    velocity = dot(values.relative_velocity, axis) +
        dot(values.separation, axis_rate)
    acceleration = dot(values.relative_acceleration, axis) +
        2 * dot(values.relative_velocity, axis_rate) +
        dot(values.separation, axis_acceleration)
    return (; position, velocity, acceleration, axis)
end

function rack_and_pinion_position(component::PlanarRackAndPinionComponent, z)
    rack = rack_translation_values(component, z)
    relative_angle = PlanarAppliedForces.marker_angle(
        component.pinion_orientation, z) -
        PlanarAppliedForces.marker_angle(
            component.pinion_carrier_orientation, z)
    rack.position + component.pitch_radius * relative_angle - component.phase
end

function rack_and_pinion_velocity(component::PlanarRackAndPinionComponent, z)
    rack = rack_translation_values(component, z)
    relative_omega = PlanarAppliedForces.marker_angular_velocity(
        component.pinion_orientation, z) -
        PlanarAppliedForces.marker_angular_velocity(
            component.pinion_carrier_orientation, z)
    rack.velocity + component.pitch_radius * relative_omega
end

function rack_and_pinion_acceleration(
        component::PlanarRackAndPinionComponent, z)
    rack = rack_translation_values(component, z)
    relative_alpha = gear_marker_alpha(component.pinion_body, z) -
        gear_marker_alpha(component.carrier_body, z)
    rack.acceleration + component.pitch_radius * relative_alpha
end

function executable_blocks(component::PlanarRackAndPinionComponent)
    dependencies = unique([
        body_dependency_indices(component.rack_body);
        body_dependency_indices(component.pinion_body);
        body_dependency_indices(component.carrier_body);
    ])
    functions = (rack_and_pinion_acceleration, rack_and_pinion_velocity,
                 rack_and_pinion_position)
    names = (:acceleration, :velocity, :position)
    rows = (component.acceleration_equation, component.velocity_equation,
            component.position_equation)
    blocks = ExecutableEquationBlock[]
    for (name, row, scalar_function) in zip(names, rows, functions)
        residual! = (equations, t, z, zdot) ->
            (equations[row] = scalar_function(component, z))
        jacobian! = function (jacobian, t, z, zdot, coefficient)
            add_scalar_finite_difference!(jacobian, row,
                trial -> scalar_function(component, trial), z, dependencies)
        end
        push!(blocks, ExecutableEquationBlock(component.name, name, [row],
            residual!, jacobian!))
    end
    blocks
end

function add_rack_and_pinion_reaction!(equations, component, z)
    rack = rack_translation_values(component, z)
    rack_kinematics = PlanarAppliedForces.point_marker_kinematics(
        component.rack_marker, z)
    pinion_kinematics = PlanarAppliedForces.point_marker_kinematics(
        component.pinion_marker, z)
    force = z[component.reaction_variable] .* rack.axis
    PlanarAppliedForces.add_body_point_force!(equations,
        component.rack_marker, rack_kinematics, force)
    PlanarAppliedForces.add_body_point_force!(equations,
        component.pinion_marker, pinion_kinematics, -force)
    return nothing
end

function equation_contributions(component::PlanarRackAndPinionComponent)
    target_rows = unique([
        collect(component.rack_marker.force_equations);
        component.rack_marker.torque_equation;
        collect(component.pinion_marker.force_equations);
        component.pinion_marker.torque_equation;
    ])
    residual! = (equations, t, z, zdot) ->
        add_rack_and_pinion_reaction!(equations, component, z)
    dependencies = unique([
        body_dependency_indices(component.rack_body);
        body_dependency_indices(component.pinion_body);
        body_dependency_indices(component.carrier_body);
        component.reaction_variable;
    ])
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        add_contribution_finite_difference!(jacobian, target_rows,
            (equations, trial) ->
                add_rack_and_pinion_reaction!(equations, component, trial),
            z, dependencies)
    end
    EquationContribution[EquationContribution(component.name,
        :reaction_to_bodies, target_rows, residual!, jacobian!)]
end

function pulley_center_kinematics(pulley::PlanarPulleyComponent, z)
    PlanarAppliedForces.point_marker_kinematics(pulley.center_marker, z)
end

"""
Return one selected common tangent of two pitch circles.

`tangent_kind` is `1` for an external tangent and `-1` for an internal
tangent. `tangent_side` selects one of the two solutions. The normals point
from each pulley center to its tangent point.
"""
function belt_tangent_geometry(center_1, radius_1, center_2, radius_2,
        tangent_kind::Integer, tangent_side::Integer)
    tangent_kind in (-1, 1) || throw(ArgumentError(
        "belt tangent kind must be -1 or 1"))
    tangent_side in (-1, 1) || throw(ArgumentError(
        "belt tangent side must be -1 or 1"))
    separation = center_2 - center_1
    center_distance = norm(separation)
    center_distance > 0 || throw(DomainError(center_distance,
        "belt pulley centers must be distinct"))
    center_unit = separation ./ center_distance
    center_perpendicular = [-center_unit[2], center_unit[1]]
    a = (radius_1 - tangent_kind * radius_2) / center_distance
    abs(a) <= 1 || throw(DomainError(a,
        "selected belt tangent is geometrically infeasible"))
    b = tangent_side * sqrt(max(zero(a), one(a) - a^2))
    normal_1 = a .* center_unit .+ b .* center_perpendicular
    normal_2 = tangent_kind .* normal_1
    point_1 = center_1 .+ radius_1 .* normal_1
    point_2 = center_2 .+ radius_2 .* normal_2
    span = point_2 - point_1
    length = norm(span)
    length > 0 || throw(DomainError(length,
        "belt tangent span must have positive length"))
    tangent = span ./ length
    beta_1 = radius_1 * dot([-normal_1[2], normal_1[1]], tangent)
    beta_2 = radius_2 * dot([-normal_2[2], normal_2[1]], tangent)
    (; point_1, point_2, normal_1, normal_2, tangent, length,
       beta_1, beta_2)
end

function belt_tangent_geometry(span::PlanarBeltSpanComponent, z)
    center_1 = pulley_center_kinematics(span.pulley_1, z)
    center_2 = pulley_center_kinematics(span.pulley_2, z)
    geometry = belt_tangent_geometry(center_1.position,
        span.pulley_1.pitch_radius, center_2.position,
        span.pulley_2.pitch_radius, span.tangent_kind, span.tangent_side)
    (; geometry..., center_1, center_2)
end

cross2(a, b) = a[1] * b[2] - a[2] * b[1]

"""Return actual tangent geometry and the span's explicit local variables."""
function belt_span_values(span::PlanarBeltSpanComponent, z)
    geometry = belt_tangent_geometry(span, z)
    point_1 = z[span.point_1_variables]
    point_2 = z[span.point_2_variables]
    tangent = z[span.tangent_variables]
    length = z[span.length_variable]
    extension = z[span.extension_variable]
    extension_rate = z[span.extension_rate_variable]
    tension = z[span.tension_variable]
    force = z[span.force_variables]
    tangent_angle_change = atan(cross2(span.reference_tangent,
        geometry.tangent), dot(span.reference_tangent, geometry.tangent))
    angle_1 = z[span.pulley_1.body.orientation_variable]
    angle_2 = z[span.pulley_2.body.orientation_variable]
    expected_extension = geometry.length - span.reference_length -
        (span.beta_2 - span.beta_1) * tangent_angle_change +
        span.beta_2 * (angle_2 - span.reference_angle_2) -
        span.beta_1 * (angle_1 - span.reference_angle_1) +
        span.initial_extension
    surface_speed_1 = dot(geometry.center_1.velocity, geometry.tangent) +
        span.beta_1 * z[span.pulley_1.body.angular_velocity_variable]
    surface_speed_2 = dot(geometry.center_2.velocity, geometry.tangent) +
        span.beta_2 * z[span.pulley_2.body.angular_velocity_variable]
    expected_extension_rate = surface_speed_2 - surface_speed_1
    (; actual_point_1 = geometry.point_1,
       actual_point_2 = geometry.point_2,
       actual_normal_1 = geometry.normal_1,
       actual_normal_2 = geometry.normal_2,
       actual_tangent = geometry.tangent,
       actual_length = geometry.length,
       center_1 = geometry.center_1, center_2 = geometry.center_2,
       point_1, point_2, tangent, length, extension,
       extension_rate, tension, force, expected_extension,
       expected_extension_rate, surface_speed_1, surface_speed_2)
end

function belt_span_geometry_residual!(equations, span, z)
    values = belt_span_values(span, z)
    rows = span.geometry_equations
    equations[rows[1:2]] .= values.point_1 .- values.actual_point_1
    equations[rows[3:4]] .= values.point_2 .- values.actual_point_2
    equations[rows[5:6]] .= values.tangent .- values.actual_tangent
    equations[rows[7]] = values.length - values.actual_length
    equations[rows[8]] = values.extension - values.expected_extension
    nothing
end

function belt_span_rate_residual!(equations, span, z)
    values = belt_span_values(span, z)
    equations[span.rate_equation] =
        values.extension_rate - values.expected_extension_rate
    nothing
end

function belt_span_load_residual!(equations, span, z)
    values = belt_span_values(span, z)
    equations[span.load_equations[1]] = values.tension -
        span.stiffness * values.extension -
        span.damping * values.extension_rate
    equations[span.load_equations[2:3]] .=
        values.force .- values.tension .* values.tangent
    nothing
end

function belt_block_jacobian!(jacobian, rows, residual!, z, dependencies)
    add_contribution_finite_difference!(jacobian, rows,
        (equations, trial) -> residual!(equations, trial), z, dependencies)
end

function executable_blocks(span::PlanarBeltSpanComponent)
    geometry_body_dependencies = unique([
        collect(span.pulley_1.body.position_variables);
        span.pulley_1.body.orientation_variable;
        collect(span.pulley_2.body.position_variables);
        span.pulley_2.body.orientation_variable])
    rate_body_dependencies = unique([geometry_body_dependencies;
        collect(span.pulley_1.body.velocity_variables);
        span.pulley_1.body.angular_velocity_variable;
        collect(span.pulley_2.body.velocity_variables);
        span.pulley_2.body.angular_velocity_variable])
    geometry_dependencies = unique([geometry_body_dependencies;
        collect(span.point_1_variables); collect(span.point_2_variables);
        collect(span.tangent_variables); span.length_variable;
        span.extension_variable])
    rate_dependencies = unique([rate_body_dependencies;
        span.extension_rate_variable])
    load_dependencies = unique([collect(span.tangent_variables);
        span.extension_variable; span.extension_rate_variable;
        span.tension_variable; collect(span.force_variables)])
    geometry! = (equations, t, z, zdot) ->
        belt_span_geometry_residual!(equations, span, z)
    rate! = (equations, t, z, zdot) ->
        belt_span_rate_residual!(equations, span, z)
    load! = (equations, t, z, zdot) ->
        belt_span_load_residual!(equations, span, z)
    geometry_jacobian! = (jacobian, t, z, zdot, coefficient) ->
        belt_block_jacobian!(jacobian, collect(span.geometry_equations),
            (equations, trial) -> belt_span_geometry_residual!(
                equations, span, trial), z, geometry_dependencies)
    rate_jacobian! = (jacobian, t, z, zdot, coefficient) ->
        belt_block_jacobian!(jacobian, [span.rate_equation],
            (equations, trial) -> belt_span_rate_residual!(
                equations, span, trial), z, rate_dependencies)
    load_jacobian! = (jacobian, t, z, zdot, coefficient) ->
        belt_block_jacobian!(jacobian, collect(span.load_equations),
            (equations, trial) -> belt_span_load_residual!(
                equations, span, trial), z, load_dependencies)
    ExecutableEquationBlock[
        ExecutableEquationBlock(span.name, :geometry,
            collect(span.geometry_equations), geometry!, geometry_jacobian!),
        ExecutableEquationBlock(span.name, :rate, [span.rate_equation],
            rate!, rate_jacobian!),
        ExecutableEquationBlock(span.name, :load,
            collect(span.load_equations), load!, load_jacobian!),
    ]
end

function add_belt_point_force!(equations, pulley, z, point, force)
    body = pulley.body
    equations[body.balance_equations[1:2]] .-= force
    radius = point - z[body.position_variables]
    equations[body.balance_equations[3]] -= cross2(radius, force)
    nothing
end

function equation_contributions(span::PlanarBeltSpanComponent)
    target_rows = unique([collect(span.pulley_1.body.balance_equations);
                          collect(span.pulley_2.body.balance_equations)])
    residual! = function (equations, t, z, zdot)
        values = belt_span_values(span, z)
        add_belt_point_force!(equations, span.pulley_1, z,
            values.point_1, values.force)
        add_belt_point_force!(equations, span.pulley_2, z,
            values.point_2, -values.force)
    end
    dependencies = unique([collect(span.pulley_1.body.position_variables);
                           collect(span.pulley_2.body.position_variables);
                           collect(span.point_1_variables);
                           collect(span.point_2_variables);
                           collect(span.force_variables)])
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        add_contribution_finite_difference!(jacobian, target_rows,
            (equations, trial) -> residual!(equations, t, trial, zdot),
            z, dependencies)
    end
    EquationContribution[EquationContribution(span.name, :load_to_bodies,
        target_rows, residual!, jacobian!)]
end

executable_blocks(::PlanarPulleyComponent) = ExecutableEquationBlock[]
equation_contributions(::PlanarPulleyComponent) = EquationContribution[]
executable_blocks(::PlanarBeltComponent) = ExecutableEquationBlock[]
equation_contributions(::PlanarBeltComponent) = EquationContribution[]

motion_marker_angle(marker, z) = PlanarAppliedForces.marker_angle(marker, z)
motion_marker_omega(marker, z) =
    PlanarAppliedForces.marker_angular_velocity(marker, z)
motion_marker_alpha(body::PlanarRigidBodyComponent, z) =
    z[body.angular_acceleration_variable]
motion_marker_alpha(::Nothing, z) = zero(eltype(z))

function motion_relative_angle(generator, z)
    angle_a = motion_marker_angle(generator.marker_a, z)
    angle_b = motion_marker_angle(generator.marker_b, z)
    x = cos(angle_b) * cos(angle_a) + sin(angle_b) * sin(angle_a)
    y = cos(angle_b) * sin(angle_a) - sin(angle_b) * cos(angle_a)
    return atan(y, x)
end

function add_motion_position_partial!(jacobian, row,
        marker::PlanarBodyOrientationMarker, coefficient)
    jacobian[row, marker.theta_variable] += coefficient
end
add_motion_position_partial!(jacobian, row,
    marker::PlanarGroundOrientationMarker, coefficient) = nothing

function add_motion_velocity_partial!(jacobian, row,
        marker::PlanarBodyOrientationMarker, coefficient)
    jacobian[row, marker.omega_variable] += coefficient
end
add_motion_velocity_partial!(jacobian, row,
    marker::PlanarGroundOrientationMarker, coefficient) = nothing

function add_motion_acceleration_partial!(jacobian, row,
        body::PlanarRigidBodyComponent, coefficient)
    jacobian[row, body.angular_acceleration_variable] += coefficient
end
add_motion_acceleration_partial!(jacobian, row, ::Nothing, coefficient) = nothing

function executable_blocks(generator::PlanarRotationalMotionGenerator)
    position! = function (equations, t, z, zdot)
        relative_angle = motion_relative_angle(generator, z)
        equations[generator.position_equations[1]] =
            atan(sin(z[generator.angle_variable] - relative_angle),
                 cos(z[generator.angle_variable] - relative_angle))
        equations[generator.position_equations[2]] =
            atan(sin(z[generator.angle_variable] - generator.motion(t)),
                 cos(z[generator.angle_variable] - generator.motion(t)))
    end
    position_jacobian! = function (jacobian, t, z, zdot, coefficient)
        relation_row, prescription_row = generator.position_equations
        jacobian[relation_row, generator.angle_variable] += 1
        add_motion_position_partial!(jacobian, relation_row,
            generator.marker_a, -1)
        add_motion_position_partial!(jacobian, relation_row,
            generator.marker_b, 1)
        jacobian[prescription_row, generator.angle_variable] += 1
    end
    velocity! = function (equations, t, z, zdot)
        relative_omega = motion_marker_omega(generator.marker_a, z) -
            motion_marker_omega(generator.marker_b, z)
        equations[generator.velocity_equations[1]] =
            z[generator.angular_velocity_variable] - relative_omega
        equations[generator.velocity_equations[2]] =
            z[generator.angular_velocity_variable] -
            generator.motion_derivative(t)
    end
    velocity_jacobian! = function (jacobian, t, z, zdot, coefficient)
        relation_row, prescription_row = generator.velocity_equations
        jacobian[relation_row, generator.angular_velocity_variable] += 1
        add_motion_velocity_partial!(jacobian, relation_row,
            generator.marker_a, -1)
        add_motion_velocity_partial!(jacobian, relation_row,
            generator.marker_b, 1)
        jacobian[prescription_row, generator.angular_velocity_variable] += 1
    end
    acceleration! = function (equations, t, z, zdot)
        relative_alpha = motion_marker_alpha(generator.body_a, z) -
            motion_marker_alpha(generator.body_b, z)
        equations[generator.acceleration_equations[1]] =
            z[generator.angular_acceleration_variable] - relative_alpha
        equations[generator.acceleration_equations[2]] =
            z[generator.angular_acceleration_variable] -
            generator.motion_second_derivative(t)
    end
    acceleration_jacobian! = function (jacobian, t, z, zdot, coefficient)
        relation_row, prescription_row = generator.acceleration_equations
        jacobian[relation_row, generator.angular_acceleration_variable] += 1
        add_motion_acceleration_partial!(jacobian, relation_row,
            generator.body_a, -1)
        add_motion_acceleration_partial!(jacobian, relation_row,
            generator.body_b, 1)
        jacobian[prescription_row,
                 generator.angular_acceleration_variable] += 1
    end
    return ExecutableEquationBlock[
        ExecutableEquationBlock(generator.name, :position,
            collect(generator.position_equations), position!,
            position_jacobian!),
        ExecutableEquationBlock(generator.name, :velocity,
            collect(generator.velocity_equations), velocity!,
            velocity_jacobian!),
        ExecutableEquationBlock(generator.name, :acceleration,
            collect(generator.acceleration_equations), acceleration!,
            acceleration_jacobian!),
    ]
end

function equation_contributions(generator::PlanarRotationalMotionGenerator)
    residual! = function (equations, t, z, zdot)
        torque = z[generator.torque_variable]
        PlanarAppliedForces.add_marker_torque!(
            equations, generator.marker_a, torque)
        PlanarAppliedForces.add_marker_torque!(
            equations, generator.marker_b, -torque)
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        PlanarAppliedForces.add_marker_torque_jacobian!(jacobian,
            generator.marker_a, generator.torque_variable, 1)
        PlanarAppliedForces.add_marker_torque_jacobian!(jacobian,
            generator.marker_b, generator.torque_variable, -1)
    end
    target_rows = Int[]
    for marker in (generator.marker_a, generator.marker_b)
        marker isa PlanarBodyOrientationMarker || continue
        push!(target_rows, marker.torque_equation)
    end
    return EquationContribution[
        EquationContribution(generator.name, :torque_to_bodies,
            unique(target_rows), residual!, jacobian!),
    ]
end

function executable_blocks(generator::PlanarTranslationalMotionGenerator)
    dependencies = inplane_dependency_indices(generator)
    scalar_functions = (inplane_position, inplane_velocity,
                        inplane_acceleration)
    laws = (generator.motion, generator.motion_derivative,
            generator.motion_second_derivative)
    variables = (generator.distance_variable, generator.velocity_variable,
                 generator.acceleration_variable)
    rows = (generator.position_equations, generator.velocity_equations,
            generator.acceleration_equations)
    names = (:position, :velocity, :acceleration)
    blocks = ExecutableEquationBlock[]
    for (name, scalar_function, law, variable, block_rows) in
            zip(names, scalar_functions, laws, variables, rows)
        relation_row, prescription_row = block_rows
        residual! = function (equations, t, z, zdot)
            equations[relation_row] =
                z[variable] - scalar_function(generator, z)
            equations[prescription_row] = z[variable] - law(t)
        end
        jacobian! = function (jacobian, t, z, zdot, coefficient)
            jacobian[relation_row, variable] += 1
            add_scalar_finite_difference!(jacobian, relation_row,
                trial -> -scalar_function(generator, trial), z, dependencies)
            jacobian[prescription_row, variable] += 1
        end
        push!(blocks, ExecutableEquationBlock(generator.name, name,
            collect(block_rows), residual!, jacobian!))
    end
    blocks
end

function equation_contributions(generator::PlanarTranslationalMotionGenerator)
    target_rows = Int[]
    geometry = directed_geometry(generator)
    for marker in (geometry.marker_i, geometry.marker_j)
        marker isa PlanarBodyPointMarker || continue
        append!(target_rows, marker.force_equations)
        push!(target_rows, marker.torque_equation)
    end
    target_rows = unique(target_rows)
    residual! = (equations, t, z, zdot) ->
        add_inplane_reaction!(equations, generator, z)
    dependencies = unique([inplane_dependency_indices(generator);
                           generator.reaction_variable])
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        add_contribution_finite_difference!(jacobian, target_rows,
            (equations, trial) ->
                add_inplane_reaction!(equations, generator, trial),
            z, dependencies)
    end
    EquationContribution[EquationContribution(generator.name,
        :reaction_to_bodies, target_rows, residual!, jacobian!)]
end

coordinate_variables(joint::PlanarRevoluteJointComponent) =
    (joint.rotation_variables[3], joint.rotation_variables[2],
     joint.rotation_variables[1])
coordinate_variables(coordinate::PlanarDistanceCoordinateComponent) =
    (coordinate.distance_variable, coordinate.velocity_variable,
     coordinate.acceleration_variable)

function executable_blocks(coordinate::PlanarDistanceCoordinateComponent)
    dependencies = inplane_dependency_indices(coordinate)
    scalar_functions = (inplane_position, inplane_velocity,
                        inplane_acceleration)
    variables = (coordinate.distance_variable, coordinate.velocity_variable,
                 coordinate.acceleration_variable)
    rows = (coordinate.position_equation, coordinate.velocity_equation,
            coordinate.acceleration_equation)
    names = (:position, :velocity, :acceleration)
    blocks = ExecutableEquationBlock[]
    for (name, scalar_function, variable, row) in
            zip(names, scalar_functions, variables, rows)
        residual! = (equations, t, z, zdot) ->
            (equations[row] = z[variable] - scalar_function(coordinate, z))
        jacobian! = function (jacobian, t, z, zdot, coefficient)
            jacobian[row, variable] += 1
            add_scalar_finite_difference!(jacobian, row,
                trial -> -scalar_function(coordinate, trial), z,
                dependencies)
        end
        push!(blocks, ExecutableEquationBlock(coordinate.name, name,
            [row], residual!, jacobian!))
    end
    state_equations = isempty(coordinate.candidate_state_equations) ?
        coordinate.state_equations : coordinate.candidate_state_equations
    if !isempty(state_equations)
        acceleration = coordinate.acceleration_variable
        velocity = coordinate.velocity_variable
        distance = coordinate.distance_variable
        state! = function (equations, t, z, zdot)
            equations[state_equations[1]] = z[acceleration] - zdot[velocity]
            equations[state_equations[2]] = z[velocity] - zdot[distance]
        end
        state_jacobian! = function (jacobian, t, z, zdot, coefficient)
            jacobian[state_equations[1], acceleration] += 1
            jacobian[state_equations[1], velocity] -= coefficient
            jacobian[state_equations[2], velocity] += 1
            jacobian[state_equations[2], distance] -= coefficient
        end
        push!(blocks, ExecutableEquationBlock(coordinate.name,
            :selected_distance_state, copy(state_equations), state!,
            state_jacobian!))
    end
    blocks
end

equation_contributions(::PlanarDistanceCoordinateComponent) =
    EquationContribution[]

function add_coordinate_reaction!(equations,
        joint::PlanarRevoluteJointComponent, generalized_force, z)
    PlanarAppliedForces.add_marker_torque!(equations,
        joint.rotation_marker_a, generalized_force)
    PlanarAppliedForces.add_marker_torque!(equations,
        joint.rotation_marker_b, -generalized_force)
end

function add_coordinate_reaction!(equations,
        coordinate::PlanarDistanceCoordinateComponent, generalized_force, z)
    values = inplane_values(coordinate, z)
    force = generalized_force .* values.direction.unit
    geometry = directed_geometry(coordinate)
    PlanarAppliedForces.add_body_point_force!(equations,
        geometry.marker_i, values.marker_i, force)
    PlanarAppliedForces.add_body_point_force!(equations,
        geometry.marker_j, values.marker_j, -force)
end

function coordinate_target_rows(joint::PlanarRevoluteJointComponent)
    rows = Int[]
    for marker in (joint.rotation_marker_a, joint.rotation_marker_b)
        marker isa PlanarBodyOrientationMarker || continue
        push!(rows, marker.torque_equation)
    end
    rows
end

function coordinate_target_rows(coordinate::PlanarDistanceCoordinateComponent)
    rows = Int[]
    geometry = directed_geometry(coordinate)
    for marker in (geometry.marker_i, geometry.marker_j)
        marker isa PlanarBodyPointMarker || continue
        append!(rows, marker.force_equations)
        push!(rows, marker.torque_equation)
    end
    rows
end

coordinate_dependency_indices(joint::PlanarRevoluteJointComponent) =
    unique([body_dependency_indices(joint.body_a);
            body_dependency_indices(joint.body_b)])
coordinate_dependency_indices(coordinate::PlanarDistanceCoordinateComponent) =
    inplane_dependency_indices(coordinate)

function executable_blocks(coupler::PlanarCoordinateCoupler)
    rows = (coupler.position_equation, coupler.velocity_equation,
            coupler.acceleration_equation)
    variable_groups = coordinate_variables.(coupler.coordinates)
    columns = ([group[1] for group in variable_groups],
               [group[2] for group in variable_groups],
               [group[3] for group in variable_groups])
    offsets = (coupler.offset, zero(coupler.offset), zero(coupler.offset))
    names = (:position, :velocity, :acceleration)
    blocks = ExecutableEquationBlock[]
    for (name, row, variables, offset) in zip(names, rows, columns, offsets)
        residual! = function (equations, t, z, zdot)
            equations[row] = dot(coupler.coefficients, z[variables]) - offset
        end
        jacobian! = function (jacobian, t, z, zdot, coefficient)
            for (variable, factor) in zip(variables, coupler.coefficients)
                jacobian[row, variable] += factor
            end
        end
        push!(blocks, ExecutableEquationBlock(coupler.name, name,
            [row], residual!, jacobian!))
    end
    blocks
end

function equation_contributions(coupler::PlanarCoordinateCoupler)
    target_rows = unique(collect(Iterators.flatten(
        coordinate_target_rows(coordinate)
        for coordinate in coupler.coordinates)))
    residual! = function (equations, t, z, zdot)
        reaction = z[coupler.reaction_variable]
        for (coordinate, factor) in
                zip(coupler.coordinates, coupler.coefficients)
            add_coordinate_reaction!(equations, coordinate,
                factor * reaction, z)
        end
    end
    jacobian! = if coupler.coordinate_kind == :rotation
        function (jacobian, t, z, zdot, coefficient)
            for (joint, factor) in
                    zip(coupler.coordinates, coupler.coefficients)
                PlanarAppliedForces.add_marker_torque_jacobian!(jacobian,
                    joint.rotation_marker_a, coupler.reaction_variable,
                    factor)
                PlanarAppliedForces.add_marker_torque_jacobian!(jacobian,
                    joint.rotation_marker_b, coupler.reaction_variable,
                    -factor)
            end
        end
    else
        dependencies = unique([coupler.reaction_variable;
            collect(Iterators.flatten(
                coordinate_dependency_indices(coordinate)
                for coordinate in coupler.coordinates))])
        function (jacobian, t, z, zdot, coefficient)
            add_contribution_finite_difference!(jacobian, target_rows,
                (equations, trial) -> residual!(equations, t, trial, zdot),
                z, dependencies)
        end
    end
    EquationContribution[EquationContribution(coupler.name,
        :reaction_to_coordinates, target_rows, residual!, jacobian!)]
end

function planar_span_values(span, z)
    marker_1 = PlanarAppliedForces.point_marker_kinematics(span.marker_1, z)
    marker_2 = PlanarAppliedForces.point_marker_kinematics(span.marker_2, z)
    length = z[span.length_variable]
    length > zero(length) || throw(DomainError(length,
        "planar span distance must be positive"))
    relative_velocity = marker_2.velocity - marker_1.velocity
    values = (; marker_1, marker_2,
       spanning = @view(z[span.spanning_variables]), length,
       unit = @view(z[span.unit_variables]),
       relative_velocity,
       length_rate = z[span.length_rate_variable])
    hasproperty(span, :force_variable) || return values
    (; values...,
       scalar_force = z[span.force_variable],
       global_force = @view(z[span.global_force_variables]))
end

spanning_values(force, z) = planar_span_values(force.element, z)

function planar_point_acceleration(body, marker_values, z)
    isnothing(body) && return zeros(eltype(z), 2)
    z[body.acceleration_variables] .+
        marker_values.d .* z[body.angular_acceleration_variable] .-
        marker_values.r .* z[body.angular_velocity_variable]^2
end

function planar_span_measure_values(measure, z)
    values = planar_span_values(measure, z)
    relative_acceleration =
        planar_point_acceleration(measure.body_2, values.marker_2, z) -
        planar_point_acceleration(measure.body_1, values.marker_1, z)
    transverse_speed_squared =
        dot(values.relative_velocity, values.relative_velocity) -
        values.length_rate^2
    (; values..., relative_acceleration,
       length_acceleration = z[measure.length_acceleration_variable],
       transverse_speed_squared)
end

function planar_span_geometry_block(span; name = nothing,
        block_name = :geometry)
    component_name = isnothing(name) ? span.name : name
    geometry_rows = [collect(span.spanning_equations);
                     span.length_equation; collect(span.unit_equations)]
    geometry! = function (equations, t, z, zdot)
        values = planar_span_values(span, z)
        equations[span.spanning_equations] .= values.spanning .-
            (values.marker_2.position - values.marker_1.position)
        equations[span.length_equation] = values.length - norm(values.spanning)
        equations[span.unit_equations] .= values.unit .-
            values.spanning ./ values.length
    end
    geometry_jacobian! = function (jacobian, t, z, zdot, coefficient)
        values = planar_span_values(span, z)
        s = span.spanning_variables
        ell = span.length_variable
        for k in 1:2
            jacobian[span.spanning_equations[k], s[k]] += 1
        end
        PlanarAppliedForces.add_spanning_marker_partials!(jacobian,
            span.spanning_equations, span.marker_1, 1, values.marker_1)
        PlanarAppliedForces.add_spanning_marker_partials!(jacobian,
            span.spanning_equations, span.marker_2, -1, values.marker_2)
        jacobian[span.length_equation, ell] += 1
        jacobian[span.length_equation, s] .-=
            values.spanning ./ norm(values.spanning)
        for k in 1:2
            row = span.unit_equations[k]
            jacobian[row, span.unit_variables[k]] += 1
            jacobian[row, s[k]] -= 1 / values.length
            jacobian[row, ell] += values.spanning[k] / values.length^2
        end
    end
    ExecutableEquationBlock(component_name, block_name, geometry_rows,
        geometry!, geometry_jacobian!)
end

function planar_span_velocity_block(span; name = nothing,
        block_name = :velocity)
    component_name = isnothing(name) ? span.name : name
    velocity! = function (equations, t, z, zdot)
        values = planar_span_values(span, z)
        equations[span.length_rate_equation] = values.length_rate -
            dot(values.unit, values.relative_velocity)
    end
    velocity_jacobian! = function (jacobian, t, z, zdot, coefficient)
        values = planar_span_values(span, z)
        row = span.length_rate_equation
        jacobian[row, span.length_rate_variable] += 1
        jacobian[row, span.unit_variables] .-= values.relative_velocity
        omega_1 = span.marker_1 isa PlanarBodyPointMarker ?
            z[span.marker_1.omega_variable] : zero(eltype(z))
        omega_2 = span.marker_2 isa PlanarBodyPointMarker ?
            z[span.marker_2.omega_variable] : zero(eltype(z))
        PlanarAppliedForces.add_length_rate_marker_partials!(jacobian, row,
            span.marker_1, -1, values.marker_1, values.unit, omega_1)
        PlanarAppliedForces.add_length_rate_marker_partials!(jacobian, row,
            span.marker_2, 1, values.marker_2, values.unit, omega_2)
    end
    ExecutableEquationBlock(component_name, block_name,
        [span.length_rate_equation], velocity!, velocity_jacobian!)
end

function executable_blocks(force::PlanarSpanningForceComponent)
    element = force.element
    load! = function (equations, t, z, zdot)
        values = spanning_values(force, z)
        equations[element.force_equation] = values.scalar_force -
            (force.active[] ? force.law(t, z) : 0.0)
        equations[element.global_force_equations] .= values.global_force .+
            values.unit .* values.scalar_force
    end
    load_jacobian! = function (jacobian, t, z, zdot, coefficient)
        values = spanning_values(force, z)
        jacobian[element.force_equation, element.force_variable] += 1
        if force.active[]
            partials = force.law.gradient(t, z)
            for (column, partial) in zip(force.law.dependencies, partials)
                jacobian[element.force_equation, column] -= partial
            end
        end
        for k in 1:2
            row = element.global_force_equations[k]
            jacobian[row, element.global_force_variables[k]] += 1
            jacobian[row, element.unit_variables[k]] += values.scalar_force
            jacobian[row, element.force_variable] += values.unit[k]
        end
    end
    return ExecutableEquationBlock[
        planar_span_geometry_block(element; name = force.name),
        planar_span_velocity_block(element; name = force.name,
            block_name = :rate),
        ExecutableEquationBlock(force.name, :load,
            [element.force_equation; collect(element.global_force_equations)],
            load!, load_jacobian!),
    ]
end

function add_span_acceleration_marker_partials!(jacobian, row, body,
        marker_values, values, relative_sign, z)
    isnothing(body) && return nothing
    unit = values.unit
    relative_velocity = values.relative_velocity
    length = values.length
    omega = z[body.angular_velocity_variable]
    alpha = z[body.angular_acceleration_variable]
    jacobian[row, body.acceleration_variables] .-= relative_sign .* unit
    jacobian[row, body.angular_acceleration_variable] -=
        relative_sign * dot(unit, marker_values.d)
    jacobian[row, body.velocity_variables] .-=
        relative_sign .* (2 .* relative_velocity ./ length)
    omega_partial = dot(unit, -2omega .* marker_values.r) +
        2dot(relative_velocity, marker_values.d) / length
    jacobian[row, body.angular_velocity_variable] -=
        relative_sign * omega_partial
    theta_partial = dot(unit,
        -alpha .* marker_values.r .- omega^2 .* marker_values.d) -
        2omega * dot(relative_velocity, marker_values.r) / length
    jacobian[row, body.orientation_variable] -=
        relative_sign * theta_partial
    nothing
end

function executable_blocks(measure::PlanarSpanMeasureComponent)
    acceleration! = function (equations, t, z, zdot)
        values = planar_span_measure_values(measure, z)
        equations[measure.length_acceleration_equation] =
            values.length_acceleration -
            dot(values.unit, values.relative_acceleration) -
            values.transverse_speed_squared / values.length
    end
    acceleration_jacobian! = function (jacobian, t, z, zdot, coefficient)
        values = planar_span_measure_values(measure, z)
        row = measure.length_acceleration_equation
        jacobian[row, measure.length_acceleration_variable] += 1
        jacobian[row, measure.unit_variables] .-=
            values.relative_acceleration
        jacobian[row, measure.length_rate_variable] +=
            2values.length_rate / values.length
        jacobian[row, measure.length_variable] +=
            values.transverse_speed_squared / values.length^2
        add_span_acceleration_marker_partials!(jacobian, row,
            measure.body_1, values.marker_1, values, -1, z)
        add_span_acceleration_marker_partials!(jacobian, row,
            measure.body_2, values.marker_2, values, 1, z)
    end
    ExecutableEquationBlock[
        planar_span_geometry_block(measure),
        planar_span_velocity_block(measure),
        ExecutableEquationBlock(measure.name, :acceleration,
            [measure.length_acceleration_equation], acceleration!,
            acceleration_jacobian!),
    ]
end

equation_contributions(::PlanarSpanMeasureComponent) = EquationContribution[]

function equation_contributions(force::PlanarSpanningForceComponent)
    element = force.element
    residual! = function (equations, t, z, zdot)
        values = spanning_values(force, z)
        PlanarAppliedForces.add_body_point_force!(equations, element.marker_1,
            values.marker_1, values.global_force)
        PlanarAppliedForces.add_body_point_force!(equations, element.marker_2,
            values.marker_2, -values.global_force)
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        values = spanning_values(force, z)
        PlanarAppliedForces.add_point_force_partials!(jacobian,
            element.marker_1, 1, values.marker_1,
            element.global_force_variables, values.global_force)
        PlanarAppliedForces.add_point_force_partials!(jacobian,
            element.marker_2, -1, values.marker_2,
            element.global_force_variables, values.global_force)
    end
    target_rows = Int[]
    for marker in (element.marker_1, element.marker_2)
        marker isa PlanarBodyPointMarker || continue
        append!(target_rows, marker.force_equations)
        push!(target_rows, marker.torque_equation)
    end
    return EquationContribution[
        EquationContribution(force.name, :load_to_bodies,
            unique(target_rows), residual!, jacobian!),
    ]
end

function executable_blocks(spring::PlanarTorsionalSpringComponent)
    element = spring.element
    residual! = function (equations, t, z, zdot)
        relative = PlanarAppliedForces.relative_rotation(element, z)
        equations[element.constitutive_equation] =
            z[element.torque_variable] +
            element.stiffness * (relative.angle - element.free_angle) +
            element.damping * relative.angular_velocity
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        row = element.constitutive_equation
        jacobian[row, element.torque_variable] += 1
        PlanarAppliedForces.add_constitutive_partials!(jacobian,
            element.marker_1, row, element.stiffness, element.damping)
        PlanarAppliedForces.add_constitutive_partials!(jacobian,
            element.marker_2, row, -element.stiffness, -element.damping)
    end
    ExecutableEquationBlock[ExecutableEquationBlock(spring.name, :load,
        [element.constitutive_equation], residual!, jacobian!)]
end

orientation_torque_rows(marker::PlanarBodyOrientationMarker) =
    [marker.torque_equation]
orientation_torque_rows(::PlanarGroundOrientationMarker) = Int[]

function equation_contributions(spring::PlanarTorsionalSpringComponent)
    element = spring.element
    residual! = function (equations, t, z, zdot)
        torque = z[element.torque_variable]
        PlanarAppliedForces.add_marker_torque!(
            equations, element.marker_1, torque)
        PlanarAppliedForces.add_marker_torque!(
            equations, element.marker_2, -torque)
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        PlanarAppliedForces.add_marker_torque_jacobian!(jacobian,
            element.marker_1, element.torque_variable, 1)
        PlanarAppliedForces.add_marker_torque_jacobian!(jacobian,
            element.marker_2, element.torque_variable, -1)
    end
    target_rows = unique([
        orientation_torque_rows(element.marker_1);
        orientation_torque_rows(element.marker_2);
    ])
    EquationContribution[EquationContribution(spring.name,
        :torque_to_bodies, target_rows, residual!, jacobian!)]
end

function bushing_values(bushing::PlanarBushingComponent, z)
    marker_1 = PlanarAppliedForces.point_marker_kinematics(bushing.marker_1.point, z)
    marker_2 = PlanarAppliedForces.point_marker_kinematics(bushing.marker_2.point, z)
    angle_2 = PlanarAppliedForces.marker_angle(bushing.marker_2.orientation, z)
    c, s = cos(angle_2), sin(angle_2)
    rotation = [c -s; s c]
    relative_position = transpose(rotation) *
        (marker_1.position - marker_2.position)
    omega_2 = PlanarAppliedForces.marker_angular_velocity(
        bushing.marker_2.orientation, z)
    relative_velocity = transpose(rotation) *
        (marker_1.velocity - marker_2.velocity) .-
        omega_2 .* [-relative_position[2], relative_position[1]]
    angle_1 = PlanarAppliedForces.marker_angle(bushing.marker_1.orientation, z)
    relative_angle = atan(sin(angle_1 - angle_2), cos(angle_1 - angle_2))
    angle_deformation = atan(sin(relative_angle - bushing.free_angle),
        cos(relative_angle - bushing.free_angle))
    relative_omega =
        PlanarAppliedForces.marker_angular_velocity(
            bushing.marker_1.orientation, z) - omega_2
    local_force = -bushing.translational_stiffness .*
        (relative_position - bushing.free_position) .-
        bushing.translational_damping .* relative_velocity
    global_force = bushing.active[] ? rotation * local_force :
        zeros(eltype(z), 2)
    torque = bushing.active[] ?
        -bushing.rotational_stiffness * angle_deformation -
            bushing.rotational_damping * relative_omega : zero(eltype(z))
    (; marker_1, marker_2, relative_position, relative_velocity,
       relative_angle, relative_omega, global_force, torque)
end

function bushing_marker_dependencies(marker::PlanarBodyPointMarker)
    [collect(marker.position_variables); marker.theta_variable;
     collect(marker.velocity_variables); marker.omega_variable]
end
bushing_marker_dependencies(::PlanarGroundPointMarker) = Int[]

function executable_blocks(bushing::PlanarBushingComponent)
    residual! = function (equations, t, z, zdot)
        values = bushing_values(bushing, z)
        equations[bushing.load_equations[1:2]] .=
            z[bushing.force_variables] .- values.global_force
        equations[bushing.load_equations[3]] =
            z[bushing.torque_variable] - values.torque
    end
    dependencies = unique([bushing_marker_dependencies(bushing.marker_1.point);
                           bushing_marker_dependencies(bushing.marker_2.point)])
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        for k in 1:2
            jacobian[bushing.load_equations[k], bushing.force_variables[k]] += 1
        end
        jacobian[bushing.load_equations[3], bushing.torque_variable] += 1
        for column in dependencies
            step = sqrt(eps(real(float(one(eltype(z)))))) *
                max(abs(z[column]), one(eltype(z)))
            plus, minus = copy(z), copy(z)
            plus[column] += step
            minus[column] -= step
            plus_values = bushing_values(bushing, plus)
            minus_values = bushing_values(bushing, minus)
            force_partial = (plus_values.global_force -
                minus_values.global_force) / (2step)
            torque_partial = (plus_values.torque - minus_values.torque) / (2step)
            jacobian[bushing.load_equations[1:2], column] .-= force_partial
            jacobian[bushing.load_equations[3], column] -= torque_partial
        end
    end
    ExecutableEquationBlock[ExecutableEquationBlock(bushing.name, :load,
        collect(bushing.load_equations), residual!, jacobian!)]
end

function equation_contributions(bushing::PlanarBushingComponent)
    residual! = function (equations, t, z, zdot)
        values = bushing_values(bushing, z)
        force = z[bushing.force_variables]
        torque = z[bushing.torque_variable]
        PlanarAppliedForces.add_body_point_force!(equations,
            bushing.marker_1.point, values.marker_1, force)
        PlanarAppliedForces.add_body_point_force!(equations,
            bushing.marker_2.point, values.marker_2, -force)
        PlanarAppliedForces.add_marker_torque!(equations,
            bushing.marker_1.orientation, torque)
        PlanarAppliedForces.add_marker_torque!(equations,
            bushing.marker_2.orientation, -torque)
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        values = bushing_values(bushing, z)
        force = z[bushing.force_variables]
        PlanarAppliedForces.add_point_force_partials!(jacobian,
            bushing.marker_1.point, 1, values.marker_1,
            bushing.force_variables, force)
        PlanarAppliedForces.add_point_force_partials!(jacobian,
            bushing.marker_2.point, -1, values.marker_2,
            bushing.force_variables, force)
        PlanarAppliedForces.add_marker_torque_jacobian!(jacobian,
            bushing.marker_1.orientation, bushing.torque_variable, 1)
        PlanarAppliedForces.add_marker_torque_jacobian!(jacobian,
            bushing.marker_2.orientation, bushing.torque_variable, -1)
    end
    target_rows = Int[]
    for marker in (bushing.marker_1, bushing.marker_2)
        marker.point isa PlanarBodyPointMarker || continue
        append!(target_rows, marker.point.force_equations)
        push!(target_rows, marker.point.torque_equation)
    end
    EquationContribution[EquationContribution(bushing.name, :load_to_bodies,
        unique(target_rows), residual!, jacobian!)]
end

function plane_contact_values(contact::PlanarPlaneContactComponent, z)
    geometry = directed_distance_values(
        contact.geometry, z; include_acceleration = false)
    marker_1, marker_2 = geometry.marker_i, geometry.marker_j
    normal = geometry.axis.unit
    gap = geometry.position - contact.radius
    gap_rate = geometry.velocity
    penetration = max(-gap, zero(gap))
    damping_multiplier = max(zero(gap),
        one(gap) - contact.damping_factor * gap_rate)
    normal_force = contact.active[] ?
        contact.stiffness * penetration * damping_multiplier : zero(gap)
    global_force = normal_force .* normal
    (; marker_1, marker_2, normal, gap, gap_rate, normal_force, global_force)
end

plane_contact_gap(contact::PlanarPlaneContactComponent, z) =
    z[contact.gap_variable]

function plane_contact_damping_surface(contact::PlanarPlaneContactComponent, z)
    contact.active[] || return one(eltype(z))
    gap = z[contact.gap_variable]
    return gap < 0 ?
        one(gap) - contact.damping_factor * z[contact.gap_rate_variable] :
        one(gap)
end

function executable_blocks(contact::PlanarPlaneContactComponent)
    residual! = function (equations, t, z, zdot)
        values = plane_contact_values(contact, z)
        rows = contact.contact_equations
        equations[rows[1]] = z[contact.gap_variable] - values.gap
        equations[rows[2]] = z[contact.gap_rate_variable] - values.gap_rate
        penetration = max(-z[contact.gap_variable], zero(eltype(z)))
        damping_multiplier = max(zero(eltype(z)), one(eltype(z)) -
            contact.damping_factor * z[contact.gap_rate_variable])
        target_force = contact.active[] ?
            contact.stiffness * penetration * damping_multiplier :
            zero(eltype(z))
        equations[rows[3]] = z[contact.normal_force_variable] - target_force
        equations[rows[4:5]] .= z[contact.global_force_variables] .-
            z[contact.normal_force_variable] .* values.normal
    end
    dependencies = directed_distance_dependencies(contact.geometry)
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        rows = contact.contact_equations
        jacobian[rows[1], contact.gap_variable] += 1
        jacobian[rows[2], contact.gap_rate_variable] += 1
        jacobian[rows[3], contact.normal_force_variable] += 1
        # Touch this structural entry in both contact states so a sparse
        # prototype assembled while separated remains valid after impact.
        penetration = max(-z[contact.gap_variable], zero(eltype(z)))
        raw_multiplier = one(eltype(z)) -
            contact.damping_factor * z[contact.gap_rate_variable]
        damping_multiplier = max(zero(eltype(z)), raw_multiplier)
        active = contact.active[] && z[contact.gap_variable] < 0 &&
            raw_multiplier > 0
        jacobian[rows[3], contact.gap_variable] +=
            active ? contact.stiffness * damping_multiplier : 0
        jacobian[rows[3], contact.gap_rate_variable] += active ?
            contact.stiffness * penetration * contact.damping_factor : 0
        for k in 1:2
            jacobian[rows[3 + k], contact.global_force_variables[k]] += 1
        end
        values = plane_contact_values(contact, z)
        jacobian[rows[4:5], contact.normal_force_variable] .-= values.normal
        for column in dependencies
            step = sqrt(eps(real(float(one(eltype(z)))))) *
                max(abs(z[column]), one(eltype(z)))
            plus, minus = copy(z), copy(z)
            plus[column] += step
            minus[column] -= step
            plus_values = plane_contact_values(contact, plus)
            minus_values = plane_contact_values(contact, minus)
            jacobian[rows[1], column] -=
                (plus_values.gap - minus_values.gap) / (2step)
            jacobian[rows[2], column] -=
                (plus_values.gap_rate - minus_values.gap_rate) / (2step)
            normal_partial = (plus_values.normal - minus_values.normal) / (2step)
            jacobian[rows[4:5], column] .-=
                z[contact.normal_force_variable] .* normal_partial
        end
    end
    ExecutableEquationBlock[ExecutableEquationBlock(contact.name, :contact,
        collect(contact.contact_equations), residual!, jacobian!)]
end

function equation_contributions(contact::PlanarPlaneContactComponent)
    residual! = function (equations, t, z, zdot)
        values = plane_contact_values(contact, z)
        force = z[contact.global_force_variables]
        PlanarAppliedForces.add_body_point_force!(equations,
            contact.marker_1.point, values.marker_1, force)
        PlanarAppliedForces.add_body_point_force!(equations,
            contact.marker_2.point, values.marker_2, -force)
    end
    jacobian! = function (jacobian, t, z, zdot, coefficient)
        values = plane_contact_values(contact, z)
        force = z[contact.global_force_variables]
        PlanarAppliedForces.add_point_force_partials!(jacobian,
            contact.marker_1.point, 1, values.marker_1,
            contact.global_force_variables, force)
        PlanarAppliedForces.add_point_force_partials!(jacobian,
            contact.marker_2.point, -1, values.marker_2,
            contact.global_force_variables, force)
    end
    target_rows = Int[]
    for marker in (contact.marker_1, contact.marker_2)
        marker.point isa PlanarBodyPointMarker || continue
        append!(target_rows, marker.point.force_equations)
        push!(target_rows, marker.point.torque_equation)
    end
    EquationContribution[EquationContribution(contact.name, :load_to_bodies,
        unique(target_rows), residual!, jacobian!)]
end

end
