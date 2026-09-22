"""
    SpatialComponentAssembly

Component-local equations for the spatial rigid-body modeler. The physical
orientation is a rotation matrix evaluated from normalized Euler parameters.
Body-fixed pseudo angles accompany selected angular-velocity states but are
not used to evaluate finite orientation.
"""
module SpatialComponentAssembly

using LinearAlgebra
using ..AutomaticAnalysis

export SpatialRigidBodyComponent, SpatialFlexibleBeamComponent,
       SpatialGravityComponent,
       spatial_body_registration, spatial_flexible_beam_registration,
       allocated_spatial_body,
       component_registration, executable_blocks, equation_contributions,
       skew, axis_angle_rotation, quaternion_rate_matrix, rotation_matrix,
       rotation_vector_jacobian, rotation_transpose_vector_jacobian,
       matrix_to_euler_parameters

"""
Allocated spatial rigid body in the unreduced canonical system.

Positions and translational vectors use global components. Angular velocity,
angular acceleration, inertia, and pseudo-angle increments use the body
reference frame. `euler_parameter_variables` hold a normalized scalar-first
quaternion that maps body components into the global frame. Acceleration stays
explicit and the three pseudo-angle state equations provide local mechanical
orientation columns without replacing the finite orientation parameters.
"""
struct SpatialRigidBodyComponent{T}
    name::Symbol
    mass::T
    inertia::Matrix{T}
    acceleration_variables::UnitRange{Int}
    angular_acceleration_variables::UnitRange{Int}
    velocity_variables::UnitRange{Int}
    angular_velocity_variables::UnitRange{Int}
    position_variables::UnitRange{Int}
    pseudo_angle_variables::UnitRange{Int}
    euler_parameter_variables::UnitRange{Int}
    balance_equations::UnitRange{Int}
    acceleration_state_equations::UnitRange{Int}
    angular_acceleration_state_equations::UnitRange{Int}
    position_state_equations::UnitRange{Int}
    pseudo_angle_state_equations::UnitRange{Int}
    orientation_equations::UnitRange{Int}
end

"""
Two-node spatial Timoshenko beam with a floating rigid reference frame.

The reference frame carries the finite translation and orientation used by a
spatial rigid body. Six small elastic coordinates describe axial extension,
two transverse deflections, torsion, and two bending rotations after the six
rigid modes have been removed from the ordinary 12-coordinate beam element.
"""
struct SpatialFlexibleBeamComponent{T}
    name::Symbol
    mass::T
    inertia::Matrix{T}
    length::T
    elastic_mass::Matrix{T}
    elastic_stiffness::Matrix{T}
    elastic_damping::Matrix{T}
    deformation_shape::Matrix{T}
    acceleration_variables::UnitRange{Int}
    angular_acceleration_variables::UnitRange{Int}
    elastic_acceleration_variables::UnitRange{Int}
    velocity_variables::UnitRange{Int}
    angular_velocity_variables::UnitRange{Int}
    elastic_velocity_variables::UnitRange{Int}
    position_variables::UnitRange{Int}
    pseudo_angle_variables::UnitRange{Int}
    euler_parameter_variables::UnitRange{Int}
    elastic_position_variables::UnitRange{Int}
    balance_equations::UnitRange{Int}
    acceleration_state_equations::UnitRange{Int}
    angular_acceleration_state_equations::UnitRange{Int}
    elastic_acceleration_state_equations::UnitRange{Int}
    position_state_equations::UnitRange{Int}
    pseudo_angle_state_equations::UnitRange{Int}
    elastic_position_state_equations::UnitRange{Int}
    orientation_equations::UnitRange{Int}
end

"""Constant global gravitational acceleration applied to one spatial body."""
struct SpatialGravityComponent{B,T}
    name::Symbol
    body::B
    acceleration::Vector{T}
end

"""Return the matrix whose product with another vector is the cross product."""
function skew(vector::AbstractVector)
    length(vector) == 3 || throw(DimensionMismatch(
        "a cross-product matrix requires three components"))
    x, y, z = vector
    [zero(x) -z y; z zero(x) -x; -y x zero(x)]
end

"""Rotation matrix for a right-handed `angle` about `axis`."""
function axis_angle_rotation(angle::Real, axis::AbstractVector)
    length(axis) == 3 || throw(DimensionMismatch(
        "a rotation axis requires three components"))
    values = promote(Float64(angle), Float64.(axis)...)
    theta = first(values)
    direction = collect(values[2:4])
    all(isfinite, values) || throw(ArgumentError(
        "axis-angle orientation must contain finite numbers"))
    magnitude = norm(direction)
    magnitude > 0 || throw(ArgumentError(
        "axis-angle orientation requires a nonzero axis"))
    unit_axis = direction ./ magnitude
    cosine, sine = cos(theta), sin(theta)
    cosine .* Matrix{Float64}(I, 3, 3) .+
        (1 - cosine) .* (unit_axis * transpose(unit_axis)) .+
        sine .* skew(unit_axis)
end

"""Matrix `Q(p)` satisfying `ṗ = Q(p)ωᵇ/2` for a unit quaternion."""
function quaternion_rate_matrix(parameters::AbstractVector)
    length(parameters) == 4 || throw(DimensionMismatch(
        "Euler parameters require four components"))
    scalar = parameters[1]
    vector = parameters[2:4]
    vcat(reshape(-vector, 1, 3),
        scalar .* Matrix{eltype(parameters)}(I, 3, 3) + skew(vector))
end

"""Rotation matrix mapping body-frame components into the global frame."""
function rotation_matrix(parameters::AbstractVector)
    length(parameters) == 4 || throw(DimensionMismatch(
        "Euler parameters require four components"))
    scalar = parameters[1]
    vector = parameters[2:4]
    identity = Matrix{eltype(parameters)}(I, 3, 3)
    (scalar^2 - dot(vector, vector)) .* identity .+
        2 .* (vector * transpose(vector)) .+ 2scalar .* skew(vector)
end

"""Partial of `rotation_matrix(parameters) * vector` with respect to `parameters`."""
function rotation_vector_jacobian(parameters::AbstractVector,
        vector::AbstractVector)
    length(parameters) == 4 || throw(DimensionMismatch(
        "Euler parameters require four components"))
    length(vector) == 3 || throw(DimensionMismatch(
        "a rotated vector requires three components"))
    scalar = parameters[1]
    quaternion_vector = parameters[2:4]
    identity = Matrix{promote_type(eltype(parameters), eltype(vector))}(I, 3, 3)
    scalar_column = 2scalar .* vector .+
        2 .* cross(quaternion_vector, vector)
    vector_columns = -2 .* vector * transpose(quaternion_vector) .+
        2 .* identity .* dot(quaternion_vector, vector) .+
        2 .* quaternion_vector * transpose(vector) .-
        2scalar .* skew(vector)
    hcat(scalar_column, vector_columns)
end

"""Partial of `transpose(rotation_matrix(parameters)) * vector`."""
function rotation_transpose_vector_jacobian(parameters::AbstractVector,
        vector::AbstractVector)
    length(parameters) == 4 || throw(DimensionMismatch(
        "Euler parameters require four components"))
    length(vector) == 3 || throw(DimensionMismatch(
        "a rotated vector requires three components"))
    scalar = parameters[1]
    quaternion_vector = parameters[2:4]
    identity = Matrix{promote_type(eltype(parameters), eltype(vector))}(I, 3, 3)
    scalar_column = 2scalar .* vector .-
        2 .* cross(quaternion_vector, vector)
    vector_columns = -2 .* vector * transpose(quaternion_vector) .+
        2 .* identity .* dot(quaternion_vector, vector) .+
        2 .* quaternion_vector * transpose(vector) .+
        2scalar .* skew(vector)
    hcat(scalar_column, vector_columns)
end

"""Convert a proper orthonormal rotation matrix to scalar-first parameters."""
function matrix_to_euler_parameters(matrix::AbstractMatrix)
    size(matrix) == (3, 3) || throw(DimensionMismatch(
        "orientation must be a 3 by 3 matrix"))
    A = Matrix{Float64}(matrix)
    norm(transpose(A) * A - I, Inf) <= 1.0e-10 || throw(ArgumentError(
        "orientation matrix must be orthonormal"))
    det(A) > 0 || throw(ArgumentError(
        "orientation matrix must be a proper rotation"))
    trace_value = tr(A)
    parameters = if trace_value > 0
        scalar = 0.5 * sqrt(1 + trace_value)
        denominator = 4scalar
        [scalar,
         (A[3, 2] - A[2, 3]) / denominator,
         (A[1, 3] - A[3, 1]) / denominator,
         (A[2, 1] - A[1, 2]) / denominator]
    else
        diagonal = diag(A)
        index = argmax(diagonal)
        if index == 1
            x = 0.5 * sqrt(1 + A[1, 1] - A[2, 2] - A[3, 3])
            denominator = 4x
            [(A[3, 2] - A[2, 3]) / denominator, x,
             (A[1, 2] + A[2, 1]) / denominator,
             (A[1, 3] + A[3, 1]) / denominator]
        elseif index == 2
            y = 0.5 * sqrt(1 - A[1, 1] + A[2, 2] - A[3, 3])
            denominator = 4y
            [(A[1, 3] - A[3, 1]) / denominator,
             (A[1, 2] + A[2, 1]) / denominator, y,
             (A[2, 3] + A[3, 2]) / denominator]
        else
            z = 0.5 * sqrt(1 - A[1, 1] - A[2, 2] + A[3, 3])
            denominator = 4z
            [(A[2, 1] - A[1, 2]) / denominator,
             (A[1, 3] + A[3, 1]) / denominator,
             (A[2, 3] + A[3, 2]) / denominator, z]
        end
    end
    parameters ./= norm(parameters)
    parameters[1] < 0 && (parameters .*= -1)
    parameters
end

function spatial_body_registration(name::Symbol)
    variables = VariableDeclaration[
        VariableDeclaration(:a_x, :acceleration, 2),
        VariableDeclaration(:a_y, :acceleration, 2),
        VariableDeclaration(:a_z, :acceleration, 2),
        VariableDeclaration(:alpha_x, :angular_acceleration, 2),
        VariableDeclaration(:alpha_y, :angular_acceleration, 2),
        VariableDeclaration(:alpha_z, :angular_acceleration, 2),
        VariableDeclaration(:V_x, :velocity, 1),
        VariableDeclaration(:V_y, :velocity, 1),
        VariableDeclaration(:V_z, :velocity, 1),
        VariableDeclaration(:omega_x, :angular_velocity, 1),
        VariableDeclaration(:omega_y, :angular_velocity, 1),
        VariableDeclaration(:omega_z, :angular_velocity, 1),
        VariableDeclaration(:R_x, :position, 0),
        VariableDeclaration(:R_y, :position, 0),
        VariableDeclaration(:R_z, :position, 0),
        VariableDeclaration(:psi_x, :orientation, 0),
        VariableDeclaration(:psi_y, :orientation, 0),
        VariableDeclaration(:psi_z, :orientation, 0),
        VariableDeclaration(:p_0, :orientation_parameter, 0),
        VariableDeclaration(:p_1, :orientation_parameter, 0),
        VariableDeclaration(:p_2, :orientation_parameter, 0),
        VariableDeclaration(:p_3, :orientation_parameter, 0),
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:balance, EquationDeclaration[
            [EquationDeclaration(Symbol(:sum_F_, axis), :balance, 2,
                :body_balance) for axis in (:x, :y, :z)];
            [EquationDeclaration(Symbol(:sum_T_, axis), :balance, 2,
                :body_balance) for axis in (:x, :y, :z)];
        ]),
        EquationBlockDeclaration(:selected_state, EquationDeclaration[
            [EquationDeclaration(Symbol(:a_, axis, :_state),
                :state_equation, 2, :translational_state)
                for axis in (:x, :y, :z)];
            [EquationDeclaration(Symbol(:alpha_, axis, :_state),
                :state_equation, 2, :angular_state)
                for axis in (:x, :y, :z)];
            [EquationDeclaration(Symbol(:V_, axis, :_state),
                :state_equation, 1, :translational_state)
                for axis in (:x, :y, :z)];
            [EquationDeclaration(Symbol(:omega_, axis, :_state),
                :state_equation, 1, :angular_state)
                for axis in (:x, :y, :z)];
        ]),
        EquationBlockDeclaration(:orientation, EquationDeclaration[
            [EquationDeclaration(Symbol(:psi_, axis, :_parameter_kinematic),
                :coordinate_relation, 1, :euler_parameter_kinematics)
                for axis in (:x, :y, :z)];
            EquationDeclaration(:parameter_normalization, :normalization, 0,
                :euler_parameter_normalization);
        ]),
    ]
    ComponentRegistration(name, variables, blocks)
end

function spatial_flexible_beam_registration(name::Symbol)
    elastic_names = (:u, :v, :w, :rx, :ry, :rz)
    variables = VariableDeclaration[
        [VariableDeclaration(Symbol(:a_, axis), :acceleration, 2)
            for axis in (:x, :y, :z)];
        [VariableDeclaration(Symbol(:alpha_, axis), :angular_acceleration, 2)
            for axis in (:x, :y, :z)];
        [VariableDeclaration(Symbol(:eta_, coordinate, :_ddot),
            :elastic_acceleration, 2) for coordinate in elastic_names];
        [VariableDeclaration(Symbol(:V_, axis), :velocity, 1)
            for axis in (:x, :y, :z)];
        [VariableDeclaration(Symbol(:omega_, axis), :angular_velocity, 1)
            for axis in (:x, :y, :z)];
        [VariableDeclaration(Symbol(:eta_, coordinate, :_dot),
            :elastic_velocity, 1) for coordinate in elastic_names];
        [VariableDeclaration(Symbol(:R_, axis), :position, 0)
            for axis in (:x, :y, :z)];
        [VariableDeclaration(Symbol(:psi_, axis), :orientation, 0)
            for axis in (:x, :y, :z)];
        [VariableDeclaration(Symbol(:p_, index), :orientation_parameter, 0)
            for index in 0:3];
        [VariableDeclaration(Symbol(:eta_, coordinate), :elastic_position, 0)
            for coordinate in elastic_names];
    ]
    blocks = EquationBlockDeclaration[
        EquationBlockDeclaration(:balance, EquationDeclaration[
            [EquationDeclaration(Symbol(:sum_F_, axis), :balance, 2,
                :beam_balance) for axis in (:x, :y, :z)];
            [EquationDeclaration(Symbol(:sum_T_, axis), :balance, 2,
                :beam_balance) for axis in (:x, :y, :z)];
            [EquationDeclaration(Symbol(:sum_Q_, coordinate), :balance, 2,
                :beam_balance) for coordinate in elastic_names];
        ]),
        EquationBlockDeclaration(:selected_state, EquationDeclaration[
            [EquationDeclaration(Symbol(:a_, axis, :_state),
                :state_equation, 2, :translational_state)
                for axis in (:x, :y, :z)];
            [EquationDeclaration(Symbol(:alpha_, axis, :_state),
                :state_equation, 2, :angular_state)
                for axis in (:x, :y, :z)];
            [EquationDeclaration(Symbol(:eta_, coordinate, :_acceleration_state),
                :state_equation, 2, :elastic_state)
                for coordinate in elastic_names];
            [EquationDeclaration(Symbol(:V_, axis, :_state),
                :state_equation, 1, :translational_state)
                for axis in (:x, :y, :z)];
            [EquationDeclaration(Symbol(:omega_, axis, :_state),
                :state_equation, 1, :angular_state)
                for axis in (:x, :y, :z)];
            [EquationDeclaration(Symbol(:eta_, coordinate, :_position_state),
                :state_equation, 1, :elastic_state)
                for coordinate in elastic_names];
        ]),
        EquationBlockDeclaration(:orientation, EquationDeclaration[
            [EquationDeclaration(Symbol(:psi_, axis, :_parameter_kinematic),
                :coordinate_relation, 1, :euler_parameter_kinematics)
                for axis in (:x, :y, :z)];
            EquationDeclaration(:parameter_normalization, :normalization, 0,
                :euler_parameter_normalization);
        ]),
    ]
    ComponentRegistration(name, variables, blocks)
end

component_registration(body::SpatialRigidBodyComponent) =
    spatial_body_registration(body.name)
component_registration(beam::SpatialFlexibleBeamComponent) =
    spatial_flexible_beam_registration(beam.name)
component_registration(::SpatialGravityComponent) = nothing

function allocated_spatial_body(layout, name, mass, inertia)
    variables = component_variable_indices(layout, name)
    state_equations = component_equation_indices(
        layout, name, :selected_state)
    orientation_equations = component_equation_indices(
        layout, name, :orientation)
    SpatialRigidBodyComponent(name, mass, Matrix(inertia),
        variables[1:3], variables[4:6], variables[7:9], variables[10:12],
        variables[13:15], variables[16:18], variables[19:22],
        component_equation_indices(layout, name, :balance),
        state_equations[1:3], state_equations[4:6],
        state_equations[7:9], state_equations[10:12],
        orientation_equations)
end

function allocated_spatial_flexible_beam(layout, name, mass, inertia, length,
        elastic_mass, elastic_stiffness, deformation_shape,
        damping_time_scale)
    variables = component_variable_indices(layout, name)
    states = component_equation_indices(layout, name, :selected_state)
    orientation = component_equation_indices(layout, name, :orientation)
    SpatialFlexibleBeamComponent(name, mass, Matrix(inertia), length,
        Matrix(elastic_mass), Matrix(elastic_stiffness),
        damping_time_scale .* Matrix(elastic_stiffness),
        Matrix(deformation_shape),
        variables[1:3], variables[4:6], variables[7:12],
        variables[13:15], variables[16:18], variables[19:24],
        variables[25:27], variables[28:30], variables[31:34],
        variables[35:40], component_equation_indices(layout, name, :balance),
        states[1:3], states[4:6], states[7:12], states[13:15],
        states[16:18], states[19:24], orientation)
end

function body_balance!(equations, z, body)
    acceleration = @view z[body.acceleration_variables]
    alpha = @view z[body.angular_acceleration_variables]
    omega = @view z[body.angular_velocity_variables]
    equations[body.balance_equations[1:3]] .= body.mass .* acceleration
    equations[body.balance_equations[4:6]] .=
        body.inertia * alpha + cross(omega, body.inertia * omega)
    nothing
end

function body_balance_jacobian!(jacobian, z, body)
    for row in 1:3, column in 1:3
        jacobian[body.balance_equations[row],
            body.acceleration_variables[column]] +=
            row == column ? body.mass : zero(body.mass)
        jacobian[body.balance_equations[row + 3],
            body.angular_acceleration_variables[column]] +=
            body.inertia[row, column]
    end
    omega = @view z[body.angular_velocity_variables]
    derivative = -skew(body.inertia * omega) + skew(omega) * body.inertia
    for row in 1:3, column in 1:3
        jacobian[body.balance_equations[row + 3],
            body.angular_velocity_variables[column]] += derivative[row, column]
    end
    nothing
end

function beam_balance!(equations, z, beam::SpatialFlexibleBeamComponent)
    body_balance!(equations, z, beam)
    rows = beam.balance_equations[7:12]
    equations[rows] .=
        beam.elastic_mass * z[beam.elastic_acceleration_variables] +
        beam.elastic_damping * z[beam.elastic_velocity_variables] +
        beam.elastic_stiffness * z[beam.elastic_position_variables]
    nothing
end

function beam_balance_jacobian!(jacobian, z,
        beam::SpatialFlexibleBeamComponent)
    body_balance_jacobian!(jacobian, z, beam)
    rows = beam.balance_equations[7:12]
    jacobian[rows, beam.elastic_acceleration_variables] .+= beam.elastic_mass
    jacobian[rows, beam.elastic_velocity_variables] .+= beam.elastic_damping
    jacobian[rows, beam.elastic_position_variables] .+= beam.elastic_stiffness
    nothing
end

function body_states!(equations, z, zdot, body)
    equations[body.acceleration_state_equations] .=
        z[body.acceleration_variables] .- zdot[body.velocity_variables]
    equations[body.angular_acceleration_state_equations] .=
        z[body.angular_acceleration_variables] .-
        zdot[body.angular_velocity_variables]
    equations[body.position_state_equations] .=
        z[body.velocity_variables] .- zdot[body.position_variables]
    equations[body.pseudo_angle_state_equations] .=
        z[body.angular_velocity_variables] .-
        zdot[body.pseudo_angle_variables]
    nothing
end

function beam_states!(equations, z, zdot, beam::SpatialFlexibleBeamComponent)
    body_states!(equations, z, zdot, beam)
    equations[beam.elastic_acceleration_state_equations] .=
        z[beam.elastic_acceleration_variables] .-
        zdot[beam.elastic_velocity_variables]
    equations[beam.elastic_position_state_equations] .=
        z[beam.elastic_velocity_variables] .-
        zdot[beam.elastic_position_variables]
    nothing
end

function body_states_jacobian!(jacobian, coefficient, body)
    pairs = (
        (body.acceleration_state_equations, body.acceleration_variables,
            body.velocity_variables),
        (body.angular_acceleration_state_equations,
            body.angular_acceleration_variables,
            body.angular_velocity_variables),
        (body.position_state_equations, body.velocity_variables,
            body.position_variables),
        (body.pseudo_angle_state_equations,
            body.angular_velocity_variables, body.pseudo_angle_variables),
    )
    for (rows, values, derivatives) in pairs, index in 1:3
        jacobian[rows[index], values[index]] += 1
        jacobian[rows[index], derivatives[index]] -= coefficient
    end
    nothing
end

function beam_states_jacobian!(jacobian, coefficient,
        beam::SpatialFlexibleBeamComponent)
    body_states_jacobian!(jacobian, coefficient, beam)
    for index in 1:6
        jacobian[beam.elastic_acceleration_state_equations[index],
            beam.elastic_acceleration_variables[index]] += 1
        jacobian[beam.elastic_acceleration_state_equations[index],
            beam.elastic_velocity_variables[index]] -= coefficient
        jacobian[beam.elastic_position_state_equations[index],
            beam.elastic_velocity_variables[index]] += 1
        jacobian[beam.elastic_position_state_equations[index],
            beam.elastic_position_variables[index]] -= coefficient
    end
    nothing
end

function orientation_equations!(equations, z, zdot, body)
    parameters = @view z[body.euler_parameter_variables]
    parameter_rates = @view zdot[body.euler_parameter_variables]
    rows = body.orientation_equations
    equations[rows[1:3]] .= zdot[body.pseudo_angle_variables] .-
        2 .* transpose(quaternion_rate_matrix(parameters)) * parameter_rates
    equations[rows[4]] = dot(parameters, parameters) - 1
    nothing
end

function orientation_jacobian!(jacobian, z, zdot, coefficient, body)
    parameters = @view z[body.euler_parameter_variables]
    parameter_rates = @view zdot[body.euler_parameter_variables]
    rows = body.orientation_equations
    for index in 1:3
        jacobian[rows[index], body.pseudo_angle_variables[index]] += coefficient
    end
    rate_scalar = parameter_rates[1]
    rate_vector = parameter_rates[2:4]
    parameter_derivative = hcat(-2 .* rate_vector,
        2rate_scalar .* Matrix{eltype(z)}(I, 3, 3) .-
            2 .* skew(rate_vector))
    derivative_derivative =
        -2coefficient .* transpose(quaternion_rate_matrix(parameters))
    for row in 1:3, column in 1:4
        jacobian[rows[row], body.euler_parameter_variables[column]] +=
            parameter_derivative[row, column] +
            derivative_derivative[row, column]
    end
    for column in 1:4
        jacobian[rows[4], body.euler_parameter_variables[column]] +=
            2parameters[column]
    end
    nothing
end

function executable_blocks(body::SpatialRigidBodyComponent)
    ExecutableEquationBlock[
        ExecutableEquationBlock(body.name, :balance,
            collect(body.balance_equations),
            (e, t, z, zd) -> body_balance!(e, z, body),
            (J, t, z, zd, c) -> body_balance_jacobian!(J, z, body)),
        ExecutableEquationBlock(body.name, :selected_state,
            collect(body.acceleration_state_equations.start:
                body.pseudo_angle_state_equations.stop),
            (e, t, z, zd) -> body_states!(e, z, zd, body),
            (J, t, z, zd, c) -> body_states_jacobian!(J, c, body)),
        ExecutableEquationBlock(body.name, :orientation,
            collect(body.orientation_equations),
            (e, t, z, zd) -> orientation_equations!(e, z, zd, body),
            (J, t, z, zd, c) -> orientation_jacobian!(J, z, zd, c, body)),
    ]
end

function executable_blocks(beam::SpatialFlexibleBeamComponent)
    ExecutableEquationBlock[
        ExecutableEquationBlock(beam.name, :balance,
            collect(beam.balance_equations),
            (e, t, z, zd) -> beam_balance!(e, z, beam),
            (J, t, z, zd, c) -> beam_balance_jacobian!(J, z, beam)),
        ExecutableEquationBlock(beam.name, :selected_state,
            collect(beam.acceleration_state_equations.start:
                beam.elastic_position_state_equations.stop),
            (e, t, z, zd) -> beam_states!(e, z, zd, beam),
            (J, t, z, zd, c) -> beam_states_jacobian!(J, c, beam)),
        ExecutableEquationBlock(beam.name, :orientation,
            collect(beam.orientation_equations),
            (e, t, z, zd) -> orientation_equations!(e, z, zd, beam),
            (J, t, z, zd, c) -> orientation_jacobian!(J, z, zd, c, beam)),
    ]
end

function equation_contributions(gravity::SpatialGravityComponent)
    rows = collect(gravity.body.balance_equations[1:3])
    residual! = function (equations, t, z, zdot)
        equations[rows] .-= gravity.body.mass .* gravity.acceleration
    end
    EquationContribution[EquationContribution(gravity.name, :body_force,
        rows, residual!, (J, t, z, zd, c) -> nothing)]
end

executable_blocks(::SpatialGravityComponent) = ExecutableEquationBlock[]
equation_contributions(::SpatialRigidBodyComponent) = EquationContribution[]
equation_contributions(::SpatialFlexibleBeamComponent) = EquationContribution[]

end
