using PracticalMechanicalSimulation
using PracticalMechanicalSimulation.SpatialComponentAssembly
using PracticalMechanicalSimulation.SpatialConstraints
using LinearAlgebra
using Test

const CV_MODEL_PATH = normpath(joinpath(
    @__DIR__, "..", "..", "models", "spatial",
    "constant-velocity-shafts.toml"))

@testset "Spatial constant-velocity phase constraint" begin
    loaded = load_spatial_model(CV_MODEL_PATH)
    joint = loaded.connections[:constant_velocity]
    @test joint isa SpatialConstantVelocityJoint
    constraint = joint.phase
    @test constraint isa SpatialCVPhaseConstraint

    # Equal shaft phases satisfy the constraint while the included shaft angle
    # changes. This is the property that distinguishes the construction from
    # a pair of ordinary projected angular coordinates.
    state = copy(loaded.initial_values)
    input_body = loaded.bodies[:input_shaft]
    output_body = loaded.bodies[:output_shaft]
    for articulation in deg2rad.((10.0, 35.0, 60.0, 100.0)),
            phase in range(-pi, pi; length = 9)
        input_orientation =
            axis_angle_rotation(-articulation / 2, [0.0, 1.0, 0.0]) *
            axis_angle_rotation(phase, [0.0, 0.0, 1.0])
        output_orientation =
            axis_angle_rotation(articulation / 2, [0.0, 1.0, 0.0]) *
            axis_angle_rotation(phase, [0.0, 0.0, 1.0])
        state[input_body.euler_parameter_variables] .=
            matrix_to_euler_parameters(input_orientation)
        state[output_body.euler_parameter_variables] .=
            matrix_to_euler_parameters(output_orientation)
        @test abs(cv_phase_position(constraint, state)) < 2.0e-14
        directions = SpatialConstraints.cv_phase_reaction_directions(
            constraint, state)
        @test norm(directions.first + directions.second) < 3.0e-14
    end

    # The differentiated equations also vanish while articulation and shaft
    # phase change simultaneously.
    for time in range(0.0, 1.7; length = 9)
        articulation = deg2rad(55.0) + 0.25sin(time)
        articulation_rate = 0.25cos(time)
        articulation_acceleration = -0.25sin(time)
        phase = 1.4time + 0.2sin(2time)
        phase_rate = 1.4 + 0.4cos(2time)
        phase_acceleration = -0.8sin(2time)
        for (body, sign) in ((input_body, -1.0), (output_body, 1.0))
            beta = sign * articulation / 2
            beta_rate = sign * articulation_rate / 2
            beta_acceleration = sign * articulation_acceleration / 2
            tilt = axis_angle_rotation(beta, [0.0, 1.0, 0.0])
            orientation = tilt *
                axis_angle_rotation(phase, [0.0, 0.0, 1.0])
            shaft_axis = tilt * [0.0, 0.0, 1.0]
            omega_global = beta_rate .* [0.0, 1.0, 0.0] .+
                phase_rate .* shaft_axis
            alpha_global = beta_acceleration .* [0.0, 1.0, 0.0] .+
                phase_acceleration .* shaft_axis .+
                phase_rate .* beta_rate .* cross(
                    [0.0, 1.0, 0.0], shaft_axis)
            state[body.euler_parameter_variables] .=
                matrix_to_euler_parameters(orientation)
            state[body.angular_velocity_variables] .=
                transpose(orientation) * omega_global
            state[body.angular_acceleration_variables] .=
                transpose(orientation) * alpha_global
        end
        @test abs(cv_phase_position(constraint, state)) < 2.0e-14
        @test abs(cv_phase_velocity(constraint, state)) < 3.0e-14
        @test abs(cv_phase_acceleration(constraint, state)) < 2.0e-13
    end


    # A common rigid-body motion must not appear as relative shaft rotation.
    articulation = deg2rad(70.0)
    phase = 0.37
    common_omega = [0.8, -1.1, 0.6]
    common_alpha = [-0.3, 0.4, 0.9]
    for (body, sign) in ((input_body, -1.0), (output_body, 1.0))
        orientation =
            axis_angle_rotation(sign * articulation / 2, [0.0, 1.0, 0.0]) *
            axis_angle_rotation(phase, [0.0, 0.0, 1.0])
        state[body.euler_parameter_variables] .=
            matrix_to_euler_parameters(orientation)
        state[body.angular_velocity_variables] .=
            transpose(orientation) * common_omega
        state[body.angular_acceleration_variables] .=
            transpose(orientation) * common_alpha
    end
    @test abs(cv_phase_velocity(constraint, state)) < 2.0e-14
    @test abs(cv_phase_acceleration(constraint, state)) < 2.0e-14

    # The constraint equations use the unnormalized shaft-axis sum, while the
    # physical reaction direction remains a unit vector.
    state[input_body.angular_velocity_variables] .= [0.2, -0.1, 1.4]
    state[output_body.angular_velocity_variables] .= [-0.4, 0.3, 0.7]
    state[input_body.angular_acceleration_variables] .= [0.5, 0.2, -0.6]
    state[output_body.angular_acceleration_variables] .= [-0.3, 0.4, 0.1]
    kinematics = SpatialConstraints.cv_phase_kinematics(constraint, state)
    first = SpatialConstraints.marker_angular_kinematics(
        constraint.marker_i, state)
    second = SpatialConstraints.marker_angular_kinematics(
        constraint.marker_j, state)
    @test cv_phase_velocity(constraint, state) ≈
        dot(first.omega - second.omega, kinematics.axis_sum)
    @test cv_phase_acceleration(constraint, state) ≈
        dot(first.alpha - second.alpha, kinematics.axis_sum) +
        dot(first.omega - second.omega, kinematics.axis_sum_velocity)
    directions = SpatialConstraints.cv_phase_reaction_directions(
        constraint, state)
    @test norm(directions.first) ≈ 1.0
    @test directions.second ≈ -directions.first

    state[input_body.euler_parameter_variables] .= [1.0, 0.0, 0.0, 0.0]
    state[output_body.euler_parameter_variables] .=
        matrix_to_euler_parameters(
            axis_angle_rotation(pi, [1.0, 0.0, 0.0]))
    @test_throws ArgumentError cv_phase_position(constraint, state)

    result = run_spatial_model(CV_MODEL_PATH; duration = 0.25, samples = 21)
    input_omega_index = loaded.connections[:input_bearing].hinge.rotation_variables[2]
    output_omega_index = loaded.connections[:output_bearing].rotation_variables[2]
    input_omega = [sample[input_omega_index] for sample in result.states]
    output_omega = [sample[output_omega_index] for sample in result.states]
    @test maximum(abs.(input_omega .- pi)) < 2.0e-7
    @test maximum(abs.(output_omega .- input_omega)) < 2.0e-7
    @test maximum(abs(cv_phase_position(constraint, sample))
        for sample in result.states) < 2.0e-8

    reaction = [sample[constraint.reaction_variable] for sample in result.states]
    @test all(isfinite, reaction)
    @test minimum(abs.(reaction)) > 0.2
end
