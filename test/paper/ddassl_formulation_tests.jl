using Test

include(joinpath(@__DIR__, "..", "..", "src", "common", "HistoricalDDASSL.jl"))
using .HistoricalDDASSL
using SciMLBase

function decay_residual!(out, t, y, yprime, parameter)
    out[1] = yprime[1] + y[1]
end

function decay_jacobian!(matrix, t, y, yprime, cj, parameter)
    matrix[1, 1] = 1 + cj
end

function semi_explicit_residual!(out, t, y, yprime, parameter)
    out[1] = yprime[1] + y[1]
    out[2] = y[2] - y[1]^2
end

@testset "HistoricalDDASSL" begin
    options = DASSLOptions{Float64}(
        rtol = 1.0e-7,
        atol = 1.0e-9,
        initial_step = 1.0e-4,
        maximum_step = 0.05,
    )
    exact = exp(-1.0)

    @test error_control_from_levels([2, 1, 0, 2], 0) ==
        Bool[false, true, true, false]
    @test error_control_from_levels([2, 1, 0, 2], 1) ==
        Bool[false, true, true, false]
    @test error_control_from_levels([2, 1, 0, 2], 2) ==
        Bool[false, false, true, false]
    @test error_control_from_levels([0, 1, 1, 2], 1;
        differential_vars = [true, true, false, false]) ==
        Bool[true, true, false, false]
    @test level_factors([0, 1, 2], 0.1) ≈ [1.0, 0.1, 0.01]
    for step in (1.0, 0.1, 0.01)
        # v - (R-Rold)/h = 0 is level one. Its coefficients become
        # independent of h after row-and-column level scaling.
        raw_matrix = [1.0 -1 / step; 0.0 1.0]
        scaled_matrix = scaled_newton_matrix(raw_matrix,
            level_factors([1, 0], step), level_factors([1, 0], step))
        @test scaled_matrix ≈ [1.0 -1.0; 0.0 1.0]
    end

    numerical = dassl(decay_residual!, [1.0], [-1.0], (0.0, 1.0);
        options)
    @test SciMLBase.successful_retcode(numerical.retcode)
    @test abs(numerical.y[end][1] - exact) < 1.0e-6
    @test maximum(numerical.orders) > 1
    @test numerical.u === numerical.y
    @test numerical.du === numerical.yprime
    @test numerical[1, end] == numerical.y[end][1]
    @test numerical[1, :] == [state[1] for state in numerical.y]
    @test abs(numerical(0.5)[1] - exp(-0.5)) < 2.0e-6
    @test abs(numerical(0.5, Val{1})[1] + exp(-0.5)) < 5.0e-6
    @test numerical.stats.nf == numerical.stats.residual_evaluations
    @test numerical.stats.naccept == numerical.stats.accepted_steps

    analytical = dassl(decay_residual!, [1.0], [-1.0], (0.0, 1.0);
        options, jacobian! = decay_jacobian!)
    @test SciMLBase.successful_retcode(analytical.retcode)
    @test abs(analytical.y[end][1] - exact) < 1.0e-6
    @test analytical.stats.residual_evaluations <
        numerical.stats.residual_evaluations

    dae = dassl(semi_explicit_residual!, [1.0, 1.0], [-1.0, 0.0],
        (0.0, 1.0); options, error_control = Bool[true, false])
    @test SciMLBase.successful_retcode(dae.retcode)
    @test abs(dae.y[end][1] - exact) < 1.0e-6
    @test abs(dae.y[end][2] - dae.y[end][1]^2) < 1.0e-12
    @test dae.error_control == Bool[true, false]

    overridden = dassl(semi_explicit_residual!, [1.0, 1.0], [-1.0, 0.0],
        (0.0, 0.0); options, variable_levels = [0, 2],
        differential_vars = [true, false],
        error_control = [false, true])
    @test overridden.error_control == Bool[false, true]

    leveled_dae = dassl(semi_explicit_residual!, [1.0, 1.0], [-1.0, 0.0],
        (0.0, 1.0); options, variable_levels = [0, 2],
        equation_levels = [1, 0], deficit = 2)
    @test SciMLBase.successful_retcode(leveled_dae.retcode)
    @test leveled_dae.variable_levels == [0, 2]
    @test leveled_dae.equation_levels == [1, 0]
    @test leveled_dae.deficit == 2
    @test abs(leveled_dae.y[end][2] - leveled_dae.y[end][1]^2) < 1.0e-12

    inconsistent = dassl(decay_residual!, [1.0], [0.0], (0.0, 1.0);
        options)
    @test inconsistent.retcode == SciMLBase.ReturnCode.InitialFailure
end

include(joinpath(@__DIR__, "..", "..", "examples", "planar",
    "implicit_acceleration_pendulum.jl"))

@testset "Deficit-zero acceleration pendulum" begin
    solution = run_acceleration_historical_ddassl(tspan = (0.0, 1.0))
    parameters = PendulumParameters()
    diagnostics = acceleration_solution_diagnostics(solution, parameters)
    comparison = compare_implicit_with_reduced(solution, parameters)

    @test SciMLBase.successful_retcode(solution.retcode)
    @test solution.variable_levels == [0, 0, 0, 1, 1, 1, 2, 2]
    @test solution.equation_levels == [2, 2, 2, 1, 1, 1, 2, 2]
    @test solution.differential_vars ==
        Bool[true, true, true, true, true, true, false, false]
    @test solution.error_control ==
        Bool[true, true, true, true, true, true, false, false]
    @test diagnostics.maximum_acceleration_constraint_error < 1.0e-7
    @test comparison.maximum_state_difference < 2.0e-4
end

include(joinpath(@__DIR__, "..", "..", "examples", "planar",
    "implicit_velocity_pendulum.jl"))

@testset "Deficit-one velocity pendulum" begin
    solution = run_velocity_historical_ddassl(tspan = (0.0, 1.0))
    parameters = PendulumParameters()
    diagnostics = velocity_solution_diagnostics(solution, parameters)
    comparison = compare_implicit_with_reduced(solution, parameters)

    @test SciMLBase.successful_retcode(solution.retcode)
    @test solution.variable_levels == [0, 0, 0, 1, 1, 1, 2, 2]
    @test solution.equation_levels == [2, 2, 2, 1, 1, 1, 1, 1]
    @test solution.differential_vars ==
        Bool[true, true, true, true, true, true, false, false]
    @test solution.error_control ==
        Bool[true, true, true, true, true, true, false, false]
    @test diagnostics.maximum_velocity_constraint_error < 1.0e-9
    @test comparison.maximum_state_difference < 3.0e-4
end

include(joinpath(@__DIR__, "..", "..", "examples", "planar",
    "gear_constraint_satisfaction_pendulum.jl"))

@testset "GearStableV constraint-satisfaction pendulum" begin
    p = PendulumParameters()
    z0, zdot0 = gear_constraint_satisfaction_initial_conditions(
        deg2rad(35.0), 0.7, p)
    equations = zeros(13)
    gear_constraint_satisfaction_pendulum!(
        equations, 0.0, z0, zdot0, p)
    @test norm(equations, Inf) < 1.0e-12

    coefficient = 2.3
    jacobian = zeros(13, 13)
    gear_constraint_satisfaction_jacobian!(
        jacobian, 0.0, z0, zdot0, coefficient, p)
    direction = collect(range(-0.3, 0.4; length = 13))
    epsilon = 1.0e-7
    perturbed = zeros(13)
    gear_constraint_satisfaction_pendulum!(perturbed, 0.0,
        z0 + epsilon * direction,
        zdot0 + coefficient * epsilon * direction, p)
    finite_difference = (perturbed - equations) / epsilon
    @test norm(jacobian * direction - finite_difference, Inf) < 1.0e-6

    solution = run_gear_constraint_satisfaction_pendulum(
        tspan = (0.0, 1.0))
    diagnostics = gear_constraint_satisfaction_diagnostics(solution, p)
    @test SciMLBase.successful_retcode(solution.retcode)
    @test solution.variable_levels == GEAR_VARIABLE_LEVELS
    @test solution.equation_levels == GEAR_EQUATION_LEVELS
    @test solution.differential_vars == GEAR_DIFFERENTIAL_VARS
    @test solution.error_control == GEAR_DIFFERENTIAL_VARS
    @test diagnostics.maximum_position_constraint_error < 1.0e-11
    @test diagnostics.maximum_velocity_constraint_error < 1.0e-9
    @test diagnostics.maximum_implicit_equation_error < 1.0e-7
    @test diagnostics.maximum_state_difference < 5.0e-4
    @test diagnostics.maximum_satisfaction_multiplier < 2.0e-4
end

@testset "GearStableA constraint-satisfaction pendulum" begin
    p = PendulumParameters()
    z0, zdot0 = complete_gear_constraint_satisfaction_initial_conditions(
        deg2rad(35.0), 0.7, p)
    equations = zeros(15)
    complete_gear_constraint_satisfaction_pendulum!(
        equations, 0.0, z0, zdot0, p)
    @test norm(equations, Inf) < 1.0e-12

    coefficient = 2.3
    jacobian = zeros(15, 15)
    complete_gear_constraint_satisfaction_jacobian!(
        jacobian, 0.0, z0, zdot0, coefficient, p)
    direction = collect(range(-0.3, 0.4; length = 15))
    epsilon = 1.0e-7
    perturbed = zeros(15)
    complete_gear_constraint_satisfaction_pendulum!(perturbed, 0.0,
        z0 + epsilon * direction,
        zdot0 + coefficient * epsilon * direction, p)
    finite_difference = (perturbed - equations) / epsilon
    @test norm(jacobian * direction - finite_difference, Inf) < 1.0e-6

    solution = run_complete_gear_constraint_satisfaction_pendulum(
        tspan = (0.0, 1.0))
    diagnostics = complete_gear_constraint_satisfaction_diagnostics(
        solution, p)
    @test SciMLBase.successful_retcode(solution.retcode)
    @test solution.variable_levels == COMPLETE_GEAR_VARIABLE_LEVELS
    @test solution.equation_levels == COMPLETE_GEAR_EQUATION_LEVELS
    @test solution.differential_vars == COMPLETE_GEAR_DIFFERENTIAL_VARS
    @test solution.error_control == COMPLETE_GEAR_DIFFERENTIAL_VARS
    @test diagnostics.maximum_position_constraint_error < 1.0e-11
    @test diagnostics.maximum_velocity_constraint_error < 1.0e-9
    @test diagnostics.maximum_acceleration_constraint_error < 1.0e-7
    @test diagnostics.maximum_implicit_equation_error < 1.0e-7
    @test diagnostics.maximum_state_difference < 5.0e-4
    @test diagnostics.maximum_acceleration_difference < 2.0e-3
    @test diagnostics.maximum_acceleration_multiplier < 5.0e-3
end

@testset "Deficit-two displacement pendulum" begin
    solution = run_displacement_historical_ddassl(tspan = (0.0, 1.0))
    parameters = PendulumParameters()
    diagnostics = implicit_solution_diagnostics(solution, parameters)
    comparison = compare_implicit_with_reduced(solution, parameters)

    @test SciMLBase.successful_retcode(solution.retcode)
    @test solution.variable_levels == [0, 0, 0, 1, 1, 1, 2, 2]
    @test solution.equation_levels == [2, 2, 2, 1, 1, 1, 0, 0]
    @test solution.differential_vars ==
        Bool[true, true, true, true, true, true, false, false]
    @test solution.error_control ==
        Bool[true, true, true, false, false, false, false, false]
    @test diagnostics.maximum_position_constraint_error < 1.0e-10
    @test comparison.maximum_state_difference < 1.0e-3
end

include(joinpath(@__DIR__, "..", "..", "examples", "planar",
    "implicit_baumgarte_pendulum.jl"))

@testset "Baumgarte-stabilized pendulum" begin
    correction_time = 0.1
    parameters = BaumgarteParameters(PendulumParameters(), correction_time)
    solution = run_baumgarte_historical_ddassl(
        tspan = (0.0, 1.0); correction_time)
    diagnostics = baumgarte_diagnostics(solution, parameters)

    @test SciMLBase.successful_retcode(solution.retcode)
    @test solution.variable_levels == [0, 0, 0, 1, 1, 1, 2, 2]
    @test solution.equation_levels == [2, 2, 2, 1, 1, 1, 2, 2]
    @test solution.error_control ==
        Bool[true, true, true, true, true, true, false, false]
    @test diagnostics.final_position_error <
        diagnostics.initial_position_error / 10
    @test diagnostics.maximum_acceleration_error ≈ 0.1 rtol = 1.0e-10
    @test diagnostics.maximum_critical_decay_difference < 1.0e-4
    @test diagnostics.maximum_implicit_equation_error < 1.0e-6
end

@testset "First-order stabilized pendulum" begin
    correction_time = 0.1
    parameters = BaumgarteParameters(PendulumParameters(), correction_time)
    solution = run_first_order_stabilized_historical_ddassl(
        tspan = (0.0, 1.0); correction_time)
    diagnostics = first_order_stabilization_diagnostics(solution, parameters)

    @test SciMLBase.successful_retcode(solution.retcode)
    @test solution.variable_levels == [0, 0, 0, 1, 1, 1, 2, 2]
    @test solution.equation_levels == [2, 2, 2, 1, 1, 1, 1, 1]
    @test solution.error_control ==
        Bool[true, true, true, true, true, true, false, false]
    @test diagnostics.final_position_error <
        diagnostics.initial_position_error / 10
    @test diagnostics.maximum_exponential_decay_difference < 1.0e-4
    @test diagnostics.maximum_implicit_equation_error < 1.0e-6
end

include(joinpath(@__DIR__, "..", "..", "examples", "planar",
    "independent_state_pendulum.jl"))

@testset "Independent-state implicit pendulum" begin
    solution = run_independent_state_pendulum(tspan = (0.0, 1.0))
    diagnostics = independent_state_diagnostics(solution, PendulumParameters())

    @test SciMLBase.successful_retcode(solution.retcode)
    @test solution.variable_levels ==
        [2, 2, 2, 1, 1, 1, 0, 0, 0, 2, 2]
    @test solution.equation_levels ==
        [2, 2, 2, 2, 2, 1, 1, 0, 0, 2, 1]
    @test solution.differential_vars ==
        Bool[false, false, false, false, false, true,
             false, false, true, false, false]
    @test solution.error_control == solution.differential_vars
    @test diagnostics.maximum_implicit_equation_error < 1.0e-7
    @test diagnostics.maximum_integrated_state_difference < 2.0e-4
    @test diagnostics.maximum_acceleration_difference < 5.0e-4
    @test diagnostics.maximum_reaction_difference < 5.0e-4
end

include(joinpath(@__DIR__, "..", "..", "examples", "planar",
    "torsional_spring_damper_pendulum.jl"))

@testset "Planar torsional spring-damper" begin
    system = TorsionalSpringPendulumSystem(
        stiffness = 3.0, damping = 0.2, free_angle = 0.1)
    z0, zdot0 = torsional_spring_initial_conditions(
        deg2rad(35.0), 0.7, system)
    relative = relative_rotation(system.spring_damper, z0)
    @test relative.angle ≈ deg2rad(35.0)
    @test relative.angular_velocity ≈ 0.7
    @test z0[12] ≈ -3.0 * (deg2rad(35.0) - 0.1) - 0.2 * 0.7

    equations = zeros(12)
    torsional_spring_damper_pendulum!(
        equations, 0.0, z0, zdot0, system)
    @test norm(equations, Inf) < 1.0e-12

    coefficient = 2.3
    jacobian = zeros(12, 12)
    torsional_spring_damper_pendulum_jacobian!(
        jacobian, 0.0, z0, zdot0, coefficient, system)
    direction = collect(range(-0.3, 0.4; length = 12))
    epsilon = 1.0e-7
    perturbed = zeros(12)
    torsional_spring_damper_pendulum!(perturbed, 0.0,
        z0 + epsilon * direction,
        zdot0 + coefficient * epsilon * direction, system)
    finite_difference = (perturbed - equations) / epsilon
    @test norm(jacobian * direction - finite_difference, Inf) < 1.0e-6

    solution, run_system = run_torsional_spring_damper_pendulum(
        stiffness = 3.0, damping = 0.2, free_angle = 0.1,
        tspan = (0.0, 1.0))
    diagnostics = torsional_spring_damper_diagnostics(solution, run_system)
    @test SciMLBase.successful_retcode(solution.retcode)
    @test solution.variable_levels == TORSIONAL_SPRING_VARIABLE_LEVELS
    @test solution.equation_levels == TORSIONAL_SPRING_EQUATION_LEVELS
    @test solution.differential_vars == TORSIONAL_SPRING_DIFFERENTIAL_VARS
    @test solution.error_control == TORSIONAL_SPRING_DIFFERENTIAL_VARS
    @test diagnostics.maximum_implicit_equation_error < 1.0e-7
    @test diagnostics.maximum_constitutive_torque_error < 1.0e-12
    @test diagnostics.maximum_selected_state_difference < 1.0e-3
    @test diagnostics.total_energy_change < 0
end

include(joinpath(@__DIR__, "..", "..", "examples", "planar",
    "spanning_spring_damper_pendulum.jl"))

@testset "Planar spanning spring-damper" begin
    system = SpanningSpringPendulumSystem()
    z0, zdot0 = spanning_spring_initial_conditions(
        deg2rad(35.0), 0.7, system)
    local_values = spanning_force_kinematics(system.spring_damper, z0)
    @test local_values.spanning ≈
        local_values.marker_2.position - local_values.marker_1.position
    @test local_values.length ≈ norm(local_values.spanning)
    @test norm(local_values.unit) ≈ 1.0
    @test local_values.length_rate ≈ dot(local_values.unit,
        local_values.marker_2.velocity - local_values.marker_1.velocity)
    @test local_values.global_force ≈
        -local_values.unit .* local_values.scalar_force

    equations = zeros(20)
    spanning_spring_damper_pendulum!(
        equations, 0.0, z0, zdot0, system)
    @test norm(equations, Inf) < 1.0e-12

    coefficient = 2.3
    jacobian = zeros(20, 20)
    spanning_spring_damper_pendulum_jacobian!(
        jacobian, 0.0, z0, zdot0, coefficient, system)
    direction = collect(range(-0.3, 0.4; length = 20))
    epsilon = 1.0e-7
    perturbed = zeros(20)
    spanning_spring_damper_pendulum!(perturbed, 0.0,
        z0 + epsilon * direction,
        zdot0 + coefficient * epsilon * direction, system)
    finite_difference = (perturbed - equations) / epsilon
    @test norm(jacobian * direction - finite_difference, Inf) < 1.0e-6

    solution, run_system = run_spanning_spring_damper_pendulum(
        tspan = (0.0, 1.0))
    diagnostics = spanning_spring_diagnostics(solution, run_system)
    @test SciMLBase.successful_retcode(solution.retcode)
    @test solution.variable_levels == SPANNING_SPRING_VARIABLE_LEVELS
    @test solution.equation_levels == SPANNING_SPRING_EQUATION_LEVELS
    @test solution.differential_vars == SPANNING_SPRING_DIFFERENTIAL_VARS
    @test solution.error_control == SPANNING_SPRING_DIFFERENTIAL_VARS
    @test diagnostics.maximum_implicit_equation_error < 1.0e-7
    @test diagnostics.maximum_selected_state_difference < 1.0e-3
    @test diagnostics.minimum_length > 0
    @test diagnostics.total_energy_change < 0
end

include(joinpath(@__DIR__, "..", "..", "examples", "planar",
    "spanning_spring_static_equilibrium.jl"))

@testset "Spanning spring static equilibrium" begin
    system = SpanningSpringPendulumSystem()
    initial = spanning_spring_static_initial_guess(deg2rad(35.0), system)
    equations = zeros(13)
    spanning_spring_static_equations!(equations, initial, system)
    jacobian = spanning_spring_static_jacobian(initial, system)
    @test size(jacobian) == (13, 13)

    direction = collect(range(-0.3, 0.4; length = 13))
    epsilon = 1.0e-7
    perturbed = zeros(13)
    spanning_spring_static_equations!(
        perturbed, initial + epsilon * direction, system)
    finite_difference = (perturbed - equations) / epsilon
    @test norm(jacobian * direction - finite_difference, Inf) < 1.0e-6

    y, solved_system, iterations =
        solve_spanning_spring_static_equilibrium()
    diagnostics = spanning_spring_static_diagnostics(
        y, solved_system, iterations)
    @test iterations <= 5
    @test diagnostics.maximum_equation_error < 1.0e-12
    @test diagnostics.equilibrium_angle_degrees ≈ 34.329079079 atol=1.0e-8
    @test diagnostics.reduced_moment_error < 1.0e-12
    @test diagnostics.spring_length > solved_system.spring_damper.free_length
    @test expand_static_spanning_state(y)[17] == 0

    y_other_damping, _, _ = solve_spanning_spring_static_equilibrium(
        damping = 50.0)
    @test y_other_damping ≈ y atol=1.0e-12
end

include(joinpath(@__DIR__, "..", "..", "examples", "planar",
    "automatic_spanning_pendulum_analyses.jl"))

@testset "Automatic spanning-pendulum analyses" begin
    catalog = SPANNING_ANALYSIS_CATALOG
    @test length(catalog.variables) == 20
    @test length(catalog.equations) == 20
    body_registration = planar_body_analysis_registration()
    pin_registration = revolute_pin_analysis_registration()
    spring_registration = spanning_spring_analysis_registration()
    builder = AnalysisCatalogBuilder()
    @test register_component_variables!(builder, body_registration) == 1:9
    @test register_component_variables!(builder, pin_registration) == 10:11
    @test register_component_variables!(builder, spring_registration) == 12:20
    @test register_component_equation_block!(
        builder, body_registration, :balance) == 1:3
    @test register_component_equation_block!(
        builder, pin_registration, :acceleration) == 4:5
    @test register_component_equation_block!(
        builder, pin_registration, :velocity) == 6:7
    @test register_component_equation_block!(
        builder, pin_registration, :position) == 8:9
    @test register_component_equation_block!(
        builder, body_registration, :selected_state) == 10:11
    @test register_component_equation_block!(
        builder, spring_registration, :geometry) == 12:16
    @test register_component_equation_block!(
        builder, spring_registration, :rate) == 17:17
    @test register_component_equation_block!(
        builder, spring_registration, :load) == 18:20
    rebuilt_catalog = finish_catalog(builder)
    @test getfield.(rebuilt_catalog.variables, :name) ==
        getfield.(catalog.variables, :name)
    @test_throws ErrorException register_component_variables!(
        builder, body_registration)
    selections = automatic_analysis_selections()
    @test selections.position.variable_indices == [7, 8, 9]
    @test selections.position.equation_indices == [8, 9]
    @test selections.velocity.variable_indices == [4, 5, 6]
    @test selections.velocity.equation_indices == [6, 7]
    @test selections.acceleration.variable_indices ==
        [1, 2, 3, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    @test selections.acceleration.equation_indices ==
        [1, 2, 3, 4, 5, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    @test selections.static.variable_indices ==
        [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20]
    @test selections.static.equation_indices ==
        [1, 2, 3, 8, 9, 12, 13, 14, 15, 16, 18, 19, 20]
    @test selections.dynamics.variable_indices == collect(1:20)
    @test selections.dynamics.equation_indices == collect(1:20)

    executable_system = SpanningSpringPendulumSystem()
    executable_model = spanning_pendulum_executable_model(executable_system)
    @test length(executable_model.equation_blocks) == 8
    @test length(executable_model.contributions) == 3
    zero_length_context = zeros(20)
    zero_length_context[7:9] .= [0.0, 0.0, deg2rad(20.0)]
    position_equations = zeros(2)
    evaluate_analysis_equations!(position_equations, executable_model,
        selections.position, 0.0, zero_length_context, zeros(20))
    @test position_equations ≈ position_constraint(
        zero_length_context[7:9], executable_system.mechanical)
    velocity_equations = zeros(2)
    evaluate_analysis_equations!(velocity_equations, executable_model,
        selections.velocity, 0.0, zero_length_context, zeros(20))
    @test velocity_equations == zeros(2)

    full_state, full_derivative = spanning_spring_initial_conditions(
        deg2rad(35.0), 0.7, executable_system)
    executable_equations = zeros(20)
    evaluate_analysis_equations!(executable_equations, executable_model,
        selections.dynamics, 0.0, full_state, full_derivative)
    monolithic_equations = zeros(20)
    spanning_spring_damper_pendulum!(monolithic_equations, 0.0,
        full_state, full_derivative, executable_system)
    @test norm(executable_equations - monolithic_equations, Inf) < 1.0e-12
    executable_jacobian = evaluate_analysis_jacobian(executable_model,
        selections.dynamics, 0.0, full_state, full_derivative, 2.3)
    monolithic_jacobian = zeros(20, 20)
    spanning_spring_damper_pendulum_jacobian!(monolithic_jacobian, 0.0,
        full_state, full_derivative, 2.3, executable_system)
    @test executable_jacobian ≈ monolithic_jacobian

    system = SpanningSpringPendulumSystem()
    position_state, _ = automatic_position_initialization(
        [0.03, -0.02, deg2rad(35.0)], system)
    @test norm(position_constraint(
        position_state[7:9], system.mechanical), Inf) < 1.0e-12
    velocity_state, _ = automatic_velocity_initialization(
        [0.1, -0.2, 0.7], position_state, system)
    @test norm(velocity_constraint(
        [velocity_state[7:9]; velocity_state[4:6]],
        system.mechanical), Inf) < 1.0e-12

    acceleration_state, acceleration_selection, iterations =
        automatic_acceleration_initialization(velocity_state, system)
    acceleration_equations = zeros(14)
    selected_spanning_equations!(acceleration_equations,
        analysis_values(acceleration_state, acceleration_selection),
        velocity_state, acceleration_selection, system)
    @test iterations <= 2
    @test norm(acceleration_equations, Inf) < 1.0e-12
end

include(joinpath(@__DIR__, "..", "..", "examples", "planar",
    "euler_parameter_independent_state_pendulum.jl"))

@testset "Euler-parameter orientation pendulum" begin
    mechanical = PendulumParameters()
    angle_orientation = run_euler_parameter_independent_state_pendulum(
        tspan = (0.0, 1.0),
        use_euler_parameters_for_orientation = false,
        control_euler_parameters = true)
    parameter_orientation = run_euler_parameter_independent_state_pendulum(
        tspan = (0.0, 1.0),
        use_euler_parameters_for_orientation = true,
        control_euler_parameters = false)
    parameter_controlled = run_euler_parameter_independent_state_pendulum(
        tspan = (0.0, 1.0),
        use_euler_parameters_for_orientation = true,
        control_euler_parameters = true,
        control_theta = false)
    angle_diagnostics = euler_parameter_pendulum_diagnostics(
        angle_orientation,
        EulerParameterPendulumParameters(mechanical, false))
    parameter_diagnostics = euler_parameter_pendulum_diagnostics(
        parameter_orientation,
        EulerParameterPendulumParameters(mechanical, true))
    parameter_controlled_diagnostics = euler_parameter_pendulum_diagnostics(
        parameter_controlled,
        EulerParameterPendulumParameters(mechanical, true))

    @test EULER_PARAMETER_DIFFERENTIAL_VARS ==
        Bool[false, false, false, false, false, true,
             false, false, true, false, false, true, true]
    @test angle_orientation.error_control ==
        EULER_PARAMETER_DIFFERENTIAL_VARS
    @test parameter_orientation.error_control ==
        EULER_PARAMETER_STATE_ERROR_CONTROL
    @test parameter_orientation.equation_levels ==
        [2, 2, 2, 2, 2, 1, 1, 0, 0, 2, 1, 1, 0]

    for (solution, diagnostics) in (
            (angle_orientation, angle_diagnostics),
            (parameter_orientation, parameter_diagnostics))
        @test SciMLBase.successful_retcode(solution.retcode)
        @test diagnostics.maximum_euler_normalization_error < 1.0e-10
        @test diagnostics.maximum_orientation_matrix_difference < 3.0e-5
        @test diagnostics.maximum_mechanical_orientation_difference < 2.0e-4
        @test diagnostics.maximum_integrated_state_difference < 3.0e-4
        @test diagnostics.maximum_acceleration_difference < 6.0e-4
    end

    @test parameter_controlled.error_control ==
        Bool[false, false, false, false, false, true,
             false, false, false, false, false, true, true]
    @test SciMLBase.successful_retcode(parameter_controlled.retcode)
    @test parameter_controlled_diagnostics.maximum_euler_normalization_error <
        1.0e-10
    @test parameter_controlled_diagnostics.maximum_orientation_matrix_difference <
        5.0e-5
    @test parameter_controlled_diagnostics.maximum_mechanical_orientation_difference <
        3.0e-4
    @test parameter_controlled_diagnostics.maximum_acceleration_difference <
        4.0e-3
end

include(joinpath(@__DIR__, "..", "..", "examples", "planar",
    "modular_sparse_two_body_pendulum.jl"))
include(joinpath(@__DIR__, "..", "..", "examples", "planar",
    "automatic_two_body_component_assembly.jl"))

include(joinpath(@__DIR__, "..", "..", "examples", "planar",
    "rotational_motion_generator_pendulum.jl"))
include(joinpath(@__DIR__, "..", "..", "examples", "planar",
    "driven_base_double_pendulum.jl"))

@testset "Rotational motion-generator pendulum" begin
    time = 0.7
    result = driven_pendulum_diagnostics(time)
    state = result.state
    assembly = result.system.model
    body = assembly.body
    motion = assembly.motion
    @test length(state) == 15
    @test length(assembly.layout.catalog.equations) == 15
    @test state[body.orientation_variable] ≈ motion.motion(time) atol=1.0e-13
    @test state[body.angular_velocity_variable] ≈
        motion.motion_derivative(time) atol=1.0e-13
    @test state[body.angular_acceleration_variable] ≈
        motion.motion_second_derivative(time) atol=1.0e-13
    marker = PlanarAppliedForces.point_marker_kinematics(
        assembly.pin.marker_a, state)
    @test marker.position ≈ assembly.pin.marker_b.position atol=1.0e-13
    @test marker.velocity ≈ zeros(2) atol=1.0e-13
    expected_acceleration = state[body.acceleration_variables] .+
        marker.d .* state[body.angular_acceleration_variable] .-
        marker.r .* state[body.angular_velocity_variable]^2
    @test expected_acceleration ≈ zeros(2) atol=1.0e-13
    @test maximum(values(result.errors)) < 1.0e-12
    position = select_analysis(assembly.layout.catalog, KinematicPosition())
    velocity = select_analysis(assembly.layout.catalog, KinematicVelocity())
    acceleration = select_analysis(
        assembly.layout.catalog, KinematicAcceleration())
    forces = select_analysis(assembly.layout.catalog, KinematicForces())
    @test length(position.variable_indices) == length(position.equation_indices) == 4
    @test length(velocity.variable_indices) == length(velocity.equation_indices) == 4
    @test length(acceleration.variable_indices) ==
        length(acceleration.equation_indices) == 4
    @test length(forces.variable_indices) == length(forces.equation_indices) == 3
    @test isfinite(state[motion.torque_variable])
    zero_gravity = driven_pendulum_diagnostics(time;
        mechanical = PendulumParameters(gravity = [0.0, 0.0]))
    kinematic_variables = [collect(1:9); collect(12:14)]
    @test zero_gravity.state[kinematic_variables] ≈
        state[kinematic_variables] atol=1.0e-13
    @test zero_gravity.state[[10, 11, 15]] != state[[10, 11, 15]]
end

@testset "Driven-base double pendulum" begin
    state_selection = driven_base_state_selection(
        deg2rad(20.0), deg2rad(-20.0))
    @test size(state_selection.D) == (5, 6)
    @test state_selection.qr_selection.rank == 5
    @test state_selection.preferred.independent == [6]
    @test norm(state_selection.D * state_selection.preferred.P, Inf) < 1.0e-12
    @test isfinite(state_selection.preferred.dependent_condition)

    solution, _, assembly = run_driven_base_double_pendulum(
        tspan = (0.0, 0.5))
    diagnostics = driven_base_double_pendulum_diagnostics(solution, assembly)
    @test SciMLBase.successful_retcode(solution.retcode)
    @test length(first(solution.u)) == 26
    @test count(solution.differential_vars) == 2
    @test findall(solution.differential_vars) ==
        [assembly.body_2.angular_velocity_variable,
         assembly.body_2.orientation_variable]
    @test solution.error_control == solution.differential_vars
    @test diagnostics.maximum_implicit_equation_error < 1.0e-10
    @test diagnostics.maximum_position_constraint_error < 1.0e-12
    @test diagnostics.maximum_prescribed_angle_error < 1.0e-12
    @test diagnostics.maximum_order >= 4
end

@testset "Modular sparse two-body pendulum" begin
    system = ModularTwoBodySystem()
    z0, zdot0 = modular_two_body_initial_conditions(
        deg2rad(35.0), 0.2, deg2rad(-20.0), -0.1, system)
    equations = zeros(22)
    modular_two_body_equations!(equations, 0.0, z0, zdot0, system)
    @test norm(equations, Inf) < 1.0e-12

    coefficient = 2.3
    jacobian = modular_two_body_sparse_jacobian(
        0.0, z0, zdot0, coefficient, system)
    @test jacobian isa SparseMatrixCSC
    @test size(jacobian) == (22, 22)
    @test nnz(jacobian) == 82
    @test nnz(jacobian) < length(jacobian) / 5

    direction = collect(range(-0.3, 0.4; length = 22))
    epsilon = 1.0e-7
    perturbed = zeros(22)
    modular_two_body_equations!(perturbed, 0.0,
        z0 + epsilon * direction,
        zdot0 + coefficient * epsilon * direction, system)
    finite_difference = (perturbed - equations) / epsilon
    @test norm(jacobian * direction - finite_difference, Inf) < 1.0e-6

    component_model = automatic_two_body_executable_model(system)
    allocated_assembly = automatic_two_body_assembly(system)
    allocated_layout = allocated_assembly.layout
    @test component_variable_indices(allocated_layout, :body_1) == 1:9
    @test component_variable_indices(allocated_layout, :body_2) == 10:18
    @test component_variable_indices(allocated_layout, :pin_1) == 19:20
    @test component_variable_indices(allocated_layout, :pin_2) == 21:22
    @test component_equation_indices(
        allocated_layout, :body_1, :balance) == 1:3
    @test component_equation_indices(
        allocated_layout, :pin_2, :position) == 17:18
    @test component_equation_indices(
        allocated_layout, :body_2, :selected_state) == 21:22
    @test length(component_model.equation_blocks) == 10
    @test length(component_model.contributions) == 4
    @test length(component_model.catalog.variables) == 22
    @test length(component_model.catalog.equations) == 22
    position_selection = select_analysis(component_model.catalog, PositionIC())
    velocity_selection = select_analysis(component_model.catalog, VelocityIC())
    acceleration_selection = select_analysis(
        component_model.catalog, AccelerationIC())
    @test position_selection.variable_indices == [7, 8, 9, 16, 17, 18]
    @test position_selection.equation_indices == [11, 12, 17, 18]
    @test velocity_selection.variable_indices == [4, 5, 6, 13, 14, 15]
    @test velocity_selection.equation_indices == [9, 10, 15, 16]
    @test acceleration_selection.variable_indices ==
        [1, 2, 3, 10, 11, 12, 19, 20, 21, 22]
    @test acceleration_selection.equation_indices ==
        [1, 2, 3, 4, 5, 6, 7, 8, 13, 14]
    component_equations = zeros(22)
    automatic_two_body_equations!(component_equations, 0.0, z0, zdot0,
        component_model)
    @test component_equations ≈ equations atol=1.0e-12
    component_jacobian = automatic_two_body_jacobian(
        0.0, z0, zdot0, coefficient, component_model)
    @test component_jacobian isa SparseMatrixCSC
    @test nnz(component_jacobian) == nnz(jacobian)
    @test component_jacobian ≈ Matrix(jacobian) atol=1.0e-12
    fixed_component_jacobian = copy(component_jacobian)
    fixed_columns = copy(fixed_component_jacobian.colptr)
    fixed_rows = copy(fixed_component_jacobian.rowval)
    automatic_two_body_sparse_jacobian!(fixed_component_jacobian,
        0.0, z0, zdot0, coefficient, component_model)
    @test fixed_component_jacobian == component_jacobian
    @test fixed_component_jacobian.colptr == fixed_columns
    @test fixed_component_jacobian.rowval == fixed_rows

    solution, run_system = run_modular_sparse_two_body_pendulum(
        tspan = (0.0, 0.25))
    diagnostics = modular_two_body_diagnostics(solution, run_system)
    @test SciMLBase.successful_retcode(solution.retcode)
    @test solution.differential_vars == MODULAR_TWO_BODY_DIFFERENTIAL_VARS
    @test solution.error_control == MODULAR_TWO_BODY_DIFFERENTIAL_VARS
    @test diagnostics.maximum_implicit_equation_error < 1.0e-7
    @test diagnostics.maximum_constraint_equation_error < 1.0e-7
    @test diagnostics.maximum_selected_state_difference < 5.0e-4
    @test diagnostics.maximum_reference_energy_drift < 1.0e-8
    @test diagnostics.maximum_bdf_energy_drift < 2.0e-3
    @test maximum(solution.orders) == 5
    @test last(solution.orders) >= 4

    component_solution, component_system, _ =
        run_automatic_two_body_component_model(tspan = (0.0, 0.25))
    component_diagnostics = modular_two_body_diagnostics(
        component_solution, component_system)
    @test SciMLBase.successful_retcode(component_solution.retcode)
    @test component_diagnostics.jacobian_stored_entries == 82
    @test component_diagnostics.maximum_implicit_equation_error < 1.0e-7
    @test last(component_solution.u) ≈ last(solution.u) atol=1.0e-11
end
