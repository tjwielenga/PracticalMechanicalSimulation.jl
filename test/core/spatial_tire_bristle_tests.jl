using Test
using LinearAlgebra
using PracticalMechanicalSimulation
using PracticalMechanicalSimulation.SpatialTires
using PracticalMechanicalSimulation.SpatialSimulationRunner
using PracticalMechanicalSimulation.AutomaticAnalysis

@testset "Spatial rolling-tire bristle law" begin
    source = read(joinpath(@__DIR__, "..", "..", "models", "spatial",
        "driven-rolling-tire.toml"), String)
    source = replace(source,
        "longitudinal_expression = \"tire.normal_force * longitudinal_stiffness_per_load * tire.slip_ratio\"" =>
            "tangential_model = \"bristle\"\n" *
            "patch_length_by_load = [[0, 0], [3000, 0.22], [6000, 0.29]]\n" *
            "cornering_stiffness_by_load = [[0, 0], [3000, 45000], [6000, 70000]]\n" *
            "longitudinal_slip_stiffness_by_load = [[0, 0], [3000, 60000], [6000, 95000]]\n" *
            "shear_release_time = 0.01",
        "lateral_expression = \"-tire.normal_force * lateral_stiffness_per_load * tire.slip_angle\"" => "")
    loaded = load_spatial_model(IOBuffer(source))
    tire = loaded.forces[:tire]
    @test tire.bristle && tire.transient
    @test SpatialTires.tire_curve(tire.lateral_curve, 4500.0) == 57500.0
    @test SpatialTires.tire_curve(tire.lateral_curve, 7500.0) == 82500.0
    @test SpatialTires.tire_bristle_properties(tire, 3000.0).stiffness_y ≈
        45000 / 0.22

    state = copy(loaded.initial_values)
    fx = tire.longitudinal_deformation_variable
    fy = tire.lateral_deformation_variable
    n = tire.normal_force_variable
    sx = tire.longitudinal_slip_velocity_variable
    sy = tire.lateral_slip_velocity_variable
    vx = tire.forward_velocity_variable
    state[n] = 3000.0
    state[vx] = 0.0
    state[sx] = state[sy] = 0.0
    state[fx] = 0.004
    state[fy] = 0.006
    @test collect(tire_deformation_rates(tire, state)) ≈ zeros(2) atol=1e-12
    trial = SpatialTires.tire_bristle_trial(tire, state, state[n])
    @test trial[1] < 0 && trial[2] < 0
    @test trial[2] ≈ -45000 / 0.22 * state[fy]

    # In the small-slip limit the steady-state slope is K*L=C, independent
    # of rolling speed; no slip-angle regularization enters these rates.
    state[fx] = state[fy] = 0.0
    for speed in (0.05, 1.0, 20.0)
        state[vx] = speed
        state[sy] = 1e-5 * speed
        length_y = SpatialTires.tire_bristle_properties(tire, state[n]).length_y
        state[fy] = 1e-5 * length_y
        @test tire_deformation_rates(tire, state)[2] ≈ 0.0 atol=1e-7
        @test SpatialTires.tire_bristle_trial(tire, state, state[n])[2] ≈
            -0.45 atol=1e-10
        state[fy] = 0.0
    end

    # The bristles remain elastic through most of the friction ellipse. At its
    # boundary, plastic release cancels only slip that would drive the elastic
    # force farther outward; reverse slip unloads it without plastic release.
    state[n] = 3000.0
    state[vx] = 10.0
    state[sx] = state[fx] = 0.0
    state[tire.deflection_rate_variable] = 0.0
    properties = SpatialTires.tire_bristle_properties(tire, state[n])
    lateral_limit = tire.mu_lateral * state[n]
    boundary_shear = lateral_limit / properties.stiffness_y
    threshold_slip = state[vx] / properties.length_y * boundary_shear
    state[fy] = 0.5 * boundary_shear
    state[sy] = 2 * threshold_slip
    @test tire_deformation_rates(tire, state)[2] > 0
    state[fy] = boundary_shear
    @test tire_deformation_rates(tire, state)[2] ≈ 0.0 atol=1e-10
    state[sy] = -threshold_slip
    @test tire_deformation_rates(tire, state)[2] < 0
    state[fy] = 1.1 * boundary_shear
    state[sy] = 0.0
    @test tire_deformation_rates(tire, state)[2] < 0

    state[sy] = 0.0
    state[vx] = 0.0
    state[fy] = 0.005
    state[tire.deflection_rate_variable] = -0.05
    @test tire_deformation_rates(tire, state)[2] < 0
    state[tire.deflection_rate_variable] = 0.05
    @test tire_deformation_rates(tire, state)[2] ≈ 0.0 atol=1e-12
    # A shrinking patch removes stored shear; growing it cannot restore
    # lost deformation without new slip.
    state[n] = 3000.0
    state[tire.deflection_rate_variable] = -0.05
    original_trial = abs(SpatialTires.tire_bristle_trial(
        tire, state, state[n])[2])
    for _ in 1:10
        state[fy] += 0.01 * tire_deformation_rates(tire, state)[2]
        state[n] -= 20.0
    end
    state[tire.deflection_rate_variable] = 0.05
    state[n] = 3000.0
    @test tire_deformation_rates(tire, state)[2] ≈ 0.0 atol=1e-12
    @test abs(SpatialTires.tire_bristle_trial(tire, state, state[n])[2]) <
        original_trial
    state[n] = 0.0
    @test SpatialTires.tire_bristle_trial(tire, state, state[n]) == (0.0, 0.0)
    @test tire_deformation_rates(tire, state)[2] ≈ -state[fy] / 0.01
    state[n] = 1e-9
    state[sy] = 1.0
    @test all(isfinite, tire_deformation_rates(tire, state))
    limited = SpatialTires.limited_tire_forces(tire, state[n], 1e6, 1e6)
    @test hypot(limited[1] / tire.mu_longitudinal,
        limited[2] / tire.mu_lateral) ≈ state[n]

    derivative = SpatialSimulationRunner.initial_spatial_derivative(
        loaded.initial_values, loaded)
    selection = AnalysisSelection(Dynamics(), loaded.active_variable_indices,
        loaded.active_equation_indices)
    residual = zeros(length(selection.equation_indices))
    evaluate_analysis_equations!(residual, loaded.model, selection, 0.0,
        loaded.initial_values, derivative)
    @test norm(residual, Inf) < 1e-9
    @test all(isfinite, Matrix(evaluate_analysis_sparse_jacobian(
        loaded.model, selection, 0.0, loaded.initial_values, derivative, 1.0)))
    # Check both constitutive and differential rows against perturbations,
    # including a nonzero shear and unloading rate away from law corners.
    probe = copy(loaded.initial_values)
    probe[n] = 3200.0
    probe[fx] = 0.002
    probe[fy] = 0.001
    probe[sx] = 0.2
    probe[sy] = 0.1
    probe[vx] = 2.0
    probe[tire.deflection_rate_variable] = -0.01
    coefficient = 3.0
    analytical = Matrix(evaluate_analysis_sparse_jacobian(loaded.model,
        selection, 0.0, probe, derivative, coefficient))
    row_map = Dict(equation => i for (i, equation) in
        enumerate(selection.equation_indices))
    column_map = Dict(variable => i for (i, variable) in
        enumerate(selection.variable_indices))
    rows = [row_map[i] for i in [tire.load_equations[2:3];
        tire.deformation_equations]]
    variables = [n, fx, fy, sx, sy, vx, tire.deflection_rate_variable,
        loaded.bodies[:wheel].angular_velocity_variables[2]]
    for variable in variables
        plus, minus = copy(probe), copy(probe)
        rate_plus, rate_minus = copy(derivative), copy(derivative)
        h = 1e-6
        plus[variable] += h
        minus[variable] -= h
        rate_plus[variable] += coefficient * h
        rate_minus[variable] -= coefficient * h
        residual_plus, residual_minus = similar(residual), similar(residual)
        evaluate_analysis_equations!(residual_plus, loaded.model, selection,
            0.0, plus, rate_plus)
        evaluate_analysis_equations!(residual_minus, loaded.model, selection,
            0.0, minus, rate_minus)
        @test analytical[rows, column_map[variable]] ≈
            (residual_plus[rows] - residual_minus[rows]) / (2h) atol=1e-5
    end
    result = run_spatial_model(IOBuffer(source); duration=0.025, samples=3)
    @test length(result.states) == 3
    @test all(isfinite, last(result.states))

    parked = run_spatial_model(joinpath(@__DIR__, "..", "..", "models",
        "spatial", "bristle-tire-side-slope.toml"))
    parked_tire = parked.loaded.forces[:tire]
    @test parked.static_relaxation_cycles >= 2
    @test first(parked.states)[parked_tire.normal_force_variable] ≈
        98.1 * cosd(10) atol=0.5
    @test abs(first(parked.states)[parked_tire.lateral_force_variable]) ≈
        98.1 * sind(10) atol=0.5
    @test last(parked.states)[parked_tire.lateral_force_variable] ≈
        first(parked.states)[parked_tire.lateral_force_variable] atol=0.1

    bounce = run_spatial_model(joinpath(@__DIR__, "..", "..", "models",
        "spatial", "bristle-tire-liftoff.toml"))
    bounce_tire = bounce.loaded.forces[:tire]
    normal = [z[bounce_tire.normal_force_variable] for z in bounce.states]
    lateral = [z[bounce_tire.lateral_force_variable] for z in bounce.states]
    shear = [z[bounce_tire.lateral_deformation_variable] for z in bounce.states]
    @test normal[1] > 0 && abs(lateral[1]) > 0
    @test all(iszero, normal[8:20])
    @test all(iszero, lateral[8:20])
    @test abs(shear[20]) < 1e-8
    touchdown = findfirst(index -> index > 20 && normal[index] > 0,
        eachindex(normal))
    @test !isnothing(touchdown)
    @test abs(shear[touchdown - 1]) < 1e-8
    @test normal[touchdown] > 0
    @test all(isfinite, reduce(vcat, bounce.states))
    @test all(eachindex(normal)) do index
        abs(lateral[index]) <= bounce_tire.mu_lateral * normal[index] + 1e-8
    end
    @test bounce.solution.stats.events_found == 0
    @test bounce.solution.stats.history_restarts == 0
    unloaded = bounce.states[10]
    unloaded_derivative = SpatialSimulationRunner.initial_spatial_derivative(
        unloaded, bounce.loaded)
    unloaded_selection = AnalysisSelection(Dynamics(),
        bounce.loaded.active_variable_indices,
        bounce.loaded.active_equation_indices)
    @test all(isfinite, Matrix(evaluate_analysis_sparse_jacobian(
        bounce.loaded.model, unloaded_selection, bounce.times[10],
        unloaded, unloaded_derivative, 1.0)))

    fore_aft = run_spatial_model(joinpath(@__DIR__, "..", "..", "models",
        "spatial", "bristle-tire-fore-aft-liftoff.toml"))
    fore_tire = fore_aft.loaded.forces[:tire]
    brake = fore_aft.loaded.forces[:brake]
    axle = fore_aft.loaded.connections[:axle]
    fore_first = first(fore_aft.states)
    fore_normal = [z[fore_tire.normal_force_variable]
        for z in fore_aft.states]
    fore_force = [z[fore_tire.longitudinal_force_variable]
        for z in fore_aft.states]
    fore_shear = [z[fore_tire.longitudinal_deformation_variable]
        for z in fore_aft.states]
    @test fore_aft.static_relaxation_cycles >= 2
    @test fore_first[fore_tire.normal_force_variable] ≈
        98.1 * cosd(10) atol=0.5
    @test fore_first[fore_tire.longitudinal_force_variable] ≈
        -98.1 * sind(10) atol=0.5
    @test abs(fore_first[fore_tire.lateral_force_variable]) < 1e-5
    @test fore_first[brake.magnitude_variable] ≈
        fore_first[fore_tire.longitudinal_force_variable] *
        (fore_tire.radius - fore_first[fore_tire.deflection_variable])
        atol=0.1
    @test fore_first[axle.hinge.rotation_variables[3]] ≈
        -fore_first[brake.magnitude_variable] / brake.stiffness atol=1e-4
    @test all(iszero, fore_normal[8:20])
    @test maximum(abs, fore_force[8:20]) < 1e-6
    @test abs(fore_shear[20]) < 1e-8
    fore_touchdown = findfirst(index -> index > 20 && fore_normal[index] > 0,
        eachindex(fore_normal))
    @test !isnothing(fore_touchdown)
    @test abs(fore_shear[fore_touchdown - 1]) < 1e-8
    @test all(eachindex(fore_normal)) do index
        abs(fore_force[index]) <= fore_tire.mu_longitudinal *
            fore_normal[index] + 1e-8
    end
    @test fore_aft.solution.stats.events_found == 0
    @test fore_aft.solution.stats.history_restarts == 0

    static_source = replace(source, "mode = \"dynamic\"" =>
        "mode = \"static\"")
    @test_throws ArgumentError load_spatial_model(IOBuffer(static_source))
    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "shear_release_time = 0.01" =>
            "shear_release_time = 0.01\nlongitudinal_expression = \"0\"")))
end
