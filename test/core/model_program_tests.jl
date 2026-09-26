using PracticalMechanicalSimulation
using PracticalMechanicalSimulation.AutomaticAnalysis
using PracticalMechanicalSimulation.PlanarComponentAssembly
using PracticalMechanicalSimulation.PlanarAppliedForces
using PracticalMechanicalSimulation.SimulationRunner
using HDF5
using LinearAlgebra
using SciMLBase

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const MODEL_DIRECTORY = joinpath(PROJECT_ROOT, "models", "planar")

@testset "Planar distance measurements" begin
    path = joinpath(MODEL_DIRECTORY, "distance-measures.toml")
    loaded = load_planar_model(path)
    measure = loaded.measures[:range]
    coordinate = loaded.connections[:height]
    initial = loaded.initial_values

    @test measure isa PlanarSpanMeasureComponent
    @test coordinate isa PlanarDistanceCoordinateComponent
    @test loaded.analysis.degrees_of_freedom == 3
    @test Symbol("range.velocity") ∉
        loaded.state_selection.selected_velocities
    @test length(loaded.active_variable_indices) ==
        length(loaded.active_equation_indices)

    separation = [1.0, 0.5]
    velocity = [0.2, -0.1]
    distance = norm(separation)
    @test initial[measure.spanning_variables] ≈ separation
    @test initial[measure.length_variable] ≈ distance
    @test initial[measure.unit_variables] ≈ separation / distance
    @test initial[measure.length_rate_variable] ≈
        dot(separation / distance, velocity)
    @test initial[measure.length_acceleration_variable] ≈
        (dot(velocity, velocity) -
         initial[measure.length_rate_variable]^2) / distance
    @test initial[coordinate.distance_variable] ≈ 0.5
    @test initial[coordinate.velocity_variable] ≈ -0.1

    force = only(loaded.forces[:measurement_probe])
    @test sort(force.magnitude.dependencies) == sort([
        measure.length_variable, coordinate.velocity_variable])

    derivative = SimulationRunner.initial_derivative(initial, loaded)
    selection = AnalysisSelection(Dynamics(), loaded.active_variable_indices,
        loaded.active_equation_indices)
    equations = zeros(length(selection.equation_indices))
    evaluate_analysis_equations!(equations, loaded.model, selection, 0.0,
        initial, derivative)
    @test norm(equations, Inf) < 1.0e-12

    coefficient = 2.1
    analytical = Matrix(evaluate_analysis_sparse_jacobian(loaded.model,
        selection, 0.0, initial, derivative, coefficient))
    numerical = similar(analytical)
    step = 1.0e-7
    plus, minus = similar(equations), similar(equations)
    for (column, variable) in enumerate(loaded.active_variable_indices)
        state_plus, state_minus = copy(initial), copy(initial)
        rate_plus, rate_minus = copy(derivative), copy(derivative)
        state_plus[variable] += step
        state_minus[variable] -= step
        rate_plus[variable] += coefficient * step
        rate_minus[variable] -= coefficient * step
        evaluate_analysis_equations!(plus, loaded.model, selection, 0.0,
            state_plus, rate_plus)
        evaluate_analysis_equations!(minus, loaded.model, selection, 0.0,
            state_minus, rate_minus)
        numerical[:, column] .= (plus .- minus) ./ (2step)
    end
    @test norm(analytical - numerical, Inf) < 3.0e-7

    result = run_planar_model(path; duration = 0.2, samples = 11)
    @test length(result.states) == 11
    monitor = result.solution.error_monitor
    active = result.loaded.active_variable_indices
    @test !monitor.mask[only(findall(==(measure.length_variable), active))]
    @test !monitor.mask[only(findall(
        ==(measure.length_rate_variable), active))]
    @test monitor.mask[only(findall(==(coordinate.velocity_variable), active))]
    for state in result.states
        marker_1 = PlanarAppliedForces.point_marker_kinematics(
            measure.marker_1, state)
        marker_2 = PlanarAppliedForces.point_marker_kinematics(
            measure.marker_2, state)
        @test state[measure.length_variable] ≈
            norm(marker_2.position - marker_1.position) atol = 2.0e-6
    end

    coincident = replace(read(path, String),
        "position = [1.0, 0.5]" => "position = [0.0, 0.0]")
    @test_throws ArgumentError load_planar_model(IOBuffer(coincident))
    preferred_span = read(path, String) * """

    [state_selection]
    method = "preferred"
    preferred_velocities = ["range.velocity"]
    """
    @test_throws ArgumentError load_planar_model(IOBuffer(preferred_span))
end

@testset "TOML model program" begin
    slider_path = joinpath(MODEL_DIRECTORY,
        "constant-speed-slider-crank.toml")
    slider = run_planar_model(slider_path; duration = 0.1, samples = 3)
    @test slider.analysis_mode == :kinematic
    @test slider.loaded.analysis.degrees_of_freedom == 0
    @test length(slider.states) == 3
    @test maximum(abs, slider.states[end] - slider.states[1]) > 0
    @test slider.loaded.analysis.formulation == :StateSelected
    @test slider.loaded.analysis.deficit == 0
    @test slider.loaded.analysis.initialization == :none
    @test slider.loaded.analysis.static_method == :newton
    @test slider.loaded.simulation.output_precision == :single
    double_precision_slider = load_planar_model(IOBuffer(replace(
        read(slider_path, String), "output_samples = 81" =>
            "output_samples = 81\noutput_precision = \"double\"")))
    @test double_precision_slider.simulation.output_precision == :double
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        read(slider_path, String), "output_samples = 81" =>
            "output_samples = 81\noutput_precision = \"extended\"")))
    degree_slider = load_planar_model(IOBuffer(replace(
        read(slider_path, String), "initial_angle = 0.0" =>
            "initial_angle = \"15°\"")))
    @test degree_slider.drivers[:crank_drive].motion(0.0) ≈ pi / 12
    @test_throws ArgumentError load_planar_model(IOBuffer(
        read(slider_path, String) * "\n[analysis]\ndeficit = 0\n"))

    static_source = read(slider_path, String) *
        "\n[analysis]\nmode = \"static\"\n"
    static_slider = run_planar_model(IOBuffer(static_source);
        duration = 0.1, samples = 3)
    @test static_slider.analysis_mode == :static
    @test length(static_slider.static_iterations) >= 3
    for (kinematic, static) in zip(slider.states, static_slider.states)
        for body in values(slider.loaded.bodies)
            @test static[body.position_variables] ≈
                kinematic[body.position_variables] atol = 2.0e-10
            @test static[body.orientation_variable] ≈
                kinematic[body.orientation_variable] atol = 2.0e-10
            @test all(iszero, static[body.velocity_variables])
            @test all(iszero, static[body.acceleration_variables])
        end
    end

    four_bar = run_planar_model(joinpath(MODEL_DIRECTORY,
        "torque-driven-four-bar.toml"); duration = 0.05, samples = 3)
    @test four_bar.analysis_mode == :dynamic
    @test four_bar.loaded.analysis.degrees_of_freedom == 1
    @test count(four_bar.solution.differential_vars) == 2

    modal = run_planar_model(joinpath(MODEL_DIRECTORY,
        "modal-pendulum.toml"))
    modal_body = modal.loaded.bodies[:pendulum]
    modal_pin = modal.loaded.connections[:pin]
    modal_spring = only(modal.loaded.forces[:pin_spring]).element
    effective_inertia = modal_body.inertia + modal_body.mass * 0.5^2
    effective_stiffness = modal_spring.stiffness + modal_body.mass * 9.81 * 0.5
    expected_decay = -modal_spring.damping / (2effective_inertia)
    expected_frequency = sqrt(effective_stiffness / effective_inertia -
        expected_decay^2)
    @test modal.analysis_mode == :modal
    @test modal.loaded.analysis.number_of_modes == 1
    @test modal.loaded.analysis.frequency_shift_hz == 0.0
    @test length(modal.eigenvalues) == 1
    @test only(modal.eigenvalues) ≈
        expected_decay + im * expected_frequency rtol = 2.0e-12
    @test only(modal.natural_frequencies_hz) ≈
        sqrt(effective_stiffness / effective_inertia) / (2pi) rtol = 2.0e-12
    @test only(modal.equation_errors) < 1.0e-12
    @test modal.sparse_factorizations == 1
    @test length(modal.differential_columns) == 2
    @test size(modal.mode_shapes) ==
        (length(modal.loaded.initial_values), 1)
    modal_theta = modal.mode_shapes[modal_pin.rotation_variables[3], 1]
    modal_torque = modal.mode_shapes[modal_spring.torque_variable, 1]
    @test modal_theta ≈ 1.0 atol = 1.0e-12
    @test modal_torque ≈
        -(modal_spring.stiffness +
          modal_spring.damping * only(modal.eigenvalues)) * modal_theta
        rtol = 2.0e-12

    modal_source = read(joinpath(MODEL_DIRECTORY,
        "modal-pendulum.toml"), String)
    instantaneous_source = replace(modal_source,
        "position = [0.0, -0.5]\norientation = 0.0" =>
        "position = [0.3535533905932738, -0.3535533905932738]\n" *
        "orientation = 0.7853981633974483")
    instantaneous = run_planar_model(IOBuffer(instantaneous_source))
    instantaneous_body = instantaneous.loaded.bodies[:pendulum]
    @test norm(instantaneous.operating_state[
        instantaneous_body.acceleration_variables]) > 1.0
    @test only(instantaneous.equation_errors) < 1.0e-12
    @test only(instantaneous.eigenvalues) != only(modal.eigenvalues)

    shifted_source = replace(read(joinpath(MODEL_DIRECTORY,
        "bushing-supported-body.toml"), String),
        "mode = \"automatic\"" =>
        "mode = \"modal\"\nmodes = 3\nfrequency_shift_hz = 0.5")
    shifted_modal = run_planar_model(IOBuffer(shifted_source))
    @test length(shifted_modal.eigenvalues) == 3
    @test shifted_modal.finite_eigenvalue_count == 6
    @test length(shifted_modal.differential_columns) == 6
    @test maximum(shifted_modal.equation_errors) < 1.0e-12

    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        modal_source, "modes = 1" => "modes = 0")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        modal_source, "frequency_shift_hz = 0.0" =>
            "frequency_shift_hz = -1.0")))
    @test_throws ArgumentError run_planar_model(IOBuffer(
        read(slider_path, String) * "\n[analysis]\nmode = \"modal\"\n"))

    bushing_path = joinpath(MODEL_DIRECTORY, "bushing-supported-body.toml")
    bushing = run_planar_model(bushing_path; duration = 0.1, samples = 3)
    body = bushing.loaded.bodies[:body]
    support = only(bushing.loaded.forces[:support])
    @test bushing.loaded.analysis.degrees_of_freedom == 3
    @test bushing.states[end][body.position_variables[2]] < -1.0
    @test bushing.states[end][support.force_variables[2]] > 0.0
    @test support.damping_time_scale == 0.1
    @test support.translational_damping == [10.0, 10.0]
    @test support.rotational_damping == 2.0
    legacy_angle_source = replace(read(bushing_path, String),
        "orientation = 0.2" => "angle = 0.2")
    legacy_angle = load_planar_model(IOBuffer(legacy_angle_source))
    @test legacy_angle.initial_values[
        legacy_angle.bodies[:body].orientation_variable] ≈ 0.2
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        read(bushing_path, String), "orientation = 0.2" =>
            "orientation = 0.2\nangle = 0.2")))

    directed_force = run_planar_model(joinpath(MODEL_DIRECTORY,
        "marker-directed-force.toml"); duration = 0.02, samples = 3)
    applied = only(directed_force.loaded.forces[:push])
    link = directed_force.loaded.bodies[:link]
    initial_directed = first(directed_force.states)
    @test applied isa PracticalMechanicalSimulation.PlanarComponentAssembly.PlanarAppliedForceComponent
    @test applied.direction_axis isa
        PracticalMechanicalSimulation.PlanarDirectedDistances.PlanarDirectedAxis
    @test isnothing(applied.reaction_marker)
    @test !haskey(directed_force.loaded.markers, Symbol("push.reaction"))
    @test initial_directed[link.angular_acceleration_variable] ≈ 4.5 atol = 1.0e-9

    no_reaction_source = replace(read(joinpath(MODEL_DIRECTORY,
        "marker-directed-force.toml"), String),
        "force = 1.0" => "expression = \"1.5*cos(t)\"",
        "markers = [\"link.tip\", \"ground.force_direction\"]" =>
            "markers = [\"link.tip\", \"link.tip\"]")
    no_reaction = run_planar_model(IOBuffer(no_reaction_source);
        duration = 0.0, samples = 1)
    no_reaction_force = only(no_reaction.loaded.forces[:push])
    @test isnothing(no_reaction_force.reaction_marker)
    @test !haskey(no_reaction.loaded.markers, Symbol("push.reaction"))
    @test no_reaction_force.magnitude(0.0) == 1.5

    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        no_reaction_source, "expression = \"1.5*cos(t)\"" =>
            "force = 1.5\nexpression = \"1.5*cos(t)\"")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        no_reaction_source, "markers = [\"link.tip\", \"link.tip\"]" =>
            "markers = [\"link.tip\"]")))

    reaction_source = read(slider_path, String) * """

    [push]
    type = "applied_force"
    markers = ["crank.end", "ground.origin"]
    force = 1.5
    reaction_body = "rod"
    """
    reaction_model = run_planar_model(IOBuffer(reaction_source);
        duration = 0.0, samples = 1)
    reaction_force = only(reaction_model.loaded.forces[:push])
    reaction_marker = reaction_model.loaded.markers[Symbol("push.reaction")]
    reaction_state = only(reaction_model.states)
    application_position = PracticalMechanicalSimulation.PlanarAppliedForces.point_marker_kinematics(
        reaction_force.application_marker, reaction_state).position
    reaction_position = PracticalMechanicalSimulation.PlanarAppliedForces.point_marker_kinematics(
        reaction_force.reaction_marker, reaction_state).position
    @test reaction_marker.owner === reaction_model.loaded.bodies[:rod]
    @test reaction_position ≈ application_position atol = 1.0e-12
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        reaction_source, "reaction_body = \"rod\"" =>
            "reaction_body = \"crank\"")))

    bushing_pendulum = run_planar_model(joinpath(MODEL_DIRECTORY,
        "bushing-pendulum.toml"); duration = 0.05, samples = 3)
    pendulum_body = bushing_pendulum.loaded.bodies[:pendulum]
    pin_bushing = only(bushing_pendulum.loaded.forces[:pin_bushing])
    pendulum_initial = first(bushing_pendulum.states)
    @test bushing_pendulum.loaded.analysis.degrees_of_freedom == 3
    @test bushing_pendulum.loaded.analysis.static_method == :dynamic_relaxation
    @test bushing_pendulum.static_relaxation_cycles == 1
    @test pin_bushing.rotational_stiffness == 0.0
    @test pin_bushing.rotational_damping == 0.0
    @test pendulum_initial[pendulum_body.position_variables] ≈
        [0.0, -0.50981] atol = 1.0e-11
    @test pendulum_initial[pendulum_body.velocity_variables] ≈
        [0.5, 0.0] atol = 1.0e-12
    @test last(bushing_pendulum.states)[pendulum_body.orientation_variable] >
        pendulum_initial[pendulum_body.orientation_variable] + 0.01
    static_relaxation_source = replace(read(joinpath(MODEL_DIRECTORY,
        "bushing-pendulum.toml"), String),
        "mode = \"dynamic\"\ninitialization = \"static_equilibrium\"" =>
            "mode = \"static\"")
    static_relaxation = run_planar_model(
        IOBuffer(static_relaxation_source); samples = 1)
    static_pendulum_body = static_relaxation.loaded.bodies[:pendulum]
    @test static_relaxation.static_relaxation_cycles == [1]
    @test only(static_relaxation.states)[static_pendulum_body.position_variables] ≈
        [0.0, -0.50981] atol = 1.0e-11

    contact_path = joinpath(MODEL_DIRECTORY, "bouncing-ball.toml")
    bounce = run_planar_model(contact_path; duration = 0.6, samples = 31)
    contact = only(bounce.loaded.forces[:floor_contact])
    @test contact.geometry isa
        PracticalMechanicalSimulation.PlanarDirectedDistances.PlanarDirectedDistance
    @test contact.radius == 0.1
    @test contact.damping_factor == 0.15
    @test bounce.loaded.analysis.degrees_of_freedom == 3
    @test minimum(state[contact.gap_variable] for state in bounce.states) < 0
    @test maximum(state[contact.normal_force_variable]
                  for state in bounce.states) > 100
    @test length(bounce.solution.events) == 2
    @test first(bounce.solution.events).time ≈ sqrt(2 * 0.9 / 9.81) atol = 2e-6
    @test bounce.solution.stats.history_restarts == 2
    maximum_force_state = argmax(state -> state[contact.normal_force_variable],
        bounce.states)
    penetration = max(-maximum_force_state[contact.gap_variable], 0.0)
    closing_speed = -maximum_force_state[contact.gap_rate_variable]
    @test maximum_force_state[contact.normal_force_variable] ≈
        contact.stiffness * penetration *
        max(0.0, 1 + contact.damping_factor * closing_speed) rtol = 2e-5
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        read(contact_path, String),
        "stiffness = 10000.0" => "stiffness = 0.0")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        read(contact_path, String), "radius = 0.1\n" => "")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        read(contact_path, String), "radius = 0.1" => "radius = 0.0")))
    elastic_contact = only(load_planar_model(IOBuffer(replace(
        read(contact_path, String), "damping_factor = 0.15\n" =>
        ""))).forces[:floor_contact])
    @test elastic_contact.damping_factor == 0.0
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        read(contact_path, String), "damping_factor = 0.15" =>
        "damping_factor = -0.1")))

    explicit_damping_source = replace(read(bushing_path, String),
        "damping_time_scale = 0.1" => """
        damping_time_scale = 0.5
        translational_damping = [3.0, 4.0]
        rotational_damping = 6.0
        """)
    explicit_damping = only(load_planar_model(
        IOBuffer(explicit_damping_source)).forces[:support])
    @test explicit_damping.damping_time_scale == 0.5
    @test explicit_damping.translational_damping == [3.0, 4.0]
    @test explicit_damping.rotational_damping == 6.0
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        read(bushing_path, String),
        "damping_time_scale = 0.1" => "damping_time_scale = -0.1")))

    static_bushing_source = replace(read(bushing_path, String),
        "mode = \"automatic\"" => "mode = \"static\"")
    equilibrium = run_planar_model(IOBuffer(static_bushing_source); samples = 1)
    equilibrium_state = only(equilibrium.states)
    @test equilibrium_state[body.position_variables] ≈
        [0.0, -1.0981] atol = 1.0e-12
    @test equilibrium_state[support.force_variables] ≈
        [0.0, 9.81] atol = 1.0e-12

    initialized_source = replace(read(bushing_path, String),
        "mode = \"automatic\"" =>
        "mode = \"dynamic\"\ninitialization = \"static_equilibrium\"")
    initialized = run_planar_model(IOBuffer(initialized_source);
        duration = 0.1, samples = 3)
    initialized_body = initialized.loaded.bodies[:body]
    @test initialized.analysis_mode == :dynamic
    @test initialized.loaded.analysis.initialization == :static_equilibrium
    @test initialized.static_initialization_iterations == 1
    @test initialized.states[1][initialized_body.position_variables] ≈
        [0.0, -1.0981] atol = 1.0e-12
    @test initialized.states[1][initialized_body.orientation_variable] ≈
        0.0 atol = 1.0e-12
    @test initialized.states[end][initialized_body.position_variables] ≈
        initialized.states[1][initialized_body.position_variables] atol = 1e-11

    moving_initialized_source = replace(initialized_source,
        "orientation = 0.2" =>
        "orientation = 0.2\nvelocity = [0.0, 0.25]\nangular_velocity = -0.1")
    moving_initialized = run_planar_model(IOBuffer(moving_initialized_source);
        duration = 0.0, samples = 1)
    moving_body = moving_initialized.loaded.bodies[:body]
    @test only(moving_initialized.states)[moving_body.position_variables] ≈
        [0.0, -1.0981] atol = 1.0e-12
    @test only(moving_initialized.states)[moving_body.velocity_variables] ==
        [0.0, 0.25]
    @test only(moving_initialized.states)[moving_body.angular_velocity_variable] ==
        -0.1

    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        read(bushing_path, String), "mode = \"automatic\"" =>
        "mode = \"dynamic\"\ninitialization = \"unknown\"")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        read(bushing_path, String), "mode = \"automatic\"" =>
        "mode = \"static\"\ninitialization = \"static_equilibrium\"")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        initialized_source, "initialization = \"static_equilibrium\"" =>
        "initialization = \"static_equilibrium\"\nstatic_method = \"unknown\"")))
    kinematic_initialization_source = read(slider_path, String) *
        "\n[analysis]\nmode = \"automatic\"\n" *
        "initialization = \"static_equilibrium\"\n"
    @test_throws ArgumentError run_planar_model(
        IOBuffer(kinematic_initialization_source); duration = 0.1, samples = 3)

    mktempdir() do directory
        result_path = joinpath(directory, "slider.simp")
        @test write_result(result_path, slider) == result_path
        stored = read_result(result_path)
        @test stored.status == :complete
        @test stored.output_precision == :single
        h5open(result_path, "r") do file
            @test eltype(file["results/values"]) == Float32
        end
        @test isempty(stored.status_message)
        @test stored.times == slider.times
        @test isapprox(stored.values,
            reduce(vcat, permutedims.(slider.states));
            rtol = 1.0e-6, atol = 1.0e-7)
        @test stored.analysis_mode == :kinematic
        @test stored.solver_statistics.events_found == 0
        @test stored.solver_statistics.symbolic_factorizations ==
            slider.solution.stats.symbolic_factorizations
        @test stored.solver_statistics.corrector_failures ==
            slider.solution.stats.corrector_failures
        @test isempty(stored.modal_eigenvalues)
        @test size(stored.modal_mode_shapes) == (0, 0)

        partial_path = joinpath(directory, "partial.simp")
        writer = begin_incremental_result(partial_path, slider.loaded,
            :kinematic)
        record_incremental_event!(writer, (; kind = :sample, time = 0.1,
            state = copy(slider.states[2])))
        finish_incremental_result!(writer, :failed;
            message = "test failure", failure_time = 0.1)
        partial = read_result(partial_path)
        @test partial.status == :failed
        @test partial.status_message == "test failure"
        @test partial.failure_time == 0.1
        @test partial.times == [0.0, 0.1]
        @test partial.output_precision == :single
        @test isapprox(partial.values[2, :], slider.states[2];
            rtol = 1.0e-6, atol = 1.0e-7)

        abandoned_path = joinpath(directory, "abandoned.simp")
        abandoned_writer = begin_incremental_result(abandoned_path,
            slider.loaded, :kinematic; flush_interval = 1)
        record_incremental_event!(abandoned_writer, (;
            kind = :sample, time = 0.1,
            state = copy(slider.states[2])))
        @test finalize_result!(abandoned_path) == abandoned_path
        abandoned = read_result(abandoned_path)
        @test abandoned.status == :interrupted
        @test abandoned.status_message == "run finalized by the user"
        @test abandoned.failure_time == 0.1
        @test abandoned.times == [0.0, 0.1]
        @test finalize_result!(abandoned_path) == abandoned_path
        @test_throws ArgumentError finalize_result!(result_path)

        collector = PracticalMechanicalSimulation.ResultIO.HealthPeakCollector()
        PracticalMechanicalSimulation.ResultIO.record_health_step!(collector,
            0.1, [1.0, 2.0], 2, 0.01, 0.5, 5.0, 1)
        PracticalMechanicalSimulation.ResultIO.record_health_step!(collector,
            0.11, [3.0, 4.0], 2, 0.01, 0.5, 30.0, 2)
        PracticalMechanicalSimulation.ResultIO.record_health_step!(collector,
            0.12, [5.0, 6.0], 2, 0.01, 0.5, 40.0, 1)
        PracticalMechanicalSimulation.ResultIO.record_health_step!(collector,
            0.13, [7.0, 8.0], 2, 0.01, 0.5, 1.0, 1)
        @test length(collector.peaks) == 1
        @test only(collector.peaks).time == 0.12
        @test only(collector.peaks).active_values == [5.0, 6.0]
        @test collector.maximum.error == 40.0

        double_path = joinpath(directory, "slider-double.simp")
        write_result(double_path, slider; output_precision = :double)
        @test read_result(double_path).output_precision == :double
        h5open(double_path, "r") do file
            @test eltype(file["results/values"]) == Float64
        end

        @test length(four_bar.solution.t) == 2
        @test length(four_bar.states) == length(four_bar.times)
        monitored_variable = first(findall(
            four_bar.solution.error_monitor.mask))
        health_time = four_bar.times[2]
        health_order = 3
        health_step_size = 0.01
        raw_peak = PracticalMechanicalSimulation.ResultIO.RawHealthPeak(
            health_time,
            four_bar.states[2][four_bar.loaded.active_variable_indices],
            30.0, 0.5, monitored_variable, health_order,
            health_step_size)
        monitored_four_bar = merge(four_bar,
            (; health_step_peaks = [raw_peak]))
        health_path = joinpath(directory, "health.simp")
        write_result(health_path, monitored_four_bar)
        stored_health = read_result(health_path)
        @test length(stored_health.health_snapshots) == 1
        health_snapshot = only(stored_health.health_snapshots)
        @test health_snapshot.time == health_time
        @test health_snapshot.physical_error == 30.0
        @test health_snapshot.controlled_error == 0.5
        @test health_snapshot.amplification == 60.0
        @test health_snapshot.severity == :check
        @test health_snapshot.order == health_order
        @test health_snapshot.step_size == health_step_size
        @test isapprox(health_snapshot.values[
            four_bar.loaded.active_variable_indices],
            four_bar.states[2][four_bar.loaded.active_variable_indices];
            rtol = 1.0e-6, atol = 1.0e-7)
        @test !isempty(health_snapshot.selected_variables)

        bounce_path = joinpath(directory, "bounce.simp")
        write_result(bounce_path, bounce)
        stored_bounce = read_result(bounce_path)
        @test stored_bounce.solver_statistics.events_found == 2
        @test stored_bounce.solver_statistics.history_restarts == 2

        initialized_path = joinpath(directory, "initialized.simp")
        write_result(initialized_path, initialized)
        stored_initialized = read_result(initialized_path)
        @test stored_initialized.static_initialization_iterations == 1
        @test stored_initialized.static_relaxation_cycles == 0

        relaxed_path = joinpath(directory, "relaxed.simp")
        write_result(relaxed_path, bushing_pendulum)
        @test read_result(relaxed_path).static_relaxation_cycles == 1

        modal_path = joinpath(directory, "modal.simp")
        write_result(modal_path, modal)
        stored_modal = read_result(modal_path)
        @test stored_modal.analysis_mode == :modal
        @test stored_modal.modal_eigenvalues ≈ modal.eigenvalues
        @test stored_modal.modal_natural_frequencies_hz ≈
            modal.natural_frequencies_hz
        @test stored_modal.modal_damped_frequencies_hz ≈
            modal.damped_frequencies_hz
        @test stored_modal.modal_damping_ratios ≈ modal.damping_ratios
        @test isapprox(stored_modal.modal_mode_shapes, modal.mode_shapes;
            rtol = 1.0e-6, atol = 1.0e-7)
        @test stored_modal.modal_equation_errors ≈ modal.equation_errors
        @test stored_modal.modal_shift ≈ modal.shift

        modal_theta_index = modal_pin.rotation_variables[3]
        modal_theta_name = string(
            stored_modal.variable_components[modal_theta_index], ".",
            stored_modal.variable_names[modal_theta_index])
        modal_csv_stem = joinpath(directory, "pendulum.csv")
        modal_csv_paths = export_result_csv(modal_csv_stem, stored_modal;
            variables = [modal_theta_name])
        @test modal_csv_paths == (
            modes = joinpath(directory, "pendulum-modes.csv"),
            mode_shapes = joinpath(directory, "pendulum-mode-shapes.csv"))

        mode_lines = readlines(modal_csv_paths.modes)
        @test first(mode_lines) ==
            "mode,eigenvalue_real,eigenvalue_imaginary,natural_frequency_hz," *
            "damped_frequency_hz,damping_ratio,equation_error"
        @test length(mode_lines) == 2
        mode_fields = split(last(mode_lines), ',')
        @test parse(Int, mode_fields[1]) == 1
        @test parse(Float64, mode_fields[2]) ≈ real(only(modal.eigenvalues))
        @test parse(Float64, mode_fields[3]) ≈ imag(only(modal.eigenvalues))
        @test parse(Float64, mode_fields[4]) ≈
            only(modal.natural_frequencies_hz)
        @test parse(Float64, mode_fields[7]) ≈ only(modal.equation_errors)

        shape_lines = readlines(modal_csv_paths.mode_shapes)
        @test first(shape_lines) ==
            "mode,component,variable,kind,real,imaginary,magnitude,phase_deg"
        @test length(shape_lines) == 2
        shape_fields = split(last(shape_lines), ',')
        modal_theta = modal.mode_shapes[modal_theta_index, 1]
        @test parse(Int, shape_fields[1]) == 1
        @test shape_fields[2] ==
            stored_modal.variable_components[modal_theta_index]
        @test shape_fields[3] ==
            stored_modal.variable_names[modal_theta_index]
        @test shape_fields[4] ==
            stored_modal.variable_kinds[modal_theta_index]
        @test parse(Float64, shape_fields[5]) ≈ real(modal_theta)
        @test parse(Float64, shape_fields[6]) ≈ imag(modal_theta)
        @test parse(Float64, shape_fields[7]) ≈ abs(modal_theta)
        @test parse(Float64, shape_fields[8]) ≈ rad2deg(angle(modal_theta))
        @test_throws ArgumentError export_result_csv(modal_csv_stem,
            stored_modal; variables = [modal_theta_name])

        static_relaxed_path = joinpath(directory, "static-relaxed.simp")
        write_result(static_relaxed_path, static_relaxation)
        @test read_result(static_relaxed_path).solver_statistics.relaxation_cycles ==
            [1]

        static_initial_path = joinpath(directory, "static-initial.simp")
        write_result(static_initial_path, equilibrium)
        transferred_model_path = joinpath(directory, "transferred.toml")
        transferred_source = replace(read(bushing_path, String),
            "orientation = 0.2" =>
                "orientation = 0.2\nvelocity = [0.25, 0.0]") * """

            [initial_conditions]
            result = "static-initial.simp"
            """
        write(transferred_model_path, transferred_source)
        transferred = run_planar_model(transferred_model_path;
            duration = 0.0, samples = 1)
        transferred_body = transferred.loaded.bodies[:body]
        transferred_state = only(transferred.states)
        @test transferred.loaded.initial_conditions.enabled
        @test transferred.loaded.initial_conditions.result_path == static_initial_path
        @test transferred.loaded.initial_conditions.sample == 1
        @test transferred.loaded.initial_conditions.source_time ==
            only(equilibrium.times)
        @test !transferred.loaded.initial_conditions.include_velocities
        @test Set(transferred.loaded.initial_conditions.transferred_variables) ==
            Set((Symbol("body.R_x"), Symbol("body.R_y"), Symbol("body.theta")))
        @test transferred_state[transferred_body.position_variables] ≈
            equilibrium_state[body.position_variables] atol = 2.0e-7
        @test transferred_state[transferred_body.orientation_variable] ≈
            equilibrium_state[body.orientation_variable] atol = 2.0e-7
        @test transferred_state[transferred_body.velocity_variables] == [0.25, 0.0]

        moving_initial_path = joinpath(directory, "moving-initial.simp")
        write_result(moving_initial_path, moving_initialized)
        velocity_model_path = joinpath(directory, "transferred-velocity.toml")
        velocity_source = read(bushing_path, String) * """

            [initial_conditions]
            result = "moving-initial.simp"
            sample = "static"
            include_velocities = true
            """
        write(velocity_model_path, velocity_source)
        velocity_transfer = run_planar_model(velocity_model_path;
            duration = 0.0, samples = 1)
        velocity_body = velocity_transfer.loaded.bodies[:body]
        velocity_state = only(velocity_transfer.states)
        moving_state = only(moving_initialized.states)
        @test velocity_transfer.loaded.initial_conditions.include_velocities
        @test velocity_transfer.loaded.initial_conditions.requested_sample ==
            "static"
        @test velocity_transfer.loaded.initial_conditions.sample == 1
        @test velocity_state[velocity_body.position_variables] ≈
            moving_state[moving_body.position_variables] atol = 2.0e-7
        @test velocity_state[velocity_body.velocity_variables] ≈
            moving_state[moving_body.velocity_variables] atol = 2.0e-7
        @test velocity_state[velocity_body.angular_velocity_variable] ≈
            moving_state[moving_body.angular_velocity_variable] atol = 2.0e-7

        incompatible_path = joinpath(directory, "incompatible.toml")
        write(incompatible_path, read(bushing_path, String) * """

            [initial_conditions]
            result = "slider.simp"
            """)
        @test_throws ArgumentError load_planar_model(incompatible_path)

        invalid_sample_path = joinpath(directory, "invalid-sample.toml")
        write(invalid_sample_path, transferred_source * "sample = 0\n")
        @test_throws ArgumentError load_planar_model(invalid_sample_path)

        ordinary_dynamic_path = joinpath(directory, "ordinary-dynamic.simp")
        write_result(ordinary_dynamic_path, bushing)
        invalid_static_path = joinpath(directory, "invalid-static.toml")
        write(invalid_static_path, read(bushing_path, String) * """

            [initial_conditions]
            result = "ordinary-dynamic.simp"
            sample = "static"
            """)
        @test_throws ArgumentError load_planar_model(invalid_static_path)

        csv_path = joinpath(directory, "slider.csv")
        @test export_result_csv(csv_path, result_path;
            variables = ["crank.theta"]) == csv_path
        @test first(readlines(csv_path)) == "time,crank.theta"

        model_path = joinpath(directory, "recovered.toml")
        @test extract_model(model_path, result_path) == model_path
        @test read(model_path, String) == read(slider_path, String)
    end

    output = IOBuffer()
    @test planar_model_main([slider_path, "0.1", "3"];
        output = output) == 0
    @test occursin("analysis: kinematic", String(take!(output)))
end

@testset "State-dependent load expressions" begin
    assembly = PracticalMechanicalSimulation.PlanarComponentAssembly

    spring_source = read(joinpath(MODEL_DIRECTORY,
        "torsional-spring-pendulum.toml"), String)
    expression_source = replace(spring_source, """
        [pin_spring]
        type = "torsional_spring_damper"
        markers = ["pendulum.pin", "ground.origin"]
        stiffness = 3.0
        damping_time_scale = 0.1
        free_angle = 0.0
        """ => """
        [pin_spring]
        type = "applied_torque"
        markers = ["pendulum.pin", "ground.origin"]
        expression = "-3.0*pin.theta-0.3*pin.omega"
        """)
    reference = run_planar_model(IOBuffer(spring_source);
        duration = 0.1, samples = 5)
    expression = run_planar_model(IOBuffer(expression_source);
        duration = 0.1, samples = 5)
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        expression_source, "type = \"applied_torque\"" =>
            "type = \"expression_torque\"")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        expression_source,
        "expression = \"-3.0*pin.theta-0.3*pin.omega\"" =>
            "torque = 1.0\nexpression = \"-3.0*pin.theta-0.3*pin.omega\"")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        expression_source,
        "expression = \"-3.0*pin.theta-0.3*pin.omega\"" => "")))
    body = expression.loaded.bodies[:pendulum]
    pin = expression.loaded.connections[:pin]
    torque = only(expression.loaded.forces[:pin_spring])
    @test torque isa assembly.PlanarAppliedTorqueComponent
    @test torque.torque_variable != 0
    @test Symbol("pin_spring.T") in Symbol.(string(variable.component, ".",
        variable.name) for variable in expression.loaded.layout.catalog.variables)
    @test torque.torque.dependencies == sort([
        pin.rotation_variables[2], pin.rotation_variables[3]])
    initial_gradient = torque.torque.gradient(0.0, first(expression.states))
    expected_gradient = [index == pin.rotation_variables[2] ? -0.3 : -3.0
        for index in torque.torque.dependencies]
    @test initial_gradient ≈ expected_gradient atol = 1.0e-14
    for (reference_state, expression_state) in
            zip(reference.states, expression.states)
        @test expression_state[body.acceleration_variables] ≈
            reference_state[body.acceleration_variables] rtol = 2.0e-7
        @test expression_state[body.velocity_variables] ≈
            reference_state[body.velocity_variables] rtol = 2.0e-7
        @test expression_state[body.position_variables] ≈
            reference_state[body.position_variables] rtol = 2.0e-7
        @test expression_state[pin.rotation_variables] ≈
            reference_state[pin.rotation_variables] rtol = 2.0e-7
        @test expression_state[torque.torque_variable] ≈
            -3.0 * expression_state[pin.rotation_variables[3]] -
            0.3 * expression_state[pin.rotation_variables[2]] atol = 1.0e-9
    end

    directed_source = read(joinpath(MODEL_DIRECTORY,
        "marker-directed-force.toml"), String)
    coordinate_source = replace(directed_source,
        "[pin_damper]" => """
        [tip_height]
        type = "distance_coordinate"
        markers = ["link.tip", "ground.force_direction"]

        [pin_damper]""",
        "force = 1.0" =>
            "expression = \"2.0-3.0*tip_height.distance-0.5*tip_height.velocity\"")
    directed = run_planar_model(IOBuffer(coordinate_source);
        duration = 0.0, samples = 1)
    force = only(directed.loaded.forces[:push])
    coordinate = directed.loaded.connections[:tip_height]
    state = only(directed.states)
    @test force.magnitude_variable != 0
    @test state[force.magnitude_variable] ≈ 2.0 atol = 1.0e-12
    @test force.magnitude.dependencies == sort([
        coordinate.distance_variable, coordinate.velocity_variable])
    @test force.magnitude.gradient(0.0, state) ≈
        [index == coordinate.distance_variable ? -3.0 : -0.5
         for index in force.magnitude.dependencies] atol = 1.0e-14

    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        coordinate_source, "tip_height.distance" => "tip_height.missing")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        coordinate_source, "tip_height.distance" =>
            "tip_height.acceleration")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        expression_source, "rotation_coordinates = true" =>
            "rotation_coordinates = false")))
end

@testset "Coordinate couplers" begin
    assembly = PracticalMechanicalSimulation.PlanarComponentAssembly

    rotational_path = joinpath(MODEL_DIRECTORY,
        "rotational-coordinate-coupler.toml")
    rotational = run_planar_model(rotational_path;
        duration = 0.5, samples = 6)
    rotational_coupler = rotational.loaded.connections[:shaft_coupler]
    input_bearing = rotational.loaded.connections[:input_bearing]
    output_bearing = rotational.loaded.connections[:output_bearing]
    @test rotational_coupler isa assembly.PlanarCoordinateCoupler
    @test rotational_coupler.coordinate_kind == :rotation
    @test rotational.loaded.analysis.degrees_of_freedom == 1
    @test rotational.loaded.state_selection.preferred_velocities ==
        [Symbol("input_bearing.omega")]
    for state in rotational.states
        input_angle = state[input_bearing.rotation_variables[3]]
        output_angle = state[output_bearing.rotation_variables[3]]
        input_rate = state[input_bearing.rotation_variables[2]]
        output_rate = state[output_bearing.rotation_variables[2]]
        @test input_angle - 2output_angle ≈ rotational_coupler.offset atol = 1e-12
        @test input_rate - 2output_rate ≈ 0.0 atol = 1e-12
        @test state[rotational_coupler.reaction_variable] ≈
            -1 / 3 atol = 1e-12
    end

    translational_path = joinpath(MODEL_DIRECTORY,
        "translational-coordinate-coupler.toml")
    translational = run_planar_model(translational_path;
        duration = 0.5, samples = 6)
    translational_coupler = translational.loaded.connections[:slider_coupler]
    upper = translational.loaded.connections[:upper_distance]
    lower = translational.loaded.connections[:lower_distance]
    @test translational_coupler isa assembly.PlanarCoordinateCoupler
    @test upper isa assembly.PlanarDistanceCoordinateComponent
    @test translational_coupler.coordinate_kind == :distance
    @test translational.loaded.analysis.degrees_of_freedom == 1
    @test translational.loaded.state_selection.preferred_velocities ==
        [Symbol("upper_distance.velocity")]
    @test length(upper.candidate_state_equations) == 2
    for state in translational.states
        @test state[upper.distance_variable] + state[lower.distance_variable] ≈
            translational_coupler.offset atol = 1e-12
        @test state[upper.velocity_variable] + state[lower.velocity_variable] ≈
            0.0 atol = 1e-12
        @test state[translational_coupler.reaction_variable] ≈
            -0.6 atol = 1e-12
    end

    screw = run_planar_model(joinpath(MODEL_DIRECTORY,
        "screw-motion-coupler.toml"); duration = 0.5, samples = 6)
    screw_coupler = screw.loaded.connections[:screw_coupler]
    screw_bearing = screw.loaded.connections[:screw_bearing]
    nut_travel = screw.loaded.connections[:nut_travel]
    @test screw_coupler.coordinate_kind == :mixed
    @test screw.loaded.analysis.degrees_of_freedom == 1
    @test screw.loaded.state_selection.preferred_velocities ==
        [Symbol("screw_bearing.omega")]
    for state in screw.states
        angle = state[screw_bearing.rotation_variables[3]]
        angular_velocity = state[screw_bearing.rotation_variables[2]]
        @test state[nut_travel.distance_variable] - 0.05angle ≈
            screw_coupler.offset atol = 2e-12
        @test state[nut_travel.velocity_variable] - 0.05angular_velocity ≈
            0.0 atol = 2e-12
        @test state[screw_coupler.reaction_variable] ≈
            10 / 81 atol = 2e-12
    end

    rotational_source = read(rotational_path, String)
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        rotational_source,
        "coefficients = [1.0, -2.0]" => "coefficients = [1.0]")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        rotational_source,
        "output_bearing.rotation" => "output_bearing.distance")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        rotational_source,
        "offset = \"initial\"" => "offset = \"unsupported\"")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        rotational_source,
        "rotation_coordinates = true" => "rotation_coordinates = false";
        count = 1)))
end

@testset "Runtime state reselection" begin
    source = read(joinpath(MODEL_DIRECTORY,
        "relative-coordinate-pendulum.toml"), String)
    source = replace(source,
        "position = [0.25, -0.4330127018922193]" =>
            "position = [0.0, -0.5]",
        "orientation = 0.5235987755982988" => "orientation = 0.0",
        "velocity = [0.5, 0.0]" =>
            "velocity = [3.132091952673165, 0.0]",
        "angular_velocity = 1.0" =>
            "angular_velocity = 6.26418390534633",
        "[pendulum.pin]" => """
            [pendulum.initial.impose]
            V_x = 3.132091952673165

            [pendulum.pin]""",
        "preferred_velocities = [\"pin.omega\"]" =>
            "preferred_velocities = [\"pendulum.V_x\"]")
    result = run_planar_model(IOBuffer(source); duration = 0.7, samples = 71)
    joint = result.loaded.connections[:pin]
    body = result.loaded.bodies[:pendulum]

    @test result.solution.stats.state_reselections == 1
    @test result.solution.stats.history_restarts == 0
    @test length(result.state_selection_changes) == 2
    change = last(result.state_selection_changes)
    @test change.reason in
        (:high_physical_error, :deteriorating_iteration_matrix)
    @test change.previous == [Symbol("pendulum.V_x")]
    @test change.selected == [Symbol("pin.omega")]
    @test 0.0 < change.time < 0.5
    @test length(body.candidate_state_equations) == 3
    @test length(joint.candidate_state_equations) == 2
    active = result.loaded.active_variable_indices
    local_index(index) = only(findall(==(index), active))
    @test result.solution.differential_vars[
        local_index(joint.rotation_variables[2])]
    @test result.solution.differential_vars[
        local_index(joint.rotation_variables[3])]
    @test !result.solution.differential_vars[
        local_index(body.velocity_variables[1])]

    mktempdir() do directory
        path = joinpath(directory, "reselected.simp")
        write_result(path, result)
        stored = read_result(path)
        @test stored.solver_statistics.state_reselections == 1
        @test length(stored.state_selection_changes) == 2
        @test last(stored.state_selection_changes).previous ==
            ["pendulum.V_x"]
        @test last(stored.state_selection_changes).selected == ["pin.omega"]
        @test length(stored.health_snapshots) ==
            length(result.health_step_peaks)
    end
end

@testset "Planar gear pair" begin
    path = joinpath(MODEL_DIRECTORY, "constant-speed-gear-pair.toml")
    result = run_planar_model(path; duration = 0.1, samples = 3)
    loaded = result.loaded
    gear_1, gear_2 = loaded.bodies[:gear1], loaded.bodies[:gear2]
    gear_pair = loaded.connections[:gear_pair]
    applied_forces = PracticalMechanicalSimulation.PlanarAppliedForces

    @test result.analysis_mode == :kinematic
    @test loaded.analysis.degrees_of_freedom == 0
    @test gear_pair.radius_1 == 0.4
    @test gear_pair.radius_2 == -0.6
    @test gear_pair.joint_1 === loaded.connections[:left_bearing]
    @test gear_pair.joint_2 === loaded.connections[:right_bearing]
    @test gear_pair.force_reference_marker ===
        loaded.connections[:left_bearing].marker_b
    @test gear_pair.contact_marker ===
        loaded.markers[Symbol("carrier.contact")].point
    @test gear_pair.marker_1 ===
        loaded.markers[Symbol("gear_pair.contact_1")].point
    @test gear_pair.marker_2 ===
        loaded.markers[Symbol("gear_pair.contact_2")].point
    @test loaded.markers[Symbol("gear_pair.contact_1")].owner === gear_1
    @test loaded.markers[Symbol("gear_pair.contact_2")].owner === gear_2
    drive = loaded.drivers[:drive]
    left_bearing = loaded.connections[:left_bearing]
    right_bearing = loaded.connections[:right_bearing]
    for state in result.states
        @test gear_pair.radius_1 * state[gear_1.orientation_variable] -
              gear_pair.radius_2 * state[gear_2.orientation_variable] ≈
            0.0 atol = 1e-12
        @test state[gear_2.angular_velocity_variable] ≈
            -(2 / 3) * state[gear_1.angular_velocity_variable] atol = 1e-12
        contact = applied_forces.point_marker_kinematics(
            loaded.markers[Symbol("gear_pair.contact_1")].point, state).position
        carrier = applied_forces.point_marker_kinematics(
            loaded.markers[Symbol("carrier.contact")].point, state).position
        @test contact ≈ carrier atol = 1e-12
        @test state[gear_pair.reaction_variable] ≈ -5 / 0.6 atol = 1e-12
        @test state[drive.torque_variable] ≈ 5 * 0.4 / 0.6 atol = 1e-12
        @test state[left_bearing.reaction_variables] ≈ [0.0, 5 / 0.6] atol = 1e-12
        @test state[right_bearing.reaction_variables] ≈ [0.0, -5 / 0.6] atol = 1e-12
        @test state[drive.torque_variable] *
              state[gear_1.angular_velocity_variable] +
              5 * state[gear_2.angular_velocity_variable] ≈ 0.0 atol = 1e-12
    end

    trial = copy(first(result.states))
    trial[gear_pair.reaction_variable] = 1.0
    equations = zeros(length(loaded.layout.catalog.equations))
    PracticalMechanicalSimulation.PlanarComponentAssembly.add_gear_reaction!(
        equations, gear_pair, trial)
    @test equations[gear_1.balance_equations[1:2]] == [0.0, -1.0]
    @test equations[gear_2.balance_equations[1:2]] == [0.0, 1.0]
    @test equations[gear_1.balance_equations[3]] ≈ -0.4
    @test equations[gear_2.balance_equations[3]] ≈ -0.6

    source = read(path, String)
    initial_phase_source = replace(source, "phase = \"initial\"\n" => "")
    initial_phase_source = replace(initial_phase_source,
        "[gear2]\ntype = \"rigid_body\"" =>
        "[gear2]\ntype = \"rigid_body\"\nangle = 0.25")
    initial_phase_model = load_planar_model(IOBuffer(initial_phase_source))
    initial_phase_pair = initial_phase_model.connections[:gear_pair]
    @test initial_phase_pair.phase ≈ 0.15 atol = 1e-15
    @test PracticalMechanicalSimulation.PlanarComponentAssembly.gear_position(
        initial_phase_pair, initial_phase_model.initial_values) ≈ 0.0 atol = 1e-15
    numeric_phase_model = load_planar_model(IOBuffer(replace(source,
        "phase = \"initial\"" => "phase = 0.125")))
    @test numeric_phase_model.connections[:gear_pair].phase == 0.125
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "phase = \"initial\"" => "phase = \"unknown\"")))

    oriented_source = replace(source,
        "[carrier.contact]\ntype = \"marker\"\nposition = [0.4, 0.0]" =>
        "[carrier.contact]\ntype = \"marker\"\nposition = [0.4, 0.0]\nangle = 1.5707963267948966")
    oriented = load_planar_model(IOBuffer(oriented_source))
    oriented_gear_pair = oriented.connections[:gear_pair]
    oriented_trial = copy(oriented.initial_values)
    oriented_trial[oriented_gear_pair.reaction_variable] = 1.0
    oriented_equations = zeros(length(oriented.layout.catalog.equations))
    PracticalMechanicalSimulation.PlanarComponentAssembly.add_gear_reaction!(
        oriented_equations, oriented_gear_pair, oriented_trial)
    @test oriented_equations[
        oriented.bodies[:gear1].balance_equations[1:2]] ≈ [0.0, -1.0] atol = 1e-15
    @test oriented_equations[
        oriented.bodies[:gear2].balance_equations[1:2]] ≈ [0.0, 1.0] atol = 1e-15
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "joints = [\"left_bearing\", \"right_bearing\"]\n" => "")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "contact_marker = \"carrier.contact\"\n" => "")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "contact_marker = \"carrier.contact\"" =>
        "contact_marker = \"gear1.rim\"")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "position = [0.4, 0.0]" => "position = [0.0, 0.0]";
        count = 1)))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "position = [0.4, 0.0]" => "position = [0.4, 0.1]";
        count = 1)))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "markers = [\"gear1.center\", \"carrier.left_bearing\"]" =>
        "markers = [\"carrier.left_bearing\", \"gear1.center\"]";
        count = 1)))

    internal_path = joinpath(MODEL_DIRECTORY,
        "constant-speed-internal-gear-pair.toml")
    internal = run_planar_model(internal_path; duration = 0.1, samples = 3)
    internal_loaded = internal.loaded
    internal_pair = internal_loaded.connections[:gear_pair]
    ring, pinion = internal_loaded.bodies[:ring], internal_loaded.bodies[:pinion]
    internal_drive = internal_loaded.drivers[:drive]
    @test internal.analysis_mode == :kinematic
    @test internal_loaded.analysis.degrees_of_freedom == 0
    @test internal_pair.radius_1 ≈ 1.4 atol = 1e-12
    @test internal_pair.radius_2 ≈ 0.4 atol = 1e-12
    @test internal_pair.force_reference_marker ===
        internal_loaded.connections[:ring_bearing].marker_b
    @test internal_loaded.markers[Symbol("gear_pair.contact_1")].owner === ring
    @test internal_loaded.markers[Symbol("gear_pair.contact_2")].owner === pinion
    for state in internal.states
        @test internal_pair.radius_1 * state[ring.orientation_variable] -
              internal_pair.radius_2 * state[pinion.orientation_variable] ≈
            0.0 atol = 1e-12
        @test state[pinion.angular_velocity_variable] ≈
            3.5 * state[ring.angular_velocity_variable] atol = 1e-12
        @test state[internal_pair.reaction_variable] ≈ 5 / 0.4 atol = 1e-12
        @test state[internal_drive.torque_variable] ≈ -5 * 1.4 / 0.4 atol = 1e-12
        @test state[internal_drive.torque_variable] *
              state[ring.angular_velocity_variable] +
              5 * state[pinion.angular_velocity_variable] ≈ 0.0 atol = 1e-12
    end
    internal_trial = copy(first(internal.states))
    internal_trial[internal_pair.reaction_variable] = 1.0
    internal_equations = zeros(length(internal_loaded.layout.catalog.equations))
    PracticalMechanicalSimulation.PlanarComponentAssembly.add_gear_reaction!(
        internal_equations, internal_pair, internal_trial)
    @test internal_equations[ring.balance_equations[1:2]] == [0.0, -1.0]
    @test internal_equations[pinion.balance_equations[1:2]] == [0.0, 1.0]
    @test internal_equations[ring.balance_equations[3]] ≈ -1.4 atol = 1e-12
    @test internal_equations[pinion.balance_equations[3]] ≈ 0.4 atol = 1e-12
end

@testset "Planetary gear pair carriers" begin
    path = joinpath(MODEL_DIRECTORY,
        "constant-speed-planetary-gear-set.toml")
    result = run_planar_model(path; duration = 0.1, samples = 3)
    loaded = result.loaded
    sun = loaded.bodies[:sun]
    planet = loaded.bodies[:planet]
    ring = loaded.bodies[:ring]
    carrier = loaded.bodies[:carrier]
    sun_planet = loaded.connections[:sun_planet]
    ring_planet = loaded.connections[:ring_planet]
    applied_forces = PracticalMechanicalSimulation.PlanarAppliedForces

    @test result.analysis_mode == :kinematic
    @test loaded.analysis.degrees_of_freedom == 0
    @test sun_planet.carrier_body === carrier
    @test sun_planet.joint_1.body_b === nothing
    @test sun_planet.radius_1 ≈ -0.4 atol = 1e-12
    @test sun_planet.radius_2 ≈ 0.3 atol = 1e-12
    @test sun_planet.force_reference_marker ===
        loaded.connections[:planet_bearing].marker_b
    @test ring_planet.radius_1 ≈ 1.0 atol = 1e-12
    @test ring_planet.radius_2 ≈ 0.3 atol = 1e-12
    @test ring_planet.force_reference_marker ===
        loaded.connections[:planet_bearing].marker_b
    @test loaded.markers[Symbol("sun_planet.contact_1")].owner === sun
    @test loaded.markers[Symbol("sun_planet.contact_2")].owner === planet
    @test loaded.markers[Symbol("ring_planet.contact_1")].owner === ring
    @test loaded.markers[Symbol("ring_planet.contact_2")].owner === planet
    @test PracticalMechanicalSimulation.PlanarComponentAssembly.gear_force_direction(
        sun_planet, first(result.states)) ≈ [0.0, -1.0] atol = 1e-12
    @test PracticalMechanicalSimulation.PlanarComponentAssembly.gear_force_direction(
        ring_planet, first(result.states)) ≈ [0.0, 1.0] atol = 1e-12

    for state in result.states
        sun_speed = state[sun.angular_velocity_variable]
        @test state[carrier.angular_velocity_variable] ≈
            (2 / 7) * sun_speed atol = 1e-11
        @test state[planet.angular_velocity_variable] ≈
            -(2 / 3) * sun_speed atol = 1e-11
        @test state[ring.angular_velocity_variable] ≈ 0.0 atol = 1e-12
        @test state[sun_planet.reaction_variable] ≈ 25 / 7 atol = 1e-8
        @test state[ring_planet.reaction_variable] ≈ -25 / 7 atol = 1e-8
        for pair_name in (:sun_planet, :ring_planet)
            pair = loaded.connections[pair_name]
            contact_position = applied_forces.point_marker_kinematics(
                pair.contact_marker, state).position
            @test applied_forces.point_marker_kinematics(
                pair.marker_1, state).position ≈ contact_position atol = 1e-12
            @test applied_forces.point_marker_kinematics(
                pair.marker_2, state).position ≈ contact_position atol = 1e-12
        end
    end
end

@testset "Planar rack and pinion" begin
    path = joinpath(MODEL_DIRECTORY,
        "constant-speed-rack-and-pinion.toml")
    result = run_planar_model(path; duration = 0.1, samples = 3)
    loaded = result.loaded
    component = loaded.connections[:rack_and_pinion]
    rack = loaded.bodies[:rack]
    pinion = loaded.bodies[:pinion]
    drive = loaded.drivers[:drive]
    assembly = PracticalMechanicalSimulation.PlanarComponentAssembly
    applied_forces = PracticalMechanicalSimulation.PlanarAppliedForces

    @test component isa assembly.PlanarRackAndPinionComponent
    @test component.translational_joint === loaded.connections[:rack_guide]
    @test component.revolute_joint === loaded.connections[:pinion_bearing]
    @test component.pitch_radius == 0.4
    @test component.phase == 0.0
    @test result.analysis_mode == :kinematic
    @test loaded.analysis.degrees_of_freedom == 0
    @test loaded.markers[Symbol("rack_and_pinion.carrier_contact")].point ===
        component.contact_marker
    @test loaded.markers[Symbol("rack_and_pinion.rack_contact")].point ===
        component.rack_marker
    @test loaded.markers[Symbol("rack_and_pinion.pinion_contact")].point ===
        component.pinion_marker
    @test component.rack_marker isa applied_forces.PlanarFloatingPointMarker
    @test component.pinion_marker isa applied_forces.PlanarFloatingPointMarker

    speed = 2pi / 5
    for state in result.states
        rack_motion = assembly.rack_translation_values(component, state)
        @test rack_motion.position +
              component.pitch_radius * state[pinion.orientation_variable] ≈
            0.0 atol = 1e-12
        @test rack_motion.velocity + component.pitch_radius * speed ≈
            0.0 atol = 1e-12
        @test state[component.reaction_variable] ≈ -9.81 atol = 1e-12
        @test state[drive.torque_variable] ≈ 9.81 * 0.4 atol = 1e-12
        carrier_contact = applied_forces.point_marker_kinematics(
            component.contact_marker, state).position
        rack_contact = applied_forces.point_marker_kinematics(
            component.rack_marker, state).position
        pinion_contact = applied_forces.point_marker_kinematics(
            component.pinion_marker, state).position
        @test rack_contact ≈ carrier_contact atol = 1e-12
        @test pinion_contact ≈ carrier_contact atol = 1e-12
        @test state[component.translational_joint.perp.reaction_variable] ≈
            0.981 atol = 1e-12
    end

    trial = copy(first(result.states))
    trial[component.reaction_variable] = 1.0
    equations = zeros(length(loaded.layout.catalog.equations))
    assembly.add_rack_and_pinion_reaction!(equations, component, trial)
    @test equations[rack.balance_equations[1:2]] ≈ [0.0, 1.0] atol = 1e-15
    @test equations[rack.balance_equations[3]] ≈ -0.1 atol = 1e-15
    @test equations[pinion.balance_equations[1:2]] ≈ [0.0, -1.0] atol = 1e-15
    @test equations[pinion.balance_equations[3]] ≈ -0.4 atol = 1e-15

    source = read(path, String)
    initial_phase_source = replace(source, "phase = \"initial\"\n" => "")
    initial_phase_source = replace(initial_phase_source,
        "[pinion]\ntype = \"rigid_body\"" =>
        "[pinion]\ntype = \"rigid_body\"\nangle = 0.25")
    initial_phase_model = load_planar_model(IOBuffer(initial_phase_source))
    initial_phase_component =
        initial_phase_model.connections[:rack_and_pinion]
    @test initial_phase_component.phase ≈ 0.1 atol = 1e-15
    @test assembly.rack_and_pinion_position(initial_phase_component,
        initial_phase_model.initial_values) ≈ 0.0 atol = 1e-15
    numeric_phase_model = load_planar_model(IOBuffer(replace(source,
        "phase = \"initial\"" => "phase = 0.125")))
    @test numeric_phase_model.connections[:rack_and_pinion].phase == 0.125
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "phase = \"initial\"" => "phase = \"unknown\"")))

    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "joints = [\"rack_guide\", \"pinion_bearing\"]\n" => "")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "pitch_radius = 0.4" => "pitch_radius = -0.4")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "joints = [\"rack_guide\", \"pinion_bearing\"]" =>
            "joints = [\"pinion_bearing\", \"rack_guide\"]")))
end

@testset "Planar pulley belt" begin
    assembly = PracticalMechanicalSimulation.PlanarComponentAssembly
    for tangent_kind in (1, -1), tangent_side in (1, -1)
        geometry = assembly.belt_tangent_geometry(
            [0.0, 0.0], 0.4, [2.0, 0.0], 0.6,
            tangent_kind, tangent_side)
        @test norm(geometry.normal_1) ≈ 1.0 atol = 1.0e-14
        @test geometry.normal_2 ≈
            tangent_kind .* geometry.normal_1 atol = 1.0e-14
        @test dot(geometry.normal_1, geometry.tangent) ≈ 0.0 atol = 1.0e-14
        @test dot(geometry.normal_2, geometry.tangent) ≈ 0.0 atol = 1.0e-14
        @test geometry.point_1 ≈ 0.4 .* geometry.normal_1 atol = 1.0e-14
        @test geometry.point_2 ≈
            [2.0, 0.0] + 0.6 .* geometry.normal_2 atol = 1.0e-14
    end
    @test_throws DomainError assembly.belt_tangent_geometry(
        [0.0, 0.0], 1.0, [1.5, 0.0], 1.0, -1, 1)

    path = joinpath(MODEL_DIRECTORY, "three-pulley-belt-tensioner.toml")
    loaded = load_planar_model(path)
    belt = loaded.connections[:belt]
    @test belt isa assembly.PlanarBeltComponent
    @test length(belt.spans) == 3
    @test all(span -> span isa assembly.PlanarBeltSpanComponent, belt.spans)
    @test all(span -> span.tangent_kind == 1, belt.spans)
    @test all(span -> span.tangent_side == 1, belt.spans)
    @test all(span -> loaded.initial_values[span.tension_variable] ≈ 100.0,
        belt.spans)
    @test loaded.analysis.degrees_of_freedom == 3

    result = run_planar_model(path; duration = 0.05, samples = 4)
    @test result.analysis_mode == :dynamic
    @test result.static_initialization_iterations > 0
    @test result.solution.retcode == SciMLBase.ReturnCode.Success
    initial_tensions = [first(result.states)[span.tension_variable]
        for span in belt.spans]
    @test initial_tensions[1] - initial_tensions[2] ≈ 2.0 atol = 2.0e-7
    @test initial_tensions[2] - initial_tensions[3] ≈ 0.0 atol = 2.0e-7
    for state in result.states, span in belt.spans
        values = assembly.belt_span_values(span, state)
        @test values.point_1 ≈ values.actual_point_1 atol = 2.0e-7
        @test values.point_2 ≈ values.actual_point_2 atol = 2.0e-7
        @test values.tangent ≈ values.actual_tangent atol = 2.0e-7
        @test values.length ≈ values.actual_length atol = 2.0e-7
        @test values.extension ≈ values.expected_extension atol = 2.0e-7
        @test values.extension_rate ≈
            values.expected_extension_rate atol = 2.0e-7
        @test values.tension ≈ span.stiffness * values.extension +
            span.damping * values.extension_rate atol = 2.0e-5
        @test values.force ≈ values.tension .* values.tangent atol = 2.0e-5
        @test values.tension > 0
    end

    source = read(path, String)
    free_length_source = replace(source, "initial_tension = 100.0" =>
        "free_length = $(belt.free_length)")
    free_length_model = load_planar_model(IOBuffer(free_length_source))
    @test free_length_model.connections[:belt].initial_tension ≈
        100.0 atol = 1.0e-8
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "initial_tension = 100.0" =>
        "initial_tension = 100.0\nfree_length = 3.0")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "near_points = [\"carrier.driver_upper_hint\", \"carrier.driven_upper_hint\"]" =>
        "near_points = [\"carrier.driver_center\", \"carrier.driven_center\"]")))
end

@testset "Perp constraint primitive" begin
    path = joinpath(MODEL_DIRECTORY, "perp-guided-slider.toml")
    result = run_planar_model(path; duration = 0.1, samples = 3)
    loaded = result.loaded
    body = loaded.bodies[:slider]
    perp = loaded.connections[:guide_orientation]
    assembly = PracticalMechanicalSimulation.PlanarComponentAssembly

    @test result.analysis_mode == :dynamic
    @test loaded.analysis.degrees_of_freedom == 1
    @test loaded.state_selection.selected_velocities == [Symbol("slider.V_x")]
    @test perp isa assembly.PlanarPerpConstraint
    for state in result.states
        @test state[body.position_variables[2]] ≈ 0.0 atol = 1e-12
        @test state[body.orientation_variable] ≈ 0.0 atol = 1e-12
        @test state[body.angular_velocity_variable] ≈ 0.0 atol = 1e-12
        @test assembly.perp_position(perp, state) ≈ 0.0 atol = 1e-12
        @test assembly.perp_velocity(perp, state) ≈ 0.0 atol = 1e-12
        @test state[perp.reaction_variable] ≈ -1.0 atol = 1e-12
        @test state[loaded.connections[:guide_position].reaction_variable] ≈
            9.81 atol = 1e-12
    end
    final = last(result.states)
    @test final[body.position_variables[1]] ≈ -0.4 atol = 1e-10
    @test final[body.velocity_variables[1]] ≈ 1.0 atol = 1e-12

    trial = copy(first(result.states))
    trial[perp.reaction_variable] = 2.0
    equations = zeros(length(loaded.layout.catalog.equations))
    assembly.add_perp_reaction!(equations, perp, trial)
    @test equations[body.balance_equations[3]] == -2.0

    source = read(path, String)
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "[guide_orientation]\ntype = \"perp\"\nmarkers = [\"slider.center\", \"ground.guide\"]" =>
            "[guide_orientation]\ntype = \"perp\"\nmarkers = [\"slider.center\"]")))
end

@testset "Translational joint composition" begin
    path = joinpath(MODEL_DIRECTORY, "translational-joint-slider.toml")
    result = run_planar_model(path; duration = 0.1, samples = 3)
    loaded = result.loaded
    joint = loaded.connections[:guide]
    assembly = PracticalMechanicalSimulation.PlanarComponentAssembly
    body = loaded.bodies[:slider]

    @test joint isa assembly.PlanarTranslationalJoint
    @test joint.inplane isa assembly.PlanarInplaneConstraint
    @test joint.perp isa assembly.PlanarPerpConstraint
    @test joint.inplane.geometry.axis.orientation === joint.perp.orientation_j
    @test result.analysis_mode == :dynamic
    @test loaded.analysis.degrees_of_freedom == 1
    @test loaded.state_selection.selected_velocities == [Symbol("slider.V_x")]
    components = Set(variable.component for variable in
        loaded.layout.catalog.variables)
    @test Symbol("guide.inplane") in components
    @test Symbol("guide.perp") in components

    for state in result.states
        @test state[body.position_variables[2]] ≈ 0.0 atol = 1e-12
        @test state[body.orientation_variable] ≈ 0.0 atol = 1e-12
        @test assembly.inplane_position(joint.inplane, state) ≈ 0.0 atol = 1e-12
        @test assembly.perp_position(joint.perp, state) ≈ 0.0 atol = 1e-12
        @test state[joint.inplane.reaction_variable] ≈ 9.81 atol = 1e-12
        @test state[joint.perp.reaction_variable] ≈ -1.0 atol = 1e-12
    end
    final = last(result.states)
    @test final[body.position_variables[1]] ≈ -0.4 atol = 1e-10
    @test final[body.velocity_variables[1]] ≈ 1.0 atol = 1e-12

    source = read(path, String)
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "markers = [\"slider.center\", \"ground.guide\"]" =>
            "markers = [\"slider.center\"]")))
end

@testset "Fixed joint composition" begin
    path = joinpath(MODEL_DIRECTORY,
        "fixed-joint-rotating-assembly.toml")
    result = run_planar_model(path; duration = 0.1, samples = 3)
    loaded = result.loaded
    fixed = loaded.connections[:rigid_connection]
    assembly = PracticalMechanicalSimulation.PlanarComponentAssembly
    carrier = loaded.bodies[:carrier]
    extension = loaded.bodies[:extension]
    applied_forces = PracticalMechanicalSimulation.PlanarAppliedForces

    @test fixed isa assembly.PlanarFixedJoint
    @test fixed.revolute isa assembly.PlanarRevoluteJointComponent
    @test fixed.perp isa assembly.PlanarPerpConstraint
    @test result.analysis_mode == :kinematic
    @test loaded.analysis.degrees_of_freedom == 0
    components = Set(variable.component for variable in
        loaded.layout.catalog.variables)
    @test Symbol("rigid_connection.revolute") in components
    @test Symbol("rigid_connection.perp") in components

    for state in result.states
        point_i = applied_forces.point_marker_kinematics(
            fixed.revolute.marker_a, state).position
        point_j = applied_forces.point_marker_kinematics(
            fixed.revolute.marker_b, state).position
        @test point_i ≈ point_j atol = 1e-11
        @test state[carrier.orientation_variable] ≈
            state[extension.orientation_variable] atol = 1e-11
        @test state[carrier.angular_velocity_variable] ≈
            state[extension.angular_velocity_variable] atol = 1e-11
    end

    initial = first(result.states)
    speed = 2pi / 5
    @test initial[fixed.revolute.reaction_variables] ≈
        [1.5speed^2, -9.81] atol = 1e-11
    @test initial[fixed.perp.reaction_variable] ≈ -4.905 atol = 1e-11

    source = read(path, String)
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "[rigid_connection]\ntype = \"fixed\"\nmarkers = [\"carrier.tip\", \"extension.base\"]" =>
            "[rigid_connection]\ntype = \"fixed\"\nmarkers = [\"carrier.tip\"]")))
end

@testset "TOML spanning force" begin
    path = joinpath(MODEL_DIRECTORY, "spanning-spring-pendulum.toml")
    result = run_planar_model(path; duration = 0.1, samples = 3)
    loaded = result.loaded
    spring = only(loaded.forces[:spring])
    element = spring.element
    applied_forces = PracticalMechanicalSimulation.PlanarAppliedForces
    assembly = PracticalMechanicalSimulation.PlanarComponentAssembly
    analysis = PracticalMechanicalSimulation.AutomaticAnalysis

    @test spring isa assembly.PlanarSpanningForceComponent
    @test loaded.analysis.degrees_of_freedom == 1
    @test loaded.state_selection.selected_velocities ==
        [Symbol("pendulum.omega")]
    @test length(analysis.component_variable_indices(
        loaded.layout, :spring)) == 9
    variable_names = Set(Symbol(variable.component, :., variable.name)
        for variable in loaded.layout.catalog.variables)
    @test Symbol("spring.length") in variable_names
    @test Symbol("spring.length_rate") in variable_names
    @test Symbol("spring.force") in variable_names
    @test Symbol("spring.ell") ∉ variable_names
    @test spring.law.dependencies ==
        [element.length_variable, element.length_rate_variable]
    @test spring.law.gradient(0.0, first(result.states)) == [-20.0, -0.5]
    for state in result.states
        marker_1 = applied_forces.point_marker_kinematics(
            element.marker_1, state)
        marker_2 = applied_forces.point_marker_kinematics(
            element.marker_2, state)
        spanning = marker_2.position - marker_1.position
        relative_velocity = marker_2.velocity - marker_1.velocity
        @test state[element.spanning_variables] ≈ spanning atol = 1e-7
        @test state[element.length_variable] ≈ norm(spanning) atol = 1e-7
        @test state[element.unit_variables] ≈
            spanning / norm(spanning) atol = 1e-7
        @test state[element.length_rate_variable] ≈
            dot(state[element.unit_variables], relative_velocity) atol = 1e-7
        @test state[element.force_variable] ≈
            -element.stiffness * (state[element.length_variable] -
                element.free_length) -
            element.damping * state[element.length_rate_variable] atol = 1e-7
        @test state[element.global_force_variables] ≈
            -state[element.unit_variables] * state[element.force_variable] atol = 1e-7
    end

    source = read(path, String)
    expression_source = replace(source, """
        stiffness = 20.0
        damping = 0.5
        free_length = 0.5
        """ => """
        expression = "-20.0*(spring.length-0.5)-0.5*spring.length_rate"
        """)
    expression_result = run_planar_model(IOBuffer(expression_source);
        duration = 0.1, samples = 3)
    expression_force = only(expression_result.loaded.forces[:spring])
    @test expression_force.law.gradient(0.0,
        first(expression_result.states)) == [-20.0, -0.5]
    for (linear_state, expression_state) in
            zip(result.states, expression_result.states)
        @test expression_state ≈ linear_state rtol = 2.0e-7
    end

    constant_source = replace(source, """
        stiffness = 20.0
        damping = 0.5
        free_length = 0.5
        """ => "force = 1.25\n")
    constant_result = run_planar_model(IOBuffer(constant_source);
        duration = 0.0, samples = 1)
    constant_force = only(constant_result.loaded.forces[:spring])
    @test only(constant_result.states)[
        constant_force.element.force_variable] ≈ 1.25 atol = 1.0e-12

    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "type = \"spanning_force\"" =>
            "type = \"spanning_spring_damper\"")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "stiffness = 20.0" => "force = 1.0\nstiffness = 20.0")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "stiffness = 20.0" => "stiffness = -1.0")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "damping = 0.5\n" => "")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "position = [0.8, -0.2]" => "position = [0.0, -1.0]")))
end

@testset "TOML torsional spring-damper" begin
    path = joinpath(MODEL_DIRECTORY, "torsional-spring-pendulum.toml")
    result = run_planar_model(path; duration = 0.1, samples = 3)
    loaded = result.loaded
    spring = only(loaded.forces[:pin_spring])
    element = spring.element
    body = loaded.bodies[:pendulum]
    assembly = PracticalMechanicalSimulation.PlanarComponentAssembly
    applied_forces = PracticalMechanicalSimulation.PlanarAppliedForces
    analysis = PracticalMechanicalSimulation.AutomaticAnalysis

    @test spring isa assembly.PlanarTorsionalSpringComponent
    @test result.analysis_mode == :dynamic
    @test loaded.analysis.degrees_of_freedom == 1
    @test loaded.state_selection.selected_velocities == [Symbol("pin.omega")]
    @test length(analysis.component_variable_indices(
        loaded.layout, :pin_spring)) == 1
    @test element.stiffness == 3.0
    @test spring.damping_time_scale == 0.1
    @test element.damping ≈ 0.3 atol = 1e-15
    @test element.free_angle == 0.0
    for state in result.states
        relative = applied_forces.relative_rotation(element, state)
        @test state[element.torque_variable] ≈
            -element.stiffness * (relative.angle - element.free_angle) -
            element.damping * relative.angular_velocity atol = 1e-7
    end

    trial = copy(first(result.states))
    trial[element.torque_variable] = 1.0
    equations = zeros(length(loaded.layout.catalog.equations))
    contribution = only(assembly.equation_contributions(spring))
    contribution.residual!(equations, 0.0, trial, zero(trial))
    @test equations[body.balance_equations[3]] == -1.0

    source = read(path, String)
    degree_source = replace(source,
        "orientation = 0.7853981633974483" => "orientation = \"45°\"",
        "free_angle = 0.0" => "free_angle = \"10 deg\"")
    degree_loaded = load_planar_model(IOBuffer(degree_source))
    @test degree_loaded.initial_values[
        degree_loaded.bodies[:pendulum].orientation_variable] ≈ pi / 4
    @test only(degree_loaded.forces[:pin_spring]).element.free_angle ≈ pi / 18
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        degree_source, "\"45°\"" => "\"forty°\"")))
    explicit = load_planar_model(IOBuffer(replace(source,
        "damping_time_scale = 0.1" =>
        "damping_time_scale = 0.1\ndamping = 0.7")))
    @test only(explicit.forces[:pin_spring]).element.damping == 0.7
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "stiffness = 3.0\n" => "")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "damping_time_scale = 0.1" => "damping_time_scale = -0.1")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "markers = [\"pendulum.pin\", \"ground.origin\"]" =>
        "markers = [\"pendulum.pin\"]"; count = 1)))
end

@testset "Revolute relative coordinate" begin
    model_path = joinpath(MODEL_DIRECTORY,
        "relative-coordinate-pendulum.toml")
    result = run_planar_model(model_path; duration = 0.05, samples = 3)
    loaded = result.loaded
    joint = loaded.connections[:pin]
    body = loaded.bodies[:pendulum]
    alpha, omega, theta = joint.rotation_variables

    @test loaded.analysis.degrees_of_freedom == 1
    @test loaded.state_selection.selected_velocities == [Symbol("pin.omega")]
    @test Symbol("pin.omega") in loaded.state_selection.qr_pivots
    @test length(joint.rotation_equations) == 3
    @test length(joint.state_equations) == 2
    @test joint.selected_state_variables == [(alpha, omega, theta)]
    @test count(result.solution.differential_vars) == 2
    @test result.static_initialization_iterations > 0
    @test first(result.states)[body.position_variables] ≈ [0.0, -0.5] atol = 1.0e-11
    @test first(result.states)[body.velocity_variables[1]] ≈
        0.5first(result.states)[body.angular_velocity_variable] atol = 1.0e-11
    @test first(result.states)[body.velocity_variables[2]] ≈ 0.0 atol = 1.0e-11
    @test first(result.states)[body.orientation_variable] ≈ 0.0 atol = 1.0e-11
    @test 0.9 < first(result.states)[body.angular_velocity_variable] < 1.0
    @test maximum(abs(state[body.orientation_variable]) for state in result.states) >
        0.01
    for state in result.states
        @test state[theta] ≈ state[body.orientation_variable] atol = 1.0e-11
        @test state[omega] ≈ state[body.angular_velocity_variable] atol = 1.0e-11
        @test state[alpha] ≈ state[body.angular_acceleration_variable] atol = 1.0e-10
    end

    default_joint = load_planar_model(joinpath(MODEL_DIRECTORY,
        "constant-speed-slider-crank.toml")).connections[:ground_pin]
    @test isempty(default_joint.rotation_variables)

    source = read(model_path, String)

    imposed_source = replace(source, "[gravity]" => """
        [pin.initial.impose]
        angle = "5 deg"
        omega = 10.0

        [gravity]""")
    imposed = load_planar_model(IOBuffer(imposed_source))
    imposed_joint = imposed.connections[:pin]
    imposed_body = imposed.bodies[:pendulum]
    imposed_state = imposed.initial_values
    _, imposed_omega, imposed_theta = imposed_joint.rotation_variables
    @test imposed_state[imposed_theta] ≈ deg2rad(5) atol = 1.0e-13
    @test imposed_state[imposed_omega] ≈ 10.0 atol = 1.0e-12
    @test imposed_state[imposed_body.orientation_variable] ≈
        deg2rad(5) atol = 1.0e-13
    @test imposed_state[imposed_body.angular_velocity_variable] ≈
        10.0 atol = 1.0e-12
    @test imposed_state[imposed_body.position_variables] ≈
        [0.5sin(deg2rad(5)), -0.5cos(deg2rad(5))] atol = 1.0e-12
    @test imposed_state[imposed_body.velocity_variables] ≈
        [5cos(deg2rad(5)), 5sin(deg2rad(5))] atol = 1.0e-11
    @test imposed.state_selection.imposed_initial_variables ==
        [Symbol("pin.omega"), Symbol("pin.theta")]

    body_imposed_source = replace(source,
        "[pendulum.pin]" => """
        [pendulum.initial.impose]
        V_x = 0.75

        [pendulum.pin]""")
    body_imposed = load_planar_model(IOBuffer(body_imposed_source))
    body_imposed_body = body_imposed.bodies[:pendulum]
    @test body_imposed.initial_values[
        body_imposed_body.velocity_variables[1]] == 0.75
    @test body_imposed.state_selection.imposed_initial_variables ==
        [Symbol("pendulum.V_x")]

    low_weight_source = replace(source, "[gravity]" => """
        [pin.initial]
        omega = 10.0
        omega_weight = 1.0

        [gravity]""")
    high_weight_source = replace(low_weight_source,
        "omega_weight = 1.0" => "omega_weight = 1000.0")
    low_weight = load_planar_model(IOBuffer(low_weight_source))
    high_weight = load_planar_model(IOBuffer(high_weight_source))
    low_joint = low_weight.connections[:pin]
    high_joint = high_weight.connections[:pin]
    low_omega = low_weight.initial_values[low_joint.rotation_variables[2]]
    high_omega = high_weight.initial_values[high_joint.rotation_variables[2]]
    @test abs(high_omega - 10.0) < abs(low_omega - 10.0)

    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "angular_velocity = 1.0" =>
            "angular_velocity = 1.0\norientation_ic_weight = 2.0")))

    two_body_source = """
        [model]
        dimension = "planar"

        [light]
        type = "rigid_body"
        mass = 1.0
        inertia = 0.1
        position = [0.0, 0.0]

        [light.point]
        type = "marker"

        [heavy]
        type = "rigid_body"
        mass = 9.0
        inertia = 0.9
        position = [1.0, 0.0]
        velocity = [1.0, 0.0]

        [heavy.point]
        type = "marker"

        [joint]
        type = "revolute"
        markers = ["light.point", "heavy.point"]
        """
    mass_weighted = load_planar_model(IOBuffer(two_body_source))
    light = mass_weighted.bodies[:light]
    heavy = mass_weighted.bodies[:heavy]
    @test mass_weighted.entered_initial_values[
        light.position_variables[1]] ≈ 0.0
    @test mass_weighted.entered_initial_values[
        heavy.position_variables[1]] ≈ 1.0
    @test mass_weighted.entered_initial_values[
        light.velocity_variables[1]] ≈ 0.0
    @test mass_weighted.entered_initial_values[
        heavy.velocity_variables[1]] ≈ 1.0
    @test mass_weighted.initial_values[light.position_variables[1]] ≈ 0.9
    @test mass_weighted.initial_values[heavy.position_variables[1]] ≈ 0.9
    @test mass_weighted.initial_values[light.velocity_variables[1]] ≈ 0.9
    @test mass_weighted.initial_values[heavy.velocity_variables[1]] ≈ 0.9
    @test mass_weighted.initial_condition_weights[:light].translational_weight == 1.0
    @test mass_weighted.initial_condition_weights[:heavy].translational_weight == 9.0

    equal_weight_source = replace(two_body_source,
        "mass = 1.0\ninertia = 0.1" =>
            "mass = 1.0\ninertia = 0.1\nic_weight_scale = 9.0")
    equally_weighted = load_planar_model(IOBuffer(equal_weight_source))
    equal_light = equally_weighted.bodies[:light]
    equal_heavy = equally_weighted.bodies[:heavy]
    @test equally_weighted.initial_values[
        equal_light.position_variables[1]] ≈ 0.5
    @test equally_weighted.initial_values[
        equal_heavy.position_variables[1]] ≈ 0.5
    @test equally_weighted.initial_values[
        equal_light.velocity_variables[1]] ≈ 0.5
    @test equally_weighted.initial_values[
        equal_heavy.velocity_variables[1]] ≈ 0.5

    massless_source = replace(two_body_source,
        "mass = 1.0\ninertia = 0.1" => "mass = 0.0\ninertia = 0.0")
    massless = load_planar_model(IOBuffer(massless_source))
    @test massless.initial_condition_weights[:light].effective_mass ≈ 9.0e-6
    @test massless.initial_condition_weights[:light].translational_weight > 0
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        two_body_source, "mass = 1.0" =>
            "mass = 1.0\nic_weight_scale = 0.0"; count = 1)))

    distance_source = """
        [model]
        dimension = "planar"

        [body]
        type = "rigid_body"
        mass = 1.0
        inertia = 1.0

        [body.point]
        type = "marker"

        [ground]
        type = "ground"

        [ground.axis]
        type = "marker"

        [travel]
        type = "distance_coordinate"
        markers = ["body.point", "ground.axis"]

        [travel.initial.impose]
        distance = 0.25
        velocity = 1.0
        """
    distance_loaded = load_planar_model(IOBuffer(distance_source))
    distance_body = distance_loaded.bodies[:body]
    distance_coordinate = distance_loaded.connections[:travel]
    @test distance_loaded.initial_values[
        distance_coordinate.distance_variable] ≈ 0.25 atol = 1.0e-13
    @test distance_loaded.initial_values[
        distance_coordinate.velocity_variable] ≈ 1.0 atol = 1.0e-13
    @test distance_loaded.initial_values[
        distance_body.position_variables[2]] ≈ 0.25 atol = 1.0e-13
    @test distance_loaded.initial_values[
        distance_body.velocity_variables[2]] ≈ 1.0 atol = 1.0e-13

    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "rotation_coordinates = true" => "rotation_coordinates = 1")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "rotation_coordinates = true\n" => "")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        imposed_source, "rotation_coordinates = true" =>
            "rotation_coordinates = false")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        imposed_source, "angle = \"5 deg\"" => "angle = \"bad deg\"")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        imposed_source, "omega = 10.0" => "acceleration = 10.0")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        body_imposed_source, "V_x = 0.75" => "speed = 0.75")))
    conflicting_source = replace(source, "[pendulum.pin]" => """
        [pendulum.initial.impose]
        R_x = 1.0
        R_y = 1.0
        theta = 0.0

        [pendulum.pin]""")
    @test_throws ArgumentError load_planar_model(IOBuffer(conflicting_source))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        low_weight_source, "omega_weight = 1.0" =>
            "omega_weight = 0.0")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        imposed_source, "[pin.initial.impose]" =>
            "[pin.initial]\nangle = \"5 deg\"\n\n[pin.initial.impose]")))

    automatic_source = replace(source, """

    [state_selection]
    method = "preferred"
    preferred_velocities = ["pin.omega"]
    allow_fallback = false
    """ => "")
    automatic = load_planar_model(IOBuffer(automatic_source))
    @test automatic.state_selection.selected_velocities ==
        [Symbol("pin.omega")]

    offset_source = replace(source,
        "[ground.origin]\ntype = \"marker\"" =>
            "[ground.origin]\ntype = \"marker\"\nangle = 0.2",
        "[pendulum.pin]\ntype = \"marker\"\nposition = [0.0, 0.5]" =>
            "[pendulum.pin]\ntype = \"marker\"\nposition = [0.0, 0.5]\nangle = -0.1")
    offset = load_planar_model(IOBuffer(offset_source))
    offset_joint = offset.connections[:pin]
    @test offset.initial_values[offset_joint.rotation_variables[3]] ≈
        offset.initial_values[offset.bodies[:pendulum].orientation_variable] - 0.3

    mktempdir() do directory
        path = joinpath(directory, "relative-coordinate.simp")
        write_result(path, result)
        stored = read_result(path)
        names = ["$(stored.variable_components[i]).$(stored.variable_names[i])"
            for i in eachindex(stored.variable_names)]
        @test all(name -> name in names,
            ("pin.theta", "pin.omega", "pin.alpha"))

        transfer_model_path = joinpath(directory, "relative-transfer.toml")
        transfer_source = replace(source,
            "initialization = \"static_equilibrium\"\n" => "",
            "position = [0.25, -0.4330127018922193]" =>
                "position = [0.0, -0.5]",
            "orientation = 0.5235987755982988" => "orientation = 0.0",
            "velocity = [0.5, 0.0]" => "velocity = [0.0, 0.0]",
            "angular_velocity = 1.0" => "angular_velocity = 0.0") * """

            [initial_conditions]
            result = "relative-coordinate.simp"
            include_velocities = true
            """
        write(transfer_model_path, transfer_source)
        transferred = run_planar_model(transfer_model_path;
            duration = 0.0, samples = 1)
        transferred_joint = transferred.loaded.connections[:pin]
        transferred_state = only(transferred.states)
        source_state = last(result.states)
        _, transferred_omega, transferred_theta =
            transferred_joint.rotation_variables
        @test Symbol("pin.theta") in
            transferred.loaded.initial_conditions.transferred_variables
        @test Symbol("pin.omega") in
            transferred.loaded.initial_conditions.transferred_variables
        @test isapprox(transferred_state[transferred_theta],
            source_state[theta]; atol = 1.0e-7)
        @test isapprox(transferred_state[transferred_omega],
            source_state[omega]; atol = 1.0e-7)
    end
end

@testset "Translational distance generator" begin
    model_path = joinpath(MODEL_DIRECTORY,
        "translational-distance-slider.toml")
    result = run_planar_model(model_path; duration = 0.1, samples = 3)
    loaded = result.loaded
    driver = loaded.drivers[:drive]
    body = loaded.bodies[:slider]
    assembly = PracticalMechanicalSimulation.PlanarComponentAssembly
    analysis = PracticalMechanicalSimulation.AutomaticAnalysis

    @test driver isa assembly.PlanarTranslationalMotionGenerator
    @test result.analysis_mode == :kinematic
    @test loaded.analysis.degrees_of_freedom == 0
    @test length(analysis.component_variable_indices(
        loaded.layout, :drive)) == 4
    @test assembly.inplane_direction(driver, first(result.states)).unit ≈
        [0.0, 1.0] atol = 1e-15
    for (time, state) in zip(result.times, result.states)
        @test state[driver.distance_variable] ≈ 0.4time atol = 1e-11
        @test state[driver.velocity_variable] ≈ 0.4 atol = 1e-12
        @test state[driver.acceleration_variable] ≈ 0.0 atol = 1e-12
        @test state[body.position_variables] ≈ [0.0, 0.4time] atol = 1e-11
        @test state[body.velocity_variables] ≈ [0.0, 0.4] atol = 1e-11
        @test state[body.acceleration_variables] ≈ [0.0, 0.0] atol = 1e-10
        @test assembly.inplane_position(driver, state) ≈
            state[driver.distance_variable] atol = 1e-11
        @test assembly.inplane_velocity(driver, state) ≈
            state[driver.velocity_variable] atol = 1e-11
        @test assembly.inplane_acceleration(driver, state) ≈
            state[driver.acceleration_variable] atol = 1e-10
        @test state[driver.reaction_variable] ≈ 9.81 atol = 1e-10
    end

    trial = copy(first(result.states))
    trial[driver.reaction_variable] = 1.0
    equations = zeros(length(loaded.layout.catalog.equations))
    assembly.add_inplane_reaction!(equations, driver, trial)
    @test equations[body.balance_equations[1:2]] ≈ [0.0, -1.0] atol = 1e-15
    @test equations[body.balance_equations[3]] ≈ 0.0 atol = 1e-15

    source = read(model_path, String)
    offset = run_planar_model(IOBuffer(replace(source,
        "initial_distance = 0.0" => "initial_distance = 0.2"));
        duration = 0.0, samples = 1)
    @test only(offset.states)[offset.loaded.bodies[:slider].position_variables] ≈
        [0.0, 0.2] atol = 1e-11

    expression_source = replace(source,
        "function = \"constant_speed\"\ninitial_distance = 0.0\nvelocity = \"drive_speed\"" =>
        "function = \"expression\"\ndistance = \"0.2*t^2\"")
    expression = run_planar_model(IOBuffer(expression_source);
        duration = 0.1, samples = 2)
    expression_driver = expression.loaded.drivers[:drive]
    expression_final = last(expression.states)
    @test expression_final[expression_driver.distance_variable] ≈
        0.002 atol = 1e-10
    @test expression_final[expression_driver.velocity_variable] ≈
        0.04 atol = 1e-10
    @test expression_final[expression_driver.acceleration_variable] ≈
        0.4 atol = 1e-10
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        expression_source, "distance = \"0.2*t^2\"" =>
            "distance = \"0.2*t^2\"\nvelocity = \"0.4*t\"")))

    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "velocity = \"drive_speed\"\n" => "")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "markers = [\"slider.center\", \"ground.drive_axis\"]" =>
        "markers = [\"slider.center\"]")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "function = \"constant_speed\"" => "function = \"unknown\"")))
end

@testset "Zero-mass planar bodies" begin
    # A massless body with no geometric radius still receives a finite,
    # model-level characteristic length for state-selection scaling.
    isolated_source = """
    [model]
    dimension = "planar"

    [body]
    type = "rigid_body"
    mass = 0.0
    inertia = 0.0
    """
    isolated = load_planar_model(IOBuffer(isolated_source))
    @test isolated.state_selection.body_characteristic_lengths[:body] == 1.0
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        isolated_source, "mass = 0.0" => "mass = -1.0")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        isolated_source, "inertia = 0.0" => "inertia = -1.0")))

    # The connecting rod is a kinematic intermediate. Its motion is determined
    # by the two pins, the in-plane constraint, and the crank generator.
    slider_source = read(joinpath(MODEL_DIRECTORY,
        "constant-speed-slider-crank.toml"), String)
    massless_slider_source = replace(slider_source,
        "mass = 2.0\ninertia = 0.16666666666666666" =>
        "mass = 0.0\ninertia = 0.0")
    slider = run_planar_model(IOBuffer(massless_slider_source);
        duration = 0.1, samples = 3)
    @test slider.analysis_mode == :kinematic
    @test slider.loaded.state_selection.selected_velocities == Symbol[]
    @test all(isfinite, reduce(vcat, slider.states))

    # In the supplied dynamic four-bar, making the coupler massless leaves the
    # automatically selected generalized coordinate on the massive crank.
    four_bar_source = read(joinpath(MODEL_DIRECTORY,
        "torque-driven-four-bar.toml"), String)
    massless_four_bar_source = replace(four_bar_source,
        "mass = 1.5\ninertia = 0.24499999999999997" =>
        "mass = 0.0\ninertia = 0.0")
    four_bar = run_planar_model(IOBuffer(massless_four_bar_source);
        duration = 0.02, samples = 2)
    @test four_bar.loaded.state_selection.selected_velocities ==
        [Symbol("crank.omega")]
    @test all(isfinite, reduce(vcat, four_bar.states))

    # With no inertia, bushing torque balance gives the first-order equation
    # c*theta_dot + k*theta = t. This checks a massless mode positioned by force
    # equilibrium rather than by an ideal constraint.
    bushing_source = read(joinpath(MODEL_DIRECTORY,
        "bushing-supported-body.toml"), String)
    massless_bushing_source = replace(bushing_source,
        "mass = 1.0\ninertia = 0.020833333333333332" =>
            "mass = 0.0\ninertia = 0.0",
        "orientation = 0.2" => "orientation = 0.0",
        "acceleration = [0.0, -9.81]" =>
            "acceleration = [0.0, 0.0]") * """

    [turning]
    type = "applied_torque"
    markers = ["body.mount", "ground.mount"]
    expression = "t"
    """
    duration = 0.1
    response = run_planar_model(IOBuffer(massless_bushing_source);
        duration, samples = 3)
    body = response.loaded.bodies[:body]
    final = last(response.states)
    stiffness, damping = 20.0, 2.0
    expected_angle = duration / stiffness -
        damping / stiffness^2 * (1 - exp(-stiffness / damping * duration))
    expected_omega = (1 - exp(-stiffness / damping * duration)) / stiffness
    @test final[body.orientation_variable] ≈ expected_angle rtol = 2.0e-6
    @test final[body.angular_velocity_variable] ≈ expected_omega rtol = 2.0e-6

    # Without inertia or damping, torque balance is purely algebraic:
    # k*theta = t. The full implicit solution drives the selected angle and
    # angular-velocity histories even though they are not inertial states.
    algebraic_source = replace(bushing_source,
        "mass = 1.0\ninertia = 0.020833333333333332" =>
            "mass = 0.0\ninertia = 0.0",
        "orientation = 0.2" =>
            "orientation = 0.0\nangular_velocity = 0.05",
        "damping_time_scale = 0.1" => "damping_time_scale = 0.0",
        "acceleration = [0.0, -9.81]" =>
            "acceleration = [0.0, 0.0]") * """

    [turning]
    type = "applied_torque"
    markers = ["body.mount", "ground.mount"]
    expression = "t"
    """
    algebraic = run_planar_model(IOBuffer(algebraic_source);
        duration = 0.1, samples = 6)
    algebraic_body = algebraic.loaded.bodies[:body]
    for (time, state) in zip(algebraic.times, algebraic.states)
        @test state[algebraic_body.orientation_variable] ≈
            time / 20 atol = 2.0e-12
        @test state[algebraic_body.angular_velocity_variable] ≈
            0.05 atol = 3.0e-12
        @test state[algebraic_body.angular_acceleration_variable] ≈
            0.0 atol = 1.0e-10
    end
end

@testset "Automatic state and redundant-row selection" begin
    source = """
    [model]
    dimension = "planar"

    [simulation]
    end_time = 0.1
    output_samples = 3

    [ground]
    type = "ground"

    [ground.origin]
    type = "marker"

    [body]
    type = "rigid_body"
    mass = 1.0
    inertia = 0.1

    [body.left]
    type = "marker"
    position = [-0.5, 0.0]

    [body.right]
    type = "marker"
    position = [0.5, 0.0]

    [left_guide]
    type = "inplane"
    markers = ["body.left", "ground.origin"]

    [right_guide]
    type = "inplane"
    markers = ["body.right", "ground.origin"]
    """
    loaded = load_planar_model(IOBuffer(source))
    @test loaded.state_selection.selected_velocities == [Symbol("body.V_x")]
    @test loaded.state_selection.qr_rank == 2
    @test loaded.state_selection.body_characteristic_lengths[:body] == 0.5

    redundant = load_planar_model(IOBuffer(source * """

    [duplicate_guide]
    type = "inplane"
    markers = ["body.left", "ground.origin"]
    """))
    @test length(redundant.state_selection.redundant_velocity_equations) == 1
    @test length(redundant.active_variable_indices) ==
        length(redundant.active_equation_indices)

    preferred = load_planar_model(IOBuffer(source * """

    [state_selection]
    method = "preferred"
    preferred_velocities = ["body.V_x"]
    allow_fallback = false
    """))
    @test preferred.state_selection.selected_velocities == [Symbol("body.V_x")]
    @test !preferred.state_selection.fallback_used
end
