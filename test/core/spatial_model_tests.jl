using PracticalMechanicalSimulation
using PracticalMechanicalSimulation.AutomaticAnalysis
using PracticalMechanicalSimulation.SpatialComponentAssembly
using PracticalMechanicalSimulation.SpatialModeling
using PracticalMechanicalSimulation.SpatialDirectedDistances
using PracticalMechanicalSimulation.SpatialSpans
using PracticalMechanicalSimulation.SpatialBelts
using PracticalMechanicalSimulation.SpatialSimulationRunner
using PracticalMechanicalSimulation.SpatialConstraints
using PracticalMechanicalSimulation.SpatialCoordinateCouplers
using PracticalMechanicalSimulation.SpatialGearPairs
using PracticalMechanicalSimulation.SpatialRackAndPinions
using PracticalMechanicalSimulation.SpatialAppliedForces
using PracticalMechanicalSimulation.SpatialBushings
using PracticalMechanicalSimulation.SpatialPlaneContacts
using PracticalMechanicalSimulation.SpatialTires
using PracticalMechanicalSimulation.SpatialMotionGenerators
using LinearAlgebra
using Test
using TOML

const SPATIAL_MODEL_DIRECTORY = normpath(joinpath(
    @__DIR__, "..", "..", "models", "spatial"))

@testset "Spatial belt tangents and drive" begin
    axis_1 = [0.0, 0.0, 1.0]
    axis_2 = [1.0, 0.0, 0.0]
    geometry = SpatialBelts.angled_tangent_geometry(
        [0.0, 0.0, 0.0], axis_1, 0.4,
        [0.4, 1.0, 0.3], axis_2, 0.3, 1, 1, 1.0e-10)
    @test geometry.point_1 ≈ [0.4, 0.0, 0.0]
    @test geometry.point_2 ≈ [0.4, 1.0, 0.0]
    @test geometry.tangent ≈ [0.0, 1.0, 0.0]
    @test geometry.beta_1 ≈ 0.4
    @test geometry.beta_2 ≈ 0.3
    @test_throws DomainError SpatialBelts.angled_tangent_geometry(
        [0.0, 0.0, 0.0], axis_1, 0.4,
        [0.0, 1.0, 0.0], axis_2, 0.3, 1, 1, 1.0e-8)
    @test_throws DomainError SpatialBelts.parallel_tangent_geometry(
        [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.4,
        [1.0, 0.0, 0.1], [0.0, 0.0, 1.0], 0.3, 1, 1, 1.0e-8)

    path = joinpath(SPATIAL_MODEL_DIRECTORY, "three-pulley-belt.toml")
    loaded = load_spatial_model(path)
    @test loaded.simulation.output_precision == :single
    double_precision = load_spatial_model(IOBuffer(replace(
        read(path, String), "output_samples = 401" =>
            "output_samples = 401\noutput_precision = \"double\"")))
    @test double_precision.simulation.output_precision == :double
    @test loaded.connections[:driver_pulley] isa SpatialPulleyComponent
    @test loaded.connections[:belt] isa SpatialBeltComponent
    @test count(force -> force isa SpatialBeltSpanComponent,
        values(loaded.forces)) == 3
    for span in loaded.connections[:belt].spans
        values = spatial_belt_span_values(span, loaded.initial_values)
        @test values.extension ≈ span.initial_extension atol = 1.0e-12
        @test values.tension ≈ span.stiffness * values.extension +
            span.damping * values.extension_rate atol = 1.0e-10
        @test values.geometry.feasibility_error < 1.0e-12
    end
    result = run_spatial_model(path; duration = 0.05, samples = 6)
    @test length(result.states) == 6
    @test result.loaded.analysis.degrees_of_freedom == 2
    @test abs(last(result.states)[
        loaded.connections[:driven_pulley].joint.hinge.rotation_variables[3]]) >
        1.0e-3

    spatial_path = joinpath(
        SPATIAL_MODEL_DIRECTORY, "out-of-plane-four-pulley-belt.toml")
    spatial = load_spatial_model(spatial_path)
    @test count(force -> force isa SpatialBeltSpanComponent,
        values(spatial.forces)) == 4
    @test all(spatial.connections[:belt].spans) do span
        values = spatial_belt_span_values(span, spatial.initial_values)
        abs(dot(values.geometry.axis_1, values.geometry.axis_2)) < 0.9 &&
            values.geometry.feasibility_error < 1.0e-10
    end
end

@testset "Spatial body reference and center-of-mass marker" begin
    source = """
        [model]
        dimension = "spatial"

        [body]
        type = "rigid_body"
        mass = 2.0
        inertia = [1.0, 2.0, 3.0]
        center_of_mass = "body.cm"
        position = [1.0, 2.0, 3.0]
        velocity = [4.0, 5.0, 6.0]
        angular_velocity = [0.0, 0.0, 2.0]
        orientation = ["90 deg", 0.0, 0.0, 1.0]

        [body.reference]
        type = "marker"

        [body.cm]
        type = "marker"
        position = [0.5, 0.0, 0.0]
        orientation = ["90 deg", 1.0, 0.0, 0.0]
        """
    loaded = load_spatial_model(IOBuffer(source))
    body = loaded.bodies[:body]
    state = loaded.initial_values
    reference = loaded.body_reference_frames[:body]
    body_orientation = rotation_matrix(
        @view state[body.euler_parameter_variables])

    @test reference.center_of_mass_marker == Symbol("body.cm")
    @test reference.center_of_mass_position == [0.5, 0.0, 0.0]
    @test state[body.position_variables] ≈ [1.0, 2.5, 3.0]
    @test state[body.velocity_variables] ≈ [3.0, 5.0, 6.0]
    @test body.inertia ≈ reference.center_of_mass_orientation *
        Diagonal([1.0, 2.0, 3.0]) *
        transpose(reference.center_of_mass_orientation)
    @test loaded.markers[Symbol("body.reference")].position_body ==
        [-0.5, 0.0, 0.0]
    @test loaded.markers[Symbol("body.cm")].position_body == zeros(3)
    @test spatial_marker_position(
        loaded.markers[Symbol("body.reference")], state) ≈ [1.0, 2.0, 3.0]
    @test spatial_marker_velocity(
        loaded.markers[Symbol("body.reference")], state) ≈ [4.0, 5.0, 6.0]
    @test spatial_marker_position(
        loaded.markers[Symbol("body.cm")], state) ≈ [1.0, 2.5, 3.0]
    @test spatial_marker_orientation(
        loaded.markers[Symbol("body.cm")], state) ≈
        body_orientation * reference.center_of_mass_orientation

    imposed_source = source * """

        [body.initial.impose]
        orientation = [0.0, 0.0, 0.0, 1.0]
        omega_z = 4.0
        """
    imposed = load_spatial_model(IOBuffer(imposed_source))
    imposed_body = imposed.bodies[:body]
    @test imposed.initial_values[imposed_body.position_variables] ≈
        [1.5, 2.0, 3.0]
    @test imposed.initial_values[imposed_body.velocity_variables] ≈
        [4.0, 7.0, 6.0]

    result = run_spatial_model(IOBuffer(source); duration = 0.01, samples = 2)
    mktempdir() do directory
        path = joinpath(directory, "offset-cm.simp")
        write_result(path, result)
        stored = read_result(path)
        @test length(stored.body_references) == 1
        stored_reference = only(stored.body_references)
        @test stored_reference.name == "body"
        @test stored_reference.center_of_mass_marker == "body.cm"
        @test stored_reference.center_of_mass_position == [0.5, 0.0, 0.0]
        @test stored_reference.center_of_mass_orientation ≈
            reference.center_of_mass_orientation
    end

    missing = replace(source, "center_of_mass = \"body.cm\"" =>
        "center_of_mass = \"body.missing\"")
    @test_throws ArgumentError load_spatial_model(IOBuffer(missing))
    wrong_type = replace(source, "center_of_mass = \"body.cm\"" =>
        "center_of_mass = [0.5, 0.0, 0.0]")
    @test_throws ArgumentError load_spatial_model(IOBuffer(wrong_type))
    wrong_owner = replace(source, "center_of_mass = \"body.cm\"" =>
        "center_of_mass = \"other.cm\"") * """

        [other]
        type = "rigid_body"
        mass = 1.0
        inertia = [1.0, 1.0, 1.0]

        [other.cm]
        type = "marker"
        """
    @test_throws ArgumentError load_spatial_model(IOBuffer(wrong_owner))
end

@testset "Spatial spur rack and pinion" begin
    path = joinpath(SPATIAL_MODEL_DIRECTORY, "spur-rack-and-pinion.toml")
    source = read(path, String)
    loaded = load_spatial_model(path)
    component = loaded.connections[:rack_and_pinion]
    initial = loaded.initial_values

    @test component isa SpatialRackAndPinion
    @test component.pitch_radius ≈ 0.4
    @test component.rolling_radius ≈ 0.4
    @test loaded.analysis.degrees_of_freedom == 1
    @test loaded.state_selection.selected_velocities ==
        [Symbol("pinion_joint.omega")]
    @test loaded.markers[Symbol("rack_and_pinion.carrier_contact")] isa
        SpatialGroundMarker
    @test loaded.markers[Symbol("rack_and_pinion.rack_contact")] isa
        SpatialFloatingMarker
    @test loaded.markers[Symbol("rack_and_pinion.pinion_contact")] isa
        SpatialFloatingMarker
    @test spatial_marker_position(component.contact_marker, initial) ≈
        [0.0, -0.4, 0.0]
    @test spatial_marker_position(component.rack_contact, initial) ≈
        spatial_marker_position(component.contact_marker, initial)
    @test spatial_marker_position(component.pinion_contact, initial) ≈
        spatial_marker_position(component.contact_marker, initial)
    @test spatial_rack_and_pinion_position(component, initial) ≈ 0.0
    @test spatial_rack_and_pinion_velocity(component, initial) ≈ 0.0

    derivative = SpatialSimulationRunner.initial_spatial_derivative(
        initial, loaded)
    selection = AnalysisSelection(Dynamics(), loaded.active_variable_indices,
        loaded.active_equation_indices)
    coefficient = 2.3
    analytical = Matrix(evaluate_analysis_sparse_jacobian(loaded.model,
        selection, 0.0, initial, derivative, coefficient))
    numerical = similar(analytical)
    step = 1.0e-7
    plus = zeros(size(analytical, 1))
    minus = similar(plus)
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
    @test norm(analytical - numerical, Inf) < 5.0e-7

    result = run_spatial_model(path; duration = 0.1, samples = 11)
    @test maximum(abs(spatial_rack_and_pinion_position(component, state))
        for state in result.states) < 1.0e-10
    @test maximum(abs(spatial_rack_and_pinion_velocity(component, state))
        for state in result.states) < 1.0e-10
    @test abs(last(result.states)[component.reaction_variable]) > 0.01

    misaligned = replace(source,
        "orientation = [\"180 deg\", 0.707106781187, 0.0, 0.707106781187]" =>
        "orientation = [\"90 deg\", 1.0, 0.0, 0.0]"; count = 1)
    @test_throws ArgumentError load_spatial_model(IOBuffer(misaligned))
end

@testset "Spatial redundant constraint removal" begin
    source = read(joinpath(
        SPATIAL_MODEL_DIRECTORY, "spherical-pendulum.toml"), String) * """

    [duplicate_pin]
    type = "spherical"
    markers = ["pendulum.pin", "ground.pin"]
    """
    loaded = load_spatial_model(IOBuffer(source))
    redundant_names = loaded.state_selection.redundant_velocity_equations

    @test loaded.analysis.degrees_of_freedom == 3
    @test length(redundant_names) == 3
    @test all(name -> startswith(String(name), "pin.") ||
        startswith(String(name), "duplicate_pin."), redundant_names)
    @test length(loaded.active_variable_indices) ==
        length(loaded.active_equation_indices)

    active_variables = Set(loaded.active_variable_indices)
    active_equations = Set(loaded.active_equation_indices)
    inactive_families = [family
        for family in values(loaded.layout.constraint_families)
        if family.reaction_variable ∉ active_variables]
    @test length(inactive_families) == 3
    @test all(family -> family.position_equation ∉ active_equations &&
        family.velocity_equation ∉ active_equations &&
        family.acceleration_equation ∉ active_equations,
        inactive_families)

    result = run_spatial_model(IOBuffer(source); duration = 0.02, samples = 3)
    @test length(result.times) == 3
    for state in result.states
        @test spatial_marker_position(
            result.loaded.connections[:pin].marker_a, state) ≈
            spatial_marker_position(
            result.loaded.connections[:pin].marker_b, state) atol = 2.0e-8
    end
    result_active_variables = Set(result.loaded.active_variable_indices)
    result_inactive_families = [family
        for family in values(result.loaded.layout.constraint_families)
        if family.reaction_variable ∉ result_active_variables]
    mktempdir() do directory
        path = joinpath(directory, "redundant-spatial.simp")
        write_result(path, result)
        stored = read_result(path)
        @test all(family -> !stored.variable_active[
                family.reaction_variable] &&
            !stored.equation_active[family.position_equation] &&
            !stored.equation_active[family.velocity_equation] &&
            !stored.equation_active[family.acceleration_equation],
            result_inactive_families)
    end

    static_source = replace(source,
        "mode = \"automatic\"" => "mode = \"static\"")
    static_source = replace(static_source,
        "position = [0.5, 0.0, 0.0]" =>
            "position = [0.0, 0.0, -0.5]\n" *
            "orientation = [\"90 deg\", 0.0, 1.0, 0.0]")
    static_result = run_spatial_model(IOBuffer(static_source);
        start_time = 0.0, end_time = 0.0, samples = 1)
    @test static_result.analysis_mode == :static
    @test length(static_result.loaded.active_variable_indices) ==
        length(static_result.loaded.active_equation_indices)

    four_bar_path = joinpath(
        SPATIAL_MODEL_DIRECTORY, "torque-driven-four-bar.toml")
    four_bar = load_spatial_model(four_bar_path)
    @test four_bar.analysis.degrees_of_freedom == 1
    @test length(four_bar.state_selection.redundant_velocity_equations) == 3
    @test length(four_bar.active_variable_indices) ==
        length(four_bar.active_equation_indices)

    four_bar_result = run_spatial_model(four_bar_path;
        duration = 0.02, samples = 3)
    @test length(four_bar_result.times) == 3
    for state in four_bar_result.states
        for joint in values(four_bar_result.loaded.connections)
            @test spatial_marker_position(joint.marker_a, state) ≈
                spatial_marker_position(joint.marker_b, state) atol = 1.0e-8
        end
    end

    modal_source = replace(read(four_bar_path, String),
        "mode = \"automatic\"" => "mode = \"modal\"\nmodes = 1")
    modal_source = replace(modal_source, "torque = 0.5" =>
        "stiffness = 2.0\ndamping_time_scale = 0.02\n" *
        "free_angle = \"initial\"")
    four_bar_modal = run_spatial_model(IOBuffer(modal_source))
    @test four_bar_modal.analysis_mode == :modal
    @test length(four_bar_modal.loaded.state_selection.
        redundant_velocity_equations) == 3
    @test length(four_bar_modal.eigenvalues) == 1
    @test only(four_bar_modal.equation_errors) < 1.0e-12
    @test four_bar_modal.sparse_factorizations == 1
end

@testset "Spatial ideal gear pair" begin
    path = joinpath(SPATIAL_MODEL_DIRECTORY, "bevel-gear-pair.toml")
    source = read(path, String)
    loaded = load_spatial_model(path)
    gear = loaded.connections[:contact]
    initial = loaded.initial_values
    document = TOML.parse(source)

    @test gear isa SpatialGearPair
    @test gear.side_1.radius ≈ 0.4
    @test gear.side_2.radius ≈ -0.2
    @test loaded.analysis.degrees_of_freedom == 1
    @test loaded.state_selection.selected_velocities ==
        [Symbol("gear1_joint.omega")]
    @test loaded.markers[Symbol("contact.contact_1")] isa
        SpatialFloatingMarker
    @test loaded.markers[Symbol("contact.contact_2")] isa
        SpatialFloatingMarker
    @test document["gear1"]["graphics"]["shape"] == "gear"
    @test document["gear1"]["graphics"]["marker"] == "gear1.pitch"
    @test document["gear2"]["graphics"]["shape"] == "gear"
    @test document["gear2"]["graphics"]["marker"] == "gear2.pitch"
    @test loaded.markers[Symbol("gear1.pitch")].body === loaded.bodies[:gear1]
    @test loaded.markers[Symbol("gear2.pitch")].body === loaded.bodies[:gear2]
    @test spatial_marker_position(gear.side_1.joint.marker_b, initial) ≈
        spatial_marker_position(gear.side_2.joint.marker_b, initial)
    @test norm(cross(
        spatial_marker_position(gear.contact_marker, initial) -
            spatial_marker_position(gear.side_1.joint.marker_b, initial),
        spatial_marker_position(gear.contact_marker, initial) -
            spatial_marker_position(gear.side_2.joint.marker_b, initial))) <
        1.0e-14
    @test spatial_gear_position(gear, initial) ≈ 0.0 atol = 1.0e-12
    @test spatial_gear_velocity(gear, initial) ≈ 0.0 atol = 1.0e-12

    derivative = SpatialSimulationRunner.initial_spatial_derivative(
        initial, loaded)
    selection = AnalysisSelection(Dynamics(), loaded.active_variable_indices,
        loaded.active_equation_indices)
    coefficient = 2.1
    analytical = Matrix(evaluate_analysis_sparse_jacobian(loaded.model,
        selection, 0.0, initial, derivative, coefficient))
    numerical = similar(analytical)
    step = 1.0e-7
    plus = zeros(size(analytical, 1))
    minus = similar(plus)
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
    @test norm(analytical - numerical, Inf) < 5.0e-7

    result = run_spatial_model(path; duration = 0.1, samples = 11)
    @test maximum(abs(spatial_gear_position(gear, state))
        for state in result.states) < 2.0e-8
    @test maximum(abs(spatial_gear_velocity(gear, state))
        for state in result.states) < 2.0e-8
    @test abs(last(result.states)[gear.reaction_variable]) > 0.05

    on_axis = replace(source, "position = [0.4, 0.0, 0.2]" =>
        "position = [0.0, 0.0, 0.2]")
    @test_throws ArgumentError load_spatial_model(IOBuffer(on_axis))

    noncollinear = replace(source, "position = [0.4, 0.0, 0.2]" =>
        "position = [0.4, 0.1, 0.2]")
    error = try
        load_spatial_model(IOBuffer(noncollinear))
        nothing
    catch caught
        caught
    end
    @test error isa ArgumentError
    @test occursin("pitch tangents are not collinear", sprint(showerror, error))
    @test occursin("degrees", sprint(showerror, error))
end

@testset "Planetary gear set with moving carrier" begin
    path = joinpath(SPATIAL_MODEL_DIRECTORY, "planetary-gear-set.toml")
    loaded = load_spatial_model(path)
    sun_planet = loaded.connections[:sun_planet]
    planet_ring = loaded.connections[:planet_ring]
    carrier = loaded.bodies[:carrier]
    initial = loaded.initial_values

    @test sun_planet isa SpatialGearPair
    @test planet_ring isa SpatialGearPair
    @test sun_planet.side_1.radius ≈ 0.4
    @test sun_planet.side_2.radius ≈ -0.3
    @test planet_ring.side_1.radius ≈ 1.0
    @test planet_ring.side_2.radius ≈ 0.3
    @test sun_planet.contact_marker.body === carrier
    @test planet_ring.contact_marker.body === carrier
    @test loaded.analysis.degrees_of_freedom == 1
    @test loaded.state_selection.selected_velocities ==
        [Symbol("sun_joint.omega")]
    @test isempty(loaded.state_selection.redundant_velocity_equations)

    result = run_spatial_model(path; duration = 0.2, samples = 11)
    states = result.states
    sun_theta = loaded.connections[:sun_joint].hinge.rotation_variables[3]
    ring_theta = loaded.connections[:ring_joint].hinge.rotation_variables[3]
    carrier_theta =
        loaded.connections[:carrier_joint].hinge.rotation_variables[3]
    planet_phase = sun_planet.position_variables[2]
    initial_planet_phase = initial[planet_phase]

    @test maximum(abs(spatial_gear_position(sun_planet, state))
        for state in states) < 2.0e-10
    @test maximum(abs(spatial_gear_position(planet_ring, state))
        for state in states) < 2.0e-10
    @test maximum(abs(state[ring_theta]) for state in states) < 1.0e-12
    @test maximum(abs(state[carrier_theta] - (2 / 7) * state[sun_theta])
        for state in states) < 2.0e-10
    @test maximum(abs((state[planet_phase] - initial_planet_phase) +
        (20 / 21) * state[sun_theta]) for state in states) < 2.0e-9
    @test norm(spatial_marker_position(sun_planet.contact_marker,
        last(states)) - spatial_marker_position(
            sun_planet.contact_marker, first(states))) > 1.0e-4
    @test abs(last(states)[sun_planet.reaction_variable]) > 1.0e-3
    @test abs(last(states)[planet_ring.reaction_variable]) > 1.0e-3
end

@testset "Spatial coordinate coupler" begin
    path = joinpath(SPATIAL_MODEL_DIRECTORY, "screw-motion-coupler.toml")
    source = read(path, String)
    loaded = load_spatial_model(path)
    coupler = loaded.connections[:screw_coupler]
    hinge = connection_hinge(loaded.connections[:screw_bearing])
    inline = connection_inline(loaded.connections[:nut_guide])
    initial = loaded.initial_values

    @test coupler isa SpatialCoordinateCoupler
    @test coupler.coordinate_kind == :mixed
    @test coupler.coefficients == [1.0, -0.08]
    @test coupler.offset == 0.0
    @test !isempty(hinge.rotation_variables)
    @test !isempty(inline.translation_variables)
    @test loaded.analysis.degrees_of_freedom == 1
    @test loaded.state_selection.selected_velocities ==
        [Symbol("screw_bearing.omega")]
    @test initial[coupler.reaction_variable] ≈ 0.1098901098901099
    @test initial[inline.translation_variables[1]] ≈
        0.08 * initial[hinge.rotation_variables[1]]

    derivative = SpatialSimulationRunner.initial_spatial_derivative(
        initial, loaded)
    selection = AnalysisSelection(Dynamics(), loaded.active_variable_indices,
        loaded.active_equation_indices)
    analytical = Matrix(evaluate_analysis_sparse_jacobian(loaded.model,
        selection, 0.0, initial, derivative, 2.1))
    numerical = similar(analytical)
    step = 1.0e-7
    plus = zeros(size(analytical, 1))
    minus = similar(plus)
    for (column, variable) in enumerate(loaded.active_variable_indices)
        state_plus, state_minus = copy(initial), copy(initial)
        rate_plus, rate_minus = copy(derivative), copy(derivative)
        state_plus[variable] += step
        state_minus[variable] -= step
        rate_plus[variable] += 2.1step
        rate_minus[variable] -= 2.1step
        evaluate_analysis_equations!(plus, loaded.model, selection, 0.0,
            state_plus, rate_plus)
        evaluate_analysis_equations!(minus, loaded.model, selection, 0.0,
            state_minus, rate_minus)
        numerical[:, column] .= (plus .- minus) ./ (2step)
    end
    @test norm(analytical - numerical, Inf) < 3.0e-7

    result = run_spatial_model(path; duration = 0.2, samples = 21)
    theta = hinge.rotation_variables[3]
    distance = inline.translation_variables[3]
    @test maximum(abs(state[distance] - 0.08state[theta])
        for state in result.states) < 1.0e-10
    @test last(result.states)[theta] > 0.05

    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "nut_guide.distance" => "nut_guide.rotation")))
    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "coefficients = [1.0, -0.08]" => "coefficients = [1.0, 0.0]")))
end

@testset "Spatial bushing" begin
    path = joinpath(SPATIAL_MODEL_DIRECTORY, "bushing-supported-body.toml")
    source = read(path, String)
    loaded = load_spatial_model(path)
    bushing = loaded.forces[:support]
    initial = loaded.initial_values

    @test bushing isa SpatialBushingComponent
    @test loaded.analysis.degrees_of_freedom == 6
    @test length(initial) == 46
    @test bushing.translational_damping ≈ [24.0, 24.0, 36.0]
    @test bushing.rotational_damping ≈ [1.2, 1.8, 1.2]
    @test initial[bushing.translation_variables] ≈ zeros(3)
    @test initial[bushing.angle_variables] ≈ zeros(3)
    @test initial[bushing.translation_rate_variables] ≈ [0.0, 0.3, 0.0]
    @test initial[bushing.angular_velocity_variables] ≈ [0.35, 0.0, 0.2]
    @test initial[bushing.local_force_variables] ≈ [0.0, -7.2, 0.0]
    @test initial[bushing.local_torque_variables] ≈ [-0.42, 0.0, -0.24]

    for axis in 1:3
        axis_vector = [axis == component ? 1.0 : 0.0 for component in 1:3]
        oriented_source = replace(source,
            "angular_velocity = [0.35, 0.0, 0.20]" =>
                "angular_velocity = [0.35, 0.0, 0.20]\n" *
                "orientation = [\"15 deg\", $(join(axis_vector, ", "))]")
        oriented = load_spatial_model(IOBuffer(oriented_source))
        oriented_bushing = oriented.forces[:support]
        expected_angles = zeros(3)
        expected_angles[axis] = pi / 12
        @test oriented.initial_values[oriented_bushing.angle_variables] ≈
            expected_angles atol = 1.0e-12
    end

    function bushing_jacobian_error(model)
        state = model.initial_values
        derivative = SpatialSimulationRunner.initial_spatial_derivative(
            state, model)
        selection = AnalysisSelection(Dynamics(),
            model.active_variable_indices, model.active_equation_indices)
        equations = zeros(length(selection.equation_indices))
        evaluate_analysis_equations!(equations, model.model, selection, 0.0,
            state, derivative)
        norm(equations, Inf) < 1.0e-11 || return Inf
        coefficient = 2.1
        analytical = Matrix(evaluate_analysis_sparse_jacobian(model.model,
            selection, 0.0, state, derivative, coefficient))
        numerical = similar(analytical)
        step = 1.0e-7
        plus, minus = similar(equations), similar(equations)
        for (column, variable) in enumerate(model.active_variable_indices)
            state_plus, state_minus = copy(state), copy(state)
            rate_plus, rate_minus = copy(derivative), copy(derivative)
            state_plus[variable] += step
            state_minus[variable] -= step
            rate_plus[variable] += coefficient * step
            rate_minus[variable] -= coefficient * step
            evaluate_analysis_equations!(plus, model.model, selection, 0.0,
                state_plus, rate_plus)
            evaluate_analysis_equations!(minus, model.model, selection, 0.0,
                state_minus, rate_minus)
            numerical[:, column] .= (plus .- minus) ./ (2step)
        end
        norm(analytical - numerical, Inf)
    end
    @test bushing_jacobian_error(loaded) < 5.0e-7

    large_z_source = replace(source,
        "angular_velocity = [0.35, 0.0, 0.20]" =>
            "angular_velocity = [0.35, 0.0, 0.20]\n" *
            "orientation = [\"120 deg\", 0.0, 0.0, 1.0]")
    large_z = load_spatial_model(IOBuffer(large_z_source))
    large_z_bushing = large_z.forces[:support]
    @test large_z.initial_values[large_z_bushing.angle_variables] ≈
        [0.0, 0.0, 2pi / 3] atol = 1.0e-12
    @test bushing_jacobian_error(large_z) < 2.0e-6

    result = run_spatial_model(path; duration = 0.1, samples = 11)
    static_state = first(result.states)
    @test result.static_initialization_iterations > 0
    @test static_state[bushing.translation_variables][3] ≈
        -9.81 / 900 atol = 2.0e-9
    @test static_state[bushing.angle_variables][2] > 0.1

    common_rotation_source = """
        [model]
        dimension = "spatial"

        [first]
        type = "rigid_body"
        mass = 1.0
        inertia = [1.0, 1.0, 1.0]
        position = [1.0, 0.0, 0.0]
        velocity = [0.0, 2.0, 0.0]
        angular_velocity = [0.0, 0.0, 2.0]

        [first.mount]
        type = "marker"

        [second]
        type = "rigid_body"
        mass = 1.0
        inertia = [1.0, 1.0, 1.0]
        angular_velocity = [0.0, 0.0, 2.0]

        [second.mount]
        type = "marker"

        [support]
        type = "bushing"
        markers = ["first.mount", "second.mount"]
        translational_stiffness = [1.0, 1.0, 1.0]
        rotational_stiffness = [1.0, 1.0, 1.0]
        """
    common_rotation = load_spatial_model(IOBuffer(common_rotation_source))
    rotating_bushing = common_rotation.forces[:support]
    @test rotating_bushing.reaction_marker isa SpatialFloatingMarker
    @test spatial_marker_position(rotating_bushing.reaction_marker,
        common_rotation.initial_values) ≈ [1.0, 0.0, 0.0]
    @test common_rotation.initial_values[
        rotating_bushing.translation_rate_variables] ≈ zeros(3)
    @test common_rotation.initial_values[
        rotating_bushing.angular_velocity_variables] ≈ zeros(3)
    @test bushing_jacobian_error(common_rotation) < 5.0e-7

    explicit_damping = replace(source,
        "damping_time_scale = 0.04" =>
            "damping_time_scale = 0.04\n" *
            "translational_damping = [1.0, 2.0, 3.0]\n" *
            "rotational_damping = [0.1, 0.2, 0.3]")
    explicit = load_spatial_model(IOBuffer(explicit_damping)).forces[:support]
    @test explicit.translational_damping == [1.0, 2.0, 3.0]
    @test explicit.rotational_damping == [0.1, 0.2, 0.3]

    static_only_source = replace(source,
        "damping_time_scale = 0.04" =>
            "damping_time_scale = 0.04\nactive_during = \"static\"")
    staged = run_spatial_model(IOBuffer(static_only_source);
        duration = 0.001, samples = 2)
    staged_bushing = staged.loaded.forces[:support]
    @test staged_bushing.active_during == (:static,)
    @test !staged_bushing.active[]
    @test first(staged.states)[staged_bushing.local_force_variables] ≈ zeros(3)
    @test first(staged.states)[staged_bushing.local_torque_variables] ≈ zeros(3)
    @test first(staged.states)[staged.loaded.bodies[:body].position_variables][3] < 0

    invalid_stage = replace(source,
        "damping_time_scale = 0.04" =>
            "damping_time_scale = 0.04\nactive_during = \"assembly\"")
    @test_throws ArgumentError load_spatial_model(IOBuffer(invalid_stage))

    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "translational_stiffness = [600.0, 600.0, 900.0]" =>
            "translational_stiffness = [600.0, -1.0, 900.0]")))
end

@testset "Spatial sphere-plane contact" begin
    source = """
        [model]
        dimension = "spatial"

        [simulation]
        end_time = 0.45
        output_samples = 46
        initial_step = 1.0e-6
        maximum_step = 0.002

        [ground]
        type = "ground"

        [ground.plane]
        type = "marker"

        [ball]
        type = "rigid_body"
        mass = 1.0
        inertia = [0.01, 0.01, 0.01]
        position = [0.0, 0.0, 0.70]

        [ball.center]
        type = "marker"

        [floor_contact]
        type = "plane_contact"
        markers = ["ball.center", "ground.plane"]
        radius = 0.50
        stiffness = 10000.0
        damping_factor = 0.15

        [gravity]
        type = "gravity"
        acceleration = [0.0, 0.0, -9.81]
        bodies = ["ball"]
        """
    loaded = load_spatial_model(IOBuffer(source))
    contact = loaded.forces[:floor_contact]
    @test contact isa SpatialPlaneContactComponent
    @test contact.geometry.axis.index == 3
    @test loaded.initial_values[contact.gap_variable] ≈ 0.20
    @test loaded.initial_values[contact.gap_rate_variable] ≈ 0.0
    @test loaded.initial_values[contact.normal_force_variable] ≈ 0.0

    lua_source = """
        local sim3d = require "sim3d"
        sim3d.model {dimension = "spatial"}
        sim3d.ground {name = "ground"}
        sim3d.marker {name = "ground.plane"}
        sim3d.rigid_body {
            name = "ball", mass = 1.0, inertia = {0.01, 0.01, 0.01},
            position = {0.0, 0.0, 0.70}
        }
        sim3d.marker {name = "ball.center"}
        sim3d.plane_contact {
            name = "floor_contact",
            markers = {"ball.center", "ground.plane"},
            radius = 0.50, stiffness = 10000.0
        }
        """
    lua_loaded = load_spatial_model(IOBuffer(lua_source); format = :lua,
        source_directory = normpath(joinpath(
            @__DIR__, "..", "..", "assemblies")))
    @test lua_loaded.forces[:floor_contact] isa SpatialPlaneContactComponent

    expression_source = replace(source,
        "stiffness = 10000.0\ndamping_factor = 0.15" =>
        "expression = \"max(0, -10000*floor_contact.gap - " *
            "2000*floor_contact.gap_rate)\"")
    expression_loaded = load_spatial_model(IOBuffer(expression_source))
    expression_contact = expression_loaded.forces[:floor_contact]
    @test expression_contact.expression
    @test expression_contact.stiffness == 0.0
    expression_state = copy(expression_loaded.initial_values)
    expression_body = expression_loaded.bodies[:ball]
    expression_state[expression_body.position_variables[3]] = 0.45
    expression_state[expression_body.velocity_variables[3]] = -1.0
    initialize_spatial_plane_contact!(expression_state, expression_contact)
    @test expression_state[expression_contact.gap_variable] ≈ -0.05
    @test expression_state[expression_contact.gap_rate_variable] ≈ -1.0
    @test expression_state[expression_contact.normal_force_variable] ≈ 2500.0

    unclamped = load_spatial_model(IOBuffer(replace(source,
        "stiffness = 10000.0\ndamping_factor = 0.15" =>
            "expression = \"-5\"")))
    @test unclamped.initial_values[
        unclamped.forces[:floor_contact].normal_force_variable] == -5.0
    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "damping_factor = 0.15" =>
            "damping_factor = 0.15\nexpression = \"0\"")))

    penetrating = copy(loaded.initial_values)
    body = loaded.bodies[:ball]
    penetrating[body.position_variables[3]] = 0.45
    penetrating[body.velocity_variables[3]] = -1.0
    initialize_spatial_plane_contact!(penetrating, contact)
    @test penetrating[contact.gap_variable] ≈ -0.05
    @test penetrating[contact.gap_rate_variable] ≈ -1.0
    @test penetrating[contact.normal_force_variable] ≈ 575.0
    @test penetrating[contact.global_force_variables] ≈ [0.0, 0.0, 575.0]
    values = spatial_plane_contact_values(contact, penetrating)
    @test values.contact_point ≈ [0.0, 0.0, 0.0]

    smooth_loaded = load_spatial_model(IOBuffer(replace(source,
        "damping_factor = 0.15" =>
            "damping_factor = 0.15\ntransition_depth = 0.001")))
    smooth_contact = smooth_loaded.forces[:floor_contact]
    smooth_state = copy(smooth_loaded.initial_values)
    smooth_state[smooth_loaded.bodies[:ball].position_variables[3]] = 0.4995
    initialize_spatial_plane_contact!(smooth_state, smooth_contact)
    @test smooth_contact.transition_depth ≈ 0.001
    @test smooth_state[smooth_contact.gap_variable] ≈ -0.0005
    @test smooth_state[smooth_contact.normal_force_variable] ≈ 0.615234375
    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "damping_factor = 0.15" =>
            "damping_factor = 0.15\ntransition_depth = -0.001")))

    derivative = SpatialSimulationRunner.initial_spatial_derivative(
        penetrating, loaded)
    selection = AnalysisSelection(Dynamics(), loaded.active_variable_indices,
        loaded.active_equation_indices)
    analytical = Matrix(evaluate_analysis_sparse_jacobian(loaded.model,
        selection, 0.0, penetrating, derivative, 1.7))
    equations = zeros(length(selection.equation_indices))
    numerical = similar(analytical)
    step = 1.0e-7
    plus, minus = similar(equations), similar(equations)
    for (column, variable) in enumerate(loaded.active_variable_indices)
        state_plus, state_minus = copy(penetrating), copy(penetrating)
        rate_plus, rate_minus = copy(derivative), copy(derivative)
        state_plus[variable] += step
        state_minus[variable] -= step
        rate_plus[variable] += 1.7step
        rate_minus[variable] -= 1.7step
        evaluate_analysis_equations!(plus, loaded.model, selection, 0.0,
            state_plus, rate_plus)
        evaluate_analysis_equations!(minus, loaded.model, selection, 0.0,
            state_minus, rate_minus)
        numerical[:, column] .= (plus .- minus) ./ (2step)
    end
    @test norm(analytical - numerical, Inf) < 2.0e-5

    expression_derivative =
        SpatialSimulationRunner.initial_spatial_derivative(
            expression_state, expression_loaded)
    expression_selection = AnalysisSelection(Dynamics(),
        expression_loaded.active_variable_indices,
        expression_loaded.active_equation_indices)
    expression_analytical = Matrix(evaluate_analysis_sparse_jacobian(
        expression_loaded.model, expression_selection, 0.0,
        expression_state, expression_derivative, 1.7))
    expression_equations = zeros(length(
        expression_selection.equation_indices))
    expression_numerical = similar(expression_analytical)
    expression_plus = similar(expression_equations)
    expression_minus = similar(expression_equations)
    for (column, variable) in enumerate(
            expression_loaded.active_variable_indices)
        state_plus, state_minus = copy(expression_state), copy(expression_state)
        rate_plus = copy(expression_derivative)
        rate_minus = copy(expression_derivative)
        state_plus[variable] += step
        state_minus[variable] -= step
        rate_plus[variable] += 1.7step
        rate_minus[variable] -= 1.7step
        evaluate_analysis_equations!(expression_plus,
            expression_loaded.model, expression_selection, 0.0,
            state_plus, rate_plus)
        evaluate_analysis_equations!(expression_minus,
            expression_loaded.model, expression_selection, 0.0,
            state_minus, rate_minus)
        expression_numerical[:, column] .=
            (expression_plus .- expression_minus) ./ (2step)
    end
    @test norm(expression_analytical - expression_numerical, Inf) < 2.0e-5

    result = run_spatial_model(IOBuffer(source))
    @test result.solution.stats.events_found == 0
    @test result.solution.stats.history_restarts == 0
    @test maximum(state[contact.normal_force_variable]
        for state in result.states) > 0

    staged = load_spatial_model(IOBuffer(replace(source,
        "damping_factor = 0.15" =>
            "damping_factor = 0.15\ninactive_during = \"static\"")))
    staged_contact = staged.forces[:floor_contact]
    state = copy(staged.initial_values)
    state[staged.bodies[:ball].position_variables[3]] = 0.45
    set_spatial_plane_contact_stage!(staged_contact, :static)
    initialize_spatial_plane_contact!(state, staged_contact)
    @test state[staged_contact.normal_force_variable] == 0.0
    set_spatial_plane_contact_stage!(staged_contact, :dynamic)
    initialize_spatial_plane_contact!(state, staged_contact)
    @test state[staged_contact.normal_force_variable] > 0.0

    moving_plane_source = """
        [model]
        dimension = "spatial"

        [plane]
        type = "rigid_body"
        mass = 2.0
        inertia = [0.1, 0.1, 0.1]

        [plane.surface]
        type = "marker"

        [sphere]
        type = "rigid_body"
        mass = 1.0
        inertia = [0.01, 0.01, 0.01]
        position = [0.0, 0.0, 0.45]

        [sphere.center]
        type = "marker"

        [contact]
        type = "plane_contact"
        markers = ["sphere.center", "plane.surface"]
        radius = 0.50
        stiffness = 10000.0
        """
    moving = load_spatial_model(IOBuffer(moving_plane_source))
    moving_contact = moving.forces[:contact]
    moving_state = copy(moving.initial_values)
    initialize_spatial_plane_contact!(moving_state, moving_contact)
    contribution = only(SpatialPlaneContacts.equation_contributions(
        moving_contact))
    contact_equations = zeros(length(moving.layout.catalog.equations))
    contribution.residual!(contact_equations, 0.0, moving_state,
        zeros(length(moving_state)))
    sphere_force = contact_equations[
        moving.bodies[:sphere].balance_equations[1:3]]
    plane_force = contact_equations[
        moving.bodies[:plane].balance_equations[1:3]]
    @test norm(sphere_force) > 0
    @test sphere_force + plane_force ≈ zeros(3)
end

@testset "Spatial rolling tire" begin
    path = joinpath(SPATIAL_MODEL_DIRECTORY, "driven-rolling-tire.toml")
    source = read(path, String)
    loaded = load_spatial_model(path)
    tire = loaded.forces[:tire]
    initial = loaded.initial_values

    @test tire isa SpatialTireComponent
    @test loaded.analysis.degrees_of_freedom == 3
    @test initial[tire.deflection_variable] ≈ 0.01 atol = 1.0e-14
    @test initial[tire.normal_force_variable] ≈ 150.0 atol = 1.0e-10
    @test initial[tire.slip_ratio_variable] == 0.0
    @test initial[tire.slip_angle_variable] > 1.0
    @test initial[tire.lateral_force_variable] ≈ -120.0 atol = 1.0e-10

    rolling = copy(initial)
    wheel = loaded.bodies[:wheel]
    spindle = loaded.bodies[:spindle]
    rolling[wheel.velocity_variables] .= [2.0, 0.0, 0.0]
    rolling[spindle.velocity_variables] .= [2.0, 0.0, 0.0]
    rolling[wheel.angular_velocity_variables] .= [0.0, 2.0 / 0.49, 0.0]
    tire_kinematics = spatial_tire_kinematics(tire, rolling)
    @test tire_kinematics.forward ≈ [1.0, 0.0, 0.0] atol = 1.0e-14
    @test tire_kinematics.lateral ≈ [0.0, 1.0, 0.0] atol = 1.0e-14
    @test tire_kinematics.longitudinal_slip_velocity ≈ 0.0 atol = 1.0e-12

    rolling[wheel.velocity_variables] .= [2.0, 0.3, 0.0]
    rolling[spindle.velocity_variables] .= [2.0, 0.3, 0.0]
    rolling[wheel.angular_velocity_variables] .= [0.0, 5.0, 0.0]
    initialize_spatial_tire!(rolling, tire, 0.0)
    @test rolling[tire.slip_ratio_variable] > 0
    @test rolling[tire.slip_angle_variable] > 0
    @test rolling[tire.longitudinal_force_variable] > 0
    @test rolling[tire.lateral_force_variable] < 0
    utilization = hypot(
        rolling[tire.longitudinal_force_variable] /
            (tire.mu_longitudinal * rolling[tire.normal_force_variable]),
        rolling[tire.lateral_force_variable] /
            (tire.mu_lateral * rolling[tire.normal_force_variable]))
    @test utilization <= 1 + 1.0e-12

    derivative = SpatialSimulationRunner.initial_spatial_derivative(
        initial, loaded)
    selection = AnalysisSelection(Dynamics(), loaded.active_variable_indices,
        loaded.active_equation_indices)
    equations = zeros(length(selection.equation_indices))
    evaluate_analysis_equations!(equations, loaded.model, selection, 0.0,
        initial, derivative)
    @test norm(equations, Inf) < 1.0e-11
    coefficient = 1.7
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
    @test norm(analytical - numerical, Inf) < 3.0e-6

    result = run_spatial_model(path; duration = 0.25, samples = 11)
    final = last(result.states)
    @test final[wheel.position_variables[1]] > 0
    @test final[wheel.angular_velocity_variables[2]] > 0
    @test final[tire.longitudinal_force_variable] > 0
    @test result.solution.stats.events_found == 0

    expression_source = replace(source,
        "normal_stiffness = 15000.0\nnormal_damping_time_scale = 0.015" =>
        "normal_expression = \"max(0, 15000*tire.deflection + " *
        "225*tire.deflection_rate)\"")
    expression_model = load_spatial_model(IOBuffer(expression_source))
    expression_tire = expression_model.forces[:tire]
    @test expression_tire.normal_expression
    @test expression_model.initial_values[
        expression_tire.normal_force_variable] ≈ 150.0 atol = 1.0e-10

    transient_source = replace(source,
        "normal_damping_time_scale = 0.015" =>
            "normal_damping_time_scale = 0.015\n" *
            "longitudinal_relaxation_length = 0.30\n" *
            "lateral_relaxation_length = 0.45",
        "tire.normal_force * longitudinal_stiffness_per_load * tire.slip_ratio" =>
            "tire.normal_force * longitudinal_stiffness_per_load / 0.30 * tire.longitudinal_deformation",
        "-tire.normal_force * lateral_stiffness_per_load * tire.slip_angle" =>
            "-tire.normal_force * lateral_stiffness_per_load / 0.45 * tire.lateral_deformation")
    transient = load_spatial_model(IOBuffer(transient_source))
    transient_tire = transient.forces[:tire]
    @test transient_tire.transient
    @test transient_tire.longitudinal_relaxation_length == 0.30
    @test transient_tire.lateral_relaxation_length == 0.45
    @test length(transient_tire.deformation_equations) == 2

    stationary = copy(transient.initial_values)
    for body in values(transient.bodies)
        stationary[body.velocity_variables] .= 0
        stationary[body.angular_velocity_variables] .= 0
    end
    stationary[transient_tire.longitudinal_deformation_variable] = 0.01
    stationary[transient_tire.lateral_deformation_variable] = -0.02
    initialize_spatial_tire!(stationary, transient_tire, 0.0)
    @test collect(tire_deformation_rates(transient_tire, stationary)) ≈
        zeros(2) atol = 1.0e-14
    @test stationary[transient_tire.longitudinal_force_variable] > 0
    @test stationary[transient_tire.lateral_force_variable] > 0

    unloaded = copy(stationary)
    unloaded[transient_tire.normal_force_variable] = 0.0
    @test unloaded[transient_tire.deflection_variable] > 0
    @test collect(tire_deformation_rates(transient_tire, unloaded)) ≈ [
        -transient_tire.regularization_speed *
            unloaded[transient_tire.longitudinal_deformation_variable] /
            transient_tire.longitudinal_relaxation_length,
        -transient_tire.regularization_speed *
            unloaded[transient_tire.lateral_deformation_variable] /
            transient_tire.lateral_relaxation_length,
    ] atol = 1.0e-14
    tiny_limited = SpatialTires.limited_tire_forces(
        transient_tire, 1.0e-300, 1000.0, -2000.0)
    @test all(isfinite, tiny_limited)
    @test hypot(tiny_limited[1] / transient_tire.mu_longitudinal,
        tiny_limited[2] / transient_tire.mu_lateral) ≈ 1.0e-300

    transient_derivative =
        SpatialSimulationRunner.initial_spatial_derivative(
            transient.initial_values, transient)
    transient_selection = AnalysisSelection(Dynamics(),
        transient.active_variable_indices, transient.active_equation_indices)
    transient_equations = zeros(length(transient_selection.equation_indices))
    evaluate_analysis_equations!(transient_equations, transient.model,
        transient_selection, 0.0, transient.initial_values,
        transient_derivative)
    @test norm(transient_equations, Inf) < 1.0e-11
    transient_coefficient = 1.7
    transient_jacobian = Matrix(evaluate_analysis_sparse_jacobian(
        transient.model, transient_selection, 0.0,
        transient.initial_values, transient_derivative,
        transient_coefficient))
    equation_location = Dict(equation => location for (location, equation) in
        enumerate(transient.active_equation_indices))
    variable_location = Dict(variable => location for (location, variable) in
        enumerate(transient.active_variable_indices))
    deformation_rows = [equation_location[equation] for equation in
        transient_tire.deformation_equations]
    deformation_columns = [transient_tire.deflection_variable,
        transient_tire.normal_force_variable,
        transient_tire.forward_velocity_variable,
        transient_tire.longitudinal_slip_velocity_variable,
        transient_tire.lateral_slip_velocity_variable,
        transient_tire.longitudinal_deformation_variable,
        transient_tire.lateral_deformation_variable]
    numerical_deformation_jacobian = zeros(2, length(deformation_columns))
    transient_plus = similar(transient_equations)
    transient_minus = similar(transient_equations)
    for (column, variable) in enumerate(deformation_columns)
        state_plus = copy(transient.initial_values)
        state_minus = copy(transient.initial_values)
        rate_plus = copy(transient_derivative)
        rate_minus = copy(transient_derivative)
        state_plus[variable] += step
        state_minus[variable] -= step
        rate_plus[variable] += transient_coefficient * step
        rate_minus[variable] -= transient_coefficient * step
        evaluate_analysis_equations!(transient_plus, transient.model,
            transient_selection, 0.0, state_plus, rate_plus)
        evaluate_analysis_equations!(transient_minus, transient.model,
            transient_selection, 0.0, state_minus, rate_minus)
        numerical_deformation_jacobian[:, column] .=
            (transient_plus[deformation_rows] .-
             transient_minus[deformation_rows]) ./ (2step)
    end
    @test transient_jacobian[deformation_rows,
        [variable_location[variable] for variable in deformation_columns]] ≈
        numerical_deformation_jacobian atol = 1.0e-8
    transient_partition =
        SpatialSimulationRunner.runtime_spatial_state_partition(transient)
    differential, error_control, _ =
        SpatialSimulationRunner.spatial_state_masks(
            transient, transient_partition)
    for variable in (transient_tire.longitudinal_deformation_variable,
            transient_tire.lateral_deformation_variable)
        @test differential[variable_location[variable]]
        @test error_control[variable_location[variable]]
    end

    modal_transient_source = replace(transient_source,
        "mode = \"dynamic\"" => "mode = \"modal\"\nmodes = 2")
    modal_transient = run_spatial_model(IOBuffer(modal_transient_source))
    @test modal_transient.analysis_mode == :modal
    @test length(modal_transient.differential_columns) ==
        2 * modal_transient.loaded.analysis.degrees_of_freedom + 2
    @test length(modal_transient.eigenvalues) == 2
    @test maximum(modal_transient.equation_errors) < 1.0e-12

    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "normal_damping_time_scale = 0.015" =>
            "normal_damping_time_scale = 0.015\n" *
            "longitudinal_relaxation_length = 0.30")))

    drop_source = """
        [model]
        dimension = "spatial"

        [simulation]
        end_time = 0.45
        output_samples = 31
        maximum_step = 0.002

        [ground]
        type = "ground"

        [ground.road]
        type = "marker"

        [wheel]
        type = "rigid_body"
        mass = 1.0
        inertia = [0.1, 0.1, 0.1]
        position = [0.0, 0.0, 0.70]

        [wheel.center]
        type = "marker"
        orientation = ["-90 deg", 1.0, 0.0, 0.0]

        [tire]
        type = "rolling_tire"
        markers = ["wheel.center", "ground.road"]
        radius = 0.50
        normal_stiffness = 5000.0
        normal_damping_time_scale = 0.02
        longitudinal_expression = "0"
        lateral_expression = "0"

        [gravity]
        type = "gravity"
        acceleration = [0.0, 0.0, -9.81]
        bodies = ["wheel"]
        """
    drop = run_spatial_model(IOBuffer(drop_source))
    drop_tire = drop.loaded.forces[:tire]
    @test drop.solution.stats.events_found == 0
    @test drop.solution.stats.history_restarts == 0
    @test maximum(state[drop_tire.normal_force_variable]
        for state in drop.states) > 0

    lateral_liftoff_source = """
        [model]
        dimension = "spatial"

        [simulation]
        end_time = 0.15
        output_samples = 31
        maximum_step = 0.002

        [ground]
        type = "ground"

        [ground.road]
        type = "marker"

        [wheel]
        type = "rigid_body"
        mass = 10.0
        inertia = [0.5, 0.5, 0.5]
        position = [0.0, 0.0, 0.49]
        velocity = [30.0, 4.0, 0.2]
        angular_velocity = [0.0, 60.0, 0.0]

        [wheel.center]
        type = "marker"
        orientation = ["-90 deg", 1.0, 0.0, 0.0]

        [tire]
        type = "rolling_tire"
        markers = ["wheel.center", "ground.road"]
        radius = 0.50
        regularization_speed = 0.1
        normal_stiffness = 5000.0
        normal_damping_time_scale = 0.02
        longitudinal_relaxation_length = 0.30
        lateral_relaxation_length = 0.45
        longitudinal_expression = "10000*tire.longitudinal_deformation"
        lateral_expression = "-10000*tire.lateral_deformation"
        friction_limit = "ellipse"
        mu_longitudinal = 1.0
        mu_lateral = 1.0
        """
    lateral_liftoff = run_spatial_model(IOBuffer(lateral_liftoff_source))
    liftoff_tire = lateral_liftoff.loaded.forces[:tire]
    liftoff_final = last(lateral_liftoff.states)
    @test lateral_liftoff.solution.stats.events_found == 0
    @test lateral_liftoff.solution.stats.history_restarts == 0
    @test liftoff_final[liftoff_tire.deflection_variable] < 0
    @test liftoff_final[liftoff_tire.normal_force_variable] == 0
    @test liftoff_final[liftoff_tire.lateral_trial_force_variable] != 0
    @test abs(liftoff_final[liftoff_tire.lateral_force_variable]) < 1.0e-40
    @test isfinite(liftoff_final[
        liftoff_tire.lateral_deformation_variable])

    bad_axis = replace(source,
        "[wheel.center]\ntype = \"marker\"\n" *
        "orientation = [\"-90 deg\", 1.0, 0.0, 0.0]" =>
        "[wheel.center]\ntype = \"marker\"")
    @test_throws ArgumentError load_spatial_model(IOBuffer(bad_axis))
end

@testset "Spatial distance measurements" begin
    path = joinpath(SPATIAL_MODEL_DIRECTORY, "distance-measures.toml")
    loaded = load_spatial_model(path)
    span_measure = loaded.measures[:range]
    directed_measure = loaded.measures[:height]
    span = span_measure.span
    initial = loaded.initial_values

    @test span_measure isa SpatialSpanMeasure
    @test directed_measure isa SpatialDirectedDistanceMeasure
    @test loaded.analysis.degrees_of_freedom == 6
    @test length(loaded.layout.catalog.variables) == 38
    @test length(loaded.layout.catalog.equations) == 38
    @test length(loaded.active_variable_indices) == 38
    @test length(loaded.active_equation_indices) == 38

    separation = [1.0, 0.5, 0.75]
    velocity = [0.2, -0.1, 0.3]
    distance = norm(separation)
    @test initial[span.separation_variables] ≈ separation
    @test initial[span.distance_variable] ≈ distance
    @test initial[span.direction_variables] ≈ separation / distance
    @test initial[span.velocity_variable] ≈
        dot(separation / distance, velocity)
    @test initial[span.acceleration_variable] ≈
        (dot(velocity, velocity) -
         initial[span.velocity_variable]^2) / distance
    @test initial[directed_measure.distance_variable] ≈ 0.75
    @test initial[directed_measure.velocity_variable] ≈ 0.3
    @test initial[directed_measure.acceleration_variable] ≈ 0.0 atol = 1.0e-13

    force = loaded.forces[:measurement_probe]
    expected_dependencies = sort([
        span.distance_variable, directed_measure.velocity_variable])
    @test sort(force.magnitude.dependencies) == expected_dependencies
    @test initial[force.magnitude_variable] == 0.0

    derivative = SpatialSimulationRunner.initial_spatial_derivative(
        initial, loaded)
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
    @test norm(analytical - numerical, Inf) < 2.0e-7

    result = run_spatial_model(path; duration = 0.2, samples = 11)
    @test length(result.states) == 11
    for state in result.states
        expected_separation = spatial_marker_position(span.marker_2, state) -
            spatial_marker_position(span.marker_1, state)
        @test state[span.distance_variable] ≈
            norm(expected_separation) atol = 2.0e-6
        marker_position = spatial_marker_position(
            directed_measure.geometry.marker_i, state)
        @test state[directed_measure.distance_variable] ≈
            marker_position[3] atol = 2.0e-8
    end

    coincident = replace(read(path, String),
        "position = [1.0, 0.5, 0.75]" => "position = [0.0, 0.0, 0.0]")
    @test_throws ArgumentError load_spatial_model(IOBuffer(coincident))
end

@testset "Spatial static analysis" begin
    static_path = joinpath(SPATIAL_MODEL_DIRECTORY,
        "quasistatic-torsional-pendulum.toml")
    loaded = load_spatial_model(static_path)
    @test loaded.analysis.mode == :static
    @test loaded.analysis.initialization == :none
    @test loaded.analysis.static_method == :newton

    static_progress = Any[]
    result = run_spatial_model(static_path; samples = 3,
        static_progress = event -> push!(static_progress, event))
    hinge = result.loaded.connections[:pin].hinge
    spherical = result.loaded.connections[:pin].spherical
    body = result.loaded.bodies[:pendulum]
    angles = [state[hinge.rotation_variables[3]] for state in result.states]
    @test result.analysis_mode == :static
    @test result.static_mass_regularized == falses(3)
    @test length(result.states) == 3
    @test result.continuation_solutions == 3
    @test static_progress[1].phase == :initial_conditions
    @test static_progress[1].status == :entered
    @test static_progress[2].phase == :initial_conditions
    @test static_progress[2].status == :consistent
    @test static_progress[1].state[body.position_variables] ≈
        result.loaded.entered_initial_values[body.position_variables]
    @test static_progress[2].state[body.position_variables] ≈
        result.loaded.initial_values[body.position_variables]
    static_solution_progress = filter(
        event -> event.phase != :initial_conditions, static_progress)
    @test first(static_solution_progress).status == :initial
    @test last(static_progress).status == :converged
    @test count(event -> event.status == :converged,
        static_progress) == 3
    @test all(event -> event.phase == :newton, static_solution_progress)
    @test all(>(0), result.static_iterations)
    @test issorted(angles)
    for (time, state, angle) in zip(result.times, result.states, angles)
        @test -8angle + 4time - 4.905cos(angle) ≈ 0.0 atol = 2.0e-10
        @test state[spherical.reaction_variables] ≈
            [0.0, 9.81, 0.0] atol = 2.0e-10
        @test norm(state[body.euler_parameter_variables]) ≈ 1.0 atol = 1.0e-13
        @test all(iszero, state[body.velocity_variables])
        @test all(iszero, state[body.angular_velocity_variables])
        @test all(iszero, state[body.acceleration_variables])
        @test all(iszero, state[body.angular_acceleration_variables])
    end

    dynamic_path = joinpath(SPATIAL_MODEL_DIRECTORY,
        "static-initialized-torsional-pendulum.toml")
    dynamic = run_spatial_model(dynamic_path; duration = 0.05, samples = 6)
    dynamic_hinge = dynamic.loaded.connections[:pin].hinge
    @test dynamic.analysis_mode == :dynamic
    @test dynamic.static_initialization_iterations == 5
    @test dynamic.static_relaxation_cycles == 0
    @test dynamic.static_mass_regularized === false
    @test first(dynamic.states)[dynamic_hinge.rotation_variables[3]] ≈
        angles[1] atol = 2.0e-10
    @test first(dynamic.states)[dynamic_hinge.rotation_variables[2]] ≈
        1.0 atol = 2.0e-12
    @test last(dynamic.states)[dynamic_hinge.rotation_variables[3]] >
        first(dynamic.states)[dynamic_hinge.rotation_variables[3]]

    mktempdir() do directory
        initialized_path = joinpath(directory, "initialized.simp")
        write_result(initialized_path, dynamic)
        @test read_result(initialized_path).static_mass_regularized === false
        base_source = replace(read(dynamic_path, String),
            "initialization = \"static_equilibrium\"" =>
                "initialization = \"none\"",
            "position = [0.4330127018922193, 0.25, 0.0]" =>
                "position = [0.1, 0.2, 0.3]",
            "orientation = [\"30 deg\", 0.0, 0.0, 1.0]" =>
                "orientation = [0.0, 0.0, 0.0, 1.0]")
        static_transfer_source = replace(base_source,
            "[pin.initial.impose]\nomega = 1.0" =>
                "[pin.initial.impose]\nomega = 0.25") * """

            [initial_conditions]
            result = "initialized.simp"
            sample = "static"
            include_velocities = false
            """
        static_transfer_path = joinpath(directory, "static-transfer.toml")
        write(static_transfer_path, static_transfer_source)
        static_transfer = run_spatial_model(static_transfer_path;
            duration = 0.0, samples = 1)
        transferred_body = static_transfer.loaded.bodies[:pendulum]
        transferred_hinge = static_transfer.loaded.connections[:pin].hinge
        transferred_state = only(static_transfer.states)
        source_body = dynamic.loaded.bodies[:pendulum]
        source_state = first(dynamic.states)
        @test static_transfer.loaded.initial_conditions.enabled
        @test static_transfer.loaded.initial_conditions.result_path ==
            initialized_path
        @test static_transfer.loaded.initial_conditions.requested_sample ==
            "static"
        @test static_transfer.loaded.initial_conditions.sample == 1
        @test static_transfer.loaded.initial_conditions.source_time ==
            first(dynamic.times)
        @test !static_transfer.loaded.initial_conditions.include_velocities
        @test transferred_state[transferred_body.position_variables] ≈
            source_state[source_body.position_variables] atol = 2.0e-7
        @test transferred_state[transferred_body.euler_parameter_variables] ≈
            source_state[source_body.euler_parameter_variables] atol = 2.0e-7
        @test norm(transferred_state[
            transferred_body.euler_parameter_variables]) ≈ 1.0
        @test transferred_state[transferred_hinge.rotation_variables[3]] ≈
            source_state[dynamic_hinge.rotation_variables[3]] atol = 2.0e-7
        @test transferred_state[transferred_hinge.rotation_variables[2]] ≈
            0.25 atol = 2.0e-12
        @test Symbol("pendulum.p_0") in
            static_transfer.loaded.initial_conditions.transferred_variables
        @test Symbol("pin.theta") in
            static_transfer.loaded.initial_conditions.transferred_variables

        moving_source = replace(base_source,
            "[pin.initial.impose]\nomega = 1.0\n" => "") * """

            [initial_conditions]
            result = "initialized.simp"
            sample = "last"
            include_velocities = true
            """
        moving_path = joinpath(directory, "moving-transfer.toml")
        write(moving_path, moving_source)
        moving_transfer = run_spatial_model(moving_path;
            duration = 0.0, samples = 1)
        moving_body = moving_transfer.loaded.bodies[:pendulum]
        moving_hinge = moving_transfer.loaded.connections[:pin].hinge
        moving_state = only(moving_transfer.states)
        final_source = last(dynamic.states)
        @test moving_transfer.loaded.initial_conditions.sample ==
            length(dynamic.times)
        @test moving_transfer.loaded.initial_conditions.include_velocities
        @test moving_state[moving_body.position_variables] ≈
            final_source[source_body.position_variables] atol = 2.0e-7
        @test moving_state[moving_body.euler_parameter_variables] ≈
            final_source[source_body.euler_parameter_variables] atol = 2.0e-7
        @test moving_state[moving_body.velocity_variables] ≈
            final_source[source_body.velocity_variables] atol = 2.0e-7
        @test moving_state[moving_body.angular_velocity_variables] ≈
            final_source[source_body.angular_velocity_variables] atol = 2.0e-7
        @test moving_state[moving_hinge.rotation_variables[3]] ≈
            final_source[dynamic_hinge.rotation_variables[3]] atol = 2.0e-7
        @test moving_state[moving_hinge.rotation_variables[2]] ≈
            final_source[dynamic_hinge.rotation_variables[2]] atol = 2.0e-7

        static_path_result = joinpath(directory, "static.simp")
        write_result(static_path_result, result)
        @test read_result(static_path_result).static_mass_regularized ==
            falses(3)
        standalone_source = replace(static_transfer_source,
            "result = \"initialized.simp\"" => "result = \"static.simp\"",
            "sample = \"static\"" => "sample = \"last\"")
        standalone_path = joinpath(directory, "standalone-transfer.toml")
        write(standalone_path, standalone_source)
        standalone_transfer = run_spatial_model(standalone_path;
            duration = 0.0, samples = 1)
        standalone_body = standalone_transfer.loaded.bodies[:pendulum]
        @test only(standalone_transfer.states)[
            standalone_body.position_variables] ≈
            last(result.states)[body.position_variables] atol = 2.0e-7
        @test standalone_transfer.loaded.initial_conditions.sample == 3

        incompatible_source = replace(static_transfer_source,
            "[pin]\ntype = \"revolute\"" =>
                "[pin]\ntype = \"hinge\"")
        incompatible_path = joinpath(directory, "incompatible.toml")
        write(incompatible_path, incompatible_source)
        @test_throws ArgumentError load_spatial_model(incompatible_path)
    end

    relaxation_source = replace(read(static_path, String),
        "mode = \"static\"" => """mode = "static"
        static_method = "dynamic_relaxation"
        relaxation_duration = 0.25
        relaxation_min_cycles = 2
        relaxation_max_cycles = 16
        handoff_acceleration = 1.0e6
        handoff_speed = 1.0e6
        handoff_correction = 1.0e6""")
    relaxation_progress = Any[]
    relaxed = run_spatial_model(IOBuffer(relaxation_source);
        duration = 0.0, samples = 1,
        static_progress = event -> push!(relaxation_progress, event))
    @test relaxed.loaded.analysis.relaxation_min_cycles == 2
    @test relaxed.loaded.analysis.relaxation_max_cycles == 16
    @test relaxed.loaded.analysis.handoff_acceleration == 1.0e6
    @test relaxed.loaded.analysis.handoff_speed == 1.0e6
    @test relaxed.loaded.analysis.handoff_correction == 1.0e6
    @test only(relaxed.static_relaxation_cycles) == 2
    @test any(event -> event.phase == :dynamic_relaxation,
        relaxation_progress)
    @test all(event -> hasproperty(event, :force_imbalance) &&
        hasproperty(event, :torque_imbalance) &&
        hasproperty(event, :constraint_error), relaxation_progress)
    @test all(event -> isfinite(event.force_imbalance) &&
        isfinite(event.torque_imbalance) &&
        isfinite(event.constraint_error), relaxation_progress)
    @test last(relaxation_progress).status == :converged
    @test last(relaxation_progress).force_imbalance <= 1.0e-8
    @test last(relaxation_progress).torque_imbalance <= 1.0e-8
    @test last(relaxation_progress).constraint_error <= 1.0e-10
    @test length(relaxed.static_progress_events) == length(relaxation_progress)
    mktempdir() do directory
        path = joinpath(directory, "static-progress.simp")
        write_result(path, relaxed)
        snapshots = read_result(path).static_snapshots
        @test length(snapshots) == length(relaxation_progress)
        @test snapshots[1].phase == :initial_conditions
        @test snapshots[1].status == :entered
        @test snapshots[2].phase == :initial_conditions
        @test snapshots[2].status == :consistent
        @test snapshots[3].phase == :dynamic_relaxation
        @test last(snapshots).phase == :newton_polish
        @test last(snapshots).status == :converged
        @test isapprox(last(snapshots).values,
            last(relaxation_progress).state; rtol = 1.0e-6, atol = 1.0e-7)
        @test last(snapshots).force_imbalance ≈
            last(relaxation_progress).force_imbalance
    end
    @test only(relaxed.states)[
        relaxed.loaded.connections[:pin].hinge.rotation_variables[3]] ≈
        angles[1] atol = 2.0e-10

    unpolished_source = replace(relaxation_source,
        "handoff_correction = 1.0e6" =>
            "handoff_correction = 1.0e6\n        relaxation_polish = false")
    unpolished_progress = Any[]
    unpolished = run_spatial_model(IOBuffer(unpolished_source);
        duration = 0.0, samples = 1,
        static_progress = event -> push!(unpolished_progress, event))
    @test !unpolished.loaded.analysis.relaxation_polish
    @test only(unpolished.static_iterations) == 0
    @test only(unpolished.static_relaxation_cycles) == 2
    @test !any(event -> event.phase == :newton_polish,
        unpolished_progress)
    @test last(unpolished_progress).phase == :dynamic_relaxation
    @test last(unpolished_progress).status == :converged

    inverted_path = joinpath(SPATIAL_MODEL_DIRECTORY,
        "static-inverted-spherical-pendulum.toml")
    inverted_progress = Any[]
    inverted = run_spatial_model(inverted_path;
        static_progress = event -> push!(inverted_progress, event))
    inverted_state = only(inverted.states)
    inverted_body = inverted.loaded.bodies[:pendulum]
    inverted_pin = inverted.loaded.connections[:pin]
    @test inverted.loaded.analysis.degrees_of_freedom == 3
    @test only(inverted.static_relaxation_cycles) > 1
    @test any(event -> event.status == :waiting, inverted_progress)
    @test any(event -> event.status == :ready, inverted_progress)
    inverted_polish = filter(event -> event.phase == :newton_polish,
        inverted_progress)
    @test maximum(event.constraint_error for event in inverted_polish) <
        1.0e-5
    @test inverted.static_mass_regularized == Bool[true]
    @test inverted_state[inverted_body.position_variables] ≈
        [0.0, 0.0, -0.5] atol = 2.0e-9
    @test inverted_state[inverted_pin.reaction_variables] ≈
        [0.0, 0.0, 9.81] atol = 2.0e-9
    @test all(iszero,
        inverted_state[inverted_body.angular_velocity_variables])

    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(
        read(static_path, String), "mode = \"static\"" =>
            "mode = \"static\"\nrelaxation_cycles = 4")))
    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(
        read(static_path, String), "mode = \"static\"" =>
            "mode = \"static\"\nrelaxation_min_cycles = 4\n" *
            "relaxation_max_cycles = 3")))
    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(
        read(static_path, String), "mode = \"static\"" =>
            "mode = \"static\"\nhandoff_acceleration = 0.0")))

    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(
        read(static_path, String), "mode = \"static\"" =>
            "mode = \"static\"\ninitialization = \"static_equilibrium\"")))
end

@testset "Spatial modal analysis" begin
    path = joinpath(SPATIAL_MODEL_DIRECTORY,
        "modal-revolute-pendulum.toml")
    source = read(path, String)
    result = run_spatial_model(path)
    body = result.loaded.bodies[:pendulum]
    hinge = result.loaded.connections[:pin].hinge
    spring = result.loaded.forces[:spring]
    effective_inertia = body.inertia[3, 3] + body.mass * 0.5^2
    expected_decay = -spring.damping / (2effective_inertia)
    expected_frequency = sqrt(spring.stiffness / effective_inertia -
        expected_decay^2)

    @test result.analysis_mode == :modal
    @test result.loaded.analysis.number_of_modes == 1
    @test length(result.eigenvalues) == 1
    @test only(result.eigenvalues) ≈
        expected_decay + im * expected_frequency rtol = 2.0e-12
    @test only(result.natural_frequencies_hz) ≈
        sqrt(spring.stiffness / effective_inertia) / (2pi) rtol = 2.0e-12
    @test only(result.equation_errors) < 1.0e-12
    @test result.sparse_factorizations == 1
    @test length(result.differential_columns) == 2
    @test size(result.state_jacobian) == (30, 30)
    @test result.mode_shapes[hinge.rotation_variables[3], 1] ≈ 1
        atol = 2.0e-12
    pseudo_angle = @view result.mode_shapes[
        body.pseudo_angle_variables, 1]
    expected_parameter_mode = 0.5 .* quaternion_rate_matrix(
        @view(result.operating_state[body.euler_parameter_variables])) *
        pseudo_angle
    @test result.mode_shapes[body.euler_parameter_variables, 1] ≈
        expected_parameter_mode atol = 2.0e-12

    mktempdir() do directory
        result_path = joinpath(directory, "modal.simp")
        write_result(result_path, result)
        stored = read_result(result_path)
        @test stored.analysis_mode == :modal
        @test stored.modal_eigenvalues ≈ result.eigenvalues
        @test isapprox(stored.modal_mode_shapes, result.mode_shapes;
            rtol = 1.0e-6, atol = 1.0e-7)
    end

    static_source = replace(read(joinpath(SPATIAL_MODEL_DIRECTORY,
        "static-initialized-torsional-pendulum.toml"), String),
        "mode = \"dynamic\"" => "mode = \"modal\"\nmodes = 1")
    static_result = run_spatial_model(IOBuffer(static_source))
    static_hinge = static_result.loaded.connections[:pin].hinge
    @test static_result.analysis_mode == :modal
    @test static_result.static_initialization_iterations == 5
    @test static_result.operating_state[
        static_hinge.rotation_variables[2]] == 0
    @test length(static_result.eigenvalues) == 1

    chain = run_spatial_model(joinpath(SPATIAL_MODEL_DIRECTORY,
        "modal-three-link-pendulum.toml"))
    mass_matrix = [7 / 3 1.5 0.5;
                   1.5 4 / 3 0.5;
                   0.5 0.5 1 / 3]
    relative_angle_map = [1.0 0 0; -1 1 0; 0 -1 1]
    stiffness_matrix = transpose(relative_angle_map) *
        Diagonal([8.0, 6.0, 4.0]) * relative_angle_map +
        Diagonal(9.81 .* [2.5, 1.5, 0.5])
    damping_matrix = transpose(relative_angle_map) *
        Diagonal(0.02 .* [8.0, 6.0, 4.0]) * relative_angle_map
    reference_matrix = vcat(
        hcat(zeros(3, 3), Matrix{Float64}(I, 3, 3)),
        hcat(-(mass_matrix \ stiffness_matrix),
            -(mass_matrix \ damping_matrix)))
    reference_values = sort([value for value in eigvals(reference_matrix)
        if imag(value) > 0]; by = abs)
    @test chain.analysis_mode == :modal
    @test chain.loaded.analysis.degrees_of_freedom == 3
    @test chain.static_initialization_iterations == 0
    @test length(chain.eigenvalues) == 3
    @test chain.eigenvalues ≈ reference_values rtol = 2.0e-12
    @test maximum(chain.equation_errors) < 2.0e-12
    @test chain.sparse_factorizations == 1

    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "modes = 1" => "modes = 0")))
    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "frequency_shift_hz = 0.0" => "frequency_shift_hz = -1.0")))
end

@testset "Spatial spanning force" begin
    path = joinpath(SPATIAL_MODEL_DIRECTORY, "spanning-spring-body.toml")
    source = read(path, String)
    loaded = load_spatial_model(path)
    force = loaded.forces[:spring]
    initial = loaded.initial_values

    @test force isa SpatialSpanningForceComponent
    @test length(loaded.layout.catalog.variables) == 35
    @test length(loaded.layout.catalog.equations) == 35
    expected_span = spatial_marker_position(force.span.marker_2, initial) -
        spatial_marker_position(force.span.marker_1, initial)
    expected_length = norm(expected_span)
    @test initial[force.span.separation_variables] ≈ expected_span
    @test initial[force.span.distance_variable] ≈ expected_length
    @test norm(initial[force.span.direction_variables]) ≈ 1.0
    @test initial[force.force_variable] ≈
        -30 * (expected_length - 0.65) -
        initial[force.span.velocity_variable]
    @test initial[force.global_force_variables] ≈
        -initial[force.force_variable] .* initial[force.span.direction_variables]
    relative_velocity =
        spatial_marker_velocity(force.span.marker_2, initial) -
        spatial_marker_velocity(force.span.marker_1, initial)
    relative_acceleration =
        spatial_marker_acceleration(force.span.marker_2, initial) -
        spatial_marker_acceleration(force.span.marker_1, initial)
    direction = initial[force.span.direction_variables]
    distance_velocity = initial[force.span.velocity_variable]
    expected_acceleration = dot(direction, relative_acceleration) +
        (dot(relative_velocity, relative_velocity) - distance_velocity^2) /
        expected_length
    @test initial[force.span.acceleration_variable] ≈ expected_acceleration
    refreshed = copy(initial)
    refreshed[force.span.acceleration_variable] = 0.0
    initialize_spatial_spanning_force!(refreshed, force, 0.0)
    @test refreshed[force.span.acceleration_variable] ≈ expected_acceleration

    derivative = SpatialSimulationRunner.initial_spatial_derivative(
        initial, loaded)
    selection = AnalysisSelection(Dynamics(), loaded.active_variable_indices,
        loaded.active_equation_indices)
    equations = zeros(length(selection.equation_indices))
    evaluate_analysis_equations!(equations, loaded.model, selection, 0.0,
        initial, derivative)
    @test norm(equations, Inf) < 1.0e-12

    coefficient = 2.3
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
    @test norm(analytical - numerical, Inf) < 2.0e-7

    result = run_spatial_model(path; duration = 0.25, samples = 11)
    @test length(result.states) == 11
    @test maximum(abs(norm(state[force.span.direction_variables]) - 1)
        for state in result.states) < 2.0e-8
    @test maximum(abs(state[force.span.distance_variable] -
        norm(state[force.span.separation_variables]))
        for state in result.states) < 2.0e-8

    expression_source = replace(source,
        "[analysis]" => "[parameters]\nk = 30.0\nc = 1.0\nl0 = 0.65\n\n[analysis]",
        "stiffness = 30.0\ndamping = 1.0\nfree_length = 0.65" =>
            "expression = \"-k * (spring.length - l0) - c * spring.length_rate\"")
    expression_model = load_spatial_model(IOBuffer(expression_source))
    expression_force = expression_model.forces[:spring]
    @test expression_model.initial_values[expression_force.force_variable] ≈
        initial[force.force_variable]

    constant_source = replace(source,
        "stiffness = 30.0\ndamping = 1.0\nfree_length = 0.65" =>
            "force = 4.0")
    constant_model = load_spatial_model(IOBuffer(constant_source))
    constant_force = constant_model.forces[:spring]
    @test constant_model.initial_values[constant_force.force_variable] == 4.0

    static_only_source = replace(source,
        "free_length = 0.65" =>
            "free_length = 0.65\nactive_during = \"static\"")
    static_only_model = load_spatial_model(IOBuffer(static_only_source))
    static_only_force = static_only_model.forces[:spring]
    staged_state = copy(static_only_model.initial_values)
    @test static_only_force.active_during == (:static,)
    @test !static_only_force.active[]
    @test staged_state[static_only_force.force_variable] == 0.0
    @test staged_state[static_only_force.global_force_variables] ≈ zeros(3)
    @test staged_state[static_only_force.span.distance_variable] > 0.0
    SpatialSimulationRunner.set_spatial_analysis_stage!(static_only_model,
        :static; state = staged_state, time = 0.0)
    @test static_only_force.active[]
    @test staged_state[static_only_force.force_variable] != 0.0
    SpatialSimulationRunner.set_spatial_analysis_stage!(static_only_model,
        :dynamic; state = staged_state, time = 0.0)
    @test !static_only_force.active[]
    @test staged_state[static_only_force.force_variable] == 0.0
    @test staged_state[static_only_force.global_force_variables] ≈ zeros(3)

    dynamic_only_source = replace(source,
        "free_length = 0.65" =>
            "free_length = 0.65\ninactive_during = \"static\"")
    dynamic_only_model = load_spatial_model(IOBuffer(dynamic_only_source))
    dynamic_only_force = dynamic_only_model.forces[:spring]
    @test dynamic_only_force.active_during == (:dynamic, :modal)
    @test dynamic_only_force.active[]
    dynamic_only_state = copy(dynamic_only_model.initial_values)
    SpatialSimulationRunner.set_spatial_analysis_stage!(dynamic_only_model,
        :static; state = dynamic_only_state, time = 0.0)
    @test !dynamic_only_force.active[]
    @test dynamic_only_state[dynamic_only_force.force_variable] == 0.0
    SpatialSimulationRunner.set_spatial_analysis_stage!(dynamic_only_model,
        :dynamic; state = dynamic_only_state, time = 0.0)
    @test dynamic_only_force.active[]

    conflicting_stages = replace(source,
        "free_length = 0.65" => "free_length = 0.65\n" *
            "active_during = \"dynamic\"\ninactive_during = \"static\"")
    @test_throws ArgumentError load_spatial_model(IOBuffer(conflicting_stages))

    bad_expression = replace(expression_source, "spring.length" =>
        "body.acceleration_x")
    @test_throws ArgumentError load_spatial_model(IOBuffer(bad_expression))
end

@testset "Spatial marker-directed applied force" begin
    path = joinpath(SPATIAL_MODEL_DIRECTORY, "directed-applied-force.toml")
    source = read(path, String)
    loaded = load_spatial_model(path)
    force = loaded.forces[:push]
    initial = loaded.initial_values
    body = loaded.bodies[:body]

    @test force isa SpatialAppliedForceComponent
    @test force.direction_axis isa SpatialDirectedAxis
    @test length(loaded.layout.catalog.variables) == 26
    @test initial[force.magnitude_variable] == 12.0
    expected_direction = spatial_marker_orientation(
        force.direction_axis.marker, initial)[:, 3]
    @test initial[force.global_force_variables] ≈ 12 .* expected_direction
    @test initial[body.acceleration_variables] ≈
        initial[force.global_force_variables] ./ body.mass .+
        [0.0, 0.0, -9.81]

    result = run_spatial_model(path; duration = 0.2, samples = 9)
    @test length(result.states) == 9
    @test maximum(abs(norm(state[force.global_force_variables]) - 12)
        for state in result.states) < 1.0e-12

    static_only_source = replace(source, "force = 12.0" =>
        "force = 12.0\nactive_during = \"static\"")
    static_only_result = run_spatial_model(IOBuffer(static_only_source);
        duration = 0.02, samples = 3)
    static_only_force = static_only_result.loaded.forces[:push]
    @test static_only_force.active_during == (:static,)
    @test !static_only_force.active[]
    @test all(state[static_only_force.magnitude_variable] == 0.0
        for state in static_only_result.states)
    @test all(state[static_only_force.global_force_variables] ≈ zeros(3)
        for state in static_only_result.states)

    staged_source = replace(source, "force = 12.0" =>
        "force = 12.0\nactive_during = [\"dynamic\", \"modal\"]")
    staged_model = load_spatial_model(IOBuffer(staged_source))
    staged_force = staged_model.forces[:push]
    @test staged_force.active_during == (:dynamic, :modal)
    @test staged_force.active[]
    staged_state = copy(staged_model.initial_values)
    SpatialSimulationRunner.set_spatial_analysis_stage!(staged_model,
        :static; state = staged_state, time = 0.0)
    @test !staged_force.active[]
    @test staged_state[staged_force.magnitude_variable] == 0.0
    @test staged_state[staged_force.global_force_variables] ≈ zeros(3)
    SpatialSimulationRunner.set_spatial_analysis_stage!(staged_model,
        :modal; state = staged_state, time = 0.0)
    @test staged_force.active[]
    @test staged_state[staged_force.magnitude_variable] == 12.0

    invalid_stage = replace(source, "force = 12.0" =>
        "force = 12.0\nactive_during = \"assembly\"")
    @test_throws ArgumentError load_spatial_model(IOBuffer(invalid_stage))
    empty_stages = replace(source, "force = 12.0" =>
        "force = 12.0\nactive_during = []")
    @test_throws ArgumentError load_spatial_model(IOBuffer(empty_stages))

    reaction_source = replace(source,
        "[analysis]" =>
            "[parameters]\nf0 = 12.0\ngain = 0.4\n\n[analysis]",
        "markers = [\"body.application\", \"ground.force_axis\"]\nforce = 12.0" =>
            "markers = [\"body.application\", \"reaction.axis\"]\n" *
            "expression = \"f0 + gain * body.V_z\"\n" *
            "reaction_body = \"reaction\"",
        "acceleration = [0.0, 0.0, -9.81]" =>
            "acceleration = [0.0, 0.0, 0.0]") * """

        [reaction]
        type = "rigid_body"
        mass = 2.0
        inertia = [0.15, 0.18, 0.20]
        position = [-0.5, 0.2, 0.1]
        orientation = ["25 deg", 0.0, 1.0, 0.0]

        [reaction.axis]
        type = "marker"
        orientation = ["15 deg", 1.0, 0.0, 0.0]
        """
    reaction_model = load_spatial_model(IOBuffer(reaction_source))
    reaction_force = reaction_model.forces[:push]
    reaction_initial = reaction_model.initial_values
    reaction_marker = reaction_force.reaction_marker
    @test reaction_marker isa SpatialFloatingMarker
    @test reaction_marker.name == Symbol("push.reaction")
    @test spatial_marker_position(reaction_marker, reaction_initial) ≈
        spatial_marker_position(reaction_force.application_marker,
            reaction_initial)
    @test reaction_initial[reaction_force.magnitude_variable] == 12.0

    derivative = SpatialSimulationRunner.initial_spatial_derivative(
        reaction_initial, reaction_model)
    selection = AnalysisSelection(Dynamics(),
        reaction_model.active_variable_indices,
        reaction_model.active_equation_indices)
    equations = zeros(length(selection.equation_indices))
    evaluate_analysis_equations!(equations, reaction_model.model, selection,
        0.0, reaction_initial, derivative)
    @test norm(equations, Inf) < 1.0e-12

    coefficient = 1.9
    analytical = Matrix(evaluate_analysis_sparse_jacobian(
        reaction_model.model, selection, 0.0, reaction_initial, derivative,
        coefficient))
    numerical = similar(analytical)
    step = 1.0e-7
    plus, minus = similar(equations), similar(equations)
    for (column, variable) in enumerate(
            reaction_model.active_variable_indices)
        state_plus, state_minus = copy(reaction_initial), copy(reaction_initial)
        rate_plus, rate_minus = copy(derivative), copy(derivative)
        state_plus[variable] += step
        state_minus[variable] -= step
        rate_plus[variable] += coefficient * step
        rate_minus[variable] -= coefficient * step
        evaluate_analysis_equations!(plus, reaction_model.model, selection,
            0.0, state_plus, rate_plus)
        evaluate_analysis_equations!(minus, reaction_model.model, selection,
            0.0, state_minus, rate_minus)
        numerical[:, column] .= (plus .- minus) ./ (2step)
    end
    @test norm(analytical - numerical, Inf) < 3.0e-7

    same_body = replace(source, "force = 12.0" =>
        "force = 12.0\nreaction_body = \"body\"")
    @test_throws ArgumentError load_spatial_model(IOBuffer(same_body))
end

@testset "Spatial joint applied torque" begin
    path = joinpath(SPATIAL_MODEL_DIRECTORY,
        "torsional-spring-pendulum.toml")
    source = read(path, String)
    loaded = load_spatial_model(path)
    torque = loaded.forces[:spring]
    joint = loaded.connections[:pin]
    body = loaded.bodies[:pendulum]
    hinge = joint.hinge
    initial = loaded.initial_values

    @test torque isa SpatialAppliedTorqueComponent
    @test !isempty(hinge.rotation_variables)
    @test loaded.state_selection.selected_velocities == [Symbol("pin.omega")]
    @test all(index -> index in loaded.active_variable_indices,
        body.pseudo_angle_variables)
    @test all(index -> index in loaded.active_equation_indices,
        body.orientation_equations)
    @test all(index -> index ∉ loaded.active_equation_indices,
        body.pseudo_angle_state_equations)
    @test torque.stiffness == 8.0
    @test torque.damping == 0.24
    @test initial[hinge.rotation_variables[3]] ≈ pi / 6
    @test initial[torque.magnitude_variable] ≈ -8pi / 6
    expected_axis = spatial_marker_orientation(hinge.marker_j, initial)[:, 3]
    @test initial[torque.global_torque_variables] ≈
        initial[torque.magnitude_variable] .* expected_axis

    static_only_source = replace(source, "free_angle = \"0 deg\"" =>
        "free_angle = \"0 deg\"\nactive_during = \"static\"")
    static_only_model = load_spatial_model(IOBuffer(static_only_source))
    static_only_torque = static_only_model.forces[:spring]
    staged_state = copy(static_only_model.initial_values)
    @test static_only_torque.active_during == (:static,)
    @test !static_only_torque.active[]
    @test staged_state[static_only_torque.magnitude_variable] == 0.0
    @test staged_state[static_only_torque.global_torque_variables] ≈ zeros(3)
    SpatialSimulationRunner.set_spatial_analysis_stage!(static_only_model,
        :static; state = staged_state, time = 0.0)
    @test static_only_torque.active[]
    @test staged_state[static_only_torque.magnitude_variable] != 0.0
    SpatialSimulationRunner.set_spatial_analysis_stage!(static_only_model,
        :modal; state = staged_state, time = 0.0)
    @test !static_only_torque.active[]
    @test staged_state[static_only_torque.magnitude_variable] == 0.0
    @test staged_state[static_only_torque.global_torque_variables] ≈ zeros(3)

    derivative = SpatialSimulationRunner.initial_spatial_derivative(
        initial, loaded)
    selection = AnalysisSelection(Dynamics(), loaded.active_variable_indices,
        loaded.active_equation_indices)
    equations = zeros(length(selection.equation_indices))
    evaluate_analysis_equations!(equations, loaded.model, selection, 0.0,
        initial, derivative)
    @test norm(equations, Inf) < 1.0e-12

    coefficient = 2.2
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

    result = run_spatial_model(path; duration = 0.5, samples = 51)
    theta = hinge.rotation_variables[3]
    @test length(result.states) == 51
    @test last(result.states)[theta] < 0
    @test maximum(abs(state[torque.magnitude_variable] -
        (-8state[theta] - 0.24state[hinge.rotation_variables[2]]))
        for state in result.states) < 1.0e-8

    expression_source = replace(source,
        "[analysis]" => "[parameters]\nk = 8.0\nc = 0.24\n\n[analysis]",
        "stiffness = 8.0\ndamping_time_scale = 0.03\nfree_angle = \"0 deg\"" =>
            "expression = \"-k * pin.theta - c * pin.omega\"")
    expression_model = load_spatial_model(IOBuffer(expression_source))
    expression_torque = expression_model.forces[:spring]
    @test expression_model.initial_values[
        expression_torque.magnitude_variable] ≈
        initial[torque.magnitude_variable]

    constant_source = replace(source,
        "stiffness = 8.0\ndamping_time_scale = 0.03\nfree_angle = \"0 deg\"" =>
            "torque = 5.0")
    constant_model = load_spatial_model(IOBuffer(constant_source))
    constant_torque = constant_model.forces[:spring]
    @test constant_model.initial_values[constant_torque.magnitude_variable] ==
        5.0

    hinge_source = replace(read(joinpath(SPATIAL_MODEL_DIRECTORY,
            "hinge-pendulum.toml"), String),
        "rotation_coordinates = true\n" => "") * """

        [brake]
        type = "applied_torque"
        joint = "axis"
        torque = -0.1
        """
    hinge_model = load_spatial_model(IOBuffer(hinge_source))
    hinge_torque = hinge_model.forces[:brake]
    @test hinge_torque.hinge === hinge_model.connections[:axis]
    @test !isempty(hinge_torque.hinge.rotation_variables)

    initial_free_source = replace(source, "free_angle = \"0 deg\"" =>
        "free_angle = \"initial\"")
    initial_free_model = load_spatial_model(IOBuffer(initial_free_source))
    initial_free_torque = initial_free_model.forces[:spring]
    @test initial_free_torque.free_angle[] ≈
        initial_free_model.initial_values[
            initial_free_model.connections[:pin].hinge.rotation_variables[3]]
    @test abs(initial_free_model.initial_values[
        initial_free_torque.magnitude_variable]) < 1.0e-14

    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "damping_time_scale = 0.03" =>
            "damping = 0.24\ndamping_time_scale = 0.03")))
    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "joint = \"pin\"" => "joint = \"missing\"")))
end

@testset "Spatial rotational motion generator" begin
    path = joinpath(SPATIAL_MODEL_DIRECTORY,
        "constant-speed-revolute-crank.toml")
    source = read(path, String)
    loaded = load_spatial_model(path)
    driver = loaded.drivers[:drive]
    joint = loaded.connections[:pin]
    hinge = joint.hinge
    initial = loaded.initial_values

    @test driver isa SpatialRotationalMotionGenerator
    @test driver.hinge === hinge
    @test !isempty(hinge.rotation_variables)
    @test loaded.analysis.degrees_of_freedom == 0
    @test isempty(loaded.state_selection.selected_velocities)
    @test initial[hinge.rotation_variables] ≈ [0.0, 2pi, 0.0]
    @test initial[driver.reaction_variable] ≈ 4.905

    derivative = SpatialSimulationRunner.initial_spatial_derivative(
        initial, loaded)
    selection = AnalysisSelection(Dynamics(), loaded.active_variable_indices,
        loaded.active_equation_indices)
    equations = zeros(length(selection.equation_indices))
    evaluate_analysis_equations!(equations, loaded.model, selection, 0.0,
        initial, derivative)
    @test norm(equations, Inf) < 1.0e-12

    jacobian_state = copy(initial)
    jacobian_state[driver.reaction_variable] = 1.3
    coefficient = 2.1
    evaluate_analysis_equations!(equations, loaded.model, selection, 0.0,
        jacobian_state, derivative)
    analytical = Matrix(evaluate_analysis_sparse_jacobian(loaded.model,
        selection, 0.0, jacobian_state, derivative, coefficient))
    numerical = similar(analytical)
    step = 1.0e-7
    plus, minus = similar(equations), similar(equations)
    for (column, variable) in enumerate(loaded.active_variable_indices)
        state_plus, state_minus = copy(jacobian_state), copy(jacobian_state)
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

    result = run_spatial_model(path; duration = 0.1, samples = 11)
    @test last(result.states)[hinge.rotation_variables] ≈
        [0.0, 2pi, 0.2pi] atol = 2.0e-9
    @test last(result.states)[loaded.bodies[:crank].position_variables] ≈
        [0.5cos(0.2pi), 0.5sin(0.2pi), 0.0] atol = 2.0e-9

    expression_source = replace(source,
        "function = \"constant_speed\"\ninitial_angle = \"0 deg\"\n" *
        "angular_velocity = \"speed\"" =>
        "function = \"expression\"\nangle = \"0.25*sin(4*t)\"")
    expression_model = load_spatial_model(IOBuffer(expression_source))
    expression_hinge = expression_model.connections[:pin].hinge
    @test expression_model.initial_values[expression_hinge.rotation_variables] ≈
        [0.0, 1.0, 0.0]
    expression_result = run_spatial_model(IOBuffer(expression_source);
        duration = 0.05, samples = 6)
    @test last(expression_result.states)[expression_hinge.rotation_variables] ≈
        [-4sin(0.2), cos(0.2), 0.25sin(0.2)] atol = 2.0e-9
    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(
        expression_source, "angle = \"0.25*sin(4*t)\"" =>
            "angle = \"0.25*sin(4*t)\"\nangular_velocity = \"cos(4*t)\"")))

    shifted_source = replace(expression_source,
        "start_time = 0.0" => "start_time = 0.2")
    shifted_model = load_spatial_model(IOBuffer(shifted_source))
    shifted_hinge = shifted_model.connections[:pin].hinge
    @test shifted_model.initial_values[shifted_hinge.rotation_variables] ≈
        [-4sin(0.8), cos(0.8), 0.25sin(0.8)] atol = 2.0e-10

    hinge_source = replace(read(joinpath(SPATIAL_MODEL_DIRECTORY,
            "hinge-pendulum.toml"), String),
        "[state_selection]\nmethod = \"preferred\"\n" *
        "preferred_velocities = [\"axis.omega\"]\n" *
        "allow_fallback = false\n\n" => "",
        "rotation_coordinates = true\n" => "") * """

        [drive]
        type = "rotational_motion"
        joint = "axis"
        function = "constant_speed"
        angular_velocity = 0.5
        """
    hinge_model = load_spatial_model(IOBuffer(hinge_source))
    @test hinge_model.drivers[:drive].hinge ===
        hinge_model.connections[:axis]
    @test hinge_model.analysis.degrees_of_freedom == 0

    duplicate_driver = source * """

        [second_drive]
        type = "rotational_motion"
        joint = "pin"
        function = "constant_speed"
        angular_velocity = 1.0
        """
    @test_throws ArgumentError load_spatial_model(IOBuffer(duplicate_driver))
    bad_joint = replace(source, "joint = \"pin\"" => "joint = \"missing\"")
    @test_throws ArgumentError load_spatial_model(IOBuffer(bad_joint))
end

@testset "Spatial translational motion generator" begin
    path = joinpath(SPATIAL_MODEL_DIRECTORY,
        "constant-speed-translational-slider.toml")
    source = read(path, String)
    loaded = load_spatial_model(path)
    driver = loaded.drivers[:drive]
    body = loaded.bodies[:slider]
    initial = loaded.initial_values

    @test driver isa SpatialTranslationalMotionGenerator
    @test driver.geometry isa SpatialDirectedDistance
    @test driver.geometry.axis.index == 3
    @test loaded.analysis.degrees_of_freedom == 0
    @test isempty(loaded.state_selection.selected_velocities)
    @test length(component_variable_indices(loaded.layout, :drive)) == 4
    @test initial[[driver.distance_variable, driver.velocity_variable,
        driver.acceleration_variable, driver.reaction_variable]] ≈
        [0.2, 0.4, 0.0, 9.81]
    @test directed_distance_direction(driver.geometry, initial) ≈
        [0.0, 0.0, 1.0]

    derivative = SpatialSimulationRunner.initial_spatial_derivative(
        initial, loaded)
    @test derivative[driver.distance_variable] == 0.4
    @test derivative[driver.velocity_variable] == 0.0
    selection = AnalysisSelection(Dynamics(), loaded.active_variable_indices,
        loaded.active_equation_indices)
    equations = zeros(length(selection.equation_indices))
    evaluate_analysis_equations!(equations, loaded.model, selection, 0.0,
        initial, derivative)
    @test norm(equations, Inf) < 1.0e-12

    coefficient = 2.4
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

    result = run_spatial_model(path; duration = 0.1, samples = 6)
    final = last(result.states)
    @test final[[driver.distance_variable, driver.velocity_variable,
        driver.acceleration_variable, driver.reaction_variable]] ≈
        [0.24, 0.4, 0.0, 9.81] atol = 2.0e-10
    @test final[body.position_variables] ≈ [0.0, 0.0, 0.24] atol = 2.0e-10

    inconsistent_source = replace(source,
        "position = [0.0, 0.0, 0.2]" =>
            "position = [0.0, 0.0, -0.1]")
    corrected = load_spatial_model(IOBuffer(inconsistent_source))
    @test corrected.initial_values[
        corrected.bodies[:slider].position_variables] ≈
        [0.0, 0.0, 0.2] atol = 1.0e-12

    expression_source = replace(source,
        "function = \"constant_speed\"\ninitial_distance = 0.2\n" *
        "velocity = \"speed\"" =>
        "function = \"expression\"\ndistance = \"0.2 + 0.1*t^2\"")
    expression = run_spatial_model(IOBuffer(expression_source);
        duration = 0.1, samples = 3)
    expression_driver = expression.loaded.drivers[:drive]
    expression_final = last(expression.states)
    @test expression_final[[expression_driver.distance_variable,
        expression_driver.velocity_variable,
        expression_driver.acceleration_variable]] ≈
        [0.201, 0.02, 0.2] atol = 2.0e-10
    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(
        expression_source, "distance = \"0.2 + 0.1*t^2\"" =>
            "distance = \"0.2 + 0.1*t^2\"\nvelocity = \"0.2*t\"")))

    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "velocity = \"speed\"\n" => "")))
    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "[drive]\ntype = \"translational_motion\"\n" *
        "markers = [\"slider.axis\", \"ground.axis\"]" =>
        "[drive]\ntype = \"translational_motion\"\n" *
        "markers = [\"slider.axis\"]")))
    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "function = \"constant_speed\"" => "function = \"unknown\"")))
end

@testset "Spatial spanning motion generator" begin
    path = joinpath(SPATIAL_MODEL_DIRECTORY,
        "constant-distance-massless-link.toml")
    source = read(path, String)
    loaded = load_spatial_model(path)
    driver = loaded.drivers[:link]
    body = loaded.bodies[:body]
    initial = loaded.initial_values

    @test driver isa SpatialSpanningMotionGenerator
    @test loaded.analysis.degrees_of_freedom == 5
    @test loaded.state_selection.selected_velocities == [
        Symbol("body.V_x"), Symbol("body.V_y"), Symbol("body.omega_x"),
        Symbol("body.omega_y"), Symbol("body.omega_z")]
    @test length(component_variable_indices(loaded.layout, :link)) == 13
    @test initial[[driver.span.distance_variable,
        driver.span.velocity_variable, driver.span.acceleration_variable,
        driver.reaction_variable]] ≈
        [1.0, 0.0, 0.0, -6.246]
    @test initial[driver.global_force_variables] ≈
        [-4.9968, 0.0, 3.7476]
    @test norm(spatial_marker_position(driver.span.marker_2, initial) -
        spatial_marker_position(driver.span.marker_1, initial)) ≈ 1.0

    derivative = SpatialSimulationRunner.initial_spatial_derivative(
        initial, loaded)
    @test derivative[driver.span.distance_variable] == 0.0
    @test derivative[driver.span.velocity_variable] == 0.0
    selection = AnalysisSelection(Dynamics(), loaded.active_variable_indices,
        loaded.active_equation_indices)
    equations = zeros(length(selection.equation_indices))
    evaluate_analysis_equations!(equations, loaded.model, selection, 0.0,
        initial, derivative)
    @test norm(equations, Inf) < 1.0e-12

    coefficient = 2.3
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

    result = run_spatial_model(path; duration = 0.2, samples = 11)
    @test maximum(abs(state[driver.span.distance_variable] - 1.0)
        for state in result.states) < 2.0e-10
    @test maximum(abs(norm(spatial_marker_position(driver.span.marker_2, state) -
        spatial_marker_position(driver.span.marker_1, state)) - 1.0)
        for state in result.states) < 5.0e-8

    speed_source = replace(source,
        "distance = 1.0" =>
        "function = \"constant_speed\"\ninitial_distance = 1.0\n" *
        "velocity = -0.1")
    speed_result = run_spatial_model(IOBuffer(speed_source);
        duration = 0.1, samples = 3)
    speed_driver = speed_result.loaded.drivers[:link]
    @test last(speed_result.states)[[speed_driver.span.distance_variable,
        speed_driver.span.velocity_variable,
        speed_driver.span.acceleration_variable]] ≈
        [0.99, -0.1, 0.0] atol = 2.0e-10

    expression_source = replace(source,
        "distance = 1.0" =>
        "distance = \"1 + 0.05*sin(2*t)\"")
    expression_result = run_spatial_model(IOBuffer(expression_source);
        duration = 0.05, samples = 3)
    expression_driver = expression_result.loaded.drivers[:link]
    expression_final = last(expression_result.states)
    @test expression_final[[expression_driver.span.distance_variable,
        expression_driver.span.velocity_variable,
        expression_driver.span.acceleration_variable]] ≈
        [1 + 0.05sin(0.1), 0.1cos(0.1), -0.2sin(0.1)] atol = 2.0e-9
    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(
        expression_source, "distance = \"1 + 0.05*sin(2*t)\"" =>
            "distance = \"1 + 0.05*sin(2*t)\"\n" *
            "acceleration = \"-0.2*sin(2*t)\"")))

    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "distance = 1.0" => "distance = 0.0")))
    coincident_source = replace(source,
        "position = [0.8, 0.0, -0.6]" => "position = [0.0, 0.0, 0.0]")
    @test_throws ArgumentError load_spatial_model(IOBuffer(coincident_source))
end

@testset "Spatial TOML free rigid body" begin
    path = joinpath(SPATIAL_MODEL_DIRECTORY, "free-rotating-body.toml")
    source = read(path, String)
    loaded = load_spatial_model(path)
    streamed = load_spatial_model(IOBuffer(source))
    @test loaded.title == "Free rotating spatial body"
    @test streamed.title == loaded.title
    @test length(loaded.layout.catalog.variables) == 22
    @test length(loaded.layout.catalog.equations) == 22
    @test loaded.analysis.degrees_of_freedom == 6
    @test loaded.analysis.formulation == :StateSelected
    @test loaded.analysis.deficit == 0
    @test Set(keys(loaded.bodies)) == Set(keys(streamed.bodies))

    body = loaded.bodies[:body]
    initial = loaded.initial_values
    initial_angle = pi / 6
    initial_orientation = axis_angle_rotation(initial_angle, [0.0, 0.0, 1.0])
    @test rotation_matrix(initial[body.euler_parameter_variables]) ≈
        initial_orientation
    @test spatial_marker_position(loaded.markers[Symbol("body.reference")], initial) ≈
        initial_orientation * [0.4, 0.0, 0.0]
    @test initial[body.acceleration_variables] == [0.0, 0.0, -9.81]
    @test initial[body.angular_acceleration_variables] == zeros(3)

    derivative = SpatialSimulationRunner.initial_spatial_derivative(
        initial, loaded)
    selection = AnalysisSelection(Dynamics(), loaded.active_variable_indices,
        loaded.active_equation_indices)
    equations = zeros(length(initial))
    evaluate_analysis_equations!(equations, loaded.model, selection, 0.0,
        initial, derivative)
    @test norm(equations, Inf) < 1.0e-13

    coefficient = 3.7
    analytical = Matrix(evaluate_analysis_sparse_jacobian(loaded.model,
        selection, 0.0, initial, derivative, coefficient))
    numerical = similar(analytical)
    step = 1.0e-7
    plus = similar(equations)
    minus = similar(equations)
    for column in eachindex(initial)
        state_plus, state_minus = copy(initial), copy(initial)
        rate_plus, rate_minus = copy(derivative), copy(derivative)
        state_plus[column] += step
        state_minus[column] -= step
        rate_plus[column] += coefficient * step
        rate_minus[column] -= coefficient * step
        evaluate_analysis_equations!(plus, loaded.model, selection, 0.0,
            state_plus, rate_plus)
        evaluate_analysis_equations!(minus, loaded.model, selection, 0.0,
            state_minus, rate_minus)
        numerical[:, column] .= (plus .- minus) ./ (2step)
    end
    @test norm(analytical - numerical, Inf) < 2.0e-8

    result = run_spatial_model(path; duration = 0.4, samples = 21)
    final_time = last(result.times)
    final = last(result.states)
    expected_position = [0.0, 0.0, 0.0] .+
        [0.8, 0.2, 4.0] .* final_time .+
        0.5 .* [0.0, 0.0, -9.81] .* final_time^2
    expected_angle = initial_angle + 4final_time
    expected_parameters = [cos(expected_angle / 2), 0.0, 0.0,
        sin(expected_angle / 2)]
    actual_parameters = final[body.euler_parameter_variables]
    parameter_error = min(norm(actual_parameters - expected_parameters, Inf),
        norm(actual_parameters + expected_parameters, Inf))
    @test final[body.position_variables] ≈ expected_position atol = 2.0e-7
    @test final[body.velocity_variables] ≈
        [0.8, 0.2, 4.0] .+ [0.0, 0.0, -9.81] .* final_time atol = 2.0e-8
    @test final[body.angular_velocity_variables] ≈ [0.0, 0.0, 4.0]
    @test final[body.pseudo_angle_variables] ≈
        [0.0, 0.0, 4final_time] atol = 2.0e-8
    @test parameter_error < 2.0e-7
    @test abs(dot(actual_parameters, actual_parameters) - 1) < 2.0e-10
    @test rotation_matrix(actual_parameters) ≈
        [cos(expected_angle) -sin(expected_angle) 0.0;
         sin(expected_angle) cos(expected_angle) 0.0;
         0.0 0.0 1.0] atol = 4.0e-7
    @test count(result.solution.differential_vars) == 16
    @test count(result.solution.error_control) == 16

    interrupted_events = Any[]
    interrupt_after_accept = (time, values, rates, order, step, statistics) ->
        throw(InterruptException())
    @test_throws InterruptException begin
        SpatialSimulationRunner.run_spatial_implicit_model(loaded,
            range(0.0, 0.1; length = 3);
            accepted_step! = interrupt_after_accept,
            sample_progress = event -> push!(interrupted_events, event))
    end
    @test first(interrupted_events).kind == :sample
    @test last(interrupted_events).kind == :failure
    @test last(interrupted_events).time > 0.0
    @test length(last(interrupted_events).state) ==
        length(loaded.layout.catalog.variables)

    mktempdir() do directory
        result_path = joinpath(directory, "free-body.simp")
        write_result(result_path, result)
        stored = read_result(result_path)
        @test stored.analysis_mode == :dynamic
        @test stored.degrees_of_freedom == 6
        @test size(stored.values) == (21, 22)
        @test length(stored.state_selection_changes) == 1
        @test stored.state_selection_changes[1].reason == :initial
        @test TOML.parse(stored.model_source)["model"]["dimension"] == "spatial"
    end

    output = IOBuffer()
    @test spatial_model_main([path, "0.01", "2"];
        output, error = IOBuffer()) == 0
    @test occursin("analysis: dynamic (6 mechanical DOF)",
        String(take!(output)))

    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "dimension = \"spatial\"" => "dimension = \"planar\"")))
    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "inertia = [0.08, 0.12, 0.16]" => "inertia = [0.08, 0.12, 0.0]")))
    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "orientation = [0.5235987755982988, 0.0, 0.0, 1.0]" =>
            "orientation = [0.5235987755982988, 0.0, 0.0, 0.0]")))
    matrix_source = replace(source,
        "orientation = [0.5235987755982988, 0.0, 0.0, 1.0]" =>
            "orientation = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], " *
            "[0.0, 0.0, 1.0]]")
    matrix_loaded = load_spatial_model(IOBuffer(matrix_source))
    matrix_body = matrix_loaded.bodies[:body]
    @test rotation_matrix(matrix_loaded.initial_values[
        matrix_body.euler_parameter_variables]) ≈
        [0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0]
    @test_throws ArgumentError load_spatial_model(IOBuffer(source *
        "\n[unsupported]\ntype = \"revolute\"\n"))
end

@testset "Spatial initial-condition weighting and imposition" begin
    two_body_source = """
        [model]
        dimension = "spatial"

        [light]
        type = "rigid_body"
        mass = 1.0
        inertia = [0.1, 0.1, 0.1]
        position = [0.0, 0.0, 0.0]

        [light.point]
        type = "marker"

        [heavy]
        type = "rigid_body"
        mass = 9.0
        inertia = [0.9, 0.9, 0.9]
        position = [1.0, 0.0, 0.0]
        velocity = [1.0, 0.0, 0.0]

        [heavy.point]
        type = "marker"

        [joint]
        type = "spherical"
        markers = ["light.point", "heavy.point"]
        """
    weighted = load_spatial_model(IOBuffer(two_body_source))
    light = weighted.bodies[:light]
    heavy = weighted.bodies[:heavy]
    @test weighted.initial_values[light.position_variables[1]] ≈ 0.9
    @test weighted.initial_values[heavy.position_variables[1]] ≈ 0.9
    @test weighted.initial_values[light.velocity_variables[1]] ≈ 0.9
    @test weighted.initial_values[heavy.velocity_variables[1]] ≈ 0.9
    @test weighted.initial_condition_weights[:light].weight == 1.0
    @test weighted.initial_condition_weights[:heavy].weight == 9.0

    equally_weighted = load_spatial_model(IOBuffer(replace(two_body_source,
        "mass = 1.0\ninertia = [0.1, 0.1, 0.1]" =>
            "mass = 1.0\ninertia = [0.1, 0.1, 0.1]\nic_weight_scale = 9.0")))
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
    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(
        two_body_source, "mass = 1.0" =>
            "mass = 1.0\nic_weight_scale = 0.0"; count = 1)))

    free_source = read(joinpath(SPATIAL_MODEL_DIRECTORY,
        "free-rotating-body.toml"), String)
    imposed_body_source = replace(free_source, "[body.reference]" => """
        [body.initial.impose]
        R_x = 2.0
        orientation = ["40 deg", 0.0, 1.0, 0.0]
        V_z = -1.25
        omega_y = 0.75

        [body.reference]""")
    imposed_body = load_spatial_model(IOBuffer(imposed_body_source))
    imposed_free = imposed_body.bodies[:body]
    @test imposed_body.initial_values[imposed_free.position_variables[1]] == 2.0
    @test imposed_body.initial_values[imposed_free.velocity_variables[3]] == -1.25
    @test imposed_body.initial_values[
        imposed_free.angular_velocity_variables[2]] == 0.75
    @test rotation_matrix(imposed_body.initial_values[
        imposed_free.euler_parameter_variables]) ≈
        axis_angle_rotation(deg2rad(40), [0.0, 1.0, 0.0])
    @test imposed_body.state_selection.imposed_initial_variables == [
        Symbol("body.R_x"), Symbol("body.V_z"),
        Symbol("body.omega_y"), Symbol("body.orientation")]

    hinge_source = read(joinpath(SPATIAL_MODEL_DIRECTORY,
        "hinge-pendulum.toml"), String)
    imposed_angle_source = replace(hinge_source, "[gravity]" => """
        [axis.initial.impose]
        angle = "5 deg"

        [gravity]""")
    imposed_angle = load_spatial_model(IOBuffer(imposed_angle_source))
    imposed_hinge = imposed_angle.connections[:axis]
    imposed_theta = imposed_hinge.rotation_variables[3]
    imposed_omega = imposed_hinge.rotation_variables[2]
    @test imposed_angle.initial_values[imposed_theta] ≈ deg2rad(5)
    @test hinge_angle(imposed_hinge, imposed_angle.initial_values) ≈
        deg2rad(5)
    @test imposed_angle.initial_values[imposed_omega] ≈
        hinge_angular_velocity(imposed_hinge, imposed_angle.initial_values)
    @test imposed_angle.initial_values[imposed_omega] != 0.4
    @test imposed_angle.state_selection.imposed_initial_variables ==
        [Symbol("axis.theta")]

    imposed_omega_source = replace(hinge_source, "[gravity]" => """
        [axis.initial.impose]
        omega = 2.0

        [gravity]""")
    imposed_omega_model = load_spatial_model(IOBuffer(imposed_omega_source))
    omega_hinge = imposed_omega_model.connections[:axis]
    omega_variable = omega_hinge.rotation_variables[2]
    @test imposed_omega_model.initial_values[omega_variable] ≈ 2.0
    @test hinge_angular_velocity(omega_hinge,
        imposed_omega_model.initial_values) ≈ 2.0
    @test imposed_omega_model.state_selection.imposed_initial_variables ==
        [Symbol("axis.omega")]
    imposed_omega_run = run_spatial_model(IOBuffer(imposed_omega_source);
        duration = 0.01, samples = 2)
    run_hinge = imposed_omega_run.loaded.connections[:axis]
    @test first(imposed_omega_run.states)[run_hinge.rotation_variables[2]] ≈
        2.0

    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(
        imposed_angle_source, "rotation_coordinates = true" =>
            "rotation_coordinates = false")))

    doubly_imposed = replace(two_body_source,
        "[light.point]" => """
        [light.initial.impose]
        R_x = 0.0

        [light.point]""",
        "[heavy.point]" => """
        [heavy.initial.impose]
        R_x = 1.0

        [heavy.point]""")
    @test_throws ArgumentError load_spatial_model(IOBuffer(doubly_imposed))

    inline_source = read(joinpath(SPATIAL_MODEL_DIRECTORY,
        "inline-pendulum.toml"), String)
    imposed_distance_source = replace(inline_source, "[gravity]" => """
        [guide.initial.impose]
        distance = 0.2
        velocity = 0.7

        [gravity]""")
    imposed_distance = load_spatial_model(IOBuffer(imposed_distance_source))
    imposed_inline = imposed_distance.connections[:guide]
    acceleration, axial_velocity, distance =
        imposed_inline.translation_variables
    @test imposed_distance.initial_values[distance] ≈ 0.2
    @test inline_distance(imposed_inline,
        imposed_distance.initial_values) ≈ 0.2
    @test imposed_distance.initial_values[axial_velocity] ≈ 0.7
    @test inline_velocity(imposed_inline,
        imposed_distance.initial_values) ≈ 0.7
end

@testset "Spatial inplane constraint" begin
    path = joinpath(SPATIAL_MODEL_DIRECTORY, "inplane-slider.toml")
    source = read(path, String)
    loaded = load_spatial_model(path)
    body = loaded.bodies[:slider]
    constraint = loaded.connections[:support]
    initial = loaded.initial_values

    @test constraint isa SpatialInplaneConstraint
    @test constraint.geometry isa SpatialDirectedDistance
    @test length(loaded.layout.catalog.variables) == 23
    @test length(loaded.layout.catalog.equations) == 25
    @test length(loaded.active_variable_indices) == 23
    @test length(loaded.active_equation_indices) == 23
    @test loaded.analysis.degrees_of_freedom == 5
    @test loaded.state_selection.selected_velocities == [
        Symbol("slider.V_x"), Symbol("slider.V_y"),
        Symbol("slider.omega_x"), Symbol("slider.omega_y"),
        Symbol("slider.omega_z")]
    @test inplane_normal(constraint, initial) == [0.0, 0.0, 1.0]
    @test abs(inplane_position(constraint, initial)) < 1.0e-14
    @test abs(inplane_velocity(constraint, initial)) < 1.0e-14
    @test abs(inplane_acceleration(constraint, initial)) < 1.0e-13
    @test initial[constraint.reaction_variable] ≈ 9.81
    @test initial[body.acceleration_variables] ≈ zeros(3) atol = 1.0e-13
    @test initial[body.angular_acceleration_variables] ≈
        zeros(3) atol = 1.0e-13

    derivative = SpatialSimulationRunner.initial_spatial_derivative(
        initial, loaded)
    selection = AnalysisSelection(Dynamics(), loaded.active_variable_indices,
        loaded.active_equation_indices)
    equations = zeros(length(loaded.active_equation_indices))
    evaluate_analysis_equations!(equations, loaded.model, selection, 0.0,
        initial, derivative)
    @test norm(equations, Inf) < 2.0e-13

    coefficient = 2.7
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
    @test norm(analytical - numerical, Inf) < 5.0e-8

    progress_events = Any[]
    result = run_spatial_model(path; duration = 0.5, samples = 21,
        result_progress = event -> push!(progress_events, event))
    @test first(progress_events).kind == :begin
    streamed = filter(event -> event.kind == :sample, progress_events)
    @test getproperty.(streamed, :time) == result.times
    @test all(streamed[index].state ≈ result.states[index]
        for index in eachindex(streamed))
    @test maximum(abs(inplane_position(constraint, state))
        for state in result.states) < 1.0e-10
    final = last(result.states)
    @test final[body.position_variables] ≈ [0.225, 0.1, 0.1]
        atol = 2.0e-8
    @test final[body.velocity_variables] ≈ [0.45, 0.20, 0.0]
        atol = 2.0e-9
    @test maximum(abs(state[constraint.reaction_variable] - 9.81)
        for state in result.states) < 2.0e-8

    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "markers = [\"slider.contact\", \"ground.plane\"]" =>
            "markers = [\"ground.plane\", \"ground.plane\"]")))

    moving_plane = load_spatial_model(IOBuffer("""
        [model]
        dimension = "spatial"

        [plane]
        type = "rigid_body"
        mass = 2.0
        inertia = [0.4, 0.5, 0.6]
        angular_velocity = [0.2, -0.1, 0.3]

        [plane.surface]
        type = "marker"

        [slider]
        type = "rigid_body"
        mass = 1.0
        inertia = [0.1, 0.2, 0.3]
        position = [0.5, 0.2, 0.0]
        velocity = [0.1, 0.1, 0.09]

        [slider.contact]
        type = "marker"

        [contact]
        type = "inplane"
        markers = ["slider.contact", "plane.surface"]

        [gravity]
        type = "gravity"
        acceleration = [0.0, 0.0, -9.81]
        bodies = ["slider"]
        """))
    moving_constraint = moving_plane.connections[:contact]
    moving_initial = moving_plane.initial_values
    @test abs(moving_initial[moving_constraint.reaction_variable]) > 1
    moving_derivative = SpatialSimulationRunner.initial_spatial_derivative(
        moving_initial, moving_plane)
    moving_selection = AnalysisSelection(Dynamics(),
        moving_plane.active_variable_indices,
        moving_plane.active_equation_indices)
    moving_equations = zeros(length(moving_plane.active_equation_indices))
    moving_analytical = Matrix(evaluate_analysis_sparse_jacobian(
        moving_plane.model, moving_selection, 0.0, moving_initial,
        moving_derivative, coefficient))
    moving_numerical = similar(moving_analytical)
    moving_plus, moving_minus = similar(moving_equations),
        similar(moving_equations)
    for (column, variable) in enumerate(moving_plane.active_variable_indices)
        state_plus, state_minus = copy(moving_initial), copy(moving_initial)
        rate_plus, rate_minus = copy(moving_derivative),
            copy(moving_derivative)
        state_plus[variable] += step
        state_minus[variable] -= step
        rate_plus[variable] += coefficient * step
        rate_minus[variable] -= coefficient * step
        evaluate_analysis_equations!(moving_plus, moving_plane.model,
            moving_selection, 0.0, state_plus, rate_plus)
        evaluate_analysis_equations!(moving_minus, moving_plane.model,
            moving_selection, 0.0, state_minus, rate_minus)
        moving_numerical[:, column] .=
            (moving_plus .- moving_minus) ./ (2step)
    end
    moving_difference = moving_analytical - moving_numerical
    @test norm(moving_difference, Inf) < 2.0e-7
end

@testset "Spatial inline constraint" begin
    path = joinpath(SPATIAL_MODEL_DIRECTORY, "inline-slider.toml")
    source = read(path, String)
    loaded = load_spatial_model(path)
    body = loaded.bodies[:slider]
    constraint = loaded.connections[:guide]
    initial = loaded.initial_values

    @test constraint isa SpatialInlineConstraint
    @test constraint.inplane_x.geometry.axis.index == 1
    @test constraint.inplane_y.geometry.axis.index == 2
    @test constraint.axial_geometry isa SpatialDirectedDistance
    @test inline_axis(constraint, initial) == [0.0, 0.0, 1.0]
    @test length(loaded.layout.catalog.variables) == 27
    @test length(loaded.layout.catalog.equations) == 33
    @test length(loaded.active_variable_indices) == 27
    @test length(loaded.active_equation_indices) == 27
    @test loaded.analysis.degrees_of_freedom == 4
    @test loaded.state_selection.selected_velocities == [
        Symbol("guide.velocity"), Symbol("slider.omega_x"),
        Symbol("slider.omega_y"), Symbol("slider.omega_z")]

    partial_preference_source = replace(source,
        "preferred_velocities = [\n" *
        "    \"guide.velocity\",\n" *
        "    \"slider.omega_x\",\n" *
        "    \"slider.omega_y\",\n" *
        "    \"slider.omega_z\",\n" *
        "]" => "preferred_velocities = [\"guide.velocity\"]")
    partial_preference = load_spatial_model(
        IOBuffer(partial_preference_source))
    @test length(partial_preference.state_selection.selected_velocities) == 4
    @test first(partial_preference.state_selection.selected_velocities) ==
        Symbol("guide.velocity")
    @test !partial_preference.state_selection.fallback_used
    @test rotation_matrix(initial[body.euler_parameter_variables]) ≈
        axis_angle_rotation(pi / 6, [0.0, 0.0, 1.0])
    @test initial[constraint.translation_variables] ≈ [-9.81, 0.0, 1.0]
    @test abs(inplane_position(constraint.inplane_x, initial)) < 1.0e-14
    @test abs(inplane_position(constraint.inplane_y, initial)) < 1.0e-14

    derivative = SpatialSimulationRunner.initial_spatial_derivative(
        initial, loaded)
    selection = AnalysisSelection(Dynamics(), loaded.active_variable_indices,
        loaded.active_equation_indices)
    equations = zeros(length(loaded.active_equation_indices))
    evaluate_analysis_equations!(equations, loaded.model, selection, 0.0,
        initial, derivative)
    @test norm(equations, Inf) < 2.0e-13

    coefficient = 2.4
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
    @test norm(analytical - numerical, Inf) < 5.0e-8

    result = run_spatial_model(path)
    final = last(result.states)
    @test final[constraint.translation_variables] ≈
        [-9.81, -9.81, 1.0 - 9.81 / 2] atol = 5.0e-7
    @test isapprox(final[body.position_variables],
        [0.0, 0.0, 1.0 - 9.81 / 2]; atol = 5.0e-7)
    @test maximum(max(abs(inplane_position(constraint.inplane_x, state)),
                      abs(inplane_position(constraint.inplane_y, state)))
        for state in result.states) < 1.0e-10

    without_coordinates = replace(source,
        "[state_selection]\nmethod = \"preferred\"\n" *
        "preferred_velocities = [\n" *
        "    \"guide.velocity\",\n" *
        "    \"slider.omega_x\",\n" *
        "    \"slider.omega_y\",\n" *
        "    \"slider.omega_z\",\n" *
        "]\nallow_fallback = false\n\n" => "",
        "translation_coordinates = true\n" => "")
    plain = load_spatial_model(IOBuffer(without_coordinates))
    @test isempty(plain.connections[:guide].translation_variables)
    @test length(plain.layout.catalog.variables) == 24
    @test length(plain.layout.catalog.equations) == 28
    @test plain.analysis.degrees_of_freedom == 4
    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "translation_coordinates = true" =>
            "translation_coordinates = \"yes\"")))
    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "\"30°\"" => "\"thirty°\"")))
end

@testset "Spatial slider-crank" begin
    path = joinpath(SPATIAL_MODEL_DIRECTORY, "slider-crank.toml")
    loaded = load_spatial_model(path)
    ground_pin = loaded.connections[:ground_pin]
    connecting_pin = loaded.connections[:connecting_pin]
    slider_guide = loaded.connections[:slider_guide]
    initial = loaded.initial_values

    @test ground_pin isa SpatialRevoluteJoint
    @test connecting_pin isa SpatialRevoluteJoint
    @test slider_guide isa SpatialInplaneConstraint
    @test length(loaded.layout.catalog.variables) == 58
    @test length(loaded.layout.catalog.equations) == 82
    @test length(loaded.active_variable_indices) == 58
    @test length(loaded.active_equation_indices) == 58
    @test loaded.analysis.degrees_of_freedom == 1
    @test loaded.state_selection.selected_velocities ==
        [Symbol("ground_pin.omega")]
    @test initial[ground_pin.hinge.rotation_variables[3]] ≈ pi / 4
    @test inplane_normal(slider_guide, initial) ≈ [0.0, 1.0, 0.0]
        atol = 1.0e-14
    @test abs(inplane_position(slider_guide, initial)) < 1.0e-13

    result = run_spatial_model(path; duration = 0.5, samples = 21)
    @test last(result.states)[ground_pin.hinge.rotation_variables[3]] < pi / 4
    @test maximum(abs(inplane_position(slider_guide, state))
        for state in result.states) < 2.0e-8
    for joint in (ground_pin, connecting_pin)
        @test maximum(norm(spatial_marker_position(
            joint.spherical.marker_a, state) - spatial_marker_position(
            joint.spherical.marker_b, state))
            for state in result.states) < 4.0e-8
        @test maximum(abs(perp_position(joint.hinge.perp_xz, state))
            for state in result.states) < 2.0e-8
        @test maximum(abs(perp_position(joint.hinge.perp_yz, state))
            for state in result.states) < 2.0e-8
    end
end

@testset "Spatial revolute joint" begin
    path = joinpath(SPATIAL_MODEL_DIRECTORY, "revolute-pendulum.toml")
    source = read(path, String)
    loaded = load_spatial_model(path)
    body = loaded.bodies[:pendulum]
    joint = loaded.connections[:pin]
    initial = loaded.initial_values

    @test joint isa SpatialRevoluteJoint
    @test connection_hinge(joint) === joint.hinge
    @test joint.spherical.marker_a === joint.hinge.marker_i
    @test joint.spherical.marker_b === joint.hinge.marker_j
    @test length(loaded.layout.catalog.variables) == 30
    @test length(loaded.layout.catalog.equations) == 42
    @test length(loaded.active_variable_indices) == 30
    @test length(loaded.active_equation_indices) == 30
    @test loaded.analysis.degrees_of_freedom == 1
    @test loaded.state_selection.selected_velocities == [Symbol("pin.omega")]
    @test length(joint.spherical.reaction_variables) == 3
    @test joint.hinge.perp_xz.reaction_variable !=
        joint.hinge.perp_yz.reaction_variable
    alpha, omega, theta = joint.hinge.rotation_variables
    @test initial[theta] ≈ 0.0 atol = 1.0e-14
    @test initial[omega] ≈ 0.4 atol = 1.0e-14
    @test initial[alpha] ≈ hinge_angular_acceleration(
        joint.hinge, initial)

    derivative = SpatialSimulationRunner.initial_spatial_derivative(
        initial, loaded)
    selection = AnalysisSelection(Dynamics(), loaded.active_variable_indices,
        loaded.active_equation_indices)
    equations = zeros(length(loaded.active_equation_indices))
    evaluate_analysis_equations!(equations, loaded.model, selection, 0.0,
        initial, derivative)
    @test norm(equations, Inf) < 2.0e-13

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
    @test norm(analytical - numerical, Inf) < 1.0e-7

    result = run_spatial_model(path; duration = 0.5, samples = 21)
    @test result.solution.stats.accepted_steps < 100
    @test maximum(norm(spatial_marker_position(
        joint.spherical.marker_a, state) - spatial_marker_position(
        joint.spherical.marker_b, state)) for state in result.states) < 2.0e-8
    @test maximum(abs(perp_position(joint.hinge.perp_xz, state))
        for state in result.states) < 2.0e-8
    @test maximum(abs(perp_position(joint.hinge.perp_yz, state))
        for state in result.states) < 2.0e-8
    @test maximum(abs(state[omega] -
        hinge_angular_velocity(joint.hinge, state))
        for state in result.states) < 1.0e-7
    @test abs(last(result.states)[theta] - initial[theta]) > 1.0

    without_coordinates = replace(source,
        "[state_selection]\nmethod = \"preferred\"\n" *
            "preferred_velocities = [\"pin.omega\"]\n" *
            "allow_fallback = false\n\n" => "",
        "rotation_coordinates = true\n" => "")
    plain = load_spatial_model(IOBuffer(without_coordinates))
    @test plain.connections[:pin] isa SpatialRevoluteJoint
    @test isempty(plain.connections[:pin].hinge.rotation_variables)
    @test length(plain.layout.catalog.variables) == 27
    @test length(plain.layout.catalog.equations) == 37
    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "rotation_coordinates = true" => "rotation_coordinates = \"yes\"")))
end


@testset "Spatial spherical joint" begin
    path = joinpath(SPATIAL_MODEL_DIRECTORY, "spherical-pendulum.toml")
    source = read(path, String)
    loaded = load_spatial_model(path)
    body = loaded.bodies[:pendulum]
    joint = loaded.connections[:pin]
    initial = loaded.initial_values

    @test joint isa SpatialSphericalJoint
    @test length(loaded.layout.catalog.variables) == 25
    @test length(loaded.layout.catalog.equations) == 31
    @test length(loaded.active_equation_indices) == 25
    @test loaded.analysis.degrees_of_freedom == 3
    @test loaded.state_selection.selected_velocities == [
        Symbol("pendulum.omega_x"), Symbol("pendulum.omega_y"),
        Symbol("pendulum.omega_z")]
    @test spatial_marker_position(joint.marker_a, initial) ≈
        spatial_marker_position(joint.marker_b, initial)
    @test initial[body.acceleration_variables] ≈ [0.0, 0.0, -7.3575]
    @test initial[body.angular_acceleration_variables] ≈
        [0.0, 14.715, 0.0]
    @test initial[joint.reaction_variables] ≈ [0.0, 0.0, 2.4525]

    derivative = SpatialSimulationRunner.initial_spatial_derivative(
        initial, loaded)
    selection = AnalysisSelection(Dynamics(), loaded.active_variable_indices,
        loaded.active_equation_indices)
    equations = zeros(length(loaded.active_equation_indices))
    evaluate_analysis_equations!(equations, loaded.model, selection, 0.0,
        initial, derivative)
    @test norm(equations, Inf) < 2.0e-13

    coefficient = 2.3
    analytical = Matrix(evaluate_analysis_sparse_jacobian(loaded.model,
        selection, 0.0, initial, derivative, coefficient))
    numerical = similar(analytical)
    step = 1.0e-7
    plus = similar(equations)
    minus = similar(equations)
    for column in eachindex(initial)
        state_plus, state_minus = copy(initial), copy(initial)
        rate_plus, rate_minus = copy(derivative), copy(derivative)
        state_plus[column] += step
        state_minus[column] -= step
        rate_plus[column] += coefficient * step
        rate_minus[column] -= coefficient * step
        evaluate_analysis_equations!(plus, loaded.model, selection, 0.0,
            state_plus, rate_plus)
        evaluate_analysis_equations!(minus, loaded.model, selection, 0.0,
            state_minus, rate_minus)
        numerical[:, column] .= (plus .- minus) ./ (2step)
    end
    @test norm(analytical - numerical, Inf) < 5.0e-8

    workspace = SpatialSimulationRunner.spatial_dynamic_jacobian_workspace(
        loaded, selection, 0.0, initial, derivative, coefficient)
    dynamic_jacobian = copy(workspace.prototype)
    SpatialSimulationRunner.evaluate_spatial_dynamic_jacobian!(
        dynamic_jacobian, workspace, loaded, selection, 0.0, initial,
        derivative, coefficient)
    local_variables = Dict(index => local_index
        for (local_index, index) in enumerate(loaded.active_variable_indices))
    local_equations = Dict(index => local_index
        for (local_index, index) in enumerate(loaded.active_equation_indices))
    parameter_columns = [local_variables[index]
        for index in body.euler_parameter_variables]
    pseudo_columns = [local_variables[index]
        for index in body.pseudo_angle_variables]
    position_rows = [local_equations[index]
        for index in joint.position_equations]
    orientation_rows = [local_equations[index]
        for index in body.orientation_equations[1:3]]
    @test norm(dynamic_jacobian[position_rows, parameter_columns], Inf) == 0
    @test norm(dynamic_jacobian[position_rows, pseudo_columns], Inf) > 0
    @test dynamic_jacobian[orientation_rows, pseudo_columns] ≈
        coefficient .* Matrix{Float64}(I, 3, 3)
    @test norm(dynamic_jacobian[orientation_rows, parameter_columns], Inf) > 0

    result = run_spatial_model(path; duration = 0.5, samples = 21)
    @test count(result.solution.differential_vars) == 10
    @test count(result.solution.error_control) == 16
    dynamic_state = last(result.states)
    dynamic_derivative = SpatialSimulationRunner.initial_spatial_derivative(
        dynamic_state, loaded)
    analytical .= Matrix(evaluate_analysis_sparse_jacobian(loaded.model,
        selection, 0.5, dynamic_state, dynamic_derivative, coefficient))
    for column in eachindex(dynamic_state)
        state_plus, state_minus = copy(dynamic_state), copy(dynamic_state)
        rate_plus, rate_minus = copy(dynamic_derivative),
            copy(dynamic_derivative)
        state_plus[column] += step
        state_minus[column] -= step
        rate_plus[column] += coefficient * step
        rate_minus[column] -= coefficient * step
        evaluate_analysis_equations!(plus, loaded.model, selection, 0.5,
            state_plus, rate_plus)
        evaluate_analysis_equations!(minus, loaded.model, selection, 0.5,
            state_minus, rate_minus)
        numerical[:, column] .= (plus .- minus) ./ (2step)
    end
    @test norm(analytical - numerical, Inf) < 1.0e-7
    initial_energy = 0.5body.mass * dot(
        initial[body.velocity_variables], initial[body.velocity_variables]) +
        0.5dot(initial[body.angular_velocity_variables],
            body.inertia * initial[body.angular_velocity_variables]) +
        body.mass * 9.81 * initial[body.position_variables[3]]
    for state in result.states
        @test spatial_marker_position(joint.marker_a, state) ≈
            spatial_marker_position(joint.marker_b, state) atol = 2.0e-8
        @test spatial_marker_velocity(joint.marker_a, state) ≈
            spatial_marker_velocity(joint.marker_b, state) atol = 3.0e-7
        @test spatial_marker_acceleration(joint.marker_a, state) ≈
            spatial_marker_acceleration(joint.marker_b, state) atol = 5.0e-6
        energy = 0.5body.mass * dot(state[body.velocity_variables],
            state[body.velocity_variables]) +
            0.5dot(state[body.angular_velocity_variables],
                body.inertia * state[body.angular_velocity_variables]) +
            body.mass * 9.81 * state[body.position_variables[3]]
        @test energy ≈ initial_energy atol = 2.0e-6
    end

    displaced = load_spatial_model(IOBuffer(replace(source,
        "position = [0.5, 0.0, 0.0]" =>
            "position = [0.6, 0.0, 0.0]")))
    displaced_joint = displaced.connections[:pin]
    @test spatial_marker_position(displaced_joint.marker_a,
        displaced.initial_values) ≈ spatial_marker_position(
            displaced_joint.marker_b, displaced.initial_values) atol = 1.0e-14
    @test displaced.initial_conditions.position_corrections == 1

    ground_joint = source * """

        [ground.other]
        type = "marker"

        [invalid]
        type = "spherical"
        markers = ["ground.pin", "ground.other"]
        """
    @test_throws ArgumentError load_spatial_model(IOBuffer(ground_joint))
end

@testset "Spatial perpendicular-axis constraint" begin
    path = joinpath(SPATIAL_MODEL_DIRECTORY,
        "spherical-perp-pendulum.toml")
    source = read(path, String)
    loaded = load_spatial_model(path)
    body = loaded.bodies[:pendulum]
    constraint = loaded.connections[:guide]
    initial = loaded.initial_values

    @test constraint isa SpatialPerpConstraint
    @test length(loaded.layout.catalog.variables) == 26
    @test length(loaded.layout.catalog.equations) == 34
    @test length(loaded.active_variable_indices) == 26
    @test length(loaded.active_equation_indices) == 26
    @test loaded.analysis.degrees_of_freedom == 2
    @test loaded.state_selection.selected_velocities == [
        Symbol("pendulum.omega_y"), Symbol("pendulum.omega_x")]
    @test isempty(loaded.state_selection.inactive_variables)
    @test perp_normal(constraint, initial) ≈ [0.0, 0.0, 1.0]
    @test abs(perp_position(constraint, initial)) < 1.0e-14
    @test abs(perp_velocity(constraint, initial)) < 1.0e-14
    @test abs(perp_acceleration(constraint, initial)) < 1.0e-13

    derivative = SpatialSimulationRunner.initial_spatial_derivative(
        initial, loaded)
    selection = AnalysisSelection(Dynamics(), loaded.active_variable_indices,
        loaded.active_equation_indices)
    equations = zeros(length(loaded.active_equation_indices))
    evaluate_analysis_equations!(equations, loaded.model, selection, 0.0,
        initial, derivative)
    @test norm(equations, Inf) < 2.0e-13

    coefficient = 2.7
    analytical = Matrix(evaluate_analysis_sparse_jacobian(loaded.model,
        selection, 0.0, initial, derivative, coefficient))
    numerical = similar(analytical)
    step = 1.0e-7
    plus = similar(equations)
    minus = similar(equations)
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
    @test norm(analytical - numerical, Inf) < 5.0e-8

    result = run_spatial_model(path; duration = 0.5, samples = 21)
    @test count(result.solution.differential_vars) == 8
    @test count(result.solution.error_control) == 15
    for state in result.states
        @test abs(perp_position(constraint, state)) < 2.0e-8
        @test abs(perp_velocity(constraint, state)) < 3.0e-8
        @test abs(perp_acceleration(constraint, state)) < 2.0e-6
        @test norm(perp_normal(constraint, state)) ≈ 1.0 atol = 2.0e-8
    end

    dynamic_state = last(result.states)
    dynamic_derivative = SpatialSimulationRunner.initial_spatial_derivative(
        dynamic_state, loaded)
    analytical .= Matrix(evaluate_analysis_sparse_jacobian(loaded.model,
        selection, 0.5, dynamic_state, dynamic_derivative, coefficient))
    for (column, variable) in enumerate(loaded.active_variable_indices)
        state_plus, state_minus = copy(dynamic_state), copy(dynamic_state)
        rate_plus, rate_minus = copy(dynamic_derivative),
            copy(dynamic_derivative)
        state_plus[variable] += step
        state_minus[variable] -= step
        rate_plus[variable] += coefficient * step
        rate_minus[variable] -= coefficient * step
        evaluate_analysis_equations!(plus, loaded.model, selection, 0.5,
            state_plus, rate_plus)
        evaluate_analysis_equations!(minus, loaded.model, selection, 0.5,
            state_minus, rate_minus)
        numerical[:, column] .= (plus .- minus) ./ (2step)
    end
    @test abs(dynamic_state[constraint.reaction_variable]) > 1.0e-4
    @test norm(analytical - numerical, Inf) < 1.0e-7

    ground_constraint = source * """

        [ground.other]
        type = "marker"

        [invalid]
        type = "perp"
        markers = ["ground.pin", "ground.other"]
        """
    @test_throws ArgumentError load_spatial_model(IOBuffer(ground_constraint))

    misaligned_source = replace(source,
        "angular_velocity = [1.5, 0.0, 0.0]" =>
            "orientation = [0.2, 0.0, 0.0, 1.0]\n" *
            "angular_velocity = [1.5, 0.0, 0.0]")
    misaligned = load_spatial_model(IOBuffer(misaligned_source))
    misaligned_body = misaligned.bodies[:pendulum]
    misaligned_constraint = misaligned.connections[:guide]
    corrected = misaligned.initial_values
    corrected_orientation = rotation_matrix(
        corrected[misaligned_body.euler_parameter_variables])
    @test misaligned.initial_conditions.position_corrections == 3
    @test abs(perp_position(misaligned_constraint, corrected)) < 1.0e-12
    @test spatial_marker_position(
        misaligned.connections[:pin].marker_a, corrected) ≈
        spatial_marker_position(
            misaligned.connections[:pin].marker_b, corrected) atol = 1.0e-12
    @test transpose(corrected_orientation) * corrected_orientation ≈ I atol = 1.0e-13
    @test det(corrected_orientation) ≈ 1.0 atol = 1.0e-13
    @test norm(corrected_orientation - axis_angle_rotation(
        0.2, [0.0, 0.0, 1.0]), Inf) > 0.1
end

@testset "Spatial hinge constraint" begin
    path = joinpath(SPATIAL_MODEL_DIRECTORY, "hinge-pendulum.toml")
    source = read(path, String)
    loaded = load_spatial_model(path)
    body = loaded.bodies[:pendulum]
    spherical = loaded.connections[:pin]
    hinge = loaded.connections[:axis]
    initial = loaded.initial_values

    @test hinge isa SpatialHingeConstraint
    @test (hinge.perp_xz.axis_i, hinge.perp_xz.axis_j) == (1, 3)
    @test (hinge.perp_yz.axis_i, hinge.perp_yz.axis_j) == (2, 3)
    @test length(loaded.layout.catalog.variables) == 30
    @test length(loaded.layout.catalog.equations) == 42
    @test length(loaded.active_variable_indices) == 30
    @test length(loaded.active_equation_indices) == 30
    @test loaded.analysis.degrees_of_freedom == 1
    @test loaded.state_selection.selected_velocities ==
        [Symbol("axis.omega")]
    @test loaded.state_selection.method == :preferred
    @test !loaded.state_selection.fallback_used
    partition = SpatialSimulationRunner.runtime_spatial_state_partition(loaded)
    @test partition.selected_names == [Symbol("axis.omega")]
    @test length(partition.equation_indices) ==
        length(loaded.active_variable_indices)
    automatic_choice, runtime_diagnostics =
        SpatialSimulationRunner.automatic_runtime_spatial_state_selection(
            loaded, initial, 0.0)
    @test automatic_choice == [Symbol("pendulum.V_y")]
    @test runtime_diagnostics.rank == loaded.state_selection.rank
    SpatialSimulationRunner.select_runtime_spatial_states!(
        partition, automatic_choice)
    @test partition.selected_names == automatic_choice
    @test length(partition.equation_indices) ==
        length(loaded.active_variable_indices)
    runtime_selection = AnalysisSelection(Dynamics(),
        loaded.active_variable_indices, partition.equation_indices)
    runtime_equations = zeros(length(partition.equation_indices))
    runtime_derivative = SpatialSimulationRunner.initial_spatial_derivative(
        initial, loaded)
    evaluate_analysis_equations!(runtime_equations, loaded.model,
        runtime_selection, 0.0, initial, runtime_derivative)
    @test norm(runtime_equations, Inf) < 2.0e-13
    runtime_prototype = evaluate_analysis_sparse_jacobian(loaded.model,
        runtime_selection, 0.0, initial, runtime_derivative, 1.0)
    @test size(runtime_prototype) ==
        (length(loaded.active_variable_indices),
         length(loaded.active_variable_indices))
    @test perp_normal(hinge.perp_xz, initial) ≈ [0.0, -1.0, 0.0]
    @test perp_normal(hinge.perp_yz, initial) ≈ [1.0, 0.0, 0.0]
    alpha, omega, theta = hinge.rotation_variables
    @test initial[theta] ≈ 0.0 atol = 1.0e-14
    @test initial[omega] ≈ hinge_angular_velocity(hinge, initial)
    @test initial[alpha] ≈ hinge_angular_acceleration(hinge, initial)

    derivative = SpatialSimulationRunner.initial_spatial_derivative(
        initial, loaded)
    selection = AnalysisSelection(Dynamics(), loaded.active_variable_indices,
        loaded.active_equation_indices)
    equations = zeros(length(loaded.active_equation_indices))
    evaluate_analysis_equations!(equations, loaded.model, selection, 0.0,
        initial, derivative)
    @test norm(equations, Inf) < 2.0e-13

    jacobian_state = copy(initial)
    jacobian_state[hinge.perp_xz.reaction_variable] = 1.2
    jacobian_state[hinge.perp_yz.reaction_variable] = -0.7
    changed_equations = similar(equations)
    evaluate_analysis_equations!(changed_equations, loaded.model, selection,
        0.0, jacobian_state, derivative)
    torque_rows = [findfirst(==(row), loaded.active_equation_indices)
        for row in body.balance_equations[4:6]]
    @test changed_equations[torque_rows] - equations[torque_rows] ≈
        [0.7, 1.2, 0.0]

    coefficient = 2.5
    analytical = Matrix(evaluate_analysis_sparse_jacobian(loaded.model,
        selection, 0.0, jacobian_state, derivative, coefficient))
    numerical = similar(analytical)
    step = 1.0e-7
    plus = similar(equations)
    minus = similar(equations)
    for (column, variable) in enumerate(loaded.active_variable_indices)
        state_plus, state_minus = copy(jacobian_state), copy(jacobian_state)
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
    @test norm(analytical - numerical, Inf) < 1.0e-7

    result = run_spatial_model(path; duration = 0.5, samples = 21)
    @test length(result.state_selection_changes) == 1
    @test result.state_selection_changes[1].reason == :initial
    @test result.state_selection_changes[1].selected == [Symbol("axis.omega")]
    @test result.solution.stats.state_reselections == 0
    @test count(result.solution.differential_vars) == 6
    @test count(result.solution.error_control) == 15
    @test maximum(abs(perp_position(hinge.perp_xz, state))
        for state in result.states) < 2.0e-8
    @test maximum(abs(perp_position(hinge.perp_yz, state))
        for state in result.states) < 2.0e-8
    @test maximum(abs(perp_velocity(hinge.perp_xz, state))
        for state in result.states) < 2.0e-8
    @test maximum(abs(perp_velocity(hinge.perp_yz, state))
        for state in result.states) < 2.0e-8
    @test maximum(norm(spatial_marker_position(spherical.marker_a, state) -
        spatial_marker_position(spherical.marker_b, state))
        for state in result.states) < 2.0e-8
    initial_energy = 0.5body.mass * dot(initial[body.velocity_variables],
        initial[body.velocity_variables]) +
        0.5dot(initial[body.angular_velocity_variables],
            body.inertia * initial[body.angular_velocity_variables]) +
        body.mass * 9.81 * initial[body.position_variables[2]]
    energy_error = maximum(abs(
        0.5body.mass * dot(state[body.velocity_variables],
            state[body.velocity_variables]) +
        0.5dot(state[body.angular_velocity_variables],
            body.inertia * state[body.angular_velocity_variables]) +
        body.mass * 9.81 * state[body.position_variables[2]] - initial_energy)
        for state in result.states)
    @test energy_error < 2.0e-6
    @test maximum(abs(state[omega] - hinge_angular_velocity(hinge, state))
        for state in result.states) < 1.0e-7
    @test maximum(abs(state[alpha] -
        hinge_angular_acceleration(hinge, state))
        for state in result.states) < 2.0e-6

    misaligned_source = replace(source,
        "angular_velocity = [0.0, 0.0, 0.4]" =>
            "orientation = [0.15, 1.0, 0.0, 0.0]\n" *
            "angular_velocity = [0.0, 0.0, 0.4]")
    misaligned = load_spatial_model(IOBuffer(misaligned_source))
    corrected_hinge = misaligned.connections[:axis]
    @test misaligned.initial_conditions.position_corrections == 3
    @test abs(perp_position(corrected_hinge.perp_xz,
        misaligned.initial_values)) < 1.0e-12
    @test abs(perp_position(corrected_hinge.perp_yz,
        misaligned.initial_values)) < 1.0e-12

    rotating_source = replace(source,
        "angular_velocity = [0.0, 0.0, 0.4]" =>
            "angular_velocity = [0.0, 0.0, 4.0]",
        "velocity = [0.0, 0.2, 0.0]" =>
            "velocity = [0.0, 2.0, 0.0]",
        "acceleration = [0.0, -9.81, 0.0]" =>
            "acceleration = [0.0, 0.0, 0.0]")
    rotating = run_spatial_model(IOBuffer(rotating_source);
        duration = 1.7, samples = 18)
    rotating_hinge = rotating.loaded.connections[:axis]
    rotating_theta = rotating_hinge.rotation_variables[3]
    @test last(rotating.states)[rotating_theta] ≈ 6.8 atol = 2.0e-6
    @test hinge_angle(rotating_hinge, last(rotating.states)) ≈
        atan(sin(6.8), cos(6.8)) atol = 2.0e-6

    ground_hinge = source * """

        [ground.other]
        type = "marker"

        [invalid]
        type = "hinge"
        markers = ["ground.pin", "ground.other"]
        """
    @test_throws ArgumentError load_spatial_model(IOBuffer(ground_hinge))
end

@testset "Spatial orient constraint" begin
    path = joinpath(SPATIAL_MODEL_DIRECTORY, "oriented-free-body.toml")
    source = read(path, String)
    loaded = load_spatial_model(path)
    constraint = loaded.connections[:orientation_lock]
    body = loaded.bodies[:locked_body]
    free_body = loaded.bodies[:free_body]
    initial = loaded.initial_values

    @test constraint isa SpatialOrientConstraint
    @test length(loaded.layout.catalog.variables) == 47
    @test length(loaded.layout.catalog.equations) == 53
    @test length(loaded.active_variable_indices) == 47
    @test length(loaded.active_equation_indices) == 47
    @test loaded.analysis.degrees_of_freedom == 9
    @test loaded.initial_conditions.position_corrections > 0
    @test loaded.initial_conditions.velocity_corrections > 0
    @test norm(spatial_marker_orientation(constraint.marker_i, initial) -
        spatial_marker_orientation(constraint.marker_j, initial)) < 2.0e-12
    @test norm(initial[body.angular_velocity_variables]) < 1.0e-13

    primitives = (constraint.hinge.perp_xz,
        constraint.hinge.perp_yz, constraint.perp_xy)
    @test maximum(abs(perp_position(item, initial))
        for item in primitives) < 5.0e-15
    @test maximum(abs(perp_velocity(item, initial))
        for item in primitives) < 5.0e-15

    derivative = SpatialSimulationRunner.initial_spatial_derivative(
        initial, loaded)
    selection = AnalysisSelection(Dynamics(), loaded.active_variable_indices,
        loaded.active_equation_indices)
    jacobian_state = copy(initial)
    for (primitive, reaction) in zip(primitives, (1.1, -0.7, 0.4))
        jacobian_state[primitive.reaction_variable] = reaction
    end
    equations = zeros(length(selection.equation_indices))
    evaluate_analysis_equations!(equations, loaded.model, selection, 0.0,
        jacobian_state, derivative)
    coefficient = 2.4
    analytical = Matrix(evaluate_analysis_sparse_jacobian(loaded.model,
        selection, 0.0, jacobian_state, derivative, coefficient))
    numerical = similar(analytical)
    step = 1.0e-7
    plus, minus = similar(equations), similar(equations)
    for (column, variable) in enumerate(loaded.active_variable_indices)
        state_plus, state_minus = copy(jacobian_state), copy(jacobian_state)
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
    @test norm(analytical - numerical, Inf) < 2.0e-7

    result = run_spatial_model(path; duration = 0.8, samples = 81)
    @test maximum(maximum(abs(perp_position(item, state))
        for item in primitives) for state in result.states) < 2.0e-8
    @test maximum(maximum(abs(perp_velocity(item, state))
        for item in primitives) for state in result.states) < 2.0e-8
    @test maximum(norm(state[body.euler_parameter_variables] -
        initial[body.euler_parameter_variables])
        for state in result.states) < 2.0e-8
    @test norm(last(result.states)[free_body.euler_parameter_variables] -
        initial[free_body.euler_parameter_variables]) > 0.2
    terminal = last(result.states)
    @test terminal[body.position_variables] ≈
        [0.8, -0.18, 0.8 - 0.5 * 9.81 * 0.8^2] atol = 2.0e-6

    ground_to_ground = replace(source,
        "markers = [\"locked_body.frame\", \"ground.reference\"]" =>
            "markers = [\"ground.world\", \"ground.reference\"]")
    @test_throws ArgumentError load_spatial_model(
        IOBuffer(ground_to_ground))
end

@testset "Spatial fixed joint" begin
    path = joinpath(SPATIAL_MODEL_DIRECTORY,
        "fixed-two-body-assembly.toml")
    source = read(path, String)
    loaded = load_spatial_model(path)
    joint = loaded.connections[:weld]
    initial = loaded.initial_values

    @test joint isa SpatialFixedJoint
    @test joint.spherical.marker_a === joint.orient.marker_i
    @test joint.spherical.marker_b === joint.orient.marker_j
    @test length(loaded.layout.catalog.variables) == 50
    @test length(loaded.layout.catalog.equations) == 62
    @test length(loaded.active_variable_indices) == 50
    @test length(loaded.active_equation_indices) == 50
    @test loaded.analysis.degrees_of_freedom == 6
    @test norm(spatial_marker_position(joint.marker_a, initial) -
        spatial_marker_position(joint.marker_b, initial)) < 1.0e-14
    @test norm(spatial_marker_orientation(joint.marker_a, initial) -
        spatial_marker_orientation(joint.marker_b, initial)) < 1.0e-14

    primitives = (joint.orient.hinge.perp_xz,
        joint.orient.hinge.perp_yz, joint.orient.perp_xy)
    @test maximum(abs(perp_position(item, initial))
        for item in primitives) < 1.0e-14
    @test maximum(abs(perp_velocity(item, initial))
        for item in primitives) < 1.0e-14

    derivative = SpatialSimulationRunner.initial_spatial_derivative(
        initial, loaded)
    selection = AnalysisSelection(Dynamics(), loaded.active_variable_indices,
        loaded.active_equation_indices)
    jacobian_state = copy(initial)
    jacobian_state[joint.spherical.reaction_variables] .= [0.8, -0.5, 0.3]
    for (primitive, reaction) in zip(primitives, (1.1, -0.7, 0.4))
        jacobian_state[primitive.reaction_variable] = reaction
    end
    equations = zeros(length(selection.equation_indices))
    evaluate_analysis_equations!(equations, loaded.model, selection, 0.0,
        jacobian_state, derivative)
    coefficient = 2.0
    analytical = Matrix(evaluate_analysis_sparse_jacobian(loaded.model,
        selection, 0.0, jacobian_state, derivative, coefficient))
    numerical = similar(analytical)
    step = 1.0e-7
    plus, minus = similar(equations), similar(equations)
    for (column, variable) in enumerate(loaded.active_variable_indices)
        state_plus, state_minus = copy(jacobian_state), copy(jacobian_state)
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

    result = run_spatial_model(path; duration = 1.2, samples = 121)
    @test maximum(norm(spatial_marker_position(joint.marker_a, state) -
        spatial_marker_position(joint.marker_b, state))
        for state in result.states) < 2.0e-8
    @test maximum(norm(spatial_marker_orientation(joint.marker_a, state) -
        spatial_marker_orientation(joint.marker_b, state))
        for state in result.states) < 2.0e-8
    @test maximum(maximum(abs(perp_position(item, state))
        for item in primitives) for state in result.states) < 2.0e-8
    @test norm(spatial_marker_orientation(joint.marker_a,
        last(result.states)) -
        spatial_marker_orientation(joint.marker_a, initial)) > 0.2

    fixed_to_ground_source = """
        [model]
        dimension = "spatial"

        [ground]
        type = "ground"

        [ground.pin]
        type = "marker"

        [body]
        type = "rigid_body"
        mass = 1.0
        inertia = [0.01, 0.08, 0.08]
        position = [0.5, 0.0, 0.0]

        [body.pin]
        type = "marker"
        position = [-0.5, 0.0, 0.0]

        [weld]
        type = "fixed"
        markers = ["body.pin", "ground.pin"]
        """
    fixed_to_ground = load_spatial_model(IOBuffer(fixed_to_ground_source))
    @test fixed_to_ground.analysis.degrees_of_freedom == 0
    fixed_to_ground_result = run_spatial_model(
        IOBuffer(fixed_to_ground_source); duration = 0.02, samples = 3)
    @test maximum(norm(state[fixed_to_ground.bodies[:body].position_variables] -
        [0.5, 0.0, 0.0]) for state in fixed_to_ground_result.states) <
        1.0e-12

    ground_to_ground = replace(source,
        "markers = [\"base.joint\", \"arm.joint\"]" =>
            "markers = [\"ground.world\", \"ground.world\"]")
    @test_throws ArgumentError load_spatial_model(
        IOBuffer(ground_to_ground))
end
