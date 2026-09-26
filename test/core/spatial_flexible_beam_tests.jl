using Test
using TOML
using LinearAlgebra
using PracticalMechanicalSimulation
using PracticalMechanicalSimulation.SpatialComponentAssembly
using PracticalMechanicalSimulation.SpatialModeling

const SPATIAL_FLEXIBLE_CANTILEVER = normpath(joinpath(@__DIR__, "..", "..",
    "models", "spatial", "flexible-cantilever.toml"))
const SPATIAL_RECTANGULAR_CANTILEVER = normpath(joinpath(@__DIR__, "..", "..",
    "models", "spatial", "rectangular-flexible-cantilever.toml"))

function tip_loaded_cantilever(path, force, orientation)
    document = TOML.parsefile(path)
    pop!(document, "gravity", nothing)
    document["analysis"]["mode"] = "static"
    document["analysis"]["static_method"] = "newton"
    document["analysis"]["static_tolerance"] = 1.0e-10
    document["ground"]["load_axis"] = Dict(
        "type" => "marker", "orientation" => orientation)
    document["tip_load"] = Dict(
        "type" => "applied_force",
        "markers" => ["beam.end_j", "ground.load_axis"],
        "force" => force)
    loaded = load_spatial_model(document; source_directory = dirname(path))
    result = run_spatial_model(loaded)
    tip = spatial_marker_position(loaded.markers[Symbol("beam.end_j")],
        last(result.states))
    (; loaded, result, tip)
end

function fixed_end_tip_compliance(beam, tip_load)
    L = beam.length
    rigid_at(offset) = [
        Matrix{Float64}(I, 3, 3) -skew(offset)
        zeros(3, 3) Matrix{Float64}(I, 3, 3)
    ]
    rigid_i = rigid_at([-L / 2, 0.0, 0.0])
    rigid_j = rigid_at([L / 2, 0.0, 0.0])
    shape_i = @view beam.deformation_shape[1:6, :]
    shape_j = @view beam.deformation_shape[7:12, :]
    reference_from_elastic = -(rigid_i \ shape_i)
    tip_from_elastic = rigid_j * reference_from_elastic + shape_j
    elastic = beam.elastic_stiffness \
        (transpose(tip_from_elastic) * tip_load)
    tip_from_elastic * elastic
end

function beam_mechanical_energy(beam, state)
    velocity = @view state[beam.velocity_variables]
    omega = @view state[beam.angular_velocity_variables]
    elastic_velocity = @view state[beam.elastic_velocity_variables]
    elastic_position = @view state[beam.elastic_position_variables]
    0.5beam.mass * dot(velocity, velocity) +
        0.5dot(omega, beam.inertia * omega) +
        0.5dot(elastic_velocity,
            beam.elastic_mass * elastic_velocity) +
        0.5dot(elastic_position,
            beam.elastic_stiffness * elastic_position)
end

@testset "Spatial floating-reference flexible beam" begin
    loaded = load_spatial_model(SPATIAL_FLEXIBLE_CANTILEVER)
    beam = loaded.bodies[:beam]
    @test beam isa SpatialFlexibleBeamComponent
    @test all(name -> haskey(loaded.markers, Symbol("beam.", name)),
        (:end_i, :cm, :end_j))
    @test loaded.analysis.degrees_of_freedom == 6
    @test Set(loaded.state_selection.selected_velocities) == Set(
        Symbol("beam.eta_$(coordinate)_dot")
        for coordinate in (:u, :v, :w, :rx, :ry, :rz))
    @test all(isposdef, (beam.elastic_mass, beam.elastic_stiffness))

    static_result = run_spatial_model(loaded)
    static_state = last(static_result.states)
    @test all(iszero, static_state[beam.elastic_velocity_variables])
    @test all(iszero, static_state[beam.elastic_acceleration_variables])
    tip = spatial_marker_position(loaded.markers[Symbol("beam.end_j")],
        static_state)
    elementary_tip_deflection = beam.mass * 9.81 * beam.length^3 /
        (8 * 2.0e7 * 8.333333333333333e-6)
    @test tip[3] < 0
    @test -tip[3] ≈ elementary_tip_deflection rtol = 0.05

    document = TOML.parsefile(SPATIAL_FLEXIBLE_CANTILEVER)
    document["analysis"]["mode"] = "dynamic"
    document["simulation"]["end_time"] = 0.2
    document["simulation"]["output_samples"] = 21
    dynamic_model = load_spatial_model(document;
        source_directory = dirname(SPATIAL_FLEXIBLE_CANTILEVER))
    motion = run_spatial_model(dynamic_model)
    tip_z = [spatial_marker_position(
        dynamic_model.markers[Symbol("beam.end_j")], state)[3]
        for state in motion.states]
    @test minimum(tip_z) < -0.01
    @test maximum(tip_z) ≈ 0.0 atol = 1.0e-12

    document = TOML.parsefile(SPATIAL_FLEXIBLE_CANTILEVER)
    document["analysis"]["mode"] = "modal"
    document["analysis"]["modes"] = 6
    modal_model = load_spatial_model(document;
        source_directory = dirname(SPATIAL_FLEXIBLE_CANTILEVER))
    modes = run_spatial_model(modal_model)
    first_euler_bernoulli_frequency = 1.875104068711961^2 / (2pi) *
        sqrt(2.0e7 * 8.333333333333333e-6 /
             (beam.mass * beam.length^3))
    @test modes.natural_frequencies_hz[1:2] ≈
        fill(first_euler_bernoulli_frequency, 2) rtol = 0.002
    @test maximum(modes.equation_errors) < 1.0e-10
end

@testset "Spatial flexible-beam analytical verification" begin
    axial_force = 10.0
    axial = tip_loaded_cantilever(SPATIAL_RECTANGULAR_CANTILEVER,
        axial_force, [pi / 2, 0.0, 1.0, 0.0])
    beam = axial.loaded.bodies[:beam]
    section = TOML.parse(axial.loaded.model_source)["beam"]
    expected_extension = axial_force * beam.length /
        (section["elastic_modulus"] * section["area"])
    @test axial.tip[1] - beam.length ≈ expected_extension rtol = 2.0e-4
    @test abs(axial.tip[2]) < 1.0e-10
    @test abs(axial.tip[3]) < 1.0e-10

    transverse_force = -0.1
    about_y = tip_loaded_cantilever(SPATIAL_RECTANGULAR_CANTILEVER,
        transverse_force, [0.0, 0.0, 0.0, 1.0])
    about_z = tip_loaded_cantilever(SPATIAL_RECTANGULAR_CANTILEVER,
        transverse_force, [-pi / 2, 1.0, 0.0, 0.0])
    E = section["elastic_modulus"]
    G = section["shear_modulus"]
    A = section["area"]
    expected_z = abs(transverse_force) * (
        beam.length^3 / (3E * section["second_moment_y"]) +
        beam.length / (section["shear_coefficient_z"] * G * A))
    expected_y = abs(transverse_force) * (
        beam.length^3 / (3E * section["second_moment_z"]) +
        beam.length / (section["shear_coefficient_y"] * G * A))
    @test -about_y.tip[3] ≈ expected_z rtol = 1.0e-3
    @test -about_z.tip[2] ≈ expected_y rtol = 1.0e-3
    @test abs(about_z.tip[2] / about_y.tip[3]) > 3.9

    torque = 0.2
    tip_motion = fixed_end_tip_compliance(beam,
        [0.0, 0.0, 0.0, torque, 0.0, 0.0])
    expected_twist = torque * beam.length /
        (G * section["torsion_constant"])
    @test tip_motion[4] ≈ expected_twist rtol = 1.0e-12
    @test norm(tip_motion[[1, 2, 3, 5, 6]]) < 1.0e-12
end

@testset "Spatial flexible-beam dynamic verification" begin
    damped_document = TOML.parsefile(SPATIAL_FLEXIBLE_CANTILEVER)
    pop!(damped_document, "gravity", nothing)
    damped_document["analysis"]["mode"] = "dynamic"
    damped_document["simulation"]["end_time"] = 0.5
    damped_document["simulation"]["output_samples"] = 101
    damped_document["simulation"]["maximum_step"] = 0.002
    damped_document["beam"]["damping_time_scale"] = 0.01
    damped_document["beam"]["elastic_position"] =
        [0.0, 0.0, 0.01, 0.0, 0.015, 0.0]
    damped = load_spatial_model(damped_document;
        source_directory = dirname(SPATIAL_FLEXIBLE_CANTILEVER))
    damped_result = run_spatial_model(damped)
    damped_beam = damped.bodies[:beam]
    energy = [beam_mechanical_energy(damped_beam, state)
        for state in damped_result.states]
    @test energy[end] < 0.4energy[1]
    @test maximum(diff(energy)) < 2.0e-5energy[1]

    free_document = TOML.parsefile(SPATIAL_FLEXIBLE_CANTILEVER)
    pop!(free_document, "gravity", nothing)
    pop!(free_document, "support", nothing)
    free_document["analysis"]["mode"] = "dynamic"
    free_document["simulation"]["end_time"] = 1.0
    free_document["simulation"]["output_samples"] = 51
    free_document["beam"]["velocity"] = [0.3, -0.2, 0.1]
    free_document["beam"]["angular_velocity"] = [2.0, 0.0, 0.0]
    free = load_spatial_model(free_document;
        source_directory = dirname(SPATIAL_FLEXIBLE_CANTILEVER))
    free_result = run_spatial_model(free)
    free_beam = free.bodies[:beam]
    maximum_elastic_motion = maximum(norm(
        @view state[free_beam.elastic_position_variables])
        for state in free_result.states)
    @test maximum_elastic_motion < 1.0e-9
    expected_center = [0.5, 0.0, 0.0] + [0.3, -0.2, 0.1]
    @test free_result.states[end][free_beam.position_variables] ≈
        expected_center atol = 2.0e-7
end

@testset "Spatial flexible-beam section conveniences" begin
    circular_document = TOML.parsefile(SPATIAL_FLEXIBLE_CANTILEVER)
    circular = circular_document["beam"]
    for field in ("mass", "area", "shear_modulus", "second_moment_y",
            "second_moment_z", "torsion_constant", "shear_coefficient_y",
            "shear_coefficient_z", "graphics")
        pop!(circular, field, nothing)
    end
    circular["density"] = 1000.0
    circular["poisson_ratio"] = 0.25
    circular["section"] = Dict("shape" => "circular", "diameter" => 0.1)
    circular_model = load_spatial_model(circular_document;
        source_directory = dirname(SPATIAL_FLEXIBLE_CANTILEVER))
    circular_beam = circular_model.bodies[:beam]
    @test circular_beam.mass ≈ 1000 * pi * 0.05^2
    @test isposdef(circular_beam.elastic_stiffness)
    expanded_circular = TOML.parse(circular_model.model_source)["beam"]
    @test expanded_circular["area"] ≈ pi * 0.05^2
    @test expanded_circular["second_moment_y"] ≈ pi * 0.05^4 / 4
    @test expanded_circular["second_moment_z"] ≈ pi * 0.05^4 / 4
    @test expanded_circular["torsion_constant"] ≈ pi * 0.05^4 / 2
    @test expanded_circular["shear_modulus"] ≈ 8.0e6
    @test expanded_circular["graphics"]["member"]["radius"] == 0.05
    circular_round_trip = load_spatial_model(
        IOBuffer(circular_model.model_source))
    @test circular_round_trip.model_source == circular_model.model_source

    rectangular_document = TOML.parsefile(SPATIAL_FLEXIBLE_CANTILEVER)
    rectangular = rectangular_document["beam"]
    for field in ("area", "second_moment_y", "second_moment_z",
            "torsion_constant", "shear_coefficient_y",
            "shear_coefficient_z", "graphics")
        pop!(rectangular, field, nothing)
    end
    rectangular["section"] = Dict(
        "shape" => "rectangular", "width" => 0.08, "height" => 0.12)
    rectangular_model = load_spatial_model(rectangular_document;
        source_directory = dirname(SPATIAL_FLEXIBLE_CANTILEVER))
    expanded_rectangle = TOML.parse(rectangular_model.model_source)["beam"]
    @test expanded_rectangle["area"] ≈ 0.08 * 0.12
    @test expanded_rectangle["second_moment_y"] ≈ 0.08 * 0.12^3 / 12
    @test expanded_rectangle["second_moment_z"] ≈ 0.12 * 0.08^3 / 12
    long_side, short_side = 0.12, 0.08
    ratio = short_side / long_side
    expected_torsion = long_side * short_side^3 *
        (1 / 3 - 0.21 * ratio * (1 - ratio^4 / 12))
    @test expanded_rectangle["torsion_constant"] ≈ expected_torsion
    @test expanded_rectangle["graphics"]["member"]["shape"] == "box"
    @test expanded_rectangle["graphics"]["member"]["width"] == 0.08
    @test expanded_rectangle["graphics"]["member"]["height"] == 0.12

    conflicting = deepcopy(circular_document)
    conflicting["beam"]["area"] = 1.0
    @test_throws ArgumentError load_spatial_model(conflicting;
        source_directory = dirname(SPATIAL_FLEXIBLE_CANTILEVER))

    api_model = PracticalMechanicalSimulation.Sim3D.Model(:section_api)
    PracticalMechanicalSimulation.Sim3D.flexible_beam!(api_model, :beam;
        length = 0.5,
        density = 2700.0,
        elastic_modulus = 7.0e10,
        poisson_ratio = 0.33,
        section = (shape = :circular, radius = 0.01))
    api_beam = load_spatial_model(api_model).bodies[:beam]
    @test api_beam.mass ≈ 2700 * pi * 0.01^2 * 0.5
end

@testset "Connected flexible-beam motion" begin
    model = PracticalMechanicalSimulation.Sim3D.Model(
        :connected_flexible_beams)
    PracticalMechanicalSimulation.Sim3D.analysis!(model; mode = :dynamic)
    PracticalMechanicalSimulation.Sim3D.simulation!(model;
        end_time = 1.0e-3, output_samples = 2,
        initial_step = 1.0e-7, maximum_step = 1.0e-3)
    ground = PracticalMechanicalSimulation.Sim3D.ground!(model, :ground)
    hub = PracticalMechanicalSimulation.Sim3D.marker!(ground, :hub)
    for (name, center) in ((:first, 0.25), (:second, 0.75))
        PracticalMechanicalSimulation.Sim3D.flexible_beam!(model, name;
            length = 0.5, density = 1000.0,
            elastic_modulus = 2.0e8, poisson_ratio = 0.3,
            damping_time_scale = 0.01,
            section = (shape = :rectangular,
                width = 0.04, height = 0.02),
            position = [center, 0.0, 0.0])
    end
    PracticalMechanicalSimulation.Sim3D.fixed!(model, :connection;
        markers = ["second.end_i", "first.end_j"])
    root = PracticalMechanicalSimulation.Sim3D.revolute!(model, :root;
        markers = ["first.end_i", hub], rotation_coordinates = true)
    PracticalMechanicalSimulation.Sim3D.rotational_motion!(model, :drive;
        joint = root, angle = "0.1*t")

    loaded = load_spatial_model(model)
    result = run_spatial_model(loaded)
    final = last(result.states)
    @test last(result.times) == 1.0e-3
    @test all(isfinite, final)
    @test spatial_marker_position(loaded.markers[Symbol("first.end_j")],
        final) ≈ spatial_marker_position(
            loaded.markers[Symbol("second.end_i")], final) atol = 1.0e-9
end
