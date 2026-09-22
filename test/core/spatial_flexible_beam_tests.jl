using Test
using TOML
using LinearAlgebra
using PracticalMechanicalSimulation
using PracticalMechanicalSimulation.SpatialComponentAssembly
using PracticalMechanicalSimulation.SpatialModeling

const SPATIAL_FLEXIBLE_CANTILEVER = normpath(joinpath(@__DIR__, "..", "..",
    "models", "spatial", "flexible-cantilever.toml"))

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
    tip = spatial_marker_position(loaded.markers[Symbol("beam.end_j")],
        last(static_result.states))
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
    @test modes.natural_frequencies_hz[1:2] ≈ [7.2176, 7.2176] rtol = 0.002
    @test maximum(modes.equation_errors) < 1.0e-10
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
