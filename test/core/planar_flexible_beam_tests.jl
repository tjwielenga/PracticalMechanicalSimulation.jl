using Test
using TOML
using LinearAlgebra
using PracticalMechanicalSimulation
using PracticalMechanicalSimulation.PlanarAppliedForces
using PracticalMechanicalSimulation.PlanarComponentAssembly

const PLANAR_FLEXIBLE_CANTILEVER = normpath(joinpath(@__DIR__, "..", "..",
    "models", "planar", "flexible-cantilever.toml"))
const PLANAR_RECTANGULAR_CANTILEVER = normpath(joinpath(@__DIR__, "..", "..",
    "models", "planar", "rectangular-flexible-cantilever.toml"))

@testset "Planar floating-reference flexible beam" begin
    loaded = load_planar_model(PLANAR_FLEXIBLE_CANTILEVER)
    beam = loaded.bodies[:beam]
    @test beam isa PlanarFlexibleBeamComponent
    @test haskey(loaded.markers, Symbol("beam.cm"))
    @test loaded.analysis.degrees_of_freedom == 3
    @test Set(loaded.state_selection.selected_velocities) == Set((
        Symbol("beam.eta_u_dot"), Symbol("beam.eta_v_dot"),
        Symbol("beam.eta_beta_dot")))
    @test all(isposdef, (beam.elastic_mass, beam.elastic_stiffness))

    static_result = run_planar_model(loaded)
    static_state = last(static_result.states)
    @test all(iszero, static_state[beam.elastic_velocity_variables])
    @test all(iszero, static_state[beam.elastic_acceleration_variables])
    tip = PlanarAppliedForces.point_marker_kinematics(
        loaded.markers[Symbol("beam.end_j")].point,
        static_state).position
    elementary_tip_deflection =
        loaded.bodies[:beam].mass * 9.81 * beam.length^3 /
        (8 * 2.0e7 * 8.333333333333333e-6)
    @test tip[2] < 0
    @test -tip[2] ≈ elementary_tip_deflection rtol = 0.02

    document = TOML.parsefile(PLANAR_FLEXIBLE_CANTILEVER)
    document["analysis"]["mode"] = "modal"
    document["analysis"]["modes"] = 3
    modal_model = load_planar_model(document;
        source_directory = dirname(PLANAR_FLEXIBLE_CANTILEVER))
    modes = run_planar_model(modal_model)
    @test modes.natural_frequencies_hz[1] ≈ 7.22 rtol = 0.02
    @test maximum(modes.equation_errors) < 1.0e-10
end

@testset "Planar flexible-beam section conveniences" begin
    rectangular_model = load_planar_model(PLANAR_RECTANGULAR_CANTILEVER)
    rectangular_beam = rectangular_model.bodies[:beam]
    expanded_rectangle = TOML.parse(rectangular_model.model_source)["beam"]
    @test rectangular_beam.mass ≈ 100.0 * 0.10 * 0.10
    @test expanded_rectangle["area"] ≈ 0.10 * 0.10
    @test expanded_rectangle["second_moment"] ≈ 0.10 * 0.10^3 / 12
    @test expanded_rectangle["shear_modulus"] ≈ 8.0e6
    @test expanded_rectangle["shear_coefficient"] ≈ 5 / 6
    @test expanded_rectangle["graphics"]["member"]["shape"] == "box"
    @test expanded_rectangle["graphics"]["member"]["width"] == 0.10
    @test expanded_rectangle["graphics"]["member"]["height"] == 0.10
    rectangular_result = run_planar_model(rectangular_model)
    @test last(rectangular_result.states)[
        rectangular_beam.elastic_position_variables[2]] < 0

    circular_document = TOML.parsefile(PLANAR_FLEXIBLE_CANTILEVER)
    circular = circular_document["beam"]
    for field in ("mass", "inertia", "area", "second_moment",
            "shear_modulus", "shear_coefficient", "graphics")
        pop!(circular, field, nothing)
    end
    circular["density"] = 1000.0
    circular["poisson_ratio"] = 0.25
    circular["section"] = Dict(
        "shape" => "circular", "diameter" => 0.10)
    circular_model = load_planar_model(circular_document;
        source_directory = dirname(PLANAR_FLEXIBLE_CANTILEVER))
    circular_beam = circular_model.bodies[:beam]
    expanded_circular = TOML.parse(circular_model.model_source)["beam"]
    @test circular_beam.mass ≈ 1000 * pi * 0.05^2
    @test expanded_circular["area"] ≈ pi * 0.05^2
    @test expanded_circular["second_moment"] ≈ pi * 0.05^4 / 4
    @test expanded_circular["shear_modulus"] ≈ 8.0e6
    @test expanded_circular["shear_coefficient"] ≈ 6 / 7
    @test expanded_circular["graphics"]["member"]["radius"] == 0.05
    circular_round_trip = load_planar_model(IOBuffer(
        circular_model.model_source))
    @test circular_round_trip.model_source == circular_model.model_source

    conflicting = deepcopy(circular_document)
    conflicting["beam"]["area"] = 1.0
    @test_throws ArgumentError load_planar_model(conflicting;
        source_directory = dirname(PLANAR_FLEXIBLE_CANTILEVER))

    api_model = PracticalMechanicalSimulation.Sim2D.Model(:section_api)
    PracticalMechanicalSimulation.Sim2D.flexible_beam!(api_model, :beam;
        length = 0.5,
        density = 2700.0,
        elastic_modulus = 7.0e10,
        poisson_ratio = 0.33,
        section = (shape = :circular, radius = 0.01))
    api_beam = load_planar_model(api_model).bodies[:beam]
    @test api_beam.mass ≈ 2700 * pi * 0.01^2 * 0.5
end
