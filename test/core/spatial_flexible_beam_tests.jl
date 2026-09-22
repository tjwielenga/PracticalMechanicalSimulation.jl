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
