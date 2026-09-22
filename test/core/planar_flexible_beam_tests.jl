using Test
using TOML
using LinearAlgebra
using PracticalMechanicalSimulation
using PracticalMechanicalSimulation.PlanarAppliedForces
using PracticalMechanicalSimulation.PlanarComponentAssembly

const PLANAR_FLEXIBLE_CANTILEVER = normpath(joinpath(@__DIR__, "..", "..",
    "models", "planar", "flexible-cantilever.toml"))

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
    tip = PlanarAppliedForces.point_marker_kinematics(
        loaded.markers[Symbol("beam.end_j")].point,
        last(static_result.states)).position
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
