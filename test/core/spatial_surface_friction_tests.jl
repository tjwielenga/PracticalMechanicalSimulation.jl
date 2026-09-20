using Test
using PracticalMechanicalSimulation
using PracticalMechanicalSimulation.SpatialFrictionForces

@testset "Spatial surface friction and sliding block" begin
    path = normpath(joinpath(@__DIR__, "..", "..", "models", "spatial",
        "sliding-block-friction.toml"))
    loaded = load_spatial_model(path)
    friction = loaded.forces[:friction]
    @test friction isa SpatialSurfaceFriction
    @test friction.contact === loaded.forces[:support]
    @test length(friction.shear_variables) == 2

    result = run_spatial_model(path)
    body = result.loaded.bodies[:block]
    friction = result.loaded.forces[:friction]
    normal = result.loaded.forces[:support]
    position = body.position_variables[1]
    velocity = body.velocity_variables[1]
    state_at(time) = result.states[argmin(abs.(result.times .- time))]
    @test abs(state_at(0.8)[velocity]) < 1.0e-3
    @test abs(state_at(0.8)[position] - state_at(1.0)[position]) < 1.0e-3
    @test state_at(1.8)[velocity] > 1.0
    @test all(result.states) do state
        hypot(state[friction.force_variables]...) <=
            friction.static_coefficient *
                max(state[normal.normal_force_variable], 0) + 1.0e-7
    end

    # The same element supplies a static load through finite shear, without
    # requiring a small nonzero sliding velocity or a temporary guide joint.
    source = read(path, String)
    source = replace(source, "mode = \"dynamic\"" => "mode = \"static\"")
    source = replace(source, "velocity = [1.0, 0.0, 0.0]" =>
        "velocity = [0.0, 0.0, 0.0]")
    source = replace(source,
        "expression = \"12*max(0,min(1,(t-1.0)/0.4))\"" =>
            "force = 2.0")
    static = run_spatial_model(IOBuffer(source); samples = 2)
    static_friction = static.loaded.forces[:friction]
    @test static.analysis_mode == :static
    @test last(static.states)[static_friction.force_variables[1]] ≈
        -2.0 atol = 1.0e-5
    @test last(static.states)[static_friction.shear_variables[1]] ≈
        2.0 / static_friction.stiffness atol = 1.0e-6

    handoff_source = replace(source, "mode = \"static\"" =>
        "mode = \"dynamic\"\ninitialization = \"static_equilibrium\"")
    handoff = run_spatial_model(IOBuffer(handoff_source);
        duration = 0.1, samples = 5)
    handoff_body = handoff.loaded.bodies[:block]
    handoff_friction = handoff.loaded.forces[:friction]
    @test first(handoff.states)[handoff_friction.shear_variables[1]] ≈
        2.0 / handoff_friction.stiffness atol = 1.0e-5
    @test abs(last(handoff.states)[handoff_body.velocity_variables[1]]) <
        1.0e-3

    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(
        read(path, String), "contact = \"support\"" =>
            "contact = \"missing\"")))
end
