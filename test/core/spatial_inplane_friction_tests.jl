using Test
using LinearAlgebra
using PracticalMechanicalSimulation
using PracticalMechanicalSimulation.SpatialFrictionForces

@testset "Spatial inplane constraint friction" begin
    path = normpath(joinpath(@__DIR__, "..", "..", "models", "spatial",
        "inplane-friction-block.toml"))
    loaded = load_spatial_model(path)
    friction = loaded.forces[:plane_friction]
    @test friction isa SpatialInplaneFriction
    @test friction.constraint === loaded.connections[:support]
    @test length(friction.shear_variables) == 2

    state = copy(loaded.initial_values)
    reaction = friction.constraint.reaction_variable
    state[reaction] = -4.0
    @test SpatialFrictionForces.inplane_normal_load(friction, state) ≈ 4.0
    state[reaction] = 0.0
    state[friction.normal_load_variable] = 0.0
    state[friction.shear_variables] .= [0.001, -0.002]
    @test SpatialFrictionForces.calculated_inplane_friction_force(
        friction, state) == [0.0, 0.0]
    @test inplane_friction_rates(friction, state) ≈
        -state[friction.shear_variables] ./ friction.release_time

    result = run_spatial_model(path)
    body = result.loaded.bodies[:block]
    friction = result.loaded.forces[:plane_friction]
    state_at(time) = result.states[argmin(abs.(result.times .- time))]
    @test maximum(abs, state_at(0.8)[body.velocity_variables[1:2]]) <
        1.0e-3
    @test norm(state_at(0.8)[body.position_variables[1:2]] -
        state_at(1.0)[body.position_variables[1:2]]) < 1.0e-3
    @test state_at(1.8)[body.velocity_variables[1]] > 1.0
    @test abs(state_at(1.8)[body.velocity_variables[2]]) < 1.0e-3
    @test all(result.states) do sample
        hypot(sample[friction.force_variables]...) <=
            friction.static_coefficient *
                sample[friction.normal_load_variable] + 1.0e-7
    end

    source = read(path, String)
    source = replace(source, "mode = \"dynamic\"" => "mode = \"static\"")
    source = replace(source, "velocity = [1.0, 0.4, 0.0]" =>
        "velocity = [0.0, 0.0, 0.0]")
    source = replace(source,
        "expression = \"12*max(0,min(1,(t-1.0)/0.4))\"" =>
            "force = 2.0")
    source *= """

        [ground.y_push]
        type = "marker"
        orientation = ["-90 deg", 1, 0, 0]

        [side_push]
        type = "applied_force"
        markers = ["block.frame", "ground.y_push"]
        force = 3.0
        """
    static = run_spatial_model(IOBuffer(source); samples = 2)
    static_friction = static.loaded.forces[:plane_friction]
    @test static.analysis_mode == :static
    @test last(static.states)[static_friction.force_variables] ≈
        [-2.0, -3.0] atol = 1.0e-6
    @test last(static.states)[static_friction.shear_variables] ≈
        [2.0, 3.0] ./ static_friction.stiffness atol = 1.0e-6

    handoff_source = replace(source, "mode = \"static\"" =>
        "mode = \"dynamic\"\ninitialization = \"static_equilibrium\"")
    handoff = run_spatial_model(IOBuffer(handoff_source);
        duration = 0.1, samples = 5)
    handoff_friction = handoff.loaded.forces[:plane_friction]
    @test first(handoff.states)[handoff_friction.shear_variables] ≈
        [2.0, 3.0] ./ handoff_friction.stiffness atol = 1.0e-5
    @test maximum(abs, last(handoff.states)[
        handoff.loaded.bodies[:block].velocity_variables[1:2]]) < 1.0e-3

    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(
        read(path, String), "constraint = \"support\"" =>
            "constraint = \"orientation\"")))
end
