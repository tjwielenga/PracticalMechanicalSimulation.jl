using Test
using PracticalMechanicalSimulation
using PracticalMechanicalSimulation.SpatialFrictionForces

@testset "Spatial translational guide friction" begin
    path = normpath(joinpath(@__DIR__, "..", "..", "models", "spatial",
        "translational-guide-friction.toml"))
    loaded = load_spatial_model(path)
    friction = loaded.forces[:guide_friction]
    @test friction isa SpatialTranslationalFriction
    @test friction.joint === loaded.connections[:guide]
    @test friction.preload == 0

    state = copy(loaded.initial_values)
    state[friction.joint.inplane_x.reaction_variable] = 3.0
    state[friction.joint.inplane_y.reaction_variable] = 4.0
    @test SpatialFrictionForces.translational_guide_load(friction, state) ≈ 5.0
    state[friction.joint.inplane_x.reaction_variable] = 0.0
    state[friction.joint.inplane_y.reaction_variable] = 0.0
    state[friction.guide_load_variable] = 0.0
    state[friction.shear_variable] = 0.001
    @test SpatialFrictionForces.calculated_translational_friction_force(
        friction, state) == 0
    @test SpatialFrictionForces.translational_friction_rate(friction, state) ≈
        -0.001 / friction.release_time

    result = run_spatial_model(path)
    friction = result.loaded.forces[:guide_friction]
    joint = result.loaded.connections[:guide]
    velocity = joint.translation_variables[2]
    distance = joint.translation_variables[3]
    state_at(time) = result.states[argmin(abs.(result.times .- time))]
    @test abs(state_at(0.8)[velocity]) < 1.0e-3
    @test abs(state_at(0.8)[distance] - state_at(1.0)[distance]) < 1.0e-3
    @test state_at(1.75)[velocity] > 1.0
    @test all(result.states) do sample
        abs(sample[friction.force_variable]) <=
            friction.static_coefficient *
                sample[friction.guide_load_variable] + 1.0e-7
    end

    source = read(path, String)
    source = replace(source, "mode = \"dynamic\"" => "mode = \"static\"")
    source = replace(source, "velocity = [0.0, 0.0, 1.0]" =>
        "velocity = [0.0, 0.0, 0.0]")
    source = replace(source,
        "expression = \"12*max(0,min(1,(t-1.0)/0.4))\"" =>
            "force = 2.0")
    static = run_spatial_model(IOBuffer(source); samples = 2)
    static_friction = static.loaded.forces[:guide_friction]
    @test static.analysis_mode == :static
    @test last(static.states)[static_friction.force_variable] ≈
        -2.0 atol = 1.0e-6
    @test last(static.states)[static_friction.shear_variable] ≈
        2.0 / static_friction.stiffness atol = 1.0e-6

    handoff_source = replace(source, "mode = \"static\"" =>
        "mode = \"dynamic\"\ninitialization = \"static_equilibrium\"")
    handoff = run_spatial_model(IOBuffer(handoff_source);
        duration = 0.1, samples = 5)
    handoff_friction = handoff.loaded.forces[:guide_friction]
    handoff_velocity =
        handoff.loaded.connections[:guide].translation_variables[2]
    @test first(handoff.states)[handoff_friction.shear_variable] ≈
        2.0 / handoff_friction.stiffness atol = 1.0e-5
    @test abs(last(handoff.states)[handoff_velocity]) < 1.0e-3

    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(
        read(path, String), "joint = \"guide\"" =>
            "joint = \"orientation\"")))
end
