using Test
using PracticalMechanicalSimulation
using PracticalMechanicalSimulation.SpatialFrictionForces

@testset "Spatial revolute bearing friction" begin
    path = normpath(joinpath(@__DIR__, "..", "..", "models", "spatial",
        "revolute-bearing-friction.toml"))
    loaded = load_spatial_model(path)
    friction = loaded.forces[:bearing_friction]
    @test friction isa SpatialRevoluteFriction
    @test friction.joint === loaded.connections[:pin]
    @test friction.preload == 0

    # Axial thrust does not contribute to the radial-bearing friction limit.
    state = copy(loaded.initial_values)
    state[friction.joint.spherical.reaction_variables] .= [0.0, 0.0, 100.0]
    @test SpatialFrictionForces.revolute_bearing_load(friction, state) < 1.0e-8
    state[friction.joint.spherical.reaction_variables] .= [3.0, 4.0, 100.0]
    @test SpatialFrictionForces.revolute_bearing_load(friction, state) ≈ 5.0

    result = run_spatial_model(path)
    friction = result.loaded.forces[:bearing_friction]
    hinge = result.loaded.connections[:pin].hinge
    omega = hinge.rotation_variables[2]
    state_at(time) = result.states[argmin(abs.(result.times .- time))]
    @test abs(state_at(1.0)[omega]) < 1.0e-2
    @test state_at(2.0)[omega] > 1.0
    @test all(result.states) do sample
        abs(sample[friction.torque_variable]) <=
            friction.static_coefficient * friction.effective_radius *
                sample[friction.bearing_load_variable] + 1.0e-7
    end

    source = read(path, String)
    source = replace(source, "mode = \"dynamic\"" => "mode = \"static\"")
    source = replace(source, "angular_velocity = [0.0, 0.0, 6.0]" =>
        "angular_velocity = [0.0, 0.0, 0.0]")
    source = replace(source,
        "expression = \"2*max(0,min(1,(t-1.0)/0.4))\"" =>
            "torque = 0.3")
    static = run_spatial_model(IOBuffer(source); samples = 2)
    static_friction = static.loaded.forces[:bearing_friction]
    @test static.analysis_mode == :static
    @test last(static.states)[static_friction.torque_variable] ≈
        -0.3 atol = 1.0e-6
    @test last(static.states)[static_friction.shear_variable] ≈
        0.3 / static_friction.stiffness atol = 1.0e-6

    handoff_source = replace(source, "mode = \"static\"" =>
        "mode = \"dynamic\"\ninitialization = \"static_equilibrium\"")
    handoff = run_spatial_model(IOBuffer(handoff_source);
        duration = 0.1, samples = 5)
    handoff_friction = handoff.loaded.forces[:bearing_friction]
    handoff_omega = handoff.loaded.connections[:pin].hinge.rotation_variables[2]
    @test first(handoff.states)[handoff_friction.shear_variable] ≈
        0.3 / handoff_friction.stiffness atol = 1.0e-5
    @test abs(last(handoff.states)[handoff_omega]) < 1.0e-3

    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(
        read(path, String), "joint = \"pin\"" =>
            "joint = \"missing\"")))
end
