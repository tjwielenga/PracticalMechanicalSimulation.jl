using Test
using PracticalMechanicalSimulation
using PracticalMechanicalSimulation.PlanarFrictionForces

planar_friction_model(name) = normpath(joinpath(@__DIR__, "..", "..",
    "models", "planar", name))

@testset "Planar revolute bearing friction" begin
    path = planar_friction_model("revolute-bearing-friction.toml")
    loaded = load_planar_model(path)
    friction = only(loaded.forces[:bearing_friction])
    @test friction isa PlanarRevoluteFriction
    @test friction.joint === loaded.connections[:pin]
    state = copy(loaded.initial_values)
    state[friction.joint.reaction_variables] .= [3.0, 4.0]
    @test revolute_bearing_load(friction, state) ≈ 5.0

    result = run_planar_model(path)
    friction = only(result.loaded.forces[:bearing_friction])
    omega = result.loaded.connections[:pin].rotation_variables[2]
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
    source = replace(source, "angular_velocity = 6.0" =>
        "angular_velocity = 0.0")
    source = replace(source,
        "expression = \"2*max(0,min(1,(t-1.0)/0.4))\"" =>
            "torque = 0.3")
    static = run_planar_model(IOBuffer(source); end_time = 0.0, samples = 1)
    static_friction = only(static.loaded.forces[:bearing_friction])
    @test only(static.states)[static_friction.torque_variable] ≈
        -0.3 atol = 1.0e-6
    @test only(static.states)[static_friction.shear_variable] ≈
        0.3 / static_friction.stiffness atol = 1.0e-6

    handoff_source = replace(source, "mode = \"static\"" =>
        "mode = \"dynamic\"\ninitialization = \"static_equilibrium\"")
    handoff = run_planar_model(IOBuffer(handoff_source);
        duration = 0.1, samples = 5)
    handoff_friction = only(handoff.loaded.forces[:bearing_friction])
    handoff_omega = handoff.loaded.connections[:pin].rotation_variables[2]
    @test first(handoff.states)[handoff_friction.shear_variable] ≈
        0.3 / handoff_friction.stiffness atol = 1.0e-5
    @test abs(last(handoff.states)[handoff_omega]) < 1.0e-3

    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        read(path, String), "joint = \"pin\"" => "joint = \"missing\"")))
end

@testset "Planar translational guide friction" begin
    path = planar_friction_model("translational-guide-friction.toml")
    loaded = load_planar_model(path)
    friction = only(loaded.forces[:guide_friction])
    @test friction isa PlanarTranslationalFriction
    @test friction.joint === loaded.connections[:guide]
    state = copy(loaded.initial_values)
    state[friction.joint.inplane.reaction_variable] = -4.0
    @test translational_guide_load(friction, state) ≈ 4.0
    state[friction.joint.inplane.reaction_variable] = 0.0
    state[friction.load_variable] = 0.0
    state[friction.shear_variable] = 0.001
    @test PlanarFrictionForces.calculated_tangential_friction_force(
        friction, state) == 0
    @test translational_friction_rate(friction, state) ≈
        -0.001 / friction.release_time

    result = run_planar_model(path)
    friction = only(result.loaded.forces[:guide_friction])
    body = result.loaded.bodies[:slider]
    state_at(time) = result.states[argmin(abs.(result.times .- time))]
    @test abs(state_at(0.8)[body.velocity_variables[1]]) < 1.0e-3
    @test abs(state_at(1.0)[body.position_variables[1]] -
        state_at(0.8)[body.position_variables[1]]) < 1.0e-3
    @test state_at(1.75)[body.velocity_variables[1]] > 1.0
    @test all(result.states) do sample
        abs(sample[friction.force_variable]) <=
            friction.static_coefficient * sample[friction.load_variable] +
                1.0e-7
    end

    source = replace(read(path, String),
        "mode = \"dynamic\"" => "mode = \"static\"",
        "velocity = [1.0, 0.0]" => "velocity = [0.0, 0.0]",
        "expression = \"12*max(0,min(1,(t-1.0)/0.4))\"" => "force = 2.0")
    static = run_planar_model(IOBuffer(source); end_time = 0.0, samples = 1)
    static_friction = only(static.loaded.forces[:guide_friction])
    @test only(static.states)[static_friction.force_variable] ≈
        -2.0 atol = 1.0e-6
    @test only(static.states)[static_friction.shear_variable] ≈
        2.0 / static_friction.stiffness atol = 1.0e-6

    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        read(path, String), "joint = \"guide\"" => "joint = \"missing\"")))
end

@testset "Planar inplane friction" begin
    path = planar_friction_model("inplane-friction-block.toml")
    loaded = load_planar_model(path)
    friction = only(loaded.forces[:plane_friction])
    @test friction isa PlanarInplaneFriction
    @test friction.constraint === loaded.connections[:support]
    state = copy(loaded.initial_values)
    state[friction.constraint.reaction_variable] = -4.0
    @test inplane_normal_load(friction, state) ≈ 4.0

    result = run_planar_model(path)
    friction = only(result.loaded.forces[:plane_friction])
    body = result.loaded.bodies[:block]
    state_at(time) = result.states[argmin(abs.(result.times .- time))]
    @test abs(state_at(0.8)[body.velocity_variables[1]]) < 1.0e-3
    @test abs(state_at(1.0)[body.position_variables[1]] -
        state_at(0.8)[body.position_variables[1]]) < 1.0e-3
    @test state_at(1.8)[body.velocity_variables[1]] > 1.0
    @test all(result.states) do sample
        abs(sample[friction.force_variable]) <=
            friction.static_coefficient * sample[friction.load_variable] +
                1.0e-7
    end

    source = replace(read(path, String),
        "mode = \"dynamic\"" => "mode = \"static\"",
        "velocity = [1.0, 0.0]" => "velocity = [0.0, 0.0]",
        "expression = \"12*max(0,min(1,(t-1.0)/0.4))\"" => "force = 2.0")
    static = run_planar_model(IOBuffer(source); end_time = 0.0, samples = 1)
    static_friction = only(static.loaded.forces[:plane_friction])
    @test only(static.states)[static_friction.force_variable] ≈
        -2.0 atol = 1.0e-6
    @test only(static.states)[static_friction.shear_variable] ≈
        2.0 / static_friction.stiffness atol = 1.0e-6

    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        read(path, String), "constraint = \"support\"" =>
            "constraint = \"missing\"")))
end

@testset "Sim2D friction builders" begin
    Sim2D = PracticalMechanicalSimulation.Sim2D
    model = Sim2D.Model(:friction_builders)
    bearing = Sim2D.revolute_friction!(model, :bearing;
        joint = :pin, effective_radius = 0.02, stiffness = 10.0,
        static_coefficient = 0.8, dynamic_coefficient = 0.6)
    guide = Sim2D.translational_friction!(model, :guide_friction;
        joint = :guide, stiffness = 100.0,
        static_coefficient = 0.8, dynamic_coefficient = 0.6)
    plane = Sim2D.inplane_friction!(model, :plane_friction;
        constraint = :support, stiffness = 100.0,
        static_coefficient = 0.8, dynamic_coefficient = 0.6)
    specification = Sim2D.document(model)
    @test specification[bearing.name]["type"] == "revolute_friction"
    @test specification[guide.name]["type"] == "translational_friction"
    @test specification[plane.name]["type"] == "inplane_friction"
end
