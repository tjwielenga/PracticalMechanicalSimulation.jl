using Test
using PracticalMechanicalSimulation
using PracticalMechanicalSimulation.PlanarComponentAssembly

const PLANAR_STAGE_MODEL_DIRECTORY = normpath(joinpath(
    @__DIR__, "..", "..", "models", "planar"))

stage_model(name) = joinpath(PLANAR_STAGE_MODEL_DIRECTORY, name)
stage_force(loaded, name) = only(loaded.forces[name])

function set_planar_stage!(loaded, stage, state; time = 0.0)
    PracticalMechanicalSimulation.SimulationRunner.set_planar_analysis_stage!(
        loaded, stage; state, time)
end

@testset "Planar stage-dependent forces" begin
    applied_source = replace(read(stage_model("marker-directed-force.toml"),
        String), "force = 1.0" =>
            "force = 1.0\nactive_during = \"static\"")
    applied = load_planar_model(IOBuffer(applied_source))
    push = stage_force(applied, :push)
    state = copy(applied.initial_values)
    @test push.active_during == (:static,)
    @test !push.active[]
    @test push.magnitude_variable != 0
    @test state[push.magnitude_variable] == 0.0
    set_planar_stage!(applied, :static, state)
    @test push.active[]
    @test state[push.magnitude_variable] == 1.0
    set_planar_stage!(applied, :dynamic, state)
    @test !push.active[]
    @test state[push.magnitude_variable] == 0.0
    applied_result = run_planar_model(IOBuffer(applied_source);
        duration = 0.01, samples = 2)
    result_push = stage_force(applied_result.loaded, :push)
    @test all(sample[result_push.magnitude_variable] == 0.0
        for sample in applied_result.states)

    torque_source = replace(read(
        stage_model("revolute-bearing-friction.toml"), String),
        "expression = \"2*max(0,min(1,(t-1.0)/0.4))\"" =>
            "torque = 2.0\nactive_during = \"static\"")
    torque_model = load_planar_model(IOBuffer(torque_source))
    torque = stage_force(torque_model, :drive)
    torque_state = copy(torque_model.initial_values)
    @test !torque.active[]
    @test torque_state[torque.torque_variable] == 0.0
    set_planar_stage!(torque_model, :static, torque_state)
    @test torque.active[]
    @test torque_state[torque.torque_variable] == 2.0
    set_planar_stage!(torque_model, :modal, torque_state)
    @test !torque.active[]
    @test torque_state[torque.torque_variable] == 0.0

    spanning_source = replace(read(
        stage_model("spanning-spring-pendulum.toml"), String),
        "free_length = 0.5" =>
            "free_length = 0.5\nactive_during = \"static\"")
    spanning_model = load_planar_model(IOBuffer(spanning_source))
    spring = stage_force(spanning_model, :spring)
    spanning_state = copy(spanning_model.initial_values)
    @test !spring.active[]
    @test spanning_state[spring.element.force_variable] == 0.0
    @test spanning_state[spring.element.global_force_variables] ≈ zeros(2)
    @test spanning_state[spring.element.length_variable] > 0.0
    set_planar_stage!(spanning_model, :static, spanning_state)
    @test spring.active[]
    @test spanning_state[spring.element.force_variable] != 0.0
    set_planar_stage!(spanning_model, :dynamic, spanning_state)
    @test spanning_state[spring.element.force_variable] == 0.0
    @test spanning_state[spring.element.global_force_variables] ≈ zeros(2)

    bushing_source = replace(read(
        stage_model("bushing-supported-body.toml"), String),
        "free_angle = 0.0" =>
            "free_angle = 0.0\nactive_during = \"static\"")
    bushing_model = load_planar_model(IOBuffer(bushing_source))
    support = stage_force(bushing_model, :support)
    bushing_state = copy(bushing_model.initial_values)
    @test !support.active[]
    set_planar_stage!(bushing_model, :static, bushing_state)
    @test support.active[]
    @test bushing_state[support.torque_variable] != 0.0
    set_planar_stage!(bushing_model, :dynamic, bushing_state)
    @test bushing_state[support.force_variables] ≈ zeros(2)
    @test bushing_state[support.torque_variable] == 0.0

    contact_source = replace(read(stage_model("bouncing-ball.toml"), String),
        "damping_factor = 0.15" =>
            "damping_factor = 0.15\ninactive_during = \"static\"")
    contact_model = load_planar_model(IOBuffer(contact_source))
    contact = stage_force(contact_model, :floor_contact)
    contact_state = copy(contact_model.initial_values)
    contact_state[contact_model.bodies[:ball].position_variables[2]] = 0.05
    set_planar_stage!(contact_model, :static, contact_state)
    @test !contact.active[]
    @test contact_state[contact.normal_force_variable] == 0.0
    set_planar_stage!(contact_model, :dynamic, contact_state)
    @test contact.active[]
    @test contact_state[contact.normal_force_variable] > 0.0

    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        applied_source, "active_during = \"static\"" =>
            "active_during = \"dynamic\"\ninactive_during = \"static\"")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        applied_source, "active_during = \"static\"" =>
            "active_during = \"assembly\"")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        applied_source, "active_during = \"static\"" =>
            "active_during = []")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(
        read(stage_model("bushing-supported-body.toml"), String),
        "acceleration = [0.0, -9.81]" =>
            "acceleration = [0.0, -9.81]\nactive_during = \"static\"")))

    staged_result = run_planar_model(
        stage_model("stage-dependent-force-drop.toml");
        duration = 0.01, samples = 2)
    @test getproperty.(staged_result.static_progress_events, :status) ==
        [:entered, :consistent, :converged]
    @test all(event.model_time == 0.0
        for event in staged_result.static_progress_events)
end
