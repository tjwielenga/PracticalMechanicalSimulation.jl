using PracticalMechanicalSimulation
using PracticalMechanicalSimulation.AutomaticAnalysis
using PracticalMechanicalSimulation.SimulationRunner:
    initial_derivative, runtime_state_partition
using Test

@testset "Planar user equation components" begin
    path = joinpath(@__DIR__, "..", "..", "models", "planar",
        "controlled-revolute-pendulum.toml")
    loaded = load_planar_model(path)
    controller = loaded.equation_components[:controller]
    @test length(controller.state_indices) == 1
    @test length(controller.algebraic_indices) == 2
    variable_index(component, name) = only(variable.index for variable in
        loaded.layout.catalog.variables if variable.component == component &&
        variable.name == name)
    angle_error = variable_index(:controller, :angle_error)
    integral_error = variable_index(:controller, :integral_error)
    control_value = variable_index(:controller, :torque)
    angle = variable_index(:pin, :theta)
    omega = variable_index(:pin, :omega)
    @test loaded.analysis.degrees_of_freedom == 2

    result = run_planar_model(path; duration = 0.05, samples = 4)
    @test length(result.times) == 4
    @test abs(last(result.states)[integral_error]) > 1.0e-4
    torque = only(result.loaded.forces[:control_torque])
    for sample in result.states
        @test sample[torque.torque_variable] ≈ sample[control_value] atol = 1.0e-7
        @test sample[angle_error] ≈ -pi / 4 - sample[angle] atol = 1.0e-7
    end

    state = first(result.states)
    derivative = initial_derivative(state, result.loaded, 0.0)
    partition = runtime_state_partition(result.loaded)
    selection = AnalysisSelection(Dynamics(),
        result.loaded.active_variable_indices, partition.equation_indices)
    jacobian = evaluate_analysis_sparse_jacobian(result.loaded.model,
        selection, 0.0, state, derivative, 2.0)
    state_row = findfirst(==(only(controller.state_equations)),
        selection.equation_indices)
    state_column = findfirst(==(integral_error), selection.variable_indices)
    @test jacobian[state_row, state_column] ≈ 2.0
    error_column = findfirst(==(angle_error), selection.variable_indices)
    @test jacobian[state_row, error_column] ≈ -1.0

    mktempdir() do directory
        filename = joinpath(directory, "controlled.simp")
        write_result(filename, result)
        stored = read_result(filename)
        @test ("controller", "integral_error") in
            collect(zip(stored.variable_components, stored.variable_names))
        restarted_source = read(path, String) * "\n[initial_conditions]\n" *
            "result = \"$filename\"\nsample = \"last\"\n"
        restarted = load_planar_model(IOBuffer(restarted_source))
        @test restarted.initial_values[
            only(restarted.equation_components[:controller].state_indices)] ≈
            last(result.states)[integral_error] atol = 1.0e-6
    end

    source = read(path, String)
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "state_equations = { integral_error = \"der(integral_error) = angle_error\" }" =>
        "state_equations = {}")))
    @test_throws ArgumentError load_planar_model(IOBuffer(replace(source,
        "inputs = { angle = \"pin.theta\", omega = \"pin.omega\" }" =>
        "inputs = { angle = \"pin.theta\", omega = \"missing.omega\" }")))

    settled = run_planar_model(path)
    final = last(settled.states)
    @test final[angle] ≈ -pi / 4 atol = 2.0e-4
    @test abs(final[omega]) < 2.0e-4
    @test final[angle_error] ≈ 0.0 atol = 2.0e-4
    @test final[control_value] ≈ 4.905 / sqrt(2) atol = 1.0e-3
end

@testset "Planar user equations in static equilibrium" begin
    source = """
        [model]
        dimension = "planar"
        [analysis]
        mode = "static"
        [ground]
        type = "ground"
        [ground.pin]
        type = "marker"
        [body]
        type = "rigid_body"
        mass = 1.0
        inertia = 1.0
        [body.pin]
        type = "marker"
        [fix]
        type = "fixed"
        markers = ["body.pin", "ground.pin"]
        [aux]
        type = "equation_component"
        states = { s = { initial = 0.0, static = "steady" } }
        variables = { y = { initial = 0.0 } }
        state_equations = { s = "1-s" }
        equations = ["y = 2*s"]
        """
    result = run_planar_model(IOBuffer(source); end_time = 0.0, samples = 1)
    component = result.loaded.equation_components[:aux]
    @test only(result.states)[only(component.state_indices)] ≈ 1.0
    @test only(result.states)[only(component.algebraic_indices)] ≈ 2.0

    held = replace(source, "static = \"steady\"" => "static = \"hold\"")
    held_result = run_planar_model(IOBuffer(held); end_time = 0.0, samples = 1)
    held_component = held_result.loaded.equation_components[:aux]
    @test only(held_result.states)[only(held_component.state_indices)] ≈ 0.0
    @test only(held_result.states)[only(held_component.algebraic_indices)] ≈ 0.0
end

@testset "Sim2D equation component builder" begin
    Sim2D = PracticalMechanicalSimulation.Sim2D
    model = Sim2D.Model(:equations)
    component = Sim2D.equation_component!(model, :aux;
        states = Dict(:s => Dict(:initial => 1.0, :static => :hold)),
        variables = Dict(:y => Dict(:initial => 0.0)),
        state_equations = Dict(:s => "der(s) = -s"),
        equations = ["y = 2*s"])
    specification = Sim2D.document(model)
    @test component.name == "aux"
    @test specification["aux"]["type"] == "equation_component"
    @test specification["aux"]["states"]["s"]["static"] == "hold"
end
