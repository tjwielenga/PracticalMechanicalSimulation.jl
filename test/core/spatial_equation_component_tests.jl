using PracticalMechanicalSimulation
using PracticalMechanicalSimulation.AutomaticAnalysis
using PracticalMechanicalSimulation.SpatialSimulationRunner:
    initial_spatial_derivative, spatial_state_masks,
    runtime_spatial_state_partition
using Test

@testset "Spatial user equation components" begin
    path = joinpath(@__DIR__, "..", "..", "models", "spatial",
        "controlled-revolute-pendulum.toml")
    loaded = load_spatial_model(path)
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
    @test loaded.initial_values[angle_error] ≈ -π / 4
    @test loaded.initial_values[control_value] ≈ -12π / 4
    torque = loaded.forces[:control_torque]
    @test loaded.initial_values[torque.magnitude_variable] ≈ -12π / 4
    partition = runtime_spatial_state_partition(loaded)
    differential, controlled, monitored = spatial_state_masks(loaded,
        partition)
    local_index = findfirst(==(only(controller.state_indices)),
        loaded.active_variable_indices)
    @test differential[local_index] && controlled[local_index] &&
        monitored[local_index]

    state = loaded.initial_values
    derivative = initial_spatial_derivative(state, loaded)
    selection = AnalysisSelection(Dynamics(), loaded.active_variable_indices,
        partition.equation_indices)
    jacobian = evaluate_analysis_sparse_jacobian(loaded.model, selection,
        0.0, state, derivative, 2.0)
    row = findfirst(==(only(controller.state_equations)),
        selection.equation_indices)
    column = findfirst(==(only(controller.state_indices)),
        selection.variable_indices)
    @test jacobian[row, column] ≈ 2.0
    omega_column = findfirst(==(omega), selection.variable_indices)
    error_column = findfirst(==(angle_error), selection.variable_indices)
    @test jacobian[row, error_column] ≈ -1.0
    @test jacobian[row, omega_column] ≈ 0.0
    row = findfirst(==(first(controller.algebraic_equations)),
        selection.equation_indices)
    angle_column = findfirst(==(angle), selection.variable_indices)
    @test jacobian[row, error_column] ≈ 1.0
    @test jacobian[row, angle_column] ≈ 1.0
    row = findfirst(==(last(controller.algebraic_equations)),
        selection.equation_indices)
    column = findfirst(==(control_value),
        selection.variable_indices)
    @test jacobian[row, column] ≈ 1.0
    @test jacobian[row, error_column] ≈ -12.0
    @test jacobian[row, findfirst(==(integral_error),
        selection.variable_indices)] ≈ -10.0
    @test jacobian[row, omega_column] ≈ 5.0

    result = run_spatial_model(path; duration = 0.05, samples = 4)
    @test length(result.times) == 4
    @test abs(last(result.states)[only(controller.state_indices)]) > 1.0e-4
    for sample in result.states
        @test sample[torque.magnitude_variable] ≈
            sample[control_value] atol = 1.0e-7
    end
    mktempdir() do directory
        filename = joinpath(directory, "controlled.simp")
        write_result(filename, result)
        stored = read_result(filename)
        @test ("controller", "integral_error") in
            collect(zip(stored.variable_components, stored.variable_names))
        restarted_source = read(path, String) * "\n[initial_conditions]\n" *
            "result = \"$filename\"\nsample = \"last\"\n"
        restarted = load_spatial_model(IOBuffer(restarted_source))
        @test restarted.initial_values[
            only(restarted.equation_components[:controller].state_indices)] ≈
            last(result.states)[only(controller.state_indices)] atol = 1.0e-6
        steady_restart = load_spatial_model(IOBuffer(replace(
            restarted_source, "static = \"hold\"" =>
                "static = \"steady\"")))
        @test steady_restart.initial_values[
            only(steady_restart.equation_components[:controller].state_indices)] ≈
            last(result.states)[only(controller.state_indices)] atol = 1.0e-6
    end

    source = read(path, String)
    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "state_equations = { integral_error = \"der(integral_error) = angle_error\" }" =>
        "state_equations = {}")))
    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "equations = [\"angle_error = target_angle - angle\", \"torque = kp*angle_error + ki*integral_error - kd*omega\"]" =>
        "equations = []")))
    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "static = \"hold\"" => "static = \"unknown\"")))
    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "state_equations = { integral_error = \"der(integral_error) = angle_error\" }" =>
        "state_equations = { integral_error = \"der(integral_error) = eval(omega)\" }")))
    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "inputs = { angle = \"pin.theta\", omega = \"pin.omega\" }" =>
        "inputs = { angle = \"pin.theta\", omega = \"missing.omega\" }")))
    @test_throws ArgumentError load_spatial_model(IOBuffer(replace(source,
        "target_angle = \"-45 deg\"" =>
        "target_angle = \"not an angle\"")))

    settled = run_spatial_model(path)
    final = last(settled.states)
    @test final[angle] ≈ -π / 4 atol = 2.0e-4
    @test abs(final[omega]) < 2.0e-4
    @test final[angle_error] ≈ 0.0 atol = 2.0e-4
    @test final[control_value] ≈ 4.905 / sqrt(2) atol = 1.0e-3

    modal_source = replace(source, "mode = \"dynamic\"" =>
        "mode = \"modal\"\nmodes = 2")
    modal = run_spatial_model(IOBuffer(modal_source))
    @test length(modal.eigenvalues) == 2
    @test length(modal.differential_columns) == 3
end

@testset "Lua-defined spatial equation component" begin
    source = """
        model { name = "lua_equations", dimension = "spatial" }
        ground { name = "ground" }
        rigid_body { name = "body", mass = 1.0,
            inertia = {1.0, 1.0, 1.0} }
        equation_component {
            name = "aux",
            states = { s = { initial = 2.0, static = "hold" } },
            variables = { y = { initial = 0.0 } },
            state_equations = { s = "der(s) = -s" },
            equations = { "y = 3*s" }
        }
        """
    loaded = load_spatial_model(IOBuffer(source); format = :lua)
    component = loaded.equation_components[:aux]
    @test loaded.initial_values[only(component.state_indices)] ≈ 2.0
    @test loaded.initial_values[only(component.algebraic_indices)] ≈ 6.0
end

@testset "User equation drives a rolling tire force" begin
    path = joinpath(@__DIR__, "..", "..", "models", "spatial",
        "driven-rolling-tire.toml")
    source = replace(read(path, String),
        "longitudinal_expression = \"tire.normal_force * longitudinal_stiffness_per_load * tire.slip_ratio\"" =>
        "longitudinal_expression = \"tire_law.fx\"") * """

        [tire_law]
        type = "equation_component"
        inputs = { normal = "tire.normal_force", slip = "tire.slip_ratio" }
        variables = { fx = { initial = 0.0, scale = 100.0 } }
        equations = ["fx = 12*normal*slip"]
        """
    result = run_spatial_model(IOBuffer(source); duration = 0.02,
        samples = 3)
    component = result.loaded.equation_components[:tire_law]
    tire = result.loaded.forces[:tire]
    for state in result.states
        # Stored output may be interpolated within a BDF step, so an
        # algebraic equality is only approximate at those samples.
        @test state[only(component.algebraic_indices)] ≈
            12 * state[tire.normal_force_variable] *
            state[tire.slip_ratio_variable] atol = 2.0e-4
    end
end

@testset "Spatial user equations in static equilibrium" begin
    source = """
        [model]
        dimension = "spatial"
        [analysis]
        mode = "static"
        [ground]
        type = "ground"
        [ground.pin]
        type = "marker"
        [body]
        type = "rigid_body"
        mass = 1.0
        inertia = [1.0, 1.0, 1.0]
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
    result = run_spatial_model(IOBuffer(source); end_time = 0.0,
        samples = 1)
    component = result.loaded.equation_components[:aux]
    @test only(result.states)[only(component.state_indices)] ≈ 1.0
    @test only(result.states)[only(component.algebraic_indices)] ≈ 2.0
    held = replace(source, "static = \"steady\"" => "static = \"hold\"")
    held_result = run_spatial_model(IOBuffer(held); end_time = 0.0,
        samples = 1)
    held_component = held_result.loaded.equation_components[:aux]
    @test only(held_result.states)[only(held_component.state_indices)] ≈ 0.0
    @test only(held_result.states)[only(held_component.algebraic_indices)] ≈ 0.0

    loop_source = replace(source,
        "states = { s = { initial = 0.0, static = \"steady\" } }" =>
            "states = {}",
        "variables = { y = { initial = 0.0 } }" =>
            "variables = { x = { initial = 0.0 }, y = { initial = 0.0 } }",
        "state_equations = { s = \"1-s\" }" =>
            "state_equations = {}",
        "equations = [\"y = 2*s\"]" =>
            "equations = [\"x + y = 1\", \"x - y = 0\"]")
    loop = load_spatial_model(IOBuffer(loop_source))
    @test loop.initial_values[loop.equation_components[:aux].algebraic_indices] ≈
        [0.5, 0.5]

    dynamic_source = replace(source, "mode = \"static\"" =>
        "mode = \"dynamic\"\ninitialization = \"static_equilibrium\"")
    dynamic = run_spatial_model(IOBuffer(dynamic_source);
        duration = 0.01, samples = 2)
    dynamic_component = dynamic.loaded.equation_components[:aux]
    @test first(dynamic.states)[only(dynamic_component.state_indices)] ≈ 1.0
    @test last(dynamic.states)[only(dynamic_component.state_indices)] ≈ 1.0
end
