using Test
using PracticalMechanicalSimulation
using PracticalMechanicalSimulation.AutomaticAnalysis
using PracticalMechanicalSimulation.PlanarCurveContacts
using PracticalMechanicalSimulation.PlanarComponentAssembly

@testset "Planar closed curves and roller contact" begin
    points = [[0.5, 0.0], [0.35, 0.35], [0.0, 0.5], [-0.35, 0.35],
              [-0.5, 0.0], [-0.35, -0.35], [0.0, -0.5], [0.35, -0.35]]
    curve = PlanarClosedCurve(points)
    start = curve_point(curve, 0.0)
    finish = curve_point(curve, curve.length)
    @test start.position ≈ finish.position atol = 1.0e-12
    @test start.derivative ≈ finish.derivative atol = 1.0e-12
    @test start.second_derivative ≈ finish.second_derivative atol = 1.0e-12

    source = """
    [model]
    name = "curve_contact_test"
    dimension = "planar"

    [simulation]
    end_time = 0.01
    output_samples = 2

    [ground]
    type = "ground"

    [ground.frame]
    type = "marker"

    [profile]
    type = "curve"
    marker = "ground.frame"
    points = [[0.5, 0.0], [0.35, 0.35], [0.0, 0.5], [-0.35, 0.35], [-0.5, 0.0], [-0.35, -0.35], [0.0, -0.5], [0.35, -0.35]]

    [roller]
    type = "rigid_body"
    mass = 1.0
    inertia = 0.01
    position = [0.0, 0.58]

    [roller.center]
    type = "marker"

    [contact]
    type = "curve_contact"
    curve = "profile"
    roller_marker = "roller.center"
    radius = 0.1
    stiffness = 10000.0
    damping_factor = 0.1
    transition_depth = 0.001
    """
    loaded = load_planar_model(IOBuffer(source))
    contact = only(loaded.forces[:contact])
    @test contact isa PlanarCurveContactComponent
    values = curve_contact_values(contact, loaded.initial_values)
    @test values.tangent_error ≈ 0.0 atol = 1.0e-10
    @test values.gap ≈ -0.02 atol = 1.0e-10
    @test values.normal ≈ [0.0, 1.0] atol = 1.0e-10
    @test loaded.initial_values[contact.normal_force_variable] > 0

    derivative = zeros(length(loaded.initial_values))
    selection = AnalysisSelection(Dynamics(),
        collect(eachindex(loaded.initial_values)),
        collect(eachindex(loaded.layout.catalog.equations)))
    residual = zeros(length(selection.equation_indices))
    evaluate_analysis_equations!(residual, loaded.model, selection, 0.0,
        loaded.initial_values, derivative)
    jacobian = evaluate_analysis_jacobian(loaded.model, selection, 0.0,
        loaded.initial_values, derivative, 0.0)
    @test all(isfinite, residual)
    @test all(isfinite, jacobian)

    expression_source = replace(source,
        "stiffness = 10000.0\ndamping_factor = 0.1\ntransition_depth = 0.001" =>
        "expression = \"max(0, -contact.gap) * 1234 * max(0, contact.normal_y)\"")
    expression_model = load_planar_model(IOBuffer(expression_source))
    expression_contact = only(expression_model.forces[:contact])
    @test expression_contact.expression
    @test expression_model.initial_values[
        expression_contact.normal_force_variable] ≈ 24.68 atol = 1.0e-8

    lua_source = """
    local sim2d = require "sim2d"
    sim2d.model {name = "lua_curve_contact", dimension = "planar"}
    sim2d.ground {name = "ground"}
    sim2d.marker {name = "ground.frame"}
    sim2d.curve {
        name = "profile", marker = "ground.frame",
        points = {{0.5, 0.0}, {0.0, 0.5}, {-0.5, 0.0}, {0.0, -0.5}}
    }
    sim2d.rigid_body {
        name = "roller", mass = 1.0, inertia = 0.01,
        position = {0.0, 0.58}
    }
    sim2d.marker {name = "roller.center"}
    sim2d.curve_contact {
        name = "contact", curve = "profile",
        roller_marker = "roller.center", radius = 0.1,
        stiffness = 10000.0
    }
    """
    lua_model = load_planar_model(IOBuffer(lua_source); format = :lua)
    @test only(lua_model.forces[:contact]) isa PlanarCurveContactComponent
end

@testset "Planar plane-contact expression parity" begin
    source = """
    [model]
    name = "plane_expression_test"
    dimension = "planar"

    [ground]
    type = "ground"

    [ground.plane]
    type = "marker"

    [ball]
    type = "rigid_body"
    mass = 1.0
    inertia = 0.01
    position = [0.0, 0.08]

    [ball.center]
    type = "marker"

    [contact]
    type = "plane_contact"
    markers = ["ball.center", "ground.plane"]
    radius = 0.1
    expression = "max(0, -contact.gap) * 2000"
    """
    loaded = load_planar_model(IOBuffer(source))
    contact = only(loaded.forces[:contact])
    @test contact.expression
    @test loaded.initial_values[contact.normal_force_variable] ≈ 40.0
    bad = replace(source,
        "expression = \"max(0, -contact.gap) * 2000\"" =>
        "expression = \"max(0, -contact.gap) * 2000\"\ndamping_factor = 0.1")
    @test_throws ArgumentError load_planar_model(IOBuffer(bad))
end
