using Test
using PracticalMechanicalSimulation
using PracticalMechanicalSimulation.AutomaticAnalysis
using PracticalMechanicalSimulation.SpatialCurveContacts

@testset "Spatial extruded curve contacts" begin
    source = """
    [model]
    name = "spatial_curve_contact_test"
    dimension = "spatial"

    [ground]
    type = "ground"

    [ground.frame]
    type = "marker"

    [profile]
    type = "curve"
    marker = "ground.frame"
    half_width = 0.08
    points = [[0.5, 0.0], [0.35, 0.35], [0.0, 0.5], [-0.35, 0.35], [-0.5, 0.0], [-0.35, -0.35], [0.0, -0.5], [0.35, -0.35]]

    [roller]
    type = "rigid_body"
    mass = 1.0
    inertia = [0.01, 0.01, 0.01]
    position = [0.0, 0.58, 0.03]

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
    loaded = load_spatial_model(IOBuffer(source))
    contact = loaded.forces[:contact]
    @test contact isa SpatialCurveContactComponent
    @test contact.curve_half_width == 0.08
    values = spatial_curve_contact_values(contact, loaded.initial_values)
    @test values.tangent_error ≈ 0.0 atol = 1.0e-10
    @test values.gap ≈ -0.02 atol = 1.0e-10
    @test values.normal ≈ [0.0, 1.0, 0.0] atol = 1.0e-10
    @test values.contact_point[3] ≈ 0.03 atol = 1.0e-12
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
    expression_model = load_spatial_model(IOBuffer(expression_source))
    expression_contact = expression_model.forces[:contact]
    @test expression_contact.expression
    @test expression_model.initial_values[
        expression_contact.normal_force_variable] ≈ 24.68 atol = 1.0e-8

    flat_source = replace(source,
        "position = [0.0, 0.58, 0.03]" =>
            "position = [0.0, 0.48, 0.03]",
        "type = \"curve_contact\"" =>
            "type = \"flat_follower_contact\"",
        "roller_marker = \"roller.center\"\nradius = 0.1" =>
            "follower_marker = \"roller.center\"")
    flat_model = load_spatial_model(IOBuffer(flat_source))
    flat_contact = flat_model.forces[:contact]
    @test flat_contact.follower_kind == :flat
    flat_values = spatial_curve_contact_values(flat_contact,
        flat_model.initial_values)
    @test flat_values.tangent_error ≈ 0.0 atol = 1.0e-10
    @test flat_values.gap ≈ -0.02 atol = 1.0e-10
    @test flat_values.normal ≈ [0.0, 1.0, 0.0] atol = 1.0e-10

    misaligned = replace(source,
        "[roller.center]\ntype = \"marker\"" =>
        "[roller.center]\ntype = \"marker\"\norientation = [\"10 deg\", 1.0, 0.0, 0.0]")
    @test_throws ArgumentError load_spatial_model(IOBuffer(misaligned))

    lua_source = """
    local sim3d = require "sim3d"
    sim3d.model {name = "lua_spatial_curve", dimension = "spatial"}
    sim3d.ground {name = "ground"}
    sim3d.marker {name = "ground.frame"}
    sim3d.curve {
        name = "profile", marker = "ground.frame", half_width = 0.08,
        points = {{0.5, 0.0}, {0.0, 0.5}, {-0.5, 0.0}, {0.0, -0.5}}
    }
    sim3d.rigid_body {
        name = "roller", mass = 1.0, inertia = {0.01, 0.01, 0.01},
        position = {0.0, 0.58, 0.0}
    }
    sim3d.marker {name = "roller.center"}
    sim3d.curve_contact {
        name = "contact", curve = "profile",
        roller_marker = "roller.center", radius = 0.1,
        stiffness = 10000.0
    }
    """
    lua_model = load_spatial_model(IOBuffer(lua_source); format = :lua,
        source_directory = normpath(joinpath(
            @__DIR__, "..", "..", "assemblies")))
    @test lua_model.forces[:contact] isa SpatialCurveContactComponent
end
