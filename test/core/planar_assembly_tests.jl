using PracticalMechanicalSimulation
using Test
using TOML

const PLANAR_ASSEMBLY_MODEL = normpath(joinpath(@__DIR__, "..", "..",
    "models", "planar", "lua-double-pendulum.lua"))

@testset "Planar Lua assemblies" begin
    loaded = load_planar_model(PLANAR_ASSEMBLY_MODEL)
    expanded = TOML.parse(loaded.model_source)

    @test loaded.title == "Reusable Sim2D Lua double pendulum"
    @test Set(keys(loaded.bodies)) == Set((
        Symbol("pendulum.first_link"), Symbol("pendulum.second_link")))
    @test haskey(loaded.markers, Symbol("pendulum.first_link.inner"))
    @test haskey(loaded.markers, Symbol("pendulum.second_link.outer"))
    @test haskey(loaded.connections, Symbol("pendulum.base_joint"))
    @test haskey(loaded.connections, Symbol("pendulum.elbow_joint"))
    @test expanded["pendulum"]["first_link"]["graphics"]["member"][
        "shape"] == "cylinder"
    @test !occursin("type = \"double_pendulum\"", loaded.model_source)

    result = run_planar_model(loaded; duration = 0.02, samples = 3)
    @test result.analysis_mode == :dynamic
    @test result.times == [0.0, 0.01, 0.02]
    @test all(isfinite, reduce(vcat, result.states))

    source = """
        local sim2d = require "sim2d"
        sim2d.model {name = "local_geometry", dimension = "planar"}
        sim2d.ground {name = "ground"}
        sim2d.rigid_body {
            name = "body", mass = 1.0, inertia = 1.0,
            position = {1.0, 2.0}, orientation = "90 deg"
        }
        sim2d.marker {
            name = "body.point",
            position = sim2d.local_point("body", {1.0, 3.0}),
            orientation = sim2d.local_orientation("body", "180 deg")
        }
        """
    geometry = load_planar_model(IOBuffer(source); format = :lua)
    geometry_document = TOML.parse(geometry.model_source)
    @test geometry_document["body"]["point"]["position"] ≈ [1.0, 0.0]
    @test geometry_document["body"]["point"]["orientation"] ≈ pi / 2

    assembly_source = read(joinpath(@__DIR__, "..", "..", "assemblies",
        "planar", "double_pendulum.lua"), String)
    assembly_error = try
        load_planar_model(IOBuffer(assembly_source); format = :lua)
        nothing
    catch error
        error
    end
    @test assembly_error isa ArgumentError
    @test occursin("returns an assembly module and does not define a model",
        sprint(showerror, assembly_error))

    output = IOBuffer()
    @test PracticalMechanicalSimulation.CommandLine.planar_model_main(
        [PLANAR_ASSEMBLY_MODEL, "0.02", "3"]; output) == 0
    @test occursin("Reusable Sim2D Lua double pendulum",
        String(take!(output)))
end
