using Test
using PracticalMechanicalSimulation

const SMOKE_PACKAGE_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

@testset "package smoke test" begin
    planar_path = joinpath(SMOKE_PACKAGE_ROOT, "models", "planar",
        "torsional-spring-pendulum.toml")
    planar_model = load_planar_model(planar_path)
    planar_result = run_planar_model(planar_model;
        duration = 0.01, samples = 2)
    @test planar_result.analysis_mode == :dynamic
    @test planar_result.times == [0.0, 0.01]
    @test all(isfinite, reduce(vcat, planar_result.states))

    spatial_path = joinpath(SMOKE_PACKAGE_ROOT, "models", "spatial",
        "revolute-pendulum.toml")
    spatial_model = load_spatial_model(spatial_path)
    spatial_result = run_spatial_model(spatial_model;
        duration = 0.01, samples = 2)
    @test spatial_result.analysis_mode == :dynamic
    @test spatial_result.times == [0.0, 0.01]
    @test all(isfinite, reduce(vcat, spatial_result.states))

    mktempdir() do directory
        result_path = joinpath(directory, "smoke.simp")
        @test write_result(result_path, spatial_result) == result_path
        stored = read_result(result_path)
        @test stored.status == :complete
        @test stored.times == [0.0, 0.01]
        @test stored.title == spatial_model.title
    end
end
