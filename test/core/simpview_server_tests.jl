module SimpViewServerTests

using HDF5
using Test

include("../../apps/SimpViewWeb/server/SimpViewServer.jl")

using .SimpViewServer

@testset "SimpView model service" begin
    model_path = joinpath(@__DIR__, "..", "..", "models", "spatial",
        "revolute-pendulum.toml")
    source = read(model_path, String)
    preview = preview_model_document(source, basename(model_path))
    @test preview["runnable"]
    @test preview["run_settings"]["end_time"] == 3.0
    @test preview["run_settings"]["frames_per_second"] == 60.0
    @test preview["run_settings"]["analysis_mode"] == "dynamic"

    planar_lua_path = joinpath(@__DIR__, "..", "..", "models", "planar",
        "lua-double-pendulum.lua")
    planar_lua_source = read(planar_lua_path, String)
    planar_lua_preview = preview_model_document(planar_lua_source,
        basename(planar_lua_path))
    @test planar_lua_preview["runnable"]
    @test planar_lua_preview["title"] ==
        "Reusable Sim2D Lua double pendulum"
    @test planar_lua_preview["run_settings"]["analysis_mode"] == "dynamic"

    initial = start_model_run(source, basename(model_path);
        analysis = "initial_conditions")
    initial_run = SimpViewServer.find_run(initial.id)
    wait(initial_run.task)
    initial_snapshot = model_run_snapshot(initial.id)
    @test initial.analysis_mode == :initial_conditions
    @test initial.samples == 1
    @test initial_snapshot.status == :complete
    @test initial_run.result.analysis_mode == :initial_conditions

    static = start_model_run(source, basename(model_path);
        analysis = "static")
    static_run = SimpViewServer.find_run(static.id)
    wait(static_run.task)
    static_snapshot = model_run_snapshot(static.id)
    @test static.analysis_mode == :static
    @test static.samples == 1
    @test static_snapshot.status == :complete

    continued = start_model_run(source, basename(model_path);
        analysis = "dynamic", initial_run_id = static.id,
        end_time = 0.01, frames_per_second = 100.0)
    continued_run = SimpViewServer.find_run(continued.id)
    wait(continued_run.task)
    continued_snapshot = model_run_snapshot(continued.id)
    @test continued_snapshot.status == :complete
    @test continued_run.starting_result.state == last(static_run.result.states)

    extended = start_model_run(source, basename(model_path);
        analysis = "dynamic", initial_run_id = continued.id,
        end_time = 0.02, frames_per_second = 100.0)
    extended_run = SimpViewServer.find_run(extended.id)
    wait(extended_run.task)
    extended_snapshot = model_run_snapshot(extended.id)
    @test extended_snapshot.status == :complete
    @test extended_run.starting_result.time == 0.01
    @test extended_run.result.times == [0.0, 0.01, 0.02]

    modal = start_model_run(source, basename(model_path); analysis = "modal",
        initial_run_id = static.id)
    modal_run = SimpViewServer.find_run(modal.id)
    wait(modal_run.task)
    modal_snapshot = model_run_snapshot(modal.id)
    @test modal_snapshot.status == :complete
    @test modal_run.starting_result.state == last(static_run.result.states)
    for body in values(modal_run.loaded.bodies)
        @test modal_run.result.operating_state[body.position_variables] ≈
            last(static_run.result.states)[body.position_variables]
        @test modal_run.result.operating_state[body.euler_parameter_variables] ≈
            last(static_run.result.states)[body.euler_parameter_variables]
    end
    @test modal_snapshot.document["choice_name"] == "Mode"
    @test length(modal_snapshot.document["choices"]) >= 1

    started = start_model_run(source, basename(model_path);
        end_time = 0.05, frames_per_second = 20.0)
    @test started.samples == 2
    run = SimpViewServer.find_run(started.id)
    wait(run.task)
    snapshot = model_run_snapshot(started.id)
    @test snapshot.status == :complete
    @test snapshot.result_ready
    @test snapshot.document["run_status"] == "complete"
    @test length(snapshot.document["choices"][1]["times"]) == 2
    h5open(run.result_path, "r") do file
        @test read(attributes(file)["status"]) == "complete"
        @test haskey(file, "graphics")
        @test eltype(file["results/values"]) == Float32
    end
    @test isnothing(model_run_snapshot(started.id, snapshot.revision))
end

end
