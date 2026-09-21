using Test
using PracticalMechanicalSimulation
using HDF5

const Sim3DAPI = PracticalMechanicalSimulation.Sim3D

function api_pendulum(; end_time = 0.02, frames_per_second = 100)
    model = Sim3DAPI.Model(:julia_api_pendulum;
        title = "Julia API pendulum")
    Sim3DAPI.analysis!(model; mode = :automatic)
    Sim3DAPI.simulation!(model;
        start_time = 0.0,
        end_time,
        frames_per_second,
        relative_tolerance = 1.0e-7,
        absolute_tolerance = 1.0e-9,
        maximum_step = 0.01)
    Sim3DAPI.graphics!(model;
        background = :white,
        body_palette = :colorblind)

    ground = Sim3DAPI.ground!(model, :ground)
    ground_pin = Sim3DAPI.marker!(ground, :pin)
    Sim3DAPI.graphics!(ground_pin;
        shape = :xy_frame,
        axis_length = 0.35,
        plane_size = 0.16)

    pendulum = Sim3DAPI.rigid_body!(model, :pendulum;
        mass = 1.0,
        inertia = [0.001, 1 / 12, 1 / 12],
        position = [0.5, 0.0, 0.0],
        velocity = [0.0, 0.2, 0.0],
        angular_velocity = [0.0, 0.0, 0.4])
    Sim3DAPI.graphics!(pendulum;
        shape = :box,
        size = [1.0, 0.08, 0.08],
        color = :steelblue)
    body_pin = Sim3DAPI.marker!(pendulum, :pin;
        position = [-0.5, 0.0, 0.0])

    pin = Sim3DAPI.revolute!(model, :pin;
        markers = [body_pin, ground_pin],
        rotation_coordinates = true)
    Sim3DAPI.gravity!(model, :gravity;
        acceleration = [0.0, -9.81, 0.0],
        bodies = [pendulum])
    Sim3DAPI.state_selection!(model;
        method = :preferred,
        preferred_velocities = [Sim3DAPI.variable(pin, :omega)],
        allow_fallback = false)
    model
end

@testset "spatial Julia API" begin
    model = api_pendulum()
    specification = Sim3DAPI.document(model)
    @test specification["pendulum"]["pin"]["type"] == "marker"
    @test specification["pin"]["markers"] ==
        ["pendulum.pin", "ground.pin"]
    @test specification["state_selection"]["preferred_velocities"] ==
        ["pin.omega"]
    @test specification["simulation"]["output_samples"] == 3

    # `document` returns a detached value rather than exposing builder state.
    specification["model"]["title"] = "changed"
    @test Sim3DAPI.document(model)["model"]["title"] ==
        "Julia API pendulum"

    loaded = load_spatial_model(model)
    @test loaded.title == "Julia API pendulum"
    @test Set(keys(loaded.bodies)) == Set([:pendulum])
    @test Set(keys(loaded.connections)) == Set([:pin])
    @test loaded.model_source ==
        load_spatial_model(Sim3DAPI.document(model)).model_source

    result = Sim3DAPI.run(model)
    @test result.analysis_mode == :dynamic
    @test result.times == [0.0, 0.01, 0.02]
    @test all(isfinite, reduce(vcat, result.states))

    mktempdir() do directory
        model_path = joinpath(directory, "pendulum.toml")
        result_path = joinpath(directory, "pendulum.simp")
        @test Sim3DAPI.write_model(model_path, model) == model_path
        reloaded = load_spatial_model(model_path)
        @test reloaded.model_source == loaded.model_source
        Sim3DAPI.save_result(result_path, result)
        stored = read_result(result_path)
        @test stored.status == :complete
        @test stored.title == "Julia API pendulum"
        @test stored.model_source == loaded.model_source
        h5open(result_path, "r") do file
            @test haskey(file, "graphics")
        end
    end

    other = Sim3DAPI.Model(:other)
    other_ground = Sim3DAPI.ground!(other, :ground)
    other_pin = Sim3DAPI.marker!(other_ground, :pin)
    @test_throws ArgumentError Sim3DAPI.revolute!(model, :bad;
        markers = [other_pin, "ground.pin"])
    @test_throws ArgumentError Sim3DAPI.rigid_body!(model, :pendulum;
        mass = 1.0, inertia = [1.0, 1.0, 1.0])

    symbolic_document = Dict(
        :model => Dict(:name => :symbolic, :dimension => :spatial),
        :ground => Dict(:type => :ground),
        :body => Dict(:type => :rigid_body, :mass => 1.0,
            :inertia => [1.0, 1.0, 1.0]))
    symbolic = load_spatial_model(symbolic_document)
    @test symbolic.title == "symbolic"
    @test symbolic.grounds == Set([:ground])

    hierarchical = Sim3DAPI.Model(:hierarchical)
    Sim3DAPI.ground!(hierarchical, :ground)
    nested_body = Sim3DAPI.rigid_body!(hierarchical, "assembly.link";
        mass = 1.0, inertia = [1.0, 1.0, 1.0])
    nested_tip = Sim3DAPI.marker!(nested_body, :tip;
        position = [1.0, 0.0, 0.0])
    nested_document = Sim3DAPI.document(hierarchical)
    @test nested_tip.name == "assembly.link.tip"
    @test nested_document["assembly"]["link"]["tip"]["type"] == "marker"
    nested_loaded = load_spatial_model(hierarchical)
    @test haskey(nested_loaded.bodies, Symbol("assembly.link"))
    @test haskey(nested_loaded.markers, Symbol("assembly.link.tip"))
end
