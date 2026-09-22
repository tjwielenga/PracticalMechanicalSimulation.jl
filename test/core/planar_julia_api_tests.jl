using Test
using PracticalMechanicalSimulation
using HDF5
using TOML

const Sim2DAPI = PracticalMechanicalSimulation.Sim2D

function api_planar_pendulum(; end_time = 0.02, frames_per_second = 100)
    model = Sim2DAPI.Model(:julia_api_planar_pendulum;
        title = "Julia API planar pendulum")
    Sim2DAPI.analysis!(model; mode = :automatic)
    Sim2DAPI.simulation!(model;
        start_time = 0.0,
        end_time,
        frames_per_second,
        relative_tolerance = 1.0e-7,
        absolute_tolerance = 1.0e-9,
        maximum_step = 0.01)
    Sim2DAPI.graphics!(model;
        background = :white,
        body_palette = :colorblind)

    ground = Sim2DAPI.ground!(model, :ground)
    ground_pin = Sim2DAPI.marker!(ground, :pin)
    pendulum = Sim2DAPI.rigid_body!(model, :pendulum;
        mass = 1.0,
        inertia = 1 / 12,
        position = [0.5, 0.0],
        velocity = [0.0, 0.2],
        angular_velocity = 0.4)
    body_pin = Sim2DAPI.marker!(pendulum, :pin;
        position = [-0.5, 0.0])
    tip = Sim2DAPI.marker!(pendulum, :tip;
        position = [0.5, 0.0])
    Sim2DAPI.graphics!(pendulum; show_default = false, color = :steelblue)
    Sim2DAPI.graphic!(pendulum, :member;
        shape = :cylinder, markers = [body_pin, tip], radius = 0.04)

    pin = Sim2DAPI.revolute!(model, :pin;
        markers = [body_pin, ground_pin],
        rotation_coordinates = true)
    Sim2DAPI.gravity!(model, :gravity;
        acceleration = [0.0, -9.81],
        bodies = [pendulum])
    Sim2DAPI.state_selection!(model;
        method = :preferred,
        preferred_velocities = [Sim2DAPI.variable(pin, :omega)],
        allow_fallback = false)
    model
end

@testset "planar Julia API" begin
    model = api_planar_pendulum()
    specification = Sim2DAPI.document(model)
    @test specification["pendulum"]["pin"]["type"] == "marker"
    @test specification["pin"]["markers"] ==
        ["pendulum.pin", "ground.pin"]
    @test specification["state_selection"]["preferred_velocities"] ==
        ["pin.omega"]
    @test specification["simulation"]["output_samples"] == 3

    specification["model"]["title"] = "changed"
    @test Sim2DAPI.document(model)["model"]["title"] ==
        "Julia API planar pendulum"

    loaded = load_planar_model(model)
    @test loaded.title == "Julia API planar pendulum"
    @test Set(keys(loaded.bodies)) == Set([:pendulum])
    @test Set(keys(loaded.connections)) == Set([:pin])
    @test loaded.model_source ==
        load_planar_model(Sim2DAPI.document(model)).model_source

    result = Sim2DAPI.run(model)
    @test result.analysis_mode == :dynamic
    @test result.times == [0.0, 0.01, 0.02]
    @test all(isfinite, reduce(vcat, result.states))

    mktempdir() do directory
        model_path = joinpath(directory, "pendulum.toml")
        result_path = joinpath(directory, "pendulum.simp")
        @test Sim2DAPI.write_model(model_path, model) == model_path
        reloaded = load_planar_model(model_path)
        @test TOML.parse(reloaded.model_source) == TOML.parse(loaded.model_source)
        Sim2DAPI.save_result(result_path, result)
        stored = read_result(result_path)
        @test stored.status == :complete
        @test stored.title == "Julia API planar pendulum"
        @test stored.model_source == loaded.model_source
        h5open(result_path, "r") do file
            @test haskey(file, "graphics")
        end
    end

    other = Sim2DAPI.Model(:other)
    other_ground = Sim2DAPI.ground!(other, :ground)
    other_pin = Sim2DAPI.marker!(other_ground, :pin)
    @test_throws ArgumentError Sim2DAPI.revolute!(model, :bad;
        markers = [other_pin, "ground.pin"])
    @test_throws ArgumentError Sim2DAPI.rigid_body!(model, :pendulum;
        mass = 1.0, inertia = 1.0)

    symbolic_document = Dict(
        :model => Dict(:name => :symbolic, :dimension => :planar),
        :ground => Dict(:type => :ground),
        :body => Dict(:type => :rigid_body, :mass => 1.0,
            :inertia => 1.0))
    symbolic = load_planar_model(symbolic_document)
    @test symbolic.title == "symbolic"

    hierarchical = Sim2DAPI.Model(:hierarchical)
    Sim2DAPI.ground!(hierarchical, :ground)
    nested_body = Sim2DAPI.rigid_body!(hierarchical, "assembly.link";
        mass = 1.0, inertia = 1.0)
    nested_tip = Sim2DAPI.marker!(nested_body, :tip;
        position = [1.0, 0.0])
    nested_document = Sim2DAPI.document(hierarchical)
    @test nested_tip.name == "assembly.link.tip"
    @test nested_document["assembly"]["link"]["tip"]["type"] == "marker"
    nested_loaded = load_planar_model(hierarchical)
    @test haskey(nested_loaded.bodies, Symbol("assembly.link"))
    @test haskey(nested_loaded.markers, Symbol("assembly.link.tip"))

    contact_model = Sim2DAPI.Model(:api_curve_contact)
    contact_ground = Sim2DAPI.ground!(contact_model, :ground)
    profile_frame = Sim2DAPI.marker!(contact_ground, :profile_frame)
    profile = Sim2DAPI.curve!(contact_model, :profile;
        marker = profile_frame,
        points = [[0.5, 0.0], [0.0, 0.5],
                  [-0.5, 0.0], [0.0, -0.5]])
    roller = Sim2DAPI.rigid_body!(contact_model, :roller;
        mass = 1.0, inertia = 0.01, position = [0.0, 0.58])
    roller_center = Sim2DAPI.marker!(roller, :center)
    Sim2DAPI.curve_contact!(contact_model, :contact;
        curve = profile, roller_marker = roller_center,
        radius = 0.1, stiffness = 10_000.0)
    contact_loaded = load_planar_model(contact_model)
    @test only(contact_loaded.forces[:contact]) isa
        PracticalMechanicalSimulation.PlanarCurveContacts.PlanarCurveContactComponent

    flat_contact_model = Sim2DAPI.Model(:api_flat_follower_contact)
    flat_ground = Sim2DAPI.ground!(flat_contact_model, :ground)
    flat_frame = Sim2DAPI.marker!(flat_ground, :profile_frame)
    flat_profile = Sim2DAPI.curve!(flat_contact_model, :profile;
        marker = flat_frame,
        points = [[0.5, 0.0], [0.0, 0.5],
                  [-0.5, 0.0], [0.0, -0.5]])
    follower = Sim2DAPI.rigid_body!(flat_contact_model, :follower;
        mass = 1.0, inertia = 0.01, position = [0.0, 0.48])
    face = Sim2DAPI.marker!(follower, :face)
    Sim2DAPI.flat_follower_contact!(flat_contact_model, :contact;
        curve = flat_profile, follower_marker = face,
        stiffness = 10_000.0)
    flat_contact_loaded = load_planar_model(flat_contact_model)
    @test only(flat_contact_loaded.forces[:contact]).follower_kind == :flat
end
