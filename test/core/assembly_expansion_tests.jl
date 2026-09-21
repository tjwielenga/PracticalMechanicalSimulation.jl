using PracticalMechanicalSimulation
using PracticalMechanicalSimulation.SpatialModeling
using Test
using TOML

const ASSEMBLY_MODEL_DIRECTORY = normpath(joinpath(
    @__DIR__, "..", "..", "models", "spatial"))

function lua_model_error(source)
    try
        PracticalMechanicalSimulation.AssemblyExpansion.parse_lua_model(source;
            element_types =
                PracticalMechanicalSimulation.SpatialModelIO.SPATIAL_ELEMENT_TYPES)
        nothing
    catch error
        error
    end
end

@testset "Modeling assemblies" begin
    @testset "shared Lua validation helpers" begin
        missing = lua_model_error("""
            local sim3d = require "sim3d"
            sim3d.required({}, "radius", "test pulley")
            """)
        @test missing isa ArgumentError
        @test occursin("test pulley requires 'radius'", sprint(showerror, missing))

        invalid_positive = lua_model_error("""
            local sim3d = require "sim3d"
            sim3d.positive({radius = 0}, "radius", "test pulley")
            """)
        @test invalid_positive isa ArgumentError
        @test occursin("test pulley 'radius' must be positive",
            sprint(showerror, invalid_positive))

        invalid_nonnegative = lua_model_error("""
            local sim3d = require "sim3d"
            sim3d.nonnegative({damping = -1}, "damping", "test bushing")
            """)
        @test invalid_nonnegative isa ArgumentError
        @test occursin("test bushing 'damping' must be nonnegative",
            sprint(showerror, invalid_nonnegative))

        assembly_source = read(joinpath(@__DIR__, "..", "..", "assemblies",
            "spatial", "gylt_245_75_r16.lua"), String)
        assembly_as_model = lua_model_error(assembly_source)
        @test assembly_as_model isa ArgumentError
        @test occursin("returns an assembly module and does not define a model",
            sprint(showerror, assembly_as_model))
    end

    path = joinpath(ASSEMBLY_MODEL_DIRECTORY,
        "leaf-spring-assembly.lua")
    loaded = load_spatial_model(path)
    expanded = TOML.parse(loaded.model_source)
    toml_loaded = load_spatial_model(joinpath(ASSEMBLY_MODEL_DIRECTORY,
        "leaf-spring-assembly.toml"))

    @test Set(keys(loaded.bodies)) == Set((:leaf,
        Symbol("rear_spring.front_link"),
        Symbol("rear_spring.rear_link"),
        Symbol("rear_spring.shackle")))
    @test haskey(loaded.markers, Symbol("leaf.rear_spring.center"))
    @test haskey(loaded.markers,
        Symbol("ground.rear_spring.front_eye"))
    @test expanded["rear_spring"]["front_link"]["graphics"]["size"][2] ≈
        0.076 / 5
    @test expanded["leaf"]["graphics"]["rear_spring"]["central_leaf"]["shape"] ==
        "box"
    @test !haskey(expanded["model"], "assemblies")
    @test !occursin("type = \"leaf_spring\"", loaded.model_source)
    @test loaded.initial_values ≈ toml_loaded.initial_values
    @test loaded.active_variable_indices == toml_loaded.active_variable_indices
    @test loaded.active_equation_indices == toml_loaded.active_equation_indices

    vehicle_path = joinpath(ASSEMBLY_MODEL_DIRECTORY,
        "large-van-static.lua")
    vehicle = load_spatial_model(vehicle_path)
    vehicle_document = TOML.parse(vehicle.model_source)
    @test length(vehicle.bodies) == 29
    @test length(vehicle.connections) == 40
    @test length(vehicle.forces) == 64
    @test haskey(vehicle.bodies, Symbol("van.left_front.spindle"))
    @test haskey(vehicle.bodies, Symbol("van.rear.axle"))
    @test haskey(vehicle.forces, Symbol("van.left_front_shock"))
    @test haskey(vehicle.connections,
        Symbol("van.front_stabilizer.left_pivot"))
    @test haskey(vehicle.connections,
        Symbol("van.steering.pitman_to_center_perp"))
    @test haskey(vehicle.drivers,
        Symbol("van.steering.steering_wheel_motion"))
    @test haskey(vehicle.connections,
        Symbol("van.steering.steering_gear"))
    @test vehicle.connections[
        Symbol("van.steering.steering_gear")].coefficients == [-1.0, 19.5]
    @test haskey(vehicle.forces,
        Symbol("van.steering.windup_on_pitman"))
    @test haskey(vehicle.forces,
        Symbol("van.steering.windup_on_shaft"))
    @test haskey(vehicle.bodies,
        Symbol("van.steering.steering_wheel"))
    @test haskey(vehicle.bodies,
        Symbol("van.steering.windup_shaft"))
    bodywork = vehicle_document["van"]["body"]["graphics"]["bodywork"]
    @test bodywork["shape"] == "surface"
    @test length(bodywork["vertices"]) == 16
    steering_rim = vehicle_document["van"]["steering"]["steering_wheel"][
        "graphics"]["rim"]
    @test steering_rim["shape"] == "surface"
    @test steering_rim["marker"] ==
        "van.steering.steering_wheel.wheel_center"
    @test length(steering_rim["vertices"]) == 32 * 8
    road_pad = vehicle_document["ground"]["road"]["graphics"][
        "asphalt_pad"]
    @test road_pad["shape"] == "surface"
    @test road_pad["include_in_fit"] == false
    @test road_pad["vertices"] == [
        [-50.0, -50.0, 0.0], [50.0, -50.0, 0.0],
        [50.0, 50.0, 0.0], [-50.0, 50.0, 0.0]]
    @test haskey(vehicle.bodies,
        Symbol("van.left_front_tire.wheel"))
    @test haskey(vehicle.forces,
        Symbol("van.left_front_tire.contact"))
    default_tire_name = Symbol("van.left_front_tire.contact")
    @test !vehicle.forces[default_tire_name].transient
    dynamic_vehicle_source = read(joinpath(ASSEMBLY_MODEL_DIRECTORY,
        "large-van.lua"), String)
    dynamic_vehicle_document =
        PracticalMechanicalSimulation.AssemblyExpansion.parse_lua_model(
            dynamic_vehicle_source;
            element_types = PracticalMechanicalSimulation.SpatialModelIO.
                SPATIAL_ELEMENT_TYPES)
    default_tire_table = dynamic_vehicle_document["van"]["left_front_tire"]["contact"]
    default_wheel_table = dynamic_vehicle_document["van"]["left_front_tire"]["wheel"]
    @test default_wheel_table["velocity"][1] < 0
    @test default_wheel_table["angular_velocity"][3] < 0
    @test !haskey(default_tire_table, "longitudinal_relaxation_length")
    @test !haskey(default_tire_table, "lateral_relaxation_length")
    @test occursin(".slip_ratio", default_tire_table["longitudinal_expression"])
    @test occursin(".lateral_slip_velocity",
        default_tire_table["lateral_expression"])
    transient_vehicle_source = replace(dynamic_vehicle_source,
        "tire_damping_time_scale = 0.01," =>
        "tire_damping_time_scale = 0.01,\n" *
        "    tire_longitudinal_relaxation_length = 0.30,\n" *
        "    tire_lateral_relaxation_length = 0.45,")
    transient_vehicle_document =
        PracticalMechanicalSimulation.AssemblyExpansion.parse_lua_model(
            transient_vehicle_source;
            element_types = PracticalMechanicalSimulation.SpatialModelIO.
                SPATIAL_ELEMENT_TYPES)
    transient_tire_table = transient_vehicle_document["van"][
        "left_front_tire"]["contact"]
    @test transient_tire_table["longitudinal_relaxation_length"] == 0.30
    @test transient_tire_table["lateral_relaxation_length"] == 0.45
    @test occursin(".longitudinal_deformation",
        transient_tire_table["longitudinal_expression"])
    @test occursin(".lateral_deformation",
        transient_tire_table["lateral_expression"])
    bristle_vehicle_source = read(joinpath(ASSEMBLY_MODEL_DIRECTORY,
        "large-van-bristle-experiment.lua"), String)
    bristle_vehicle_document =
        PracticalMechanicalSimulation.AssemblyExpansion.parse_lua_model(
            bristle_vehicle_source;
            element_types = PracticalMechanicalSimulation.SpatialModelIO.
                SPATIAL_ELEMENT_TYPES)
    bristle_tire_table = bristle_vehicle_document["van"][
        "left_front_tire"]["contact"]
    @test bristle_tire_table["tangential_model"] == "bristle"
    @test bristle_tire_table["patch_length_by_load"] == [
        [0, 0], [2000, 0.15], [4000, 0.21], [6000, 0.26],
        [8000, 0.30], [10000, 0.33]]
    @test bristle_tire_table["cornering_stiffness_by_load"][4] ==
        [6000, 23000]
    @test bristle_tire_table["longitudinal_slip_stiffness_by_load"][4] ==
        [6000, 4300]
    @test bristle_tire_table["shear_release_time"] == 0.01
    @test bristle_tire_table["lateral_relaxation_fraction"] == 0.5
    @test !haskey(bristle_tire_table, "longitudinal_expression")
    @test !haskey(bristle_tire_table, "lateral_expression")
    left_front_bending = vehicle.forces[
        Symbol("van.rear.left.front_bending")]
    left_rear_bending = vehicle.forces[
        Symbol("van.rear.left.rear_bending")]
    @test vehicle.initial_values[
        left_front_bending.local_torque_variables[3]] ≈ 0.5 *
        (947.0*9.8-(0.08*1947.0)*9.8) * 0.75 *
        0.7318749021052059 / 2 rtol = 1.0e-10
    @test vehicle.initial_values[
        left_rear_bending.local_torque_variables[3]] ≈ -0.5 *
        (947.0*9.8-(0.08*1947.0)*9.8) * 0.75 *
        0.6809977036169096 / 2 rtol = 1.0e-10
    @test vehicle.forces[Symbol("van.left_front_jounce")].active_during ==
        (:dynamic, :modal)
    @test vehicle.forces[Symbol("van.left_front_rebound")].active_during ==
        (:dynamic, :modal)
    @test vehicle.forces[Symbol("van.left_rear_jounce")].active_during ==
        (:dynamic, :modal)
    plane_contacts = PracticalMechanicalSimulation.SpatialPlaneContacts
    plane_contact_type = plane_contacts.SpatialPlaneContactComponent
    for contact_name in ("van.left_front_jounce",
            "van.right_front_jounce", "van.left_front_rebound",
            "van.right_front_rebound", "van.left_rear_jounce",
            "van.right_rear_jounce")
        contact = vehicle.forces[Symbol(contact_name)]
        @test contact isa plane_contact_type
        @test contact.radius ≈ 0.03
        @test !contact.expression
        @test contact.stiffness == 4.0e6
        @test contact.damping_factor == 0.15
        @test contact.transition_depth == 0.001
        values = plane_contacts.spatial_plane_contact_values(
            contact, vehicle.initial_values)
        @test values.gap > 0
        @test values.normal_force == 0
    end
    @test vehicle.forces[Symbol("van.left_front_jounce")].sphere_marker.body.name ==
        Symbol("van.left_front.lower_control_arm")
    @test vehicle.forces[Symbol("van.left_front_jounce")].plane_marker.body.name ==
        Symbol("van.body")
    @test vehicle.forces[Symbol("van.left_front_rebound")].sphere_marker.body.name ==
        Symbol("van.left_front.upper_control_arm")
    @test vehicle.forces[Symbol("van.left_front_rebound")].plane_marker.body.name ==
        Symbol("van.body")
    @test vehicle.forces[Symbol("van.left_rear_jounce")].sphere_marker.body.name ==
        Symbol("van.body")
    @test vehicle.forces[Symbol("van.left_rear_jounce")].plane_marker.body.name ==
        Symbol("van.rear.axle")
    @test vehicle.forces[Symbol("van.right_rear_jounce")].sphere_marker.body.name ==
        Symbol("van.body")
    @test vehicle.forces[Symbol("van.right_rear_jounce")].plane_marker.body.name ==
        Symbol("van.rear.axle")
    for corner in ("front_left", "rear_left", "front_right", "rear_right")
        contact_name = "van.roof_$(corner)_bumper"
        contact = vehicle.forces[Symbol(contact_name)]
        @test contact isa plane_contact_type
        @test contact.active_during == (:dynamic, :modal)
        @test contact.radius ≈ 0.04
        @test contact.stiffness == 4.0e6
        @test contact.damping_factor == 0.05
        @test contact.transition_depth == 0.005
        @test contact.sphere_marker.body.name == Symbol("van.body")
        @test contact.plane_marker.name == Symbol("ground.road")
    end
    rear_mount = vehicle.forces[Symbol("van.rear.left.front_mount")]
    @test rear_mount.translational_stiffness == fill(1.0e6, 3)
    rear_shock = vehicle.forces[Symbol("van.left_rear_shock")]
    extension_state = copy(vehicle.initial_values)
    extension_state[rear_shock.span.velocity_variable] = 0.12
    @test rear_shock.law(0.0, extension_state) ≈ -420.0
    extension_state[rear_shock.span.velocity_variable] = -0.12
    @test rear_shock.law(0.0, extension_state) ≈ 144.0
    left_front_spring = vehicle.forces[Symbol("van.left_front_spring")]
    right_front_spring = vehicle.forces[Symbol("van.right_front_spring")]
    @test left_front_spring.stiffness ≈ 206719.93827462336
    @test right_front_spring.stiffness ≈ 196636.55316553017
    @test vehicle.initial_values[left_front_spring.force_variable] ≈
        9524.458934373224
    @test vehicle.initial_values[right_front_spring.force_variable] ≈
        9289.262954029411

    surface_model = load_spatial_model(joinpath(ASSEMBLY_MODEL_DIRECTORY,
        "simple-body-surface.lua"))
    surface_document = TOML.parse(surface_model.model_source)
    surface = surface_document["vehicle"]["graphics"]["bodywork"]
    @test surface["shape"] == "surface"
    @test length(surface["vertices"]) == 16
    @test Set(keys(surface["patches"])) == Set(("body", "glass"))
    @test length(surface["patches"]["body"]["faces"]) == 10
    @test length(surface["patches"]["glass"]["faces"]) == 4

    mktempdir() do directory
        marker_path = joinpath(directory, "located_marker.toml")
        wrapper_path = joinpath(directory, "marker_group.toml")
        model_path = joinpath(directory, "nested_model.toml")
        write(marker_path, """
            [assembly]
            name = "located_marker"
            dimension = "spatial"

            [interface.owner]
            kind = "body_or_ground"

            [parameters]
            point = { kind = "vector", required = true }

            [owner.tip]
            type = "marker"
            position = "=local_point(owner,point)"
            """)
        write(wrapper_path, """
            [assembly]
            name = "marker_group"
            dimension = "spatial"
            assemblies = ["located_marker.toml"]

            [interface.owner]
            kind = "body_or_ground"

            [parameters]
            point = { kind = "vector", required = true }

            [child]
            type = "located_marker"
            owner = "@owner"
            point = "=point"
            """)
        write(model_path, """
            [model]
            dimension = "spatial"
            assemblies = ["marker_group.toml"]

            [ground]
            type = "ground"

            [vehicle.body]
            type = "rigid_body"
            mass = 1.0
            inertia = [1.0, 1.0, 1.0]

            [vehicle.first]
            type = "marker_group"
            owner = "vehicle.body"
            point = [1.0, 2.0, 3.0]

            [vehicle.second]
            type = "marker_group"
            owner = "vehicle.body"
            point = [-1.0, -2.0, -3.0]
            """)
        nested = load_spatial_model(model_path)
        marker = nested.markers[Symbol("vehicle.body.first.child.tip")]
        @test marker.position_body == [1.0, 2.0, 3.0]
        second = nested.markers[Symbol("vehicle.body.second.child.tip")]
        @test second.position_body == [-1.0, -2.0, -3.0]

        optional_path = joinpath(directory, "optional_marker.lua")
        optional_model_path = joinpath(directory, "optional_model.toml")
        write(optional_path, """
            return {
                name = "optional_marker",
                dimension = "spatial",
                interfaces = {
                    owner = {kind = "body"}
                },
                parameters = {
                    point = {kind = "vector", required = true},
                    enabled = {kind = "boolean", default = true}
                },
                build = function(p)
                    local result = {}
                    if p.enabled then
                        result.owner = {
                            tip = {
                                type = "marker",
                                position = local_point(p.owner, p.point)
                            }
                        }
                    end
                    return result
                end
            }
            """)
        write(optional_model_path, """
            [model]
            dimension = "spatial"
            assemblies = ["optional_marker.lua"]

            [body]
            type = "rigid_body"
            mass = 1.0
            inertia = [1.0, 1.0, 1.0]

            [shown]
            type = "optional_marker"
            owner = "body"
            point = [1.0, 0.0, 0.0]

            [hidden]
            type = "optional_marker"
            owner = "body"
            point = [-1.0, 0.0, 0.0]
            enabled = false
            """)
        optional = load_spatial_model(optional_model_path)
        @test haskey(optional.markers, Symbol("body.shown.tip"))
        @test !haskey(optional.markers, Symbol("body.hidden.tip"))

        unknown_statement_path = joinpath(directory, "unknown_statement.lua")
        write(unknown_statement_path, """
            model {
                name = "unknown_statement",
                dimension = "spatial"
            }

            optional_marker {
                name = "tip",
                owner = "body",
                point = {1.0, 0.0, 0.0}
            }
            """)
        unknown_statement_error = try
            load_spatial_model(unknown_statement_path)
            nothing
        catch error
            error
        end
        @test unknown_statement_error isa ArgumentError
        @test occursin("unknown model statement 'optional_marker'",
            sprint(showerror, unknown_statement_error))

        recursive_path = joinpath(directory, "recursive.toml")
        recursive_model_path = joinpath(directory, "recursive_model.toml")
        write(recursive_path, """
            [assembly]
            name = "recursive"
            dimension = "spatial"

            [again]
            type = "recursive"
            """)
        write(recursive_model_path, """
            [model]
            dimension = "spatial"
            assemblies = ["recursive.toml"]

            [item]
            type = "recursive"
            """)
        @test_throws ArgumentError load_spatial_model(recursive_model_path)
    end
end
