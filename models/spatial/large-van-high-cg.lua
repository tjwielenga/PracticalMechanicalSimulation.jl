-- The same van and maneuver as large-van.lua, with only the sprung-body
-- center of mass raised 4 cm in its body-reference frame.
local sim3d = require "sim3d"
local model, analysis, simulation, graphics, state_selection,
    initial_conditions, ground, marker, gravity =
    sim3d.model, sim3d.analysis, sim3d.simulation, sim3d.graphics,
    sim3d.state_selection, sim3d.initial_conditions,
    sim3d.ground, sim3d.marker, sim3d.gravity
local vector, cross, unit, frame =
    sim3d.vector, sim3d.cross, sim3d.unit, sim3d.frame
local large_van = require "spatial.large_van"
local ground_pad = require "spatial.ground_pad"

model {
    name = "large_van_high_cg",
    title = "Large Van - High CG",
    dimension = "spatial"
}

analysis {
    mode = "dynamic",
    modes = 9,
    frequency_shift_hz = 2.0
}

initial_conditions {
    result = "../../results/examples/spatial/large-van-static.simp",
    sample = "last",
    include_velocities = false
}

simulation {
    start_time = 0.0,
    end_time = 6.0,
    output_samples = 361,
    relative_tolerance = 1.0e-5,
    absolute_tolerance = 1.0e-7,
    initial_step = 1.0e-7,
    maximum_step = 0.1
}

state_selection {
    method = "preferred",
    preferred_velocities = {
        "van.left_front_tire.axle.omega",
        "van.right_front_tire.axle.omega",
        "van.left_rear_tire.axle.omega",
        "van.right_rear_tire.axle.omega",
        "van.left_front.upper_hinge.omega",
        "van.right_front.upper_hinge.omega",
        "van.rear.axle.V_z",
        "van.rear.axle.omega_x",
        "van.body.V_x", "van.body.V_y", "van.body.V_z",
        "van.body.omega_x", "van.body.omega_y", "van.body.omega_z"
    },
    allow_fallback = false
}

graphics {
    background = "white",
    body_palette = "colorblind"
}

ground {name = "ground"}
local road_z = unit(vector {-0.0010453444, -0.0020088528, 1.0})
local road_x = unit(cross(vector {0.0, 1.0, 0.0}, road_z))
local road_y = cross(road_z, road_x)
marker {
    name = "ground.road",
    position = {0.0, 0.0, -0.0295190411},
    orientation = frame(road_x, road_y, road_z),
    graphics = {
        shape = "xy_frame", axis_length = 0.5, plane_size = 0.5,
        opacity = 0.0,
        label = "road"
    }
}
ground_pad {
    name = "asphalt_pad", marker = "ground.road",
    length = 360.0, width = 360.0, color = "gray25"
}

local u = "min(max(t,0),1)"
local van = large_van {
    name = "van",
    road_marker = "ground.road",
    steering_wheel_motion_angle =
        "-pi*(10*" .. u .. "^3-15*" .. u .. "^4+6*" .. u .. "^5)",
    tire_damping_time_scale = 0.01,
    forward_speed = 30.0,
    sprung_center = {2.102, -0.020, 0.767},
    roof_contacts = true,
    colors = {
        body = "gray70",
        control_arm = "steelblue",
        spindle = "darkorange",
        stabilizer = "firebrick",
        stabilizer_link = "orange",
        steering_link = "gray35",
        tie_rod = "seagreen",
        steering_wheel = "black",
        axle = "gray35",
        leaf = "steelblue",
        shackle = "darkorange",
        tire = "gray20"
    }
}

gravity {
    name = "gravity",
    acceleration = {0.0, 0.0, -9.81},
    bodies = van.bodies
}
