local sim3d = require "sim3d"
local model, analysis, simulation, graphics, initial_conditions,
    ground, marker, gravity =
    sim3d.model, sim3d.analysis, sim3d.simulation, sim3d.graphics,
    sim3d.initial_conditions, sim3d.ground, sim3d.marker, sim3d.gravity
local vector, cross, unit, frame =
    sim3d.vector, sim3d.cross, sim3d.unit, sim3d.frame
local large_van = require "spatial.large_van"
local ground_pad = require "spatial.ground_pad"

model {
    name = "large_van_modal",
    title = "Modes of the Large Van suspension (transient tires)",
    dimension = "spatial"
}

analysis {
    mode = "modal",
    modes = 9,
    -- Avoid factoring directly at the free vehicle's road-plane rigid modes.
    frequency_shift_hz = 2.0
}

-- Transfer the settled configuration, but build the modal system without the
-- temporary static guide or any other static-only applied force.
initial_conditions {
    result = "../../results/examples/spatial/large-van-static.simp",
    sample = "last",
    include_velocities = false
}

simulation {
    start_time = 0.0,
    relative_tolerance = 1.0e-5,
    absolute_tolerance = 1.0e-7
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
    length = 100.0, width = 100.0, color = "gray25"
}

local van = large_van {
    name = "van",
    road_marker = "ground.road",
    -- Hold the steering input at the static value. No forward velocity or
    -- wheel spin is imposed at the stationary modal operating point.
    steering_wheel_motion_angle = 0.0,
    tire_damping_time_scale = 0.01,
    -- This separate linearization example retains the transient tire states.
    tire_longitudinal_relaxation_length = 0.30,
    tire_lateral_relaxation_length = 0.45,
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
