local sim3d = require "sim3d"
local model, analysis, simulation, graphics, ground, marker, gravity =
    sim3d.model, sim3d.analysis, sim3d.simulation, sim3d.graphics,
    sim3d.ground, sim3d.marker, sim3d.gravity
local vector, cross, unit, frame =
    sim3d.vector, sim3d.cross, sim3d.unit, sim3d.frame
local large_van = require "spatial.large_van"
local ground_pad = require "spatial.ground_pad"

model {
    name = "large_van_full_preview",
    title = "Complete Large Van — entered model",
    dimension = "spatial"
}

-- This entry exists for inspecting the complete assembly before a successful
-- static result is available. The development static and dynamic entries are
-- used for actual analysis.
analysis {mode = "dynamic"}

simulation {
    start_time = 0.0,
    end_time = 0.1,
    output_samples = 7,
    relative_tolerance = 1.0e-5,
    absolute_tolerance = 1.0e-7,
    initial_step = 1.0e-7,
    maximum_step = 0.001
}

graphics {background = "white", body_palette = "colorblind"}

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
    tire_damping_time_scale = 0.01,
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
