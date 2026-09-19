local sim3d = require "sim3d"
local model, analysis, simulation, graphics, ground, marker, gravity,
    bushing =
    sim3d.model, sim3d.analysis, sim3d.simulation, sim3d.graphics,
    sim3d.ground, sim3d.marker, sim3d.gravity,
    sim3d.bushing
local vector, cross, unit, frame =
    sim3d.vector, sim3d.cross, sim3d.unit, sim3d.frame
local large_van = require "spatial.large_van"
local ground_pad = require "spatial.ground_pad"

model {
    name = "large_van_static",
    title = "Large Van static assembly",
    dimension = "spatial"
}

analysis {
    mode = "static",
    static_method = "newton"
}

simulation {
    start_time = 0.0,
    end_time = 0.0,
    output_samples = 1,
    relative_tolerance = 1.0e-5,
    absolute_tolerance = 1.0e-5,
    initial_step = 1.0e-7,
    maximum_step = 0.001
}

graphics {
    background = "white",
    body_palette = "colorblind"
}

ground {name = "ground"}
-- The measured wheel centers lie within a few millimeters of this plane.
-- Using it as the nominal road makes the four starting tire loads close to
-- their measured front and rear axle loads instead of multiplying those
-- small survey-height differences by the tire stiffness.
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
marker {name = "ground.static_guide"}

local van = large_van {
    name = "van",
    road_marker = "ground.road",
    -- Hold the steering wheel straight while the complete steering column,
    -- gear coupler, and windup compliance settle.
    steering_wheel_motion_angle = 0.0,
    tire_damping_time_scale = 0.01,
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

-- A weak static-only bushing removes the neutral road-plane translations and
-- heading. It has no heave, roll, or pitch stiffness, so those coordinates
-- remain determined by suspension and tire force balance. The dynamic model
-- does not contain this temporary element.
marker {name = van.body .. ".static_guide"}
bushing {
    name = "static_guide",
    markers = {van.body .. ".static_guide", "ground.static_guide"},
    translational_stiffness = {5000.0, 5000.0, 0.0},
    rotational_stiffness = {0.0, 0.0, 5000.0},
    damping_time_scale = 0.10,
    active_during = "static"
}

gravity {
    name = "gravity",
    acceleration = {0.0, 0.0, -9.81},
    bodies = van.bodies
}
