local sim2d = require "sim2d"
local model, analysis, simulation, graphics, state_selection =
    sim2d.model, sim2d.analysis, sim2d.simulation, sim2d.graphics,
    sim2d.state_selection
local ground, marker = sim2d.ground, sim2d.marker
local double_pendulum = require "planar.double_pendulum"

model {
    name = "lua_double_pendulum",
    title = "Reusable Sim2D Lua double pendulum",
    dimension = "planar"
}

analysis {mode = "automatic"}

simulation {
    start_time = 0.0,
    end_time = 4.0,
    output_samples = 241,
    relative_tolerance = 1.0e-7,
    absolute_tolerance = 1.0e-9,
    maximum_step = 0.01
}

graphics {background = "white", body_palette = "colorblind"}

ground {name = "ground"}
marker {name = "ground.pin", position = {0.0, 0.0}}

local pendulum = double_pendulum {
    name = "pendulum",
    ground_pin = "ground.pin",
    pivot = {0.0, 0.0},
    first_length = 1.0,
    second_length = 0.8,
    first_mass = 1.0,
    second_mass = 0.7
}

state_selection {
    method = "preferred",
    preferred_velocities = {
        pendulum.base_joint .. ".omega",
        pendulum.elbow_joint .. ".omega"
    },
    allow_fallback = false
}
