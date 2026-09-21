-- Reusable planar double pendulum built from two reusable link constructions.

local sim2d = require "sim2d"
local rigid_body, marker, revolute, gravity =
    sim2d.rigid_body, sim2d.marker, sim2d.revolute, sim2d.gravity
local required, positive, vector =
    sim2d.required, sim2d.positive, sim2d.vector
local context = "double_pendulum"

local function direction(angle)
    return vector({math.cos(angle), math.sin(angle)})
end

local function link(name, inner_position, angle, length, mass, width, color)
    local axis = direction(angle)
    local center = inner_position + 0.5*length*axis
    local inner = name .. ".inner"
    local outer = name .. ".outer"

    rigid_body {
        name = name,
        mass = mass,
        inertia = mass*(length*length + width*width)/12,
        position = center,
        orientation = angle,
        graphics = {
            show_default = false,
            color = color,
            member = {
                shape = "cylinder",
                markers = {inner, outer},
                radius = width/2,
                color = color
            }
        }
    }
    marker {name = inner, position = {-length/2, 0.0}}
    marker {name = outer, position = {length/2, 0.0}}

    return {
        body = name,
        inner = inner,
        outer = outer,
        outer_position = inner_position + length*axis
    }
end

local function double_pendulum(p)
    local name = required(p, "name", context)
    local ground_pin = required(p, "ground_pin", context)
    local pivot = vector(p.pivot or {0.0, 0.0})
    local first_length = positive(p, "first_length", context)
    local second_length = positive(p, "second_length", context)
    local first_mass = positive(p, "first_mass", context)
    local second_mass = positive(p, "second_mass", context)
    local first_angle = p.first_angle or 0.0
    local second_angle = p.second_angle or 0.0
    local width = p.width or 0.08
    local root = name .. "."

    local first = link(root .. "first_link", pivot, first_angle,
        first_length, first_mass, width, p.first_color or "steelblue")
    local second = link(root .. "second_link", first.outer_position,
        second_angle, second_length, second_mass, width,
        p.second_color or "darkorange")
    local base_joint = root .. "base_joint"
    local elbow_joint = root .. "elbow_joint"

    revolute {
        name = base_joint,
        markers = {first.inner, ground_pin},
        rotation_coordinates = true
    }
    revolute {
        name = elbow_joint,
        markers = {second.inner, first.outer},
        rotation_coordinates = true
    }
    gravity {
        name = root .. "gravity",
        acceleration = p.gravity or {0.0, -9.81},
        bodies = {first.body, second.body}
    }

    return {
        bodies = {first.body, second.body},
        base_joint = base_joint,
        elbow_joint = elbow_joint,
        first = first,
        second = second
    }
end

return double_pendulum
