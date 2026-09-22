-- Small vehicle force constructions used by converted historical assemblies.
-- Geometry inputs are initial global points and all units are SI.

local sim3d = require "sim3d"
local marker, spanning_force, plane_contact =
    sim3d.marker, sim3d.spanning_force, sim3d.plane_contact
local vector, dot, cross, norm, unit, frame =
    sim3d.vector, sim3d.dot, sim3d.cross, sim3d.norm,
    sim3d.unit, sim3d.frame
local local_point, local_frame = sim3d.local_point, sim3d.local_frame
local required = sim3d.required

local function viewer_graphics(p)
    if p.assembly == nil then
        return p.graphics
    end
    local result = {}
    for key, value in pairs(p.graphics or {}) do
        result[key] = value
    end
    result.assembly = p.assembly
    return result
end

local function external_marker(body, name, point, orientation)
    local specification = {
        name = body .. "." .. name,
        position = local_point(body, vector(point))
    }
    if orientation ~= nil then
        specification.orientation = local_frame(body, orientation)
    end
    marker(specification)
    return body .. "." .. name
end

local function frame_with_z(z_direction)
    local z = unit(z_direction)
    local hint = vector {1.0, 0.0, 0.0}
    local x = hint-dot(hint, z)*z
    if norm(x) <= 1.0e-12 then
        hint = vector {0.0, 1.0, 0.0}
        x = hint-dot(hint, z)*z
    end
    x = unit(x)
    return frame(x, cross(z, x), z)
end

local function spring(p)
    local name = required(p, "name", "vehicle spring")
    local first_body = required(p, "first_body", "vehicle spring")
    local second_body = required(p, "second_body", "vehicle spring")
    local first = external_marker(first_body, name .. ".first",
        required(p, "first_point", "vehicle spring"))
    local second = external_marker(second_body, name .. ".second",
        required(p, "second_point", "vehicle spring"))
    spanning_force {
        name = name,
        markers = {first, second},
        stiffness = required(p, "stiffness", "vehicle spring"),
        damping = p.damping or 0.0,
        free_length = required(p, "free_length", "vehicle spring"),
        graphics = viewer_graphics(p)
    }
end

local function number_text(value)
    return string.format("%.17g", value)
end

-- Return a restricted Sim3D expression for a linearly extrapolated table.
-- The hinge representation is exact at every supplied knot.  Slope changes
-- are introduced with max(x-x_k,0), which the force-expression reader and
-- its automatic differentiation already support.
local function piecewise_linear_expression(argument, points)
    if type(points) ~= "table" or #points < 2 then
        error("tabulated damper requires at least two [velocity, force] pairs")
    end
    local x0, y0 = points[1][1], points[1][2]
    local previous_slope =
        (points[2][2]-points[1][2])/(points[2][1]-points[1][1])
    local terms = {
        number_text(y0),
        number_text(previous_slope) .. "*((" .. argument .. ")-(" ..
            number_text(x0) .. "))"
    }
    for index = 2, #points-1 do
        local left, right = points[index], points[index+1]
        if right[1] <= left[1] then
            error("tabulated damper velocities must increase strictly")
        end
        local slope = (right[2]-left[2])/(right[1]-left[1])
        local slope_change = slope-previous_slope
        if math.abs(slope_change) > 1.0e-14 then
            table.insert(terms, number_text(slope_change) .. "*max((" ..
                argument .. ")-(" .. number_text(left[1]) .. "),0)")
        end
        previous_slope = slope
    end
    return table.concat(terms, "+")
end

local function tabulated_damper(p)
    local name = required(p, "name", "tabulated damper")
    local first_body = required(p, "first_body", "tabulated damper")
    local second_body = required(p, "second_body", "tabulated damper")
    local first = external_marker(first_body, name .. ".first",
        required(p, "first_point", "tabulated damper"))
    local second = external_marker(second_body, name .. ".second",
        required(p, "second_point", "tabulated damper"))
    local argument = p.reverse_velocity == false and
        name .. ".length_rate" or "-(" .. name .. ".length_rate)"
    -- The table returns the ADAMS force on the first marker. Preserve its
    -- historical compression-velocity argument and use its force sign
    -- directly.
    local expression = piecewise_linear_expression(argument,
        required(p, "table", "tabulated damper"))
    spanning_force {
        name = name,
        markers = {first, second},
        expression = expression,
        graphics = viewer_graphics(p)
    }
end

local function bumper(p)
    local name = required(p, "name", "vehicle bumper")
    local lower_body = required(p, "lower_body", "vehicle bumper")
    local upper_body = required(p, "upper_body", "vehicle bumper")
    local lower_point = vector(required(p, "lower_point", "vehicle bumper"))
    local upper_point = vector(required(p, "upper_point", "vehicle bumper"))
    local kind = required(p, "kind", "vehicle bumper")
    if kind ~= "jounce" and kind ~= "rebound" then
        error("vehicle bumper kind must be 'jounce' or 'rebound'")
    end
    local radius = p.radius or 0.03
    local stiffness = required(p, "stiffness", "vehicle bumper")
    local damping_factor = p.damping_factor or 0.15
    local transition_depth = p.transition_depth or 0.001
    if radius <= 0 then
        error("vehicle bumper radius must be positive")
    end
    if stiffness < 0 then
        error("vehicle bumper stiffness must be nonnegative")
    end
    if damping_factor < 0 then
        error("vehicle bumper damping_factor must be nonnegative")
    end
    if transition_depth < 0 then
        error("vehicle bumper transition_depth must be nonnegative")
    end
    if norm(upper_point-lower_point) <= 1.0e-12 then
        error("vehicle bumper points must be separated")
    end

    local sphere_on = p.sphere_on or (kind == "jounce" and "lower" or "upper")
    if sphere_on ~= "lower" and sphere_on ~= "upper" then
        error("vehicle bumper sphere_on must be 'lower' or 'upper'")
    end
    local sphere_body, sphere_point, sphere_role
    local plane_body, plane_point, plane_role
    if sphere_on == "lower" then
        sphere_body, sphere_point, sphere_role = lower_body, lower_point, "lower"
        plane_body, plane_point, plane_role = upper_body, upper_point, "upper"
    else
        sphere_body, sphere_point, sphere_role = upper_body, upper_point, "upper"
        plane_body, plane_point, plane_role = lower_body, lower_point, "lower"
    end
    -- The plane normal follows the former span line but points from the plane
    -- toward the sphere, so positive contact force separates the two bodies.
    local plane_orientation = frame_with_z(sphere_point-plane_point)
    local sphere = external_marker(sphere_body,
        name .. "." .. sphere_role, sphere_point)
    local plane = external_marker(plane_body,
        name .. "." .. plane_role, plane_point, plane_orientation)
    local active_during = p.active_during
    local inactive_during = p.inactive_during
    if active_during == nil and inactive_during == nil then
        inactive_during = "static"
    end
    plane_contact {
        name = name,
        markers = {sphere, plane},
        radius = radius,
        stiffness = stiffness,
        damping_factor = damping_factor,
        transition_depth = transition_depth,
        active_during = active_during,
        inactive_during = inactive_during,
        graphics = viewer_graphics(p)
    }
end

return {
    spring = spring,
    tabulated_damper = tabulated_damper,
    bumper = bumper,
    piecewise_linear_expression = piecewise_linear_expression
}
