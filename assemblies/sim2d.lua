---@diagnostic disable: undefined-global

-- Public Lua modeling interface supplied by Sim2D. The host creates these
-- functions before a model is read; this module gives them an explicit Lua
-- namespace that editors and reusable assembly modules can understand.
local function validation_context(context)
    if type(context) ~= "string" or context == "" then
        error("validation context must be a nonempty string")
    end
    return context
end

local function required(parameters, field, context)
    context = validation_context(context)
    if type(parameters) ~= "table" then
        error(context .. " parameters must be a table")
    end
    if type(field) ~= "string" or field == "" then
        error(context .. " required field name must be a nonempty string")
    end
    local value = parameters[field]
    if value == nil then
        error(context .. " requires '" .. field .. "'")
    end
    return value
end

local function positive(parameters, field, context)
    local value = required(parameters, field, context)
    if type(value) ~= "number" or value <= 0 then
        error(context .. " '" .. field .. "' must be positive")
    end
    return value
end

local function nonnegative(parameters, field, context)
    local value = required(parameters, field, context)
    if type(value) ~= "number" or value < 0 then
        error(context .. " '" .. field .. "' must be nonnegative")
    end
    return value
end

local function local_point(body, point)
    return {
        __simp_deferred = "planar_local_point",
        body = body,
        value = point
    }
end

local function local_orientation(body, angle)
    return {
        __simp_deferred = "planar_local_orientation",
        body = body,
        value = angle
    }
end

local sim2d = {
    model = model,
    analysis = analysis,
    simulation = simulation,
    parameters = parameters,
    graphics = graphics,
    state_selection = state_selection,
    initial_conditions = initial_conditions,

    ground = ground,
    rigid_body = rigid_body,
    marker = marker,
    floating_marker = floating_marker,
    revolute = revolute,
    inplane = inplane,
    perp = perp,
    translational = translational,
    fixed = fixed,
    gear_pair = gear_pair,
    rack_and_pinion = rack_and_pinion,
    distance_coordinate = distance_coordinate,
    coupler = coupler,
    span = span,
    pulley = pulley,
    belt = belt,
    belt_span = belt_span,
    rotational_motion = rotational_motion,
    translational_motion = translational_motion,
    gravity = gravity,
    applied_force = applied_force,
    applied_torque = applied_torque,
    bushing = bushing,
    curve = curve,
    plane_contact = plane_contact,
    curve_contact = curve_contact,
    flat_follower_contact = flat_follower_contact,
    spanning_force = spanning_force,
    torsional_spring_damper = torsional_spring_damper,
    surface_friction = surface_friction,
    revolute_friction = revolute_friction,
    translational_friction = translational_friction,
    inplane_friction = inplane_friction,
    equation_component = equation_component,

    model_table = model_table,
    required = required,
    positive = positive,
    nonnegative = nonnegative,
    vector = vector,
    dot = dot,
    norm = norm,
    unit = unit,
    local_point = local_point,
    local_orientation = local_orientation
}

return sim2d
