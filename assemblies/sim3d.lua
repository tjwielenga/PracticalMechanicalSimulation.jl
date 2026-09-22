---@diagnostic disable: undefined-global

-- Public Lua modeling interface supplied by Sim3D. The host creates these
-- functions before a model is read; this module gives them an explicit Lua
-- namespace that editors and assembly modules can understand.
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

local sim3d = {
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
    gravity = gravity,
    applied_force = applied_force,
    applied_torque = applied_torque,
    spanning_force = spanning_force,
    spherical = spherical,
    perp = perp,
    inplane = inplane,
    inline = inline,
    hinge = hinge,
    orient = orient,
    revolute = revolute,
    fixed = fixed,
    rotational_motion = rotational_motion,
    translational_motion = translational_motion,
    spanning_motion = spanning_motion,
    span = span,
    directed_distance = directed_distance,
    bushing = bushing,
    curve = curve,
    curve_contact = curve_contact,
    flat_follower_contact = flat_follower_contact,
    plane_contact = plane_contact,
    rolling_tire = rolling_tire,
    coupler = coupler,
    gear_pair = gear_pair,
    rack_and_pinion = rack_and_pinion,
    pulley = pulley,
    belt = belt,
    belt_span = belt_span,

    model_table = model_table,
    required = required,
    positive = positive,
    nonnegative = nonnegative,
    vector = vector,
    dot = dot,
    cross = cross,
    norm = norm,
    unit = unit,
    frame = frame,
    link_frame = link_frame,
    rotation = rotation,
    compose = compose,
    box_inertia = box_inertia,
    local_point = local_point,
    local_frame = local_frame
}

return sim3d
