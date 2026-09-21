-- Goodyear LT245/75R16 wheel and tire converted from gylt24575r16_um.py
-- and its available analytical ADAMS tire property files.
-- Geometry inputs are initial global points and all units are SI.

local sim3d = require "sim3d"
local rigid_body, marker, revolute, bushing, rolling_tire =
    sim3d.rigid_body, sim3d.marker, sim3d.revolute, sim3d.bushing,
    sim3d.rolling_tire
local vector, cross, norm, unit, frame =
    sim3d.vector, sim3d.cross, sim3d.norm, sim3d.unit, sim3d.frame
local local_point, local_frame = sim3d.local_point, sim3d.local_frame
local required = sim3d.required
local context = "gylt_245_75_r16"

local function frame_with_z(z_direction, x_hint)
    local z = unit(z_direction)
    local x = x_hint-z*(x_hint[1]*z[1]+x_hint[2]*z[2]+x_hint[3]*z[3])
    if norm(x) <= 1.0e-12 then
        error("gylt_245_75_r16 forward direction is parallel to its axle")
    end
    x = unit(x)
    return frame(x, cross(z, x), z)
end

local function gylt_245_75_r16(p)
    local name = required(p, "name", context)
    local spindle = required(p, "spindle", context)
    local center = vector(required(p, "center", context))
    local road_marker = required(p, "road_marker", context)
    local pressure = p.pressure_psi or 60
    if pressure ~= 60 and pressure ~= 88 then
        error("gylt_245_75_r16 pressure_psi must be 60 or 88")
    end

    local axle = vector(p.axle_direction or {0.0, 1.0, 0.0})
    local forward = vector(p.forward_direction or {1.0, 0.0, 0.0})
    local wheel_frame = frame_with_z(axle, forward)
    local width = p.width or 0.245
    local rim_radius = (p.rim_diameter or 0.406)/2
    local outside_radius = p.outside_radius or
        rim_radius+(p.aspect_ratio or 0.75)*width
    local rolling_radius = p.rolling_radius or
        (pressure == 60 and 14.98*0.0254 or 15.1*0.0254)
    local mass = p.mass or 39.04
    local normal_stiffness = p.normal_stiffness or
        (pressure == 60 and 2378*175.126835 or 3021*175.126835)
    local longitudinal_stiffness = p.longitudinal_stiffness or 4000.0
    local cornering_stiffness = p.cornering_stiffness or
        162*4.4482216152605*180/math.pi
    local camber_stiffness = p.camber_stiffness or
        13.5*4.4482216152605*180/math.pi
    local friction = p.friction or 1.0
    local normal_damping_time_scale = p.normal_damping_time_scale or 0.0
    local regularization_speed = p.regularization_speed or 0.1
    local tangential_model = p.tangential_model or "expression"
    if tangential_model ~= "expression" and tangential_model ~= "bristle" then
        error(name .. " tangential_model must be expression or bristle")
    end
    local longitudinal_relaxation_length = p.longitudinal_relaxation_length
    local lateral_relaxation_length = p.lateral_relaxation_length
    if (longitudinal_relaxation_length == nil) ~=
            (lateral_relaxation_length == nil) then
        error(name .. " requires both tire relaxation lengths or neither")
    end
    local transient = longitudinal_relaxation_length ~= nil
    if tangential_model == "bristle" and transient then
        error(name .. " bristle model cannot use expression relaxation lengths")
    end
    if transient then
        sim3d.positive(p, "longitudinal_relaxation_length", name)
        sim3d.positive(p, "lateral_relaxation_length", name)
    end
    local forward_speed = p.forward_speed or 0.0
    local static_spin_stiffness = p.static_spin_stiffness or 0.0
    local static_spin_damping_time_scale =
        p.static_spin_damping_time_scale or 0.10
    sim3d.nonnegative(
        {static_spin_stiffness = static_spin_stiffness},
        "static_spin_stiffness", name)
    sim3d.nonnegative(
        {static_spin_damping_time_scale = static_spin_damping_time_scale},
        "static_spin_damping_time_scale", name)
    local color = p.color or "gray20"
    local visible = p.visible
    if visible == nil then visible = true end
    if type(visible) ~= "boolean" then
        error("gylt_245_75_r16 visible must be true or false")
    end
    local show_default = p.show_default
    if show_default == nil then show_default = true end
    if type(show_default) ~= "boolean" then
        error("gylt_245_75_r16 show_default must be true or false")
    end

    local axial_inertia = 0.5*mass*(outside_radius^2+rim_radius^2)
    local transverse_inertia = mass*
        (3*(outside_radius^2+rim_radius^2)+width^2)/12
    local wheel = name .. ".wheel"
    local axle_joint = name .. ".axle"
    local contact = name .. ".contact"
    local wheel_definition = {
        name = wheel,
        mass = mass,
        inertia = {transverse_inertia, transverse_inertia, axial_inertia},
        center_of_mass = wheel .. ".center",
        position = center,
        orientation = wheel_frame,
        graphics = {
            shape = "cylinder", marker = wheel .. ".center", axis = "z",
            radius = outside_radius, length = width,
            color = color, opacity = p.opacity or 0.82,
            visible = visible, show_default = show_default
        }
    }
    if forward_speed ~= 0.0 then
        wheel_definition.velocity = forward*forward_speed
        wheel_definition.angular_velocity = {
            0.0, 0.0, -forward_speed/rolling_radius
        }
    end
    rigid_body(wheel_definition)
    marker {name = wheel .. ".center"}
    marker {name = wheel .. ".hub"}
    marker {
        name = spindle .. "." .. name .. ".hub",
        position = local_point(spindle, center),
        orientation = local_frame(spindle, wheel_frame)
    }
    revolute {
        name = axle_joint,
        markers = {wheel .. ".hub", spindle .. "." .. name .. ".hub"},
        rotation_coordinates = true
    }
    local static_spin_hold = nil
    if static_spin_stiffness > 0.0 then
        static_spin_hold = name .. ".static_spin_hold"
        bushing {
            name = static_spin_hold,
            markers = {
                wheel .. ".hub", spindle .. "." .. name .. ".hub"
            },
            translational_stiffness = {0.0, 0.0, 0.0},
            rotational_stiffness = {0.0, 0.0, static_spin_stiffness},
            damping_time_scale = static_spin_damping_time_scale,
            active_during = "static"
        }
    end
    local tire_definition = {
        name = contact,
        markers = {wheel .. ".center", road_marker},
        radius = rolling_radius,
        regularization_speed = regularization_speed,
        normal_stiffness = normal_stiffness,
        normal_damping_time_scale = normal_damping_time_scale,
        friction_limit = "ellipse",
        mu_longitudinal = friction,
        mu_lateral = friction
    }
    if tangential_model == "bristle" then
        tire_definition.tangential_model = "bristle"
        tire_definition.patch_length_by_load = required(
            p, "patch_length_by_load", name)
        tire_definition.cornering_stiffness_by_load = required(
            p, "cornering_stiffness_by_load", name)
        tire_definition.longitudinal_slip_stiffness_by_load = required(
            p, "longitudinal_slip_stiffness_by_load", name)
        tire_definition.longitudinal_relaxation_fraction =
            p.longitudinal_relaxation_fraction or 1.0
        tire_definition.lateral_relaxation_fraction =
            p.lateral_relaxation_fraction or 1.0
        tire_definition.shear_release_time = p.shear_release_time or 0.01
    elseif transient then
        tire_definition.longitudinal_relaxation_length =
            longitudinal_relaxation_length
        tire_definition.lateral_relaxation_length = lateral_relaxation_length
        tire_definition.longitudinal_expression = string.format(
            "%.17g*%s.longitudinal_deformation",
            longitudinal_stiffness/longitudinal_relaxation_length, contact)
        tire_definition.lateral_expression = string.format(
            "-%.17g*%s.lateral_deformation-%.17g*%s.camber_angle",
            cornering_stiffness/lateral_relaxation_length, contact,
            camber_stiffness, contact)
    else
        tire_definition.longitudinal_expression = string.format(
            "%.17g*%s.slip_ratio", longitudinal_stiffness, contact)
        tire_definition.lateral_expression = string.format(
            "-%.17g*%s.lateral_slip_velocity/sqrt(" ..
            "%s.forward_velocity^2+%.17g^2)-%.17g*%s.camber_angle",
            cornering_stiffness, contact, contact, regularization_speed,
            camber_stiffness, contact)
    end
    rolling_tire(tire_definition)

    return {
        body = wheel,
        axle = axle_joint,
        static_spin_hold = static_spin_hold,
        tire = contact,
        mass = mass,
        outside_radius = outside_radius,
        rolling_radius = rolling_radius,
        initial_spin = -forward_speed/rolling_radius,
        rolling_resistance_coefficient = 0.015,
        omitted_outputs = {"rolling resistance", "aligning torque",
            "overturning moment", "separate rollover Pogo contact"}
    }
end

return gylt_245_75_r16
