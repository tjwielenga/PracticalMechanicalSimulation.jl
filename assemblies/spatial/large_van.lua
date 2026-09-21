-- Large-van chassis and suspension adapted from the author's earlier model.
-- Geometry inputs are initial global points and all units are SI.

local sim3d = require "sim3d"
local rigid_body, marker, spanning_motion, model_table, plane_contact =
    sim3d.rigid_body, sim3d.marker, sim3d.spanning_motion,
    sim3d.model_table, sim3d.plane_contact
local vector, dot, cross, norm, unit, link_frame =
    sim3d.vector, sim3d.dot, sim3d.cross, sim3d.norm, sim3d.unit,
    sim3d.link_frame
local sla_suspension = require "spatial.sla_suspension"
local stabilizer_bar = require "spatial.stabilizer_bar"
local steering_linkage = require "spatial.steering_linkage"
local hotchkiss_suspension = require "spatial.hotchkiss_suspension"
local gylt_245_75_r16 = require "spatial.gylt_245_75_r16"
local vehicle_forces = require "spatial.vehicle_force_elements"
local simple_body = require "spatial.simple_body"

local function large_van(p)
    p = p or {}
    local name = p.name or "large_van"
    local root = name .. "."
    local body = root .. "body"
    local up = vector {0.0, 0.0, 1.0}
    -- The historical van hardpoints place the front axle at the smaller x
    -- coordinate, so the vehicle's natural forward direction is -x.
    local forward = vector(p.forward_direction or {-1.0, 0.0, 0.0})
    local forward_speed = p.forward_speed or 0.0
    local initial_velocity = forward*forward_speed

    local colors = p.colors or {}
    local steel_density = p.steel_density or 7850.0
    local tire_mass = p.tire_mass or 39.04
    local wheel_mass = p.wheel_mass or (p.road_marker and tire_mass or 10.0)
    local curb_mass = p.curb_mass or 1947.0
    -- These body properties were recorded with the earlier van model; its
    -- original script still contained a 100 kg placeholder.
    local sprung_mass = p.sprung_mass or 1649.1
    local sprung_inertia = p.sprung_inertia or {607.0, 2990.0, 3089.0}
    local sprung_center = p.sprung_center or {2.102, -0.020, 0.727}

    rigid_body {
        name = body,
        mass = sprung_mass,
        inertia = sprung_inertia,
        center_of_mass = body .. ".cm",
        graphics = {show_default = false, color = colors.body or "gray60"}
    }
    marker {name = body .. ".cm", position = sprung_center}
    local body_graphic = nil
    if p.body_graphics ~= false then
        local dimensions = type(p.body_graphics) == "table" and
            p.body_graphics or {}
        body_graphic = simple_body {
            name = "bodywork",
            body = body,
            origin = dimensions.origin or {0.0, 0.0, 0.08},
            hood_length = dimensions.hood_length or 1.05,
            roof_front = dimensions.roof_front or 1.45,
            roof_back = dimensions.roof_back or 4.35,
            window_back = dimensions.window_back or 4.70,
            length = dimensions.length or 5.20,
            front_hood_height = dimensions.front_hood_height or 0.78,
            windshield_height = dimensions.windshield_height or 1.08,
            roof_height = dimensions.roof_height or 1.86,
            back_glass_height = dimensions.back_glass_height or 1.16,
            trunk_height = dimensions.trunk_height or 0.78,
            rocker_width = dimensions.rocker_width or 1.45,
            body_width = dimensions.body_width or 1.84,
            roof_width = dimensions.roof_width or 1.64,
            body_color = colors.body or "gray60",
            body_opacity = p.body_opacity or 0.72,
            glass_color = colors.glass or "gray20",
            glass_opacity = p.glass_opacity or 0.38,
            draw_edges = true,
            edge_color = colors.body_edge or "gray35",
            edge_width = 1.0
        }
    end

    -- Optional roof-to-road contacts share the exact roof corners used by
    -- the body graphic. They are inactive during static assembly and become
    -- physical sphere/plane bumpers if the vehicle rolls over.
    if p.roof_contacts then
        if not p.road_marker or not body_graphic then
            error("large_van roof_contacts require road_marker and body graphics")
        end
        local options = type(p.roof_contacts) == "table" and p.roof_contacts or {}
        for _, corner in ipairs({"front_left", "rear_left", "front_right", "rear_right"}) do
            local contact_name = root .. "roof_" .. corner .. "_bumper"
            local sphere_marker = body .. ".roof_" .. corner .. "_bumper"
            marker {name = sphere_marker, position = body_graphic.roof_corners[corner]}
            plane_contact {
                name = contact_name,
                markers = {sphere_marker, p.road_marker},
                radius = options.radius or 0.04,
                stiffness = options.stiffness or 4.0e6,
                damping_factor = options.damping_factor or 0.05,
                transition_depth = options.transition_depth or 0.005,
                inactive_during = "static"
            }
        end
    end

    -- Wheel centers are retained for future wheel and tire assemblies.
    local left_front_wheel = vector {0.759, -0.874, 0.338}
    local right_front_wheel = vector {0.771, 0.874, 0.342}
    local left_rear_wheel = vector {4.474, -0.866, 0.343}
    local right_rear_wheel = vector {4.475, 0.866, 0.346}

    -- Front SLA hardpoints.
    local left_lower_front = vector {0.571, -0.312, 0.293}
    local left_lower_rear = vector {0.945, -0.316, 0.305}
    local left_lower_ball = vector {0.757, -0.751, 0.236}
    local right_lower_front = vector {0.575, 0.318, 0.294}
    local right_lower_rear = vector {0.947, 0.316, 0.304}
    local right_lower_ball = vector {0.765, 0.754, 0.239}
    local left_front = sla_suspension {
        name = root .. "left_front",
        frame = body,
        spindle_mass = (14.6+0.045*curb_mass)/2-wheel_mass,
        upper_front = {0.707, -0.461, 0.571},
        upper_rear = {0.925, -0.480, 0.542},
        upper_ball = {0.769, -0.710, 0.534},
        lower_front = left_lower_front,
        lower_rear = left_lower_rear,
        lower_ball = left_lower_ball,
        damping_time_scale = 0.01,
        spindle_diameter = 0.03,
        arm_color = colors.control_arm,
        spindle_color = colors.spindle
    }
    local right_front = sla_suspension {
        name = root .. "right_front",
        frame = body,
        spindle_mass = (14.6+0.045*curb_mass)/2-wheel_mass,
        upper_front = {0.708, 0.462, 0.571},
        upper_rear = {0.925, 0.474, 0.542},
        upper_ball = {0.777, 0.710, 0.530},
        lower_front = right_lower_front,
        lower_rear = right_lower_rear,
        lower_ball = right_lower_ball,
        damping_time_scale = 0.01,
        spindle_diameter = 0.03,
        arm_color = colors.control_arm,
        spindle_color = colors.spindle
    }

    -- Front coil springs.  The free lengths reproduce the curb-load
    -- calculation in the Python model.
    local left_spring_bottom = vector {0.789, -0.527, 0.211}
    local left_spring_top = vector {0.821, -0.501, 0.502}
    local right_spring_bottom = vector {0.788, 0.532, 0.212}
    local right_spring_top = vector {0.813, 0.511, 0.500}
    local front_effective_wheel_stiffness =
        p.front_spring_stiffness or 44100.0
    local front_tire_force = p.front_tire_force or 1000.0*9.8
    local spindle_mass = (14.6+0.045*curb_mass)/2-wheel_mass
    local front_sprung_wheel_force =
        front_tire_force/2-(spindle_mass+wheel_mass)*9.8
    local function coil_parameters(lower_front, lower_rear, lower_ball,
            spring_bottom, spring_top)
        local hinge_axis = unit(lower_front-lower_rear)
        local spring_axis = unit(spring_bottom-spring_top)
        local wheel_lever = math.abs(dot(hinge_axis,
            cross(lower_ball-lower_front, up)))
        local spring_lever = math.abs(dot(hinge_axis,
            cross(spring_bottom-lower_front, spring_axis)))
        if wheel_lever <= 1.0e-12 or spring_lever <= 1.0e-12 then
            error("front spring hardpoints give a zero suspension lever arm")
        end
        local motion_ratio = spring_lever/wheel_lever
        -- k_wheel = k_coil*(ds/dz_wheel)^2. The initial coil compression is
        -- the required coil force divided by that converted coil rate.
        local stiffness = front_effective_wheel_stiffness/motion_ratio^2
        local deflection = front_sprung_wheel_force*motion_ratio/
            front_effective_wheel_stiffness
        return stiffness, deflection
    end
    local left_spring_stiffness, left_spring_deflection = coil_parameters(
        left_lower_front, left_lower_rear, left_lower_ball,
        left_spring_bottom, left_spring_top)
    local right_spring_stiffness, right_spring_deflection = coil_parameters(
        right_lower_front, right_lower_rear, right_lower_ball,
        right_spring_bottom, right_spring_top)
    vehicle_forces.spring {
        name = root .. "left_front_spring",
        first_body = body, first_point = left_spring_top,
        second_body = left_front.lower_control_arm,
        second_point = left_spring_bottom,
        free_length = norm(left_spring_bottom-left_spring_top)+
            left_spring_deflection,
        stiffness = left_spring_stiffness
    }
    vehicle_forces.spring {
        name = root .. "right_front_spring",
        first_body = body, first_point = right_spring_top,
        second_body = right_front.lower_control_arm,
        second_point = right_spring_bottom,
        free_length = norm(right_spring_bottom-right_spring_top)+
            right_spring_deflection,
        stiffness = right_spring_stiffness
    }

    local stabilizer = nil
    if p.include_front_stabilizer ~= false then
        stabilizer = stabilizer_bar {
            name = root .. "front_stabilizer",
            torsional_stiffness = p.stabilizer_stiffness or 3200.0,
            center_body = body,
            left_center_point = {0.305, -0.383, 0.388},
            right_center_point = {0.305, 0.397, 0.388},
            left_body = left_front.lower_control_arm,
            left_body_point = {0.623, -0.531, 0.388},
            left_offset_point = {0.623, -0.531, 0.288},
            right_body = right_front.lower_control_arm,
            right_body_point = {0.628, 0.536, 0.389},
            right_offset_point = {0.628, 0.536, 0.289},
            diameter = 0.027,
            damping_time_scale = p.stabilizer_damping_time_scale or 0.0,
            bar_color = colors.stabilizer,
            link_color = colors.stabilizer_link
        }
    end

    local steering = nil
    if p.include_steering ~= false then
        steering = steering_linkage {
            name = root .. "steering",
            frame = body,
            left_spindle = left_front.spindle,
            right_spindle = right_front.spindle,
            pitman_pivot = {0.499, -0.452, 0.545},
            pitman_axis_point = {0.341, -0.455, 0.541},
            center_to_pitman = {0.547, -0.449, 0.392},
            center_to_idler = {0.526, 0.487, 0.396},
            idler_pivot = {0.603, 0.485, 0.561},
            idler_axis_point = {0.496, 0.473, 0.553},
            left_inner_tie = {0.577, -0.348, 0.402},
            left_knuckle = {0.596, -0.756, 0.371},
            right_inner_tie = {0.591, 0.358, 0.399},
            right_knuckle = {0.607, 0.762, 0.361},
            steering_wheel = {1.727, -0.367, 1.133},
            steering_column_end = {0.390, -0.367, 0.478},
            diameter = 0.029,
            wheel_diameter = 0.30,
            free_play = p.steering_free_play or 0.127*2,
            rotational_damping = p.steering_rotational_damping or 100.0,
            pitman_motion_angle = p.steering_pitman_motion_angle,
            wheel_motion_angle = p.steering_wheel_motion_angle,
            -- Sim3D's revolute-coordinate orientation already accounts for
            -- the rear-facing pitman arm.  A positive ratio therefore gives
            -- the intended steering-wheel-to-road-wheel direction.
            gear_ratio = p.steering_gear_ratio or 19.5,
            windup_stiffness = p.steering_windup_stiffness or 27000.0,
            windup_damping = p.steering_windup_damping or 0.0,
            link_color = colors.steering_link,
            tie_color = colors.tie_rod,
            wheel_color = colors.steering_wheel
        }
    else
        local function fixed_toe_link(side, spindle, inner_point, outer_point)
            local inner = body .. "." .. name .. "." .. side .. "_toe_inner"
            local outer = spindle .. ".fixed_toe_outer"
            marker {name = inner, position = local_point(body, inner_point)}
            marker {name = outer, position = local_point(spindle, outer_point)}
            spanning_motion {
                name = root .. side .. "_fixed_toe",
                markers = {outer, inner},
                distance = norm(vector(outer_point)-vector(inner_point))
            }
        end
        fixed_toe_link("left", left_front.spindle,
            {0.577, -0.348, 0.402}, {0.596, -0.756, 0.371})
        fixed_toe_link("right", right_front.spindle,
            {0.591, 0.358, 0.399}, {0.607, 0.762, 0.361})
    end

    local front_shock_table = {
        {-1.92, -8256.0}, {-0.96, -4128.0}, {-0.48, -2064.0},
        {-0.24, -1032.0}, {-0.12, -516.0}, {0.0, 0.0},
        {0.12, 228.0}, {0.24, 456.0}, {0.48, 912.0},
        {0.96, 1824.0}, {1.92, 3648.0}
    }
    vehicle_forces.tabulated_damper {
        name = root .. "left_front_shock",
        first_body = body, first_point = {1.040, -0.477, 0.530},
        second_body = left_front.lower_control_arm,
        second_point = {0.910, -0.556, 0.261},
        table = front_shock_table
    }
    vehicle_forces.tabulated_damper {
        name = root .. "right_front_shock",
        first_body = body, first_point = {1.038, 0.477, 0.532},
        second_body = right_front.lower_control_arm,
        second_point = {0.922, 0.564, 0.257},
        table = front_shock_table
    }

    local bumper_stiffness = p.bump_stop_stiffness or 4.0e6
    local bumper_radius = p.bump_stop_radius or 0.030
    local bumper_damping_factor = p.bump_stop_damping_factor or 0.15
    local bumper_transition_depth = p.bump_stop_transition_depth or 0.001
    if p.include_bump_stops ~= false then
        local front_bumper_height = 0.075
        local rebound_height = 0.030
        vehicle_forces.bumper {
        name = root .. "left_front_jounce",
        kind = "jounce",
        lower_body = left_front.lower_control_arm,
        lower_point = {0.699, -0.642, 0.377-front_bumper_height},
        upper_body = body, upper_point = {0.740, -0.628, 0.411},
        radius = bumper_radius, stiffness = bumper_stiffness,
        damping_factor = bumper_damping_factor,
        transition_depth = bumper_transition_depth
        }
        vehicle_forces.bumper {
        name = root .. "right_front_jounce",
        kind = "jounce",
        lower_body = right_front.lower_control_arm,
        lower_point = {0.755, 0.648, 0.381-front_bumper_height},
        upper_body = body, upper_point = {0.750, 0.635, 0.415},
        radius = bumper_radius, stiffness = bumper_stiffness,
        damping_factor = bumper_damping_factor,
        transition_depth = bumper_transition_depth
        }
        local function rebound_base(top, bottom)
        local top_point, bottom_point = vector(top), vector(bottom)
        return bottom_point-(top_point-bottom_point)/norm(top_point-bottom_point)*
            rebound_height
        end
        vehicle_forces.bumper {
        name = root .. "left_front_rebound",
        kind = "rebound",
        upper_body = left_front.upper_control_arm,
        upper_point = {0.805, -0.610, 0.491},
        lower_body = body,
        lower_point = rebound_base(
            {0.805, -0.610, 0.491}, {0.813, -0.604, 0.453}),
        radius = bumper_radius, stiffness = bumper_stiffness,
        damping_factor = bumper_damping_factor,
        transition_depth = bumper_transition_depth
        }
        vehicle_forces.bumper {
        name = root .. "right_front_rebound",
        kind = "rebound",
        upper_body = right_front.upper_control_arm,
        upper_point = {0.815, 0.606, 0.491},
        lower_body = body,
        lower_point = rebound_base(
            {0.815, 0.606, 0.491}, {0.824, 0.602, 0.460}),
        radius = bumper_radius, stiffness = bumper_stiffness,
        damping_factor = bumper_damping_factor,
        transition_depth = bumper_transition_depth
        }
    end

    -- Rear axle and Hotchkiss leaf-spring suspension.
    local rear_axle_mass = 0.08*curb_mass-2*wheel_mass
    local rear_tire_force = p.rear_tire_force or 947.0*9.8
    local rear_leaf_nominal_load = p.rear_leaf_nominal_load or
        (rear_tire_force-(rear_axle_mass+2*wheel_mass)*9.8)/2
    local axle_diameter = 0.076
    local axle_length = norm(right_rear_wheel-left_rear_wheel)
    local axle_radius = axle_diameter/2
    local rear_axle_inertia = {
        0.5*rear_axle_mass*axle_radius^2,
        rear_axle_mass*(3*axle_radius^2+axle_length^2)/12,
        rear_axle_mass*(3*axle_radius^2+axle_length^2)/12
    }
    local rear = hotchkiss_suspension {
        name = root .. "rear",
        frame = body,
        axle_mass = rear_axle_mass,
        axle_inertia = rear_axle_inertia,
        left_axle_end = left_rear_wheel,
        right_axle_end = right_rear_wheel,
        axle_diameter = axle_diameter,
        differential_diameter = 0.200,
        left_front_eye = {3.742, -0.628, 0.494},
        left_leaf_center = {4.474, -0.650, 0.506},
        left_rear_shackle_eye = {5.152, -0.673, 0.572},
        left_rear_frame_eye = {5.133, -0.663, 0.684},
        right_front_eye = {3.742, 0.637, 0.499},
        right_leaf_center = {4.474, 0.661, 0.512},
        right_rear_shackle_eye = {5.150, 0.685, 0.575},
        right_rear_frame_eye = {5.130, 0.678, 0.687},
        up = up,
        leaf_stiffness = p.rear_leaf_stiffness or 38600.0,
        twist_stiffness = p.rear_leaf_twist_stiffness or 500.0,
        leaf_width = 0.064,
        -- The historical value 1000 is an automotive bushing rate in N/mm.
        -- Sim3D uses metres, so the translational rate is 1e6 N/m.
        mount_stiffness = p.rear_mount_stiffness or
            {1.0e6, 1.0e6, 1.0e6},
        mount_rotational_stiffness = {300.0, 300.0, 300.0},
        damping_time_scale = 0.01,
        nominal_load = rear_leaf_nominal_load,
        density = steel_density,
        axle_color = colors.axle,
        leaf_color = colors.leaf,
        shackle_color = colors.shackle
    }

    local rear_shock_table = {
        {-1.92, -6720.0}, {-0.96, -3360.0}, {-0.48, -1680.0},
        {-0.24, -840.0}, {-0.12, -420.0}, {0.0, 0.0},
        {0.12, 144.0}, {0.24, 288.0}, {0.48, 576.0},
        {0.96, 1152.0}, {1.92, 2304.0}
    }
    vehicle_forces.tabulated_damper {
        name = root .. "left_rear_shock",
        first_body = body, first_point = {4.726, -0.392, 0.687},
        second_body = rear.axle, second_point = {4.512, -0.429, 0.225},
        table = rear_shock_table
    }
    vehicle_forces.tabulated_damper {
        name = root .. "right_rear_shock",
        first_body = body, first_point = {4.168, 0.425, 0.680},
        second_body = rear.axle, second_point = {4.440, 0.440, 0.227},
        table = rear_shock_table
    }
    if p.include_bump_stops ~= false then
        local rear_bumper_height = 0.080
        vehicle_forces.bumper {
            name = root .. "left_rear_jounce",
            kind = "jounce",
            sphere_on = "upper",
            lower_body = rear.axle, lower_point = {4.466, -0.515, 0.385},
            upper_body = body,
            upper_point = {4.466, -0.522, 0.488+rear_bumper_height},
            radius = bumper_radius, stiffness = bumper_stiffness,
            damping_factor = bumper_damping_factor,
            transition_depth = bumper_transition_depth
        }
        vehicle_forces.bumper {
            name = root .. "right_rear_jounce",
            kind = "jounce",
            sphere_on = "upper",
            lower_body = rear.axle, lower_point = {4.483, 0.522, 0.387},
            upper_body = body,
            upper_point = {4.460, 0.525, 0.493+rear_bumper_height},
            radius = bumper_radius, stiffness = bumper_stiffness,
            damping_factor = bumper_damping_factor,
            transition_depth = bumper_transition_depth
        }
    end

    local tires = nil
    if p.road_marker then
        local tire_color = colors.tire or "gray20"
        local common_tire = {
            road_marker = p.road_marker,
            mass = tire_mass,
            axle_direction = p.wheel_axle_direction or {0.0, 1.0, 0.0},
            forward_direction = forward,
            forward_speed = forward_speed,
            normal_damping_time_scale = p.tire_damping_time_scale or 0.0,
            regularization_speed = p.tire_regularization_speed or 0.1,
            tangential_model = p.tire_tangential_model or "expression",
            patch_length_by_load = p.tire_patch_length_by_load,
            cornering_stiffness_by_load =
                p.tire_cornering_stiffness_by_load,
            longitudinal_slip_stiffness_by_load =
                p.tire_longitudinal_slip_stiffness_by_load,
            longitudinal_relaxation_fraction =
                p.tire_longitudinal_relaxation_fraction,
            lateral_relaxation_fraction =
                p.tire_lateral_relaxation_fraction,
            shear_release_time = p.tire_shear_release_time,
            longitudinal_relaxation_length =
                p.tire_longitudinal_relaxation_length,
            lateral_relaxation_length = p.tire_lateral_relaxation_length,
            friction = p.tire_friction or 1.0,
            color = tire_color
        }
        local function add_tire(role, spindle, center, pressure)
            return gylt_245_75_r16 {
                name = root .. role,
                spindle = spindle,
                center = center,
                road_marker = common_tire.road_marker,
                mass = common_tire.mass,
                pressure_psi = pressure,
                axle_direction = common_tire.axle_direction,
                forward_direction = common_tire.forward_direction,
                forward_speed = common_tire.forward_speed,
                normal_damping_time_scale =
                    common_tire.normal_damping_time_scale,
                regularization_speed = common_tire.regularization_speed,
                tangential_model = common_tire.tangential_model,
                patch_length_by_load = common_tire.patch_length_by_load,
                cornering_stiffness_by_load =
                    common_tire.cornering_stiffness_by_load,
                longitudinal_slip_stiffness_by_load =
                    common_tire.longitudinal_slip_stiffness_by_load,
                longitudinal_relaxation_fraction =
                    common_tire.longitudinal_relaxation_fraction,
                lateral_relaxation_fraction =
                    common_tire.lateral_relaxation_fraction,
                shear_release_time = common_tire.shear_release_time,
                longitudinal_relaxation_length =
                    common_tire.longitudinal_relaxation_length,
                lateral_relaxation_length =
                    common_tire.lateral_relaxation_length,
                friction = common_tire.friction,
                color = common_tire.color
            }
        end
        local front_pressure = p.front_tire_pressure_psi or 60
        local rear_pressure = p.rear_tire_pressure_psi or 60
        tires = {
            left_front = add_tire("left_front_tire", left_front.spindle,
                left_front_wheel, front_pressure),
            right_front = add_tire("right_front_tire", right_front.spindle,
                right_front_wheel, front_pressure),
            left_rear = add_tire("left_rear_tire", rear.axle,
                left_rear_wheel, rear_pressure),
            right_rear = add_tire("right_rear_tire", rear.axle,
                right_rear_wheel, rear_pressure)
        }
    end

    local bodies = {body}
    local function append_bodies(items)
        for _, item in ipairs(items or {}) do
            table.insert(bodies, item)
        end
    end
    append_bodies(left_front.bodies)
    append_bodies(right_front.bodies)
    if stabilizer then append_bodies(stabilizer.bodies) end
    if steering then append_bodies(steering.bodies) end
    append_bodies(rear.bodies)
    if tires then
        table.insert(bodies, tires.left_front.body)
        table.insert(bodies, tires.right_front.body)
        table.insert(bodies, tires.left_rear.body)
        table.insert(bodies, tires.right_rear.body)
    end

    -- A statically initialized dynamic analysis holds these declared
    -- velocities aside during equilibrium and restores them before the
    -- dynamic consistency solve. Giving every body the common vehicle
    -- velocity avoids asking the velocity projection to manufacture the
    -- initial translation from the chassis alone. Each tire assembly adds
    -- the corresponding free-rolling wheel spin separately.
    if forward_speed ~= 0.0 then
        local tire_bodies = {}
        if tires then
            tire_bodies[tires.left_front.body] = true
            tire_bodies[tires.right_front.body] = true
            tire_bodies[tires.left_rear.body] = true
            tire_bodies[tires.right_rear.body] = true
        end
        for _, moving_body in ipairs(bodies) do
            if not tire_bodies[moving_body] then
                model_table {name = moving_body, velocity = initial_velocity}
            end
        end
    end

    return {
        body = body,
        body_graphic = body_graphic,
        left_front = left_front,
        right_front = right_front,
        steering = steering,
        rear = rear,
        tires = tires,
        bodies = bodies,
        wheel_centers = {
            left_front = left_front_wheel,
            right_front = right_front_wheel,
            left_rear = left_rear_wheel,
            right_rear = right_rear_wheel
        },
        limitations = {
            "rolling resistance and tire aligning/overturning moments are not yet applied",
            "the separate rollover Pogo contact is not yet represented"
        }
    }
end

return large_van
