-- Simple vehicle-body surface converted from the author's earlier model.
-- Dimensions and vertex coordinates are expressed in the body's reference
-- frame, or in the optional marker frame when marker is supplied.

local sim3d = require "sim3d"
local model_table, vector = sim3d.model_table, sim3d.vector
local required, positive = sim3d.required, sim3d.positive
local context = "simple_body"

local function simple_body(p)
    local name = required(p, "name", context)
    local body = required(p, "body", context)
    local origin = vector(p.origin or {0.0, 0.0, 0.0})

    local hood_length = positive(p, "hood_length", context)
    local roof_front = positive(p, "roof_front", context)
    local roof_back = positive(p, "roof_back", context)
    local window_back = positive(p, "window_back", context)
    local length = positive(p, "length", context)
    local front_hood_height = positive(p, "front_hood_height", context)
    local windshield_height = positive(p, "windshield_height", context)
    local roof_height = positive(p, "roof_height", context)
    local back_glass_height = positive(p, "back_glass_height", context)
    local trunk_height = positive(p, "trunk_height", context)
    local rocker_width = positive(p, "rocker_width", context)
    local body_width = positive(p, "body_width", context)
    local roof_width = positive(p, "roof_width", context)

    if not (hood_length < roof_front and roof_front < roof_back and
            roof_back < window_back and window_back < length) then
        error("simple_body longitudinal stations must increase from " ..
            "hood_length through window_back and remain below length")
    end
    if rocker_width > body_width or roof_width > body_width then
        error("simple_body rocker_width and roof_width must not exceed " ..
            "body_width")
    end

    local body_color = p.body_color or "steelblue"
    local glass_color = p.glass_color or "gray15"
    local glass_opacity = p.glass_opacity or 0.55
    if type(glass_opacity) ~= "number" or
            glass_opacity < 0 or glass_opacity > 1 then
        error("simple_body glass_opacity must be between zero and one")
    end

    local half_body_width = body_width/2
    local half_roof_width = roof_width/2
    local half_rocker_width = rocker_width/2

    -- Left-side vertices retain the names used in the historical routine.
    local L1 = origin+vector({0.0, -half_rocker_width, 0.0})
    local L2 = origin+vector({0.0, -half_body_width, front_hood_height})
    local L3 = origin+vector({hood_length, -half_body_width,
        windshield_height})
    local L4 = origin+vector({window_back, -half_body_width,
        back_glass_height})
    local L5 = origin+vector({length, -half_body_width, trunk_height})
    local L6 = origin+vector({length, -half_rocker_width, 0.0})
    local L7 = origin+vector({roof_front, -half_roof_width, roof_height})
    local L8 = origin+vector({roof_back, -half_roof_width, roof_height})

    -- Right-side vertices are the corresponding positive-y points.
    local R1 = origin+vector({0.0, half_rocker_width, 0.0})
    local R2 = origin+vector({0.0, half_body_width, front_hood_height})
    local R3 = origin+vector({hood_length, half_body_width,
        windshield_height})
    local R4 = origin+vector({window_back, half_body_width,
        back_glass_height})
    local R5 = origin+vector({length, half_body_width, trunk_height})
    local R6 = origin+vector({length, half_rocker_width, 0.0})
    local R7 = origin+vector({roof_front, half_roof_width, roof_height})
    local R8 = origin+vector({roof_back, half_roof_width, roof_height})

    -- Vertex indices are one-based. The winding follows the original ADAMS
    -- outlines and gives outward-facing normals on the left and right sides.
    local vertices = {
        L1, L2, L3, L4, L5, L6, L7, L8,
        R1, R2, R3, R4, R5, R6, R7, R8
    }
    local body_faces = {
        {6, 5, 2, 1},       -- left lower body
        {5, 4, 3, 2},       -- left middle body
        {9, 10, 13, 14},    -- right lower body
        {10, 11, 12, 13},   -- right middle body
        {1, 2, 10, 9},      -- grille
        {2, 3, 11, 10},     -- hood
        {7, 8, 16, 15},     -- roof
        {4, 5, 13, 12},     -- trunk lid
        {5, 6, 14, 13},     -- rear panel
        {6, 1, 9, 14}       -- underside
    }
    local glass_faces = {
        {4, 8, 7, 3},       -- left side glass
        {11, 15, 16, 12},   -- right side glass
        {3, 7, 15, 11},     -- windshield
        {8, 4, 12, 16}      -- rear window
    }

    local surface = {
        name = body .. ".graphics." .. name,
        shape = "surface",
        vertices = vertices,
        draw_edges = p.draw_edges ~= false,
        edge_color = p.edge_color or "gray25",
        edge_width = p.edge_width or 1.0,
        patches = {
            body = {
                color = body_color,
                opacity = p.body_opacity or 1.0,
                faces = body_faces
            },
            glass = {
                color = glass_color,
                opacity = glass_opacity,
                faces = glass_faces
            }
        }
    }
    if p.marker ~= nil then
        surface.marker = p.marker
    end
    model_table(surface)

    return {
        graphic = surface.name,
        vertex_count = #vertices,
        face_count = #body_faces+#glass_faces,
        roof_corners = {
            front_left = L7,
            rear_left = L8,
            front_right = R7,
            rear_right = R8
        }
    }
end

return simple_body
