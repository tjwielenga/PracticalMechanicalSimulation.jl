module ViewerData

export BodyTrajectory, GearTrajectory, PulleyTrajectory, BeltSpanTrajectory,
       BeltWrapTrajectory, GuideTrajectory, ConnectorTrajectory,
       SpanMeasurementTrajectory, DirectedDistanceMeasurementTrajectory,
       TorsionalSpringTrajectory, PlaneTrajectory, SphereTrajectory,
       JointTrajectory, ForceArrowTrajectory, TorqueArrowTrajectory,
       GraphicCylinderTrajectory, GraphicFrustumTrajectory,
       GraphicSurfacePatch, GraphicSurfaceTrajectory,
       GraphicMarkerTrajectory, XYFrameTrajectory,
       GraphicStyle, ViewerAppearance, ViewerBookmark, MechanismResult

"""Sampled motion and simple default appearance for one rigid body."""
struct BodyTrajectory{T}
    name::Symbol
    point_a::Matrix{T}
    point_b::Matrix{T}
    center::Matrix{T}
    angle::Vector{T}
    radius::T
    ellipsoid_axes::NTuple{3,T}
    show_default::Bool
    reference_length::T
end

BodyTrajectory(name::Symbol, point_a::Matrix{T}, point_b::Matrix{T},
        center::Matrix{T}, angle::Vector{T}, radius::T,
        ellipsoid_axes::NTuple{3,T}) where {T} =
    BodyTrajectory(name, point_a, point_b, center, angle, radius,
        ellipsoid_axes, true, zero(T))

"""Sampled center and orientation for a simple planar gear disc."""
struct GearTrajectory{T}
    name::Symbol
    center::Matrix{T}
    angle::Vector{T}
    radius::T
    half_width::T
    internal::Bool
end

"""Sampled center and orientation for a planar belt pulley."""
struct PulleyTrajectory{T}
    name::Symbol
    center::Matrix{T}
    angle::Vector{T}
    radius::T
    half_width::T
end

"""Sampled endpoints of one straight belt span."""
struct BeltSpanTrajectory{T}
    name::Symbol
    point_1::Matrix{T}
    point_2::Matrix{T}
end

"""Sampled entry and exit points for one directed pulley wrap."""
struct BeltWrapTrajectory{T}
    name::Symbol
    center::Matrix{T}
    point_1::Matrix{T}
    point_2::Matrix{T}
    axis::Matrix{T}
    radius::T
    direction::Int
end

function BeltWrapTrajectory(name::Symbol, center::Matrix{T},
        point_1::Matrix{T}, point_2::Matrix{T}, radius::T,
        direction::Int) where {T}
    axis = zeros(T, size(center))
    axis[:, 3] .= one(T)
    BeltWrapTrajectory(name, center, point_1, point_2, axis, radius,
        direction)
end

"""Sampled endpoints and radius for a simple guide or ground line."""
struct GuideTrajectory{T}
    name::Symbol
    point_a::Matrix{T}
    point_b::Matrix{T}
    radius::T
    kind::Symbol
end

GuideTrajectory(name::Symbol, point_a::Matrix{T}, point_b::Matrix{T},
        radius::T) where {T} =
    GuideTrajectory(name, point_a, point_b, radius, :guide)

"""Sampled endpoints and radius for a compliant force connector."""
struct ConnectorTrajectory{T}
    name::Symbol
    point_a::Matrix{T}
    point_b::Matrix{T}
    radius::T
    kind::Symbol
end

ConnectorTrajectory(name::Symbol, point_a::Matrix{T}, point_b::Matrix{T},
        radius::T) where {T} =
    ConnectorTrajectory(name, point_a, point_b, radius, :connector)

"""Sampled endpoints of a positive marker-to-marker distance measurement."""
struct SpanMeasurementTrajectory{T}
    name::Symbol
    point_1::Matrix{T}
    point_2::Matrix{T}
end

"""Sampled geometry of a signed point-to-plane distance measurement."""
struct DirectedDistanceMeasurementTrajectory{T}
    name::Symbol
    measured_point::Matrix{T}
    projected_point::Matrix{T}
    normal::Matrix{T}
    plane_x::Matrix{T}
    plane_y::Matrix{T}
end

"""Sampled center and endpoint orientations for a torsional spring symbol."""
struct TorsionalSpringTrajectory{T}
    name::Symbol
    center::Matrix{T}
    angle_1::Vector{T}
    angle_2::Vector{T}
    radius::T
end

"""Sampled tangent endpoints and out-of-plane half-width for a constraint plane."""
struct PlaneTrajectory{T}
    name::Symbol
    point_a::Matrix{T}
    point_b::Matrix{T}
    half_width::T
end

"""Sampled center and fixed radius for spherical contact geometry."""
struct SphereTrajectory{T}
    name::Symbol
    center::Matrix{T}
    radius::T
end

"""Sampled position of a named joint or contact point."""
struct JointTrajectory{T}
    name::Symbol
    position::Matrix{T}
    diameter::T
end

JointTrajectory(name::Symbol, position::Matrix{T}) where {T} =
    JointTrajectory(name, position, T(0.075))

"""Sampled application point and physical force vector for a load arrow."""
struct ForceArrowTrajectory{T}
    name::Symbol
    position::Matrix{T}
    force::Matrix{T}
    category::Symbol
    on_ground::Bool
end

ForceArrowTrajectory(name::Symbol, position::Matrix{T}, force::Matrix{T},
    category::Symbol) where {T} =
    ForceArrowTrajectory(name, position, force, category, false)

"""Sampled application point, magnitude, and axis for a torque symbol."""
struct TorqueArrowTrajectory{T}
    name::Symbol
    position::Matrix{T}
    torque::Vector{T}
    axis::Matrix{T}
    category::Symbol
    on_ground::Bool
end

TorqueArrowTrajectory(name::Symbol, position::Matrix{T}, torque::Vector{T},
        category::Symbol, on_ground::Bool = false) where {T} = begin
    axis = zeros(T, length(torque), 3)
    axis[:, 3] .= one(T)
    TorqueArrowTrajectory(name, position, torque, axis, category, on_ground)
end

"""Appearance overrides for an inferred graphical element."""
struct GraphicStyle
    visible::Bool
    show_default::Bool
    color::Union{Nothing,String}
    opacity::Union{Nothing,Float64}
end

GraphicStyle() = GraphicStyle(true, true, nothing, nothing)

"""Global appearance settings reconstructed from the embedded model."""
struct ViewerAppearance
    background::String
    body_palette::Vector{String}
    styles::Dict{Symbol,GraphicStyle}
    reaction_color::String
    applied_color::String
    show_reactions::Bool
    show_applied_loads::Bool
    show_torques::Bool
    show_ground_loads::Bool
    assembly_paths::Dict{Symbol,Vector{String}}
end

ViewerAppearance(background, body_palette, styles, reaction_color,
        applied_color, show_reactions, show_applied_loads, show_torques) =
    ViewerAppearance(background, body_palette, styles, reaction_color,
        applied_color, show_reactions, show_applied_loads, show_torques,
        false, Dict{Symbol,Vector{String}}())

ViewerAppearance(background, body_palette, styles, reaction_color,
        applied_color, show_reactions, show_applied_loads, show_torques,
        show_ground_loads) = ViewerAppearance(background, body_palette,
    styles, reaction_color, applied_color, show_reactions,
    show_applied_loads, show_torques, show_ground_loads,
    Dict{Symbol,Vector{String}}())

ViewerAppearance() = ViewerAppearance("white",
    ["steelblue", "darkorange", "seagreen", "orchid"],
    Dict{Symbol,GraphicStyle}(), "gold2", "darkorange2",
    true, true, true, false, Dict{Symbol,Vector{String}}())

"""Named sample to which the viewer can jump directly."""
struct ViewerBookmark
    label::String
    sample::Int
end

"""Sampled endpoints for an explicit cylindrical graphical primitive."""
struct GraphicCylinderTrajectory{T}
    name::Symbol
    point_a::Matrix{T}
    point_b::Matrix{T}
    radius::T
    color::String
    opacity::T
    deformation_group::String
    reference_point_a::NTuple{3,T}
    reference_point_b::NTuple{3,T}
    orientation_direction::Matrix{T}
    show_orientation_line::Bool
    graphic_shape::Symbol
    cross_section_size::NTuple{2,T}
end

GraphicCylinderTrajectory(name::Symbol, point_a::Matrix{T},
        point_b::Matrix{T}, radius::T, color::String, opacity::T) where {T} =
    GraphicCylinderTrajectory(name, point_a, point_b, radius, color, opacity,
        "", (zero(T), zero(T), zero(T)), (zero(T), zero(T), zero(T)),
        zeros(T, 0, 3), false, :cylinder, (2 * radius, 2 * radius))

GraphicCylinderTrajectory(name::Symbol, point_a::Matrix{T},
        point_b::Matrix{T}, radius::T, color::String, opacity::T,
        deformation_group::String, reference_point_a::NTuple{3,T},
        reference_point_b::NTuple{3,T}) where {T} =
    GraphicCylinderTrajectory(name, point_a, point_b, radius, color, opacity,
        deformation_group, reference_point_a, reference_point_b,
        zeros(T, 0, 3), false, :cylinder, (2 * radius, 2 * radius))

GraphicCylinderTrajectory(name::Symbol, point_a::Matrix{T},
        point_b::Matrix{T}, radius::T, color::String, opacity::T,
        deformation_group::String, reference_point_a::NTuple{3,T},
        reference_point_b::NTuple{3,T}, orientation_direction::Matrix{T},
        show_orientation_line::Bool) where {T} =
    GraphicCylinderTrajectory(name, point_a, point_b, radius, color, opacity,
        deformation_group, reference_point_a, reference_point_b,
        orientation_direction, show_orientation_line, :cylinder,
        (2 * radius, 2 * radius))

"""Sampled axis endpoints and end radii for a gear or conical frustum."""
struct GraphicFrustumTrajectory{T}
    name::Symbol
    kind::Symbol
    point_a::Matrix{T}
    point_b::Matrix{T}
    radius_a::T
    radius_b::T
    color::String
    opacity::T
end

"""One colored triangular patch within an indexed surface graphic."""
struct GraphicSurfacePatch
    name::Symbol
    faces::Vector{NTuple{3,Int}}
    color::String
    opacity::Float64
end

Base.:(==)(first::GraphicSurfacePatch, second::GraphicSurfacePatch) =
    first.name == second.name && first.faces == second.faces &&
    first.color == second.color && first.opacity == second.opacity
Base.isequal(first::GraphicSurfacePatch, second::GraphicSurfacePatch) =
    first == second
Base.hash(patch::GraphicSurfacePatch, seed::UInt) = hash(
    (patch.name, patch.faces, patch.color, patch.opacity), seed)

"""Sampled world positions and fixed topology for an indexed surface."""
struct GraphicSurfaceTrajectory{T}
    name::Symbol
    vertices::Array{T,3}
    patches::Vector{GraphicSurfacePatch}
    edges::Vector{NTuple{2,Int}}
    edge_color::String
    edge_width::T
    category::Symbol
    include_in_fit::Bool
end

GraphicSurfaceTrajectory(name::Symbol, vertices::Array{T,3},
        patches::Vector{GraphicSurfacePatch}, edges::Vector{NTuple{2,Int}},
        edge_color::String, edge_width::T) where {T} =
    GraphicSurfaceTrajectory(name, vertices, patches, edges, edge_color,
        edge_width, :graphic, true)

GraphicSurfaceTrajectory(name::Symbol, vertices::Array{T,3},
        patches::Vector{GraphicSurfacePatch}, edges::Vector{NTuple{2,Int}},
        edge_color::String, edge_width::T, category::Symbol) where {T} =
    GraphicSurfaceTrajectory(name, vertices, patches, edges, edge_color,
        edge_width, category, true)

"""Sampled pose and dimensions for an explicit marker-based graphic."""
struct GraphicMarkerTrajectory{T}
    name::Symbol
    shape::Symbol
    center::Matrix{T}
    angle::Vector{T}
    size::NTuple{3,T}
    color::String
    opacity::T
    category::Symbol
end

GraphicMarkerTrajectory(name::Symbol, shape::Symbol, center::Matrix{T},
        angle::Vector{T}, size::NTuple{3,T}, color::String,
        opacity::T) where {T} = GraphicMarkerTrajectory(name, shape, center,
    angle, size, color, opacity, :graphic)

"""Sampled pose and appearance of a marker-attached x-y frame."""
struct XYFrameTrajectory{T}
    name::Symbol
    origin::Matrix{T}
    x_direction::Matrix{T}
    y_direction::Matrix{T}
    z_direction::Matrix{T}
    axis_length::T
    plane_size::T
    plane_color::String
    plane_opacity::T
    label::String
    category::Symbol
end

XYFrameTrajectory(name::Symbol, origin::Matrix{T}, x_direction::Matrix{T},
        y_direction::Matrix{T}, z_direction::Matrix{T}, axis_length::T,
        plane_size::T, plane_color::String, plane_opacity::T,
        label::String) where {T} = XYFrameTrajectory(name, origin,
    x_direction, y_direction, z_direction, axis_length, plane_size,
    plane_color, plane_opacity, label, :graphic)

"""Completed simulation data needed by the provisional mechanism viewer."""
struct MechanismResult{T}
    title::String
    dimension::Symbol
    times::Vector{T}
    bodies::Vector{BodyTrajectory{T}}
    gears::Vector{GearTrajectory{T}}
    pulleys::Vector{PulleyTrajectory{T}}
    belt_spans::Vector{BeltSpanTrajectory{T}}
    belt_wraps::Vector{BeltWrapTrajectory{T}}
    joints::Vector{JointTrajectory{T}}
    force_arrows::Vector{ForceArrowTrajectory{T}}
    torque_arrows::Vector{TorqueArrowTrajectory{T}}
    guides::Vector{GuideTrajectory{T}}
    connectors::Vector{ConnectorTrajectory{T}}
    span_measurements::Vector{SpanMeasurementTrajectory{T}}
    directed_distance_measurements::Vector{
        DirectedDistanceMeasurementTrajectory{T}}
    torsional_springs::Vector{TorsionalSpringTrajectory{T}}
    planes::Vector{PlaneTrajectory{T}}
    spheres::Vector{SphereTrajectory{T}}
    graphic_cylinders::Vector{GraphicCylinderTrajectory{T}}
    graphic_frustums::Vector{GraphicFrustumTrajectory{T}}
    graphic_surfaces::Vector{GraphicSurfaceTrajectory{T}}
    graphic_markers::Vector{GraphicMarkerTrajectory{T}}
    xy_frames::Vector{XYFrameTrajectory{T}}
    appearance::ViewerAppearance
    signals::Vector{Pair{String,Vector{T}}}
    bookmarks::Vector{ViewerBookmark}
end

function MechanismResult(title::String, times::Vector{T},
        bodies::Vector{BodyTrajectory{T}}, gears::Vector{GearTrajectory{T}},
        pulleys::Vector{PulleyTrajectory{T}},
        belt_spans::Vector{BeltSpanTrajectory{T}},
        belt_wraps::Vector{BeltWrapTrajectory{T}},
        joints::Vector{Matrix{T}}, guides::Vector{GuideTrajectory{T}},
        connectors::Vector{ConnectorTrajectory{T}},
        torsional_springs::Vector{TorsionalSpringTrajectory{T}},
        planes::Vector{PlaneTrajectory{T}}, spheres::Vector{SphereTrajectory{T}},
        signals::Vector{Pair{String,Vector{T}}}) where {T}
    named_joints = [JointTrajectory(Symbol("joint$index"), position)
        for (index, position) in enumerate(joints)]
    MechanismResult(title, :planar, times, bodies, gears, pulleys, belt_spans,
        belt_wraps, named_joints, ForceArrowTrajectory{T}[],
        TorqueArrowTrajectory{T}[], guides,
        connectors, SpanMeasurementTrajectory{T}[],
        DirectedDistanceMeasurementTrajectory{T}[], torsional_springs, planes,
        spheres, GraphicCylinderTrajectory{T}[],
        GraphicFrustumTrajectory{T}[],
        GraphicSurfaceTrajectory{T}[],
        GraphicMarkerTrajectory{T}[], XYFrameTrajectory{T}[],
        ViewerAppearance(), signals,
        ViewerBookmark[])
end

function MechanismResult(title::String, times::Vector{T},
        bodies::Vector{BodyTrajectory{T}}, joints::Vector{Matrix{T}},
        signals::Vector{Pair{String,Vector{T}}}) where {T}
    return MechanismResult(title, times, bodies, GearTrajectory{T}[],
        PulleyTrajectory{T}[], BeltSpanTrajectory{T}[], BeltWrapTrajectory{T}[], joints,
        GuideTrajectory{T}[], ConnectorTrajectory{T}[],
        TorsionalSpringTrajectory{T}[],
        PlaneTrajectory{T}[], SphereTrajectory{T}[], signals)
end

function MechanismResult(title::String, times::Vector{T},
        bodies::Vector{BodyTrajectory{T}}, joints::Vector{Matrix{T}},
        guides::Vector{GuideTrajectory{T}},
        signals::Vector{Pair{String,Vector{T}}}) where {T}
    return MechanismResult(title, times, bodies, GearTrajectory{T}[],
        PulleyTrajectory{T}[], BeltSpanTrajectory{T}[], BeltWrapTrajectory{T}[], joints, guides,
        ConnectorTrajectory{T}[], TorsionalSpringTrajectory{T}[],
        PlaneTrajectory{T}[],
        SphereTrajectory{T}[], signals)
end

function MechanismResult(title::String, times::Vector{T},
        bodies::Vector{BodyTrajectory{T}}, joints::Vector{Matrix{T}},
        guides::Vector{GuideTrajectory{T}},
        connectors::Vector{ConnectorTrajectory{T}},
        signals::Vector{Pair{String,Vector{T}}}) where {T}
    return MechanismResult(title, times, bodies, GearTrajectory{T}[],
        PulleyTrajectory{T}[], BeltSpanTrajectory{T}[], BeltWrapTrajectory{T}[], joints,
        guides, connectors, TorsionalSpringTrajectory{T}[],
        PlaneTrajectory{T}[], SphereTrajectory{T}[], signals)
end

end
