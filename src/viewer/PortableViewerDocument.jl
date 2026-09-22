module PortableViewerDocument

using JSON
using HDF5
using LinearAlgebra
using ..ViewerData

export viewer_document, write_viewer_document, write_graphics

const FORMAT_NAME = "SimpView"
const FORMAT_VERSION = 4
const GRAPHICS_FORMAT = "SimpGraphics"
const GRAPHICS_FORMAT_VERSION = 1

const SCENE_FIELDS = (
    :bodies,
    :gears,
    :pulleys,
    :belt_spans,
    :belt_wraps,
    :joints,
    :force_arrows,
    :torque_arrows,
    :guides,
    :connectors,
    :span_measurements,
    :directed_distance_measurements,
    :torsional_springs,
    :planes,
    :spheres,
    :graphic_cylinders,
    :graphic_frustums,
    :graphic_surfaces,
    :graphic_markers,
    :xy_frames,
)

portable_value(value::Nothing) = nothing
portable_value(value::Union{Bool,Integer,AbstractString}) = value
portable_value(value::AbstractFloat) = begin
    isfinite(value) || throw(ArgumentError(
        "portable viewer documents require finite numeric values"))
    value
end
portable_value(value::Symbol) = String(value)
portable_value(value::Tuple) = portable_value.(collect(value))
portable_value(value::Pair) = Dict(
    "name" => portable_value(first(value)),
    "values" => portable_value(last(value)),
)
portable_value(value::AbstractDict) = Dict(
    string(key) => portable_value(item) for (key, item) in value)
portable_value(value::AbstractVector) = portable_value.(value)

# JSON.jl follows Julia's column-major iteration for multidimensional arrays.
# Write rows explicitly because trajectories are stored as samples by values.
portable_value(value::AbstractMatrix) = [
    [portable_value(value[row, column]) for column in axes(value, 2)]
    for row in axes(value, 1)
]

portable_value(value::AbstractArray{<:Any,3}) = [
    [[portable_value(value[first, second, third])
      for third in axes(value, 3)]
     for second in axes(value, 2)]
    for first in axes(value, 1)
]

function portable_struct(value)
    type = typeof(value)
    isstructtype(type) || throw(ArgumentError(
        "cannot write $(type) to a portable viewer document"))
    Dict(string(field) => portable_value(getfield(value, field))
        for field in fieldnames(type))
end

portable_value(value) = portable_struct(value)

function surface_center(vertices)
    vec(sum(vertices; dims = 1)) ./ size(vertices, 1)
end

function surface_basis(vertices, center, first_index, second_index)
    first_direction = collect(@view(vertices[first_index, :])) .- center
    first_length = norm(first_direction)
    first_length > 1.0e-12 || return nothing
    x_axis = first_direction ./ first_length
    second_direction = collect(@view(vertices[second_index, :])) .- center
    second_direction .-= dot(second_direction, x_axis) .* x_axis
    second_length = norm(second_direction)
    second_length > 1.0e-12 || return nothing
    y_axis = second_direction ./ second_length
    z_axis = cross(x_axis, y_axis)
    hcat(x_axis, y_axis, z_axis)
end

function surface_reference(vertices)
    center = surface_center(vertices)
    offsets = [collect(@view(vertices[index, :])) .- center
        for index in axes(vertices, 1)]
    first_index = argmax(norm.(offsets))
    first_direction = offsets[first_index]
    first_length = norm(first_direction)
    first_length > 1.0e-12 || return nothing
    x_axis = first_direction ./ first_length
    second_index = argmax([norm(cross(x_axis, offset))
        for offset in offsets])
    basis = surface_basis(vertices, center, first_index, second_index)
    isnothing(basis) && return nothing
    local_vertices = (Matrix(vertices) .- reshape(center, 1, :)) * basis
    (; center, basis, local_vertices, first_index, second_index)
end

function rotation_quaternion(rotation)
    trace_value = tr(rotation)
    x, y, z, w = if trace_value > 0
        scale = 2sqrt(trace_value + 1)
        ((rotation[3, 2] - rotation[2, 3]) / scale,
         (rotation[1, 3] - rotation[3, 1]) / scale,
         (rotation[2, 1] - rotation[1, 2]) / scale,
         scale / 4)
    elseif rotation[1, 1] > rotation[2, 2] &&
            rotation[1, 1] > rotation[3, 3]
        scale = 2sqrt(1 + rotation[1, 1] - rotation[2, 2] - rotation[3, 3])
        (scale / 4,
         (rotation[1, 2] + rotation[2, 1]) / scale,
         (rotation[1, 3] + rotation[3, 1]) / scale,
         (rotation[3, 2] - rotation[2, 3]) / scale)
    elseif rotation[2, 2] > rotation[3, 3]
        scale = 2sqrt(1 + rotation[2, 2] - rotation[1, 1] - rotation[3, 3])
        ((rotation[1, 2] + rotation[2, 1]) / scale,
         scale / 4,
         (rotation[2, 3] + rotation[3, 2]) / scale,
         (rotation[1, 3] - rotation[3, 1]) / scale)
    else
        scale = 2sqrt(1 + rotation[3, 3] - rotation[1, 1] - rotation[2, 2])
        ((rotation[1, 3] + rotation[3, 1]) / scale,
         (rotation[2, 3] + rotation[3, 2]) / scale,
         scale / 4,
         (rotation[2, 1] - rotation[1, 2]) / scale)
    end
    quaternion = [x, y, z, w]
    quaternion ./= norm(quaternion)
    quaternion[4] < 0 && (quaternion .*= -1)
    quaternion
end

function compact_surface(surface::GraphicSurfaceTrajectory)
    reference_vertices = @view surface.vertices[1, :, :]
    reference = surface_reference(reference_vertices)
    isnothing(reference) && return nothing
    samples = size(surface.vertices, 1)
    positions = zeros(Float64, samples, 3)
    quaternions = zeros(Float64, samples, 4)
    scale = max(maximum(norm(@view(reference.local_vertices[index, :]))
        for index in axes(reference.local_vertices, 1)), 1.0)
    tolerance = 1.0e-8 * scale
    for sample in 1:samples
        vertices = @view surface.vertices[sample, :, :]
        center = surface_center(vertices)
        basis = surface_basis(vertices, center, reference.first_index,
            reference.second_index)
        isnothing(basis) && return nothing
        for vertex in axes(vertices, 1)
            reconstructed = center + basis *
                @view(reference.local_vertices[vertex, :])
            norm(reconstructed - @view(vertices[vertex, :])) <= tolerance ||
                return nothing
        end
        positions[sample, :] .= center
        quaternions[sample, :] .= rotation_quaternion(basis)
    end
    Dict(
        "name" => String(surface.name),
        "encoding" => "rigid_pose",
        "vertices" => portable_value(reference.local_vertices),
        "position" => portable_value(positions),
        "quaternion" => portable_value(quaternions),
        "patches" => portable_value(surface.patches),
        "edges" => portable_value(surface.edges),
        "edge_color" => surface.edge_color,
        "edge_width" => surface.edge_width,
        "category" => String(surface.category),
        "include_in_fit" => surface.include_in_fit,
    )
end

function continuous_quaternion(rotation, previous = nothing)
    quaternion = rotation_quaternion(rotation)
    !isnothing(previous) && dot(quaternion, previous) < 0 &&
        (quaternion .*= -1)
    quaternion
end

function direction_basis(direction)
    y_axis = collect(direction)
    length = norm(y_axis)
    length > 1.0e-12 || return (Matrix{Float64}(I, 3, 3), 0.0)
    y_axis ./= length
    reference = abs(y_axis[3]) < 0.9 ? [0.0, 0.0, 1.0] : [1.0, 0.0, 0.0]
    x_axis = normalize(cross(y_axis, reference))
    z_axis = cross(x_axis, y_axis)
    (hcat(x_axis, y_axis, z_axis), length)
end

function cylinder_pose(point_a, point_b, radius_a, radius_b = radius_a)
    samples = size(point_a, 1)
    positions = zeros(Float64, samples, 3)
    quaternions = zeros(Float64, samples, 4)
    scales = zeros(Float64, samples, 3)
    previous = nothing
    for sample in 1:samples
        first = collect(@view(point_a[sample, :]))
        second = collect(@view(point_b[sample, :]))
        basis, length = direction_basis(second - first)
        quaternion = continuous_quaternion(basis, previous)
        positions[sample, :] .= (first + second) ./ 2
        quaternions[sample, :] .= quaternion
        scales[sample, :] .= (radius_a, length, radius_b)
        previous = quaternion
    end
    (; positions, quaternions, scales)
end

function vector_pose(position, values; magnitude_scale = 1.0)
    samples = size(position, 1)
    positions = Matrix{Float64}(position)
    quaternions = zeros(Float64, samples, 4)
    scales = zeros(Float64, samples, 3)
    previous = nothing
    for sample in 1:samples
        basis, magnitude = direction_basis(@view(values[sample, :]))
        quaternion = continuous_quaternion(basis, previous)
        quaternions[sample, :] .= quaternion
        display_magnitude = magnitude_scale * magnitude
        scales[sample, :] .= (display_magnitude, display_magnitude,
            display_magnitude)
        previous = quaternion
    end
    (; positions, quaternions, scales)
end

function frame_pose(frame::XYFrameTrajectory)
    samples = size(frame.origin, 1)
    quaternions = zeros(Float64, samples, 4)
    scales = ones(Float64, samples, 3)
    previous = nothing
    for sample in 1:samples
        basis = hcat(collect(@view(frame.x_direction[sample, :])),
            collect(@view(frame.y_direction[sample, :])),
            collect(@view(frame.z_direction[sample, :])))
        quaternion = continuous_quaternion(basis, previous)
        quaternions[sample, :] .= quaternion
        previous = quaternion
    end
    (; positions = frame.origin, quaternions, scales)
end

function planar_pose(center, angle; scale = (1.0, 1.0, 1.0))
    samples = size(center, 1)
    quaternions = zeros(Float64, samples, 4)
    scales = repeat(reshape(collect(scale), 1, 3), samples, 1)
    for sample in 1:samples
        half_angle = angle[sample] / 2
        quaternions[sample, :] .= (0.0, 0.0, sin(half_angle), cos(half_angle))
    end
    (; positions = center, quaternions, scales)
end

function characteristic_graphic_length(result::MechanismResult)
    points = Vector{Vector{Float64}}()
    for body in result.bodies
        center = collect(@view(body.center[1, :]))
        half_axes = collect(body.ellipsoid_axes) ./ 2
        push!(points, center .- half_axes)
        push!(points, center .+ half_axes)
        for field in (:point_a, :point_b)
            point = collect(@view(getfield(body, field)[1, :]))
            push!(points, point .- body.radius)
            push!(points, point .+ body.radius)
        end
    end
    for surface in result.graphic_surfaces
        surface.include_in_fit || continue
        for vertex in axes(surface.vertices, 2)
            push!(points, collect(@view(surface.vertices[1, vertex, :])))
        end
    end
    for cylinder in result.graphic_cylinders, field in (:point_a, :point_b)
        point = collect(@view(getfield(cylinder, field)[1, :]))
        push!(points, point .- cylinder.radius)
        push!(points, point .+ cylinder.radius)
    end
    for marker in result.graphic_markers
        center = collect(@view(marker.center[1, :]))
        half_size = collect(marker.size) ./ 2
        push!(points, center .- half_size)
        push!(points, center .+ half_size)
    end
    for joint in result.joints
        push!(points, collect(@view(joint.position[1, :])))
    end
    isempty(points) && return 1.0
    lower = reduce((a, b) -> min.(a, b), points)
    upper = reduce((a, b) -> max.(a, b), points)
    max(norm(upper - lower), 1.0e-6)
end

function add_track!(tracks, id, pose; deformation = nothing)
    compact_rows(values) = size(values, 1) > 1 && all(
        @view(values[row, :]) == @view(values[1, :])
        for row in 2:size(values, 1)) ? values[1:1, :] : values
    track = Dict(
        "id" => id,
        "position" => portable_value(compact_rows(pose.positions)),
        "quaternion" => portable_value(compact_rows(pose.quaternions)),
        "scale" => portable_value(compact_rows(pose.scales)),
    )
    isnothing(deformation) || (track["deformation"] = deformation)
    push!(tracks, track)
    id
end

function add_instance!(instances, name, mesh, track, category, color,
        opacity; local_position = [0.0, 0.0, 0.0],
        local_quaternion = [0.0, 0.0, 0.0, 1.0],
        local_scale = [1.0, 1.0, 1.0],
        include_in_fit = true,
        path = [titlecase(String(category)); split(String(name), '.')])
    push!(instances, Dict(
        "name" => String(name),
        "mesh" => mesh,
        "track" => track,
        "category" => String(category),
        "path" => String.(path),
        "color" => color,
        "opacity" => opacity,
        "local_position" => local_position,
        "local_quaternion" => local_quaternion,
        "local_scale" => local_scale,
        "include_in_fit" => include_in_fit,
    ))
end

function quaternion_rotation(quaternion)
    x, y, z, w = quaternion
    [1 - 2(y^2 + z^2)  2(x * y - z * w)  2(x * z + y * w);
     2(x * y + z * w)  1 - 2(x^2 + z^2)  2(y * z - x * w);
     2(x * z - y * w)  2(y * z + x * w)  1 - 2(x^2 + y^2)]
end

track_row(track, field, sample) = begin
    rows = track[field]
    Float64.(rows[min(sample, length(rows))])
end

function relative_track_transform(track, body_track)
    sample_count = maximum(length(track[field])
        for field in ("position", "quaternion", "scale"))
    body_samples = maximum(length(body_track[field])
        for field in ("position", "quaternion", "scale"))
    sample_count == body_samples || sample_count == 1 || return nothing
    all(all(isapprox.(track_row(body_track, "scale", sample), 1.0;
        atol = 1.0e-10)) for sample in 1:body_samples) || return nothing

    body_position = track_row(body_track, "position", 1)
    body_rotation = quaternion_rotation(
        track_row(body_track, "quaternion", 1))
    position = body_rotation' *
        (track_row(track, "position", 1) - body_position)
    rotation = body_rotation' * quaternion_rotation(
        track_row(track, "quaternion", 1))
    scale = track_row(track, "scale", 1)
    length_scale = max(norm(body_position), norm(position), 1.0)
    position_tolerance = 2.0e-7 * length_scale
    rotation_tolerance = 2.0e-7
    scale_tolerance = 2.0e-7 * max(norm(scale), 1.0)

    for sample in 1:max(sample_count, body_samples)
        current_body_position = track_row(body_track, "position", sample)
        current_body_rotation = quaternion_rotation(
            track_row(body_track, "quaternion", sample))
        current_position = track_row(track, "position", sample)
        current_rotation = quaternion_rotation(
            track_row(track, "quaternion", sample))
        current_scale = track_row(track, "scale", sample)
        norm(current_body_position + current_body_rotation * position -
            current_position) <= position_tolerance || return nothing
        norm(current_body_rotation * rotation - current_rotation) <=
            rotation_tolerance || return nothing
        norm(current_scale - scale) <= scale_tolerance || return nothing
    end
    (; position, rotation, scale)
end

function combine_local_transform!(instance, relative)
    local_position = Float64.(instance["local_position"])
    local_rotation = quaternion_rotation(
        Float64.(instance["local_quaternion"]))
    local_scale = Float64.(instance["local_scale"])
    instance["local_position"] = relative.position +
        relative.rotation * (relative.scale .* local_position)
    instance["local_quaternion"] = rotation_quaternion(
        relative.rotation * local_rotation)
    instance["local_scale"] = relative.scale .* local_scale
end

"""Replace rigidly repeated child histories by a body track and local pose."""
function share_body_tracks!(tracks, instances, body_tracks)
    isempty(body_tracks) && return tracks
    track_lookup = Dict(String(track["id"]) => track for track in tracks)
    users = Dict{String,Vector{Int}}()
    for (index, instance) in enumerate(instances)
        push!(get!(users, String(instance["track"]), Int[]), index)
    end
    protected = Set(values(body_tracks))
    for (track_id, indices) in users
        track_id in protected && continue
        haskey(track_lookup[track_id], "deformation") && continue
        owners = Set{String}()
        for index in indices
            owner = matching_body(instances[index]["name"], body_tracks)
            isnothing(owner) || push!(owners, owner)
        end
        length(owners) == 1 || continue
        owner = only(owners)
        all(matching_body(instances[index]["name"], body_tracks) == owner
            for index in indices) || continue
        relative = relative_track_transform(track_lookup[track_id],
            track_lookup[body_tracks[owner]])
        isnothing(relative) && continue
        for index in indices
            combine_local_transform!(instances[index], relative)
            instances[index]["track"] = body_tracks[owner]
        end
    end
    used = Set(String(instance["track"]) for instance in instances)
    union!(used, protected)
    filter(track -> String(track["id"]) in used, tracks)
end

graphic_path(group, name, subgroups...) =
    [String(group); collect(String.(subgroups)); split(String(name), '.')]

graphic_item_path(group, name, suffixes...) =
    [String(group); split(String(name), '.'); collect(String.(suffixes))]

function owned_graphic_path(owner, group, name = owner, suffixes...)
    owner_parts = split(String(owner), '.')
    name_parts = split(String(name), '.')
    remainder = name_parts[1:min(length(owner_parts), length(name_parts))] ==
        owner_parts ? name_parts[length(owner_parts) + 1:end] : name_parts
    !isempty(remainder) && first(remainder) == "graphics" && popfirst!(remainder)
    assembly = owner_parts[1:end-1]
    body = last(owner_parts)
    ["Model"; assembly; "Bodies"; body; String(group); remainder;
     collect(String.(suffixes))]
end

function element_graphic_path(group, name, suffixes...;
        assembly_paths = Dict{Symbol,Vector{String}}())
    parts = split(String(name), '.')
    matches = [String(element) for element in keys(assembly_paths)
        if String(name) == String(element) ||
           startswith(String(name), String(element) * ".")]
    if !isempty(matches)
        element = matches[argmax(length.(matches))]
        element_parts = split(element, '.')
        remainder = parts[length(element_parts) + 1:end]
        return ["Model"; assembly_paths[Symbol(element)]; String(group);
                last(element_parts); remainder; collect(String.(suffixes))]
    end
    ["Model"; parts[1:end-1]; String(group); last(parts);
     collect(String.(suffixes))]
end

function embedded_assembly_element(name, assembly_paths)
    parts = split(String(name), '.')
    matches = NamedTuple[]
    for element in keys(assembly_paths)
        element_parts = split(String(element), '.')
        length(element_parts) > length(parts) && continue
        for first in 1:length(parts)-length(element_parts)+1
            last = first + length(element_parts) - 1
            parts[first:last] == element_parts || continue
            push!(matches, (; element, element_parts, last))
        end
    end
    isempty(matches) && return nothing
    match = matches[argmax(length(item.element_parts) for item in matches)]
    (; assembly = assembly_paths[match.element],
       element = last(match.element_parts),
       remainder = parts[match.last + 1:end])
end

function categorized_graphic_path(group, name, body_tracks, suffixes...;
        assembly_paths = Dict{Symbol,Vector{String}}())
    declared = embedded_assembly_element(name, assembly_paths)
    if !isnothing(declared)
        return ["Model"; declared.assembly; String(group); declared.element;
                declared.remainder; collect(String.(suffixes))]
    end
    owner = matching_body(name, body_tracks)
    text = String(name)
    isnothing(owner) && (text == "ground" || startswith(text, "ground.")) &&
        (owner = "ground")
    if owner == "ground"
        parts = split(text, '.')
        return ["Model", "ground", String(group), parts[2:end]...,
                String.(suffixes)...]
    end
    isnothing(owner) ? element_graphic_path(group, name, suffixes...;
        assembly_paths) :
        owned_graphic_path(owner, group, name, suffixes...)
end

function force_graphic_path(name, kind;
        assembly_paths = Dict{Symbol,Vector{String}}())
    parts = split(String(name), '.')
    generated = length(parts) > 1 &&
        (occursin(r"^\d+$", last(parts)) ||
         occursin(r"^coordinate_\d+$", last(parts))) ? [pop!(parts)] : String[]
    element_graphic_path("Forces", join(parts, '.'), kind, generated...;
        assembly_paths)
end

function embedded_named_element(name, elements)
    parts = split(String(name), '.')
    matches = NamedTuple[]
    for element in elements
        element_parts = split(String(element), '.')
        length(element_parts) > length(parts) && continue
        for first in 1:length(parts)-length(element_parts)+1
            last = first + length(element_parts) - 1
            parts[first:last] == element_parts || continue
            push!(matches, (; element, element_parts))
        end
    end
    isempty(matches) && return nothing
    matches[argmax(length(match.element_parts) for match in matches)].element
end

function default_force_graphic_path(name, force_elements;
        assembly_paths = Dict{Symbol,Vector{String}}())
    element = embedded_named_element(name, force_elements)
    isnothing(element) ? nothing : force_graphic_path(
        element, "Default graphic"; assembly_paths)
end

joint_default_graphic_path(name;
        assembly_paths = Dict{Symbol,Vector{String}}()) =
    element_graphic_path("Joints", name, "Default graphic"; assembly_paths)

function surface_pose(surface)
    reference = surface_reference(@view(surface.vertices[1, :, :]))
    isnothing(reference) && return nothing
    samples = size(surface.vertices, 1)
    positions = zeros(Float64, samples, 3)
    quaternions = zeros(Float64, samples, 4)
    scales = ones(Float64, samples, 3)
    bases = Matrix{Float64}[]
    previous = nothing
    tolerance = 1.0e-8 * max(maximum(norm(
        @view(reference.local_vertices[index, :]))
        for index in axes(reference.local_vertices, 1)), 1.0)
    for sample in 1:samples
        vertices = @view surface.vertices[sample, :, :]
        center = surface_center(vertices)
        centered = Matrix(vertices) .- reshape(center, 1, :)
        factorization = svd(reference.local_vertices' * centered)
        basis = factorization.V * Diagonal([
            1.0, 1.0,
            sign(det(factorization.V * factorization.U'))]) *
            factorization.U'
        all(norm(center + basis * @view(reference.local_vertices[index, :]) -
            @view(vertices[index, :])) <= tolerance
            for index in axes(vertices, 1)) || return nothing
        quaternion = continuous_quaternion(basis, previous)
        positions[sample, :] .= center
        quaternions[sample, :] .= quaternion
        push!(bases, basis)
        previous = quaternion
    end
    (; positions, quaternions, scales, bases,
       local_vertices = reference.local_vertices)
end

function body_name(surface::GraphicSurfaceTrajectory)
    suffix = ".inertia_ellipsoid"
    name = String(surface.name)
    endswith(name, suffix) ? name[1:end-length(suffix)] : nothing
end

function matching_body(name, bodies)
    text = String(name)
    matches = [body for body in keys(bodies)
        if text == body || startswith(text, body * ".")]
    isempty(matches) ? nothing : matches[argmax(length.(matches))]
end

function body_local_cylinder(point_a, point_b, body_pose)
    local_a = body_pose.bases[1]' *
        (collect(@view(point_a[1, :])) - collect(@view(body_pose.positions[1, :])))
    local_b = body_pose.bases[1]' *
        (collect(@view(point_b[1, :])) - collect(@view(body_pose.positions[1, :])))
    scale = max(norm(local_b - local_a), 1.0)
    tolerance = 1.0e-8 * scale
    for sample in axes(point_a, 1)
        center = collect(@view(body_pose.positions[sample, :]))
        basis = body_pose.bases[sample]
        norm(center + basis * local_a - @view(point_a[sample, :])) <= tolerance ||
            return nothing
        norm(center + basis * local_b - @view(point_b[sample, :])) <= tolerance ||
            return nothing
    end
    basis, length = direction_basis(local_b - local_a)
    (; position = (local_a + local_b) ./ 2,
       quaternion = rotation_quaternion(basis), length)
end

function beam_deformation_reference(body, body_track, radius)
    Dict(
        "group" => String(body.name),
        "reference_track" => body_track,
        "local_position" => [0.0, 0.0, 0.0],
        "local_quaternion" => rotation_quaternion(
            first(direction_basis([1.0, 0.0, 0.0]))),
        "local_scale" => [radius, body.reference_length, radius],
    )
end

function surface_mesh!(meshes, surface, pose, mesh_prefix)
    mesh_ids = String[]
    for (index, patch) in enumerate(surface.patches)
        id = "$mesh_prefix.patch_$index"
        meshes[id] = Dict(
            "kind" => "surface",
            "vertices" => portable_value(pose.local_vertices),
            "faces" => portable_value(patch.faces),
        )
        push!(mesh_ids, id)
    end
    mesh_ids
end

"""Convert sampled graphics to reusable meshes, pose tracks, and instances."""
function scene_graph(result::MechanismResult)
    appearance = result.appearance
    assembly_paths = appearance.assembly_paths
    palette = appearance.body_palette
    meshes = Dict{String,Any}(
        "unit_cylinder" => Dict("kind" => "cylinder"),
        "unit_sphere" => Dict("kind" => "sphere"),
        "unit_box" => Dict("kind" => "box"),
        "unit_arrow" => Dict("kind" => "arrow"),
        "unit_torque" => Dict("kind" => "torque"),
    )
    tracks = Any[]
    instances = Any[]
    body_poses = Dict{String,Any}()
    body_tracks = Dict{String,String}()
    dynamic_surfaces = Any[]
    force_elements = Set{Symbol}()
    foreach(arrow -> push!(force_elements, arrow.name), result.force_arrows)
    foreach(arrow -> push!(force_elements, arrow.name), result.torque_arrows)
    foreach(item -> push!(force_elements, item.name), result.connectors)
    foreach(item -> push!(force_elements, item.name), result.belt_spans)
    foreach(item -> push!(force_elements, item.name), result.torsional_springs)
    for (index, body) in enumerate(result.bodies)
        name = String(body.name)
        color = palette[mod1(index, length(palette))]
        body_track = "body_frame:$name"
        add_track!(tracks, body_track, planar_pose(body.center, body.angle;
            scale = (1.0, 1.0, 1.0)))
        if body.show_default
            cylinder_track = "body:$name"
            deformation = body.reference_length > 0 ?
                beam_deformation_reference(body, body_track, body.radius) :
                nothing
            add_track!(tracks, cylinder_track, cylinder_pose(body.point_a,
                body.point_b, body.radius); deformation)
            add_instance!(instances, body.name, "unit_cylinder", cylinder_track,
                "geometry", color, 1.0;
                path = owned_graphic_path(name, "Geometry", name,
                    "Default body"))
        end
        add_instance!(instances, "$(body.name).inertia", "unit_sphere",
            body_track, "inertia", color, 0.35;
            local_scale = collect(body.ellipsoid_axes),
            path = owned_graphic_path(name, "Inertia", name))
        body_tracks[name] = body_track
    end

    for (index, gear) in enumerate(result.gears)
        track = "gear:$(gear.name)"
        point_a = copy(gear.center)
        point_b = copy(gear.center)
        point_a[:, 3] .-= gear.half_width
        point_b[:, 3] .+= gear.half_width
        add_track!(tracks, track, cylinder_pose(point_a, point_b, gear.radius))
        add_instance!(instances, gear.name, "unit_cylinder", track,
            "geometry", palette[mod1(index, length(palette))],
            gear.internal ? 0.22 : 0.72;
            path = categorized_graphic_path(
                "Geometry", gear.name, body_tracks; assembly_paths))
    end
    for (index, pulley) in enumerate(result.pulleys)
        track = "pulley:$(pulley.name)"
        point_a = copy(pulley.center)
        point_b = copy(pulley.center)
        point_a[:, 3] .-= pulley.half_width
        point_b[:, 3] .+= pulley.half_width
        add_track!(tracks, track,
            cylinder_pose(point_a, point_b, pulley.radius))
        add_instance!(instances, pulley.name, "unit_cylinder", track,
            "geometry", palette[mod1(index, length(palette))], 0.68;
            path = categorized_graphic_path(
                "Geometry", pulley.name, body_tracks; assembly_paths))
    end

    # Inertia surfaces provide an unambiguous sampled pose for spatial bodies.
    for surface in result.graphic_surfaces
        pose = surface_pose(surface)
        if isnothing(pose)
            push!(dynamic_surfaces, portable_struct(surface))
            continue
        end
        body = body_name(surface)
        track = "surface:$(surface.name)"
        add_track!(tracks, track, pose)
        if !isnothing(body)
            body_poses[body] = pose
            body_tracks[body] = track
        end
        mesh_ids = surface_mesh!(meshes, surface, pose, "surface:$(surface.name)")
        category = surface.category == :inertia ? "inertia" : "geometry"
        group = surface.category == :inertia ? "Inertia" : "Geometry"
        path_name = surface.category == :inertia && !isnothing(body) ?
            body : surface.name
        for (mesh, patch) in zip(mesh_ids, surface.patches)
            add_instance!(instances, "$(surface.name).$(patch.name)", mesh,
                track, category, patch.color, patch.opacity;
                include_in_fit = surface.include_in_fit,
                path = categorized_graphic_path(
                    group, path_name, body_tracks; assembly_paths))
        end
        if !isempty(surface.edges)
            mesh = "surface:$(surface.name).edges"
            meshes[mesh] = Dict("kind" => "lines",
                "vertices" => portable_value(pose.local_vertices),
                "edges" => portable_value(surface.edges))
            add_instance!(instances, "$(surface.name).edges", mesh, track,
                category, surface.edge_color, 1.0;
                include_in_fit = surface.include_in_fit,
                path = categorized_graphic_path(
                    group, path_name, body_tracks; assembly_paths))
        end
    end

    for cylinder in result.graphic_cylinders
        body = matching_body(cylinder.name, body_poses)
        local_pose = isnothing(body) ? nothing : body_local_cylinder(
            cylinder.point_a, cylinder.point_b, body_poses[body])
        if isnothing(local_pose)
            track = "cylinder:$(cylinder.name)"
            deformation = if !isempty(cylinder.deformation_group)
                reference_a = reshape(collect(cylinder.reference_point_a), 1, 3)
                reference_b = reshape(collect(cylinder.reference_point_b), 1, 3)
                reference_pose = cylinder_pose(reference_a, reference_b,
                    cylinder.radius)
                Dict(
                    "group" => cylinder.deformation_group,
                    "reference_track" => body_tracks[
                        cylinder.deformation_group],
                    "local_position" => only(portable_value(
                        reference_pose.positions)),
                    "local_quaternion" => only(portable_value(
                        reference_pose.quaternions)),
                    "local_scale" => only(portable_value(
                        reference_pose.scales)),
                )
            else
                nothing
            end
            add_track!(tracks, track, cylinder_pose(cylinder.point_a,
                cylinder.point_b, cylinder.radius); deformation)
            default_path = default_force_graphic_path(cylinder.name,
                force_elements; assembly_paths)
            add_instance!(instances, cylinder.name, "unit_cylinder", track,
                "geometry", cylinder.color, cylinder.opacity;
                path = isnothing(default_path) ? categorized_graphic_path(
                    "Geometry", cylinder.name, body_tracks; assembly_paths) :
                    default_path)
        else
            default_path = default_force_graphic_path(cylinder.name,
                force_elements; assembly_paths)
            add_instance!(instances, cylinder.name, "unit_cylinder",
                body_tracks[body], "geometry", cylinder.color,
                cylinder.opacity; local_position = local_pose.position,
                local_quaternion = local_pose.quaternion,
                local_scale = [cylinder.radius, local_pose.length,
                    cylinder.radius],
                path = isnothing(default_path) ? categorized_graphic_path(
                    "Geometry", cylinder.name, body_tracks; assembly_paths) :
                    default_path)
        end
    end


    for frustum in result.graphic_frustums
        mesh = "frustum:$(frustum.name)"
        meshes[mesh] = Dict("kind" => "frustum",
            "radius_a" => frustum.radius_a,
            "radius_b" => frustum.radius_b,
            "segments" => frustum.kind == :gear ? 32 : 24)
        track = "frustum:$(frustum.name)"
        add_track!(tracks, track, cylinder_pose(frustum.point_a,
            frustum.point_b, 1.0))
        add_instance!(instances, frustum.name, mesh, track, "geometry",
            frustum.color, frustum.opacity;
            path = categorized_graphic_path(
                "Geometry", frustum.name, body_tracks; assembly_paths))
    end

    for marker in result.graphic_markers
        mesh = marker.shape == :box ? "unit_box" : "unit_sphere"
        track = "marker:$(marker.name)"
        add_track!(tracks, track, planar_pose(marker.center, marker.angle;
            scale = marker.size))
        category = marker.category == :marker ? "markers" : "geometry"
        default_path = marker.category == :marker ? nothing :
            default_force_graphic_path(marker.name, force_elements;
                assembly_paths)
        add_instance!(instances, marker.name, mesh, track, category,
            marker.color, marker.opacity;
            path = isnothing(default_path) ? categorized_graphic_path(
                category == "markers" ? "Markers" : "Geometry",
                marker.name, body_tracks; assembly_paths) : default_path)
    end

    axis_quaternions = Dict(
        :x => rotation_quaternion(first(direction_basis([1.0, 0.0, 0.0]))),
        :y => [0.0, 0.0, 0.0, 1.0],
        :z => rotation_quaternion(first(direction_basis([0.0, 0.0, 1.0]))),
    )
    for frame in result.xy_frames
        track = "frame:$(frame.name)"
        add_track!(tracks, track, frame_pose(frame))
        is_joint = frame.category == :perp
        category = frame.category == :marker ? "markers" :
            is_joint ? "joints" : "frames"
        axes = is_joint ? ((:x, "#9a9a9a"), (:y, "#9a9a9a")) :
            ((:x, "#bd3131"), (:y, "#27914d"), (:z, "#315dbd"))
        for (axis, color) in axes
            add_instance!(instances, "$(frame.name).$axis", "unit_arrow",
                track, category, color, 1.0;
                local_quaternion = axis_quaternions[axis],
                local_scale = fill(frame.axis_length, 3),
                path = is_joint ? joint_default_graphic_path(frame.name;
                    assembly_paths) :
                    categorized_graphic_path(
                        category == "markers" ? "Markers" : "Frames",
                        frame.name, body_tracks; assembly_paths))
        end
        if frame.plane_size > 0
            add_instance!(instances, "$(frame.name).plane", "unit_box",
                track, category, frame.plane_color, frame.plane_opacity;
                local_scale = [frame.plane_size, frame.plane_size,
                    max(frame.plane_size * 0.006, 1.0e-6)],
                path = is_joint ? joint_default_graphic_path(frame.name;
                    assembly_paths) :
                    categorized_graphic_path(
                        category == "markers" ? "Markers" : "Frames",
                        frame.name, body_tracks; assembly_paths))
        end
    end

    for (field, category, color) in (
            (:guides, "joints", "#aeb6be"),
            (:connectors, "geometry", "#25282b"))
        for cylinder in getfield(result, field)
            track = "$(field):$(cylinder.name)"
            add_track!(tracks, track, cylinder_pose(cylinder.point_a,
                cylinder.point_b, cylinder.radius))
            path = if field == :guides
                joint_default_graphic_path(cylinder.name; assembly_paths)
            else
                force_graphic_path(cylinder.name, "Default graphic";
                    assembly_paths)
            end
            add_instance!(instances, cylinder.name, "unit_cylinder", track,
                category, color, 1.0;
                path)
        end
    end
    for span in result.belt_spans
        track = "belt_span:$(span.name)"
        add_track!(tracks, track, cylinder_pose(span.point_1, span.point_2,
            0.008))
        add_instance!(instances, span.name, "unit_cylinder", track,
            "geometry", "#25282b", 1.0;
            path = force_graphic_path(span.name, "Default graphic";
                assembly_paths))
    end

    for joint in result.joints
        samples = size(joint.position, 1)
        pose = (; positions = joint.position,
            quaternions = repeat([0.0 0.0 0.0 1.0], samples, 1),
            scales = repeat(reshape(fill(joint.diameter, 3), 1, 3), samples, 1))
        track = "joint:$(joint.name)"
        add_track!(tracks, track, pose)
        is_contact = any(sphere -> sphere.name == joint.name,
            result.spheres) || any(plane -> plane.name == joint.name,
            result.planes)
        add_instance!(instances, joint.name, "unit_sphere", track,
            "joints", "#b4bbc2", 1.0;
            path = is_contact ? force_graphic_path(joint.name,
                "Default graphic"; assembly_paths) :
                joint_default_graphic_path(joint.name; assembly_paths))
    end
    for sphere in result.spheres
        samples = size(sphere.center, 1)
        pose = (; positions = sphere.center,
            quaternions = repeat([0.0 0.0 0.0 1.0], samples, 1),
            scales = repeat(reshape(fill(2sphere.radius, 3), 1, 3), samples, 1))
        track = "sphere:$(sphere.name)"
        add_track!(tracks, track, pose)
        add_instance!(instances, sphere.name, "unit_sphere", track,
            "geometry", "#b4bbc2", 1.0;
            path = force_graphic_path(sphere.name, "Default graphic";
                assembly_paths))
    end

    graphic_length = characteristic_graphic_length(result)
    force_maximum = maximum((norm(@view(arrow.force[sample, :]))
        for arrow in result.force_arrows
        for sample in axes(arrow.force, 1)); init = 0.0)
    torque_maximum = maximum((abs(value) for arrow in result.torque_arrows
        for value in arrow.torque); init = 0.0)
    force_scale = force_maximum > 0 ? 0.2graphic_length / force_maximum : 1.0
    torque_scale = torque_maximum > 0 ? 0.16graphic_length / torque_maximum : 1.0
    for arrow in result.force_arrows
        track = "force:$(arrow.name):$(length(tracks))"
        add_track!(tracks, track, vector_pose(arrow.position, arrow.force;
            magnitude_scale = force_scale))
        color = arrow.category == :reaction ? appearance.reaction_color :
            appearance.applied_color
        add_instance!(instances, arrow.name, "unit_arrow", track, "loads",
            color, 1.0;
            path = force_graphic_path(arrow.name,
                arrow.category == :reaction ? "Reaction" : "Applied";
                assembly_paths))
    end
    for arrow in result.torque_arrows
        vectors = arrow.axis .* reshape(arrow.torque, :, 1)
        track = "torque:$(arrow.name):$(length(tracks))"
        add_track!(tracks, track, vector_pose(arrow.position, vectors;
            magnitude_scale = torque_scale))
        color = arrow.category == :reaction ? appearance.reaction_color :
            appearance.applied_color
        add_instance!(instances, arrow.name, "unit_torque", track,
            "torques", color, 1.0;
            path = force_graphic_path(arrow.name,
                arrow.category == :reaction ? "Reaction torque" :
                "Applied torque"; assembly_paths))
    end

    tracks = share_body_tracks!(tracks, instances, body_tracks)
    follow_targets = [Dict("name" => name, "track" => body_tracks[name])
        for name in sort!(collect(keys(body_tracks)))]
    Dict("encoding" => "mesh_instances", "meshes" => meshes,
        "tracks" => tracks, "instances" => instances,
        "follow_targets" => follow_targets,
        "dynamic_surfaces" => dynamic_surfaces)
end

function portable_value(surface::GraphicSurfaceTrajectory)
    compact = compact_surface(surface)
    isnothing(compact) ? portable_struct(surface) : compact
end

function portable_appearance(appearance::ViewerAppearance)
    Dict(
        "background" => appearance.background,
        "body_palette" => appearance.body_palette,
        "styles" => portable_value(appearance.styles),
        "reaction_color" => appearance.reaction_color,
        "applied_color" => appearance.applied_color,
        "show_reactions" => appearance.show_reactions,
        "show_applied_loads" => appearance.show_applied_loads,
        "show_torques" => appearance.show_torques,
        "show_ground_loads" => appearance.show_ground_loads,
    )
end

function viewer_choice(result::MechanismResult, label; include_signals = true)
    Dict(
        "label" => String(label),
        "times" => portable_value(result.times),
        "scene" => scene_graph(result),
        "signals" => include_signals ? portable_value(result.signals) : Any[],
        "bookmarks" => portable_value(result.bookmarks),
    )
end

"""Create a versioned, presentation-independent document for SimpView."""
function viewer_document(results::AbstractVector{<:MechanismResult};
        labels = ["Result $index" for index in eachindex(results)],
        choice_name = "Result", include_signals = true)
    isempty(results) && throw(ArgumentError(
        "at least one mechanism result is required"))
    length(labels) == length(results) || throw(DimensionMismatch(
        "viewer labels must match the number of results"))
    dimension = first(results).dimension
    all(result -> result.dimension == dimension, results) ||
        throw(ArgumentError(
            "all results in a viewer document must have one dimension"))
    Dict(
        "format" => FORMAT_NAME,
        "version" => FORMAT_VERSION,
        "title" => first(results).title,
        "dimension" => String(dimension),
        "choice_name" => String(choice_name),
        "appearance" => portable_appearance(first(results).appearance),
        "choices" => [viewer_choice(result, label; include_signals)
            for (result, label) in zip(results, labels)],
    )
end

viewer_document(result::MechanismResult; label = "Result") =
    viewer_document([result]; labels = [label])

"""Write a portable SimpView JSON document and return its absolute path."""
function write_viewer_document(path::AbstractString,
        results::AbstractVector{<:MechanismResult}; kwargs...)
    document = viewer_document(results; kwargs...)
    open(path, "w") do io
        JSON.print(io, document)
    end
    abspath(path)
end

function write_viewer_document(path::AbstractString,
        result::MechanismResult; kwargs...)
    document = viewer_document(result; kwargs...)
    open(path, "w") do io
        JSON.print(io, document)
    end
    abspath(path)
end

function viewer_document_bytes(document)
    io = IOBuffer()
    JSON.print(io, document)
    take!(io)
end

function row_matrix(rows, columns; T = Float64)
    matrix = Matrix{T}(undef, length(rows), columns)
    for (row, values) in enumerate(rows)
        length(values) == columns || throw(DimensionMismatch(
            "graphics rows must all have $columns entries"))
        matrix[row, :] .= values
    end
    matrix
end

function triangle_matrix(faces)
    triangles = NTuple{3,Int}[]
    for face in faces
        length(face) >= 3 || throw(ArgumentError(
            "a graphics face must contain at least three vertices"))
        for index in 2:(length(face) - 1)
            push!(triangles, (Int(face[1]), Int(face[index]),
                Int(face[index + 1])))
        end
    end
    row_matrix(triangles, 3; T = Int32)
end

function write_graphics_meshes(parent, meshes)
    group = create_group(parent, "meshes")
    ids = sort!(collect(keys(meshes)))
    group["id"] = ids
    for (index, id) in enumerate(ids)
        mesh = meshes[id]
        item = create_group(group, string(index))
        item["kind"] = String(mesh["kind"])
        haskey(mesh, "vertices") &&
            (item["vertices"] = row_matrix(mesh["vertices"], 3;
                T = Float32))
        haskey(mesh, "faces") &&
            (item["triangles"] = triangle_matrix(mesh["faces"]))
        haskey(mesh, "edges") &&
            (item["edges"] = row_matrix(mesh["edges"], 2; T = Int32))
        for name in ("radius_a", "radius_b")
            haskey(mesh, name) && (item[name] = Float32(mesh[name]))
        end
        haskey(mesh, "segments") &&
            (item["segments"] = Int32(mesh["segments"]))
    end
end

function packed_track_field(tracks, field, columns)
    offsets = Vector{Int32}(undef, length(tracks))
    counts = Vector{Int32}(undef, length(tracks))
    rows = Vector{Any}()
    sizehint!(rows, sum(length(track[field]) for track in tracks))
    for (index, track) in enumerate(tracks)
        offsets[index] = Int32(length(rows))
        counts[index] = Int32(length(track[field]))
        append!(rows, track[field])
    end
    row_matrix(rows, columns; T = Float32), offsets, counts
end

function write_compressed_graphics_matrix(group, name, values;
        compression = 3)
    isempty(values) && return group[name] = values
    chunk_rows = min(size(values, 1), 1024)
    group[name, chunk = (chunk_rows, size(values, 2)), shuffle = (),
        deflate = compression] = values
end

function write_graphics_tracks(parent, tracks)
    group = create_group(parent, "tracks")
    group["id"] = String[String(track["id"]) for track in tracks]
    for (field, columns) in (("position", 3), ("quaternion", 4),
            ("scale", 3))
        values, offsets, counts = packed_track_field(tracks, field, columns)
        write_compressed_graphics_matrix(group, field, values)
        group["$(field)_offset"] = offsets
        group["$(field)_count"] = counts
    end
    group["deformation_group"] = String[
        haskey(track, "deformation") ?
            String(track["deformation"]["group"]) : "" for track in tracks]
    group["deformation_reference_track"] = String[
        haskey(track, "deformation") ?
            String(track["deformation"]["reference_track"]) : ""
        for track in tracks]
    group["deformation_local_position"] = row_matrix([
        haskey(track, "deformation") ?
            track["deformation"]["local_position"] : [0.0, 0.0, 0.0]
        for track in tracks], 3; T = Float32)
    group["deformation_local_quaternion"] = row_matrix([
        haskey(track, "deformation") ?
            track["deformation"]["local_quaternion"] : [0.0, 0.0, 0.0, 1.0]
        for track in tracks], 4; T = Float32)
    group["deformation_local_scale"] = row_matrix([
        haskey(track, "deformation") ?
            track["deformation"]["local_scale"] : [1.0, 1.0, 1.0]
        for track in tracks], 3; T = Float32)
end

function write_graphics_instances(parent, instances)
    group = create_group(parent, "instances")
    group["name"] = String[String(item["name"]) for item in instances]
    group["mesh"] = String[String(item["mesh"]) for item in instances]
    group["track"] = String[String(item["track"]) for item in instances]
    group["category"] = String[String(item["category"]) for item in instances]
    group["color"] = String[String(item["color"]) for item in instances]
    group["opacity"] = Float32[item["opacity"] for item in instances]
    group["include_in_fit"] = UInt8[
        get(item, "include_in_fit", true) for item in instances]
    group["local_position"] = row_matrix(
        [item["local_position"] for item in instances], 3; T = Float32)
    group["local_quaternion"] = row_matrix(
        [item["local_quaternion"] for item in instances], 4; T = Float32)
    group["local_scale"] = row_matrix(
        [item["local_scale"] for item in instances], 3; T = Float32)
    group["path"] = String[join(String.(item["path"]), '\n')
        for item in instances]
end

function write_graphics_follow_targets(parent, targets)
    group = create_group(parent, "follow_targets")
    group["name"] = String[String(item["name"]) for item in targets]
    group["track"] = String[String(item["track"]) for item in targets]
end

function surface_vertex_array(vertices)
    samples = length(vertices)
    vertex_count = samples == 0 ? 0 : length(first(vertices))
    values = Array{Float32,3}(undef, samples, vertex_count, 3)
    for sample in 1:samples, vertex in 1:vertex_count
        length(vertices[sample]) == vertex_count || throw(DimensionMismatch(
            "a dynamic surface changed its vertex count"))
        length(vertices[sample][vertex]) == 3 || throw(DimensionMismatch(
            "surface vertices must have three coordinates"))
        values[sample, vertex, :] .= vertices[sample][vertex]
    end
    values
end

function write_dynamic_surfaces(parent, surfaces)
    group = create_group(parent, "dynamic_surfaces")
    group["count"] = length(surfaces)
    for (index, surface) in enumerate(surfaces)
        item = create_group(group, string(index))
        for name in ("name", "category", "edge_color")
            item[name] = String(surface[name])
        end
        item["edge_width"] = Float32(surface["edge_width"])
        item["include_in_fit"] = UInt8(surface["include_in_fit"])
        item["vertices"] = surface_vertex_array(surface["vertices"])
        item["edges"] = row_matrix(surface["edges"], 2; T = Int)
        patches = create_group(item, "patches")
        patches["count"] = length(surface["patches"])
        for (patch_index, patch) in enumerate(surface["patches"])
            patch_group = create_group(patches, string(patch_index))
            patch_group["name"] = String(patch["name"])
            patch_group["color"] = String(patch["color"])
            patch_group["opacity"] = Float32(patch["opacity"])
            patch_group["triangles"] = triangle_matrix(patch["faces"])
        end
    end
end

function write_graphics_scene(parent, scene)
    scene["encoding"] == "mesh_instances" || throw(ArgumentError(
        "native HDF5 graphics require mesh-instance encoding"))
    write_graphics_meshes(parent, scene["meshes"])
    write_graphics_tracks(parent, scene["tracks"])
    write_graphics_instances(parent, scene["instances"])
    write_graphics_follow_targets(parent, scene["follow_targets"])
    write_dynamic_surfaces(parent, scene["dynamic_surfaces"])
end

function write_graphics_signals(parent, signals)
    group = create_group(parent, "signals")
    group["name"] = String[String(signal["name"]) for signal in signals]
    isempty(signals) && return
    samples = length(first(signals)["values"])
    values = Matrix{Float64}(undef, samples, length(signals))
    for (column, signal) in enumerate(signals)
        length(signal["values"]) == samples || throw(DimensionMismatch(
            "viewer signals must have equal sample counts"))
        values[:, column] .= signal["values"]
    end
    group["values"] = values
end

function write_graphics_group(file, document)
    haskey(file, "viewer") && delete_object(file, "viewer")
    haskey(file, "graphics") && delete_object(file, "graphics")
    graphics = create_group(file, "graphics")
    try
        attributes(graphics)["format"] = GRAPHICS_FORMAT
        attributes(graphics)["format_version"] = GRAPHICS_FORMAT_VERSION
        graphics["title"] = String(document["title"])
        graphics["dimension"] = String(document["dimension"])
        graphics["choice_name"] = String(document["choice_name"])
        appearance = create_group(graphics, "appearance")
        for name in ("background", "reaction_color", "applied_color")
            appearance[name] = String(document["appearance"][name])
        end
        appearance["body_palette"] = String.(
            document["appearance"]["body_palette"])
        for name in ("show_reactions", "show_applied_loads", "show_torques",
                "show_ground_loads")
            appearance[name] = UInt8(document["appearance"][name])
        end
        choices = create_group(graphics, "choices")
        choices["label"] = String[String(choice["label"])
            for choice in document["choices"]]
        for (index, choice) in enumerate(document["choices"])
            item = create_group(choices, string(index))
            item["time"] = Float64.(choice["times"])
            write_graphics_scene(item, choice["scene"])
            write_graphics_signals(item, choice["signals"])
        end
    finally
        close(graphics)
    end
end

function copy_result_without_graphics(source_path, destination_path)
    h5open(source_path, "r") do source
        h5open(destination_path, "w") do destination
            source_attributes = attributes(source)
            destination_attributes = attributes(destination)
            for name in keys(source_attributes)
                destination_attributes[name] = read(source_attributes[name])
            end
            for name in keys(source)
                name in ("graphics", "viewer") && continue
                copy_object(source, name, destination, name)
            end
        end
    end
end

"""Write normalized graphics as native HDF5 datasets in a `.simp` result."""
function write_graphics(path::AbstractString,
        results::AbstractVector{<:MechanismResult};
        labels = ["Result $index" for index in eachindex(results)],
        choice_name = "Result", include_signals = length(results) > 1)
    document = viewer_document(results; labels, choice_name, include_signals)
    existing_graphics = h5open(path, "r") do file
        haskey(file, "graphics") || haskey(file, "viewer")
    end
    if existing_graphics
        permissions = stat(path).mode & 0o777
        temporary_path, io = mktemp(dirname(abspath(path)); cleanup = false)
        close(io)
        try
            copy_result_without_graphics(path, temporary_path)
            h5open(temporary_path, "r+") do file
                write_graphics_group(file, document)
            end
            chmod(temporary_path, permissions)
            mv(temporary_path, path; force = true)
        finally
            isfile(temporary_path) && rm(temporary_path)
        end
    else
        h5open(path, "r+") do file
            write_graphics_group(file, document)
        end
    end
    abspath(path)
end

function write_graphics(path::AbstractString, result::MechanismResult;
        label = "Result", kwargs...)
    write_graphics(path, [result]; labels = [label], kwargs...)
end

end
