module StoredResultViewer

using LinearAlgebra
using TOML
using PracticalMechanicalSimulation
using PracticalMechanicalSimulation.PlanarAppliedForces
using PracticalMechanicalSimulation.PlanarComponentAssembly
using PracticalMechanicalSimulation.SpatialComponentAssembly
using PracticalMechanicalSimulation.SpatialModeling
using PracticalMechanicalSimulation.SpatialDirectedDistances
using PracticalMechanicalSimulation.SpatialConstraints
using PracticalMechanicalSimulation.SpatialCoordinateCouplers
using PracticalMechanicalSimulation.SpatialGearPairs
using PracticalMechanicalSimulation.SpatialRackAndPinions
using PracticalMechanicalSimulation.SpatialSpans
using PracticalMechanicalSimulation.SpatialBelts
using PracticalMechanicalSimulation.SpatialAppliedForces
using PracticalMechanicalSimulation.SpatialBushings
using PracticalMechanicalSimulation.SpatialPlaneContacts
using PracticalMechanicalSimulation.SpatialFrictionForces
using PracticalMechanicalSimulation.SpatialTires
using PracticalMechanicalSimulation.SpatialMotionGenerators
using ..ViewerData

export planar_viewer_context, planar_mechanism_result,
       spatial_viewer_context, spatial_mechanism_result,
       load_viewer_model, model_mechanism_result,
       simulation_mechanism_result, stored_mechanism_result

const GRAPHIC_PALETTES = Dict(
    "default" => ["steelblue", "darkorange", "seagreen", "orchid"],
    "colorblind" => ["#0072B2", "#E69F00", "#009E73", "#CC79A7",
                     "#56B4E9", "#D55E00", "#F0E442", "black"])

# Full color parsing belongs to the renderer. Reconstruction only requires a
# nonempty portable color specification.
function validate_graphic_color(specification::AbstractString)
    isempty(strip(specification)) && throw(ArgumentError(
        "graphic color specifications must not be empty"))
    String(specification)
end

function is_angular_displacement_name(name)
    value = lowercase(String(name))
    value == "theta" || startswith(value, "theta_") ||
        startswith(value, "psi_") || occursin("angle", value)
end

function viewer_signal(component, name, values; prefix = "")
    label = "$prefix$component.$name"
    if is_angular_displacement_name(name)
        return "$label (deg)" => rad2deg.(collect(values))
    end
    label => collect(values)
end

function equivalent_pair_spacing(rotational_stiffness,
        translational_stiffness)
    rotational_stiffness > 0 && translational_stiffness > 0 || return nothing
    value = 2sqrt(rotational_stiffness / translational_stiffness)
    isfinite(value) && value > 0 ? value : nothing
end

function bounded_connection_length(candidates, mechanism_span;
        fraction = 0.025)
    nominal = fraction * mechanism_span
    inferred = if isempty(candidates)
        nominal
    else
        exp(sum(log, candidates) / length(candidates))
    end
    clamp(inferred, 0.65nominal, 2.5nominal)
end

function spatial_bushing_symbol_size(component, mechanism_span)
    candidates = Float64[]
    for (rotational, translational) in (
            (component.rotational_stiffness[1],
                component.translational_stiffness[2]),
            (component.rotational_stiffness[2],
                component.translational_stiffness[1]))
        spacing = equivalent_pair_spacing(rotational, translational)
        isnothing(spacing) || push!(candidates, spacing)
    end
    length = bounded_connection_length(candidates, mechanism_span)
    radius = clamp(0.24length, 0.005mechanism_span,
        0.012mechanism_span)
    (; length, radius)
end

function planar_bushing_symbol_size(component, mechanism_span)
    spacing = equivalent_pair_spacing(component.rotational_stiffness,
        component.translational_stiffness[2])
    candidates = isnothing(spacing) ? Float64[] : [spacing]
    length = bounded_connection_length(candidates, mechanism_span;
        fraction = 0.04)
    radius = clamp(0.22length, 0.006mechanism_span,
        0.014mechanism_span)
    (; length, radius)
end

function collect_model_tables!(result, table, path = String[])
    for (key, value) in table
        value isa AbstractDict || continue
        key == "graphics" && continue
        element_path = [path; key]
        haskey(value, "type") &&
            (result[Symbol(join(element_path, "."))] = value)
        collect_model_tables!(result, value, element_path)
    end
    result
end

function viewer_model_source(document)
    reconstruction = deepcopy(document)
    # Saved-result initialization has already done its job. Replaying it while
    # reconstructing geometry would make viewing depend on the source result
    # still existing and on its relative path resolving from the viewer's
    # current directory.
    delete!(reconstruction, "initial_conditions")
    io = IOBuffer()
    TOML.print(io, reconstruction)
    String(take!(io))
end

function require_graphic_boolean(table, field, default, label)
    value = get(table, field, default)
    value isa Bool || throw(ArgumentError("$label.$field must be Boolean"))
    value
end

function optional_graphic_opacity(table, label)
    haskey(table, "opacity") || return nothing
    opacity = table["opacity"]
    opacity isa Number || throw(ArgumentError(
        "$label.opacity must be numeric"))
    value = Float64(opacity)
    0 <= value <= 1 || throw(ArgumentError(
        "$label.opacity must be between zero and one"))
    value
end

function resolve_graphic_color(specification, aliases, label,
        resolving = Set{String}())
    specification isa AbstractString || throw(ArgumentError(
        "$label color must be a string"))
    name = String(specification)
    if haskey(aliases, name)
        name in resolving && throw(ArgumentError(
            "$label contains a cyclic color definition involving '$name'"))
        push!(resolving, name)
        resolved = resolve_graphic_color(aliases[name], aliases, label,
            resolving)
        delete!(resolving, name)
        return resolved
    end
    validate_graphic_color(name)
end

function parse_graphic_style(table, aliases, label)
    visible = require_graphic_boolean(table, "visible", true, label)
    show_default = require_graphic_boolean(
        table, "show_default", true, label)
    color = haskey(table, "color") ?
        resolve_graphic_color(table["color"], aliases, label) : nothing
    GraphicStyle(visible, show_default, color,
        optional_graphic_opacity(table, label))
end

function parse_viewer_appearance(document, element_tables)
    table = get(document, "graphics", Dict{String,Any}())
    table isa AbstractDict || throw(ArgumentError(
        "graphics must be a TOML table"))
    aliases = get(table, "colors", Dict{String,Any}())
    aliases isa AbstractDict || throw(ArgumentError(
        "graphics.colors must be a TOML table"))
    all(value -> value isa AbstractString, values(aliases)) ||
        throw(ArgumentError("graphics.colors values must be strings"))
    background = resolve_graphic_color(
        get(table, "background", "white"), aliases, "graphics.background")
    palette_specification = get(table, "body_palette", "default")
    palette = if palette_specification isa AbstractString
        haskey(GRAPHIC_PALETTES, palette_specification) ||
            throw(ArgumentError(
                "graphics.body_palette must be 'default', 'colorblind', or an array of colors"))
        [resolve_graphic_color(color, aliases, "graphics.body_palette")
         for color in GRAPHIC_PALETTES[palette_specification]]
    elseif palette_specification isa Vector
        isempty(palette_specification) && throw(ArgumentError(
            "graphics.body_palette must not be empty"))
        [resolve_graphic_color(color, aliases, "graphics.body_palette")
         for color in palette_specification]
    else
        throw(ArgumentError(
            "graphics.body_palette must be a palette name or array of colors"))
    end
    styles = Dict{Symbol,GraphicStyle}()
    for (name, element) in element_tables
        haskey(element, "graphics") || continue
        graphics = element["graphics"]
        graphics isa AbstractDict || throw(ArgumentError(
            "$name.graphics must be a TOML table"))
        styles[name] = parse_graphic_style(graphics, aliases,
            "$name.graphics")
    end
    loads = get(table, "loads", Dict{String,Any}())
    loads isa AbstractDict || throw(ArgumentError(
        "graphics.loads must be a TOML table"))
    reaction_color = resolve_graphic_color(
        get(loads, "reaction_color", "gold2"), aliases,
        "graphics.loads.reaction_color")
    applied_color = resolve_graphic_color(
        get(loads, "applied_color", "darkorange2"), aliases,
        "graphics.loads.applied_color")
    show_reactions = require_graphic_boolean(loads,
        "show_reactions", true, "graphics.loads")
    show_applied_loads = require_graphic_boolean(loads,
        "show_applied_loads", true, "graphics.loads")
    show_torques = require_graphic_boolean(loads,
        "show_torques", true, "graphics.loads")
    show_ground_loads = require_graphic_boolean(loads,
        "show_ground_loads", false, "graphics.loads")
    ViewerAppearance(background, palette, styles,
        reaction_color, applied_color, show_reactions,
        show_applied_loads, show_torques, show_ground_loads), aliases
end

function validate_catalog(stored, loaded)
    variables = loaded.layout.catalog.variables
    names = string.(getproperty.(variables, :name))
    components = string.(getproperty.(variables, :component))
    names == stored.variable_names && components == stored.variable_components ||
        throw(ArgumentError("stored variable catalog does not match the embedded model"))
    size(stored.values, 2) == length(variables) ||
        throw(DimensionMismatch("stored history has the wrong number of variables"))
end

function marker_history(marker, values)
    history = zeros(size(values, 1), 3)
    for sample in axes(values, 1)
        history[sample, 1:2] .=
            PlanarAppliedForces.point_marker_kinematics(
                marker.point, @view(values[sample, :])).position
    end
    history
end

function orientation_history(marker, values)
    [PlanarAppliedForces.marker_angle(
        marker.orientation, @view(values[sample, :]))
     for sample in axes(values, 1)]
end

function collect_graphic_shapes!(result, table, path = String[])
    for (key, value) in table
        value isa AbstractDict || continue
        shape_path = [path; key]
        haskey(value, "shape") && push!(result, (shape_path, value))
        collect_graphic_shapes!(result, value, shape_path)
    end
    result
end

function required_graphic_marker(loaded, specification, label)
    specification isa AbstractString || throw(ArgumentError(
        "$label marker name must be a string"))
    name = Symbol(specification)
    haskey(loaded.markers, name) || throw(ArgumentError(
        "$label references unknown marker '$name'"))
    loaded.markers[name]
end

function check_graphic_marker_owner(marker, owner_name, owner_type, label)
    valid = if owner_type == "rigid_body"
        !isnothing(marker.owner) && marker.owner.name == owner_name
    else
        isnothing(marker.owner)
    end
    valid || throw(ArgumentError(
        "$label marker '$(marker.name)' does not belong to '$owner_name'"))
    marker
end

function required_positive_graphic_number(table, field, label)
    haskey(table, field) || throw(ArgumentError("$label requires $field"))
    specification = table[field]
    specification isa Number || throw(ArgumentError(
        "$label.$field must be numeric"))
    value = Float64(specification)
    value > 0 || throw(ArgumentError("$label.$field must be positive"))
    value
end

function optional_positive_graphic_number(table, field, default, label)
    haskey(table, field) || return Float64(default)
    specification = table[field]
    specification isa Number || throw(ArgumentError(
        "$label.$field must be numeric"))
    value = Float64(specification)
    value > 0 || throw(ArgumentError("$label.$field must be positive"))
    value
end

function required_graphic_size(table, label)
    haskey(table, "size") || throw(ArgumentError("$label requires size"))
    specification = table["size"]
    specification isa Vector && length(specification) == 3 &&
        all(value -> value isa Number, specification) ||
        throw(ArgumentError("$label.size must contain three numbers"))
    values = Float64.(specification)
    all(>(0), values) || throw(ArgumentError(
        "$label.size values must be positive"))
    Tuple(values)
end

function centered_cylinder_history(marker, values, axis, cylinder_length)
    centers = marker_history(marker, values)
    angles = orientation_history(marker, values)
    point_a, point_b = similar(centers), similar(centers)
    for sample in axes(centers, 1)
        angle = angles[sample]
        direction = if axis == "x"
            [cos(angle), sin(angle), 0.0]
        elseif axis == "y"
            [-sin(angle), cos(angle), 0.0]
        else
            [0.0, 0.0, 1.0]
        end
        point_a[sample, :] .= centers[sample, :] .-
            cylinder_length / 2 .* direction
        point_b[sample, :] .= centers[sample, :] .+
            cylinder_length / 2 .* direction
    end
    point_a, point_b
end

function explicit_graphics(loaded, values, marker_histories, element_tables,
        appearance, aliases)
    cylinders = GraphicCylinderTrajectory{Float64}[]
    markers = GraphicMarkerTrajectory{Float64}[]
    body_names = sort!(collect(keys(loaded.bodies)))
    body_colors = Dict(name => appearance.body_palette[
        mod1(index, length(appearance.body_palette))]
        for (index, name) in enumerate(body_names))
    for owner_name in sort!(collect(keys(element_tables)))
        element = element_tables[owner_name]
        owner_type = string(get(element, "type", ""))
        owner_type in ("rigid_body", "ground", "marker") || continue
        haskey(element, "graphics") || continue
        graphics = element["graphics"]
        shapes = collect_graphic_shapes!(Tuple{Vector{String},Any}[], graphics)
        parent_style = get(appearance.styles, owner_name, GraphicStyle())
        parent_style.visible || continue
        default_color = owner_type == "rigid_body" ?
            body_colors[owner_name] : "dimgray"
        for (path, table) in shapes
            label = "$(owner_name).graphics.$(join(path, '.'))"
            visible = require_graphic_boolean(table, "visible", true, label)
            visible || continue
            color = haskey(table, "color") ?
                resolve_graphic_color(table["color"], aliases, label) :
                something(parent_style.color, default_color)
            opacity = something(optional_graphic_opacity(table, label),
                parent_style.opacity, 1.0)
            shape_specification = table["shape"]
            shape_specification isa AbstractString || throw(ArgumentError(
                "$label.shape must be a string"))
            shape = lowercase(String(shape_specification))
            name = Symbol(label)
            if shape == "cylinder"
                radius = required_positive_graphic_number(
                    table, "radius", label)
                point_a, point_b = if haskey(table, "markers")
                    endpoints = table["markers"]
                    endpoints isa Vector && length(endpoints) == 2 ||
                        throw(ArgumentError(
                            "$label.markers must contain two marker names"))
                    first_marker = check_graphic_marker_owner(
                        required_graphic_marker(loaded, endpoints[1], label),
                        owner_name, owner_type, label)
                    second_marker = check_graphic_marker_owner(
                        required_graphic_marker(loaded, endpoints[2], label),
                        owner_name, owner_type, label)
                    (marker_histories[first_marker.name],
                     marker_histories[second_marker.name])
                else
                    marker = check_graphic_marker_owner(required_graphic_marker(
                        loaded, get(table, "marker", nothing), label),
                        owner_name, owner_type, label)
                    axis = lowercase(string(get(table, "axis", "z")))
                    axis in ("x", "y", "z") || throw(ArgumentError(
                        "$label.axis must be 'x', 'y', or 'z'"))
                    cylinder_length = required_positive_graphic_number(
                        table, "length", label)
                    centered_cylinder_history(
                        marker, values, axis, cylinder_length)
                end
                push!(cylinders, GraphicCylinderTrajectory(name,
                    point_a, point_b, radius, color, opacity))
            elseif shape in ("sphere", "ellipsoid", "box")
                marker = check_graphic_marker_owner(required_graphic_marker(
                    loaded, get(table, "marker", nothing), label),
                    owner_name, owner_type, label)
                size = shape == "sphere" ? begin
                    radius = required_positive_graphic_number(
                        table, "radius", label)
                    (2radius, 2radius, 2radius)
                end : required_graphic_size(table, label)
                push!(markers, GraphicMarkerTrajectory(name, Symbol(shape),
                    marker_histories[marker.name],
                    orientation_history(marker, values), size, color, opacity))
            else
                throw(ArgumentError(
                    "$label has unsupported shape '$shape'"))
            end
        end
    end
    cylinders, markers
end

function body_endpoint_markers(loaded, body)
    candidates = [marker for marker in values(loaded.markers)
                  if !isnothing(marker.owner) && marker.owner.name == body.name &&
                     marker.point isa PlanarBodyPointMarker]
    isempty(candidates) && return nothing, nothing
    length(candidates) == 1 && return only(candidates), nothing
    best_pair = (candidates[1], candidates[2])
    best_distance = -Inf
    for first_index in 1:(length(candidates) - 1)
        for second_index in (first_index + 1):length(candidates)
            distance = norm(candidates[first_index].point.r_body -
                            candidates[second_index].point.r_body)
            if distance > best_distance
                best_distance = distance
                best_pair = (candidates[first_index], candidates[second_index])
            end
        end
    end
    best_pair
end

function body_trajectory(loaded, body, values, marker_histories)
    samples = size(values, 1)
    center = zeros(samples, 3)
    center[:, 1:2] .= values[:, body.position_variables]
    angle = collect(values[:, body.orientation_variable])
    marker_a, marker_b = body_endpoint_markers(loaded, body)
    if isnothing(marker_a)
        half_length = max(sqrt(body.inertia / body.mass), 0.1)
        point_a, point_b = zeros(samples, 3), zeros(samples, 3)
        for sample in 1:samples
            rotation = [cos(angle[sample]) -sin(angle[sample]);
                        sin(angle[sample])  cos(angle[sample])]
            point_a[sample, 1:2] .= center[sample, 1:2] +
                rotation * [-half_length, 0.0]
            point_b[sample, 1:2] .= center[sample, 1:2] +
                rotation * [half_length, 0.0]
        end
    elseif isnothing(marker_b)
        point_a, point_b = marker_histories[marker_a.name], copy(center)
    else
        point_a = marker_histories[marker_a.name]
        point_b = marker_histories[marker_b.name]
    end
    length_scale = maximum(norm(point_b[row, :] - point_a[row, :])
                           for row in 1:samples)
    radius = max(0.02, 0.06 * length_scale)
    BodyTrajectory(body.name, point_a, point_b, center, angle, radius,
        (max(0.22 * length_scale, radius), 2radius, 2radius))
end

function marker_name_for_point(loaded, point)
    for (name, marker) in loaded.markers
        marker.point === point && return name
    end
    throw(ArgumentError("connection refers to an unknown marker"))
end

point_is_on_ground(point) = point isa PlanarGroundPointMarker
marker_is_on_ground(marker) = isnothing(marker.owner)
orientation_is_on_ground(marker) = marker isa PlanarGroundOrientationMarker

function orientation_position_history(loaded, orientation, history_values,
        all_bodies)
    if orientation isa PlanarAppliedForces.PlanarBodyOrientationMarker
        for body in values(loaded.bodies)
            body.orientation_variable == orientation.theta_variable &&
                return all_bodies[body.name].center
        end
        throw(ArgumentError("orientation marker refers to an unknown body"))
    end
    zeros(size(history_values, 1), 3)
end

function vector_history(values, variables; sign = 1.0)
    history = zeros(size(values, 1), 3)
    history[:, 1:2] .= sign .* values[:, variables]
    history
end

scalar_history(values, variable; sign = 1.0) =
    sign .* collect(values[:, variable])

function separated_positions(first_position, second_position)
    maximum(norm(@view(first_position[sample, :]) .-
                 @view(second_position[sample, :]))
        for sample in axes(first_position, 1)) > 1.0e-10
end

function gear_center_point(gear, joint, body)
    joint.body_a === body && return joint.marker_a
    joint.body_b === body && return joint.marker_b
    throw(ArgumentError(
        "'$(gear.name)' references a joint not attached to '$(body.name)'"))
end

function time_matches(first_time, second_time)
    tolerance = 100eps(Float64) * max(
        abs(first_time), abs(second_time), 1.0)
    abs(first_time - second_time) <= tolerance
end

function simulation_viewer_history(stored)
    if stored.analysis_mode == :static &&
            hasproperty(stored, :static_snapshots) &&
            !isempty(stored.static_snapshots)
        snapshots = stored.static_snapshots
        times = Float64.(0:(length(snapshots) - 1))
        values = Matrix{Float64}(undef, length(snapshots),
            size(stored.values, 2))
        for (index, snapshot) in enumerate(snapshots)
            values[index, :] .= snapshot.values
        end
        bookmarks = ViewerBookmark[]
        for (sample, snapshot) in enumerate(snapshots)
            initial_configuration = snapshot.phase == :initial_conditions
            cycle_end = snapshot.phase == :dynamic_relaxation &&
                snapshot.status in (:waiting, :ready)
            newton_end = snapshot.phase in (:newton, :newton_polish) &&
                snapshot.status in (:converged, :failed)
            (initial_configuration || cycle_end || newton_end) || continue
            if initial_configuration
                label = snapshot.status == :entered ? "Model input" :
                    "Consistent initial conditions"
                push!(bookmarks, ViewerBookmark(label, sample))
                continue
            end
            phase = snapshot.phase == :dynamic_relaxation ?
                "relaxation" : replace(String(snapshot.phase), '_' => ' ')
            label = "$phase cycle $(snapshot.relaxation_cycle), " *
                "iteration $(snapshot.iteration): $(snapshot.status)"
            push!(bookmarks, ViewerBookmark(label, sample))
        end
        extra_signals = Pair{String,Vector{Float64}}[
            "Static force imbalance (N)" =>
                getproperty.(snapshots, :force_imbalance),
            "Static torque imbalance (N m)" =>
                getproperty.(snapshots, :torque_imbalance),
            "Static constraint error" =>
                getproperty.(snapshots, :constraint_error),
            "Static equivalent acceleration (m/s^2)" =>
                getproperty.(snapshots, :equivalent_acceleration)]
        return (; title = stored.title * " — static convergence",
            times, values, modal = false, bookmarks, extra_signals)
    end
    times = copy(stored.times)
    rows = [collect(@view(stored.values[index, :]))
        for index in axes(stored.values, 1)]
    for snapshot in sort(stored.health_snapshots; by = snapshot -> snapshot.time)
        position = searchsortedfirst(times, snapshot.time)
        if position <= length(times) &&
                time_matches(times[position], snapshot.time)
            times[position] = snapshot.time
            rows[position] = copy(snapshot.values)
        else
            insert!(times, position, snapshot.time)
            insert!(rows, position, copy(snapshot.values))
        end
    end
    values = Matrix{Float64}(undef, length(rows), size(stored.values, 2))
    for (index, row) in enumerate(rows)
        values[index, :] .= row
    end
    bookmarks = ViewerBookmark[]
    for (number, snapshot) in enumerate(stored.health_snapshots)
        sample = findfirst(time -> time_matches(time, snapshot.time), times)
        isnothing(sample) && continue
        severity = snapshot.severity == :check ? "Check" : "Warning"
        label = "$number $severity: $(snapshot.dominant_variable), E=" *
            string(round(snapshot.physical_error; sigdigits = 4)) *
            " at t=" * string(round(snapshot.time; sigdigits = 6))
        push!(bookmarks, ViewerBookmark(label, sample))
    end
    (; title = stored.title, times, values, modal = false, bookmarks,
       extra_signals = Pair{String,Vector{Float64}}[])
end

function viewer_history(stored; mode = 1, animation_samples = 121,
        animation_amplitude = 0.1)
    stored.analysis_mode == :modal ||
        return simulation_viewer_history(stored)
    1 <= mode <= length(stored.modal_eigenvalues) || throw(ArgumentError(
        "modal result has no mode $mode"))
    animation_samples >= 3 || throw(ArgumentError(
        "modal animation requires at least three samples"))
    animation_amplitude > 0 || throw(ArgumentError(
        "modal animation amplitude must be positive"))
    size(stored.values, 1) == 1 || throw(DimensionMismatch(
        "modal result must store one operating point"))
    size(stored.modal_mode_shapes, 1) == size(stored.values, 2) ||
        throw(DimensionMismatch(
            "modal mode shapes do not match the variable catalog"))

    eigenvalue = stored.modal_eigenvalues[mode]
    display_rate = max(abs(imag(eigenvalue)), abs(eigenvalue), eps(Float64))
    period = 2pi / display_rate
    times = collect(range(0.0, period; length = animation_samples))
    values = repeat(stored.values, animation_samples, 1)
    shape = @view stored.modal_mode_shapes[:, mode]
    for sample in eachindex(times)
        phase = 2pi * (sample - 1) / (animation_samples - 1)
        values[sample, :] .+= animation_amplitude .* real.(shape .* cis(phase))
    end
    title = stored.title * " — mode $mode, " *
        string(round(stored.modal_natural_frequencies_hz[mode]; sigdigits = 6)) *
        " Hz"
    (; title, times, values, modal = true, bookmarks = ViewerBookmark[],
       extra_signals = Pair{String,Vector{Float64}}[])
end

function spatial_marker_history(marker, values)
    history = zeros(size(values, 1), 3)
    for sample in axes(values, 1)
        history[sample, :] .= spatial_marker_position(
            marker, @view(values[sample, :]))
    end
    history
end

function spatial_frame_history(marker, values)
    samples = size(values, 1)
    origin = zeros(samples, 3)
    directions = ntuple(_ -> zeros(samples, 3), 3)
    for sample in axes(values, 1)
        state = @view values[sample, :]
        origin[sample, :] .= spatial_marker_position(marker, state)
        orientation = spatial_marker_orientation(marker, state)
        for axis in 1:3
            directions[axis][sample, :] .= @view orientation[:, axis]
        end
    end
    origin, directions
end

function unit_sphere_surface(longitudes = 18, latitudes = 10)
    vertices = NTuple{3,Float64}[(0.0, 0.0, 1.0)]
    for latitude in 1:(latitudes - 1)
        polar = pi * latitude / latitudes
        for longitude in 0:(longitudes - 1)
            azimuth = 2pi * longitude / longitudes
            push!(vertices, (sin(polar) * cos(azimuth),
                sin(polar) * sin(azimuth), cos(polar)))
        end
    end
    south = length(vertices) + 1
    push!(vertices, (0.0, 0.0, -1.0))
    ring(latitude, longitude) = 2 +
        (latitude - 1) * longitudes + mod(longitude, longitudes)
    faces = NTuple{3,Int}[]
    for longitude in 0:(longitudes - 1)
        push!(faces, (1, ring(1, longitude), ring(1, longitude + 1)))
    end
    for latitude in 1:(latitudes - 2), longitude in 0:(longitudes - 1)
        first = ring(latitude, longitude)
        second = ring(latitude, longitude + 1)
        third = ring(latitude + 1, longitude)
        fourth = ring(latitude + 1, longitude + 1)
        push!(faces, (first, third, fourth), (first, fourth, second))
    end
    for longitude in 0:(longitudes - 1)
        push!(faces, (south, ring(latitudes - 1, longitude + 1),
            ring(latitudes - 1, longitude)))
    end
    reduce(vcat, (reshape(collect(vertex), 1, 3) for vertex in vertices)),
        faces
end

function spatial_inertia_ellipsoids(loaded, values, marker_histories,
        appearance)
    isempty(loaded.bodies) && return GraphicSurfaceTrajectory{Float64}[]
    first_points = [collect(@view history[1, :])
        for history in Base.values(marker_histories)]
    for body in Base.values(loaded.bodies)
        push!(first_points,
            collect(@view values[1, body.position_variables]))
    end
    low = reduce((first, second) -> min.(first, second), first_points)
    high = reduce((first, second) -> max.(first, second), first_points)
    mechanism_span = max(maximum(high - low), 1.0)
    maximum_mass = maximum(body.mass for body in Base.values(loaded.bodies))
    unit_vertices, faces = unit_sphere_surface()
    trajectories = GraphicSurfaceTrajectory{Float64}[]
    body_names = sort!(collect(keys(loaded.bodies)))
    for (palette_index, name) in enumerate(body_names)
        body = loaded.bodies[name]
        style = get(appearance.styles, name, GraphicStyle())
        style.visible || continue
        principal = eigen(Symmetric(body.inertia))
        moments = principal.values
        axis_squared = [
            moments[2] + moments[3] - moments[1],
            moments[1] + moments[3] - moments[2],
            moments[1] + moments[2] - moments[3]
        ]
        largest = maximum(axis_squared)
        floor_value = max(1.0e-6 * largest, eps(Float64))
        shape = sqrt.(max.(axis_squared, floor_value))
        shape ./= cbrt(prod(shape))
        mass_radius = 0.055 * mechanism_span * cbrt(body.mass / maximum_mass)
        semiaxes = mass_radius .* shape
        vertices = zeros(Float64, size(values, 1), size(unit_vertices, 1), 3)
        for sample in axes(values, 1)
            center = @view values[sample, body.position_variables]
            body_orientation = rotation_matrix(
                @view values[sample, body.euler_parameter_variables])
            principal_orientation = body_orientation * principal.vectors
            for vertex in axes(unit_vertices, 1)
                vertices[sample, vertex, :] .= center .+
                    principal_orientation *
                    (semiaxes .* @view(unit_vertices[vertex, :]))
            end
        end
        color = something(style.color,
            appearance.body_palette[mod1(palette_index,
                length(appearance.body_palette))])
        opacity = min(something(style.opacity, 1.0), 0.24)
        ellipsoid_name = Symbol(name, ".inertia_ellipsoid")
        patch = GraphicSurfacePatch(ellipsoid_name, faces, color, opacity)
        push!(trajectories, GraphicSurfaceTrajectory(ellipsoid_name,
            vertices, [patch], NTuple{2,Int}[], "gray25", 1.0, :inertia))
    end
    trajectories
end

function spatial_torque_magnitude_axis(torque_vectors, fallback_axis)
    samples = size(torque_vectors, 1)
    magnitude = zeros(Float64, samples)
    axis = Matrix{Float64}(undef, samples, 3)
    for sample in 1:samples
        vector = @view torque_vectors[sample, :]
        magnitude[sample] = norm(vector)
        if magnitude[sample] > eps(Float64)
            axis[sample, :] .= vector ./ magnitude[sample]
        else
            axis[sample, :] .= @view fallback_axis[sample, :]
        end
    end
    magnitude, axis
end

function spatial_xy_frames(loaded, values, element_tables, aliases)
    trajectories = XYFrameTrajectory{Float64}[]
    for owner_name in sort!(collect(keys(element_tables)))
        element = element_tables[owner_name]
        owner_type = string(get(element, "type", ""))
        haskey(element, "graphics") || continue
        graphics = element["graphics"]
        graphics isa AbstractDict || throw(ArgumentError(
            "$owner_name.graphics must be a TOML table"))
        require_graphic_boolean(
            graphics, "visible", true, "$owner_name.graphics") || continue
        shapes = Tuple{Vector{String},Any}[]
        haskey(graphics, "shape") && push!(shapes, (String[], graphics))
        collect_graphic_shapes!(shapes, graphics)
        if owner_type != "marker"
            any(shape -> begin
                specification = get(last(shape), "shape", nothing)
                specification isa AbstractString &&
                    lowercase(String(specification)) == "xy_frame"
            end, shapes) && throw(ArgumentError(
                "an xy_frame graphic must be placed beneath a marker"))
            continue
        end
        for (path, table) in shapes
            label = isempty(path) ? "$(owner_name).graphics" :
                "$(owner_name).graphics.$(join(path, '.'))"
            shape_specification = table["shape"]
            shape_specification isa AbstractString || throw(ArgumentError(
                "$label.shape must be a string"))
            lowercase(String(shape_specification)) == "xy_frame" || continue
            require_graphic_boolean(table, "visible", true, label) || continue
            haskey(table, "marker") && throw(ArgumentError(
                "$label must be placed under its marker and must not " *
                "specify a separate marker"))
            haskey(loaded.markers, owner_name) || throw(ArgumentError(
                "$label does not belong to a model marker"))
            marker = loaded.markers[owner_name]
            axis_length = optional_positive_graphic_number(
                table, "axis_length", 0.35, label)
            plane_size = optional_positive_graphic_number(
                table, "plane_size", 0.45axis_length, label)
            plane_color = haskey(table, "plane_color") ?
                resolve_graphic_color(table["plane_color"], aliases,
                    "$label.plane_color") : "gray65"
            plane_opacity = something(
                optional_graphic_opacity(table, label), 0.18)
            frame_label = get(table, "label", "")
            frame_label isa AbstractString || throw(ArgumentError(
                "$label.label must be a string"))
            origin, directions = spatial_frame_history(marker, values)
            push!(trajectories, XYFrameTrajectory(Symbol(label), origin,
                directions..., axis_length, plane_size, plane_color,
                plane_opacity, String(frame_label)))
        end
    end
    trajectories
end


function spatial_body_graphic_frame(body, loaded, values, graphics, label)
    marker = if haskey(graphics, "marker")
        spatial_body_graphic_marker(body, loaded, graphics, label)
    else
        center_offset = loaded.body_reference_frames[
            body.name].center_of_mass_position
        SpatialBodyMarker(Symbol(body.name, ".__reference__"), body,
            -center_offset, Matrix{Float64}(I, 3, 3))
    end
    spatial_frame_history(marker, values)
end

function spatial_body_trajectories(body, loaded, values, table, appearance,
        palette_index)
    graphics = get(table, "graphics", Dict{String,Any}())
    graphics isa AbstractDict || throw(ArgumentError(
        "$(body.name).graphics must be a TOML table"))
    shape = lowercase(String(get(graphics, "shape", "box")))
    shape in ("box", "cylinder") || throw(ArgumentError(
        "the spatial viewer supports body.graphics.shape = 'box' or 'cylinder'"))
    box_size = haskey(graphics, "size") ?
        required_graphic_size(graphics, "$(body.name).graphics") :
        (0.6, 0.3, 0.2)
    style = get(appearance.styles, body.name, GraphicStyle())
    style.visible && style.show_default ||
        return GraphicCylinderTrajectory{Float64}[]
    color = something(style.color,
        appearance.body_palette[mod1(palette_index,
            length(appearance.body_palette))])
    opacity = something(style.opacity, 1.0)
    label = "$(body.name).graphics"
    center, directions = spatial_body_graphic_frame(
        body, loaded, values, graphics, label)
    if shape == "cylinder"
        radius = required_positive_graphic_number(
            graphics, "radius", "$(body.name).graphics")
        length_value = required_positive_graphic_number(
            graphics, "length", "$(body.name).graphics")
        axis_name = lowercase(string(get(graphics, "axis", "z")))
        axis_name in ("x", "y", "z") || throw(ArgumentError(
            "$(body.name).graphics.axis must be 'x', 'y', or 'z'"))
        axis_index = findfirst(==(axis_name), ("x", "y", "z"))
        point_a = zeros(size(values, 1), 3)
        point_b = similar(point_a)
        for sample in axes(values, 1)
            marker_axis = @view directions[axis_index][sample, :]
            half_axis = 0.5length_value .* marker_axis
            point_a[sample, :] .= @view(center[sample, :]) .- half_axis
            point_b[sample, :] .= @view(center[sample, :]) .+ half_axis
        end
        return GraphicCylinderTrajectory{Float64}[
            GraphicCylinderTrajectory(Symbol(body.name, ".cylinder"),
                point_a, point_b, radius, color, opacity)]
    end
    half = collect(box_size) ./ 2
    corners = [[sx * half[1], sy * half[2], sz * half[3]]
        for sx in (-1, 1), sy in (-1, 1), sz in (-1, 1)]
    corner_list = vec(corners)
    edge_pairs = Tuple{Int,Int}[]
    for first_index in eachindex(corner_list)
        for second_index in (first_index + 1):length(corner_list)
            difference_count = count(value -> !iszero(value),
                corner_list[first_index] .- corner_list[second_index])
            difference_count == 1 && push!(edge_pairs,
                (first_index, second_index))
        end
    end
    radius = max(0.0125 * maximum(box_size), 1.0e-4)
    trajectories = GraphicCylinderTrajectory{Float64}[]
    for (edge, (first_index, second_index)) in enumerate(edge_pairs)
        point_a = zeros(size(values, 1), 3)
        point_b = similar(point_a)
        for sample in axes(values, 1)
            orientation = hcat((@view(directions[axis][sample, :])
                for axis in 1:3)...)
            point_a[sample, :] .= @view(center[sample, :]) .+
                orientation * corner_list[first_index]
            point_b[sample, :] .= @view(center[sample, :]) .+
                orientation * corner_list[second_index]
        end
        push!(trajectories, GraphicCylinderTrajectory(
            Symbol(body.name, ".edge_", edge), point_a, point_b,
            radius, color, opacity))
    end
    trajectories
end

function spatial_body_graphic_marker(body, loaded, graphics, label)
    haskey(graphics, "marker") || throw(ArgumentError(
        "$label requires marker"))
    marker = required_graphic_marker(loaded, graphics["marker"], label)
    marker isa SpatialBodyMarker && marker.body.name == body.name ||
        throw(ArgumentError(
            "$label marker '$(marker.name)' does not belong to '$(body.name)'"))
    marker
end

function spatial_body_graphic_appearance(body, appearance, palette_index)
    style = get(appearance.styles, body.name, GraphicStyle())
    style.visible && style.show_default || return nothing
    color = something(style.color,
        appearance.body_palette[mod1(palette_index,
            length(appearance.body_palette))])
    opacity = something(style.opacity, 1.0)
    (; color, opacity)
end

function spatial_graphic_marker_for_owner(loaded, specification, owner_name,
        owner_type, label)
    marker = required_graphic_marker(loaded, specification, label)
    valid = if owner_type == "rigid_body"
        (marker isa SpatialBodyMarker || marker isa SpatialFloatingMarker) &&
            marker.body.name == owner_name
    elseif owner_type == "ground"
        marker isa SpatialGroundMarker
    elseif owner_type == "marker"
        marker.name == owner_name
    else
        false
    end
    valid || throw(ArgumentError(
        "$label marker '$(marker.name)' does not belong to '$owner_name'"))
    marker
end

function spatial_cylinder_trajectories(loaded, values, element_tables,
        appearance, aliases)
    trajectories = GraphicCylinderTrajectory{Float64}[]
    body_names = sort!(collect(keys(loaded.bodies)))
    body_colors = Dict(name => appearance.body_palette[
        mod1(index, length(appearance.body_palette))]
        for (index, name) in enumerate(body_names))
    for owner_name in sort!(collect(keys(element_tables)))
        element = element_tables[owner_name]
        owner_type = string(get(element, "type", ""))
        owner_type in ("rigid_body", "ground", "marker") || continue
        graphics = get(element, "graphics", nothing)
        graphics isa AbstractDict || continue
        parent_style = get(appearance.styles, owner_name, GraphicStyle())
        parent_style.visible || continue
        shapes = Tuple{Vector{String},Any}[]
        haskey(graphics, "shape") && push!(shapes, (String[], graphics))
        collect_graphic_shapes!(shapes, graphics)
        for (path, table) in shapes
            shape = get(table, "shape", nothing)
            shape isa AbstractString && lowercase(String(shape)) == "cylinder" ||
                continue
            # A body's root cylinder is already handled by
            # spatial_body_trajectories. This pass adds nested primitives.
            isempty(path) && owner_type == "rigid_body" && continue
            label = isempty(path) ? "$(owner_name).graphics" :
                "$(owner_name).graphics.$(join(path, '.'))"
            require_graphic_boolean(table, "visible", true, label) || continue
            radius = required_positive_graphic_number(table, "radius", label)
            point_a, point_b = if haskey(table, "markers")
                endpoints = table["markers"]
                endpoints isa Vector && length(endpoints) == 2 ||
                    throw(ArgumentError(
                        "$label.markers must contain two marker names"))
                first_marker = spatial_graphic_marker_for_owner(loaded,
                    endpoints[1], owner_name, owner_type, label)
                second_marker = spatial_graphic_marker_for_owner(loaded,
                    endpoints[2], owner_name, owner_type, label)
                (spatial_marker_history(first_marker, values),
                 spatial_marker_history(second_marker, values))
            else
                marker = if owner_type == "marker" &&
                        !haskey(table, "marker")
                    loaded.markers[owner_name]
                else
                    spatial_graphic_marker_for_owner(loaded,
                        get(table, "marker", nothing), owner_name,
                        owner_type, label)
                end
                length_value = required_positive_graphic_number(
                    table, "length", label)
                axis_name = lowercase(string(get(table, "axis", "z")))
                axis_index = findfirst(==(axis_name), ("x", "y", "z"))
                isnothing(axis_index) && throw(ArgumentError(
                    "$label.axis must be 'x', 'y', or 'z'"))
                center, directions = spatial_frame_history(marker, values)
                (center .- length_value / 2 .* directions[axis_index],
                 center .+ length_value / 2 .* directions[axis_index])
            end
            default_color = owner_type == "rigid_body" ?
                body_colors[owner_name] : "dimgray"
            color = haskey(table, "color") ? resolve_graphic_color(
                table["color"], aliases, label) :
                something(parent_style.color, default_color)
            opacity = something(optional_graphic_opacity(table, label),
                parent_style.opacity, 1.0)
            push!(trajectories, GraphicCylinderTrajectory(Symbol(label),
                point_a, point_b, radius, color, opacity))
        end
    end
    trajectories
end

function spatial_sphere_trajectories(loaded, values, element_tables,
        appearance, aliases)
    trajectories = GraphicMarkerTrajectory{Float64}[]
    body_names = sort!(collect(keys(loaded.bodies)))
    body_colors = Dict(name => appearance.body_palette[
        mod1(index, length(appearance.body_palette))]
        for (index, name) in enumerate(body_names))
    for owner_name in sort!(collect(keys(element_tables)))
        element = element_tables[owner_name]
        owner_type = string(get(element, "type", ""))
        owner_type in ("rigid_body", "ground", "marker") || continue
        graphics = get(element, "graphics", nothing)
        graphics isa AbstractDict || continue
        parent_style = get(appearance.styles, owner_name, GraphicStyle())
        parent_style.visible || continue
        shapes = Tuple{Vector{String},Any}[]
        haskey(graphics, "shape") && push!(shapes, (String[], graphics))
        collect_graphic_shapes!(shapes, graphics)
        for (path, table) in shapes
            shape = get(table, "shape", nothing)
            shape isa AbstractString && lowercase(String(shape)) == "sphere" ||
                continue
            label = isempty(path) ? "$(owner_name).graphics" :
                "$(owner_name).graphics.$(join(path, '.'))"
            require_graphic_boolean(table, "visible", true, label) || continue
            marker = if owner_type == "marker" &&
                    !haskey(table, "marker")
                loaded.markers[owner_name]
            else
                spatial_graphic_marker_for_owner(loaded,
                    get(table, "marker", nothing), owner_name,
                    owner_type, label)
            end
            radius = required_positive_graphic_number(table, "radius", label)
            default_color = owner_type == "rigid_body" ?
                body_colors[owner_name] : "dimgray"
            color = haskey(table, "color") ? resolve_graphic_color(
                table["color"], aliases, label) :
                something(parent_style.color, default_color)
            opacity = something(optional_graphic_opacity(table, label),
                parent_style.opacity, 1.0)
            center = spatial_marker_history(marker, values)
            angle = zeros(size(values, 1))
            push!(trajectories, GraphicMarkerTrajectory(Symbol(label),
                :sphere, center, angle, (2radius, 2radius, 2radius), color,
                opacity))
        end
    end
    trajectories
end

function spatial_gear_trajectory(body, loaded, values, table, appearance,
        palette_index)
    graphics = table["graphics"]
    label = "$(body.name).graphics"
    marker = spatial_body_graphic_marker(body, loaded, graphics, label)
    pitch_radius = required_positive_graphic_number(
        graphics, "pitch_radius", label)
    width = required_positive_graphic_number(graphics, "width", label)
    cone_height = haskey(graphics, "cone_height") ?
        required_positive_graphic_number(graphics, "cone_height", label) :
        nothing
    if !isnothing(cone_height) && cone_height <= width / 2
        throw(ArgumentError(
            "$label.cone_height must be greater than half its width"))
    end
    graphic_appearance = spatial_body_graphic_appearance(
        body, appearance, palette_index)
    isnothing(graphic_appearance) && return nothing
    center, directions = spatial_frame_history(marker, values)
    point_a = center .- width / 2 .* directions[3]
    point_b = center .+ width / 2 .* directions[3]
    radius_a, radius_b = if isnothing(cone_height)
        (pitch_radius, pitch_radius)
    else
        (pitch_radius * (1 - width / (2cone_height)),
         pitch_radius * (1 + width / (2cone_height)))
    end
    GraphicFrustumTrajectory(Symbol(body.name, ".gear"), :gear,
        point_a, point_b, radius_a, radius_b, graphic_appearance.color,
        graphic_appearance.opacity)
end

function spatial_frustum_trajectory(body, loaded, values, table, appearance,
        palette_index)
    graphics = table["graphics"]
    label = "$(body.name).graphics"
    marker = spatial_body_graphic_marker(body, loaded, graphics, label)
    length_value = required_positive_graphic_number(graphics, "length", label)
    radius_1 = required_positive_graphic_number(graphics, "radius_1", label)
    radius_2 = required_positive_graphic_number(graphics, "radius_2", label)
    graphic_appearance = spatial_body_graphic_appearance(
        body, appearance, palette_index)
    isnothing(graphic_appearance) && return nothing
    center, directions = spatial_frame_history(marker, values)
    point_a = center .- length_value / 2 .* directions[3]
    point_b = center .+ length_value / 2 .* directions[3]
    GraphicFrustumTrajectory(Symbol(body.name, ".frustum"), :frustum,
        point_a, point_b, radius_1, radius_2, graphic_appearance.color,
        graphic_appearance.opacity)
end

function required_surface_vertices(table, label)
    specification = get(table, "vertices", nothing)
    specification isa AbstractVector && length(specification) >= 3 ||
        throw(ArgumentError(
            "$label.vertices must contain at least three vertices"))
    vertices = zeros(Float64, length(specification), 3)
    for (vertex_number, vertex) in enumerate(specification)
        vertex isa AbstractVector && length(vertex) == 3 &&
            all(value -> value isa Number && !(value isa Bool), vertex) ||
            throw(ArgumentError(
                "$label.vertices[$vertex_number] must contain three numbers"))
        vertices[vertex_number, :] .= Float64.(vertex)
        all(isfinite, @view vertices[vertex_number, :]) ||
            throw(ArgumentError(
                "$label.vertices[$vertex_number] must be finite"))
    end
    vertices
end

function triangulate_surface_faces(specification, vertex_count, label)
    specification isa AbstractVector && !isempty(specification) ||
        throw(ArgumentError("$label.faces must contain at least one face"))
    triangles = NTuple{3,Int}[]
    edges = Set{NTuple{2,Int}}()
    for (face_number, face) in enumerate(specification)
        face isa AbstractVector && length(face) >= 3 &&
            all(index -> index isa Integer && !(index isa Bool), face) ||
            throw(ArgumentError(
                "$label.faces[$face_number] must contain at least three " *
                "integer vertex indices"))
        indices = Int.(face)
        all(index -> index in 1:vertex_count, indices) ||
            throw(ArgumentError(
                "$label.faces[$face_number] contains a vertex index " *
                "outside 1:$vertex_count"))
        length(unique(indices)) == length(indices) || throw(ArgumentError(
            "$label.faces[$face_number] must not repeat a vertex"))
        for index in 2:(length(indices) - 1)
            push!(triangles, (indices[1], indices[index], indices[index + 1]))
        end
        for index in eachindex(indices)
            endpoints = (indices[index], indices[mod1(index + 1,
                length(indices))])
            push!(edges, minmax(endpoints...))
        end
    end
    triangles, edges
end

function spatial_surface_frame(owner_name, owner_type, loaded, values,
        table, label)
    if owner_type == "rigid_body"
        return spatial_body_graphic_frame(
            loaded.bodies[owner_name], loaded, values, table, label)
    elseif owner_type == "marker"
        haskey(table, "marker") && throw(ArgumentError(
            "$label is already beneath a marker and must not specify marker"))
        haskey(loaded.markers, owner_name) || throw(ArgumentError(
            "$label does not belong to a model marker"))
        return spatial_frame_history(loaded.markers[owner_name], values)
    elseif owner_type == "ground"
        if haskey(table, "marker")
            marker = check_graphic_marker_owner(required_graphic_marker(
                loaded, table["marker"], label), owner_name, owner_type, label)
            return spatial_frame_history(marker, values)
        end
        samples = size(values, 1)
        origin = zeros(Float64, samples, 3)
        directions = ntuple(3) do axis
            values = zeros(Float64, samples, 3)
            values[:, axis] .= 1.0
            values
        end
        return origin, directions
    end
    throw(ArgumentError(
        "$label surface must be beneath a rigid body, marker, or ground"))
end

function spatial_surface_trajectories(loaded, values, element_tables,
        appearance, aliases)
    trajectories = GraphicSurfaceTrajectory{Float64}[]
    body_names = sort!(collect(keys(loaded.bodies)))
    body_colors = Dict(name => appearance.body_palette[
        mod1(index, length(appearance.body_palette))]
        for (index, name) in enumerate(body_names))
    for owner_name in sort!(collect(keys(element_tables)))
        element = element_tables[owner_name]
        owner_type = string(get(element, "type", ""))
        owner_type in ("rigid_body", "ground", "marker") || continue
        graphics = get(element, "graphics", nothing)
        graphics isa AbstractDict || continue
        parent_style = get(appearance.styles, owner_name, GraphicStyle())
        parent_style.visible || continue
        shapes = Tuple{Vector{String},Any}[]
        haskey(graphics, "shape") && push!(shapes, (String[], graphics))
        collect_graphic_shapes!(shapes, graphics)
        for (path, table) in shapes
            shape = get(table, "shape", nothing)
            shape isa AbstractString && lowercase(String(shape)) == "surface" ||
                continue
            label = isempty(path) ? "$(owner_name).graphics" :
                "$(owner_name).graphics.$(join(path, '.'))"
            require_graphic_boolean(table, "visible", true, label) || continue
            local_vertices = required_surface_vertices(table, label)
            default_color = owner_type == "rigid_body" ?
                body_colors[owner_name] : "dimgray"
            surface_color = haskey(table, "color") ?
                resolve_graphic_color(table["color"], aliases, label) :
                something(parent_style.color, default_color)
            surface_opacity = something(optional_graphic_opacity(table, label),
                parent_style.opacity, 1.0)
            patch_specifications = get(table, "patches", nothing)
            patch_tables = if isnothing(patch_specifications)
                haskey(table, "faces") || throw(ArgumentError(
                    "$label requires faces or a patches table"))
                ["surface" => table]
            else
                patch_specifications isa AbstractDict &&
                    !isempty(patch_specifications) || throw(ArgumentError(
                        "$label.patches must be a nonempty table"))
                haskey(table, "faces") && throw(ArgumentError(
                    "$label must use either faces or patches, not both"))
                [String(name) => patch_specifications[name]
                    for name in sort!(collect(keys(patch_specifications));
                        by = string)]
            end
            patches = GraphicSurfacePatch[]
            surface_edges = Set{NTuple{2,Int}}()
            for (patch_name, patch_table) in patch_tables
                patch_table isa AbstractDict || throw(ArgumentError(
                    "$label.patches.$patch_name must be a table"))
                patch_label = isnothing(patch_specifications) ? label :
                    "$label.patches.$patch_name"
                triangles, edges = triangulate_surface_faces(
                    get(patch_table, "faces", nothing),
                    size(local_vertices, 1), patch_label)
                union!(surface_edges, edges)
                color = haskey(patch_table, "color") ?
                    resolve_graphic_color(
                        patch_table["color"], aliases, patch_label) :
                    surface_color
                opacity = something(
                    optional_graphic_opacity(patch_table, patch_label),
                    surface_opacity)
                push!(patches, GraphicSurfacePatch(Symbol(patch_label),
                    triangles, color, opacity))
            end
            draw_edges = require_graphic_boolean(
                table, "draw_edges", false, label)
            edges = draw_edges ? sort!(collect(surface_edges)) :
                NTuple{2,Int}[]
            edge_color = resolve_graphic_color(
                get(table, "edge_color", "gray25"), aliases,
                "$label.edge_color")
            edge_width = optional_positive_graphic_number(
                table, "edge_width", 1.0, label)
            include_in_fit = require_graphic_boolean(
                table, "include_in_fit", true, label)
            center, directions = spatial_surface_frame(owner_name,
                owner_type, loaded, values, table, label)
            vertices = zeros(Float64, size(values, 1),
                size(local_vertices, 1), 3)
            for sample in axes(values, 1)
                orientation = hcat((@view(directions[axis][sample, :])
                    for axis in 1:3)...)
                for vertex in axes(local_vertices, 1)
                    vertices[sample, vertex, :] .=
                        @view(center[sample, :]) .+
                        orientation * @view(local_vertices[vertex, :])
                end
            end
            push!(trajectories, GraphicSurfaceTrajectory(Symbol(label),
                vertices, patches, edges, edge_color, edge_width, :graphic,
                include_in_fit))
        end
    end
    trajectories
end


function stored_spatial_mechanism_result(stored, document, element_tables,
        appearance, aliases; mode = 1, animation_samples = 121,
        animation_amplitude = 0.1, loaded_model = nothing,
        include_signals = true)
    stored.analysis_mode != :modal && mode != 1 && throw(ArgumentError(
        "mode selection is only available for modal results"))
    loaded = isnothing(loaded_model) ?
        load_spatial_model(IOBuffer(viewer_model_source(document))) :
        loaded_model
    isnothing(loaded_model) && validate_catalog(stored, loaded)
    history = viewer_history(stored; mode, animation_samples,
        animation_amplitude)
    # Output interpolation can leave Euler parameters a few parts per million
    # away from unit length.  Normalize a private display copy so every body
    # graphic remains exactly rigid.  This also gives an arbitrary but harmless
    # orientation within any plane whose two principal inertia directions are
    # equal.
    values = copy(history.values)
    times = history.times
    for body in Base.values(loaded.bodies), sample in axes(values, 1)
        parameters = @view values[sample, body.euler_parameter_variables]
        magnitude = norm(parameters)
        magnitude > eps(Float64) && (parameters ./= magnitude)
    end
    marker_histories = Dict(name => spatial_marker_history(marker, values)
        for (name, marker) in loaded.markers)
    first_marker_points = [collect(@view history[1, :])
        for history in Base.values(marker_histories)]
    marker_low = reduce((first, second) -> min.(first, second),
        first_marker_points)
    marker_high = reduce((first, second) -> max.(first, second),
        first_marker_points)
    mechanism_span = max(maximum(marker_high - marker_low), 1.0e-3)
    # Use one nominal point-constraint size throughout a spatial model. The
    # cap prevents a vehicle's overall dimensions from making its ball joints
    # much larger than the spheres inside its revolute joints.
    joint_sphere_diameter = 0.04min(mechanism_span, 1.0)
    graphic_cylinders = spatial_cylinder_trajectories(loaded, values,
        element_tables, appearance, aliases)
    graphic_markers = spatial_sphere_trajectories(loaded, values,
        element_tables, appearance, aliases)
    graphic_frustums = GraphicFrustumTrajectory{Float64}[]
    graphic_surfaces = spatial_surface_trajectories(loaded, values,
        element_tables, appearance, aliases)
    append!(graphic_surfaces, spatial_inertia_ellipsoids(
        loaded, values, marker_histories, appearance))
    for (palette_index, name) in
            enumerate(sort!(collect(keys(loaded.bodies))))
        table = element_tables[name]
        graphics = get(table, "graphics", Dict{String,Any}())
        shape = lowercase(String(get(graphics, "shape", "box")))
        if shape in ("gear", "pulley")
            shape == "pulley" && haskey(graphics, "cone_height") &&
                throw(ArgumentError(
                    "$(name).graphics pulley shape does not use cone_height"))
            trajectory = spatial_gear_trajectory(loaded.bodies[name], loaded,
                values, table, appearance, palette_index)
            isnothing(trajectory) || push!(graphic_frustums, trajectory)
        elseif shape == "frustum"
            trajectory = spatial_frustum_trajectory(loaded.bodies[name],
                loaded, values, table, appearance, palette_index)
            isnothing(trajectory) || push!(graphic_frustums, trajectory)
        elseif shape in ("surface", "sphere")
            nothing
        else
            append!(graphic_cylinders, spatial_body_trajectories(
                loaded.bodies[name], loaded, values, table, appearance,
                palette_index))
        end
    end
    xy_frames = spatial_xy_frames(loaded, values, element_tables, aliases)
    joint_marker_names = Set{Symbol}()
    joints = JointTrajectory{Float64}[]
    guides = GuideTrajectory{Float64}[]
    connectors = ConnectorTrajectory{Float64}[]
    span_measurements = SpanMeasurementTrajectory{Float64}[]
    directed_distance_measurements =
        DirectedDistanceMeasurementTrajectory{Float64}[]
    belt_spans = BeltSpanTrajectory{Float64}[]
    belt_wraps = BeltWrapTrajectory{Float64}[]
    force_arrows = ForceArrowTrajectory{Float64}[]
    torque_arrows = TorqueArrowTrajectory{Float64}[]
    rack_guide_names = Set(connection.inline.name
        for connection in Base.values(loaded.connections)
        if connection isa SpatialRackAndPinion)
    function add_spatial_inplane_arrow!(constraint)
        marker_i = constraint.geometry.marker_i
        point = marker_histories[marker_i.name]
        reaction = Matrix{Float64}(undef, length(times), 3)
        for sample in axes(values, 1)
            state = @view values[sample, :]
            reaction[sample, :] .= state[constraint.reaction_variable] .*
                inplane_normal(constraint, state)
        end
        push!(force_arrows, ForceArrowTrajectory(constraint.name, point,
            reaction, :reaction,
            marker_i isa SpatialGroundMarker))
    end
    function add_spatial_hinge_symbol!(name, hinge)
        origin, directions = spatial_frame_history(hinge.marker_j, values)
        half_length = 1.875joint_sphere_diameter
        point_a = origin .- half_length .* directions[3]
        point_b = origin .+ half_length .* directions[3]
        pin_radius = joint_sphere_diameter / 2.4
        push!(guides, GuideTrajectory(name, point_a, point_b,
            pin_radius, :hinge))
        pin_radius
    end
    function add_spatial_perp_symbol!(name, perp)
        origin = marker_histories[perp.marker_i.name]
        _, directions_i = spatial_frame_history(perp.marker_i, values)
        _, directions_j = spatial_frame_history(perp.marker_j, values)
        direction_i = directions_i[perp.axis_i]
        direction_j = directions_j[perp.axis_j]
        normal = similar(direction_i)
        for sample in axes(normal, 1)
            value = cross(@view(direction_i[sample, :]),
                @view(direction_j[sample, :]))
            magnitude = norm(value)
            normal[sample, :] .= magnitude > 0 ? value ./ magnitude :
                directions_j[3][sample, :]
        end
        push!(xy_frames, XYFrameTrajectory(name, origin, direction_i,
            direction_j, normal, 2.5joint_sphere_diameter,
            0.875joint_sphere_diameter, "gray35", 0.0, "", :perp))
        push!(joint_marker_names, perp.marker_i.name)
        push!(joint_marker_names, perp.marker_j.name)
    end
    for name in sort!(collect(keys(loaded.connections)))
        connection = loaded.connections[name]
        if connection isa SpatialBeltComponent
            for index in eachindex(connection.spans)
                incoming = connection.spans[index]
                outgoing = connection.spans[mod1(index + 1,
                    length(connection.spans))]
                pulley = incoming.pulley_2
                center = marker_histories[pulley.center_marker.name]
                _, directions = spatial_frame_history(
                    pulley.center_marker, values)
                point_1 = Matrix(values[:, incoming.point_2_variables])
                point_2 = Matrix(values[:, outgoing.point_1_variables])
                push!(belt_wraps, BeltWrapTrajectory(
                    Symbol(connection.name, ".", pulley.name), center,
                    point_1, point_2, directions[3], pulley.pitch_radius,
                    Int(sign(incoming.beta_2))))
            end
            continue
        elseif connection isa SpatialPulleyComponent
            continue
        elseif connection isa SpatialGearPair
            push!(joint_marker_names, connection.contact_marker.name)
            point = marker_histories[connection.contact_marker.name]
            force = Matrix{Float64}(undef, length(times), 3)
            for sample in axes(values, 1)
                state = @view values[sample, :]
                force[sample, :] .= state[connection.reaction_variable] .*
                    spatial_gear_contact_direction(connection, state)
            end
            # The first gear receives the action force and the second gear
            # receives its equal and opposite reaction at the same contact
            # point.  Their colors therefore also identify their directions.
            push!(force_arrows, ForceArrowTrajectory(name, point, force,
                :applied, false))
            push!(force_arrows, ForceArrowTrajectory(name, point, -force,
                :reaction, false))
            continue
        elseif connection isa SpatialRackAndPinion
            push!(joint_marker_names, connection.contact_marker.name)
            point = marker_histories[connection.contact_marker.name]
            force = Matrix{Float64}(undef, length(times), 3)
            for sample in axes(values, 1)
                state = @view values[sample, :]
                force[sample, :] .= state[connection.reaction_variable] .*
                    spatial_rack_contact_direction(connection, state)
            end
            # The rack receives the action force and the pinion receives its
            # equal and opposite reaction at the generated pitch contact.
            push!(force_arrows, ForceArrowTrajectory(name, point, force,
                :applied, false))
            push!(force_arrows, ForceArrowTrajectory(name, point, -force,
                :reaction, false))
            continue
        elseif connection isa SpatialCoordinateCoupler
            reaction = collect(values[:, connection.reaction_variable])
            for (index, (coordinate, factor)) in enumerate(zip(
                    connection.coordinates, connection.coefficients))
                arrow_name = Symbol(name, ".coordinate_", index)
                magnitude = factor .* reaction
                if coordinate isa SpatialHingeConstraint
                    position = marker_histories[coordinate.marker_i.name]
                    _, directions = spatial_frame_history(
                        coordinate.marker_j, values)
                    push!(torque_arrows, TorqueArrowTrajectory(arrow_name,
                        position, magnitude, directions[3], :reaction,
                        coordinate.marker_i isa SpatialGroundMarker))
                else
                    marker_i = coordinate.axial_geometry.marker_i
                    marker_j = coordinate.axial_geometry.marker_j
                    position = marker_histories[marker_i.name]
                    _, directions = spatial_frame_history(marker_j, values)
                    force = directions[3] .* reshape(magnitude, :, 1)
                    push!(force_arrows, ForceArrowTrajectory(arrow_name,
                        position, force, :reaction,
                        marker_i isa SpatialGroundMarker))
                end
            end
            continue
        elseif connection isa SpatialInplaneConstraint
            marker_i = connection.geometry.marker_i
            marker_j = connection.geometry.marker_j
            push!(joint_marker_names, marker_i.name)
            push!(joint_marker_names, marker_j.name)
            point = marker_histories[marker_i.name]
            push!(joints, JointTrajectory(name, point,
                joint_sphere_diameter))
            add_spatial_inplane_arrow!(connection)
            continue
        elseif connection isa SpatialInlineConstraint
            push!(joint_marker_names, connection.marker_i.name)
            push!(joint_marker_names, connection.marker_j.name)
            point = marker_histories[connection.marker_i.name]
            name in rack_guide_names ||
                push!(joints, JointTrajectory(name, point,
                    joint_sphere_diameter))
            add_spatial_inplane_arrow!(connection.inplane_x)
            add_spatial_inplane_arrow!(connection.inplane_y)
            origin, directions = spatial_frame_history(
                connection.marker_j, values)
            axial_distance = [inline_distance(connection,
                @view(values[sample, :])) for sample in axes(values, 1)]
            half_length = max(maximum(abs, axial_distance), 0.5) + 0.25
            point_a = origin .- half_length .* directions[3]
            point_b = origin .+ half_length .* directions[3]
            push!(guides, GuideTrajectory(name, point_a, point_b,
                0.0125 * half_length))
            continue
        elseif connection isa SpatialPerpConstraint
            add_spatial_perp_symbol!(name, connection)
            continue
        elseif connection isa SpatialHingeConstraint
            push!(joint_marker_names, connection.marker_i.name)
            push!(joint_marker_names, connection.marker_j.name)
            add_spatial_hinge_symbol!(name, connection)
            continue
        end
        spherical = if connection isa SpatialSphericalJoint
            connection
        elseif connection isa SpatialRevoluteJoint
            add_spatial_hinge_symbol!(name, connection.hinge)
            connection.spherical
        elseif connection isa SpatialFixedJoint
            connection.spherical
        else
            continue
        end
        push!(joint_marker_names, spherical.marker_a.name)
        push!(joint_marker_names, spherical.marker_b.name)
        point = marker_histories[spherical.marker_a.name]
        push!(joints, JointTrajectory(name, point,
            joint_sphere_diameter))
        reaction = values[:, spherical.reaction_variables]
        push!(force_arrows, ForceArrowTrajectory(name, point,
            Matrix(reaction), :reaction,
            spherical.marker_a isa SpatialGroundMarker))
    end
    for name in sort!(collect(keys(loaded.markers)))
        name in joint_marker_names && continue
        loaded.markers[name] isa SpatialFloatingMarker && continue
        origin, directions = spatial_frame_history(loaded.markers[name], values)
        push!(xy_frames, XYFrameTrajectory(name, origin, directions...,
            0.08, 0.028, "gray70", 0.12, "", :marker))
    end
    for name in sort!(collect(keys(loaded.measures)))
        measure = loaded.measures[name]
        if measure isa SpatialSpanMeasure
            point_1 = marker_histories[measure.span.marker_1.name]
            point_2 = marker_histories[measure.span.marker_2.name]
            push!(span_measurements, SpanMeasurementTrajectory(
                name, point_1, point_2))
        elseif measure isa SpatialDirectedDistanceMeasure
            geometry = measure.geometry
            measured = marker_histories[geometry.marker_i.name]
            reference = marker_histories[geometry.marker_j.name]
            _, directions = spatial_frame_history(geometry.marker_j, values)
            normal = directions[3]
            signed_distance = vec(sum(
                (measured .- reference) .* normal; dims = 2))
            projected = measured .-
                reshape(signed_distance, :, 1) .* normal
            push!(directed_distance_measurements,
                DirectedDistanceMeasurementTrajectory(name, measured,
                    projected, normal, directions[1], directions[2]))
        end
    end
    for name in sort!(collect(keys(loaded.forces)))
        component = loaded.forces[name]
        if component isa SpatialBeltSpanComponent
            point_1 = Matrix(values[:, component.point_1_variables])
            point_2 = Matrix(values[:, component.point_2_variables])
            push!(belt_spans, BeltSpanTrajectory(name, point_1, point_2))
            force = Matrix(values[:, component.force_variables])
            push!(force_arrows, ForceArrowTrajectory(name, point_1, force,
                :applied, false))
            push!(force_arrows, ForceArrowTrajectory(name, point_2, -force,
                :reaction, false))
        elseif component isa SpatialGravityComponent
            center = zeros(length(times), 3)
            center[:, :] .= values[:, component.body.position_variables]
            force = repeat(reshape(
                component.body.mass .* component.acceleration, 1, 3),
                length(times), 1)
            push!(force_arrows, ForceArrowTrajectory(name, center, force,
                :applied, false))
        elseif component isa SpatialSpanningForceComponent
            point_1 = marker_histories[component.span.marker_1.name]
            point_2 = marker_histories[component.span.marker_2.name]
            push!(connectors, ConnectorTrajectory(name, point_1, point_2,
                0.012))
            force = Matrix(values[:, component.global_force_variables])
            push!(force_arrows, ForceArrowTrajectory(name, point_1, force,
                :applied, component.span.marker_1 isa SpatialGroundMarker))
            push!(force_arrows, ForceArrowTrajectory(name, point_2, -force,
                :reaction, component.span.marker_2 isa SpatialGroundMarker))
        elseif component isa SpatialAppliedForceComponent
            point = marker_histories[component.application_marker.name]
            force = Matrix(values[:, component.global_force_variables])
            push!(force_arrows, ForceArrowTrajectory(name, point, force,
                :applied, false))
            if !isnothing(component.reaction_marker)
                reaction = marker_histories[component.reaction_marker.name]
                push!(force_arrows, ForceArrowTrajectory(name, reaction,
                    -force, :reaction, false))
            end
        elseif component isa SpatialAppliedTorqueComponent
            marker_i = component.hinge.marker_i
            marker_j = component.hinge.marker_j
            position_i = marker_histories[marker_i.name]
            position_j = marker_histories[marker_j.name]
            _, directions = spatial_frame_history(marker_j, values)
            axis = directions[3]
            magnitude = collect(values[:, component.magnitude_variable])
            push!(torque_arrows, TorqueArrowTrajectory(name, position_i,
                magnitude, axis, :applied,
                marker_i isa SpatialGroundMarker))
            push!(torque_arrows, TorqueArrowTrajectory(name, position_j,
                -magnitude, axis, :reaction,
                marker_j isa SpatialGroundMarker))
        elseif component isa SpatialRevoluteFriction
            marker_i = component.joint.hinge.marker_i
            marker_j = component.joint.hinge.marker_j
            position_i = marker_histories[marker_i.name]
            position_j = marker_histories[marker_j.name]
            _, directions = spatial_frame_history(marker_j, values)
            axis = directions[3]
            magnitude = collect(values[:, component.torque_variable])
            push!(torque_arrows, TorqueArrowTrajectory(name, position_i,
                magnitude, axis, :applied,
                marker_i isa SpatialGroundMarker))
            push!(torque_arrows, TorqueArrowTrajectory(name, position_j,
                -magnitude, axis, :reaction,
                marker_j isa SpatialGroundMarker))
        elseif component isa SpatialBushingComponent
            marker_i = component.marker_i
            marker_j = component.marker_j
            point_i = marker_histories[marker_i.name]
            point_j = marker_histories[marker_j.name]
            center = (point_i .+ point_j) ./ 2
            _, frame_j = spatial_frame_history(marker_j, values)
            symbol = spatial_bushing_symbol_size(component, mechanism_span)
            point_a = center .- 0.5 * symbol.length .* frame_j[3]
            point_b = center .+ 0.5 * symbol.length .* frame_j[3]
            push!(connectors, ConnectorTrajectory(name, point_a, point_b,
                symbol.radius, :bushing))
            force = Matrix(values[:, component.global_force_variables])
            push!(force_arrows, ForceArrowTrajectory(name, point_i, force,
                :applied, marker_i isa SpatialGroundMarker))
            push!(force_arrows, ForceArrowTrajectory(name, point_i, -force,
                :reaction, marker_j isa SpatialGroundMarker))
            torque = Matrix(values[:, component.global_torque_variables])
            magnitude, axis = spatial_torque_magnitude_axis(torque,
                frame_j[3])
            push!(torque_arrows, TorqueArrowTrajectory(name, point_i,
                magnitude, axis, :applied,
                marker_i isa SpatialGroundMarker))
            push!(torque_arrows, TorqueArrowTrajectory(name, point_i,
                -magnitude, axis, :reaction,
                marker_j isa SpatialGroundMarker))
        elseif component isa SpatialPlaneContactComponent
            style = get(appearance.styles, name, GraphicStyle())
            if style.visible && style.show_default
                sphere_center = marker_histories[component.sphere_marker.name]
                diameter = 2component.radius
                push!(graphic_markers, GraphicMarkerTrajectory(
                    Symbol(name, ".sphere"), :sphere, sphere_center,
                    zeros(length(times)), (diameter, diameter, diameter),
                    "gray10", 0.82))
                plane_center, plane_directions = spatial_frame_history(
                    component.plane_marker, values)
                plane_thickness = max(0.03component.radius,
                    1.0e-5 * mechanism_span)
                plane_normal = plane_directions[3]
                point_a = plane_center .-
                    0.5plane_thickness .* plane_normal
                point_b = plane_center .+
                    0.5plane_thickness .* plane_normal
                push!(graphic_cylinders, GraphicCylinderTrajectory(
                    Symbol(name, ".plane"), point_a, point_b,
                    3component.radius, "gray60", 0.72))
            end
            contact = zeros(length(times), 3)
            for sample in eachindex(times)
                contact[sample, :] .= spatial_plane_contact_values(component,
                    @view(values[sample, :])).contact_point
            end
            force = Matrix(values[:, component.global_force_variables])
            push!(force_arrows, ForceArrowTrajectory(name, contact, force,
                :applied, component.sphere_marker isa SpatialGroundMarker))
            push!(force_arrows, ForceArrowTrajectory(name, contact, -force,
                :reaction, component.plane_marker isa SpatialGroundMarker))
        elseif component isa SpatialSurfaceFriction
            contact = component.contact
            points = zeros(length(times), 3)
            for sample in eachindex(times)
                points[sample, :] .= spatial_plane_contact_values(contact,
                    @view(values[sample, :])).contact_point
            end
            force = Matrix(values[:, component.global_force_variables])
            push!(force_arrows, ForceArrowTrajectory(name, points, force,
                :applied, contact.sphere_marker isa SpatialGroundMarker))
            push!(force_arrows, ForceArrowTrajectory(name, points, -force,
                :reaction, contact.plane_marker isa SpatialGroundMarker))
        elseif component isa SpatialTranslationalFriction
            marker_i = component.joint.marker_i
            marker_j = component.joint.marker_j
            point = marker_histories[marker_i.name]
            force = Matrix(values[:, component.global_force_variables])
            push!(force_arrows, ForceArrowTrajectory(name, point, force,
                :applied, marker_i isa SpatialGroundMarker))
            push!(force_arrows, ForceArrowTrajectory(name, point, -force,
                :reaction, marker_j isa SpatialGroundMarker))
        elseif component isa SpatialInplaneFriction
            marker_i = component.constraint.geometry.marker_i
            marker_j = component.constraint.geometry.marker_j
            point = marker_histories[marker_i.name]
            force = Matrix(values[:, component.global_force_variables])
            push!(force_arrows, ForceArrowTrajectory(name, point, force,
                :applied, marker_i isa SpatialGroundMarker))
            push!(force_arrows, ForceArrowTrajectory(name, point, -force,
                :reaction, marker_j isa SpatialGroundMarker))
        elseif component isa SpatialTireComponent
            contact = zeros(length(times), 3)
            for sample in eachindex(times)
                contact[sample, :] .= spatial_tire_kinematics(component,
                    @view(values[sample, :])).contact_point
            end
            force = Matrix(values[:, component.global_force_variables])
            push!(force_arrows, ForceArrowTrajectory(name, contact, force,
                :applied, false))
            push!(force_arrows, ForceArrowTrajectory(name, contact, -force,
                :reaction, component.road_marker isa SpatialGroundMarker))
        end
    end
    for name in sort!(collect(keys(loaded.drivers)))
        driver = loaded.drivers[name]
        if driver isa SpatialRotationalMotionGenerator
            marker_i = driver.hinge.marker_i
            marker_j = driver.hinge.marker_j
            position_i = marker_histories[marker_i.name]
            position_j = marker_histories[marker_j.name]
            _, directions = spatial_frame_history(marker_j, values)
            axis = directions[3]
            magnitude = collect(values[:, driver.reaction_variable])
            push!(torque_arrows, TorqueArrowTrajectory(name, position_i,
                magnitude, axis, :applied,
                marker_i isa SpatialGroundMarker))
            push!(torque_arrows, TorqueArrowTrajectory(name, position_j,
                -magnitude, axis, :reaction,
                marker_j isa SpatialGroundMarker))
        elseif driver isa SpatialTranslationalMotionGenerator
            marker_i = driver.geometry.marker_i
            marker_j = driver.geometry.marker_j
            position = marker_histories[marker_i.name]
            _, directions = spatial_frame_history(marker_j, values)
            axis = directions[3]
            magnitude = collect(values[:, driver.reaction_variable])
            force = axis .* reshape(magnitude, :, 1)
            push!(force_arrows, ForceArrowTrajectory(name, position,
                force, :applied, marker_i isa SpatialGroundMarker))
            push!(force_arrows, ForceArrowTrajectory(name, position,
                -force, :reaction, marker_j isa SpatialGroundMarker))
        elseif driver isa SpatialSpanningMotionGenerator
            marker_1 = driver.span.marker_1
            marker_2 = driver.span.marker_2
            point_1 = marker_histories[marker_1.name]
            point_2 = marker_histories[marker_2.name]
            push!(connectors, ConnectorTrajectory(name, point_1, point_2,
                0.012))
            force = Matrix(values[:, driver.global_force_variables])
            push!(force_arrows, ForceArrowTrajectory(name, point_1,
                force, :applied, marker_1 isa SpatialGroundMarker))
            push!(force_arrows, ForceArrowTrajectory(name, point_2,
                -force, :reaction, marker_2 isa SpatialGroundMarker))
        end
    end
    signals = Pair{String,Vector{Float64}}[]
    if include_signals
        for column in axes(values, 2)
            prefix = history.modal ? "Δ" : ""
            signal = history.modal ?
                collect(values[:, column] .- stored.values[1, column]) :
                collect(values[:, column])
            push!(signals, viewer_signal(stored.variable_components[column],
                stored.variable_names[column], signal; prefix))
        end
        append!(signals, history.extra_signals)
    end
    MechanismResult(history.title, :spatial, times,
        BodyTrajectory{Float64}[], GearTrajectory{Float64}[],
        PulleyTrajectory{Float64}[], belt_spans,
        belt_wraps, joints,
        force_arrows, torque_arrows,
        guides, connectors,
        span_measurements, directed_distance_measurements,
        TorsionalSpringTrajectory{Float64}[], PlaneTrajectory{Float64}[],
        SphereTrajectory{Float64}[], graphic_cylinders, graphic_frustums,
        graphic_surfaces,
        graphic_markers,
        xy_frames, appearance, signals, history.bookmarks)
end

"""Prepare model information reused while converting live spatial states."""
function spatial_viewer_context(loaded)
    document = TOML.parse(loaded.model_source)
    element_tables = collect_model_tables!(Dict{Symbol,Any}(), document)
    appearance, color_aliases =
        parse_viewer_appearance(document, element_tables)
    variables = loaded.layout.catalog.variables
    (; loaded, document, element_tables, appearance, color_aliases,
       variable_names = string.(getproperty.(variables, :name)),
       variable_components = string.(getproperty.(variables, :component)))
end

"""Convert in-memory spatial states to renderer-neutral viewer data."""
function spatial_mechanism_result(context::NamedTuple,
        values::AbstractMatrix, times::AbstractVector;
        title = context.loaded.title, include_signals = true)
    size(values, 1) == length(times) || throw(DimensionMismatch(
        "spatial viewer state rows must match the time samples"))
    loaded = context.loaded
    variables = loaded.layout.catalog.variables
    size(values, 2) == length(variables) || throw(DimensionMismatch(
        "spatial viewer states must match the model variable catalog"))
    stored = (;
        title = String(title),
        times = Float64.(times),
        values = Matrix{Float64}(values),
        variable_names = context.variable_names,
        variable_components = context.variable_components,
        analysis_mode = :static,
        health_snapshots = Any[])
    stored_spatial_mechanism_result(stored, context.document,
        context.element_tables, context.appearance, context.color_aliases;
        loaded_model = loaded, include_signals)
end

spatial_mechanism_result(loaded, values::AbstractMatrix,
        times::AbstractVector; kwargs...) = spatial_mechanism_result(
    spatial_viewer_context(loaded), values, times; kwargs...)

spatial_mechanism_result(model_or_context, state::AbstractVector;
        time = 0.0, kwargs...) = spatial_mechanism_result(model_or_context,
    reshape(collect(state), 1, :), [time]; kwargs...)

"""Convert planar states into renderer-neutral viewer data."""
function planar_mechanism_result(stored, document, element_tables,
        appearance, color_aliases; loaded_model = nothing, mode = 1,
        animation_samples = 121, animation_amplitude = 0.1,
        include_signals = true)
    loaded = isnothing(loaded_model) ?
        load_planar_model(IOBuffer(viewer_model_source(document))) :
        loaded_model
    isnothing(loaded_model) && validate_catalog(stored, loaded)
    history = viewer_history(stored; mode, animation_samples,
        animation_amplitude)
    history_values, times = history.values, history.times
    marker_histories = Dict(name => marker_history(marker, history_values)
        for (name, marker) in loaded.markers)
    graphic_cylinders, graphic_markers = explicit_graphics(loaded,
        history_values, marker_histories, element_tables, appearance,
        color_aliases)
    all_bodies = Dict(name => body_trajectory(loaded, loaded.bodies[name],
        history_values, marker_histories)
        for name in sort!(collect(keys(loaded.bodies))))

    gear_geometry = Dict{Symbol,NamedTuple}()
    pulley_geometry = Dict{Symbol,PlanarPulleyComponent}()
    for connection in values(loaded.connections)
        if connection isa PlanarGearPairComponent
            internal = connection.radius_1 * connection.radius_2 > 0
            gear_geometry[connection.body_1.name] = (
                radius = abs(connection.radius_1),
                center = gear_center_point(connection, connection.joint_1,
                    connection.body_1),
                internal = internal &&
                    abs(connection.radius_1) > abs(connection.radius_2))
            gear_geometry[connection.body_2.name] = (
                radius = abs(connection.radius_2),
                center = gear_center_point(connection, connection.joint_2,
                    connection.body_2),
                internal = internal &&
                    abs(connection.radius_2) > abs(connection.radius_1))
        elseif connection isa PlanarRackAndPinionComponent
            gear_geometry[connection.pinion_body.name] = (
                radius = connection.pitch_radius,
                center = gear_center_point(connection,
                    connection.revolute_joint, connection.pinion_body),
                internal = false)
        elseif connection isa PlanarPulleyComponent
            pulley_geometry[connection.body.name] = connection
        end
    end
    gears = GearTrajectory{Float64}[]
    for name in sort!(collect(keys(gear_geometry)))
        body = all_bodies[name]
        geometry = gear_geometry[name]
        center_name = marker_name_for_point(loaded, geometry.center)
        push!(gears, GearTrajectory(name, marker_histories[center_name],
            body.angle, geometry.radius, max(0.04, 0.12 * geometry.radius),
            geometry.internal))
    end
    pulleys = PulleyTrajectory{Float64}[]
    for name in sort!(collect(keys(pulley_geometry)))
        body = all_bodies[name]
        pulley = pulley_geometry[name]
        center_name = marker_name_for_point(loaded, pulley.center_marker)
        push!(pulleys, PulleyTrajectory(pulley.name,
            marker_histories[center_name], body.angle, pulley.pitch_radius,
            max(0.04, 0.12 * pulley.pitch_radius)))
    end
    bodies = [all_bodies[name] for name in sort!(collect(keys(all_bodies)))
              if !haskey(gear_geometry, name) &&
                 !haskey(pulley_geometry, name)]

    joints = JointTrajectory{Float64}[]
    force_arrows = ForceArrowTrajectory{Float64}[]
    torque_arrows = TorqueArrowTrajectory{Float64}[]
    guides = GuideTrajectory{Float64}[]
    connectors = ConnectorTrajectory{Float64}[]
    span_measurements = SpanMeasurementTrajectory{Float64}[]
    directed_distance_measurements =
        DirectedDistanceMeasurementTrajectory{Float64}[]
    torsional_springs = TorsionalSpringTrajectory{Float64}[]
    planes = PlaneTrajectory{Float64}[]
    spheres = SphereTrajectory{Float64}[]
    belt_spans = BeltSpanTrajectory{Float64}[]
    belt_wraps = BeltWrapTrajectory{Float64}[]
    body_spans = [maximum(norm(body.point_b[row, :] - body.point_a[row, :])
        for row in axes(body.point_a, 1)) for body in bodies]
    mechanism_span = maximum([body_spans;
        2 .* [geometry.radius for geometry in values(gear_geometry)];
        2 .* [pulley.pitch_radius for pulley in values(pulley_geometry)]; 1.0])
    for name in sort!(collect(keys(loaded.measures)))
        table = element_tables[name]
        endpoints = Symbol.(table["markers"])
        push!(span_measurements, SpanMeasurementTrajectory(name,
            marker_histories[endpoints[1]], marker_histories[endpoints[2]]))
    end
    for name in sort!(collect(keys(loaded.connections)))
        connection = loaded.connections[name]
        connection isa PlanarDistanceCoordinateComponent || continue
        endpoints = Symbol.(element_tables[name]["markers"])
        measured = marker_histories[endpoints[1]]
        reference = marker_histories[endpoints[2]]
        angles = orientation_history(loaded.markers[endpoints[2]],
            history_values)
        normal = zeros(length(times), 3)
        plane_x = zeros(length(times), 3)
        plane_y = zeros(length(times), 3)
        plane_y[:, 3] .= 1.0
        for sample in eachindex(times)
            normal[sample, 1:2] .=
                [-sin(angles[sample]), cos(angles[sample])]
            plane_x[sample, 1:2] .=
                [cos(angles[sample]), sin(angles[sample])]
        end
        signed_distance = vec(sum(
            (measured .- reference) .* normal; dims = 2))
        projected = measured .-
            reshape(signed_distance, :, 1) .* normal
        push!(directed_distance_measurements,
            DirectedDistanceMeasurementTrajectory(name, measured,
                projected, normal, plane_x, plane_y))
    end
    for components in values(loaded.forces), component in components
        component isa PlanarBeltSpanComponent || continue
        point_1 = zeros(length(times), 3)
        point_2 = zeros(length(times), 3)
        point_1[:, 1:2] .= history_values[:, component.point_1_variables]
        point_2[:, 1:2] .= history_values[:, component.point_2_variables]
        push!(belt_spans, BeltSpanTrajectory(
            component.name, point_1, point_2))
    end
    for connection in values(loaded.connections)
        connection isa PlanarBeltComponent || continue
        for index in eachindex(connection.spans)
            incoming = connection.spans[index]
            outgoing = connection.spans[mod1(index + 1,
                length(connection.spans))]
            pulley = incoming.pulley_2
            center_name = marker_name_for_point(loaded, pulley.center_marker)
            point_1 = zeros(length(times), 3)
            point_2 = zeros(length(times), 3)
            point_1[:, 1:2] .= history_values[:, incoming.point_2_variables]
            point_2[:, 1:2] .= history_values[:, outgoing.point_1_variables]
            push!(belt_wraps, BeltWrapTrajectory(
                Symbol(connection.name, ".", pulley.name),
                marker_histories[center_name], point_1, point_2,
                pulley.pitch_radius, Int(sign(incoming.beta_2))))
        end
    end
    for name in sort!(collect(keys(loaded.connections)))
        connection = loaded.connections[name]
        if connection isa PlanarRevoluteJointComponent
            marker_name = marker_name_for_point(loaded, connection.marker_a)
            position = marker_histories[marker_name]
            half_length = 0.04 * mechanism_span
            point_a, point_b = copy(position), copy(position)
            point_a[:, 3] .-= half_length
            point_b[:, 3] .+= half_length
            push!(guides, GuideTrajectory(name, point_a, point_b,
                0.012 * mechanism_span, :revolute))
            push!(force_arrows, ForceArrowTrajectory(name,
                position, vector_history(history_values,
                    connection.reaction_variables), :reaction,
                isnothing(connection.body_a)))
        elseif connection isa PlanarFixedJoint
            marker_name = marker_name_for_point(
                loaded, connection.revolute.marker_a)
            push!(joints, JointTrajectory(name,
                marker_histories[marker_name]))
            push!(force_arrows, ForceArrowTrajectory(name,
                marker_histories[marker_name], vector_history(history_values,
                    connection.revolute.reaction_variables), :reaction,
                isnothing(connection.revolute.body_a)))
            push!(torque_arrows, TorqueArrowTrajectory(name,
                marker_histories[marker_name], scalar_history(history_values,
                    connection.perp.reaction_variable), :reaction,
                isnothing(connection.perp.body_i)))
        elseif connection isa PlanarGearPairComponent
            marker_name = marker_name_for_point(loaded, connection.marker_1)
            position = marker_histories[marker_name]
            force = zeros(length(times), 3)
            for sample in eachindex(times)
                state = @view(history_values[sample, :])
                direction = PlanarComponentAssembly.gear_force_direction(
                    connection, state)
                reaction = history_values[sample,
                    connection.reaction_variable]
                force[sample, 1] = reaction * direction[1]
                force[sample, 2] = reaction * direction[2]
            end
            push!(force_arrows,
                ForceArrowTrajectory(name, position, force, :reaction))
        elseif connection isa PlanarRackAndPinionComponent
            marker_name = marker_name_for_point(
                loaded, connection.contact_marker)
            push!(joints, JointTrajectory(name,
                marker_histories[marker_name]))
            position_name = marker_name_for_point(loaded,
                connection.rack_marker)
            force = zeros(length(times), 3)
            for sample in eachindex(times)
                rack = PlanarComponentAssembly.rack_translation_values(
                    connection, @view(history_values[sample, :]))
                force[sample, 1:2] .= history_values[sample,
                    connection.reaction_variable] .* rack.axis
            end
            push!(force_arrows, ForceArrowTrajectory(name,
                marker_histories[position_name], force, :reaction))
        elseif connection isa PlanarCoordinateCoupler
            reaction = scalar_history(history_values,
                connection.reaction_variable)
            for (index, (coordinate, factor)) in enumerate(zip(
                    connection.coordinates, connection.coefficients))
                arrow_name = Symbol(name, ".coordinate_", index)
                generalized = factor .* reaction
                if coordinate isa PlanarRevoluteJointComponent
                    push!(torque_arrows, TorqueArrowTrajectory(arrow_name,
                        orientation_position_history(loaded,
                            coordinate.rotation_marker_a, history_values,
                            all_bodies), generalized, :reaction,
                        isnothing(coordinate.body_a)))
                else
                    geometry = coordinate.geometry
                    marker_name = marker_name_for_point(
                        loaded, geometry.marker_i)
                    force = zeros(length(times), 3)
                    for sample in eachindex(times)
                        direction = PlanarComponentAssembly.inplane_direction(
                            coordinate, @view(history_values[sample, :]))
                        force[sample, 1:2] .=
                            generalized[sample] .* direction.unit
                    end
                    push!(force_arrows, ForceArrowTrajectory(arrow_name,
                        marker_histories[marker_name], force, :reaction,
                        isnothing(geometry.body_i)))
                end
            end
        elseif connection isa PlanarInplaneConstraint ||
                connection isa PlanarTranslationalJoint
            inplane = connection isa PlanarTranslationalJoint ?
                connection.inplane : connection
            geometry = inplane.geometry
            marker_i_name = marker_name_for_point(loaded, geometry.marker_i)
            marker_j_name = marker_name_for_point(loaded, geometry.marker_j)
            push!(joints, JointTrajectory(name,
                marker_histories[marker_i_name]))
            origin = marker_histories[marker_j_name]
            follower = marker_histories[marker_i_name]
            samples = length(times)
            point_a, point_b = zeros(samples, 3), zeros(samples, 3)
            reaction_force = zeros(samples, 3)
            half_span = max(1.5 * mechanism_span, 0.5)
            for sample in 1:samples
                direction = PlanarComponentAssembly.inplane_direction(
                    inplane, @view(history_values[sample, :]))
                # `unit` is marker j's local y axis and therefore the plane
                # normal. Its perpendicular is marker j's local x direction,
                # which spans the plane together with the out-of-plane axis.
                tangent = direction.normal
                point_a[sample, 1:2] .= origin[sample, 1:2] .-
                    half_span .* tangent
                point_b[sample, 1:2] .= origin[sample, 1:2] .+
                    half_span .* tangent
                reaction_force[sample, 1:2] .= history_values[sample,
                    inplane.reaction_variable] .* direction.unit
            end
            push!(planes, PlaneTrajectory(name, point_a, point_b,
                max(0.6 * half_span, 0.1)))
            push!(force_arrows, ForceArrowTrajectory(name, follower,
                reaction_force, :reaction, isnothing(geometry.body_i)))
            if connection isa PlanarTranslationalJoint
                push!(torque_arrows, TorqueArrowTrajectory(name, follower,
                    scalar_history(history_values,
                        connection.perp.reaction_variable), :reaction,
                    isnothing(connection.perp.body_i)))
            end
        elseif connection isa PlanarPerpConstraint
            push!(torque_arrows, TorqueArrowTrajectory(name,
                orientation_position_history(loaded,
                    connection.orientation_i, history_values, all_bodies),
                scalar_history(history_values,
                    connection.reaction_variable), :reaction,
                isnothing(connection.body_i)))
        end
    end

    for name in sort!(collect(keys(loaded.drivers)))
        driver = loaded.drivers[name]
        if driver isa PlanarRotationalMotionGenerator
            push!(torque_arrows, TorqueArrowTrajectory(name,
                orientation_position_history(loaded, driver.marker_a,
                    history_values, all_bodies),
                scalar_history(history_values, driver.torque_variable),
                :reaction, isnothing(driver.body_a)))
            continue
        end
        driver isa PlanarTranslationalMotionGenerator || continue
        geometry = driver.geometry
        marker_i_name = marker_name_for_point(loaded, geometry.marker_i)
        marker_j_name = marker_name_for_point(loaded, geometry.marker_j)
        follower = marker_histories[marker_i_name]
        origin = marker_histories[marker_j_name]
        samples = length(times)
        point_a, point_b = zeros(samples, 3), zeros(samples, 3)
        projection = zeros(samples, 3)
        tangents = Vector{Vector{Float64}}(undef, samples)
        tangent_coordinates = zeros(samples)
        for sample in 1:samples
            direction = PlanarComponentAssembly.inplane_direction(
                driver, @view(history_values[sample, :]))
            tangents[sample] = direction.normal
            separation = follower[sample, 1:2] - origin[sample, 1:2]
            tangent_coordinates[sample] = dot(separation, direction.normal)
            distance = dot(separation, direction.unit)
            projection[sample, 1:2] .=
                follower[sample, 1:2] .- distance .* direction.unit
        end
        lower, upper = extrema(tangent_coordinates)
        margin = max(0.25 * mechanism_span, 0.1)
        for sample in 1:samples
            point_a[sample, 1:2] .= origin[sample, 1:2] .+
                (lower - margin) .* tangents[sample]
            point_b[sample, 1:2] .= origin[sample, 1:2] .+
                (upper + margin) .* tangents[sample]
        end
        half_span = (upper - lower + 2margin) / 2
        push!(planes, PlaneTrajectory(name, point_a, point_b,
            max(0.6 * half_span, 0.1)))
        push!(connectors, ConnectorTrajectory(name, projection, follower,
            max(0.012, 0.012 * mechanism_span)))
        force = zeros(samples, 3)
        for sample in 1:samples
            direction = PlanarComponentAssembly.inplane_direction(
                driver, @view(history_values[sample, :]))
            force[sample, 1:2] .= history_values[sample,
                driver.reaction_variable] .* direction.unit
        end
        push!(force_arrows, ForceArrowTrajectory(name, follower,
            force, :reaction, isnothing(geometry.body_i)))
    end

    for components in values(loaded.forces)
        for component in components
            if component isa PlanarGravityComponent
                position = all_bodies[component.body.name].center
                force = zeros(length(times), 3)
                force[:, 1:2] .= transpose(
                    component.body.mass .* component.acceleration)
                push!(force_arrows, ForceArrowTrajectory(component.name,
                    position, force, :applied))
            elseif component isa PlanarAppliedForceComponent
                application_name = marker_name_for_point(
                    loaded, component.application_marker)
                samples = length(times)
                force = zeros(samples, 3)
                for sample in 1:samples
                    force[sample, 1:2] .=
                        PlanarComponentAssembly.applied_force_value(component,
                            times[sample], @view(history_values[sample, :]))
                end
                push!(force_arrows, ForceArrowTrajectory(component.name,
                    marker_histories[application_name], force, :applied))
                if !isnothing(component.reaction_marker)
                    reaction_name = marker_name_for_point(
                        loaded, component.reaction_marker)
                    push!(force_arrows, ForceArrowTrajectory(component.name,
                        marker_histories[reaction_name], -force, :reaction,
                        point_is_on_ground(component.reaction_marker)))
                end
            elseif component isa PlanarConstantTorqueComponent ||
                    component isa PlanarAppliedTorqueComponent
                torque = component isa PlanarConstantTorqueComponent ?
                    fill(Float64(component.torque), length(times)) :
                    [Float64(PlanarComponentAssembly.applied_torque_value(
                        component, time, @view(history_values[sample, :])))
                     for (sample, time) in enumerate(times)]
                point_a = component isa PlanarAppliedTorqueComponent ?
                    component.point_a : nothing
                point_b = component isa PlanarAppliedTorqueComponent ?
                    component.point_b : nothing
                position_a = isnothing(point_a) ?
                    orientation_position_history(loaded, component.marker_a,
                        history_values, all_bodies) :
                    marker_histories[marker_name_for_point(loaded, point_a)]
                position_b = isnothing(point_b) ?
                    orientation_position_history(loaded, component.marker_b,
                        history_values, all_bodies) :
                    marker_histories[marker_name_for_point(loaded, point_b)]
                push!(torque_arrows, TorqueArrowTrajectory(component.name,
                    position_a, torque, :applied,
                    orientation_is_on_ground(component.marker_a)))
                push!(torque_arrows, TorqueArrowTrajectory(component.name,
                    position_b, -torque, :reaction,
                    orientation_is_on_ground(component.marker_b)))
            elseif component isa PlanarBeltSpanComponent
                point_1 = zeros(length(times), 3)
                point_2 = zeros(length(times), 3)
                point_1[:, 1:2] .= history_values[:,
                    component.point_1_variables]
                point_2[:, 1:2] .= history_values[:,
                    component.point_2_variables]
                force = vector_history(history_values,
                    component.force_variables)
                push!(force_arrows, ForceArrowTrajectory(component.name,
                    point_1, force, :applied))
                push!(force_arrows, ForceArrowTrajectory(component.name,
                    point_2, -force, :reaction))
            elseif component isa PlanarSpanningForceComponent
                marker_1_name = marker_name_for_point(
                    loaded, component.element.marker_1)
                marker_2_name = marker_name_for_point(
                    loaded, component.element.marker_2)
                push!(connectors, ConnectorTrajectory(component.name,
                    marker_histories[marker_1_name],
                    marker_histories[marker_2_name],
                    max(0.012, 0.012 * mechanism_span)))
                force = vector_history(history_values,
                    component.element.global_force_variables)
                push!(force_arrows, ForceArrowTrajectory(component.name,
                    marker_histories[marker_1_name], force, :applied,
                    marker_is_on_ground(loaded.markers[marker_1_name])))
                push!(force_arrows, ForceArrowTrajectory(component.name,
                    marker_histories[marker_2_name], -force, :reaction,
                    marker_is_on_ground(loaded.markers[marker_2_name])))
            elseif component isa PlanarTorsionalSpringComponent
                point_1 = marker_histories[component.marker_1.name]
                point_2 = marker_histories[component.marker_2.name]
                center = (point_1 + point_2) ./ 2
                push!(torsional_springs, TorsionalSpringTrajectory(
                    component.name, center,
                    orientation_history(component.marker_1, history_values),
                    orientation_history(component.marker_2, history_values),
                    max(0.08, 0.08 * mechanism_span)))
                torque = scalar_history(history_values,
                    component.element.torque_variable)
                push!(torque_arrows, TorqueArrowTrajectory(component.name,
                    point_1, torque, :applied,
                    marker_is_on_ground(component.marker_1)))
                push!(torque_arrows, TorqueArrowTrajectory(component.name,
                    point_2, -torque, :reaction,
                    marker_is_on_ground(component.marker_2)))
            elseif component isa PlanarBushingComponent
                point_a = marker_histories[component.marker_1.name]
                point_b = marker_histories[component.marker_2.name]
                center = (point_a .+ point_b) ./ 2
                angle = orientation_history(component.marker_2,
                    history_values)
                axis = zeros(length(times), 3)
                axis[:, 1] .= cos.(angle)
                axis[:, 2] .= sin.(angle)
                symbol = planar_bushing_symbol_size(component,
                    mechanism_span)
                symbol_a = center .- 0.5 * symbol.length .* axis
                symbol_b = center .+ 0.5 * symbol.length .* axis
                push!(connectors, ConnectorTrajectory(component.name,
                    symbol_a, symbol_b, symbol.radius, :bushing))
                force = vector_history(history_values,
                    component.force_variables)
                torque = scalar_history(history_values,
                    component.torque_variable)
                push!(force_arrows, ForceArrowTrajectory(component.name,
                    point_a, force, :applied,
                    marker_is_on_ground(component.marker_1)))
                push!(force_arrows, ForceArrowTrajectory(component.name,
                    point_b, -force, :reaction,
                    marker_is_on_ground(component.marker_2)))
                push!(torque_arrows, TorqueArrowTrajectory(component.name,
                    point_a, torque, :applied,
                    marker_is_on_ground(component.marker_1)))
                push!(torque_arrows, TorqueArrowTrajectory(component.name,
                    point_b, -torque, :reaction,
                    marker_is_on_ground(component.marker_2)))
            elseif component isa PlanarPlaneContactComponent
                sphere_center = marker_histories[component.marker_1.name]
                origin = marker_histories[component.marker_2.name]
                push!(spheres, SphereTrajectory(component.name,
                    sphere_center, component.radius))
                samples = length(times)
                point_a, point_b = zeros(samples, 3), zeros(samples, 3)
                surface_point = zeros(samples, 3)
                half_span = max(1.5 * mechanism_span, 0.5)
                for sample in 1:samples
                    angle = PlanarAppliedForces.marker_angle(
                        component.marker_2.orientation,
                        @view(history_values[sample, :]))
                    tangent = [cos(angle), sin(angle)]
                    normal = [-sin(angle), cos(angle)]
                    point_a[sample, 1:2] .= origin[sample, 1:2] .-
                        half_span .* tangent
                    point_b[sample, 1:2] .= origin[sample, 1:2] .+
                        half_span .* tangent
                    surface_point[sample, 1:2] .=
                        sphere_center[sample, 1:2] .-
                        component.radius .* normal
                end
                push!(joints, JointTrajectory(component.name, surface_point))
                push!(planes, PlaneTrajectory(component.name, point_a, point_b,
                    max(0.6 * half_span, 2component.radius)))
                force = vector_history(history_values,
                    component.global_force_variables)
                second_point = marker_histories[component.marker_2.name]
                push!(force_arrows, ForceArrowTrajectory(component.name,
                    sphere_center, force, :applied,
                    marker_is_on_ground(component.marker_1)))
                push!(force_arrows, ForceArrowTrajectory(component.name,
                    second_point, -force, :reaction,
                    marker_is_on_ground(component.marker_2)))
            end
        end
    end

    signals = Pair{String,Vector{Float64}}[]
    kind_priority = Dict("position" => 1, "orientation" => 2,
        "relative_position" => 2,
        "velocity" => 3, "angular_velocity" => 4,
        "relative_velocity" => 4,
        "acceleration" => 5, "angular_acceleration" => 6,
        "relative_acceleration" => 6,
        "reaction" => 7, "applied_geometry" => 8,
        "applied_rate" => 9, "applied_load" => 10)
    if include_signals
        columns = sort(collect(axes(history_values, 2)); by = column ->
            (get(kind_priority, stored.variable_kinds[column], 99), column))
        for column in columns
            prefix = history.modal ? "Δ" : ""
            signal = history.modal ?
                collect(history_values[:, column] .- stored.values[1, column]) :
                collect(history_values[:, column])
            push!(signals, viewer_signal(stored.variable_components[column],
                stored.variable_names[column], signal; prefix))
        end
        append!(signals, history.extra_signals)
    end
    dimension = Symbol(document["model"]["dimension"])
    MechanismResult(history.title, dimension, times,
        bodies, gears, pulleys,
        belt_spans, belt_wraps, joints, force_arrows, torque_arrows,
        guides, connectors,
        span_measurements, directed_distance_measurements,
        torsional_springs, planes, spheres, graphic_cylinders,
        GraphicFrustumTrajectory{Float64}[],
        GraphicSurfaceTrajectory{Float64}[], graphic_markers,
        XYFrameTrajectory{Float64}[], appearance, signals,
        history.bookmarks)
end

"""Prepare planar model information reused while converting in-memory states."""
function planar_viewer_context(loaded)
    document = TOML.parse(loaded.model_source)
    element_tables = collect_model_tables!(Dict{Symbol,Any}(), document)
    appearance, color_aliases =
        parse_viewer_appearance(document, element_tables)
    variables = loaded.layout.catalog.variables
    (; loaded, document, element_tables, appearance, color_aliases,
       variable_names = string.(getproperty.(variables, :name)),
       variable_components = string.(getproperty.(variables, :component)),
       variable_kinds = string.(getproperty.(variables, :kind)))
end

"""Convert in-memory planar states to renderer-neutral viewer data."""
function planar_mechanism_result(context::NamedTuple,
        values::AbstractMatrix, times::AbstractVector;
        title = context.loaded.title, include_signals = true)
    size(values, 1) == length(times) || throw(DimensionMismatch(
        "planar viewer state rows must match the time samples"))
    loaded = context.loaded
    variables = loaded.layout.catalog.variables
    size(values, 2) == length(variables) || throw(DimensionMismatch(
        "planar viewer states must match the model variable catalog"))
    stored = (;
        title = String(title),
        times = Float64.(times),
        values = Matrix{Float64}(values),
        variable_names = context.variable_names,
        variable_components = context.variable_components,
        variable_kinds = context.variable_kinds,
        analysis_mode = :dynamic,
        health_snapshots = Any[])
    planar_mechanism_result(stored, context.document,
        context.element_tables, context.appearance, context.color_aliases;
        loaded_model = loaded, include_signals)
end

planar_mechanism_result(loaded, values::AbstractMatrix,
        times::AbstractVector; kwargs...) = planar_mechanism_result(
    planar_viewer_context(loaded), values, times; kwargs...)

planar_mechanism_result(model_or_context, state::AbstractVector;
        time = 0.0, kwargs...) = planar_mechanism_result(model_or_context,
    reshape(collect(state), 1, :), [time]; kwargs...)

"""
    load_viewer_model(path)

Load a planar or spatial TOML model, or a spatial Lua assembly model, for
previewing in the common viewer.
"""
function load_viewer_model(path::AbstractString)
    extension = lowercase(splitext(path)[2])
    if extension == ".lua"
        return load_spatial_model(path)
    elseif extension == ".toml"
        document = TOML.parsefile(path)
        model = get(document, "model", nothing)
        model isa AbstractDict || throw(ArgumentError(
            "model preview requires a [model] table"))
        dimension = String(get(model, "dimension", ""))
        return dimension == "planar" ? load_planar_model(path) :
            dimension == "spatial" ? load_spatial_model(path) :
            throw(ArgumentError(
                "model.dimension must be 'planar' or 'spatial'"))
    else
        throw(ArgumentError(
            "model preview requires a .toml or .lua file"))
    end
end

"""
    model_mechanism_result(model; configuration = :entered)

Build a one-sample viewer result from a loaded model. `:entered` preserves the
configuration supplied by the model; `:consistent` shows the values after
initial-condition correction.
"""
function model_mechanism_result(loaded::Union{
        LoadedPlanarModel,LoadedSpatialModel}; configuration = :entered)
    configuration in (:entered, :consistent) || throw(ArgumentError(
        "model preview configuration must be :entered or :consistent"))
    state = configuration == :entered ? loaded.entered_initial_values :
        loaded.initial_values
    label = configuration == :entered ? "model input" :
        "consistent initial conditions"
    title = loaded.title * " — " * label
    time = loaded.simulation.start_time
    if loaded isa LoadedSpatialModel
        return spatial_mechanism_result(loaded, state; time, title,
            include_signals = false)
    end
    planar_mechanism_result(loaded, state; time, title,
        include_signals = false)
end

function model_mechanism_result(path::AbstractString;
        configuration = :entered)
    model_mechanism_result(load_viewer_model(path); configuration)
end

function in_memory_static_snapshots(result)
    hasproperty(result, :static_progress_events) ||
        return StoredStaticSnapshot[]
    [StoredStaticSnapshot(Symbol(event.phase), Symbol(event.status),
         Float64(event.model_time),
         isnothing(event.pseudo_time) ? nothing : Float64(event.pseudo_time),
         Int(event.iteration), Int(event.relaxation_cycle), copy(event.state),
         Float64(event.force_imbalance), Float64(event.torque_imbalance),
         Float64(event.constraint_error),
         Float64(event.equivalent_acceleration),
         Float64(event.reciprocal_condition), Bool(event.mass_regularized),
         String(event.message))
     for event in result.static_progress_events]
end

"""Normalize an in-memory solver result for the common viewer conversion."""
function in_memory_viewer_storage(result)
    variables = result.loaded.layout.catalog.variables
    modal = result.analysis_mode == :modal
    health = hasproperty(result, :solution) ?
        PracticalMechanicalSimulation.ResultIO.health_snapshots(result) :
        StoredHealthSnapshot[]
    (; title = result.loaded.title,
       times = Float64.(result.times),
       values = PracticalMechanicalSimulation.ResultIO.history_matrix(
           result.states),
       variable_names = string.(getproperty.(variables, :name)),
       variable_components = string.(getproperty.(variables, :component)),
       variable_kinds = string.(getproperty.(variables, :kind)),
       analysis_mode = result.analysis_mode,
       health_snapshots = health,
       static_snapshots = in_memory_static_snapshots(result),
       modal_eigenvalues = modal ? ComplexF64.(result.eigenvalues) :
           ComplexF64[],
       modal_natural_frequencies_hz = modal ?
           Float64.(result.natural_frequencies_hz) : Float64[],
       modal_damped_frequencies_hz = modal ?
           Float64.(result.damped_frequencies_hz) : Float64[],
       modal_damping_ratios = modal ? Float64.(result.damping_ratios) :
           Float64[],
       modal_mode_shapes = modal ? ComplexF64.(result.mode_shapes) :
           Matrix{ComplexF64}(undef, 0, 0),
       modal_equation_errors = modal ? Float64.(result.equation_errors) :
           Float64[],
       modal_shift = modal ? ComplexF64(result.shift) : 0.0 + 0.0im)
end

"""Convert a completed in-memory analysis into the common viewer data model."""
function simulation_mechanism_result(result; mode = 1,
        animation_samples = 121, animation_amplitude = 0.1)
    stored = in_memory_viewer_storage(result)
    context = result.loaded isa LoadedSpatialModel ?
        spatial_viewer_context(result.loaded) :
        planar_viewer_context(result.loaded)
    if result.loaded isa LoadedSpatialModel
        return stored_spatial_mechanism_result(stored, context.document,
            context.element_tables, context.appearance,
            context.color_aliases; mode, animation_samples,
            animation_amplitude, loaded_model = result.loaded)
    end
    planar_mechanism_result(stored, context.document,
        context.element_tables, context.appearance,
        context.color_aliases; mode, animation_samples,
        animation_amplitude, loaded_model = result.loaded)
end

"""Convert a stored planar or spatial run into the viewer data model."""
function stored_mechanism_result(path::AbstractString; mode = 1,
        animation_samples = 121, animation_amplitude = 0.1)
    stored = read_result(path)
    isempty(stored.model_source) &&
        throw(ArgumentError("result does not contain an embedded TOML model"))
    document = TOML.parse(stored.model_source)
    element_tables = collect_model_tables!(Dict{Symbol,Any}(), document)
    appearance, color_aliases =
        parse_viewer_appearance(document, element_tables)
    dimension = String(document["model"]["dimension"])
    dimension == "spatial" && return stored_spatial_mechanism_result(
        stored, document, element_tables, appearance, color_aliases; mode,
        animation_samples, animation_amplitude)
    dimension == "planar" || throw(ArgumentError(
        "stored model dimension must be 'planar' or 'spatial'"))
    planar_mechanism_result(stored, document, element_tables, appearance,
        color_aliases; mode, animation_samples, animation_amplitude)
end

end
