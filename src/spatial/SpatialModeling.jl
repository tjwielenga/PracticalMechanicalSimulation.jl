"""Geometry and allocation helpers for the spatial modeler."""
module SpatialModeling

using LinearAlgebra
using StaticArrays: SMatrix, SVector
using ..SpatialComponentAssembly

export SpatialGroundMarker, SpatialBodyMarker, SpatialFlexibleBeamMarker,
       SpatialFloatingMarker,
       spatial_marker_position, spatial_marker_orientation,
       spatial_marker_velocity, spatial_marker_acceleration,
       spatial_marker_dependency_indices, is_flexible_marker,
       spatial_timoshenko_matrices, allocated_spatial_flexible_beam,
       spatial_flexible_beam_marker, allocated_spatial_body

"""Oriented marker fixed in the global frame."""
struct SpatialGroundMarker{T}
    name::Symbol
    position::Vector{T}
    orientation::Matrix{T}
end

"""
Oriented marker fixed in a body's reference frame.

Both `position_body` and `orientation_body` are expressed from that reference
frame, not necessarily from the body's center of mass. Body-reference to CM
offsets are handled when the model is loaded.
"""
struct SpatialBodyMarker{B,T}
    name::Symbol
    body::B
    position_body::Vector{T}
    orientation_body::Matrix{T}
end

"""Generated end or center marker on a floating-reference flexible beam."""
struct SpatialFlexibleBeamMarker{B,T}
    name::Symbol
    body::B
    position_body::Vector{T}
    orientation_body::Matrix{T}
    translation_shape::Matrix{T}
    orientation_shape::Matrix{T}
end

"""
Point owned by one body while following another marker's global point.

Its position, velocity, and acceleration follow `follower`, but applied force
and moment contributions go to `body`. Gear, rack, and belt elements use this
to place a contact load on a moving body without adding a kinematic constraint.
"""
struct SpatialFloatingMarker{B,F}
    name::Symbol
    body::B
    follower::F
end

is_flexible_marker(::Any) = false
is_flexible_marker(::SpatialFlexibleBeamMarker) = true

spatial_marker_dependency_indices(::SpatialGroundMarker) = Int[]
function spatial_marker_dependency_indices(marker::SpatialBodyMarker)
    body = marker.body
    [collect(body.acceleration_variables);
     collect(body.angular_acceleration_variables);
     collect(body.velocity_variables);
     collect(body.angular_velocity_variables);
     collect(body.position_variables);
     collect(body.euler_parameter_variables)]
end
function spatial_marker_dependency_indices(marker::SpatialFlexibleBeamMarker)
    body = marker.body
    [collect(body.acceleration_variables);
     collect(body.angular_acceleration_variables);
     collect(body.elastic_acceleration_variables);
     collect(body.velocity_variables);
     collect(body.angular_velocity_variables);
     collect(body.elastic_velocity_variables);
     collect(body.position_variables);
     collect(body.euler_parameter_variables);
     collect(body.elastic_position_variables)]
end
spatial_marker_dependency_indices(marker::SpatialFloatingMarker) =
    spatial_marker_dependency_indices(marker.follower)

function spatial_marker_position(marker::SpatialGroundMarker, z)
    marker.position
end

function spatial_marker_position(marker::SpatialBodyMarker, z)
    parameters = SVector{4}(@view z[marker.body.euler_parameter_variables])
    position = SVector{3}(@view z[marker.body.position_variables])
    position + rotation_matrix(parameters) * SVector{3}(marker.position_body)
end


function spatial_marker_position(marker::SpatialFlexibleBeamMarker, z)
    body = marker.body
    parameters = @view z[body.euler_parameter_variables]
    elastic = @view z[body.elastic_position_variables]
    offset = marker.position_body + marker.translation_shape * elastic
    z[body.position_variables] .+ rotation_matrix(parameters) * offset
end

spatial_marker_position(marker::SpatialFloatingMarker, z) =
    spatial_marker_position(marker.follower, z)

spatial_marker_velocity(::SpatialGroundMarker, z) =
    zero(SVector{3,eltype(z)})

function spatial_marker_velocity(marker::SpatialBodyMarker, z)
    body = marker.body
    parameters = SVector{4}(@view z[body.euler_parameter_variables])
    omega = SVector{3}(@view z[body.angular_velocity_variables])
    velocity = SVector{3}(@view z[body.velocity_variables])
    velocity + rotation_matrix(parameters) *
        cross(omega, SVector{3}(marker.position_body))
end


function spatial_marker_velocity(marker::SpatialFlexibleBeamMarker, z)
    body = marker.body
    parameters = @view z[body.euler_parameter_variables]
    omega = @view z[body.angular_velocity_variables]
    elastic = @view z[body.elastic_position_variables]
    elastic_rate = @view z[body.elastic_velocity_variables]
    offset = marker.position_body + marker.translation_shape * elastic
    local_rate = marker.translation_shape * elastic_rate
    z[body.velocity_variables] .+ rotation_matrix(parameters) *
        (cross(omega, offset) + local_rate)
end

spatial_marker_velocity(marker::SpatialFloatingMarker, z) =
    spatial_marker_velocity(marker.follower, z)

spatial_marker_acceleration(::SpatialGroundMarker, z) =
    zero(SVector{3,eltype(z)})

function spatial_marker_acceleration(marker::SpatialBodyMarker, z)
    body = marker.body
    parameters = SVector{4}(@view z[body.euler_parameter_variables])
    omega = SVector{3}(@view z[body.angular_velocity_variables])
    alpha = SVector{3}(@view z[body.angular_acceleration_variables])
    offset = SVector{3}(marker.position_body)
    local_acceleration = cross(alpha, offset) .+
        cross(omega, cross(omega, offset))
    SVector{3}(@view z[body.acceleration_variables]) .+
        rotation_matrix(parameters) * local_acceleration
end


function spatial_marker_acceleration(marker::SpatialFlexibleBeamMarker, z)
    body = marker.body
    parameters = @view z[body.euler_parameter_variables]
    omega = @view z[body.angular_velocity_variables]
    alpha = @view z[body.angular_acceleration_variables]
    elastic = @view z[body.elastic_position_variables]
    elastic_rate = @view z[body.elastic_velocity_variables]
    elastic_acceleration = @view z[body.elastic_acceleration_variables]
    offset = marker.position_body + marker.translation_shape * elastic
    local_rate = marker.translation_shape * elastic_rate
    local_acceleration = cross(alpha, offset) +
        cross(omega, cross(omega, offset)) +
        2 .* cross(omega, local_rate) +
        marker.translation_shape * elastic_acceleration
    z[body.acceleration_variables] .+
        rotation_matrix(parameters) * local_acceleration
end

spatial_marker_acceleration(marker::SpatialFloatingMarker, z) =
    spatial_marker_acceleration(marker.follower, z)

spatial_marker_orientation(marker::SpatialGroundMarker, z) = marker.orientation

function spatial_marker_orientation(marker::SpatialBodyMarker, z)
    parameters = SVector{4}(@view z[marker.body.euler_parameter_variables])
    rotation_matrix(parameters) * SMatrix{3,3}(marker.orientation_body)
end


function spatial_marker_orientation(marker::SpatialFlexibleBeamMarker, z)
    body = marker.body
    parameters = @view z[body.euler_parameter_variables]
    elastic = @view z[body.elastic_position_variables]
    elastic_rotation = marker.orientation_shape * elastic
    rotation_matrix(parameters) *
        (Matrix{eltype(z)}(I, 3, 3) + skew(elastic_rotation)) *
        marker.orientation_body
end

function spatial_marker_orientation(marker::SpatialFloatingMarker, z)
    parameters = SVector{4}(@view z[marker.body.euler_parameter_variables])
    rotation_matrix(parameters)
end


"""Return full and reduced matrices for a straight spatial Timoshenko beam."""
function spatial_timoshenko_matrices(mass, length, area, elastic_modulus,
        shear_modulus, second_moment_y, second_moment_z, torsion_constant,
        shear_coefficient_y, shear_coefficient_z)
    L = Float64(length)
    M = zeros(Float64, 12, 12)
    K = zeros(Float64, 12, 12)

    function add_two_node!(matrix, indices, block)
        matrix[indices, indices] .+= block
    end

    axial = elastic_modulus * area / L
    add_two_node!(K, [1, 7], axial .* [1.0 -1.0; -1.0 1.0])
    add_two_node!(M, [1, 7], mass / 6 .* [2.0 1.0; 1.0 2.0])

    torsion = shear_modulus * torsion_constant / L
    add_two_node!(K, [4, 10], torsion .* [1.0 -1.0; -1.0 1.0])
    rotary_mass = mass * torsion_constant / (6area)
    add_two_node!(M, [4, 10], rotary_mass .* [2.0 1.0; 1.0 2.0])

    function bending_blocks(second_moment, shear_coefficient)
        EI = elastic_modulus * second_moment
        phi = 12EI / (shear_coefficient * shear_modulus * area * L^2)
        k1 = 12EI / (L^3 * (1 + phi))
        k2 = 6EI / (L^2 * (1 + phi))
        k3 = (4 + phi) * EI / (L * (1 + phi))
        k4 = (2 - phi) * EI / (L * (1 + phi))
        stiffness = [k1 k2 -k1 k2;
                     k2 k3 -k2 k4;
                     -k1 -k2 k1 -k2;
                     k2 k4 -k2 k3]
        inertia = mass / 420 .* [
            156 22L 54 -13L;
            22L 4L^2 13L -3L^2;
            54 13L 156 -22L;
            -13L -3L^2 -22L 4L^2]
        stiffness, inertia
    end

    # Local v bends about z. Local w bends about y and therefore uses the
    # opposite nodal-rotation sign under the right-handed convention.
    K_v, M_v = bending_blocks(second_moment_z, shear_coefficient_y)
    add_two_node!(K, [2, 6, 8, 12], K_v)
    add_two_node!(M, [2, 6, 8, 12], M_v)
    K_w, M_w = bending_blocks(second_moment_y, shear_coefficient_z)
    sign_change = Diagonal([1.0, -1.0, 1.0, -1.0])
    add_two_node!(K, [3, 5, 9, 11], sign_change * K_w * sign_change)
    add_two_node!(M, [3, 5, 9, 11], sign_change * M_w * sign_change)

    rigid = zeros(Float64, 12, 6)
    for (node, x) in ((1, -L / 2), (2, L / 2))
        translation = node == 1 ? (1:3) : (7:9)
        rotation = node == 1 ? (4:6) : (10:12)
        offset = [x, 0.0, 0.0]
        rigid[translation, 1:3] .= Matrix{Float64}(I, 3, 3)
        rigid[translation, 4:6] .= -skew(offset)
        rigid[rotation, 4:6] .= Matrix{Float64}(I, 3, 3)
    end
    seed = zeros(Float64, 12, 6)
    seed[1:6, :] .= -0.5 .* Matrix{Float64}(I, 6, 6)
    seed[7:12, :] .= 0.5 .* Matrix{Float64}(I, 6, 6)
    shape = seed - rigid * ((transpose(rigid) * M * rigid) \
        (transpose(rigid) * M * seed))
    elastic_mass_raw = transpose(shape) * M * shape
    elastic_stiffness_raw = transpose(shape) * K * shape
    elastic_mass = (elastic_mass_raw + transpose(elastic_mass_raw)) / 2
    elastic_stiffness =
        (elastic_stiffness_raw + transpose(elastic_stiffness_raw)) / 2
    rigid_mass = transpose(rigid) * M * rigid
    (; shape, elastic_mass, elastic_stiffness, rigid_mass)
end

"""Allocate a six-mode floating-reference spatial Timoshenko beam."""
function allocated_spatial_flexible_beam(layout, name, mass, length, area,
        elastic_modulus, shear_modulus, second_moment_y, second_moment_z,
        torsion_constant, shear_coefficient_y, shear_coefficient_z,
        damping_time_scale; inertia = nothing)
    matrices = spatial_timoshenko_matrices(mass, length, area,
        elastic_modulus, shear_modulus, second_moment_y, second_moment_z,
        torsion_constant, shear_coefficient_y, shear_coefficient_z)
    default_inertia = Matrix(matrices.rigid_mass[4:6, 4:6])
    body_inertia = isnothing(inertia) ? default_inertia : Matrix(inertia)
    SpatialComponentAssembly.allocated_spatial_flexible_beam(layout, name,
        mass, body_inertia, length, matrices.elastic_mass,
        matrices.elastic_stiffness, matrices.shape, damping_time_scale)
end

"""Construct the generated `end_i`, `cm`, or `end_j` beam marker."""
function spatial_flexible_beam_marker(beam, node::Symbol)
    node in (:end_i, :cm, :end_j) || throw(ArgumentError(
        "flexible beam marker must be end_i, cm, or end_j"))
    if node == :cm
        position = zeros(Float64, 3)
        translation_shape = zeros(Float64, 3, 6)
        orientation_shape = zeros(Float64, 3, 6)
    else
        rows = node == :end_i ? (1:6) : (7:12)
        position = [node == :end_i ? -beam.length / 2 : beam.length / 2,
                    0.0, 0.0]
        translation_shape = Matrix(beam.deformation_shape[rows[1:3], :])
        orientation_shape = Matrix(beam.deformation_shape[rows[4:6], :])
    end
    SpatialFlexibleBeamMarker(Symbol(beam.name, :., node), beam, position,
        Matrix{Float64}(I, 3, 3), translation_shape, orientation_shape)
end

end
