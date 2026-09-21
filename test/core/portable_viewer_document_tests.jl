module PortableViewerDocumentTests

using JSON
using HDF5
using LinearAlgebra
using Test

include("../../src/viewer/ViewerData.jl")
include("../../src/viewer/PortableViewerDocument.jl")

using .ViewerData
using .PortableViewerDocument

@testset "portable viewer document" begin
    @test PortableViewerDocument.graphic_path(
        "Forces", Symbol("van.spring"), "Reactions") ==
        ["Forces", "Reactions", "van", "spring"]
    @test PortableViewerDocument.graphic_item_path(
        "Frames", Symbol("van.chassis.cm")) ==
        ["Frames", "van", "chassis", "cm"]

    times = [0.0, 0.5]
    centers = [0.0 0.0 0.0; 1.0 2.0 3.0]
    body = BodyTrajectory(:link, centers, centers, centers, zeros(2),
        0.1, (0.2, 0.3, 0.4))
    result = MechanismResult("Portable test", times, [body],
        Matrix{Float64}[], ["link.R_x" => [0.0, 1.0]])

    document = viewer_document(result)
    @test PortableViewerDocument.characteristic_graphic_length(result) > 0.4
    @test document["format"] == "SimpView"
    @test document["version"] == 4
    @test document["dimension"] == "planar"
    @test document["choices"][1]["scene"]["encoding"] == "mesh_instances"
    @test document["choices"][1]["scene"]["follow_targets"] ==
        [Dict("name" => "link", "track" => "body_frame:link")]
    body_track = only(filter(track -> track["id"] == "body_frame:link",
        document["choices"][1]["scene"]["tracks"]))
    @test body_track["position"] ==
        [[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]]
    @test length(document["choices"][1]["scene"]["tracks"]) == 1
    @test all(instance["track"] == "body_frame:link" for instance in
        document["choices"][1]["scene"]["instances"])
    @test !haskey(document["choices"][1]["scene"], "legacy")
    @test any(instance -> instance["name"] == "link",
        document["choices"][1]["scene"]["instances"])

    temporary = tempname() * ".simpview.json"
    try
        @test write_viewer_document(temporary, result) == abspath(temporary)
        parsed = JSON.parsefile(temporary)
        @test parsed["choices"][1]["signals"][1]["name"] == "link.R_x"
        @test parsed["choices"][1]["signals"][1]["values"] == [0.0, 1.0]
    finally
        isfile(temporary) && rm(temporary)
    end

    native = tempname() * ".simp"
    try
        h5open(native, "w") do file
            file["sentinel"] = [1]
        end
        @test write_graphics(native, result) == abspath(native)
        h5open(native, "r") do file
            @test haskey(file, "graphics/choices/1/meshes")
            @test haskey(file, "graphics/choices/1/tracks")
            @test haskey(file, "graphics/choices/1/instances")
            @test haskey(file,
                "graphics/choices/1/tracks/position_offset")
            @test !haskey(file, "graphics/choices/1/tracks/1")
            @test eltype(read(file[
                "graphics/choices/1/tracks/position"])) == Float32
            @test eltype(read(file[
                "graphics/choices/1/instances/local_position"])) == Float32
            @test !haskey(file, "viewer")
            @test read(attributes(file["graphics"])["format"]) ==
                "SimpGraphics"
            @test read(attributes(file["graphics"])["format_version"]) == 1
        end
        first_size = filesize(native)
        chmod(native, 0o644)
        @test write_graphics(native, result) == abspath(native)
        @test filesize(native) <= first_size
        @test stat(native).mode & 0o777 == 0o644
        h5open(native, "r") do file
            @test read(file["sentinel"]) == [1]
        end
    finally
        isfile(native) && rm(native)
    end

    vertices = zeros(2, 3, 3)
    vertices[1, :, :] .= [0.0 0.0 0.0; 1.0 0.0 0.0; 0.0 1.0 0.0]
    rotation = [0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0]
    for vertex in axes(vertices, 2)
        vertices[2, vertex, :] .= [2.0, 3.0, 4.0] +
            rotation * @view(vertices[1, vertex, :])
    end
    surface = GraphicSurfaceTrajectory(:panel, vertices,
        [GraphicSurfacePatch(:face, [(1, 2, 3)], "steelblue", 1.0)],
        [(1, 2), (2, 3), (1, 3)], "gray25", 1.0)
    compact = PortableViewerDocument.portable_value(surface)
    @test compact["encoding"] == "rigid_pose"
    @test compact["include_in_fit"] == true
    @test size(reduce(hcat, compact["vertices"])) == (3, 3)
    @test all(isapprox(norm(quaternion), 1.0) for quaternion in
        compact["quaternion"])

    background_vertices = 1000 .* vertices
    background_surface = GraphicSurfaceTrajectory(:terrain,
        background_vertices,
        [GraphicSurfacePatch(:face, [(1, 2, 3)], "gray25", 1.0)],
        Tuple{Int,Int}[], "gray25", 1.0, :graphic, false)
    @test PortableViewerDocument.portable_value(background_surface)[
        "include_in_fit"] == false
    characteristic_length = PortableViewerDocument.characteristic_graphic_length(
        result)
    push!(result.graphic_surfaces, background_surface)
    @test PortableViewerDocument.characteristic_graphic_length(result) ==
        characteristic_length

    cylinder_center = zeros(2, 3)
    push!(result.graphic_cylinders, GraphicCylinderTrajectory(:contact_plane,
        cylinder_center, cylinder_center, 0.3, "gray60", 1.0))
    cylinder_length = PortableViewerDocument.characteristic_graphic_length(
        result)
    @test cylinder_length > 1.0

    sphere_center = [0.0 0.0 1.0; 0.0 0.0 1.0]
    push!(result.graphic_markers, GraphicMarkerTrajectory(:contact_sphere,
        :sphere, sphere_center, zeros(2), (0.2, 0.2, 0.2), "gray10", 1.0))
    @test PortableViewerDocument.characteristic_graphic_length(result) >
        cylinder_length
end

end
