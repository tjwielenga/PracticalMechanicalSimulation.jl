#!/usr/bin/env julia

using PracticalMechanicalSimulation

const Sim3D = PracticalMechanicalSimulation.Sim3D

"""Add one rigid pendulum link and return its useful handles."""
function pendulum_link!(model, name;
        length, mass, center, color, width = 0.08)
    transverse_inertia = mass * (length^2 + width^2) / 12
    axial_inertia = mass * width^2 / 6
    body = Sim3D.rigid_body!(model, name;
        mass,
        inertia = [axial_inertia, transverse_inertia, transverse_inertia],
        position = center)
    Sim3D.graphics!(body;
        shape = :box,
        size = [length, width, width],
        color)
    inner = Sim3D.marker!(body, :inner;
        position = [-length / 2, 0.0, 0.0])
    outer = Sim3D.marker!(body, :outer;
        position = [length / 2, 0.0, 0.0])
    (; body, inner, outer)
end

"""Build a complete double-pendulum assembly below one hierarchical name."""
function double_pendulum!(model, name, ground_pin;
        first_length = 1.0,
        second_length = 0.8,
        first_mass = 1.0,
        second_mass = 0.7)
    prefix = string(name)
    first = pendulum_link!(model, "$prefix.first_link";
        length = first_length,
        mass = first_mass,
        center = [first_length / 2, 0.0, 0.0],
        color = :steelblue)
    second = pendulum_link!(model, "$prefix.second_link";
        length = second_length,
        mass = second_mass,
        center = [first_length + second_length / 2, 0.0, 0.0],
        color = :darkorange)

    base_joint = Sim3D.revolute!(model, "$prefix.base_joint";
        markers = [first.inner, ground_pin],
        rotation_coordinates = true)
    elbow_joint = Sim3D.revolute!(model, "$prefix.elbow_joint";
        markers = [second.inner, first.outer],
        rotation_coordinates = true)
    Sim3D.gravity!(model, "$prefix.gravity";
        acceleration = [0.0, -9.81, 0.0],
        bodies = [first.body, second.body])

    (; bodies = [first.body, second.body], base_joint, elbow_joint,
       first, second)
end

model = Sim3D.Model(:hierarchical_double_pendulum;
    title = "Hierarchical Julia API double pendulum")
Sim3D.analysis!(model; mode = :automatic)
Sim3D.simulation!(model;
    end_time = 4.0,
    frames_per_second = 60,
    relative_tolerance = 1.0e-7,
    absolute_tolerance = 1.0e-9,
    maximum_step = 0.01)
Sim3D.graphics!(model;
    background = :white,
    body_palette = :colorblind)

ground = Sim3D.ground!(model, :ground)
ground_pin = Sim3D.marker!(ground, :pin)
Sim3D.graphics!(ground_pin;
    shape = :xy_frame,
    axis_length = 0.35,
    plane_size = 0.16,
    plane_color = :gray65,
    label = "ground.pin")

pendulum = double_pendulum!(model, :pendulum, ground_pin)
Sim3D.state_selection!(model;
    method = :preferred,
    preferred_velocities = [
        Sim3D.variable(pendulum.base_joint, :omega),
        Sim3D.variable(pendulum.elbow_joint, :omega),
    ],
    allow_fallback = false)

root = normpath(joinpath(@__DIR__, "..", ".."))
output = isempty(ARGS) ?
    joinpath(root, "results", "examples", "spatial",
        "julia-api-double-pendulum.simp") : abspath(ARGS[1])
mkpath(dirname(output))

result = Sim3D.run(model)
Sim3D.save_result(output, result; overwrite = true)

println("Wrote $output")
println("Open it with bin/simpView")
