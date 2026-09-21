#!/usr/bin/env julia

using PracticalMechanicalSimulation

const Sim2D = PracticalMechanicalSimulation.Sim2D

"""Add one rigid planar link and return its useful handles."""
function pendulum_link!(model, name;
        length, mass, center, color, width = 0.08)
    body = Sim2D.rigid_body!(model, name;
        mass,
        inertia = mass * (length^2 + width^2) / 12,
        position = center)
    inner = Sim2D.marker!(body, :inner;
        position = [-length / 2, 0.0])
    outer = Sim2D.marker!(body, :outer;
        position = [length / 2, 0.0])
    Sim2D.graphics!(body; show_default = false, color)
    Sim2D.graphic!(body, :member;
        shape = :cylinder,
        markers = [inner, outer],
        radius = width / 2)
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
        center = [first_length / 2, 0.0],
        color = :steelblue)
    second = pendulum_link!(model, "$prefix.second_link";
        length = second_length,
        mass = second_mass,
        center = [first_length + second_length / 2, 0.0],
        color = :darkorange)

    base_joint = Sim2D.revolute!(model, "$prefix.base_joint";
        markers = [first.inner, ground_pin],
        rotation_coordinates = true)
    elbow_joint = Sim2D.revolute!(model, "$prefix.elbow_joint";
        markers = [second.inner, first.outer],
        rotation_coordinates = true)
    Sim2D.gravity!(model, "$prefix.gravity";
        acceleration = [0.0, -9.81],
        bodies = [first.body, second.body])

    (; bodies = [first.body, second.body], base_joint, elbow_joint,
       first, second)
end


model = Sim2D.Model(:hierarchical_double_pendulum;
    title = "Hierarchical Sim2D Julia API double pendulum")
Sim2D.analysis!(model; mode = :automatic)
Sim2D.simulation!(model;
    end_time = 4.0,
    frames_per_second = 60,
    relative_tolerance = 1.0e-7,
    absolute_tolerance = 1.0e-9,
    maximum_step = 0.01)
Sim2D.graphics!(model;
    background = :white,
    body_palette = :colorblind)

ground = Sim2D.ground!(model, :ground)
ground_pin = Sim2D.marker!(ground, :pin)

pendulum = double_pendulum!(model, :pendulum, ground_pin)
Sim2D.state_selection!(model;
    method = :preferred,
    preferred_velocities = [
        Sim2D.variable(pendulum.base_joint, :omega),
        Sim2D.variable(pendulum.elbow_joint, :omega),
    ],
    allow_fallback = false)

root = normpath(joinpath(@__DIR__, "..", ".."))
output = isempty(ARGS) ?
    joinpath(root, "results", "examples", "planar",
        "julia-api-double-pendulum.simp") : abspath(ARGS[1])
mkpath(dirname(output))

result = Sim2D.run(model)
Sim2D.save_result(output, result; overwrite = true)

println("Wrote $output")
println("Open it with bin/simpview-web")
