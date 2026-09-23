#!/usr/bin/env julia

using PracticalMechanicalSimulation

const Sim3D = PracticalMechanicalSimulation.Sim3D

model = Sim3D.Model(:julia_api_pendulum;
    title = "Pendulum built with the Julia Sim3D API")
Sim3D.analysis!(model; mode = :automatic)
Sim3D.simulation!(model;
    end_time = 3.0,
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
    label = "revolute axis")

pendulum = Sim3D.rigid_body!(model, :pendulum;
    mass = 1.0,
    inertia = [0.001, 1 / 12, 1 / 12],
    position = [0.5, 0.0, 0.0],
    velocity = [0.0, 0.2, 0.0],
    angular_velocity = [0.0, 0.0, 0.4])
Sim3D.graphics!(pendulum;
    shape = :box,
    size = [1.0, 0.08, 0.08],
    color = :steelblue)
body_pin = Sim3D.marker!(pendulum, :pin;
    position = [-0.5, 0.0, 0.0])

pin = Sim3D.revolute!(model, :pin;
    markers = [body_pin, ground_pin],
    rotation_coordinates = true)
Sim3D.gravity!(model, :gravity;
    acceleration = [0.0, -9.81, 0.0],
    bodies = [pendulum])
Sim3D.state_selection!(model;
    method = :preferred,
    preferred_velocities = [Sim3D.variable(pin, :omega)],
    allow_fallback = false)

root = normpath(joinpath(@__DIR__, "..", ".."))
output = isempty(ARGS) ?
    joinpath(root, "results", "examples", "spatial",
        "julia-api-pendulum.simp") : abspath(ARGS[1])
mkpath(dirname(output))

result = Sim3D.run(model)
Sim3D.save_result(output, result; overwrite = true)

println("Wrote $output")
println("Open it with bin/simpView")
