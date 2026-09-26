#!/usr/bin/env julia

using Printf
using PracticalMechanicalSimulation

const Sim3D = PracticalMechanicalSimulation.Sim3D
const SpatialSimulationRunner =
    PracticalMechanicalSimulation.SpatialSimulationRunner
const SpatialConstraints = PracticalMechanicalSimulation.SpatialConstraints

const TOTAL_LENGTH = 4.0
const CHORD = 0.20
const THICKNESS = 0.020
const DENSITY = 1000.0
const ELASTIC_MODULUS = 2.0e10
const POISSON_RATIO = 0.30
const DAMPING_TIME_SCALE = 0.03
const RAMP_TIME = 4.0
const END_TIME = 6.0
const MAXIMUM_SPEED = 7.5

segment_name(index) = "blade.segment_" * lpad(index, 2, '0')

"""Smooth angle history whose speed rises from zero to `target_speed`."""
function spinup_angle_expression()
    u = "min(max(t/ramp_time,0),1)"
    "target_speed*(t-ramp_time*($u)+" *
        "ramp_time*(2.5*($u)^4-3*($u)^5+($u)^6))"
end

"""Build one segmented floating-reference blade model."""
function rotating_blade_model(segment_count;
        target_speed = MAXIMUM_SPEED,
        end_time = END_TIME,
        frames_per_second = 60)
    segment_count >= 1 || throw(ArgumentError(
        "segment_count must be positive"))
    segment_length = TOTAL_LENGTH / segment_count

    model = Sim3D.Model(:rotating_flexible_blade;
        title = "Spin stiffening of a segmented flexible blade")
    Sim3D.parameters!(model;
        target_speed,
        ramp_time = RAMP_TIME)
    Sim3D.analysis!(model;
        mode = :dynamic,
        initialization = :static_equilibrium,
        static_method = :newton,
        static_tolerance = 1.0e-7,
        modes = 12,
        frequency_shift_hz = 0.0)
    Sim3D.simulation!(model;
        start_time = 0.0,
        end_time,
        frames_per_second,
        relative_tolerance = 1.0e-6,
        absolute_tolerance = 1.0e-8,
        initial_step = 1.0e-7,
        maximum_step = 0.01)
    Sim3D.graphics!(model;
        background = :white,
        body_palette = :colorblind,
        characteristic_length = TOTAL_LENGTH,
        characteristic_mass = DENSITY * CHORD * THICKNESS * TOTAL_LENGTH,
        characteristic_velocity = MAXIMUM_SPEED * TOTAL_LENGTH)

    ground = Sim3D.ground!(model, :ground)
    hub = Sim3D.marker!(ground, :hub)
    Sim3D.graphics!(hub;
        shape = :xy_frame,
        axis_length = 0.45,
        plane_size = 0.20,
        plane_color = :gray65,
        label = "rotor axis")

    beams = Sim3D.ElementRef[]
    for index in 1:segment_count
        name = segment_name(index)
        beam = Sim3D.flexible_beam!(model, name;
            length = segment_length,
            density = DENSITY,
            elastic_modulus = ELASTIC_MODULUS,
            poisson_ratio = POISSON_RATIO,
            damping_time_scale = DAMPING_TIME_SCALE,
            section = (shape = :rectangular,
                width = CHORD,
                height = THICKNESS),
            position = [(index - 0.5) * segment_length, 0.0, 0.0])
        push!(beams, beam)
        index == 1 && continue
        Sim3D.fixed!(model, "blade.connection_$(lpad(index - 1, 2, '0'))";
            markers = ["$name.end_i", "$(segment_name(index - 1)).end_j"])
    end

    root_joint = Sim3D.revolute!(model, "blade.root_joint";
        markers = ["$(segment_name(1)).end_i", hub],
        rotation_coordinates = true)
    Sim3D.rotational_motion!(model, "blade.drive";
        joint = root_joint,
        angle = spinup_angle_expression())
    Sim3D.gravity!(model, "blade.gravity";
        acceleration = [0.0, 0.0, -9.80665],
        bodies = beams)

    if segment_count == 4
        # The four-segment demonstration has a compact, durable body-local
        # partition. Longer discretizations deliberately retain automatic QR
        # selection so the condition monitor can adapt the partition as the
        # blade turns.
        preferred_velocities = Sim3D.VariableRef[]
        elastic_rates = (:eta_u_dot, :eta_v_dot, :eta_w_dot,
            :eta_rx_dot, :eta_ry_dot, :eta_rz_dot)
        for beam in beams[2:end], rate in elastic_rates
            push!(preferred_velocities, Sim3D.variable(beam, rate))
        end
        append!(preferred_velocities, [
            Sim3D.variable(first(beams), :eta_v_dot),
            Sim3D.variable(first(beams), :eta_w_dot),
            Sim3D.variable(last(beams), :V_x),
            Sim3D.variable(last(beams), :V_y),
            Sim3D.variable(last(beams), :V_z),
            Sim3D.variable(last(beams), :omega_x),
        ])
        Sim3D.state_selection!(model;
            method = :preferred,
            preferred_velocities,
            allow_fallback = false)
    end
    model
end

function variable_index(loaded, component, name)
    matches = [variable.index for variable in loaded.layout.catalog.variables
        if variable.component == Symbol(component) &&
           variable.name == Symbol(name)]
    length(matches) == 1 || error(
        "expected one variable named '$component.$name', " *
        "found $(length(matches))")
    only(matches)
end

function blade_tip_height(loaded, state, segment_count)
    marker = loaded.markers[
        Symbol(segment_name(segment_count) * ".end_j")]
    SpatialConstraints.spatial_marker_position(marker, state)[3]
end

"""Run one spin-up case and linearize its final rotating state."""
function run_case(segment_count, target_speed; save_dynamic = false,
        output_path = nothing, report_progress = false,
        end_time = END_TIME)
    model = rotating_blade_model(segment_count; target_speed, end_time)
    loaded = Sim3D.load_model(model)
    last_reported_second = Ref(-1)
    progress = if report_progress
        function (event)
            if event.kind == :sample
                second = floor(Int, event.time)
                if second > last_reported_second[]
                    last_reported_second[] = second
                    println("  reached t = ", round(event.time; digits = 3), " s")
                end
            elseif event.kind == :failure
                println("  integration stopped at t = ",
                    round(event.time; digits = 6), " s")
            end
        end
    else
        nothing
    end
    dynamic = PracticalMechanicalSimulation.run_spatial_model(loaded;
        result_progress = progress)
    time = last(dynamic.times)
    state = last(dynamic.states)
    modal = SpatialSimulationRunner.run_spatial_modal_model(
        loaded, time; initial_state = state)

    if save_dynamic
        isnothing(output_path) && error(
            "output_path is required when save_dynamic is true")
        Sim3D.save_result(output_path, dynamic; overwrite = true)
    end

    omega_index = variable_index(loaded, "blade.root_joint", "omega")
    (; dynamic, modal,
       speed = state[omega_index],
       initial_tip_z = blade_tip_height(
           loaded, first(dynamic.states), segment_count),
       final_tip_z = blade_tip_height(loaded, state, segment_count),
       first_frequency = first(modal.natural_frequencies_hz),
       equation_error = maximum(modal.equation_errors))
end

root = normpath(joinpath(@__DIR__, "..", ".."))
if !isempty(ARGS) && first(ARGS) == "--single"
    length(ARGS) in (3, 4, 5) || error(
        "usage: rotating_flexible_blade.jl --single SEGMENTS SPEED " *
        "[OUTPUT [END_TIME]]")
    count = parse(Int, ARGS[2])
    speed = parse(Float64, ARGS[3])
    single_output = length(ARGS) >= 4 ? abspath(ARGS[4]) :
        joinpath(root, "results", "scratch", "rotating-flexible-blade.simp")
    single_end_time = length(ARGS) == 5 ? parse(Float64, ARGS[5]) : END_TIME
    mkpath(dirname(single_output))
    result = run_case(count, speed;
        save_dynamic = true, output_path = single_output,
        report_progress = true, end_time = single_end_time)
    println("target speed: ", speed, " rad/s")
    println("actual speed: ", result.speed, " rad/s")
    println("static/final tip z: ", result.initial_tip_z, " / ",
        result.final_tip_z, " m")
    println("modal frequencies (Hz): ",
        join(round.(result.modal.natural_frequencies_hz; digits = 6), ", "))
    println("modal damping ratios: ",
        join(round.(result.modal.damping_ratios; digits = 6), ", "))
    println("modal equation errors: ",
        join(map(value -> @sprintf("%.3e", value),
            result.modal.equation_errors), ", "))
    println("Wrote $single_output")
    exit()
end

output = isempty(ARGS) ?
    joinpath(root, "results", "examples", "spatial",
        "rotating-flexible-blade.simp") : abspath(ARGS[1])
output_stem, output_extension = splitext(basename(output))
eight_segment_output = joinpath(dirname(output),
    output_stem * "-eight-segment" * output_extension)
mkpath(dirname(output))

speeds = collect(range(0.0, MAXIMUM_SPEED; length = 4))
speed_results = NamedTuple[]
println("Four-segment blade speed sweep")
println("| target rad/s | actual rad/s | tip z static m | tip z rotating m | " *
    "first frequency Hz | max modal error |")
println("|--:|--:|--:|--:|--:|--:|")
for speed in speeds
    result = run_case(4, speed;
        save_dynamic = speed == MAXIMUM_SPEED,
        output_path = output)
    push!(speed_results, result)
    @printf("| %.3f | %.6f | %.6f | %.6f | %.6f | %.3e |\n",
        speed, result.speed, result.initial_tip_z, result.final_tip_z,
        result.first_frequency, result.equation_error)
end

println()
println("Maximum-speed segment comparison")
println("| segments | tip z static m | tip z rotating m | " *
    "first frequency Hz | max modal error |")
println("|--:|--:|--:|--:|--:|")
eight_segment_result = run_case(8, MAXIMUM_SPEED;
    save_dynamic = true, output_path = eight_segment_output)
for (count, result) in ((4, last(speed_results)),
        (8, eight_segment_result))
    @printf("| %d | %.6f | %.6f | %.6f | %.3e |\n",
        count, result.initial_tip_z, result.final_tip_z,
        result.first_frequency, result.equation_error)
end

println()
println("Wrote $output")
println("Wrote $eight_segment_output")
println("Open either result with bin/simpView")
