#!/usr/bin/env julia

using LinearAlgebra

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUTPUT_DIRECTORY = joinpath(@__DIR__, "0.2.0")

include(joinpath(ROOT, "examples", "planar", "implicit_displacement_pendulum.jl"))
include(joinpath(ROOT, "examples", "planar", "implicit_acceleration_pendulum.jl"))
include(joinpath(ROOT, "examples", "planar", "implicit_velocity_pendulum.jl"))
include(joinpath(ROOT, "examples", "planar", "implicit_baumgarte_pendulum.jl"))
include(joinpath(ROOT, "examples", "planar", "gear_constraint_satisfaction_pendulum.jl"))
include(joinpath(ROOT, "examples", "planar", "independent_state_pendulum.jl"))

function csv_field(value)
    text = value isa AbstractFloat ? repr(value) : string(value)
    occursin(r"[\",\n\r]", text) ?
        "\"" * replace(text, "\"" => "\"\"") * "\"" : text
end

function write_csv(path, rows)
    columns = collect(keys(first(rows)))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((csv_field(getproperty(row, column))
                for column in columns), ','))
        end
    end
end

function single_level_row(name, solution, diagnostics, comparison)
    (; formulation = name,
       accepted_steps = solution.stats.accepted_steps,
       rejected_steps = solution.stats.rejected_steps,
       maximum_position_constraint_error = diagnostics.maximum_position_constraint_error,
       maximum_velocity_constraint_error = diagnostics.maximum_velocity_constraint_error,
       maximum_acceleration_constraint_error = diagnostics.maximum_acceleration_constraint_error,
       maximum_state_difference = comparison.maximum_state_difference,
       maximum_reaction_difference = comparison.maximum_reaction_difference,
       maximum_energy_error = diagnostics.maximum_energy_error)
end

function independent_constraint_errors(solution, p)
    position = 0.0
    velocity = 0.0
    acceleration = 0.0
    for z in solution.u
        r_global, d_global = marker_vectors(z[9], p)
        position = max(position, norm(z[7:8] + r_global - p.pin, Inf))
        velocity = max(velocity, norm(z[4:5] + d_global .* z[6], Inf))
        acceleration = max(acceleration,
            norm(z[1:2] + d_global .* z[3] - r_global .* z[6]^2, Inf))
    end
    (; position, velocity, acceleration)
end

function multiple_level_row(name, variables, controlled, levels, solution;
        position, velocity, acceleration, state, reaction, energy)
    (; formulation = name,
       simultaneous_variables = variables,
       error_controlled_variables = controlled,
       retained_constraint_levels = levels,
       accepted_steps = solution.stats.accepted_steps,
       rejected_steps = solution.stats.rejected_steps,
       maximum_position_constraint_error = position,
       maximum_velocity_constraint_error = velocity,
       maximum_acceleration_constraint_error = acceleration,
       maximum_state_difference = state,
       maximum_reaction_difference = reaction,
       maximum_energy_error = energy)
end

function main()
    mkpath(OUTPUT_DIRECTORY)
    p = PendulumParameters()

    acceleration_solution = run_acceleration_historical_ddassl()
    acceleration_diagnostics = acceleration_solution_diagnostics(acceleration_solution, p)
    acceleration_comparison = compare_implicit_with_reduced(acceleration_solution, p)

    velocity_solution = run_velocity_historical_ddassl()
    velocity_diagnostics = velocity_solution_diagnostics(velocity_solution, p)
    velocity_comparison = compare_implicit_with_reduced(velocity_solution, p)

    position_solution = run_displacement_historical_ddassl()
    position_diagnostics = implicit_solution_diagnostics(position_solution, p)
    position_comparison = compare_implicit_with_reduced(position_solution, p)

    single_level = [
        single_level_row("Acceleration, Deficit Zero", acceleration_solution,
            acceleration_diagnostics, acceleration_comparison),
        single_level_row("Velocity, Deficit One", velocity_solution,
            velocity_diagnostics, velocity_comparison),
        single_level_row("Position, Deficit Two", position_solution,
            position_diagnostics, position_comparison),
    ]
    write_csv(joinpath(OUTPUT_DIRECTORY, "formulation_single_level.csv"),
        single_level)

    baumgarte = NamedTuple[]
    for correction_time in (0.20, 0.10, 0.05)
        parameters = BaumgarteParameters(p, correction_time)
        solution = run_baumgarte_historical_ddassl(; correction_time)
        diagnostics = baumgarte_diagnostics(solution, parameters)
        push!(baumgarte, (;
            correction_time_s = correction_time,
            accepted_steps = solution.stats.accepted_steps,
            rejected_steps = solution.stats.rejected_steps,
            final_position_constraint_error = diagnostics.final_position_error,
            final_velocity_constraint_error = diagnostics.final_velocity_error,
            maximum_acceleration_constraint_error =
                diagnostics.maximum_acceleration_error,
            maximum_critical_decay_difference = diagnostics.maximum_critical_decay_difference,
            maximum_order = maximum(solution.orders)))
    end
    write_csv(joinpath(OUTPUT_DIRECTORY, "formulation_baumgarte.csv"),
        baumgarte)

    gear_v_solution = run_gear_constraint_satisfaction_pendulum()
    gear_v = gear_constraint_satisfaction_diagnostics(gear_v_solution, p)
    gear_a_solution = run_complete_gear_constraint_satisfaction_pendulum()
    gear_a = complete_gear_constraint_satisfaction_diagnostics(gear_a_solution, p)
    independent_solution = run_independent_state_pendulum()
    independent = independent_state_diagnostics(independent_solution, p)
    constraints = independent_constraint_errors(independent_solution, p)

    multiple_level = [
        multiple_level_row("Velocity constraint", 8, 6, "velocity",
            velocity_solution;
            position = velocity_diagnostics.maximum_position_constraint_error,
            velocity = velocity_diagnostics.maximum_velocity_constraint_error,
            acceleration = velocity_diagnostics.maximum_acceleration_constraint_error,
            state = velocity_comparison.maximum_state_difference,
            reaction = velocity_comparison.maximum_reaction_difference,
            energy = velocity_diagnostics.maximum_energy_error),
        multiple_level_row("GearStableV", 13, 6, "position; velocity",
            gear_v_solution;
            position = gear_v.maximum_position_constraint_error,
            velocity = gear_v.maximum_velocity_constraint_error,
            acceleration = gear_v.maximum_acceleration_constraint_error,
            state = gear_v.maximum_state_difference,
            reaction = gear_v.maximum_reaction_difference,
            energy = gear_v.maximum_energy_error),
        multiple_level_row("GearStableA", 15, 6,
            "position; velocity; acceleration", gear_a_solution;
            position = gear_a.maximum_position_constraint_error,
            velocity = gear_a.maximum_velocity_constraint_error,
            acceleration = gear_a.maximum_acceleration_constraint_error,
            state = gear_a.maximum_state_difference,
            reaction = gear_a.maximum_reaction_difference,
            energy = gear_a.maximum_energy_error),
        multiple_level_row("Fully Consistent", 11, 2,
            "position; velocity; acceleration", independent_solution;
            position = constraints.position,
            velocity = constraints.velocity,
            acceleration = constraints.acceleration,
            state = independent.maximum_integrated_state_difference,
            reaction = independent.maximum_reaction_difference,
            energy = independent.maximum_energy_error),
    ]
    write_csv(joinpath(OUTPUT_DIRECTORY, "formulation_multiple_level.csv"),
        multiple_level)
    println("Historical formulation evidence written to $OUTPUT_DIRECTORY")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
