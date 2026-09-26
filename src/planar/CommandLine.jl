module CommandLine

using ..SimulationRunner
using ..ResultIO

export planar_model_main

function print_usage(io)
    println(io, "usage: simp2d [MODEL.toml|MODEL.lua|-] [duration] [samples]",
        " [--output RESULT.simp] [--overwrite]")
end

function parse_command_line(args)
    positional = String[]
    output_path = nothing
    overwrite = false
    index = 1
    while index <= length(args)
        argument = args[index]
        if argument == "--output" || argument == "-o"
            index == length(args) && throw(ArgumentError("$argument requires a filename"))
            output_path = args[index + 1]
            index += 2
        elseif argument == "--overwrite"
            overwrite = true
            index += 1
        else
            push!(positional, argument)
            index += 1
        end
    end
    (; positional, output_path, overwrite)
end

"""
    planar_model_main(args=ARGS; input=stdin, output=stdout, error=stderr)

Command-line entry point for planar TOML or Lua analysis. A missing model name
or `-` reads TOML from `input`. Optional positional duration and sample count
override the model; `--output` stores HDF5 and `--overwrite` permits replacement.
Returns zero on success and two for a usage error.
"""
function planar_model_main(args = ARGS; input = stdin, output = stdout, error = stderr)
    options = parse_command_line(args)
    length(options.positional) <= 3 || begin
        print_usage(error)
        return 2
    end
    positional = options.positional
    source_name = isempty(positional) ? "-" : positional[1]
    duration = length(positional) >= 2 ? parse(Float64, positional[2]) : nothing
    samples = length(positional) >= 3 ? parse(Int, positional[3]) : nothing
    source = source_name == "-" ? input : source_name
    writer = Ref{Union{Nothing,IncrementalResultWriter}}(nothing)
    failure_event = Ref{Any}(nothing)
    result_progress = if isnothing(options.output_path)
        nothing
    else
        function (event)
            if event.kind == :begin
                writer[] = begin_incremental_result(options.output_path,
                    event.loaded, event.analysis_mode;
                    overwrite = options.overwrite,
                    initial_time = first(event.times))
            else
                isnothing(writer[]) || record_incremental_event!(writer[], event)
                event.kind == :failure && (failure_event[] = event)
            end
        end
    end
    result = try
        run_planar_model(source; duration, samples, result_progress)
    catch exception
        if !isnothing(writer[])
            failure_time = isnothing(failure_event[]) ?
                writer[].last_sample_time : failure_event[].time
            statistics = isnothing(failure_event[]) ? nothing :
                (hasproperty(failure_event[], :statistics) ?
                    failure_event[].statistics : nothing)
            status = exception isa InterruptException ? :interrupted : :failed
            message = if !isnothing(failure_event[]) &&
                    hasproperty(failure_event[], :message)
                failure_event[].message
            else
                sprint(showerror, exception)
            end
            finish_incremental_result!(writer[], status;
                message, failure_time, statistics)
        end
        rethrow()
    end
    if !isnothing(options.output_path)
        isnothing(writer[]) || finish_incremental_result!(writer[], :complete)
        write_result(options.output_path, result; overwrite = true)
    end
    println(output, result.loaded.title)
    state_word = result.loaded.analysis.degrees_of_freedom == 1 ? "state" : "states"
    println(output, "  analysis: ", result.analysis_mode,
        " (", result.loaded.analysis.degrees_of_freedom, " ", state_word, ")")
    println(output, "  variables: ", length(result.loaded.layout.catalog.variables))
    println(output, "  equations: ", length(result.loaded.layout.catalog.equations))
    if result.analysis_mode == :modal
        println(output, "  modes: ", length(result.eigenvalues))
        for index in eachindex(result.eigenvalues)
            println(output, "    ", index, ": ",
                round(result.natural_frequencies_hz[index]; sigdigits = 7),
                " Hz, damping ratio ",
                round(result.damping_ratios[index]; sigdigits = 5))
        end
        println(output, "  sparse modal factorizations: ",
            result.sparse_factorizations)
        println(output, "  maximum modal equation error: ",
            maximum(result.equation_errors))
    else
        println(output, "  samples: ", length(result.times))
    end
    println(output, "  initial consistency corrections: ",
        result.initial_consistency_iterations)
    if result.loaded.initial_conditions.enabled
        saved = result.loaded.initial_conditions
        velocity_note = saved.include_velocities ?
            "configuration and velocities" : "configuration"
        source_note = saved.requested_sample == "static" ?
            "static initialization, sample $(saved.sample)" :
            "sample $(saved.sample)"
        println(output, "  saved initial conditions: ", velocity_note,
            " from ", source_note, " (t=", saved.source_time,") of ",
            saved.result_path)
    end
    if hasproperty(result, :static_initialization_iterations) &&
            !isnothing(result.static_initialization_iterations)
        println(output, "  static equilibrium initialization corrections: ",
            result.static_initialization_iterations)
    end
    if hasproperty(result, :static_relaxation_cycles) &&
            !isnothing(result.static_relaxation_cycles) &&
            result.static_relaxation_cycles isa Integer &&
            result.static_relaxation_cycles > 0
        println(output, "  dynamic relaxation cycles: ",
            result.static_relaxation_cycles)
    end
    if hasproperty(result, :solution)
        println(output, "  accepted/rejected integration steps: ",
            result.solution.stats.accepted_steps, "/",
            result.solution.stats.rejected_steps)
        if hasproperty(result, :physical_error_peak) &&
                !isnothing(result.physical_error_peak)
            peak = result.physical_error_peak
            canonical_index = result.loaded.active_variable_indices[
                peak.dominant_local_index]
            variable = result.loaded.layout.catalog.variables[canonical_index]
            println(output, "  maximum physical predictor error: ",
                round(peak.error; sigdigits = 6), " in ",
                variable.component, ".", variable.name,
                " at t=", round(peak.time; sigdigits = 8))
        else
            monitor = result.solution.error_monitor
            isnothing(monitor) || begin
                finite_steps = findall(isfinite, monitor.maximum_errors)
                if !isempty(finite_steps)
                    peak_step = finite_steps[argmax(
                        monitor.maximum_errors[finite_steps])]
                    local_index = monitor.maximum_error_indices[peak_step]
                    canonical_index = result.loaded.active_variable_indices[
                        local_index]
                    variable = result.loaded.layout.catalog.variables[
                        canonical_index]
                    println(output, "  maximum physical predictor error: ",
                        round(monitor.maximum_errors[peak_step]; sigdigits = 6),
                        " in ", variable.component, ".", variable.name,
                        " at t=", round(result.solution.t[peak_step];
                            sigdigits = 8))
                end
            end
        end
        result.solution.stats.events_found > 0 &&
            println(output, "  located force-law transitions: ",
                result.solution.stats.events_found)
        result.solution.stats.state_reselections > 0 &&
            println(output, "  runtime state reselections: ",
                result.solution.stats.state_reselections)
    elseif hasproperty(result, :static_iterations)
        println(output, "  static equilibrium solutions: ",
            result.continuation_solutions)
        println(output, "  static Newton iterations: ",
            sum(result.static_iterations))
        hasproperty(result, :static_relaxation_cycles) &&
            sum(result.static_relaxation_cycles) > 0 &&
            println(output, "  dynamic relaxation cycles: ",
                sum(result.static_relaxation_cycles))
    end
    !isnothing(options.output_path) &&
        println(output, "  result file: ", options.output_path)
    return 0
end

end
