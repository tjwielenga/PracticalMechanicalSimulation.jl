module SpatialCommandLine

using ..SpatialSimulationRunner
using ..ResultIO

export spatial_model_main

function print_usage(io)
    println(io, "usage: run_spatial_model.jl [MODEL.toml|MODEL.lua|-] [duration] [samples]",
        " [--output RESULT.simp] [--overwrite]")
end

function parse_command_line(args)
    positional = String[]
    output_path = nothing
    overwrite = false
    index = 1
    while index <= length(args)
        argument = args[index]
        if argument in ("--output", "-o")
            index == length(args) && throw(ArgumentError(
                "$argument requires a filename"))
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

function spatial_model_main(args = ARGS; input = stdin, output = stdout,
        error = stderr, static_progress = nothing)
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
        run_spatial_model(source; duration, samples, static_progress,
            result_progress)
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
    println(output, "  analysis: ", result.analysis_mode, " (",
        result.loaded.analysis.degrees_of_freedom, " mechanical DOF)")
    println(output, "  variables: ",
        length(result.loaded.layout.catalog.variables))
    println(output, "  equations: ",
        length(result.loaded.layout.catalog.equations))
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
    if hasproperty(result, :static_initialization_iterations) &&
            !isnothing(result.static_initialization_iterations)
        println(output, "  static equilibrium initialization corrections: ",
            result.static_initialization_iterations)
        result.static_relaxation_cycles > 0 &&
            println(output, "  dynamic relaxation cycles: ",
                result.static_relaxation_cycles)
        hasproperty(result, :static_mass_regularized) &&
            result.static_mass_regularized === true &&
            println(output, "  mass-regularized static Jacobian: yes")
    end
    if hasproperty(result, :solution)
        println(output, "  accepted/rejected integration steps: ",
            result.solution.stats.accepted_steps, "/",
            result.solution.stats.rejected_steps)
    elseif hasproperty(result, :static_iterations)
        println(output, "  static equilibrium solutions: ",
            result.continuation_solutions)
        println(output, "  static Newton iterations: ",
            sum(result.static_iterations))
        sum(result.static_relaxation_cycles) > 0 &&
            println(output, "  dynamic relaxation cycles: ",
                sum(result.static_relaxation_cycles))
        hasproperty(result, :static_mass_regularized) &&
            any(result.static_mass_regularized) &&
            println(output, "  mass-regularized static solutions: ",
                count(identity, result.static_mass_regularized))
    end
    !isnothing(options.output_path) &&
        println(output, "  result file: ", options.output_path)
    0
end

end
