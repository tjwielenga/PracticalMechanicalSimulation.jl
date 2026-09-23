module CSVExport

using ..ResultIO

export export_result_csv, export_modal_csv, export_result_main

qualified_variable_names(result::StoredSimulationResult) =
    ["$(result.variable_components[index]).$(result.variable_names[index])"
     for index in eachindex(result.variable_names)]

function selected_columns(result, requested)
    names = qualified_variable_names(result)
    isnothing(requested) && return collect(eachindex(names)), names
    columns = Int[]
    for requested_name in requested
        requested_name == "time" &&
            throw(ArgumentError("time is included automatically"))
        column = findfirst(==(requested_name), names)
        isnothing(column) &&
            throw(ArgumentError("unknown result variable '$requested_name'"))
        push!(columns, column)
    end
    columns, names[columns]
end

function csv_field(value::AbstractString)
    if occursin(r"[\",\r\n]", value)
        return "\"$(replace(value, "\"" => "\"\""))\""
    end
    value
end

function check_modal_dimensions(result::StoredSimulationResult)
    modes = length(result.modal_eigenvalues)
    modes > 0 || throw(ArgumentError("modal result contains no modes"))
    length(result.modal_natural_frequencies_hz) == modes ||
        throw(DimensionMismatch(
            "modal natural-frequency count does not match the eigenvalues"))
    length(result.modal_damped_frequencies_hz) == modes ||
        throw(DimensionMismatch(
            "modal damped-frequency count does not match the eigenvalues"))
    length(result.modal_damping_ratios) == modes ||
        throw(DimensionMismatch(
            "modal damping-ratio count does not match the eigenvalues"))
    length(result.modal_equation_errors) == modes ||
        throw(DimensionMismatch(
            "modal equation-error count does not match the eigenvalues"))
    size(result.modal_mode_shapes) == (length(result.variable_names), modes) ||
        throw(DimensionMismatch(
            "modal mode-shape dimensions do not match the variable and mode catalogs"))
    modes
end

"""Write one summary row per retained mode to a CSV stream."""
function write_modal_summary_csv(output::IO, result::StoredSimulationResult)
    modes = check_modal_dimensions(result)
    println(output, join(("mode", "eigenvalue_real", "eigenvalue_imaginary",
        "natural_frequency_hz", "damped_frequency_hz", "damping_ratio",
        "equation_error"), ','))
    for mode in 1:modes
        eigenvalue = result.modal_eigenvalues[mode]
        println(output, join(repr.((mode, real(eigenvalue), imag(eigenvalue),
            result.modal_natural_frequencies_hz[mode],
            result.modal_damped_frequencies_hz[mode],
            result.modal_damping_ratios[mode],
            result.modal_equation_errors[mode])), ','))
    end
    output
end

"""Write complex canonical mode shapes in long form to a CSV stream."""
function write_modal_shapes_csv(output::IO, result::StoredSimulationResult;
        variables = nothing)
    modes = check_modal_dimensions(result)
    requested = isnothing(variables) ? nothing : String.(variables)
    columns, _ = selected_columns(result, requested)
    println(output, join(("mode", "component", "variable", "kind", "real",
        "imaginary", "magnitude", "phase_deg"), ','))
    for mode in 1:modes, column in columns
        value = result.modal_mode_shapes[column, mode]
        fields = (repr(mode), csv_field(result.variable_components[column]),
            csv_field(result.variable_names[column]),
            csv_field(result.variable_kinds[column]), repr(real(value)),
            repr(imag(value)), repr(abs(value)), repr(rad2deg(angle(value))))
        println(output, join(fields, ','))
    end
    output
end

"""
    export_modal_csv(summary_output, shapes_output, result; variables=nothing)

Write a modal summary and long-form complex mode shapes to two CSV streams.
Optional qualified `component.variable` names restrict the shape rows without
changing the mode summary.
"""
function export_modal_csv(summary_output::IO, shapes_output::IO,
        result::StoredSimulationResult; variables = nothing)
    result.analysis_mode == :modal ||
        throw(ArgumentError("modal CSV export requires a modal result"))
    check_modal_dimensions(result)
    requested = isnothing(variables) ? nothing : String.(variables)
    selected_columns(result, requested)
    write_modal_summary_csv(summary_output, result)
    write_modal_shapes_csv(shapes_output, result; variables = requested)
    summary_output, shapes_output
end

"""
    export_result_csv(output::IO, result; variables=nothing)

Write time and selected canonical histories to a CSV stream. Variable names are
qualified as `component.variable`; omit `variables` to export all histories.
Time is always the first column and must not be requested explicitly.
"""
function export_result_csv(output::IO, result::StoredSimulationResult;
        variables = nothing)
    result.analysis_mode == :modal && throw(ArgumentError(
        "modal export requires two streams; use export_modal_csv or the " *
        "path form of export_result_csv"))
    requested = isnothing(variables) ? nothing : String.(variables)
    columns, names = selected_columns(result, requested)
    println(output, join(csv_field.(["time"; names]), ','))
    for row in eachindex(result.times)
        print(output, repr(result.times[row]))
        for column in columns
            print(output, ',', repr(result.values[row, column]))
        end
        println(output)
    end
    output
end

export_result_csv(output::IO, input_path::AbstractString; kwargs...) =
    export_result_csv(output, read_result(input_path); kwargs...)

"""
    export_result_csv(output_path, result; overwrite=false, variables=nothing)
    export_result_csv(output_path, input_path; overwrite=false, variables=nothing)

Export an in-memory or stored `.simp` result to CSV. Existing output files
require `overwrite=true`. A time-domain result returns `output_path`. For a
modal result, `output_path` is treated as a stem and the function returns the
generated `-modes.csv` and `-mode-shapes.csv` paths.
"""
function modal_csv_paths(output_path::AbstractString)
    stem, extension = splitext(output_path)
    base = lowercase(extension) == ".csv" ? stem : output_path
    (modes = base * "-modes.csv", mode_shapes = base * "-mode-shapes.csv")
end

function check_output_path(path, overwrite)
    isfile(path) && !overwrite &&
        throw(ArgumentError("CSV file already exists: $path"))
end

"""
    export_modal_csv(summary_path, shapes_path, result;
        overwrite=false, variables=nothing)

Export modal properties and long-form complex mode shapes to separate CSV
files. Returns the two paths as `(modes=summary_path,
mode_shapes=shapes_path)`.
"""
function export_modal_csv(summary_path::AbstractString,
        shapes_path::AbstractString, result::StoredSimulationResult;
        overwrite = false, variables = nothing)
    result.analysis_mode == :modal ||
        throw(ArgumentError("modal CSV export requires a modal result"))
    check_modal_dimensions(result)
    requested = isnothing(variables) ? nothing : String.(variables)
    selected_columns(result, requested)
    check_output_path(summary_path, overwrite)
    check_output_path(shapes_path, overwrite)
    open(summary_path, "w") do summary_output
        open(shapes_path, "w") do shapes_output
            export_modal_csv(summary_output, shapes_output, result;
                variables = requested)
        end
    end
    (modes = String(summary_path), mode_shapes = String(shapes_path))
end

function export_result_csv(output_path::AbstractString,
        result::StoredSimulationResult; overwrite = false, kwargs...)
    if result.analysis_mode == :modal
        paths = modal_csv_paths(output_path)
        return export_modal_csv(paths.modes, paths.mode_shapes, result;
            overwrite, kwargs...)
    end
    check_output_path(output_path, overwrite)
    open(output_path, "w") do output
        export_result_csv(output, result; kwargs...)
    end
    output_path
end

function export_result_csv(output_path::AbstractString,
        input_path::AbstractString; kwargs...)
    export_result_csv(output_path, read_result(input_path); kwargs...)
end

function print_usage(io)
    println(io, "usage: simpCSV RESULT.simp OUTPUT.csv",
        " [VARIABLE ...] [--overwrite]")
end

"""Command-line entry point for `.simp` to CSV conversion."""
function export_result_main(args = ARGS; output = stdout, error = stderr)
    overwrite = "--overwrite" in args
    positional = filter(!=("--overwrite"), args)
    length(positional) >= 2 || begin
        print_usage(error)
        return 2
    end
    input_path, output_path = positional[1:2]
    variables = length(positional) > 2 ? positional[3:end] : nothing
    written = export_result_csv(output_path, input_path; variables, overwrite)
    if written isa NamedTuple
        println(output, "Wrote ", written.modes)
        println(output, "Wrote ", written.mode_shapes)
    else
        println(output, "Wrote ", written)
    end
    return 0
end

end
