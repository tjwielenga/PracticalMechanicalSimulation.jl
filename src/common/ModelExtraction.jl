module ModelExtraction

using ..ResultIO

export extract_model, extract_model_main

function require_embedded_model(result::StoredSimulationResult)
    isempty(result.model_source) &&
        throw(ArgumentError("result file does not contain an embedded TOML model"))
    result.model_source
end

"""
    extract_model(output::IO, result)

Write the verbatim TOML embedded in a stored result to `output`. An empty model
record is rejected. The stream remains open and is returned.
"""
function extract_model(output::IO, result::StoredSimulationResult)
    write(output, require_embedded_model(result))
    output
end

extract_model(output::IO, input_path::AbstractString) =
    extract_model(output, read_result(input_path))

"""
    extract_model(output_path, result; overwrite=false)
    extract_model(output_path, input_path; overwrite=false)

Extract embedded TOML from an in-memory or on-disk result. Existing output
files require `overwrite=true`. Returns `output_path`.
"""
function extract_model(output_path::AbstractString,
        result::StoredSimulationResult; overwrite = false)
    isfile(output_path) && !overwrite &&
        throw(ArgumentError("model file already exists: $output_path"))
    open(output_path, "w") do output
        extract_model(output, result)
    end
    output_path
end

function extract_model(output_path::AbstractString,
        input_path::AbstractString; kwargs...)
    extract_model(output_path, read_result(input_path); kwargs...)
end

function print_usage(io)
    println(io, "usage: simpExtract RESULT.simp MODEL.toml|- [--overwrite]")
end

"""Command-line entry point for extracting a result's embedded TOML model."""
function extract_model_main(args = ARGS; output = stdout, error = stderr)
    overwrite = "--overwrite" in args
    positional = filter(!=("--overwrite"), args)
    length(positional) == 2 || begin
        print_usage(error)
        return 2
    end
    input_path, output_path = positional
    if output_path == "-"
        extract_model(output, input_path)
    else
        extract_model(output_path, input_path; overwrite)
        println(output, "Wrote ", output_path)
    end
    return 0
end

end
