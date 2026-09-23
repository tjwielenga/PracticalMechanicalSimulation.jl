using PracticalMechanicalSimulation

isdefined(Main, :ViewerData) ||
    include(joinpath(@__DIR__, "ViewerData.jl"))
isdefined(Main, :StoredResultViewer) ||
    include(joinpath(@__DIR__, "StoredResultViewer.jl"))
isdefined(Main, :PortableViewerDocument) ||
    include(joinpath(@__DIR__, "PortableViewerDocument.jl"))

using .StoredResultViewer
using .PortableViewerDocument

function reconstructed_viewer_results(input; mode = nothing)
    stored = read_result(input)
    if stored.analysis_mode == :modal && isnothing(mode)
        mode_count = length(stored.modal_eigenvalues)
        results = [stored_mechanism_result(input; mode = mode_number)
            for mode_number in 1:mode_count]
        labels = ["Mode $mode_number — " *
            string(round(stored.modal_natural_frequencies_hz[mode_number];
                sigdigits = 6)) * " Hz" for mode_number in 1:mode_count]
        return (; results, labels, choice_name = "Mode")
    end
    if stored.analysis_mode in (:dynamic, :kinematic) &&
            !isempty(stored.static_snapshots)
        results = [
            stored_mechanism_result(input; analysis_mode = :static),
            stored_mechanism_result(input),
        ]
        return (; results,
            labels = ["Static initialization", "Dynamic"],
            choice_name = "Analysis")
    end
    selected_mode = isnothing(mode) ? 1 : mode
    result = stored_mechanism_result(input; mode = selected_mode)
    label = stored.analysis_mode == :modal ? "Mode $selected_mode" : "Result"
    (; results = [result], labels = [label], choice_name = "Result")
end

function write_stored_result_graphics(input; mode = nothing)
    reconstructed = reconstructed_viewer_results(input; mode)
    write_graphics(input, reconstructed.results;
        labels = reconstructed.labels,
        choice_name = reconstructed.choice_name)
end
