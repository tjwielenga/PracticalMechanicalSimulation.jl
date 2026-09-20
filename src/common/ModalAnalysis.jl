"""
    ModalAnalysis

Small-signal modal analysis of the complete active implicit model. Algebraic
definitions, ideal constraints, reactions, and applied-force laws stay in the
sparse descriptor equations. Sparse shift-invert solves eliminate them
implicitly, leaving only a dense eigenproblem whose order is twice the number
of selected mechanical degrees of freedom plus any independent first-order
states owned by elements.
"""
module ModalAnalysis

using LinearAlgebra
using SparseArrays
using ..AutomaticAnalysis

export solve_modal_system

const DISPLACEMENT_KINDS = Set((:position, :orientation,
                                :relative_position))

"""
    solve_modal_system(loaded, state, derivative, time;
                       number_of_modes, frequency_shift_hz,
                       tolerance)

Linearize the complete active implicit equation set and calculate the finite
eigenpairs nearest `frequency_shift_hz`. Algebraic definitions and force laws
remain in the sparse descriptor pencil; only its low-rank derivative part is
reduced after one sparse shift-invert factorization. Dimension-specific
callers may supply reduced canonical index sets and a `jacobian_evaluator`
that returns their tangent-coordinate Jacobian for a derivative coefficient.
"""
function solve_modal_system(loaded, state, derivative, time;
        number_of_modes = loaded.analysis.number_of_modes,
        frequency_shift_hz = loaded.analysis.frequency_shift_hz,
        tolerance = loaded.analysis.modal_tolerance,
        variable_indices = loaded.active_variable_indices,
        equation_indices = loaded.active_equation_indices,
        jacobian_evaluator = nothing)
    number_of_modes >= 1 || throw(ArgumentError(
        "number_of_modes must be positive"))
    frequency_shift_hz >= 0 || throw(ArgumentError(
        "frequency_shift_hz must be nonnegative"))
    tolerance > 0 || throw(ArgumentError("modal tolerance must be positive"))

    state_jacobian, combined_jacobian = if isnothing(jacobian_evaluator)
        selection = AnalysisSelection(Dynamics(), variable_indices,
            equation_indices)
        (evaluate_analysis_sparse_jacobian(loaded.model, selection, time,
             state, derivative, 0.0),
         evaluate_analysis_sparse_jacobian(loaded.model, selection, time,
             state, derivative, 1.0))
    else
        (jacobian_evaluator(0.0), jacobian_evaluator(1.0))
    end
    size(state_jacobian) == (length(equation_indices),
        length(variable_indices)) || throw(DimensionMismatch(
        "modal state Jacobian does not match its equation and variable sets"))
    size(combined_jacobian) == size(state_jacobian) ||
        throw(DimensionMismatch("modal Jacobians must have equal dimensions"))
    derivative_jacobian = combined_jacobian - state_jacobian
    dropzeros!(derivative_jacobian)

    _, derivative_columns, _ = findnz(derivative_jacobian)
    differential_columns = sort!(unique(derivative_columns))
    isempty(differential_columns) && throw(ArgumentError(
        "modal analysis requires at least one independent state"))
    variable_declarations = loaded.layout.catalog.variables
    internal_columns = [column for column in differential_columns
        if variable_declarations[variable_indices[column]].kind in
           (:internal_state, :user_state_hold, :user_state_steady)]
    expected_columns = 2 * loaded.analysis.degrees_of_freedom +
        length(internal_columns)
    length(differential_columns) == expected_columns || throw(ArgumentError(
        "the implicit derivative matrix has $(length(differential_columns)) " *
        "differential columns; expected $expected_columns from state " *
        "selection and independent element states"))

    target_rate = 2pi * frequency_shift_hz
    guard = sqrt(eps(Float64)) * max(1.0, target_rate)
    shift = iszero(target_rate) ? -guard : complex(-guard, target_rate)
    shifted_jacobian = state_jacobian + shift * derivative_jacobian
    factorization = lu(shifted_jacobian; check = true)

    # E has columns only for the selected positions and velocities and any
    # independent first-order element states. If
    # U = (J + shift*E)^(-1) E[:, differential_columns], the nonzero
    # eigenvalues of the full shift-invert operator are exactly those of
    # U[differential_columns, :]. This requires sparse solves with the full
    # implicit equations but only a small dense eigenproblem of order 2*DOF
    # plus the number of first-order element states.
    right_hand_sides = Matrix{eltype(shifted_jacobian)}(
        derivative_jacobian[:, differential_columns])
    transformed_columns = ldiv!(factorization, right_hand_sides)
    reduced_operator = Matrix(transformed_columns[differential_columns, :])
    reduced_eigen = eigen(reduced_operator)

    transformed_values = reduced_eigen.values
    finite_indices = [index for index in eachindex(transformed_values)
        if isfinite(transformed_values[index]) &&
           abs(transformed_values[index]) > tolerance]
    isempty(finite_indices) && throw(ArgumentError(
        "the linearized implicit system has no finite dynamic eigenvalues"))

    all_eigenvalues = ComplexF64[
        shift - inv(transformed_values[index]) for index in finite_indices]
    all_local_modes = Matrix{ComplexF64}(undef,
        length(variable_indices), length(finite_indices))
    for (column, reduced_index) in enumerate(finite_indices)
        theta = transformed_values[reduced_index]
        all_local_modes[:, column] .= transformed_columns *
            reduced_eigen.vectors[:, reduced_index] / theta
    end

    candidate_indices = [index for index in eachindex(all_eigenvalues)
        if imag(all_eigenvalues[index]) >=
           -tolerance * max(1.0, abs(all_eigenvalues[index]))]
    sort!(candidate_indices; by = index -> (
        abs(all_eigenvalues[index] - complex(0.0, target_rate)),
        abs(all_eigenvalues[index]), real(all_eigenvalues[index])))

    selected_indices = Int[]
    for candidate in candidate_indices
        push!(selected_indices, candidate)
        length(selected_indices) == min(number_of_modes,
            length(candidate_indices)) && break
    end
    isempty(selected_indices) && throw(ArgumentError(
        "no modal eigenpairs could be selected from the finite spectrum"))

    active_lookup = Dict(canonical => local_index
        for (local_index, canonical) in
        enumerate(variable_indices))
    displacement_variables = [variable for variable in
        loaded.layout.catalog.variables
        if variable.kind in DISPLACEMENT_KINDS &&
           haskey(active_lookup, variable.index)]
    isempty(displacement_variables) && throw(ArgumentError(
        "modal analysis found no physical displacement variables"))

    mode_shapes = zeros(ComplexF64, length(loaded.initial_values),
        length(selected_indices))
    equation_errors = zeros(Float64, length(selected_indices))
    for (mode_number, spectrum_index) in enumerate(selected_indices)
        local_mode = copy(all_local_modes[:, spectrum_index])
        scaled_amplitudes = [abs(local_mode[active_lookup[variable.index]]) /
            variable.scale for variable in displacement_variables]
        pivot_number = argmax(scaled_amplitudes)
        amplitude = scaled_amplitudes[pivot_number]
        amplitude > tolerance || throw(ArgumentError(
            "modal eigenvector has no measurable physical displacement"))
        pivot = active_lookup[displacement_variables[pivot_number].index]
        local_mode .*= cis(-angle(local_mode[pivot])) / amplitude
        mode_shapes[variable_indices, mode_number] .= local_mode

        value = all_eigenvalues[spectrum_index]
        equation_error = state_jacobian * local_mode +
            value * (derivative_jacobian * local_mode)
        denominator = max(norm(state_jacobian * local_mode) +
            abs(value) * norm(derivative_jacobian * local_mode), eps(Float64))
        equation_errors[mode_number] = norm(equation_error) / denominator
    end

    eigenvalues = all_eigenvalues[selected_indices]
    natural_frequencies_hz = abs.(eigenvalues) ./ (2pi)
    damped_frequencies_hz = abs.(imag.(eigenvalues)) ./ (2pi)
    damping_ratios = [abs(value) <= tolerance ? 0.0 :
        -real(value) / abs(value) for value in eigenvalues]

    (; eigenvalues, natural_frequencies_hz, damped_frequencies_hz,
       damping_ratios, mode_shapes, equation_errors,
       state_jacobian, derivative_jacobian, differential_columns,
       shift, finite_eigenvalue_count = length(finite_indices),
       sparse_factorizations = 1)
end

end
