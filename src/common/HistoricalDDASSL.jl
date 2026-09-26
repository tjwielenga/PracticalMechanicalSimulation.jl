"""
    HistoricalDDASSL

Variable-step, variable-order BDF integrator for square implicit systems. It
is a Julia reconstruction of the DDASSL method adapted to this project's
explicit equation levels, sparse analytical Jacobians, monitored physical
errors, root location, and equivalent-equation reconfiguration.

The name separates this maintained solver from paper reproductions of older
algorithms. It is the production integrator used by the planar modeler.
"""
module HistoricalDDASSL

using LinearAlgebra
using SciMLBase
using SparseArrays

export DASSLEvent, DASSLErrorMonitor, DASSLOptions, DASSLResult, DASSLStats,
    dassl, error_control_from_levels, level_factors, scaled_newton_matrix

# SuiteSparse numbers UMFPACK_INFO entries from zero; Julia arrays from one.
const UMFPACK_RCOND_INFO_INDEX = 68

"""
    DASSLOptions{T}

Numerical controls for [`dassl`](@ref). `rtol` and `atol` define the weighted
local-error test. Step limits and `maximum_order` control the variable-step BDF
history. The Newton tolerances govern correction convergence independently of
the integration-error test.

`maximum_factorization_steps` and `factorization_coefficient_tolerance`
control reuse of the numerical Newton factorization. Symbolic sparse analysis
is retained until repeated corrector failures reach
`symbolic_reanalysis_failures`. These reuse controls affect cost, not the
equations selected by the mechanical model.
"""
Base.@kwdef struct DASSLOptions{T<:AbstractFloat}
    rtol::T = T(1.0e-6)
    atol::T = T(1.0e-8)
    initial_step::T = zero(T)
    minimum_step::T = zero(T)
    maximum_step::T = T(Inf)
    maximum_order::Int = 5
    maximum_steps::Int = 100_000
    maximum_newton_iterations::Int = 8
    maximum_factorization_steps::Int = 5
    factorization_coefficient_tolerance::T = T(0.25)
    symbolic_reanalysis_failures::Int = 2
    newton_rtol::T = T(1.0e-8)
    newton_atol::T = T(1.0e-10)
    residual_tolerance::T = T(1.0e-8)
    safety::T = T(0.85)
    minimum_factor::T = T(0.2)
    maximum_factor::T = T(2.0)
    order_change_advantage::T = T(1.2)
    minimum_order_steps::Int = 2
end

"""
Work counters returned with a DASSL solution.

`factorizations` counts numerical factorizations and
`symbolic_factorizations` counts new sparse orderings. `history_restarts`
includes event restarts, while `state_reselections` counts successful
equivalent-equation changes requested by the mechanical runner.
"""
mutable struct DASSLStats
    accepted_steps::Int
    rejected_steps::Int
    residual_evaluations::Int
    jacobian_evaluations::Int
    factorizations::Int
    symbolic_factorizations::Int
    newton_iterations::Int
    corrector_failures::Int
    root_evaluations::Int
    events_found::Int
    history_restarts::Int
    state_reselections::Int
end

DASSLStats() = DASSLStats(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

"""A located zero of one or more user-supplied discontinuity functions."""
struct DASSLEvent{T,V}
    time::T
    indices::Vector{Int}
    directions::Vector{Int}
    y::V
    yprime::V
end

"""Normalized predictor-error history for a separately monitored set of variables."""
struct DASSLErrorMonitor{T}
    mask::BitVector
    controlled_errors::Vector{T}
    rms_errors::Vector{T}
    maximum_errors::Vector{T}
    maximum_error_indices::Vector{Int}
end

"""Solution and diagnostic information returned by [`dassl`](@ref)."""
struct DASSLResult{T,V}
    t::Vector{T}
    y::Vector{V}
    yprime::Vector{V}
    orders::Vector{Int}
    steps::Vector{T}
    stats::DASSLStats
    retcode::SciMLBase.ReturnCode.T
    message::String
    variable_levels::Union{Nothing,Vector{Int}}
    equation_levels::Union{Nothing,Vector{Int}}
    differential_vars::Union{Nothing,BitVector}
    error_control::BitVector
    error_monitor::Union{Nothing,DASSLErrorMonitor{T}}
    deficit::Int
    events::Vector{DASSLEvent{T,V}}
end

# SciML names the saved state and derivative arrays `u` and `du`. Retain the
# mechanically familiar `y` and `yprime` fields while exposing both spellings.
function Base.getproperty(solution::DASSLResult, name::Symbol)
    name === :u && return getfield(solution, :y)
    name === :du && return getfield(solution, :yprime)
    return getfield(solution, name)
end

function Base.propertynames(solution::DASSLResult, private::Bool = false)
    names = fieldnames(typeof(solution))
    return (names..., :u, :du)
end

# Familiar SciML DEStats names, without hiding the explicit local names.
function Base.getproperty(stats::DASSLStats, name::Symbol)
    name === :nf && return getfield(stats, :residual_evaluations)
    name === :njacs && return getfield(stats, :jacobian_evaluations)
    name === :nw && return getfield(stats, :factorizations)
    name === :nnonliniter && return getfield(stats, :newton_iterations)
    name === :naccept && return getfield(stats, :accepted_steps)
    name === :nreject && return getfield(stats, :rejected_steps)
    return getfield(stats, name)
end

function Base.propertynames(stats::DASSLStats, private::Bool = false)
    names = fieldnames(typeof(stats))
    return (names..., :nf, :njacs, :nw, :nnonliniter, :naccept, :nreject)
end

Base.length(solution::DASSLResult) = length(solution.t)
Base.firstindex(solution::DASSLResult) = firstindex(solution.t)
Base.lastindex(solution::DASSLResult) = lastindex(solution.t)
Base.size(solution::DASSLResult) = (length(first(solution.y)), length(solution.t))
Base.size(solution::DASSLResult, dimension::Int) = size(solution)[dimension]
Base.axes(solution::DASSLResult) = (Base.OneTo(size(solution, 1)),
    Base.OneTo(size(solution, 2)))
Base.axes(solution::DASSLResult, dimension::Int) = axes(solution)[dimension]
Base.firstindex(solution::DASSLResult, dimension::Int) =
    first(axes(solution, dimension))
Base.lastindex(solution::DASSLResult, dimension::Int) =
    last(axes(solution, dimension))
Base.getindex(solution::DASSLResult, time_index::Int) = solution.y[time_index]
Base.getindex(solution::DASSLResult, variable, time_index::Int) =
    solution.y[time_index][variable]
Base.getindex(solution::DASSLResult, variable, ::Colon) =
    [state[variable] for state in solution.y]
Base.iterate(solution::DASSLResult, state = 1) = state > length(solution) ?
    nothing : (solution[state], state + 1)

function interpolation_interval(solution::DASSLResult, time)
    first_time, last_time = first(solution.t), last(solution.t)
    lower, upper = minmax(first_time, last_time)
    lower <= time <= upper || throw(DomainError(time,
        "interpolation time lies outside the solution interval"))
    length(solution.t) > 1 || return 1
    if first_time <= last_time
        return clamp(searchsortedlast(solution.t, time), 1,
            length(solution.t) - 1)
    end
    return clamp(searchsortedlast(-solution.t, -time), 1,
        length(solution.t) - 1)
end

select_components(value, ::Nothing) = value
select_components(value, ::Colon) = value
select_components(value, indices) = value[indices]

"""
Evaluate the active BDF history polynomial over an accepted step.

The order saved at the right endpoint identifies the corrected BDF polynomial
for that step. Its right endpoint and preceding accepted history values provide
the same interpolation nodes that formed the BDF derivative relation.
"""
function interpolate_state(solution::DASSLResult, time)
    if length(solution.t) == 1
        return copy(first(solution.y)), copy(first(solution.yprime))
    end
    interval = interpolation_interval(solution, time)
    T = eltype(solution.t)
    interpolation_time = T(time)

    # At a saved endpoint, retain the corrected state and its corresponding
    # implicit derivative exactly. Adjacent BDF polynomials can have slightly
    # different derivatives there after a step-size or order change.
    if interpolation_time == solution.t[interval]
        return copy(solution.y[interval]), copy(solution.yprime[interval])
    elseif interpolation_time == solution.t[interval + 1]
        return copy(solution.y[interval + 1]),
            copy(solution.yprime[interval + 1])
    end

    right = interval + 1
    order = min(solution.orders[right], right - 1)
    node_indices = right:-1:(right - order)
    nodes = solution.t[node_indices]
    values = solution.y[node_indices]
    return interpolation_state(interpolation_time, nodes, values)
end

function (solution::DASSLResult)(time::Real; idxs = nothing)
    value, _ = interpolate_state(solution, time)
    return select_components(value, idxs)
end

function (solution::DASSLResult)(time::Real, ::Type{Val{1}}; idxs = nothing)
    _, derivative = interpolate_state(solution, time)
    return select_components(derivative, idxs)
end

function (solution::DASSLResult)(times::AbstractVector; idxs = nothing)
    [solution(time; idxs) for time in times]
end

function (solution::DASSLResult)(times::AbstractVector, ::Type{Val{1}};
        idxs = nothing)
    [solution(time, Val{1}; idxs) for time in times]
end

scaled_rms(v, scale, mask) = sqrt(sum(abs2(v[i] / scale[i])
    for i in eachindex(v) if mask[i]) / count(mask))

function monitored_error_summary(y, predictor, options, variable_factors,
        mask, order)
    multiplier = inv(eltype(y)(order + 1))
    sum_squares = zero(eltype(y))
    maximum_error = zero(eltype(y))
    maximum_index = 0
    for index in eachindex(y)
        mask[index] || continue
        factor = variable_factors[index]
        scaled_y = factor * y[index]
        scaled_predictor = factor * predictor[index]
        scale = options.atol + options.rtol * max(
            abs(scaled_y), abs(scaled_predictor))
        error = multiplier * abs((scaled_y - scaled_predictor) / scale)
        sum_squares += abs2(error)
        if maximum_index == 0 || error > maximum_error
            maximum_error = error
            maximum_index = index
        end
    end
    return sqrt(sum_squares / count(mask)), maximum_error, maximum_index
end

"""Construct the default mechanical integration-error control mask."""
function error_control_from_levels(variable_levels, deficit::Integer;
        differential_vars = trues(length(variable_levels)))
    deficit >= 0 || throw(ArgumentError("deficit must be nonnegative"))
    length(differential_vars) == length(variable_levels) ||
        throw(DimensionMismatch(
            "differential_vars must match variable_levels"))
    maximum_controlled_level = max(0, min(1, 2 - deficit))
    return BitVector(differential && level <= maximum_controlled_level
        for (level, differential) in zip(variable_levels, differential_vars))
end

"""Factors that convert level-`l` quantities toward level zero."""
function level_factors(levels, step::T) where T
    scale_step = max(abs(step), sqrt(eps(T)))
    return T[scale_step^level for level in levels]
end

"""
Scale a Newton matrix by equation rows and variable columns.

For equation level `e[i]` and variable level `v[j]`, the resulting entry is
`abs(step)^(e[i]-v[j]) * matrix[i,j]`. The scaled correction is related to
the physical correction by `scaled_correction[j] = abs(step)^v[j] *
correction[j]`.
"""
function scaled_newton_matrix(matrix, equation_factors, variable_factors)
    return Diagonal(equation_factors) * matrix *
        Diagonal(inv.(variable_factors))
end

function component_scale(y, reference, options)
    options.atol .+ options.rtol .* max.(abs.(y), abs.(reference))
end

"""
Return derivative weights at `x0` for a polynomial through `nodes`.

For BDF order `k`, `nodes` contains the new time followed by `k` past
times. The result satisfies

```
yprime_new = sum(weights[j] * y_at_node[j])
```

exactly for polynomials through degree `k`.
"""
function derivative_weights(x0::T, nodes::AbstractVector{T}) where T
    n = length(nodes)
    offsets = nodes .- x0
    moments = Matrix{T}(undef, n, n)
    for degree in 0:n-1, j in 1:n
        moments[degree + 1, j] = offsets[j]^degree
    end
    rhs = zeros(T, n)
    rhs[2] = one(T)
    return moments \ rhs
end

"""Evaluate the Lagrange interpolant through `values` at `x`."""
function extrapolate(x::T, nodes::AbstractVector{T}, values) where T
    result = zeros(T, length(first(values)))
    for j in eachindex(nodes)
        coefficient = one(T)
        for m in eachindex(nodes)
            m == j && continue
            coefficient *= (x - nodes[m]) / (nodes[j] - nodes[m])
        end
        result .+= coefficient .* values[j]
    end
    return result
end

function interpolation_state(time, nodes, values)
    value = extrapolate(time, nodes, values)
    weights = derivative_weights(time, nodes)
    derivative = zeros(eltype(value), length(value))
    for index in eachindex(nodes)
        derivative .+= weights[index] .* values[index]
    end
    return value, derivative
end

function evaluate_roots!(root!, roots, time, y, yprime, parameter, stats)
    root!(roots, time, y, yprime, parameter)
    stats.root_evaluations += 1
    all(isfinite, roots) || throw(ArgumentError(
        "discontinuity functions must return finite values"))
    return roots
end

function crossing_direction(left, right)
    right > left && return 1
    right < left && return -1
    return 0
end

function crossing_detected(left, right, requested_direction)
    crossed = (left < 0 && right >= 0) || (left > 0 && right <= 0)
    crossed || return false
    direction = crossing_direction(left, right)
    return requested_direction == 0 || direction == requested_direction
end

function locate_root(root!, root_index, left_time, right_time,
        left_value, right_value, nodes, values, parameter, stats, root_count,
        time_tolerance)
    T = eltype(nodes)
    roots = zeros(T, root_count)
    left, right = left_time, right_time
    gleft, gright = left_value, right_value
    if iszero(gright)
        y, yprime = interpolation_state(right, nodes, values)
        return right, y, yprime
    end
    for _ in 1:100
        abs(right - left) <= time_tolerance && break
        # Bisection is deliberately used here: the discontinuity functions can
        # be much less smooth than the BDF interpolation polynomial.
        middle = (left + right) / 2
        y, yprime = interpolation_state(middle, nodes, values)
        evaluate_roots!(root!, roots, middle, y, yprime, parameter, stats)
        gmiddle = roots[root_index]
        if iszero(gmiddle)
            return middle, y, yprime
        elseif signbit(gmiddle) == signbit(gleft)
            left, gleft = middle, gmiddle
        else
            right, gright = middle, gmiddle
        end
    end
    time = (left + right) / 2
    y, yprime = interpolation_state(time, nodes, values)
    return time, y, yprime
end

function find_first_event(root!, left_time, right_time, left_roots,
        right_roots, root_directions, nodes, values, parameter, stats,
        time_tolerance)
    candidates = Int[index for index in eachindex(left_roots)
        if crossing_detected(left_roots[index], right_roots[index],
            root_directions[index])]
    isempty(candidates) && return nothing
    located = [locate_root(root!, index, left_time, right_time,
        left_roots[index], right_roots[index], nodes, values, parameter,
        stats, length(left_roots), time_tolerance) for index in candidates]
    direction = sign(right_time - left_time)
    first_location = argmin(direction * located[index][1]
        for index in eachindex(located))
    event_time, event_y, event_yprime = located[first_location]
    simultaneous = Int[]
    directions = Int[]
    for (candidate, location) in zip(candidates, located)
        abs(location[1] - event_time) <= time_tolerance || continue
        push!(simultaneous, candidate)
        push!(directions, crossing_direction(left_roots[candidate],
            right_roots[candidate]))
    end
    return DASSLEvent(event_time, simultaneous, directions,
        event_y, event_yprime)
end

function numerical_iteration_matrix!(matrix, residual!, residual, trial,
        t, y, yprime, cj, parameter, stats)
    T = eltype(y)
    perturbed_y = copy(y)
    perturbed_yprime = copy(yprime)
    square_epsilon = sqrt(eps(T))
    for j in eachindex(y)
        delta = square_epsilon * max(abs(y[j]), abs(yprime[j] / cj), one(T))
        delta = max(delta, eps(T))
        perturbed_y[j] = y[j] + delta
        perturbed_yprime[j] = yprime[j] + cj * delta
        residual!(trial, t, perturbed_y, perturbed_yprime, parameter)
        stats.residual_evaluations += 1
        @views matrix[:, j] .= (trial .- residual) ./ delta
        perturbed_y[j] = y[j]
        perturbed_yprime[j] = yprime[j]
    end
    stats.jacobian_evaluations += 1
    return matrix
end

function evaluate_iteration_matrix!(matrix, jacobian!, residual!, residual,
        trial, t, y, yprime, cj, parameter, stats)
    if isnothing(jacobian!)
        numerical_iteration_matrix!(matrix, residual!, residual, trial,
            t, y, yprime, cj, parameter, stats)
    else
        jacobian!(matrix, t, y, yprime, cj, parameter)
        stats.jacobian_evaluations += 1
    end
    return matrix
end

"""Internal outcome of one predicted and corrected BDF step attempt."""
struct StepAttempt{T,V}
    converged::Bool
    y::V
    yprime::V
    error::T
    residual_norm::T
    iterations::Int
    reason::Symbol
    order_errors::Vector{T}
    monitored_rms_error::T
    monitored_maximum_error::T
    monitored_maximum_index::Int
end

"""
Reusable Newton-matrix factorization.

Sparse matrices retain UMFPACK's symbolic ordering while numerical values are
refactored. The numerical factors themselves may serve several steps and all
Newton corrections on those steps until the scaled BDF coefficient changes
too much or convergence shows that a refresh is needed.
"""
mutable struct LinearFactorizationCache{T}
    factorization::Any
    numeric_valid::Bool
    scaled_coefficient::T
    steps::Int
    reciprocal_condition::T
    condition_reference::T
end

LinearFactorizationCache(::Type{T}) where T =
    LinearFactorizationCache{T}(nothing, false, zero(T), 0, T(NaN), T(NaN))

function invalidate_numeric_factorization!(cache::LinearFactorizationCache)
    cache.numeric_valid = false
    cache.steps = 0
    cache
end

function invalidate_symbolic_factorization!(cache::LinearFactorizationCache)
    cache.factorization = nothing
    cache.numeric_valid = false
    cache.steps = 0
    cache.reciprocal_condition = typeof(cache.reciprocal_condition)(NaN)
    cache.condition_reference = typeof(cache.condition_reference)(NaN)
    cache
end

function update_factorization_condition!(cache::LinearFactorizationCache)
    factorization = cache.factorization
    condition = hasproperty(factorization, :info) &&
        length(factorization.info) >= UMFPACK_RCOND_INFO_INDEX ?
        typeof(cache.reciprocal_condition)(
            factorization.info[UMFPACK_RCOND_INFO_INDEX]) :
        typeof(cache.reciprocal_condition)(NaN)
    cache.reciprocal_condition = condition
    if isfinite(condition) && condition > zero(condition) &&
            (!isfinite(cache.condition_reference) ||
             condition > cache.condition_reference)
        cache.condition_reference = condition
    end
    cache
end

function conditioning_warning(cache::LinearFactorizationCache)
    current = cache.reciprocal_condition
    reference = cache.condition_reference
    isfinite(current) && current > zero(current) &&
        isfinite(reference) && reference > zero(reference) &&
        current < reference / 4
end

function factor_iteration_matrix!(cache::LinearFactorizationCache,
        matrix, scaled_coefficient, stats)
    if issparse(matrix)
        if isnothing(cache.factorization)
            cache.factorization = lu(matrix; check = true)
            stats.symbolic_factorizations += 1
        else
            cache.factorization = lu!(cache.factorization, matrix;
                check = true, reuse_symbolic = true)
        end
    else
        cache.factorization = lu!(matrix; check = true)
    end
    cache.numeric_valid = true
    cache.scaled_coefficient = scaled_coefficient
    cache.steps = 1
    update_factorization_condition!(cache)
    stats.factorizations += 1
    cache.factorization
end

function coefficient_requires_refactorization(cache::LinearFactorizationCache,
        scaled_coefficient, tolerance)
    cache.numeric_valid || return true
    iszero(cache.scaled_coefficient) && return true
    ratio = abs(scaled_coefficient / cache.scaled_coefficient)
    lower = (one(tolerance) - tolerance) / (one(tolerance) + tolerance)
    upper = inv(lower)
    return !isfinite(ratio) || ratio < lower || ratio > upper
end

function proposed_step_factor(error, order, options; upper = options.maximum_factor)
    !isfinite(error) && return options.minimum_factor
    error <= eps(typeof(error)) && return upper
    return clamp(options.safety * error^(-inv(typeof(error)(order + 1))),
        options.minimum_factor, upper)
end

function candidate_order_errors(y, tnew, history_t, history_y,
        options, variable_factors, error_mask)
    T = eltype(y)
    errors = fill(T(Inf), options.maximum_order)
    maximum_candidate = min(options.maximum_order, length(history_t) - 1)
    scaled_y = variable_factors .* y
    for candidate_order in 1:maximum_candidate
        points = candidate_order + 1
        predictor = extrapolate(tnew, history_t[1:points], history_y[1:points])
        scaled_predictor = variable_factors .* predictor
        scale = component_scale(scaled_y, scaled_predictor, options)
        errors[candidate_order] = inv(T(candidate_order + 1)) * scaled_rms(
            scaled_y - scaled_predictor, scale, error_mask)
    end
    return errors
end

function attempt_step(residual!, jacobian!, jacobian_prototype, parameter,
        history_t, history_y,
        h, order, options, error_mask, monitor_mask, variable_levels,
        equation_levels, stats, factorization_cache, predictor!)
    T = eltype(first(history_y))
    tnew = first(history_t) + h
    usable_order = min(order, length(history_t), options.maximum_order)

    # Use one more history point for extrapolation when available. The
    # predictor-corrector difference then has the local-error order k+1.
    predictor_points = min(usable_order + 1, length(history_t))
    ypred = extrapolate(tnew, history_t[1:predictor_points],
        history_y[1:predictor_points])

    nodes = [tnew; history_t[1:usable_order]]
    weights = derivative_weights(tnew, nodes)
    cj = first(weights)
    history_derivative = zeros(T, length(ypred))
    for j in 1:usable_order
        history_derivative .+= weights[j + 1] .* history_y[j]
    end
    isnothing(predictor!) || predictor!(ypred, tnew, cj,
        history_derivative, parameter)

    residual = similar(ypred)
    trial = similar(ypred)
    matrix = isnothing(jacobian_prototype) ?
        Matrix{T}(undef, length(ypred), length(ypred)) :
        copy(jacobian_prototype)
    variable_factors = isnothing(variable_levels) ? ones(T, length(ypred)) :
        level_factors(variable_levels, h)
    equation_factors = isnothing(equation_levels) ? ones(T, length(ypred)) :
        level_factors(equation_levels, h)
    scaled_coefficient = isnothing(variable_levels) ||
        isnothing(equation_levels) ? abs(cj) : abs(cj * h)
    if coefficient_requires_refactorization(factorization_cache,
            scaled_coefficient,
            options.factorization_coefficient_tolerance)
        invalidate_numeric_factorization!(factorization_cache)
    end
    if factorization_cache.numeric_valid &&
            factorization_cache.steps >= options.maximum_factorization_steps
        invalidate_numeric_factorization!(factorization_cache)
    end

    # A cached matrix is used as a chord-iteration matrix. If it cannot
    # converge the corrector, restart from the predictor once with a newly
    # evaluated and numerically factored matrix before rejecting the step.
    reuse_first_matrix = factorization_cache.numeric_valid
    phases = reuse_first_matrix ? 2 : 1
    total_iterations = 0
    last_y = copy(ypred)
    last_yprime = cj .* last_y .+ history_derivative
    last_residual_norm = T(Inf)
    last_reason = :newton_iteration_limit

    for phase in 1:phases
        y = copy(ypred)
        yprime = cj .* y .+ history_derivative
        previous_correction_norm = T(Inf)
        factorization = factorization_cache.numeric_valid ?
            factorization_cache.factorization : nothing
        using_reused_matrix = phase == 1 && reuse_first_matrix
        using_reused_matrix && (factorization_cache.steps += 1)
        correction_multiplier = using_reused_matrix ?
            T(2) / (one(T) + scaled_coefficient /
                factorization_cache.scaled_coefficient) : one(T)
        # A chord iteration is inexpensive, so carry it two decimal decades
        # beyond the nominal test and retain fresh-Newton equation accuracy.
        convergence_limit = using_reused_matrix ? T(0.01) : one(T)
        failure_reason = :newton_iteration_limit
        residual_norm = T(Inf)

        for iteration in 1:options.maximum_newton_iterations
            residual!(residual, tnew, y, yprime, parameter)
            stats.residual_evaluations += 1
            scaled_residual = equation_factors .* residual
            residual_scale = options.residual_tolerance .* max.(
                variable_factors .* abs.(y), one(T))
            residual_norm = sqrt(sum(abs2,
                scaled_residual ./ residual_scale) / length(y))

            if isnothing(factorization)
                evaluate_iteration_matrix!(matrix, jacobian!, residual!,
                    residual, trial, tnew, y, yprime, cj, parameter, stats)
                try
                    scaled_matrix = scaled_newton_matrix(copy(matrix),
                        equation_factors, variable_factors)
                    factorization = factor_iteration_matrix!(
                        factorization_cache, scaled_matrix,
                        scaled_coefficient, stats)
                catch exception
                    exception isa SingularException || rethrow()
                    invalidate_numeric_factorization!(factorization_cache)
                    return StepAttempt(false, y, yprime, T(Inf),
                        residual_norm, total_iterations + iteration,
                        :singular_iteration_matrix, T[], T(NaN), T(NaN), 0)
                end
            end

            scaled_correction = try
                factorization \ scaled_residual
            catch exception
                exception isa SingularException || rethrow()
                invalidate_numeric_factorization!(factorization_cache)
                return StepAttempt(false, y, yprime, T(Inf),
                    residual_norm, total_iterations + iteration,
                    :singular_iteration_matrix, T[], T(NaN), T(NaN), 0)
            end
            scaled_correction .*= correction_multiplier
            correction = scaled_correction ./ variable_factors
            stats.newton_iterations += 1
            total_iterations += 1
            y .-= correction
            yprime .-= cj .* correction

            correction_scale = options.newton_atol .+
                options.newton_rtol .* max.(
                    variable_factors .* abs.(y), one(T))
            correction_norm = sqrt(sum(abs2,
                scaled_correction ./ correction_scale) / length(correction))

            if correction_norm <= convergence_limit
                residual!(residual, tnew, y, yprime, parameter)
                stats.residual_evaluations += 1
                scaled_residual = equation_factors .* residual
                residual_scale = options.residual_tolerance .* max.(
                    variable_factors .* abs.(y), one(T))
                residual_norm = sqrt(sum(abs2,
                    scaled_residual ./ residual_scale) / length(y))
                if residual_norm <= convergence_limit
                    order_errors = candidate_order_errors(y, tnew,
                        history_t, history_y, options, variable_factors,
                        error_mask)
                    if !isfinite(order_errors[usable_order])
                        scaled_y = variable_factors .* y
                        scaled_predictor = variable_factors .* ypred
                        scale = component_scale(scaled_y, scaled_predictor,
                            options)
                        order_errors[usable_order] =
                            inv(T(usable_order + 1)) * scaled_rms(
                                scaled_y - scaled_predictor, scale, error_mask)
                    end
                    error = order_errors[usable_order]
                    monitored_rms_error, monitored_maximum_error,
                        monitored_maximum_index = isnothing(monitor_mask) ?
                        (T(NaN), T(NaN), 0) : monitored_error_summary(
                            y, ypred, options, variable_factors,
                            monitor_mask, usable_order)
                    return StepAttempt(true, y, yprime, error, residual_norm,
                        total_iterations, :success, order_errors,
                        monitored_rms_error, monitored_maximum_error,
                        monitored_maximum_index)
                end
            end

            if iteration > 1
                rate = correction_norm / previous_correction_norm
                if !isfinite(rate) || rate > T(0.9)
                    failure_reason = :newton_stagnation
                    break
                end
            end
            previous_correction_norm = correction_norm
        end

        last_y = y
        last_yprime = yprime
        last_residual_norm = residual_norm
        last_reason = failure_reason
        if using_reused_matrix
            invalidate_numeric_factorization!(factorization_cache)
        else
            return StepAttempt(false, last_y, last_yprime, T(Inf),
                last_residual_norm, total_iterations, last_reason, T[],
                T(NaN), T(NaN), 0)
        end
    end

    return StepAttempt(false, last_y, last_yprime, T(Inf),
        last_residual_norm, total_iterations, last_reason, T[],
        T(NaN), T(NaN), 0)
end

function initial_step(t0, tf, y, yprime, options)
    options.initial_step > 0 && return min(options.initial_step,
        options.maximum_step, abs(tf - t0))
    scale = component_scale(y, y, options)
    derivative_norm = sqrt(sum(abs2, yprime ./ scale) / length(y))
    estimate = derivative_norm > 0 ? 0.01 / derivative_norm : 0.001 * abs(tf - t0)
    return min(max(estimate, 100eps(eltype(y)) * max(abs(t0), one(t0))),
        options.maximum_step, abs(tf - t0))
end

"""
    dassl(residual!, y0, yprime0, (t0, tf); kwargs...)

Integrate the square implicit system `G(t, y, yprime) = 0` using a compact,
modern reconstruction of DDASSL's variable-step BDF method.

The callback is

```
residual!(out, t, y, yprime, parameter)
```

An optional analytical iteration-matrix callback has the form

```
jacobian!(matrix, t, y, yprime, cj, parameter)
```

and must return `G_y + cj * G_yprime`. If it is omitted, the combined
matrix is formed by directional finite differences, as in DDASSL.

Set `jacobian_prototype` to a matrix with the required storage type and
structural pattern. The integrator copies this prototype for each step and
passes it to `jacobian!`. A sparse prototype therefore selects sparse numerical
factorization while preserving the model's supplied nonzero structure.
The numerical factors are also retained across steps and used for modified
Newton corrections while the level-scaled BDF coefficient remains nearby.
`factorization_coefficient_tolerance` controls that neighborhood; its default
of 0.25 gives coefficient-ratio bounds of 0.6 and 5/3, following DDASSL. If an
old matrix has served `maximum_factorization_steps` attempts, five by default,
it is refreshed. If it cannot converge a corrector, the step is retried once
from its predictor with a new numerical factorization before the step is rejected.
UMFPACK's symbolic analysis remains available for that refactorization while
the sparsity structure is unchanged. During one call, the symbolic analysis is
discarded after `symbolic_reanalysis_failures` consecutive corrector failures,
whose default is two.

Set `differential_vars` to classify variables whose derivatives appear in the
implicit equations. Together with `variable_levels` and `deficit`, it defines
the default integration-error mask. An explicit Boolean `error_control` vector
overrides that default while retaining every variable in the Newton solve.

An optional Boolean `error_monitor` vector records a separate predictor-error
history without affecting step acceptance. Its RMS and maximum component
values use the same tolerances and level scaling as the integration-error test.
For example, a level-one velocity difference is multiplied by the current step
size so it measures the implied position error over that step.

Set `root!` and `number_of_roots` to locate discontinuity surfaces. The callback
has the same arguments as `residual!` and writes one scalar function per
surface. `root_directions` may contain `-1`, `0`, or `1` to select decreasing,
either, or increasing crossings along the direction of integration.
`root_restart` selects one response to every located root: `:hard` discards
history, returns to order one, and quarters the next step; `:soft` retains
history, caps the order at two, and halves the next step; and `:none` retains
both history and order. Every policy refreshes the numeric iteration matrix.
The conservative generic default is `:hard`; callers that know their state and
force remain continuous may select `:soft`.

An optional `reconfigure!` callback can replace an equivalent set of implicit
equations without changing the variables or system size. It is called after a
severe monitored-error episode, a singular iteration matrix, or repeated
corrector failures. The callback receives
`(reason, t, y, yprime, error_control, differential_vars, parameter)`. It may
modify the two masks and return a named tuple containing a replacement
`jacobian_prototype` and `equation_levels`. Returning `nothing` leaves the
formulation unchanged. A proactive change made after an accepted step preserves
the complete BDF history. A recovery change after a failed attempted step
returns to the last accepted state, discards the now-questionable polynomial
history, and resumes at order one. Both force a new symbolic and numerical
factorization.

An optional `predictor!` callback may adjust auxiliary predicted values before
the corrector. It receives `(ypred, tnew, cj, history_derivative, parameter)`;
after it returns, the integrator evaluates `yprime = cj * ypred +
history_derivative`. This is intended for redundant internal coordinates whose
predicted values may be changed without changing the physical state.

An optional `accepted_step!` callback receives
`(time, y, yprime, order, step, stats)` after each accepted integration step.
It is intended for progress displays and diagnostics. The arrays are the
integrator's current values; a callback that retains them must make copies.

An optional `accepted_interval!` callback receives
`(left_time, right_time, nodes, values, order, step, stats)` after each
accepted step. `nodes` and `values` define the corrected BDF history
polynomial over that step. This permits output samples to be written as soon
as their times are crossed, using the same interpolation as the completed
solution.

An optional `accepted_step_filter!` callback receives the same arguments just
before an accepted state is stored. It may modify `y` and `yprime`. This is a
specialized hook for first-order dynamic relaxation; altering an accepted state
while retaining a higher-order history is not generally valid.

Set `save_everystep=false` to retain only the initial and latest accepted
states in the returned result. The internal BDF history remains unchanged and
bounded by `maximum_order + 1`. This mode is intended for callers that consume
requested output through `accepted_interval!`; dense interpolation over the
discarded integration interval is then unavailable from the returned result.

An optional `accepted_diagnostics!` callback receives `(time, y, yprime,
order, step, stats, controlled_error, monitored_rms_error,
monitored_maximum_error, monitored_maximum_index)` after each accepted step.
It permits bounded diagnostic collectors to retain selected events without
keeping the complete integration trace.
"""
function dassl(residual!, y0::AbstractVector{T},
        yprime0::AbstractVector{T}, tspan::Tuple{T,T};
        parameter = nothing, jacobian! = nothing, jacobian_prototype = nothing,
        options = DASSLOptions{T}(),
        error_control = nothing, error_monitor = nothing, variable_levels = nothing,
        equation_levels = nothing, differential_vars = nothing,
        deficit::Integer = 0, root! = nothing, number_of_roots::Integer = 0,
        root_directions = nothing, root_restart::Symbol = :hard,
        reconfigure! = nothing, predictor! = nothing,
        accepted_step! = nothing,
        accepted_interval! = nothing,
        accepted_step_filter! = nothing,
        accepted_diagnostics! = nothing,
        save_everystep::Bool = true) where {T<:AbstractFloat}
    length(y0) == length(yprime0) ||
        throw(DimensionMismatch("y0 and yprime0 must have equal lengths"))
    if !isnothing(jacobian_prototype)
        size(jacobian_prototype) == (length(y0), length(y0)) ||
            throw(DimensionMismatch(
                "jacobian_prototype must match the implicit system size"))
        isnothing(jacobian!) && throw(ArgumentError(
            "jacobian_prototype requires an analytical jacobian! callback"))
    end
    deficit >= 0 || throw(ArgumentError("deficit must be nonnegative"))
    if !isnothing(variable_levels)
        length(variable_levels) == length(y0) || throw(DimensionMismatch(
            "variable_levels must match y0"))
        all(level -> level >= 0, variable_levels) || throw(ArgumentError(
            "variable levels must be nonnegative"))
        variable_levels = collect(Int, variable_levels)
    end
    if !isnothing(equation_levels)
        length(equation_levels) == length(y0) || throw(DimensionMismatch(
            "equation_levels must match the residual equation count"))
        all(level -> level >= 0, equation_levels) || throw(ArgumentError(
            "equation levels must be nonnegative"))
        equation_levels = collect(Int, equation_levels)
    end
    if !isnothing(differential_vars)
        length(differential_vars) == length(y0) || throw(DimensionMismatch(
            "differential_vars must match y0"))
        differential_vars = BitVector(differential_vars)
    end
    if isnothing(error_control)
        if isnothing(variable_levels)
            error_control = isnothing(differential_vars) ? trues(length(y0)) :
                copy(differential_vars)
        else
            default_differential = isnothing(differential_vars) ?
                trues(length(y0)) : differential_vars
            error_control = error_control_from_levels(variable_levels, deficit;
                differential_vars = default_differential)
        end
    else
        length(error_control) == length(y0) ||
            throw(DimensionMismatch("error_control must match y0"))
        error_control = BitVector(error_control)
    end
    any(error_control) || throw(ArgumentError(
        "at least one variable must participate in integration error control"))
    if !isnothing(error_monitor)
        length(error_monitor) == length(y0) || throw(DimensionMismatch(
            "error_monitor must match y0"))
        error_monitor = BitVector(error_monitor)
        any(error_monitor) || throw(ArgumentError(
            "at least one variable must participate in error monitoring"))
    end
    1 <= options.maximum_order <= 5 || throw(ArgumentError(
        "maximum_order must lie between one and five"))
    options.order_change_advantage >= one(T) || throw(ArgumentError(
        "order_change_advantage must be at least one"))
    options.minimum_order_steps >= 1 || throw(ArgumentError(
        "minimum_order_steps must be positive"))
    options.maximum_factorization_steps >= 1 || throw(ArgumentError(
        "maximum_factorization_steps must be positive"))
    zero(T) <= options.factorization_coefficient_tolerance < one(T) ||
        throw(ArgumentError(
            "factorization_coefficient_tolerance must lie in [0, 1)"))
    options.symbolic_reanalysis_failures >= 1 || throw(ArgumentError(
        "symbolic_reanalysis_failures must be positive"))
    number_of_roots >= 0 || throw(ArgumentError(
        "number_of_roots must be nonnegative"))
    root_restart in (:none, :soft, :hard) || throw(ArgumentError(
        "root_restart must be :none, :soft, or :hard"))
    isnothing(root!) == iszero(number_of_roots) || throw(ArgumentError(
        "root! and a positive number_of_roots must be supplied together"))
    if isnothing(root_directions)
        root_directions = zeros(Int, number_of_roots)
    else
        length(root_directions) == number_of_roots || throw(DimensionMismatch(
            "root_directions must match number_of_roots"))
        all(direction -> direction in (-1, 0, 1), root_directions) ||
            throw(ArgumentError("root directions must be -1, 0, or 1"))
        root_directions = collect(Int, root_directions)
    end

    t0, tf = tspan
    direction = sign(tf - t0)
    events = DASSLEvent{T,Vector{T}}[]
    monitor = isnothing(error_monitor) ? nothing : DASSLErrorMonitor(
        error_monitor, T[T(NaN)], T[T(NaN)], T[T(NaN)], Int[0])
    direction == 0 && return DASSLResult(T[t0], [copy(y0)],
        [copy(yprime0)], Int[1], T[zero(T)], DASSLStats(),
        SciMLBase.ReturnCode.Success, "initial and final times are equal",
        variable_levels, equation_levels, differential_vars, error_control,
        monitor, Int(deficit), events)

    stats = DASSLStats()
    t_values = T[t0]
    y_values = [copy(y0)]
    yprime_values = [copy(yprime0)]
    orders = Int[1]
    steps = T[zero(T)]
    history_t = T[t0]
    history_y = [copy(y0)]
    order = 1
    steps_at_order = 0
    previous_accepted_order = 0
    previous_accepted_step = zero(T)
    consecutive_failures = 0
    consecutive_corrector_failures = 0
    health_episode_active = false
    factorization_cache = LinearFactorizationCache(T)
    h = direction * initial_step(t0, tf, y0, yprime0, options)
    minimum_step = options.minimum_step > 0 ? options.minimum_step :
        100eps(T) * max(abs(t0), abs(tf), one(T))
    root_time_tolerance = 100eps(T) * max(abs(t0), abs(tf), one(T))
    left_roots = zeros(T, number_of_roots)
    !isnothing(root!) && evaluate_roots!(root!, left_roots, t0, y0, yprime0,
        parameter, stats)

    residual_check = similar(y0)
    residual!(residual_check, t0, y0, yprime0, parameter)
    stats.residual_evaluations += 1
    initial_residual_scale = options.residual_tolerance .* max.(abs.(y0), one(T))
    initial_residual_norm = sqrt(sum(abs2,
        residual_check ./ initial_residual_scale) / length(y0))
    if initial_residual_norm > one(T)
        return DASSLResult(t_values, y_values, yprime_values, orders, steps,
            stats, SciMLBase.ReturnCode.InitialFailure,
            "the supplied initial y and yprime do not satisfy the residual",
            variable_levels, equation_levels, differential_vars,
            error_control, monitor, Int(deficit), events)
    end

    function reconfigure_system!(reason, time, values, derivatives)
        isnothing(reconfigure!) && return false
        update = reconfigure!(reason, time, values, derivatives,
            error_control, differential_vars, parameter)
        isnothing(update) && return false
        hasproperty(update, :jacobian_prototype) || throw(ArgumentError(
            "reconfigure! must return jacobian_prototype"))
        hasproperty(update, :equation_levels) || throw(ArgumentError(
            "reconfigure! must return equation_levels"))
        replacement = update.jacobian_prototype
        size(replacement) == (length(y0), length(y0)) ||
            throw(DimensionMismatch(
                "replacement jacobian_prototype has the wrong size"))
        levels = collect(Int, update.equation_levels)
        length(levels) == length(y0) || throw(DimensionMismatch(
            "replacement equation_levels have the wrong length"))
        jacobian_prototype = replacement
        equation_levels = levels
        invalidate_symbolic_factorization!(factorization_cache)
        stats.state_reselections += 1
        true
    end

    while direction * (tf - first(history_t)) > zero(T)
        stats.accepted_steps + stats.rejected_steps < options.maximum_steps ||
            return DASSLResult(t_values, y_values, yprime_values, orders,
                steps, stats, SciMLBase.ReturnCode.MaxIters,
                "maximum step count reached", variable_levels,
                equation_levels, differential_vars, error_control, monitor,
                Int(deficit), events)

        remaining = tf - first(history_t)
        h = direction * min(abs(h), abs(remaining), options.maximum_step)
        abs(h) >= minimum_step || return DASSLResult(t_values, y_values,
            yprime_values, orders, steps, stats,
            SciMLBase.ReturnCode.DtLessThanMin,
            "step size fell below the minimum", variable_levels,
            equation_levels, differential_vars, error_control, monitor,
            Int(deficit), events)

        attempt = attempt_step(residual!, jacobian!, jacobian_prototype,
            parameter, history_t,
            history_y, h, order, options, error_control, error_monitor,
            variable_levels, equation_levels, stats, factorization_cache,
            predictor!)

        if attempt.converged && attempt.error <= one(T)
            tnew = first(history_t) + h
            event = nothing
            right_roots = similar(left_roots)
            if !isnothing(root!)
                evaluate_roots!(root!, right_roots, tnew, attempt.y,
                    attempt.yprime, parameter, stats)
                usable_order = min(order, length(history_t),
                    options.maximum_order)
                interpolation_nodes = [tnew; history_t[1:usable_order]]
                interpolation_values = vcat([copy(attempt.y)],
                    history_y[1:usable_order])
                event = find_first_event(root!, first(history_t), tnew,
                    left_roots, right_roots, root_directions,
                    interpolation_nodes, interpolation_values, parameter,
                    stats, root_time_tolerance)
            end

            accepted_time = isnothing(event) ? tnew : event.time
            accepted_y = isnothing(event) ? attempt.y : event.y
            accepted_yprime = isnothing(event) ? attempt.yprime : event.yprime
            accepted_step = accepted_time - first(history_t)
            isnothing(accepted_step_filter!) || accepted_step_filter!(
                accepted_time, accepted_y, accepted_yprime, order,
                accepted_step, stats)
            if save_everystep || length(t_values) == 1
                push!(t_values, accepted_time)
                push!(y_values, copy(accepted_y))
                push!(yprime_values, copy(accepted_yprime))
                push!(orders, order)
                push!(steps, accepted_step)
            else
                t_values[end] = accepted_time
                y_values[end] .= accepted_y
                yprime_values[end] .= accepted_yprime
                orders[end] = order
                steps[end] = accepted_step
            end
            controlled_error = isnothing(event) ? attempt.error : T(NaN)
            monitored_rms_error = isnothing(event) ?
                attempt.monitored_rms_error : T(NaN)
            monitored_maximum_error = isnothing(event) ?
                attempt.monitored_maximum_error : T(NaN)
            monitored_maximum_index = isnothing(event) ?
                attempt.monitored_maximum_index : 0
            if !isnothing(monitor)
                if save_everystep || length(monitor.controlled_errors) == 1
                    push!(monitor.controlled_errors, controlled_error)
                    push!(monitor.rms_errors, monitored_rms_error)
                    push!(monitor.maximum_errors, monitored_maximum_error)
                    push!(monitor.maximum_error_indices,
                        monitored_maximum_index)
                else
                    monitor.controlled_errors[end] = controlled_error
                    monitor.rms_errors[end] = monitored_rms_error
                    monitor.maximum_errors[end] = monitored_maximum_error
                    monitor.maximum_error_indices[end] =
                        monitored_maximum_index
                end
            end
            stats.accepted_steps += 1
            isnothing(accepted_step!) || accepted_step!(accepted_time,
                accepted_y, accepted_yprime, order, accepted_step, stats)
            isnothing(accepted_diagnostics!) || accepted_diagnostics!(
                accepted_time, accepted_y, accepted_yprime, order,
                accepted_step, stats, controlled_error,
                monitored_rms_error, monitored_maximum_error,
                monitored_maximum_index)
            if !isnothing(accepted_interval!)
                usable_order = min(order, length(history_t),
                    options.maximum_order)
                interpolation_nodes = [accepted_time;
                    history_t[1:usable_order]]
                interpolation_values = vcat([copy(accepted_y)],
                    history_y[1:usable_order])
                accepted_interval!(first(history_t), accepted_time,
                    interpolation_nodes, interpolation_values, order,
                    accepted_step, stats)
            end
            consecutive_failures = 0
            consecutive_corrector_failures = 0

            # The sparse factorization already estimates the reciprocal
            # condition of the scaled iteration matrix. A substantial loss
            # of condition is a cheap warning that an equivalent mechanical
            # state partition may be deteriorating. Let the model compare
            # alternatives while this accepted state is still accurate.
            if isnothing(event) && conditioning_warning(factorization_cache)
                changed = reconfigure_system!(
                    :deteriorating_iteration_matrix, accepted_time,
                    accepted_y, accepted_yprime)
                if changed
                    # Every physical variable retains its own accepted BDF
                    # history. Only the equivalent state-equation rows have
                    # changed, so no coordinate transfer or history restart
                    # is required.
                    health_episode_active = false
                else
                    # No better partition was available at this state. Use
                    # the present condition as the reference so another QR
                    # comparison requires a further substantial deterioration.
                    factorization_cache.condition_reference =
                        factorization_cache.reciprocal_condition
                end
            end

            if !isnothing(event)
                push!(events, event)
                stats.events_found += 1
                health_episode_active = false
                previous_accepted_order = 0
                previous_accepted_step = zero(T)
                # A root commonly marks a discontinuous force or stiffness.
                # Preserve the sparse ordering but refresh its numeric values.
                invalidate_numeric_factorization!(factorization_cache)
                # Keep the sign on the far side of the crossing. This avoids
                # reporting the same zero again when the next step starts at it.
                left_roots .= right_roots
                if root_restart == :hard
                    history_t = T[accepted_time]
                    history_y = [copy(accepted_y)]
                    order = 1
                    steps_at_order = 0
                    stats.history_restarts += 1
                    h *= T(0.25)
                else
                    pushfirst!(history_t, accepted_time)
                    pushfirst!(history_y, copy(accepted_y))
                    resize!(history_t,
                        min(length(history_t), options.maximum_order + 1))
                    resize!(history_y,
                        min(length(history_y), options.maximum_order + 1))
                end
                if root_restart == :soft
                    # A compliant-force transition leaves the state and force
                    # continuous while changing a force derivative. Preserve
                    # useful history, but do not carry a high-order formula
                    # across that loss of smoothness.
                    order = min(order, 2)
                    steps_at_order = 0
                    stats.history_restarts += 1
                    h *= T(0.5)
                end
                continue
            end

            pushfirst!(history_t, tnew)
            pushfirst!(history_y, copy(attempt.y))
            resize!(history_t, min(length(history_t), options.maximum_order + 1))
            resize!(history_y, min(length(history_y), options.maximum_order + 1))
            !isnothing(root!) && (left_roots .= right_roots)
            steps_at_order += 1

            stable_formula = previous_accepted_order != 0 && begin
                previous_step = abs(previous_accepted_step)
                ratio = iszero(previous_step) ? one(T) :
                    abs(accepted_step) / previous_step
                order == previous_accepted_order &&
                    T(2 / 3) <= ratio <= T(3 / 2)
            end
            amplification = attempt.monitored_maximum_error /
                max(attempt.error, T(0.1))
            unhealthy = stable_formula &&
                attempt.monitored_maximum_error > T(25) &&
                amplification > T(10)
            if unhealthy && !health_episode_active
                reconfigure_system!(:high_physical_error, tnew,
                    attempt.y, attempt.yprime)
                health_episode_active = true
            elseif !unhealthy
                health_episode_active = false
            end
            previous_accepted_order = order
            previous_accepted_step = accepted_step

            current_order = order
            current_factor = proposed_step_factor(
                attempt.order_errors[current_order], current_order, options)
            selected_order = current_order
            selected_factor = current_factor
            can_change_order = steps_at_order >=
                max(options.minimum_order_steps, current_order + 1)
            if can_change_order
                candidate_orders = UnitRange(max(1, current_order - 1),
                    min(options.maximum_order, current_order + 1))
                for candidate in candidate_orders
                    candidate == current_order && continue
                    candidate > length(attempt.order_errors) && continue
                    candidate_error = attempt.order_errors[candidate]
                    isfinite(candidate_error) || continue
                    candidate_factor = proposed_step_factor(
                        candidate_error, candidate, options)
                    if candidate_factor >
                            options.order_change_advantage * selected_factor
                        selected_order = candidate
                        selected_factor = candidate_factor
                    end
                end
            end
            if selected_order != current_order
                order = selected_order
                steps_at_order = 0
            end
            factor = selected_factor
            h *= factor
        else
            stats.rejected_steps += 1
            consecutive_failures += 1
            if attempt.converged
                consecutive_corrector_failures = 0
            else
                invalidate_numeric_factorization!(factorization_cache)
                stats.corrector_failures += 1
                consecutive_corrector_failures += 1
                recovery_reason = attempt.reason == :singular_iteration_matrix ?
                    :singular_iteration_matrix : :repeated_corrector_failure
                should_reconfigure =
                    attempt.reason == :singular_iteration_matrix ||
                    consecutive_corrector_failures >=
                        options.symbolic_reanalysis_failures
                if should_reconfigure && reconfigure_system!(recovery_reason,
                        first(history_t), first(history_y),
                        last(yprime_values))
                    resize!(history_t, 1)
                    resize!(history_y, 1)
                    order = 1
                    steps_at_order = 0
                    stats.history_restarts += 1
                    h *= T(0.25)
                    consecutive_failures = 0
                    consecutive_corrector_failures = 0
                    continue
                elseif consecutive_corrector_failures >=
                        options.symbolic_reanalysis_failures
                    invalidate_symbolic_factorization!(factorization_cache)
                    consecutive_corrector_failures = 0
                end
            end
            factor = attempt.converged && isfinite(attempt.error) ?
                proposed_step_factor(attempt.error, order, options;
                    upper = T(0.8)) : T(0.25)
            if consecutive_failures >= 2 && order > 1
                order -= 1
                steps_at_order = 0
                consecutive_failures = 0
            end
            h *= factor
        end
    end

    return DASSLResult(t_values, y_values, yprime_values, orders, steps,
        stats, SciMLBase.ReturnCode.Success,
        "integration reached the final time", variable_levels,
        equation_levels, differential_vars, error_control, monitor,
        Int(deficit), events)
end

end # module
