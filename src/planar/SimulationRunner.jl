"""
    SimulationRunner

Analysis dispatcher for a loaded planar model. It provides simultaneous
consistent initialization, dynamic or kinematic implicit integration, static
equilibrium and quasi-static continuation, optional dynamic relaxation, and
modal operating-point preparation.

Dynamic integration always advances the full active canonical vector. A
`RuntimeStatePartition` changes only which candidate coordinate equations and
error-control entries define the minimal state set; therefore a recovery
selection can preserve the accepted BDF history and all physical variables.
"""
module SimulationRunner

using LinearAlgebra
using SciMLBase
using ..AutomaticAnalysis
using ..HistoricalDDASSL
using ..ModalAnalysis
using ..PlanarComponentAssembly
using ..PlanarEquationComponents
using ..PlanarModelIO

export run_planar_model

struct StaticConvergenceError <: Exception
    message::String
end

Base.showerror(io::IO, error::StaticConvergenceError) = print(io, error.message)

"""One acceleration, velocity, position triple eligible as a physical state."""
struct RuntimeStateCandidate
    name::Symbol
    acceleration::Int
    velocity::Int
    position::Int
    equations::Vector{Int}
end

"""
Mutable equation partition used for runtime state recovery.

`base_equations` never changes. `equation_indices` appends the defining
equations for `selected_names` and is rebuilt in place when QR chooses an
equivalent state set.
"""
mutable struct RuntimeStatePartition
    candidates::Dict{Symbol,RuntimeStateCandidate}
    base_equations::Vector{Int}
    equation_indices::Vector{Int}
    selected_names::Vector{Symbol}
end

"""Collect every allocated body or relative-coordinate state candidate."""
function runtime_state_candidates(loaded)
    candidates = RuntimeStateCandidate[]
    for body in values(loaded.bodies)
        for kind in (:V_x, :V_y, :omega)
            haskey(body.candidate_state_equations, kind) || continue
            acceleration, velocity, position =
                body.candidate_state_variables[kind]
            push!(candidates, RuntimeStateCandidate(
                Symbol(body.name, :., kind), acceleration, velocity, position,
                collect(body.candidate_state_equations[kind])))
        end
    end
    for connection in values(loaded.connections)
        hasproperty(connection, :candidate_state_equations) || continue
        isempty(connection.candidate_state_equations) && continue
        if connection isa PlanarRevoluteJointComponent
            acceleration, velocity, position = connection.rotation_variables
            name = Symbol(connection.name, :., :omega)
        elseif connection isa PlanarDistanceCoordinateComponent
            acceleration = connection.acceleration_variable
            velocity = connection.velocity_variable
            position = connection.distance_variable
            name = Symbol(connection.name, :., :velocity)
        else
            continue
        end
        push!(candidates, RuntimeStateCandidate(name, acceleration, velocity,
            position, copy(connection.candidate_state_equations)))
    end
    sort!(candidates; by = candidate -> candidate.velocity)
    Dict(candidate.name => candidate for candidate in candidates)
end

"""
Build the initial dynamic row partition from the loader's selected velocities.
Rows belonging to dormant candidates are omitted but remain allocated in the
canonical model.
"""
function runtime_state_partition(loaded)
    candidates = runtime_state_candidates(loaded)
    candidate_equations = Set(Iterators.flatten(
        candidate.equations for candidate in values(candidates)))
    base_equations = [index for index in loaded.active_equation_indices
        if index ∉ candidate_equations]
    selected_names = copy(loaded.state_selection.selected_velocities)
    selected_equations = collect(Iterators.flatten(
        candidates[name].equations for name in selected_names))
    RuntimeStatePartition(candidates, base_equations,
        [base_equations; selected_equations], selected_names)
end

"""Replace candidate state rows without changing the number of freedoms."""
function select_runtime_states!(partition, names)
    length(names) == length(partition.selected_names) ||
        throw(DimensionMismatch("runtime state selection changed the number of states"))
    all(name -> haskey(partition.candidates, name), names) ||
        throw(ArgumentError("runtime state selection named an unknown candidate"))
    partition.selected_names = collect(names)
    empty!(partition.equation_indices)
    append!(partition.equation_indices, partition.base_equations)
    for name in partition.selected_names
        append!(partition.equation_indices,
            partition.candidates[name].equations)
    end
    partition
end

function runtime_state_variables(partition)
    selected = [partition.candidates[name] for name in partition.selected_names]
    getproperty.(selected, :position), getproperty.(selected, :velocity)
end

function automatic_runtime_state_selection(loaded, state, time)
    body_names = sort!(collect(keys(loaded.bodies));
        by = name -> loaded.bodies[name].velocity_variables[1])
    connection_names = collect(keys(loaded.connections))
    PlanarModelIO.automatic_velocity_selection(loaded.model, loaded.bodies,
        body_names, loaded.markers, loaded.connections, connection_names,
        loaded.state_selection.body_characteristic_lengths, state, time)
end

function selected_state_variables(loaded)
    orientations = Int[]
    velocities = Int[]
    for components in (values(loaded.bodies), values(loaded.connections))
        for component in components
            hasproperty(component, :selected_state_variables) || continue
            for (_, velocity, position) in component.selected_state_variables
                push!(orientations, position)
                push!(velocities, velocity)
            end
        end
    end
    orientations, velocities
end

"""
Construct `zdot` from the explicit position/velocity/acceleration levels.

This is not a separate dynamics calculation. Acceleration initialization has
already made the supplied canonical state consistent, and this routine only
expresses its level relationships in the form required by DASSL.
"""
function initial_derivative(state, loaded,
        time = loaded.simulation.start_time)
    derivative = zeros(length(state))
    for body in values(loaded.bodies)
        derivative[body.position_variables] .= state[body.velocity_variables]
        derivative[body.orientation_variable] = state[body.angular_velocity_variable]
        derivative[body.velocity_variables] .= state[body.acceleration_variables]
        derivative[body.angular_velocity_variable] =
            state[body.angular_acceleration_variable]
    end
    for driver in values(loaded.drivers)
        if driver isa PlanarRotationalMotionGenerator
            derivative[driver.angle_variable] =
                state[driver.angular_velocity_variable]
            derivative[driver.angular_velocity_variable] =
                state[driver.angular_acceleration_variable]
        else
            derivative[driver.distance_variable] = state[driver.velocity_variable]
            derivative[driver.velocity_variable] =
                state[driver.acceleration_variable]
        end
    end
    for connection in values(loaded.connections)
        if connection isa PlanarRevoluteJointComponent &&
                !isempty(connection.rotation_variables)
            alpha, omega, theta = connection.rotation_variables
            derivative[theta] = state[omega]
            derivative[omega] = state[alpha]
        elseif connection isa PlanarDistanceCoordinateComponent
            derivative[connection.distance_variable] =
                state[connection.velocity_variable]
            derivative[connection.velocity_variable] =
                state[connection.acceleration_variable]
        end
    end
    for component in values(loaded.equation_components)
        planar_equation_state_rates!(derivative, component,
            time, state)
    end
    derivative
end

"""Simultaneously correct all initial variables while holding selected states."""
function initialize_implicit_model!(state, loaded, time;
        tolerance = 1.0e-10, maximum_iterations = 20)
    selection = AnalysisSelection(Dynamics(), loaded.active_variable_indices,
        loaded.active_equation_indices)
    state_positions, state_velocities = selected_state_variables(loaded)
    user_states = collect(Iterators.flatten(component.state_indices
        for component in values(loaded.equation_components)))
    fixed = Set([state_positions; state_velocities; user_states])
    algebraic = [index for index in loaded.active_variable_indices
                 if index ∉ fixed]
    differential = [state_positions; state_velocities; user_states]
    derivative = initial_derivative(state, loaded, time)
    equations = zeros(eltype(state), length(loaded.active_equation_indices))
    for iteration in 0:maximum_iterations
        evaluate_analysis_equations!(equations, loaded.model, selection,
            time, state, derivative)
        norm(equations, Inf) <= tolerance &&
            return initial_derivative(state, loaded, time), iteration
        iteration == maximum_iterations && break
        algebraic_selection = AnalysisSelection(Dynamics(), algebraic,
            loaded.active_equation_indices)
        J0 = evaluate_analysis_sparse_jacobian(loaded.model, algebraic_selection,
            time, state, derivative, 0.0)
        matrix = if isempty(differential)
            J0
        else
            differential_selection = AnalysisSelection(Dynamics(), differential,
                loaded.active_equation_indices)
            Jdot0 = evaluate_analysis_sparse_jacobian(loaded.model,
                differential_selection, time, state, derivative, 0.0)
            J1 = evaluate_analysis_sparse_jacobian(loaded.model,
                differential_selection,
                time, state, derivative, 1.0)
            hcat(J0, J1 - Jdot0)
        end
        correction = matrix \ equations
        state[algebraic] .-= correction[1:length(algebraic)]
        isempty(differential) ||
            (derivative[differential] .-= correction[(length(algebraic) + 1):end])
    end
    error("initial simultaneous implicit solve failed to converge")
end

const DYNAMIC_VELOCITY_KINDS = Set((:velocity, :angular_velocity,
                                    :relative_velocity))
const DYNAMIC_ACCELERATION_KINDS = Set((:acceleration, :angular_acceleration,
                                        :relative_acceleration))

function set_planar_analysis_stage!(loaded, stage)
    stage in (:static, :dynamic, :modal) || throw(ArgumentError(
        "unknown planar analysis stage '$stage'"))
    for component in values(loaded.equation_components)
        set_planar_equation_stage!(component, stage)
    end
    loaded
end

"""Find static equilibrium, then restore declared velocities for dynamics."""
function statically_initialized_state(loaded, time)
    set_planar_analysis_stage!(loaded, :static)
    selection = static_selection(loaded)
    static_values, iterations, relaxation_cycles = solve_static_equilibrium(
        loaded, selection, time,
        loaded.initial_values[selection.variable_indices])
    state = copy(loaded.initial_values)
    state[selection.variable_indices] .= static_values
    # Static force and reaction values are useful Newton guesses, but model-
    # declared velocities remain the dynamic initial conditions.
    for variable in loaded.layout.catalog.variables
        variable.kind in DYNAMIC_VELOCITY_KINDS || continue
        state[variable.index] = loaded.initial_values[variable.index]
    end
    PlanarModelIO.correct_initial_velocities!(state, loaded.model,
        loaded.initial_variable_weights,
        Set(loaded.imposed_initial_variable_indices), time)
    set_planar_analysis_stage!(loaded, :dynamic)
    return state, iterations, relaxation_cycles
end

"""
Advance one active planar canonical system with the BDF integrator.

The Newton unknown contains every active canonical variable. The differential
and error-control masks identify only the minimal physical state set; they do
not remove algebraic variables from the corrector. Runtime reselection changes
the candidate state-equation rows and masks, then reuses the accepted history.
"""
function run_implicit_model(loaded, times, analysis_mode;
        initial_state = nothing, fixed_model_time = nothing,
        sample_progress = nothing, analysis_stage = :dynamic)
    set_planar_analysis_stage!(loaded, analysis_stage)
    # Establish one complete consistent canonical state before reducing the
    # arrays passed to DASSL to their active rows and columns.
    state, static_initialization_iterations, static_relaxation_cycles =
        if !isnothing(initial_state)
            (copy(initial_state), nothing, nothing)
        elseif loaded.analysis.initialization == :static_equilibrium
            statically_initialized_state(loaded, first(times))
        else
            (copy(loaded.initial_values), nothing, nothing)
        end
    model_time(t) = isnothing(fixed_model_time) ? t : fixed_model_time
    derivative, initial_iterations =
        initialize_implicit_model!(state, loaded, model_time(first(times)))
    isnothing(sample_progress) || sample_progress((;
        kind = :sample, time = Float64(first(times)), state = copy(state)))
    active_variables = loaded.active_variable_indices
    partition = runtime_state_partition(loaded)
    selection = AnalysisSelection(Dynamics(), active_variables,
        partition.equation_indices)

    # These closures translate between DASSL's compact active vectors and the
    # canonical indexing used by every component callback.
    canonical = copy(state)
    canonical_derivative = zeros(length(state))
    residual! = function (equations, t, z, zdot, parameter)
        canonical[active_variables] .= z
        canonical_derivative .= 0
        canonical_derivative[active_variables] .= zdot
        evaluate_analysis_equations!(equations, loaded.model, selection,
            model_time(t), canonical, canonical_derivative)
    end
    jacobian! = function (matrix, t, z, zdot, coefficient, parameter)
        canonical[active_variables] .= z
        canonical_derivative .= 0
        canonical_derivative[active_variables] .= zdot
        evaluate_analysis_sparse_jacobian!(matrix, loaded.model, selection,
            model_time(t), canonical, canonical_derivative, coefficient)
    end
    prototype = evaluate_analysis_sparse_jacobian(loaded.model, selection,
        model_time(first(times)), state, derivative, 1.0)
    differential = falses(length(active_variables))
    positions, velocities = runtime_state_variables(partition)
    active_lookup = Dict(variable => local_index
        for (local_index, variable) in enumerate(active_variables))
    differential[[active_lookup[index] for index in positions]] .= true
    differential[[active_lookup[index] for index in velocities]] .= true
    for component in values(loaded.equation_components)
        for index in component.state_indices
            haskey(active_lookup, index) || continue
            differential[active_lookup[index]] = true
        end
    end

    # Measurements are useful output but should not cause step rejection.
    # Physical positions and velocities provide a broader health monitor than
    # the selected states alone.
    measurement_variables = Set{Int}()
    for name in keys(loaded.measures)
        union!(measurement_variables,
            component_variable_indices(loaded.layout, name))
    end
    physical_error_variables = BitVector(
        variable.index ∉ measurement_variables &&
        variable.kind in (:position, :orientation, :velocity,
                          :angular_velocity, :relative_position,
                          :relative_velocity, :user_state_hold,
                          :user_state_steady)
        for variable in loaded.layout.catalog.variables[active_variables])
    error_control = if any(differential)
        copy(differential)
    else
        copy(physical_error_variables)
    end
    state_selection_changes = NamedTuple[(;
        time = Float64(first(times)), reason = :initial,
        previous = Symbol[], selected = copy(partition.selected_names))]

    # A new state set is an equivalent row partition of the same canonical
    # equations. DASSL keeps its polynomial history but must rebuild the sparse
    # pattern and factorization for the replacement rows.
    function reselect_states!(reason, t, z, zdot, internal_error_control,
            internal_differential, parameter)
        canonical[active_variables] .= z
        canonical_derivative .= 0
        canonical_derivative[active_variables] .= zdot
        choice, diagnostics = automatic_runtime_state_selection(
            loaded, canonical, model_time(t))
        length(choice) == length(partition.selected_names) || return nothing
        Set(choice) == Set(partition.selected_names) && return nothing
        previous = copy(partition.selected_names)
        select_runtime_states!(partition, choice)
        new_positions, new_velocities = runtime_state_variables(partition)
        fill!(internal_differential, false)
        internal_differential[[active_lookup[index]
            for index in new_positions]] .= true
        internal_differential[[active_lookup[index]
            for index in new_velocities]] .= true
        internal_error_control .= internal_differential
        push!(state_selection_changes, (; time = Float64(t), reason,
            previous, selected = copy(partition.selected_names),
            rank = diagnostics.rank))
        new_prototype = evaluate_analysis_sparse_jacobian(loaded.model,
            selection, model_time(t), canonical, canonical_derivative, 1.0)
        new_equation_levels = getproperty.(
            loaded.layout.catalog.equations[partition.equation_indices], :level)
        (; jacobian_prototype = new_prototype,
           equation_levels = new_equation_levels)
    end
    settings = loaded.simulation
    options = HistoricalDDASSL.DASSLOptions{Float64}(
        atol = settings.absolute_tolerance,
        rtol = settings.relative_tolerance,
        initial_step = settings.initial_step,
        maximum_step = settings.maximum_step)
    contacts = PlanarPlaneContactComponent[]
    for components in values(loaded.forces), component in components
        component isa PlanarPlaneContactComponent && push!(contacts, component)
    end
    root_count = sum(contact -> contact.damping_factor > 0 ? 2 : 1, contacts;
        init = 0)
    root_callback = if isempty(contacts)
        nothing
    else
        function (roots, t, z, zdot, parameter)
            canonical[active_variables] .= z
            index = 1
            for contact in contacts
                roots[index] = plane_contact_gap(contact, canonical)
                index += 1
                if contact.damping_factor > 0
                    roots[index] = plane_contact_damping_surface(
                        contact, canonical)
                    index += 1
                end
            end
        end
    end

    # Output frames use the accepted BDF history polynomial. They therefore do
    # not force integration steps at every requested display time.
    next_output = Ref(2)
    accepted_interval! = if isnothing(sample_progress)
        nothing
    else
        function (left_time, right_time, nodes, history_values, order,
                step, statistics)
            tolerance = 16eps(Float64) * max(abs(left_time),
                abs(right_time), 1.0)
            while next_output[] <= length(times) &&
                    times[next_output[]] <= right_time + tolerance
                output_time = Float64(times[next_output[]])
                if output_time > left_time + tolerance
                    values, _ = HistoricalDDASSL.interpolation_state(
                        output_time, nodes, history_values)
                    expanded = copy(state)
                    expanded[active_variables] .= values
                    sample_progress((; kind = :sample,
                        time = output_time, state = expanded))
                end
                next_output[] += 1
            end
        end
    end
    last_accepted_time = Ref(Float64(first(times)))
    last_accepted_values = copy(state[active_variables])
    last_accepted_statistics = Ref{Any}(nothing)
    tracked_accepted_step! = if isnothing(sample_progress)
        nothing
    else
        function (time, values, rates, order, step, statistics)
            last_accepted_values .= values
            last_accepted_time[] = Float64(time)
            last_accepted_statistics[] = statistics
        end
    end
    solution = try
        HistoricalDDASSL.dassl(residual!, state[active_variables],
            derivative[active_variables],
            (first(times), last(times)); jacobian!, jacobian_prototype = prototype,
            options, variable_levels = getproperty.(
                loaded.layout.catalog.variables[active_variables], :level),
            equation_levels = getproperty.(
                loaded.layout.catalog.equations[partition.equation_indices], :level),
            differential_vars = differential, error_control,
            error_monitor = physical_error_variables,
            deficit = loaded.analysis.deficit, root! = root_callback,
            number_of_roots = root_count, root_restart = :soft,
            reconfigure! = reselect_states!,
            accepted_step! = tracked_accepted_step!, accepted_interval!)
    catch exception
        if !isnothing(sample_progress)
            expanded = copy(state)
            expanded[active_variables] .= last_accepted_values
            sample_progress((; kind = :failure,
                time = last_accepted_time[], state = expanded,
                message = exception isa InterruptException ?
                    "planar dynamic integration interrupted" :
                    "planar dynamic integration stopped: " *
                        sprint(showerror, exception),
                statistics = last_accepted_statistics[]))
        end
        rethrow()
    end
    if !SciMLBase.successful_retcode(solution.retcode)
        expanded = copy(state)
        expanded[active_variables] .= last(solution.y)
        isnothing(sample_progress) || sample_progress((;
            kind = :failure, time = Float64(last(solution.t)),
            state = expanded,
            message = "dynamic integration failed: " * solution.message,
            statistics = solution.stats))
        error("dynamic integration failed: $(solution.message)")
    end
    active_states = solution(collect(times))
    states = [begin
        expanded = copy(state)
        expanded[active_variables] .= active_state
        expanded
    end for active_state in active_states]
    (; times = collect(times), states,
       initial_consistency_iterations = initial_iterations, loaded,
       static_initialization_iterations, static_relaxation_cycles,
       analysis_mode, solution, state_selection_changes)
end

"""Select the active position, reaction, and load equations for statics."""
function static_selection(loaded)
    complete = select_analysis(loaded.layout.catalog, StaticEQ())
    active_variables = Set(loaded.active_variable_indices)
    active_equations = Set(loaded.active_equation_indices)
    selection = AnalysisSelection(StaticEQ(),
        [index for index in complete.variable_indices if index in active_variables],
        [index for index in complete.equation_indices if index in active_equations])
    length(selection.variable_indices) == length(selection.equation_indices) ||
        throw(ArgumentError("static system has $(length(selection.variable_indices)) " *
            "unknowns and $(length(selection.equation_indices)) equations"))
    selection
end

"""
Solve one planar static operating point by Newton correction.

No velocity, acceleration, or inertia terms are retained. The result is an
exact equilibrium of the selected static equations to the requested Newton
tolerances.
"""
function static_equilibrium(loaded, selection, time, initial;
        tolerance = loaded.simulation.absolute_tolerance,
        maximum_iterations = 60)
    y = copy(initial)
    canonical = zeros(eltype(y), length(loaded.layout.catalog.variables))
    derivative = zeros(eltype(y), length(canonical))
    equations = zeros(eltype(y), length(selection.equation_indices))
    function populate_canonical!(source)
        canonical .= 0
        for component in Base.values(loaded.equation_components)
            for (index, mode) in zip(component.state_indices,
                    component.state_modes)
                mode == :hold &&
                    (canonical[index] = loaded.initial_values[index])
            end
        end
        canonical[selection.variable_indices] .= source
    end
    function evaluate!(target, source)
        populate_canonical!(source)
        evaluate_analysis_equations!(target, loaded.model, selection,
            time, canonical, derivative)
    end
    for iteration in 0:maximum_iterations
        try
            evaluate!(equations, y)
        catch error
            error isa DomainError || rethrow()
            throw(StaticConvergenceError(
                "static equations are undefined at time $time"))
        end
        norm(equations, Inf) <= tolerance && return y, iteration
        iteration == maximum_iterations && break
        populate_canonical!(y)
        jacobian = evaluate_analysis_sparse_jacobian(loaded.model, selection,
            time, canonical, derivative, 0.0)
        correction = try
            -(jacobian \ equations)
        catch error
            error isa LinearAlgebra.SingularException || rethrow()
            # An ideal mechanism may have a singular first static tangent when
            # its reaction estimates are zero. A pivoted least-squares step
            # establishes those reactions; later Newton tangents then recover
            # the configuration stiffness supplied by the reactions.
            -(qr(Matrix(jacobian), ColumnNorm()) \ equations)
        end
        all(isfinite, correction) || throw(StaticConvergenceError(
            "static Jacobian is singular at time $time"))
        initial_norm = norm(equations, Inf)
        factor = 1.0
        accepted = false
        while factor >= 1 / 1024
            trial = y + factor .* correction
            trial_equations = similar(equations)
            try
                evaluate!(trial_equations, trial)
                if norm(trial_equations, Inf) < initial_norm
                    y = trial
                    accepted = true
                    break
                end
            catch error
                error isa DomainError || rethrow()
            end
            factor /= 2
        end
        accepted || throw(StaticConvergenceError(
            "static Newton line search failed at time $time"))
    end
    throw(StaticConvergenceError(
        "static Newton iteration did not converge at time $time"))
end

"""
Approach planar equilibrium through repeated dynamic pseudo-time intervals.

Velocities and accelerations are depleted between intervals. This adds
controllable damping while preserving the ordinary dynamic equations. The
relaxed state is subsequently polished by the static Newton solve.
"""
function dynamic_relaxation_equilibrium(loaded, selection, time, initial)
    state = copy(loaded.initial_values)
    state[selection.variable_indices] .= initial
    for variable in loaded.layout.catalog.variables
        variable.kind in DYNAMIC_VELOCITY_KINDS ||
            variable.kind in DYNAMIC_ACCELERATION_KINDS || continue
        state[variable.index] = 0.0
    end
    settings = loaded.analysis
    pseudo_times = [0.0, settings.relaxation_duration]
    for cycle in 1:settings.relaxation_cycles
        relaxation = run_implicit_model(loaded, pseudo_times, :dynamic;
            initial_state = state, fixed_model_time = time,
            analysis_stage = :static)
        state = last(relaxation.states)
        for variable in loaded.layout.catalog.variables
            variable.kind in DYNAMIC_VELOCITY_KINDS ||
                variable.kind in DYNAMIC_ACCELERATION_KINDS || continue
            state[variable.index] *= settings.relaxation_reduction_factor
        end
        try
            values, iterations = static_equilibrium(loaded, selection, time,
                state[selection.variable_indices]; maximum_iterations = 12)
            return values, iterations, cycle
        catch error
            error isa StaticConvergenceError || rethrow()
        end
    end
    throw(StaticConvergenceError(
        "dynamic relaxation did not reach static equilibrium at time $time " *
        "after $(settings.relaxation_cycles) cycles"))
end

function solve_static_equilibrium(loaded, selection, time, initial)
    set_planar_analysis_stage!(loaded, :static)
    if loaded.analysis.static_method == :newton
        values, iterations = static_equilibrium(
            loaded, selection, time, initial)
        return values, iterations, 0
    end
    dynamic_relaxation_equilibrium(loaded, selection, time, initial)
end

"""Predict the next quasi-static configuration from up to three solved points."""
function static_predict(times, states, target)
    length(states) < 2 && return copy(last(states))
    t0, t1 = times[end - 1], times[end]
    t1 == t0 && return copy(last(states))
    last(states) .+ ((target - t1) / (t1 - t0)) .* (states[end] .- states[end - 1])
end

"""Advance one quasi-static target, bisecting model time after a failed solve."""
function advance_static!(accepted_times, accepted_states, target, loaded,
        selection, iterations, relaxation_cycles; minimum_step)
    prediction = static_predict(accepted_times, accepted_states, target)
    try
        state, count, cycles = solve_static_equilibrium(
            loaded, selection, target, prediction)
        push!(accepted_times, target)
        push!(accepted_states, state)
        push!(iterations, count)
        push!(relaxation_cycles, cycles)
        return state
    catch error
        error isa StaticConvergenceError || rethrow()
        previous = last(accepted_times)
        abs(target - previous) > minimum_step || rethrow()
        midpoint = previous + (target - previous) / 2
        advance_static!(accepted_times, accepted_states, midpoint, loaded,
            selection, iterations, relaxation_cycles; minimum_step)
        return advance_static!(accepted_times, accepted_states, target, loaded,
            selection, iterations, relaxation_cycles; minimum_step)
    end
end

function run_static_model(loaded, times)
    selection = static_selection(loaded)
    initial = loaded.initial_values[selection.variable_indices]
    first_state, first_iterations, first_relaxation_cycles = solve_static_equilibrium(
        loaded, selection, first(times), initial)
    accepted_times = Float64[first(times)]
    accepted_states = [first_state]
    all_iterations = Int[first_iterations]
    all_relaxation_cycles = Int[first_relaxation_cycles]
    output_states = Vector{Vector{Float64}}()
    function expand_static(values)
        state = zeros(length(loaded.initial_values))
        for component in Base.values(loaded.equation_components)
            for (index, mode) in zip(component.state_indices,
                    component.state_modes)
                mode == :hold &&
                    (state[index] = loaded.initial_values[index])
            end
        end
        state[selection.variable_indices] .= values
        state
    end
    push!(output_states, expand_static(first_state))
    interval = max(abs(last(times) - first(times)), 1.0)
    minimum_step = sqrt(eps(Float64)) * interval
    for target in Iterators.drop(times, 1)
        state = advance_static!(accepted_times, accepted_states, target,
            loaded, selection, all_iterations, all_relaxation_cycles;
            minimum_step)
        push!(output_states, expand_static(state))
    end
    (; times = collect(times), states = output_states,
       initial_consistency_iterations = first_iterations, loaded,
       analysis_mode = :static, static_iterations = all_iterations,
       static_relaxation_cycles = all_relaxation_cycles,
       static_variable_indices = selection.variable_indices,
       static_equation_indices = selection.equation_indices,
       continuation_solutions = length(accepted_times))
end

function modal_operating_point(loaded, time; initial_state = nothing)
    if !isnothing(initial_state)
        state = copy(initial_state)
        derivative, initial_iterations =
            initialize_implicit_model!(state, loaded, time)
        return (; state, derivative, initial_iterations,
            static_initialization_iterations = nothing,
            static_relaxation_cycles = nothing)
    end
    state, static_initialization_iterations, static_relaxation_cycles =
        if loaded.analysis.initialization == :static_equilibrium
            values = statically_initialized_state(loaded, time)
            (values[1], values[2], values[3])
        else
            (copy(loaded.initial_values), nothing, nothing)
        end
    if loaded.analysis.initialization == :static_equilibrium
        for variable in loaded.layout.catalog.variables
            variable.kind in DYNAMIC_VELOCITY_KINDS ||
                variable.kind in DYNAMIC_ACCELERATION_KINDS || continue
            state[variable.index] = 0.0
        end
    end
    derivative, initial_iterations =
        initialize_implicit_model!(state, loaded, time)
    (; state, derivative, initial_iterations,
       static_initialization_iterations, static_relaxation_cycles)
end

function run_modal_model(loaded, time; initial_state = nothing)
    set_planar_analysis_stage!(loaded, :modal)
    operating_point = modal_operating_point(loaded, time; initial_state)
    set_planar_analysis_stage!(loaded, :modal)
    modes = solve_modal_system(loaded, operating_point.state,
        operating_point.derivative, time)
    (; times = [Float64(time)], states = [operating_point.state],
       operating_state = operating_point.state,
       operating_derivative = operating_point.derivative,
       initial_consistency_iterations = operating_point.initial_iterations,
       loaded, analysis_mode = :modal,
       static_initialization_iterations =
           operating_point.static_initialization_iterations,
       static_relaxation_cycles = operating_point.static_relaxation_cycles,
       modes...)
end

"""
    run_planar_model(source; start_time=nothing, end_time=nothing,
                     duration=nothing, samples=nothing)

Load and run a planar TOML model from a path or `IO` stream.

The model's simulation table supplies values for omitted keywords. `duration`
takes precedence over `end_time` and is measured from the effective start time.
Automatic analysis uses dynamic integration when independent states exist and
kinematic continuation otherwise; an explicit static mode performs one or more
equilibrium solutions. Modal mode linearizes the complete implicit equations at
the consistent state at `start_time`. A dynamic or modal model may request one
static-equilibrium configuration solve first. Static equilibrium may use direct
Newton iteration or dynamic relaxation at fixed model time followed by a Newton
correction. The returned named tuple contains canonical sampled states or a
modal operating point, times, the loaded model, initialization diagnostics,
and analysis-specific solver statistics.
"""
function run_planar_model(source; start_time = nothing, end_time = nothing,
        duration = nothing, samples = nothing, result_progress = nothing)
    loaded = source isa LoadedPlanarModel ? source : load_planar_model(source)
    start = isnothing(start_time) ? loaded.simulation.start_time : Float64(start_time)
    finish = if !isnothing(duration)
        Float64(duration) >= 0 || throw(ArgumentError("duration must be nonnegative"))
        start + Float64(duration)
    else
        isnothing(end_time) ? loaded.simulation.end_time : Float64(end_time)
    end
    finish >= start || throw(ArgumentError("end_time must not precede start_time"))
    count = isnothing(samples) ? loaded.simulation.output_samples : Int(samples)
    count >= 1 || throw(ArgumentError("samples must be positive"))
    times = count == 1 ? [start] : range(start, finish; length = count)
    mode = loaded.analysis.mode == :automatic ?
        (loaded.analysis.degrees_of_freedom == 0 ? :kinematic : :dynamic) :
        loaded.analysis.mode
    isnothing(result_progress) || result_progress((; kind = :begin, loaded,
        analysis_mode = mode, times = collect(times)))
    mode == :kinematic && loaded.analysis.degrees_of_freedom > 0 &&
        throw(ArgumentError("kinematic analysis requires a mechanism with no states"))
    mode == :modal && loaded.analysis.degrees_of_freedom == 0 &&
        throw(ArgumentError("modal analysis requires at least one independent state"))
    loaded.analysis.initialization == :static_equilibrium &&
            mode ∉ (:dynamic, :modal) &&
        throw(ArgumentError(
            "static_equilibrium initialization requires dynamic or modal analysis"))
    mode == :modal ? run_modal_model(loaded, start) :
        mode == :static ? run_static_model(loaded, times) :
        run_implicit_model(loaded, times, mode;
            sample_progress = result_progress)
end

end
