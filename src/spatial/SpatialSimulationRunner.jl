"""
    SpatialSimulationRunner

Analysis dispatcher for a loaded spatial model. It provides weighted
consistent initialization, static Newton and dynamic-relaxation solutions,
quasi-static continuation, fully implicit BDF dynamics, runtime state
reselection, and sparse modal linearization.

Finite orientation remains in normalized Euler parameters. Mechanical
orientation columns in the dynamic Newton matrix are transformed to three
body-fixed pseudo angles; only the Euler-parameter kinematic and normalization
equations retain parameter columns.
"""
module SpatialSimulationRunner

using LinearAlgebra
using SciMLBase
using SparseArrays
using ..AutomaticAnalysis
using ..HistoricalDDASSL
using ..ModalAnalysis
using ..SpatialComponentAssembly
using ..SpatialConstraints: connection_hinge, connection_inline,
    initialize_hinge_coordinates!, initialize_inline_coordinates!
using ..SpatialGearPairs: SpatialGearPair,
    initialize_spatial_gear_coordinates!
using ..SpatialSpans: SpatialSpanMeasure,
    initialize_spatial_span_measure!, initialize_spatial_span!
using ..SpatialBelts: SpatialBeltSpanComponent,
    initialize_spatial_belt_span!
using ..SpatialDirectedDistances: SpatialDirectedDistanceMeasure,
    initialize_spatial_directed_distance_measure!
using ..SpatialAppliedForces: SpatialSpanningForceComponent,
    SpatialAppliedForceComponent, SpatialAppliedTorqueComponent,
    initialize_spatial_spanning_force!, initialize_spatial_applied_force!,
    initialize_spatial_applied_torque!, set_spatial_applied_force_stage!,
    set_spatial_applied_torque_stage!, set_spatial_spanning_force_stage!
using ..SpatialBushings: SpatialBushingComponent, initialize_spatial_bushing!,
    set_spatial_bushing_stage!
using ..SpatialPlaneContacts: SpatialPlaneContactComponent,
    initialize_spatial_plane_contact!, set_spatial_plane_contact_stage!
using ..SpatialFrictionForces: SpatialSurfaceFriction,
    initialize_spatial_surface_friction!, set_spatial_surface_friction_stage!,
    surface_friction_rates, SpatialRevoluteFriction,
    initialize_spatial_revolute_friction!,
    set_spatial_revolute_friction_stage!, revolute_friction_rate,
    SpatialTranslationalFriction,
    initialize_spatial_translational_friction!,
    set_spatial_translational_friction_stage!, translational_friction_rate,
    SpatialInplaneFriction, initialize_spatial_inplane_friction!,
    set_spatial_inplane_friction_stage!, inplane_friction_rates
using ..SpatialTires: SpatialTireComponent, initialize_spatial_tire!,
    tire_deformation_rates
using ..SpatialEquationComponents: spatial_equation_state_rates!,
    set_spatial_equation_stage!
using ..SpatialMotionGenerators: SpatialRotationalMotionGenerator,
    SpatialTranslationalMotionGenerator, SpatialSpanningMotionGenerator
using ..SpatialModelIO

export run_spatial_model

struct SpatialStaticConvergenceError <: Exception
    message::String
end

Base.showerror(io::IO, error::SpatialStaticConvergenceError) =
    print(io, error.message)

const SPATIAL_STATIC_RATE_KINDS = Set((
    :acceleration, :angular_acceleration, :relative_acceleration,
    :velocity, :angular_velocity, :relative_velocity, :applied_rate))
const STATIC_RECIPROCAL_CONDITION_LIMIT = 1.0e-10
const STATE_RESELECTION_IMPROVEMENT = 2.0
# SuiteSparse numbers UMFPACK_INFO entries from zero. Julia's Vector uses one.
const UMFPACK_RCOND_INFO_INDEX = 68

emit_static_progress(progress, event) =
    isnothing(progress) ? nothing : progress(event)

"""One spatial velocity and coordinate pair eligible for state selection."""
struct RuntimeSpatialStateCandidate
    name::Symbol
    velocity::Int
    coordinate::Int
    equations::Vector{Int}
end

"""
Active spatial state-equation rows and the dormant equivalent alternatives.

`base_equations` is invariant during integration. Reselection replaces only
the candidate rows appended to it, without changing system size.
"""
mutable struct RuntimeSpatialStatePartition
    candidates::Dict{Symbol,RuntimeSpatialStateCandidate}
    base_equations::Vector{Int}
    equation_indices::Vector{Int}
    selected_names::Vector{Symbol}
end

function runtime_spatial_state_partition(loaded)
    state_selection = loaded.state_selection
    candidates = Dict{Symbol,RuntimeSpatialStateCandidate}()
    for (name, velocity) in zip(state_selection.candidate_names,
            state_selection.velocity_indices)
        candidates[name] = RuntimeSpatialStateCandidate(name, velocity,
            state_selection.coordinate_indices[velocity],
            copy(state_selection.state_equations[velocity]))
    end
    candidate_equations = Set(Iterators.flatten(
        candidate.equations for candidate in values(candidates)))
    base_equations = [index for index in loaded.active_equation_indices
        if index ∉ candidate_equations]
    selected_names = copy(state_selection.selected_velocities)
    selected_equations = collect(Iterators.flatten(
        candidates[name].equations for name in selected_names))
    RuntimeSpatialStatePartition(candidates, base_equations,
        [base_equations; selected_equations], selected_names)
end

function select_runtime_spatial_states!(partition, names)
    length(names) == length(partition.selected_names) ||
        throw(DimensionMismatch(
            "runtime spatial state selection changed the number of states"))
    all(name -> haskey(partition.candidates, name), names) ||
        throw(ArgumentError(
            "runtime spatial state selection named an unknown candidate"))
    partition.selected_names = collect(names)
    empty!(partition.equation_indices)
    append!(partition.equation_indices, partition.base_equations)
    for name in partition.selected_names
        append!(partition.equation_indices,
            partition.candidates[name].equations)
    end
    partition
end

function automatic_runtime_spatial_state_selection(loaded, state, time)
    selection = loaded.state_selection
    indices, diagnostics =
        SpatialModelIO.automatic_spatial_velocity_selection(
            loaded.model, selection.velocity_indices,
            selection.velocity_rows, selection.column_scales, state, time)
    names = [Symbol(loaded.layout.catalog.variables[index].component, ".",
        loaded.layout.catalog.variables[index].name) for index in indices]
    names, diagnostics
end

"""
Construct canonical rates from explicit kinematic levels and internal states.

Euler-parameter rates are calculated from body angular velocity. Hinge,
inline, span-generator, and transient-tire states obtain their rates from the
corresponding explicit velocity or constitutive relation.
"""
function initial_spatial_derivative(state, loaded,
        time = loaded.simulation.start_time)
    derivative = zeros(eltype(state), length(state))
    for body in values(loaded.bodies)
        derivative[body.position_variables] .= state[body.velocity_variables]
        derivative[body.pseudo_angle_variables] .=
            state[body.angular_velocity_variables]
        derivative[body.velocity_variables] .=
            state[body.acceleration_variables]
        derivative[body.angular_velocity_variables] .=
            state[body.angular_acceleration_variables]
        parameters = @view state[body.euler_parameter_variables]
        omega = @view state[body.angular_velocity_variables]
        derivative[body.euler_parameter_variables] .=
            0.5 .* quaternion_rate_matrix(parameters) * omega
    end
    for connection in values(loaded.connections)
        hinge = connection_hinge(connection)
        if !isnothing(hinge) && !isempty(hinge.rotation_variables)
            alpha, omega, theta = hinge.rotation_variables
            derivative[theta] = state[omega]
            derivative[omega] = state[alpha]
        end
        inline = connection_inline(connection)
        if !isnothing(inline) && !isempty(inline.translation_variables)
            acceleration, velocity, distance = inline.translation_variables
            derivative[distance] = state[velocity]
            derivative[velocity] = state[acceleration]
        end
    end
    for driver in values(loaded.drivers)
        driver isa Union{SpatialTranslationalMotionGenerator,
            SpatialSpanningMotionGenerator} || continue
        if driver isa SpatialSpanningMotionGenerator
            derivative[driver.span.distance_variable] =
                state[driver.span.velocity_variable]
            derivative[driver.span.velocity_variable] =
                state[driver.span.acceleration_variable]
        else
            derivative[driver.distance_variable] =
                state[driver.velocity_variable]
            derivative[driver.velocity_variable] =
                state[driver.acceleration_variable]
        end
    end
    for force in values(loaded.forces)
        force isa SpatialTireComponent && force.transient || continue
        longitudinal_rate, lateral_rate =
            tire_deformation_rates(force, state)
        derivative[force.longitudinal_deformation_variable] =
            longitudinal_rate
        derivative[force.lateral_deformation_variable] = lateral_rate
    end
    for force in values(loaded.forces)
        force isa SpatialSurfaceFriction || continue
        derivative[force.shear_variables] .= surface_friction_rates(force, state)
    end
    for force in values(loaded.forces)
        force isa SpatialRevoluteFriction || continue
        derivative[force.shear_variable] = revolute_friction_rate(force, state)
    end
    for force in values(loaded.forces)
        force isa SpatialTranslationalFriction || continue
        derivative[force.shear_variable] =
            translational_friction_rate(force, state)
    end
    for force in values(loaded.forces)
        force isa SpatialInplaneFriction || continue
        derivative[force.shear_variables] .= inplane_friction_rates(force, state)
    end
    for component in values(loaded.equation_components)
        spatial_equation_state_rates!(derivative, component,
            time, state)
    end
    derivative
end

"""
Return DASSL differential, controlled-error, and monitored-error masks.

All body positions, Euler parameters, and physical velocities are monitored.
Only the selected independent coordinates and true internal states are marked
differential; algebraic reactions and force definitions remain in the Newton
solve but outside integration error control.
"""
function spatial_state_masks(loaded, partition)
    differential = falses(length(loaded.initial_values))
    error_control = falses(length(loaded.initial_values))
    physical_monitor = falses(length(loaded.initial_values))
    for body in values(loaded.bodies)
        for indices in (body.position_variables,
                body.euler_parameter_variables, body.velocity_variables,
                body.angular_velocity_variables)
            error_control[indices] .= true
            physical_monitor[indices] .= true
        end
        differential[body.euler_parameter_variables] .= true
        for (variable, equation) in zip(body.pseudo_angle_variables,
                body.pseudo_angle_state_equations)
            equation in partition.base_equations &&
                (differential[variable] = true)
        end
    end
    for name in partition.selected_names
        candidate = partition.candidates[name]
        differential[candidate.velocity] = true
        error_control[candidate.velocity] = true
        physical_monitor[candidate.velocity] = true
        candidate.coordinate in loaded.active_variable_indices &&
            (differential[candidate.coordinate] = true)
        candidate.coordinate in loaded.active_variable_indices &&
            (error_control[candidate.coordinate] = true)
        candidate.coordinate in loaded.active_variable_indices &&
            (physical_monitor[candidate.coordinate] = true)
    end
    for force in values(loaded.forces)
        force isa SpatialTireComponent && force.transient || continue
        for variable in (force.longitudinal_deformation_variable,
                force.lateral_deformation_variable)
            differential[variable] = true
            error_control[variable] = true
            physical_monitor[variable] = true
        end
    end
    for force in values(loaded.forces)
        force isa SpatialSurfaceFriction || continue
        differential[force.shear_variables] .= true
        error_control[force.shear_variables] .= true
        physical_monitor[force.shear_variables] .= true
    end
    for force in values(loaded.forces)
        force isa SpatialRevoluteFriction || continue
        differential[force.shear_variable] = true
        error_control[force.shear_variable] = true
        physical_monitor[force.shear_variable] = true
    end
    for force in values(loaded.forces)
        force isa SpatialTranslationalFriction || continue
        differential[force.shear_variable] = true
        error_control[force.shear_variable] = true
        physical_monitor[force.shear_variable] = true
    end
    for force in values(loaded.forces)
        force isa SpatialInplaneFriction || continue
        differential[force.shear_variables] .= true
        error_control[force.shear_variables] .= true
        physical_monitor[force.shear_variables] .= true
    end
    for component in values(loaded.equation_components)
        for variable in component.state_indices
            differential[variable] = true
            error_control[variable] = true
            physical_monitor[variable] = true
        end
    end
    active = loaded.active_variable_indices
    local_differential = differential[active]
    local_differential, error_control[active], physical_monitor[active]
end

"""
Mapping from the canonical Euler-parameter Jacobian assembled by the spatial
components to the tangent-coordinate Jacobian used by spatial dynamics.

Component equations evaluate finite orientation from Euler parameters. Their
mechanical orientation partials are re-expressed in the three body-fixed
pseudo-angle columns. Only the Euler-parameter kinematic and normalization
rows retain Euler-parameter columns.
"""
struct SpatialDynamicJacobianWorkspace{T,Ti}
    raw_matrix::SparseMatrixCSC{T,Ti}
    prototype::SparseMatrixCSC{T,Ti}
    direct_locations::Vector{Int}
    mapped_bodies::Vector{Int}
    parameter_components::Vector{Int}
    pseudo_locations::Vector{NTuple{3,Int}}
    bodies::Vector
end

function sparse_location_map(matrix::SparseMatrixCSC)
    locations = Dict{Tuple{Int,Int},Int}()
    for column in axes(matrix, 2)
        for location in matrix.colptr[column]:(matrix.colptr[column + 1] - 1)
            locations[(matrix.rowval[location], column)] = location
        end
    end
    locations
end

function spatial_dynamic_jacobian_workspace(loaded, selection, time, state,
        derivative, coefficient)
    raw = evaluate_analysis_sparse_jacobian(loaded.model, selection, time,
        state, derivative, coefficient)
    variable_columns = Dict(canonical => local_index
        for (local_index, canonical) in
        enumerate(selection.variable_indices))
    orientation_rows = Set(Iterators.flatten(
        body.orientation_equations for body in values(loaded.bodies)))
    bodies = [loaded.bodies[name] for name in sort!(collect(keys(loaded.bodies)))]
    parameter_lookup = Dict{Int,Tuple{Int,Int,NTuple{3,Int}}}()
    for (body_index, body) in enumerate(bodies)
        pseudo_columns = Tuple(variable_columns[index]
            for index in body.pseudo_angle_variables)
        for (component, parameter) in enumerate(body.euler_parameter_variables)
            parameter_lookup[variable_columns[parameter]] =
                (body_index, component, pseudo_columns)
        end
    end

    pattern_rows = Int[]
    pattern_columns = Int[]
    for column in axes(raw, 2)
        parameter = get(parameter_lookup, column, nothing)
        for location in raw.colptr[column]:(raw.colptr[column + 1] - 1)
            row = raw.rowval[location]
            canonical_row = selection.equation_indices[row]
            if !isnothing(parameter) && canonical_row ∉ orientation_rows
                for pseudo_column in parameter[3]
                    push!(pattern_rows, row)
                    push!(pattern_columns, pseudo_column)
                end
            else
                push!(pattern_rows, row)
                push!(pattern_columns, column)
            end
        end
    end
    prototype = sparse(pattern_rows, pattern_columns,
        ones(eltype(raw), length(pattern_rows)), size(raw)...)
    fill!(prototype.nzval, zero(eltype(raw)))
    output_locations = sparse_location_map(prototype)
    direct_locations = zeros(Int, nnz(raw))
    mapped_bodies = zeros(Int, nnz(raw))
    parameter_components = zeros(Int, nnz(raw))
    pseudo_locations = fill((0, 0, 0), nnz(raw))
    for column in axes(raw, 2)
        parameter = get(parameter_lookup, column, nothing)
        for location in raw.colptr[column]:(raw.colptr[column + 1] - 1)
            row = raw.rowval[location]
            canonical_row = selection.equation_indices[row]
            if !isnothing(parameter) && canonical_row ∉ orientation_rows
                body_index, component, pseudo_columns = parameter
                mapped_bodies[location] = body_index
                parameter_components[location] = component
                pseudo_locations[location] = Tuple(
                    output_locations[(row, pseudo_column)]
                    for pseudo_column in pseudo_columns)
            else
                direct_locations[location] = output_locations[(row, column)]
            end
        end
    end
    SpatialDynamicJacobianWorkspace(raw, prototype, direct_locations,
        mapped_bodies, parameter_components, pseudo_locations, bodies)
end

function evaluate_spatial_dynamic_jacobian!(matrix, workspace, loaded,
        selection, time, state, derivative, coefficient)
    evaluate_analysis_sparse_jacobian!(workspace.raw_matrix, loaded.model,
        selection, time, state, derivative, coefficient)
    fill!(matrix.nzval, zero(eltype(matrix)))
    parameter_maps = [0.5 .* quaternion_rate_matrix(
        @view state[body.euler_parameter_variables]) for body in workspace.bodies]
    for location in eachindex(workspace.raw_matrix.nzval)
        value = workspace.raw_matrix.nzval[location]
        body_index = workspace.mapped_bodies[location]
        if iszero(body_index)
            matrix.nzval[workspace.direct_locations[location]] += value
        else
            component = workspace.parameter_components[location]
            targets = workspace.pseudo_locations[location]
            transform = parameter_maps[body_index]
            for axis in 1:3
                matrix.nzval[targets[axis]] += value * transform[component, axis]
            end
        end
    end
    matrix
end

"""
Normalize each predicted Euler parameter and choose the corresponding
nonphysical pseudo-angle value so that
`psi_dot = 2Q(p)' * p_dot` for the current BDF history.
"""
function project_spatial_orientation_predictor!(predicted,
        history_derivative, coefficient, predictor_indices)
    for indices in predictor_indices
        p0_index, p1_index, p2_index, p3_index = indices.parameters
        p0 = predicted[p0_index]
        p1 = predicted[p1_index]
        p2 = predicted[p2_index]
        p3 = predicted[p3_index]
        magnitude = sqrt(p0^2 + p1^2 + p2^2 + p3^2)
        p0 /= magnitude
        p1 /= magnitude
        p2 /= magnitude
        p3 /= magnitude
        predicted[p0_index] = p0
        predicted[p1_index] = p1
        predicted[p2_index] = p2
        predicted[p3_index] = p3

        d0 = coefficient * p0 + history_derivative[p0_index]
        d1 = coefficient * p1 + history_derivative[p1_index]
        d2 = coefficient * p2 + history_derivative[p2_index]
        d3 = coefficient * p3 + history_derivative[p3_index]
        psi1_index, psi2_index, psi3_index = indices.pseudo_angles
        psi1_rate = 2(-p1 * d0 + p0 * d1 + p3 * d2 - p2 * d3)
        psi2_rate = 2(-p2 * d0 - p3 * d1 + p0 * d2 + p1 * d3)
        psi3_rate = 2(-p3 * d0 + p2 * d1 - p1 * d2 + p0 * d3)
        predicted[psi1_index] =
            (psi1_rate - history_derivative[psi1_index]) / coefficient
        predicted[psi2_index] =
            (psi2_rate - history_derivative[psi2_index]) / coefficient
        predicted[psi3_index] =
            (psi3_rate - history_derivative[psi3_index]) / coefficient
    end
    nothing
end

"""
Select a square static system and its three-component body rotation updates.

Euler parameters evaluate the operating orientation but are not four
independent static unknowns. Each body instead contributes three local
rotation corrections, giving a nonsingular tangent-coordinate Newton system.
"""
function spatial_static_system(loaded)
    complete = select_analysis(loaded.layout.catalog, StaticEQ())
    active_variables = Set(loaded.active_variable_indices)
    active_equations = Set(loaded.active_equation_indices)
    algebraic_variables = [index for index in complete.variable_indices
        if index in active_variables &&
           loaded.layout.catalog.variables[index].kind != :orientation]
    equation_indices = [index for index in complete.equation_indices
        if index in active_equations]
    body_names = sort!(collect(keys(loaded.bodies)))
    parameter_variables = reduce(vcat,
        [collect(loaded.bodies[name].euler_parameter_variables)
         for name in body_names]; init = Int[])
    jacobian_variables = [algebraic_variables; parameter_variables]
    unknown_count = length(algebraic_variables) + 3length(body_names)
    unknown_count == length(equation_indices) || throw(ArgumentError(
        "spatial static system has $unknown_count unknowns and " *
        "$(length(equation_indices)) equations"))
    system = (; algebraic_variables, equation_indices, body_names,
       parameter_variables, jacobian_variables,
       body_lengths = SpatialModelIO.spatial_body_lengths(
           loaded.bodies, loaded.markers),
       result_variables = sort!(unique!([
           algebraic_variables; parameter_variables])))
    merge(system, (;
        mass_regularization = spatial_static_mass_regularization(
            loaded, system)))
end

"""
Mass and body-frame inertia mapped onto static correction coordinates.

This matrix regularizes a numerically singular static Jacobian. It is not an
inertia force in the equilibrium equations: it is added only to the Newton
matrix and only after the unregularized factorization is judged unhealthy.
"""
function spatial_static_mass_regularization(loaded, system;
        time_scale = 1.0)
    time_scale > 0 || throw(ArgumentError(
        "static regularization time scale must be positive"))
    equation_rows = Dict(equation => row for (row, equation) in
        enumerate(system.equation_indices))
    variable_columns = Dict(variable => column for (column, variable) in
        enumerate(system.algebraic_variables))
    rows = Int[]
    columns = Int[]
    coefficients = Float64[]
    scale = inv(Float64(time_scale)^2)
    angle_offset = length(system.algebraic_variables)
    for (body_number, name) in enumerate(system.body_names)
        body = loaded.bodies[name]
        for axis in 1:3
            row = get(equation_rows, body.balance_equations[axis], 0)
            column = get(variable_columns, body.position_variables[axis], 0)
            iszero(row) || iszero(column) || begin
                push!(rows, row)
                push!(columns, column)
                push!(coefficients, scale * body.mass)
            end
        end
        for row_axis in 1:3, column_axis in 1:3
            row = get(equation_rows,
                body.balance_equations[row_axis + 3], 0)
            column = angle_offset + 3(body_number - 1) + column_axis
            coefficient = scale * body.inertia[row_axis, column_axis]
            iszero(row) || iszero(coefficient) || begin
                push!(rows, row)
                push!(columns, column)
                push!(coefficients, coefficient)
            end
        end
    end
    count = length(system.equation_indices)
    sparse(rows, columns, coefficients, count, count)
end

function static_sparse_factorization(jacobian)
    factorization = lu(jacobian; check = false)
    reciprocal_condition = length(factorization.info) >=
        UMFPACK_RCOND_INFO_INDEX ?
        factorization.info[UMFPACK_RCOND_INFO_INDEX] : NaN
    healthy = factorization.status == 0 &&
        isfinite(reciprocal_condition) &&
        reciprocal_condition >= STATIC_RECIPROCAL_CONDITION_LIMIT
    factorization, reciprocal_condition, healthy
end

function normalize_spatial_orientations!(state, loaded)
    for body in values(loaded.bodies)
        indices = body.euler_parameter_variables
        parameters = @view state[indices]
        magnitude = norm(parameters)
        magnitude > 0 || throw(SpatialStaticConvergenceError(
            "spatial static orientation became undefined"))
        parameters ./= magnitude
    end
    state
end

function zero_spatial_static_rates!(state, loaded)
    for variable in loaded.layout.catalog.variables
        variable.kind in SPATIAL_STATIC_RATE_KINDS || continue
        state[variable.index] = 0.0
    end
    state
end

function spatial_static_equation_diagnostics(loaded, system, equations)
    equation_rows = Dict(index => row for (row, index) in
        enumerate(system.equation_indices))
    force_imbalance = 0.0
    torque_imbalance = 0.0
    equivalent_acceleration = 0.0
    for name in system.body_names
        body = loaded.bodies[name]
        force_rows = [get(equation_rows, index, 0)
            for index in body.balance_equations[1:3]]
        torque_rows = [get(equation_rows, index, 0)
            for index in body.balance_equations[4:6]]
        filter!(!iszero, force_rows)
        filter!(!iszero, torque_rows)
        if !isempty(force_rows)
            body_force = equations[force_rows]
            force_imbalance = max(force_imbalance, norm(body_force))
            body.mass > 0 && (equivalent_acceleration = max(
                equivalent_acceleration, norm(body_force) / body.mass))
        end
        if !isempty(torque_rows)
            body_torque = equations[torque_rows]
            torque_imbalance = max(torque_imbalance, norm(body_torque))
            if body.mass > 0
                angular_acceleration = body.inertia \ body_torque
                equivalent_acceleration = max(equivalent_acceleration,
                    system.body_lengths[name] * norm(angular_acceleration))
            end
        end
    end
    constraint_rows = [row for (row, index) in
        enumerate(system.equation_indices) if begin
            equation = loaded.layout.catalog.equations[index]
            equation.kind == :constraint && equation.level == 0
        end]
    constraint_error = isempty(constraint_rows) ? 0.0 :
        norm(equations[constraint_rows], Inf)
    (; force_imbalance, torque_imbalance, constraint_error,
       equivalent_acceleration)
end

function spatial_static_diagnostics(loaded, system, time, state)
    static_state = copy(state)
    zero_spatial_static_rates!(static_state, loaded)
    normalize_spatial_orientations!(static_state, loaded)
    initialize_spatial_measurements_and_forces!(static_state, loaded, time)
    selection = AnalysisSelection(StaticEQ(), system.jacobian_variables,
        system.equation_indices)
    derivative = zeros(eltype(static_state), length(static_state))
    equations = zeros(eltype(static_state), length(system.equation_indices))
    evaluate_analysis_equations!(equations, loaded.model, selection, time,
        static_state, derivative)
    merge((; residual = norm(equations, Inf)),
        spatial_static_equation_diagnostics(loaded, system, equations))
end

"""Publish the entered and constraint-consistent configurations before static motion."""
function emit_spatial_initial_configurations(progress, loaded, system, time)
    isnothing(progress) && return nothing
    consistency_corrections =
        loaded.initial_conditions.position_corrections +
        loaded.initial_conditions.velocity_corrections +
        loaded.initial_conditions.acceleration_corrections
    configurations = (
        (:entered, loaded.entered_initial_values, 0,
            "body positions and orientations as supplied by the model"),
        (:consistent, loaded.initial_values, consistency_corrections,
            "after initial-condition consistency correction"))
    for (status, source, iteration, message) in configurations
        state = copy(source)
        zero_spatial_static_rates!(state, loaded)
        normalize_spatial_orientations!(state, loaded)
        initialize_spatial_measurements_and_forces!(state, loaded, time)
        diagnostics = spatial_static_diagnostics(
            loaded, system, time, state)
        emit_static_progress(progress, (;
            phase = :initial_conditions, status, loaded, state,
            model_time = Float64(time), pseudo_time = nothing,
            iteration, relaxation_cycle = 0, diagnostics...,
            reciprocal_condition = NaN, mass_regularized = false,
            message))
    end
    nothing
end

function spatial_static_equivalent_speed(state, loaded, system)
    maximum((max(norm(state[loaded.bodies[name].velocity_variables]),
        system.body_lengths[name] * norm(
            state[loaded.bodies[name].angular_velocity_variables]))
        for name in system.body_names); init = 0.0)
end

function spatial_static_configuration_correction(loaded, system, time, state)
    static_state = copy(state)
    zero_spatial_static_rates!(static_state, loaded)
    normalize_spatial_orientations!(static_state, loaded)
    initialize_spatial_measurements_and_forces!(static_state, loaded, time)
    selection = AnalysisSelection(StaticEQ(), system.jacobian_variables,
        system.equation_indices)
    derivative = zeros(eltype(static_state), length(static_state))
    equations = zeros(eltype(static_state), length(system.equation_indices))
    evaluate_analysis_equations!(equations, loaded.model, selection, time,
        static_state, derivative)
    canonical_jacobian = evaluate_analysis_sparse_jacobian(
        loaded.model, selection, time, static_state, derivative, 0.0)
    jacobian = canonical_jacobian *
        spatial_static_correction_map(static_state, loaded, system)
    factorization, _, healthy = static_sparse_factorization(jacobian)
    solve_jacobian = healthy ? jacobian : jacobian + system.mass_regularization
    solve_factorization = healthy ? factorization :
        first(static_sparse_factorization(solve_jacobian))
    correction = solve_factorization.status == 0 ?
        -(solve_factorization \ equations) :
        -(qr(Matrix(solve_jacobian), ColumnNorm()) \ equations)
    all(isfinite, correction) || return Inf
    algebraic_columns = Dict(variable => column for (column, variable) in
        enumerate(system.algebraic_variables))
    configuration_correction = 0.0
    for (body_number, name) in enumerate(system.body_names)
        body = loaded.bodies[name]
        position_columns = [algebraic_columns[variable]
            for variable in body.position_variables]
        translation = norm(correction[position_columns]) /
            system.body_lengths[name]
        rotation_start = length(system.algebraic_variables) +
            3(body_number - 1) + 1
        rotation = norm(correction[rotation_start:(rotation_start + 2)])
        configuration_correction = max(
            configuration_correction, translation, rotation)
    end
    configuration_correction
end

function spatial_static_relaxation_readiness(loaded, system, time, state)
    diagnostics = spatial_static_diagnostics(loaded, system, time, state)
    equivalent_speed = spatial_static_equivalent_speed(state, loaded, system)
    configuration_correction =
        spatial_static_configuration_correction(loaded, system, time, state)
    constraint_limit = max(1.0e-8,
        sqrt(loaded.simulation.absolute_tolerance))
    ready = diagnostics.equivalent_acceleration <=
            loaded.analysis.handoff_acceleration &&
        equivalent_speed <= loaded.analysis.handoff_speed &&
        configuration_correction <= loaded.analysis.handoff_correction &&
        diagnostics.constraint_error <= constraint_limit
    merge(diagnostics, (; equivalent_speed, configuration_correction,
        constraint_limit, ready))
end

function spatial_static_correction_map(state, loaded, system)
    row_for_variable = Dict(index => row for (row, index) in
        enumerate(system.jacobian_variables))
    column_count = length(system.algebraic_variables) +
        3length(system.body_names)
    rows = Int[]
    columns = Int[]
    values = eltype(state)[]
    for (column, variable) in enumerate(system.algebraic_variables)
        push!(rows, row_for_variable[variable])
        push!(columns, column)
        push!(values, one(eltype(state)))
    end
    column = length(system.algebraic_variables)
    for name in system.body_names
        body = loaded.bodies[name]
        parameters = @view state[body.euler_parameter_variables]
        parameter_rows = [row_for_variable[index]
            for index in body.euler_parameter_variables]
        rate_matrix = 0.5 .* quaternion_rate_matrix(parameters)
        for axis in 1:3
            column += 1
            for local_row in 1:4
                value = rate_matrix[local_row, axis]
                iszero(value) && continue
                push!(rows, parameter_rows[local_row])
                push!(columns, column)
                push!(values, value)
            end
        end
    end
    sparse(rows, columns, values, length(system.jacobian_variables),
        column_count)
end

function apply_spatial_static_correction(state, correction, factor,
        loaded, system)
    trial = copy(state)
    algebraic_count = length(system.algebraic_variables)
    trial[system.algebraic_variables] .+=
        factor .* correction[1:algebraic_count]
    offset = algebraic_count
    for name in system.body_names
        body = loaded.bodies[name]
        rotation = factor .* correction[(offset + 1):(offset + 3)]
        offset += 3
        angle = norm(rotation)
        iszero(angle) && continue
        parameters = @view trial[body.euler_parameter_variables]
        orientation = rotation_matrix(parameters)
        increment = axis_angle_rotation(angle, rotation)
        parameters .= matrix_to_euler_parameters(orientation * increment)
    end
    normalize_spatial_orientations!(trial, loaded)
end

"""
Solve one spatial static operating point with a damped Newton iteration.

The line search accepts only corrections that reduce the complete static
equation norm. If the unregularized sparse Jacobian is singular or poorly
conditioned, a mass-and-inertia matrix remains active for subsequent Newton
corrections in that solve. Progress events retain both equation imbalance and
whether regularization was used.
"""
function spatial_static_equilibrium(loaded, system, time, initial;
        tolerance = isnothing(loaded.analysis.static_tolerance) ?
            loaded.simulation.absolute_tolerance :
            loaded.analysis.static_tolerance,
        maximum_iterations = 60, progress = nothing,
        phase = :newton, relaxation_cycle = 0)
    state = copy(initial)
    zero_spatial_static_rates!(state, loaded)
    normalize_spatial_orientations!(state, loaded)
    selection = AnalysisSelection(StaticEQ(), system.jacobian_variables,
        system.equation_indices)
    derivative = zeros(eltype(state), length(state))
    equations = zeros(eltype(state), length(system.equation_indices))
    mass_regularization_active = false
    reciprocal_condition = NaN
    equation_summary = function (values; count = 6)
        indices = partialsortperm(abs.(values),
            1:min(count, length(values)); rev = true)
        join([begin
            declaration = loaded.layout.catalog.equations[
                system.equation_indices[index]]
            string(declaration.component, ".", declaration.name, "=",
                values[index])
        end for index in indices], ", ")
    end
    for iteration in 0:maximum_iterations
        try
            evaluate_analysis_equations!(equations, loaded.model, selection,
                time, state, derivative)
        catch error
            error isa DomainError || rethrow()
            emit_static_progress(progress, (;
                phase, status = :failed, loaded, state = copy(state),
                model_time = Float64(time), pseudo_time = nothing,
                iteration, relaxation_cycle,
                residual = Inf, reciprocal_condition,
                force_imbalance = Inf, torque_imbalance = Inf,
                constraint_error = Inf, equivalent_acceleration = Inf,
                mass_regularized = mass_regularization_active,
                message = "static equations are undefined"))
            throw(SpatialStaticConvergenceError(
                "spatial static equations are undefined at time $time"))
        end
        residual = norm(equations, Inf)
        diagnostics = spatial_static_equation_diagnostics(
            loaded, system, equations)
        converged = residual <= tolerance
        emit_static_progress(progress, (;
            phase, status = converged ? :converged :
                (iteration == 0 ? :initial : :accepted),
            loaded, state = copy(state), model_time = Float64(time),
            pseudo_time = nothing, iteration, relaxation_cycle,
            residual, diagnostics..., reciprocal_condition,
            mass_regularized = mass_regularization_active,
            message = ""))
        converged && return state, iteration, mass_regularization_active
        iteration == maximum_iterations && break
        canonical_jacobian = evaluate_analysis_sparse_jacobian(
            loaded.model, selection, time, state, derivative, 0.0)
        jacobian = canonical_jacobian *
            spatial_static_correction_map(state, loaded, system)
        factorization, reciprocal_condition, healthy =
            static_sparse_factorization(jacobian)
        mass_regularization_active |= !healthy
        solve_jacobian = mass_regularization_active ?
            jacobian + system.mass_regularization : jacobian
        solve_factorization = mass_regularization_active ?
            first(static_sparse_factorization(solve_jacobian)) : factorization
        correction = if solve_factorization.status == 0
            -(solve_factorization \ equations)
        else
            -(qr(Matrix(solve_jacobian), ColumnNorm()) \ equations)
        end
        all(isfinite, correction) || throw(SpatialStaticConvergenceError(
            "spatial static Jacobian is singular at time $time"))
        initial_norm = norm(equations, Inf)
        factor = 1.0
        accepted = false
        while factor >= 1 / 1024
            trial = apply_spatial_static_correction(state, correction,
                factor, loaded, system)
            trial_equations = similar(equations)
            try
                evaluate_analysis_equations!(trial_equations, loaded.model,
                    selection, time, trial, derivative)
                if norm(trial_equations, Inf) < initial_norm
                    state = trial
                    accepted = true
                    break
                end
            catch error
                error isa DomainError || rethrow()
            end
            factor /= 2
        end
        if !accepted
            message = "Newton line search failed; largest equations: " *
                equation_summary(equations)
            emit_static_progress(progress, (;
                phase, status = :failed, loaded, state = copy(state),
                model_time = Float64(time), pseudo_time = nothing,
                iteration, relaxation_cycle, residual,
                diagnostics..., reciprocal_condition,
                mass_regularized = mass_regularization_active, message))
            throw(SpatialStaticConvergenceError(
                "spatial static Newton line search failed at time $time " *
                "after iteration $iteration; largest equations: " *
                equation_summary(equations)))
        end
    end
    emit_static_progress(progress, (;
        phase, status = :failed, loaded, state = copy(state),
        model_time = Float64(time), pseudo_time = nothing,
        iteration = maximum_iterations, relaxation_cycle,
        residual = norm(equations, Inf),
        spatial_static_equation_diagnostics(
            loaded, system, equations)...,
        reciprocal_condition,
        mass_regularized = mass_regularization_active,
        message = "maximum static Newton iterations reached"))
    throw(SpatialStaticConvergenceError(
        "spatial static Newton iteration did not converge at time $time; " *
        "largest equations: " * equation_summary(equations)))
end

function set_spatial_driver_rates!(state, loaded, time)
    for driver in values(loaded.drivers)
        if driver isa SpatialRotationalMotionGenerator
            alpha, omega, theta = driver.hinge.rotation_variables
            state[theta] = driver.motion(time)
            state[omega] = driver.motion_derivative(time)
            state[alpha] = driver.motion_second_derivative(time)
        elseif driver isa SpatialTranslationalMotionGenerator
            state[driver.distance_variable] = driver.motion(time)
            state[driver.velocity_variable] = driver.motion_derivative(time)
            state[driver.acceleration_variable] =
                driver.motion_second_derivative(time)
        elseif driver isa SpatialSpanningMotionGenerator
            state[driver.span.distance_variable] = driver.motion(time)
            state[driver.span.velocity_variable] =
                driver.motion_derivative(time)
            state[driver.span.acceleration_variable] =
                driver.motion_second_derivative(time)
        end
    end
    state
end

function initialize_spatial_measurements_and_forces!(state, loaded, time)
    for measure in values(loaded.measures)
        if measure isa SpatialSpanMeasure
            initialize_spatial_span_measure!(state, measure)
        elseif measure isa SpatialDirectedDistanceMeasure
            initialize_spatial_directed_distance_measure!(state, measure)
        end
    end
    for force in values(loaded.forces)
        if force isa SpatialBeltSpanComponent
            initialize_spatial_belt_span!(state, force)
        elseif force isa SpatialSpanningForceComponent
            initialize_spatial_spanning_force!(state, force, time)
        elseif force isa SpatialAppliedForceComponent
            initialize_spatial_applied_force!(state, force, time)
        elseif force isa SpatialAppliedTorqueComponent
            initialize_spatial_applied_torque!(state, force, time)
        elseif force isa SpatialBushingComponent
            initialize_spatial_bushing!(state, force)
        elseif force isa SpatialPlaneContactComponent
            initialize_spatial_plane_contact!(state, force)
        elseif force isa SpatialSurfaceFriction
            initialize_spatial_surface_friction!(state, force)
        elseif force isa SpatialRevoluteFriction
            initialize_spatial_revolute_friction!(state, force)
        elseif force isa SpatialTranslationalFriction
            initialize_spatial_translational_friction!(state, force)
        elseif force isa SpatialInplaneFriction
            initialize_spatial_inplane_friction!(state, force)
        elseif force isa SpatialTireComponent
            initialize_spatial_tire!(state, force, time)
        end
    end
    state
end

"""
Set the analysis stage used by loads with stage-dependent activation.

When `state` is supplied, all affected load and measurement variables are
reinitialized immediately so the canonical vector agrees with the new stage.
"""
function set_spatial_analysis_stage!(loaded, stage; state = nothing,
        time = loaded.simulation.start_time)
    stage in (:static, :dynamic, :modal) ||
        throw(ArgumentError("unknown spatial analysis stage '$stage'"))
    for force in values(loaded.forces)
        if force isa SpatialBushingComponent
            set_spatial_bushing_stage!(force, stage)
        elseif force isa SpatialAppliedForceComponent
            set_spatial_applied_force_stage!(force, stage)
        elseif force isa SpatialAppliedTorqueComponent
            set_spatial_applied_torque_stage!(force, stage)
        elseif force isa SpatialSpanningForceComponent
            set_spatial_spanning_force_stage!(force, stage)
        elseif force isa SpatialPlaneContactComponent
            set_spatial_plane_contact_stage!(force, stage)
        elseif force isa SpatialSurfaceFriction
            set_spatial_surface_friction_stage!(force, stage)
        elseif force isa SpatialRevoluteFriction
            set_spatial_revolute_friction_stage!(force, stage)
        elseif force isa SpatialTranslationalFriction
            set_spatial_translational_friction_stage!(force, stage)
        elseif force isa SpatialInplaneFriction
            set_spatial_inplane_friction_stage!(force, stage)
        end
    end
    for component in values(loaded.equation_components)
        set_spatial_equation_stage!(component, stage)
    end
    isnothing(state) ||
        initialize_spatial_measurements_and_forces!(state, loaded, time)
    loaded
end

"""
Prepare a position-consistent state for dynamic or modal equations.

The routine optionally restores declared velocities after a static solution,
projects velocity constraints, refreshes relative coordinates and measures,
solves acceleration/reaction consistency, and finally refreshes all applied
loads. It modifies and returns `state`.
"""
function prepare_spatial_dynamic_state!(state, loaded, time;
        restore_declared_velocities = false)
    if restore_declared_velocities
        for variable in loaded.layout.catalog.variables
            variable.kind in (:velocity, :angular_velocity,
                :relative_velocity) || continue
            state[variable.index] = loaded.initial_values[variable.index]
        end
    end
    set_spatial_driver_rates!(state, loaded, time)
    SpatialModelIO.weighted_constraint_projection!(state, loaded.model,
        loaded.state_selection.velocity_indices,
        loaded.state_selection.velocity_rows,
        loaded.initial_variable_weights,
        Set(loaded.imposed_initial_variable_indices), loaded.layout,
        "spatial velocity after static equilibrium"; time)
    for connection in values(loaded.connections)
        hinge = connection_hinge(connection)
        isnothing(hinge) || initialize_hinge_coordinates!(state, hinge, 1)
        inline = connection_inline(connection)
        isnothing(inline) || initialize_inline_coordinates!(state, inline, 1)
        connection isa SpatialGearPair &&
            initialize_spatial_gear_coordinates!(state, connection, 1)
    end
    initialize_spatial_measurements_and_forces!(state, loaded, time)
    inactive_variables = setdiff(
        collect(eachindex(loaded.layout.catalog.variables)),
        loaded.active_variable_indices)
    inactive_equations = setdiff(
        collect(eachindex(loaded.layout.catalog.equations)),
        loaded.active_equation_indices)
    SpatialModelIO.initialize_spatial_accelerations!(
        state, loaded.model, loaded.layout, time;
        inactive_variables, inactive_equations)
    for connection in values(loaded.connections)
        hinge = connection_hinge(connection)
        isnothing(hinge) || initialize_hinge_coordinates!(state, hinge, 2)
        inline = connection_inline(connection)
        isnothing(inline) || initialize_inline_coordinates!(state, inline, 2)
        connection isa SpatialGearPair &&
            initialize_spatial_gear_coordinates!(state, connection, 2)
    end
    initialize_spatial_measurements_and_forces!(state, loaded, time)
    state
end

"""
Advance the active spatial canonical system with the common BDF integrator.

All active algebraic and differential variables remain Newton unknowns.
Mechanical orientation partials are mapped from four Euler-parameter columns
to three pseudo-angle columns before sparse factorization. Accepted output
times are interpolated from the corrected BDF history and may be reported
incrementally through `sample_progress`.
"""
function run_spatial_implicit_model(loaded, times;
        initial_state = nothing, fixed_model_time = nothing,
        static_initialization_iterations = nothing,
        static_relaxation_cycles = nothing,
        static_mass_regularized = nothing,
        accepted_step! = nothing, accepted_step_filter! = nothing,
        sample_progress = nothing,
        maximum_step = nothing, maximum_order = nothing,
        relative_tolerance = nothing, absolute_tolerance = nothing,
        error_control_derivatives = true)
    # The model is consistent in canonical storage before the compact active
    # system and its state-equation partition are passed to DASSL.
    start, finish = first(times), last(times)
    model_time(time) = isnothing(fixed_model_time) ? time : fixed_model_time
    state = isnothing(initial_state) ? copy(loaded.initial_values) :
        copy(initial_state)
    isnothing(sample_progress) || sample_progress((;
        kind = :sample, time = Float64(start), state = copy(state)))
    derivative = initial_spatial_derivative(state, loaded, model_time(start))
    active_variables = loaded.active_variable_indices
    partition = runtime_spatial_state_partition(loaded)
    selection = AnalysisSelection(Dynamics(), loaded.active_variable_indices,
        partition.equation_indices)
    initial_equations = zeros(length(partition.equation_indices))
    evaluate_analysis_equations!(initial_equations, loaded.model, selection,
        model_time(start), state, derivative)
    consistency_limit = max(1.0e-8,
        isnothing(absolute_tolerance) ? loaded.simulation.absolute_tolerance :
        Float64(absolute_tolerance))
    norm(initial_equations, Inf) <= consistency_limit || error(
        "spatial initial values are not consistent: maximum equation error " *
        string(norm(initial_equations, Inf)))

    # Component callbacks continue to use canonical indices. These closures
    # unpack DASSL's active vectors and map orientation partials into the
    # pseudo-angle Newton columns.
    canonical = copy(state)
    canonical_derivative = zeros(length(state))
    residual! = function (equations, time, values, rates, parameter)
        canonical[active_variables] .= values
        canonical_derivative .= 0
        canonical_derivative[active_variables] .= rates
        evaluate_analysis_equations!(equations, loaded.model, selection,
            model_time(time), canonical, canonical_derivative)
    end
    workspace = spatial_dynamic_jacobian_workspace(loaded, selection,
        model_time(start), state, derivative, 1.0)
    jacobian! = function (matrix, time, values, rates, coefficient, parameter)
        canonical[active_variables] .= values
        canonical_derivative .= 0
        canonical_derivative[active_variables] .= rates
        evaluate_spatial_dynamic_jacobian!(matrix, workspace, loaded,
            selection, model_time(time), canonical, canonical_derivative,
            coefficient)
    end
    prototype = copy(workspace.prototype)
    evaluate_spatial_dynamic_jacobian!(prototype, workspace, loaded,
        selection, model_time(start), state, derivative, 1.0)
    local_variable = Dict(canonical => local_index
        for (local_index, canonical) in enumerate(active_variables))
    predictor_indices = [(;
        parameters = Tuple(local_variable[index]
            for index in body.euler_parameter_variables),
        pseudo_angles = Tuple(local_variable[index]
            for index in body.pseudo_angle_variables))
        for body in values(loaded.bodies)]
    predictor! = function (predicted, time, coefficient,
            history_derivative, parameter)
        project_spatial_orientation_predictor!(predicted,
            history_derivative, coefficient, predictor_indices)
    end
    differential, error_control, physical_monitor =
        spatial_state_masks(loaded, partition)
    if !error_control_derivatives
        for (local_index, canonical_index) in enumerate(active_variables)
            kind = loaded.layout.catalog.variables[canonical_index].kind
            kind in SPATIAL_STATIC_RATE_KINDS || continue
            error_control[local_index] = false
            physical_monitor[local_index] = false
        end
    end
    state_selection_changes = NamedTuple[(;
        time = Float64(start), reason = :initial,
        previous = Symbol[], selected = copy(partition.selected_names))]

    function partition_condition(selected_names, diagnostics)
        selected_velocities = Set(partition.candidates[name].velocity
            for name in selected_names)
        dependent_columns = findall(index ->
            !(index in selected_velocities),
            loaded.state_selection.velocity_indices)
        length(dependent_columns) == diagnostics.rank || return 0.0
        singular_values = svdvals(
            diagnostics.scaled_partial[:, dependent_columns])
        isempty(singular_values) && return Inf
        first(singular_values) > 0 || return 0.0
        last(singular_values) / first(singular_values)
    end

    # Reselection changes equivalent state-equation rows only. The complete
    # physical state and BDF history stay intact; the changed sparse pattern
    # receives a new symbolic factorization.
    function reselect_states!(reason, time, values, rates,
            internal_error_control, internal_differential, parameter)
        canonical[active_variables] .= values
        canonical_derivative .= 0
        canonical_derivative[active_variables] .= rates
        choice, diagnostics = automatic_runtime_spatial_state_selection(
            loaded, canonical, model_time(time))
        length(choice) == length(partition.selected_names) || return nothing
        Set(choice) == Set(partition.selected_names) && return nothing
        previous = copy(partition.selected_names)
        current_condition = partition_condition(previous, diagnostics)
        proposed_condition = partition_condition(choice, diagnostics)
        reason != :singular_iteration_matrix &&
            proposed_condition <
                STATE_RESELECTION_IMPROVEMENT * current_condition &&
            return nothing
        select_runtime_spatial_states!(partition, choice)
        new_differential, new_error_control, _ =
            spatial_state_masks(loaded, partition)
        internal_differential .= new_differential
        internal_error_control .= new_error_control
        push!(state_selection_changes, (; time = Float64(time), reason,
            previous, selected = copy(partition.selected_names),
            rank = diagnostics.rank, current_condition,
            selected_condition = proposed_condition))
        workspace = spatial_dynamic_jacobian_workspace(loaded, selection,
            model_time(time), canonical, canonical_derivative, 1.0)
        new_prototype = copy(workspace.prototype)
        evaluate_spatial_dynamic_jacobian!(new_prototype, workspace, loaded,
            selection, model_time(time), canonical, canonical_derivative, 1.0)
        new_equation_levels = getproperty.(
            loaded.layout.catalog.equations[partition.equation_indices],
            :level)
        (; jacobian_prototype = new_prototype,
           equation_levels = new_equation_levels)
    end
    settings = loaded.simulation
    integration_maximum_step = isnothing(maximum_step) ?
        settings.maximum_step : Float64(maximum_step)
    integration_maximum_order = isnothing(maximum_order) ? 5 :
        Int(maximum_order)
    integration_relative_tolerance = isnothing(relative_tolerance) ?
        settings.relative_tolerance : Float64(relative_tolerance)
    integration_absolute_tolerance = isnothing(absolute_tolerance) ?
        settings.absolute_tolerance : Float64(absolute_tolerance)
    options = HistoricalDDASSL.DASSLOptions{Float64}(
        atol = integration_absolute_tolerance,
        rtol = integration_relative_tolerance,
        initial_step = settings.initial_step,
        maximum_step = integration_maximum_step,
        maximum_order = integration_maximum_order)
    # Tire and compliant-contact forces are continuous at their switching
    # surfaces. Let the ordinary BDF corrector and error controller cross
    # those constitutive-law corners without locating an event or forcing a
    # history restart. This is especially important when several stiff vehicle
    # contacts switch in close succession.
    root_count = 0
    root_callback = nothing

    # Requested result frames are interpolated from the corrected BDF history
    # instead of constraining the internal step sequence.
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
    last_accepted_time = Ref(Float64(start))
    last_accepted_values = copy(state[active_variables])
    last_accepted_statistics = Ref{Any}(nothing)
    tracked_accepted_step! = if isnothing(sample_progress) &&
            isnothing(accepted_step!)
        nothing
    else
        function (time, values, rates, order, step, statistics)
            if !isnothing(sample_progress)
                last_accepted_values .= values
                last_accepted_time[] = Float64(time)
                last_accepted_statistics[] = statistics
            end
            isnothing(accepted_step!) || accepted_step!(time, values, rates,
                order, step, statistics)
        end
    end
    solution = try
        HistoricalDDASSL.dassl(residual!, state[active_variables],
            derivative[active_variables],
            (start, finish); jacobian!, jacobian_prototype = prototype, options,
            variable_levels = getproperty.(
                loaded.layout.catalog.variables[active_variables], :level),
            equation_levels = getproperty.(
                loaded.layout.catalog.equations[partition.equation_indices], :level),
            differential_vars = differential, error_control,
            error_monitor = physical_monitor, deficit = loaded.analysis.deficit,
            root! = root_callback, number_of_roots = root_count,
            root_restart = :soft,
            reconfigure! = reselect_states!, predictor!,
            accepted_step! = tracked_accepted_step!, accepted_interval!,
            accepted_step_filter!)
    catch exception
        if !isnothing(sample_progress)
            expanded = copy(state)
            expanded[active_variables] .= last_accepted_values
            sample_progress((; kind = :failure,
                time = last_accepted_time[], state = expanded,
                message = exception isa InterruptException ?
                    "spatial dynamic integration interrupted" :
                    "spatial dynamic integration stopped: " *
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
            message = "spatial dynamic integration failed: " *
                solution.message,
            statistics = solution.stats))
        error("spatial dynamic integration failed: $(solution.message)")
    end
    active_states = solution(times)
    states = [begin
        expanded = copy(state)
        expanded[active_variables] .= values
        expanded
    end for values in active_states]
    initial_consistency_iterations =
        loaded.initial_conditions.position_corrections +
        loaded.initial_conditions.velocity_corrections +
        loaded.initial_conditions.acceleration_corrections
    (; times = collect(times), states, initial_consistency_iterations, loaded,
       static_initialization_iterations,
       static_relaxation_cycles, static_mass_regularized,
       analysis_mode = :dynamic, solution, state_selection_changes)
end

"""
Approach spatial equilibrium by repeatedly integrating first-order BDF
pseudo-time cycles while depleting rates after every accepted step.

The physical model time remains fixed. Handoff to the exact static Newton
solve occurs only after the configured minimum cycles and force, acceleration,
speed, and configuration-correction tests are satisfied.
"""
function dynamic_relaxation_spatial_equilibrium(loaded, system, time,
        initial; progress = nothing)
    state = copy(initial)
    zero_spatial_static_rates!(state, loaded)
    settings = loaded.analysis
    pseudo_times = [0.0, settings.relaxation_duration]
    for cycle in 1:settings.relaxation_max_cycles
        prepare_spatial_dynamic_state!(state, loaded, time)
        diagnostics = spatial_static_diagnostics(loaded, system, time, state)
        emit_static_progress(progress, (;
            phase = :dynamic_relaxation, status = :initial,
            loaded, state = copy(state), model_time = Float64(time),
            pseudo_time = 0.0, iteration = 0, relaxation_cycle = cycle,
            diagnostics..., reciprocal_condition = NaN,
            mass_regularized = false,
            message = "starting relaxation cycle"))
        last_display_time = Ref(-Inf)
        local_derivative_variables = [local_index for
            (local_index, canonical_index) in
                enumerate(loaded.active_variable_indices)
            if loaded.layout.catalog.variables[canonical_index].kind in
                SPATIAL_STATIC_RATE_KINDS]
        accepted_step_filter! = function (pseudo_time, values, rates,
                order, step, statistics)
            step_fraction = abs(step) / settings.relaxation_duration
            step_reduction = settings.relaxation_reduction_factor^
                step_fraction
            values[local_derivative_variables] .*= step_reduction
            rates .*= step_reduction
            nothing
        end
        accepted_step! = if isnothing(progress)
            nothing
        else
            function (pseudo_time, values, rates, order, step, statistics)
                display_interval = 1 / 30
                pseudo_time + 100eps(Float64) < settings.relaxation_duration &&
                    pseudo_time - last_display_time[] < display_interval &&
                    return
                last_display_time[] = pseudo_time
                snapshot = copy(state)
                snapshot[loaded.active_variable_indices] .= values
                diagnostics = spatial_static_diagnostics(
                    loaded, system, time, snapshot)
                emit_static_progress(progress, (;
                    phase = :dynamic_relaxation, status = :accepted,
                    loaded, state = snapshot, model_time = Float64(time),
                    pseudo_time = Float64(pseudo_time), iteration =
                        statistics.accepted_steps,
                    relaxation_cycle = cycle, diagnostics...,
                    reciprocal_condition = NaN, mass_regularized = false,
                    message = "BDF order $order, step $(Float64(step))"))
            end
        end
        relaxation = run_spatial_implicit_model(loaded, pseudo_times;
            initial_state = state, fixed_model_time = time, accepted_step!,
            accepted_step_filter!, maximum_step = Inf, maximum_order = 1,
            relative_tolerance = max(
                loaded.simulation.relative_tolerance, 1.0e-3),
            absolute_tolerance = max(
                loaded.simulation.absolute_tolerance, 1.0e-6),
            error_control_derivatives = false)
        state = last(relaxation.states)
        readiness = spatial_static_relaxation_readiness(
            loaded, system, time, state)
        acceleration_text = string(round(
            readiness.equivalent_acceleration; sigdigits = 4))
        speed_text = string(round(readiness.equivalent_speed; sigdigits = 4))
        correction_text = string(round(
            readiness.configuration_correction; sigdigits = 4))
        minimum_cycles_complete = cycle >= settings.relaxation_min_cycles
        handoff_ready = minimum_cycles_complete && readiness.ready
        minimum_text = minimum_cycles_complete ? "" :
            ", minimum cycle $cycle of $(settings.relaxation_min_cycles)"
        readiness_message = "handoff " *
            (handoff_ready ? "ready" : "waiting") *
            ": equivalent acceleration $acceleration_text m/s^2, " *
            "speed $speed_text m/s, configuration correction " *
            correction_text * minimum_text
        emit_static_progress(progress, (;
            phase = :dynamic_relaxation,
            status = handoff_ready ? :ready : :waiting,
            loaded, state = copy(state), model_time = Float64(time),
            pseudo_time = settings.relaxation_duration,
            iteration = relaxation.solution.stats.accepted_steps,
            relaxation_cycle = cycle,
            residual = readiness.residual,
            force_imbalance = readiness.force_imbalance,
            torque_imbalance = readiness.torque_imbalance,
            constraint_error = readiness.constraint_error,
            equivalent_acceleration = readiness.equivalent_acceleration,
            reciprocal_condition = NaN, mass_regularized = false,
            message = readiness_message))
        handoff_ready || continue
        if !settings.relaxation_polish
            emit_static_progress(progress, (;
                phase = :dynamic_relaxation, status = :converged,
                loaded, state = copy(state), model_time = Float64(time),
                pseudo_time = settings.relaxation_duration,
                iteration = relaxation.solution.stats.accepted_steps,
                relaxation_cycle = cycle,
                residual = readiness.residual,
                force_imbalance = readiness.force_imbalance,
                torque_imbalance = readiness.torque_imbalance,
                constraint_error = readiness.constraint_error,
                equivalent_acceleration = readiness.equivalent_acceleration,
                reciprocal_condition = NaN, mass_regularized = false,
                message = "relaxed solution accepted without Newton polish"))
            return state, 0, cycle, false
        end
        try
            equilibrium, iterations, mass_regularized = spatial_static_equilibrium(
                loaded, system, time, state; maximum_iterations = 12,
                progress, phase = :newton_polish, relaxation_cycle = cycle)
            return equilibrium, iterations, cycle, mass_regularized
        catch error
            error isa SpatialStaticConvergenceError || rethrow()
        end
    end
    throw(SpatialStaticConvergenceError(
        "spatial dynamic relaxation did not reach static equilibrium at " *
        "time $time after $(settings.relaxation_max_cycles) cycles"))
end

function solve_spatial_static_equilibrium(loaded, system, time, initial;
        progress = nothing)
    if loaded.analysis.static_method == :newton
        state, iterations, mass_regularized = spatial_static_equilibrium(
            loaded, system, time, initial; progress)
        return state, iterations, 0, mass_regularized
    end
    dynamic_relaxation_spatial_equilibrium(
        loaded, system, time, initial; progress)
end

function spatial_statically_initialized_state(loaded, time;
        progress = nothing)
    set_spatial_analysis_stage!(loaded, :static;
        state = loaded.initial_values, time)
    system = spatial_static_system(loaded)
    emit_spatial_initial_configurations(progress, loaded, system, time)
    state, iterations, relaxation_cycles, mass_regularized =
        solve_spatial_static_equilibrium(
            loaded, system, time, loaded.initial_values; progress)
    set_spatial_analysis_stage!(loaded, :dynamic)
    prepare_spatial_dynamic_state!(state, loaded, time;
        restore_declared_velocities = true)
    state, iterations, relaxation_cycles, mass_regularized
end

function spatial_modal_operating_point(loaded, time; progress = nothing,
        initial_state = nothing)
    if !isnothing(initial_state)
        state = copy(initial_state)
        set_spatial_analysis_stage!(loaded, :modal; state, time)
        prepare_spatial_dynamic_state!(state, loaded, time)
        derivative = initial_spatial_derivative(state, loaded, time)
        return (; state, derivative, initial_iterations = 0,
            static_initialization_iterations = nothing,
            static_relaxation_cycles = nothing,
            static_mass_regularized = nothing)
    end
    set_spatial_analysis_stage!(loaded,
        loaded.analysis.initialization == :static_equilibrium ?
            :static : :modal; state = loaded.initial_values, time)
    state, static_initialization_iterations, static_relaxation_cycles,
        static_mass_regularized =
        if loaded.analysis.initialization == :static_equilibrium
            system = spatial_static_system(loaded)
            emit_spatial_initial_configurations(
                progress, loaded, system, time)
            equilibrium, iterations, cycles, mass_regularized =
                solve_spatial_static_equilibrium(
                    loaded, system, time, loaded.initial_values; progress)
            (equilibrium, iterations, cycles, mass_regularized)
        else
            (copy(loaded.initial_values), nothing, nothing, nothing)
        end
    set_spatial_analysis_stage!(loaded, :modal)
    prepare_spatial_dynamic_state!(state, loaded, time)
    derivative = initial_spatial_derivative(state, loaded, time)
    initial_iterations = loaded.initial_conditions.position_corrections +
        loaded.initial_conditions.velocity_corrections +
        loaded.initial_conditions.acceleration_corrections
    (; state, derivative, initial_iterations,
       static_initialization_iterations, static_relaxation_cycles,
       static_mass_regularized)
end

"""
Build the square tangent-coordinate Jacobian evaluator used for modal analysis.

Euler-parameter kinematic equations and their four variables are removed from
the eigenproblem. Three pseudo-angle columns per body retain the mechanical
orientation partials at the supplied finite operating orientation.
"""
function spatial_modal_linearization(loaded, state, derivative, time)
    # Euler parameters define the finite operating orientation, but their four
    # perturbations are replaced by the three body-fixed pseudo-angle columns.
    parameter_variables = Set(Iterators.flatten(
        body.euler_parameter_variables for body in values(loaded.bodies)))
    orientation_equations = Set(Iterators.flatten(
        body.orientation_equations for body in values(loaded.bodies)))
    variable_indices = [index for index in loaded.active_variable_indices
        if index ∉ parameter_variables]
    equation_indices = [index for index in loaded.active_equation_indices
        if index ∉ orientation_equations]
    length(variable_indices) == length(equation_indices) ||
        throw(DimensionMismatch(
            "spatial modal variable and equation sets must be square"))

    full_selection = AnalysisSelection(Dynamics(),
        loaded.active_variable_indices, loaded.active_equation_indices)
    workspace = spatial_dynamic_jacobian_workspace(loaded, full_selection,
        time, state, derivative, 1.0)
    matrix = copy(workspace.prototype)
    full_variable_lookup = Dict(canonical => local_index
        for (local_index, canonical) in
        enumerate(full_selection.variable_indices))
    full_equation_lookup = Dict(canonical => local_index
        for (local_index, canonical) in
        enumerate(full_selection.equation_indices))
    local_columns = [full_variable_lookup[index] for index in variable_indices]
    local_rows = [full_equation_lookup[index] for index in equation_indices]
    jacobian_evaluator = function (coefficient)
        evaluate_spatial_dynamic_jacobian!(matrix, workspace, loaded,
            full_selection, time, state, derivative, coefficient)
        copy(matrix[local_rows, local_columns])
    end
    (; variable_indices, equation_indices, jacobian_evaluator)
end

function run_spatial_modal_model(loaded, time; static_progress = nothing,
        initial_state = nothing)
    operating_point = spatial_modal_operating_point(loaded, time;
        progress = static_progress, initial_state)
    linearization = spatial_modal_linearization(loaded,
        operating_point.state, operating_point.derivative, time)
    modes = solve_modal_system(loaded, operating_point.state,
        operating_point.derivative, time;
        variable_indices = linearization.variable_indices,
        equation_indices = linearization.equation_indices,
        jacobian_evaluator = linearization.jacobian_evaluator)

    mode_shapes = copy(modes.mode_shapes)
    # Store the tangent Euler-parameter perturbation so the ordinary spatial
    # viewer can animate the mode. It was not an eigenproblem coordinate.
    for body in values(loaded.bodies), mode in axes(mode_shapes, 2)
        pseudo_angle = @view mode_shapes[body.pseudo_angle_variables, mode]
        parameters = @view operating_point.state[
            body.euler_parameter_variables]
        mode_shapes[body.euler_parameter_variables, mode] .=
            0.5 .* quaternion_rate_matrix(parameters) * pseudo_angle
    end
    modes = merge(modes, (; mode_shapes))
    (; times = [Float64(time)], states = [operating_point.state],
       operating_state = operating_point.state,
       operating_derivative = operating_point.derivative,
       initial_consistency_iterations = operating_point.initial_iterations,
       loaded, analysis_mode = :modal,
       static_initialization_iterations =
           operating_point.static_initialization_iterations,
       static_relaxation_cycles =
           operating_point.static_relaxation_cycles,
       static_mass_regularized = operating_point.static_mass_regularized,
       modes...)
end

"""Linearly predict a quasi-static target and renormalize every orientation."""
function spatial_static_predict(times, states, target, loaded)
    length(states) < 2 && return copy(last(states))
    first_time, last_time = times[end - 1], times[end]
    last_time == first_time && return copy(last(states))
    prediction = last(states) .+ ((target - last_time) /
        (last_time - first_time)) .* (states[end] .- states[end - 1])
    normalize_spatial_orientations!(prediction, loaded)
end

"""Advance one quasi-static target, bisecting model time after a failed solve."""
function advance_spatial_static!(accepted_times, accepted_states, target,
        loaded, system, iterations, relaxation_cycles, mass_regularized;
        minimum_step, progress = nothing)
    prediction = spatial_static_predict(
        accepted_times, accepted_states, target, loaded)
    try
        state, count, cycles, regularized = solve_spatial_static_equilibrium(
            loaded, system, target, prediction; progress)
        push!(accepted_times, target)
        push!(accepted_states, state)
        push!(iterations, count)
        push!(relaxation_cycles, cycles)
        push!(mass_regularized, regularized)
        return state
    catch error
        error isa SpatialStaticConvergenceError || rethrow()
        previous = last(accepted_times)
        abs(target - previous) > minimum_step || rethrow()
        midpoint = previous + (target - previous) / 2
        advance_spatial_static!(accepted_times, accepted_states, midpoint,
            loaded, system, iterations, relaxation_cycles, mass_regularized;
            minimum_step, progress)
        advance_spatial_static!(accepted_times, accepted_states, target,
            loaded, system, iterations, relaxation_cycles, mass_regularized;
            minimum_step, progress)
    end
end

function run_spatial_static_model(loaded, times; static_progress = nothing)
    set_spatial_analysis_stage!(loaded, :static;
        state = loaded.initial_values, time = first(times))
    system = spatial_static_system(loaded)
    emit_spatial_initial_configurations(
        static_progress, loaded, system, first(times))
    first_state, first_iterations, first_relaxation_cycles,
        first_mass_regularized =
        solve_spatial_static_equilibrium(
            loaded, system, first(times), loaded.initial_values;
            progress = static_progress)
    accepted_times = Float64[first(times)]
    accepted_states = [first_state]
    all_iterations = Int[first_iterations]
    all_relaxation_cycles = Int[first_relaxation_cycles]
    all_mass_regularized = Bool[first_mass_regularized]
    output_states = [copy(first_state)]
    interval = max(abs(last(times) - first(times)), 1.0)
    minimum_step = sqrt(eps(Float64)) * interval
    for target in Iterators.drop(times, 1)
        state = advance_spatial_static!(accepted_times, accepted_states,
            target, loaded, system, all_iterations,
            all_relaxation_cycles, all_mass_regularized; minimum_step,
            progress = static_progress)
        push!(output_states, copy(state))
    end
    (; times = collect(times), states = output_states,
       initial_consistency_iterations = first_iterations, loaded,
       analysis_mode = :static, static_iterations = all_iterations,
       static_relaxation_cycles = all_relaxation_cycles,
       static_mass_regularized = all_mass_regularized,
       static_variable_indices = system.result_variables,
       static_equation_indices = system.equation_indices,
       continuation_solutions = length(accepted_times))
end

"""
    run_spatial_model(source; start_time=nothing, end_time=nothing,
                      duration=nothing, samples=nothing)

Load and run a spatial TOML model. Automatic analysis selects dynamic
integration. `analysis.mode = "static"` solves one equilibrium or a
time-stepped sequence of equilibria. Modal analysis linearizes the sparse
implicit equations at the consistent operating point. A dynamic or modal model
may use `analysis.initialization = "static_equilibrium"`. Dynamics then restores
declared velocities before BDF integration; modal analysis retains the
stationary equilibrium.
"""
function run_spatial_model(source; start_time = nothing, end_time = nothing,
        duration = nothing, samples = nothing, static_progress = nothing,
        result_progress = nothing)
    loaded = source isa LoadedSpatialModel ? source : load_spatial_model(source)
    static_progress_events = Any[]
    progress = function (event)
        push!(static_progress_events, event)
        isnothing(static_progress) || static_progress(event)
        isnothing(result_progress) || result_progress(
            merge((; kind = :static), event))
    end
    with_static_progress(result) =
        merge(result, (; static_progress_events))
    settings = loaded.simulation
    start = isnothing(start_time) ? settings.start_time : Float64(start_time)
    finish = if !isnothing(duration)
        Float64(duration) >= 0 || throw(ArgumentError(
            "duration must be nonnegative"))
        start + Float64(duration)
    else
        isnothing(end_time) ? settings.end_time : Float64(end_time)
    end
    finish >= start || throw(ArgumentError(
        "end_time must not precede start_time"))
    count = isnothing(samples) ? settings.output_samples : Int(samples)
    count >= 1 || throw(ArgumentError("samples must be positive"))
    times = count == 1 ? [start] :
        collect(range(start, finish; length = count))
    mode = loaded.analysis.mode == :automatic ? :dynamic : loaded.analysis.mode
    isnothing(result_progress) || result_progress((; kind = :begin, loaded,
        analysis_mode = mode, times = collect(times)))
    mode == :modal && loaded.analysis.degrees_of_freedom == 0 &&
        throw(ArgumentError(
            "modal analysis requires at least one independent state"))
    mode == :modal && return with_static_progress(
        run_spatial_modal_model(loaded, start; static_progress = progress))
    if mode == :static
        return with_static_progress(
            run_spatial_static_model(loaded, times;
                static_progress = progress))
    end
    if loaded.analysis.initialization != :static_equilibrium
        set_spatial_analysis_stage!(loaded, :dynamic)
        return with_static_progress(run_spatial_implicit_model(loaded, times;
            sample_progress = result_progress))
    end
    state, iterations, relaxation_cycles, mass_regularized =
        spatial_statically_initialized_state(loaded, start;
            progress)
    with_static_progress(run_spatial_implicit_model(loaded, times;
        initial_state = state,
        static_initialization_iterations = iterations,
        static_relaxation_cycles = relaxation_cycles,
        static_mass_regularized = mass_regularized,
        sample_progress = result_progress))
end

end
