"""
    AutomaticAnalysis

Canonical allocation and analysis selection for component-local implicit
equations. Components declare variables and equation blocks without knowing
global indices. This module assigns those indices, records engineering
metadata, selects the variables and equations required by an analysis, and
assembles dense or sparse Jacobians from local callbacks.

Every canonical equation has exactly one owning `ExecutableEquationBlock`.
Loads and connections may additionally contribute terms to equations owned by
another component through `EquationContribution`. That distinction lets a
body own its force balance while each attached element remains local.
"""
module AutomaticAnalysis

using SparseArrays

export VariableMetadata, EquationMetadata, AnalysisCatalog,
       PositionIC, VelocityIC, AccelerationIC, StaticEQ, Dynamics,
       KinematicPosition, KinematicVelocity, KinematicAcceleration,
       KinematicForces,
       AnalysisSelection, select_analysis, analysis_values,
       expand_analysis_values!, VariableDeclaration, EquationDeclaration,
       EquationBlockDeclaration, ConstraintFamilyDeclaration,
       AllocatedConstraintFamily, ComponentRegistration,
       AnalysisCatalogBuilder, register_component_variables!,
       register_component_equation_block!, finish_catalog,
       ModelLayoutBuilder, ModelLayout, allocate_component_variables!,
       allocate_component_equation_block!, finish_layout,
       component_variable_indices, component_equation_indices,
       ExecutableEquationBlock, EquationContribution,
       ExecutableAnalysisModel, AnalysisEquationWorkspace,
       analysis_equation_workspace, configure_analysis_equation_workspace!,
       evaluate_analysis_equations!,
       evaluate_analysis_jacobian, evaluate_analysis_sparse_jacobian,
       evaluate_analysis_sparse_jacobian!

"""Rule selecting one square analysis from the unreduced canonical system."""
abstract type AnalysisPolicy end

"""Correct configuration variables against level-zero constraint equations."""
struct PositionIC <: AnalysisPolicy end

"""Correct velocity variables against level-one constraint equations."""
struct VelocityIC <: AnalysisPolicy end

"""Solve accelerations, reactions, and applied-load definitions initially."""
struct AccelerationIC <: AnalysisPolicy end

"""Solve positions, reactions, and static load definitions with inertia omitted."""
struct StaticEQ <: AnalysisPolicy end

"""Retain every active canonical variable and equation for implicit dynamics."""
struct Dynamics <: AnalysisPolicy end

"""Solve the position level of a fully prescribed kinematic mechanism."""
struct KinematicPosition <: AnalysisPolicy end

"""Solve the velocity level of a fully prescribed kinematic mechanism."""
struct KinematicVelocity <: AnalysisPolicy end

"""Solve the acceleration level of a fully prescribed kinematic mechanism."""
struct KinematicAcceleration <: AnalysisPolicy end

"""Recover reactions and applied loads after prescribed kinematics are known."""
struct KinematicForces <: AnalysisPolicy end

"""Engineering identity and scaling of one canonical scalar variable."""
struct VariableMetadata
    index::Int
    name::Symbol
    kind::Symbol
    level::Int
    component::Symbol
    scale::Float64
end

"""Engineering identity and constraint-family membership of one equation."""
struct EquationMetadata
    index::Int
    name::Symbol
    kind::Symbol
    level::Int
    component::Symbol
    family::Symbol
end

"""Component-local declaration of a scalar variable before allocation."""
struct VariableDeclaration
    name::Symbol
    kind::Symbol
    level::Int
    scale::Float64
end

VariableDeclaration(name, kind, level) =
    VariableDeclaration(name, kind, level, 1.0)

"""Component-local declaration of a scalar implicit equation before allocation."""
struct EquationDeclaration
    name::Symbol
    kind::Symbol
    level::Int
    family::Symbol
end

"""Named group of equations owned and evaluated by one component callback."""
struct EquationBlockDeclaration
    name::Symbol
    equations::Vector{EquationDeclaration}
end

"""Names linking one scalar ideal constraint across its three levels and reaction."""
struct ConstraintFamilyDeclaration
    id::Symbol
    reaction_variable::Symbol
    position_equation::Symbol
    velocity_equation::Symbol
    acceleration_equation::Symbol
end

"""Canonical indices belonging to one allocated scalar ideal constraint."""
struct AllocatedConstraintFamily
    component::Symbol
    id::Symbol
    reaction_variable::Int
    position_equation::Int
    velocity_equation::Int
    acceleration_equation::Int
end

"""
Declaration of one component's canonical variables, owned equation blocks,
and ideal-constraint families. Registration contains names and metadata only;
the corresponding executable callbacks are supplied after allocation.
"""
struct ComponentRegistration
    name::Symbol
    variables::Vector{VariableDeclaration}
    equation_blocks::Dict{Symbol,EquationBlockDeclaration}
    constraint_families::Vector{ConstraintFamilyDeclaration}
end

function ComponentRegistration(name::Symbol,
        variables::AbstractVector{<:VariableDeclaration},
        equation_blocks::AbstractVector{<:EquationBlockDeclaration},
        constraint_families::AbstractVector{<:ConstraintFamilyDeclaration} =
            ConstraintFamilyDeclaration[])
    blocks = Dict(block.name => block for block in equation_blocks)
    length(blocks) == length(equation_blocks) ||
        error("Component $name has duplicate equation-block names")
    family_ids = getproperty.(constraint_families, :id)
    length(unique(family_ids)) == length(family_ids) ||
        error("Component $name has duplicate constraint-family identifiers")
    return ComponentRegistration(name, collect(variables), blocks,
        collect(constraint_families))
end

"""First-pass allocator that assigns consecutive canonical scalar indices."""
mutable struct AnalysisCatalogBuilder
    variables::Vector{VariableMetadata}
    equations::Vector{EquationMetadata}
    registered_variable_components::Set{Symbol}
    registered_equation_blocks::Set{Tuple{Symbol,Symbol}}
end

AnalysisCatalogBuilder() = AnalysisCatalogBuilder(
    VariableMetadata[], EquationMetadata[], Set{Symbol}(),
    Set{Tuple{Symbol,Symbol}}())

"""Allocate all variables of `component` as one consecutive canonical range."""
function register_component_variables!(builder::AnalysisCatalogBuilder,
                                       component::ComponentRegistration)
    component.name in builder.registered_variable_components &&
        error("Variables for component $(component.name) are already registered")
    first_index = length(builder.variables) + 1
    for declaration in component.variables
        push!(builder.variables, VariableMetadata(
            length(builder.variables) + 1, declaration.name,
            declaration.kind, declaration.level, component.name,
            declaration.scale))
    end
    push!(builder.registered_variable_components, component.name)
    return first_index:length(builder.variables)
end

"""Allocate one owned equation block as a consecutive canonical range."""
function register_component_equation_block!(builder::AnalysisCatalogBuilder,
        component::ComponentRegistration, block_name::Symbol)
    key = (component.name, block_name)
    key in builder.registered_equation_blocks &&
        error("Equation block $key is already registered")
    haskey(component.equation_blocks, block_name) ||
        error("Component $(component.name) has no equation block $block_name")
    block = component.equation_blocks[block_name]
    first_index = length(builder.equations) + 1
    for declaration in block.equations
        push!(builder.equations, EquationMetadata(
            length(builder.equations) + 1, declaration.name,
            declaration.kind, declaration.level, component.name,
            declaration.family))
    end
    push!(builder.registered_equation_blocks, key)
    return first_index:length(builder.equations)
end

finish_catalog(builder::AnalysisCatalogBuilder) =
    AnalysisCatalog(builder.variables, builder.equations)

"""Collect canonical ranges by component name while constructing a model."""
mutable struct ModelLayoutBuilder
    catalog_builder::AnalysisCatalogBuilder
    variable_indices::Dict{Symbol,UnitRange{Int}}
    equation_indices::Dict{Tuple{Symbol,Symbol},UnitRange{Int}}
end

ModelLayoutBuilder() = ModelLayoutBuilder(AnalysisCatalogBuilder(),
    Dict{Symbol,UnitRange{Int}}(),
    Dict{Tuple{Symbol,Symbol},UnitRange{Int}}())

function allocate_component_variables!(builder::ModelLayoutBuilder,
                                       component::ComponentRegistration)
    indices = register_component_variables!(builder.catalog_builder, component)
    builder.variable_indices[component.name] = indices
    return indices
end

function allocate_component_equation_block!(builder::ModelLayoutBuilder,
        component::ComponentRegistration, block_name::Symbol)
    indices = register_component_equation_block!(
        builder.catalog_builder, component, block_name)
    builder.equation_indices[(component.name, block_name)] = indices
    return indices
end

"""
Completed mapping from component-local declarations to canonical indices.

The canonical catalog is deliberately larger than any one analysis. It keeps
all levels of variables and equations, including inactive candidate state
equations, so an analysis or runtime state change can select a different square
subsystem without reallocating model data.
"""
struct ModelLayout{C}
    catalog::C
    variable_indices::Dict{Symbol,UnitRange{Int}}
    equation_indices::Dict{Tuple{Symbol,Symbol},UnitRange{Int}}
    constraint_families::Dict{Tuple{Symbol,Symbol},AllocatedConstraintFamily}
end

"""Freeze a completed layout and resolve every ideal constraint family."""
function finish_layout(builder::ModelLayoutBuilder,
        registrations = ComponentRegistration[])
    catalog = finish_catalog(builder.catalog_builder)
    families = Dict{Tuple{Symbol,Symbol},AllocatedConstraintFamily}()
    for registration in registrations
        for declaration in registration.constraint_families
            variable_matches = [variable for variable in catalog.variables
                if variable.component == registration.name &&
                   variable.name == declaration.reaction_variable]
            length(variable_matches) == 1 || error(
                "Constraint family $(registration.name).$(declaration.id) " *
                "must name one reaction variable")
            variable_matches[1].kind == :reaction || error(
                "Constraint family reaction must have kind :reaction")
            function equation_index(name, level)
                matches = [equation for equation in catalog.equations
                    if equation.component == registration.name &&
                       equation.name == name]
                length(matches) == 1 || error(
                    "Constraint family $(registration.name).$(declaration.id) " *
                    "must name one $level-level equation")
                matches[1].kind == :constraint && matches[1].level == level ||
                    error("Constraint family equation $name has the wrong kind or level")
                matches[1].index
            end
            family = AllocatedConstraintFamily(registration.name,
                declaration.id, variable_matches[1].index,
                equation_index(declaration.position_equation, 0),
                equation_index(declaration.velocity_equation, 1),
                equation_index(declaration.acceleration_equation, 2))
            families[(registration.name, declaration.id)] = family
        end
    end
    ModelLayout(catalog, copy(builder.variable_indices),
        copy(builder.equation_indices), families)
end

component_variable_indices(layout::ModelLayout, component::Symbol) =
    layout.variable_indices[component]
component_equation_indices(layout::ModelLayout, component::Symbol,
                           block::Symbol) =
    layout.equation_indices[(component, block)]

"""Complete ordered metadata for the unreduced variables and equations."""
struct AnalysisCatalog
    variables::Vector{VariableMetadata}
    equations::Vector{EquationMetadata}
    function AnalysisCatalog(variables, equations)
        variable_indices = sort([variable.index for variable in variables])
        equation_indices = sort([equation.index for equation in equations])
        variable_indices == collect(1:length(variables)) ||
            error("Variable metadata must cover consecutive canonical indices")
        equation_indices == collect(1:length(equations)) ||
            error("Equation metadata must cover consecutive canonical indices")
        new(collect(variables), collect(equations))
    end
end

"""One ordered square view of the canonical variables and equations."""
struct AnalysisSelection{P<:AnalysisPolicy}
    policy::P
    variable_indices::Vector{Int}
    equation_indices::Vector{Int}
end

"""
An equation block owned by one component.

`residual!` and `jacobian!` write using canonical indices. The block owns each
index in `equation_indices`; other components may add contributions to those
same equations but may not own them.
"""
struct ExecutableEquationBlock{R,J}
    component::Symbol
    name::Symbol
    equation_indices::Vector{Int}
    residual!::R
    jacobian!::J
end


"""
Component-local additive terms applied to equations owned elsewhere.

Typical examples are a spring force added to two body balances and a joint
reaction added to the balances of its connected bodies.
"""
struct EquationContribution{R,J}
    component::Symbol
    name::Symbol
    target_equation_indices::Vector{Int}
    residual!::R
    jacobian!::J
end


"""
Executable canonical model made from owned equation blocks and additive
contributions. Construction verifies that every canonical equation is owned
once and only once.
"""
struct ExecutableAnalysisModel
    catalog::AnalysisCatalog
    equation_blocks::Vector{ExecutableEquationBlock}
    contributions::Vector{EquationContribution}
    function ExecutableAnalysisModel(catalog, equation_blocks, contributions)
        owned = sort(vcat((block.equation_indices
                           for block in equation_blocks)...))
        owned == collect(1:length(catalog.equations)) ||
            error("Executable equation blocks must own every equation once")
        new(catalog, collect(equation_blocks), collect(contributions))
    end
end

select_variable(variable, ::PositionIC) =
    variable.kind in (:position, :orientation, :elastic_position)

select_equation(equation, ::PositionIC) =
    equation.level == 0 &&
    equation.kind in (:constraint, :normalization, :coordinate_relation)

select_variable(variable, ::VelocityIC) =
    variable.kind in (:velocity, :angular_velocity, :elastic_velocity)

select_equation(equation, ::VelocityIC) =
    equation.level == 1 &&
    equation.kind in (:constraint, :normalization, :coordinate_relation)

select_variable(variable, ::AccelerationIC) =
    variable.kind in (:acceleration, :angular_acceleration, :reaction,
                      :elastic_acceleration,
                      :applied_geometry, :applied_rate, :applied_load,
                      :user_algebraic)

select_equation(equation, ::AccelerationIC) =
    equation.kind == :balance ||
    (equation.kind == :constraint && equation.level == 2) ||
    equation.kind in (:applied_definition, :user_algebraic)

select_variable(variable, ::StaticEQ) =
    variable.kind in (:position, :orientation, :relative_position, :reaction,
                      :elastic_position,
                      :applied_geometry, :applied_load, :user_algebraic,
                      :user_state_steady)

select_equation(equation, ::StaticEQ) =
    equation.kind == :balance ||
    (equation.kind == :constraint && equation.level == 0) ||
    (equation.kind in (:coordinate_relation, :motion) && equation.level == 0) ||
    (equation.kind == :applied_definition && equation.level != 1) ||
    equation.kind in (:user_algebraic, :user_differential_steady)

select_variable(variable, ::Dynamics) = true
select_equation(equation, ::Dynamics) = true

select_variable(variable, ::KinematicPosition) =
    variable.kind in (:position, :orientation, :elastic_position,
                      :relative_position)
select_equation(equation, ::KinematicPosition) =
    equation.level == 0 &&
    equation.kind in (:constraint, :coordinate_relation, :motion)

select_variable(variable, ::KinematicVelocity) =
    variable.kind in (:velocity, :angular_velocity, :elastic_velocity,
                      :relative_velocity)
select_equation(equation, ::KinematicVelocity) =
    equation.level == 1 &&
    equation.kind in (:constraint, :coordinate_relation, :motion)

select_variable(variable, ::KinematicAcceleration) =
    variable.kind in (:acceleration, :angular_acceleration,
                      :elastic_acceleration,
                      :relative_acceleration)
select_equation(equation, ::KinematicAcceleration) =
    equation.level == 2 &&
    equation.kind in (:constraint, :coordinate_relation, :motion)

select_variable(variable, ::KinematicForces) =
    variable.kind in (:reaction, :applied_geometry, :applied_rate,
                      :applied_load, :user_algebraic)
select_equation(equation, ::KinematicForces) =
    equation.kind in (:balance, :applied_definition, :user_algebraic)

"""
    select_analysis(catalog, policy) -> AnalysisSelection

Select the ordered canonical rows and columns required by `policy`. This
metadata selection does not test squareness or numerical rank; the model
loader performs those checks after redundant families and state equations are
known.
"""
function select_analysis(catalog::AnalysisCatalog, policy::AnalysisPolicy)
    variable_indices = [variable.index for variable in catalog.variables
                        if select_variable(variable, policy)]
    equation_indices = [equation.index for equation in catalog.equations
                        if select_equation(equation, policy)]
    return AnalysisSelection(policy, variable_indices, equation_indices)
end

analysis_values(canonical, selection::AnalysisSelection) =
    canonical[selection.variable_indices]

function expand_analysis_values!(canonical, analysis,
                                 selection::AnalysisSelection)
    length(analysis) == length(selection.variable_indices) ||
        throw(DimensionMismatch("Analysis vector has the wrong length"))
    canonical[selection.variable_indices] .= analysis
    return canonical
end

function selected_block_active(indices, selected::Set{Int})
    return any(index in selected for index in indices)
end

"""Reusable canonical storage and active-component lists for residual assembly."""
mutable struct AnalysisEquationWorkspace{T}
    full_equations::Vector{T}
    selected_equations::BitVector
    active_blocks::Vector{Int}
    active_contributions::Vector{Int}
end

"""Refresh a residual workspace after the selected equation rows change."""
function configure_analysis_equation_workspace!(
        workspace::AnalysisEquationWorkspace,
        model::ExecutableAnalysisModel, selection::AnalysisSelection)
    fill!(workspace.selected_equations, false)
    workspace.selected_equations[selection.equation_indices] .= true
    empty!(workspace.active_blocks)
    for (index, block) in enumerate(model.equation_blocks)
        any(@view workspace.selected_equations[block.equation_indices]) &&
            push!(workspace.active_blocks, index)
    end
    empty!(workspace.active_contributions)
    for (index, contribution) in enumerate(model.contributions)
        any(@view workspace.selected_equations[
            contribution.target_equation_indices]) &&
            push!(workspace.active_contributions, index)
    end
    workspace
end

"""Allocate residual assembly storage for one model and analysis selection."""
function analysis_equation_workspace(model::ExecutableAnalysisModel,
        selection::AnalysisSelection, ::Type{T} = Float64) where {T}
    workspace = AnalysisEquationWorkspace(zeros(T,
            length(model.catalog.equations)),
        falses(length(model.catalog.equations)), Int[], Int[])
    configure_analysis_equation_workspace!(workspace, model, selection)
end

"""
Evaluate the selected implicit equations from their owning blocks and then add
all active component contributions. Component callbacks continue to address
canonical indices; only the returned vector is packed in selection order.
"""
function evaluate_analysis_equations!(equations, model::ExecutableAnalysisModel,
        selection::AnalysisSelection, t, canonical, canonical_derivative)
    workspace = analysis_equation_workspace(model, selection,
        eltype(canonical))
    evaluate_analysis_equations!(equations, model, selection, t, canonical,
        canonical_derivative, workspace)
end

"""Evaluate selected equations using preallocated canonical storage."""
function evaluate_analysis_equations!(equations, model::ExecutableAnalysisModel,
        selection::AnalysisSelection, t, canonical, canonical_derivative,
        workspace::AnalysisEquationWorkspace)
    length(equations) == length(selection.equation_indices) ||
        throw(DimensionMismatch("Selected equation vector has the wrong length"))
    length(workspace.full_equations) == length(model.catalog.equations) ||
        throw(DimensionMismatch("Equation workspace has the wrong length"))
    full_equations = workspace.full_equations
    fill!(full_equations, zero(eltype(full_equations)))
    for index in workspace.active_blocks
        model.equation_blocks[index].residual!(
            full_equations, t, canonical, canonical_derivative)
    end
    for index in workspace.active_contributions
        model.contributions[index].residual!(
            full_equations, t, canonical, canonical_derivative)
    end
    for (local_index, canonical_index) in
            enumerate(selection.equation_indices)
        equations[local_index] = full_equations[canonical_index]
    end
    return equations
end

"""Assemble a dense selected Newton matrix for diagnostics or small solves."""
function evaluate_analysis_jacobian(model::ExecutableAnalysisModel,
        selection::AnalysisSelection, t, canonical, canonical_derivative,
        derivative_coefficient)
    selected = Set(selection.equation_indices)
    full_jacobian = zeros(eltype(canonical),
        length(model.catalog.equations), length(model.catalog.variables))
    for block in model.equation_blocks
        selected_block_active(block.equation_indices, selected) || continue
        block.jacobian!(full_jacobian, t, canonical, canonical_derivative,
                        derivative_coefficient)
    end
    for contribution in model.contributions
        selected_block_active(contribution.target_equation_indices, selected) ||
            continue
        contribution.jacobian!(
            full_jacobian, t, canonical, canonical_derivative,
            derivative_coefficient)
    end
    return full_jacobian[
        selection.equation_indices, selection.variable_indices]
end

"""
Dictionary-backed matrix used once to discover the selected Jacobian pattern.
Explicitly assigned zeros are retained because they may be nonzero at another
mechanism configuration.
"""
mutable struct SparseJacobianAccumulator{T} <: AbstractMatrix{T}
    entries::Dict{Tuple{Int,Int},T}
    dimensions::Tuple{Int,Int}
end

SparseJacobianAccumulator{T}(rows::Integer, columns::Integer) where {T} =
    SparseJacobianAccumulator(Dict{Tuple{Int,Int},T}(), (rows, columns))

Base.size(accumulator::SparseJacobianAccumulator) = accumulator.dimensions
Base.IndexStyle(::Type{<:SparseJacobianAccumulator}) = IndexCartesian()
Base.getindex(accumulator::SparseJacobianAccumulator{T}, row::Int,
              column::Int) where {T} =
    get(accumulator.entries, (row, column), zero(T))

function Base.setindex!(accumulator::SparseJacobianAccumulator{T}, value,
                        row::Int, column::Int) where {T}
    converted = convert(T, value)
    # Retain an explicitly assigned zero. It represents a structural entry
    # whose value may become nonzero at another configuration.
    accumulator.entries[(row, column)] = converted
    return converted
end

"""Assemble selected Jacobian entries without materializing a dense matrix."""
function evaluate_analysis_sparse_jacobian(model::ExecutableAnalysisModel,
        selection::AnalysisSelection, t, canonical, canonical_derivative,
        derivative_coefficient)
    selected_equations = Set(selection.equation_indices)
    accumulator = SparseJacobianAccumulator{eltype(canonical)}(
        length(model.catalog.equations), length(model.catalog.variables))
    for block in model.equation_blocks
        selected_block_active(block.equation_indices, selected_equations) ||
            continue
        block.jacobian!(accumulator, t, canonical, canonical_derivative,
                        derivative_coefficient)
    end
    for contribution in model.contributions
        selected_block_active(contribution.target_equation_indices,
                              selected_equations) || continue
        contribution.jacobian!(accumulator, t, canonical,
            canonical_derivative, derivative_coefficient)
    end

    local_rows = Dict(index => local_index for (local_index, index) in
                      enumerate(selection.equation_indices))
    local_columns = Dict(index => local_index for (local_index, index) in
                         enumerate(selection.variable_indices))
    rows = Int[]
    columns = Int[]
    values = eltype(canonical)[]
    for ((row, column), value) in accumulator.entries
        haskey(local_rows, row) && haskey(local_columns, column) || continue
        push!(rows, local_rows[row])
        push!(columns, local_columns[column])
        push!(values, value)
    end
    return sparse(rows, columns, values,
        length(selection.equation_indices), length(selection.variable_indices))
end

"""
Matrix-like canonical-index view of an existing selected CSC matrix.

After pattern discovery, component callbacks write numerical values directly
into fixed sparse storage. Writes to canonical variables outside the current
analysis are ignored; writes outside the established selected pattern are an
error.
"""
struct FixedSparseJacobianAccumulator{T,Ti} <: AbstractMatrix{T}
    matrix::SparseMatrixCSC{T,Ti}
    canonical_to_local_row::Vector{Int}
    canonical_to_local_column::Vector{Int}
end

Base.size(accumulator::FixedSparseJacobianAccumulator) = (
    length(accumulator.canonical_to_local_row),
    length(accumulator.canonical_to_local_column))
Base.IndexStyle(::Type{<:FixedSparseJacobianAccumulator}) = IndexCartesian()

function fixed_sparse_location(accumulator::FixedSparseJacobianAccumulator,
                               row::Int, column::Int)
    local_row = accumulator.canonical_to_local_row[row]
    local_column = accumulator.canonical_to_local_column[column]
    (iszero(local_row) || iszero(local_column)) && return 0
    matrix = accumulator.matrix
    for location in matrix.colptr[local_column]:(
            matrix.colptr[local_column + 1] - 1)
        stored_row = matrix.rowval[location]
        stored_row == local_row && return location
        stored_row > local_row && break
    end
    return -1
end

function Base.getindex(accumulator::FixedSparseJacobianAccumulator{T},
                       row::Int, column::Int) where {T}
    location = fixed_sparse_location(accumulator, row, column)
    location > 0 && return accumulator.matrix.nzval[location]
    return zero(T)
end

function Base.setindex!(accumulator::FixedSparseJacobianAccumulator{T}, value,
                        row::Int, column::Int) where {T}
    converted = convert(T, value)
    location = fixed_sparse_location(accumulator, row, column)
    # A callback may touch a canonical column omitted by this analysis.
    iszero(location) && return converted
    location > 0 || error(
        "Jacobian callback wrote outside the established sparse pattern at " *
        "($row, $column)")
    accumulator.matrix.nzval[location] = converted
    return converted
end

function selection_map(canonical_size, selected_indices)
    mapping = zeros(Int, canonical_size)
    for (local_index, canonical_index) in enumerate(selected_indices)
        mapping[canonical_index] = local_index
    end
    return mapping
end

"""Accumulate into an existing selected CSC matrix with a fixed pattern."""
function evaluate_analysis_sparse_jacobian!(matrix::SparseMatrixCSC,
        model::ExecutableAnalysisModel, selection::AnalysisSelection,
        t, canonical, canonical_derivative, derivative_coefficient)
    size(matrix) == (length(selection.equation_indices),
                     length(selection.variable_indices)) ||
        throw(DimensionMismatch("Sparse Jacobian has the wrong dimensions"))
    fill!(matrix.nzval, zero(eltype(matrix)))
    accumulator = FixedSparseJacobianAccumulator(matrix,
        selection_map(length(model.catalog.equations),
                      selection.equation_indices),
        selection_map(length(model.catalog.variables),
                      selection.variable_indices))
    selected_equations = Set(selection.equation_indices)
    for block in model.equation_blocks
        selected_block_active(block.equation_indices, selected_equations) ||
            continue
        block.jacobian!(accumulator, t, canonical, canonical_derivative,
                        derivative_coefficient)
    end
    for contribution in model.contributions
        selected_block_active(contribution.target_equation_indices,
                              selected_equations) || continue
        contribution.jacobian!(accumulator, t, canonical,
            canonical_derivative, derivative_coefficient)
    end
    return matrix
end

end
