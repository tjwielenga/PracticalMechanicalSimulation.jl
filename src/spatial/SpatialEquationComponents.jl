"""User-defined, component-local scalar algebraic and differential equations."""
module SpatialEquationComponents

using ..AutomaticAnalysis
using ..InputUnits: angle_value
using ..ScalarExpressions
using ..SpatialComponentAssembly

import ..SpatialComponentAssembly: executable_blocks, equation_contributions

export SpatialEquationComponent, spatial_equation_registration,
    allocated_spatial_equation_component, initialize_spatial_equation_component!,
    spatial_equation_state_rates!, set_spatial_equation_stage!

struct SpatialEquationComponent
    name::Symbol
    state_indices::Vector{Int}
    state_modes::Vector{Symbol}
    state_equations::Vector{Int}
    state_laws::Vector{ScalarLaw}
    algebraic_indices::Vector{Int}
    algebraic_equations::Vector{Int}
    algebraic_laws::Vector{Tuple{ScalarLaw,ScalarLaw}}
    stage::Base.RefValue{Symbol}
end

function equation_component_tables(name, table)
    states = get(table, "states", Dict{String,Any}())
    variables = get(table, "variables", Dict{String,Any}())
    rates = get(table, "state_equations", Dict{String,Any}())
    equations = get(table, "equations", String[])
    inputs = get(table, "inputs", Dict{String,Any}())
    parameters = get(table, "parameters", Dict{String,Any}())
    for (label, value) in (("states", states), ("variables", variables),
            ("state_equations", rates), ("inputs", inputs),
            ("parameters", parameters))
        value isa AbstractDict || throw(ArgumentError(
            "equation component '$name'.$label must be a table"))
    end
    equations isa Vector && all(item -> item isa AbstractString, equations) ||
        throw(ArgumentError(
            "equation component '$name'.equations must be an array of strings"))
    state_names = sort!(Symbol.(collect(keys(states))); by = string)
    variable_names = sort!(Symbol.(collect(keys(variables))); by = string)
    isempty(state_names) && isempty(variable_names) && throw(ArgumentError(
        "equation component '$name' must declare a state or variable"))
    names = [state_names; variable_names]
    length(unique(names)) == length(names) || throw(ArgumentError(
        "equation component '$name' has a duplicate variable name"))
    Set(Symbol.(keys(rates))) == Set(state_names) || throw(ArgumentError(
        "equation component '$name' needs one state_equations entry per state"))
    length(equations) == length(variable_names) || throw(ArgumentError(
        "equation component '$name' needs one algebraic equation per variable"))
    reserved = Set((:t, :pi, :e, :der))
    union!(reserved, keys(ScalarExpressions.EXPRESSION_FUNCTIONS))
    aliases = [names; Symbol.(keys(inputs)); Symbol.(keys(parameters))]
    all(item -> Base.isidentifier(String(item)), aliases) ||
        throw(ArgumentError(
            "equation component '$name' names must be simple identifiers"))
    length(unique(aliases)) == length(aliases) &&
        isempty(intersect(Set(aliases), reserved)) || throw(ArgumentError(
            "equation component '$name' has a duplicate or reserved name"))
    states, variables, rates, equations, inputs, parameters,
        state_names, variable_names
end

function variable_specification(name, variable_name, value; state = false)
    value isa AbstractDict || throw(ArgumentError(
        "equation component '$name.$variable_name' must be a table"))
    initial = get(value, "initial", 0.0)
    scale = get(value, "scale", 1.0)
    initial isa Number && isfinite(initial) || throw(ArgumentError(
        "equation component '$name.$variable_name'.initial must be finite"))
    scale isa Number && isfinite(scale) && scale > 0 || throw(ArgumentError(
        "equation component '$name.$variable_name'.scale must be positive"))
    raw_mode = get(value, "static", "hold")
    raw_mode isa AbstractString || throw(ArgumentError(
        "equation component '$name.$variable_name'.static must be a string"))
    mode = Symbol(raw_mode)
    state && mode ∉ (:hold, :steady) && throw(ArgumentError(
        "equation component '$name.$variable_name'.static must be hold or steady"))
    !state && haskey(value, "static") && throw(ArgumentError(
        "equation component '$name.$variable_name'.static is only for states"))
    (; initial = Float64(initial), scale = Float64(scale), mode)
end

function spatial_equation_registration(name::Symbol, table)
    states, variables, _, equations, _, _, state_names, variable_names =
        equation_component_tables(name, table)
    declarations = VariableDeclaration[]
    state_rows = EquationDeclaration[]
    for variable in state_names
        specification = variable_specification(name, variable,
            states[String(variable)]; state = true)
        kind = specification.mode == :steady ?
            :user_state_steady : :user_state_hold
        push!(declarations, VariableDeclaration(variable, kind, 0,
            specification.scale))
        push!(state_rows, EquationDeclaration(Symbol(:rate_, variable),
            specification.mode == :steady ?
                :user_differential_steady : :user_differential_hold,
            0, :user_equation))
    end
    algebraic_rows = EquationDeclaration[]
    for (number, variable) in enumerate(variable_names)
        specification = variable_specification(name, variable,
            variables[String(variable)])
        push!(declarations, VariableDeclaration(variable,
            :user_algebraic, 2, specification.scale))
        push!(algebraic_rows, EquationDeclaration(Symbol(:equation_, number),
            :user_algebraic, 2, :user_equation))
    end
    blocks = EquationBlockDeclaration[]
    isempty(state_rows) || push!(blocks,
        EquationBlockDeclaration(:differential, state_rows))
    isempty(equations) || push!(blocks,
        EquationBlockDeclaration(:algebraic, algebraic_rows))
    ComponentRegistration(name, declarations, blocks)
end

"""Replace local aliases by qualified canonical names before safe compilation."""
function expand_equation_aliases(node, aliases)
    if node isa Symbol
        return get(aliases, node, node)
    elseif node isa Expr && node.head == :call
        return Expr(:call, node.args[1],
            (expand_equation_aliases(argument, aliases)
                for argument in node.args[2:end])...)
    end
    node
end

function compile_component_law(source, parameters, layout, aliases)
    source isa AbstractString || throw(ArgumentError(
        "user equation expression must be a string"))
    syntax = expand_equation_aliases(Meta.parse(source), aliases)
    ScalarExpressions.compile_model_expression(syntax, parameters, layout,
        variable -> true)
end

function compile_component_equality(source, parameters, layout, aliases)
    parts = split(source, '='; limit = 2)
    length(parts) == 2 && !isempty(strip(parts[1])) &&
        !isempty(strip(parts[2])) || throw(ArgumentError(
            "user algebraic equation must have the form 'left = right'"))
    left = compile_component_law(strip(parts[1]), parameters, layout,
        aliases)
    right = compile_component_law(strip(parts[2]), parameters, layout,
        aliases)
    (left, right)
end

function compile_state_rate(source, variable, parameters, layout, aliases)
    source isa AbstractString || throw(ArgumentError(
        "state equation for '$variable' must be a string"))
    parts = split(source, '='; limit = 2)
    if length(parts) == 2
        left = Meta.parse(strip(parts[1]))
        left isa Expr && left.head == :call && length(left.args) == 2 &&
            left.args[1] == :der && left.args[2] == variable ||
            throw(ArgumentError(
                "state equation for '$variable' must begin with der($variable)"))
        return compile_component_law(strip(parts[2]), parameters, layout,
            aliases)
    end
    compile_component_law(source, parameters, layout, aliases)
end

function allocated_spatial_equation_component(layout, name, table,
        global_parameters)
    states, variables, rates, equations, inputs, local_parameters,
        state_names, variable_names = equation_component_tables(name, table)
    parameters = Dict{String,Float64}(string(key) => Float64(value)
        for (key, value) in global_parameters)
    for (key, value) in local_parameters
        converted = if value isa Number
            Float64(value)
        elseif value isa AbstractString
            angle_value(value, "equation component '$name' parameter '$key'")
        else
            throw(ArgumentError(
                "equation component '$name' parameter '$key' must be " *
                "numeric or a quoted degree angle"))
        end
        isfinite(converted) || throw(ArgumentError(
            "equation component '$name' parameter '$key' must be finite"))
        haskey(parameters, String(key)) && throw(ArgumentError(
            "equation component '$name' parameter '$key' duplicates a model parameter"))
        parameters[String(key)] = converted
    end
    for alias in [state_names; variable_names; Symbol.(keys(inputs))]
        haskey(global_parameters, String(alias)) && throw(ArgumentError(
            "equation component '$name' name '$alias' conflicts with a model parameter"))
    end
    aliases = Dict{Symbol,Any}()
    for variable in [state_names; variable_names]
        aliases[variable] = Meta.parse(string(name, ".", variable))
    end
    for (alias, target) in inputs
        target isa AbstractString || throw(ArgumentError(
            "equation component '$name' input '$alias' must name a model variable"))
        syntax = Meta.parse(String(target))
        qualified = ScalarExpressions.qualified_expression_name(syntax)
        qualified !== nothing && occursin('.', String(target)) &&
            any(variable -> Symbol(variable.component, :., variable.name) ==
                qualified, layout.catalog.variables) || throw(ArgumentError(
            "equation component '$name' input '$alias' names unknown " *
            "model variable '$target'"))
        aliases[Symbol(alias)] = syntax
    end
    indices = component_variable_indices(layout, name)
    nstates = length(state_names)
    state_indices = collect(first(indices):(first(indices) + nstates - 1))
    algebraic_indices = collect((first(indices) + nstates):last(indices))
    state_modes = [variable_specification(name, variable,
        states[String(variable)]; state = true).mode for variable in state_names]
    state_laws = [compile_state_rate(rates[String(variable)], variable,
        parameters, layout, aliases) for variable in state_names]
    algebraic_laws = Tuple{ScalarLaw,ScalarLaw}[
        compile_component_equality(source, parameters, layout, aliases)
        for source in equations]
    state_rows = isempty(state_names) ? Int[] : collect(
        component_equation_indices(layout, name, :differential))
    algebraic_rows = isempty(variable_names) ? Int[] : collect(
        component_equation_indices(layout, name, :algebraic))
    SpatialEquationComponent(name, state_indices, state_modes, state_rows,
        state_laws, algebraic_indices, algebraic_rows, algebraic_laws,
        Ref(:dynamic))
end

function initialize_spatial_equation_component!(initial, component, table)
    states, variables, _, _, _, _, state_names, variable_names =
        equation_component_tables(component.name, table)
    for (index, name) in zip(component.state_indices, state_names)
        initial[index] = variable_specification(component.name, name,
            states[String(name)]; state = true).initial
    end
    for (index, name) in zip(component.algebraic_indices, variable_names)
        initial[index] = variable_specification(component.name, name,
            variables[String(name)]).initial
    end
    initial
end

function set_spatial_equation_stage!(component, stage)
    component.stage[] = stage
    component
end

function spatial_equation_state_rates!(derivative, component, t, z)
    for (index, mode, law) in zip(component.state_indices,
            component.state_modes, component.state_laws)
        derivative[index] = component.stage[] == :static && mode == :hold ?
            0.0 : law(t, z)
    end
    derivative
end

function executable_blocks(component::SpatialEquationComponent)
    blocks = ExecutableEquationBlock[]
    if !isempty(component.state_indices)
        residual! = function (equations, t, z, zdot)
            for (row, index, mode, law) in zip(component.state_equations,
                    component.state_indices, component.state_modes,
                    component.state_laws)
                held = component.stage[] == :static && mode == :hold
                equations[row] = zdot[index] - (held ? 0.0 : law(t, z))
            end
        end
        jacobian! = function (jacobian, t, z, zdot, coefficient)
            for (row, index, mode, law) in zip(component.state_equations,
                    component.state_indices, component.state_modes,
                    component.state_laws)
                jacobian[row, index] += coefficient
                component.stage[] == :static && mode == :hold && continue
                for (column, partial) in zip(law.dependencies,
                        law.gradient(t, z))
                    jacobian[row, column] -= partial
                end
            end
        end
        push!(blocks, ExecutableEquationBlock(component.name, :differential,
            component.state_equations, residual!, jacobian!))
    end
    if !isempty(component.algebraic_indices)
        residual! = function (equations, t, z, zdot)
            for (row, (left, right)) in zip(component.algebraic_equations,
                    component.algebraic_laws)
                equations[row] = left(t, z) - right(t, z)
            end
        end
        jacobian! = function (jacobian, t, z, zdot, coefficient)
            for (row, (left, right)) in zip(component.algebraic_equations,
                    component.algebraic_laws)
                for (column, partial) in zip(left.dependencies,
                        left.gradient(t, z))
                    jacobian[row, column] += partial
                end
                for (column, partial) in zip(right.dependencies,
                        right.gradient(t, z))
                    jacobian[row, column] -= partial
                end
            end
        end
        push!(blocks, ExecutableEquationBlock(component.name, :algebraic,
            component.algebraic_equations, residual!, jacobian!))
    end
    blocks
end

equation_contributions(::SpatialEquationComponent) = EquationContribution[]

end
