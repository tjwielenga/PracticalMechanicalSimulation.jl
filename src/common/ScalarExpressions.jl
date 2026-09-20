"""Restricted scalar expressions shared by planar and spatial model input."""
module ScalarExpressions

using ForwardDiff

export ScalarLaw, compile_time_expression, compile_motion_expression,
       compile_model_expression, constant_scalar_law,
       expression_uses_model_variables

const EXPRESSION_FUNCTIONS = Dict{Symbol,Function}(
    :sin => sin, :cos => cos, :tan => tan, :asin => asin, :acos => acos,
    :atan => atan, :tanh => tanh, :exp => exp, :log => log, :sqrt => sqrt,
    :abs => abs, :min => min, :max => max)
const EXPRESSION_OPERATORS = Dict{Symbol,Function}(
    :+ => +, :- => -, :* => *, :/ => /, :^ => ^)

"""A scalar law with its canonical dependencies and local gradient."""
struct ScalarLaw{F,G}
    value::F
    gradient::G
    dependencies::Vector{Int}
end

(law::ScalarLaw)(t, z) = law.value(t, z)
function (law::ScalarLaw)(t)
    isempty(law.dependencies) || throw(ArgumentError(
        "a state-dependent scalar law requires the current model values"))
    law.value(t, Float64[])
end

function qualified_expression_name(node)
    node isa Symbol && return node
    node isa Expr && node.head == :. && length(node.args) == 2 || return nothing
    property = node.args[2]
    property isa QuoteNode && property.value isa Symbol || return nothing
    parent = qualified_expression_name(node.args[1])
    isnothing(parent) && return nothing
    Symbol(parent, :., property.value)
end

function evaluate_scalar_expression(node, environment)
    node isa Number && return node
    if node isa Symbol
        haskey(environment, node) || throw(ArgumentError(
            "unknown expression name '$node'"))
        return environment[node]
    end
    qualified = qualified_expression_name(node)
    if !isnothing(qualified)
        haskey(environment, qualified) || throw(ArgumentError(
            "unknown expression variable '$qualified'"))
        return environment[qualified]
    end
    node isa Expr && node.head == :call || throw(ArgumentError(
        "unsupported expression syntax: $node"))
    operation = node.args[1]
    arguments = [evaluate_scalar_expression(argument, environment)
                 for argument in node.args[2:end]]
    haskey(EXPRESSION_OPERATORS, operation) &&
        return EXPRESSION_OPERATORS[operation](arguments...)
    haskey(EXPRESSION_FUNCTIONS, operation) &&
        return EXPRESSION_FUNCTIONS[operation](arguments...)
    throw(ArgumentError("function '$operation' is not allowed in expressions"))
end

function collect_expression_names!(bare_names, qualified_names, node)
    node isa Number && return nothing
    if node isa Symbol
        push!(bare_names, node)
        return nothing
    end
    qualified = qualified_expression_name(node)
    if !isnothing(qualified)
        push!(qualified_names, qualified)
        return nothing
    end
    node isa Expr && node.head == :call || throw(ArgumentError(
        "unsupported expression syntax: $node"))
    operation = node.args[1]
    (haskey(EXPRESSION_OPERATORS, operation) ||
     haskey(EXPRESSION_FUNCTIONS, operation)) || throw(ArgumentError(
        "function '$operation' is not allowed in expressions"))
    for argument in node.args[2:end]
        collect_expression_names!(bare_names, qualified_names, argument)
    end
    nothing
end

function expression_uses_model_variables(source::AbstractString)
    bare_names, qualified_names = Set{Symbol}(), Set{Symbol}()
    collect_expression_names!(bare_names, qualified_names, Meta.parse(source))
    !isempty(qualified_names)
end

"""Compile a restricted expression depending only on time and parameters."""
function compile_time_expression(source::AbstractString, parameters)
    syntax = Meta.parse(source)
    constants = Dict{Symbol,Float64}(:pi => pi, :e => MathConstants.e)
    for (name, value) in parameters
        value isa Number || throw(ArgumentError(
            "parameter '$name' must be numeric"))
        constants[Symbol(name)] = Float64(value)
    end
    expression = function (t)
        environment = Dict{Symbol,Any}(constants)
        environment[:t] = t
        value = evaluate_scalar_expression(syntax, environment)
        value isa Number || throw(ArgumentError("expression is not scalar"))
        value
    end
    expression(0.0)
    expression
end

"""
Compile a prescribed coordinate expression and obtain its first two time
derivatives with nested forward-mode differentiation.
"""
function compile_motion_expression(source::AbstractString, parameters)
    coordinate = compile_time_expression(source, parameters)
    velocity = t -> ForwardDiff.derivative(coordinate, t)
    acceleration = t -> ForwardDiff.derivative(velocity, t)
    velocity(0.0)
    acceleration(0.0)
    coordinate, velocity, acceleration
end

"""Compile a scalar law depending on supported canonical model variables."""
function compile_model_expression(source, parameters, layout,
        variable_supported)
    syntax = source isa AbstractString ? Meta.parse(source) : source
    bare_names, qualified_names = Set{Symbol}(), Set{Symbol}()
    collect_expression_names!(bare_names, qualified_names, syntax)

    constants = Dict{Symbol,Float64}(:pi => pi, :e => MathConstants.e)
    for (name, value) in parameters
        value isa Number || throw(ArgumentError(
            "parameter '$name' must be numeric"))
        constants[Symbol(name)] = Float64(value)
    end
    allowed_bare = Set([:t; collect(keys(constants))])
    unknown_bare = sort!(collect(setdiff(bare_names, allowed_bare)); by = string)
    isempty(unknown_bare) || throw(ArgumentError(
        "unknown expression name '$(first(unknown_bare))'"))

    variables = Dict(Symbol(variable.component, :., variable.name) => variable
        for variable in layout.catalog.variables)
    ordered_names = sort!(collect(qualified_names); by = string)
    dependencies = Int[]
    for name in ordered_names
        haskey(variables, name) || throw(ArgumentError(
            "unknown expression variable '$name'"))
        variable = variables[name]
        variable_supported(variable) || throw(ArgumentError(
            "expression variable '$name' has unsupported kind '$(variable.kind)'"))
        push!(dependencies, variable.index)
    end

    function evaluate_inputs(t, inputs)
        environment = Dict{Symbol,Any}(constants)
        environment[:t] = t
        for (name, value) in zip(ordered_names, inputs)
            environment[name] = value
        end
        value = evaluate_scalar_expression(syntax, environment)
        value isa Number || throw(ArgumentError("expression is not scalar"))
        value
    end
    value = (t, z) -> evaluate_inputs(t, @view(z[dependencies]))
    gradient = isempty(dependencies) ?
        ((t, z) -> Float64[]) :
        ((t, z) -> ForwardDiff.gradient(
            inputs -> evaluate_inputs(t, inputs), z[dependencies]))
    ScalarLaw(value, gradient, dependencies)
end

constant_scalar_law(value::Number) = ScalarLaw(
    (t, z) -> Float64(value), (t, z) -> Float64[], Int[])

end
