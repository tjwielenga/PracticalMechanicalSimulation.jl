"""Build-time expansion of reusable TOML modeling assemblies."""
module AssemblyExpansion

using LinearAlgebra
using LuaCall
using TOML
using ..InputUnits: angle_value

export expand_model_assemblies, parse_lua_model

const RESERVED_ASSEMBLY_TABLES = Set((
    "assembly", "interface", "parameters", "calculate"))
const REFERENCE_FIELDS = Set((
    "marker", "markers", "body", "bodies", "joint", "joints",
    "center_of_mass", "contact_marker", "span", "spans", "pulley",
    "pulleys"))
const DEFAULT_ASSEMBLY_DIRECTORY = normpath(joinpath(
    @__DIR__, "..", "..", "assemblies"))

struct AssemblyDefinition
    name::String
    dimension::String
    document::Dict{String,Any}
    path::String
end

struct BodyReference
    name::String
end

const LUA_ASSEMBLY_PRELUDE = raw"""
local vector_mt = {}

function vector(values)
    return setmetatable(values, vector_mt)
end

local function vector_result(a, b, operation)
    if #a ~= #b then
        error("assembly vectors must have the same length")
    end
    local result = {}
    for i = 1, #a do
        result[i] = operation(a[i], b[i])
    end
    return vector(result)
end

function vector_mt.__add(a, b)
    return vector_result(a, b, function(x, y) return x + y end)
end

function vector_mt.__sub(a, b)
    return vector_result(a, b, function(x, y) return x - y end)
end

function vector_mt.__unm(a)
    local result = {}
    for i = 1, #a do
        result[i] = -a[i]
    end
    return vector(result)
end

function vector_mt.__mul(a, b)
    local result = {}
    if type(a) == "number" then
        for i = 1, #b do result[i] = a*b[i] end
        return vector(result)
    elseif type(b) == "number" then
        for i = 1, #a do result[i] = a[i]*b end
        return vector(result)
    end
    error("assembly vectors may only be multiplied by scalars")
end

function vector_mt.__div(a, b)
    if type(b) ~= "number" then
        error("assembly vectors may only be divided by scalars")
    end
    local result = {}
    for i = 1, #a do result[i] = a[i]/b end
    return vector(result)
end

function dot(a, b)
    if #a ~= #b then
        error("assembly vectors must have the same length")
    end
    local result = 0
    for i = 1, #a do result = result+a[i]*b[i] end
    return result
end

function cross(a, b)
    if #a ~= 3 or #b ~= 3 then
        error("cross products require three-vectors")
    end
    return vector({
        a[2]*b[3] - a[3]*b[2],
        a[3]*b[1] - a[1]*b[3],
        a[1]*b[2] - a[2]*b[1]
    })
end

function norm(a)
    return math.sqrt(dot(a, a))
end

function unit(a)
    local length = norm(a)
    if length <= 1.0e-14 then
        error("cannot form a unit vector from a zero-length vector")
    end
    return a / length
end

abs = math.abs
sqrt = math.sqrt
sin = math.sin
cos = math.cos
tan = math.tan
asin = math.asin
acos = math.acos
atan = math.atan
min = math.min
max = math.max

function frame(x, y, z)
    return {
        {x[1], y[1], z[1]},
        {x[2], y[2], z[2]},
        {x[3], y[3], z[3]}
    }
end

local function transpose(a)
    return {
        {a[1][1], a[2][1], a[3][1]},
        {a[1][2], a[2][2], a[3][2]},
        {a[1][3], a[2][3], a[3][3]}
    }
end

local function matrix_vector(a, x)
    return vector({dot(vector(a[1]), x), dot(vector(a[2]), x),
        dot(vector(a[3]), x)})
end

function compose(a, b)
    local bt = transpose(b)
    local result = {{}, {}, {}}
    for i = 1, 3 do
        for j = 1, 3 do
            result[i][j] = dot(vector(a[i]), vector(bt[j]))
        end
    end
    return result
end

function rotation(angle, axis)
    if type(angle) == "string" then
        local number, unit_name = string.match(angle,
            "^%s*([%+%-]?[%d%.]+)%s*([%a°]*)%s*$")
        if number == nil or (unit_name ~= "deg" and unit_name ~= "°") then
            error("quoted assembly angles must end in 'deg' or '°'")
        end
        angle = tonumber(number)*math.pi/180
    end
    local a = unit(vector(axis))
    local x, y, z = a[1], a[2], a[3]
    local c, s, d = math.cos(angle), math.sin(angle), 1-math.cos(angle)
    return {
        {c+x*x*d, x*y*d-z*s, x*z*d+y*s},
        {y*x*d+z*s, c+y*y*d, y*z*d-x*s},
        {z*x*d-y*s, z*y*d+x*s, c+z*z*d}
    }
end

function link_frame(point_a, point_b, side)
    local x = unit(point_b-point_a)
    local z = unit(side-dot(side, x)*x)
    return frame(x, cross(z, x), z)
end

function box_inertia(mass, size)
    local x, y, z = size[1], size[2], size[3]
    return vector({
        mass*(y*y+z*z)/12,
        mass*(x*x+z*z)/12,
        mass*(x*x+y*y)/12
    })
end

function local_point(body, point)
    if type(body) == "string" then
        return {
            __simp_deferred = "local_point",
            body = body,
            value = point
        }
    end
    return matrix_vector(transpose(body.orientation), point-body.position)
end

function local_frame(body, orientation)
    if type(body) == "string" then
        return {
            __simp_deferred = "local_frame",
            body = body,
            value = orientation
        }
    end
    return compose(transpose(body.orientation), orientation)
end

local function prepare_parameters(definition, p)
    for name, specification in pairs(definition.parameters or {}) do
        if specification.kind == "vector" and p[name] ~= nil then
            p[name] = vector(p[name])
        end
    end
    for name, _ in pairs(definition.interfaces or {}) do
        local body = p[name]
        body.position = vector(body.position)
    end
    return p
end

-- Assembly files are trusted model code, but they do not need file, process,
-- package, or Julia access while constructing their result table.
julia = nil
os = nil
io = nil
package = nil
require = nil
dofile = nil
loadfile = nil
load = nil
debug = nil
"""

const LUA_MODEL_RECORDER = raw"""
local simp_model_entries = {}

local function record_section(name)
    return function(values)
        if type(values) ~= "table" then
            error(name .. " must be called with one table")
        end
        table.insert(simp_model_entries, {
            section = name,
            values = values
        })
    end
end

model = record_section("model")
analysis = record_section("analysis")
simulation = record_section("simulation")
parameters = record_section("parameters")
graphics = record_section("graphics")
state_selection = record_section("state_selection")
initial_conditions = record_section("initial_conditions")

local function register_element_type(element_type)
    rawset(_ENV, element_type, function(values)
        if type(values) ~= "table" then
            error(element_type .. " must be called with one table")
        end
        if values.name == nil then
            error(element_type .. " requires name")
        end
        table.insert(simp_model_entries, {
            element_type = element_type,
            values = values
        })
    end)
end

for _, element_type in ipairs(simp_builtin_element_types) do
    register_element_type(element_type)
end

function model_table(values)
    if type(values) ~= "table" then
        error("model_table must be called with one table")
    end
    if type(values.name) ~= "string" or values.name == "" then
        error("model_table requires name")
    end
    table.insert(simp_model_entries, {
        table_path = values.name,
        values = values
    })
end

require = simp_safe_require

setmetatable(_ENV, {
    __index = function(_, name)
        error("unknown model statement '" .. tostring(name) .. "'")
    end
})
"""

path_text(path) = join(path, ".")

function table_at(document, path; create = false)
    table = document
    for key in path
        if create
            child = get!(table, key, Dict{String,Any}())
        else
            haskey(table, key) || return nothing
            child = table[key]
        end
        child isa AbstractDict || throw(ArgumentError(
            "TOML path '$(path_text(path))' conflicts with a value"))
        table = child
    end
    table
end

function remove_table!(document, path)
    isempty(path) && return document
    parent = table_at(document, path[1:end-1])
    isnothing(parent) || delete!(parent, path[end])
    document
end

function common_prefix_length(first, second)
    count = 0
    for (a, b) in zip(first, second)
        a == b || break
        count += 1
    end
    count
end

function imported_assembly_paths(table, label)
    value = get(table, "assemblies", String[])
    value isa Vector && all(item -> item isa AbstractString, value) ||
        throw(ArgumentError("$label.assemblies must be an array of paths"))
    String.(value)
end

function lua_stack_value(state, index)
    type_index = LuaCall.LibLua.lua_type(state, index)
    if type_index == LuaCall.LibLua.LUA_TBOOLEAN
        return !iszero(LuaCall.LibLua.lua_toboolean(state, index))
    end
    LuaCall.getstack(state, index, type_index)
end

function lua_to_julia(value)
    value isa LuaCall.LuaTable || return value
    entries = Pair{Any,Any}[]
    state = LuaCall.LS(value)
    LuaCall.checkstack(state, 2)
    push!(state, nothing)
    while !iszero(LuaCall.LibLua.lua_next(state, LuaCall.idx(value)))
        key = lua_stack_value(state, -2)
        item = lua_stack_value(state, -1)
        push!(entries, lua_to_julia(key) => lua_to_julia(item))
        pop!(state, 1)
    end
    integer_keys = [key for (key, _) in entries if key isa Integer]
    if length(integer_keys) == length(entries) &&
            sort!(Int.(integer_keys)) == collect(1:length(entries))
        result = Vector{Any}(undef, length(entries))
        for (key, item) in entries
            result[Int(key)] = item
        end
        return result
    end
    all(pair -> first(pair) isa AbstractString, entries) ||
        throw(ArgumentError("Lua assembly tables must use string keys or consecutive integer keys"))
    Dict{String,Any}(String(key) => item for (key, item) in entries)
end

function lua_assembly_document(path)
    source = read(path, String)
    state = LuaCall.LuaState(nothing)
    try
        LuaCall.@luascope state begin
            chunk = LuaCall.lualoadstring(state,
                LUA_ASSEMBLY_PRELUDE * "\nlocal definition = (function()\n" *
                source * "\nend)()\nreturn definition")
            definition = chunk()
            definition isa LuaCall.LuaTable || throw(ArgumentError(
                "Lua assembly '$path' must return a table"))
            name = definition["name"]
            dimension = definition["dimension"]
            interfaces = definition["interfaces"]
            parameters = definition["parameters"]
            imports = definition["assemblies"]
            name isa AbstractString || throw(ArgumentError(
                "Lua assembly '$path' requires a string name"))
            dimension isa AbstractString || throw(ArgumentError(
                "Lua assembly '$path' requires a string dimension"))
            interface_document = isnothing(interfaces) ? Dict{String,Any}() :
                lua_to_julia(interfaces)
            parameter_document = isnothing(parameters) ? Dict{String,Any}() :
                lua_to_julia(parameters)
            interface_document isa AbstractVector && isempty(interface_document) &&
                (interface_document = Dict{String,Any}())
            parameter_document isa AbstractVector && isempty(parameter_document) &&
                (parameter_document = Dict{String,Any}())
            import_document = isnothing(imports) ? String[] :
                lua_to_julia(imports)
            Dict{String,Any}(
                "assembly" => Dict{String,Any}(
                    "name" => String(name),
                    "dimension" => String(dimension),
                    "assemblies" => import_document),
                "interface" => interface_document,
                "parameters" => parameter_document)
        end
    catch error
        error isa ArgumentError && rethrow()
        throw(ArgumentError(
            "could not read Lua assembly '$path': $(sprint(showerror, error))"))
    end
end

function merge_lua_model_table!(document, path, values, label)
    output = table_at(document, path; create = true)
    for (key, value) in values
        haskey(output, key) && throw(ArgumentError(
            "$label defines '$(path_text(path)).$key' more than once"))
        output[String(key)] = value
    end
    document
end


function lua_string_literal(value)
    '"' * replace(String(value), '\\' => "\\\\", '"' => "\\\"") * '"'
end

function lua_require_setup(source_directory)
    directories = isnothing(source_directory) ? String[] :
        [abspath(String(source_directory))]
    push!(directories, DEFAULT_ASSEMBLY_DIRECTORY)
    paths = String[]
    for directory in unique(directories)
        push!(paths, joinpath(directory, "?.lua"))
        push!(paths, joinpath(directory, "?", "init.lua"))
    end
    "local simp_safe_require = require\n" *
        "package.path = " * lua_string_literal(join(paths, ";")) * "\n" *
        "package.cpath = \"\"\n" *
        "package.searchers[3] = nil\n" *
        "package.searchers[4] = nil\n"
end

function resolve_lua_deferred_value(value, document)
    if value isa AbstractDict
        operation = get(value, "__simp_deferred", nothing)
        if operation isa AbstractString
            body = get(value, "body", nothing)
            body isa AbstractString || throw(ArgumentError(
                "deferred Lua geometry requires a body name"))
            input = get(value, "value", nothing)
            if operation == "planar_local_point"
                position, orientation = planar_body_pose(
                    document, BodyReference(body))
                point = numeric_planar_vector(
                    input, "planar_local_point value")
                c, s = cos(orientation), sin(orientation)
                offset = point - position
                return [c*offset[1] + s*offset[2],
                    -s*offset[1] + c*offset[2]]
            elseif operation == "planar_local_orientation"
                _, orientation = planar_body_pose(
                    document, BodyReference(body))
                return angle_value(input,
                    "planar_local_orientation value") - orientation
            end
            position, orientation = body_pose(document, BodyReference(body))
            if operation == "local_point"
                point = numeric_vector(input, "local_point value")
                return transpose(orientation) * (point - position)
            elseif operation == "local_frame"
                frame_value = numeric_matrix(input, "local_frame value")
                local_value = transpose(orientation) * frame_value
                return [collect(row) for row in eachrow(local_value)]
            end
            throw(ArgumentError(
                "unknown deferred Lua geometry operation '$operation'"))
        end
        for key in collect(keys(value))
            value[key] = resolve_lua_deferred_value(value[key], document)
        end
    elseif value isa AbstractVector
        for index in eachindex(value)
            value[index] = resolve_lua_deferred_value(value[index], document)
        end
    end
    value
end

"""Convert a Lua data-description model into the ordinary model document."""
function parse_lua_model(source::AbstractString; label = "Lua model",
        source_directory = nothing, element_types = String[])
    builtin_types = "local simp_builtin_element_types = {" *
        join(lua_string_literal.(sort!(String.(collect(element_types)))), ",") * "}\n"
    program = lua_require_setup(source_directory) * LUA_ASSEMBLY_PRELUDE *
        "\n" * builtin_types * LUA_MODEL_RECORDER *
        "\nlocal simp_source_result = (function()\n" * source *
        "\nend)()\n" *
        "if #simp_model_entries == 0 and simp_source_result ~= nil then\n" *
        "    error(\"this Lua file returns an assembly module and does not " *
        "define a model; open a model file that calls the assembly\")\n" *
        "end\n" *
        "return simp_model_entries"
    state = LuaCall.LuaState(nothing)
    try
        LuaCall.@luascope state begin
            chunk = LuaCall.lualoadstring(state, program)
            entries = chunk()
            entries isa LuaCall.LuaTable || throw(ArgumentError(
                "$label did not produce a model entry list"))
            converted = lua_to_julia(entries)
            converted isa AbstractVector || throw(ArgumentError(
                "$label model entries are not an array"))
            document = Dict{String,Any}()
            for entry in converted
                entry isa AbstractDict || throw(ArgumentError(
                    "$label contains an invalid model entry"))
                values = get(entry, "values", nothing)
                values isa AbstractDict || throw(ArgumentError(
                    "$label entries must contain named fields"))
                fields = Dict{String,Any}(values)
                if haskey(entry, "table_path")
                    name = pop!(fields, "name")
                    merge_lua_model_table!(document,
                        split(String(name), '.'), fields, label)
                    continue
                end
                if haskey(entry, "section")
                    section = String(entry["section"])
                    merge_lua_model_table!(document, [section], fields, label)
                    continue
                end
                kind = String(get(entry, "element_type", ""))
                haskey(fields, "name") || throw(ArgumentError(
                    "$label $kind entry requires name"))
                name = fields["name"]
                name isa AbstractString || throw(ArgumentError(
                    "$label $kind entry name must be a string"))
                delete!(fields, "name")
                haskey(fields, "type") && throw(ArgumentError(
                    "$label entry '$name' must not specify type; its function name is the type"))
                fields["type"] = kind
                merge_lua_model_table!(document, split(String(name), '.'),
                    fields, label)
            end
            resolve_lua_deferred_value(document, document)
            document
        end
    catch error
        error isa ArgumentError && rethrow()
        throw(ArgumentError(
            "could not read $label: $(sprint(showerror, error))"))
    end
end

function load_assembly_definition!(registry, path, dimension, loading)
    full_path = abspath(path)
    full_path in loading && throw(ArgumentError(
        "cyclic assembly-file import involving '$full_path'"))
    isfile(full_path) || throw(ArgumentError(
        "assembly definition '$full_path' does not exist"))
    push!(loading, full_path)
    document = endswith(lowercase(full_path), ".lua") ?
        lua_assembly_document(full_path) : TOML.parsefile(full_path)
    header = get(document, "assembly", nothing)
    header isa AbstractDict || throw(ArgumentError(
        "assembly definition '$full_path' requires an [assembly] table"))
    haskey(header, "name") || throw(ArgumentError(
        "assembly definition '$full_path' requires assembly.name"))
    name = String(header["name"])
    assembly_dimension = String(get(header, "dimension", dimension))
    assembly_dimension == dimension || throw(ArgumentError(
        "assembly '$name' has dimension '$assembly_dimension', expected '$dimension'"))
    definition = AssemblyDefinition(name, assembly_dimension, document, full_path)
    if haskey(registry, name)
        registry[name].path == full_path || throw(ArgumentError(
            "assembly type '$name' is defined by more than one file"))
    else
        registry[name] = definition
    end
    for imported in imported_assembly_paths(header, "assembly '$name'")
        load_assembly_definition!(registry,
            joinpath(dirname(full_path), imported), dimension, loading)
    end
    delete!(loading, full_path)
    definition
end

function load_assembly_registry(document, source_directory, dimension)
    model = get(document, "model", Dict{String,Any}())
    registry = Dict{String,AssemblyDefinition}()
    for imported in imported_assembly_paths(model, "model")
        load_assembly_definition!(registry,
            joinpath(source_directory, imported), dimension, Set{String}())
    end
    registry
end

function assembly_instance_paths!(paths, table, registry, path = String[])
    for (key, value) in table
        value isa AbstractDict || continue
        child_path = [path; key]
        kind = get(value, "type", nothing)
        if kind isa AbstractString && haskey(registry, String(kind))
            push!(paths, child_path)
        else
            assembly_instance_paths!(paths, value, registry, child_path)
        end
    end
    paths
end

function finite_build_value(value, label)
    if value isa Number
        isfinite(value) || throw(ArgumentError("$label must be finite"))
    elseif value isa AbstractArray
        all(item -> item isa Number && isfinite(item), value) ||
            throw(ArgumentError("$label must contain finite numbers"))
    end
    value
end

function numeric_vector(value, label)
    value isa AbstractVector && length(value) == 3 &&
        all(item -> item isa Number && isfinite(item), value) ||
        throw(ArgumentError("$label must be a finite three-vector"))
    Float64.(value)
end

function numeric_planar_vector(value, label)
    value isa AbstractVector && length(value) == 2 &&
        all(item -> item isa Number && isfinite(item), value) ||
        throw(ArgumentError("$label must be a finite two-vector"))
    Float64.(value)
end

function numeric_matrix(value, label)
    if value isa AbstractMatrix && size(value) == (3, 3)
        result = Float64.(value)
    elseif value isa AbstractVector && length(value) == 3 &&
            all(row -> row isa AbstractVector && length(row) == 3 &&
                all(item -> item isa Number, row), value)
        result = reduce(vcat,
            permutedims.([Float64.(row) for row in value]))
    else
        throw(ArgumentError("$label must be a 3 by 3 matrix"))
    end
    all(isfinite, result) || throw(ArgumentError("$label must be finite"))
    result
end

function axis_angle_rotation(angle, axis)
    direction = numeric_vector(axis, "rotation axis")
    magnitude = norm(direction)
    magnitude > sqrt(eps(Float64)) || throw(ArgumentError(
        "rotation axis must have nonzero length"))
    x, y, z = direction / magnitude
    theta = angle_value(angle, "assembly rotation angle")
    c, s, d = cos(theta), sin(theta), 1 - cos(theta)
    [c + x*x*d x*y*d-z*s x*z*d+y*s;
     y*x*d+z*s c+y*y*d y*z*d-x*s;
     z*x*d-y*s z*y*d+x*s c+z*z*d]
end

function body_pose(document, reference::BodyReference)
    path = split(reference.name, '.')
    table = table_at(document, path)
    table isa AbstractDict || throw(ArgumentError(
        "assembly interface references unknown body or ground '$(reference.name)'"))
    kind = get(table, "type", nothing)
    kind in ("rigid_body", "ground") || throw(ArgumentError(
        "assembly interface '$(reference.name)' must name a rigid body or ground"))
    position = kind == "ground" ? zeros(3) :
        numeric_vector(get(table, "position", zeros(3)),
            "body '$(reference.name)'.position")
    orientation_value = get(table, "orientation", nothing)
    orientation = if isnothing(orientation_value)
        Matrix{Float64}(I, 3, 3)
    elseif orientation_value isa AbstractVector && length(orientation_value) == 4
        axis_angle_rotation(orientation_value[1], orientation_value[2:4])
    else
        numeric_matrix(orientation_value, "body '$(reference.name)'.orientation")
    end
    position, orientation
end

function planar_body_pose(document, reference::BodyReference)
    path = split(reference.name, '.')
    table = table_at(document, path)
    table isa AbstractDict || throw(ArgumentError(
        "assembly interface references unknown body or ground '$(reference.name)'"))
    kind = get(table, "type", nothing)
    kind in ("rigid_body", "ground") || throw(ArgumentError(
        "assembly interface '$(reference.name)' must name a rigid body or ground"))
    position = kind == "ground" ? zeros(2) :
        numeric_planar_vector(get(table, "position", zeros(2)),
            "body '$(reference.name)'.position")
    orientation = kind == "ground" ? 0.0 : angle_value(
        get(table, "orientation", 0.0),
        "body '$(reference.name)'.orientation")
    position, orientation
end

function lua_quoted_string(value)
    escaped = replace(String(value),
        '\\' => "\\\\", '"' => "\\\"", '\n' => "\\n",
        '\r' => "\\r", '\t' => "\\t")
    "\"" * escaped * "\""
end

function lua_literal(value, document)
    value isa Bool && return value ? "true" : "false"
    value isa Integer && return string(value)
    value isa AbstractFloat && begin
        isfinite(value) || throw(ArgumentError(
            "Lua assembly inputs must be finite"))
        return repr(Float64(value))
    end
    value isa Real && begin
        converted = Float64(value)
        isfinite(converted) || throw(ArgumentError(
            "Lua assembly inputs must be finite"))
        return repr(converted)
    end
    value isa AbstractString && return lua_quoted_string(value)
    if value isa BodyReference
        position, orientation = body_pose(document, value)
        return "{name=" * lua_quoted_string(value.name) *
            ",position=" * lua_literal(position, document) *
            ",orientation=" * lua_literal(orientation, document) * "}"
    elseif value isa AbstractMatrix
        rows = [collect(row) for row in eachrow(value)]
        return lua_literal(rows, document)
    elseif value isa AbstractVector
        return "{" * join((lua_literal(item, document)
            for item in value), ",") * "}"
    elseif value isa AbstractDict
        fields = sort!(collect(keys(value)); by = string)
        return "{" * join(("[" * lua_quoted_string(string(key)) * "]=" *
            lua_literal(value[key], document) for key in fields), ",") * "}"
    end
    throw(ArgumentError(
        "unsupported Lua assembly input value of type $(typeof(value))"))
end

function lua_assembly_template(definition, environment, document)
    source = read(definition.path, String)
    parameters = Dict(String(name) => value for (name, value) in environment)
    input = lua_literal(parameters, document)
    program = LUA_ASSEMBLY_PRELUDE *
        "\nlocal definition = (function()\n" * source *
        "\nend)()\nlocal parameters = " * input *
        "\nprepare_parameters(definition, parameters)" *
        "\nreturn definition.build(parameters)"
    state = LuaCall.LuaState(nothing)
    try
        LuaCall.@luascope state begin
            chunk = LuaCall.lualoadstring(state, program)
            template = chunk()
            template isa LuaCall.LuaTable || throw(ArgumentError(
                "Lua assembly '$(definition.name)' build function must return a table"))
            converted = lua_to_julia(template)
            converted isa AbstractVector && isempty(converted) &&
                (converted = Dict{String,Any}())
            converted isa AbstractDict || throw(ArgumentError(
                "Lua assembly '$(definition.name)' build result must be a table with named fields"))
            Dict{String,Any}(converted)
        end
    catch error
        error isa ArgumentError && rethrow()
        throw(ArgumentError(
            "could not build Lua assembly '$(definition.name)': " *
            sprint(showerror, error)))
    end
end

function build_function(name, arguments, document)
    name == :dot && return dot(arguments...)
    name == :cross && return cross(arguments...)
    name == :norm && return norm(arguments...)
    if name == :unit
        vector = numeric_vector(only(arguments), "unit argument")
        magnitude = norm(vector)
        magnitude > sqrt(eps(Float64)) || throw(ArgumentError(
            "unit argument must have nonzero length"))
        return vector / magnitude
    elseif name == :frame
        length(arguments) == 3 || throw(ArgumentError(
            "frame requires x, y, and z vectors"))
        return hcat((numeric_vector(argument, "frame axis")
            for argument in arguments)...)
    elseif name == :link_frame
        length(arguments) == 3 || throw(ArgumentError(
            "link_frame requires two points and a side direction"))
        point_a = numeric_vector(arguments[1], "link_frame first point")
        point_b = numeric_vector(arguments[2], "link_frame second point")
        x_axis = point_b - point_a
        x_axis /= norm(x_axis)
        side = numeric_vector(arguments[3], "link_frame side direction")
        z_axis = side - dot(side, x_axis) * x_axis
        norm(z_axis) > sqrt(eps(Float64)) || throw(ArgumentError(
            "link_frame side direction must not be parallel to its span"))
        z_axis /= norm(z_axis)
        return hcat(x_axis, cross(z_axis, x_axis), z_axis)
    elseif name == :rotation
        length(arguments) == 2 || throw(ArgumentError(
            "rotation requires an angle and an axis"))
        return axis_angle_rotation(arguments...)
    elseif name == :compose
        return *(arguments...)
    elseif name == :box_inertia
        length(arguments) == 2 || throw(ArgumentError(
            "box_inertia requires mass and size"))
        mass = arguments[1]
        size = numeric_vector(arguments[2], "box_inertia size")
        mass isa Number && mass >= 0 || throw(ArgumentError(
            "box_inertia mass must be nonnegative"))
        x, y, z = size
        return Float64(mass) .* [y^2 + z^2, x^2 + z^2, x^2 + y^2] ./ 12
    elseif name == :local_point
        length(arguments) == 2 && arguments[1] isa BodyReference ||
            throw(ArgumentError("local_point requires a body interface and point"))
        position, orientation = body_pose(document, arguments[1])
        return transpose(orientation) *
            (numeric_vector(arguments[2], "local_point point") - position)
    elseif name == :local_frame
        length(arguments) == 2 && arguments[1] isa BodyReference ||
            throw(ArgumentError("local_frame requires a body interface and frame"))
        _, orientation = body_pose(document, arguments[1])
        return transpose(orientation) *
            numeric_matrix(arguments[2], "local_frame frame")
    end
    scalar_functions = Dict{Symbol,Function}(
        :abs => abs, :sqrt => sqrt, :sin => sin, :cos => cos, :tan => tan,
        :asin => asin, :acos => acos, :atan => atan, :min => min, :max => max)
    haskey(scalar_functions, name) && return scalar_functions[name](arguments...)
    throw(ArgumentError("function '$name' is not allowed in assembly expressions"))
end

function evaluate_build_expression(node, environment, document)
    node isa Number && return node
    node isa AbstractString && return node
    if node isa Symbol
        haskey(environment, node) || throw(ArgumentError(
            "unknown assembly expression name '$node'"))
        return environment[node]
    end
    node isa Expr || throw(ArgumentError(
        "unsupported assembly expression syntax: $node"))
    node.head == :vect && return [evaluate_build_expression(
        argument, environment, document) for argument in node.args]
    node.head == :call || throw(ArgumentError(
        "unsupported assembly expression syntax: $node"))
    operation = node.args[1]
    arguments = [evaluate_build_expression(argument, environment, document)
        for argument in node.args[2:end]]
    operations = Dict{Symbol,Function}(
        :+ => +, :- => -, :* => *, :/ => /, :^ => ^)
    haskey(operations, operation) && return operations[operation](arguments...)
    operation isa Symbol || throw(ArgumentError(
        "unsupported assembly expression operation '$operation'"))
    build_function(operation, arguments, document)
end

function evaluate_build_source(source, environment, document, label)
    syntax = try
        Meta.parse(source)
    catch error
        throw(ArgumentError("invalid $label: $(sprint(showerror, error))"))
    end
    finite_build_value(evaluate_build_expression(syntax, environment, document), label)
end

function materialize_value(value, environment, document, label)
    if value isa AbstractString && startswith(value, "=")
        return evaluate_build_source(strip(value[2:end]), environment, document,
            label)
    elseif value isa AbstractVector
        return [materialize_value(item, environment, document, label)
            for item in value]
    end
    value
end

function toml_value(value)
    value isa BodyReference && return value.name
    value isa AbstractMatrix && return [collect(row) for row in eachrow(value)]
    value isa AbstractVector && return [toml_value(item) for item in value]
    value
end

function parameter_environment(definition, instance, interfaces, document)
    environment = Dict{Symbol,Any}(:pi => pi, :e => MathConstants.e)
    for (name, binding) in interfaces
        environment[Symbol(name)] = BodyReference(binding)
    end
    descriptors = get(definition.document, "parameters", Dict{String,Any}())
    descriptors isa AbstractDict || throw(ArgumentError(
        "assembly '$(definition.name)'.parameters must be a table"))
    allowed = Set(["type"; collect(keys(interfaces)); collect(keys(descriptors))])
    unknown = sort!(collect(setdiff(Set(String.(keys(instance))), allowed)))
    isempty(unknown) || throw(ArgumentError(
        "assembly instance has unknown field '$(first(unknown))'"))
    remaining_parameters = Dict{String,Any}()
    for (name, descriptor) in descriptors
        descriptor isa AbstractDict || throw(ArgumentError(
            "assembly '$(definition.name)' parameter '$name' must be an inline table"))
        supplied = haskey(instance, name)
        required = Bool(get(descriptor, "required", false))
        supplied || haskey(descriptor, "default") || !required ||
            throw(ArgumentError("assembly '$(definition.name)' requires parameter '$name'"))
        supplied || haskey(descriptor, "default") || continue
        remaining_parameters[String(name)] =
            (descriptor, supplied ? instance[name] : descriptor["default"])
    end
    while !isempty(remaining_parameters)
        progress = false
        for name in sort!(collect(keys(remaining_parameters)))
            descriptor, raw = remaining_parameters[name]
            value = try
                materialize_value(raw, environment, document,
                    "assembly parameter '$name'")
            catch error
                message = sprint(showerror, error)
                startswith(message,
                    "ArgumentError: unknown assembly expression name") && continue
                rethrow()
            end
            kind = String(get(descriptor, "kind", "number"))
            if kind == "number"
                value isa Number && !(value isa Bool) && isfinite(value) || throw(ArgumentError(
                    "assembly parameter '$name' must be numeric"))
                value = Float64(value)
            elseif kind == "vector"
                value = numeric_vector(value, "assembly parameter '$name'")
            elseif kind == "string"
                value isa AbstractString || throw(ArgumentError(
                    "assembly parameter '$name' must be a string"))
                value = String(value)
            elseif kind == "boolean"
                value isa Bool || throw(ArgumentError(
                    "assembly parameter '$name' must be true or false"))
            else
                throw(ArgumentError(
                    "assembly parameter '$name' has unsupported kind '$kind'"))
            end
            environment[Symbol(name)] = value
            delete!(remaining_parameters, name)
            progress = true
        end
        progress || throw(ArgumentError(
            "assembly parameters contain unknown names or a dependency cycle: " *
            join(sort!(collect(keys(remaining_parameters))), ", ")))
    end
    calculations = get(definition.document, "calculate", Dict{String,Any}())
    calculations isa AbstractDict || throw(ArgumentError(
        "assembly '$(definition.name)'.calculate must be a table"))
    remaining = Dict(String(name) => source for (name, source) in calculations)
    while !isempty(remaining)
        progress = false
        for name in sort!(collect(keys(remaining)))
            source = remaining[name]
            source isa AbstractString && startswith(source, "=") ||
                throw(ArgumentError("assembly calculation '$name' must be an expression beginning with '='"))
            try
                environment[Symbol(name)] = evaluate_build_source(
                    strip(source[2:end]), environment, document,
                    "assembly calculation '$name'")
                delete!(remaining, name)
                progress = true
            catch error
                message = sprint(showerror, error)
                startswith(message, "ArgumentError: unknown assembly expression name") || rethrow()
            end
        end
        progress || throw(ArgumentError(
            "assembly calculations contain unknown names or a dependency cycle: " *
            join(sort!(collect(keys(remaining))), ", ")))
    end
    environment
end

function interface_bindings(definition, instance, document)
    declarations = get(definition.document, "interface", Dict{String,Any}())
    declarations isa AbstractDict || throw(ArgumentError(
        "assembly '$(definition.name)'.interface must be a table"))
    bindings = Dict{String,String}()
    for (name, declaration) in declarations
        declaration isa AbstractDict || throw(ArgumentError(
            "assembly interface '$name' must be a table"))
        haskey(instance, name) || throw(ArgumentError(
            "assembly '$(definition.name)' requires interface '$name'"))
        binding = instance[name]
        binding isa AbstractString || throw(ArgumentError(
            "assembly interface '$name' must name a body or ground"))
        kind = String(get(declaration, "kind", "body_or_ground"))
        body_table = table_at(document, split(String(binding), '.'))
        body_table isa AbstractDict || throw(ArgumentError(
            "assembly interface '$name' references unknown body '$(binding)'"))
        body_kind = get(body_table, "type", nothing)
        valid = kind == "body" ? body_kind == "rigid_body" :
            kind == "body_or_ground" ? body_kind in ("rigid_body", "ground") : false
        valid || throw(ArgumentError(
            "assembly interface '$name' must reference a $kind"))
        bindings[String(name)] = String(binding)
    end
    bindings
end

function mapped_path(local_path, instance_path, interfaces)
    isempty(local_path) && return instance_path
    port = local_path[1]
    haskey(interfaces, port) || return [instance_path; local_path]
    owner = split(interfaces[port], '.')
    common = common_prefix_length(owner, instance_path)
    suffix = instance_path[common+1:end]
    if length(local_path) >= 2 && local_path[2] == "graphics"
        return [owner; "graphics"; suffix; local_path[3:end]]
    end
    [owner; suffix; local_path[2:end]]
end

function mapped_reference(reference, instance_path, interfaces)
    text = String(reference)
    if startswith(text, "@")
        name = text[2:end]
        haskey(interfaces, name) || throw(ArgumentError(
            "unknown assembly interface '$name'"))
        return interfaces[name]
    end
    mapped_path(split(text, '.'), instance_path, interfaces) |> path_text
end

function transform_field(key, value, environment, document,
        instance_path, interfaces, label)
    materialized = materialize_value(value, environment, document, label)
    if key in REFERENCE_FIELDS
        if materialized isa AbstractString
            materialized = mapped_reference(materialized, instance_path, interfaces)
        elseif materialized isa AbstractVector &&
                all(item -> item isa AbstractString, materialized)
            materialized = [mapped_reference(item, instance_path, interfaces)
                for item in materialized]
        end
    end
    toml_value(materialized)
end

function child_tables(table)
    [(String(key), value) for (key, value) in table if value isa AbstractDict]
end

function emit_template_table!(document, definition, registry, local_path,
        table, instance_path, interfaces, environment, stack)
    kind = get(table, "type", nothing)
    if kind isa AbstractString && haskey(registry, String(kind))
        child_instance = Dict{String,Any}()
        for (key, value) in table
            value isa AbstractDict && continue
            if key == "type"
                child_instance[key] = value
            elseif key in keys(get(registry[String(kind)].document,
                    "interface", Dict{String,Any}()))
                value isa AbstractString || throw(ArgumentError(
                    "nested assembly interface '$key' must be a name"))
                child_instance[key] = mapped_reference(value,
                    instance_path, interfaces)
            else
                child_instance[key] = materialize_value(value, environment,
                    document, "nested assembly field '$key'")
            end
        end
        child_path = mapped_path(local_path, instance_path, interfaces)
        instantiate!(document, registry[String(kind)], registry,
            child_path, child_instance, stack)
        return
    end

    output_path = mapped_path(local_path, instance_path, interfaces)
    output = table_at(document, output_path; create = true)
    for (key, value) in table
        value isa AbstractDict && continue
        haskey(output, key) && throw(ArgumentError(
            "assembly output '$(path_text(output_path)).$key' already exists"))
        output[key] = transform_field(String(key), value, environment,
            document, instance_path, interfaces,
            "assembly table '$(path_text(local_path))'.$key")
    end
    for (key, child) in child_tables(table)
        emit_template_table!(document, definition, registry,
            [local_path; key], child, instance_path, interfaces,
            environment, stack)
    end
end

function instantiate!(document, definition, registry, instance_path,
        instance, stack)
    definition.name in stack && throw(ArgumentError(
        "recursive assembly instantiation: " *
        join([stack; definition.name], " -> ")))
    interfaces = interface_bindings(definition, instance, document)
    environment = parameter_environment(
        definition, instance, interfaces, document)
    next_stack = [stack; definition.name]
    template = endswith(lowercase(definition.path), ".lua") ?
        lua_assembly_template(definition, environment, document) :
        definition.document
    for (key, table) in template
        key in RESERVED_ASSEMBLY_TABLES && continue
        table isa AbstractDict || throw(ArgumentError(
            "assembly template value '$key' must be a table"))
        emit_template_table!(document, definition, registry, [String(key)],
            table, instance_path, interfaces, environment, next_stack)
    end
    document
end

"""
Expand custom assembly instances into ordinary TOML model elements.

Assembly definition paths are relative to the model containing
`model.assemblies`. The returned document no longer depends on those files.
"""
function expand_model_assemblies(document, source_directory;
        dimension = "spatial")
    registry = load_assembly_registry(document, source_directory, dimension)
    isempty(registry) && return document
    instances = assembly_instance_paths!(Vector{Vector{String}}(),
        document, registry)
    for path in instances
        instance = table_at(document, path)
        definition = registry[String(instance["type"])]
        instance_copy = Dict{String,Any}(instance)
        remove_table!(document, path)
        instantiate!(document, definition, registry, path,
            instance_copy, String[])
    end
    model = get(document, "model", nothing)
    model isa AbstractDict && delete!(model, "assemblies")
    document
end

end
