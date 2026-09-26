using Documenter
using PracticalMechanicalSimulation

const REPOSITORY_URL =
    "https://github.com/tjwielenga/PracticalMechanicalSimulation.jl"
const REPOSITORY_ROOT = dirname(@__DIR__)
const PUBLISHED_PAGES = [
    "index.md",
    "getting-started.md",
    "common/README.md",
    "common/result-viewer.md",
    "common/numerical-stiffness.md",
    "planar/README.md",
    "planar/using-planar-modeler.md",
    "planar/toml-reference.md",
    "planar/julia-api.md",
    "planar/modeling-assemblies.md",
    "planar/library-api.md",
    "spatial/README.md",
    "spatial/toml-reference.md",
    "spatial/julia-api.md",
    "spatial/modeling-assemblies.md",
    "spatial/flexible-beam-verification.md",
    "spatial/rotating-flexible-blade.md",
]

function published_link(source_file, destination)
    any(prefix -> startswith(destination, prefix),
        ("#", "http:", "https:", "mailto:")) &&
        return destination
    parts = split(destination, '#'; limit = 2)
    path = first(parts)
    isempty(path) && return destination
    target = normpath(joinpath(@__DIR__, dirname(source_file), path))
    inside_documentation = startswith(target, normpath(@__DIR__) * Base.Filesystem.path_separator)
    target_relative_to_docs = inside_documentation ? relpath(target, @__DIR__) : ""
    if inside_documentation &&
            (target_relative_to_docs in PUBLISHED_PAGES || !endswith(path, ".md"))
        return destination
    end
    startswith(target, REPOSITORY_ROOT * Base.Filesystem.path_separator) ||
        return destination
    repository_path = replace(relpath(target, REPOSITORY_ROOT), '\\' => '/')
    kind = isdir(target) ? "tree" : "blob"
    anchor = length(parts) == 2 ? "#" * parts[2] : ""
    "$REPOSITORY_URL/$kind/main/$repository_path$anchor"
end

function documenter_math(text)
    output = String[]
    in_code_fence = false
    in_display_math = false
    for line in split(text, '\n'; keepempty = true)
        stripped = strip(line)
        if !in_code_fence && stripped == "\$\$"
            push!(output, in_display_math ? "```" : "```math")
            in_display_math = !in_display_math
            continue
        end
        if !in_display_math &&
                (startswith(stripped, "```") || startswith(stripped, "~~~"))
            in_code_fence = !in_code_fence
            push!(output, line)
            continue
        end
        if !in_code_fence && !in_display_math
            line = replace(line, r"(?<!\\)\$([^$\n]+?)(?<!\\)\$" => s"``\1``")
        end
        push!(output, line)
    end
    join(output, '\n')
end

function prepare_documentation_source()
    source = mktempdir()
    for relative_path in PUBLISHED_PAGES
        input_path = joinpath(@__DIR__, relative_path)
        output_path = joinpath(source, relative_path)
        mkpath(dirname(output_path))
        text = documenter_math(read(input_path, String))
        text = replace(text, r"(?<=\]\()([^)]+)(?=\))" =>
            destination -> published_link(relative_path, destination))
        write(output_path, text)
    end
    assets = joinpath(@__DIR__, "common", "assets")
    isdir(assets) && cp(assets, joinpath(source, "common", "assets"))
    source
end

const DOCUMENTATION_SOURCE = prepare_documentation_source()

makedocs(
    modules = [PracticalMechanicalSimulation],
    sitename = "Practical Mechanical Simulation",
    source = DOCUMENTATION_SOURCE,
    build = joinpath(@__DIR__, "build"),
    clean = true,
    remotes = nothing,
    checkdocs = :none,
    warnonly = [:cross_references, :missing_docs],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        repolink = REPOSITORY_URL,
        edit_link = nothing,
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting-started.md",
        "Common" => [
            "Overview" => "common/README.md",
            "SimpView and Results" => "common/result-viewer.md",
            "Numerical Stiffness" => "common/numerical-stiffness.md",
        ],
        "Planar Modeler" => [
            "User's Guide" => "planar/README.md",
            "Operating Workflow" => "planar/using-planar-modeler.md",
            "TOML Reference" => "planar/toml-reference.md",
            "Julia API" => "planar/julia-api.md",
            "Lua Assemblies" => "planar/modeling-assemblies.md",
            "Library API" => "planar/library-api.md",
        ],
        "Spatial Modeler" => [
            "User's Guide" => "spatial/README.md",
            "TOML Reference" => "spatial/toml-reference.md",
            "Julia API" => "spatial/julia-api.md",
            "Lua Assemblies" => "spatial/modeling-assemblies.md",
            "Flexible Beam Verification" =>
                "spatial/flexible-beam-verification.md",
            "Rotating Flexible Blade" =>
                "spatial/rotating-flexible-blade.md",
        ],
    ],
)

deploydocs(
    repo = "github.com/tjwielenga/PracticalMechanicalSimulation.jl.git",
    devbranch = "main",
)
