# Practical Mechanical Simulation

[![CI](https://github.com/tjwielenga/PracticalMechanicalSimulation.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/tjwielenga/PracticalMechanicalSimulation.jl/actions/workflows/ci.yml)
[![Documentation](https://github.com/tjwielenga/PracticalMechanicalSimulation.jl/actions/workflows/documentation.yml/badge.svg)](https://github.com/tjwielenga/PracticalMechanicalSimulation.jl/actions/workflows/documentation.yml)
[![codecov](https://codecov.io/gh/tjwielenga/PracticalMechanicalSimulation.jl/graph/badge.svg)](https://codecov.io/gh/tjwielenga/PracticalMechanicalSimulation.jl)

Practical Mechanical Simulation is a Julia package for component-based planar
and spatial mechanical-system simulation. Components contribute local
equations and Jacobian terms to a complete sparse system in which acceleration,
velocity, position, reaction, and applied-force variables remain explicit.

The package performs initial-condition, kinematic, dynamic, static,
quasi-static, and modal analyses as applicable to a model. Models may be built
with the Julia API, read from TOML, or assembled hierarchically with Lua.
Results are stored in portable `.simp` files that can be inspected with the
browser-based SimpView application.

## Installation

Install the registered package from Julia:

```julia
using Pkg
Pkg.add("PracticalMechanicalSimulation")
```

Then load it with:

```julia
using PracticalMechanicalSimulation
```

The source checkout additionally contains the command-line programs,
SimpView, executable models, benchmarks, paper studies, and the complete test
suite. See [Getting Started](docs/getting-started.md) when working from the
repository.

## A first simulation from Julia

The following example loads an included planar pendulum, runs its configured
analysis, and writes a result file:

```julia
using PracticalMechanicalSimulation

model_file = joinpath(
    pkgdir(PracticalMechanicalSimulation),
    "models", "planar", "torsional-spring-pendulum.toml")

model = load_planar_model(model_file)
result = run_planar_model(model)
write_result("torsional-spring-pendulum.simp", result; overwrite = true)
```

To view the result, start SimpView from the root of a source checkout
(Node.js 20.19 or newer is required):

```bash
bin/simpView
```

When SimpView opens in the browser, choose **Open** and select
`torsional-spring-pendulum.simp`. The result can then be animated and its
stored variables plotted.

Planar and spatial models can instead be constructed directly with the
[`Sim2D`](docs/planar/julia-api.md) and
[`Sim3D`](docs/spatial/julia-api.md) builder APIs. Both APIs use the same
loaders, validation, initial-condition correction, sparse assembly, and
analysis code as file-based models.

## Documentation

- [Documentation index](docs/index.md)
- [Planar Modeler User's Guide](docs/planar/README.md)
- [Spatial Modeler User's Guide](docs/spatial/README.md)
- [SimpView and stored results](docs/common/result-viewer.md)
- [Technical Manual](architecture/README.md)
- [Project Summary](PROJECT_SUMMARY.md)
- [Contributing and focused tests](CONTRIBUTING.md)

The hosted documentation is built with Documenter.jl from the same Markdown
manuals kept in this repository.

## Development disclosure

This project was developed through an extended collaboration between Thomas
Wielenga and OpenAI Codex. Wielenga supplied the mechanical formulations,
modeling conventions, engineering judgment, review criteria, and direction of
the project. Codex assisted substantially with code implementation, tests,
documentation, and repository maintenance. The human maintainer reviews
changes through their source, analytical checks, regression tests, and model
behavior, and is responsible for the published package.

## License and citation

The software and documentation are released under the [MIT License](LICENSE).
Citation metadata for software releases is provided in
[`CITATION.cff`](CITATION.cff). The public repository is
[tjwielenga/PracticalMechanicalSimulation.jl](https://github.com/tjwielenga/PracticalMechanicalSimulation.jl).
