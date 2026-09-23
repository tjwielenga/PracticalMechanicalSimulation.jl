# Getting Started

## Install the Julia package

From the Julia package prompt, install and load the package:

```julia
pkg> add PracticalMechanicalSimulation

julia> using PracticalMechanicalSimulation
```

A complete planar example using `load_planar_model`, `run_planar_model`, and
`write_result` is shown in the
[repository README](https://github.com/tjwielenga/PracticalMechanicalSimulation.jl#readme).

For programmatic model construction, continue with the
[Sim2D Julia API](planar/julia-api.md) or
[Sim3D Julia API](spatial/julia-api.md).

## Work from the source repository

The command-line programs, browser viewer, models, paper studies, and complete
test suite are distributed with the source repository:

```bash
git clone https://github.com/tjwielenga/PracticalMechanicalSimulation.jl.git
cd PracticalMechanicalSimulation.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Run an included planar model:

```bash
./bin/simp2d models/planar/constant-speed-slider-crank.toml \
    --output results/examples/planar/constant-speed-slider-crank.simp \
    --overwrite
```

Run a spatial model:

```bash
./bin/simp3d models/spatial/free-rotating-body.toml \
    --output results/examples/spatial/free-rotating-body.simp \
    --overwrite
```

## SimpView

SimpView requires Node.js 20.19 or newer. Start its local model service and
browser application with:

```bash
bin/simpView
```

Use **Open** to inspect a TOML or Lua model before simulation, or to load a
stored `.simp` result. See the [Result Viewer guide](common/result-viewer.md)
for animation, plotting, hierarchy, scaling, and camera controls.

## Tests

Run the supported package test suite with:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Focused test groups and the separately maintained paper verification programs
are described in
[Contributing](https://github.com/tjwielenga/PracticalMechanicalSimulation.jl/blob/main/CONTRIBUTING.md).
