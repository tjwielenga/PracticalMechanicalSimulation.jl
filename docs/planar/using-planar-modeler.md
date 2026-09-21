# Using the Planar Modeler

This guide follows one model from TOML input through simulation, graphical
inspection, data export, and a later run. The detailed meaning of every field
is in the [TOML User's Guide](toml-reference.md).

## 1. Begin with a model

The easiest starting point is to copy a maintained model from `models/planar`.
The
[constant-speed slider-crank](../../models/planar/constant-speed-slider-crank.toml) is a
fully driven kinematic example. The
[torque-driven four-bar](../../models/planar/torque-driven-four-bar.toml) is a dynamic
example with an automatically selected state.

Every file begins by declaring a planar model:

```toml
[model]
name = "my_mechanism"
title = "My mechanism"
dimension = "planar"
```

Named element tables then define bodies, markers, connections, forces, and
motion generators. A nested name expresses ownership. For example,
`[crank.end]` is a marker owned by `[crank]`; mechanical connections still
reference their marker endpoints explicitly.

Simulation controls normally remain in the model:

```toml
[simulation]
start_time = 0.0
end_time = 5.0
output_samples = 301
relative_tolerance = 1.0e-7
absolute_tolerance = 1.0e-9
```

## 2. Run and save the analysis

From the repository root:

```bash
./bin/simp2d \
    models/planar/constant-speed-slider-crank.toml \
    --output results/examples/planar/constant-speed-slider-crank.simp --overwrite
```

The executable selects the project relative to its own location. It can
therefore also be invoked by absolute path or through a symbolic link placed on
the shell `PATH`.

The command reports the resolved analysis type, mobility, canonical variable
and equation counts, initialization work, and solver statistics. With
`analysis.mode = "automatic"`, a model with independent states runs
dynamically and a completely driven model runs kinematically.

The optional positional arguments override the duration and number of output
samples without changing the TOML:

```bash
./bin/simp2d MODEL.toml 2.5 501
```

The same command accepts a Lua model. Lua is useful when reusable assemblies
need calculations, loops, or optional construction. See
[Sim2D Lua Assemblies](modeling-assemblies.md):

```bash
./bin/simp2d models/planar/lua-double-pendulum.lua 2.5 151
```

`--output` is optional. Without it the analysis runs and prints its summary but
does not retain a result. Existing result files are protected unless
`--overwrite` is present.

## 3. Inspect the result

Open any saved time-domain or modal result with the same viewer:

```bash
bin/simpview-web
```

The viewer reconstructs the mechanism from the TOML embedded in the result.
Its menus can plot every stored canonical variable. Force and torque controls
show applied and reaction loads, and modal results provide analysis and mode
selectors. The original analysis need not be rerun.

A `.simp` file contains the requested output samples, not every internal
integration step. The command-line runner appends those samples while the
analysis proceeds. If integration fails, the same file remains readable
through its last accepted state and records the failure. It also retains the
complete variable catalog, inactive-variable status, diagnostics, embedded
model, and native viewer graphics in compressed HDF5 form.

## 4. Export or recover data

Export all histories to CSV:

```bash
julia --project=. bin/export_results.jl \
    results/examples/planar/constant-speed-slider-crank.simp \
    results/examples/planar/slider-crank.csv
```

List qualified variable names after the output filename to export only those
columns:

```bash
julia --project=. bin/export_results.jl \
    results/examples/planar/constant-speed-slider-crank.simp \
    results/examples/planar/angles.csv crank.theta rod.theta
```

Time is always included. The `.simp` file remains the complete result record;
CSV is intended for spreadsheets and other exchange tools.

For a modal result, the same command creates two files. An output name of
`pendulum.csv` produces `pendulum-modes.csv`, with one row per mode, and
`pendulum-mode-shapes.csv`, with one row per canonical variable per mode.
Optional qualified variable names filter the mode-shape rows.

Recover the exact TOML stored with a result:

```bash
julia --project=. bin/extract_model.jl \
    results/examples/planar/constant-speed-slider-crank.simp \
    results/examples/planar/recovered.toml
```

This is useful when the original model has moved or when a previous run should
be used as the starting point for a variation.
For a Lua-authored model, extraction returns the expanded primitive TOML that
was actually simulated rather than the original Lua source.

## 5. Choose an analysis sequence

### Static equilibrium followed by dynamics

One model can find equilibrium and then begin dynamic integration:

```toml
[analysis]
mode = "dynamic"
initialization = "static_equilibrium"
```

The equilibrium positions seed the dynamic run. Body velocities declared in
the model are restored afterward, so a pendulum or vehicle can begin moving
from its statically loaded configuration.

Use `static_method = "dynamic_relaxation"` when direct static Newton iteration
does not have a good enough starting configuration. Relaxation approaches the
equilibrium, and a final Newton correction produces the accurate static state.

### Quasi-static continuation

Use `mode = "static"` with multiple output times to solve a sequence of static
positions while time-dependent loads or motion generators change. Inertial
forces are absent. Each converged position predicts the next one.

### Modal analysis

Use `mode = "modal"` to linearize the complete implicit equations at the
starting operating point. Add `initialization = "static_equilibrium"` for modes
about equilibrium. The viewer can select modes from the resulting `.simp`
file.

### Continue from an earlier result

A second model can take its starting configuration from a saved analysis:

```toml
[initial_conditions]
result = "first-run.simp"
sample = "last"
include_velocities = true
```

Matching uses element, variable, and type names rather than column positions.
This permits the second model to change its forces while reusing the bodies and
connections from the first. Set `include_velocities = false` to transfer only
configuration.

For a dynamic result that began with static-equilibrium initialization,
`sample = "static"` recovers the saved equilibrium configuration while
discarding the dynamic initial velocities. For a standalone static result, use
`"first"`, `"last"`, or a one-based sample number.

## 6. Use standard input when convenient

A missing model filename or `-` reads TOML from standard input:

```bash
./bin/simp2d - \
    --output results/examples/planar/run.simp --overwrite < model.toml
```

Relative filenames in `[initial_conditions]` are resolved from the current
working directory for streamed input and from the model's directory for a
named TOML file.

## 7. When a model does not solve

First check marker order, marker orientation, and whether the initial geometry
can satisfy every ideal connection. A singular system usually means that a
motion is not determined or that constraints are redundant in a way that does
not form identifiable scalar families. Preferred states must equal the model
mobility; allowing fallback lets the QR selection choose a usable set.

Zero-mass bodies are supported when connections or force balance determine
their motion. An algebraically force-balanced body must start at a consistent
position or use static-equilibrium initialization to find one.

The command summary, state-selection diagnostics, health snapshots, and stored
active/inactive status provide the first evidence to examine before changing
solver tolerances.
