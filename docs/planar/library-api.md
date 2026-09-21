# Julia Library API

Use the separate [Sim2D Julia API](julia-api.md) when constructing a planar
model programmatically. This page documents the lower-level run and result
operations that work with both TOML and Julia-built models.

The small public API supports running TOML models and processing stored results.
Lower-level assembly modules are available to the paper examples but are still
internal extension points and may change as the spatial design develops.

```julia
using PracticalMechanicalSimulation
```

## Load and run

```julia
loaded = load_planar_model("models/planar/torque-driven-four-bar.toml")
result = run_planar_model("models/planar/torque-driven-four-bar.toml")
```

`load_planar_model` accepts a file path, an `IO` containing TOML, an
`AbstractDict` model document, or a `Sim2D.Model`, and returns a
`LoadedPlanarModel`. Its useful inspection fields include `title`, `layout`,
`analysis`, `simulation`, `state_selection`, `initial_conditions`,
`initial_values`, and the active canonical variable and equation indices.
`initial_condition_weights`, `initial_variable_weights`, and
`imposed_initial_variable_indices` expose the resolved IC correction policy;
`state_selection.imposed_initial_variables` gives the corresponding qualified
names. Each body entry in `initial_condition_weights` records its TOML scale,
effective mass, characteristic length, translational weight $s\,m$, and
rotational weight $s\,m\,L^2$.
When a TOML file imports a saved configuration, `initial_conditions` records
the resolved result path, requested selector, selected sample and time, whether
velocities were included, and the qualified variables transferred.

`run_planar_model` accepts the same source forms and supports keyword overrides:

```julia
result = run_planar_model(path;
    start_time = 0.0,
    end_time = 2.0,
    duration = nothing,
    samples = 401)
```

`duration`, when supplied, takes precedence over `end_time` and is measured
from the effective start time. A result contains sampled `times`, canonical
`states`, the `loaded` model, `analysis_mode`, initialization diagnostics, and
either dynamic solver statistics or static continuation diagnostics.

For `mode = "modal"`, the result contains one operating-point state plus
`eigenvalues`, `natural_frequencies_hz`, `damped_frequencies_hz`,
`damping_ratios`, complex canonical `mode_shapes`, and `equation_errors`. The
modal solver uses one sparse factorization. `solve_modal_system` is also
exported for linearizing a prepared consistent state and derivative directly.

When `[analysis]` requests `initialization = "static_equilibrium"`, the dynamic
result's `static_initialization_iterations` records the static Newton
corrections separately from `initial_consistency_iterations`. It is `nothing`
for an ordinary dynamic run. The first sampled state is the dynamically
consistent state at the statically determined configuration.
When dynamic relaxation is selected, `static_relaxation_cycles` records its
number of pseudo-time intervals; `static_initialization_iterations` records
only the final exact Newton corrections.

For a dynamic model containing `plane_contact` forces,
`result.solution.events` records each located force-law transition. Every event
contains its time, one or more root-function indices and crossing directions,
and the interpolated canonical state and derivative. The solver statistics also
report `root_evaluations`, `events_found`, and `history_restarts`.
The default root response is a soft history restart for continuous compliant
force transitions: recent values are retained, BDF order is capped at two,
the numerical iteration matrix is refreshed, and the next step is halved.
The modeling runner selects this policy from the contact element semantics;
there is no corresponding TOML setting. The lower-level DDASSL API offers the
single `root_restart` policy with `:hard`, `:soft`, and `:none` values and uses
the conservative `:hard` default.

Dynamic results also provide `state_selection_changes`. The initial record
contains the starting velocity coordinates. A recovery record identifies a
high physical predictor error, singular iteration matrix, or repeated
corrector failure and gives the previous and replacement coordinates. Runtime
reselection retains the BDF history but rebuilds the sparse Jacobian pattern
and its UMFPACK symbolic factorization.

For sparse dynamic systems, `factorizations` counts numerical Newton-matrix
factorizations and `symbolic_factorizations` counts UMFPACK analyses of the
sparse pattern. The level-scaled numerical factors are reused across accepted
step attempts as a modified-Newton matrix while the scaled BDF coefficient
remains within the configured ratio window, for at most
`maximum_factorization_steps` attempts (five by default). All Newton
corrections in an attempt use the same factors. If an old matrix cannot
converge, the solver restarts that corrector from its predictor with a freshly
evaluated and numerically factored matrix before reducing the step. A located
event also refreshes the numerical values because it may mark a force or
stiffness discontinuity.

The symbolic result is normally reused for the whole analysis, including these
numeric refreshes. It is rebuilt after repeated corrector failures, and a new
analysis following a change of selected states necessarily begins with a new
pattern analysis. `corrector_failures` records only failures remaining after
the fresh-matrix retry; a stale-matrix refresh that recovers immediately is not
counted as a rejected corrector.

`compile_time_expression(text, parameters)` compiles the restricted expression
language documented in the [TOML User's Guide](toml-reference.md). It returns a
callable scalar function of time and does not evaluate arbitrary Julia code.

Saved-result initialization is selected in the TOML rather than by a run-time
API keyword. This keeps the initialization source and transfer policy in the
model's reproducibility record. Relative result paths passed to
`load_planar_model(path)` or `run_planar_model(path)` are resolved beside that
model file; stream inputs use the current working directory.

## Store and read

```julia
write_result("run.simp", result; overwrite = false, compression = 3)
stored = read_result("run.simp")
```

State histories use the model's `simulation.output_precision`, which defaults
to `"single"`. `write_result(...; output_precision = :double)` can override it
for a particular API write. Calculations and arrays returned by `read_result`
remain `Float64` in either case.

`read_result` returns a `StoredSimulationResult` containing only portable data:
times, the history matrix, canonical variable and equation metadata,
active/inactive flags, embedded TOML, analysis metadata, and solver statistics.
It does not rebuild a live model or rerun the analysis.

`stored.health_snapshots` contains one `StoredHealthSnapshot` for the peak of
each elevated physical-error episode. A snapshot records the retained canonical
values, time, BDF step and order, dominant variable, selected variables,
physical and selected-state errors, amplification, and warning/check severity.
The ordinary stored history still contains only the requested output samples.

`stored.state_selection_changes` contains `StoredStateSelectionChange`
records. The first record gives the initial velocity selection. Later records
give the time, recovery reason, previous selection, and replacement selection.
The corresponding solver count is
`stored.solver_statistics.state_reselections`.

For a modal result, the additional `modal_*` fields contain the complex
eigenvalues and canonical mode shapes, frequencies, damping ratios, equation
errors, and shift. Ordinary time-domain results leave these arrays empty.

History column `i` is described by `stored.variable_names[i]`,
`stored.variable_components[i]`, `stored.variable_kinds[i]`,
`stored.variable_levels[i]`, and `stored.variable_active[i]`.

## CSV and model extraction

```julia
export_result_csv("run.csv", stored; overwrite = true)
export_result_csv("motion.csv", "run.simp";
    variables = ["crank.theta", "coupler.theta"], overwrite = true)

extract_model("recovered.toml", stored; overwrite = true)
```

Stream forms are also supported:

```julia
export_result_csv(stdout, stored; variables = ["crank.theta"])
extract_model(stdout, stored)
```

CSV always begins with `time`. Requested variable names are qualified as
`component.variable`; requesting `time` explicitly is an error.

For a modal stored result, the path form treats the requested filename as a
stem and returns a named tuple containing the generated `-modes.csv` and
`-mode-shapes.csv` paths. To write modal data to streams, use
`export_modal_csv(summary_output, shapes_output, stored)`. Qualified variable
names filter the long-form mode-shape rows without changing the mode summary.

## Command entry points

The exported `planar_model_main`, `export_result_main`, and
`extract_model_main` functions implement the wrappers under `bin/`. They
accept an argument vector and return a process-style status code, which makes
them directly testable without starting a subprocess.

## Internal extension boundary

`AutomaticAnalysis`, `PlanarComponentAssembly`, `PlanarModeling`, and
`HistoricalDDASSL` are submodules rather than top-level public API. The paper
examples intentionally use them to expose formulation details. New program
components currently require coordinated changes at this internal boundary;
follow [Adding a planar component](../../architecture/planar/adding-planar-component.md).
