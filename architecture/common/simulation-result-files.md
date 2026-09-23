# Simulation Result Files

Status: shared executable HDF5 format, version 2

## 1. Purpose

A completed run can be retained independently of the solver and viewed or
exported later without repeating the simulation. The public `write_result` and
`read_result` functions use HDF5 under the `.simp` filename extension. The
extension identifies the mechanical-result schema; the underlying file remains
readable by ordinary HDF5 software.

The reader validates the format and version stored inside the HDF5 file rather
than relying only on its `.simp` filename extension.

The format stores plain arrays, strings, and scalar attributes rather than
serialized Julia objects. This keeps the persistent representation independent
of Julia type layouts.

## 2. Version 2 layout

```text
attributes
    format = "PracticalMechanicalSimulation"
    format_version = 2
    title
    status = "running" | "complete" | "failed" | "interrupted"
    status_message
    output_precision = "single" | "double"
    failure_time                                      # optional
model/toml
model/body_references/name                        # spatial only
model/body_references/center_of_mass_marker
model/body_references/center_of_mass_position
model/body_references/center_of_mass_orientation
simulation/time
results/values
variables/name
variables/component
variables/kind
variables/level
variables/scale
variables/active
equations/name
equations/component
equations/kind
equations/level
equations/active
diagnostics/initial_consistency_iterations
diagnostics/static_initialization_iterations  # optional
diagnostics/static_relaxation_cycles           # optional
diagnostics/health/time                         # unhealthy episodes only
diagnostics/health/values
diagnostics/health/physical_error
diagnostics/health/controlled_error
diagnostics/health/amplification
diagnostics/health/dominant_variable
diagnostics/health/dominant_variable_index
diagnostics/health/order
diagnostics/health/step_size
diagnostics/health/severity
diagnostics/health/selected_variables
diagnostics/state_selection/time                # initial choice and changes
diagnostics/state_selection/reason
diagnostics/state_selection/previous
diagnostics/state_selection/selected
diagnostics/static_progress/...                       # retained corrections
graphics/appearance
graphics/choices/.../meshes
graphics/choices/.../tracks
graphics/choices/.../instances
graphics/choices/.../dynamic_surfaces
graphics/choices/.../signals
modal/eigenvalue_real                           # modal only
modal/eigenvalue_imaginary                      # modal only
modal/natural_frequency_hz                      # modal only
modal/damped_frequency_hz                       # modal only
modal/damping_ratio                             # modal only
modal/mode_shape_real                           # modal only
modal/mode_shape_imaginary                      # modal only
modal/equation_error                            # modal only
```

`results/values` has one row per requested output time and one column per
canonical variable. The arrays under `variables` describe those columns in
canonical order. Every acceleration, velocity, position, orientation, reaction,
and applied-force variable retained by the unreduced system is therefore
available for later plotting.

The `active` arrays identify the variables and equations included in the
packed implicit solve. They retain the same canonical ordering as their
metadata groups. A reaction belonging to a redundant constraint remains in the
stored history with value zero, but its inactive flag distinguishes that
suppressed value from a calculated zero reaction. The corresponding inactive
position-, velocity-, and acceleration-equation entries identify the complete
suppressed scalar constraint family.

The original TOML text is embedded verbatim. It provides model topology,
marker geometry, declared simulation settings, and a reproducibility record.
The actual saved time vector records any command-line override.

Spatial results retain CM-based canonical histories. The constant
`model/body_references` arrays store each modeled reference frame's CM marker,
position offset, and orientation offset. Reference-frame and marker motion can
therefore be reconstructed without duplicating their time histories. These
datasets are an optional schema section; older and planar result files
return an empty body-reference list.

`initial_consistency_iterations` records the number of simultaneous Newton
corrections used to establish the starting dynamic solution. When a dynamic
run first solves static equilibrium, the optional
`static_initialization_iterations` dataset records that separate Newton count.
The optional `static_relaxation_cycles` dataset distinguishes a dynamically
relaxed equilibrium from a direct Newton equilibrium. Standalone static runs
store the corresponding per-solution arrays under `diagnostics/static`.
All continued results record the actual analysis mode, number of independent
states, and DDASSL
counts for accepted and rejected steps, equation and Jacobian evaluations,
numerical and symbolic factorizations, corrector failures, and Newton
iterations. Runs with discontinuity surfaces also record root-function
evaluations, located transitions, and BDF-history restarts.

Static relaxation and Newton-correction states are stored under
`diagnostics/static_progress`, including force and torque imbalance,
constraint error, condition estimate, phase, and iteration number.

The command-line runners create the requested `.simp` file before analysis
starts and mark it `running`. Requested dynamic frames are appended in small
batches as DDASSL crosses their times. They are evaluated from the current BDF
history polynomial, the same polynomial used by the completed solution. A
normal finish replaces the working record with the complete result. On solver
failure, the last accepted state, failure time, message, solver counts, and all
previously written frames remain in the same file; no separate partial file is
used.

Dynamic runs also monitor the normalized predictor error in every physical
position and velocity without using those additional variables to accept or
reject a step. Velocity errors are multiplied by the current step size, making
them comparable with position errors. A physical error above 10 accompanied by
at least tenfold amplification over the selected-state error defines an
elevated-health episode. The optional `diagnostics/health` group retains the
accepted canonical state at the peak of each episode. Errors above 25
are marked `check`; lower peaks are marked `warning`. Each record includes the
dominant variable, selected variables, step size, and BDF order.

The first crossing of the check threshold in an episode invokes the pivoted-QR
selection at the accepted configuration. If QR chooses a different valid
partition, the solver replaces the closing state-equation rows, changes the
differential and error-control masks, and refactors the sparse Jacobian. The
complete BDF history is retained. A singular iteration matrix and repeated
corrector failures request the same recovery operation before reducing the
step. `diagnostics/state_selection` records the initial choice and every actual
change, including its cause and the old and new velocity names. The health
snapshot at a change retains the old selection so the analyst can inspect the
configuration that caused it.

A modal result stores one operating-point row in `results/values`. The optional
`modal` group contains one complex canonical mode-shape column per retained
mode, split into portable real and imaginary datasets, together with its
eigenvalue, frequencies, damping ratio, and original-equation error. This is an
optional schema section; older time-domain files remain valid and produce
empty modal arrays when read.

The `/graphics` group contains renderer-independent native HDF5 data for
SimpView. Fixed meshes are shared by compact pose tracks and instances.
Rigidly attached graphics use a fixed local transformation from a shared body
track. Viewer geometry and packed, compressed pose tracks use `Float32`.
Genuinely deforming surfaces
retain sampled vertices. The browser reads these datasets directly; `.simp`
files do not contain an opaque JSON viewer attachment.

## 3. Size policy

Only requested output samples and the last accepted failure state are written
to the ordinary result history; internal solver work is not saved. The few
peak states retained for unhealthy
episodes are stored separately and do not change the requested history.
The solver and its in-memory results always use `Float64`. By default, state
histories, retained static and health states, and modal shapes are converted to
`Float32` when written. A model can set
`simulation.output_precision = "double"` to retain those arrays as `Float64`.
Readers return `Float64` arrays in either case, and saved-result initialization
recalculates consistency. The history matrix is chunked by time and uses the
portable HDF5 shuffle and deflate filters at a moderate compression level.
Existing files are not overwritten unless `overwrite=true` is explicitly
requested or `--overwrite` is passed on the command line.

## 4. Command-line use

```bash
./bin/simp2d model.toml \
    --output results/examples/run.simp
```

For standard input:

```bash
./bin/simp2d - \
    --output results/examples/run.simp < model.toml
```

## 5. Stored-result viewer

SimpView reads the result file without rerunning the simulation. Start it,
then open the desired `.simp` file:

```bash
bin/simpView
```

For a modal result, the viewer provides a mode selector above the history plot.

The viewer generates one display cycle from each stored complex mode shape and
the operating point. Changing the selector switches the animation, title,
history plot, and playback period. Those animation samples are not stored in
the `.simp` file.

When health snapshots are present, the time-domain viewer inserts those exact
states into its in-memory timeline and provides a Health selector above the
plot. Selecting an episode moves the animation, time slider, and plot cursor to
the unhealthy configuration. The inserted display samples remain separate
from the ordinary stored history and CSV export.

It validates the stored canonical catalog against the embedded model, rebuilds
body and marker geometry, and supplies every canonical history to the plot
menus. Body endpoints are inferred from their most widely separated markers.
In-plane constraints use the second marker's oriented axis and the recorded
travel to construct rectangular surfaces extending in the out-of-plane
direction. A plane-contact force uses the same surface representation and adds
its marker-centered contact sphere and surface contact point to the scene.

The result reader is also the boundary used by the CSV exporter. Clients use
the stored metadata rather than depending on a live simulation.

## 6. CSV export

The schema-aware converter writes time followed by fully qualified canonical
variable names:

```bash
bin/simpCSV \
    results/examples/run.simp results/examples/run.csv
```

List variables after the output filename to create a smaller table:

```bash
bin/simpCSV \
    results/examples/run.simp results/examples/motion.csv \
    crank.theta rod.theta rod.R_x
```

`time` is always the first column and need not be requested. Numeric values use
Julia's shortest round-trip `Float64` representation, preserving the stored
values when the CSV is read back. Existing CSV files are protected unless
`--overwrite` is supplied. CSV is an exchange format; the compressed `.simp`
file remains the complete result record.

For a modal result, the requested output name is used as a stem:

```bash
bin/simpCSV \
    results/examples/modal-run.simp results/examples/modal-run.csv
```

This creates `modal-run-modes.csv` and `modal-run-mode-shapes.csv`. The first
file has one row per mode with its complex eigenvalue, natural and damped
frequencies, damping ratio, and equation error. The second uses a long form
with one row per canonical variable per mode. It records the component,
variable, kind, real and imaginary values, magnitude, and phase in degrees.
Qualified variable names supplied on the command line restrict the mode-shape
rows but do not change the mode summary.

## 7. Model extraction

The embedded model can be recovered verbatim for inspection, editing, or a new
run:

```bash
bin/simpExtract \
    results/examples/run.simp results/examples/recovered.toml
```

An existing output is protected unless `--overwrite` is supplied. Use `-` as
the output name to write only the TOML text to standard output:

```bash
bin/simpExtract results/examples/run.simp -
```

This permits pipelines such as redirecting the model to a temporary file or
another TOML-aware tool without mixing status text into the extracted model.

## 8. Reusing a saved configuration

A later TOML model can use a stored sample as its initial configuration:

```toml
[initial_conditions]
result = "static-equilibrium.simp"
sample = "static"
include_velocities = false
```

This reads the portable history and metadata directly; it does not rerun the
embedded model. Matching is by qualified component name, element type,
canonical variable name, and variable kind. The transfer is deliberately
limited to configuration and, when requested, velocity variables. The new
model recomputes accelerations, reactions, and force variables during initial
consistency correction.

For a dynamic source run that used static-equilibrium initialization,
`sample = "static"` selects the first stored configuration after confirming the
corresponding initialization diagnostic is present. With velocity transfer
disabled, this recovers the equilibrium position while ignoring the dynamic
velocities restored before the source run began.
