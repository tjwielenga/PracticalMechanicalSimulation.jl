# Practical Mechanical Simulation

Practical Mechanical Simulation is a Julia program and written treatment for
component-based mechanical-system simulation. It reads planar and spatial
models from TOML, and reusable planar and spatial assemblies from Lua, then
assembles their complete unreduced implicit equations into sparse systems. The program
performs initial-condition, kinematic, dynamic, static, quasi-static, and modal
analysis as applicable to the model.

The implementation keeps acceleration, velocity, position, reaction, and
applied-force variables explicit. Components contribute local equations and
Jacobian terms to a sparse global system. This is also the working code base for
a paper comparing alternative formulations, so the repository intentionally
contains both the supported modeling program and smaller executable analysis
studies.

The long-term audiences, publication routes, and documentation plan are
recorded in [Project Goals and Distribution Plan](PROJECT_GOALS.md). The
[Project Summary](PROJECT_SUMMARY.md) is the current map of the method,
programs, documentation, verification, limitations, and priorities.

The planar rigid-mechanism modeler is mature for its intended scope. The
spatial modeler has a broad element library and can run substantial models,
including a preliminary full vehicle, but it has had less use and remains
under validation. Their source, models, examples, and documentation are kept
separate so the two levels of maturity remain clear.

## Quick start

The project targets Julia 1.12. SimpView Web additionally requires Node.js
20.19 or newer; Node.js 24 is used by continuous integration. Obtain the
source and enter the repository:

```bash
git clone https://github.com/tjwielenga/PracticalMechanicalSimulation.git
cd PracticalMechanicalSimulation
```

Install the declared Julia packages once:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Run the constant-speed slider-crank and save a compact HDF5 result:

```bash
./bin/simp2d \
    models/planar/constant-speed-slider-crank.toml \
    --output results/examples/planar/constant-speed-slider-crank.simp --overwrite
```

Start the browser-based SimpView application, then use **Open** to select the
stored result:

```bash
bin/simpview-web
```

SimpView can also open a planar or spatial TOML or Lua model to inspect its
entered or consistent initial configuration and
run its configured analysis:

```bash
bin/simpview-web
```

Open `models/spatial/hinge-pendulum.toml`, inspect its two initial
configurations, and use **Run** to follow the configured analysis. **End time**
and **Frames/s** can be changed before starting.

The model can instead be piped through standard input:

```bash
./bin/simp2d - \
    --output results/examples/planar/constant-speed-slider-crank.simp --overwrite \
    < models/planar/constant-speed-slider-crank.toml
```

See [Using the Planar Modeler](docs/planar/using-planar-modeler.md) for the complete
operating workflow. The [documentation index](docs/README.md) links the TOML
User's Guide, result tools, Technical Manual, and paper studies.
Julia users can construct planar models directly through the
[Sim2D Julia API](docs/planar/julia-api.md), without writing TOML.
Calculated and hierarchical planar constructions can instead be packaged as
[Sim2D Lua assemblies](docs/planar/modeling-assemblies.md).

The first spatial model is run separately:

```bash
./bin/simp3d models/spatial/free-rotating-body.toml \
    --output results/examples/spatial/free-rotating-body.simp --overwrite
bin/simpview-web
```

See the [Spatial Modeler User's Guide](docs/spatial/README.md) for the current
capability, conventions, and limitations of the 3D implementation.
Julia users can construct the same models programmatically through the
[Sim3D Julia API](docs/spatial/julia-api.md), without writing TOML or Lua.

## Selected included models

Spatial models:

- [`free-rotating-body.toml`](models/spatial/free-rotating-body.toml) checks
  ballistic translation and finite 3D rotation from axis-angle input.
- [`spherical-pendulum.toml`](models/spatial/spherical-pendulum.toml) constrains
  one body point to ground while leaving all three rotations free.
- [`slider-crank.toml`](models/spatial/slider-crank.toml) combines two spatial
  revolute joints and an inplane primitive into a one-freedom mechanism that
  falls under gravity.
- [`inline-slider.toml`](models/spatial/inline-slider.toml) constrains a body
  point to a marker's $z$-axis and exposes the relative translation as a state
  candidate.
- [`inline-pendulum.toml`](models/spatial/inline-pendulum.toml) suspends a
  vertical pendulum from a horizontal inline so its support translates while
  the body swings under gravity.
- [`oriented-free-body.toml`](models/spatial/oriented-free-body.toml) compares
  a tumbling free body with an identical translating body whose three rotations
  are fixed to a reference marker.
- [`fixed-two-body-assembly.toml`](models/spatial/fixed-two-body-assembly.toml)
  joins two bodies into a freely falling and tumbling rigid assembly.
- [`torsional-spring-pendulum.toml`](models/spatial/torsional-spring-pendulum.toml)
  applies a damped torsional spring through a spatial revolute joint.
- [`controlled-revolute-pendulum.toml`](models/spatial/controlled-revolute-pendulum.toml)
  uses a PID equation component to bring a revolute pendulum to rest at
  45 degrees from vertical.
- [`constant-speed-revolute-crank.toml`](models/spatial/constant-speed-revolute-crank.toml)
  drives a spatial revolute joint through one complete turn and records the
  required drive torque.
- [`constant-speed-translational-slider.toml`](models/spatial/constant-speed-translational-slider.toml)
  prescribes motion along a marker's $z$-axis and records the drive force
  required to lift a guided slider.
- [`constant-distance-massless-link.toml`](models/spatial/constant-distance-massless-link.toml)
  uses a positive spanning-distance constraint as an ideal massless link with
  spherical ends.
- [`screw-motion-coupler.toml`](models/spatial/screw-motion-coupler.toml)
  couples a revolute angle to an inline distance and maps the scalar coupler
  reaction back to torque and force.
- [`bevel-gear-pair.toml`](models/spatial/bevel-gear-pair.toml) couples two
  right-angle revolute axes using a carrier contact point, automatically
  checked collinear pitch tangents, and signed effective radii.
- [`quasistatic-torsional-pendulum.toml`](models/spatial/quasistatic-torsional-pendulum.toml)
  follows the equilibrium of a gravity-loaded pendulum while an applied torque
  changes with model time.
- [`static-initialized-torsional-pendulum.toml`](models/spatial/static-initialized-torsional-pendulum.toml)
  finds the pendulum's loaded equilibrium before restoring an imposed joint
  velocity and beginning dynamic integration.
- [`static-inverted-spherical-pendulum.toml`](models/spatial/static-inverted-spherical-pendulum.toml)
  uses dynamic relaxation to bring a three-rotation spherical pendulum from
  nearly inverted to its hanging equilibrium.
- [`bushing-supported-body.toml`](models/spatial/bushing-supported-body.toml)
  finds the loaded equilibrium of a free body supported only by a spatial
  six-component bushing, then begins damped dynamic motion.
- [`bushing-pendulum.toml`](models/spatial/bushing-pendulum.toml) removes the
  bushing's local $z$ rotational stiffness so a compliant pendulum can swing
  under gravity.
- [`driven-rolling-tire.toml`](models/spatial/driven-rolling-tire.toml) gives a
  wheel lateral velocity and axle torque, then calculates normal,
  longitudinal-slip, lateral-slip, and combined-friction tire forces.
- [`steered-tire-test-rig.toml`](models/spatial/steered-tire-test-rig.toml)
  translates a braked tire while sweeping its spindle through
  $\mathord{\pm}30^\circ$, demonstrating how longitudinal and lateral forces
  share the available friction ellipse.

Planar models:

- [`stage-dependent-force-drop.toml`](models/planar/stage-dependent-force-drop.toml)
  balances a body with static-only supports, then removes them for dynamics.
- [`sliding-block-surface-friction.toml`](models/planar/sliding-block-surface-friction.toml)
  demonstrates one-sided planar contact friction through sliding, sticking,
  and force-driven breakaway.
- [`translational-guide-friction.toml`](models/planar/translational-guide-friction.toml)
  demonstrates compliant sliding, sticking, and force-driven breakaway in a
  planar guide.
- [`controlled-revolute-pendulum.toml`](models/planar/controlled-revolute-pendulum.toml)
  uses a user-defined differential state and algebraic equations to bring a
  gravity-loaded pendulum to a specified angle with PID control.
- [`modal-pendulum.toml`](models/planar/modal-pendulum.toml) linearizes the complete
  implicit equations of a damped torsional-spring pendulum and stores its
  complex mode shape.
- [`constant-speed-slider-crank.toml`](models/planar/constant-speed-slider-crank.toml)
  is a zero-state kinematic mechanism driven through one crank rotation and
  demonstrates marker-attached graphical primitives and named colors.
- [`constant-speed-gear-pair.toml`](models/planar/constant-speed-gear-pair.toml) uses
  carrier-located floating markers and an ideal external gear-pair constraint.
- [`constant-speed-internal-gear-pair.toml`](models/planar/constant-speed-internal-gear-pair.toml)
  derives an internal gear ratio and rotation sign from its carrier geometry.
- [`constant-speed-rack-and-pinion.toml`](models/planar/constant-speed-rack-and-pinion.toml)
  couples a translational rack to a revolute pinion using a pitch radius.
- [`three-pulley-belt-tensioner.toml`](models/planar/three-pulley-belt-tensioner.toml)
  uses stiff elastic tangent spans and a spring-loaded moving tensioner.
- [`perp-guided-slider.toml`](models/planar/perp-guided-slider.toml) combines
  `inplane` and `perp` primitives to form a planar prismatic guide.
- [`translational-joint-slider.toml`](models/planar/translational-joint-slider.toml)
  groups those primitives into a one-freedom translational joint.
- [`translational-distance-slider.toml`](models/planar/translational-distance-slider.toml)
  prescribes one signed marker distance along a reference marker's y-axis.
- [`rotational-coordinate-coupler.toml`](models/planar/rotational-coordinate-coupler.toml)
  couples two optional revolute-joint rotation coordinates at a fixed ratio.
- [`translational-coordinate-coupler.toml`](models/planar/translational-coordinate-coupler.toml)
  couples two unprescribed marker-distance coordinates and uses one as the
  preferred state.
- [`screw-motion-coupler.toml`](models/planar/screw-motion-coupler.toml) mixes a
  rotation coordinate and a distance coordinate to convert screw rotation
  into nut travel and axial force.
- [`fixed-joint-rotating-assembly.toml`](models/planar/fixed-joint-rotating-assembly.toml)
  constructs a fixed joint from revolute and perp primitives.
- [`spanning-spring-pendulum.toml`](models/planar/spanning-spring-pendulum.toml)
  uses the fully expanded nine-equation spanning force with its predefined
  linear spring-damper law.
- [`nonlinear-spanning-force-pendulum.toml`](models/planar/nonlinear-spanning-force-pendulum.toml)
  uses the same component with a length- and length-rate-dependent expression.
- [`torsional-spring-pendulum.toml`](models/planar/torsional-spring-pendulum.toml)
  uses a marker-to-marker torsional spring with a local damping time scale.
- [`nonlinear-expression-pendulum.toml`](models/planar/nonlinear-expression-pendulum.toml)
  defines torque from a revolute angle and angular velocity and differentiates
  that local constitutive expression with dual numbers.
- [`marker-directed-force.toml`](models/planar/marker-directed-force.toml) applies a
  force at a link tip along a fixed ground marker's y axis.
- [`torque-driven-four-bar.toml`](models/planar/torque-driven-four-bar.toml) is a
  one-degree-of-freedom dynamic mechanism with automatic state selection.
- [`bushing-supported-body.toml`](models/planar/bushing-supported-body.toml) exercises
  an unconstrained dynamic body supported by a compliant planar bushing.
- [`relative-coordinate-pendulum.toml`](models/planar/relative-coordinate-pendulum.toml)
  selects an optional revolute-joint angle and rate as its dynamic states.
- [`bushing-pendulum.toml`](models/planar/bushing-pendulum.toml) replaces that ideal
  joint with a translational bushing having zero rotational coefficients.
- [`bouncing-ball.toml`](models/planar/bouncing-ball.toml) exercises DDASSL
  discontinuity detection with a compliant one-sided plane contact.
- [`rotating-cam-follower.toml`](models/planar/rotating-cam-follower.toml)
  follows a smooth marker-fixed cam profile with an explicit roller-contact
  station and a compliant normal force.
- [`rotating-cam-flat-follower.toml`](models/planar/rotating-cam-flat-follower.toml)
  uses the same periodic profile with a marker-oriented flat plate, explicit
  tangency, separation, and recontact.

Simulation settings normally belong in the TOML file. Optional positional
`duration` and `samples` arguments override its end time and output count:

```bash
./bin/simp2d MODEL.toml 2.5 501
```

## Stored results

The `.simp` extension stands for **Simulation Made Practical**. A `.simp` file
is a versioned, compressed HDF5 file containing requested output samples, the
complete canonical variable catalog, active/inactive status, solver
diagnostics, native viewer graphics, and the original TOML model. The command
line creates this file before analysis and appends requested frames as they are
reached. If analysis fails, the same file remains readable and contains the
history through the last accepted state together with the failure message.
Internal integration steps are not stored.

```bash
# Convert all histories to CSV.
julia --project=. bin/export_results.jl \
    results/examples/run.simp results/examples/run.csv

# Export selected histories.
julia --project=. bin/export_results.jl \
    results/examples/run.simp results/examples/motion.csv \
    crank.theta rod.theta rod.R_x

# Recover the embedded model.
julia --project=. bin/extract_model.jl \
    results/examples/run.simp results/examples/recovered.toml
```

Existing output files are protected unless `--overwrite` is supplied.

Regenerate the complete supported planar and spatial example-result set with:

```bash
bin/regenerate-example-results
```

The generated files are separated under `results/examples/planar/` and
`results/examples/spatial/`. They are ignored by Git and can be recreated at
any time.

In the result viewer, left-drag rotates the mechanism and Shift+left-drag
translates it sideways or vertically. Right-drag also translates the view,
the scroll wheel zooms, and Control+left-click restores the original view.
Angular-displacement plots are labeled and displayed in degrees; the `.simp`
file continues to store the underlying values in radians. The independent X
and Y selectors support both time histories and variable-versus-variable
plots. See the [Result Viewer guide](docs/common/result-viewer.md) for all
controls.

Run and view the modal example with:

```bash
./bin/simp2d models/planar/modal-pendulum.toml \
    --output results/examples/planar/modal-pendulum.simp --overwrite
bin/simpview-web
```

## Repository organization

- The [`src/` source guide](src/README.md) explains the execution path through
  `src/common`, `src/planar`, `src/spatial`, and `src/viewer` and gives a
  recommended reading order for the numerical code.
- `apps/SimpViewWeb` is the browser-based viewer and local Julia model service.
- `bin/` contains the `simp2d` and `simp3d` executables and thin Julia launchers for
  simulation, viewing, CSV export, and embedded-model extraction.
- `models/planar` and `models/spatial` contain supported executable TOML and
  Lua models.
- `examples/planar/` contains the planar paper studies and explanatory
  examples; `examples/spatial/` contains the spatial examples.
- `results/` separates ignored generated output into example, benchmark, and
  scratch subdirectories.
- The [`architecture/` Technical Manual](architecture/README.md) separates the
  current mathematical and software design from the earlier design and
  verification record.
- `test/core/` checks the supported program; `test/paper/` preserves executable
  evidence for the alternative analysis formulations.
- `paper/` contains the planar-methods consolidation, claims-to-evidence
  matrix, outline, and working Fully Consistent methods manuscript.

Run the supported program test suite with

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Focused `planar`, `spatial`, `viewer`, and `ddassl` test groups can be selected
through `test/runtests.jl`; see [Contributing](CONTRIBUTING.md).

Run the larger paper verification suite separately with

```bash
julia --project=test/paper -e 'using Pkg; Pkg.instantiate()'
julia --project=test/paper test/paper/run_all.jl
```

Known limitations and current priorities are recorded in the
[Project Summary](PROJECT_SUMMARY.md). The [release checklist](RELEASE_CHECKLIST.md)
records the remaining publication steps for the first public version.

## License and citation

The software and documentation in this repository are released under the
[MIT License](LICENSE). Citation metadata for software releases is in
[`CITATION.cff`](CITATION.cff); the methods paper can be cited separately when
it is published. The public source repository is
[tjwielenga/PracticalMechanicalSimulation](https://github.com/tjwielenga/PracticalMechanicalSimulation).
