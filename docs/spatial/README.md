# Spatial Modeler User's Guide

The spatial modeler  has an extensive set of modeling elements. It reads
a spatial TOML model, assembles the Fully Consistent implicit equations for one
or more rigid bodies, performs dynamic, static, or quasi-static analysis,
performs sparse modal analysis, stores a `.simp` result, and reconstructs the
result in the shared viewer.

The modeling element set includes:

Bodies and geometry:
- rigid bodies with six degrees of freedom;
- flexible body Timoshenko beams with six elastic coordinates;
- ground and body-fixed oriented markers;

Joints and constraints:
- spherical joints;
- perpendicular-axis constraints;
- inplane point-to-plane constraints;
- inline point-to-line constraints;
- hinge orientation constraints;
- orient constraints that fix all three relative rotations;
- revolute joints;
- fixed joints that fix two bodies together;
- cylindrical joints with free translation and rotation about one axis;
- translational joints with one free translation;

Couplers and contacts:
- coordinate couplers that relate rotations and translations;
- ideal gear pairs between parallel or nonparallel revolute axes;
- rack and pinion coupling;
- bumpers - one-sided compliant sphere to plane contacts;
- cams with circular roller or flat-plate followers;
- rolling tires with various options;

Applied forces
- marker-directed applied forces with optional reaction bodies;
- marker-directed applied torques with optional reaction bodies;
- six-component spatial bushings with marker-defined unloaded geometry;
- joint-based applied torques with constant, expression, or torsional
  spring-damper laws;
- marker-to-marker spanning forces with constant, expression, or linear
  spring-damper laws; and
- gravity.

Motion generators:
- translational motion generators;
- spanning motion generators, including constant-distance massless links;
- rotational motion generators for hinge and revolute joints;

Measurements
- reaction-free span and directed-distance measurements;
- joint translations and rotations;

Expressions and state equations:
- expressions can be defined for many elements that are functions of
  model states and measurements;
- user-defined differential and algebraic equations and first-order states can
  be defined and drive existing force and torque elements;

The complete planar program remains separate
under [`docs/planar`](../planar/README.md).

## Running and viewing models

Run any spatial TOML or Lua model from the repository root with:

```bash
./bin/simp3d MODEL --output RESULT --overwrite
```

For example:

```bash
./bin/simp3d models/spatial/revolute-pendulum.toml \
    --output results/examples/spatial/revolute-pendulum.simp --overwrite
```

The model's end time and output rate are normally defined in the input file.
They may be overridden by supplying an end time and sample count after the
model name:

```bash
./bin/simp3d models/spatial/revolute-pendulum.toml 2.0 121
```

A TOML model can also be supplied through standard input:

```bash
./bin/simp3d - < models/spatial/revolute-pendulum.toml
```

Start SimpView and open either the model or its stored result:

```bash
bin/simpView
```

SimpView can run an opened model, display its entered and consistent initial
configurations, follow static-equilibrium progress, and display dynamic,
quasi-static, or modal results. A failed solve retains every accepted result
sample or static iteration that was written before the failure. A completed
`.simp` file can also be converted to CSV or used to recover its embedded
model.

## Example models

The maintained spatial models are in
[`models/spatial`](../../models/spatial/). Run any model with the command
shown above. Representative models include:

Bodies, joints, and flexible members:

- [`free-rotating-body.toml`](../../models/spatial/free-rotating-body.toml),
  an unconstrained rigid body;
- [`offset-cm-spherical-pendulum.toml`](../../models/spatial/offset-cm-spherical-pendulum.toml),
  a body whose modeling reference and center of mass are different;
- [`spherical-pendulum.toml`](../../models/spatial/spherical-pendulum.toml),
  a body with all three relative rotations free;
- [`spherical-perp-pendulum.toml`](../../models/spatial/spherical-perp-pendulum.toml),
  a spherical joint with one rotation removed;
- [`hinge-pendulum.toml`](../../models/spatial/hinge-pendulum.toml) and
  [`revolute-pendulum.toml`](../../models/spatial/revolute-pendulum.toml),
  the primitive and packaged forms of a one-axis rotational joint;
- [`inplane-slider.toml`](../../models/spatial/inplane-slider.toml),
  [`inline-slider.toml`](../../models/spatial/inline-slider.toml), and
  [`inline-pendulum.toml`](../../models/spatial/inline-pendulum.toml),
  examples of the point-to-plane and point-to-line constraints;
- [`cylindrical-joint.toml`](../../models/spatial/cylindrical-joint.toml) and
  [`translational-joint.toml`](../../models/spatial/translational-joint.toml),
  the corresponding user-facing spatial joints;
- [`oriented-free-body.toml`](../../models/spatial/oriented-free-body.toml)
  and [`fixed-two-body-assembly.toml`](../../models/spatial/fixed-two-body-assembly.toml),
  examples of fixed relative orientation and a fixed joint;
- [`slider-crank.toml`](../../models/spatial/slider-crank.toml), a spatial
  mechanism assembled from revolute and inplane constraints;
- [`leaf-spring-assembly.lua`](../../models/spatial/leaf-spring-assembly.lua),
  a hierarchical Lua assembly;
- [`simple-body-surface.lua`](../../models/spatial/simple-body-surface.lua),
  a Lua-built vehicle body with shared surface vertices and separate painted
  and glass patches;
- [`flexible-cantilever.toml`](../../models/spatial/flexible-cantilever.toml)
  and [`rectangular-flexible-cantilever.toml`](../../models/spatial/rectangular-flexible-cantilever.toml),
  floating-reference flexible beams with circular and rectangular sections.

For Lua input conventions, see [Sim3D Lua Assemblies](modeling-assemblies.md).
For beam properties and limitations, see the
[flexible-beam reference](toml-reference.md#flexible-beam) and
[verification notes](flexible-beam-verification.md).

Forces, controls, and measurements:

- [`directed-applied-force.toml`](../../models/spatial/directed-applied-force.toml),
  an off-center marker-directed force;
- [`directed-applied-torque.toml`](../../models/spatial/directed-applied-torque.toml),
  a pure torque along an independently oriented marker axis;
- [`spanning-spring-body.toml`](../../models/spatial/spanning-spring-body.toml),
  a marker-to-marker spring-damper;
- [`bushing-supported-body.toml`](../../models/spatial/bushing-supported-body.toml)
  and [`bushing-pendulum.toml`](../../models/spatial/bushing-pendulum.toml),
  six-component compliant connections;
- [`torsional-spring-pendulum.toml`](../../models/spatial/torsional-spring-pendulum.toml),
  an applied torque using a spring-damper law;
- [`controlled-revolute-pendulum.toml`](../../models/spatial/controlled-revolute-pendulum.toml),
  a PID controller built from user-defined differential and algebraic
  equations;
- [`distance-measures.toml`](../../models/spatial/distance-measures.toml),
  reaction-free span and directed-distance measurements.

The controlled-pendulum equations are described in the
[equation-component reference](toml-reference.md#user-defined-equation-component).
The distance example records range, signed height, and their first two time
derivatives for plotting in SimpView.

Motion generators and coordinate couplers:

- [`constant-speed-revolute-crank.toml`](../../models/spatial/constant-speed-revolute-crank.toml),
  prescribed revolute motion;
- [`constant-speed-translational-slider.toml`](../../models/spatial/constant-speed-translational-slider.toml),
  prescribed translation along a guide;
- [`constant-distance-massless-link.toml`](../../models/spatial/constant-distance-massless-link.toml),
  a spanning motion generator used as a massless link;
- [`screw-motion-coupler.toml`](../../models/spatial/screw-motion-coupler.toml),
  coupled revolute rotation and inline translation;
- [`bevel-gear-pair.toml`](../../models/spatial/bevel-gear-pair.toml),
  an ideal gear pair with perpendicular axes;
- [`planetary-gear-set.toml`](../../models/spatial/planetary-gear-set.toml),
  external and internal gear pairs on a moving carrier;
- [`spur-rack-and-pinion.toml`](../../models/spatial/spur-rack-and-pinion.toml),
  coupled rack translation and pinion rotation.

Contact and tire models:

- [`bouncing-ball.toml`](../../models/spatial/bouncing-ball.toml), compliant
  sphere-to-plane contact;
- the [spatial cam-follower examples](../../examples/spatial/rotating-cam-followers.md),
  covering roller and flat followers on translating and rocking bodies;
- [`driven-rolling-tire.toml`](../../models/spatial/driven-rolling-tire.toml),
  combined longitudinal and lateral tire force;
- [`steered-tire-test-rig.toml`](../../models/spatial/steered-tire-test-rig.toml),
  prescribed steering with braking slip and combined-friction limiting;
- the [single-wheel bristle experiments](bristle-tire-liftoff-experiment.md),
  which examine parked holding, lift-off, shear release, and re-contact;
- [`large-van.lua`](../../models/spatial/large-van.lua) and
  [`large-van-high-cg.lua`](../../models/spatial/large-van-high-cg.lua),
  hierarchical full-vehicle examples.

Static, quasi-static, and modal analyses:

- [`quasistatic-torsional-pendulum.toml`](../../models/spatial/quasistatic-torsional-pendulum.toml),
  consecutive equilibrium solutions under a time-varying torque;
- [`static-initialized-torsional-pendulum.toml`](../../models/spatial/static-initialized-torsional-pendulum.toml),
  static equilibrium followed by dynamic motion;
- [`static-inverted-spherical-pendulum.toml`](../../models/spatial/static-inverted-spherical-pendulum.toml),
  a difficult equilibrium found by dynamic relaxation and Newton correction;
- [`modal-revolute-pendulum.toml`](../../models/spatial/modal-revolute-pendulum.toml),
  a one-freedom analytical modal check;
- [`modal-three-link-pendulum.toml`](../../models/spatial/modal-three-link-pendulum.toml),
  a three-mode coupled check.

## Additional modeling notes

A spatial model may reuse a stored configuration through its
`[initial_conditions]` table. It can select the first or last sample, a sample
number, or the static configuration saved during a dynamically initialized
run. Matching positions, Euler parameters, and relative coordinates transfer
by component name and type; velocities are optional. The new model
recalculates its accelerations, reactions, and forces.

A bushing with stiff translation in all three directions, stiff rotation
around the second marker's $x$ and $y$ axes, and free rotation around its
$z$ axis acts as a compliant revolute joint. It retains six relative freedoms
instead of adding five ideal constraints. This representation includes joint
flexibility directly and can open a closed constraint loop, although its large
coefficients make the system numerically stiff. Use only enough stiffness to
obtain the required deflection accuracy.

Spatial models may contain redundant ideal constraints. After assembling a
consistent initial configuration, the program uses pivoted QR to retain an
independent set of velocity-constraint rows. For every redundant row it
deactivates the corresponding reaction variable and the complete position,
velocity, and acceleration constraint family. The inactive status is retained
in the `.simp` result.

The independent set is not unique. QR may remove a different collection of
scalar constraints than an analyst would choose by inspection while retaining
the same rank and number of freedoms. Reactions associated with redundant
ideal constraints are not individually unique; the reported reactions
correspond to the retained independent equations. Use compliant connections
when the physical distribution of load among nominally redundant connections
matters.

## Model input

Spatial and planar files both use `.toml` for readable model input and `.simp`
for stored results. They are distinguished by
`[model].dimension = "spatial"` and by their separate model directories and
executables.

The [Spatial TOML Reference](toml-reference.md) describes the currently
accepted tables and coordinate conventions. The working example is
[`free-rotating-body.toml`](../../models/spatial/free-rotating-body.toml).

## Julia API

Julia programs can either load an existing model or construct one without
authoring TOML. The complete builder interface is described in the
[Sim3D Julia API](julia-api.md).

The programmatic loading entry points mirror the planar interface:

```julia
using PracticalMechanicalSimulation

loaded = load_spatial_model("models/spatial/free-rotating-body.toml")
result = run_spatial_model("models/spatial/free-rotating-body.toml")
write_result("results/examples/spatial/free-rotating-body.simp", result;
    overwrite = true)
```

`load_spatial_model` accepts either a file path or an `IO` stream. The returned
`LoadedSpatialModel` contains the allocated variable catalog, executable
implicit model, bodies, markers, forces, initial values, and analysis metadata.

Run the Julia-built pendulum with:

```text
julia --project=. examples/spatial/julia_api_pendulum.jl
```

Its result can be opened in SimpView in the same way as a file-based model.

The
[`julia_api_double_pendulum.jl`](../../examples/spatial/julia_api_double_pendulum.jl)
example demonstrates reusable nested assembly functions and the resulting
model and SimpView hierarchy.

## Rotation convention

The physical orientation of a body is its rotation matrix from body components
to global components. TOML accepts the shorter axis-angle form
`[angle, axis_x, axis_y, axis_z]` or an explicit 3 by 3 rotation matrix.
Internally, four normalized Euler parameters provide a nonsingular way to
integrate that matrix. Angular velocity and angular acceleration are resolved
in the body frame. Position, velocity, acceleration, and gravity are resolved
in the global frame.

Three body-fixed pseudo angles provide local orientation columns for the Fully
Consistent formulation. When an angular velocity is selected as a state, its
corresponding pseudo angle is integrated with it. All three pseudo angles stay
in the implicit system so constraints can determine dependent orientation
corrections. Their accumulated values do not define finite orientation and
should not be interpreted as Euler angles.
