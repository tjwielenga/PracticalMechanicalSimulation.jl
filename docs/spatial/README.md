# Spatial Modeler User's Guide

The spatial modeler now has a small but complete first vertical slice. It reads
a spatial TOML model, assembles the Fully Consistent implicit equations for one
or more rigid bodies, performs dynamic, static, or quasi-static analysis,
performs sparse modal analysis, stores a `.simp` result, and reconstructs the
result in the shared viewer.

The present element set is deliberately limited to:

- free rigid bodies with optional body-reference and CM-marker separation;
- straight floating-reference Timoshenko beams with six elastic coordinates;
- ground and body-fixed oriented markers;
- spherical joints;
- perpendicular-axis constraints;
- inplane point-to-plane constraints;
- inline point-to-line constraints with optional translation coordinates;
- hinge orientation constraints;
- orient constraints fixing all three relative rotations;
- revolute joints;
- fixed joints;
- linear coordinate couplers among hinge or revolute rotations and inline
  translations;
- ideal gear pairs between parallel or nonparallel revolute axes, with
  carrier contact geometry and generated floating contact markers;
- reaction-free span and directed-distance measurements;
- marker-directed applied forces with optional reaction bodies;
- six-component spatial bushings with marker-defined unloaded geometry;
- one-sided compliant sphere-plane contacts;
- marker-fixed planar cam profiles extruded along their local $z$ axes, with
  circular roller or flat-plate followers;
- rolling tires with normal compliance, expression forces or optional
  load-dependent bristle states, and combined-slip limiting;
- user-defined scalar algebraic equations and first-order states that can
  drive existing force and torque elements;
- joint-based applied torques with constant, expression, or torsional
  spring-damper laws;
- marker-based translational motion generators;
- spanning motion generators, including constant-distance massless links;
- rotational motion generators for hinge and revolute joints;
- marker-to-marker spanning forces with constant, expression, or linear
  spring-damper laws; and
- gravity.

Arbitrary surface-to-surface contact and general floating torques are not yet
supported.
The complete planar program remains separate
under [`docs/planar`](../planar/README.md).

## Running the example

From the repository root:

```bash
./bin/simp3d models/spatial/free-rotating-body.toml \
    --output results/examples/spatial/free-rotating-body.simp --overwrite
```

The duration and sample count may be overridden in the same way as for
`simp2d`:

```bash
./bin/simp3d models/spatial/free-rotating-body.toml 2.0 401
```

A model can also be supplied through standard input:

```bash
./bin/simp3d - < models/spatial/free-rotating-body.toml
```

Start `bin/simpView`, open a model, and select **Static equilibrium** to
watch accepted Newton corrections and DDASSL pseudo-time steps as they are
calculated. Before static motion begins, the viewer records the entered model
configuration and the configuration after initial-condition consistency. A
failed solve retains those configurations and every accepted static step for
inspection.

View the stored motion by starting SimpView and opening the `.simp` file:

```bash
bin/simpView
```

The example body appears as an oriented box, its body-fixed marker appears as a
sphere, and gravity appears as an applied-force arrow. The `.simp` file can
also be converted to CSV or used to recover the embedded TOML model with the
same common result tools used by planar analyses.

The offset-CM spherical pendulum defines its marker geometry from the pivot
reference frame and names a separate CM marker. The loader retains CM-based
equations and result variables while the viewer reconstructs the modeled
reference geometry:

```bash
./bin/simp3d models/spatial/offset-cm-spherical-pendulum.toml \
    --output results/examples/spatial/offset-cm-spherical-pendulum.simp --overwrite
bin/simpView
```

The leaf spring is the first Lua-authored hierarchical model. Its module is
loaded with Lua's `require` function and the returned `leaf_spring` function is
called directly with a table. The resulting ordinary Sim3D model is stored in
the result:

```bash
./bin/simp3d models/spatial/leaf-spring-assembly.lua \
    --output results/examples/spatial/leaf-spring-assembly.simp --overwrite
bin/simpView
```

See [`modeling-assemblies.md`](modeling-assemblies.md) for the Lua input and
assembly-definition conventions.

The spatial flexible cantilever fixes one generated beam end marker to ground.
Its six elastic coordinates produce two equal first bending frequencies, while
SimpView can amplify the three-dimensional centerline deformation:

```bash
./bin/simp3d models/spatial/flexible-cantilever.toml \
    --output results/examples/spatial/flexible-cantilever.simp --overwrite
bin/simpView
```

Beam input may give the canonical section properties directly or use a solid
`circular` or `rectangular` section preset. The presets derive area, bending
inertias, torsion constant, and shear coefficients from ordinary dimensions;
mass and shear modulus can additionally be derived from density and Poisson's
ratio. See the [TOML reference](toml-reference.md#flexible-beam).
The [beam verification notes](flexible-beam-verification.md) record the
analytical static, modal, damped-motion, and rigid-motion checks.

The [controlled pendulum](../../models/spatial/controlled-revolute-pendulum.toml)
uses a PID equation component to bring a revolute pendulum to rest at
45 degrees from downward vertical. Run it with:

```bash
./bin/simp3d models/spatial/controlled-revolute-pendulum.toml \
    --output results/examples/spatial/controlled-revolute-pendulum.simp --overwrite
```

See the [equation-component reference](toml-reference.md#user-defined-equation-component)
for its equations, initial values, and static behavior.

The [single-wheel bristle experiment](bristle-tire-liftoff-experiment.md)
checks parked side-slope holding, lift-off, shear release, and re-contact.
Its fore–aft companion uses a wheel brake torque to hold the slope.
Its tire curves are illustrative rather than measured.

The measurement example records both the marker-to-marker range and the
signed height along the ground marker's $z$-axis:

```bash
./bin/simp3d models/spatial/distance-measures.toml \
    --output results/examples/spatial/distance-measures.simp --overwrite
bin/simpView
```

Select `range.distance`, `range.velocity`, `range.acceleration`, or the
corresponding `height` variables in the viewer plot.

The spatial cam examples use the same closed profile with circular and flat
followers. The profile lies in its marker's local $x$-$y$ plane and is
extruded along local $z$ for contact and display:

See the [spatial cam-follower example notes](../../examples/spatial/rotating-cam-followers.md)
for the marker and joint conventions.

```bash
./bin/simp3d models/spatial/rotating-cam-roller-follower.toml \
    --output results/examples/spatial/rotating-cam-roller-follower.simp --overwrite
./bin/simp3d models/spatial/rotating-cam-flat-follower.toml \
    --output results/examples/spatial/rotating-cam-flat-follower.simp --overwrite
bin/simpView
```

The quasi-static pendulum solves a sequence of equilibria while its applied
torque increases with model time:

```bash
./bin/simp3d models/spatial/quasistatic-torsional-pendulum.toml \
    --output results/examples/spatial/quasistatic-torsional-pendulum.simp --overwrite
bin/simpView
```

The static-initialized pendulum first finds the equilibrium under gravity and
its torsional spring. It then restores an imposed relative angular velocity
and begins dynamic integration:

```bash
./bin/simp3d models/spatial/static-initialized-torsional-pendulum.toml \
    --output results/examples/spatial/static-initialized-torsional-pendulum.simp --overwrite
bin/simpView
```

A separate spatial model may reuse a stored configuration through its
`[initial_conditions]` table. It can select the first or last sample, a sample
number, or the static configuration from a dynamically initialized result.
Matching positions, Euler parameters, and relative coordinates transfer by
component name and type; velocities are optional. The new model recalculates
its accelerations, reactions, and forces.

The inverted spherical pendulum is a more difficult static problem. Its
spherical joint leaves all three rotations free, and its center of mass starts
five degrees away from the exact unstable inverted equilibrium. Dynamic
relaxation lets it travel almost a full link length before the final static
Newton correction:

```bash
./bin/simp3d models/spatial/static-inverted-spherical-pendulum.toml \
    --output results/examples/spatial/static-inverted-spherical-pendulum.simp --overwrite
bin/simpView
```

The spherical-joint pendulum is the first constrained spatial model:

```bash
./bin/simp3d models/spatial/spherical-pendulum.toml \
    --output results/examples/spatial/spherical-pendulum.simp --overwrite
bin/simpView
```

Its viewer shows the joint, its reaction force on the first marker, gravity,
and the body motion.

The spherical-plus-perp example removes one of the spherical joint's three
relative rotations:

```bash
./bin/simp3d models/spatial/spherical-perp-pendulum.toml \
    --output results/examples/spatial/spherical-perp-pendulum.simp --overwrite
bin/simpView
```

The inplane example lets a body slide freely while keeping one of its marker
points in a plane:

```bash
./bin/simp3d models/spatial/inplane-slider.toml \
    --output results/examples/spatial/inplane-slider.simp --overwrite
bin/simpView
```

The inline example constrains a body point to the second marker's $z$-axis and
uses the relative axial velocity as a preferred state:

```bash
./bin/simp3d models/spatial/inline-slider.toml \
    --output results/examples/spatial/inline-slider.simp --overwrite
bin/simpView
```

The inline-pendulum example gives a vertical pendulum both motion along a
horizontal guide and transverse motion that starts it swinging:

```bash
./bin/simp3d models/spatial/inline-pendulum.toml \
    --output results/examples/spatial/inline-pendulum.simp --overwrite
bin/simpView
```

The spatial slider-crank combines two revolute joints with an inplane
primitive. It begins with its crank at $45$ degrees and collapses under
gravity:

```bash
./bin/simp3d models/spatial/slider-crank.toml \
    --output results/examples/spatial/slider-crank.simp --overwrite
bin/simpView
```

The hinge pendulum combines a spherical joint with two perpendicular-axis
constraints, leaving only rotation about the ground marker's $z$-axis:

```bash
./bin/simp3d models/spatial/hinge-pendulum.toml \
    --output results/examples/spatial/hinge-pendulum.simp --overwrite
bin/simpView
```

The revolute-pendulum example packages those same five constraint families
into one revolute joint and exposes its optional relative angle:

```bash
./bin/simp3d models/spatial/revolute-pendulum.toml \
    --output results/examples/spatial/revolute-pendulum.simp --overwrite
bin/simpView
```

The orient example compares two translating bodies. The orange body tumbles,
while the blue body has all three rotations fixed to a ground marker and
retains that orientation throughout the motion:

```bash
./bin/simp3d models/spatial/oriented-free-body.toml \
    --output results/examples/spatial/oriented-free-body.simp --overwrite
bin/simpView
```

The fixed-joint example joins two bodies into one freely moving L-shaped
assembly. The assembly falls and tumbles while preserving the marker
coincidence and relative orientation:

```bash
./bin/simp3d models/spatial/fixed-two-body-assembly.toml \
    --output results/examples/spatial/fixed-two-body-assembly.simp --overwrite
bin/simpView
```

The spanning-force example suspends a free body from a ground marker with a
linear spring-damper. Its transverse initial velocity makes the body swing
through three-dimensional space:

```bash
./bin/simp3d models/spatial/spanning-spring-body.toml \
    --output results/examples/spatial/spanning-spring-body.simp --overwrite
bin/simpView
```

The bushing example finds the gravity-loaded static position of a free body,
restores its declared velocity, and then integrates the damped six-freedom
motion:

```bash
./bin/simp3d models/spatial/bushing-supported-body.toml \
    --output results/examples/spatial/bushing-supported-body.simp --overwrite
bin/simpView
```

The plots include the three local translations and Bryant angles, their
rates, and the local and global bushing loads.

The bushing-pendulum example sets the local $z$ rotational stiffness to zero.
It starts horizontally and swings under gravity while the translational part
of the bushing remains compliant:

```bash
./bin/simp3d models/spatial/bushing-pendulum.toml \
    --output results/examples/spatial/bushing-pendulum.simp --overwrite
bin/simpView
```

This is also a useful general modeling pattern. A bushing with stiff
translation in all three directions, stiff rotation around the second
marker's $x$ and $y$ axes, and free rotation around its $z$ axis acts as a
compliant revolute joint. It retains six relative freedoms rather than adding
five ideal constraints. Force balance keeps five motions small and leaves the
sixth free. This representation can include joint flexibility directly and
can open a closed constraint loop, although its large coefficients make the
system numerically stiff. Use only enough stiffness to obtain the needed
deflection accuracy.

The driven-tire example mounts a wheel on a revolute carried by a spindle that
can translate in the road plane. It starts with lateral velocity, and an
applied axle torque creates longitudinal slip. The tire then develops both
cornering and traction forces:

```bash
./bin/simp3d models/spatial/driven-rolling-tire.toml \
    --output results/examples/spatial/driven-rolling-tire.simp --overwrite
bin/simpView
```

The tire plots include normal deflection and force, transport velocities,
unnormalized slip velocities, slip ratio, slip angle, trial forces, and the
friction-limited forces. The blue applied-force arrow begins at the calculated
road contact point. Its ground reaction is hidden unless ground loads are
enabled.

The steered test rig prescribes the carriage travel and a complete sinusoidal
steering cycle through $+30^\circ$ and $-30^\circ$. The wheel starts at its
free-rolling speed and a constant resisting axle torque then creates braking
slip. Its smooth lateral-force law reaches about 85 percent of its independent
limit at $10^\circ$ slip angle and 95 percent at $15^\circ$. The friction
ellipse reduces that force when longitudinal and lateral demand occur
together:

```bash
./bin/simp3d models/spatial/steered-tire-test-rig.toml \
    --output results/examples/spatial/steered-tire-test-rig.simp --overwrite
bin/simpView
```

Select `steering_joint.theta (deg)`, `tire.slip_angle (deg)`,
`tire.longitudinal_force`, or `tire.lateral_force` in the viewer plot to
compare the commanded steering, resulting slip angle, braking force, and
combined-friction-limited cornering force.

The applied-force example uses the $z$-axis of an oriented ground marker to
direct an off-center force on a free body:

```bash
./bin/simp3d models/spatial/directed-applied-force.toml \
    --output results/examples/spatial/directed-applied-force.simp --overwrite
bin/simpView
```

The torsional-spring pendulum applies a spring-damper torque through a
revolute joint. The same element also accepts a constant torque or an
expression:

```bash
./bin/simp3d models/spatial/torsional-spring-pendulum.toml \
    --output results/examples/spatial/torsional-spring-pendulum.simp --overwrite
bin/simpView
```

The modal version of the revolute pendulum provides the first analytical
spatial modal check. It has one rotational freedom, and its frequency and
damping follow directly from the spring, damper, and effective inertia about
the pin:

```bash
./bin/simp3d models/spatial/modal-revolute-pendulum.toml \
    --output results/examples/spatial/modal-revolute-pendulum.simp --overwrite
bin/simpView
```

The viewer supplies a mode selector when more than one mode is stored.

The three-link modal pendulum checks several coupled modes against an
independent three-coordinate mass, damping, and stiffness calculation. It
first confirms the hanging static equilibrium and then stores three modes:

```bash
./bin/simp3d models/spatial/modal-three-link-pendulum.toml \
    --output results/examples/spatial/modal-three-link-pendulum.simp --overwrite
bin/simpView
```

The constant-speed crank uses a rotational motion generator to prescribe one
complete revolution while gravity determines the required drive torque:

```bash
./bin/simp3d models/spatial/constant-speed-revolute-crank.toml \
    --output results/examples/spatial/constant-speed-revolute-crank.simp --overwrite
bin/simpView
```

The constant-speed translational slider combines an inline and orient
constraint into a prismatic guide, then prescribes distance along the ground
marker's $z$-axis. Its drive force balances gravity:

```bash
./bin/simp3d models/spatial/constant-speed-translational-slider.toml \
    --output results/examples/spatial/constant-speed-translational-slider.simp --overwrite
bin/simpView
```

The massless-link example uses a constant-distance spanning motion to suspend
a freely rotating body. The link removes only radial motion, so the body can
swing in three dimensions:

```bash
./bin/simp3d models/spatial/constant-distance-massless-link.toml \
    --output results/examples/spatial/constant-distance-massless-link.simp --overwrite
bin/simpView
```

The screw example uses one coordinate coupler to convert revolute rotation
into inline translation. The coupler automatically enables both relative
coordinates and maps its scalar reaction back to a torque on the screw and a
force on the nut:

```bash
./bin/simp3d models/spatial/screw-motion-coupler.toml \
    --output results/examples/spatial/screw-motion-coupler.simp --overwrite
bin/simpView
```

The bevel-gear example couples two revolute axes at right angles. Their
centers share an apex, so the contact direction must come from the two pitch
tangents rather than the cross product of the two center-to-contact vectors:

```bash
./bin/simp3d models/spatial/bevel-gear-pair.toml \
    --output results/examples/spatial/bevel-gear-pair.simp --overwrite
bin/simpView
```

The blue action arrow is the force applied to the first gear named by the gear
pair. The red reaction arrow is the equal and opposite force on the second
gear. Each gear body explicitly owns a `shape = "gear"` graphic positioned and
oriented by its referenced pitch marker. Supplying `cone_height` gives the
current smooth conical-frustum representation; omitting it gives a cylinder.
The gear-pair element does not infer either graphic.

The planetary example uses two ideal gear pairs at carrier-fixed contact
markers. The sun-to-planet pair is external, the planet-to-ring pair is
internal, and the planet bearing and both contact points move with the
carrier:

```bash
./bin/simp3d models/spatial/planetary-gear-set.toml \
    --output results/examples/spatial/planetary-gear-set.simp --overwrite
bin/simpView
```

The spur rack-and-pinion example uses the inline base marker as its pitch
frame. A separate orient constraint completes the rack guide, while the
rack-and-pinion element couples the remaining translation to the pinion
rotation and generates the contact-force markers:

```bash
./bin/simp3d models/spatial/spur-rack-and-pinion.toml \
    --output results/examples/spatial/spur-rack-and-pinion.simp --overwrite
bin/simpView
```

Spatial models may contain redundant ideal constraints. After assembling a
consistent initial configuration, the program uses pivoted QR to retain an
independent set of velocity-constraint rows. For every redundant row it
deactivates the corresponding reaction variable and the complete position,
velocity, and acceleration constraint family. The inactive status is retained
in the `.simp` result.

The independent set is not unique. QR may remove a different collection of
scalar constraints than an analyst would choose by inspection, while retaining
the same constraint rank and number of freedoms. For example, a spatial model
of a planar four-bar might retain every joint-center constraint in the
$z$ direction and remove redundant perpendicular-axis constraints around the
closed loop. This does not add compliance or another freedom: the retained
equations still imply the omitted constraints on the modeled branch.

Reactions associated with redundant ideal constraints are not individually
unique. The reported reactions correspond to the independent equation set
chosen by QR. A compliant model, such as one using bushings, is needed when the
physical distribution of load among nominally redundant connections matters.

## Cylindrical and translational joints

Two familiar joints can be constructed directly from the spatial primitives,
so they do not need separate element types.

A cylindrical joint combines an `inline` and a `hinge` using the same ordered
marker pair:

```toml
[cylinder.translation]
type = "inline"
markers = ["shaft.axis", "housing.axis"]
translation_coordinates = true

[cylinder.rotation]
type = "hinge"
markers = ["shaft.axis", "housing.axis"]
rotation_coordinates = true
```

It leaves two relative degrees of freedom: translation along and rotation
about the second marker's $z$-axis.

A translational joint replaces the hinge with an `orient`:

```toml
[slider.translation]
type = "inline"
markers = ["slider.axis", "rail.axis"]
translation_coordinates = true

[slider.orientation]
type = "orient"
markers = ["slider.axis", "rail.axis"]
```

The orient fixes all relative rotation, leaving only translation along the
second marker's $z$-axis. The coordinate options are needed only when the
relative displacement, velocity, acceleration, or rotation variables are
wanted for output, expressions, initial conditions, or state selection.

## Model input

Spatial and planar files both use `.toml` for readable model input and `.simp`
for stored results. They are distinguished by
`[model].dimension = "spatial"` and by their separate model directories and
executables.

The [Spatial TOML Reference](toml-reference.md) describes the currently
accepted tables and coordinate conventions. The working example is
[`free-rotating-body.toml`](../../models/spatial/free-rotating-body.toml).

The Lua assembly example
[`simple-body-surface.lua`](../../models/spatial/simple-body-surface.lua)
builds a vehicle body from shared indexed vertices with separate painted and
glass patches. Run and view it with:

```text
./bin/simp3d models/spatial/simple-body-surface.lua \
    --output results/examples/spatial/simple-body-surface.simp --overwrite
bin/simpView
```

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

The Julia-built pendulum can be run and opened directly in SimpView with:

```text
julia --project=. examples/spatial/julia_api_pendulum.jl
bin/simpView
```

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
