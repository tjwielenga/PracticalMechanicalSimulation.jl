# Spatial TOML Reference

This reference describes only the currently implemented spatial vertical
slice. Fields or element types documented for the planar modeler are not
automatically available in spatial models.

## Model

Every spatial file requires:

```toml
[model]
name = "free_body"
title = "Free spatial body"
dimension = "spatial"
```

`name` and `title` are optional labels. `dimension = "spatial"` is required.

## Analysis and simulation

The analysis table accepts the following fields:

| Field | Required | Default | Accepted values |
| --- | --- | --- | --- |
| `mode` | no | `"automatic"` | `"automatic"`, `"dynamic"`, `"static"`, or `"modal"` |
| `initialization` | no | `"none"` | `"none"` or `"static_equilibrium"` |
| `static_method` | no | `"newton"` | `"newton"` or `"dynamic_relaxation"` |
| `relaxation_duration` | no | `0.5` | positive pseudo-time interval |
| `relaxation_reduction_factor` | no | `0.25` | at least zero and less than one |
| `relaxation_min_cycles` | no | `1` | positive integer |
| `relaxation_max_cycles` | no | `20` | integer at least as large as `relaxation_min_cycles` |
| `relaxation_polish` | no | `true` | `true` or `false` |
| `handoff_acceleration` | no | `0.01` | positive equivalent body-acceleration limit, m/s² |
| `handoff_speed` | no | `0.01` | positive equivalent body-speed limit, m/s |
| `handoff_correction` | no | `0.01` | positive normalized Newton-correction limit |
| `modes` | no | `10` | positive number of modes to retain |
| `frequency_shift_hz` | no | `0.0` | nonnegative frequency around which modes are ordered |
| `modal_tolerance` | no | `1.0e-9` | positive eigenvalue filtering tolerance |

`automatic` currently selects dynamic analysis, as does an explicit
`dynamic`:

```toml
[analysis]
mode = "automatic"       # or "dynamic"
```

`mode = "static"` solves force and moment balance together with the
position-level constraints. Velocities and accelerations are zero, so inertial
and damping terms do not contribute. One output sample gives one equilibrium.
Multiple output times give quasi-static continuation: time-dependent forces
and motion generators are evaluated at each time, and the preceding
equilibrium predicts the next one. A failed continuation interval is divided
automatically.

A dynamic analysis may find a static equilibrium before integration:

```toml
[analysis]
mode = "dynamic"
initialization = "static_equilibrium"
```

The equilibrium positions and orientations become the dynamic starting
configuration. Body and relative velocities declared in the model are then
restored, and the ordinary dynamic initializer solves consistent velocities,
accelerations, reactions, and force variables. The same forces are active in
both stages.

`mode = "modal"` linearizes the complete implicit system at
`simulation.start_time`. The declared configuration and velocity are made
consistent first, but the operating point need not be a static equilibrium.
Use `initialization = "static_equilibrium"` for conventional modes about a
stationary equilibrium. Algebraic definitions, force laws, constraints, and
reactions remain in the sparse eigenproblem and in the stored mode shapes.
`modes` limits the number retained, and `frequency_shift_hz` orders modes near
the requested frequency. End time and output sample count do not affect a
modal calculation.

The default `static_method = "newton"` solves the static equations directly.
Each ordinary sparse factorization supplies an inexpensive reciprocal-condition
estimate. If the static Jacobian is singular or poorly conditioned, the solver
adds the bodies' mass and body-frame inertia to the correction Jacobian with a
one-second pseudo-time scale. It retains that regularization for the rest of
the equilibrium solve. Mass and inertia do not enter the static equations, and
no velocities or accelerations are solved. They only keep neutral correction
directions regular and tend to preserve the model's supplied position and
orientation. The result records which static solutions used this correction.
For a difficult starting configuration,
`static_method = "dynamic_relaxation"` advances the dynamic equations over a
short pseudo-time interval while physical model time is held fixed. Relaxation
uses first-order BDF and reduces the velocity and acceleration variables after
every accepted step. `relaxation_reduction_factor` specifies the cumulative
reduction over one `relaxation_duration`, so the result is independent of how
many adaptive steps DDASSL takes. One complete interval is always performed
before a Newton handoff is considered. The handoff waits until the mass- and
inertia-scaled static imbalance, the corresponding body speeds, and a trial
Newton configuration correction are all small. This avoids large Newton
motions that could abruptly activate a tire, bumper, or another stiff
one-sided force. `relaxation_min_cycles` requires a minimum number of complete
intervals so an unstable starting equilibrium has time to depart.
`relaxation_max_cycles` limits the total number of intervals.
By default, a ready relaxed configuration is polished to the static-equation
tolerance with Newton iteration. Setting `relaxation_polish = false` accepts
the relaxed configuration as soon as all handoff limits are met. This is
useful for diagnosing a difficult model and for systems with one-sided forces
that should not be crossed by a Newton correction. The result is limited by
the handoff criteria rather than the tighter static-equation tolerance.

The default handoff limits are deliberately conservative. The largest
force-over-mass or inertia-scaled torque imbalance must be below
`handoff_acceleration`; the largest translational or
characteristic-length-scaled angular speed must be below `handoff_speed`; and
the trial Newton translation divided by body characteristic length, or trial
rotation in radians, must be below `handoff_correction`. Absolute force and
torque limits are not used because the same numerical values would have very
different meanings for a small linkage and a vehicle.

DDASSL controls its own pseudo-time step during relaxation; the physical
simulation's `maximum_step` does not limit it. Relaxation uses a deliberately
looser integration-error tolerance because its transient path is not a result
of the analysis. Its implicit corrector and constraint equations retain their
normal accuracy. Selected relaxation states and the Newton corrections are
stored in the result for later inspection. They are diagnostic histories; the
final equilibrium remains the result of the static analysis.

The simulation table is optional. Its defaults are shown here:

```toml
[simulation]
start_time = 0.0
end_time = 1.0
output_samples = 201
output_precision = "single"
relative_tolerance = 1.0e-7
absolute_tolerance = 1.0e-9
initial_step = 1.0e-6
maximum_step = 0.01
```

`output_precision` may be `"single"` or `"double"`. The default `"single"`
stores canonical state histories, static-progress states, health snapshots,
and modal shapes as `Float32` to limit file size. All equations, Newton
iterations, Jacobians, and integration calculations remain `Float64`.
`read_result` converts the stored arrays back to `Float64`; saved-result
initialization then recalculates consistent initial conditions. Select
`"double"` when unusually accurate numerical reuse of the stored values is
more important than file size.

### Saved-result initialization

A spatial model can start from a sample in an existing spatial `.simp` result:

```toml
[initial_conditions]
result = "static-assembly.simp"
sample = "static"
include_velocities = false
```

`sample` may be `"first"`, `"last"`, `"static"`, or a one-based stored sample
number. `"static"` selects the first configuration from a dynamic result that
used static-equilibrium initialization. For a standalone static result, use
`"first"`, `"last"`, or a sample number. A relative result filename is resolved
from the model file's directory; an `IO` model uses the current directory.

Matching uses qualified component names, element types, canonical variable
names, and variable kinds. Body positions, normalized Euler parameters, and
available hinge, revolute, or inline coordinates, and named user-defined
differential states are transferred.
`include_velocities = true` also transfers matching linear, angular, and
relative velocities. Pseudo angles, accelerations, reactions, forces, and
other work variables are not transferred. The new model then performs its
ordinary position, velocity, and acceleration consistency calculations.

Explicit values under `initial.impose` in the new model take precedence over
saved values. Motion generators likewise impose their new model laws at the
new start time. This permits a static assembly model and a later dynamic or
modal model to use different forces while sharing named bodies and joints.

## Parameters

Numeric parameters can be named once and used in force expressions:

```toml
[parameters]
k = 30.0
c = 1.0
free = 0.65
```

Parameter names are local to the model file.

## User-defined equation component

An `equation_component` adds named scalar variables and first-order states to
Sim3D's implicit equation system. It is useful for a controller, tire law,
hydraulic element, or other auxiliary model. For example:

```toml
[controller]
type = "equation_component"
inputs = { angle = "pin.theta", omega = "pin.omega" }
parameters = { target_angle = "-45 deg", kp = 12.0, ki = 10.0, kd = 5.0 }
states = { integral_error = { initial = 0.0, scale = 1.0, static = "hold" } }
variables = { angle_error = { initial = 0.0, scale = 1.0 }, torque = { initial = 0.0, scale = 10.0 } }
state_equations = { integral_error = "der(integral_error) = angle_error" }
equations = ["angle_error = target_angle - angle", "torque = kp*angle_error + ki*integral_error - kd*omega"]

[control_torque]
type = "applied_torque"
joint = "pin"
expression = "controller.torque"
```

`inputs` gives short names to existing model variables. Local state and
algebraic-variable names, local parameters, model parameters, and `t` may be
used in the equations. Fully qualified model-variable names also work.
Local parameters may be numbers or quoted degree angles such as `"-45 deg"`;
the latter are converted to radians. In this example the pendulum starts
horizontal and the negative target places it 45 degrees to the right of
downward vertical. The proportional and derivative terms bring it toward
rest; the integral term supplies the torque needed to hold it against gravity.
`state_equations` must contain one equation per state, written either as
`"der(name) = expression"` or just its right-hand expression. `equations`
contains one scalar equality per declared algebraic variable; either side
may contain variables, so algebraic loops are solved with the model rather
than ordered as assignments. The count is checked when loading the model.

`initial` is the starting value for a differential state and an initial guess
for an algebraic variable. `scale` is a positive characteristic magnitude for
numerical scaling. A state may set
`static = "steady"`, which solves its equation with `der(state) = 0` in static
equilibrium, or `static = "hold"` (the default), which retains its entered
value during static analysis. Algebraic equations remain active in both
cases. Differential states are included in BDF error control; algebraic
variables are not selected as mechanical states by QR. Both are saved in the
`.simp` result, and state values can be transferred by saved-result
initialization even when the new model chooses a different static policy.

Existing applied-force, applied-torque, and spanning-force expressions can
reference an equation component's named variables to apply its loads. Lua
assemblies may construct the same `equation_component` table. A rolling tire's
`longitudinal_expression` and `lateral_expression` can likewise use such
variables, so an auxiliary tire law can use the tire's contact measurements
and supply its trial forces. Equations use the same restricted arithmetic and
elementary functions as force expressions; they cannot execute arbitrary
Julia or Lua during the solve. This first version handles continuous scalar
equations and explicit first-order state rates. It does not yet provide
general effort/flow connectors, discrete events, automatic unit checking,
or a Modelica importer. Use SI units consistently in entered equations.

The complete example is
[`controlled-revolute-pendulum.toml`](../../models/spatial/controlled-revolute-pendulum.toml).

## Modeling assemblies

A modeling assembly is a reusable TOML definition that expands into ordinary
bodies, markers, joints, forces, and graphics before the model is allocated.
The simulation equations do not contain a special assembly object.

List the definition files in the model table:

```toml
[model]
dimension = "spatial"
assemblies = ["../../assemblies/spatial/leaf_spring.toml"]
```

The path is relative to the model file. An instance uses the assembly name as
its `type`:

```toml
[rear_left_spring]
type = "leaf_spring"
frame = "chassis"
leaf = "rear_axle"
front_eye = [-0.64, 0.75, -0.07]
leaf_center = [0.0, 0.75, 0.0]
rear_shackle_eye = [0.70, 0.75, 0.07]
rear_frame_eye = [0.71, 0.75, 0.20]
leaf_stiffness = 36532.0
width = 0.076
mount_stiffness = [7.0e6, 7.0e6, 1.0e6]
mount_rotational_stiffness = [1.0e4, 1.0e4, 81.4]
```

`frame` and `leaf` are interfaces declared by this assembly. They name bodies
already present in the enclosing model. Other fields are checked against the
assembly's parameter declarations. Unknown fields, missing required values,
wrong dimensions, and wrong interface types are reported while the model is
read.

Internal names are qualified by the instance name. The example produces
`rear_left_spring.front_link`, `rear_left_spring.rear_link`, and
`rear_left_spring.shackle`. A marker contributed to the existing axle is named
under both owners, such as `rear_axle.rear_left_spring.front_joint`. This keeps
the marker owned by the axle without risking a name collision with another
spring.

Assemblies may instantiate other imported assemblies. The nested names retain
the full instance hierarchy, and recursive assembly cycles are rejected. The
expanded ordinary TOML is stored in the `.simp` result, so model extraction and
rerunning do not require the original assembly definition files.

The supplied Lua leaf-spring assembly, its declarative TOML counterpart, and
its executable example are:

- [`assemblies/spatial/leaf_spring.lua`](../../assemblies/spatial/leaf_spring.lua)
- [`assemblies/spatial/leaf_spring.toml`](../../assemblies/spatial/leaf_spring.toml)
- [`models/spatial/leaf-spring-assembly.lua`](../../models/spatial/leaf-spring-assembly.lua)
- [`models/spatial/leaf-spring-assembly.toml`](../../models/spatial/leaf-spring-assembly.toml)

The assembly-file format and its build-time expressions are described in
[`modeling-assemblies.md`](modeling-assemblies.md).

## Ground

A ground entry provides an owner for fixed spatial markers:

```toml
[ground]
type = "ground"
```

Ground is optional for a model containing only free bodies.

## Rigid body

Each free body requires positive mass and positive-definite inertia:

```toml
[body]
type = "rigid_body"
mass = 2.0
inertia = [0.08, 0.12, 0.16]
position = [0.0, 0.0, 0.0]
velocity = [0.8, 0.2, 4.0]
angular_velocity = [0.0, 0.0, 4.0]
orientation = [0.5235987755982988, 0.0, 0.0, 1.0]
ic_weight_scale = 1.0
```

`center_of_mass` optionally names a marker belonging to the body. All body
marker positions are entered from the body reference-frame origin. The named
marker gives the CM position and the axes in which `inertia` is entered.
`inertia` may be three principal moments or a symmetric positive-definite 3 by
3 matrix about the CM. The reader transforms it into body-reference
components when the CM marker is oriented.

`position` and `velocity` are global components for the body reference-frame
origin. `angular_velocity` is resolved in the body reference frame. The reader
converts position and velocity to the internal CM state. These three motion
fields default to zero. If `center_of_mass` is omitted, the reference origin
and axes are also the CM origin and inertia axes, preserving the original
model convention.

For example, a link may be described naturally from its first joint while its
mass center lies partway along the link:

```toml
[link]
type = "rigid_body"
mass = 4.0
inertia = [0.12, 0.40, 0.44]
center_of_mass = "link.cm"

[link.joint]
type = "marker"

[link.cm]
type = "marker"
position = [0.4, 0.05, 0.0]

[link.tip]
type = "marker"
position = [1.0, 0.0, 0.0]
```

The CM marker is otherwise an ordinary oriented body marker. It may be used
by elements or graphics. It must belong to the body that names it.

The short orientation form is

```text
[angle, axis_x, axis_y, axis_z]
```

The angle is in radians and follows the right-hand rule. The reader normalizes
the axis, which must be nonzero. The rotation axis has the same components in
the frame before and after the rotation because rotation about an axis leaves
that axis unchanged.

A quoted angle ending in a degree symbol or `deg` is converted to radians:

```toml
orientation = ["30°", 0.0, 0.0, 1.0]
orientation = ["30 deg", 0.0, 0.0, 1.0]
```

An orientation may instead be supplied as a proper orthonormal rotation matrix
that maps body components to global components:

```toml
orientation = [
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
]
```

Orientation defaults to the identity matrix. Euler parameters are an internal
integration representation and are not entered in the model.

### Initial-condition assembly

The body fields above are initial guesses. Before selecting states, the loader
corrects body positions, orientations, velocities, and angular velocities so
that the joint equations are consistent.

The translational correction weight is the body's mass. A more massive body
therefore moves less than a lighter body during assembly. For rotation, the
program uses

$$
w_{\mathrm{rotation}}=mL^2,
$$

where $L$ is the body's characteristic length. Thus $mL^2$ acts like a scalar
rotational inertia. It has the correct dimensions but is not the body's full
inertia tensor. The same weights are used for position and velocity
correction.

The characteristic length is normally the greatest distance from the body
origin to one of its markers. If the body has no offset marker, the program
estimates $L$ from its inertia and mass. `ic_weight_scale` is an optional
positive scalar that multiplies both the translational weight $m$ and the
rotational weight $mL^2$. Its default is `1.0`.

An exact body initial condition is placed under `initial.impose`:

```toml
[body.initial.impose]
R_x = 0.0
R_y = 0.0
R_z = 1.0
orientation = ["30 deg", 0.0, 1.0, 0.0]
V_x = 0.5
omega_z = 1.0
```

The available component names are `R_x`, `R_y`, `R_z`, `V_x`, `V_y`, `V_z`,
`omega_x`, `omega_y`, and `omega_z`. Any subset may be imposed. Position and
velocity conditions are independent. Orientation is imposed as one complete
axis-angle or matrix value; individual Euler parameters cannot be imposed.
The explicit `R` and `V` names refer to the canonical CM position and velocity;
the body-level `position` and `velocity` fields above locate the reference
origin and are converted before initial-condition assembly.

An imposed value is removed from the corresponding consistency correction.
The remaining bodies and components must be able to satisfy the model
equations or loading fails with an initial-condition error. Imposition applies
only at the initial time; it does not constrain subsequent motion.

## Marker

A marker must be nested under its owning body or ground:

```toml
[body.tip]
type = "marker"
position = [0.4, 0.0, 0.0]

[ground.reference]
type = "marker"
position = [0.0, 0.0, 0.0]
```

The position of a body marker is expressed from the body reference-frame
origin. The reader subtracts the optional CM offset before assembling the
CM-based equations. The position of a ground marker is expressed globally. A
marker may also have an `orientation` matrix. It maps marker components into
its owner's reference-frame components and defaults to the identity matrix.

## Measurements

Measurements report kinematic quantities without applying a force, adding a
constraint, or changing the model's mobility. They can be plotted from the
stored result like any other canonical variable.

A `span` measures the positive straight-line distance between two ordered
marker points:

```toml
[range]
type = "span"
markers = ["ground.origin", "body.center"]
```

The markers must initially be separated. The element reports `range.distance`,
`range.velocity`, and `range.acceleration`. It also retains the three
separation components and the three components of the instantaneous unit
direction as explicit outputs.

A `directed_distance` measures the signed distance of the first marker from
the plane through the second marker, along the second marker's oriented local
$z$-axis:

```toml
[height]
type = "directed_distance"
markers = ["body.center", "ground.origin"]
```

It reports `height.distance`, `height.velocity`, and `height.acceleration`.
Positive distance places the first marker on the positive-$z$ side of the
second marker's plane. Motion and rotation of either marker are included in
the velocity and acceleration.

Measurement distance and velocity may be used in force expressions, for
example:

```toml
expression = "k*range.distance + c*range.velocity"
```

Acceleration is currently an output, not an allowed force-expression input.
Measurements are not state candidates and do not participate in
initial-condition correction.

The viewer draws a span as a dotted line with arrowheads pointing outward
toward its two markers. A directed distance is a dotted line from the
projected point on the second marker's plane to the measured marker. Its
single arrow always points along the second marker's positive $z$-axis. Thus,
for a negative distance, the arrowhead lies near the projected point and
points away from the measured marker. A small translucent square at the
projected point shows the measurement plane. Arrowheads stop slightly short
of marker symbols so that the markers do not hide them.

The viewer's `Measurements` control shows or hides these symbols. Their color,
opacity, and initial visibility may be overridden like other inferred element
graphics:

```toml
[range.graphics]
color = "darkorange"
opacity = 0.8
visible = true
```

## Applied force

An applied force uses two ordered markers:

```toml
[push]
type = "applied_force"
markers = ["body.application", "ground.force_axis"]
force = 12.0
```

The first marker is the application point. The force direction is the second
marker's oriented local $z$-axis. The direction marker may belong to any body,
including the application or reaction body, or it may belong to ground. A
positive force acts along the marker's positive $z$ direction.

An expression may replace the constant magnitude:

```toml
expression = "f0 * sin(2*pi*t) + gain * body.V_z"
```

Exactly one of `force` or `expression` is required. Expressions use the same
parameters, model variables, mathematical functions, and dual-number
derivatives described for the spanning force below.

`active_during` may restrict the force to one analysis stage or a list of
stages:

```toml
active_during = "static"
```

```toml
active_during = ["dynamic", "modal"]
```

The available stages are `"static"`, `"dynamic"`, and `"modal"`. Omitting
the field, or using `"always"`, activates the force in all three. Static
includes both dynamic relaxation and Newton polish. When inactive, the
element retains its output equations but reports zero force and contributes
nothing to the body balances.

When it is shorter to name the excluded stages, use `inactive_during`
instead:

```toml
inactive_during = "static"
```

This is equivalent to `active_during = ["dynamic", "modal"]`. The two fields
cannot be used together.

`reaction_body` is optional:

```toml
reaction_body = "actuator"
```

When it is present, the element creates a floating marker named
`element-name.reaction`. This marker stays coincident with the application
point but is owned by the reaction body, where it applies the opposite force
and its corresponding moment. The reaction body must differ from the
application body. When `reaction_body` is omitted, ground is assumed and no
reaction is entered in a body balance or drawn by the viewer.

The element contributes the output variables `force`, `F_x`, `F_y`, and
`F_z`. The viewer draws the first force in the applied-force color and an
optional opposite force in the reaction color.

## Applied torque

An applied torque references one hinge or revolute joint:

```toml
[drive]
type = "applied_torque"
joint = "pin"
torque = 5.0
```

A positive torque acts on the first side of the joint around the positive
$z$-axis of its second marker. The second side receives the equal and opposite
torque. A ground-side torque is included physically but hidden by the viewer
unless ground loads are enabled.

Exactly one of three torque definitions is required. In addition to a
constant `torque`, an expression may use `t`, parameters, and supported model
variables:

```toml
expression = "drive_torque - c * pin.omega"
```

Referencing a joint automatically enables its relative `theta`, `omega`, and
`alpha` definitions, so they are available to expressions, output, and state
selection.

The built-in torsional spring-damper form is:

```toml
stiffness = 8.0
damping_time_scale = 0.03
free_angle = "0 deg"
```

It applies

$$
T=-k(\theta-\theta_0)-c\omega.
$$

`stiffness` and the resulting damping must be nonnegative. Specify either
`damping` directly or `damping_time_scale`, for which
$c=k\,\mathtt{damping\_time\_scale}$. If neither is present, damping is zero.
The free angle defaults to zero and accepts radians, a quoted degree value, or
`"initial"`. The last choice records the joint angle after consistent position
initialization, so the spring begins without angular preload.

`active_during` may be `"static"`, `"dynamic"`, `"modal"`, or an array of
those stage names. Omitting it, or specifying `"always"`, keeps the torque
active in every stage. An inactive torque reports zero magnitude and zero
global torque and contributes nothing to either body's balance equations.
`inactive_during` may be used instead to name the excluded stages; the two
fields cannot be combined.

The element contributes `torque`, `T_x`, `T_y`, and `T_z`. The viewer uses a
square-shaft, double-cone symbol aligned with the joint axis. The first-side
torque uses the applied-load color and the second-side torque uses the reaction
color.

## Bushing

A spatial bushing is a six-component spring-damper between two ordered
markers:

```toml
[support]
type = "bushing"
markers = ["body.mount", "ground.mount"]
translational_stiffness = [600.0, 600.0, 900.0]
rotational_stiffness = [30.0, 45.0, 30.0]
damping_time_scale = 0.04
active_during = "always"
```

The first marker receives the bushing force and torque. The second marker
defines the local $x$, $y$, and $z$ directions used by the coefficients. The
opposite force and torque act on the second body at a floating point coincident
with the first marker. A ground-side load is included physically and hidden by
the viewer unless ground loads are enabled.

`translational_stiffness` is required and contains three nonnegative values.
`rotational_stiffness` is optional and defaults to three zeros. Damping may be
specified directly with three-component `translational_damping` and
`rotational_damping` fields. When either is omitted, it is estimated from

$$
c=k\,\mathtt{damping\_time\_scale}.
$$

The time scale is local to this bushing and defaults to zero. An explicit
damping vector overrides the estimate for that part of the bushing.

`active_during` optionally limits the bushing to one analysis stage or an
array of stages. Its values are `"always"` (the default), `"static"`,
`"dynamic"`, and `"modal"`; `"always"` cannot be combined with other values.
Alternatively, `inactive_during` names the excluded stage or stages. For
example, `inactive_during = "static"` leaves a bushing active in dynamic and
modal analyses. The two fields cannot be used together.
A static-only bushing remains active throughout dynamic relaxation and Newton
polish, then becomes unloaded before dynamic velocities, accelerations, and
reactions are initialized. This is useful as a weak numerical hold on neutral
static coordinates without leaving a joint or restraint in the subsequent
simulation. Model time does not distinguish these stages: it is held fixed
while dynamic relaxation advances in pseudo-time.

Marker placement and orientation define the unloaded bushing. The marker
points coincide at zero translational deformation, and their corresponding
axes align at zero rotational deformation. There are no separate free-position
or free-angle fields.

The rotational deformation uses one $z$-$y$-$x$ Bryant-angle set. Large
relative rotation should be around the second marker's local $z$-axis. The
local $x$ and $y$ rotations should remain moderate; the angle set becomes
singular if the $y$ rotation reaches $90$ degrees. The same definition works
normally when all three rotations are small. The periodic angle equations
allow the reported `angle_z` to remain continuous through complete turns.

The bushing contributes no constraint. It leaves all six relative freedoms
and adds the result variables `r_x`, `r_y`, `r_z`, `v_x`, `v_y`, `v_z`,
`angle_x`, `angle_y`, `angle_z`, `omega_x`, `omega_y`, `omega_z`, local force
and torque components `f_*` and `tau_*`, and global components `F_*` and `T_*`.
The viewer draws a connector and the applied and reaction force and torque
symbols.

### Compliant revolute joint

A bushing can be used in place of an ideal revolute joint. Give it high
translational stiffness in all three directions, high rotational stiffness
around the second marker's $x$ and $y$ axes, and zero rotational stiffness
and damping around its $z$ axis:

```toml
[compliant_pin]
type = "bushing"
markers = ["link.pin", "ground.pin"]
translational_stiffness = [6000.0, 6000.0, 9000.0]
rotational_stiffness = [300.0, 450.0, 0.0]
damping_time_scale = 0.04
```

The five stiff directions approximate the five constraints of a revolute
joint, while rotation around $z$ remains free. Unlike an ideal revolute, the
bushing does not remove any degrees of freedom. Its restrained motions have
small finite deflections determined by force balance, and the large
coefficients make the dynamic equations numerically stiff. This can be useful
for representing joint flexibility or for opening a closed constraint loop.
The stiff integrator is intended to handle this kind of model, but coefficients
should be only as large as the required deflection accuracy demands.

## Plane contact

A plane contact applies a one-sided compliant force between a
marker-centered sphere and an oriented plane:

```toml
[floor_contact]
type = "plane_contact"
markers = ["ball.center", "ground.plane"]
radius = 0.1
stiffness = 10000.0
damping_factor = 0.15
transition_depth = 0.001
```

| Field | Required | Default |
| --- | --- | --- |
| `markers` | yes, exactly two | — |
| `radius` | yes, positive | — |
| `stiffness` | exactly one of `stiffness` or `expression` | — |
| `expression` | exactly one of `stiffness` or `expression` | — |
| `damping_factor` | no, nonnegative, s/m | `0.0` |
| `transition_depth` | no, nonnegative length | `0.0` |
| `active_during` | no | all stages |
| `inactive_during` | no | none |

The first marker is the center of a sphere with the specified radius. The
second marker defines an infinite plane through its point, with its local
$z$-axis as the outward normal. Positive `gap` means separation and negative
`gap` means penetration. This is a compliant force, not a constraint.

By default, the viewer draws a dark rubber-colored sphere on the first marker
and a gray circular plane on the second marker. The plane radius is three
times the contact radius. These inferred graphics use the same radius as the
contact equations, so force activation agrees with the visible geometry.

`transition_depth` smooths contact engagement. Over this first amount of
penetration, the tangent stiffness rises continuously from zero to the
specified `stiffness`; deeper contact retains that linear stiffness. A zero
value gives an immediate stiffness change. A small positive value is
recommended for stiff contacts because BDF correctors, like other implicit
integrators, work better without a sharp stiffness corner.

`damping_factor` increases the force during closing motion and reduces it
during rebound without permitting a tensile contact force. `active_during`
and `inactive_during` use the same `static`, `dynamic`, and `modal` stage names
as other spatial force elements. A vehicle jounce or rebound bumper can
therefore specify `inactive_during = "static"` while remaining active during
the subsequent dynamic analysis.

An expression replaces the complete built-in normal-force law. It may use the
contact's own `gap` and `gap_rate`, as well as the other model variables and
functions available to spatial force expressions. For example, this law uses
a one-micrometre activation ramp, a linear stiffness, and a constant viscous
coefficient:

```toml
[parameters]
contact_stiffness = 4.0e6
contact_damping = 2.0e4
activation_depth = 1.0e-6

[stop]
type = "plane_contact"
markers = ["arm.stop", "body.stop"]
radius = 0.03
expression = "min(1,max(0,-stop.gap/activation_depth))*max(0,-contact_stiffness*stop.gap-contact_damping*stop.gap_rate)"
```

The element applies the expression's value exactly. It does not automatically
set an expression force to zero during separation and does not clamp a
negative value. Contact gating, smoothing, and any no-tension rule therefore
belong in the expression. `damping_factor` and `transition_depth` describe
only the built-in stiffness law and cannot be combined with `expression`.

The component outputs `gap`, `gap_rate`, `normal_force`, `F_x`, `F_y`, and
`F_z`. The global force acts on the sphere body along the plane normal. A
non-ground plane body receives the equal-and-opposite force at the same
projected contact point. The force is continuous at contact entry, exit, and
damping-branch changes. The dynamic integrator crosses these points through
its ordinary BDF correction and error control without forcing a history
restart.

See [Plane contact in the Technical
Manual](../../architecture/spatial/spatial-element-formulations.md#plane-contact)
for the equations and load-transfer convention.

## Rolling tire

A rolling tire uses a wheel-center marker and a planar-road marker:

```toml
[tire]
type = "rolling_tire"
markers = ["wheel.center", "ground.road"]
radius = 0.50
regularization_speed = 0.10
normal_stiffness = 15000.0
normal_damping_time_scale = 0.015
longitudinal_expression = "4000*tire.slip_ratio"
lateral_expression = "-40000*tire.slip_angle"
friction_limit = "ellipse"
mu_longitudinal = 0.90
mu_lateral = 0.80
```

The first marker must belong to the wheel body. Its local $z$-axis defines
the axle, but its $x$ and $y$ directions may spin with the wheel. The second
marker defines an infinite road plane through its point with its local
$z$-axis as the outward normal. The wheel axle must not be parallel to the
road normal. The tire does not refer to a joint, so the wheel may be supported
by a revolute, a bushing, or another connection.

`radius` is required and positive. `regularization_speed` is positive and
defaults to 0.1. It keeps the normalized slip quantities finite when the
forward speed approaches zero. The tire also reports unnormalized slip
velocities, which are often more useful in a low-speed force law.

The normal force is zero when the tire is clear of the plane and cannot become
tensile. A linear normal law uses positive `normal_stiffness` and either
nonnegative `normal_damping` or `normal_damping_time_scale`:

$$
F_n=\max(0,k_n\delta+c_n\dot\delta),
\qquad
c_n=k_n\,\mathtt{normal\_damping\_time\_scale}.
$$

If neither damping field is present, normal damping is zero. A nonlinear law
may replace the stiffness fields with `normal_expression`. It may use
`tire.deflection`, `tire.deflection_rate`, model parameters, and the other
ordinary expression inputs. The tire still makes the force zero outside
contact and clips a negative expression value to zero.

`longitudinal_expression` and `lateral_expression` are required. They define
trial forces and may use the tire measurements and `tire.normal_force`. A
positive longitudinal force acts in the tire's calculated forward direction.
A positive lateral force acts toward the axle's positive side projected into
the road plane. These are ordinary restricted expressions, so their partials
are calculated with dual numbers.

The default is **no relaxation**: force expressions use the current slip
without adding tread-deformation states. The Goodyear LT245/75R16 Lua assembly
used by the large van has the same default. To enable relaxation in that
assembly, supply both `longitudinal_relaxation_length` and
`lateral_relaxation_length`; the large-van wrapper exposes them as
`tire_longitudinal_relaxation_length` and
`tire_lateral_relaxation_length`.

`friction_limit` may be `"none"`, which is the default, or `"ellipse"`. The
ellipse requires positive `mu_longitudinal` and `mu_lateral`. When the two
trial forces lie outside

$$
\left(\frac{F_x}{\mu_xF_n}\right)^2+
\left(\frac{F_y}{\mu_yF_n}\right)^2=1,
$$

they are scaled together to lie on the ellipse. This retains their trial-force
direction and couples braking or traction to cornering capacity.

When both relaxation lengths are supplied, the tire adds two independent
first-order internal equations:

```toml
longitudinal_relaxation_length = 0.30
lateral_relaxation_length = 0.45
longitudinal_expression = "tire.normal_force * 12 / 0.30 * tire.longitudinal_deformation"
lateral_expression = "-tire.normal_force * 8 / 0.45 * tire.lateral_deformation"
```

$$
\dot\delta_x=-v_{sx}-\frac{|v_x|}{L_x}\delta_x,
\qquad
\dot\delta_y=v_{sy}-\frac{|v_x|}{L_y}\delta_y.
$$

The implementation uses a smooth approximation to $|v_x|$ near zero. The
deformations are ordinary expression inputs named `longitudinal_deformation`
and `lateral_deformation`. They are integrated and kept under error control,
but they are not candidates in the mechanical QR state selection. At zero
transport speed a loaded contact retains its deformation and can therefore
carry a static tangential force. An unloaded tire releases stored deformation.

Both positive relaxation lengths must be present to use the transient model.
If both are omitted, the tire has no internal deformation states and the
expressions may use `slip_ratio`, `slip_angle`, or slip velocities to define
the earlier instantaneous model.

The deformation states supply forces; they are not kinematic constraints. A
free vehicle can therefore still have neutral longitudinal, lateral, or yaw
motions during a direct static solution. A direct static solve retains the
supplied or transferred deformation because its zero-speed differential
equations do not select a unique value. Dynamic relaxation can develop the
deformation when brakes or other loads oppose tire motion.

The outputs are `deflection`, `deflection_rate`, `forward_velocity`,
`lateral_velocity`, `longitudinal_slip_velocity`,
`lateral_slip_velocity`, `slip_ratio`, `slip_angle`, `camber_angle`,
`normal_force`, `longitudinal_trial_force`, `lateral_trial_force`,
`longitudinal_force`, `lateral_force`, and the global components `F_x`,
`F_y`, and `F_z`. A transient tire also outputs
`longitudinal_deformation` and `lateral_deformation`. The viewer draws the
tire force at the projected contact
point. Its road reaction follows the usual rule for ground-load visibility.

## Translational motion

A translational motion generator prescribes the signed distance of one marker
from the plane through a second marker:

```toml
[drive]
type = "translational_motion"
markers = ["slider.axis", "ground.axis"]
function = "constant_speed"
initial_distance = 0.2
velocity = "speed"
```

The ordered markers define

$$
d=(P_i-P_j)^T\hat z_j,
$$

where $\hat z_j$ is the second marker's oriented local $z$-axis. The second
marker therefore defines both the reference plane and the positive direction.
This is one primitive constraint: it does not constrain relative motion within
the plane or any relative rotation. If the second marker moves or rotates, the
velocity and acceleration equations include the corresponding motion of its
point and $z$-axis.

For `function = "constant_speed"`, `velocity` is required and may be a number
or a numeric parameter/expression evaluated at the initial time.
`initial_distance` defaults to zero.

A general prescribed history supplies one distance expression:

```toml
[drive]
type = "translational_motion"
markers = ["slider.axis", "ground.axis"]
function = "expression"
distance = "0.2 + 0.1*t^2"
```

The program uses automatic differentiation to obtain velocity and acceleration
from the distance expression. The expression must therefore be twice
differentiable over the simulation interval. Omitting `function` selects
`expression`. Do not supply separate `velocity` or `acceleration` fields with
an expression.

The component stores `distance_m`, `velocity_m`, `acceleration_m`, and
`force_m`. A positive drive force acts on the first marker along $\hat z_j$.
The second marker's body receives the equal-and-opposite force and its
corresponding moment; a ground-side reaction is hidden by the viewer unless
ground loads are enabled.

## Spanning motion

A spanning motion generator prescribes the positive distance between two
marker points:

$$
\ell=\lVert P_2-P_1\rVert>0.
$$

A constant-distance generator is an ideal massless rigid link with spherical
ends:

```toml
[link]
type = "spanning_motion"
markers = ["body.link", "ground.anchor"]
distance = 1.0
```

The supplied distance must be positive. The two initial marker points must be
separated because their initial line establishes the direction used during
consistent position initialization.

The `distance` field may be a positive number or a time expression. The program
automatically differentiates it to obtain the prescribed velocity and
acceleration. The `constant_speed` form remains available and uses
`initial_distance` and `velocity`, just like the translational generator. The
prescribed distance is checked whenever the equations are evaluated and the
simulation stops if it becomes zero or negative.

The velocity and acceleration definitions are

$$
\dot\ell=\hat u^T(v_2-v_1),
$$

$$
\ddot\ell=\hat u^T(a_2-a_1)+
\frac{(v_2-v_1)^T(v_2-v_1)-\dot\ell^2}{\ell},
$$

where $\hat u=(P_2-P_1)/\ell$. The second term in the acceleration equation
accounts for motion perpendicular to the link.

The component's internal span stores `s_x`, `s_y`, `s_z`, `distance_m`,
`u_x`, `u_y`, `u_z`, `velocity_m`, and `acceleration_m`. The generator also
stores `force_m` and the global components `F_x`, `F_y`, and `F_z`. The stored
reaction is the force on the first marker along the line from the second marker
toward the first. Positive `force_m` is compression and negative `force_m` is
tension. The viewer draws the link and both endpoint forces, hiding a
ground-side force by default.

## Rotational motion

A rotational motion generator prescribes the relative rotation of a hinge or
revolute joint:

```toml
[drive]
type = "rotational_motion"
joint = "pin"
function = "constant_speed"
initial_angle = "0 deg"
angular_velocity = "speed"
```

The joint is required. Referencing it automatically enables its relative
`theta`, `omega`, and `alpha` variables. Positive motion rotates the first side
of the joint around the positive $z$-axis of its second marker. The generator
adds one drive-torque variable around that axis and prescribes the relative
angle, angular velocity, and angular acceleration at their corresponding
equation levels. It therefore removes the joint's remaining rotational
freedom.

For `function = "constant_speed"`, `angular_velocity` is required and may be a
number or a numeric parameter/expression evaluated at the initial time.
`initial_angle` defaults to zero and accepts radians or a quoted degree value.

A general prescribed history supplies one angle expression:

```toml
[drive]
type = "rotational_motion"
joint = "pin"
function = "expression"
angle = "0.25*sin(4*t)"
```

The program obtains angular velocity and acceleration by automatically
differentiating the angle expression, which must be twice differentiable over
the simulation interval. Omitting `function` selects `expression`. Do not
supply separate `angular_velocity` or `angular_acceleration` fields. A joint
may have only one rotational motion, and a driven joint cannot also have an
`[initial]` table. The drive torque is stored as `element-name.torque`. The
viewer draws it in the applied-load color on the first side of the joint and
draws its opposite in the reaction color on the second side. A ground-side
reaction is hidden unless ground loads are enabled.

## Spanning force

A spanning force acts along the straight line between two ordered markers:

```toml
[spring]
type = "spanning_force"
markers = ["body.spring_point", "ground.anchor"]
stiffness = 30.0
damping = 1.0
free_length = 0.65
```

The scalar force is the force on the first marker along the line from the
second marker toward the first. Positive force is compression and repels the
markers; negative force is tension and attracts them. The second marker
receives the opposite force. A body marker also receives the moment caused by
applying the force away from its body origin. A ground-side force is not added
to the equations and its viewer arrow is hidden by default.

Exactly one of three force definitions is required. A constant force uses:

```toml
force = 5.0
```

The linear spring-damper form requires all three fields `stiffness`, `damping`,
and `free_length`. Their values must be nonnegative. Its force is

$$
f=-k(\ell-\ell_0)-c\dot\ell.
$$

An expression can use `t`, numeric values from `[parameters]`, supported body
or relative coordinates and velocities, and the spanning force's own
`length` and `length_rate`:

```toml
expression = "-k * (spring.length - free) - c * spring.length_rate"
```

The available mathematical functions and expression restrictions are the
same as in the planar modeler. Expressions cannot assign values or call Julia
code. Their partial derivatives are evaluated with dual numbers for the
analytical system Jacobian.

`active_during` has the same meaning as for an applied force. It may select
one or more of `"static"`, `"dynamic"`, and `"modal"`; omission or `"always"`
selects all three. When inactive, the span measurements remain available but
the scalar and global force outputs are zero and no force enters either body
balance. `inactive_during` may instead name one or more excluded stages.

The element contributes the output variables `s_x`, `s_y`, `s_z`, `length`,
`u_x`, `u_y`, `u_z`, `length_rate`, `length_acceleration`, `force`, `F_x`,
`F_y`, and `F_z`. The `s_*` and `u_*` vectors point from the first marker to the
second; the `F_*` components are the force on the first marker and therefore
equal `-force * u_*`. The initial marker points must not coincide. The viewer
draws the connector and both force arrows, subject to the normal rule that
ground arrows are hidden.

## Pulley and belt

A pulley identifies the pitch circle carried by the first side of a spatial
revolute joint:

```toml
[driver_pulley]
type = "pulley"
body = "driver"
joint = "driver_joint"
pitch_radius = 0.4
```

The named body must own the revolute joint's first marker. Referencing the
joint automatically enables its continuous rotation coordinate. The pitch
radius must be positive.

A belt is an ordered closed loop of two or more straight tangent spans:

```toml
[belt]
type = "belt"
spans = ["belt.upper", "belt.return"]
stiffness = 50000.0
damping_time_scale = 0.001
initial_tension = 40.0

[belt.upper]
type = "belt_span"
pulleys = ["driver_pulley", "driven_pulley"]
near_points = ["ground.driver_upper", "ground.driven_upper"]
```

Each pair of near-point markers selects the intended common tangent at the
initial configuration. The pulley order must continue from one span to the
next and close back on itself. Stiffness, damping, and damping time scale may
be set on the belt and overridden on an individual span. Specify either a
nonnegative `initial_tension` or a positive `free_length`, but not both.

Pulley axes may be parallel, antiparallel, or angled. For angled axes, the
locations, axes, and pitch radii must admit an exact straight line tangent to
both pitch circles. The model is rejected when no such tangent exists. The
optional positive `feasibility_tolerance`, normally `1.0e-5`, controls this
geometric check.

The current belt is massless, elastic, bilateral, and has no slip. Each span
stores its tangent points, direction, length, extension, extension rate,
tension, and global force. The viewer draws the straight spans, pulley wraps,
and the action and reaction force arrows. Pulley body graphics remain explicit
body graphics; they are not generated by the belt element.

## Gravity

Gravity names one or more bodies and uses global acceleration components:

```toml
[gravity]
type = "gravity"
acceleration = [0.0, 0.0, -9.81]
bodies = ["body"]
```

Multiple gravity entries are allowed and their accelerations add on a body.

## Spherical joint

A spherical joint makes two marker points coincident while leaving their three
relative rotations free:

```toml
[pin]
type = "spherical"
markers = ["pendulum.pin", "ground.pin"]
```

`markers` is required and ordered. At position, velocity, and acceleration
levels the first marker point equals the second marker point. The joint has a
three-component global reaction. Positive reaction components are the force on
the first marker; the opposite force acts on the second marker. Marker
orientations have no effect on a spherical joint.

The loader makes supplied body positions, orientations, and velocities
consistent before it selects states, then solves the initial acceleration and
reaction equations. Position assembly applies finite body-frame rotation
corrections, so the resulting orientation remains a proper rotation matrix.
After assembly, pivoted QR identifies any redundant scalar ideal constraints.
The program deactivates each redundant reaction and its complete position,
velocity, and acceleration equation family before selecting states.

## Perpendicular-axis constraint

The `perp` primitive makes the $x$-axis of its first marker perpendicular to
the $y$-axis of its second marker:

```toml
[guide]
type = "perp"
markers = ["pendulum.pin", "ground.pin"]
```

`markers` is required and ordered. The constraint removes one relative
rotation. Its scalar reaction is a torque on the first marker in the
instantaneous direction

$$
\hat n=\hat x_i\times\hat y_j.
$$

The opposite torque acts on the second marker. The direction changes with
either marker. When the constraint is satisfied, the two marker axes are unit
and perpendicular, so $\hat n$ is also a unit vector.

## Inplane constraint

The `inplane` primitive keeps the first marker point in the plane defined by
the second marker:

```toml
[support]
type = "inplane"
markers = ["slider.contact", "ground.plane"]
```

`markers` is required and ordered. The plane passes through the second marker
and its normal is the second marker's oriented local $z$-axis. The element
contributes one scalar constraint family at position, velocity, and
acceleration levels. Its scalar reaction acts along the plane normal at the
first marker point. The opposite force acts at the coincident floating point
on the second marker's body; no opposite force is stored when that marker is
on ground.

## Inline constraint

The `inline` primitive keeps the first marker point on the line through the
second marker along its oriented local $z$-axis:

```toml
[guide]
type = "inline"
markers = ["slider.center", "ground.axis"]
translation_coordinates = true
```

Internally it applies two inplane constraints with normals $\hat x_j$ and
$\hat y_j$. It therefore contributes two transverse reaction components but
no reaction along $\hat z_j$. It does not constrain relative orientation.

With `translation_coordinates = true`, the element adds the signed relative
variables `guide.distance`, `guide.velocity`, and `guide.acceleration` along
$+\hat z_j$. `guide.velocity` is eligible for automatic or preferred state
selection and the coordinates can be used by coordinate couplers. The default
is `false` unless another element references the coordinate.

The optional relative initial conditions use scalar values and weights:

```toml
[guide.initial]
distance = 0.2
velocity = 0.5
distance_weight = 10.0
velocity_weight = 2.0
```

For an exact value instead:

```toml
[guide.initial.impose]
distance = 0.2
```

Do not specify the same quantity as both a weighted guess and an imposed
value. Distance and velocity are independent. An `initial` table requires
`translation_coordinates = true`. Without an initial specification, the
relative variables simply follow the assembled body motion and do not
influence it.

## Coordinate coupler

A `coupler` imposes one linear relation among two or more scalar relative
coordinates:

```toml
[screw_coupler]
type = "coupler"
coordinates = ["nut_guide.distance", "screw_bearing.rotation"]
coefficients = [1.0, -0.08]
offset = "initial"
```

Each coordinate name ends in `.rotation` for a `hinge` or `revolute`, or
`.distance` for an `inline`. Referencing a coordinate automatically enables
its angle, distance, velocity, and acceleration variables. The example above
imposes

$$
q_{\mathrm{nut}}-0.08\theta_{\mathrm{screw}}=q_0.
$$

`coordinates` must contain at least two unique ports. `coefficients` supplies
one finite number per port and at least two coefficients must be nonzero. In a
mixed rotation-translation relation, the coefficients carry the units needed
to make every term consistent.

The default `offset = "initial"` evaluates the left side from the supplied
initial coordinates, so adding the coupler does not change their initial
relation. A finite numeric offset may instead establish an absolute relation.
The same coefficients couple relative velocities and accelerations. The
element contributes one scalar reaction, `screw_coupler.lambda`, which maps
back to force or torque at each referenced coordinate.

## Ideal gear pair

An ideal spatial gear pair references two revolute joints and one contact
marker fixed to the carrier:

```toml
[contact]
type = "gear_pair"
joints = ["gear1_joint", "gear2_joint"]
contact_marker = "carrier.contact"
phase = "initial"
```

The joint order identifies the first and second gear. Each revolute must list
its gear marker first and its base marker second. At least one of those base
markers must have the same owner as the contact marker. The joint axes may be
parallel, as in spur and internal gears, or nonparallel, as in a bevel pair.

The program calculates the pitch tangent for each gear from its joint axis,
joint center, and the carrier contact point. The contact point must not lie on
either joint axis. The two calculated tangents must be parallel or
antiparallel; otherwise the supplied geometry cannot describe one ideal point
of contact and the model is rejected. The error message reports the angular
mismatch in degrees.

The element calculates signed effective pitch radii from the same geometry.
It then constrains the two gear rotations at position, velocity, and
acceleration levels. `phase = "initial"`, which is the default, retains the
supplied initial gear orientations. A numeric phase instead specifies an
absolute angular registration in radians.

The element creates floating markers named `contact.contact_1` and
`contact.contact_2`. They belong to the corresponding gear bodies while
remaining coincident with `carrier.contact`. The scalar reaction applies
opposite forces there along the common calculated tangent. The viewer draws
the force on the first gear as the blue action arrow and the equal and opposite
force on the second gear as the red reaction arrow. Gear graphics are specified
separately on their bodies; the gear-pair element does not infer or alter body
graphics. This first implementation represents zero pressure angle; a
separately specified helical or spiral-bevel contact direction is not yet
included.

## Rack and pinion

An ideal spatial spur rack and pinion references an inline constraint and a
revolute joint:

```toml
[rack_and_pinion]
type = "rack_and_pinion"
joints = ["rack_guide", "pinion_bearing"]
pitch_radius = 0.4
phase = "initial"
```

`joints` is required and ordered. The inline constraint must list its
rack-body marker first, and the revolute must list its pinion-body marker
first. Their second markers must belong to the same carrier. Referencing the
elements automatically enables the inline distance and continuous revolute
angle coordinates.

The inline base marker is the rack-and-pinion pitch frame. Its $x$-axis must
be parallel or antiparallel to the pinion revolute axis, its $y$-axis points
from the pinion center toward the rack contact, and its $z$-axis defines rack
travel and contact-force direction. The required `pitch_radius` is positive.
The optional positive `alignment_tolerance` is an angular tolerance in
radians and defaults to $10^{-5}$.

The default `phase = "initial"` preserves the supplied rack position and
pinion angle. A finite numeric phase instead establishes an absolute
registration in the model's length units. The same relation is imposed on
the relative velocities and accelerations.

The element generates `rack_and_pinion.carrier_contact`,
`rack_and_pinion.rack_contact`, and `rack_and_pinion.pinion_contact`. The last
two are body-owned floating markers that remain coincident with the
carrier-fixed pitch contact. The scalar reaction applies a force along the
inline base marker's $z$-axis to the rack and the opposite force to the pinion.
The viewer shows these as the blue action and red reaction arrows. Rack and
pinion body graphics remain explicit and are not generated by the mechanical
element.

This element deliberately represents an ordinary zero-pressure-angle spur
rack and pinion. Skew racks and helical tooth geometry are not included.

## Hinge constraint

The `hinge` orientation constraint leaves relative rotation only around the
$z$-axis of its second marker:

```toml
[axis]
type = "hinge"
markers = ["pendulum.pin", "ground.pin"]
```

`markers` is required and ordered. Internally, the hinge applies two
perpendicular-axis constraints:

$$
\hat x_i\mathbin{\cdot}\hat z_j=0,
\qquad
\hat y_i\mathbin{\cdot}\hat z_j=0.
$$

The corresponding reaction-torque directions are
$\hat x_i\times\hat z_j$ and $\hat y_i\times\hat z_j$. Together they can
transmit torque perpendicular to the hinge axis but not around it.

A hinge does not make its marker points coincident. Combining a hinge with a
spherical joint at the same markers gives a spatial revolute joint with one
relative rotational freedom.

Set `rotation_coordinates = true` to add the continuous relative variables
`theta`, `omega`, and `alpha`:

```toml
[axis]
type = "hinge"
markers = ["pendulum.pin", "ground.pin"]
rotation_coordinates = true
```

Positive `theta` rotates the first marker's $x$-axis from the second marker's
$x$-axis according to the right-hand rule around the second marker's $z$-axis.
Its initial value is found from the supplied marker orientations. Unlike an
angle calculated only for plotting, `theta` is continuous through complete
revolutions. `omega` and `alpha` are its relative angular velocity and angular
acceleration.

The relative angular velocity also becomes a state-selection candidate. It
can be preferred with:

```toml
[state_selection]
method = "preferred"
preferred_velocities = ["axis.omega"]
allow_fallback = false
```

The list may contain all of the required states or only the choices that matter
to the analyst. For a partial list, pivoted QR completes the state set while
retaining every named velocity. The loader checks that the requested
velocities can coexist in an independent state set. If they cannot,
`allow_fallback = true` permits a wholly automatic selection; `false` reports
the incompatible preference as an input error.

A preferred selection establishes the initial states. During simulation the
solver may replace it with an automatic QR selection if the current state set
causes a singular iteration matrix, repeated corrector failures, or an
unhealthy physical predictor error. A change does not restart the model from
its initial conditions.

`rotation_coordinates` defaults to `false`, so these additional variables and
equations are only allocated when they are useful to the model.

Relative rotation accepts the same initial-condition forms:

```toml
[axis.initial]
angle = "5 deg"
omega = 1.0
angle_weight = 10.0
omega_weight = 2.0
```

For an exact value instead:

```toml
[axis.initial.impose]
angle = "5 deg"
```

Angle and angular velocity are independent. An `initial` table requires
`rotation_coordinates = true`. A specified value is a weighted guess unless
it appears under `initial.impose`. With no initial specification, the relative
variables follow the assembled body orientations and angular velocities
without influencing them.

## Orient constraint

The `orient` constraint fixes all three relative rotations between two
markers while leaving their positions independent:

```toml
[orientation_lock]
type = "orient"
markers = ["body.frame", "ground.reference"]
```

It combines the two perpendicular-axis constraints of a hinge with a third
constraint:

$$
\hat x_i\mathbin{\cdot}\hat z_j=0,
\qquad
\hat y_i\mathbin{\cdot}\hat z_j=0,
\qquad
\hat x_i\mathbin{\cdot}\hat y_j=0.
$$

These are three fixed rotations and provide three independent reaction-torque
directions near the aligned configuration. Consistent initial-condition
assembly moves the body orientations to the nearest configuration satisfying
them and removes inconsistent relative angular velocity.

An `orient` does not constrain either marker position. Combining it with a
spherical joint at the same markers fixes both translation and rotation.

## Revolute joint

A `revolute` is the convenient combination of a spherical joint and a hinge
constraint acting between the same two markers:

```toml
[pin]
type = "revolute"
markers = ["pendulum.pin", "ground.pin"]
rotation_coordinates = true
```

The spherical part makes the marker points coincident. The hinge part aligns
their $z$-axes and leaves rotation only around the second marker's $z$-axis.
The joint therefore has one relative degree of freedom, a three-component
force reaction, and two torque reactions perpendicular to its free axis.

`rotation_coordinates` has the same optional meaning and sign convention as
it does for a hinge. When enabled, this example provides `pin.theta`,
`pin.omega`, and `pin.alpha`, and `pin.omega` is available for state
selection. It defaults to `false`.

## Fixed joint

A `fixed` joint combines a spherical joint and an orient constraint at the
same two markers:

```toml
[weld]
type = "fixed"
markers = ["base.joint", "arm.joint"]
```

The spherical part makes the marker points coincident. The orient part fixes
all three relative rotations. The joint therefore removes all six relative
degrees of freedom and carries a three-component force reaction together with
three independent torque reactions.

Marker locations and orientations establish the fixed relative pose between
the two bodies. The body frames do not need to have the same orientation; only
the two joint-marker frames are assembled into coincidence and alignment.

## Graphics

The shared viewer understands the common top-level graphics settings. A named
`xy_frame` graphic owned by ground draws red, green, and blue arrows along the
positive $x$, $y$, and $z$ directions. It also draws a small translucent patch
of the positive $x$-$y$ plane:

```toml
[ground]
type = "ground"

[ground.origin]
type = "marker"

[ground.origin.graphics]
shape = "xy_frame"
axis_length = 0.35
plane_size = 0.16
plane_color = "gray65"
label = "ground.origin"
```

The graphic must be placed beneath a marker. That marker automatically
supplies its position and orientation; there is no separate `marker` field.
The optional `label` is drawn beside the positive $z$ arrow and can identify
particular frames in a crowded model.

`axis_length` defaults to `0.35`. `plane_size` defaults to 45 percent of the
axis length, `plane_color` defaults to `"gray65"`, and `opacity` defaults to
`0.18`. Setting `visible = false` omits the graphic. For example, a body-marker
frame can be added with:

```toml
[body.tip]
type = "marker"
position = [0.4, 0.0, 0.0]

[body.tip.graphics]
shape = "xy_frame"
axis_length = 0.18
plane_size = 0.07
plane_color = "lightskyblue"
opacity = 0.16
label = "body.tip"
```

The initial spatial viewer also recognizes a body-attached box:

```toml
[graphics]
background = "white"
body_palette = "colorblind"

[body.graphics]
shape = "box"
size = [0.8, 0.4, 0.2]
color = "steelblue"
```

Its `size` values are the marker- or body-reference-frame lengths along $x$,
$y$, and $z$. With no `marker` field, a box or cylinder is centered and
oriented on the body reference frame. An optional marker belonging to the body
can position and orient either shape elsewhere:

```toml
[body.graphics]
shape = "cylinder"
marker = "body.cm"
axis = "z"
radius = 0.08
length = 1.0
```

A physical sphere is specified by its marker and radius:

```toml
[ball.graphics]
shape = "sphere"
marker = "ball.center"
radius = 0.1
color = "steelblue"
```

This is distinct from the automatically displayed inertia ellipsoid. The
inertia ellipsoid shows relative mass and principal-inertia proportions and
is deliberately scaled to the mechanism; it does not define physical body or
contact geometry.

A gear graphic is owned by a body and positioned by a marker belonging to
that body:

```toml
[gear.pitch]
type = "marker"

[gear.graphics]
shape = "gear"
marker = "gear.pitch"
pitch_radius = 0.4
width = 0.06
cone_height = 0.2
color = "steelblue"
```

The marker origin is the center of the pitch circle and its oriented local
$z$-axis is the gear axis. `pitch_radius` and `width` must be positive. If
`cone_height` is omitted, the present viewer draws a cylindrical gear. If it
is supplied, the viewer draws a conical frustum whose virtual apex is

$$
P_{\mathrm{apex}}=P_{\mathrm{marker}}-
    h\hat z_{\mathrm{marker}}.
$$

The pitch radius is measured in the marker plane at the axial center of the
gear. The two face radii are

$$
r_- = r_p\left(1-\frac{w}{2h}\right),\qquad
r_+ = r_p\left(1+\frac{w}{2h}\right).
$$

Consequently, `cone_height` must be greater than half the width. The smooth
cylinder or frustum is the initial representation of the `gear` shape; teeth
can be added later without changing the model description.

The `pulley` shape uses the same marker, `pitch_radius`, `width`, and color
fields but does not accept `cone_height`. It currently draws the pulley as a
smooth cylinder centered on its pitch circle:

```toml
[pulley.graphics]
shape = "pulley"
marker = "pulley.center"
pitch_radius = 0.4
width = 0.12
color = "steelblue"
```

A general conical frustum uses the same marker convention but specifies its
two face radii directly:

```toml
[body.graphics]
shape = "frustum"
marker = "body.frame"
length = 0.3
radius_1 = 0.12
radius_2 = 0.20
color = "darkorange"
```

`radius_1` is at the face displaced by `-length / 2` along the marker's local
$z$-axis. `radius_2` is at the face displaced by `+length / 2`. Both radii and
the length must be positive. Unlike `gear`, `frustum` has no pitch-radius or
cone-apex interpretation.

### User-defined surfaces

A `surface` graphic uses a shared array of vertices and one or more named
patches. Vertices are given in the body reference frame. If `marker` is
specified, they are instead given in that marker's frame. Vertex indices in
the faces are one-based:

```toml
[body.graphics.bodywork]
shape = "surface"
vertices = [
    [-1.0, -0.5, 0.0],
    [ 1.0, -0.5, 0.0],
    [ 1.0,  0.5, 0.0],
    [-1.0,  0.5, 0.0],
    [-0.6, -0.4, 0.7],
    [ 0.6, -0.4, 0.7],
    [ 0.6,  0.4, 0.7],
    [-0.6,  0.4, 0.7]
]
draw_edges = true
edge_color = "gray25"
edge_width = 1.0

[body.graphics.bodywork.patches.paint]
color = "steelblue"
opacity = 1.0
faces = [
    [1, 2, 3, 4],
    [5, 8, 7, 6]
]

[body.graphics.bodywork.patches.glass]
color = "gray15"
opacity = 0.5
faces = [
    [1, 5, 6, 2],
    [4, 3, 7, 8]
]
```

Each face is a planar convex polygon containing at least three vertex indices
listed consecutively around its boundary. Their winding determines the outward
normal used for shading. The viewer divides polygons into triangles without
adding diagonals to the optional edge drawing. Every patch may specify its own
`color` and `opacity`; omitted values
come from the surface and then its owning body. A surface with only one
material may put `faces`, `color`, and `opacity` directly in the surface table
instead of defining `patches`.

Surfaces may be attached directly to a rigid body, to ground, or beneath a
marker. A body or ground surface may name a `marker` explicitly. A surface
placed beneath a marker uses that marker automatically and must not name a
second one. `draw_edges` defaults to `false`, `edge_color` to `"gray25"`, and
`edge_width` to `1.0`. `include_in_fit = false` keeps a large background
surface, such as a road or terrain patch, from determining the initial **Fit
view** limits. The surface remains visible and may still be scaled or hidden.
