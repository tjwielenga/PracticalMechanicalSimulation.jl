# Planar TOML User's Guide

Status: supported executable format

This is the user-facing reference for the current planar model reader. It
explains what can be entered in a TOML model and how each element behaves. The
[Technical Manual](../../architecture/README.md) gives the detailed equations and
program design.

## File structure and naming

A model consists of reserved configuration tables and named element tables.
Every element table has a `type`. Its TOML path is its qualified name:

```toml
[crank]
type = "rigid_body"

[crank.end]
type = "marker"
position = [0.25, 0.0]
```

Here `crank.end` is contained by `crank`. Containment gives the marker its
owner; it does not make a mechanical connection. Joints, forces, and drivers
refer to qualified marker names explicitly.

The reserved top-level tables are `model`, `parameters`, `analysis`,
`simulation`, `initial_conditions`, `state_selection`, and `graphics`. Element
names may otherwise be chosen freely. Elements are sorted by qualified name
when the canonical system is allocated, so file order does not determine
variable indices.

All current values use a consistent user-chosen unit system. The examples use
SI units, radians, and seconds. A quoted scalar angle ending in `°` or `deg`
is converted from degrees when the model is read.

## Element index

Each supported element is summarized below. This table is a quick map
from a modeling purpose to its TOML `type`.

| Purpose | Type | Main references |
| --- | --- | --- |
| Fixed reference | `ground` | owned markers and graphics |
| Planar inertial or algebraic body | `rigid_body` | mass, inertia, initial motion |
| Two-node floating-reference flexible member | `flexible_beam` | section properties and generated end markers |
| Body- or ground-fixed frame | `marker` | local position and orientation |
| Body-owned point following another marker | `floating_marker` | owner and `follows` |
| Pin connection | `revolute` | two markers; optional rotation coordinate |
| One normal translation constraint | `inplane` | two ordered markers |
| One relative-rotation constraint | `perp` | two ordered markers |
| One-freedom prismatic joint | `translational` | inplane and perp primitives |
| Rigid connection | `fixed` | revolute and perp primitives |
| Measured translation | `distance_coordinate` | two markers; state candidate |
| Marker-to-marker measurement | `span` | two markers; reaction-free output |
| Linear coordinate relation | `coupler` | coordinate ports and coefficients |
| Ideal gears | `gear_pair` | two revolutes and carrier contact marker |
| Ideal rack and pinion | `rack_and_pinion` | translational and revolute joints |
| Belt pitch surface | `pulley` | body, revolute joint, pitch radius |
| Closed elastic belt | `belt` | ordered belt spans |
| Tangent elastic segment | `belt_span` | pulleys and near-point markers |
| Uniform body force | `gravity` | body list and acceleration vector |
| Marker-directed force | `applied_force` | application and direction markers |
| Applied marker torque | `applied_torque` | two markers; torque or expression |
| Rotational elastic element | `torsional_spring_damper` | two markers and coefficients |
| Marker-to-marker axial force | `spanning_force` | two points and a scalar force law |
| Planar six-coefficient support | `bushing` | two markers and local coefficients |
| Compliant sphere-plane contact | `plane_contact` | sphere and plane markers |
| Smooth closed planar profile | `curve` | marker and local profile points |
| Circular roller-profile contact | `curve_contact` | curve and roller marker |
| Flat follower-profile contact | `flat_follower_contact` | curve and face marker |
| Tangential contact friction | `surface_friction` | sphere-plane contact |
| Revolute bearing friction | `revolute_friction` | revolute joint and effective radius |
| Translational guide friction | `translational_friction` | translational joint |
| Bilateral plane-guide friction | `inplane_friction` | standalone inplane constraint |
| Prescribed translation | `translational_motion` | two markers and motion law |
| Prescribed rotation | `rotational_motion` | two markers and motion law |
| Auxiliary equations and states | `equation_component` | inputs, variables, states, and equations |

For an end-to-end operating sequence, see
[Using the Planar Modeler](using-planar-modeler.md).

## Configuration tables

### `model`

| Field | Required | Default | Meaning |
| --- | --- | --- | --- |
| `dimension` | yes | — | Must be `"planar"`. |
| `name` | no | — | Short model identifier; also used as the title fallback. |
| `title` | no | `name`, then `"Planar model"` | Human-readable result title. |

### `parameters`

Each entry used by an expression must be numeric. Parameters can appear in
force and torque expressions, restricted prescribed-motion expressions, and
the `angular_velocity` of a `constant_speed` rotational motion:

```toml
[parameters]
speed = 1.2566370614359172
forcing_frequency = 4.0
```

### User-defined equation component

An `equation_component` adds named scalar algebraic variables and first-order
states to the same implicit equation system as the mechanism. It can represent
a controller, actuator law, hydraulic state, or another auxiliary model. For
example:

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
markers = ["pendulum.pin", "ground.pin"]
expression = "controller.torque"
```

`inputs` assigns short local names to existing model variables. Equations may
also use local states, algebraic variables, local or model parameters, `t`, and
fully qualified model-variable names. A local parameter may be numeric or a
quoted degree angle. Each state needs one `state_equations` entry, written as
either `"der(name) = expression"` or only the rate expression. The number of
algebraic `equations` must equal the number of declared `variables`. These are
implicit equalities, not ordered assignments, so algebraic loops are allowed.

`initial` supplies a differential-state starting value or an algebraic initial
guess. `scale` is a positive characteristic magnitude used for numerical
scaling. `static = "steady"` solves a state's rate equation with its derivative
set to zero during static equilibrium. `static = "hold"`, the default, retains
the entered state value during statics. Algebraic equations remain active in
both cases. Differential states participate in BDF error control and are not
mechanical candidates for QR state selection. User states and variables are
stored by name in the `.simp` result; saved-result initialization transfers
matching state values and recalculates the algebraic variables.

Applied-force, applied-torque, and spanning-force expressions may reference
the component's qualified variables. Equations use the same restricted
arithmetic and elementary functions as force expressions. They cannot execute
arbitrary Julia code during the solve. The current component supports
continuous scalar equations and explicit first-order state rates; it does not
yet provide discrete events, unit checking, or effort/flow connectors.

The complete example is
[`controlled-revolute-pendulum.toml`](../../models/planar/controlled-revolute-pendulum.toml).

### `graphics`

Graphics are optional presentation data. They contribute no mass, inertia,
force, constraint, variable, or equation. The simulation preserves them in the
TOML embedded in a `.simp` result, and the stored-result viewer reconstructs
them from the saved marker histories.

The top-level table selects the viewer background and automatic body palette:

```toml
[graphics]
background = "white"
body_palette = "colorblind"

[graphics.colors]
moving_body = "steelblue"
second_body = "darkorange"
support = "dimgray"
```

`body_palette` accepts `"default"`, `"colorblind"`, or a nonempty array of
colors. A color can be a name from `graphics.colors`, a standard CSS/SVG color
name, or a CSS color specification such as `"#4E79A7"`,
`"rgb(78,121,167)"`, or `"hsl(210,40%,48%)"`. Model-specific color names are
resolved before standard names.

Every mechanical element may have an optional appearance table:

```toml
[spring.graphics]
color = "second_body"
opacity = 0.8
visible = true
show_default = true
```

| Field | Default | Meaning |
| --- | --- | --- |
| `color` | viewer default | Named or explicit color for inferred graphics and contained primitives. |
| `opacity` | viewer default | Number from zero through one. |
| `visible` | `true` | Display this element and its explicit primitives. |
| `show_default` | `true` | Display the viewer's inferred representation of this element. |

The optional load-graphics table sets the initial force and torque display:

```toml
[graphics.loads]
show_reactions = true
show_applied_loads = true
show_torques = true
show_ground_loads = false
reaction_color = "gold2"
applied_color = "darkorange2"
```

The colors may use names defined in `graphics.colors`. The viewer provides
the same visibility controls while it is running, together with force
and torque size controls. All force arrows use one common linear scale. A
force arrow has a round shaft and one conical head. A torque symbol has a
square shaft and two consecutive cones pointing in the same direction. The
torque-size control scales the complete symbol without changing its shape.
Cone area is proportional to torque magnitude. Zero and unrenderably small
forces or torques are not drawn.

Gold reaction graphics and orange applied-load graphics are the defaults. For
an ideal connection whose equal-and-opposite forces act at one point, only
the force on the first marker or body is shown. For every two-ended applied
force or torque element, the load on the first marker uses the applied-load
color and the equal-and-opposite load on the second marker uses the reaction
color. This includes spanning springs, torsional springs, bushings, contact,
belts, and explicitly applied forces and torques.

An element's own `graphics.color` controls its body, connector, spring, or
contact symbol; it does not override this load-side color convention. Change
`applied_color` and `reaction_color` in `[graphics.loads]` when a different
load color pair is wanted.

Forces and torques acting on ground are hidden by default and are excluded
from the initial automatic load scale. Set `show_ground_loads = true` or use
the viewer's **Ground loads** checkbox to display them. This visibility setting
does not remove the ground load from the mechanical equations.

Named graphical primitives may be nested under `graphics` on a `rigid_body` or
`ground`. They use `shape`, not the mechanical `type` field:

```toml
[crank.graphics]
show_default = false
color = "moving_body"

[crank.graphics.member]
shape = "cylinder"
markers = ["crank.center", "crank.end"]
radius = 0.04

[crank.graphics.hub]
shape = "cylinder"
marker = "crank.center"
axis = "z"
length = 0.08
radius = 0.07

[crank.graphics.pin]
shape = "sphere"
marker = "crank.end"
radius = 0.055
```

The first cylinder form spans two markers. The second is centered at one
marker and uses its local `"x"`, `"y"`, or `"z"` axis; the default axis is
`"z"`. The marker supplies the primitive's position and orientation.

Marker-centered boxes and ellipsoids use full three-dimensional sizes even in
a planar model, allowing a useful out-of-plane thickness:

```toml
[body.graphics.mass]
shape = "ellipsoid"
marker = "body.cm"
size = [0.40, 0.12, 0.10]
opacity = 0.5

[ground.graphics.base]
shape = "box"
marker = "ground.origin"
size = [0.20, 0.12, 0.08]
```

Each referenced marker must belong to the rigid body or ground element that
contains the primitive. A primitive may override its parent's `color`,
`opacity`, or `visible` setting. Explicit graphics describe appearance only;
contact and other physical geometry must still be declared through modeling
elements.

### `analysis`

| Field | Required | Default | Accepted values |
| --- | --- | --- | --- |
| `mode` | no | `"automatic"` | `"automatic"`, `"dynamic"`, `"kinematic"`, `"static"`, or `"modal"` |
| `initialization` | no | `"none"` | `"none"` or `"static_equilibrium"` |
| `static_method` | no | `"newton"` | `"newton"` or `"dynamic_relaxation"` |
| `relaxation_duration` | no | `0.5` | positive pseudo-time interval |
| `relaxation_reduction_factor` | no | `0.25` | at least zero and less than one |
| `relaxation_cycles` | no | `8` | positive integer |
| `modes` | no | `10` | positive number of modes to retain |
| `frequency_shift_hz` | no | `0.0` | nonnegative frequency around which modes are ordered |
| `modal_tolerance` | no | `1.0e-9` | positive eigenvalue filtering and pairing tolerance |

With `automatic`, a model with selected independent states is dynamic; a
zero-state fully driven mechanism is kinematic. `dynamic` selects the same
implicit integration path explicitly, including for a zero-state mechanism.
`kinematic` requires zero independent states. `static` solves equilibrium with
all velocities and accelerations zero. Multiple static output times give a
quasi-static continuation through time-dependent forces and prescribed motion.

For a dynamic model, `initialization = "static_equilibrium"` performs one
static solution at `simulation.start_time` before integration. Static positions,
orientations, and relative positions become the starting configuration. The
velocities declared by the bodies and motion generators are then restored, and
the ordinary simultaneous dynamic initializer solves consistent accelerations,
reactions, and applied-force variables. Static reaction and force values serve
only as Newton initial guesses. Forces with stage-dependent activation are
switched between the static and dynamic solves. The option is rejected when
the resolved analysis is kinematic or static.

`static_method = "dynamic_relaxation"` applies to either a standalone static
analysis or static initialization of dynamics. The program starts the
relaxation from rest, advances the existing implicit dynamic equations through
`relaxation_duration` of pseudo-time while holding model time fixed, and
multiplies the velocity and acceleration variables by
`relaxation_reduction_factor` between cycles. After every cycle it attempts a
short static Newton correction. That correction supplies the final accurate
equilibrium; relaxation only
brings a difficult starting configuration into its convergence region. At most
`relaxation_cycles` pseudo-time intervals are attempted. These controls are
ignored by the default `newton` method.

`modal` linearizes the complete sparse implicit equation set at
`simulation.start_time`. The declared initial position and velocity are first
made dynamically consistent, so the operating point must satisfy the
constraints but need not be an equilibrium. Algebraic definitions, force laws,
and reactions remain in the eigenproblem and in the complex mode shapes.
`initialization = "static_equilibrium"` instead calculates conventional modes
about a stationary equilibrium. `modes` limits the retained eigenpairs, while
`frequency_shift_hz` controls their ordering. One representative of an
oscillatory conjugate pair is retained; distinct real decay roots are retained
individually. Modal analysis requires at least one selected state. The
simulation end time and output sample count do not affect a modal calculation.

The modeling program currently implements only the `StateSelected` dynamic
formulation, whose constraint-derivative deficit is zero. Deficit is derived
internally and is not a TOML option. Alternative deficit and stabilization
formulations remain executable paper studies but are not selectable through
the model file.

### `simulation`

| Field | Default | Validation |
| --- | ---: | --- |
| `start_time` | `0.0` | numeric model time |
| `end_time` | `5.0` | not less than `start_time` |
| `output_samples` | `81` | integer, at least 1 |
| `output_precision` | `"single"` | `"single"` or `"double"` |
| `relative_tolerance` | `1.0e-7` | positive |
| `absolute_tolerance` | `1.0e-9` | positive |
| `initial_step` | `1.0e-6` | positive |
| `maximum_step` | `0.005` | positive |

Only requested output times are retained in a `.simp` result; internal solver
steps do not control result-file size. `output_precision = "single"` stores
canonical state histories and modal shapes as `Float32`, although the analysis
itself always uses `Float64`. The reader converts stored values back to
`Float64`, and saved-result initialization performs its ordinary consistency
calculation. Use `"double"` when the saved values will be used for sensitive
numerical postprocessing or a particularly accurate restart.

### `initial_conditions`

This optional table starts a model from a sample in an existing `.simp`
result:

```toml
[initial_conditions]
result = "assembled-vehicle.simp"
sample = "static"
include_velocities = false
```

| Field | Required | Default | Meaning |
| --- | --- | --- | --- |
| `result` | yes | — | Source `.simp` filename. |
| `sample` | no | `"last"` | `"first"`, `"last"`, `"static"`, or a one-based stored sample number. |
| `include_velocities` | no | `false` | Also transfer matching physical and relative velocity variables. |

A relative result filename is resolved from the directory containing the TOML
file. For an `IO` model, including standard input, it is resolved from the
current working directory.

`sample = "static"` is specifically for a dynamic result whose analysis began
with `initialization = "static_equilibrium"`. It selects the first saved
configuration and verifies from the result diagnostics that static
initialization actually occurred. The first sample contains that equilibrium
configuration even though its reactions, forces, and accelerations have
already been recalculated for dynamics. Using `include_velocities = false`
therefore extracts the static position without importing the subsequently
assigned dynamic velocities. A source that did not perform static
initialization is rejected for this selector.

The loader compares qualified component names and element types, then compares
canonical variable names and kinds. Matching positions, orientations, and
relative positions are transferred. When `include_velocities = true`, matching
linear, angular, and relative velocities are transferred as well. Components
that exist in only one model are left alone, which permits the second model to
add or remove forces and other elements. A same-named component with a changed
type is rejected rather than copied accidentally. At least one configuration
variable must match.

Accelerations, reactions, applied-force variables, and other work variables
are never transferred. They are recalculated by the ordinary simultaneous
initializer. Motion-generator values also follow the laws in the new model at
its start time. The saved configuration is applied before position correction
and automatic state selection. Consequently, when velocities are included,
the velocities chosen as states in the new model retain their transferred
values while dependent velocities are made consistent.

This supports a two-file sequence in which the first model performs assembly
or static equilibrium and the second changes its forces and runs dynamics. If
the second model also requests `analysis.initialization = "static_equilibrium"`,
the transferred configuration is the starting estimate for that new static
solution.

### `state_selection`

The table is normally omitted. The program corrects the supplied configuration
onto the position constraints, forms the scaled velocity-constraint partial
matrix, and uses column-pivoted QR to choose individual physical velocities.

```toml
[state_selection]
method = "preferred"
preferred_velocities = ["crank.omega"]
allow_fallback = false
maximum_condition_number = 6.7108864e7
```

| Field | Default | Meaning |
| --- | --- | --- |
| `method` | `"automatic"` | Use `"preferred"` to request named physical velocities. |
| `preferred_velocities` | `[]` | Qualified body `V_x`, `V_y`, or `omega`; revolute-coordinate `omega`; or distance-coordinate `velocity` names. Required for `preferred`. |
| `allow_fallback` | `true` | Use automatic QR selection if the preferred dependent block is unusable. |
| `maximum_condition_number` | $1/\sqrt{\epsilon}$ | Largest accepted condition number of the scaled dependent block. |

The number of preferred velocities must equal the mechanism mobility. Position
partners or body-fixed pseudo angles are inferred. Redundant scalar constraint
families are detected separately by QR of the transpose and are suppressed from
the active solve while retaining inactive canonical result entries.

Preferred velocities determine the initial partition. During dynamic analysis,
the recovery logic may replace them if physical predictor errors become much
larger than the selected-state error, the iteration matrix is singular, or
corrector failures repeat. The new selection is recorded in the `.simp` result.

## Geometric and structural elements

This guide describes the fields an analyst supplies. The corresponding
equations and assembly details are in the [Planar Element Formulations chapter
of the Technical Manual](../../architecture/planar/planar-element-formulations.md).

### Ground

```toml
[ground]
type = "ground"
```

Ground owns no variables. An ordinary marker must be nested beneath a ground
or rigid-body element. A flexible beam supplies its own end markers.

### Rigid body

```toml
[link]
type = "rigid_body"
mass = 2.0
inertia = 0.16666666666666666
position = [1.0, 0.0]
orientation = 0.0
velocity = [0.0, 0.0]
angular_velocity = 0.0
```

| Field | Required | Default | Meaning |
| --- | --- | --- | --- |
| `mass` | yes | — | Nonnegative planar body mass. Zero denotes no translational inertia. |
| `inertia` | yes | — | Nonnegative center-of-mass moment about the out-of-plane axis. Zero denotes no rotational inertia. |
| `position` | no | `[0.0, 0.0]` | Initial global center-of-mass position. |
| `orientation` | no | `0.0` | Initial body orientation about the positive $z$-axis. Numbers are radians; quoted values may use `"45°"` or `"45 deg"`. |
| `velocity` | no | `[0.0, 0.0]` | Initial global translational velocity. |
| `angular_velocity` | no | `0.0` | Initial angular velocity. |
| `characteristic_length` | no | inferred | Positive length used for rotational IC weighting and automatic state-selection scaling. |
| `ic_weight_scale` | no | `1.0` | Positive scalar multiplying the body's mass-based IC weight. |

The supplied positions and velocities are initial guesses. Before selecting
states, the program makes them consistent with the position- and
velocity-level equations. Translation uses the body weight

$$
w_{\mathrm{translation}}=s\,m,
$$

and rotation uses

$$
w_{\mathrm{rotation}}=s\,m\,L^2,
$$

where $m$ is body mass, $L$ is its characteristic length, and $s$ is
`ic_weight_scale`. Thus $mL^2$ acts like a scalar rotational inertia. It has
the correct dimensions but is not necessarily the body's specified moment of
inertia. The same two weights are used for position and velocity correction.
A larger `ic_weight_scale` makes the entire body resist correction more
strongly.

The characteristic length is normally the greatest distance from the body
origin to one of its markers. If no marker supplies a positive length, it is
estimated from inertia and mass. A positive `characteristic_length` overrides
the inferred value.

Individual body values can instead be imposed exactly during initialization:

```toml
[link.initial.impose]
R_x = 1.0
theta = "30 deg"
V_x = 0.5
omega = 2.0
```

The accepted names are `R_x`, `R_y`, `theta`, `V_x`, `V_y`, and `omega`.
Imposed values are temporary initial conditions, not motion generators. They
are held while the remaining coordinates are corrected. An incompatible set
of imposed values causes initialization to stop with an error.

When dynamic analysis requests static-equilibrium initialization, the static
solve may change imposed positions as it establishes equilibrium. Supplied
velocities are applied afterward and the weighted velocity projection is
repeated at the equilibrium configuration, so imposed velocities remain exact
at the beginning of dynamics.

The inferred characteristic length normally comes from the body's marker
locations. The program supplies reasonable fallbacks when the body has no
noncentral marker.

Zero mass and inertia are supported for intermediate bodies. The balance
equations for a zero-inertia direction are algebraic. Models must still supply
enough connections or force laws to determine that motion; otherwise
initialization or integration reports a singular system.

For IC correction only, a zero-mass body receives a small positive effective
mass equal to $10^{-6}$ times the smallest positive body mass in the model. If
the model has no positive mass, one model mass unit supplies the reference.
This avoids a singular numerical weight while allowing the massless body to
move readily during assembly. `ic_weight_scale` can increase or decrease that
default resistance.

Automatic state selection still follows kinematic mobility, so a velocity on a
massless body may be retained as a state coordinate. This does not add inertia.
An ordinary dynamic run must begin at a force-consistent algebraic position;
`initialization = "static_equilibrium"` can find one when needed.

### Flexible beam

```toml
[beam]
type = "flexible_beam"
length = 1.0
area = 0.01
second_moment = 8.333333333333333e-6
elastic_modulus = 2.0e7
shear_modulus = 8.0e6
mass = 1.0
position = [0.5, 0.0]
damping_time_scale = 0.002
```

The planar flexible beam is a two-node Timoshenko member carried by a floating
reference frame. Its reference center can translate and rotate through large
motions while its three elastic coordinates remain small. It owns mass and
center-of-mass inertia like a rigid body and adds axial, transverse, and
relative-rotation deformation.

| Field | Required | Default | Meaning |
| --- | --- | --- | --- |
| `length` | yes | — | Undeformed end-to-end length. |
| `area` | yes | — | Cross-sectional area. |
| `second_moment` | yes | — | Area second moment for in-plane bending. |
| `elastic_modulus` | yes | — | Young's modulus $E$. |
| `shear_modulus` | yes | — | Shear modulus $G$. |
| `shear_coefficient` | no | $5/6$ | Effective shear-area coefficient. |
| `mass` | yes | — | Positive total beam mass. |
| `inertia` | no | $mL^2/12$ | Reference-frame center-of-mass inertia. |
| `damping_time_scale` | no | `0.0` | Multiplies the elastic stiffness matrix to obtain viscous damping. |
| `position`, `orientation` | no | zero | Initial pose of the undeformed center frame. |
| `velocity`, `angular_velocity` | no | zero | Initial reference-frame motion. |
| `elastic_position` | no | `[0, 0, 0]` | Initial axial, transverse, and relative-rotation elastic coordinates. |
| `elastic_velocity` | no | `[0, 0, 0]` | Initial elastic-coordinate rates. |

Every beam automatically creates `beam.end_i`, `beam.cm`, and `beam.end_j`.
The center marker follows the floating reference frame; the end markers also
include elastic translation and rotation. They may be used by ordinary joints,
constraints, forces, and graphics. End loads are projected into the beam's
three rigid and three elastic balance equations. The elastic velocities also
participate in automatic state selection, so a fixed cantilever naturally
selects only its three elastic states.

See [`flexible-cantilever.toml`](../../models/planar/flexible-cantilever.toml)
for a complete static example.

### Marker

```toml
[link.end]
type = "marker"
position = [0.5, 0.0]
orientation = 0.0
```

Both fields are optional and default to zero. On a body they are expressed in
the body frame. On ground they are global. A marker combines position and
orientation; its orientation is available even when a particular joint uses
only its point. Because planar rotation is always about the positive $z$-axis,
the orientation is one scalar angle. A number is in radians; a quoted value
such as `"45°"` or `"45 deg"` is converted from degrees. The earlier field
name `angle` remains accepted for existing models, but a table must not specify
both names.

### Floating marker

```toml
[gear.contact]
type = "floating_marker"
follows = "carrier.contact"
orientation = 0.0
```

A floating marker must be nested beneath a rigid body. Its global point follows
the position and velocity of the ordinary marker named by `follows`, while its
orientation and force ownership remain with the containing body. Its optional
`orientation` is expressed in the containing body frame and defaults to zero. The
current use is a gear contact point that is located by a carrier but transmits
force to a gear. A floating marker cannot be used as an endpoint of a revolute
or in-plane joint.

## Ideal connections

The complete position, velocity, and acceleration equations and their reaction
conventions are given under [Ideal connections in the Technical
Manual](../../architecture/planar/planar-element-formulations.md#ideal-connections).

### Revolute joint

```toml
[pin]
type = "revolute"
markers = ["crank.end", "rod.end"]
rotation_coordinates = true
```

`markers` is required and must name exactly two markers. Their points coincide.
The joint contributes two reaction variables and two scalar constraint families
at position, velocity, and acceleration levels. Marker orientations remain
free. `rotation_coordinates` is optional and defaults to `false`. When true,
the joint also owns continuous `theta`, `omega`, and `alpha` variables defined
as the first marker's orientation, angular velocity, and angular acceleration
relative to the second marker. These are ordinary stored outputs, and `omega`
also becomes a state-selection candidate paired with `theta` and `alpha`.

To prefer this relative coordinate as the independent state, use its qualified
velocity name:

```toml
[state_selection]
method = "preferred"
preferred_velocities = ["pin.omega"]
allow_fallback = false
```

The relative angle and angular velocity normally come from the supplied body
and marker values. Weighted relative guesses may instead be supplied with an
`initial` table:

```toml
[pin.initial]
angle = "5 deg"
omega = 10.0
angle_weight = 100.0
omega_weight = 10.0
```

Larger weights resist correction. Exact relative initial conditions use the
nested `impose` table:

```toml
[pin.initial.impose]
angle = "5 deg"
omega = 10.0
```

Weighted and imposed fields may be mixed, but the same field must not appear
in both tables. A revolute initial table requires
`rotation_coordinates = true`.

### Span measurement

```toml
[range]
type = "span"
markers = ["ground.origin", "body.center"]
```

A span reports the positive straight-line distance between two ordered marker
points. It owns the separation components `s_x` and `s_y`, `distance`, unit
direction components `u_x` and `u_y`, `velocity`, and `acceleration`. The
velocity and acceleration are the first and second time derivatives of the
distance.

The span is a reaction-free measurement. It does not constrain either marker
and is not a state-selection candidate. Its `distance` and `velocity` may be
used in force expressions as `range.distance` and `range.velocity`. Because it
cannot become a replacement state, its predictor error does not trigger the
state-health monitor. The two markers must not initially coincide, and
floating markers are not currently accepted.

When **Measurements** is enabled, the viewer shows a dotted line between the
markers with arrowheads pointing outward toward both ends.

### Distance coordinate

```toml
[slider_distance]
type = "distance_coordinate"
markers = ["slider.center", "ground.measurement_axis"]
```

A distance coordinate measures, but does not prescribe or constrain, the
signed separation of the first marker from the second marker along the second
marker's local $y$-axis. It owns `distance`, `velocity`, and `acceleration`
variables. These are available as stored outputs, and `velocity` is an
automatic state-selection candidate paired with `distance` and `acceleration`.
For example, `slider_distance.velocity` may be named in
`state_selection.preferred_velocities`.

Initial distance and velocity values use the same weighted and imposed forms
as revolute coordinates:

```toml
[slider_distance.initial]
velocity = 1.0
velocity_weight = 10.0

[slider_distance.initial.impose]
distance = 0.25
```

In this example distance is imposed exactly and therefore must be omitted from
the parent `initial` table; velocity remains a weighted guess.

The coordinate contributes no force by itself. It is intended for output,
state selection, and use by a coordinate coupler. Floating markers are not
accepted. When **Measurements** is enabled, the viewer shows the signed
distance with a dotted line, a small plane at the projection point, and an
arrow pointing along the second marker's positive local $y$-axis.

### Coordinate coupler

```toml
[shaft_coupler]
type = "coupler"
coordinates = ["input_pin.rotation", "output_pin.rotation"]
coefficients = [1.0, -2.0]
offset = "initial"
```

A coupler supplies one ideal linear relation among two or more coordinate
ports. `coordinates` and `coefficients` are ordered arrays of equal length. A
rotation port is available as `joint.rotation` when that revolute joint has
`rotation_coordinates = true`. A distance port is available as
`coordinate.distance` on a `distance_coordinate` element. Rotation and
distance ports may be mixed, with the coefficients supplying any required
unit conversion. For example, a screw with travel $p$ per radian uses

```toml
coordinates = ["nut_travel.distance", "screw_bearing.rotation"]
coefficients = [1.0, -0.05]
```

to impose $x-0.05\theta=b$.

`offset` may be numeric. Its default value, `"initial"`, evaluates the weighted
coordinate sum from the supplied initial configuration, allowing the coupler
to preserve the starting phase without requiring the user to calculate it.
Use a numeric zero when absolute phasing is wanted.

The coupler owns one reaction `lambda`. The viewer shows the resulting force or
torque at each coupled coordinate using the reaction color.

One coupler always supplies one relation. Coupling $N$ coordinates completely
to a master therefore normally requires $N-1$ independent couplers.

### Inplane constraint

```toml
[slider]
type = "inplane"
markers = ["rod.slider_end", "ground.guide"]
```

`markers` is required and ordered. The first marker point is constrained to the
plane defined by the second marker. The plane passes through the second marker
and its normal is the second marker's oriented local $y$-axis. The element
contributes one reaction along that normal and one scalar constraint family at
all three levels.

### Perp constraint

```toml
[guide_orientation]
type = "perp"
markers = ["slider.center", "ground.guide"]
```

`markers` is required and ordered. The constraint makes the first marker's
local $x$-axis perpendicular to the second marker's local $y$-axis, thereby
preventing relative rotation. It owns one equal-and-opposite reaction torque.
A `perp` and an `inplane` constraint between the same markers form the two
primitive constraints of a planar prismatic joint.

### Translational joint

```toml
[guide]
type = "translational"
markers = ["slider.center", "ground.guide"]
```

`markers` is required and ordered. A translational joint is a convenience
composition of an `inplane` primitive and a `perp` primitive. The second
marker's local $y$-axis defines the constrained direction, so relative
translation remains free along its local $x$-axis. Relative rotation is
constrained.

The generated components are named `guide.inplane` and `guide.perp` in this
example. They retain their individual scalar reaction force and reaction
torque in stored results. Floating markers are not accepted.

### Fixed joint

```toml
[rigid_connection]
type = "fixed"
markers = ["body_1.end", "body_2.end"]
```

`markers` is required and ordered. A fixed joint is a convenience composition,
not a separate constraint formulation. The loader constructs a revolute
primitive and a perp primitive between the same markers. The revolute removes
two relative translations and the perp removes relative rotation.

The primitive reactions remain individually visible in stored results. For the
example above they are named `rigid_connection.revolute.lambda_x`,
`rigid_connection.revolute.lambda_y`, and
`rigid_connection.perp.lambda`. Floating markers are not accepted because the
revolute primitive requires body-fixed or ground-fixed points.

### Ideal gear pair

```toml
[gear_pair]
type = "gear_pair"
joints = ["left_bearing", "right_bearing"]
contact_marker = "carrier.contact"
phase = "initial"
```

The ordered `joints` must be two revolute joints whose first marker belongs to
the corresponding gear body. `contact_marker` denotes the location of the
contact point between the two gears. It must be fixed to the base body of at
least one of the two revolute joints.

The two joint centers and contact point must be collinear, and the distances
from the centers to the contact point become the pitch radii. A contact point
between the centers makes an external pair whose gears rotate in opposite
directions. A contact point outside the centers makes an internal pair whose
gears rotate in the same direction.

The gear pair creates `gear_pair.contact_1` and `gear_pair.contact_2` as
body-owned floating markers that follow this contact point. `phase` may be a
number or `"initial"`. It defaults to `"initial"`, which preserves the initial
gear angles. A number specifies the phase directly.

The contact marker's body is the geometric carrier. It need not support both
joints, which allows a planet-to-carrier joint to be paired with a
sun-to-ground or ring-to-ground joint. The supporting joints must maintain the
required pitch geometry as the mechanism moves.

The gear-pair reaction force acts on the floating marker belonging to the
first gear named by `joints[1]`, along the calculated tangent. The
equal-and-opposite force acts on the floating marker belonging to the second
gear named by `joints[2]`. The forces at their coincident contact point produce
the corresponding gear moments and bearing loads. The viewer's gear-contact
arrow shows the force on the first gear. The contact marker's `orientation` has no
effect on the gear pair.

See [Ideal gear pair in the Technical
Manual](../../architecture/planar/planar-element-formulations.md#ideal-gear-pair) for the
signed-radius constraint and reaction-force construction.

### Rack and pinion

```toml
[rack_and_pinion]
type = "rack_and_pinion"
joints = ["rack_guide", "pinion_bearing"]
pitch_radius = 0.4
phase = "initial"
```

`joints` is required and ordered. Its first entry must be a `translational`
joint whose first marker belongs to the rack and whose second marker belongs
to the common carrier. Its second entry must be a `revolute` joint connecting
the pinion to that carrier. `pitch_radius` is required and positive. `phase`
may be a number or `"initial"`. It defaults to `"initial"`, which preserves the
supplied rack position and pinion angle. A numeric phase is an equivalent rack
displacement in the model's length units.

The carrier-side translational marker's $x$-axis defines rack motion and force.
Its $y$-axis selects which side of the pinion contacts the rack. The pitch
radius locates that contact point from the carrier-side revolute marker.

The element generates three markers:
`rack_and_pinion.carrier_contact`, `rack_and_pinion.rack_contact`, and
`rack_and_pinion.pinion_contact`. The latter two are body-owned floating
markers following the carrier-fixed contact point. Equal-and-opposite forces
act there along the carrier-side translational marker's $x$-axis, retaining the
correct moment arms on both bodies.

See [Rack and pinion in the Technical
Manual](../../architecture/planar/planar-element-formulations.md#rack-and-pinion) for the
constraint and force construction.

### Pulley and elastic belt

A pulley names its body, supporting revolute joint, and pitch radius:

```toml
[driver_pulley]
type = "pulley"
body = "driver"
joint = "driver_bearing"
pitch_radius = 0.35
```

The named body must be attached to the revolute joint. A pulley nested beneath
its body may omit `body`. Supporting joints need not share a common carrier.

An ordered parent belt owns two or more tangent spans:

```toml
[belt]
type = "belt"
spans = ["belt.driver_to_driven", "belt.driven_to_tensioner",
         "belt.tensioner_to_driver"]
stiffness = 100000.0
damping_time_scale = 0.002
initial_tension = 100.0

[belt.driver_to_driven]
type = "belt_span"
pulleys = ["driver_pulley", "driven_pulley"]
near_points = ["carrier.driver_upper_hint", "carrier.driven_upper_hint"]
```

The ordered spans must form a closed loop: the second pulley of each span is
the first pulley of the next span. The two ordinary `near_points` select one of
the four feasible common tangents at assembly. The selected external or
internal branch remains fixed as the mechanism moves. The loader rejects an
infeasible tangent, ambiguous hints, or adjacent spans that imply incompatible
wrap directions.

`stiffness`, `damping`, and `damping_time_scale` may be specified on the belt
and overridden on a span. Stiffness is required and positive. When `damping`
is omitted, $c=k\tau_d$. A belt may specify either nonnegative
`initial_tension` or positive `free_length`, but not both; omitted preload
defaults to zero. A free length longer than the assembled pitch path produces
negative initial tension and a warning because the present force law is
bilateral.

Each span stores its tangent points, length, extension, extension rate, tension,
and force. Pulley rotation transfers belt material between adjacent spans, and
the assembled pulley angles are used as the extension reference. Positive
tension pulls the pulleys together at their tangent points.

The current belt is massless, no-slip, and bilateral. Negative tension is
reported rather than replaced by slack-belt behavior. See [Pulley, belt, and
belt span in the Technical
Manual](../../architecture/planar/planar-element-formulations.md#pulley-belt-and-belt-span)
for the tangent and constitutive formulation.

## Forces

“Force” is used broadly here and includes torques. Forces affect balances and
reactions but do not by themselves reduce mobility.

The force laws and their balance contributions are described under [Force
elements in the Technical
Manual](../../architecture/planar/planar-element-formulations.md#force-elements).

Applied forces, applied torques, spanning forces, bushings, plane contacts,
and both curve-contact followers may be limited to selected analysis stages.
The available stages are
`"static"`, `"dynamic"`, and `"modal"`:

```toml
active_during = "static"
```

```toml
active_during = ["dynamic", "modal"]
```

Omitting `active_during`, or setting it to `"always"`, activates the force in
all three stages. `inactive_during` may instead name the excluded stage or
stages. For example,

```toml
inactive_during = "static"
```

is equivalent to `active_during = ["dynamic", "modal"]`. The two fields
cannot be used together. Static includes dynamic relaxation and Newton polish.
An inactive element retains its kinematic geometry and equations but reports
zero force or torque and contributes nothing to body balances. SimpView hides
its load and connector graphics during that stage. This makes a static-only
assembly support or a dynamic-only bumper explicit in one model file. See
[`stage-dependent-force-drop.toml`](../../models/planar/stage-dependent-force-drop.toml).
When a dynamic result includes static initialization, SimpView presents the
saved static snapshots and dynamic time history as separate choices under
**Analysis**.

### Gravity

```toml
[gravity]
type = "gravity"
acceleration = [0.0, -9.81]
bodies = ["crank", "rod"]
```

`acceleration` is required. `bodies` defaults to all rigid bodies. Each listed
name must identify a body.

### Applied force

```toml
[push]
type = "applied_force"
markers = ["link.tip", "actuator.axis"]
force = 100.0
reaction_body = "actuator"
```

The force is applied at the first marker in the direction of the second
marker's local y axis. The direction marker may belong to any body, including
the application body or reaction body, or it may belong to ground. A positive
`force` acts in the positive y direction of that marker.

`reaction_body` is optional. When it is present, the element creates a
floating marker named `element-name.reaction`. This marker remains coincident
with the application point, is owned by the reaction body, and receives the
equal-and-opposite force. When `reaction_body` is omitted, ground is assumed
and no reaction is entered in a body balance. The reaction body must differ
from the application body.

An expression may be used in place of a constant magnitude:

```toml
[push]
type = "applied_force"
markers = ["link.tip", "actuator.axis"]
expression = "100.0*sin(2*pi*t)"
```

Exactly one of `force` or `expression` is required. In addition to time and
parameters, the expression may use the named model coordinates described for
applied torque below.

### Applied torque

```toml
[drive]
type = "applied_torque"
markers = ["crank.ground_end", "ground.origin"]
torque = 10.0
```

The two marker orientations are required. Positive torque acts on the first
marker and the equal-and-opposite torque acts on the second. An expression may
be used in place of the constant magnitude:

```toml
[drive]
type = "applied_torque"
markers = ["crank.ground_end", "ground.origin"]
expression = "-k*pin.theta-c*pin.omega"
```

Exactly one of `torque` or `expression` is required.

The expression may use numeric parameters, numeric literals, `t`, `pi`, `e`,
arithmetic `+`, `-`, `*`, `/`, and `^`, and `sin`, `cos`, `tan`, `asin`,
`acos`, `atan`, `tanh`, `exp`, `log`, `sqrt`, `abs`, `min`, and `max`. It may
also refer directly to
these qualified kinematic variables:

| Source | Available names |
| --- | --- |
| Rigid body | `body.R_x`, `body.R_y`, `body.theta`, `body.V_x`, `body.V_y`, `body.omega` |
| Revolute with `rotation_coordinates = true` | `joint.theta`, `joint.omega` |
| `distance_coordinate` | `coordinate.distance`, `coordinate.velocity` |
| `span` | `measure.distance`, `measure.velocity` |
| `spanning_force` | `element.length`, `element.length_rate` |

Other registered position and velocity coordinates, including prescribed
motion coordinates, follow the names shown in the result-variable catalog.
Accelerations, reactions, and applied-load variables are not accepted in a
load expression.

A state-dependent expression contributes an explicit local force or torque
variable and the implicit equation

$$
L-f(q,\dot q,t)=0.
$$

The load is consequently stored in the result file. Its partial derivatives
are evaluated locally with forward-mode automatic differentiation and entered
in the sparse Jacobian.

The dotted syntax is only a lookup of a registered scalar model variable; it
does not perform Julia property access. Assignment, arbitrary function calls,
and general Julia evaluation are rejected. A missing coordinate is diagnosed
while loading the model. For example, using `pin.theta` requires
`rotation_coordinates = true` on `pin`.

Expressions should presently be smooth. Use a dedicated contact element for a
force law that switches branches and needs event location or a BDF restart.

### Torsional spring-damper

```toml
[pin_spring]
type = "torsional_spring_damper"
markers = ["pendulum.pin", "ground.origin"]
stiffness = 3.0
damping_time_scale = 0.1
free_angle = 0.0
```

| Field | Required | Default |
| --- | --- | --- |
| `markers` | yes, exactly two | — |
| `stiffness` | yes, nonnegative | — |
| `damping_time_scale` | no, nonnegative | `0.0` |
| `damping` | no, nonnegative | `damping_time_scale * stiffness` |
| `free_angle` | no | `0.0` radians; quoted degrees are accepted |

The marker positions locate the viewer symbol and their orientations define the
relative angle. The stored torque is the torque applied to the first marker:

$$
T=-k(\theta-\theta_0)-c\omega.
$$

The second marker receives the opposite torque. Positive stiffness and damping
therefore oppose relative displacement and velocity. An explicit `damping`
overrides the value estimated from the local time scale. The torque is stored
as `element-name.T`. The viewer uses the applied-load color on the first marker
and the reaction color on the second; ground-side symbols are hidden by
default.

### Spanning force

```toml
[spring]
type = "spanning_force"
markers = ["body.tip", "ground.anchor"]
stiffness = 20.0
damping = 0.5
free_length = 0.5
```

The predefined linear spring-damper law requires `stiffness`, `damping`, and
`free_length`. All three values must be nonnegative. A general expression may
be supplied instead:

```toml
[spring]
type = "spanning_force"
markers = ["body.tip", "ground.anchor"]
expression = "-k1*(spring.length-l0)-k3*(spring.length-l0)^3-c*spring.length_rate"
```

A constant signed magnitude uses `force = value`. Exactly one of `force`,
`expression`, or the spring-damper coefficient set is required. The stored
magnitude is the force on the first marker in the direction from the second
marker toward the first. A positive magnitude is compression and a negative
magnitude is tension. Expressions may also reference the body and coordinate
variables listed above and obtain their local partial derivatives with dual
numbers.

`markers` is always required. The two marker points must not initially
coincide, and floating markers are not currently accepted.

The stored spanning vector and unit direction run from the first marker to the
second. Consequently, the global force on the first marker is the negative of
`force * [u_x, u_y]`. The element retains its spanning vector, length, unit
direction, length rate, scalar force, and first-marker global force components
as the nine stored variables `s_x`, `s_y`, `length`, `u_x`, `u_y`,
`length_rate`, `force`, `F_x`, and `F_y`. The viewer reconstructs the connector
from them.

### Planar bushing

```toml
[support]
type = "bushing"
markers = ["body.mount", "ground.mount"]
translational_stiffness = [100.0, 100.0]
rotational_stiffness = 20.0
damping_time_scale = 0.1
free_position = [0.0, -1.0]
free_angle = 0.0
```

| Field | Required | Default |
| --- | --- | --- |
| `markers` | yes, exactly two | — |
| `translational_stiffness` | yes, two nonnegative values | — |
| `damping_time_scale` | no, nonnegative | `0.0` |
| `translational_damping` | no, two nonnegative values | `damping_time_scale * translational_stiffness` |
| `rotational_stiffness` | no, nonnegative | `0.0` |
| `rotational_damping` | no, nonnegative | `damping_time_scale * rotational_stiffness` |
| `free_position` | no | `[0.0, 0.0]` |
| `free_angle` | no | `0.0` radians; quoted degrees are accepted |

Translational values are diagonal coefficients and an unloaded displacement in
the second marker's coordinate system. The reported `F_x`, `F_y`, and `T` are
global loads on the first marker; equal-and-opposite loads act on the second.
The bushing is an ordinary compliant force and contributes no constraint.

`damping_time_scale` is local to this bushing and has units of time. For
$\tau=$ `damping_time_scale`, omitted coefficients are estimated from

$$
C_t=\tau K_t,
\qquad
c_r=\tau k_r.
$$

An explicitly supplied `translational_damping` replaces the translational
estimate, and an explicitly supplied `rotational_damping` independently
replaces the rotational estimate. With no scale or explicit coefficient, the
corresponding damping is zero.

### Plane contact

```toml
[floor_contact]
type = "plane_contact"
markers = ["ball.center", "ground.plane"]
radius = 0.1
stiffness = 10000.0
damping_factor = 0.15
```

| Field | Required | Default |
| --- | --- | --- |
| `markers` | yes, exactly two | — |
| `radius` | yes, positive | — |
| `stiffness` | exactly one of `stiffness` or `expression` | — |
| `expression` | exactly one of `stiffness` or `expression` | — |
| `damping_factor` | no, nonnegative, s/m | `0.0` |
| `transition_depth` | no, nonnegative length | `0.0` |

The first marker is the center of a sphere with the specified radius. The
second marker defines the plane and its local y axis is the outward normal.
Positive gap means separation and negative gap means penetration. This is a
one-sided compliant force, not a constraint. `damping_factor` increases the
force during closing motion and reduces it during rebound without allowing a
tensile contact force.

`transition_depth` makes the onset of the built-in stiffness smooth over the
specified penetration. It is useful when an abrupt stiffness change would
otherwise cause repeated corrector failures. An `expression` may use
`floor_contact.gap` and `floor_contact.gap_rate` and is applied exactly as
written. It is not automatically turned off during separation and is not
clamped to a compressive force. Because that behavior belongs in the
expression, `damping_factor` and `transition_depth` cannot accompany it.

The component owns the output variables `gap`, `gap_rate`,
`normal_force`, `F_x`, and `F_y`. The dynamic integrator locates contact entry,
exit, and damping-branch changes so it can restart promptly when the stiffness
changes. The exact force law and restart policy are given under [Plane contact
in the Technical
Manual](../../architecture/planar/planar-element-formulations.md#plane-contact).

The force is normal to the plane. Applying it at the sphere center is
mechanically equivalent to applying it at the surface contact point because
the center-to-contact offset is parallel to the force and adds no moment.

### Smooth curve

```toml
[cam_profile]
type = "curve"
marker = "cam.profile_frame"
closed = true
points = [
    [0.20, 0.00],
    [0.14, 0.18],
    [0.00, 0.26],
    [-0.14, 0.18],
    [-0.20, 0.00],
    [-0.14, -0.18],
    [0.00, -0.26],
    [0.14, -0.18],
]
```

| Field | Required | Default |
| --- | --- | --- |
| `marker` | yes | — |
| `points` | yes, at least four finite `[x, y]` pairs | — |
| `closed` | no | `true` |

The points are expressed in the named marker frame. The program constructs a
periodic cubic profile that passes through them and is continuous through its
second derivative. Consecutive points must differ and the control polygon must
enclose a nonzero area. Only closed profiles are supported in this release.
The supplied order establishes the outward side; clockwise and
counterclockwise lists are both accepted.

### Curve contact

```toml
[roller_contact]
type = "curve_contact"
curve = "cam_profile"
roller_marker = "follower.center"
radius = 0.05
side = "outside"
stiffness = 50000.0
damping_factor = 0.02
transition_depth = 0.0005
```

| Field | Required | Default |
| --- | --- | --- |
| `curve` | yes, names a `curve` | — |
| `roller_marker` | yes | — |
| `radius` | yes, positive | — |
| `side` | no, `"outside"` or `"inside"` | `"outside"` |
| `stiffness` | exactly one of `stiffness` or `expression` | — |
| `expression` | exactly one of `stiffness` or `expression` | — |
| `damping_factor` | no, nonnegative, s/m | `0.0` |
| `transition_depth` | no, nonnegative length | `0.0` |
| `initial_station` | no | nearest profile point |

The roller center and curve marker must belong to different bodies, although
either may be on ground. The contact station is an explicit local unknown.
The program solves for the point at which the center-to-profile vector is
normal to the curve; it does not hide a nearest-point search inside the force
law. The force acts on the roller along the selected normal and the opposite
force acts at the profile point.

The outputs are `station`, `station_rate`, signed `curvature`, `gap`,
`gap_rate`, `normal_force`, `contact_x`, `contact_y`, `normal_x`, `normal_y`,
`F_x`, and `F_y`. These geometry and rate names may be used in a contact
`expression`. The expression follows the same exact-law convention as
`plane_contact`.

The station is allowed to continue below zero or above the profile length as
the contact travels through the closing point. Its value modulo the profile
length identifies the geometric location. This unwrapped value avoids an
artificial jump in the solver history. A nonconvex curve can have more than
one tangent contact location; `initial_station`, or otherwise the nearest
initial profile point, selects the branch that is then followed continuously.

The complete
[`rotating-cam-follower.toml`](../../models/planar/rotating-cam-follower.toml)
example drives an elliptical cam beneath a guided circular follower.

### Flat follower contact

```toml
[flat_contact]
type = "flat_follower_contact"
curve = "cam_profile"
follower_marker = "follower.face"
stiffness = 50000.0
damping_factor = 10.0
transition_depth = 0.00002
```

| Field | Required | Default |
| --- | --- | --- |
| `curve` | yes, names a `curve` | — |
| `follower_marker` | yes | — |
| `stiffness` | exactly one of `stiffness` or `expression` | — |
| `expression` | exactly one of `stiffness` or `expression` | — |
| `damping_factor` | no, nonnegative, s/m | `0.0` |
| `transition_depth` | no, nonnegative length | `0.0` |
| `initial_station` | no | supporting tangent point |

The follower marker defines the flat face. Its local $x$-axis lies along the
face and its local $y$-axis points in the positive force direction on the
follower. The contact station makes the curve tangent parallel to the face.
The gap is positive when the face separates in its positive $y$ direction and
negative when the cam penetrates it. `side` and `radius` do not apply because
the follower marker defines both the face direction and its location.

The force acts on the follower at the curve contact point, so an offset contact
point produces the correct moment on a follower that is free to rotate. The
outputs, built-in force law, optional exact expression, staging fields, and
unwrapped periodic station are the same as for `curve_contact`. Without an
explicit `initial_station`, initialization chooses the profile point having
the greatest projection along the follower marker's positive $y$-axis and
then corrects it to exact tangency.

The complete
[`rotating-cam-flat-follower.toml`](../../models/planar/rotating-cam-flat-follower.toml)
example drives a guided plate rapidly enough to demonstrate separation and
recontact. The
[`rotating-cam-rocker-follower.toml`](../../models/planar/rotating-cam-rocker-follower.toml)
example puts the same face on a remotely pivoted rocker and demonstrates the
contact moment produced about its revolute joint.

## Friction

The four planar friction elements use the same compliant stick-slip law. They
carry a scalar shear state, use a higher friction coefficient near zero slip,
and limit the resulting force or torque by a normal load inferred from the
referenced constraint reaction. Common fields are:

| Field | Required | Default |
| --- | --- | --- |
| `stiffness` | yes, positive | — |
| `damping` | no, nonnegative | `0.0` |
| `static_coefficient` | yes, at least `dynamic_coefficient` | — |
| `dynamic_coefficient` | yes, nonnegative | — |
| `transition_speed` | no, positive | `0.1` rad/s for revolute; `0.01` m/s otherwise |
| `preload` | joint friction only; nonnegative normal load | `0.0` |
| `release_time` | no, positive seconds | `0.01` |

At zero slip, the stored shear can support a load below the static limit. If
the normal capacity disappears, the force becomes zero and the shear decays
with `release_time`. Static analysis measures shear from the corrected initial
position, so friction can balance a load without artificial sliding. The
static shear is retained when dynamics begins.

### Surface friction

```toml
[friction]
type = "surface_friction"
contact = "support"
stiffness = 2000.0
damping = 90.0
static_coefficient = 0.7
dynamic_coefficient = 0.5
transition_speed = 0.025
```

`contact` must name a planar `plane_contact`. Friction acts along the plane
marker's local $x$-axis at the projected contact point. The sphere's surface
velocity includes its rotation, and motion of a body-owned plane is also
included. The contact's positive normal force sets the friction capacity;
while separated, the tangential force is zero and stored shear decays with
`release_time`. `preload` is not accepted because it would create friction
without contact.

The element reports `shear`, signed tangential `slip`, scalar `friction`,
`F_x`, and `F_y`. Equal-and-opposite forces act at the same projected contact
point. See
[`sliding-block-surface-friction.toml`](../../models/planar/sliding-block-surface-friction.toml).

### Revolute bearing friction

```toml
[bearing_friction]
type = "revolute_friction"
joint = "pin"
effective_radius = 0.08
stiffness = 15.0
damping = 1.6
static_coefficient = 1.0
dynamic_coefficient = 0.8
transition_speed = 0.2
```

`joint` must name a planar `revolute`. The positive `effective_radius` converts
the magnitude of its two-component point-force reaction into a friction-torque
limit. The element reports `shear`, relative angular `slip`, `bearing_load`,
and `torque`. Equal-and-opposite torques act on the two joined bodies. See
[`revolute-bearing-friction.toml`](../../models/planar/revolute-bearing-friction.toml).

### Translational guide friction

```toml
[guide_friction]
type = "translational_friction"
joint = "guide"
stiffness = 2000.0
damping = 100.0
static_coefficient = 0.8
dynamic_coefficient = 0.6
transition_speed = 0.025
```

`joint` must name a planar `translational` joint. Friction acts along the
second marker's local $x$-axis. The magnitude of the joint's transverse
inplane reaction, plus optional preload, is the guide load. The perpendicular
constraint's reaction torque is not converted into contact pressure. Outputs
are `shear`, axial `slip`, `guide_load`, scalar `force`, `F_x`, and `F_y`.
See [`translational-guide-friction.toml`](../../models/planar/translational-guide-friction.toml).

### Inplane constraint friction

```toml
[plane_friction]
type = "inplane_friction"
constraint = "support"
stiffness = 2000.0
damping = 100.0
static_coefficient = 0.8
dynamic_coefficient = 0.6
transition_speed = 0.025
```

`constraint` must name a standalone planar `inplane`. Friction acts along the
second marker's local $x$-axis and uses the magnitude of the signed normal
reaction plus optional preload. It reports `shear`, tangential `slip`,
`normal_load`, scalar `force`, `F_x`, and `F_y`. An ideal inplane is bilateral,
so this element permits friction for either reaction sign. It does not detect
separation. Use `plane_contact` with `surface_friction` when the bodies must be
able to separate. See
[`inplane-friction-block.toml`](../../models/planar/inplane-friction-block.toml).

## Motion generator

The differentiated constraint equations and reaction conventions are given
under [Motion generators in the Technical
Manual](../../architecture/planar/planar-element-formulations.md#motion-generators).

### Translational distance

For constant speed:

```toml
[drive]
type = "translational_motion"
markers = ["slider.center", "ground.drive_axis"]
function = "constant_speed"
initial_distance = 0.0
velocity = "drive_speed"
```

The ordered markers define the signed distance

$$
d=(P_i-P_j)^T\hat y_j,
$$

where $\hat y_j$ is the second marker's oriented local $y$-axis. Thus the
second marker defines a reference plane normal to that axis, and the first
marker is prescribed to remain at distance $d(t)$ from the plane. This is one
primitive translational constraint; it does not constrain relative motion
along the plane or relative rotation.

`velocity` is required and may be a number or a numeric parameter/expression
evaluated at the initial time. `initial_distance` defaults to zero. The loader
supplies consistent distance, velocity, and acceleration laws.

For a general prescribed history:

```toml
[drive]
type = "translational_motion"
markers = ["slider.center", "ground.drive_axis"]
function = "expression"
distance = "0.1*sin(2*t)"
```

The program obtains velocity and acceleration by automatically differentiating
the distance expression, which must be twice differentiable over the
simulation interval. Do not supply separate `velocity` or `acceleration`
fields. The component retains `distance_m`, `velocity_m`,
`acceleration_m`, and the scalar reaction `force_m` as explicit variables.
Positive reaction acts on the first marker along $\hat y_j$, with the
equal-and-opposite force on the second. Floating markers are not currently
accepted. The rack-and-pinion constraint uses the corresponding signed
translation along a marker axis.

### Rotational motion

For constant speed:

```toml
[crank_drive]
type = "rotational_motion"
markers = ["crank.ground_end", "ground.origin"]
function = "constant_speed"
initial_angle = 0.0
angular_velocity = "speed"
```

`markers` is required and ordered. `angular_velocity` is required and may be a
number or a numeric parameter/expression evaluated at the initial time.
`initial_angle` defaults to zero. It may be a number in radians or a quoted
degree value. The loader supplies consistent angle, angular-velocity, and
angular-acceleration laws.

For a general prescribed history:

```toml
[crank_drive]
type = "rotational_motion"
markers = ["crank.ground_end", "ground.origin"]
function = "expression"
angle = "0.25*sin(4*t)"
```

The program obtains angular velocity and acceleration by automatically
differentiating the angle expression. It must be twice differentiable and use
the restricted expression language listed above. Do not supply separate
`angular_velocity` or `angular_acceleration` fields. Omitting `function`
selects `expression`.

## Command line

```text
simp2d [MODEL.toml|-] [duration] [samples]
                    [--output RESULT.simp] [--overwrite]
```

A missing model name or `-` reads TOML from standard input. `duration` is
measured from the chosen start time and overrides `simulation.end_time`.
`samples` overrides `simulation.output_samples`. `--output` writes a result;
without it, the simulation runs and prints a summary but is not retained.

```bash
./bin/simp2d model.toml
./bin/simp2d model.toml 10.0 1001 \
    --output results/examples/planar/run.simp
```

The [result-file description](../../architecture/common/simulation-result-files.md)
documents viewing, CSV conversion, and extraction of the embedded TOML.

## Current boundaries

The executable format does not yet include spatial bodies, units, reusable
element definitions, conditional element behavior, or compound elements with
an implemented interface contract.
Qualified nested tables leave room for compound elements without committing to
that future design now.
