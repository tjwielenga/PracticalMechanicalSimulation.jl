# Spatial Modeling Assemblies

A modeling assembly packages a useful mechanical construction without adding
a new solver element. It creates ordinary bodies, markers, joints, forces, and
graphics already understood by Sim3D.

Assemblies may be written as ordinary Lua modules or as declarative TOML
templates. Lua is preferable when construction requires optional values,
conditions, loops, or substantial vector calculations. TOML remains useful
for assemblies that are mostly fixed collections of elements. Both forms
produce the same ordinary primitive model tables.

## Lua model file

A Lua-authored model loads an assembly with Lua's normal `require` function.
The returned function is called with one table, and `name` identifies the
assembly instance. Lua loads and caches the module once; the returned function
may then be called for any number of instances:

```lua
local sim3d = require "sim3d"
local model, ground, rigid_body =
    sim3d.model, sim3d.ground, sim3d.rigid_body
local leaf_spring = require "spatial.leaf_spring"

model {
    name = "leaf_spring_test",
    dimension = "spatial"
}

ground {
    name = "ground"
}

rigid_body {
    name = "axle",
    mass = 40.0,
    inertia = {1.0, 2.0, 2.0}
}

leaf_spring {
    name = "rear_left",
    frame = "ground",
    leaf = "axle",
    front_eye = {-0.64, 0.75, -0.07},
    leaf_center = {0.0, 0.75, 0.0},
    rear_shackle_eye = {0.70, 0.75, 0.07},
    rear_frame_eye = {0.71, 0.75, 0.20},
    leaf_stiffness = 36532.0,
    nominal_load = 3800.0,
    width = 0.076,
    mount_stiffness = {7.0e6, 7.0e6, 1.0e6},
    mount_rotational_stiffness = {1.0e4, 1.0e4, 81.4}
}
```

`nominal_load` is the vertical load carried by one complete leaf spring at
the entered geometry. The assembly converts half of that load at each inner
leaf joint into unloaded bending angles, so the front and rear bending
bushings are already preloaded after initial-condition assembly. It defaults
to zero. `hotchkiss_suspension` accepts the same per-leaf value and passes it
to both leaf springs.

`model`, `analysis`, `simulation`, `parameters`, `graphics`,
`state_selection`, and `initial_conditions` describe the corresponding
single model tables. Built-in elements are also functions and require `name`.
For example:

```lua
marker {
    name = "ground.origin",
    position = {0.0, 0.0, 0.0}
}
```

Lua does not have a separate `from ... import ...` statement. Assigning module
members to local names is its normal equivalent. Only the names used by a
model need to be listed. They can instead remain qualified when that is
clearer:

```lua
local sim3d = require "sim3d"

sim3d.marker {
    name = "ground.origin"
}
```

Both forms give the Lua extension explicit definitions for completion and
diagnostics.

The assembly function runs when it is called and adds ordinary primitive
tables to the same model entry list. After the Lua file finishes, the reader
resolves names and reference-frame conversions throughout the completed set
of tables. Declaration order therefore does not control reference resolution.

Run a Lua-authored spatial model exactly like a TOML model:

```bash
./bin/simp3d models/spatial/leaf-spring-assembly.lua \
    --output results/examples/spatial/leaf-spring-assembly.simp --overwrite
```

## Converted vehicle assemblies

The Lua assembly directory also contains a staged conversion of the author's
earlier large-van vehicle model:

- `sla_suspension.lua` constructs one short-long-arm front suspension;
- `stabilizer_bar.lua` constructs the segmented front anti-roll bar;
- `steering_linkage.lua` constructs the pitman-arm steering linkage;
- `leaf_spring.lua` and `hotchkiss_suspension.lua` construct the solid rear
  axle and its two leaf springs;
- `gylt_245_75_r16.lua` constructs the wheel bodies, axle joints, graphics,
  and Sim3D rolling tires from the available LT245/75R16 data;
- `vehicle_force_elements.lua` supplies linear springs, piecewise-linear
  tabulated dampers, and one-sided sphere-plane bump stops; and
- `large_van.lua` supplies the original hardpoints, adds a simple painted
  and glazed body surface, and assembles those pieces into one vehicle.

The pitman-arm connection to the steering cross link is a spherical joint plus
a Perp primitive. The pitman marker's first axis lies along the arm. The cross-
link marker's second axis is the fore-aft direction projected perpendicular to
that arm. The spherical joint transmits the steering motion while the Perp
prevents the cross link from freely rolling about its connections. The optional
assembly parameter `fore_aft` supplies the vehicle fore-aft direction and
defaults to the global $x$ direction; its sign does not affect the constraint.

The examples `models/spatial/large-van.lua` and
`models/spatial/large-van-high-cg.lua` instantiate the same assembly and
30 m/s right-steer maneuver. The high-CG version raises only the body center
of mass by 40 mm. Both include four roof-corner sphere-plane contacts against
the road, inactive during static assembly. Run `large-van-static.lua` first;
both dynamic models read its settled position and then establish their own
30 m/s consistent velocities. All three runs use the same steady-slip tire
model without relaxation states. Direct, mass-regularized static Newton is
used to settle the vehicle before either maneuver. The `normalSteer.py`
dependency stops after constructing steering
bodies, but the older `normalSteer2.pl` supplies the missing steering-gear
construction. The Lua assembly follows that formulation: a coordinate coupler
imposes the gear ratio, a separate rotating shaft permits windup, and a
torsional spring-damper connects that shaft to the pitman arm. When free play
is nonzero, the spring torque instead uses the Perl model's gear-side dead
band. The optional pitman rotational damper is also retained.

```bash
bin/simp3d models/spatial/large-van-static.lua \
    --output results/examples/spatial/large-van-static.simp --overwrite
bin/simp3d models/spatial/large-van.lua \
    --output results/examples/spatial/large-van.simp --overwrite
bin/simp3d models/spatial/large-van-high-cg.lua \
    --output results/examples/spatial/large-van-high-cg.simp --overwrite
bin/simpView
```

`large-van-modal.lua` is a separate linearization example. It explicitly uses
transient tire states, but neither dynamic rollover run depends on it.

The steering column is part of the steering-wheel body. Its outer endpoint
owns an oriented marker carrying a ring-shaped steering-wheel graphic. The
column revolute, gear coupler, and any rotational generator therefore move the
column and steering wheel together.

The large-van front jounce stops use a 30 mm sphere on the lower suspension
marker and a plane on the upper marker. The front rebound stops reverse those
roles. At the rear, the jounce sphere is fixed to the van body and its plane is
fixed to the axle. Each plane normal follows the former marker-to-marker span
line and points toward the sphere. Each stop uses a plane-contact expression
that supplies the historical 4 MN/m stiffness and 20 kN s/m viscous damping.
The expression owns its one-sided clamp and uses a one-micrometre activation
ramp so the full damping becomes available almost immediately after contact.
The assembly parameters `bump_stop_stiffness`, `bump_stop_damping`, and
`bump_stop_activation_depth` override these values. Stops remain inactive
during static equilibrium.

`wheel_motion_angle` may be a number or time expression. It places a
rotational generator on the steering-column revolute while retaining the gear
coupler, windup compliance, and pitman damper. The large-van wrapper exposes it
as `steering_wheel_motion_angle`. The diagnostic `pitman_motion_angle` option
instead drives the pitman directly and omits that upstream steering gear; the
two motion options cannot be used together.

The large-van assembly estimates each rear leaf's nominal load from the
historical measured rear axle load after subtracting the rear unsprung mass.
It can be overridden with `rear_leaf_nominal_load`. This restores the preload
calculation that was present but commented out in the earlier script; the old
positional Python call did not intentionally pass that value into the leaf
assembly.

The tire wrapper and its 60 and 88 psi property files supply wheel-and-tire
mass, section and rim dimensions, rolling radius, normal stiffness,
longitudinal slip stiffness, cornering stiffness, camber stiffness, and the
friction limit. By default, the wrapper uses current slip without tire
relaxation. Supplying both `longitudinal_relaxation_length` and
`lateral_relaxation_length` enables two first-order tread-deformation states
that retain tangential loading at zero transport speed. The large-van
assembly passes these through as `tire_longitudinal_relaxation_length` and
`tire_lateral_relaxation_length`. The current element does not yet apply
rolling resistance, aligning torque, or overturning moment.

The tire assembly option `static_spin_stiffness` adds a torsional bushing in
parallel with the axle revolute. It acts only during static analysis and
removes the otherwise neutral wheel-spin coordinate. Its optional
`static_spin_damping_time_scale` defaults to 0.10 seconds. The bushing is
absent when the stiffness is zero and is automatically unloaded before a
dynamic or modal stage. `show_default = false` hides the tire cylinder while
leaving the body's inertia ellipsoid available in the viewer; `visible =
false` hides both.

The historical Pogo tire was an additional contact model used when a rolling
tire overturned or struck the ground outside its normal operating attitude.
It was not a replacement for the ADAMS tire. Sim3D does not yet have that
separate three-dimensional rollover contact, so the converted assembly records
it as an explicit remaining limitation.

The shock tables are preserved as linearly extrapolated, piecewise-linear
force expressions. Bump stops use spatial plane contacts whose expressions
define their stiffness, damping, contact gating, and no-tension behavior.
These use ordinary Sim3D elements and remain visible in the expanded model and stored results. The
jounce and rebound bumpers default to `inactive_during = "static"`, so a
temporary contact cannot interfere with settling the suspension. They become
active for the subsequent dynamic simulation.

## Lua assembly module

A Lua assembly module returns an ordinary function. That function calculates
its geometry and calls the same primitive model functions used by the main
model:

```lua
local sim3d = require "sim3d"
local rigid_body, marker = sim3d.rigid_body, sim3d.marker
local vector, norm, link_frame, box_inertia =
    sim3d.vector, sim3d.norm, sim3d.link_frame, sim3d.box_inertia

local function example_linkage(p)
    local first = vector(assert(p.first_point, "first_point is required"))
    local second = vector(assert(p.second_point, "second_point is required"))
    local side = vector(p.side or {0.0, 0.0, 1.0})
    local thickness = p.thickness or p.width/5
    local length = norm(second-first)
    local center = 0.5*(first+second)
    local size = {length, thickness, p.width}
    local mass = (p.density or 7850.0)*length*thickness*p.width
    local body = p.name .. ".link"

    rigid_body {
        name = body,
        mass = mass,
        inertia = box_inertia(mass, size),
        center_of_mass = body .. ".cm",
        position = center,
        graphics = {
            shape = "box",
            marker = body .. ".cm",
            size = size,
            color = p.color or "steelblue"
        }
    }

    marker {
        name = body .. ".cm",
        orientation = link_frame(first, second, side)
    }
end

return example_linkage
```

The model loads and calls it in the usual Lua form:

```lua
local sim3d = require "sim3d"
local model = sim3d.model
local example_linkage = require "spatial.example_linkage"

example_linkage {
    name = "steering_link",
    first_point = {0.0, 0.0, 0.0},
    second_point = {0.5, 0.2, 0.0},
    width = 0.03
}
```

Input vectors are ordinary Lua tables. Call `vector` before using vector
arithmetic on them. Lua's `nil`, Boolean operations, `if` statements, and
loops provide defaults and optional construction directly.

Lua assembly code runs through an embedded Lua 5.4 runtime. File, process,
and Julia access remain unavailable. `require` is restricted to Lua modules
in the model directory and the Sim3D assembly directory; native-library
loading is disabled. Assembly files are nevertheless model code and should
come from a trusted source.

## Declarative TOML definition

An assembly definition starts with its name and dimension:

```toml
[assembly]
name = "example_linkage"
dimension = "spatial"
```

Interfaces identify bodies supplied by the enclosing model:

```toml
[interface.frame]
kind = "body_or_ground"

[interface.moving_body]
kind = "body"
```

The implemented interface kinds are `body` and `body_or_ground`.

Parameters are typed and may be required or have defaults:

```toml
[parameters]
first_point = { kind = "vector", required = true }
second_point = { kind = "vector", required = true }
width = { kind = "number", required = true }
thickness = { kind = "number", default = "=width/5" }
color = { kind = "string", default = "steelblue" }
```

The implemented parameter kinds are `number`, three-component `vector`,
`string`, and `boolean`. A string beginning with `=` is a build-time
expression. An ordinary string remains an ordinary TOML value; this
distinction leaves runtime force expressions unchanged.

Intermediate calculations keep repeated geometry out of the element tables:

```toml
[calculate]
axis = "=unit(second_point-first_point)"
length = "=norm(second_point-first_point)"
center = "=0.5*(first_point+second_point)"
size = "=[length,thickness,width]"
mass = "=density*length*thickness*width"
inertia = "=box_inertia(mass,size)"
```

Calculations are resolved by their dependencies rather than their order in the
file. A dependency cycle or unknown name is an input error.

The remainder of the file looks like an ordinary spatial model. Its fields may
use calculated values:

```toml
[link]
type = "rigid_body"
mass = "=mass"
inertia = "=inertia"
center_of_mass = "link.cm"
position = "=center"

[link.cm]
type = "marker"
orientation = "=link_frame(first_point,second_point,side)"
```

## Geometry functions

Lua builders and declarative TOML assembly expressions share these geometry
functions:

- `abs`, `sqrt`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `min`, and
  `max`;
- `dot`, `cross`, `norm`, and `unit`;
- `frame(x, y, z)` to construct an orientation from three column axes;
- `link_frame(point_a, point_b, side)` to construct link principal axes;
- `rotation(angle, axis)` and `compose(A, B)` for orientations;
- `box_inertia(mass, size)` for the three principal moments of a uniform box;
- `local_point(body_name, global_point)` and
  `local_frame(body_name, global_frame)` to express global assembly geometry
  in a named body's reference frame. These conversions are completed in the
  second pass after all body tables are available.

Angles passed to `rotation` may be radians or quoted degree values such as
`"5 deg"` and `"5°"`.

`model_table` adds an untyped table when a Lua assembly needs to contribute
graphics or other data below an existing model object:

```lua
model_table {
    name = p.body .. ".graphics." .. p.name,
    shape = "box",
    marker = p.body .. "." .. p.name .. ".center",
    size = size
}
```

The same method can add an indexed `surface` graphic. The assembly supplies
one `vertices` array and named material patches whose `faces` contain one-based
indices into that array. This keeps geometry calculations readable in Lua and
avoids repeating shared vertices. See `assemblies/spatial/simple_body.lua` for
a bodywork example with separate painted and glass patches.

`assemblies/spatial/ground_pad.lua` provides a simpler example. It attaches a
rectangular surface to an oriented ground marker. The dynamic large-van models
use it for a 360 m by 360 m dark-gray asphalt pad. The pad sets
`include_in_fit = false`, so the viewer initially frames the vehicle rather
than the entire road.

## Declarative TOML names and interfaces

An internal name is prefixed with the instance name. For an instance named
`rear_left`, `link.cm` becomes `rear_left.link.cm`.

A table whose first name is an interface contributes an item to the supplied
body. For example:

```toml
[moving_body.connection]
type = "marker"
position = "=local_point(moving_body,connection_point)"
```

If `moving_body = "axle"` on instance `rear_left`, the expanded marker is
`axle.rear_left.connection`. Graphics contributed through an interface are
placed in the supplied body's `graphics` hierarchy in the same way.

Use `@interface_name` when a nested assembly must receive one of the parent
assembly's interface bodies directly:

```toml
[subassembly]
type = "another_assembly"
frame = "@frame"
```

An assembly definition imports definitions needed by its nested instances in
its own header:

```toml
[assembly]
name = "vehicle_corner"
dimension = "spatial"
assemblies = ["leaf_spring.toml", "damper.toml"]
```

This allows a vehicle assembly to contain suspension assemblies, which may in
turn contain still smaller assemblies. Expansion stops at the ordinary
primitive model elements.

## Building up a vehicle model

The vehicle examples begin with a rigid chassis and four production tire
assemblies before adding suspension subsystems. This gives each new layer a
small, understandable test model.

First solve and save the static configuration:

```sh
bin/simp3d models/spatial/vehicle-development/rigid-four-tire-static.lua \
  --output results/examples/spatial/rigid-four-tire-static.simp --overwrite
```

The vehicle has neutral road-plane translations and heading. The static solver
detects the resulting singular Jacobian and retains mass and inertia in its
correction matrix. This leaves those neutral coordinates at their supplied
values while heave, roll, and pitch settle from force balance. The tire
deformation states supply tangential forces but are not constraints. The
companion dynamic model reads the saved equilibrium and assigns the chassis
velocity and free-rolling wheel speeds:

```sh
bin/simp3d models/spatial/vehicle-development/rigid-four-tire-dynamic.lua \
  --output results/examples/spatial/rigid-four-tire-dynamic.simp --overwrite
```

For the nominal symmetric model, the static result carries 4.405 kN on each
tire. The next development stages can therefore add the front suspension,
rear suspension, and steering separately without hiding a basic tire or
static-transfer problem inside the complete vehicle.

The next stage replaces the two front rigid wheel attachments with independent
short-long-arm (SLA) suspensions. `assemblies/spatial/sla_suspension.lua`
builds each corner from ordinary Sim3D elements: hinged upper and lower
control arms, a spindle, and two spherical ball joints. The vehicle model adds
a spring-damper from the chassis to the lower arm and a constant-length toe
link from the chassis to the spindle. Keeping those elements outside the SLA
assembly allows a later vehicle model to substitute a shock absorber and a
steering linkage without changing the suspension geometry.

The earlier `assemblies/spatial/macpherson_suspension.lua` remains available
as a separate assembly, but it is no longer used by this vehicle-development
model.

During static assembly, a static-only bushing weakly holds the chassis in the
otherwise neutral road-plane $x$, $y$, and yaw directions. It has no heave,
roll, or pitch stiffness, so ride height and attitude remain determined by the
suspension and tires. `active_during = "static"` keeps the hold active through
dynamic relaxation and Newton polish but removes its loads before dynamic
initialization. The dynamic model reads the named body and joint coordinates
from the static result and starts with a completely free chassis.
Run the pair in order:

```sh
bin/simp3d models/spatial/vehicle-development/front-sla-static.lua \
  --output results/examples/spatial/front-sla-static.simp --overwrite

bin/simp3d models/spatial/vehicle-development/front-sla-dynamic.lua \
  --output results/examples/spatial/front-sla-dynamic.simp --overwrite
```

The model retains the rigid rear wheel attachments so that the new front
suspension can be judged separately. In the nominal symmetric static result,
the SLA springs have a stiffness of 175 kN/m. Their initial free lengths retain
the original 4.2 kN preload, but the increased stiffness accounts for the
motion ratio between the lower-arm spring attachment and the ball joint. The
chassis settles about 36 mm from its entered position. The left and right front
tire loads are both about 4.576 kN and the rear loads are about 4.234 kN.
Static equilibrium takes nine relaxation cycles followed by one Newton
correction. The one-second dynamic run uses the transient tire model and 61
output samples.
