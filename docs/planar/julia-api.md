# Sim2D Julia API

The Julia API builds the same hierarchical planar model document accepted from
TOML. It then uses the ordinary Sim2D loader, validation, initial-condition
correction, state selection, sparse assembly, and analysis code. Julia input is
not a separate solver path.

```julia
using PracticalMechanicalSimulation

const Sim2D = PracticalMechanicalSimulation.Sim2D
```

## A complete model

```julia
model = Sim2D.Model(:pendulum;
    title = "Pendulum built in Julia")

Sim2D.analysis!(model; mode = :automatic)
Sim2D.simulation!(model;
    end_time = 3.0,
    frames_per_second = 60,
    maximum_step = 0.01)

ground = Sim2D.ground!(model, :ground)
ground_pin = Sim2D.marker!(ground, :pin)

pendulum = Sim2D.rigid_body!(model, :pendulum;
    mass = 1.0,
    inertia = 1 / 12,
    position = [0.5, 0.0])
body_pin = Sim2D.marker!(pendulum, :pin;
    position = [-0.5, 0.0])
tip = Sim2D.marker!(pendulum, :tip;
    position = [0.5, 0.0])

Sim2D.graphics!(pendulum; show_default = false, color = :steelblue)
Sim2D.graphic!(pendulum, :member;
    shape = :cylinder,
    markers = [body_pin, tip],
    radius = 0.04)

pin = Sim2D.revolute!(model, :pin;
    markers = [body_pin, ground_pin],
    rotation_coordinates = true)
Sim2D.gravity!(model, :gravity;
    acceleration = [0.0, -9.81],
    bodies = [pendulum])

Sim2D.state_selection!(model;
    method = :preferred,
    preferred_velocities = [Sim2D.variable(pin, :omega)],
    allow_fallback = false)

result = Sim2D.run(model)
Sim2D.save_result("pendulum.simp", result; overwrite = true)
```

`save_result` writes the ordinary portable `.simp` history and adds the native
graphics section read directly by SimpView. Use the general `write_result`
function instead when the graphics section is not wanted.

## Handles and reusable assemblies

Builder functions return `ElementRef` handles. A marker builder takes its body
or ground handle, so ownership and its qualified name do not need to be
repeated. Fields that expect element names accept handles directly. A handle
from another model is rejected.

Ordinary Julia functions can build reusable assemblies. They receive a model
and name prefix, add their elements, and return the handles that callers may
need. Dotted names produce the corresponding hierarchy in the generated model
and in SimpView.

The executable
[`julia_api_double_pendulum.jl`](../../examples/planar/julia_api_double_pendulum.jl)
uses a `double_pendulum!` assembly that calls `pendulum_link!` twice:

```text
pendulum
├── first_link
│   ├── inner
│   └── outer
├── second_link
│   ├── inner
│   └── outer
├── base_joint
├── elbow_joint
└── gravity
```

Run it with:

```bash
julia --project=. examples/planar/julia_api_double_pendulum.jl
bin/simpView
```

## Model settings and elements

The following functions update their corresponding top-level model tables:

- `analysis!`
- `simulation!`
- `state_selection!`
- `initial_conditions!`
- `parameters!`
- `graphics!`

Their keywords are the fields documented in the
[Planar TOML User's Guide](toml-reference.md). Symbols are converted to
strings. `simulation!` additionally accepts `frames_per_second` and derives
the output sample count from the configured start and end times.

The API provides named builders for the current planar library: bodies,
flexible beams, markers, joint primitives and compound joints, measurements, gears,
rack-and-pinion sets, belts, forces, bushings, contact, motion generators, and
user equation components. The friction builders are `surface_friction!`,
`revolute_friction!`, `translational_friction!`, and `inplane_friction!`. For
example:

```julia
beam = Sim2D.flexible_beam!(model, :beam;
    length = 1.0, area = 0.01,
    second_moment = 8.333333333333333e-6,
    elastic_modulus = 2.0e7, shear_modulus = 8.0e6,
    mass = 1.0, position = [0.5, 0.0])
beam_root = Sim2D.beam_marker(beam, :end_i)
beam_tip = Sim2D.beam_marker(beam, :end_j)
```

The section conveniences use the same fields as TOML. For example, a solid
circular beam can derive its properties from material and section data:

```julia
beam = Sim2D.flexible_beam!(model, :beam;
    length = 1.0,
    density = 7800.0,
    elastic_modulus = 2.0e11,
    poisson_ratio = 0.30,
    section = (shape = :circular, diameter = 0.025),
    position = [0.5, 0.0])
```

Use `section = (shape = :rectangular, width = ..., height = ...)` for a
rectangular member. The expanded canonical properties and generated member
graphic are retained in the model stored with the result.

`beam_marker` returns handles for the beam's generated `end_i`, `cm`, and
`end_j` markers without adding duplicate marker tables to the model.
The planar beam is the same simplified, small-deformation element documented
in the [TOML flexible-beam reference](toml-reference.md#flexible-beam); the
Julia interface does not add deformation-dependent inertia or nonlinear beam
effects.

For example:

```julia
profile = Sim2D.curve!(model, :cam_profile;
    marker = profile_frame,
    points = [[0.20, 0.0], [0.0, 0.26],
              [-0.20, 0.0], [0.0, -0.26]])
Sim2D.curve_contact!(model, :roller_contact;
    curve = profile, roller_marker = roller_center,
    radius = 0.05, stiffness = 50_000.0)
Sim2D.flat_follower_contact!(model, :flat_contact;
    curve = profile, follower_marker = follower_face,
    stiffness = 50_000.0)
```

The curve and contact fields have the same meanings as in the TOML guide.
For example, a Julia vector comprehension may generate any desired set of
profile points before calling `curve!`.

User equation components use the same builder path:

```julia
Sim2D.equation_component!(model, :controller;
    inputs = Dict(:error => "pin.theta", :omega => "pin.omega"),
    parameters = Dict(:kp => 12.0, :ki => 10.0, :kd => 5.0),
    states = Dict(:integral_error => Dict(:initial => 0.0, :static => :hold)),
    variables = Dict(:torque => Dict(:initial => 0.0)),
    state_equations = Dict(:integral_error => "-error"),
    equations = ["torque = -kp*error - ki*integral_error - kd*omega"])
```

Builder names follow element types with a trailing `!`. `element!` is the
generic escape hatch for a newly introduced type before it receives a named
convenience function. All generated fields still pass through the ordinary
planar loader's validation.

`graphics!(element; ...)` attaches the element's primary graphics table.
`graphic!(element, name; ...)` adds one named graphic below it.

## Load, run, and inspect

The package entry points accept a Julia-built model directly:

```julia
loaded = load_planar_model(model)
result = run_planar_model(model; duration = 2.0, samples = 121)
```

`Sim2D.load_model(model)`, `Sim2D.simulate(model)`, and `Sim2D.run(model)` are
shorter equivalents. The returned objects are the same types used for TOML
models.

`Sim2D.document(model)` returns a detached dictionary. The ordinary planar
loader also accepts an `AbstractDict` directly and normalizes symbolic keys and
values to their portable strings.

```julia
specification = Sim2D.document(model)
loaded = load_planar_model(specification)
Sim2D.write_model("generated-model.toml", model; overwrite = true)
```

The generated document is embedded in every stored result. Model extraction,
saved-result initialization, CSV conversion, and SimpView therefore behave the
same way for Julia-built, TOML, and Lua models.

## Current boundary

The API builds model specifications; allocated variable and equation indices
remain internal. Force fields accept the same portable expression strings as
TOML. Arbitrary Julia closures are not accepted as force laws because they
could not be reconstructed from a saved result.
