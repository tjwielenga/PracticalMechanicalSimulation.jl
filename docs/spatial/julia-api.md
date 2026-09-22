# Sim3D Julia API

The Julia API builds the same hierarchical model document accepted from TOML
and produced by Lua assemblies. It then uses the ordinary Sim3D loader,
validation, initial-condition correction, state selection, sparse assembly,
and analysis code. Julia input is not a separate solver path.

```julia
using PracticalMechanicalSimulation

const Sim3D = PracticalMechanicalSimulation.Sim3D
```

## A complete model

```julia
model = Sim3D.Model(:pendulum;
    title = "Pendulum built in Julia")

Sim3D.analysis!(model; mode = :automatic)
Sim3D.simulation!(model;
    end_time = 3.0,
    frames_per_second = 60,
    maximum_step = 0.01)
Sim3D.graphics!(model;
    background = :white,
    body_palette = :colorblind)

ground = Sim3D.ground!(model, :ground)
ground_pin = Sim3D.marker!(ground, :pin)

pendulum = Sim3D.rigid_body!(model, :pendulum;
    mass = 1.0,
    inertia = [0.001, 1 / 12, 1 / 12],
    position = [0.5, 0.0, 0.0],
    velocity = [0.0, 0.2, 0.0],
    angular_velocity = [0.0, 0.0, 0.4])
Sim3D.graphics!(pendulum;
    shape = :box,
    size = [1.0, 0.08, 0.08],
    color = :steelblue)
body_pin = Sim3D.marker!(pendulum, :pin;
    position = [-0.5, 0.0, 0.0])

pin = Sim3D.revolute!(model, :pin;
    markers = [body_pin, ground_pin],
    rotation_coordinates = true)
Sim3D.gravity!(model, :gravity;
    acceleration = [0.0, -9.81, 0.0],
    bodies = [pendulum])

Sim3D.state_selection!(model;
    method = :preferred,
    preferred_velocities = [Sim3D.variable(pin, :omega)],
    allow_fallback = false)

result = Sim3D.run(model)
Sim3D.save_result("pendulum.simp", result; overwrite = true)
```

`save_result` writes the ordinary portable `.simp` history and adds its native
graphics section. The result can be opened directly by SimpView. Use the
general `write_result` function instead when the graphics section is not
wanted.

The executable version of this example is
[`julia_api_pendulum.jl`](../../examples/spatial/julia_api_pendulum.jl).

## Handles and references

Builder functions return `ElementRef` objects. A marker builder takes its body
or ground handle, so ownership and its qualified name do not need to be
repeated:

```julia
body = Sim3D.rigid_body!(model, :body; mass = 5.0,
    inertia = [1.0, 2.0, 2.0])
tip = Sim3D.marker!(body, :tip; position = [1.0, 0.0, 0.0])
```

Here `tip` refers to `body.tip`. Element fields that expect names accept these
handles directly:

```julia
joint = Sim3D.spherical!(model, :joint;
    markers = [tip, ground_point])
```

`Sim3D.variable(element, name)` returns a qualified `VariableRef` for state
preferences, expressions constructed by higher-level tools, and result
selection. A reference from a different model is rejected rather than silently
turning into an invalid name.

Julia functions can create reusable assemblies by receiving a model and a
name prefix, adding elements, and returning their handles:

```julia
function link!(model, name; length, mass)
    body = Sim3D.rigid_body!(model, name;
        mass,
        inertia = [mass * length^2 / 100,
                   mass * length^2 / 12,
                   mass * length^2 / 12])
    first = Sim3D.marker!(body, :first;
        position = [-length / 2, 0.0, 0.0])
    second = Sim3D.marker!(body, :second;
        position = [length / 2, 0.0, 0.0])
    (; body, first, second)
end
```

The function is ordinary Julia and may use loops, arrays, calculations, and
other assembly functions. A dotted body name such as `"vehicle.front.link"`
creates corresponding hierarchy in the portable document.

The executable
[`julia_api_double_pendulum.jl`](../../examples/spatial/julia_api_double_pendulum.jl)
shows the full pattern. Its `double_pendulum!` assembly calls
`pendulum_link!` twice and places all generated elements below `pendulum`:

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

The assembly returns handles for its bodies, joints, and subassemblies, so the
calling model can refer to their variables without reconstructing qualified
names. Run and view it with:

```bash
julia --project=. examples/spatial/julia_api_double_pendulum.jl
bin/simpview-web
```

A floating-reference spatial beam uses the same fields as TOML:

```julia
beam = Sim3D.flexible_beam!(model, :beam;
    mass = 1.0, length = 1.0, area = 0.01,
    elastic_modulus = 2.0e7, shear_modulus = 8.0e6,
    second_moment_y = 8.333333333333333e-6,
    second_moment_z = 8.333333333333333e-6,
    torsion_constant = 1.6666666666666667e-5,
    damping_time_scale = 0.002,
    position = [0.5, 0.0, 0.0])
```

The loader creates the qualified markers `beam.end_i`, `beam.cm`, and
`beam.end_j`. They may be referenced by name in builders until dedicated
generated-marker handles are added to the Julia interface.

## Model and analysis settings

The following functions update the corresponding top-level model table:

- `analysis!`
- `simulation!`
- `state_selection!`
- `initial_conditions!`
- `parameters!`
- `graphics!`

Their keywords use the fields documented in the
[Spatial TOML Reference](toml-reference.md). Symbols are converted to strings,
so `mode = :dynamic` and `mode = "dynamic"` are equivalent. Values set by a
later call replace the same fields but retain unrelated settings.

`simulation!` additionally accepts `frames_per_second`. It calculates
`output_samples` from the configured start and end times. The stored model
still contains the resulting sample count and therefore remains independent
of this Julia convenience.

## Elements and graphics

The API has named builders for the current spatial library, including bodies,
markers, joint primitives, compound joints, bushings, contacts, tires, forces,
motion generators, couplers, gears, rack-and-pinion sets, belts, measurements,
friction elements, and user equation components. Builder names follow element
types with a trailing `!`; for example:

```julia
Sim3D.bushing!(model, :mount; markers = [first, second], ...)
Sim3D.rolling_tire!(model, :tire; markers = [center, road], ...)
profile = Sim3D.curve!(model, :profile;
    marker = cam_frame, points = profile_points, half_width = 0.05)
Sim3D.curve_contact!(model, :roller_contact;
    curve = profile, roller_marker = roller_center,
    radius = 0.05, stiffness = 5.0e4)
Sim3D.equation_component!(model, :controller; states = ..., equations = ...)
```

`Sim3D.element!` is the generic form and permits a newly introduced element
type to be used before it receives a named convenience function:

```julia
element = Sim3D.element!(model, :revolute, :pin;
    markers = [body_pin, ground_pin],
    rotation_coordinates = true)
```

All fields still pass through ordinary loader validation.

`Sim3D.graphics!(element; ...)` attaches a simple graphics table.
`Sim3D.graphic!(element, name; ...)` adds one named graphic beneath it, which
allows an element or body to carry several shapes.

## Load, run, and inspect

The exported package entry points accept a Julia-built model directly:

```julia
loaded = load_spatial_model(model)
result = run_spatial_model(model; duration = 2.0, samples = 121)
```

`Sim3D.load_model(model)`, `Sim3D.simulate(model)`, and `Sim3D.run(model)` are
shorter equivalents. The returned `LoadedSpatialModel` and analysis result are
the same types returned for TOML or Lua input.

`Sim3D.document(model)` returns a detached dictionary containing the generated
model document. `load_spatial_model` also accepts an `AbstractDict` directly;
symbolic keys and values are normalized to their portable string form.

```julia
specification = Sim3D.document(model)
loaded = load_spatial_model(specification)
```

The builder can write an inspectable or reusable TOML form without requiring
the user to author TOML:

```julia
Sim3D.write_model("generated-model.toml", model; overwrite = true)
```

The expanded document is embedded in every stored result. Model extraction,
saved-result initialization, CSV conversion, and SimpView therefore behave the
same way for Julia, TOML, and Lua models.

## Current API boundary

The public builder creates model specifications; it deliberately does not
expose allocated equation indices or internal component constructors. Those
remain implementation details and can change without changing the builder.

Force and equation fields currently accept the same portable expression
strings as TOML. A future Julia expression macro can translate Julia syntax
into that representation while retaining automatic differentiation and
result-file reproducibility. Arbitrary Julia closures are not yet accepted as
force laws because they cannot be reconstructed from a saved result.
