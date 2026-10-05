# Python Modeling API

The Python interface constructs the same hierarchical documents accepted from
TOML and produced by the Julia and Lua builders. It then writes portable TOML
and invokes the existing `simp2d` or `simp3d` program. Python input is therefore
not a separate solver path: model validation, initial-condition correction,
state selection, sparse assembly, integration, and result storage remain in
the Julia implementation.

The modeling layer is independent of Blender and other graphics programs.
Graphics entries are ordinary model data. SimpView reads their stored form from
the resulting `.simp` file.

## Install from a source checkout

The initial Python package is distributed with the source repository. It has no
third-party Python runtime dependencies and supports Python 3.10 or newer:

```bash
python3 -m pip install ./python
```

The solver bridge locates `bin/simp2d` and `bin/simp3d` automatically in a
source checkout. An installed interface may instead use executables on `PATH`,
or the environment variables `SIMP2D_EXECUTABLE`, `SIMP3D_EXECUTABLE`, or
`PRACTICAL_MECHANICAL_SIMULATION_ROOT`.

## A planar model

The capitalized module aliases match the Julia `Sim2D` and `Sim3D` namespaces.
Python function names omit Julia's `!` suffix, while model fields retain their
TOML names.

```python
from practical_mechanical_simulation import Sim2D

model = Sim2D.Model("pendulum", title="Pendulum built in Python")
Sim2D.analysis(model, mode="automatic")
Sim2D.simulation(
    model,
    end_time=3.0,
    frames_per_second=60,
    maximum_step=0.01,
)

ground = Sim2D.ground(model, "ground")
ground_pin = Sim2D.marker(ground, "pin")

pendulum = Sim2D.rigid_body(
    model,
    "pendulum",
    mass=1.0,
    inertia=1 / 12,
    position=[0.5, 0.0],
)
body_pin = Sim2D.marker(pendulum, "pin", position=[-0.5, 0.0])
tip = Sim2D.marker(pendulum, "tip", position=[0.5, 0.0])
Sim2D.graphics(pendulum, show_default=False, color="steelblue")
Sim2D.graphic(
    pendulum,
    "member",
    shape="cylinder",
    markers=[body_pin, tip],
    radius=0.04,
)

pin = Sim2D.revolute(
    model,
    "pin",
    markers=[body_pin, ground_pin],
    rotation_coordinates=True,
)
Sim2D.gravity(
    model,
    "gravity",
    acceleration=[0.0, -9.81],
    bodies=[pendulum],
)
Sim2D.state_selection(
    model,
    method="preferred",
    preferred_velocities=[Sim2D.variable(pin, "omega")],
    allow_fallback=False,
)

run = Sim2D.run(model, output="pendulum.simp", overwrite=True)
print(run.stdout)
```

`Sim2D.run` first writes a temporary TOML file beside the model's
`source_directory`, so paths to model resources keep the same meaning. The
temporary file is removed after the solver exits. Pass `model_path=...` to
retain the generated TOML, or call `Sim2D.write_model(...)` without running it.

## Spatial models and reusable assemblies

The spatial interface follows the same pattern:

```python
from practical_mechanical_simulation import Sim3D

model = Sim3D.Model("free_body")
Sim3D.simulation(model, end_time=1.0, frames_per_second=60)
body = Sim3D.rigid_body(
    model,
    "body",
    mass=1.0,
    inertia=[0.1, 0.2, 0.3],
    position=[0.0, 0.0, 1.0],
)
Sim3D.graphics(body, shape="box", size=[0.8, 0.4, 0.2])
Sim3D.run(model, output="free-body.simp", overwrite=True)
```

Ordinary Python functions provide the assembly mechanism. They receive a model
and name prefix, add elements with dotted names, and return whatever handles a
caller needs. Element references passed as keyword values are converted to
qualified TOML names. References from a different model are rejected.

`marker(owner, name, ...)` accepts ground, rigid-body, and flexible-beam
owners. A flexible beam already generates its `end_i`, `cm`, and `end_j`
markers, so `marker(beam, "end_i")` returns that marker's handle without adding
a duplicate table. The equivalent explicit spelling is
`beam_marker(beam, "end_i")`. Arbitrary additional beam stations are not yet
part of the simplified flexible-beam element.

The equivalent executable double-pendulum examples are:

- [`examples/planar/python_api_double_pendulum.py`](../../examples/planar/python_api_double_pendulum.py)
- [`examples/spatial/python_api_double_pendulum.py`](../../examples/spatial/python_api_double_pendulum.py)

Run them from the repository root with:

```bash
python3 examples/planar/python_api_double_pendulum.py
python3 examples/spatial/python_api_double_pendulum.py
bin/simpView
```

The examples can also run directly from a checkout before editable
installation; each adds the local `python/src` directory only when the package
cannot already be imported.

## Builder correspondence

The shared model operations are:

| Purpose | Python | Julia |
|---|---|---|
| Create a model | `Sim2D.Model(...)` | `Sim2D.Model(...)` |
| Add an element | `Sim2D.element(...)` | `Sim2D.element!(...)` |
| Add a marker | `Sim2D.marker(...)` | `Sim2D.marker!(...)` |
| Attach one graphic | `Sim2D.graphic(...)` | `Sim2D.graphic!(...)` |
| Set model or element graphics | `Sim2D.graphics(...)` | `Sim2D.graphics!(...)` |
| Change an element | `Sim2D.set_properties(...)` | `Sim2D.set_properties!(...)` |
| Select a variable | `Sim2D.variable(...)` | `Sim2D.variable(...)` |
| Write TOML | `Sim2D.write_model(...)` | `Sim2D.write_model(...)` |
| Run the model | `Sim2D.run(...)` | `Sim2D.run(...)` |

The named element builders cover the element libraries exported by the Julia
`Sim2D` and `Sim3D` modules. Their keyword fields are intentionally not
reimplemented or interpreted in Python. Use the
[planar TOML reference](../planar/toml-reference.md) and
[spatial TOML reference](../spatial/toml-reference.md) for the authoritative
element definitions; the ordinary loader reports the same validation errors
for Python-generated models as it does for hand-written TOML.

`Sim2D.run` and `Sim3D.run` return a `SolverRun` record containing the command,
exit status, captured standard output and error, and output path. A failed
process raises `SolverError` by default. Use `check=False` when a diagnostic
program needs to inspect an expected failed run itself.
