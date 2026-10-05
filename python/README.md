# Practical Mechanical Simulation Python API

This package builds planar Sim2D and spatial Sim3D models with ordinary Python
functions. It writes the same portable TOML consumed by the existing programs
and invokes `simp2d` or `simp3d`; it does not duplicate either solver.

From a source checkout, install the interface with:

```bash
python3 -m pip install ./python
```

Then import the dimension-specific namespace:

```python
from practical_mechanical_simulation import Sim2D, Sim3D
```

The public conventions mirror the Julia builders. Python omits Julia's `!`
suffix: `Sim2D.rigid_body(...)`, `Sim2D.marker(...)`, and
`Sim2D.simulation(...)` correspond to `Sim2D.rigid_body!(...)`,
`Sim2D.marker!(...)`, and `Sim2D.simulation!(...)`. Element keywords are the
same fields documented in the TOML guides.

See [`docs/common/python-api.md`](../docs/common/python-api.md) for the complete
introduction and the executable examples under `examples/planar` and
`examples/spatial`.

The Python modeling package has no graphics-system dependency. In particular,
it does not import or require Blender. Geometry descriptions remain model data
and are interpreted by SimpView from the resulting `.simp` file.
