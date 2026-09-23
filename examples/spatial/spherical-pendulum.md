# Spatial Spherical-Joint Pendulum

[`spherical-pendulum.toml`](../../models/spatial/spherical-pendulum.toml)
connects one horizontal body to ground with a spherical joint. The joint makes
the two marker points coincident while leaving all three relative rotations
free. Gravity starts the body swinging in the global $x$-$z$ plane.
The named ground `xy_frame` graphic is attached to `ground.pin`. A second,
smaller frame follows `pendulum.tip`, demonstrating the position and
orientation of a body-fixed marker.

Run and store it with:

```bash
./bin/simp3d models/spatial/spherical-pendulum.toml \
    --output results/examples/spatial/spherical-pendulum.simp --overwrite
```

View the result with:

```bash
bin/simpView
```
