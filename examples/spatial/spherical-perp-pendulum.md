# Spatial Spherical-and-Perp Pendulum

[`spherical-perp-pendulum.toml`](../../models/spatial/spherical-perp-pendulum.toml)
combines a spherical joint with one perpendicular-axis primitive. The
spherical joint makes the two pin markers coincident. The `perp` constraint
makes the pendulum marker's $x$-axis perpendicular to the ground marker's
$y$-axis, leaving two relative rotational freedoms.

The pendulum begins with angular velocity about its own $x$-axis. The long
axis can swing in the global $x$-$z$ plane while the body also spins around
that axis. The frame at `pendulum.tip` makes the second motion visible.

Run and store it with:

```bash
./bin/simp3d models/spatial/spherical-perp-pendulum.toml \
    --output results/examples/spatial/spherical-perp-pendulum.simp --overwrite
```

View the result with:

```bash
bin/simpView
```
