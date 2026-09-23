# Spatial Hinge Pendulum

[`hinge-pendulum.toml`](../../models/spatial/hinge-pendulum.toml) combines a
spherical joint with a hinge orientation constraint. The spherical joint
makes the two pin markers coincident. The hinge makes the pendulum marker's
$x$- and $y$-axes perpendicular to the ground marker's $z$-axis. Together,
they leave only rotation around the ground marker's $z$-axis.

The hinge requests its optional rotation coordinates. The preferred-state
setting selects `axis.omega`, so `axis.theta` is the continuous integrated
angle and `axis.alpha` is its angular acceleration.

The link starts horizontally with a small positive angular velocity. Gravity
acts in the negative global $y$ direction, so the link swings in the global
$x$-$y$ plane. The two displayed marker frames make the retained relative
rotation visible.

Run and store it with:

```bash
./bin/simp3d models/spatial/hinge-pendulum.toml \
    --output results/examples/spatial/hinge-pendulum.simp --overwrite
```

View the result with:

```bash
bin/simpView
```
