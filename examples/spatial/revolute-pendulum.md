# Spatial Revolute Pendulum

[`revolute-pendulum.toml`](../../models/spatial/revolute-pendulum.toml) uses
one `revolute` element in place of separate spherical and hinge constraints.
Internally, its spherical primitive makes the two pin points coincident and
its two perpendicular-axis primitives leave only rotation around the second
marker's $z$-axis.

The revolute requests its optional rotation coordinates. The preferred state
is `pin.omega`, while `pin.theta` is the continuous integrated joint angle and
`pin.alpha` is its angular acceleration. The link starts horizontally with a
small positive angular velocity and swings in the global $x$-$y$ plane.

Run and store it with:

```bash
./bin/simp3d models/spatial/revolute-pendulum.toml \
    --output results/examples/spatial/revolute-pendulum.simp --overwrite
```

View the result with:

```bash
bin/simpView
```
