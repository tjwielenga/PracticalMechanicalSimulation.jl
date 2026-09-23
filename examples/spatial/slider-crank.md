# Spatial Slider-Crank

[`slider-crank.toml`](../../models/spatial/slider-crank.toml) builds a
one-degree-of-freedom slider-crank from two revolute joints and one `inplane`
primitive. The crank is one unit long and the connecting rod is two units
long.

Both revolute axes are parallel to the global $z$-axis, so the links move in
the $x$-$y$ plane. The second marker of `slider_guide` is rotated so that its
$z$-axis points along global $y$. Its plane is therefore the global $x$-$z$
plane. The rod's slider marker is already confined to $z=0$ by the revolute
joints, so the additional inplane equation makes it move along the global
$x$-axis.

The crank starts at $45$ degrees with zero velocity. Gravity acts in the
negative global $y$ direction and causes the mechanism to collapse from its
initial position.

Run and store it with:

```bash
./bin/simp3d models/spatial/slider-crank.toml \
    --output results/examples/spatial/slider-crank.simp --overwrite
```

View the result with:

```bash
bin/simpView
```
