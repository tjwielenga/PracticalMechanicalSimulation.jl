# Spatial Inplane Slider

[`inplane-slider.toml`](../../models/spatial/inplane-slider.toml) constrains
the `slider.contact` point to the plane through `ground.plane`. The plane
normal is the ground marker's $z$-axis. The body begins with velocity in the
plane and gravity presses it against the plane.

The large translucent marker frame shows the constraint plane. The viewer
also shows the scalar reaction as a force arrow at the constrained point.

Run and store it with:

```bash
./bin/simp3d models/spatial/inplane-slider.toml \
    --output results/examples/spatial/inplane-slider.simp --overwrite
```

View the result with:

```bash
bin/simpView
```
