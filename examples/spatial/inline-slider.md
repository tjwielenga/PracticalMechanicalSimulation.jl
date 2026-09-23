# Spatial Inline Slider

[`inline-slider.toml`](../../models/spatial/inline-slider.toml) constrains the
body's center to the line through `ground.axis` along that marker's $z$-axis.
The body retains its three rotational freedoms while gravity accelerates its
one translational freedom down the line.

The optional translation coordinates are enabled. The preferred physical
states are `guide.velocity` and the body's three angular velocities. The body
also uses the quoted degree form

```toml
orientation = ["30°", 0.0, 0.0, 1.0]
```

which the loader converts to radians.

Run and store it with:

```bash
./bin/simp3d models/spatial/inline-slider.toml \
    --output results/examples/spatial/inline-slider.simp --overwrite
```

View the result with:

```bash
bin/simpView
```
