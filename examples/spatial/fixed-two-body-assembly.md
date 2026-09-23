# Spatial fixed joint

The `fixed` element combines the three translational constraints of a
spherical joint with the three rotational constraints of an `orient`.

The blue and orange bodies begin as an L-shaped assembly. Their center
velocities and body-frame angular velocities describe one consistent rigid
motion. The complete assembly falls and tumbles, while the marker points remain
coincident and the relative orientation between the bodies remains fixed.

Run and view it with:

```bash
./bin/simp3d models/spatial/fixed-two-body-assembly.toml \
    --output results/examples/spatial/fixed-two-body-assembly.simp --overwrite
bin/simpView
```
