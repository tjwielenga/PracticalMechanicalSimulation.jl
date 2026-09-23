# Pendulum Suspended from an Inline

[`inline-pendulum.toml`](../../models/spatial/inline-pendulum.toml) constrains
the pendulum's suspension point to the horizontal line through
`ground.guide`. The guide marker is rotated so its local $z$-axis points along
global $x$.

The pendulum begins vertical. Only its center-of-mass velocity is supplied in
the model, with one component along the guide and one transverse component.
This is not initially consistent because the suspension point would also move
transversely. Initial-condition projection adjusts the translational and
angular velocities to satisfy the Inline constraint. The suspension point
then travels along the guide while the pendulum swings in the $y$-$z$ plane
under gravity.

Run and store the result with:

```bash
./bin/simp3d models/spatial/inline-pendulum.toml \
    --output results/examples/spatial/inline-pendulum.simp --overwrite
```

View the result with:

```bash
bin/simpView
```
