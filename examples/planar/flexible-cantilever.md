# Floating-reference flexible cantilever

[`flexible-cantilever.toml`](../../models/planar/flexible-cantilever.toml)
defines one planar two-node Timoshenko beam. Its left generated end marker is
fixed to ground and uniform gravity bends the beam downward. Automatic QR state
selection retains the three elastic velocities and removes the reference-frame
velocities constrained by the support.

Run and save the static result with:

```bash
bin/simp2d models/planar/flexible-cantilever.toml \
  --output results/examples/flexible-cantilever.simp --overwrite
```

Then start SimpView and open the saved result:

```bash
bin/simpView
```

The member appears under **Model → Bodies → beam → Geometry**. Select
**Model → Bodies → beam → Deformation** and move the logarithmic amplification slider
to make the small elastic curve easier to see without changing its rigid
motion or thickness.

For the supplied properties, the computed free-end displacement is about
7.43 mm downward. Elementary uniform-load cantilever theory gives 7.36 mm.
The same model's first linear natural frequency is about 7.22 Hz.
