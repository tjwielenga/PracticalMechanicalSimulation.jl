# Translational-Distance Slider

The maintained model in
[`translational-distance-slider.toml`](../../models/planar/translational-distance-slider.toml)
demonstrates the one-direction translational distance generator. For ordered
markers $i$ and $j$, it prescribes

$$
d=(P_i-P_j)^T\hat y_j,
$$

where the second marker's local $y$-axis defines the positive distance
direction. The generator does not constrain motion parallel to its reference
plane or relative marker orientation.

In this example the drive-axis marker's $y$-axis points vertically. An ordinary
`inplane` primitive uses a second marker rotated by $90^\circ$ to prevent
horizontal motion, while a `perp` primitive prevents rotation. Together with
the distance generator, these three scalar constraints prescribe the slider
completely. The generator moves it upward at $0.4\ \mathrm{m/s}$:

$$
d(t)=0.4t,\qquad \dot d=0.4,\qquad \ddot d=0.
$$

Gravity provides the opposing load, so the generator reports a constant
positive reaction `drive.force_m` of $9.81\ \mathrm{N}$. The explicit
coordinate histories are `drive.distance_m`, `drive.velocity_m`, and
`drive.acceleration_m`.

Run the model and save its result with

```bash
./bin/simp2d \
    models/planar/translational-distance-slider.toml \
    --output results/examples/planar/translational-distance-slider.simp --overwrite
```

Then view the slider, the generator's reference plane and distance indicator,
and the selectable histories with

```bash
bin/simpView
```

The same signed-distance kinematics can later provide the translational
coordinate used by a rack-and-pinion constraint.
