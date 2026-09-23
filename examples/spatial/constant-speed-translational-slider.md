# Spatial Constant-Speed Translational Slider

This example combines an inline constraint and an orient constraint to make a
spatial translational joint. A translational motion generator then drives its
remaining freedom upward at constant speed while gravity acts downward.

The ordered generator markers define the signed distance

$$
d=(P_i-P_j)^T\hat z_j.
$$

The first marker therefore moves normal to the plane through the second
marker. The second marker's local $z$-axis gives the positive direction.

```toml
[drive]
type = "translational_motion"
markers = ["slider.axis", "ground.axis"]
function = "constant_speed"
initial_distance = 0.2
velocity = "speed"
```

The complete mechanism has no selected states. The generator supplies the
position, velocity, and acceleration at every time, while its solved drive
force balances gravity. The viewer shows that force on the slider and hides
the equal-and-opposite ground reaction unless ground loads are enabled.

Run and view the example from the repository root:

```bash
./bin/simp3d models/spatial/constant-speed-translational-slider.toml \
    --output results/examples/spatial/constant-speed-translational-slider.simp --overwrite
bin/simpView
```
