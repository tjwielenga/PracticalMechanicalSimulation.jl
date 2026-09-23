# Revolute relative-coordinate pendulum

This model enables the optional rotation coordinate on the ground pin. The
joint then owns the continuous histories `pin.theta`, `pin.omega`, and
`pin.alpha`, defined from the orientation of its first marker relative to its
second marker. The model also requests `pin.omega` as its preferred independent
velocity, so the joint angle and rate are the two integrated physical states.

Before integration, the model performs a static-equilibrium initialization.
Gravity brings the pendulum to its vertically downward configuration. The
declared center-of-mass velocity `[0.5, 0.0]` m/s is perpendicular to gravity,
and `angular_velocity = 1.0` rad/s makes that velocity consistent with the
half-meter distance from the pin to the center of mass. These velocities are
restored after the static solution, so the dynamically consistent pendulum
starts moving horizontally rather than remaining at rest.

Run the model and retain its result with

```bash
./bin/simp2d \
    models/planar/relative-coordinate-pendulum.toml \
    --output results/examples/planar/relative-coordinate-pendulum.simp --overwrite
```

View the mechanism and select any of the joint histories from the plot menus:

```bash
bin/simpView
```

The marker angle offsets participate in the coordinate definition. For the
ordered marker pair $(a,b)$,

$$
\theta_r = \theta_a-\theta_b,\qquad
\omega_r = \omega_a-\omega_b,\qquad
\alpha_r = \alpha_a-\alpha_b.
$$

The relative angle is not wrapped, so it remains a usable integrated coordinate
when a joint passes through a complete revolution.
