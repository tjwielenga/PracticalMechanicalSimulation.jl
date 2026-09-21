# PID-controlled planar pendulum

[`controlled-revolute-pendulum.toml`](../../models/planar/controlled-revolute-pendulum.toml)
demonstrates a user-defined differential and algebraic equation component in
Sim2D. The pendulum starts horizontally. A PID controller moves it to 45
degrees to the right of downward vertical and brings it to rest there.

The controller adds one differential state, the integral of angle error, and
two algebraic variables for angle error and commanded torque. Its torque is
applied by the ordinary planar `applied_torque` element. Controller equations,
mechanism equations, and the force law are solved together in the same sparse
implicit system.

Run and view the example with:

```bash
./bin/simp2d models/planar/controlled-revolute-pendulum.toml \
    --output results/examples/planar/controlled-revolute-pendulum.simp \
    --overwrite
bin/simpview-web
```

Useful plotted quantities include `pin.theta`, `pin.omega`,
`controller.angle_error`, `controller.integral_error`, and
`controller.torque`.
