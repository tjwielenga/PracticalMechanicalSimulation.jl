# Nonlinear expression pendulum

[`nonlinear-expression-pendulum.toml`](../../models/planar/nonlinear-expression-pendulum.toml)
uses a revolute joint's explicit relative angle and angular velocity in an
applied torque expression:

$$
T=-k_1(\theta-\theta_0)-k_3(\theta-\theta_0)^3-c\omega.
$$

The joint has `rotation_coordinates = true`, so `pin.theta` is continuous and
`pin.omega` is available both to the force law and as the preferred state. The
applied torque owns an explicit `nonlinear_spring.T` variable and its local
constitutive implicit equation. Forward-mode dual numbers supply the two
partials used in that equation's sparse Jacobian row.

Run and view it with:

```sh
./bin/simp2d models/planar/nonlinear-expression-pendulum.toml \
  --output results/examples/planar/nonlinear-expression-pendulum.simp --overwrite
bin/simpView
```

The model includes a link graphic and enables the applied torque symbol. The
equal-and-opposite ground torque is hidden by the default ground-load policy.
