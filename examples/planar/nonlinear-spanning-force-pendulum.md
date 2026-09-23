# Nonlinear Spanning-Force Pendulum

[`nonlinear-spanning-force-pendulum.toml`](../../models/planar/nonlinear-spanning-force-pendulum.toml)
uses the same expanded marker-to-marker geometry as the linear spanning-spring
example, but supplies its scalar force with

$$
f=-k_1(\ell-\ell_0)-k_3(\ell-\ell_0)^3-c\dot\ell.
$$

The TOML expression refers to `spring.length` and `spring.length_rate`. The
component owns the local `spring.force` variable, and dual numbers provide the
constitutive partials for its sparse Jacobian row. The stored force is the load
on the first marker in the $J$-to-$I$ direction. Tension is therefore negative,
and all three terms oppose extension or extension rate.

Run and view the model with:

```bash
./bin/simp2d models/planar/nonlinear-spanning-force-pendulum.toml \
  --output results/examples/planar/nonlinear-spanning-force-pendulum.simp --overwrite
bin/simpView
```

The viewer draws the pendulum, the spanning connector, and the applied and
reaction force arrows. All nine local spanning-force histories are available
in the plot menu.
