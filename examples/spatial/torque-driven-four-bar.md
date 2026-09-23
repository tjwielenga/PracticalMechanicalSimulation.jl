# Torque-Driven Spatial Four-Bar

[`torque-driven-four-bar.toml`](../../models/spatial/torque-driven-four-bar.toml)
is a deliberately overconstrained spatial model of a planar four-bar. Three
moving bodies contribute 18 body degrees of freedom. Four parallel revolute
joints contribute 20 nominal scalar constraints, although their closed-loop
rank is only 17. The automatic row analysis therefore removes three complete
scalar constraint families and leaves one state.

For the present ordering and scaling, QR deactivates
`coupler_rocker.Phi_xz`, `coupler_rocker.Phi_yz`, and
`crank_ground.Phi_xz`. An analyst might instead expect the program to remove
one $z$-position closure equation and two perpendicular-axis equations at one
joint. These are different valid bases for the same rank-17 constraint system;
the particular QR choice is not part of the model interface and may change if
the equation ordering or scaling changes. The simulated body origins and axes
remain in the original plane.

Because ideal redundant reactions are not unique, the reaction components in
this example are those of the retained independent equation set. A model with
finite connection stiffness would be required to predict how load is
physically shared among the redundant connections.

A constant torque acts on the crank-ground revolute. Gravity is omitted so
the example isolates the closed-loop constraint behavior. The ten-second run
carries the crank through several complete revolutions.

Run and store the result with:

```bash
./bin/simp3d models/spatial/torque-driven-four-bar.toml \
    --output results/examples/spatial/torque-driven-four-bar.simp --overwrite
```

View the result with:

```bash
bin/simpView
```
