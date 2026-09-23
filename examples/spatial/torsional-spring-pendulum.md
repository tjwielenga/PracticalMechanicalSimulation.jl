# Spatial torsional-spring pendulum

This example uses one `applied_torque` element on a revolute joint. Its
built-in spring-damper law returns the pendulum from an initial angle of
$30^\circ$ toward a zero free angle. The damping coefficient is inferred from
the local damping time scale: $c=k\tau=8(0.03)=0.24$ N m s/rad.

The torque reference automatically enables `pin.theta`, `pin.omega`, and
`pin.alpha`. The viewer draws the torque on the first joint side in the
applied-load color. Its opposite ground-side torque is hidden by default.

Run and view it with:

```bash
./bin/simp3d models/spatial/torsional-spring-pendulum.toml \
    --output results/examples/spatial/torsional-spring-pendulum.simp --overwrite
bin/simpView
```
