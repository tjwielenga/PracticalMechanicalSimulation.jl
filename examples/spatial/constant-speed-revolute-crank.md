# Spatial Constant-Speed Revolute Crank

This example adds a rotational motion generator to a spatial revolute joint.
The generator drives the crank through one complete turn in one second while
gravity acts on it. The mechanism has no remaining state variables: its
position, velocity, and acceleration are determined by the revolute joint and
the prescribed motion at every time.

The generator references the joint rather than repeating its markers:

```toml
[drive]
type = "rotational_motion"
joint = "pin"
function = "constant_speed"
initial_angle = "0 deg"
angular_velocity = "speed"
```

Referencing `pin` automatically enables its relative `theta`, `omega`, and
`alpha` variables. The generator prescribes all three consistently and solves
the required drive torque around the second marker's $z$-axis. In this example
that torque varies through the turn as it balances gravity.

Run and view the example from the repository root:

```bash
./bin/simp3d models/spatial/constant-speed-revolute-crank.toml \
    --output results/examples/spatial/constant-speed-revolute-crank.simp --overwrite
bin/simpView
```

The viewer plots every canonical variable, including `pin.theta`, `pin.omega`,
`pin.alpha`, and `drive.torque`. It draws the drive torque on the crank and,
when ground loads are enabled, the equal-and-opposite reaction on ground.
