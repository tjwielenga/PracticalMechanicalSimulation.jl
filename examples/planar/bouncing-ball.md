# Bouncing Ball with Compliant Contact

Status: executable TOML contact and DDASSL-discontinuity example

The model drops a planar rigid body under gravity. A sphere of radius $r=0.1$
m is centered on the body's `ball.center` marker and interacts with a
horizontal plane through an ordinary one-sided compliant force,

$$
F_n=k\delta\max(0,1+d v_c),
$$

where

$$
g=(P_1-P_2)^T\hat y_2-r
$$

is positive when the sphere is separated from the plane,
$\delta=\max(-g,0)$, and $v_c=-\dot g$ is positive during approach. The model
uses `damping_factor = 0.15` s/m. Omitting the field or setting it to zero gives
an elastic rebound. No contact
constraint, complementarity condition, or slack variable is introduced. The
body remains part of the same unreduced implicit dynamic system in both
force-law branches.

At every sign change of $g$, the integrator locates the crossing with its BDF
interpolation polynomial, accepts a shortened step ending at the crossing, and
restarts its history at order one. This does not create an impulsive velocity
change: the compliant force is continuous at $g=0$. The restart simply prevents
the high-order history from straddling an abrupt change in stiffness.

Run and view the example with

```bash
./bin/simp2d \
    models/planar/bouncing-ball.toml \
    --output results/examples/planar/bouncing-ball.simp --overwrite
bin/simpView
```

The viewer draws the contact sphere, a rectangular contact plane extending in
the out-of-plane direction, and the surface contact point. Its variable
selectors include `floor_contact.gap`, `floor_contact.gap_rate`, and
`floor_contact.normal_force`, so the force transition can be inspected beside
the animation.
