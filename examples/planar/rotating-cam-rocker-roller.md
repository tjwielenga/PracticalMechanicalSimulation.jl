# Rotating cam and rocker roller follower

[`rotating-cam-rocker-roller.toml`](../../models/planar/rotating-cam-rocker-roller.toml)
places a circular roller at the end of a rocker whose remote end is joined to
ground by a revolute joint. Gravity supplies the closing moment. Static
equilibrium establishes the initial contact force before the cam begins its
prescribed rotation.

Together with the translating roller, translating flat, and rocker flat
examples, this completes the four basic planar cam-follower arrangements. It
checks that the roller contact force is applied at the moving profile point
and produces the proper moment about the rocker pivot.

Run and view it with:

```bash
bin/simp2d models/planar/rotating-cam-rocker-roller.toml \
    --output results/examples/planar/rotating-cam-rocker-roller.simp --overwrite
bin/simpview-web
```

Useful plots include `roller_contact.station`, `roller_contact.gap`,
`roller_contact.normal_force`, `rocker_pivot.theta`, and
`rocker_pivot.omega`.
