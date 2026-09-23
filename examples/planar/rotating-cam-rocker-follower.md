# Rotating cam and revolute rocker follower

[`rotating-cam-rocker-follower.toml`](../../models/planar/rotating-cam-rocker-follower.toml)
places the flat follower face at the end of a rocker whose remote end is joined
to ground by a revolute joint. Gravity supplies the closing moment. Static
equilibrium establishes the initial contact force before the cam begins its
prescribed rotation.

This example checks behavior that the translating follower does not exercise:
the follower face rotates, the tangent point moves along the plate, and the
normal force produces a changing moment about the rocker pivot. The contact
station remains an explicit, unwrapped profile coordinate.

Run and view it with:

```bash
bin/simp2d models/planar/rotating-cam-rocker-follower.toml \
    --output results/examples/planar/rotating-cam-rocker-follower.simp --overwrite
bin/simpView
```

Useful plots include `flat_contact.station`, `flat_contact.gap`,
`flat_contact.normal_force`, `rocker_pivot.theta`, and
`rocker_pivot.omega`.
