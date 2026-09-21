# Rotating cam and circular roller follower

[`rotating-cam-follower.toml`](../../models/planar/rotating-cam-follower.toml)
drives a smooth closed cam beneath a vertically guided circular roller. The
cam profile is a periodic cubic spline through marker-local points. The roller
contact adds an explicit contact station, differentiates its tangency equation
for station rate, and applies a compliant normal force. The cam turns at
$10\ \mathrm{rad/s}$ for one and one-half revolutions. The follower cannot
remain on the retreating profile throughout the cycle, so the example includes
free flight, loss of contact force, and recontact.

Run and view it with:

```bash
bin/simp2d models/planar/rotating-cam-follower.toml \
    --output results/examples/planar/rotating-cam-follower.simp --overwrite
bin/simpview-web
```

Open the result in SimpView. The default contact graphics show the profile,
roller, contact point, and equal-and-opposite normal loads. Useful plots
include `roller_contact.station`, `roller_contact.gap`,
`roller_contact.normal_force`, and `follower.R_y`.

The example uses the built-in force law. Replacing `stiffness`,
`damping_factor`, and `transition_depth` by an `expression` makes the complete
normal law user-defined. Such an expression must explicitly provide any
separation test or force clamp it needs.
