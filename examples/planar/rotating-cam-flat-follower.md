# Rotating cam and flat plate follower

[`rotating-cam-flat-follower.toml`](../../models/planar/rotating-cam-flat-follower.toml)
drives the same smooth closed cam used by the circular-roller example beneath
a vertically guided flat plate. The plate is defined by `follower.face`: its
local $x$-axis lies along the face and its local $y$-axis is the positive force
direction on the follower.

The cam turns at $12\ \mathrm{rad/s}$ for one and one-half revolutions. The
contact station solves the curve-to-plate tangency equation explicitly and is
allowed to pass continuously through the periodic curve boundary. At this
speed the follower loses contact, enters free flight, and later strikes the
cam again.

Run and view it with:

```bash
bin/simp2d models/planar/rotating-cam-flat-follower.toml \
    --output results/examples/planar/rotating-cam-flat-follower.simp --overwrite
bin/simpview-web
```

Useful plots include `flat_contact.station`, `flat_contact.gap`,
`flat_contact.normal_force`, and `follower.R_y`.
