# Spatial cam followers

The spatial cam examples retain the planar profile construction while giving
the contacting geometry a physical width. The profile is fixed in the cam
marker's local $x$-$y$ plane and extruded along local $z$.

[`rotating-cam-roller-follower.toml`](../../models/spatial/rotating-cam-roller-follower.toml)
uses a cylindrical roller. Its marker point is the cylinder center and its
local $z$-axis is the roller axis.

[`rotating-cam-flat-follower.toml`](../../models/spatial/rotating-cam-flat-follower.toml)
uses a plate whose face is the marker's local $x$-$z$ plane. The marker's
positive local $y$-axis gives the contact-force direction.

[`rotating-cam-rocker-roller.toml`](../../models/spatial/rotating-cam-rocker-roller.toml)
places the cylindrical roller at the end of a revolute rocker.
[`rotating-cam-rocker-flat-follower.toml`](../../models/spatial/rotating-cam-rocker-flat-follower.toml)
places the flat plate at the end of the corresponding rocker. These complete
the translating/rocker and roller/flat set of spatial examples.

Both examples use an inline and an orient constraint to leave only vertical
translation. Those joints also preserve the required parallel cam and
follower $z$-axes; the contact elements themselves remain compliant force
elements and do not add an alignment constraint.

Run and view them with:

```bash
bin/simp3d models/spatial/rotating-cam-roller-follower.toml \
    --output results/examples/spatial/rotating-cam-roller-follower.simp --overwrite
bin/simp3d models/spatial/rotating-cam-flat-follower.toml \
    --output results/examples/spatial/rotating-cam-flat-follower.simp --overwrite
bin/simp3d models/spatial/rotating-cam-rocker-roller.toml \
    --output results/examples/spatial/rotating-cam-rocker-roller.simp --overwrite
bin/simp3d models/spatial/rotating-cam-rocker-flat-follower.toml \
    --output results/examples/spatial/rotating-cam-rocker-flat-follower.simp --overwrite
bin/simpView
```

Useful plots include `roller_contact.station`, `roller_contact.gap`,
`roller_contact.normal_force`, and the corresponding `flat_contact`
variables.
