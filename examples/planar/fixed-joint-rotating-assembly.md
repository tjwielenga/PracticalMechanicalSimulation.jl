# Fixed-Joint Rotating Assembly

The model in
[`fixed-joint-rotating-assembly.toml`](../../models/planar/fixed-joint-rotating-assembly.toml)
connects two rigid bodies with

```toml
[rigid_connection]
type = "fixed"
markers = ["carrier.tip", "extension.base"]
```

This user-facing joint introduces no new constraint equations. The loader
constructs two named primitives:

- `rigid_connection.revolute` makes the marker points coincide and owns the
  two force reactions.
- `rigid_connection.perp` makes the marker orientations agree and owns the
  torque reaction.

The first body is connected to ground by a revolute joint and driven through
one complete rotation in five seconds. The fixed connection makes the second
body rotate with it as one rigid assembly. Gravity and centripetal acceleration
produce nonzero primitive reactions, all of which remain available in the
stored result.

Run and view the example with

```bash
./bin/simp2d \
    models/planar/fixed-joint-rotating-assembly.toml \
    --output results/examples/planar/fixed-joint-rotating-assembly.simp --overwrite

bin/simpView
```
