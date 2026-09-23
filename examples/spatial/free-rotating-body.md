# Free Rotating Spatial Body

[`free-rotating-body.toml`](../../models/spatial/free-rotating-body.toml) is the
first complete spatial model. A free rigid body translates under gravity while
rotating about its body-fixed $z$ principal axis.
Its named ground `xy_frame` graphic is attached to `ground.origin` and
establishes the spatial orientation.

The exact motion is ballistic translation with constant angular speed. The
physical orientation is the rotation matrix evaluated from four normalized
Euler parameters. Three body-fixed pseudo angles accompany the selected
angular-velocity states but do not define finite orientation.

Run and store the model with:

```bash
./bin/simp3d models/spatial/free-rotating-body.toml \
  --output results/examples/spatial/free-rotating-body.simp --overwrite
```

View the stored result with:

```bash
bin/simpView
```
