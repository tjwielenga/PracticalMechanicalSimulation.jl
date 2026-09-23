# Spatial spanning spring-damper

This example suspends a free spatial body from a ground marker with a
spanning spring-damper. The body begins away from the vertical plane and with
a transverse velocity, so the spring both changes length and swings through
three-dimensional space.

Run the model from the project root:

```bash
./bin/simp3d models/spatial/spanning-spring-body.toml \
    --output results/examples/spatial/spanning-spring-body.simp --overwrite
```

View the stored result:

```bash
bin/simpView
```

The scalar is the force on the first marker along the line from the second
marker toward the first. Positive force is compression and negative force is
tension. The applied-force arrow is drawn at the first marker and the
reaction-color arrow at the second marker. Because the second marker belongs
to ground, its arrow is hidden by default.
