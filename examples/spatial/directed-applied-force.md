# Spatial marker-directed applied force

[`directed-applied-force.toml`](../../models/spatial/directed-applied-force.toml)
applies a constant force at an offset body marker. The force direction is the
oriented $z$-axis of `ground.force_axis`. The offset force translates and
rotates the free body while gravity acts at its center of mass.

Run and store the example with:

```bash
./bin/simp3d models/spatial/directed-applied-force.toml \
    --output results/examples/spatial/directed-applied-force.simp --overwrite
```

View the stored result with:

```bash
bin/simpView
```

The named frame graphic shows the direction marker's orientation. The blue
force arrow follows its positive $z$-axis. No reaction body is specified, so
the assumed ground reaction is not drawn.
