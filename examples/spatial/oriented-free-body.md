# Spatial orient constraint

The `orient` element fixes all three relative rotations between two markers
without constraining the distance between them. Two identical bodies make the
effect visible. Both translate and fall together, but the orange free body
tumbles while the blue constrained body remains aligned with the ground
reference frame.

The constrained body's supplied orientation and angular velocity are
intentionally inconsistent. Initial-condition projection aligns its frame with
the reference and removes its relative angular velocity. The free body begins
parallel with the same supplied angular velocity and retains it.

Run and view it with:

```bash
./bin/simp3d models/spatial/oriented-free-body.toml \
    --output results/examples/spatial/oriented-free-body.simp --overwrite
bin/simpView
```
