# Spatial Constant-Distance Massless Link

This example suspends a freely rotating body from ground with a spanning
motion generator. The generator prescribes the positive distance between its
two markers:

$$
\ell=\lVert P_2-P_1\rVert=1.
$$

```toml
[link]
type = "spanning_motion"
markers = ["body.link", "ground.anchor"]
distance = 1.0
```

Because only the distance is constrained, each end behaves like a spherical
joint. With constant distance the generator is an ideal massless rigid link.
The body starts with velocity transverse to the link, so it swings in three
dimensions while also rotating freely about its center of mass.

Run and view the example from the repository root:

```bash
./bin/simp3d models/spatial/constant-distance-massless-link.toml \
    --output results/examples/spatial/constant-distance-massless-link.simp --overwrite
bin/simpView
```

The viewer draws the massless link and the force it applies to the body. The
equal-and-opposite ground reaction remains hidden unless ground loads are
enabled.
