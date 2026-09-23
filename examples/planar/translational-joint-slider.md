# Translational-Joint Slider

The model in
[`translational-joint-slider.toml`](../../models/planar/translational-joint-slider.toml)
uses one `translational` joint to give a rigid body one translational degree of
freedom. Its ordered markers are $i$ on the slider and $j$ on ground.

The joint is a convenience composition of two primitive constraints. Its
generated `guide.inplane` component imposes

$$
(P_i-P_j)^T\hat y_j=0,
$$

and `guide.perp` makes $\hat x_i$ perpendicular to $\hat y_j$. The first
constraint removes relative translation along marker $j$'s local $y$-axis;
the second removes relative rotation. Translation along marker $j$'s local
$x$-axis remains free.

The two primitives retain their separate equations and reactions. In this
example, `guide.inplane.lambda` supports the slider against gravity and
`guide.perp.lambda` balances the applied torque. Automatic state selection
chooses `slider.V_x` and its corresponding position as the single state pair.

The stored-result viewer draws the guide plane and the sliding body. Run the
example with

```bash
./bin/simp2d \
    models/planar/translational-joint-slider.toml \
    --output results/examples/planar/translational-joint-slider.simp --overwrite

bin/simpView
```
