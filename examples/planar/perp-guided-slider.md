# Perp-Guided Slider

The model in [`perp-guided-slider.toml`](../../models/planar/perp-guided-slider.toml)
demonstrates a planar prismatic joint assembled from two primitive
constraints. The ground guide marker's $x$-axis points horizontally. An
`inplane` constraint keeps the slider center on that axis, while a `perp`
constraint keeps the slider orientation aligned with it.

Let the marker angles be $\theta_i$ and $\theta_j$. The position equation for
the `perp` primitive is

$$
\Phi=(\hat x_i^g)^T\hat y_j^g
    =\sin(\theta_i-\theta_j)=0.
$$

On the perpendicular branch, its equivalent velocity and acceleration
equations are

$$
\Phi_v=\omega_{iz}-\omega_{jz}=0,
\qquad
\Phi_a=\alpha_{iz}-\alpha_{jz}=0.
$$

One scalar reaction torque acts on the first marker and the opposite torque
acts on the second. Together, the primitives constrain vertical translation
and rotation while leaving horizontal translation free. Automatic state
selection therefore chooses `slider.V_x` as the model's single independent
velocity.

The slider starts at $x=-0.5$ with unit positive velocity. It has no applied
force along the guide, so it moves uniformly to $x=0.5$ in one second. Gravity
acts normal to the guide and a constant $+1\ \mathrm{N\,m}$ torque acts on the
slider. The `inplane` reaction supports the slider with $9.81\ \mathrm{N}$,
while the `perp` reaction is $-1\ \mathrm{N\,m}$. Thus the example checks both
primitive reactions without changing the one allowed translation.

Run and view the example with

```bash
./bin/simp2d \
    models/planar/perp-guided-slider.toml \
    --output results/examples/planar/perp-guided-slider.simp --overwrite

bin/simpView
```
