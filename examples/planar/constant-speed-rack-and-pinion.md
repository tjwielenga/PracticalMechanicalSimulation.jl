# Constant-Speed Rack and Pinion

The model in
[`constant-speed-rack-and-pinion.toml`](../../models/planar/constant-speed-rack-and-pinion.toml)
couples a pinion supported by a revolute joint to a rack supported by a
translational joint. Both joints connect their bodies to the same carrier.

The carrier-side marker of the translational joint defines the orthonormal
directions $\hat x_c$ and $\hat y_c$. Rack translation and contact force act
along $\hat x_c$. The positive $\hat y_c$ direction points from the pinion
center $C$ toward the pitch contact, so the specified pitch radius constructs

$$
P_c=C+r_p\hat y_c.
$$

The loader creates a carrier-fixed contact marker and both rack-owned and
pinion-owned floating markers at $P_c$. The two floating markers preserve the
correct contact-force moment arms while the material points at contact change
as the rack translates and the pinion rotates. This remains correct when the
translational joint's marker line is offset from the rack pitch line.

Let $q$ be rack displacement along $\hat x_c$ and let $\theta$ be pinion
rotation relative to the carrier. The ideal no-slip constraint is

$$
q+r_p\theta-\phi=0.
$$

Here $r_p=0.4\ \mathrm{m}$. The model uses `phase = "initial"`; its zero initial
rack displacement and pinion angle therefore give $\phi=0$. The pinion turns
once in five seconds, so the rack travels one pitch circumference in the
opposite signed direction. The rack guide is vertical and gravity loads the
rack downward, providing a nonzero contact reaction and drive torque for
inspection.

Run and view the model with

```bash
./bin/simp2d \
    models/planar/constant-speed-rack-and-pinion.toml \
    --output results/examples/planar/constant-speed-rack-and-pinion.simp --overwrite

bin/simpView
```
