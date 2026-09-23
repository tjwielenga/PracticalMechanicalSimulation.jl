# Constant-Speed External Gear Pair

The model in
[`constant-speed-gear-pair.toml`](../../models/planar/constant-speed-gear-pair.toml)
contains two gears supported by revolute bearings on a fixed carrier. The
bearing centers are one unit apart. Their pitch radii are $r_1=0.4$ and
$r_2=0.6$, as calculated from the two base-side bearing markers to the
carrier contact marker at $(0.4,0)$. Because the contact lies between the
bearings, the signed radii used in the constraint have opposite signs. The
loader therefore identifies an external pair.

The gear-pair element generates one floating contact marker on each gear. The
two floating points follow the carrier contact marker, while their orientations
and applied-force ownership remain with their respective gears. The element
names both revolute joints and the carrier contact marker. It imposes

$$
\rho_1(\theta_1-\theta_c)-\rho_2(\theta_2-\theta_c)=0,
\qquad \rho_1=0.4,\quad \rho_2=-0.6.
$$

The carrier is ground in this model, and a constant-speed generator drives
gear 1 through one complete revolution in five seconds. A constant
$+5\ \mathrm{N\,m}$ torque acts on gear 2 relative to the carrier. Since gear 2
rotates in the negative direction, this is a resisting load. Consequently,

$$
\theta_2=-\frac{r_1}{r_2}\theta_1=-\frac{2}{3}\theta_1.
$$

The gear pair owns one scalar reaction. Its force on the first gear is directed
along the tangent calculated from the bearing and contact positions; the second
gear receives the equal-and-opposite force. Applying these forces at the common
contact point produces the gear moments and bearing reactions. The contact
marker therefore needs no separately specified force direction.

The constant-speed motion has zero angular acceleration, so moment balance on
gear 2 gives the expected gear-pair reaction magnitude

$$
|\lambda|=\frac{5}{0.6}=8.333\ldots\ \mathrm{N}.
$$

The required driving torque on gear 1 is therefore

$$
T_d=0.4|\lambda|=3.333\ldots\ \mathrm{N\,m},
$$

which also satisfies power balance:
$T_d\omega_1+5\omega_2=0$.

Run the model and save its result with

```bash
./bin/simp2d \
    models/planar/constant-speed-gear-pair.toml \
    --output results/examples/planar/constant-speed-gear-pair.simp --overwrite
```

Then open the two gear discs, rotating spokes, bearing points, contact point,
and selectable histories with

```bash
bin/simpView
```
