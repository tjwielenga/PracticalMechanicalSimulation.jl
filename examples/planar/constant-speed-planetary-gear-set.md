# Constant-Speed Planetary Gear Set

The model in
[`constant-speed-planetary-gear-set.toml`](../../models/planar/constant-speed-planetary-gear-set.toml)
demonstrates that a gear pair's geometric carrier does not have to support both
gear bearings. The sun and ring bearings are fixed to ground. The planet
bearing is fixed to a moving carrier. Both contact markers are fixed to that
carrier and therefore follow the planet as it revolves around the common sun
and ring center.

The pitch radii are

$$
r_s=0.4,\qquad r_p=0.3,\qquad r_r=1.0.
$$

The sun--planet contact lies between their centers and gives the external
pair. The tangent constructed from the planet bearing gives two negative
constraint coefficients. Equivalently, its projected signed radii are
$\rho_s=-r_s$ and $\rho_p=r_p$, so its constraint is

$$
-r_s(\theta_s-\theta_c)-r_p(\theta_p-\theta_c)=0.
$$

The ring--planet contact lies outside both centers and gives the internal
pair. Its projected signed radii have the same sign:

$$
r_r(\theta_r-\theta_c)-r_p(\theta_p-\theta_c)=0.
$$

The ring is held fixed, the sun is driven at constant speed, and a $5$ N m
resisting torque is applied to the carrier. The resulting carrier and planet
speeds are

$$
\omega_c=\frac{r_s}{r_s+r_r}\omega_s=\frac{2}{7}\omega_s,
\qquad
\omega_p=-\frac{2}{3}\omega_s.
$$

Each gear-pair element generates one floating contact marker on each of its
two gears. The planet therefore owns one generated marker for the sun contact
and another for the ring contact. Equal-and-opposite reactions are applied at
the moving contact locations. The carrier receives the resulting loads
through the planet bearing rather than directly from the gear-pair elements.
The viewer draws the reaction on the first gear as an arrow tangent to the
pitch circles. The direction is calculated from the planet-bearing and contact
positions. Arrow lengths use one common linear scale, so their relative lengths
show their relative force magnitudes.

The simulation covers one complete carrier revolution. Run it with

```bash
./bin/simp2d \
    models/planar/constant-speed-planetary-gear-set.toml \
    --output results/examples/planar/constant-speed-planetary-gear-set.simp --overwrite
```

Then view the moving carrier, external sun--planet contact, internal
ring--planet contact, and selectable histories with

```bash
bin/simpView
```
