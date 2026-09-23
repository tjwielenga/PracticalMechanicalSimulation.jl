# Constant-Speed Internal Gear Pair

The model in
[`constant-speed-internal-gear-pair.toml`](../../models/planar/constant-speed-internal-gear-pair.toml)
places the ring-gear bearing at $(0,0)$, the pinion bearing at $(1,0)$, and the
carrier contact marker at $(1.4,0)$. The contact therefore lies outside the
two bearing centers. The loader obtains physical pitch radii of $1.4$ and
$0.4$ from the carrier geometry. The equal-and-opposite contact forces give
the signed radii

$$
\rho_1=1.4,\qquad \rho_2=0.4.
$$

Their matching signs identify an internal pair. Its constraint is

$$
\rho_1(\theta_1-\theta_c)-\rho_2(\theta_2-\theta_c)=0.
$$

The gear pair generates one floating contact marker on each gear. Both markers
follow the carrier contact point while applying the equal-and-opposite reaction
forces to their respective gear bodies. Their direction is tangent to the pitch
circles, as calculated from the bearing and contact positions.

The fixed carrier gives $\omega_2=3.5\omega_1$: the ring and pinion rotate in
the same direction. A generator drives the ring at $0.6\ \mathrm{rad/s}$ and a
constant $+5\ \mathrm{N\,m}$ torque acts on the pinion. The resulting generator
torque is $-17.5\ \mathrm{N\,m}$, satisfying

$$
(-17.5)(0.6)+(5)(2.1)=0.
$$

Run the model and save its result with

```bash
./bin/simp2d \
    models/planar/constant-speed-internal-gear-pair.toml \
    --output results/examples/planar/constant-speed-internal-gear-pair.simp --overwrite
```

View the ring, pinion, bearings, contact point, and result histories with

```bash
bin/simpView
```

The general stored-result viewer draws the larger internal member with a
lightly shaded, outlined disc so that the enclosed pinion remains visible.
