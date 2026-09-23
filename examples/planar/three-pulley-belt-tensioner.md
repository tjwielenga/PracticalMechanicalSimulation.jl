# Three-Pulley Belt with a Tensioner

The model in
[`three-pulley-belt-tensioner.toml`](../../models/planar/three-pulley-belt-tensioner.toml)
connects driver and driven pulleys through a third pulley carried by a moving
tensioner arm. The driver and driven bearings connect directly to ground. The
tensioner bearing connects its pulley to the arm, so its two belt spans verify
that tangent geometry does not require a common carrier.

Three ordered elastic spans form the closed belt path. Near-point markers
select the upper driver-to-driven tangent and the two lower tangents adjoining
the tensioner. At each step the span recalculates its tangent points on the
three pitch circles and applies equal-and-opposite tensile forces there.

Straight-span distance alone cannot transmit rotation between fixed-center
pulleys. For each contact, the program therefore projects the pulley surface
velocity onto the span. If $u_1$ and $u_2$ are the two signed surface speeds,

$$
\dot e=u_2-u_1,
\qquad
T=ke+c\dot e.
$$

The assembled geometry provides the material reference, and the belt begins
with a common $100\ \mathrm{N}$ tension. Initial pulley angular velocities are
consistent with the prescribed $2\ \mathrm{rad/s}$ driver speed, avoiding an
artificial damping impulse at the start. Static initialization lets the
tensioner arm settle against its torsional spring and lets the driven pulley
rotate slightly to balance its $1\ \mathrm{N\,m}$ resisting torque. The normal
static minimization retains nearby pulley phases without adding temporary
constraints.

Run and view the model with

```bash
./bin/simp2d \
    models/planar/three-pulley-belt-tensioner.toml \
    --output results/examples/planar/three-pulley-belt-tensioner.simp --overwrite

bin/simpView
```

The viewer reconstructs the pulley pitch discs, three straight spans, and the
three wrap arcs from the stored result. The plot menus include every span's
extension, extension rate, tension, tangent points, and force components.
