# Bushing-supported pendulum

Status: executable TOML static-to-dynamic example

This model replaces the ideal revolute joint in the relative-coordinate
pendulum with an ordinary planar bushing. The bushing acts between the body pin
marker and the ground origin. Its translational coefficients are

$$
K_t=\operatorname{diag}(1000,1000)\ \mathrm{N/m},\qquad
C_t=\operatorname{diag}(20,20)\ \mathrm{N\,s/m}.
$$

Both rotational coefficients are explicitly zero:

```toml
rotational_stiffness = 0.0
rotational_damping = 0.0
```

Consequently the bushing locates the pin compliantly but supplies no direct
restoring or damping torque. Unlike an ideal revolute joint, it adds no
constraint and the planar body retains three independent coordinates.

The analysis selects `static_method = "dynamic_relaxation"`. It advances the
same dynamic equations with model time held at zero, starting with zero
velocity. For this model one default pseudo-time interval moves the initially
angled body into the static Newton corrector's convergence region; five final
Newton corrections produce the equilibrium. The relaxation history is not
part of the physical simulation or its result history.

The static initialization aligns the pendulum with gravity. At equilibrium the
vertical bushing extension is

$$
\delta_y=\frac{mg}{k_y}=\frac{9.81}{1000}=0.00981\ \mathrm{m},
$$

so the pin lies slightly below the ground origin. The declared horizontal
center-of-mass velocity of $0.5$ m/s and angular velocity of $1$ rad/s make the
pin velocity initially zero. They are restored after static equilibrium and
launch the pendulum into its dynamic motion.

Run and view the example with

```bash
./bin/simp2d \
    models/planar/bushing-pendulum.toml \
    --output results/examples/planar/bushing-pendulum.simp --overwrite
bin/simpView
```
