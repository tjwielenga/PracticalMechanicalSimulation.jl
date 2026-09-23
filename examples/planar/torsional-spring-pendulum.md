# Torsional-Spring Pendulum

The maintained model in
[`torsional-spring-pendulum.toml`](../../models/planar/torsional-spring-pendulum.toml)
adds a torsional spring-damper between the two orientation markers of the
pendulum's revolute joint. The pendulum starts at $45^\circ$ with zero angular
velocity and oscillates under gravity and the spring torque.

For relative marker angle $\theta_{12}$ and relative angular velocity
$\omega_{12}$, the element defines

$$
T=-k(\theta_{12}-\theta_0)-c\omega_{12}.
$$

Here $k=3\ \mathrm{N\,m/rad}$, $\theta_0=0$, and the local damping time scale
$\tau=0.1\ \mathrm{s}$ supplies

$$
c=\tau k=0.3\ \mathrm{N\,m\,s/rad}.
$$

The stored scalar $T$ is the actual torque on the first marker; the second
marker receives the opposite torque. Thus positive torque follows the first
marker's positive rotational direction, matching the ADAMS convention, while
the minus signs make this built-in law restoring. `pin_spring.T` is an explicit
level-two variable with its own local implicit equation and is stored with the
other result histories.

Run the model and save its result with

```bash
./bin/simp2d \
    models/planar/torsional-spring-pendulum.toml \
    --output results/examples/planar/torsional-spring-pendulum.simp --overwrite
```

Then open the animated pendulum, torsional-spring symbol, and selectable result
histories with

```bash
bin/simpView
```
