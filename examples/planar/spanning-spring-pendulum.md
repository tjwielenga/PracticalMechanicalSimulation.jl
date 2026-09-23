# Spanning-Spring Pendulum

The maintained model in
[`spanning-spring-pendulum.toml`](../../models/planar/spanning-spring-pendulum.toml)
uses the fully expanded `spanning_force` component with its predefined linear
spring-damper law. It joins the pendulum tip to an offset ground marker.

For marker positions $P_1$ and $P_2$, the element defines

$$
s=P_2-P_1,
\qquad
\ell=\|s\|,
\qquad
u=\frac{s}{\ell},
\qquad
\dot\ell=u^T(V_2-V_1).
$$

The scalar and global forces are

$$
f=-k(\ell-\ell_0)-c\dot\ell,
\qquad
F_1=-fu.
$$

The stored variables `s_x`, `s_y`, `length`, `u_x`, `u_y`, `length_rate`,
`force`, `F_x`, and `F_y` each have a corresponding local implicit equation.
The scalar is the force on the first marker along the line from the second
marker toward the first. Positive force is compression and negative force is
tension. The global components are the force on the first marker, and their
opposite is applied to ground.

Run and view the maintained model with

```bash
./bin/simp2d \
    models/planar/spanning-spring-pendulum.toml \
    --output results/examples/planar/spanning-spring-pendulum.simp --overwrite

bin/simpView
```

The viewer draws the force element between its two markers, and its nine local
histories are available in the plot menus and CSV output.
