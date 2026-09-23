# Planar Bushing

Status: executable TOML force-element example

The planar bushing joins two oriented markers without adding a kinematic
constraint. It supplies two translational force components and one torque. Its
translational deformation and rate are measured in the second marker's frame.
The rotational deformation is the relative marker angle.

For diagonal translational stiffness and damping, the force on marker 1 is

$$
F_1^g=A_2\left[-K_t(d^2-d_0^2)-C_t\dot d^2\right],
$$

and its applied couple is

$$
T_1=-k_r(\theta_1-\theta_2-\theta_0)
    -c_r(\omega_1-\omega_2).
$$

Marker 2 receives the equal-and-opposite force and torque. The component owns
three explicit applied-load variables, `F_x`, `F_y`, and `T`, and three local
constitutive equations. It participates normally in dynamic, acceleration
initialization, and static-equilibrium analyses.

The example specifies the bushing-local value
`damping_time_scale = 0.1`. Therefore its omitted translational and rotational
damping coefficients are calculated from $C_t=0.1K_t$ and $c_r=0.1k_r$,
giving `[10.0, 10.0]` N s/m and `2.0` N m s/rad. Either damping field can still
be supplied explicitly to override its corresponding estimate.

The example model suspends a free planar body from ground. It starts at the
unloaded translational position with an angular displacement, so gravity and
the rotational spring produce damped vertical and rotational responses. Run
and view it with

```bash
./bin/simp2d \
    models/planar/bushing-supported-body.toml \
    --output results/examples/planar/bushing-supported-body.simp --overwrite
bin/simpView
```

Changing `analysis.mode` to `"static"` gives the equilibrium displacement. For
the supplied mass, gravity, and vertical stiffness, the exact center position
is

$$
y=-1-\frac{9.81}{100}=-1.0981\ \mathrm{m},
$$

and the bushing vertical force is $9.81$ N.

The same equilibrium can initialize a dynamic run without changing analysis
time:

```toml
[analysis]
mode = "dynamic"
initialization = "static_equilibrium"
```

With the example's zero declared velocities, the body begins at
$y=-1.0981$ m and remains at rest. If body velocities are supplied, the static
configuration is retained but those velocities are restored before the
dynamic consistency solve and integration.
