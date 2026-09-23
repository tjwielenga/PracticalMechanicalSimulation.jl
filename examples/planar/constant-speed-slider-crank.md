# Constant-Speed Slider-Crank Mechanism

Status: verified kinematic mechanism

## 1. Mechanical model

The slider-crank contains two planar rigid bodies. A crank of length $r$ is
connected to ground by a revolute joint. A connecting rod of length

$$
\ell=2r
$$

is connected to the crank by a second revolute joint. The marker at the other
end of the rod is connected to ground by an in-plane constraint.

The ground point and the in-plane reference point both lie at the crank's fixed
pivot. The ground-fixed in-plane direction is

$$
\hat u^g=\begin{bmatrix}0&1\end{bmatrix}^T.
$$

Consequently, the constraint suppresses vertical displacement of the rod-end
marker while allowing it to slide horizontally along the line through the crank
pivot. Its scalar reaction is vertical; no horizontal slider reaction is
introduced.

The default dimensions are $r=0.5$ m and $\ell=1.0$ m.

## 2. Constant-speed drive

A rotational motion generator between the crank and ground prescribes

$$
\theta_c(t)=\theta_{c0}+\omega_c t,
\qquad
\dot\theta_c=\omega_c,
\qquad
\ddot\theta_c=0.
$$

The default angular velocity is $2\pi/5$ rad/s, so the five-second example
contains one complete crank revolution. The generator's reaction variable is
the driving torque required to overcome inertia and gravity.

## 3. Equation inventory

The canonical model contains 27 variables and 27 available equations:

- 18 rigid-body variables and 6 body-balance equations;
- 4 revolute-joint reactions and 12 revolute-constraint equations;
- 1 in-plane reaction and 3 in-plane constraint equations; and
- 4 motion-generator variables and 6 motion-generator equations.

The metadata policies select the following square systems:

| Analysis | Unknowns | Equations |
|---|---:|---:|
| Kinematic position | 7 | 7 |
| Kinematic velocity | 7 | 7 |
| Kinematic acceleration | 7 | 7 |
| Kinematic forces | 6 | 6 |

The force analysis determines the two ground-pin reactions, two crank-pin
reactions, the in-plane reaction, and the required driving torque from the six
body force and torque balances.

## 4. Predictor and correction

Only the first configuration is seeded from the elementary slider-crank
geometry. Subsequent configurations use the preceding positions, velocities,
and accelerations to form the second-order prediction

$$
R_{n+1}^{(0)}=R_n+hV_n+\frac{h^2}{2}a_n,
$$

$$
\theta_{n+1}^{(0)}=\theta_n+h\omega_n+\frac{h^2}{2}\alpha_n.
$$

Position analysis corrects this prediction onto the assembled constraint
manifold. Over the complete default cycle, no sampled position requires more
than two Newton corrections.

## 5. Analytical check

For the right-hand assembly branch, the slider position is

$$
x_s=r\cos\theta_c+
\sqrt{\ell^2-r^2\sin^2\theta_c}.
$$

The regression test compares this expression against the component-assembled
solution at 81 positions around the complete revolution. It also checks the rod
length and zero vertical slider displacement directly.

Run the analysis with

```bash
julia --project=. examples/planar/constant_speed_slider_crank.jl
```

Save and view the TOML-defined simulation with

```bash
./bin/simp2d \
    models/planar/constant-speed-slider-crank.toml --output results/examples/planar/constant-speed-slider-crank.simp
bin/simpView
```

The viewer reads only the saved result. It reconstructs a ground guide, the two
moving links, and the three connection locations. Every stored canonical
variable is available in the history menus.
