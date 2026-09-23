# Torque-Driven Planar Four-Bar Linkage

Status: verified one-degree-of-freedom closed-loop dynamics

## 1. From prescribed motion to applied load

The earlier driven four-bar prescribed the crank angle, angular velocity, and
angular acceleration. Its motion generator removed the mechanism's remaining
degree of freedom, and a force analysis recovered the torque required to produce
that motion.

This example removes the motion generator and applies a constant torque of

$$
T_a=10\ \mathrm{N\,m}
$$

between ground and the crank. The crank motion is now an outcome of the body
inertia, gravity, applied torque, and ideal joint constraints.

## 2. Constant-torque component

`PlanarConstantTorqueComponent` is an explicit applied-force component. Like the
torsional spring-damper, it acts between two orientation markers and contributes
equal-and-opposite moments. Positive $T_a$ is applied to the first marker and
$-T_a$ to the second marker.

Unlike the spring-damper, the constant torque has no deformation-dependent
constitutive law. Its magnitude is a parameter, so it owns no solution variable
and no equation row. With the second marker fixed to ground, only the crank's
torque balance receives a contribution.

## 3. Dynamic equation selection

Removing the motion generator removes four variables and six kinematic
equations. The three bodies and four revolute joints supply 35 variables but only
33 balance and constraint equations. The remaining two equations are the
selected crank state equations

$$
\alpha_c-\dot\omega_c=0,
\qquad
\omega_c-\dot\theta_c=0.
$$

The resulting 35-equation system retains all three levels of the four ideal-loop
constraints. Only $\omega_c$ and $\theta_c$ are differential variables and
participate in integration-error control. Dependent positions and velocities,
all accelerations, and all eight joint-reaction components remain simultaneous
unknowns in the unreduced solve.

This is a deficit-zero independent-state formulation. The independent state has
two scalar entries, but the Newton system still contains the complete mechanical
solution.

## 4. Initial conditions

The initial crank angle and angular velocity are specified. Consistent initial
conditions are found in three stages:

1. solve the eight position constraints for the other eight configuration
   variables;
2. solve the eight velocity constraints for the other eight velocity variables;
3. solve nine body balances and eight acceleration constraints simultaneously
   for nine accelerations and eight joint reactions.

The applied torque participates in the third stage. It does not affect the
position or velocity consistency calculations.

## 5. Work and energy check

The mechanical energy is

$$
E=\sum_b\left(
\frac12m_b V_b^T V_b+
\frac12J_b\omega_b^2-
m_b(g^g)^T R_b^g
\right).
$$

For constant applied torque, the supplied work is

$$
W_a=T_a\left(\theta_c-\theta_{c0}\right).
$$

The diagnostic monitors

$$
E(t)-E(0)-W_a(t).
$$

Over the default one-second calculation, its maximum magnitude is approximately
$3\times10^{-4}$ J while the crank completes slightly more than one revolution. The
position-constraint error remains near machine precision. The shorter display
interval keeps the rapidly accelerating mechanism visually readable.

Run the dynamic example with

```bash
julia --project=. examples/planar/torque_driven_planar_four_bar.jl
```

Run the current TOML model and common stored-result viewer with

```bash
./bin/simp2d models/planar/torque-driven-four-bar.toml \
    --output results/examples/planar/torque-driven-four-bar.simp --overwrite
bin/simpView
```

The viewer plots the actual integrated motion and makes the stored body,
velocity, applied-load, and reaction variables available as histories.
