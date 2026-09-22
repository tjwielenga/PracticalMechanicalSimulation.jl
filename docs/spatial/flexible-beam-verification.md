# Flexible-Beam Verification

The floating-reference beam is checked against elementary beam results and
against properties that follow directly from its stated assumptions. These
checks are maintained in
[`spatial_flexible_beam_tests.jl`](../../test/core/spatial_flexible_beam_tests.jl).
They verify the compact beam that Sim3D implements; they do not extend it into
a nonlinear flexible-body formulation.

## Static compliance

The rectangular test beam is one meter long, with local width $b=0.04$ m and
height $h=0.08$ m. Its section properties are

$$
A=bh,qquad I_y=\frac{bh^3}{12},qquad I_z=\frac{hb^3}{12}.
$$

For a cantilever carrying an axial tip force $P$, the expected extension is

$$
u=\frac{PL}{EA}.
$$

With $P=10$ N, the expected extension is $1.5625\times10^{-4}$ m. The test
applies the force through the generated `beam.end_j` marker and compares the
complete static solution with this value.

For a transverse tip force, the Timoshenko tip displacement is

$$
\delta=\frac{PL^3}{3EI}+\frac{PL}{\kappa GA}.
$$

The same $0.1$ N force is applied separately along the local $y$ and $z$
directions. The expected magnitudes are $3.9109375\times10^{-3}$ m about the
weak section axis and $9.8125\times10^{-4}$ m about the strong axis. Checking
both directions detects an interchange of the beam's local section axes.

A component-level torsion check applies an end couple $T$ to the reduced beam
coordinates. The expected twist is

$$
\phi=\frac{TL}{GJ}.
$$

For $T=0.2$ N m, the rectangular-section test gives
$\phi=2.1333485\times10^{-2}$ rad.

## Modal behavior

The circular-section cantilever checks its two equal first bending
frequencies against the Euler--Bernoulli reference

$$
f_1=\frac{1.875104^2}{2\pi}
\sqrt{\frac{EI}{mL^3}}.
$$

The reference value is 7.2243 Hz. The Sim3D frequencies are approximately
7.2176 Hz. The small difference is consistent with the Timoshenko shear term
and the single two-node element. Modal equation errors must also remain below
$10^{-10}$.

## Dynamic behavior

Two dynamic properties are checked:

1. A displaced, damped cantilever must lose mechanical energy. The test
   calculates rigid translational, rigid rotational, elastic kinetic, and
   elastic strain energy from the canonical variables. Energy must decrease
   smoothly and fall below 40 percent of its initial value over the test.
2. A free beam translates while rotating two radians about its local $x$
   axis. Its elastic coordinates must remain below $10^{-9}$. This checks that
   a large rigid motion does not create false elastic deformation in the
   simplified floating-reference formulation.

The existing gravity-loaded cantilever remains as a distributed-load check.
Its static tip displacement is compared with the uniform-load cantilever
result, and a short dynamic release checks the direction and scale of motion.

## Defects found by verification

Adding these cases found two implementation errors:

- the spatial loader rejected generated flexible-beam markers as application
  points for an applied force, even though the force projection already
  supported them; and
- spatial static analysis did not zero elastic acceleration variables before
  evaluating equilibrium. Part of a static load could therefore be balanced
  by a fictitious elastic acceleration.

The same rate classification was corrected in Sim2D's static initialization
and dynamic-relaxation bookkeeping. These cases remain automated regression
tests.

