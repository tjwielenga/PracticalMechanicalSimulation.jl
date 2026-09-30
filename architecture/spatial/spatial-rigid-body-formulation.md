# Spatial Rigid-Body Formulation

## Scope

This document describes the first executable spatial component. It is the
smallest complete test of the intended architecture: a free rigid body under
gravity, represented by local implicit equations and integrated through the
same sparse BDF machinery used by the planar modeler.

## Frames and variables

The body rotation matrix $A^{gb}$ maps body-reference-frame components into
global components. Translation variables are resolved globally. Angular
variables are resolved in the body reference frame. The internal translational
origin is always the center of mass, even when model geometry is entered from
a different body reference origin.

Model input may define this matrix directly or by a right-handed angle
$\theta$ and nonzero axis $u$. With $\hat u=u/\|u\|$, the reader constructs

$$
A=\cos\theta\,I+(1-\cos\theta)\hat u\hat u^T+
\sin\theta\,\widetilde{\hat u}.
$$

The axis has identical components in the before and after frames because it is
unchanged by its own rotation.

Each body allocates 22 canonical variables:

| Level | Variables | Meaning |
|---:|---|---|
| 2 | $a^g$ | translational acceleration |
| 2 | $\alpha^b$ | body-fixed angular acceleration |
| 1 | $V^g$ | translational velocity |
| 1 | $\omega^b$ | body-fixed angular velocity |
| 0 | $R^g$ | center-of-mass position |
| 0 | $\vartheta^b$ | body-fixed pseudo angles |
| 0 | $p=(p_0,e^T)^T$ | scalar-first Euler parameters |

The six selected physical velocities are the components of $V^g$ and
$\omega^b$. Their corresponding integrated coordinates are $R^g$ and
$\vartheta^b$. The Euler parameters separately carry finite orientation and add four
differential variables plus one normalization equation. The pseudo angles are
not finite rotation coordinates. For a selected angular state,

$$
\dot\vartheta^b=\omega^b.
$$

All three pseudo angles remain in the system whether or not their angular
velocities are selected as states. This gives the Newton matrix three local
orientation columns for every body.

## Rotation matrix and orientation rate

For scalar-first Euler parameters $p=(p_0,e^T)^T$, the physical rotation
matrix is

$$
A^{gb}=(p_0^2-e^Te)I+2ee^T+2p_0\widetilde e.
$$

Define

$$
Q(p)=
\begin{bmatrix}
-e^T\\
p_0I+\widetilde e
\end{bmatrix}.
$$

The orientation bridge equations are

$$
\dot\vartheta^b-2Q(p)^T\dot p=0,
\qquad
p^Tp-1=0.
$$

The bridge carries a local pseudo-angle correction into the Euler parameters.
For a selected angular state, the separate equation
$\omega^b-\dot\vartheta^b=0$ connects angular velocity to the bridge. For a
dependent angular component, the position and velocity constraints determine
the pseudo-angle correction and angular velocity instead. The accumulated
value of $\vartheta^b$ has no physical meaning and is never used to calculate the
rotation matrix.

The normalization equation is solved as part of the implicit system. Euler
parameters are also included in BDF error control. The BDF predictor is
normalized before correction so a harmless radial prediction error does not
enter the mechanical equations. The normalization equation remains active in
the corrector.

The component equations still evaluate their finite geometry from $p$. In the
dynamic Newton matrix, their Euler-parameter orientation partials are moved to
the three pseudo-angle columns. At the leading level the coordinate map is

$$
\Delta p=\frac{1}{2}Q(p)\Delta\vartheta^b.
$$

The complete lower-order parameter-rate term remains in the orientation
bridge. Only the orientation bridge and normalization equations retain
Euler-parameter columns. After normalizing $p$, the predictor resets the
nonphysical pseudo-angle prediction so the bridge rate is initially
consistent. The separated matrix then gives the same leading physical
correction as the expanded Euler-parameter matrix.

During initial-condition assembly, orientation is corrected through a
three-component body-frame rotation $\Delta\theta^b$. Let
$\delta=\|\Delta\theta^b\|$ and
$K=\widetilde{\Delta\theta^b}$. The exact incremental rotation is

$$
\Delta A=I+\frac{\sin\delta}{\delta}K+
\frac{1-\cos\delta}{\delta^2}K^2,
$$

and the body orientation is updated on the right:

$$
A_{\mathit{new}}^{gb}=A_{\mathit{old}}^{gb}\Delta A.
$$

The updated matrix is converted back to normalized Euler parameters. Position
assembly therefore uses three independent rotational corrections, preserves a
proper rotation after every Newton step, and does not treat the four Euler
parameters as independent corrections. Body characteristic length scales the
rotational corrections in the weighted minimum-norm solve.

## Balance and state equations

With applied global force $F^g$ and applied body-frame torque $T^b$, the body
balance equations are

$$
m a^g-F^g=0,
$$

$$
J^b\alpha^b+\omega^b\times(J^b\omega^b)-T^b=0.
$$

The state equations retain acceleration as an explicit unknown:

$$
a^g-\dot V^g=0,
\qquad
\alpha^b-\dot\omega^b=0,
$$

$$
V^g-\dot R^g=0,
\qquad
\omega^b-\dot\vartheta^b=0.
$$

Together with the four orientation equations, these give 22 equations for the
22 canonical variables. A free body has six selected physical states: the
three components of $V^g$ and the three components of $\omega^b$.

Gravity contributes $-m g^g$ to the translational balance equations and no
torque. Initial translational acceleration is the sum of the specified gravity
accelerations. Initial angular acceleration is found from the free Euler
equation,

$$
\alpha^b=-{J^b}^{-1}\left(\omega^b\times J^b\omega^b\right).
$$

## Markers

Let $c^b$ be the position of the optional CM marker measured from the modeled
body reference origin. A marker entered at $r_m^b$ is converted once during
model loading to the CM-relative offset

$$
\bar r_m^b=r_m^b-c^b.
$$

If the input body position and velocity describe the reference origin, their
internal CM values are

$$
R_C^g=R_B^g+A^{gb}c^b,
$$

$$
V_C^g=V_B^g+A^{gb}(\omega^b\times c^b).
$$

A body-fixed marker with CM-relative position $\bar r_m^b$ and local
orientation $A^{bm}$ then has

$$
P^g=R_C^g+A^{gb}\bar r_m^b,
\qquad
A^{gm}=A^{gb}A^{bm}.
$$

The CM marker orientation $A^{bc}$ may also define the axes of the entered
CM inertia. The body-balance equations use the constant transformed tensor

$$
J^b=A^{bc}J^c{A^{bc}}^T.
$$

Without a named CM marker, $c^b=0$ and $A^{bc}=I$, which is exactly the
original CM-frame input convention.

The present vertical slice evaluates marker geometry for output and graphics.
These same definitions provide the interfaces for the next spatial constraint
components.

## Sparse assembly and verification

Each body contributes three executable equation blocks: balance, candidate
state equations, and orientation equations. State selection activates only
the state-equation pairs that are needed. Gravity contributes locally to the
body's translational balance rows. The component Jacobian has the BDF form

$$
\frac{\partial G}{\partial z}
+c\frac{\partial G}{\partial\dot z}.
$$

The exact component Jacobian is checked against central finite differences for
all 22 columns. Spatial dynamics then makes the orientation-column
transformation described above. The free-body example is also checked against
exact ballistic translation and constant principal-axis rotation, including
the rotation matrix and Euler-parameter normalization.
