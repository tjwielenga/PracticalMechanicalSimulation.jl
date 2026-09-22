# Spatial Floating-Reference Beam

The spatial beam is a straight, two-node Timoshenko member whose local $x$
axis follows its undeformed centerline. Each node begins with the six ordinary
small beam coordinates

$$
d_i=[u_i,v_i,w_i,r_{xi},r_{yi},r_{zi}]^T,
\qquad
d_j=[u_j,v_j,w_j,r_{xj},r_{yj},r_{zj}]^T.
$$

The conventional 12 by 12 element matrices contain six rigid modes. A matrix
$R$ maps the translation and small rotation of the center reference frame into
the two nodal coordinate sets. Starting from six antisymmetric nodal shapes
$S_0$, the implemented elastic shape matrix is

$$
S=S_0-R(R^TMR)^{-1}R^TMS_0.
$$

Thus $R^TMS=0$: rigid reference motion and elastic deformation cannot describe
the same motion. The reduced elastic matrices are

$$
M_e=S^TMS,
\qquad
K_e=S^TKS,
\qquad
C_e=t_dK_e.
$$

Let $R_t$ denote the three rigid-translation columns of $R$. A uniform
gravitational acceleration $g$ has the consistent element load $MR_tg$.
Its direct projection onto the elastic coordinates is

$$
Q_{e,g}=S^TMR_tg=0
$$

because $R^TMS=0$. Gravity therefore appears as the resultant $mg$ in the
floating-reference translational balance rather than as a separate term in
the six elastic equations. Constraints still couple those balances and
produce the correct distributed self-weight deformation. The single gravity
arrow drawn at the reference center is only a viewer convention; it does not
represent a concentrated center load.

The six generalized elastic coordinates are ordered as axial deformation,
two transverse deformations, torsion, and two bending rotations. Their balance
equations are

$$
M_e\ddot\eta+C_e\dot\eta+K_e\eta-Q=0.
$$

The floating reference retains the ordinary spatial body translation, body
angular velocity, body-fixed pseudo angles, and normalized Euler parameters.
Consequently large rigid motion uses the same formulation as every rigid body,
while only the deformation relative to that frame is assumed small.

For an end marker with reference offset $r_0$, translational shape $N_t$, and
rotational shape $N_r$,

$$
r=r_0+N_t\eta,
$$

and its global position is

$$
P=R+A(p)r.
$$

Velocity and acceleration include rigid rotation, elastic motion, and the
Coriolis term. The marker orientation uses the small elastic rotation
$N_r\eta$ relative to the finite floating-frame orientation. A global end
force $f$ contributes its ordinary force and moment to the reference balance
and the generalized elastic load

$$
Q_f=N_t^TA^Tf.
$$

An end torque $m$ similarly contributes

$$
Q_m=N_r^TA^Tm.
$$

Generated `end_i`, `cm`, and `end_j` markers allow fixed, spherical, and
perpendicular-axis constraints and point loads to connect to the member.
Compound joints made from those primitives use the same marker mechanics.
Constraint and reaction Jacobians involving a flexible marker are
differentiated locally with ForwardDiff. Rigid-only models retain their
existing analytical sparse Jacobian path.
