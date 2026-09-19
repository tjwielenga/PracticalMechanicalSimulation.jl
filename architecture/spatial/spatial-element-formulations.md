# Spatial Element Formulations

## Initial-condition correction

Before state selection, the supplied body configuration and velocity guesses
are projected onto the position- and velocity-level constraint equations. The
configuration correction for body $b$ consists of a global translation
$\Delta R_b$ and a small body-frame rotation $\Delta\phi_b$. The rotation is
applied by right multiplication,

$$
A_b^+ = A_b\exp([\Delta\phi_b]_\times),
$$

and converted back to normalized Euler parameters. Consequently, the
correction never treats the four Euler parameters as independent physical
coordinates.

Let $L_b$ be the body's characteristic length and $s_b$ its user-supplied
initial-condition scale. With positive mass, the default body weight is

$$
w_b=s_bm_b.
$$

The position projection chooses the Newton correction having minimum weighted
size

$$
\sum_b w_b\left(
\Delta R_b^T\Delta R_b+
L_b^2\Delta\phi_b^T\Delta\phi_b
\right).
$$

The coefficient of $\Delta R_b^T\Delta R_b$ is therefore $s_bm_b$, while the
coefficient of $\Delta\phi_b^T\Delta\phi_b$ is
$s_bm_bL_b^2$. The latter is a scalar approximation to rotational inertia; it
is not the body's full inertia tensor. The velocity projection uses the same
coefficients for $\Delta V_b$ and $\Delta\omega_b$.

$L_b$ is normally the largest nonzero body-frame marker offset. If the body
has no offset marker, it is estimated from the largest diagonal inertia and
mass as $\sqrt{J_{\max}/m_b}$, with a fallback lower bound of one model length
unit. This mass weighting makes lighter bodies absorb more of an otherwise
ambiguous assembly correction. A future zero-mass spatial body uses a small
positive effective mass based on the smallest positive body mass in the
model, preventing a singular numerical weight while leaving it easy to move.

Imposed translational components are omitted from the correction columns. A
complete imposed orientation omits all three small-rotation columns for that
body. Relative-coordinate definitions are normally evaluated after body
assembly. They join the weighted projection only when the user supplies a
relative initial value, weight, or imposed value.

## Directed distance

Several spatial elements use the same marker geometry. Select one of the
second marker's oriented axes and call it $\hat n_j$. The signed distance from
the plane through the second marker to the first marker is

$$
d=(P_i-P_j)^T\hat n_j.
$$

Its derivatives are

$$
\dot d=(V_i-V_j)^T\hat n_j+(P_i-P_j)^T\dot{\hat n}_j,
$$

$$
\ddot d=(a_i-a_j)^T\hat n_j+
2(V_i-V_j)^T\dot{\hat n}_j+
(P_i-P_j)^T\ddot{\hat n}_j.
$$

The Inplane uses the $z_j$ directed distance as a constraint. Translational
motion prescribes the same distance. The Inline uses the $x_j$ and $y_j$
distances as constraints and, when requested, reports the $z_j$ distance as
its free relative coordinate. An applied force uses only the directed axis.
These elements share the axis, distance, derivative, reaction, and analytical
Jacobian calculations.

A standalone `directed_distance` measurement uses the second marker's
$z_j$-axis and adds three local definitions,

$$
d_m-d=0,\qquad v_m-\dot d=0,\qquad a_m-\ddot d=0.
$$

It contributes no body load or reaction. Its three variables and equations
therefore leave the model mobility unchanged. The distance and velocity are
available to scalar force expressions, while all three levels are stored as
outputs. Measurement variables are deliberately excluded from state selection
and initial-condition correction.

## Applied force

The applied force references an application marker $a$ and a direction marker
$d$. Its global direction is the direction marker's oriented $z$-axis,

$$
\hat d^g=A_d^g\hat e_z.
$$

The element retains four explicit load variables and equations,

$$
f-F(t,z)=0,
\qquad
F^g-f\hat d^g=0.
$$

Positive $f$ therefore acts along $+\hat d^g$. The global force $F^g$ is
applied at the first marker. If its body-fixed offset is $r_a^b$, the force
also gives the body-frame moment

$$
T_a^b=r_a^b\times (A^{gb})^TF^g.
$$

The direction marker may be fixed to ground or any moving body. Its position
does not affect the force; only its orientation is used.

An optional reaction body receives $-F^g$ at a generated floating marker that
remains coincident with the application point. For reaction body $r$, its
instantaneous body-frame lever is

$$
r_f^r=(A_r^{gr})^T(P_a^g-R_r^g).
$$

The reaction moment is $r_f^r\times (A_r^{gr})^T(-F^g)$. Thus the two forces
are equal, opposite, and coincident even though the reaction point moves in
the reaction body's frame. Without a reaction body, ground is assumed and no
opposite body force is assembled.

The scalar law may be constant or a restricted model expression. The
constitutive, direction, application-point, and floating-reaction partials
are all supplied to the analytical sparse Jacobian.

## Spatial bushing

The spatial bushing joins oriented markers $i$ and $j$ without adding a
constraint. Marker $j$ defines the frame used for its three translational and
three rotational laws. Its local translation is

$$
r^j=A_j^T(P_i-P_j).
$$

Because the components are measured in a rotating frame, their exact rate is

$$
v^j=A_j^T(V_i-V_j)-\omega_j^j\times r^j,
$$

where $\omega_j^j=A_j^T\Omega_j$. The cross-product term is required. Without
it, a common rigid rotation of the two markers would appear to change the
bushing deformation.

Let the relative orientation be

$$
A_r=A_j^TA_i=
R_z(\alpha_z)R_y(\alpha_y)R_x(\alpha_x).
$$

The bushing uses the corresponding $z$-$y$-$x$ Bryant angles,

$$
\alpha_x=\operatorname{atan2}(A_{r,32},A_{r,33}),
$$

$$
\alpha_y=\operatorname{atan2}\left(
-A_{r,31},\sqrt{A_{r,32}^2+A_{r,33}^2}\right),
$$

$$
\alpha_z=\operatorname{atan2}(A_{r,21},A_{r,11}).
$$

This sequence allows large rotation around the second marker's $z$-axis while
the $x$ and $y$ rotations remain moderate. A pure $z$ rotation leaves
$\alpha_x=\alpha_y=0$ through a complete turn. The coordinate singularity is
at $\alpha_y=\pm90^\circ$, outside the intended bushing range. When every
rotation is small, the three angles reduce to the ordinary small rotation
components.

Each angle equation uses a periodic difference between its explicit variable
and the angle calculated from $A_r$. The same construction is used for the
hinge angle. The integration history and predictor can therefore retain a
continuous $\alpha_z$ through complete turns instead of forcing it to jump at
$\pm\pi$.

Rotational damping uses the relative angular velocity resolved in the same
frame,

$$
\omega_{ij}^j=A_j^T(\Omega_i-\Omega_j).
$$

For diagonal stiffness and damping vectors, the local constitutive laws are

$$
f^j=-K_t r^j-C_t v^j,
\qquad
\tau^j=-K_r\alpha-C_r\omega_{ij}^j.
$$

The global loads on marker $i$ are

$$
F^g=A_jf^j,
\qquad
T^g=A_j\tau^j.
$$

The second body receives the opposite wrench at a generated floating marker
that follows $P_i$. The force pair is therefore equal, opposite, and
coincident. Together with the opposite torque pair, this introduces no
unintended net wrench.

Translation, translation rate, Bryant angles, relative angular velocity,
local loads, and global loads are explicit component variables. They are tied
to the marker kinematics and constitutive laws by 24 local implicit equations.
Marker-point and marker-axis partials provide the kinematic derivatives. Those
partials and the load-transfer derivatives are inserted directly into the
sparse system Jacobian.

Setting $K_{r,z}=C_{r,z}=0$ while making the other five stiffness directions
large gives a compliant approximation to a revolute joint around $z_j$. It
does not reduce the model mobility: force and moment balance determine the five
small restrained deflections. This is useful when physical joint compliance is
wanted or when an ideal closed loop is opened by a stiff element. The required
stiffness should be chosen from the acceptable joint deflection rather than
made arbitrarily large.

## Plane contact

The plane-contact element places a sphere of radius $r$ at marker $i$ and
uses the local $z$-axis of marker $j$ as the outward normal $\hat n$ of a
plane through $P_j$. The signed gap and its rate are

$$
g=(P_i-P_j)^T\hat n-r,
$$

$$
\dot g=(V_i-V_j)^T\hat n+(P_i-P_j)^T\dot{\hat n}.
$$

Positive gap denotes separation. With penetration

$$
\delta=\max(-g,0).
$$

For the normalized penetration $x=\delta/e$, use the smooth stiffness ramp

$$
S(x)=126x^5-420x^6+540x^7-315x^8+70x^9
$$

and its integral

$$
I(x)=21x^6-60x^7+\frac{135}{2}x^8-35x^9+7x^{10}.
$$

With optional transition depth $e$, define

$$
\delta_e=
\begin{cases}
\delta, & e=0,\\
eI(\delta/e), & 0<\delta<e,\\
\delta-\dfrac{e}{2}, & \delta\ge e.
\end{cases}
$$

The compliant normal law is

$$
F_n=k\delta_e\max(0,1-d\dot g).
$$

For positive $e$, force and its first five penetration derivatives are
continuous at both ends of the transition. The tangent rises monotonically
from zero to $k$ and remains $k$ afterward. Setting $e=0$ recovers the
unsmoothed linear law.

Thus closing motion increases the contact force, rebound reduces it, and the
element cannot produce tension. The global force on the sphere body is

$$
F^g=F_n\hat n.
$$

Alternatively, an input expression may define $F_n$ from the registered model
variables, including this element's $g$ and $\dot g$. In that form the value of
the expression is used without modification. The plane-contact element adds
no separation test, nonnegative clamp, or smoothing; those are part of the
analyst's expression when wanted. This permits ordinary linear viscous
damping, nonlinear bump-stop curves, or another contact law without attaching
external Julia code. Forward-mode differentiation supplies its local
Jacobian partials.

The point

$$
Q=P_i-\left((P_i-P_j)^T\hat n\right)\hat n
$$

is the projection of the sphere center onto the plane. The sphere body
receives $F^g$ at $Q$, and a non-ground plane body receives $-F^g$ at the same
point. Moving the sphere-side application point from $Q$ to the sphere center
would give the same moment because the displacement is parallel to the force;
using the common point makes the equal-and-opposite pair explicitly
coincident.

Gap, gap rate, normal force, and the three global-force components are
explicit variables with six component-local implicit equations. Local
forward-mode differentiation supplies the point, normal, and load-transfer
partials. Because the compliant force remains continuous at its switching
surfaces, DDASSL crosses them through ordinary correction and error control.
No event location or forced history restart is used.

## Rolling tire

The rolling tire is a force element rather than a no-slip constraint. Marker
$w$ lies at the wheel center and its local $z$-axis gives the axle direction
$\hat a$. Marker $r$ defines a planar road with normal $\hat n$. The contact
frame is

$$
\hat e_x=\frac{\hat a\times\hat n}
{\|\hat a\times\hat n\|},
\qquad
\hat e_y=\hat n\times\hat e_x.
$$

The cross product must be nonzero. Its order establishes the positive rolling
direction from marker orientation alone. For wheel center $C$, road point $P$,
and unloaded radius $R$,

$$
h=(C-P)^T\hat n,
\qquad
\delta=R-h,
\qquad
Q=C-h\hat n.
$$

$Q$ is the projection of the center onto the instantaneous road plane and is
the force application point. Because the road normal may move,

$$
\dot h=(V_C-V_P)^T\hat n+(C-P)^T\dot{\hat n},
\qquad
\dot\delta=-\dot h.
$$

Let $V_r(Q)$ be the velocity of the road-body point instantaneously
coincident with $Q$. The transport and tread-point relative velocities are

$$
v_c=V_C-V_r(Q),
$$

$$
v_q=V_C+\Omega_w\times(Q-C)-V_r(Q).
$$

The first velocity describes where the wheel is traveling. The second
describes slip of the wheel material at the contact point. They give

$$
V_x=v_c^T\hat e_x,\qquad V_y=v_c^T\hat e_y,
$$

$$
s_x=v_q^T\hat e_x,\qquad s_y=v_q^T\hat e_y.
$$

For regularization speed $V_0>0$,

$$
\kappa=-\frac{s_x}{\sqrt{V_x^2+V_0^2}},
\qquad
\alpha=\operatorname{atan2}
\left(V_y,\sqrt{V_x^2+V_0^2}\right).
$$

Thus pure rolling gives $s_x=0$. A driven wheel whose tread moves backward
relative to the road has positive $\kappa$. The camber output is

$$
\gamma=\operatorname{atan2}
\left(\hat a^T\hat n,\|\hat a\times\hat n\|\right).
$$

The one-sided normal law first calculates a candidate $F_n^*$. The applied
normal force is

$$
F_n=\begin{cases}
\max(0,F_n^*),&\delta>0,\\
0,&\delta\leq0.
\end{cases}
$$

The longitudinal and lateral expressions calculate trial forces $F_x^*$ and
$F_y^*$. They are zero when $F_n=0$. With the optional friction ellipse,

$$
u=\sqrt{
\left(\frac{F_x^*}{\mu_xF_n}\right)^2+
\left(\frac{F_y^*}{\mu_yF_n}\right)^2},
$$

and both trial forces are divided by $u$ when $u>1$. The global force on the
wheel is

$$
F^g=F_x\hat e_x+F_y\hat e_y+F_n\hat n.
$$

It acts at $Q$. A non-ground road body receives the opposite force at the same
point. The force pair is therefore coincident and introduces no artificial
net moment.

The tire owns nine kinematic variables, five scalar load variables, three
global-force variables, and the same number of local implicit equations.
Without explicit relaxation lengths, its force expressions use instantaneous
slip and there are no tread-deformation states. Supplying both relaxation
lengths adds one longitudinal and one lateral first-order deformation state.
Component-local forward-mode differentiation supplies the marker-kinematic,
contact-frame, and load-transfer partials. The three scalar laws use the
ordinary restricted-expression gradients. Tire lift-off is handled by the
ordinary BDF corrector and error controller without event location or a
forced history restart.

## Span coordinates

Spanning forces and spanning motions share an internal span coordinate. For
ordered marker points $P_1$ and $P_2$, it allocates the explicit variables

$$
s=P_2-P_1,\qquad
\ell=\|s\|,\qquad
\hat u=\frac{s}{\ell},
$$

$$
\dot\ell=\hat u^T(V_2-V_1),
$$

and

$$
\ddot\ell=\hat u^T(a_2-a_1)+
\frac{(V_2-V_1)^T(V_2-V_1)-\dot\ell^2}{\ell}.
$$

These definitions appear as nine local implicit equations for the three
components of $s$, the distance, the three components of $\hat u$, the
distance velocity, and the distance acceleration. The distance must be
positive. Keeping these intermediate quantities explicit makes them available
as outputs and keeps the element Jacobian small and local. The transverse
velocity term in $\ddot\ell$ accounts for a span whose direction changes.

A standalone `span` measurement consists only of these nine definitions. It
has no scalar-force equation, reaction, or prescribed-motion equation, so it
does not change mobility. Its distance and velocity may be used by scalar
force expressions and all nine variables are retained as outputs. Like the
directed-distance measurement, it is excluded from state selection and
initial-condition correction.

## Spanning force

A spanning force adds the definitions

$$
f=F(t,z),\qquad
F_I^g=-f\hat u
$$

to the shared span. These are four more local implicit equations: one for the
scalar force and three for its global components.

The span direction $\hat u$ points from the first marker toward the second.
The stored $f$ is the force on the first marker in the opposite line-of-sight
direction. Positive $f$ is compression and negative $f$ is tension. The force
on the second marker is $-F_I^g$. For a body marker at body-fixed offset
$r^b$, the corresponding body-frame moment is

$$
T^b=r^b\times (A^{gb})^T F^g,
$$

with the appropriate sign for that end. A ground end contributes no body
balance equation.

The built-in linear law is

$$
f=-k(\ell-\ell_0)-c\dot\ell.
$$

A constant or restricted scalar expression may replace it. Expression
dependencies are named canonical variables. Forward-mode dual numbers supply
the local partials, which are inserted into the same analytical sparse
Jacobian as the geometry and body-force contributions.

## Spatial pulley belts

A pulley pitch circle has center $C_i$, unit revolute axis $a_i$, and positive
radius $r_i$. A straight span joins tangent points

$$
P_i=C_i+r_i n_i,
$$

where $n_i$ is a unit radial vector perpendicular to $a_i$. Its unit direction
from the first pulley to the second is

$$
t=\frac{P_2-P_1}{\|P_2-P_1\|}.
$$

Tangency requires

$$
t^Tn_1=t^Tn_2=t^Ta_1=t^Ta_2=0.
$$

For parallel axes, the calculation is the ordinary planar common-tangent
construction performed in their common pitch-circle plane. It provides the
two external and, when the radii and center distance allow them, two crossed
tangents.

For nonparallel axes, the possible line direction is fixed apart from sign:

$$
t_0=\frac{a_1\times a_2}{\|a_1\times a_2\|}.
$$

The four radial candidates are obtained from

$$
n_i=\sigma_i(t_0\times a_i),\qquad \sigma_i\in\{-1,1\}.
$$

After forming $P_1$ and $P_2$, the program tests the four tangency dot
products above. This is also an existence test: arbitrary nonparallel pitch
circles need not have a common straight tangent. Initial near-point markers
select the feasible branch closest to the intended contact locations.

Define the signed rolling radius at each end by

$$
\beta_i=r_i(a_i\times n_i)^Tt.
$$

Its magnitude is the pitch radius and its sign describes how positive pulley
rotation feeds belt into the span. Let $\theta_i$ be the continuous revolute
coordinate and let $\gamma_i$ be the angle of $n_i$ around the base marker's
axis. The span extension is retained as an algebraic definition,

$$
e=\ell-\ell_0+
\beta_2[(\theta_2-\gamma_2)-(\theta_{2,0}-\gamma_{2,0})]-
\beta_1[(\theta_1-\gamma_1)-(\theta_{1,0}-\gamma_{1,0})]+e_0.
$$

This is the spatial form of the planar length, tangent-migration, and pulley
rotation relation. It does not add an integrated belt state. The instantaneous
surface velocity at contact is

$$
v_i^c=V_{C_i}+\omega_i\times(P_i-C_i),
$$

so

$$
\dot e=t^T(v_2^c-v_1^c).
$$

Each span uses the bilateral spring-damper law

$$
T=ke+c\dot e,
$$

and applies $Tt$ to the first pulley and $-Tt$ to the second at the calculated
tangent points. Ordered spans form a closed loop. The initial path length is
the sum of the straight lengths and the circular wrap lengths between the
incoming and outgoing contact points on each pulley. Initial tension is
distributed among the spans according to their compliance.

The first implementation assumes a massless belt, no slip, circular pitch
surfaces, and a geometrically feasible straight span. It deliberately does not
model slack, belt mass, axial walking, or a freely twisted belt between
incompatible pulley planes.

## Spanning motion

A spanning motion adds prescribed distance, velocity, and acceleration
equations to the shared span. The user supplies only a positive distance
expression $\ell_p(t)$. Forward-mode differentiation supplies
$\dot\ell_p(t)$ and $\ddot\ell_p(t)$, so the three prescribed equations remain
consistent by construction. A scalar reaction produces equal-and-opposite
forces at the two marker points. The reaction scalar is the force on the first
marker along the $J$-to-$I$ line of sight, so the first-marker force is
$-f\hat u$. A constant distance is an ideal massless link with spherical ends.

## Spherical joint

A spherical joint connects two marker points but does not constrain their
orientations. Let the marker positions be

$$
P_i^g=R_i^g+A^{gi}r_i^i,
\qquad
P_j^g=R_j^g+A^{gj}r_j^j.
$$

The position equation is

$$
\Phi=P_i^g-P_j^g=0.
$$

Its velocity equation is

$$
V_i^g+A^{gi}(\omega_i^i\times r_i^i)
-V_j^g-A^{gj}(\omega_j^j\times r_j^j)=0.
$$

Its acceleration equation is

$$
a_i^g+A^{gi}\left(
\alpha_i^i\times r_i^i+
\omega_i^i\times(\omega_i^i\times r_i^i)
\right)
-a_j^g-A^{gj}\left(
\alpha_j^j\times r_j^j+
\omega_j^j\times(\omega_j^j\times r_j^j)
\right)=0.
$$

A ground marker contributes its fixed position and zero velocity and
acceleration.

The joint allocates the global reaction vector $\lambda^g$. Positive
$\lambda^g$ is the force applied to the first marker, and $-\lambda^g$ is
applied to the second marker. For a body-side marker, this gives the body-frame
moment

$$
T^b=r^b\times\left((A^{gb})^T\lambda^g\right).
$$

The reaction force and moment enter the ordinary body balance equations. No
special joint balance equation is needed.

## Perpendicular-axis constraint

Let $\hat x_i$ be the first marker's oriented $x$-axis and $\hat y_j$ the
second marker's oriented $y$-axis. The position equation is

$$
\Phi=\hat x_i\mathbin{\cdot}\hat y_j=0.
$$

Both axes move with their owners. Their derivatives give the exact velocity
and acceleration equations

$$
\dot\Phi=
\dot{\hat x}_i\mathbin{\cdot}\hat y_j+
\hat x_i\mathbin{\cdot}\dot{\hat y}_j=0,
$$

$$
\ddot\Phi=
\ddot{\hat x}_i\mathbin{\cdot}\hat y_j+
2\dot{\hat x}_i\mathbin{\cdot}\dot{\hat y}_j+
\hat x_i\mathbin{\cdot}\ddot{\hat y}_j=0.
$$

For a body-fixed unit vector $u^b$,

$$
\dot{\hat u}^{g}=A^{gb}(\omega^b\times u^b),
$$

$$
\ddot{\hat u}^{g}=A^{gb}\left[
\alpha^b\times u^b+
\omega^b\times(\omega^b\times u^b)\right].
$$

The scalar reaction $\lambda$ applies the global torque

$$
T_i^g=\lambda\hat n,
\qquad
T_j^g=-\lambda\hat n,
\qquad
\hat n=\hat x_i\times\hat y_j.
$$

This normal is configuration dependent. The analytical Jacobian therefore
includes its derivatives with respect to both marker orientations.

## Inplane constraint

The Inplane uses the common directed distance with the second marker's
oriented $z$-axis. Its position equation is

$$
\Phi=(P_i-P_j)^T\hat z_j=0.
$$

With $d=P_i-P_j$, the differentiated equations are

$$
\dot\Phi=(V_i-V_j)^T\hat z_j+d^T\dot{\hat z}_j=0,
$$

$$
\ddot\Phi=(a_i-a_j)^T\hat z_j+
2(V_i-V_j)^T\dot{\hat z}_j+d^T\ddot{\hat z}_j=0.
$$

The direction derivatives use the same body-fixed-vector formulas given for
the perpendicular-axis constraint. The analytical Jacobian includes both the
marker-point motion and the changing plane normal.

The scalar reaction applies $\lambda\hat z_j$ at $P_i$. Its opposite acts on
the second marker's owner at a floating point coincident with $P_i$. Thus the
two forces have no net force or moment on the complete system even when
$P_i$ is not coincident with $P_j$. If the second marker belongs to ground,
only the force on the first body enters the stored model equations.

## Inline constraint

The inline primitive is the intersection of two inplane constraints. With
$d=P_i-P_j$, its position equations are

$$
d^T\hat x_j=0,
\qquad
d^T\hat y_j=0.
$$

The two equations use the inplane velocity and acceleration derivatives with
normal axes $\hat x_j$ and $\hat y_j$. Their scalar reactions produce the
transverse force

$$
F_i=\lambda_x\hat x_j+\lambda_y\hat y_j,
$$

with the opposite force applied to the second body at a floating point
coincident with $P_i$. There is no axial reaction and no constraint on
relative orientation.

When translation coordinates are requested, the signed distance is

$$
s=d^T\hat z_j.
$$

Its exact derivatives are

$$
v=(V_i-V_j)^T\hat z_j+d^T\dot{\hat z}_j,
$$

$$
a=(a_i-a_j)^T\hat z_j+
2(V_i-V_j)^T\dot{\hat z}_j+d^T\ddot{\hat z}_j.
$$

The coordinate state equations are

$$
a-\dot v=0,
\qquad
v-\dot s=0.
$$

The velocity $v$ joins the body and hinge-coordinate velocities considered by
pivoted QR. All coordinate definitions have analytical Jacobian
contributions.

## Coordinate coupler

Let $q_k$, $v_k$, and $a_k$ be the position, velocity, and acceleration of a
hinge rotation or inline translation coordinate. One coupler adds

$$
\Phi=\sum_k c_kq_k-q_0=0,
$$

$$
\dot\Phi=\sum_k c_kv_k=0,
\qquad
\ddot\Phi=\sum_k c_ka_k=0.
$$

These are component-local implicit equations with direct Jacobian entries
$c_k$. A single reaction variable $\lambda$ supplies the generalized reaction
$c_k\lambda$ at coordinate $k$.

For a hinge coordinate, $c_k\lambda$ becomes equal and opposite torques about
the second marker's moving $z$-axis. For an inline coordinate, it becomes
equal and opposite forces along that axis, applied at the first marker point
and its coincident floating point on the second body. These are the same
reaction mappings used by the corresponding rotational and translational
motion generators. Their derivatives include the changing marker axis and,
for translation, the changing moment arms.

An initial offset records $\sum_k c_kq_k$ after the relative coordinates have
been initialized from the supplied model configuration and any relative
initial conditions. A numeric offset instead imposes an absolute coordinate
relation. The coupler equation participates in the velocity-constraint partial
matrix used for QR state selection, so any of its referenced velocities can
be retained when the partition permits it.

## Ideal gear pair

For side $i$, let $C_i$ be the base-side marker of the referenced revolute,
$\hat a_i$ its $z$-axis, and $P$ the contact point fixed to the carrier. Define

$$
r_i=P-C_i,
\qquad
t_i=\hat a_i\times r_i.
$$

Each $t_i$ must be nonzero. The normalized directions must also be collinear,
allowing either sign:

$$
\left|\hat t_1^T\hat t_2\right|\simeq 1.
$$

The implementation rejects an angular mismatch greater than $10^{-5}$
radians and reports it in degrees. This check works for both parallel-axis
gears and common-apex bevel gears. In particular, the two radial vectors of a
common-apex bevel pair can be identical, so $r_1\times r_2$ is not a useful
contact-direction test.

Choose $\hat t=\hat t_1$. The signed effective pitch radius is the scalar
triple product

$$
\rho_i=\hat a_i^T(r_i\times\hat t).
$$

Its magnitude is the moment arm of the contact force about the joint axis.
Its sign includes whether the second pitch tangent is parallel or
antiparallel to the chosen contact direction. The gear equation is

$$
\Phi=\rho_1\theta_1-\rho_2\theta_2-\phi_0=0.
$$

The angles $\theta_i$ measure the first marker of each revolute relative to
the carrier frame about the corresponding carrier-stored axis. Their
definitions are periodic implicit equations, so their reported values can
remain continuous through complete rotations. For `phase = "initial"`,

$$
\phi_0=\rho_1\theta_{1,0}-\rho_2\theta_{2,0}.
$$

The velocity and acceleration constraints are evaluated directly from the
carrier-relative body angular velocities and accelerations:

$$
\dot\Phi=\rho_1\omega_1-\rho_2\omega_2=0,
\qquad
\ddot\Phi=\rho_1\alpha_1-\rho_2\alpha_2=0.
$$

Writing these two equations directly in terms of body motion lets the QR
state-selection partial matrix see the physical coupling. The element still
retains $\theta_i$, $\omega_i$, and $\alpha_i$ as named output variables.

One scalar reaction $\lambda$ applies $+\lambda\hat t$ to the first gear and
$-\lambda\hat t$ to the second at generated floating contact markers. The
axial moments are therefore $\rho_1\lambda$ and $-\rho_2\lambda$, exactly the
generalized reactions of the scalar constraint. This automatic formulation
uses the pitch tangent as a zero-pressure-angle force direction. A later
extension can supply a different marker-defined contact direction for helical
or spiral-bevel geometry.

The mechanical element does not create graphics. A body may separately own a
`gear` graphic positioned and oriented by a referenced body marker. This keeps
the pitch geometry used by the equations independent of its visual
representation.

The carrier need not be the base body of both revolute joints. In a planetary
set, a contact marker fixed to the moving planet carrier supplies the
carrier-relative radial direction and angular reference even when the sun or
ring bearing is fixed to ground. An external sun-to-planet pair and internal
planet-to-ring pair can therefore share the same planet revolute while using
two different carrier contact markers.

## Rack and pinion

The spatial rack and pinion is limited to an ordinary spur arrangement. Let
the carrier-side marker of the inline constraint define the orthonormal pitch
frame $(\hat x_c,\hat y_c,\hat z_c)$. Its axes have the following meanings:

- $\hat x_c$ is parallel or antiparallel to the pinion revolute axis
  $\hat a$;
- $\hat y_c$ points from the pinion center toward the pitch contact; and
- $\hat z_c$ is the rack translation and contact-force direction.

For pinion center $C$ and positive pitch radius $r_p$, the carrier-fixed
contact point is

$$
P=C+r_p\hat y_c.
$$

The signed rolling radius is the scalar triple-product reduction

$$
\rho=r_p\hat a^T(\hat y_c\times\hat z_c)
     =r_p\hat a^T\hat x_c.
$$

Thus $\rho$ is positive when the pinion axis follows $\hat x_c$ and negative
when it is antiparallel. Let $q$ be the inline distance along $\hat z_c$ and
$\theta$ the continuous revolute angle about $\hat a$. The ideal no-slip
constraint is

$$
\Phi=q-\rho\theta-\phi_0=0,
$$

with corresponding velocity and acceleration equations

$$
\dot\Phi=\dot q-\rho\omega=0,
\qquad
\ddot\Phi=\ddot q-\rho\alpha=0.
$$

For `phase = "initial"`, $\phi_0=q_0-\rho\theta_0$, so adding the
element does not alter the supplied initial registration. A numeric phase has
units of length.

One scalar reaction $\lambda$ applies $+\lambda\hat z_c$ to the rack and
$-\lambda\hat z_c$ to the pinion at generated floating markers coincident
with $P$. The pinion moment about its joint is $-\rho\lambda$, which is the
generalized reaction obtained from the constraint. The complete force and
moment partials, including a moving carrier frame and contact point, are
inserted in the sparse Jacobian by local forward differentiation.

The loader requires the inline and revolute base markers to share one carrier
and rejects a pinion axis whose angular mismatch from $\pm\hat x_c$ exceeds
the specified tolerance. Skew and helical rack geometry are intentionally
outside this element.

## Hinge constraint

The hinge is a pair of perpendicular-axis primitives. Its ordered markers
supply the axes in

$$
\Phi_1=\hat x_i\mathbin{\cdot}\hat z_j=0,
\qquad
\Phi_2=\hat y_i\mathbin{\cdot}\hat z_j=0.
$$

Because $\hat x_i$ and $\hat y_i$ span the plane normal to $\hat z_i$, these
equations make $\hat z_i$ and $\hat z_j$ collinear. The initial configuration
selects the parallel or antiparallel branch. Relative rotation around the
common axis is not constrained.

Each scalar equation uses the velocity and acceleration derivatives given for
the perpendicular-axis primitive. With two scalar reactions, the torque on
the first marker is

$$
T_i^g=\lambda_1(\hat x_i\times\hat z_j)+
      \lambda_2(\hat y_i\times\hat z_j),
$$

and the opposite torque acts on the second marker. Thus the hinge transmits
the two torque components perpendicular to its axis but no axial torque. It
does not constrain translation. A spherical joint and hinge acting at the
same ordered markers form a spatial revolute joint.

### Optional hinge rotation coordinate

Let

$$
c=\hat x_j\mathbin{\cdot}\hat x_i,
\qquad
s=\hat z_j\mathbin{\cdot}(\hat x_j\times\hat x_i).
$$

When requested, the hinge allocates the relative variables $\theta$, $\omega$,
and $\alpha$. The initial angle is $\operatorname{atan2}(s,c)$. During motion,
the periodic implicit equation

$$
\operatorname{atan2}
\left(
\sin\theta\,c-\cos\theta\,s,
\cos\theta\,c+\sin\theta\,s
\right)=0
$$

relates the orientation to an unwrapped $\theta$. Because the equation uses
only the periodic functions of $\theta$, the integrated coordinate may pass
continuously through any number of complete revolutions without encountering
an artificial discontinuity at $\pm\pi$.

With owner angular quantities resolved globally, the rate definitions are

$$
\omega=(\omega_i^g-\omega_j^g)\mathbin{\cdot}\hat z_j,
$$

$$
\alpha=(\alpha_i^g-\alpha_j^g)\mathbin{\cdot}\hat z_j+
(\omega_i^g-\omega_j^g)\mathbin{\cdot}\dot{\hat z}_j,
$$

where

$$
\dot{\hat z}_j=\omega_j^g\times\hat z_j.
$$

The two optional state equations are

$$
\alpha-\dot\omega=0,
\qquad
\omega-\dot\theta=0.
$$

Pivoted QR considers $\omega$ alongside the body velocity components. If it
is selected, these two state equations replace the equivalent body state
equations. All coordinate definitions and their analytical Jacobian
contributions remain local to the hinge.

## Revolute joint

The spatial revolute joint is a composite of a spherical joint and a hinge
between the same ordered markers. It contributes five scalar constraint
families: three make the marker points coincident and two remove relative
rotation perpendicular to the second marker's $z$-axis. Its reaction variables
are the spherical joint's three global force components and the hinge's two
scalar torque components.

The composite retains the primitive equation implementations and analytical
Jacobian contributions. It only combines their registrations under one
component name. If `rotation_coordinates = true`, the hinge part also supplies
the relative `theta`, `omega`, and `alpha` definitions and the two candidate
state equations described above.

## State equations and system size

The scalar velocity constraints remove independent velocity components.
Pivoted QR of their partial matrix chooses the remaining physical velocities.
For a dependent translation, its acceleration-to-velocity and
velocity-to-position state equations are inactive. For a dependent angular
velocity, its acceleration-to-velocity and velocity-to-pseudo-angle state
equations are inactive. All three pseudo angles remain active, together with
the three equations that carry their rates into the Euler parameters. The
position constraints determine dependent pseudo-angle corrections. The active
system therefore remains square and keeps only the selected state equations.

Finite geometry is evaluated from Euler parameters. In the dynamic Newton
matrix, mechanical orientation partials are placed in the pseudo-angle columns
instead of the Euler-parameter columns. The orientation bridge converts those
corrections to Euler-parameter corrections. Before a corrector begins, the BDF
predictor adjusts the nonphysical pseudo-angle values so this bridge is
consistent. The pseudo-angle values themselves are not physical outputs.

For a single body connected to ground, the usual selection retains the three
body-fixed angular velocities. Translation is then determined by the spherical
constraint.

Adding a hinge at the spherical joint removes two angular-velocity states and
leaves only the angular velocity around the hinge axis. For the aligned
single-body example, QR selects `omega_z`.

State selection is not frozen after initialization. DDASSL can request a new
selection after a singular iteration matrix, repeated corrector failures, or
an unhealthy physical predictor error. The spatial runner evaluates the
velocity-constraint partial matrix at the current configuration and repeats
the pivoted QR selection. If the selected set changes, it replaces only the
candidate state equations and the differential and error-control masks. It
then rebuilds the sparse Jacobian pattern and symbolic factorization. The
canonical variables and accepted BDF history are retained.

## Initial conditions

Initial consistency is established before state selection:

1. Body origins and orientations are corrected together so the position
   equations hold. Each orientation correction is a three-component rotation
   in the current body frame rather than an additive change to the four Euler
   parameters.
2. Translational and angular velocities are projected onto the velocity
   constraints.
3. Body accelerations and joint reactions are solved together from the body
   balances and acceleration constraints.

The spatial implementation uses pivoted QR of the transpose of the scaled
velocity-constraint partial matrix to choose independent equation rows. A
dependent row identifies one scalar ideal-constraint family. Its reaction
variable and its position-, velocity-, and acceleration-level equations are
all deactivated. Initial position and velocity projection may temporarily use
a rank-revealing independent subset so a consistent redundant model can be
assembled before this removal. This is separate from the column-pivoted QR
used to choose state variables.

The row basis chosen by pivoted QR is not unique. Dependencies can couple
translation and orientation equations, so QR need not deactivate the scalar
family that appears most obviously redundant from the mechanism geometry. In
the parallel-axis spatial four-bar, for example, an analyst might remove one
$z$-closure equation and two axis-alignment equations at the closing joint.
The numerical factorization can instead retain all four $z$-closure equations
and remove three perpendicular-axis families distributed around the loop.
Both choices have rank 17 and leave the same single mechanism freedom.

This equivalence applies to the local constraint manifold represented by the
retained independent rows; it does not make the omitted constraint into a
compliance. Nevertheless, the individual multipliers of an ideal redundant
model are indeterminate before a row basis is chosen. The calculated reactions
therefore belong to the selected independent system and must not be interpreted
as a unique physical load distribution among redundant ideal joints. Replacing
the relevant ideal connections with finite-stiffness elements makes that load
distribution determinate when it is important.

### Possible preferred constraint removals

Experience with particular classes of mechanisms may show that some valid row
bases work better than others. The program may eventually allow the analyst to
identify constraint families that are preferred for removal, in the same way
that states can now be preferred. Automatic QR would remain the default.

A preferred-removal option should be a request rather than an unchecked
deactivation. The loader would verify that each name refers to a complete
scalar constraint family, that the requested removals do not lower the
required rank, and that the remaining equations are adequately conditioned.
The analyst could specify only part of the redundant set and allow QR to choose
the remainder. An `allow_fallback` option could either permit a fully automatic
choice or make an unusable preference an input error.

This option should not be added until experience shows what makes one valid
choice better than another. Possible reasons include better conditioning over
finite motion, more useful reaction reporting, and removal choices that are
easier for an analyst to understand. Preferred removal would not make the
reactions of an ideal redundant model unique.

## Static analysis

Static analysis uses the same component equations and variables as dynamics,
but selects only force and moment balance, position-level constraints and
motion prescriptions, and the static definitions needed by the active force
elements. All velocity and acceleration variables are set to zero. This
removes inertia and velocity-dependent damping without requiring separate
static versions of the components.

A spatial body has three rotational degrees of freedom even though its
orientation is stored in four Euler parameters. The static Newton system
therefore uses a three-component incremental rotation $\delta\phi$ in the
current body frame. At the linearization point, this maps into the canonical
Euler-parameter columns as

$$
\delta p=\frac{1}{2}Q(p)\delta\phi.
$$

The Newton correction is applied multiplicatively,

$$
A_{\mathrm{new}}=A\,\exp([\delta\phi]_\times),
$$

and converted back to normalized Euler parameters. Thus the static system has
three orientation unknowns for each body, remains square, and every Newton
iterate stays on the rotation manifold. Translational positions, relative
coordinates, reactions, and applied-force variables receive ordinary additive
corrections. A backtracking line search rejects corrections that do not reduce
the largest implicit-equation error.

With more than one requested output time, each converged equilibrium predicts
the next one by linear extrapolation. If Newton iteration cannot cross an
interval, the interval is divided and the intermediate equilibrium becomes a
new continuation point. These intermediate points help convergence but are not
included among the requested output samples.

Dynamic relaxation is an optional way to approach a difficult equilibrium.
The ordinary spatial dynamic equations advance through pseudo-time while every
component sees a fixed physical model time. Velocity and acceleration are
reduced between relaxation cycles, and an exact static Newton correction is
tried after every cycle. Only that final static solution is retained.

For static initialization of dynamics, the static position and orientation are
retained, the velocities declared in the model are restored, and the velocity
constraints are projected. The dynamic initializer then solves consistent
accelerations, reactions, and applied-force variables before BDF integration
begins.

## Jacobian

Each component supplies analytical derivatives of all three constraint levels
and its reaction contributions. Quaternion derivatives are taken through the
physical rotation matrix. The complete active BDF Jacobian

$$
\frac{\partial G}{\partial z}
+c\frac{\partial G}{\partial\dot z}
$$

is checked against central finite differences in the examples.
