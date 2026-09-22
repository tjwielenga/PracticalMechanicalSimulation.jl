# Planar Element Formulations

This chapter of the Technical Manual describes the equations assembled by the
current planar modeling elements. The [TOML User's
Guide](../../docs/planar/toml-reference.md) describes how to enter these elements. This
chapter explains what the program does with them.

## Common notation

A body retains center-of-mass position $R$, orientation $\theta$, translational
velocity $V$, angular velocity $\omega$, translational acceleration $A$, and
angular acceleration $\alpha$. For a body-fixed marker with local position
$r^b$,

$$
P=R+C(\theta)r^b,
$$

where $C$ is the planar rotation matrix. Its velocity and acceleration are

$$
\dot P=V+\omega\,\hat z\times C r^b,
$$

$$
\ddot P=A+\alpha\,\hat z\times C r^b-\omega^2 C r^b.
$$

Ground markers have fixed global position and orientation. A floating marker
uses the point motion of another marker while retaining force and torque
ownership on its own body.

Each ideal connection contributes position-, velocity-, and acceleration-level
implicit equations. Its reaction variables enter the body balances. Compliant
elements contribute forces and torques without reducing kinematic mobility.

### Directed distance

Several planar elements use the same marker geometry. The second marker's
oriented local $y$-axis is the direction

$$
\hat n=\hat y_j,
$$

and the signed distance from its plane to the first marker is

$$
d=(P_i-P_j)^T\hat n.
$$

Its first two derivatives are

$$
\dot d=(\dot P_i-\dot P_j)^T\hat n+(P_i-P_j)^T\dot{\hat n},
$$

$$
\ddot d=(\ddot P_i-\ddot P_j)^T\hat n
       +2(\dot P_i-\dot P_j)^T\dot{\hat n}
       +(P_i-P_j)^T\ddot{\hat n}.
$$

The inplane constraint sets $d=0$. A translational motion generator prescribes
$d(t)$. A distance coordinate reports $d$, $\dot d$, and $\ddot d$. Plane
contact uses $d-r$ as the gap for a sphere of radius $r$. The applied-force
element needs only the common directed axis $\hat n$. These elements share one
internal calculation of the axis and directed-distance kinematics.

### Marker-to-marker spans

For two marker points, define

$$
s=P_2-P_1,\qquad \ell=\|s\|,\qquad \hat u=\frac{s}{\ell}.
$$

The distance derivatives are

$$
\dot\ell=\hat u^T(\dot P_2-\dot P_1),
$$

$$
\ddot\ell=\hat u^T(\ddot P_2-\ddot P_1)
 +\frac{\|\dot P_2-\dot P_1\|^2-\dot\ell^2}{\ell}.
$$

The reaction-free `span` element reports these quantities. The
`spanning_force` element uses the same geometry and velocity equations before
applying its scalar constitutive force along $\hat u$. A zero initial length is
rejected because the unit direction and the derivative equations would be
undefined.

## Bodies

A planar rigid body contributes three balance equations:

$$
mA_x-\sum F_x=0,\qquad
mA_y-\sum F_y=0,\qquad
I\alpha-\sum T=0.
$$

Mass and inertia may be zero. The corresponding balance equation is then
algebraic. Connections or force laws must still determine that direction, or
the assembled system is singular. A massless body's velocity may remain in the
selected state history; this does not add inertia because its algebraic balance
continues to determine the motion.

The characteristic length used in automatic state selection is normally the
largest distance from the body center to one of its markers. If no useful
marker radius exists, the loader tries the radius of gyration $\sqrt{I/m}$ and
then a model-wide fallback length.

## Ideal connections

### Revolute joint

For marker points $P_i$ and $P_j$, a revolute joint imposes

$$
\Phi=P_i-P_j=0.
$$

Its velocity and acceleration equations are obtained by differentiating this
equation. The two reaction components act at the coincident marker points with
equal and opposite signs.

When rotation coordinates are requested, the joint also defines the continuous
relative quantities

$$
\theta_r=\theta_i-\theta_j,\qquad
\omega_r=\omega_i-\omega_j,\qquad
\alpha_r=\alpha_i-\alpha_j.
$$

These definitions add explicit variables and equations. They do not add a
constraint. The relative velocity becomes a state-selection candidate.

### Distance coordinate

For the second marker's oriented local $y$-axis $\hat y_j$, the measured
distance is

$$
q=(P_i-P_j)^T\hat y_j.
$$

Differentiating this definition supplies its velocity and acceleration
variables. The coordinate itself contributes no reaction. It provides output,
a possible state coordinate, and an interface for a coordinate coupler.

### Coordinate coupler

A coupler imposes one scalar relation among coordinate ports:

$$
\Phi=\sum_k c_k q_k-b=0.
$$

The same coefficients multiply the coordinate velocities and accelerations.
The generalized reaction on coordinate $k$ is $c_k\lambda$. A rotation port
turns this into equal-and-opposite marker torques. A distance port turns it
into equal-and-opposite forces along the measurement axis.

For `offset = "initial"`, the loader evaluates $\sum c_kq_k$ from the initial
configuration. A numeric offset instead establishes an absolute phase.

### Inplane constraint

The inplane constraint is

$$
\Phi=(P_i-P_j)^T\hat y_j=0.
$$

It constrains the first point to the plane through the second marker whose
normal is $\hat y_j$. Differentiating $\Phi$ supplies the velocity and
acceleration equations, including the motion of the normal when marker $j$
rotates. One scalar reaction acts along $\hat y_j$.

### Perp constraint

The perp constraint makes the first marker's local $x$-axis perpendicular to
the second marker's local $y$-axis:

$$
\Phi=(\hat x_i)^T\hat y_j=0.
$$

On the constraint surface, its differentiated equations reduce to

$$
\omega_i-\omega_j=0,\qquad
\alpha_i-\alpha_j=0.
$$

Its scalar reaction is an equal-and-opposite marker torque.

### Translational and fixed joints

A translational joint combines one inplane and one perp primitive. Translation
remains free along the second marker's local $x$-axis. The two primitive
reactions remain separate.

A fixed joint combines a revolute and a perp primitive. Its two force reactions
come from the revolute and its torque reaction comes from the perp constraint.

### Ideal gear pair

An ideal gear pair references two revolute joints and a carrier contact point.
The gear center for each side is the base-side marker of its revolute joint.
At least one joint base must own the contact marker. That joint supplies the
reference colinear unit vector

$$
d_b=P-C_b, \qquad
\hat d_0 = \frac{d_b}{\|d_b\|}
$$

and tangent to the pitch circles

$$
\hat t=\hat z_b \times \hat d.
$$

Let $d_1$ and $d_2$ run from the two base-side joint markers to the contact
point. The signed pitch radii are the dot products of these vectors with the colinear unit vector .

$$
\rho_1=d_1 ・\hat d_0,\qquad
\rho_2=d_2 ・ \hat d_0.
$$

Their magnitudes are the physical pitch radii. They have the same sign for an
internal pair and opposite signs for an external pair. The constraint is

$$
\rho_1(\theta_1-\theta_c) -
\rho_2(\theta_2-\theta_c)-\phi_0=0.
$$

The same signed radii multiply the relative angular velocities and
accelerations.

For `phase = "initial"`,

$$
\phi_0=
\rho_1(\theta_{1,0}-\theta_{c,0})-
\rho_2(\theta_{2,0}-\theta_{c,0}),
$$

so the gear constraint does not alter the supplied initial angles. A numeric
phase establishes absolute tooth registration.

The scalar reaction applies $\lambda\hat t$ to the first gear and
$-\lambda\hat t$ to the second at generated floating contact markers. These
forces produce contact forces on the gears which result in gear moments
around their joints.

### Rack and pinion

The carrier-side translational marker supplies axes $\hat x_c,\hat y_c$. If
$C$ is the carrier-side revolute marker and $r_p$ is the pitch radius, the
contact point is

$$
P_c=C+r_p\hat y_c.
$$

The ideal constraint is

$$
q+r_p(\theta_p-\theta_c)-\phi_0=0,
$$

where $q$ is rack displacement along $\hat x_c$. The element creates one
carrier contact marker and floating markers on the rack and pinion.
Equal-and-opposite forces along $\hat x_c$ generate the rack force and pinion
moment.

### Pulley, belt, and belt span

Each span selects a common tangent between two pitch circles. The supplied near
points choose one of the four possible tangent branches during assembly. That branch is
then followed as the mechanism moves.

For a span, the program retains tangent points $P_1,P_2$, tangent direction
$t$, straight length $\ell$, extension $e$, extension rate $\dot e$, tension
$T$, and global force. Pulley surface velocity contributes to $\dot e$, so
fixed-center pulley rotation transfers belt material between adjacent spans.
The assembled pulley angles define the extension reference and prevent
arbitrary initial angular phase from creating tension.

The constitutive law is

$$
T=ke+c\dot e.
$$

When damping is inferred from a local time scale, $c=k\tau_d$. Positive tension
pulls both pulleys toward the opposite tangent point. The forces are applied at
the tangent points and therefore produce the pulley moments. The current belt
is massless, bilateral, and no-slip; negative tension is retained rather than
being replaced by a slack-belt model.

## Force elements

Applied force, applied torque, spanning force, bushing, plane-contact, and
curve-contact components carry an analysis-stage activation set. Their constitutive
equations set all explicit load variables to zero while inactive; geometry and
rate definitions remain present. The body-balance contributions continue to
read those explicit variables, so switching stages changes the load without
changing the allocated sparse system. Static-to-dynamic handoff reinitializes
the load variables immediately after changing the active stage. SimpView does
not draw the inactive element's load or connector symbol.

### Gravity and applied loads

Gravity contributes $m g$ to each selected body at its center of mass.

A marker-directed applied force uses the second marker's local $y$-axis as its
direction and the first marker as its application point. When a reaction body
is supplied, a generated floating marker applies the equal-and-opposite force
to that body.

Applied torque acts on the first marker orientation and applies its
equal-and-opposite to the second marker orientation. Constant and
time-expression versions use the same mechanics.

An applied force or applied torque may depend on named scalar configuration
and velocity variables. For a state-dependent load law $f(x,t)$, the element
owns a load variable $L$ and the local implicit equation

$$
L-f(x,t)=0.
$$

The body balances depend on $L$ and the marker geometry. The constitutive row
depends only on $L$ and the qualified kinematic variables found while parsing
the expression. This preserves component locality and gives the sparse matrix
its exact structural columns before the analysis starts.

The expression is evaluated with multidirectional dual numbers. If its named
inputs are $x_1,\ldots,x_n$, their independent dual components return all
$\partial f/\partial x_i$ in one local forward-mode evaluation. The
constitutive Jacobian entries are

$$
\frac{\partial G_L}{\partial L}=1,\qquad
\frac{\partial G_L}{\partial x_i}=-\frac{\partial f}{\partial x_i}.
$$

Acceleration, reaction, and applied-load variables are excluded. These would
change the character of the force law or permit implicit load cycles and can
be considered separately if a physical need appears.

### Torsional spring-damper

For principal relative angle $\theta_{12}$,

$$
T=-k(\theta_{12}-\theta_0)-c(\omega_1-\omega_2).
$$

The stored torque $T$ is the load on the first marker. The second marker
receives $-T$. Positive stiffness and damping therefore oppose positive
relative displacement and velocity. If the damping coefficient is omitted,
$c=k\tau_d$.

### Spanning force

Let $s=P_2-P_1$, $\ell=\|s\|$, and $u=s/\ell$. The vector $u$ points from the
first marker toward the second. The stored scalar $f$ is the force on the first
marker in the opposite, $J$-to-$I$, direction. The global force on the first
marker is

$$
F_1=-fu,
$$

and $F_2=-F_1$. A positive $f$ is compression and repels the markers; a
negative $f$ is tension and attracts them. The scalar law may be a constant,
the predefined spring-damper

$$
f=-k(\ell-\ell_0)-c\dot\ell,
$$

or an expression using `element.length`, `element.length_rate`, and other
named configuration and velocity variables. The linear partials are supplied
analytically, while expression partials are calculated with dual numbers.

The program retains the nine explicit local variables `s_x`, `s_y`, `length`,
`u_x`, `u_y`, `length_rate`, `force`, `F_x`, and `F_y`. Their definition
equations are part of the unreduced system rather than calculations hidden
inside one force callback.

### Planar bushing

The bushing measures translation in the second marker's frame and relative
rotation between the marker frames. Diagonal translational stiffness and
damping produce the two local force components; rotational stiffness and
damping produce the torque. The resulting loads are transformed to global
coordinates and applied with equal-and-opposite signs.

For local damping time scale $\tau$,

$$
C_t=\tau K_t,\qquad c_r=\tau k_r.
$$

Explicit damping coefficients independently override these estimates.

### Plane contact

The first marker is the center of a sphere of radius $r$. The second marker
defines a plane with outward normal $\hat y_2$. The signed gap is

$$
g=(P_1-P_2)^T\hat y_2-r.
$$

With penetration $\delta=\max(-g,0)$ and closing speed $v_c=-\dot g$, the
normal force is

$$
F_n=k\delta\max(0,1+d v_c).
$$

An optional transition depth replaces $\delta$ near first contact by an
integrated smootherstep stiffness. The effective penetration and its first
five derivatives join smoothly at contact entry and at the end of the
transition. A user expression is an alternative to the built-in law. It is
evaluated from the explicit `gap` and `gap_rate` variables and is applied
without an implicit separation test or force clamp.

The law is compliant rather than an inequality constraint. It is continuous at
first contact, increases damping during closing motion, reduces force during
rebound, and cannot pull the sphere toward the plane.

The integrator locates contact entry and exit and ends the current step at the
crossing. It retains recent BDF history, limits the order to two, refreshes the
numerical iteration matrix, and reduces the next step. When damping is active,
it also locates a zero of $1+d v_c$ because that changes the force-law branch.
Expression contacts do not request inferred root events because their branch
structure is not known to the element.

### Closed curve and circular roller contact

A `curve` is a periodic cubic spline $q(s)$ in a marker frame. Chord length of
the supplied control polygon defines its station coordinates. A cyclic spline
solve establishes second derivatives so $q$, $q'$, and $q''$ are continuous
where the profile closes.

For roller center $C$ and global profile point $Q(s)$, the contact owns the
station $s$ and enforces

$$
(C-Q)^T\hat t=0,
$$

where $\hat t=Q'/\lVert Q'\rVert$. The differentiated equation solves the
explicit `station_rate`. The control-polygon orientation selects a consistent
outward normal $\hat n$; `side = "inside"` reverses it. For roller radius $r$,

$$
g=(C-Q)^T\hat n-r.
$$

Both $\dot g$ and the tangency-rate equation include body motion, curve-frame
rotation, and motion of $Q$ along the curve. The signed curvature is

$$
\kappa=\frac{d\hat t}{d\ell}\mathbin{\cdot}\hat n.
$$

The normal law is the same built-in or expression law used by plane contact.
The roller receives $F_n \hat n$ and the curve body receives the opposite force
at $Q$. The local variables retain the station, station rate, curvature, gap,
gap rate, contact point, normal, scalar force, and global force in the
unreduced equation system. This makes the geometry available for plotting and
for later force extensions without repeating the contact search.

### Closed curve and flat follower contact

The flat follower reuses the periodic curve and the same explicit station,
rate, curvature, gap, and force variables. Let $P_f$ be its marker point and
let $\hat x_f$ and $\hat y_f$ be the marker axes. The local $x$-axis lies in
the plate and the local $y$-axis is the positive force direction. Tangency is
the scalar equation

$$
\hat t(s)^T\hat y_f=0.
$$

Its time derivative contains the curve-frame angular velocity, station rate,
and follower angular velocity. It therefore determines the explicit station
rate without a search inside the force law. The signed gap is

$$
g=(P_f-Q(s))^T\hat y_f.
$$

Positive gap is separation in the positive follower-normal direction;
negative gap is compliant penetration. The follower receives
$F_n\hat y_f$ at $Q(s)$ and the curve body receives its opposite. Applying the
force at $Q$ rather than at $P_f$ retains the moment from a contact point that
is offset along the face.

Unless the user supplies `initial_station`, initialization samples the profile
for the largest value of $Q^T\hat y_f$ and Newton-corrects that station to the
tangency equation. The selected station then remains unwrapped and continuous
through the periodic boundary. The built-in normal law, expression law,
staging, and BDF transition handling are identical to circular roller contact.

## Friction

The planar surface, revolute, translational, and inplane friction elements
share one scalar bristle law. For slip $v$, static and dynamic coefficients
$\mu_s\geq\mu_d$, and transition speed $v_s$,

$$
\mu(v)=\mu_d+(\mu_s-\mu_d)e^{-v^2/v_s^2},\qquad G=\mu(v)N.
$$

Here $N$ is the effective normal load and $G$ is a force capacity. For a
revolute, $G$ is additionally multiplied by the effective bearing radius and
is therefore a torque capacity. With shear state $s$, stiffness $k$, damping
$c$, and nonzero capacity,

$$
\dot s=v-\frac{k|v|}{G}s,\qquad f=\operatorname{clamp}(-ks-c\dot s,-G,G).
$$

When $G$ is negligible, $f=0$ and $\dot s=-s/t_r$, so stored shear is released
without division by zero. In statics the differential equation is replaced by
$s=q-q_0$, where $q$ is the associated relative angle or tangential distance.
The anchor $q_0$ is established after initial-position correction and includes
any shear transferred from a saved result.

For `surface_friction`, let $P_i$ be the sphere center, $P_j$ the plane
marker, $\hat n$ the plane normal, and $\hat t$ its local $x$-axis. The plane
projection of the sphere center is

$$
Q=P_i-[(P_i-P_j)^T\hat n]\hat n.
$$

With sphere radius $r$, the signed slip includes motion of both contacting
surfaces:

$$
v=\hat t^T\left[V_i+\omega_i\mathbin{\times}(-r\hat n)
-V_j-\omega_j\mathbin{\times}(Q-P_j)\right].
$$

The planar cross products denote the in-plane velocity produced by the scalar
angular velocity. The contact's one-sided normal force gives
$G=\mu(v)\max(F_n,0)$. The tangential force and its opposite act at $Q$; when
contact opens, the force is zero and stored shear decays with $t_r$.

For a planar revolute with point-reaction vector $R$ and effective radius
$r_b$,

$$
N=\|R\|+N_0,\qquad G=\mu(\omega)r_bN.
$$

The slip is the first marker's angular velocity relative to the second. The
friction torque acts on the first marker and its opposite on the second.

For a translational joint or standalone inplane constraint, let $\lambda$ be
the inplane primitive's signed normal reaction. Then

$$
N=|\lambda|+N_0.
$$

The second marker's local $x$-axis is the tangential direction. Its rotation is
included when differentiating the tangential coordinate. The resulting global
force and its opposite act at the two marker points. The translational joint's
perpendicular reaction is a torque and is deliberately not included in this
normal-load estimate. Because the underlying inplane constraint is bilateral,
its friction law is also bilateral; it is not a substitute for one-sided
contact detection.

## Motion generators

### Translational motion

The prescribed coordinate is

$$
d=(P_i-P_j)^T\hat y_j.
$$

The position equation sets this coordinate equal to the requested distance
history. Its differentiated equations prescribe velocity and acceleration.
The reaction is an equal-and-opposite marker force along $\hat y_j$. Because
only one scalar distance is prescribed, translation along the plane and
relative rotation remain free.

### Rotational motion

The generator prescribes the first marker's angle relative to the second
marker, together with consistent angular velocity and acceleration histories.
Its scalar reaction is an equal-and-opposite marker torque.

Constant-speed generators construct all three histories from an initial value
and speed. Expression generators accept the three histories separately; the
analyst is responsible for making them mutually consistent.

## Assembly and Jacobians

Each element registers the variables and equations it owns before allocation.
After allocation it supplies callbacks for its implicit equations, balance
contributions, and Jacobian entries. Component contributions are accumulated
into the unreduced sparse system. See [Modular sparse
assembly](../design-record/modular-sparse-assembly.md) for the assembly contract and [Adding a
planar component](adding-planar-component.md) for a complete example.
