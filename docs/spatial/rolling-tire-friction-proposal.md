# Rolling-tire friction proposal

This records the design discussion preceding the first experimental bristle
tire option. The current implementation and its limitations are described in
the [spatial technical manual](../../architecture/spatial/spatial-element-formulations.md#optional-load-dependent-tire-bristles)
and [TOML reference](toml-reference.md#rolling-tire). Some alternatives below
remain proposals, not claims about the implemented element.

## Why use the new friction law?

The present tire calculates wheel-road slip, a one-sided normal force, and
longitudinal and lateral forces from user expressions. It can limit the two
forces together with a friction ellipse. By default it has no tangential
deformation states. Optional relaxation lengths add two states, but a direct
static solve cannot develop their values from a stationary initial position
alone. See [the current tire formulation](../../architecture/spatial/spatial-element-formulations.md#rolling-tire).

Our new [surface-friction law](../../architecture/spatial/spatial-element-formulations.md#tangential-surface-friction)
has two shear states. They let a stationary contact resist a force below its
friction limit, and they release when contact force disappears. This is useful
for a tire, particularly when settling a vehicle on the road.

We should **not** add a `surface_friction` element alongside a `rolling_tire`.
Both would apply tangential force at the same contact. Nor should we simply
replace the tire's force expressions with the existing bristle equations. A
rolling tire continuously replaces the tread in its contact patch; a stationary
block does not. With no slip, the existing bristle law retains its shear even
while the tire rolls. Tire shear must instead relax as tread passes through the
patch.

## Proposed first model

Keep the tire's present wheel marker, road marker, contact point, slip
kinematics, one-sided normal force, and force application. Add an optional
*tire bristle* tangential law in place of the longitudinal and lateral force
expressions. It would own two deformation states, $u_x$ and $u_y$, measured in
meters, in the tire's forward and lateral directions. The two forces would
still be reported as `longitudinal_force` and `lateral_force` and would still
act at the existing contact point.

Let $s_x$ and $s_y$ be the existing longitudinal and lateral **tread slip
velocities**. They include wheel rotation. Let $V_x$ be the existing forward
**transport velocity** of the wheel center relative to the road. Let the
positive tangential stiffnesses $K_x(F_n),K_y(F_n)$ and relaxation lengths
$L_x(F_n),L_y(F_n)$ take their values at the current normal load. Suppressing
that argument below, a useful starting law is

$$
\dot u_x=s_x-\frac{|V_x|}{L_x}u_x-pu_x,
\qquad
\dot u_y=s_y-\frac{|V_x|}{L_y}u_y-pu_y,
$$

The $|V_x|/L$ terms replace tread as the tire rolls. The $p$ term releases
excess elastic shear when the contact slides. It is the tire counterpart of
the sliding term in our new friction law. The sign convention makes a positive
slip velocity produce a negative force on the wheel. At $V_x=0$ the rolling
replacement term is zero, so a parked tire can retain shear.

For the first version, use the elastic trial force

$$
F_x^*=-K_xu_x,
\qquad
F_y^*=-K_yu_y.
$$

At small, steady slip and nonzero rolling speed, the sliding term is of higher
order in slip and these equations give approximately

$$
F_x\simeq-K_xL_x\frac{s_x}{|V_x|},
\qquad
F_y\simeq-K_yL_y\frac{s_y}{|V_x|}.
$$

For small physical slip angle $\alpha\simeq s_y/|V_x|$, the steady cornering
stiffness is therefore $C_\alpha=K_yL_y$, independent of speed. The lateral
force takes time $L_y/|V_x|$ to relax, or travel distance $L_y$. The same
calibration applies longitudinally with small slip ratio. This gives the
analyst a way to choose $K_y=C_\alpha/L_y$ from tire data. This rolling
formula is not used at exactly zero speed: there, stored $u_y$ supplies the
static holding force.

This relationship is also used in lumped tire models and experimental
relaxation-length measurements; see [this test-rig study](https://link.springer.com/article/10.1007/s11012-023-01684-z).

The tire's current reported `slip_angle` uses a regularization speed. That is
useful to avoid a division at rest, but reusing the current regularized
transport speed in the new shear equation would make the predicted cornering
stiffness vary at very low speed. The bristle law must use the actual $|V_x|$
and tread slip velocities, with explicit handling at zero speed. If a narrow
numerical smoothing band proves necessary, it must lie below the speeds used
to measure cornering stiffness, and its effect on the measured slope must be
checked. Low-speed calibration should use the geometric slip angle, not the
regularized output near its cutoff.

An added term $-C_ys_y$ would contribute $C_y|V_x|$ to the *steady*
cornering stiffness. We should therefore omit this damping term from the
first version rather than give the model an unintended speed dependence.
Transient damping, if needed, should be designed and tested separately.

## Cornering stiffness and vertical load

Cornering stiffness must depend on normal load. It approaches zero as
$F_n\to0$, but over the heavily loaded range it generally grows less than
proportionally to $F_n$. This load sensitivity is important to vehicle
handling and belongs in the first tire-bristle model, not a later refinement.
If measured values are available, interpolate a user-supplied
$C_\alpha(F_n)$ curve. A simple concave law could be a fallback, but its
shape should be checked against tire measurements. The longitudinal
small-slip stiffness needs its own load curve; it should not silently inherit
the lateral one.

At each positive load, once an effective relaxation length $L_y(F_n)$ is
chosen, set

$$
K_y(F_n)=\frac{C_\alpha(F_n)}{L_y(F_n)}.
$$

This makes the desired force-versus-slip slope the calibration target.
If the relaxation length follows the shrinking contact patch, $K_y$ must
change with load as well; changing patch length alone would impose an
arbitrary cornering-stiffness curve. At zero load, use the no-contact branch
rather than evaluate $0/0$.

For two tires carrying $F_0-\Delta F$ and $F_0+\Delta F$, a concave
$C_\alpha(F_n)$ gives less combined cornering stiffness than two tires each
carrying $F_0$. The heavily loaded tire does not make up for what the
lightly loaded tire loses. This is why the load curve affects the Large Van's
turning and rollover response even when total axle load is unchanged.
The friction limits $\mu_xF_n$ and $\mu_yF_n$ are separate calibrations;
cornering stiffness alone does not determine peak force.

## Friction limit and sliding

Use the tire's calculated normal force $F_n\geq0$, not a joint reaction. For
positive longitudinal and lateral coefficients $\mu_x,\mu_y$, define

$$
G_x=\mu_xF_n,\qquad G_y=\mu_yF_n.
$$

The trial forces above would be projected onto the current friction ellipse
when necessary:

$$
H=\sqrt{(F_x^*/G_x)^2+(F_y^*/G_y)^2},
\qquad
(F_x,F_y)=\frac{(F_x^*,F_y^*)}{\max(1,H)}.
$$

This retains the existing combined-braking-and-cornering behavior. The
internal sliding rate $p$ must also prevent unbounded deformation when a
stationary wheel slides. For a loaded contact, one candidate is

$$
p=\sqrt{(K_xs_x/G_x)^2+(K_ys_y/G_y)^2}.
$$

With no rolling transport, this has the same steady sliding limit as the
new two-direction bristle law. With rolling transport, it also allows small
slip forces below the limit. This is a **candidate**, not a settled numerical
formula: its denominators require careful low-load treatment, as discussed
below. The force ellipse must be enforced even during a transient in which
the deformation state has not yet adjusted.

For the first version I suggest using the existing $\mu_x$ and $\mu_y$ as
constant limiting coefficients. The surface-friction law's separate static
and dynamic coefficients could be added later, but making that transition a
function of slip *speed* alone would move the peak to a different slip ratio
when vehicle speed changes. If a distinct peak and limiting friction are
needed, we should define that variation against tire slip and load, then
compare it with measured force curves. This simple model is not the full
contact-patch brush model or an aligning-moment model.

## Lift-off and very small normal loads

At $F_n=0$, the tire must apply **exactly zero** tangential force. Its stored
deformations should decay toward zero over a specified release time so an
old force does not reappear at the next contact. All divisions by $G_x$ or
$G_y$ must be bypassed in this unloaded branch.

The zero-load branch alone is not enough. Immediately before lift-off,
$G_x,G_y$ are small but positive, so the candidate $p$ above could become
very large and make integration difficult. The implementation should use a
small-load transition or regularization that keeps state rates and Jacobian
coefficients bounded while the **applied** force still tends exactly to zero
with $F_n$. We should choose and test this transition in the single-wheel
prototype rather than copy the present surface-friction threshold unchanged.
The tangential states should remain under error control.

The effective relaxation length need not equal the geometric patch length:
carcass compliance can contribute to the measured value. A future
load-dependent length may decrease as the patch shrinks, but it must not be
inserted as a literal zero into $|V_x|/L$. At lift-off, the zero-force and
shear-release behavior takes over. A distributed brush model would handle
patch size and tread replacement more directly, at greater complexity.

Shrinking contact must also remove stored shear. In the two-state equations
above, $p=0$ when slip is zero. A falling $F_n$ then reduces the force limit
but leaves $u$ unchanged, although some previously loaded tread has left the
road. Merely clipping the force would allow old shear to reappear if the patch
grows again. The first implementation needs an unloading rule, not just a
zero-load rule.

One possible lumped interpretation is to carry the *integrated elastic shear
force* of the active patch. If an effective loaded patch length $a(F_n)$
decreases and shear is roughly uniform, that stored force would decrease in
proportion to the surviving patch length. Growing the patch would not restore
the lost force: newly contacting tread starts unstrained. This is only an
approximation, since pressure and shear are not uniform across a real patch.
It is not yet a selected update formula. We must choose a consistent state
definition and unloading rule before implementation; we should not both scale
the state and scale a load-dependent stiffness for the same lost tread. A
small distributed brush may be preferable if the lumped rule cannot pass the
unload/reload tests.

## Static equilibrium and saved states

At the supplied initial configuration, zero shear means an unstressed tire
contact. A saved result may instead provide nonzero shear. During a
relaxation-based static solve, the changing slip can build shear and balance
a force without prescribing an artificial horizontal tire force. Preserve the
resulting shear when dynamics starts; do not reset it at the handoff.

A direct Newton static solve needs more thought. At zero velocity, a
differential shear equation says only $\dot u=0$; it does not determine $u$.
Also, the shear caused by rolling depends on the path taken to the current
pose, not just on the final wheel-center position. We should not claim that
adding two states by itself fixes direct static equilibrium. My recommended
first implementation would support relaxation-based static initialization
and its transfer to dynamics. For Newton polishing, use an incremental
no-slip relation anchored at the last accepted relaxation pose, updating the
anchor when a step slides. If that cannot be made reliable in the first
implementation, disable polish for this tire law rather than silently freeze
or arbitrarily solve its shear states. Direct Newton static from an arbitrary
pose can be a later extension.

The tire's forward and lateral axes change as the wheel steers. Stored shear
must be transported into the new contact frame; otherwise steering could
rotate a pre-existing force without slip. In local coordinates this requires
the frame-rotation term $-E^T\dot E\,u$, where $E$ contains the forward and
lateral unit vectors. The same issue arises if the road plane moves. This is
part of the first implementation, not an optional refinement.

## Possible TOML interface

This is illustrative input, **not valid in the current reader**:

```toml
[tire]
type = "rolling_tire"
markers = ["wheel.center", "ground.road"]
radius = 0.50
normal_stiffness = 15000.0
normal_damping_time_scale = 0.015

tangential_model = "bristle"
# Each row is [normal load (N), measured value]. Illustrative data only.
patch_length_by_load = [[0.0, 0.0], [3000.0, 0.22], [6000.0, 0.29]]
cornering_stiffness_by_load = [[0.0, 0.0], [3000.0, 45000.0], [6000.0, 70000.0]]
longitudinal_slip_stiffness_by_load = [[0.0, 0.0], [3000.0, 60000.0], [6000.0, 95000.0]]
lateral_relaxation_fraction = 1.0       # L_y / patch length; calibrate
longitudinal_relaxation_fraction = 1.0  # L_x / patch length; calibrate
mu_longitudinal = 0.9
mu_lateral = 0.8
shear_release_time = 0.01         # s, when unloaded
```

`tangential_model = "expression"` would mean the current behavior and remain
the default. In `"bristle"` mode the user would not also supply
`longitudinal_expression` or `lateral_expression`; mixing the two would be an
input error. The names and defaults above are proposals to review before
changing the reader. All numerical load curves are illustrative, not fitted
to the Large Van tires. The effective relaxation lengths would be the patch
length times their calibrated fractions. The $K_x(F_n)$ and $K_y(F_n)$ used
in the shear equations would be derived from the specified small-slip
stiffness curves and those lengths. An eventual reader must define
interpolation and out-of-range behavior explicitly.

## Tests before trying the van

1. A parked vehicle holds on a side slope by stored lateral shear. With the
   wheels braked, it also holds on a forward slope within the friction limit.
   These cases must work after relaxation-based static initialization.
2. Measure the small-slip lateral slope at several rolling speeds, including
   the lowest intended test-rig speed. It should remain $K_yL_y$. Check the
   longitudinal slope similarly. Under combined slip, both forces stay
   inside the ellipse.
3. Repeat the small-slip tests at several normal loads. The lateral slope
   must match $C_\alpha(F_n)$, approach zero at lift-off, and reproduce the
   reduction in combined axle stiffness under left-right load transfer.
4. At pure rolling after an earlier slip, stored shear decays as tread passes
   through the patch. Steering does not rotate that shear spuriously.
5. On lift-off, both forces reach zero and the shear state releases without
   a singular Jacobian or excessive step rejection. Re-contact starts without
   an old force spike.
6. With zero slip, lower the normal load without fully lifting the tire, then
   raise it again. Stored shear must decrease as contact shrinks and must not
   reappear when new tread enters the growing patch.
7. A relaxation-based static solve transfers its tire forces and shear states
   to a subsequent dynamic run. Saved-result initialization gives the same
   handoff.
8. Compare this model with the current expression tire on the single-wheel
   rig at several speeds and normal loads. Only then try it on the Large Van.

This gives us a small, inspectable first tire law. More detailed contact-patch
pressure, load-dependent peak slip, camber thrust, and aligning moment can be
considered separately; they need not be mixed into this first change.
