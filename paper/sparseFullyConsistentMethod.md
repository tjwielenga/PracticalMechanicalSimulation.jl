# Sparse Fully Consistent Modeling Method

Thomas Wielenga

Wielenga Innovation Foundation, Inc.

4218 N Surf Rd, Hollywood, FL  33019

thomas@wielenga.org

ASME Member:  20001319227

## ABSTRACT

This paper presents the Sparse Fully Consistent Modeling Method, in which all
variables satisfy all levels of the constraint equations (displacement,
velocity and acceleration) at all times. The equations are formulated in
implicit form and integrated by a stiff (BDF) integrator. The equations are
produced component by component, are relatively simple, and allow the
formulation of the partial derivatives needed for a Jacobian in simple terms.
Assembly of these partials into a Jacobian produces a large but very sparse
Jacobian which can be solved efficiently. Comparisons are made to other
existing methods and model formulations. Planar examples and benchmarks are
examined, and the results show that the large unreduced system can be solved
efficiently. Both a working 2D mechanism simulation program and a 3D mechanism
simulation program are freely available for close inspection and reuse. They
include an extensive set of modeling elements and methods for calculating
consistent initial conditions, static equilibrium, and modal analysis. In
addition, a graphical interface allows animation and plotting of system
variables, and a Lua modeling interface allows the building of complex
hierarchical models.


## 1. INTRODUCTION

Real mechanical systems often contain stiff bushings, dampers, contact forces,
flexible parts, and driveline components. These components introduce high
frequencies that can decay rapidly. A "non-stiff" integrator must take steps
small enough for the highest frequency or fastest mode, even after its motion
has decayed. A stiff integrator can take much larger steps (orders of magnitude
larger), but it requires a Newton-like correction and an efficient solution of
a matrix of partial derivatives (a Jacobian). For efficiency reasons, a 
practical general mechanism
simulation program requires the use of a stiff integrator.

Many integrators require the user to explicitly solve for the derivative of the state
variables. Other integrators instead formulate the
equations in implicit form. These implicit integrators can even solve equations that include
algebraic relationships in them. These kinds of equations are called
differential-algebraic equations (DAEs) and they must be solved by a stiff integrator.

Some modeling methods reduce the system to a smaller set of coordinates and
then solve for the highest derivatives in order to use an explicit integrator.
This produces fewer equations, but solving them requires a factorization of a
complicated mass matrix. If the system is stiff, a Jacobian is required and
this process produces an inordinately complicated Jacobian. Jacobian-free
methods can avoid forming the entire Jacobian, but they may require a
preconditioning matrix that is still difficult to create.

Using a stiff implicit DAE solver is the method of choice for practical
simulation programs and has been used for many years in the popular ADAMS
simulation program. But even here, there are many choices as to which equations
and variables to use in the formulation. We will explore some of these choices
and compare their advantages and disadvantages. We will present the Sparse
Fully Consistent Modeling Method. It keeps the component equations in 
implicit form, and each component produces a small number of variables,
equations, and Jacobian partials. These contributions are assembled into a
large sparse system. The sparse matrix solver finds an efficient order for
solving for the variables, eliminating any need for special solution paths for
open chains, tree topologies, or closed loops.

Fully Consistent means the position, velocity, and acceleration constraint
equations are present and satisfied during the entire solution process, not
just at output times. The method keeps all the physical variables including
forces, accelerations, velocities and positions in the solution set, and adds
differential equations only for the states needed to represent the degrees of
freedom. This produces a solvable set of equations at each solution step.


## 2. EQUATION FORMULATIONS

To illustrate the different variable and equations choices, we will use the planar pendulum.  The pendulum has a mass and rotational inertia and it is acted on by gravity. It is attached to ground by a revolute joint that allows only rotation about its pivot.

**Figure 1: Pendulum under gravity**


There are variables for the position of the mass, its velocity and its acceleration.  There is an angle, an angular velocity, and an angular acceleration.  Finally, there is an angle that can be defined in the revolute joint and its derivatives.  Because there is a constraint, there are also constraint forces that enforce the constraint.  


**Table 1. Candidate variables for the planar pendulum**

| Variable | Symbol |Count |
|---|---|---:|
| Center-of-mass position| $R^g=[R_x\;\;R_y]^T$ | 2 |
|  Translational velocity| $V^g=[V_x\;\;V_y]^T$ | 2 |
|  Translational acceleration| $a^g=[a_x\;\;a_y]^T$ | 2 |
|  Body angle, velocity, acceleration | $\theta_b$, $\omega_b$, $\alpha_b$  | 3 |
|  Pin angle velocity, acceleration| $\theta_p$, $\omega_p$, $\alpha_p$|3| 
|  Pin reaction| $\lambda^g=[\lambda_x\;\;\lambda_y]^T$ | 2 | 
| **Total variables** | | **14** |


**Table 2.  Candidate equations**

| Kind | Equation block | Count |
|---|---|---:|
|  Force balance| $m a^g-\lambda^g-mg^g=0$ | 2 |
|  Moment balance| $J\alpha-(d^g)^T\lambda^g=0$ | 1 |
|  Acceleration pin constraint| $a^g+d^g\alpha-r^g\omega^2=0$ | 2 |
|  Velocity pin constraint| $V^g+d^g\omega=0$ | 2 |
|  Position pin constraint| $R^g+r^g-p_0^g=0$ | 2 |
|   **Total equations**  || **9** |

In each method the mechanism remains the same but the equations used change. The following
formulations will be compared:

- reduced coordinates;
- displacement constraints with Deficit Two;
- velocity constraints with Deficit One;
- acceleration constraints with Deficit Zero;
- Baumgarte stabilization;
- GearStableV;
- GearStableA; and
- the Fully Consistent method.

An ideal position constraint is written

$$
\Phi(q,t)=0.
$$

The velocity and acceleration constraint equations are obtained by
differentiating the position constraint:

$$
\dot\Phi=D(q,t)v+\phi(q,t)=0,
$$

$$
\ddot\Phi=D(q,t)a+\gamma(q,v,t)=0.
$$

Constraint Deficit is the number of additional differentiations needed before
the accelerations and reactions can be solved. A formulation using only the
position constraint has Deficit Two. The velocity formulation has Deficit One,
and the acceleration formulation has Deficit Zero. This definition is based on
the mechanical quantities of interest. It does not increase the deficit to
account for derivatives of algebraic variables that are not needed by the
mechanical solution.

The Index of a set of differential-algebraic equations is related to
Constraint Deficit, but it depends on the mathematical definition of the DAE
and on which variables are considered solvable. Constraint Deficit gives the
analyst a direct description of how much constraint differentiation is needed
to obtain accelerations and reactions.

GearStableV retains the position and velocity constraint equations and adds
one set of satisfaction multipliers. GearStableA also retains the acceleration
constraint equations and uses two sets of satisfaction multipliers. These
methods keep the indicated constraint levels satisfied, but allow small
differences between some physical variables and their time derivatives.

The Fully Consistent method retains all three levels of the constraint
equations. It adds differential equations for a minimal set of physical states
instead of adding satisfaction multipliers.


## 3. Fully Consistent equations

The variables used by the program can be grouped as

$$
x=(q,v,a,\lambda,\xi).
$$

The variables $q$, $v$, and $a$ are the positions, velocities, and
accelerations. The variable $\lambda$ contains the constraint reactions. The
variable $\xi$ contains other component variables used to define geometry,
rates, forces, and constitutive equations.

The complete set of implicit equations can be written as

$$
F=
\begin{bmatrix}
F_B\\
F_{Kq}\\
F_{Kv}\\
F_C\\
F_F\\
F_S
\end{bmatrix}=0.
$$

The equations $F_B$ are the force and moment balances. The equations $F_{Kq}$
and $F_{Kv}$ relate positions, velocities, accelerations, and their time
derivatives. The equations $F_C$ are the active constraint equations at the
position, velocity, and acceleration levels. The equations $F_F$ define forces
and other constitutive relationships. The equations $F_S$ are the differential
equations for the selected physical states.

Acceleration is an explicit unknown. Constraint reactions are applied to the
balance equations using the transpose of the virtual-power relationships.
Geometric and force quantities can also remain as explicit unknowns. They do
not have to be eliminated before the equations are assembled.

Redundant constraint equations make the equation set singular. The program
finds redundant constraint rows and deactivates the entire equation family
associated with each redundant row. All independent levels of the active
constraints remain in the equation set.

**Table 2, pendulum equations and variables:** The final paper should include
an equation-and-variable table for the Fully Consistent pendulum. The rows
should show the force balance, moment balance, three levels of the pin
constraint, and the two selected state equations. The columns should show the
eleven variables $a^g$, $\alpha$, $V^g$, $\omega$, $R^g$, $\theta$, and
$\lambda^g$. This table will show which variables occur in each equation and
will give a specific example of the general notation above. It can be adapted
from the full equation inventory in
[`planar-rigid-body-pendulum.md`](../examples/planar/planar-rigid-body-pendulum.md)
using the eleven-variable selection in
[`planar-pendulum-independent-state-implicit-system.md`](../examples/planar/planar-pendulum-independent-state-implicit-system.md).

**Figure 2:** The sparsity pattern of the pendulum or four-bar Jacobian. The
rows are labeled by equation type and the columns by variable type.

## 4. Preparing the equations for solution

The supplied model is an initial estimate. Its positions must satisfy the
position constraint equations before states are selected. The program first
corrects the body positions and orientations. The correction is weighted so
that a massive body normally moves less than a light body. A translational
weight is based on mass. A rotational weight uses $mL^2$, where $L$ is a
characteristic body length. The velocities are corrected in the same way when
consistent velocities are needed.

The velocity constraint matrix is then evaluated at the consistent
configuration. It shows how the physical velocities are related:

$$
D=\frac{\partial\dot\Phi}{\partial v}.
$$

The translational velocity columns are scaled using a characteristic length
for each body. Angular velocity columns have a scale of one. The nonzero rows
are normalized after column scaling. This keeps the rank test and state choice
from being determined only by the units used for the equations.

A QR decomposition of the transpose is used when the constraint matrix does
not have full row rank. It identifies redundant constraint rows. When one row
is redundant, the complete position, velocity, and acceleration equation
family belonging to that row is deactivated. The remaining constraint levels
stay in the equation set.

A separate column-pivoted QR decomposition is used to select states:

$$
DP=QR.
$$

The pivot columns identify velocity components that are dependent on the
constraints. The remaining velocity components are independent and are
candidates for states. The matching positions or angular coordinates are
selected with them. The number of selected state pairs is equal to the number
of degrees of freedom. The other physical variables remain in the equation set
and are solved along with the states.

The user can request preferred states. A preferred state is used only if the
constraint matrix retains the required rank. The chosen active equation
families and selected state equations close the square system. The remaining
accelerations, reactions, and force variables are then initialized together.

State selection and redundant-equation detection are different calculations
even though both use QR decomposition.

**Table 3:** Compare automatic state selection, body angular velocities, and
relative joint angular velocities for the ten-second pendulum simulation.

**Figure 3:** The steps used to prepare the equations for analysis.

```text
model -> allocate the system -> correct positions and velocities
      -> find active constraints -> choose states -> form the square system
      -> initialize remaining variables -> analyze
```

## 5. Scaling

The equations contain translations, rotations, forces, moments, and
dimensionless constraint equations. The model therefore needs physical scales
as well as derivative-level scales. Characteristic length, mass, and velocity
are quantities that an analyst can estimate for a machine or vehicle. They
provide the related time, force, and moment scales. Body marker locations and
inertia properties can provide local characteristic lengths when the user does
not supply them.

The position, velocity, and acceleration variables and equations also have
different derivative levels. For a step of size $h$, let $e_i$ be the level of
equation $i$ and $v_j$ the level of variable $j$. A Jacobian entry is scaled as

$$
\widehat J_{ij}=|h|^{e_i-v_j}J_{ij}.
$$

A coefficient produced by a BDF derivative contains the corresponding inverse
power of $h$. The level scale cancels that power. The leading coefficients in
the scaled Jacobian therefore remain approximately independent of step size.
This improves conditioning and makes it practical to reuse a Jacobian
factorization as the integrator changes its step size.

The paper should show this first for

$$
\dot R-V=0
$$

and then for

$$
m\dot V-F=0.
$$

This simple example makes the purpose of equation and variable levels clear.

## 6. Integration

At each correction, the BDF integrator solves a linearized set of the implicit
equations using

$$
J=F_x+c_jF_{\dot x}.
$$

The variable-step, variable-order BDF method integrates the selected physical
states. All variables are stored in the history polynomial so positions,
velocities, accelerations, reactions, and force variables can be predicted at
the next step and corrected together. Physical positions and velocities can be
included in the integration-error test even when they were not selected as
states. Reactions and other algebraic variables are normally excluded from
error control, but they remain in the Newton correction and in the stored
history.

The sparsity pattern remains fixed while the selected states and active
equations remain fixed. The symbolic factorization can therefore be reused.
Numerical factors can also be reused for several Newton corrections when the
scaled Jacobian has not changed greatly. If correction fails with the old
factors, a new Jacobian and numerical factorization are tried before the step
is rejected. Repeated correction failures can cause the symbolic factorization
to be repeated. A change in state selection or active equations also requires
the affected factorization work to be repeated.

When a continuous force law changes stiffness abruptly, as at contact, a soft
restart retains useful history while reducing the prediction order. A state
choice that becomes unhealthy does not require an immediate repartition at
every step. Singular factorization, repeated correction failures, or an error
that does not respond normally to step reduction can request a new state
selection. The complete variable history is retained, but the closing state
equations and sparse factorization are replaced.

## 7. Other analyses

The initialized equation assembly supports several kinds of analysis:

- A kinematic mechanism is solved as a complete implicit system at each time.
- A dynamic analysis integrates the selected physical states.
- Static equilibrium can be found directly by Newton iteration or approached
  by dynamic relaxation before a final Newton correction.
- Consecutive static solutions can follow a slowly changing input without
  inertial forces.
- A static solution can be used as the initial configuration for a dynamic or
  modal analysis.
- Values from a saved result can be transferred to a compatible model when
  component names and types agree.

### Modal analysis

The complete implicit equations can be linearized at an operating point. The
linear modal equations have the form

$$
(J+sE)\hat x=0.
$$

The calculation uses the shifted matrix

$$
U=(J+\sigma E)^{-1}E.
$$

Only columns of $E$ associated with differential states are nonzero. A single
sparse factorization reduces the remaining calculation to a dense eigenvalue
problem based on those states. The large algebraic part of the model is
handled by sparse solves and does not add spurious finite modes. The original
implicit equations are used to check each calculated mode.

**Table 4:** Give the natural frequencies and original-equation errors for a
model with several degrees of freedom. The one-degree pendulum is useful as a
unit test. A rotor train or multi-link spatial pendulum is a better paper
example.

The important point is that these analyses reuse the component equations. They
do not require separate kinematic, dynamic, static, and modal versions of every
body, joint, or force.

## 8. Planar and spatial mechanics

The planar model provides the simplest way to explain the equations. A planar
body has two translational components and one angular component. The spatial
model uses three translations, three angular velocities, and three angular
accelerations. Euler parameters represent the finite orientation, while
body-fixed pseudo angles provide the infinitesimal rotation columns used in
the Jacobian and in state selection. The kinematic differential equations
update the Euler parameters from angular velocity.

Markers give positions and orientations on bodies. Joints and force elements
are defined between markers rather than being written directly for each pair
of bodies. Virtual-power relationships apply scalar and vector reactions to
the force and moment balances. The same relationships provide compact local
Jacobian contributions.

The program now implements the method in both planar and spatial models. The
spatial examples exercise consistent orientation correction, redundant
constraints, state selection and reselection, static equilibrium, dynamic
integration, modal analysis, stiff bushings, contact, tires, gears, belts, and
saved-result initialization. The paper should use only enough of this material
to show that the method is not limited to planar equations. It should not
become a catalog of model elements.

## 9. Examples and results

The program contains enough components to test the method on useful
mechanisms. Revolute, in-plane, perpendicular, fixed, translational, spherical,
hinge, inline, and orient constraints provide ideal constraint equations.
Bushings, torsional spring-dampers, spanning forces, tires, and contact provide
stiff force elements. Gear pairs, rack-and-pinion elements, couplers, and
pulley belts add internal variables and moving geometry.

Organize the examples by the question they answer:

- The pendulum compares equation formulations and accuracy.
- Serial pendulum chains show sparse scaling in open systems.
- Four-bars and parallelogram chains exercise closed loops, redundant
  constraints, and state choices.
- Rotor trains provide a smooth stiff problem with an analytical modal
  reference.
- Bouncing balls exercise a continuous force whose stiffness changes at
  contact.
- Static and modal examples show reuse of the complete equation assembly.
- A concise spatial mechanism shows that the formulation carries into three
  dimensions.
- The large van shows the breadth of a practical spatial model containing
  suspensions, steering, bushings, tires, static equilibrium, dynamics, and
  modal analysis. It is evidence of breadth rather than the principal
  performance benchmark.

These elements and examples demonstrate component equations, ideal
constraints, stiff forces, internal variables, and changes in force stiffness.
The paper does not need to describe every model input field. The complete
models can be supplied so that the results can be run and viewed.

Each timing result must identify the Git commit, Julia version, computer,
operating system, model, tolerances, warm-up procedure, number of runs, and
reported statistic. Model loading, initialization, and integration time should
be reported separately.

The results should include the number of model variables, selected states, and
nonzero Jacobian entries. They should also include memory allocation, accepted
and rejected steps, Newton corrections, symbolic factorizations, numerical
factorizations, and a measure of solution error. These values help explain the
timing results.

**Table 5, open chains:** Use the 10-, 25-, and 50-link pendulums. Show the
growth in model size and run time. Show separately the improvement from using
UMFPACK, reusing the symbolic factorization, and reusing numerical factors.

**Table 6, ideal closed loops:** Use the parallelogram chains. Show the run time
and the growth of the dense QR work used for initialization and state
selection.

**Table 7, a smooth stiff system:** Use the rotor train. Include solution
errors, modal errors, and the range of natural frequencies.

**Table 8, changing force stiffness:** Use the bouncing-ball bank to compare
hard and soft integrator restarts. Report rejected steps separately from Newton
correction failures.

**Figures 4--6:** Plot run time and memory allocation against the number of
model variables. Plot the number of symbolic and numerical factorizations
against the attempted steps. Plot the hard and soft restart results with their
rejected-step counts.

## 10. What the results show

The results gathered so far show that a large unreduced system does not imply
a dense solution. The canonical size and sparse nonzero count grow nearly
linearly for the serial chains. Replacing generic factorization with UMFPACK
changed the 50-link pendulum from an impractical calculation into a short one.
Reusing symbolic analysis and then reusing numerical factors reduced the work
further.

The dense QR calculations used for initialization, redundant-row detection,
and state selection have not been the principal cost for the tested models.
Their growth is more noticeable in the closed-loop chains and remains a limit
to watch. The state comparisons also show that a physically convenient state
is not automatically the fastest state, and a state that is valid initially
may become unhealthy later.

The rotor train shows that the method can span a widening range of smooth
frequencies while retaining the requested accuracy. The bouncing-ball bank
shows that compliant contact is more difficult, but also shows that a soft
restart can materially reduce rejected work. Planar and spatial models use the
same component-local assembly and sparse solution strategy.

The final text in this section should follow the regenerated publication
tables. It should not claim more accuracy or scaling than those tables show.

The Fully Consistent method uses the same component description for open
systems, closed loops, and loops connected by stiff forces. Constraint
reactions and other physical results remain directly available. The selected
states are recognizable physical variables. The same model supports
kinematic, dynamic, static, and modal analysis.

## 11. Limitations and future work

The method solves a larger equation set than a reduced-coordinate method. Its
performance depends on sparse ordering and factorization. Dense QR will become
more expensive for very large systems with many closed loops. Sparse
rank-revealing methods may eventually be needed.

Runtime state reselection is presently a recovery method rather than a
continuous optimization of the state partition. Better health measures and
better rules for choosing a replacement state remain useful work. Contact and
friction still cause rejected steps and can make the complete model difficult
to solve. Massless bodies are allowed, but their automatic state selection and
general behavior need more testing.

Automatic differentiation can help form component Jacobian partials.
Jacobian-free Newton--Krylov methods can avoid assembling the complete
Jacobian, but they still require a good preconditioner. Flexible bodies and
more extensive spatial validation are also future work. These additions can
use the existing component assembly without changing the central method.

## 12. Conclusions

The examples and benchmarks should determine the final conclusion. The work
done so far shows that a mechanism model can retain its physical variables and
all useful levels of its constraint equations without using a dense matrix. A minimal
set of physical states can be integrated while the complete sparse system is
solved. Component equations, equation and variable scaling, and sparse
factorization make the large system practical for the planar and spatial
mechanisms tested so far.

The conclusion should return directly to the main question. A mechanical
system can keep its positions, velocities, accelerations, reactions,
constraint equations, and force definitions consistent, integrate only the
states representing its degrees of freedom, and still be solved efficiently.

## Reproducing the results

The paper should give the commands used to run the program tests, paper tests,
and benchmarks. The models and the data used to produce the tables should be
stored with the paper. The Git commit, Julia version, computer, solver settings,
and tolerances should be recorded. Examples using the supported program should
be distinguished from older examples that compare equation formulations.

## Work still needed

1. Freeze or configure the historical integrator used by the formulation
   comparison tests and restore a passing paper-test baseline.
2. Store the benchmark results as data and provide one command to regenerate
   the publication tables.
3. Choose and verify the model used for the modal-analysis table.
4. Produce the equation table and figures identified in the outline.
5. Review published work before making claims about what is new.
6. Reconcile the abstract with the final evidence without changing its writing
   style.
7. Prepare the final manuscript using the ASME journal format.
