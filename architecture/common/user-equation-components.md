# User-defined equation components

Sim2D and Sim3D share one implementation for scalar algebraic variables and
first-order differential states declared by an `equation_component`. The
component uses the same canonical allocation and owned equation blocks as the
built-in mechanical elements; it is not a separate solver.

For states $s_i$, algebraic variables $a_j$, and named inputs from the
mechanical system, the entered equations have the form

$$
\dot s_i=f_i(t,s,a,z), \qquad
L_j(t,s,a,z)=R_j(t,s,a,z).
$$

An input alias is resolved to a canonical variable index before integration.
Algebraic equations become simultaneous implicit equations $L_j-R_j=0$ and
are not executed in source order. They may therefore form algebraic loops with
one another and with the mechanical equations. The number of algebraic
equations must equal the number of declared algebraic variables; the normal
Jacobian solve determines whether those equations are independent.

Each state contributes the implicit BDF equation

$$
F_i=\dot s_i-f_i(t,s,a,z)=0.
$$

If $c$ is the current BDF derivative coefficient, its Jacobian row is

$$
\frac{\partial F_i}{\partial z_k}
=c\delta_{ik}-\frac{\partial f_i}{\partial z_k}.
$$

The restricted expression compiler discovers the local dependencies and uses
dual-number automatic differentiation for their partial derivatives. The
component writes only those entries into the sparse system Jacobian. Load
expressions can then refer to the component's named variables, and the normal
load elements contribute their forces and torques to body balance equations.

Declared differential states are marked differential, included in BDF error
control, and stored by qualified name in `.simp`. They add physical states to
the model but are not candidates for QR selection of mechanical coordinates.
Algebraic variables and their equations participate in initial conditions,
statics, dynamics, and modal linearization. Acceleration initialization solves
the algebraic equations with the mechanical balances, after which the state
rates are evaluated from their equations.

During statics, `static = "steady"` adds a state to the static unknowns and
solves $f_i=0$. `static = "hold"` retains its starting value and temporarily
uses $\dot s_i=0$ during dynamic relaxation. The normal dynamic equation is
restored for dynamics and modal analysis. Saved-result initialization
transfers matching user states, while algebraic variables are recalculated.

The current implementation deliberately supports continuous scalar equations
and one explicit first-order rate law per state. It does not yet implement
arbitrary expressions containing `der(...)`, effort/flow connectors, discrete
events, or automatic unit checking. These can be added later without changing
the component-local allocation and sparse assembly.
