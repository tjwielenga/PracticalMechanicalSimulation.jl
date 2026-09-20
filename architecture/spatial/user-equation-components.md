# User-defined equation components

Sim3D can now allocate scalar algebraic variables and first-order differential
states for an `equation_component`. The user-facing syntax is in the
[spatial TOML reference](../../docs/spatial/toml-reference.md#user-defined-equation-component).
This component uses the same canonical allocation and owned equation blocks
as the built-in tire and load components; it is not a separate solver.

For state $s_i$ and algebraic variables $a_j$, the entered equations have
the form

$$
\dot s_i=f_i(t,s,a,z), \qquad
L_j(t,s,a,z)=R_j(t,s,a,z),
$$

where $z$ includes measurements and other named model variables. An input
alias resolves to a canonical variable index before integration. Algebraic
equations are not executed in source order. They become implicit equations
$L_j-R_j=0$ and can form loops with other model equations. The count of
algebraic equations must equal the count of declared algebraic variables;
the usual Jacobian solve determines whether the equations are independent.

Each state equation contributes the implicit BDF row

$$
F_i=\dot s_i-f_i(t,s,a,z)=0.
$$

The local Jacobian contribution at a BDF step is

$$
\frac{\partial F_i}{\partial z_k}
=\alpha\delta_{ik}-\frac{\partial f_i}{\partial z_k},
$$

where $\alpha$ is the BDF derivative coefficient. The restricted expression
compiler identifies dependencies and differentiates each small local
expression with `ForwardDiff`. It only writes those dependencies to the
sparse Jacobian. Existing applied-force, applied-torque, and spanning-force
elements can refer to a component's variables in their magnitude expressions;
their existing contributions then place the load in body balance equations.

All declared differential states are marked differential, included in BDF
error control, and stored by qualified name in `.simp`. They are not
candidates for QR selection of mechanical coordinates. An algebraic variable
and its equation participate in acceleration initialization, statics, and
dynamics. The acceleration initialization solves these algebraic equations
together with loads and mechanical balance. The initial derivative of each
state is then evaluated from its state equation.

For a static solution, `static = "steady"` adds the state to the static
unknowns and solves $f_i=0$. `static = "hold"` leaves that state at its
starting value; during dynamic relaxation its equation is temporarily
$\dot s_i=0$. The normal dynamic equation is restored for dynamics and
modal linearization. The `.simp` saved-initial-condition transfer accepts
matching user states, while algebraic values are recalculated.

The first version intentionally requires one explicit first-order rate law
per state and supports continuous scalar equations only. It does not
differentiate an arbitrary expression containing `der(...)`, create
Modelica-style effort/flow connectors, or implement events and discrete
state changes. These can be added while retaining the same component-local
allocation and sparse assembly.
