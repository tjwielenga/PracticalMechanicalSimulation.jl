# State Selection from Velocity Constraints

Status: common method, currently implemented in the planar modeler

## 1. Purpose

This document records a procedure for selecting individual mechanical variables as integration states. The objective is not to construct an optimal linear combination of coordinates. It is to select actual translational-velocity and body-angular-velocity components, pair them with corresponding position and pseudo-angle variables, and add state equations only for that independent set.

The common case is a holonomically constrained mechanical system. Velocity constraints are used for selection because their translational- and angular-velocity columns already have the geometric form needed for mechanical Jacobians. Position constraints remain in the implicit equation system and continue to maintain assembly.

The procedure should:

- detect the numerical rank of the constraint equations;
- identify redundant constraint rows only when redundancy exists;
- select individual dependent and independent velocity components;
- select corresponding translational positions and body-fixed pseudo angles;
- avoid using Euler-parameter components as independent states;
- preserve useful user-preferred coordinates when they are well-conditioned; and
- require only one pivoted QR factorization in the ordinary full-row-rank case.

## 2. Velocity constraint partial

Collect the body velocity coordinates in

$$
\nu=
\begin{bmatrix}
V^g\\
\omega^b
\end{bmatrix}.
$$

For an assembled system, write the velocity constraints as

$$
\dot\Phi(q,\nu,t)
=D(q,t)\nu+c(q,t)
=0,
$$

where

$$
D=\frac{\partial\dot\Phi}{\partial\nu}
\in\mathbb R^{m\times n}.
$$

The columns of $D$ correspond to individual physical velocity components. Translational columns correspond to components of $V^g$. Rotational columns correspond to components of body-referenced angular velocity $\omega^b$.

For marker constraints, the coefficient of $\omega^b$ has the same infinitesimal-rotation geometry as the position-constraint partial with respect to a body-fixed pseudo angle. This permits the angular-velocity selection to be transferred directly to pseudo-angle state selection.

## 3. Scaling before rank and coordinate selection

Numerical rank and pivot choices depend on scaling. Introduce diagonal row and column scales

$$
S_r=\operatorname{diag}(s_{r,1},\ldots,s_{r,m}),
$$

$$
S_c=\operatorname{diag}(s_{c,1},\ldots,s_{c,n}),
$$

and form

$$
\overline D=S_rDS_c.
$$

Row scaling should account for the units and characteristic magnitudes of the constraint equations. Column scaling should account for the characteristic magnitudes of translational and angular velocities. Without suitable scaling, QR may prefer one type of coordinate merely because its numerical coefficients are larger.

For each body $i$, define a characteristic length $L_i$ as the maximum
distance from its center of mass to any marker fixed on that body. If there is
no noncentral marker, use the radius of gyration $\sqrt{I_i/m_i}$. A model-level
positive floor handles a pathologically small inferred value, and the model
description may provide an explicit positive override. In the planar ordering
$[V_x,V_y,\omega]$, use column scales $[L_i,L_i,1]$. After applying these
column scales, normalize every nonzero row by its Euclidean norm. Structurally
zero rows remain zero.

The resulting $\overline D$ is used consistently for numerical rank, column
pivots, transpose-QR redundant-row detection, and preferred dependent-block
rank and conditioning. The chosen lengths and both diagonal scales are kept in
the initialization diagnostics so the decision can be reproduced.

Scaling affects the partition selected by QR but does not introduce linear combinations into the final state set. The permutation still identifies individual columns of the original physical variable vector.

## 4. Ordinary full-row-rank procedure

Compute a column-pivoted QR factorization

$$
\overline D\Pi_c
=Q
\begin{bmatrix}
R_{11}&R_{12}
\end{bmatrix}.
$$

The diagonal entries of $R$ determine the numerical rank $r$. A practical rank test compares each candidate pivot with a tolerance based on matrix size, floating-point precision, and the leading pivot. The exact tolerance is an implementation choice and should be reported in initialization diagnostics.

In the ordinary case,

$$
r=m.
$$

All velocity-constraint rows are then independent. No factorization of $D^T$ is required.

The first $m$ pivot columns identify dependent velocity components:

$$
\nu_d=\nu[\mathcal J_d].
$$

The remaining $f=n-m$ columns identify independent velocity components:

$$
\nu_i=\nu[\mathcal J_i].
$$

After applying the column permutation, the velocity constraint is

$$
D_d\nu_d+D_i\nu_i+c=0,
$$

with nonsingular $D_d$. Therefore,

$$
\nu_d=-D_d^{-1}(D_i\nu_i+c).
$$

The corresponding tangent mapping is

$$
\nu
=H\nu_i+\nu_0,
$$

where, in permuted ordering,

$$
H=
\begin{bmatrix}
-D_d^{-1}D_i\\
I
\end{bmatrix}.
$$

This mapping is useful for analysis and initialization. The implicit mechanical system need not explicitly eliminate the dependent variables; it can retain the full solution vector and add state equations only for $\nu_i$.

## 5. Selecting individual position states

Each independent velocity component selects a corresponding position-level state.

For an independent translational velocity component $V_k^g$, introduce or select the corresponding translational position $R_k^g$ and add

$$
V_k^g-\dot R_k^g=0.
$$

For an independent body-angular-velocity component $\omega_k^b$, introduce or select the corresponding body-fixed pseudo angle $\vartheta_k$ and add

$$
\omega_k^b-\dot\vartheta_k=0.
$$

The pseudo angle is an integration and infinitesimal-rotation coordinate. Its accumulated value is not used to evaluate orientation. Euler parameters retain that responsibility. The spatial orientation bridge is

$$
A=A(p),
$$

$$
\dot{\vartheta}^b=2Q(p)^T\dot p,
$$

$$
p^Tp-1=0.
$$

The selected translational positions and pseudo angles are differential variables and participate in integration-error control. The independent velocities are also differential states when their derivatives are used by the kinetic equations. All pseudo angles remain in the implicit system, but dependent pseudo angles are not included in error control. Euler parameters remain differential variables and are included in spatial orientation error control.

## 6. Why pseudo angles are preferred over Euler-parameter states

Individual Euler-parameter components are poor candidates for independent mechanical states because:

- they are coupled by the normalization constraint;
- their sensitivities to physical rotation vary with orientation;
- the pairs $p$ and $-p$ describe the same orientation;
- componentwise parameter tolerances are not direct angular-error tolerances; and
- mechanical position Jacobians are more complicated with respect to $p$ than with respect to an infinitesimal rotation.

Body-fixed pseudo angles are unconstrained tangent coordinates. The position partial with respect to a pseudo-angle component has the same geometric coefficient as the velocity partial with respect to the corresponding body-angular-velocity component. This permits reuse of component-local Jacobian information.

For an ideal revolute joint with body-fixed unit hinge direction $a^b$, the relative angular velocity is

$$
\omega_{\mathrm{rel}}^b=a^bu.
$$

Selecting the largest component of $a^b$ gives

$$
\max_k|a_k^b|\geq\frac{1}{\sqrt{3}}.
$$

Because $a^b$ is fixed in the joint-marker frame, a good initial pseudo-angle selection remains well-conditioned as the body moves globally.

## 7. Redundant constraint rows

The first QR factorization also reports the numerical rank. If

$$
r<m,
$$

the velocity-constraint set contains redundant rows. The $Q$ factor identifies independent linear combinations of constraints, but it does not identify a subset of actual component equations to retain.

Let $\mathcal J$ be the first $r$ velocity-column pivots from the original factorization and define the thin full-column-rank matrix

$$
D_J=\overline D[:,\mathcal J]
\in\mathbb R^{m\times r}.
$$

Perform a second column-pivoted QR factorization only on its transpose:

$$
D_J^T\Pi_r=Q_rR_r.
$$

The first $r$ row pivots identify a set $\mathcal I$ of actual constraint rows such that

$$
\overline D[\mathcal I,\mathcal J]
$$

is nonsingular. Retain those rows and report the remaining rows as redundant.

Thus:

- the ordinary case uses one QR factorization of $D$;
- the transpose QR is performed only when redundant rows are detected; and
- even then, it is applied to the smaller $r\times m$ transpose of an already-selected thin block, not to the entire original matrix.

This supporting rank check belongs to initialization and repartitioning. It should not dominate the time-stepping formulation or the presentation of the mechanical method.

## 8. Preferred coordinates

Mechanical knowledge may provide a more understandable state selection than an unrestricted QR pivot order. Examples include a known revolute-joint rate or a body-fixed pseudo angle aligned with a joint axis.

The implementation should allow a preferred independent set $\mathcal J_i^{\mathrm{preferred}}$. Its complementary dependent block should first be tested for rank and conditioning. If the block is acceptable, retain the preferred selection. If it is singular or poorly conditioned, fall back to pivoted QR.

Column scaling can also bias QR toward mechanically preferred variables, but an explicit preference followed by a conditioning check is easier to explain and reproduce.

## 9. Position-level compatibility check

For a holonomic constraint

$$
\Phi(q,t)=0,
$$

the velocity constraint is obtained by differentiation. With consistent translational and body-fixed infinitesimal-rotation coordinates, the position and velocity partials describe the same local constrained directions.

The velocity-selected partition should nevertheless be checked during initialization against the position constraints. The dependent position and pseudo-angle columns must provide the rank needed to solve the position assembly problem. This is a verification of the selected partition, not normally a second coordinate-selection procedure.

If this compatibility check fails, the implementation may perform separate QR factorizations at the position and velocity levels. That is a fallback rather than the default method.

## 10. Nonholonomic extension

For a genuinely nonholonomic constraint,

$$
D_v(q,t)\nu+c(q,t)=0,
$$

there may be no corresponding position constraint. The counts of independent positions and velocities can then differ:

$$
n_{q,\mathrm{ind}}=n-\operatorname{rank}(C_q),
$$

$$
n_{v,\mathrm{ind}}=n-\operatorname{rank}(D_v).
$$

In that case, separate position and velocity selections are appropriate. The position states are selected from the holonomic position constraints, while the independent velocity states are selected from the complete velocity-constraint set.

Because purely nonholonomic constraints are uncommon in the initial class of examples, this separate selection should be implemented only when required.

## 11. Conditioning and repartitioning

Initialization should record:

- the estimated rank;
- the selected independent constraint rows;
- the dependent and independent velocity indices;
- the selected translational-position and pseudo-angle states;
- the diagonal pivots of $R$; and
- a condition estimate for the dependent block.

The dependent velocity block is not normally available as a separate matrix during integration. Its coefficients are distributed through the assembled implicit-system Jacobian. Reassembling and factoring the velocity-constraint partial at every step would add expense without improving an already satisfactory state selection.

Instead, the normal integration process should monitor the conditioning and performance of the complete Newton solution matrix. Warning indicators include:

- a poor reciprocal-condition estimate from the existing numerical factorization;
- small or unstable pivots;
- excessive pivot growth or iterative-refinement corrections;
- unusually large Newton corrections; and
- deteriorating Newton convergence that is not corrected by a smaller time step.

Poor conditioning of the complete matrix does not prove that the state selection is responsible. Other causes include stiff force elements, poor equation scaling, redundant equations, and a physical mechanism singularity. If the complete solution matrix raises a warning, the initialization analysis is repeated as a diagnostic operation:

1. Reassemble the velocity-constraint partial $D$.
2. Test the current dependent block and preferred alternative selections.
3. Use pivoted QR if a satisfactory preferred alternative is unavailable.
4. Distinguish a poor state partition from constraint redundancy or a physical singularity.
5. Replace the selected state equations only when a better valid partition is found.

Repartitioning is therefore exceptional rather than part of ordinary time stepping. Hysteresis should be used to prevent repeated switching between similarly conditioned partitions.

The implemented recovery test uses the separately monitored physical predictor
errors. A stable-formula step requests QR when the maximum normalized physical
error exceeds 25 and is more than ten times the selected-state error. Only the
first such step in a contiguous episode requests selection, which supplies the
hysteresis. A singular iteration matrix or repeated corrector failures also
requests QR. If QR returns the existing selection, integration continues
without altering the implicit equations.

For an ideal revolute joint expressed in its body-fixed marker frame, the selected hinge component remains fixed and this repartitioning is not expected.

### 11.1 Localizing state equations in the implicit system

Place the selected state equations after the permanent mechanical equations. The complete implicit equation set then has the block form

$$
F(y,\dot y,t)=
\begin{bmatrix}
F_m(y,\dot y,t)\\
F_s(y,\dot y,t;\mathcal S)
\end{bmatrix}
=0,
$$

where $F_m$ contains the force, torque, constraint, orientation, and other permanent equations. The closing block $F_s$ contains the state equations selected by the current index set $\mathcal S$.

The BDF Newton matrix has the corresponding form

$$
J=
\frac{\partial F}{\partial y}
+c_j\frac{\partial F}{\partial\dot y}
=
\begin{bmatrix}
J_m\\
J_s(\mathcal S)
\end{bmatrix},
$$

where $c_j$ is the current BDF leading derivative coefficient. When the selected states change, the permanent block $J_m$ is unchanged. Only the final state-equation rows and their entries in the existing variable columns are replaced.

For example, selecting an angular state adds

$$
\omega_k^b-\dot\theta_k^\ast=0.
$$

Its Newton row has entries only in the columns for $\omega_k^b$ and $\theta_k^\ast$. A translational state equation has the same local form in the corresponding velocity and position columns.


### 11.2 Preserving the BDF history

The unreduced formulation carries accepted solution history for all variables, including variables not currently selected for integration-error control. Consequently, changing the selected state equations does not require restarting the BDF method or discarding its history.

At a repartitioning point:

- the stored values of all variables remain unchanged;
- BDF predictions and numerical derivatives remain available for all variables;
- the final state-equation rows are replaced;
- the differential-variable and error-control selections are updated; and
- the implicit equations, right-hand side, and numerical Jacobian are reevaluated.

All candidate state-equation pairs are allocated when the model is loaded, but
only the chosen pairs are active. They occupy the closing rows of the dynamic
system. A change substitutes new candidate pairs in those same local rows and
therefore does not change the system size. The new sparse pattern receives a
fresh UMFPACK symbolic and numerical factorization.

There is no physical discontinuity and no coordinate reconstruction. The selection changes which existing relationships close the implicit system and which existing variables control integration error; it does not change the mechanical configuration represented by the solution history.

## 12. Recommended default procedure

The recommended state-selection procedure is:

1. Assemble and scale the velocity-constraint partial $D$.
2. Apply column-pivoted QR to $D$.
3. Determine numerical rank from $R$.
4. If the row rank is full, retain all constraint rows and perform no transpose factorization.
5. If redundant rows exist, apply pivoted QR to the transpose of the thin selected column block and retain the identified actual constraint rows.
6. Use the pivot velocity columns as dependent components.
7. Use the remaining velocity columns as individual independent components.
8. Pair independent translational velocities with translational positions.
9. Pair independent body-angular velocities with body-fixed pseudo angles.
10. Add state equations and integration-error control for the selected pairs.
11. Retain Euler parameters for orientation evaluation, kinematic consistency, and normalization, but do not select their individual components as mechanical states.
12. Verify the velocity-selected partition against the position constraints.
13. Use separate position and velocity selections only for nonholonomic constraints or a failed compatibility check.
14. During integration, use the complete Newton matrix to trigger exceptional state-selection diagnostics.
15. If repartitioning is necessary, replace only the closing state equations, update error-control selections, and continue with the complete existing BDF history.

This procedure selects understandable physical state variables, preserves the simple component Jacobians associated with infinitesimal rotations, and keeps redundant-row processing out of the ordinary initialization path.
