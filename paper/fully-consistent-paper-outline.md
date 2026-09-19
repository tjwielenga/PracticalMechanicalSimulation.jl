# Sparse Fully Consistent Modeling Method

Status: working paper outline

Target journal: ASME *Journal of Computational and Nonlinear Dynamics*

## 1. The problem

- Practical mechanical systems are often numerically stiff.
- Stiff integration requires Newton-like iteration and an efficient Jacobian
  solution.
- Coordinate reduction can make the equations and Jacobian complicated.
- Component-local implicit equations are simple to formulate, but produce a
  large system.
- Sparse matrix methods make that large system practical.

## 2. Ways to formulate a mechanism

Use the planar pendulum as the common example.

- Reduced coordinates
- Displacement constraints: Deficit Two
- Velocity constraints: Deficit One
- Acceleration constraints: Deficit Zero
- Baumgarte stabilization
- GearStableV
- GearStableA
- Fully Consistent formulation

Explain the strengths and weaknesses of each method. Introduce Constraint
Deficit and distinguish it from mathematical DAE index.

## 3. The Fully Consistent equations

- Positions, velocities, accelerations, reactions, and component variables
  remain explicit.
- Position, velocity, and acceleration constraint equations are all retained.
- Only the independent physical states receive differential equations.
- Acceleration remains an unknown.
- Bodies, joints, and forces contribute local equations and Jacobian partials.
- Assemble those contributions into one sparse system.

Include the pendulum equation-and-variable table here.

## 4. Preparing the equations for solution

- Correct the supplied positions to satisfy constraints.
- Correct velocities when necessary.
- Select independent physical states using column-pivoted QR of the velocity
  constraint matrix.
- Detect redundant constraint equations using QR of the transposed constraint
  matrix.
- Apply preferred state choices when they preserve rank.
- Construct the active square equation system.
- Solve for accelerations, reactions, and force variables.

## 5. Scaling

- Characteristic length, mass, and velocity
- Scaling translational and rotational quantities
- Equation and variable levels
- Step-size-dependent BDF coefficients
- Keeping the scaled Jacobian from changing unnecessarily as the step size
  changes

Explain why scaling helps both conditioning and factorization reuse.

## 6. Integration

- Variable-step, variable-order BDF
- Prediction from the history polynomial
- Newton correction of the complete system
- Error control on physical positions and velocities
- Sparse symbolic and numerical factorization
- Reusing numerical factors during modified Newton iteration
- Rejected steps and restarts when force stiffness changes
- State reselection as a recovery method

## 7. Other analyses

Show how the same assembled equations support:

- consistent initial conditions;
- kinematic analysis;
- dynamic analysis;
- static equilibrium;
- dynamic relaxation;
- consecutive static solutions;
- modal analysis; and
- continuation from saved results.

The important point is that these analyses do not require separate component
formulations.

## 8. Planar and spatial mechanics

- Planar coordinates and rotations
- Spatial orientation using Euler parameters
- Pseudo-angle Jacobian columns
- Angular velocities and accelerations
- Marker frames
- Virtual-power application of forces and reactions
- The same component assembly in two and three dimensions

Demonstrate that the method is not inherently planar without turning this
section into a catalog of spatial elements.

## 9. Examples and results

Organize the evidence by what it tests:

- Pendulum: equation formulations and accuracy
- Serial pendulum chains: open-system sparse scaling
- Parallelogram chains and four-bars: closed loops and redundant constraints
- Rotor trains: smooth numerical stiffness
- Bouncing balls: changing force stiffness
- Static and modal examples
- A spatial mechanism
- The large van: a practical spatial model containing stiff forces,
  tires, suspensions, steering, static equilibrium, dynamics, and modal
  analysis

The large van should be evidence of breadth rather than the principal
performance benchmark.

## 10. What the results show

- A large unreduced system need not mean a dense solution.
- Sparse factorization is essential.
- Symbolic and numerical reuse materially reduce work.
- Dense QR is inexpensive for the models tested, although it may eventually
  limit very large closed-loop systems.
- Physical state choices affect performance.
- Compliant contact remains difficult but can be handled.
- Component-local equations work for both planar and spatial systems.

## 11. Limitations and future work

- Large systems with many closed loops
- Better sparse rank detection
- Improved criteria for state reselection
- More testing of massless bodies
- Discontinuous contact and friction
- Automatic differentiation and Jacobian-free alternatives
- Flexible bodies
- More complete spatial validation

## 12. Conclusions

Return directly to the main question: a mechanical system can retain all
physical variables and all useful constraint levels, integrate only a minimal
set of states, and still be solved efficiently using sparse methods.

The paper is primarily about the method. The simulation program and its models
provide evidence rather than turning the paper into a software manual.
