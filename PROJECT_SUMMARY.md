# Practical Mechanical Simulation Project Summary

Status: living project map

Updated: September 23, 2026

## Current state

Practical Mechanical Simulation is an open mechanical-system modeling and
simulation project written in Julia. It now includes:

- a planar rigid-mechanism modeler;
- a spatial rigid-body modeler with a broad element library;
- the Sparse Fully Consistent Modeling Method;
- a variable-step BDF integrator with sparse analytical Jacobians;
- initial-condition, kinematic, dynamic, static, quasi-static, and modal
  analyses;
- automatic state selection and redundant-constraint removal;
- TOML model input and Lua hierarchical modeling assemblies;
- portable `.simp` result files; and
- the browser-based SimpView model and result viewer.

The public source release is available, and registration in Julia's General
registry is awaiting review. The immediate priorities are maintaining the
documentation, completing package registration, and making the program easier
to discover. The methods paper remains important, but it is not a gate for
making the program available. The agreed goals and sequence are recorded in
[Project Goals and Distribution Plan](PROJECT_GOALS.md).

## 1. Purpose

The project has two connected purposes:

1. Provide a practical program for modeling and simulating mechanical systems.
2. Explain and test a formulation in which component-local implicit equations
   are assembled into a sparse, fully consistent system.

The intended readers are engineers, analysts, researchers, developers, the
author, and electronic intelligences helping any of those people. The project
therefore favors explicit terminology, plain-text documentation, executable
examples, and tests that state expected behavior.

The main public entry points are intended to be a Julia package, a public
repository, a Foundation website, SimpView, and a methods paper. These serve
different purposes and should support one another rather than being forced
into one document.

## 2. Sparse Fully Consistent Modeling Method

The supported formulation is called the **Sparse Fully Consistent Modeling
Method**, or the **Fully Consistent method** for short.

Each body, constraint, force, measurement, and generator owns local variables,
implicit equations, and Jacobian contributions. These are allocated into one
canonical unreduced system. Position, velocity, acceleration, reaction, force,
measurement, and auxiliary quantities remain explicit where their elements
need them. Sparse assembly and factorization exploit the locality without
requiring the model author to derive reduced equations or traversal rules.

The program first makes the supplied configuration and velocities consistent.
It then uses pivoted QR of the velocity-constraint partial matrix to select a
minimal set of individual physical velocity states. Corresponding positions,
body-fixed pseudo angles, or relative coordinates complete the state pairs. A
separate QR of the transpose finds redundant constraint rows. The full
canonical history is retained even though only a square active equation set is
solved.

During spatial dynamics, normalized Euler parameters represent finite body
orientation. Mechanical orientation partials are assembled with respect to
body-fixed pseudo angles. Euler-parameter kinematic equations connect those
local angular corrections to angular velocity and finite orientation.

The supported model-file formulation has constraint-derivative deficit zero.
The project uses **deficit** to mean the number of further constraint
differentiations needed to determine the highest mechanical derivatives and
reactions of interest. Alternative deficit formulations and Gear stabilization
methods remain executable paper studies, but are not options in the supported
modeling program.

The mathematical foundation is described in
[Mathematical Architecture](architecture/common/mathematical-architecture.md).
State and redundancy selection are described in
[State Selection from Velocity Constraints](architecture/common/state-selection-from-velocity-constraints.md).

## 3. Planar and spatial programs

The two modelers share numerical and result machinery but keep their elements,
input documentation, and examples separate.

| Program | Input | Command | Present role |
| --- | --- | --- | --- |
| Planar | TOML, Lua, or Julia | `bin/simp2d` | 2D mechanism program |
| Spatial | TOML, Lua, or Julia | `bin/simp3d` | 3D mechanism program |

The starting points are:

- [Planar Modeler User's Guide](docs/planar/README.md)
- [Spatial Modeler User's Guide](docs/spatial/README.md)
- [Common User Information](docs/common/README.md)
- [Technical Manual](architecture/README.md)
- [Source Guide](src/README.md)

## 4. Analyses

The program supports the following analyses.

### Initial conditions

Position and velocity consistency are established before state selection.
Weighted minimum-motion corrections normally favor moving lighter bodies over
heavier bodies. A body weight is based on mass; rotational weighting also uses
its characteristic length. Relative coordinates or named body quantities can
be imposed when the model requires a particular initial condition.

Accelerations, reactions, applied loads, and other algebraic quantities are
then solved simultaneously. A saved `.simp` configuration can be applied to a
new model by matching component names, types, variables, and variable kinds.

### Kinematic analysis

A fully driven planar mechanism with no independent states is solved at each
requested time. Motion generators provide their positions; their velocities
and accelerations are derived by automatic differentiation. An explicit
dynamic analysis may also be requested for a zero-state mechanism.

### Dynamic analysis

The dynamic program retains a minimal state set for integration but solves the
complete active implicit system at every solution step. Variable-step BDF handles
the stiff force elements that motivate the project. State selection is usually
held fixed during a successful run. Singular iteration matrices, repeated
corrector failures, or unhealthy physical-variable errors can request a new
QR selection without discarding the accepted BDF history.  Preferred states can also be specified.

### Static and quasi-static analysis

Direct static analysis sets velocity and acceleration to zero and solves force,
moment, and position-constraint equations. A mass-and-inertia regularized
correction Jacobian is available when neutral directions make the ordinary
static Jacobian singular or poorly conditioned. Mass and inertia do not enter
the static equilibrium equations themselves.

Dynamic relaxation is available for difficult starting configurations. It
uses first-order BDF pseudo-time steps, reduces derivatives after every
accepted step, and hands off to a static Newton polish only after force,
acceleration, speed, and configuration-correction tests are satisfied. The
accepted relaxation and Newton steps are stored for inspection.

Static analysis at several model times gives quasi-static continuation. The
previous equilibrium predicts the next one while time-dependent forces or
generators change.

A dynamic or modal analysis can request static-equilibrium initialization in
the same run. Declared dynamic velocities are restored after the equilibrium
configuration is found.

### Modal analysis

Sparse modal analysis linearizes the complete implicit equation set at a
consistent operating point. Algebraic definitions, force laws, constraints,
and reactions remain in the linear system. Static initialization gives the
usual stationary operating point, but equilibrium is not required when the
analyst intentionally wants a linearization about another consistent state.

The method and stored modal data are described in
[Sparse Modal Linear Analysis](architecture/common/modal-linear-analysis.md).

## 5. Modeling elements

The tables below are an index, not a replacement for the two model references.

### Planar elements

| Category | Implemented elements |
| --- | --- |
| Structure | ground, rigid body, floating-reference flexible beam, oriented marker, floating marker |
| Measurements and coordinates | span, directed distance, optional revolute angle |
| Constraint primitives | inplane, perp |
| Joints | revolute, translational, fixed |
| Ideal transmissions | coordinate coupler, gear pair, rack and pinion, elastic pulley belt |
| Forces | gravity, directed force, applied torque, torsional spring-damper, spanning force, bushing, sphere-plane contact, cam contact with rollers or flat followers, surface and joint friction |
| Motion | rotational and translational-distance generators |

The complete fields, conventions, and examples are in the
[Planar TOML User's Guide](docs/planar/toml-reference.md). The equations are in
[Planar Element Formulations](architecture/planar/planar-element-formulations.md).

### Spatial elements

| Category | Implemented elements |
| --- | --- |
| Structure | ground, rigid body, floating-reference flexible beam, oriented marker, generated floating marker |
| Measurements and coordinates | span, directed distance, optional hinge/revolute angle, inline translation |
| Constraint primitives | spherical, perp, inplane, inline, hinge, orient |
| Joints | revolute, fixed, cylindrical, translational |
| Ideal transmissions | coordinate coupler, gear pairs, rack and pinion, pulley belt |
| Forces | gravity, directed force, directed torque, joint torque, spanning force, bushing, sphere-plane contact, cam contacts with roller or flat followers, surface and joint friction, rolling tire |
| Motion | rotational, translational, and spanning motion generators |

The spatial fields and present boundaries are in the
[Spatial TOML Reference](docs/spatial/toml-reference.md). The derivations are
in [Spatial Element Formulations](architecture/spatial/spatial-element-formulations.md)
and [Spatial Rigid-Body Formulation](architecture/spatial/spatial-rigid-body-formulation.md).

## 6. Model input and hierarchy

### TOML

TOML is the direct model format. Dotted table names provide readable ownership
and hierarchy. A marker normally appears below its ground or body table, while
connections, forces, measurements, and generators reference qualified marker
or joint names. Missing optional marker position and orientation fields use
documented defaults. Degree strings such as `"30 deg"` are accepted where
angles would otherwise be awkward to enter in radians.

Models may declare parameters near the beginning of the file and reuse them in
supported numeric fields. Applied torque, spanning force, tire, and contact
laws can use scalar expressions. Expressions can reference time and exposed
element coordinates such as angle, angular velocity, span, span rate, gap, or
gap rate. Dual-number differentiation supplies their local Jacobian partials.

### Lua assemblies

Lua is used when a model needs reusable hierarchical construction rather than
one direct TOML hierarchy. An assembly module is an ordinary Lua module that
returns a function. Calling the function with a table creates bodies, markers,
elements, subassemblies, and graphics through the Sim2D or Sim3D interface.
Assemblies can call other assemblies. For example, a planar mechanism can
contain reusable linkages, while a 3D vehicle can contain suspension,
leaf-spring, stabilizer-bar, wheel, and tire assemblies.

Lua expansion produces the same ordinary planar or spatial model used by the
TOML reader. It is a model-construction layer, not a second solver. See
[Sim2D Lua Assemblies](docs/planar/modeling-assemblies.md) and
[Spatial Modeling Assemblies](docs/spatial/modeling-assemblies.md).

## 7. Numerical implementation

The central calculation path is:

```text
TOML or Lua model
    -> model validation and component allocation
    -> consistent position and velocity initial conditions
    -> redundant-row and independent-state selection
    -> active sparse implicit system
    -> dynamic, static, quasi-static, kinematic, or modal solver
    -> .simp result
    -> SimpView, CSV export, restart, or model extraction
```

The Julia DDASSL conversion uses variable-step, variable-order BDF through
order five. Its Newton matrix is level-scaled so the leading coefficients of
an equation's highest-level variables remain approximately independent of
step size. Physical integrated states control local error by default, while
positions and velocities outside the selected state set provide a health
monitor. Reaction forces are not ordinary error-control variables.

Float64 sparse CSC Jacobians are factored with UMFPACK. A symbolic
factorization is reused while the sparsity pattern remains valid. Numerical
factors are also reused for modified-Newton corrections and refreshed when
step changes, convergence behavior, events, or state reselection require it.

The implementation is described in
[Julia DDASSL Conversion](architecture/common/julia-ddassl-conversion.md).

## 8. Results and SimpView

`.simp` is the common result extension for planar and spatial models.
The extension stands for "simulation made practical." It is an
HDF5 file containing the source model, canonical catalog, body references,
requested state history, diagnostics, static convergence, state-selection
changes, optional modal data, and renderer-independent viewer graphics.

Calculations use Float64. Result histories and viewer data use Float32 by
default to control file size, but a model can request double-precision output.
Initializing a model from a saved result recalculates initial conditions.

Simulation results can be viewed in the browser-based SimpView. It supports
model inspection, consistent-IC display, static, dynamic, and modal result
display, hierarchical graphics selection and emphasis, force and reaction
emphasis, body following, plots against time or another variable, and running
an existing model and saving its results in a `.simp` file.

The format is documented in
[Simulation Result Files](architecture/common/simulation-result-files.md), and
viewer controls are in [SimpView and the Result Viewer](docs/common/result-viewer.md).

A command-line run of `simp2d` or `simp3d` with the `--output` argument can be
stopped with a single Control-C. The last accepted integrator state and pending
static progress are saved in the same file with status `interrupted`. If a
hard kill (a second Control-C) leaves a valid file marked `running`, the
`bin/simpFinalize` utility preserves the last history already flushed to disk
and marks it interrupted. These files can be viewed by SimpView to help
diagnose model problems.

## 9. Repository map

| Path | Purpose |
| --- | --- |
| `src/common` | BDF, automatic analysis, expressions, modal analysis, saved results, CSV and extraction |
| `src/planar` | Planar bodies, elements, TOML/Lua loading, assembly, analysis, and command line |
| `src/spatial` | Spatial bodies, elements, TOML/Lua loading, assembly, analysis, and command line |
| `apps/SimpViewWeb` | Browser SimpView client |
| `bin` | User command wrappers and result utilities |
| `models/planar` | Maintained planar TOML and Lua input models |
| `models/spatial` | Maintained spatial TOML, Lua, and vehicle-development models |
| `assemblies/planar`, `assemblies/spatial` | Reusable Lua modeling assemblies |
| `examples` | Explanations and historical executable studies |
| `docs` | User guides and workflow documentation |
| `architecture` | Technical manual, equations, implementation, and design record |
| `benchmark` | Generated performance models and recorded measurements |
| `test/core` | Supported-program verification |
| `test/paper` | Alternative formulations and paper evidence |
| `paper` | Methods-paper plan, evidence, and working manuscript |
| `results` | Regenerable examples, benchmarks, and scratch results |

The [Source Guide](src/README.md) gives the recommended reading order and
explains ownership of canonical equations and contributions.

## 10. Important examples and benchmarks

Good first planar examples are:

- [`constant-speed-slider-crank.toml`](models/planar/constant-speed-slider-crank.toml)
  for a fully driven mechanism;
- [`torque-driven-four-bar.toml`](models/planar/torque-driven-four-bar.toml)
  for closed-loop dynamics;
- [`nonlinear-expression-pendulum.toml`](models/planar/nonlinear-expression-pendulum.toml)
  for an expression torque;
- [`bouncing-ball.toml`](models/planar/bouncing-ball.toml) for one-sided contact;
  and
- [`modal-pendulum.toml`](models/planar/modal-pendulum.toml) for modal analysis.

Good first spatial examples are:

- [`free-rotating-body.toml`](models/spatial/free-rotating-body.toml) for finite
  orientation and free flight;
- [`revolute-pendulum.toml`](models/spatial/revolute-pendulum.toml) for the
  basic joint and relative angle;
- [`static-inverted-spherical-pendulum.toml`](models/spatial/static-inverted-spherical-pendulum.toml)
  for dynamic relaxation;
- [`bevel-gear-pair.toml`](models/spatial/bevel-gear-pair.toml) for nonparallel
  revolute axes;
- [`out-of-plane-four-pulley-belt.toml`](models/spatial/out-of-plane-four-pulley-belt.toml)
  for spatial belt geometry;
- [`steered-tire-test-rig.toml`](models/spatial/steered-tire-test-rig.toml) for
  combined tire forces; and
- [`large-van.lua`](models/spatial/large-van.lua) and
  [`large-van-high-cg.lua`](models/spatial/large-van-high-cg.lua) for a
  full-vehicle rollover comparison with roof-ground contacts.

The [Planar Performance Benchmarks](benchmark/README.md) cover open pendulum
chains, closed parallelogram chains, impact banks, and rotor trains. They record
the change from generic sparse LU to UMFPACK, symbolic and numerical
factorization reuse, state-coordinate comparisons, initialization cost, and
scaling with model size. The measurements are development-machine baselines,
not universal performance claims.

## 11. Common commands

From the repository root, run and save one model:

```bash
./bin/simp2d models/planar/torque-driven-four-bar.toml \
    --output results/examples/planar/torque-driven-four-bar.simp --overwrite

./bin/simp3d models/spatial/revolute-pendulum.toml \
    --output results/examples/spatial/revolute-pendulum.simp --overwrite
```

Start SimpView and open a `.simp`, `.toml`, or `.lua` file:

```bash
bin/simpView
```

Export selected histories to CSV:

```bash
bin/simpCSV RESULT.simp RESULT.csv \
    body.R_x body.V_x
```

Recover the embedded model:

```bash
bin/simpExtract RESULT.simp recovered.toml
```

Finalize a valid result left `running` by a hard termination:

```bash
bin/simpFinalize RESULT.simp
```

Run focused verification:

```bash
julia --project=. test/runtests.jl planar
julia --project=. test/runtests.jl spatial
julia --project=. test/runtests.jl ddassl
```

Run the default maintained suite with:

```bash
julia --project=. test/runtests.jl
```

The separate paper studies are run with:

```bash
julia --project=. test/paper/run_all.jl
```

See [Test Suites](test/README.md) before using paper-study results as current
publication evidence.

## 12. Verification status

The most recent focused maintained-suite checks completed as follows on the
development machine:

| Suite | Checkpoint | Checks | Result |
| --- | --- | ---: | --- |
| Planar | September 22, 2026 | 1,181 | passed |
| Spatial, including assembly expansion | September 23, 2026 | 1,582 | passed |

These counts record the latest work, not a permanent release qualification.
The viewer and DDASSL groups were not rerun as part of that focused checkpoint.
A public-release audit should run the complete maintained suite, record Julia
and dependency versions, and reproduce a selected benchmark set.

The paper suite is intentionally separate. It preserves older formulations and
comparisons, and should not force the supported source to remain compatible
with every historical implementation. Its numerical baseline must be settled
before final publication comparisons are claimed.

## 13. Known limitations

- Planar and spatial floating-reference beams provide compact linear
  Timoshenko flexibility. They assume small elastic deformation and fixed
  reference mass and inertia; general nonlinear flexible bodies are not
  implemented.
- General spatial surface-to-surface contact and automatic contact search are
  not implemented. Current contact elements include compliant sphere-to-plane
  and cam-follower contact; the rolling tire supplies its own longitudinal and
  lateral force model.
- The preliminary vehicle model still needs validation through tire lift-off,
  stiff bumper engagement, rollover, and long dynamic runs.
- Gear pairs are ideal kinematic transmissions. Pressure angle, tooth
  compliance, backlash, and load-dependent contact-side switching are
  deferred.
- Model quantities normally use consistent SI units. There is no general unit
  algebra system, although degree strings are supported for angle entry.
- Lua assemblies are useful and composable. Their planar interface is new and
  needs use on larger assemblies; both dimensions should retain the same
  validation conventions as their libraries grow.
- Dense pivoted QR is used for initial state and redundant-row selection. It
  has not been a bottleneck in present models, but very large closed-loop
  systems may eventually justify a sparse rank-revealing method.
- Allocation in long dynamic simulations remains substantial even though
  sparse factorization performance is good.
- Substantial speedups may still be available with further benchmarking.
- Registration in Julia's General registry is awaiting review. TagBot will
  create the release tag after registration is accepted.

## 14. Publication and distribution

The intended distribution has several parts:

- an open Julia package for installation and reuse;
- a public source repository for development and verification;
- a Foundation website as a stable public introduction;
- SimpView for direct inspection of models and results; and
- an ASME methods paper explaining and evaluating the Fully Consistent method.

The intended paper venue is the ASME *Journal of Computational and Nonlinear
Dynamics*. The working manuscript is
[Sparse Fully Consistent Modeling Method](paper/sparseFullyConsistentMethod.md).
The [paper workspace](paper/README.md) also contains the outline, claims and
evidence map, and planar-methods consolidation.

The working manuscript is Markdown. Word and PDF are review or submission
artifacts produced at checkpoints. The author prefers plain, direct technical
writing rather than an academic style.

## 15. Current priorities and deferred work

The current sequence is:

1. Maintain this project summary as the central map.
2. Complete registration in Julia's General registry and create the release
   tag through TagBot.
3. Continue improving installation, first-run instructions, examples, and
   reproducible verification as users encounter them.
4. Add a Foundation website page linking the program, documentation,
   repository, examples, and SimpView.
5. Learn from initial users and incorporate their experience.
6. Finish the methods paper using the released program and reproducible
   evidence.

Deferred modeling work includes further vehicle diagnosis, more realistic gear
contact, general spatial contact, nonlinear flexible bodies, and performance
work justified by larger models. These are valuable extensions, but they
should not prevent documentation and public release of the program that
already exists.

## 16. Glossary

**Active equation or variable**

A quantity included in the solution system for the current
analysis. Inactive quantities remain named in the catalog and result format.

**Canonical system**

The complete collection of component variables and implicit
equations. Equations and variables are selected from the canonical set
for each particular analysis.

**Constraint-derivative deficit**

The number of further constraint differentiations needed to determine the
highest mechanical derivatives and reactions of interest. The supported Fully
Consistent model-file formulation has deficit zero.

**Consistent initial conditions**

A configuration, velocity set, acceleration set, reactions, and force variables
that satisfy the applicable implicit equations at the starting time.

**Dynamic relaxation**

A static-equilibrium method that uses dynamics equations but damps
velocities and accelerations and advances using first-order BDF pseudo-time
steps but fixed physical model time.

**Floating marker**

A marker owned by one body whose instantaneous location is defined by geometry
on another body or carrier. Gear, rack, belt, and reaction-force elements use
floating markers to apply loads at the correct moving point.

**Fully Consistent method**

Short name for the Sparse Fully Consistent Modeling Method.

**Implicit equation**

An equation which is defined to be equal to zero.

**Measurement**

A reaction-free element that defines useful output such as span distance,
directed distance, velocity, or acceleration.

**Physical state**

An individually selected velocity and its position level partner.
They correspond to mechanism degrees of freedom.
Their derivatives are included in the solution set and are
linked to solution variables to make a square solvable system.
The BDF integrator estimates error using these variables.

**Pseudo angles**

Body-fixed local angular coordinates used for spatial Newton partials, state
selection, and integration coupling. Their derivatives correspond to the
body-fixed components of angular velocity, but their accumulated values are
not interpreted as the body's finite physical orientation. Normalized Euler
parameters represent that orientation.

**Relative coordinate**

An optional angle or distance owned by a joint or measurement and made
available for output, expressions, couplers, initial conditions, and state
selection.

**SimpView**

The browser-based application for opening models and `.simp` results, running
analyses, animating motion or modes, controlling graphics, and plotting
canonical variables.

**`.simp` result**

The HDF5 result file shared by planar and spatial programs. It contains the
model, catalog, histories, diagnostics, modes when present, and native viewer
data.
