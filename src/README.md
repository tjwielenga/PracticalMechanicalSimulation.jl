# Source Guide

The source is organized around the path followed by a model rather than around
one large mechanism data structure:

1. `PlanarModelIO` or `SpatialModelIO` reads the model, validates references,
   allocates canonical variables and equations, and constructs the components.
2. Initial position and velocity equations are solved before redundant
   constraints and independent states are selected.
3. `SimulationRunner` or `SpatialSimulationRunner` selects the square system
   required by the requested analysis.
4. Components evaluate their own implicit equations and add their forces or
   reactions to body balance equations.
5. `HistoricalDDASSL`, the static solver, or `ModalAnalysis` solves that
   selected system.
6. `ResultIO` stores requested samples and diagnostics in a `.simp` file.

The public entry points are re-exported by
`PracticalMechanicalSimulation.jl`. Most model users need only
`load_planar_model`, `run_planar_model`, `load_spatial_model`,
`run_spatial_model`, `write_result`, and `read_result`.
The public `Sim3D` submodule constructs the same spatial model document
directly from Julia and retains the ordinary loader and result pipeline.

## Canonical system

The canonical system contains every allocated variable and equation, including
acceleration, velocity, position, reaction, applied-load, normalization, and
dormant candidate-state quantities. A component does not know where another
component was allocated. It retains only its own canonical indices and marker
references.

`AutomaticAnalysis.jl` provides three separate layers:

- declarations describe component-local variables and equation blocks;
- `ModelLayout` maps those declarations into canonical indices; and
- `AnalysisSelection` supplies an ordered square view for one calculation.

Every canonical equation has one owner. A rigid body owns its balance and
state equations. Joints and force elements use `EquationContribution` objects
to add reactions or applied loads to those owned balances. This keeps element
formulations local without obscuring ownership of the assembled equation.

The common names used in equation callbacks are:

- `z`: the complete canonical variable vector;
- `zdot`: its BDF derivative vector;
- `equations`: the complete or selected implicit-equation vector; and
- `coefficient`: the BDF value multiplying partials with respect to `zdot` in
  the Newton matrix `G_z + coefficient*G_zdot`.

## Initial conditions and state selection

Model loading retains both the entered and corrected initial configurations.
Position and velocity consistency are established before state selection,
because the velocity-constraint partial matrix is meaningful only at a
consistent configuration.

Pivoted QR of the velocity-constraint partial matrix selects individual
physical velocities. A separate QR of its transpose identifies redundant
constraint rows. Redundant ideal constraints are removed as complete
position, velocity, acceleration, and reaction families. Candidate state
equations remain allocated, so the dynamic runner can replace a poor state
partition without rebuilding the model or discarding accepted BDF history.

## Planar and spatial rotation

Planar bodies use one physical angle. Spatial bodies use normalized
scalar-first Euler parameters to evaluate finite orientation. Mechanical
orientation partials are assembled in body-fixed pseudo-angle columns. The
Euler-parameter kinematic equations map angular velocity through to the finite
orientation, while the pseudo angles provide the local coordinates used by
state selection and the Newton correction. Their accumulated values are not
interpreted as a finite orientation.

## Reading order

For the planar method, begin with:

1. `common/AutomaticAnalysis.jl`
2. `planar/PlanarComponentAssembly.jl`
3. `planar/PlanarModelIO.jl`
4. `planar/SimulationRunner.jl`
5. `common/HistoricalDDASSL.jl`

For the spatial extension, next read:

1. `spatial/SpatialComponentAssembly.jl`
2. `spatial/SpatialModeling.jl`
3. `spatial/SpatialConstraints.jl`
4. `spatial/SpatialModelIO.jl`
5. `spatial/SpatialSimulationRunner.jl`

The mathematical derivations and element sign conventions belong in the
Technical Manual under `architecture/`. Source comments explain how those
formulations are represented and solved; they should not become a second copy
of the derivations.
