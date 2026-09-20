# Spatial Technical Manual

The spatial modeler is implemented independently from the mature planar
element library. It contains free rigid bodies, oriented markers, constraint
primitives and composite joints, gravity, applied forces and torques, spanning
forces, bushings, sphere-plane contacts, rolling tires, coordinate couplers,
gear pairs, rack-and-pinion sets, and planar or out-of-plane pulley belts. It
also contains TOML and Lua model loading, sparse analytical Jacobians,
mass-weighted initial-condition correction with optional imposed values, BDF
integration, direct and relaxed static equilibrium, quasi-static continuation,
static initialization of dynamics, saved-result initialization, sparse modal
analysis, `.simp` storage, and shared result viewing.
Consistent redundant ideal constraints are detected by pivoted QR and removed
as complete three-level constraint families before state selection.

The [Spatial Rigid-Body Formulation](spatial-rigid-body-formulation.md)
documents the variables, frames, equations, and orientation representation now
implemented. The general [Mathematical Architecture](../common/mathematical-architecture.md)
remains the design basis for future spatial joints and forces.
The [Spatial Element Formulations](spatial-element-formulations.md) records the
joint and spanning-force equations, composite construction, and reaction
conventions.
The [User-defined Equation Components](user-equation-components.md) chapter
describes how auxiliary algebraic equations and first-order states join the
same sparse system.

The implemented source path is:

```text
spatial TOML model
   -> SpatialModelIO
   -> AutomaticAnalysis allocation
   -> SavedInitialConditions result transfer when requested
   -> SpatialModeling, SpatialComponentAssembly, SpatialDirectedDistances,
      SpatialConstraints, SpatialSpans, SpatialBelts,
      SpatialRackAndPinions, and SpatialAppliedForces
   -> SpatialSimulationRunner, HistoricalDDASSL, and ModalAnalysis
   -> ResultIO
```

The common `SavedInitialConditions.jl` module performs the same result-sample,
name, type, and variable matching for both modelers. Spatial configurations
transfer physical orientation through normalized Euler parameters; body-fixed
pseudo angles are deliberately recalculated rather than transferred.

The spatial program can run substantial mechanisms and preliminary vehicle
models. It has had less validation than the planar program, especially for
difficult contact and tire-lift-off behavior. Each new joint or force must be
documented here when its local equations and tests are added.
