# Practical Mechanical Simulation

Practical Mechanical Simulation is a component-based Julia package for planar
and spatial mechanical-system simulation. It assembles component-local
implicit equations into complete sparse systems and supports initial-condition,
kinematic, dynamic, static, quasi-static, and modal analysis.

## Entry points

```@docs
PracticalMechanicalSimulation
PracticalMechanicalSimulation.PlanarModelIO.load_planar_model
PracticalMechanicalSimulation.SimulationRunner.run_planar_model
PracticalMechanicalSimulation.SpatialModelIO.load_spatial_model
PracticalMechanicalSimulation.SpatialSimulationRunner.run_spatial_model
PracticalMechanicalSimulation.ResultIO.write_result
PracticalMechanicalSimulation.ResultIO.read_result
Sim2D
Sim3D
```

## Reading paths

- [Getting Started](getting-started.md) covers package and source-checkout use.
- [Common User Information](common/README.md) describes stored results and
  SimpView.
- [Planar Modeler User's Guide](planar/README.md) covers the mature 2D
  implementation.
- [Spatial Modeler User's Guide](spatial/README.md) covers the 3D
  implementation and its current limitations.
- The [Technical Manual](https://github.com/tjwielenga/PracticalMechanicalSimulation.jl/tree/main/architecture)
  records equations, implementation details, and the design history.
