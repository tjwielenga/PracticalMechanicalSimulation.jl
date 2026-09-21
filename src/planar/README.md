# Planar Source

This directory contains the implemented two-dimensional modeler: planar body
and force definitions, component assembly, TOML loading, analysis execution,
and command-line handling. The result viewer is shared with spatial models and
therefore lives under `src/viewer`.

These files remain direct submodules of `PracticalMechanicalSimulation`, and
the public planar API is unchanged by the directory organization.

`PlanarModelIO.jl` is the model-building coordinator. It performs consistency
correction, redundant-row detection, and state selection before it returns a
`LoadedPlanarModel`. `SimulationRunner.jl` does not rebuild that model. It
chooses active canonical rows and columns for kinematic, dynamic, static,
quasi-static, or modal analysis and supplies them to the common solvers.

`PlanarComponentAssembly.jl` is the central element implementation file. A
component normally has four pieces: a declaration of its local variables and
owned equations, an allocated structure holding canonical indices, executable
callbacks for its owned equations, and additive contributions to body
balances. Small construction helpers in `PlanarModeling.jl` connect those
pieces after the complete layout is known.

`PlanarDirectedDistances.jl` contains the common directed-axis and signed
marker-distance kinematics used by inplane constraints, distance coordinates,
translational motion, plane contact, and marker-directed applied forces.

`PlanarFrictionForces.jl` adds carried-shear friction to revolute joints,
translational joints, and standalone inplane constraints. Its scalar bristle
law is shared with Sim3D through `common/FrictionLaws.jl`.

Marker-to-marker span geometry and velocity equations are shared by the
reaction-free `span` measurement and the `spanning_force` element in
`PlanarComponentAssembly.jl`.
