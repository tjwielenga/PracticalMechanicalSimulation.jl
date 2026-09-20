# Spatial Source

This directory contains the three-dimensional modeler. It remains separate
from `src/planar` so planar and spatial rotation conventions cannot be
mistaken for one another. It currently provides:

- a free rigid body with explicit translational and angular accelerations;
- physical rotation matrices evaluated from normalized Euler parameters;
- body-fixed pseudo angles accompanying selected angular-velocity states;
- ground and body-fixed oriented markers;
- spherical joints with global reaction vectors;
- perpendicular-axis constraints with configuration-dependent reaction
  directions;
- inplane constraints with moving plane normals and point reactions;
- inline constraints composed from two inplane primitives, with optional
  axial distance, velocity, and acceleration coordinates;
- hinge orientation constraints assembled from two perpendicular-axis
  primitives;
- revolute joints assembled from spherical and hinge primitives;
- optional continuous hinge angle, angular-velocity, and angular-acceleration
  coordinates eligible for state selection;
- linear coordinate couplers among two or more hinge rotations and inline
  translations, including mixed screw relations and physical reactions;
- ideal parallel-axis and bevel gear pairs with automatic pitch-tangent
  geometry, signed effective radii, and floating contact reactions;
- ideal spur rack-and-pinion sets coupling an inline distance to a revolute
  angle with generated floating contact reactions;
- mass-weighted position and velocity consistency correction, scalar body
  weight scales, and optional imposed body or relative initial conditions;
- runtime QR state reselection requested by DDASSL recovery diagnostics;
- pivoted-QR detection and complete three-level removal of redundant scalar
  ideal-constraint families;
- reaction-free span and directed-distance measurements with explicit
  distance, velocity, and acceleration outputs;
- marker-directed applied forces with constant or expression magnitudes and
  optional floating reaction markers;
- six-component spatial bushings with $z$-$y$-$x$ Bryant angles, exact
  rotating-frame rates, and coincident floating reactions;
- tangential bristle friction on compliant sphere-plane contact, with static
  anchoring, slip-dependent force capacity, and two carried shear states;
- revolute bearing friction driven by the joint's radial reaction, with one
  carried angular shear state and equal-and-opposite friction torques;
- axial translational friction driven by an inline guide's two transverse
  reactions, with one carried shear state and equal-and-opposite forces;
- two-direction tangential friction on a bilateral inplane primitive, using
  its signed normal reaction's magnitude and two carried shear states;
- steady-state rolling tires with planar contact kinematics, expression-based
  longitudinal and lateral slip forces, and a combined-slip friction ellipse;
- marker-to-marker spanning forces with explicit geometry, rate, scalar-force,
  and global-force variables;
- gravity;
- direct static equilibrium, quasi-static continuation, dynamic relaxation,
  and static initialization of dynamics;
- a spatial TOML reader and static and dynamic runners; and
- sparse analytical Jacobian contributions.

`SpatialModelIO.jl` owns the spatial TOML vocabulary. `SpatialComponentAssembly.jl`
defines body variables and implicit equations, `SpatialModeling.jl` defines
marker geometry, `SpatialDirectedDistances.jl` supplies the common oriented-axis
and point-to-plane kinematics, `SpatialSpans.jl` supplies shared
marker-to-marker distance kinematics, and `SpatialConstraints.jl` contains
ideal joints.
`SpatialAppliedForces.jl` contains applied-force elements,
`SpatialCoordinateCouplers.jl` contains scalar relative-coordinate couplers,
`SpatialGearPairs.jl` contains ideal spatial gear geometry and reactions,
`SpatialRackAndPinions.jl` contains ordinary spur rack-and-pinion geometry
and reactions,
`SpatialBelts.jl` contains oriented pulley pitch circles, common-tangent
geometry, elastic no-slip spans, and pulley loads,
`SpatialBushings.jl` contains the spatial bushing,
`SpatialPlaneContacts.jl` contains one-sided compliant sphere-plane contact,
`SpatialFrictionForces.jl` adds tangential friction to that contact,
and revolute-bearing friction to ideal revolute joints,
and axial guide friction to spatial inline constraints,
and tangential friction to spatial inplane constraints,
and `SpatialTires.jl` contains rolling-tire contact kinematics, load equations,
and event surfaces.
`SpatialSimulationRunner.jl` provides the static solver and connects dynamic
models to the shared DDASSL implementation. It also converts mechanical
orientation partials from Euler-parameter columns to local pseudo-angle
columns for the dynamic Newton matrix. Other spatial joints and contacts with
general surface geometry remain to be implemented.
