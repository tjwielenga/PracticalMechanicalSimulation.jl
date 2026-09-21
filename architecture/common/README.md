# Common Technical Manual

These documents describe mathematical and numerical machinery shared by the
planar modeler and the initial spatial vertical slice. Documents state when a
capability is implemented only for one dimension.

- [Mathematical architecture](mathematical-architecture.md) defines the
  component-local unreduced implicit system and its intended spatial form.
- [Julia DDASSL conversion](julia-ddassl-conversion.md) documents the shared
  variable-step BDF integrator, scaling, sparse Newton solution, error control,
  and event handling.
- [State selection from velocity constraints](state-selection-from-velocity-constraints.md)
  describes initial QR selection, redundancy detection, health monitoring, and
  runtime recovery.
- [Sparse modal linear analysis](modal-linear-analysis.md) develops
  linearization of the full implicit equations and the sparse shift-invert
  solution.
- [Simulation result files](simulation-result-files.md) defines the portable
  `.simp` HDF5 schema.
- [User-defined equation components](user-equation-components.md) explains how
  auxiliary algebraic variables and first-order states join either modeler's
  component-local sparse implicit system.
- [ADAMS compatibility specification](adams-compatibility-specification.md)
  defines marker order, relative-coordinate, load-sign, inertia, and unit
  conventions and records the current compatibility audit.

Dimension-specific body, marker, and element equations do not belong here.
They are documented in the [planar](../planar/README.md) and
[spatial](../spatial/README.md) sections.
