# Changelog

This project uses semantic version numbers for public releases. Development
history before the first release remains available in Git.

## Unreleased

### Added

- Julia model-building APIs for planar and spatial models, including
  hierarchical construction examples.
- User-defined algebraic and differential equation components for both
  modelers.
- Planar and spatial friction elements for surfaces, revolute joints, and
  translational joints, together with inplane friction.
- Planar and spatial cam contacts with roller and flat followers.
- Simplified floating-reference Timoshenko beams for planar and spatial
  models.
- Planar Lua assemblies and stage-dependent planar forces.
- Spatial cylindrical and translational joints and a marker-directed torque.

### Changed

- Improved the spatial rolling-tire model with load-dependent bristle
  behavior.
- Reorganized SimpView graphics controls, improved plotting and camera
  following, and added flexible-deformation display controls.
- Standardized the command names under `bin` and added package documentation,
  selective CI workflows, and coverage reporting.

## 0.1.0 - 2026-09-16

First public release of Practical Mechanical Simulation.

### Included

- Planar and spatial rigid-mechanism modeling with component-local implicit
  equations assembled into a sparse, unreduced system.
- Initial-condition, kinematic, dynamic, static, quasi-static, and sparse
  modal analyses using the Sparse Fully Consistent Modeling Method.
- A variable-step BDF integrator with analytical sparse Jacobians, automatic
  state selection, redundant-constraint removal, restart handling, and
  readable partial results after interrupted or failed runs.
- TOML model input, reusable Lua spatial assemblies, command-line programs,
  and a broad planar and spatial element library.
- Versioned compressed `.simp` HDF5 results with single-precision output by
  default, CSV and embedded-model extraction, and saved-result initialization.
- SimpView for inspecting models, running analyses, viewing time and modal
  animation, plotting signals, following bodies, and scaling or hiding
  graphical categories.
- MIT licensing, contribution and security guidance, software citation
  metadata, and automated Julia and web-viewer checks.
- A separate paper-test environment for DASSL, OrdinaryDiffEq, and Sundials
  comparisons, keeping research-only packages out of normal installation.

### Known limitations

- The included large van is a preliminary integration and performance
  model, not a validated vehicle prediction.
- Flexible bodies, detailed gear contact, and general contact search are not
  implemented in this release.
- SimpView is served locally by Julia and a web browser; it is not yet a
  packaged desktop application or hosted website.
- The methods manuscript and its larger research-verification suite remain
  works in progress and are not required to use the simulation program.
