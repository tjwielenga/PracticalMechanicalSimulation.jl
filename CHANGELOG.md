# Changelog

This project uses semantic version numbers for public releases. Development
history before the first release remains available in Git.

## Unreleased

No changes recorded yet.

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

- The planar modeler is the more mature program. The spatial modeler has a
  broad working library but has had less use and remains under validation.
- The included large van is a preliminary integration and performance
  model, not a validated vehicle prediction.
- Flexible bodies, detailed gear contact, and general contact search are not
  implemented in this release.
- SimpView is served locally by Julia and a web browser; it is not yet a
  packaged desktop application or hosted website.
- The methods manuscript and its larger research-verification suite remain
  works in progress and are not required to use the simulation program.
