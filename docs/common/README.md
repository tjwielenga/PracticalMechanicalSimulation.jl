# Common User Information

This section is for user-facing behavior shared by the planar and spatial
modelers.

The `.simp` result format identifies the model dimension inside the file. The
same viewing, CSV-export, and model-extraction workflow therefore serves both
modelers. Its detailed current schema is in
[Simulation Result Files](../../architecture/common/simulation-result-files.md).

Current command wrappers are:

```text
bin/simpView
bin/simpCSV
bin/simpExtract
```

The [Result Viewer](result-viewer.md) documents animation navigation, load
visibility, result and health selection, angle display, and time-history or
variable-versus-variable plots.

[Numerical Stiffness in Mechanical-System Simulation](numerical-stiffness.md)
explains why bushings, contact, tires, and other fast local behavior motivate
the shared implicit BDF solver. It develops accuracy and stability in the
$h\lambda$ plane and relates the historical analysis to the current sparse
implementation.

Dimension-specific model input and analysis instructions remain in the
[planar](../planar/README.md) and [spatial](../spatial/README.md) guides.
