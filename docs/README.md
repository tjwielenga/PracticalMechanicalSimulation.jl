# User Documentation

The user documentation is separated by model dimension. A document about
bodies, markers, elements, or model input belongs under `planar` or `spatial`;
only behavior shared by both modelers belongs under `common`.

## Planar modeler

The [Planar Modeler User's Guide](planar/README.md) documents the completed and
supported two-dimensional program. It includes the complete run workflow, TOML
reference, Julia API, and links to maintained planar examples.

Run a planar model with:

```bash
./bin/simp2d models/planar/MODEL.toml
```

## Common tools

[Common user information](common/README.md) covers facilities used by both
dimensions, including the [result viewer](common/result-viewer.md), stored
`.simp` results, the command-line result tools, and the background chapter on
[numerical stiffness](common/numerical-stiffness.md).

## Spatial modeler

The [Spatial Modeler User's Guide](spatial/README.md) documents the working
rigid-body modeler, its analysis workflows, TOML interface, Lua assemblies,
and broad element library. Its limits are stated explicitly because it has had
less validation than the planar program.
Julia programs can construct spatial models directly through the
[Sim3D Julia API](spatial/julia-api.md).

## Other reading paths

- The [Technical Manual](../architecture/README.md) contains equations,
  implementation details, and the design record.
- The [development roadmap](roadmap.md) records the public-release phase and
  deliberately deferred ideas.
- The [`paper` workspace](../paper/README.md) contains the working manuscript
  and evidence for the Fully Consistent methods paper.
