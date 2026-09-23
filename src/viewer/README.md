# Result viewer source

`ViewerData.jl` defines the presentation-independent trajectories used by both
viewer paths. `StoredResultViewer.jl` reads the model dimension stored in a
`.simp` result and reconstructs planar or spatial geometry.

`PortableViewerDocument.jl` converts the reconstructed data to the normalized
mesh, pose, instance, signal, and appearance representation consumed by the
current browser viewer in
[`apps/SimpViewWeb`](../../apps/SimpViewWeb/README.md). It also writes that
representation into native `/graphics` datasets in a `.simp` result.

SimpView is the supported viewer. This directory contains only
renderer-neutral data conversion and storage support.
