# SimpView

SimpView is the browser-based model and result viewer.

The viewer reads the native `/graphics` datasets in `.simp` HDF5 results.
Reusable meshes, sampled rigid poses, instances, signals, result choices, and
appearance are stored as ordinary HDF5 arrays rather than an opaque JSON
attachment. The representation is independent of Three.js, so another
renderer can read the same data.

Results written by the `simp2d` and `simp3d` executables receive these native
graphics automatically.

For `.simp` results only, start the development viewer. Node.js 20.19 or newer
is required; continuous integration uses Node.js 24.

```sh
cd apps/SimpViewWeb
npm install
npm run dev
```

Open the address printed by Vite, then select or drop the `.simp` file itself.

To preview results and models with one command, run:

```sh
bin/simpView
```

This starts the persistent Julia model service, waits for it to become ready,
and opens the web interface. For development, the two processes can instead be
started separately. Start the model service in one terminal:

```sh
julia --threads=auto --project=. apps/SimpViewWeb/server/start.jl
```

Then start the web interface in a second terminal as shown above. The same
**Open** button accepts `.simp`, `.toml`, and `.lua`. Model previews contain a
configuration selector for the supplied model input and the corrected,
consistent initial conditions. Lua `require` calls resolve modules from the
repository's standard `assemblies` directory. Because browsers do not reveal a
selected file's absolute path, a model that imports private modules stored only
beside that model will need the later path-aware file browser.

Opening a model also displays an **Analysis** selector with **Initial
conditions**, **Static equilibrium**, **Dynamic simulation**, and **Modal
analysis** choices. Initial conditions stop after model consistency. Static
equilibrium starts from those initial conditions and stops at equilibrium.
Dynamic simulation starts from consistent initial conditions without silently
performing static equilibrium first. If a static analysis has just completed,
its equilibrium configuration is reused when dynamic simulation is selected;
the model's declared initial velocities are then restored and made consistent.
Repeating dynamic simulation continues from the last time and state and retains
the earlier time history in the new result. The requested end time must be
later than the last simulated time.

After modal analysis completes, SimpView switches from the preceding time
history to the first mode and exposes the **Mode** selector. Modal animation
uses normalized pseudo-time and plays one full oscillation every two seconds;
it does not use the physical modal frequency as wall-clock playback speed. If
a static or dynamic run has just completed, modal analysis uses that final
configuration as its operating point. The logarithmic **Mode scale** control
enlarges or reduces the displayed mode shape and its plotted perturbations
without changing the saved eigenvector.

The default display rate is 60 frames per simulated second. The analysis runs
on a background Julia thread while the browser remains responsive. Spatial
static corrections and requested dynamic output frames appear as they are
accepted; if the slider is at the newest frame, it continues following the
analysis. Moving the slider back lets an earlier configuration be inspected
without being pulled forward by later updates.

The service writes a normal incremental `.simp` result while the analysis is
running. A solver failure therefore leaves the accepted motion and static
iterations available in the viewer. When the run completes or fails,
**Save result** downloads that self-contained file, including its native
graphics. Supported browsers open a **Save As** dialog so its directory can be
chosen; other browsers use their normal download directory. The temporary
server copy is not placed in the project tree.

Current scope:

- local file opening and drag-and-drop;
- TOML and Lua model preview through the loopback-only Julia model service;
- background model execution with end-time and frame-rate controls;
- live spatial-static corrections and dynamic output frames;
- completed and failed `.simp` result download;
- Three.js scene navigation and fit-to-model;
- a planar camera that initially faces the $x$-$y$ plane with $+y$ up, and a
  spatial $+z$-up camera that starts from an oblique view;
- continuously looping animation with play, pause, and position-slider
  controls; display frames are interpolated between stored result samples;
- a searchable assembly/component tree for selecting plot variables;
- independent X and Y variables, with time as the default X axis;
- result or mode selection when the exported document contains alternatives;
- a hierarchical graphics selector with independent visibility and
  multiplicative logarithmic scaling from `0.001x` to `1000x`; applied and
  reaction forces and torques are separate branches, and `Reset` affects only
  the selected branch while `Reset all` restores every scale to `1x`;
- a **Follow** branch in the graphics selector that lists bodies by assembly;
  selecting a body follows its translation while leaving camera rotation and
  zoom under user control. **Rotate with followed body** also carries the
  camera orientation with that body, and **Fixed scene** stops following;

The simulation arrays within the `.simp` file remain authoritative. Native
graphics are derived from them when an analysis completes or stops with a
solver failure.
