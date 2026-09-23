# SimpView and the Result Viewer

The browser-based SimpView is the current interface. Start it from the
repository root with:

```bash
bin/simpView
```

Its **Open** button accepts `.simp`, `.toml`, and `.lua` files. A model can be
inspected before analysis and run through initial conditions, static
equilibrium, dynamics, or modal analysis. Static does not continue
automatically into dynamics. A later dynamic selection reuses a completed
static configuration when one is available, restores the model's declared
initial velocities, and makes them consistent before integration. Otherwise,
dynamics starts from the model's consistent initial conditions without a
static solve. A repeated dynamic run continues from its last accepted time and
state, preserving the earlier samples in the resulting file. Modal completion
selects the first mode and exposes the mode selector instead of retaining the
preceding time-history display. Mode shapes use normalized pseudo-time and play
one complete oscillation every two seconds, independent of modal frequency.
The viewer interpolates positions, orientations, and graphical scales between
stored samples so browser refresh frames do not step through the modal cycle.
If static or dynamic analysis has already completed in the current model
session, modal analysis uses that final configuration as its operating point.
The logarithmic **Mode scale** slider enlarges or reduces the displayed mode
shape and its plotted perturbations. It changes only the display and does not
alter the calculated or stored eigenvector.

Dynamic runs accept **End time** and **Frames/s** overrides and can be followed
while output frames are produced. Spatial static corrections are displayed as
they are accepted. **Save result** downloads the completed or failed
self-contained `.simp` file. It opens a **Save As** dialog in browsers that
provide one and otherwise uses the browser's download folder.

For a command-line run with `--output`, pressing **Control-C once** requests a
cooperative stop. The integrator saves its exact last accepted state, pending
static progress is flushed, and the same `.simp` file is finalized with status
`interrupted`. It can then be opened normally to examine the motion leading to
the stop. A second interrupt or an external hard kill cannot run this cleanup.
If such a kill leaves a valid file marked `running`, finalize the samples that
had already reached disk with:

```bash
bin/simpFinalize RESULT.simp
```

That repair command cannot recover a state that existed only in process
memory. The ordinary cooperative stop is therefore preferable.

The model, analysis histories, native graphics, and variable catalog are read
directly from the HDF5 result file. A failed command-line analysis leaves the
same requested `.simp` file in place. SimpView labels it as failed and displays
the requested frames through the last accepted state, which makes the motion
leading to a failure available for diagnosis.

When an analysis is started from SimpView, accepted output frames and
spatial static corrections are displayed while they are calculated. Completed
spatial static results retain these convergence states, so reopening the
`.simp` file replays the relaxation and Newton iterations rather than showing
only the final equilibrium. The first two samples are **Model input**, showing
the body positions and orientations supplied by the model, and **Consistent
initial conditions**, showing the configuration after the joint equations have
been satisfied. They are published before static equilibrium begins, so they
remain available if the later static solution fails.

Dashed vertical lines in the combined convergence plot separate these initial
configurations, the relaxation cycles, and Newton polish. The live title
identifies each phase separately. It also reports the current static body-force
imbalance, body-torque imbalance, and position-constraint error, so the
relaxation's progress toward equilibrium can be distinguished from DDASSL's
integration accuracy.

## Animation controls

| Action | Control |
|---|---|
| Rotate the mechanism | Left-drag in the animation; horizontal motion orbits around global $z$ |
| Translate it sideways or vertically | Shift+left-drag, or right-drag |
| Zoom | Scroll in the animation |
| Restore the original camera | Use **Fit view** |
| Start or stop time playback | **Play/Pause** |
| Inspect one sample | Move the time slider |
| Keep a translating mechanism in view | Select its main body under **Follow** |

Planar results initially face the global $x$-$y$ plane. Spatial results start
from an oblique three-dimensional view. Camera changes affect only the display.
SimpView uses global $+y$ as the initial upward screen direction for a
planar model and global $+z$ as up for a spatial model.
The graphics hierarchy begins with a **Follow** branch listing bodies by assembly.
It normally starts at **Ground**, which keeps the camera in the fixed global
frame. Selecting a moving body translates the camera
and its orbit center with that body's center of mass while preserving camera
rotation and zoom. It is useful for a vehicle or other mechanism that travels
far from its initial position; it does not rotate the camera with the body.

The remaining selector follows the model rather than the renderer. **Model**
expands first into assemblies and subassemblies using their hierarchical names.
Within each assembly, **Bodies** contains body geometry, inertia displays,
marker frames, and flexible deformation controls; **Joints** contains ideal
connections and their primitive symbols; and **Forces** contains applied loads,
reactions, torques, bushings, and other force-element graphics. Selecting an
assembly therefore lets the common logarithmic scale control resize everything
owned by that assembly, and its checkbox shows or hides the whole assembly.
Parent and child scales multiply.

All branches start collapsed. Expand only the assembly being inspected. Lua
assemblies can place an element whose public name lies outside its natural
namespace by setting `graphics.assembly` to the owning assembly name. SimpView
then places all graphics for that element under the declared assembly instead
of guessing from the element name.

A flexible body has a **Deformation** leaf. Selecting it changes the same slider
to **Amplification**. This magnifies only elastic displacement and rotation
relative to the body's floating reference frame. It does not magnify rigid-body
motion or member thickness, and **Reset** restores the physical $1\times$
shape.

The checkboxes show or hide reaction loads, applied loads, torques, loads on
ground, measurement symbols, ordinary marker frames, joint graphics, and
inertia ellipsoids.
**Markers** does not hide joints; **Joints** controls the mechanical joint and
bushing symbols separately. A perpendicular-axis primitive shows its two
constrained directions as a heavy light-gray two-axis frame. A hinge is drawn
as a light-gray pin along the second marker's $z$-axis. Every spatial
point-constraint sphere uses the same nominal diameter, whether it belongs to
a spherical, fixed, inline, inplane, or revolute connection. This common
sphere diameter is the base dimension for every inferred ideal-joint symbol.
Each Perp axis is 2.5 sphere diameters long. The hinge pin is 3.75 sphere
diameters long and its diameter is the sphere diameter divided by 1.2. A
revolute combines the hinge pin and spherical point symbols because it
combines those two constraints; its sphere is 20 percent larger than the pin
diameter so both parts remain visible. A
bushing is drawn as a contrasting black cylinder along its second marker's
$z$-axis. An
inferred spatial marker is drawn as a small oriented $x$-$y$ frame rather
than a sphere. **Force scale** uniformly scales the length, shaft diameter, and
head of every force arrow. **Torque size** uniformly scales the square shaft
and both cones. Both symbols therefore retain one fixed shape at every
magnitude. Their linear dimensions are proportional to the square root of the
force or torque magnitude, so the area of an arrowhead represents the load.
This square-root scale keeps small loads visible when the same model also
contains much larger loads.
Each joint or force element that supplies an inferred symbol has a
**Default graphic** item below its name in the graphics tree. It controls the
joint pin, point-constraint sphere, bushing cylinder, spring connector, or
contact geometry independently of the element's applied and reaction arrows.
For example, a plane-contact force's sphere and plane share one **Default
graphic** item, so they can be enlarged or hidden even when the contact force
is zero. The tree controls use logarithmic scaling from $0.001$ through $1000$
times nominal size, allowing millimeter-scale vehicle connections to be
inspected without hiding other model graphics. None of these controls changes
the calculated values.

When possible, the bushing cylinder length is estimated from its bending and
translational stiffnesses. Two translational springs, each with half of the
combined stiffness $k$ and separated by $L$, give

$$
k_\theta=\frac{kL^2}{4},\qquad
L=2\sqrt{\frac{k_\theta}{k}}.
$$

For a cylinder along local $z$, the viewer uses the available
$(k_{\theta x},k_y)$ and $(k_{\theta y},k_x)$ pairs. It combines compatible
values and bounds the displayed length relative to the model size. A bushing
without usable translational and bending stiffness, such as a rotational-only
leaf bending bushing, receives a reasonable default symbol.

Each visible spatial body also has a translucent inertia ellipsoid centered at
its center of mass. **Inertia ellipsoids** turns all of them on or off without
hiding the bodies. The eigenvectors of the body inertia tensor orient the
ellipsoid. If its principal moments are $J_1$, $J_2$, and $J_3$, the unscaled
axis shapes are obtained from

$$
a_1^2=J_2+J_3-J_1,\qquad
a_2^2=J_1+J_3-J_2,\qquad
a_3^2=J_1+J_2-J_3.
$$

This complementary-moment form makes the long axis of a slender link lie
along the link rather than perpendicular to it. The axis shape is normalized,
then the overall linear size is made proportional to the cube root of body
mass. Thus ellipsoid volume represents mass while shape and orientation
represent the inertia tensor. It is a mass-property symbol, not the physical
surface of the body. Setting a body's graphic `visible = false` also hides its
inertia ellipsoid.

When a file contains more than one analysis result or modal mode, the selector
at the upper-left of the animation chooses the result being animated. Keeping
this selector with the animation separates it from the hierarchical X/Y plot
browser. An **Event** selector, when present, moves directly to a stored initial
configuration, static phase boundary, or solver-health event.

## Plotting histories

The **X** and **Y** buttons show the variables currently assigned to the plot
axes. **X** defaults to `time (s)`, producing the usual history plot. Press
either button to choose the axis to change. The shared browser then opens at
the selected variable's component. Entries ending in `/` are assemblies,
subassemblies, bodies, or elements; select them to descend through the model
name hierarchy. **Up** moves to the parent and **Root** returns to the top.
Entries without `/` are variables and assign that history to the chosen axis.
`time (s)` is available at the root when choosing X. The red point marks the
animation's current sample.

For example, the steered tire test rig can show its combined-friction path by
using the browser to select:

```text
X: tire.longitudinal_force
Y: tire.lateral_force
```

Changing either selector rescales both plot axes to the selected data. A
spatial static result that contains convergence history can plot the force and
torque imbalance, constraint error, and equivalent acceleration through the
stored iteration sequence.

Named angular displacements are labeled with `(deg)` and displayed in degrees.
This includes body and joint angles, tire slip and camber angles, and bushing
angles. Angular velocity and acceleration retain their stored units. The
underlying `.simp` histories remain in radians, so CSV export and model restart
are unaffected by this display conversion.

## Force and torque colors

Applied forces and torques are orange by default. Their equal-and-opposite
reactions are gold. The load colors are distinct from the red, green, and blue
marker axes while identifying the two sides of an action/reaction pair. A
model can override either color through `graphics.loads`. Ground reactions are
hidden by default and can be enabled with **Ground loads**.
