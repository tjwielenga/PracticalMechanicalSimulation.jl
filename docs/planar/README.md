# Planar Modeler User's Guide

This directory contains the user documentation for the supported two-dimensional
modeling program. It does not describe a spatial modeler.

## Using the planar modeling program

1. Start with the [project overview](../../README.md) and run one supplied model.
2. Follow [Using the Planar Modeler](using-planar-modeler.md) for a complete
   run, view, export, and restart workflow.
3. Use the [TOML User's Guide](toml-reference.md) while building or editing a
   mechanism.
4. Read the [result-file description](../../architecture/common/simulation-result-files.md)
   for HDF5 storage, viewing, CSV conversion, and model extraction.
5. Use [Sim2D Lua Assemblies](modeling-assemblies.md) to package calculated or
   hierarchical constructions for reuse.
6. Use the [Sim2D Julia model-building API](julia-api.md) to construct models
   directly in Julia, and the [Julia library API](library-api.md) for loading,
   running, and processing existing models.
7. Read [Sparse modal linear analysis](../../architecture/common/modal-linear-analysis.md)
   for operating-point linearization, the descriptor eigenproblem, and modal
   result storage.

The maintained TOML and Lua models are in
[`models/planar`](../../models/planar/). Their companion
explanations include:
- [constant-speed slider-crank](../../examples/planar/constant-speed-slider-crank.md), a kinematic slider crank mechanism;
- [torque-driven four-bar](../../examples/planar/torque-driven-planar-four-bar.md), a four bar mechanism driven by a torque on the driver link;
- [planar bushing](../../examples/planar/planar-bushing.md), an example of using a bushing;
- [relative-coordinate pendulum](../../examples/planar/relative-coordinate-pendulum.md)
demonstrates making a revolute-owned angle the preferred state;
- [bushing-supported pendulum](../../examples/planar/bushing-pendulum.md) replaces a revolute
joint with a bushing;
- [constant-speed external gear pair](../../examples/planar/constant-speed-gear-pair.md)
introduces generated floating markers and an ideal gear-pair constraint.
- [internal gear pair](../../examples/planar/constant-speed-internal-gear-pair.md)
demonstrates an internal gear set;
- [planetary gear set](../../examples/planar/constant-speed-planetary-gear-set.md) demonstrates
how a moving carrier can produce a planetary gear set;
- [constant-speed rack and pinion](../../examples/planar/constant-speed-rack-and-pinion.md)
couples a translational joint to a revolute joint;
- [three-pulley belt](../../examples/planar/three-pulley-belt-tensioner.md) models a three pulley
system with elastic tangent spans and a spring-loaded tensioner;
- [perp-guided slider](../../examples/planar/perp-guided-slider.md) combines the
`inplane` and `perp` primitives to make a prismatic guide;
- [translational-joint slider](../../examples/planar/translational-joint-slider.md)
uses the corresponding user-facing one-freedom joint;
- [translational-distance slider](../../examples/planar/translational-distance-slider.md)
adds a one-direction prescribed distance, velocity, acceleration, and
reaction force;
- [`rotational-coordinate-coupler.toml`](../../models/planar/rotational-coordinate-coupler.toml) 
couples two revolute joints together;
- [`translational-coordinate-coupler.toml`](../../models/planar/translational-coordinate-coupler.toml)
couples measured distance ports;
- [`screw-motion-coupler.toml`](../../models/planar/screw-motion-coupler.toml) couples a rotation to a translation;
- [fixed-joint rotating assembly](../../examples/planar/fixed-joint-rotating-assembly.md)
shows a convenience joint assembled from revolute and perp primitives;
- [spanning-spring pendulum](../../examples/planar/spanning-spring-pendulum.md) uses the
spanning force element;
- [nonlinear spanning-force pendulum](../../examples/planar/nonlinear-spanning-force-pendulum.md)
uses a spanning force with an expression based on its length and length rate;
- [`distance-measures.toml`](../../models/planar/distance-measures.toml) model
shows reaction-free straight-line and directed-distance measurements and their
viewer symbols;
- [torsional-spring pendulum](../../examples/planar/torsional-spring-pendulum.md) uses a
 rotational spring-damper and displays its explicit torque;
- [nonlinear expression pendulum](../../examples/planar/nonlinear-expression-pendulum.md)
uses a revolute angle and angular velocity in a user-written constitutive
torque law;
- [PID-controlled pendulum](../../examples/planar/controlled-revolute-pendulum.md)
adds an integral state and algebraic control variables to the mechanism's
implicit system with a user-defined equation component.
- [bouncing ball](../../examples/planar/bouncing-ball.md) demonstrates compliant
one-sided contact and BDF history restart at force-law transitions;
- [rotating cam and follower](../../examples/planar/rotating-cam-follower.md)
uses a cubic profile cam and a circular roller with a damped normal force;
- [rotating cam and flat follower](../../examples/planar/rotating-cam-flat-follower.md)
uses a cubic profile cam a flat plate with moving contact point;
- [rocker follower](../../examples/planar/rotating-cam-rocker-follower.md)
cam with rocking flat plate follower;
- [rocker roller follower](../../examples/planar/rotating-cam-rocker-roller.md)
cam with rocking roller follower;
- [planar friction examples](../../examples/planar/friction-elements.md)
demonstrate one-sided surface friction, revolute bearing friction,
translational guide friction, and bilateral inplane friction through sliding,
sticking, separation behavior, and breakaway;
- [stage-dependent force example](../../examples/planar/stage-dependent-forces.md)
uses static-only supports that are removed automatically when dynamics begins;
- [`modal-pendulum.toml`](../../models/planar/modal-pendulum.toml) model provides the
small analytical reference case for modal linear analysis.

The executables:
- [`julia_api_double_pendulum.jl`](../../examples/planar/julia_api_double_pendulum.jl)
demonstrates reusable nested Sim2D assembly functions and the resulting model
and SimpView hierarchy;
- [`lua-double-pendulum.lua`](../../models/planar/lua-double-pendulum.lua)
constructs the same kind of hierarchy from a reusable Lua assembly module.

## Related material

The [planar Technical Manual](../../architecture/planar/README.md) explains the
implementation and equations. Shared integration, state-selection, result-file,
and modal methods are in the [common Technical Manual](../../architecture/common/README.md).
The separate [`paper` workspace](../../paper/README.md) retains the analysis
studies and evidence used for the Fully Consistent methods paper.
