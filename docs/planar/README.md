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

The planar feature boundary, current documentation work, deferred ideas, and
the next 3D phase are recorded in the [development roadmap](../roadmap.md).

The [`paper` workspace](../../paper/README.md) preserves the planar-methods
consolidation, maps prospective claims to executable evidence, and provides a
section-by-section scaffold for the Fully Consistent methods paper.

The maintained TOML and Lua models are in
[`models/planar`](../../models/planar/). Their companion
explanations include the [constant-speed slider-crank](../../examples/planar/constant-speed-slider-crank.md),
[torque-driven four-bar](../../examples/planar/torque-driven-planar-four-bar.md), and
[planar bushing](../../examples/planar/planar-bushing.md). The
[relative-coordinate pendulum](../../examples/planar/relative-coordinate-pendulum.md)
demonstrates a revolute-owned angle used as the preferred state, while the
[bushing-supported pendulum](../../examples/planar/bushing-pendulum.md) replaces that
ideal joint with a compliant translational support.
The [constant-speed external gear pair](../../examples/planar/constant-speed-gear-pair.md)
introduces generated floating markers and a scalar ideal gear-pair constraint.
The [internal gear pair](../../examples/planar/constant-speed-internal-gear-pair.md)
uses the same element with its contact point outside the two bearing centers.
The [planetary gear set](../../examples/planar/constant-speed-planetary-gear-set.md)
uses two gear pairs whose contact geometry follows a moving carrier even though
the sun and ring bearings are supported by ground.
The [constant-speed rack and pinion](../../examples/planar/constant-speed-rack-and-pinion.md)
uses generated floating contact markers to couple a translational joint to a
revolute joint.
The [three-pulley belt](../../examples/planar/three-pulley-belt-tensioner.md) uses
rotation-dependent elastic tangent spans and a spring-loaded moving tensioner.
The [perp-guided slider](../../examples/planar/perp-guided-slider.md) combines the
`inplane` and `perp` primitives to make a prismatic guide.
The [translational-joint slider](../../examples/planar/translational-joint-slider.md)
uses the corresponding user-facing one-freedom joint.
The [translational-distance slider](../../examples/planar/translational-distance-slider.md)
adds a one-direction prescribed distance, velocity, acceleration, and
reaction force.
The [`rotational-coordinate-coupler.toml`](../../models/planar/rotational-coordinate-coupler.toml)
and [`translational-coordinate-coupler.toml`](../../models/planar/translational-coordinate-coupler.toml)
models apply one general linear coupler to revolute rotation ports and measured
distance ports, respectively.
The [`screw-motion-coupler.toml`](../../models/planar/screw-motion-coupler.toml) model
mixes those port types to convert screw rotation into nut translation.
The [fixed-joint rotating assembly](../../examples/planar/fixed-joint-rotating-assembly.md)
shows a convenience joint assembled from revolute and perp primitives.
The [spanning-spring pendulum](../../examples/planar/spanning-spring-pendulum.md) uses the
fully expanded axial force element. The
[nonlinear spanning-force pendulum](../../examples/planar/nonlinear-spanning-force-pendulum.md)
uses the same geometry with a force expression based on its length and length
rate. The
[`distance-measures.toml`](../../models/planar/distance-measures.toml) model
shows reaction-free straight-line and directed-distance measurements and their
viewer symbols. The
[torsional-spring pendulum](../../examples/planar/torsional-spring-pendulum.md) uses the
corresponding rotational spring-damper and displays its explicit torque.
The [nonlinear expression pendulum](../../examples/planar/nonlinear-expression-pendulum.md)
uses a revolute angle and angular velocity in a user-written constitutive
torque law and obtains its sparse Jacobian partials with dual numbers.
The [PID-controlled pendulum](../../examples/planar/controlled-revolute-pendulum.md)
adds an integral state and algebraic control variables to the mechanism's
implicit system with a user-defined equation component.
The [bouncing ball](../../examples/planar/bouncing-ball.md) demonstrates compliant
one-sided contact and BDF history restart at force-law transitions.
The [rotating cam and follower](../../examples/planar/rotating-cam-follower.md)
uses a periodic cubic profile, an explicit contact station, and a circular
roller with a compliant normal force.
The [rotating cam and flat follower](../../examples/planar/rotating-cam-flat-follower.md)
uses the same profile with a marker-oriented plate and applies the force at
the moving tangent point.
The [rocker follower](../../examples/planar/rotating-cam-rocker-follower.md)
mounts the plate on a remote revolute so its orientation and moment arm both
change during contact.
The [planar friction examples](../../examples/planar/friction-elements.md)
demonstrate one-sided surface friction, revolute bearing friction,
translational guide friction, and bilateral inplane friction through sliding,
sticking, separation behavior, and breakaway.
The [stage-dependent force example](../../examples/planar/stage-dependent-forces.md)
uses static-only supports that are removed automatically when dynamics begins.
The [`modal-pendulum.toml`](../../models/planar/modal-pendulum.toml) model provides the
small analytical reference case for modal linear analysis.

The executable
[`julia_api_double_pendulum.jl`](../../examples/planar/julia_api_double_pendulum.jl)
demonstrates reusable nested Sim2D assembly functions and the resulting model
and SimpView hierarchy.
The executable
[`lua-double-pendulum.lua`](../../models/planar/lua-double-pendulum.lua)
constructs the same kind of hierarchy from a reusable Lua assembly module.

## Related material

The [planar Technical Manual](../../architecture/planar/README.md) explains the
implementation and equations. Shared integration, state-selection, result-file,
and modal methods are in the [common Technical Manual](../../architecture/common/README.md).
The separate [`paper` workspace](../../paper/README.md) retains the analysis
studies and evidence used for the Fully Consistent methods paper.
