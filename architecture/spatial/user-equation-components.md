# User-defined equation components

Sim3D uses the dimension-independent equation-component implementation
described in the
[common technical chapter](../common/user-equation-components.md). Its states
and algebraic variables are allocated in the spatial canonical system and are
solved with the body, joint, force, and tire equations. The spatial TOML syntax
and controlled-pendulum example are documented in the
[Spatial TOML User's Guide](../../docs/spatial/toml-reference.md#user-defined-equation-component).

Spatial loads may use qualified equation-component variables in their
expressions. In particular, applied forces, applied torques, spanning forces,
and rolling-tire trial-force expressions can turn an auxiliary equation model
into mechanical forces and torques. Saved-result initialization, statics,
dynamics, modal analysis, and `.simp` storage otherwise follow the common
behavior.
