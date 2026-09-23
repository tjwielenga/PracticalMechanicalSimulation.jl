# Stage-dependent planar forces

Force elements may be limited to selected analysis stages. The
[`stage-dependent-force-drop.toml`](../../models/planar/stage-dependent-force-drop.toml)
example balances gravity with an upward applied force during static
equilibrium. A static-only bushing makes that equilibrium unique without
moving the entered position. Both supports are inactive during dynamics, so
the body subsequently falls under gravity.

```toml
[static_support]
type = "applied_force"
markers = ["body.center", "ground.up"]
force = 9.81
active_during = "static"
```

Run and view it with:

```bash
./bin/simp2d models/planar/stage-dependent-force-drop.toml \
    --output results/examples/planar/stage-dependent-force-drop.simp \
    --overwrite
bin/simpView
```

The saved result contains separate **Static initialization** and **Dynamic**
choices under SimpView's **Analysis** selector. The static choice shows the
entered model, consistent initial conditions, and converged equilibrium with
the support active. The dynamic choice starts from that equilibrium with both
static-only support elements removed.

The same fields apply to planar applied forces, applied torques, spanning
forces, bushings, and plane contacts. An inactive force retains its allocated
equations and kinematic geometry but reports zero force or torque. SimpView
does not draw its load or connector symbol during an inactive stage.
