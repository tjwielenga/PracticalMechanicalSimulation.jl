# Sim2D Lua Assemblies

A Lua assembly is an ordinary function that adds a useful collection of
Sim2D bodies, markers, joints, forces, and graphics to a model. It does not
add a new kind of solver element. The Lua reader expands the construction into
the same hierarchical model tables used by TOML and the Julia API.

Lua is useful when an assembly needs calculations, optional values, loops, or
nested reusable assemblies. TOML remains the simpler input for a model that is
written out directly.

## A Lua model

The `sim2d` module supplies the model tables and all planar element functions:

```lua
local sim2d = require "sim2d"
local model, ground, marker =
    sim2d.model, sim2d.ground, sim2d.marker
local double_pendulum = require "planar.double_pendulum"

model {
    name = "lua_double_pendulum",
    dimension = "planar"
}

ground {name = "ground"}
marker {name = "ground.pin", position = {0.0, 0.0}}

local pendulum = double_pendulum {
    name = "pendulum",
    ground_pin = "ground.pin",
    pivot = {0.0, 0.0},
    first_length = 1.0,
    second_length = 0.8,
    first_mass = 1.0,
    second_mass = 0.7
}
```

The complete executable model is
[`lua-double-pendulum.lua`](../../models/planar/lua-double-pendulum.lua).
It uses the reusable module
[`double_pendulum.lua`](../../assemblies/planar/double_pendulum.lua).

Run it exactly like a TOML model:

```bash
./bin/simp2d models/planar/lua-double-pendulum.lua \
    --output results/examples/planar/lua-double-pendulum.simp --overwrite
```

SimpView can also open the Lua model directly, display its entered and
consistent initial configurations, and run the analysis.

## An assembly module

An assembly module returns one ordinary Lua function:

```lua
local sim2d = require "sim2d"
local rigid_body, marker, revolute =
    sim2d.rigid_body, sim2d.marker, sim2d.revolute

local function pendulum(p)
    local name = sim2d.required(p, "name", "pendulum")
    local body = name .. ".body"
    local pin = body .. ".pin"

    rigid_body {
        name = body,
        mass = sim2d.positive(p, "mass", "pendulum"),
        inertia = p.inertia,
        position = p.center
    }
    marker {name = pin, position = p.pin_offset}
    revolute {
        name = name .. ".joint",
        markers = {pin, p.ground_pin}
    }

    return {body = body, joint = name .. ".joint"}
end

return pendulum
```

The instance `name` becomes the hierarchy prefix. An outer assembly may call
this function more than once and may itself be called by another assembly.
The returned table can expose body, marker, joint, or coordinate names needed
by the caller.

`require` searches beside the model and in the repository `assemblies`
directory. A module stored as `assemblies/planar/linkage.lua` is therefore
loaded as `require "planar.linkage"`. Lua caches the module, but its returned
function may be called for any number of instances.

## Supplied helpers

The `sim2d` module provides:

- the top-level `model`, `analysis`, `simulation`, `parameters`, `graphics`,
  `state_selection`, and `initial_conditions` functions;
- one function for every Sim2D element type;
- `required`, `positive`, and `nonnegative` parameter checks;
- two-dimensional `vector` arithmetic, `dot`, `norm`, and `unit`; and
- `local_point(body, global_point)` and
  `local_orientation(body, global_angle)` for geometry specified in the
  initial global frame.

The local-geometry helpers are resolved after the complete Lua model has been
built, so the referenced body may appear before or after the assembly call.
Angles use the same numeric-radian or quoted-degree syntax accepted by TOML.

Lua model code runs through the embedded Lua runtime. It does not receive
Julia, process, operating-system, or unrestricted file access. It constructs
portable model data; it does not become part of the numerical solver.

## Stored results

Before analysis, a Lua model is expanded into ordinary primitive model tables.
That expanded TOML representation is embedded in the `.simp` result. The
result can therefore be viewed, restarted, or extracted without the original
Lua modules. Extracting the model recovers the expanded TOML, not the original
Lua source.
