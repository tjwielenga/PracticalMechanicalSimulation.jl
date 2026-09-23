# Generated Results

This directory keeps generated analysis output separate from models, examples,
and source code.

- `examples/planar/` and `examples/spatial/` contain results and exports
  produced while running documented models. Regenerate the complete supported
  set with `tools/regenerate-example-results`.
- `benchmarks/` contains retained benchmark runs used for timing or viewing.
- `scratch/` contains short-lived experiments.

The contents of these subdirectories are ignored by Git. Only their empty
structure is retained. A result needed permanently by an automated test belongs
in `test/fixtures/` instead.

The regeneration script runs every top-level planar TOML model and every
top-level spatial TOML model. It also runs the supported Lua assembly examples,
with dependent large-van analyses ordered after the static solution. Preliminary
models under `models/spatial/vehicle-development` are intentionally excluded.
The ignored `results/examples/regeneration-report.txt` records successes and
failures from the latest batch. The normal-CG and high-CG large-van models
share the same 30 m/s steering maneuver and roof-ground contacts, making the
change in rollover response easy to compare in SimpView.
