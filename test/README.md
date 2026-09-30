# Test Suites

The default development check is a short smoke suite. It loads and runs one
small planar model and one small spatial model, then writes and reads a `.simp`
result:

```bash
julia --project=. test/runtests.jl
```

This is the normal check for a small change. It is deliberately kept short.

For element or solver development, run the affected focused suite:

```bash
julia --project=. test/runtests.jl spatial
julia --project=. test/runtests.jl planar
julia --project=. test/runtests.jl ddassl
```

More than one name may be supplied when a change crosses boundaries. `all`
and `core` run the complete supported-program suite, including the smoke test.
The complete suite covers representative kinematic, dynamic, static, bushing,
state-selection, result-file, CSV, model-extraction, command-line, and DDASSL
behavior:

```bash
julia --project=. test/runtests.jl all
```

Run it at major checkpoints or before a release rather than after every small
change.

The paper environment preserves the alternative formulations and the
comparison packages used by the written treatment. The current publication
tables are regenerated through the evidence runner:

```bash
julia --project=test/paper -e 'using Pkg; Pkg.instantiate()'
julia --project=. paper/evidence/regenerate.jl
```

Its environment contains DASSL, OrdinaryDiffEq, and Sundials. Those comparison
packages are intentionally not installed with the supported modeling program.
The separate manifest pins their versions so the formulation comparisons can
be reproduced without constraining dependencies of the registered package.

The broader historical suite can still be run with
`julia --project=test/paper test/paper/run_all.jl`. It exercises executable analysis studies in
`examples/planar/`, including the deficit formulations, stabilization methods,
reduced-coordinate comparisons, Euler parameters, force-element development,
and intermediate assembly designs. Some of those development examples use
superseded internal component APIs and are not part of the release gate. A
formulation should not be removed merely because the current TOML program
supersedes it.

The paper suite is operationally separate from the default suite, but it still
imports `src/common/HistoricalDDASSL.jl`. Publication evidence therefore names
the package release and source commit as well as retaining this environment's
manifest. The [paper evidence directory](../paper/evidence/README.md) records
the regeneration procedure and the resulting data.
