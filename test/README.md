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

The paper verification suite preserves the alternative formulations and the
development comparisons used to support the written treatment:

```bash
julia --project=test/paper -e 'using Pkg; Pkg.instantiate()'
julia --project=test/paper test/paper/run_all.jl
```

Its environment contains DASSL, OrdinaryDiffEq, and Sundials. Those comparison
packages are intentionally not installed with the supported modeling program.

These checks intentionally exercise the executable analysis studies in
`examples/planar/`, including the deficit formulations, stabilization methods,
reduced-coordinate comparisons, Euler parameters, force-element development,
and intermediate assembly designs. A formulation should not be removed from
this suite merely because the current TOML program supersedes it.

The paper suite is operationally separate from the default suite, but it still
imports the evolving `src/common/HistoricalDDASSL.jl`; it is not a numerically frozen
historical baseline. The [planar consolidation](../paper/planar-methods-consolidation.md)
records the resulting verification issue and the choice that must be made
before publication comparisons are finalized.
