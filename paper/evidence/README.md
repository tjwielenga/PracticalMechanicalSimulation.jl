# Paper evidence

This directory contains the reproducible measurements used by
`paper/sparseFullyConsistentMethod.md`. The evidence for each release is kept
in a separate subdirectory so that later performance changes do not silently
alter the measurements quoted by the paper.

Run the complete current-release measurement set from the repository root:

```bash
julia --project=. paper/evidence/regenerate.jl
```

The runner warms each timed path, repeats timed calculations three times, and
records the median-time run. Use `--repeats=N` to choose a different positive
repetition count. `--quick` runs one repetition and the smallest scalable
case; it checks the machinery but is not publication evidence.

The generated release directory contains:

- `metadata.toml`: release, source, Julia, operating-system, processor, and
  numerical-protocol information;
- one CSV file for every reproducible table or table section;
- `tables.md`: a compact Markdown view generated from those CSV data; and
- `scope.md`: an explicit distinction between current-release measurements
  and historical implementation comparisons.

The runner verifies that `src/` and `Project.toml` still match the named
release commit. Changes to the manuscript, evidence runner, or examples do
not invalidate that check because they do not change the simulation program.

Tables that compare obsolete implementations cannot honestly be regenerated
from one current checkout. In particular, the intermediate rows in the sparse
factorization and Large Van optimization tables remain historical measurements.
The runner measures the final 0.2.0 implementation and labels the older rows
as historical rather than presenting them as new 0.2.0 results.
