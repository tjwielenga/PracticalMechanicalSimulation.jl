# Fully Consistent Methods Paper Workspace

Status: working methods manuscript

The intended publication venue is the ASME *Journal of Computational and
Nonlinear Dynamics*.

The intended access plan is normal no-fee journal publication accompanied by
a freely available preprint, open source code and supporting results, and
self-archiving of the accepted manuscript after ASME's embargo period.

This directory captures the methods, evidence, paper structure, and working
manuscript. It is deliberately a working record, not a polished manuscript.

- [Planar methods consolidation](planar-methods-consolidation.md) defines the
  present Fully Consistent method, records the development sequence, and identifies
  strengths, weaknesses, and unresolved questions.
- [Paper outline](fully-consistent-paper-outline.md) defines the intended
  argument and section order before the manuscript is expanded.
- [Claims and evidence](claims-and-evidence.md) connects prospective paper
  claims to executable tests, models, and measured benchmarks.
- [Working manuscript](sparseFullyConsistentMethod.md) is the authoritative
  draft of the methods paper and will be converted to PDF for review.
- [JCND submission checklist](jcnd-submission-checklist.md) records the current
  journal requirements and the corresponding project preparation steps.

## References

Zotero remains the reference library for the paper. The recommended connection
to the Markdown draft is the Better BibTeX extension for Zotero:

1. Create a Zotero collection for the paper.
2. Export that collection to `references.bib` in this directory using **Better
   BibLaTeX**.
3. Select **Keep updated** when exporting. Zotero will then update the file
   whenever an item in the collection changes.
4. Use the exported citation keys in the manuscript.

Pandoc citation forms used in the manuscript include:

```markdown
One reference [@wielenga1986numerical].

Several references [@gear1971; @petzold1982].

A reference with page numbers [@wielenga1986numerical, pp. 4–6].
```

The manuscript metadata names `references.bib` and the included
`american-society-of-mechanical-engineers.csl` style. Pandoc therefore formats
the citations and builds the reference list in ASME style. To make a Word
review copy, install Pandoc if necessary and run:

```sh
make -C paper docx
```

The generated Word citations are formatted text rather than live Zotero
fields. Zotero and the Markdown manuscript therefore remain the authoritative
sources until the manuscript is ready for final submission.

The editable pendulum illustration is retained as
`figures/planar-pendulum.svg`. The manuscript uses its rendered PNG counterpart
so that Word output does not require an SVG conversion utility.

The supported program and the paper-verification studies have separate test
paths. See the verification note in the consolidation document before using
the historical paper suite as publication evidence; old check counts and
intermediate commits are not publication baselines.

The author's local historical archive is not part of the public repository.
It provides private historical traceability, but the new paper should
distinguish original terminology from the present reconstruction and should
not incorporate external proprietary material.

The public user documentation contains a newly written background chapter on
[numerical stiffness](../docs/common/numerical-stiffness.md). The close
transcription of the author's 1986 conference paper remains in the ignored
local historical archive and is not part of the MIT-licensed repository.
