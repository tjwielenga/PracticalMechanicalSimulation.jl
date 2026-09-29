# JCND Submission Checklist

Checked September 27, 2026 against the current [Journal of Computational and
Nonlinear Dynamics author guidelines](https://asme.onlinelibrary.wiley.com/hub/journal/15551423/author-guidelines)
and [Wiley Authors submission help](https://authors.wiley.com/help/submitting-your-manuscript.html).

## Manuscript

- Submit through the Wiley Authors Dashboard.
- An editable Word `.docx` or LaTeX source is accepted.
- JCND uses Free Format Submission. Headings, citations, and references must be
  internally consistent, but they do not need to be converted to the final
  journal typography before review.
- The abstract must be self-contained, contain no citations, and be no more
  than 300 words.
- Include keywords.
- The title page must include the author name, affiliation, email address, and
  the submitting author's ORCID iD.
- The suggested order is title, abstract, main text, references, tables,
  figures, figure legends, and appendices.

## References

- References may use any consistent format during review.
- The bibliography must include enough information to identify each source.
- A research manuscript should cite more than five references, and at least
  half of the references should be journal articles.
- Zotero remains the reference library. Better BibTeX maintains
  `references.bib`, and Pandoc renders the citations and reference list into
  the Word manuscript. The Word submission should contain the rendered
  reference list; the `.bib` file remains available as a separate source file
  or for a later LaTeX submission.

## Figures and tables

- Figures may be embedded in the review manuscript or uploaded separately.
- Keep every figure as a separate, highest-resolution source file even when it
  is also embedded in the manuscript.
- Number figures with Arabic numerals in order of appearance.
- Every figure needs a legend that identifies the figure and explains its
  content, symbols, and abbreviations.
- Each table needs a title and any necessary footnotes.
- The general Wiley submission system permits up to 500 MB combined across all
  files. Follow any smaller per-file limit displayed by the JCND submission
  dashboard.

## Required statements and submission information

- Acknowledgments, if applicable.
- Funding statement.
- Conflict-of-interest disclosure.
- Data availability statement.
- Data citations in the reference list when published data are used.
- Ethics statements, if applicable.
- Disclosure describing any use of generative artificial intelligence in
  preparing the manuscript.
- Permission and attribution for any copyrighted material owned by others.

## Publication considerations

- A cover letter is optional.
- Preprints are permitted.
- Open access is optional and may require an article publication charge.
- The first 12 typeset pages have no page charge. Each additional typeset page
  is currently listed at $200.
- Color is free online, but color figures printed in color are currently listed
  at $500 per figure.

## Project preparation

- Keep `sparseFullyConsistentMethod.md` as the authoritative manuscript.
- Generate the review Word file with `make -C paper docx`.
- Keep the pendulum SVG as the editable vector master and its PNG as the Word
  rendering.
- Keep the original SimpView screenshot as a separate PNG.
- Before submission, prepare a shorter journal manuscript if the estimated
  typeset length is substantially greater than 12 pages. Detailed equation
  catalogs, derivations, and verification records can remain in appendices,
  supporting information, or the public repository as appropriate.
