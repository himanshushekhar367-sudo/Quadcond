# Manuscript, software release v0.5.1

`manuscript_source.md` is the text. `build_manuscript.py` renders it to Markdown
and .docx, filling six figures and four generated tables from frozen inputs. It
fits nothing and recomputes no performance figure; every number it inserts is
copied from a result file under `inputs/`.

    python build_manuscript.py     # -> QuadCond_AENNA_NAR_manuscript.{md,docx}
    python verify_numbers.py       # checks every quantitative claim in the prose

`verify_numbers.py` is the reason to trust the prose. It holds a list of claims
that appear in the rendered manuscript and recomputes each one from
`inputs/benchmark_results.csv`, `inputs/scanner_benchmark.csv`,
`inputs/posctrl_summary.csv`, `inputs/genomewide_summary.json` and
`inputs/genomewide_sensitivity.json`. A claim that cannot be recomputed from
those files does not belong in the manuscript; add the source or drop the claim.
Run it after every edit that touches a number.

`inputs/` holds frozen copies of the result files, so the manuscript builds from
this directory alone -- `gw_out/` is run output and is not in the repository. If
an analysis is rerun, copy its result file here deliberately and rerun both
scripts; nothing updates silently.

Still outstanding before submission: the public server URL, the archived release
DOI, the corresponding-author contact line, funding, CRediT roles and the
conflict declaration.

Three citations remain marked CITATION TO ADD in the text, each for a different
reason:

- **G4detector** -- journal and DOI verified (IEEE/ACM Trans. Comput. Biol.
  Bioinform., `10.1109/TCBB.2021.3073595`), author list not confirmed. Crossref,
  PubMed and Europe PMC all rate-limited the lookup, and IEEE returned 403.
  Pull the author list from the article page and add the record to
  `verified_references.json` in the same shape as the others.
- **Peptide-directed RNA i-motif folding** (`10.1039/D6SC01203E`) -- cited in the
  limitations for the point that a binding partner shifted folding without
  shifting Tm. Not yet verified against the publisher record.

`g4mismatch` and `g4snvhunter` are done: both records were read from the
publisher's own citation block and are in `verified_references.json`.

The previous draft for v0.4.9 is kept unchanged under
`release-v0.4.9/manuscript/`.
