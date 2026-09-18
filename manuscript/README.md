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
DOI, the corresponding-author contact line, funding, CRediT roles, the conflict
declaration, and bibliography entries for the tools marked CITATION TO ADD
(G4detector, G4mismatch, G4SNVHunter). The previous draft for v0.4.9 is kept
unchanged under `release-v0.4.9/manuscript/`.
