# External validation panels

Everything in `../published_tools` compares QuadCond with other tools on
collections QuadCond was trained on, using grouped out-of-fold predictions. That
is a fair internal protocol. It is not external validation, and the manuscript
says so.

A panel here is a set of measurements published *after* the model artifact was
frozen, scored by that artifact with no refitting. It is the only evidence that
answers "does this work on data nobody could have tuned it to".

    ..\..\.venv\Scripts\python.exe run_panel.py panels\<panel>.tsv

Use the checkout's own `.venv`, not the ambient shell. It is a Python 3.11
environment built from `requirements-model.txt` (numpy, scipy, scikit-learn and
xgboost pinned to the versions the artifact was saved under) with `quadcond`
installed editable. The model artifact will not load anywhere else, and a
conda `base` or a WSL shell that happens to have pandas is not the same thing.
`run_panel.py` says so rather than reporting a bare missing-module error.

Before extracting real sequences, check the harness runs:

    ..\..\.venv\Scripts\python.exe selftest.py

`run_panel.py` refuses to score a panel that has no publication date in its
provenance block, and refuses to score one whose declared model hash does not
match the artifact actually loaded. Both refusals exist because the value of a
panel is entirely in its provenance: an undated panel scored by an unidentified
model is a number with no claim attached.

## Current panels

| Panel | Source | What it tests | Status |
| --- | --- | --- | --- |
| `plasma_cfdna_2026` | bioRxiv 10.64898/2026.09.05.747829, 10 Sep 2026 | G4 folding and topology under a near-physiological mixed Na+/K+/Mg2+ buffer, on 50 individually CD-tested oligonucleotides | sequences to extract from the supplement |
| `fam230_kcl_titration` | Chem. Sci. 10.1039/D6SC04173F, 3 Sep 2026 | whether predicted stability tracks a 0-100 mM K+ titration on two cfDNA-derived sequences and their G-disrupted mutants | sequences to extract from the SI |

Both were published after `quadcond_model.joblib`
(`493640193ad70974d35582c55b4a721095bfcba50edc6827d9636d2f7b7be0c7`) was frozen,
which is what makes them usable this way and what the provenance block records.

## The rule that makes this directory worth having

Nothing in here is ever consulted while making a modelling decision. Not to pick
a threshold, not to choose between two artifacts, not to decide whether a head
is ready. The moment a panel informs a choice it stops being held out, and the
honest response is to move it into `../published_tools` and relabel it, not to
keep calling it external. There is no way to detect this from the files, so it
has to be a discipline.

## Adding a panel

Copy the provenance block from an existing panel, fill in the source, the
publication date and the model hash the panel is held out against, and add one
row per measurement with its own conditions. Conditions are per row and are
never assumed: this model takes buffer as an input, so scoring a physiological
mixed-cation panel at a defaulted 100 mM K+ measures nothing at all.

Keep proxy evidence out. Genomic enrichment intervals, antibody peaks and
fragmentomic scores are not folding measurements, and a panel that mixes them
with CD calls cannot be read as either.
