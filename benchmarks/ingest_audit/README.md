# Ingest audit

Run before a candidate dataset is allowed anywhere near training.

    python g4all_overlap.py /path/to/G4All.xlsx

Nothing in here merges, converts or retrains. It answers one question — how much
of a candidate collection is actually new — and writes a quarantine table of the
rows that are.

## Why this exists

The atlas already holds the G4STAB stability collection, the G4ShapePredictor
topology collection and the iM-Seeker transitional-pH collection, and the
published-tool benchmark evaluates QuadCond *and every competitor* on those same
collections. A curated compilation of the published literature — which is what
G4All is — will contain many of the same measurements under different source
records.

Two things go wrong if it is ingested without checking, and neither announces
itself:

**The benchmark quietly becomes partly in-sample.** Every headline number in the
manuscript is a grouped out-of-fold figure computed over those evaluation sets.
Sequences re-entering training through a differently-named file turn that into
a partially in-sample comparison, and no test in the pipeline would fail.

**Grouped folds leak.** The same melting temperature arriving under two source
DOIs is assigned to two different folds, and cross-validated performance rises
for a reason that has nothing to do with the model.

The script measures both. It prints a verdict, not just counts, because
"3,895 records" and "400 sequences the atlas did not already have" are the same
file described two ways and only one of them belongs in a paper.

## Reading the output

`fraction_new` is how the addition should be described. `benchmark_overlap` is
the number that decides whether retraining is safe: anything above a few percent
means the benchmark must either exclude those sequences from training or be
rerun on sets that exclude them, with both reported.

Compare against the full `data/atlas.db`. Against the trimmed `atlas_core.db`
the overlap is a lower bound, and the script says so on stderr.

## After the audit

The quarantine table is a starting point, not an ingest. Before anything is
merged the conditions embedded in free text need parsing into real columns,
unstable and uncertain folding calls need preserving rather than collapsing into
a binary, and every row needs its primary-paper DOI so grouping can use it.
