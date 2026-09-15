"""What each head is allowed to claim.

One place decides this, because the alternative is what v0.4.0 shipped: the
report called two heads proxies, the model metadata called them
``grounded_in_measurements=true``, and the model card printed
"measurement-grounded: yes". All three were reading the same field and meaning
different things by it.

The field was the problem. ``grounded_in_measurements`` collapsed two questions
that have to stay apart:

**Was anything measured?** For an iMab CUT&Tag peak, yes. An antibody was
raised, a library was sequenced, and peaks were called from real signal.

**Was folding measured?** No. CUT&Tag reports antibody occupancy at a locus in
chromatin, which is a different observable from a melting transition in a
defined buffer -- and one whose interpretation is contested, since iMab
specificity and the accessible-chromatin contribution to targeted CUT&Tag
signal are both active questions in the literature. A score built on it cannot
license a statement about folding, at any confidence.

So three fields replace the one:

``has_experimental_observation``
    Something was measured somewhere upstream of this label.
``target_semantics``
    What the label actually is: ``biophysical``, ``genomic_proxy``,
    ``derived``, ``predicted`` -- or ``synthetic``, which was added after a
    branch of this project shipped a generator that registered made-up rows as
    ``evidence_tier="experimental"`` and trained into the canonical artifact
    path. Had it run, every field here would have certified a fit to a
    generated CSV as measurement-grounded. ``synthetic`` is checked first and
    poisons the other two fields, so no downstream formatting can undo it.
``biophysically_grounded``
    The label is a physical measurement of the structure this head predicts.
    Only this field licenses a folding or stability claim.

The v0.4.0 headline "7 measurement-grounded heads of 11" is replaced by

    7 biophysically anchored · 3 genomic-proxy · 2 derived/predicted auxiliary

which is the same twelve heads described so that the reader cannot mistake the
middle two for the first seven.
"""
from __future__ import annotations

BIOPHYSICAL = "biophysical"
GENOMIC_PROXY = "genomic_proxy"
DERIVED = "derived"
PREDICTED = "predicted"
#: Nothing was measured anywhere upstream of this label. See below.
SYNTHETIC = "synthetic"

SEMANTICS_LABEL = {
    BIOPHYSICAL: "biophysically anchored",
    GENOMIC_PROXY: "genomic proxy",
    DERIVED: "derived / control",
    PREDICTED: "prior-model output",
    SYNTHETIC: "SYNTHETIC — no measurement upstream",
}

# Sources whose positives are antibody occupancy at a genomic locus rather than
# a physical measurement of the structure.
_PROXY_SOURCE_PREFIXES = ("gse220882",)

# Claim strings that v0.4.0 got wrong. The originals said the heads predict
# "P(sequence forms an i-motif in cells)", which antibody occupancy cannot
# license. These replace them at metadata-correction time.
CORRECTED_CLAIMS = {
    "im_fold_genomic": (
        "GENOMIC-PROXY SCORE, NOT A FOLDING PROBABILITY. Calibrated score for "
        "separating canonical i-motif motifs found in reproducible iMab CUT&Tag "
        "peaks (live HEK293T, GSE220882) from composition-matched motifs taken "
        "from dinucleotide shuffles of the same peak windows. The positive label "
        "is antibody occupancy at a locus, not a measured folding event, and no "
        "buffer was measured -- the attached condition is a nominal guess with "
        "every field flagged imputed. Read it as 'does this look like the motifs "
        "enriched under iMab peaks', and never as P(folds). Its value over "
        "im_fold is sequence diversity: thousands of independent loci instead of "
        "variants of a handful of designed constructs."
    ),
    "g4_fold_genomic": (
        "GENOMIC-PROXY SCORE, NOT A FOLDING PROBABILITY. The BG4 counterpart of "
        "im_fold_genomic, on the same peaks and the same matched-shuffle control. "
        "It exists mainly as a control: the two heads share a pipeline, so a gap "
        "between them says something about the antibodies and the underlying "
        "motif grammar rather than about the model."
    ),
}

# What the calibration number does and does not cover, by head class. ECE is
# measured on the task the head was trained on, and for every folding head that
# task is a constructed, roughly balanced discrimination against generated
# negatives -- not a sample of any population a user will actually query.
CALIBRATION_SCOPE = {
    "binary": (
        "ECE and Brier are measured on the constructed, approximately balanced "
        "positive-versus-matched-shuffle task this head was trained on. The "
        "negatives are generated, not measured non-folders, so these are not "
        "absolute probabilities for an arbitrary genomic, transcriptomic or "
        "aptamer population. Recalibrate against the prevalence of your own "
        "candidate set before treating a number as a probability."
    ),
    "multiclass": (
        "ECE is measured over the three topology classes on sequences already "
        "known to form a G4, in K+ buffer. It says nothing about whether a "
        "sequence folds, and nothing about topology in Na+."
    ),
    "regression": (
        "Regression heads carry an OOF-residual interval, not a conformal "
        "guarantee: the half-width and its reported coverage come from the same "
        "pool of out-of-fold residuals, so the coverage figure is in-sample to "
        "that pool. Treat the width as a typical error scale."
    ),
}


def target_semantics(name: str, training_meta: dict) -> str:
    """Classify what a head's label actually is.

    Derived from the head's own recorded sources and tiers rather than a
    hand-kept list, so a new head cannot be added without being classified.
    """
    sources = tuple(training_meta.get("sources") or ())
    tiers = tuple(training_meta.get("tiers") or ())
    label_classes = tuple(training_meta.get("label_classes") or ())

    if any(s.startswith(_PROXY_SOURCE_PREFIXES) or
           s.startswith("shuffled::" + _PROXY_SOURCE_PREFIXES[0]) for s in sources) \
       or GENOMIC_PROXY in label_classes:
        return GENOMIC_PROXY
    # Checked before everything else, because a synthetic corpus contaminates
    # any classification it is mixed into. A branch of this project shipped a
    # generator that registered its rows as `evidence_tier="experimental"` and
    # trained into the canonical artifact path; every field below would then
    # have certified a fit to a generated CSV as measurement-grounded, and the
    # model card, /info and the viewer's badges would all have agreed. The tier
    # is what stops that, and it stops it here rather than in the one script
    # that happens to know it made the data up.
    if SYNTHETIC in tiers:
        return SYNTHETIC
    if tiers == (PREDICTED,):
        return PREDICTED
    if "experimental" in tiers:
        return BIOPHYSICAL
    return DERIVED


def semantics(name: str, training_meta: dict, task: str) -> dict:
    """The three claim-basis fields plus the calibration scope, for one head."""
    ts = target_semantics(name, training_meta)
    tiers = tuple(training_meta.get("tiers") or ())
    # Synthetic contamination poisons the claim rather than diluting it. A
    # corpus tagged ("synthetic", "experimental") still lit this field up, which
    # is exactly the shape the branch generator produced -- it registered
    # generated rows as experimental, so a head trained on a mix would have
    # reported that something was measured upstream of its label. For that head
    # nobody can say which rows drove which prediction, so the answer to "was
    # anything measured upstream of this?" is not yes.
    measured = "experimental" in tiers and SYNTHETIC not in tiers
    return {
        "has_experimental_observation": measured,
        "target_semantics": ts,
        "target_semantics_label": SEMANTICS_LABEL[ts],
        "biophysically_grounded": ts == BIOPHYSICAL,
        "calibration_scope": CALIBRATION_SCOPE.get(task, ""),
    }


def headline(heads: dict | list) -> str:
    """The 7 / 2 / 2 sentence, counted from the model rather than typed."""
    counts = tally(heads)
    line = (f"{counts[BIOPHYSICAL]} biophysically anchored · "
            f"{counts[GENOMIC_PROXY]} genomic-proxy · "
            f"{counts[DERIVED] + counts[PREDICTED]} derived/predicted auxiliary")
    # Leads, rather than appearing at the end of the list, because a reader who
    # stops after the first clause has to have been told.
    if counts[SYNTHETIC]:
        return (f"SYNTHETIC MODEL — {counts[SYNTHETIC]} of {sum(counts.values())} "
                f"heads trained on generated data, not measurements. "
                f"Not for scientific use. ({line})")
    return line


def tally(heads: dict | list) -> dict[str, int]:
    """Count heads by target semantics. Accepts Head objects or metadata dicts."""
    counts = {BIOPHYSICAL: 0, GENOMIC_PROXY: 0, DERIVED: 0, PREDICTED: 0, SYNTHETIC: 0}
    items = heads.values() if isinstance(heads, dict) else heads
    for h in items:
        meta = getattr(h, "training_meta", None)
        if meta is None:
            meta = h.get("training_meta", h)
        name = getattr(h, "name", None) or (h.get("name") if isinstance(h, dict) else "")
        counts[target_semantics(name, meta)] += 1
    return counts
