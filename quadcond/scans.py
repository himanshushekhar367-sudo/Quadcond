"""Comparisons: one sequence against its mutants, one buffer against a range.

Both are the same shape of question -- *what changes* -- and both were being
answered badly, in ways that a single-prediction endpoint cannot fix from the
outside.

A branch of the viewer shipped a mutation explorer that called the local
heuristic energy function, ranked substitutions by its dG, and labelled
anything past a fixed magnitude "significant". Three separate problems: the
heuristic is the function that cannot tell K+ from Na+ and whose R2 against
measured folding free energy is 0.074; a difference between two structures'
energies is not a property of either; and a magnitude threshold is not a
significance test, it is a threshold with the word attached.

This module answers the difference question with the head that actually covers
the quantity, and it computes the uncertainty **on the difference** rather than
inheriting it from the two marginal predictions.

Why that distinction is the whole point
---------------------------------------

Every regression head is an ensemble. ``head.predict`` averages the members and
reports the spread as ``ensemble_std``; the interval beside a value is a
quantile of the out-of-fold residual pool, a typical error scale for *one*
prediction. Neither is the uncertainty of a difference.

Wild type and mutant differ by one base. The ensemble members largely agree
about what that base does even when they disagree about the absolute melting
temperature, so the two errors are strongly correlated and mostly cancel.
Adding the marginal half-widths in quadrature would report a difference far
noisier than it is; quoting one marginal half-width would report one far
tighter. The honest statistic is computed **paired**: ask each ensemble member
for both sequences, take that member's difference, and report the spread of
those differences. That is what ``_paired_delta`` does, and it is the only
uncertainty this module attaches to a delta.

It is still an ensemble-disagreement statistic, not a validated interval on the
difference. Mutation-effect prediction has to be validated as its own task
against measured wild-type/mutant panels; the accuracy of absolute Tm says
nothing about it. Every delta this module emits carries that sentence.
"""
from __future__ import annotations

from typing import Any, Iterable, Sequence

import numpy as np

from . import capabilities, claims
from .conditions import Condition
from .features import featurize
from .motifs import clean, find_g4, find_im, revcomp
from .schema import PREDICTION_SCHEMA_VERSION

#: Substitutions only. Indels change length, and length is an applicability
#: axis for several heads, so an indel scan silently mixes "this base matters"
#: with "this sequence is now outside the training range". A separate feature.
BASES = ("A", "C", "G", "T")

#: A scan is quadratic in nothing, but it is 3n predictions, and every head runs
#: on every one. Bounded so a pasted chromosome cannot wedge the service.
MAX_SCAN_LENGTH = 120
MAX_SEQUENCE_PER_SCAN = 5000
MAX_CONDITION_POINTS = 64
#: points x responsive heads, because that product is the actual cost.
MAX_CONDITION_PREDICTIONS = 512
MAX_BATCH_SEQUENCES = 256

DELTA_CAVEAT = (
    "Spread of the paired per-estimator differences, not a validated interval "
    "on the difference. Mutation-effect prediction is a different task from "
    "absolute prediction and has not been validated here against measured "
    "wild-type/mutant panels; this head's error on absolute values does not "
    "transfer to its error on differences."
)

NO_SIGNIFICANCE_NOTE = (
    "Ranked by magnitude of predicted change. No significance test is applied "
    "and no threshold separates 'significant' from 'not': there is no null "
    "distribution here to test against. Use the paired spread and the "
    "applicability state to judge which rows are worth measuring."
)


# --------------------------------------------------------------------- helpers
def _regression_heads(model, heads: Iterable[str] | None) -> list[str]:
    wanted = set(heads) if heads else set(model.heads)
    return [n for n, h in model.heads.items()
            if n in wanted and h.task == "regression"]


def _paired_deltas(head, wt_seq: str, mut_seqs: Sequence[str],
                   cond: Condition) -> list[dict]:
    """Per-estimator differences for many mutants at once, against one wild type.

    The pairing is the point: member *k* predicts both sequences, and the
    difference of those two numbers is member *k*'s opinion about the mutation.
    The mean of those opinions is the reported effect and their spread is the
    disagreement about it -- both computed on the difference, where the shared
    error that dominates each absolute prediction has already cancelled.

    Batched because the per-pair version was the scan's whole cost. It
    featurized two sequences and ran the estimator ensemble once **per mutant
    per head** -- on top of the `Predictor.predict` call that had already
    featurized the same mutant -- so a 22-mer cost several hundred round trips
    through sklearn. Measured in a browser competing with software WebGL for the
    same CPU, that was 78 seconds; the request was still in flight when the page
    gave up on it, and the panel reset to its start state with no error, because
    the only thing that had gone wrong was time.
    """
    X = featurize([wt_seq, *mut_seqs], [cond] * (len(mut_seqs) + 1),
                  kind=head.kind, use_conditions=head.use_conditions)
    raw = np.asarray(head._raw(X), dtype=float)   # (1 + n_mutants, n_estimators)
    base = raw[0]
    out: list[dict] = []
    for row in raw[1:]:
        per_member = row - base
        out.append({
            "delta": round(float(per_member.mean()), 3),
            "delta_paired_sd": round(float(per_member.std()), 3),
            "n_estimators": int(per_member.size),
        })
    return out


COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}


def _g_richness(seq: str) -> int:
    """G count -- the tie-break the locus ingest used to canonicalise windows."""
    return seq.upper().count("G")


def _strand_for_kind(kind: str, forward: str, reverse: str) -> tuple[str, str]:
    """Which strand carries this head's motif, and its sequence.

    A G-quadruplex and an i-motif at one duplex position sit on **opposite
    strands**. `im_pht` is a C-tract model; evaluating it on a G-rich strand is
    the category error this project has spent four releases removing from the
    evidence endpoint, and a mutation scan that ignores strand would have
    reintroduced it at a new one -- a first pass here did exactly that, asking
    three i-motif heads about Tel22's G-rich strand and getting confident
    numbers back.

    Ties and empty cases resolve to the strand richer in the relevant base
    rather than to the typed strand, because the typed strand is an accident of
    which one the user happened to paste.
    """
    # The locus head is not a motif model. Its training windows were
    # canonicalised to the G-richer strand at ingest, so it must be asked in
    # that orientation -- routing it through the i-motif finder (which is what
    # an `else` on two kinds did) asked a 201-nt window model about whichever
    # strand happened to carry more C-tracts.
    if kind == "locus":
        return ("+", forward) if _g_richness(forward) >= _g_richness(reverse) else ("-", reverse)
    finder = find_g4 if kind == "G4" else find_im
    f_hits, r_hits = len(finder(forward)), len(finder(reverse))
    if f_hits or r_hits:
        return ("+", forward) if f_hits >= r_hits else ("-", reverse)
    base = "G" if kind == "G4" else "C"
    return ("+", forward) if forward.count(base) >= reverse.count(base) else ("-", reverse)


def _motif_state(wt: str, mut: str, kind: str) -> dict:
    """Did the substitution keep, break or create the motif this head is about?

    An explicit state rather than a null. "No prediction" and "the motif is
    gone" look identical in a table of blanks, and they are opposite findings:
    the second is the strongest result a scan can produce.
    """
    finder = find_g4 if kind == "G4" else find_im
    before, after = len(finder(wt)), len(finder(mut))
    if not before and not after:
        # Distinct from "retained", which a naive equality test called this and
        # which reads as reassurance. There was never a motif here for this head
        # to be about, on this strand -- the head is not applicable at all, and
        # a delta computed anyway would be a number about nothing.
        return {"state": "no_motif", "before": 0, "after": 0,
                "note": "Neither the wild type nor the mutant carries this "
                        "head's canonical motif on this strand. The head is not "
                        "applicable, before or after."}
    if before and not after:
        return {"state": "motif_lost", "before": before, "after": after,
                "note": "The substitution removes the canonical motif. The head "
                        "is not applicable to the mutant, and that is the result."}
    if after and not before:
        return {"state": "motif_gained", "before": before, "after": after,
                "note": "The substitution creates a canonical motif absent from "
                        "the wild type."}
    if after != before:
        return {"state": "motif_count_changed", "before": before, "after": after}
    return {"state": "motif_retained", "before": before, "after": after}


def _applicability_state(entry: dict | None) -> dict:
    """Flatten a prediction's applicability into one reportable state."""
    if entry is None:
        return {"state": "head_absent"}
    ap = entry.get("applicability") or {}
    if entry.get("refused") or ap.get("refused"):
        return {"state": "refused",
                "reason": entry.get("refusal_reason") or ap.get("refusal_reason", "")}
    if not ap.get("in_domain", True):
        return {"state": "out_of_domain", "warnings": ap.get("warnings", [])}
    return {"state": "in_domain", "warnings": ap.get("warnings", [])}


CONDITION_RESPONSE_NOTE = (
    "Per head, what the training rows can support for each condition knob. "
    "'varied' -- the rows span more than one value; this does not demonstrate "
    "that the fitted model learned an effect. 'fixed' -- the training rows "
    "share one value or the head uses no condition features; a condition "
    "response is not supported by this record. 'predicts' -- the knob is this head's target, "
    "not an input. 'unknown' -- the artifact records no range for it. "
    "Unknown fields cannot establish either responsiveness or invariance."
)

#: A head whose target IS a condition variable does not take it as an input.
_TARGET_CONDITION_FIELD = {"tm": "temperature", "ph_t": "ph"}


def condition_responsiveness(model, names) -> dict:
    """Which condition knobs each head can actually respond to.

    Read off the artifact's own applicability record rather than assumed from
    the head's name: `im_pht` sounds condition-aware and is not -- every one of
    its 160 training rows sits at one buffer, so each field reports n_unique 1
    and the head can only return the same number whatever the conditions say.
    Describing it alongside `im_pht_condition` as condition-responsive was the
    overclaim this function exists to stop.
    """
    from .conditions import CONDITION_FIELDS

    out: dict[str, dict[str, str]] = {}
    for name in names:
        head = getattr(model, "heads", {}).get(name)
        if head is None:
            continue
        target_field = _TARGET_CONDITION_FIELD.get(getattr(head, "target", None))
        per: dict[str, str] = {}
        for field in CONDITION_FIELDS:
            if field == target_field:
                per[field] = "predicts"
                continue
            record = (head.applicability or {}).get(field)
            n_unique = record.get("n_unique") if isinstance(record, dict) else None
            if not getattr(head, "use_conditions", True):
                # No condition features in the model at all: what the rows did
                # is irrelevant, the knob is not wired to anything.
                per[field] = "fixed"
            elif n_unique is None:
                per[field] = "unknown"
            else:
                per[field] = "varied" if n_unique > 1 else "fixed"
        out[name] = per
    return out


def _run_record(pred, *, heads=None, **extra) -> dict:
    """The identity block every comparison workflow carries.

    One function rather than three hand-written dictionaries, because the batch
    endpoint grew a model artifact hash and an atlas fingerprint that the
    mutation and condition responses never got -- so two of the three exports
    could not say *which* model produced them, which is the whole point of
    exporting provenance at all.
    """
    from . import provenance as _prov

    model = pred.model
    # The artifact hash falls back to the hash of the file that was opened.
    # An empty string here used to satisfy the field and identify nothing.
    sha = (getattr(model, "artifact_sha256", None)
           or getattr(pred, "model_sha256", None) or None)
    rec = {
        "quadcond_model_version": getattr(model, "version", None),
        "model_artifact_sha256": sha,
        "atlas_fingerprint": getattr(model, "dataset_fingerprint", None) or None,
        "prediction_schema_version": PREDICTION_SCHEMA_VERSION,
        "provenance": _prov.run_provenance(pred),
    }
    if heads:
        rec["condition_responsiveness"] = condition_responsiveness(model, heads)
        rec["condition_note"] = CONDITION_RESPONSE_NOTE
    rec.update(extra)
    return rec


def _head_provenance(name: str, head) -> dict:
    sem = claims.semantics(name, head.training_meta, head.task)
    return {
        "head": name,
        "kind": head.kind,
        "task": head.task,
        "target": head.target,
        "units": {"tm": "degC", "ph_t": "pH units"}.get(head.target, ""),
        "claim": head.training_meta.get("claim", ""),
        "target_semantics": sem["target_semantics"],
        "biophysically_grounded": sem["biophysically_grounded"],
        "n_training_rows": head.training_meta.get("n_rows"),
    }


# ------------------------------------------------------------- mutation scan
def mutation_scan(
    pred,
    sequence: str,
    condition: Condition | None = None,
    *,
    heads: Sequence[str] | None = None,
    positions: Sequence[int] | None = None,
) -> dict:
    """Every single substitution, compared with a frozen wild type.

    The wild type and the conditions are computed once, at the top, and every
    mutant is compared with that one baseline. The panel this replaces
    recomputed the displayed baseline when a mutant was selected while keeping
    the earlier scan's rows on screen, so the table and the reference it was
    supposedly relative to could disagree without saying so.
    """
    wt = clean(sequence)
    if not wt:
        raise ValueError("empty sequence")
    if len(wt) > MAX_SCAN_LENGTH:
        raise ValueError(
            f"sequence is {len(wt)} nt; the mutation scan is capped at "
            f"{MAX_SCAN_LENGTH} nt ({3 * MAX_SCAN_LENGTH} predictions per head)")
    cond = condition or Condition()

    names = _regression_heads(pred.model, heads)
    if not names:
        raise ValueError(
            "a mutation scan compares a measurable quantity; none of the "
            "requested heads is a regression head. Classification heads change "
            "class, not magnitude, and a difference of two calibrated "
            "probabilities is not a physical effect size.")

    sites = list(positions) if positions is not None else list(range(len(wt)))
    if len(sites) > MAX_SCAN_LENGTH:
        raise ValueError(f"too many mutation positions; capped at {MAX_SCAN_LENGTH}")
    if any(isinstance(i, bool) or not isinstance(i, int) for i in sites):
        raise ValueError("mutation positions must be integer indices")
    if len(set(sites)) != len(sites):
        raise ValueError("duplicate mutation positions are not allowed")
    for i in sites:
        if not 0 <= i < len(wt):
            raise ValueError(f"position {i} is outside the sequence")

    comp = revcomp(wt)

    # Each head is evaluated on the strand that carries its motif. A
    # substitution the user specifies at forward position `i` is one edit to a
    # duplex position: base `b` on the forward strand and its complement at
    # position len-1-i on the reverse. Mutating only the typed strand would ask
    # the i-motif heads about a sequence whose complement no longer matches it.
    strand_of: dict[str, tuple[str, str]] = {
        name: _strand_for_kind(pred.model.heads[name].kind, wt, comp)
        for name in names
    }

    # Enumerate the edits once, then do all the model work in batches. The
    # per-substitution version issued two `Predictor.predict` calls and a
    # featurize-plus-ensemble pass per head for every single mutant; this issues
    # one predict call per strand and one ensemble pass per head, whatever the
    # length of the sequence.
    edits: list[dict] = []
    for i in sites:
        for base in BASES:
            if base == wt[i]:
                continue
            j = len(wt) - 1 - i
            edits.append({
                "position": i,
                "position_1based": i + 1,
                "wild_type_base": wt[i],
                "mutant_base": base,
                "label": f"{wt[i]}{i + 1}{base}",
                "sequence": wt[:i] + base + wt[i + 1:],
                "complement": comp[:j] + COMPLEMENT[base] + comp[j + 1:],
            })
    if not edits:
        raise ValueError("no substitutions to score")

    # Keyed by strand identity -- "+" and "-" -- and never by the sequence
    # string. A self-reverse-complementary input has ``wt == comp``, so a dict
    # keyed by the sequence collapsed the two strands into one entry: the "-"
    # list overwrote the "+" list, and every forward edit was then scored on its
    # reverse complement while the cell still reported strand "+". Those are two
    # different mutants of the same duplex position, and the label pointed at
    # the wrong one -- silently, and only for the palindromic sequences that a
    # G4/i-motif tool sees constantly.
    strand_sequence = {"+": wt, "-": comp}
    by_strand = {"+": [e["sequence"] for e in edits],
                 "-": [e["complement"] for e in edits]}
    used_strands = {strand for strand, _ in strand_of.values()}

    wt_by_strand: dict[str, dict] = {}
    mut_by_strand: dict[str, list[dict]] = {}
    for strand in used_strands:
        wt_by_strand[strand] = pred.predict(
            strand_sequence[strand], cond, heads=names,
            n_neighbours=0)[0]["predictions"]
        mut_by_strand[strand] = [
            r["predictions"] for r in
            pred.predict(by_strand[strand], cond, heads=names, n_neighbours=0)
        ]

    # One ensemble pass per head over every mutant on that head's strand.
    deltas: dict[str, list[dict]] = {}
    for name in names:
        strand, _ = strand_of[name]
        deltas[name] = _paired_deltas(
            pred.model.heads[name], strand_sequence[strand],
            by_strand[strand], cond)

    rows: list[dict] = []
    for k, edit in enumerate(edits):
        per_head: dict[str, Any] = {}
        for name in names:
            head = pred.model.heads[name]
            strand, strand_seq = strand_of[name]
            mut_seq = by_strand[strand][k]
            wt_entry = wt_by_strand[strand].get(name)
            mut_entry = mut_by_strand[strand][k].get(name)
            motif = _motif_state(strand_seq, mut_seq, head.kind)
            other_strand = "-" if strand == "+" else "+"
            other = _motif_state(strand_sequence[other_strand],
                                 by_strand[other_strand][k], head.kind)
            cell: dict[str, Any] = {
                # Which strand this head was asked about. Without it a reader
                # comparing a G4 row and an i-motif row is comparing two
                # sequences. The strand *sequences* are not repeated per cell --
                # they are the same two strings for the whole scan and sit at
                # the top level, where 330 copies of them used to be.
                "strand": strand,
                "motif": motif,
                "wild_type": _applicability_state(wt_entry),
                "mutant": _applicability_state(mut_entry),
            }
            # The routing above chose one strand from the wild type. A
            # substitution that builds this head's motif on the other one
            # cannot move that choice, so without this the scan reports
            # "no motif" for a variant the prescreen flagged as a possible
            # gain. Reported, never scored: the head was not asked about this
            # strand, and inventing a value for it would be the strand
            # confusion the routing exists to prevent.
            if other["state"] == "motif_gained":
                cell["other_strand"] = {
                    "strand": other_strand,
                    "state": "motif_gained",
                    "before": other["before"],
                    "after": other["after"],
                    "scored": False,
                    "note": ("This substitution creates this head's canonical "
                             "motif on the " + other_strand + " strand, which "
                             "this scan did not score: the head was routed to "
                             "the " + strand + " strand from the wild-type "
                             "sequence. This is a sequence-pattern candidate, "
                             "not a structural prediction. Scoring it requires "
                             "an explicit strand-selection workflow; reversing "
                             "the input alone does not guarantee that selection."),
                }
            # A delta is emitted only when both sides carry a value and the head
            # is applicable to both. A mutant whose motif is gone, or whose query
            # the head refused, produces a state and no number -- there is no key
            # to read a zero out of.
            if (wt_entry and "value" in wt_entry
                    and mut_entry and "value" in mut_entry
                    and motif["state"] not in ("motif_lost", "no_motif")):
                cell["wild_type_value"] = wt_entry["value"]
                cell["mutant_value"] = mut_entry["value"]
                cell.update(deltas[name][k])
            per_head[name] = cell
        rows.append({**edit, "heads": per_head})

    return {
        "schema_version": PREDICTION_SCHEMA_VERSION,
        "scan": "mutation",
        "model_version": getattr(pred.model, "version", None),
        "atlas_fingerprint": getattr(pred.model, "dataset_fingerprint", None) or None,
        # Same contract as /batch. Two of the three comparison workflows used to
        # ship without the model artifact hash, so an exported table could name
        # a model version but not the file that produced it.
        "run_record": _run_record(
            pred,
            heads=sorted(names),
            workflow="mutation_scan",
            heads_requested=sorted(names),
            n_substitutions=len(rows),
            condition=cond.to_dict(),
            condition_imputed_fields=list(cond.imputed),
            strand_by_head={n: strand_of[n][0] for n in names},
            delta_semantics=DELTA_CAVEAT,
        ),
        # The frozen baseline travels with the rows. A consumer that stores the
        # table cannot later pair it with a different wild type or buffer
        # without the mismatch being visible in the file.
        "wild_type": {
            "sequence": wt,
            "complement": comp,
            "length": len(wt),
            "condition": cond.label(),
            "condition_detail": cond.to_dict(),
            "condition_imputed_fields": list(cond.imputed),
            "predictions": {
                n: wt_by_strand[strand_of[n][0]].get(n) for n in names},
            "strand_by_head": {n: strand_of[n][0] for n in names},
            # The two strand sequences, once, instead of on every cell.
            "strand_sequence": {"+": wt, "-": comp},
        },
        "heads": {n: {**_head_provenance(n, pred.model.heads[n]),
                      **capabilities.head_capabilities(pred.model.heads[n])}
                  for n in names},
        "substitutions": rows,
        # Said once. It was attached to every cell, which put 330 copies of a
        # 250-character paragraph into a 377 KB response for a 22-mer -- most of
        # the payload, and no more informative for being repeated.
        "delta_note": DELTA_CAVEAT,
        "ranking_note": NO_SIGNIFICANCE_NOTE,
    }


# ------------------------------------------------------------ condition scan
def condition_scan(
    pred,
    sequence: str,
    axis: str,
    values: Sequence[float],
    condition: Condition | None = None,
    *,
    heads: Sequence[str] | None = None,
    n_neighbours: int = 3,
) -> dict:
    """One duplex position across a range of one condition axis.

    Replaces moving a slider and watching a number, which answers the question
    one point at a time and gives the reader no way to see where the supported
    region ends.

    Heads that did not learn a response to ``axis`` are not silently included
    with a flat line -- a flat line is a claim of insensitivity. They are
    reported in ``inert_heads`` with the reason, and no series is drawn.

    **Each head is asked about its own strand**, which is the correction that
    matters most here. Until v0.4.8 this function sent whichever strand the user
    typed to every responsive head, while ``evidence_payload`` and
    ``mutation_scan`` already routed by motif -- so sweeping K+ over the G-rich
    strand of a telomeric repeat had ``im_tm_condition`` refuse at every point,
    and pasting the C-rich complement of the *same duplex position* made those
    same queries in-domain. The available result changed with which strand the
    user happened to paste, which is not a property of the molecule.
    """
    seq = clean(sequence)
    if not seq:
        raise ValueError("empty sequence")
    if axis not in capabilities.CONDITION_AXES:
        raise ValueError(
            f"unknown condition axis {axis!r}; expected one of "
            f"{', '.join(capabilities.CONDITION_AXES)}")
    axis_values = [float(v) for v in values]
    if not axis_values:
        raise ValueError("no values to scan")
    if len(axis_values) > MAX_CONDITION_POINTS:
        raise ValueError(
            f"{len(axis_values)} points; capped at {MAX_CONDITION_POINTS}")
    if len(seq) > MAX_SEQUENCE_PER_SCAN:
        raise ValueError(
            f"sequence is {len(seq)} nt; a condition scan is capped at "
            f"{MAX_SEQUENCE_PER_SCAN} nt")

    comp = revcomp(seq)
    base = condition or Condition()
    wanted = set(heads) if heads else set(pred.model.heads)

    responsive, inert = [], []
    for name, head in pred.model.heads.items():
        if name not in wanted:
            continue
        state = capabilities.axis_response(head)[axis]
        (responsive if state["state"] == "varied" else inert).append((name, state))

    # The total work, not just the point count. A per-axis cap bounds one
    # dimension of a request whose cost is points x heads x strands; a public
    # deployment needs the product bounded, because that is what it pays for.
    if len(axis_values) * max(len(responsive), 1) > MAX_CONDITION_PREDICTIONS:
        raise ValueError(
            f"{len(axis_values)} points x {len(responsive)} responsive heads "
            f"exceeds the {MAX_CONDITION_PREDICTIONS}-prediction budget for one "
            f"condition scan; narrow the head list or the range")

    # Strand routing, identical in kind to the mutation scan's. `strand_of`
    # holds the head's strand and that strand's sequence; predictions are then
    # issued one call per strand rather than one per head.
    strand_of: dict[str, tuple[str, str]] = {
        name: _strand_for_kind(pred.model.heads[name].kind, seq, comp)
        for name, _ in responsive
    }
    strand_sequence = {"+": seq, "-": comp}
    used_strands = sorted({s for s, _ in strand_of.values()})
    heads_by_strand = {
        s: [n for n, _ in responsive if strand_of[n][0] == s] for s in used_strands
    }

    series: dict[str, list[dict]] = {n: [] for n, _ in responsive}
    for value in axis_values:
        cond = base.replace(**{axis: value})
        for strand in used_strands:
            here = heads_by_strand[strand]
            if not here:
                continue
            out = pred.predict(strand_sequence[strand], cond,
                               heads=here, n_neighbours=0)[0]
            for name in here:
                entry = out["predictions"].get(name)
                point: dict[str, Any] = {
                    "value_of_axis": value,
                    "strand": strand,
                    "applicability": _applicability_state(entry),
                }
                if entry and "value" in entry:
                    point["value"] = entry["value"]
                    if "interval" in entry:
                        point["interval"] = entry["interval"]
                    if "folded_fraction_at_condition" in entry:
                        point["folded_fraction"] = entry["folded_fraction_at_condition"]
                        point["transition_model"] = entry.get("transition_model")
                elif entry and entry.get("probability") is not None:
                    point["probability"] = entry["probability"]
                series[name].append(point)

    # Did the prediction actually move -- across the part of the range this head
    # supports?
    #
    # `capabilities.responds_to` reports that an axis varied in training, which
    # licenses a dependence but does not demonstrate one. A composition-matched
    # shuffle inherits its parent row's buffer, so a classification head can
    # show wide training variation in an axis its label never depended on -- and
    # a flat line drawn without comment reads as "this head says the structure
    # is insensitive to salt", which is a scientific claim nobody made.
    #
    # Only in-domain points enter the summary. Including extrapolated ones let
    # an unsupported tail of the sweep decide the headline verdict for the
    # supported part, which is the opposite of what a domain check is for.
    #
    # The yardstick is the head's own error scale: the half-width of its
    # interval, which is the typical size of its mistakes. This is a descriptive
    # comparison of a range with a marginal residual half-width -- not a test,
    # and not a statement that a difference is unresolvable.
    observed: dict[str, dict] = {}
    for name, _ in responsive:
        head_points = series[name]
        supported = [q for q in head_points
                     if q["applicability"]["state"] == "in_domain"]
        census = {
            "n_points": len(head_points),
            "n_supported": len(supported),
            "n_out_of_domain": sum(
                1 for q in head_points
                if q["applicability"]["state"] == "out_of_domain"),
            "n_refused": sum(
                1 for q in head_points if q["applicability"]["state"] == "refused"),
            "summary_basis": "in_domain_points_only",
        }
        key = "value" if any("value" in q for q in head_points) else "probability"
        vals = [q[key] for q in supported if key in q]
        if len(vals) < 2:
            states = {q["applicability"]["state"] for q in head_points}
            if states and states <= {"refused"}:
                observed[name] = {
                    **census, "quantity": key,
                    "verdict": "refused_throughout",
                    "note": "This head refused the query at every point on the "
                            "axis, so there is no response to report. The reason "
                            "is on each point's applicability.",
                }
            elif census["n_out_of_domain"] and census["n_supported"] < 2:
                observed[name] = {
                    **census, "quantity": key,
                    "verdict": "insufficient_supported_points",
                    "note": (
                        f"{census['n_supported']} of {census['n_points']} points "
                        f"are inside this head's applicability domain, which is "
                        f"too few to describe a response. The extrapolated points "
                        f"are returned and drawn as gaps; they are deliberately "
                        f"not summarised."),
                }
            else:
                observed[name] = {**census, "quantity": key,
                                  "verdict": "insufficient_points"}
            continue
        spread = max(vals) - min(vals)
        widths = [abs(q["interval"][1] - q["interval"][0]) / 2
                  for q in supported if "interval" in q]
        scale = min(widths) if widths else None
        entry = {
            **census,
            "quantity": key,
            "min": round(min(vals), 4),
            "max": round(max(vals), 4),
            "spread": round(spread, 4),
            "error_scale": round(scale, 3) if scale is not None else None,
        }
        if scale is None:
            entry["verdict"] = "unknown_scale"
            entry["note"] = (
                "This head reports no interval, so there is no error scale to "
                "judge the spread against. The numbers above say how far the "
                f"{key} moved across the supported range; whether that is more "
                "than this head can resolve is not something this scan can tell "
                "you.")
        elif spread < scale:
            entry["verdict"] = "below_error_scale"
            entry["note"] = (
                f"Across the supported range the prediction moves {spread:.3g}, "
                f"which is smaller than this head's own typical error "
                f"({scale:.3g}). That comparison is descriptive -- a range set "
                f"beside a marginal residual half-width -- and does not test "
                f"whether the difference is resolvable. It is not a measured "
                f"insensitivity to {capabilities.AXIS_LABEL[axis]}.")
        else:
            entry["verdict"] = "responds"
            entry["note"] = (
                f"Across the supported range the prediction moves {spread:.3g}, "
                f"against a typical error of {scale:.3g}. Descriptive, not a "
                f"significance test.")
        observed[name] = entry

    # Measured neighbours means measured. `Atlas.neighbours` defaults to
    # ("experimental", "derived"), and the derived tier includes distilled and
    # shuffled-control rows -- a composition-matched shuffle is not a
    # measurement of anything, and listing one under "measured neighbours" is
    # the claim-basis error this project spends its whole apparatus preventing.
    neighbours = []
    if pred.atlas is not None and n_neighbours:
        try:
            neighbours = pred.atlas.neighbours(
                seq, base, n=n_neighbours, tiers=("experimental",))
        except Exception:                                    # noqa: BLE001
            neighbours = []

    return {
        "schema_version": PREDICTION_SCHEMA_VERSION,
        "scan": "condition",
        "model_version": getattr(pred.model, "version", None),
        "atlas_fingerprint": getattr(pred.model, "dataset_fingerprint", None) or None,
        "run_record": _run_record(
            pred,
            heads=sorted(n for n, _ in responsive),
            workflow="condition_scan",
            axis=axis,
            n_points=len(axis_values),
            heads_requested=sorted(n for n, _ in responsive),
            base_condition=base.to_dict(),
            condition_imputed_fields=list(base.imputed),
            strand_by_head={n: strand_of[n][0] for n, _ in responsive},
        ),
        "sequence": seq,
        "complement": comp,
        # The two strand strings once, so a series can name its strand without
        # repeating a sequence on every point.
        "strand_sequence": strand_sequence,
        "axis": axis,
        "axis_label": capabilities.AXIS_LABEL[axis],
        # The axis values, and only the axis values. This key used to be
        # overwritten with one head's list of point dictionaries by a loop
        # variable that reused the name, so a client declaring `number[]`
        # received objects.
        "values": axis_values,
        "base_condition": base.to_dict(),
        "condition_imputed_fields": list(base.imputed),
        "series": {
            name: {
                **_head_provenance(name, pred.model.heads[name]),
                # Which strand this head was asked about, beside the numbers it
                # produced. Without it, two series in one chart can be about two
                # different molecules.
                "strand": strand_of[name][0],
                "strand_sequence": strand_of[name][1],
                "training_span": {k: state[k] for k in ("min", "max", "n_unique")
                                  if k in state},
                "observed_response": observed[name],
                "points": series[name],
            }
            for name, state in responsive
        },
        # Named, with the reason, rather than omitted. A head missing from a
        # chart reads as "not applicable"; a head listed as inert says which of
        # the two very different reasons applies.
        "inert_heads": [
            {**_head_provenance(name, pred.model.heads[name]),
             "axis_state": state["state"], "reason": state["note"]}
            for name, state in inert
        ],
        "measured_neighbours": neighbours,
        "neighbour_tiers": ["experimental"],
        "neighbour_note": (
            "Nearby rows from the experimental tier of the atlas, shown for "
            "context. Derived rows -- distillation targets and composition-"
            "matched shuffles -- are excluded, because a shuffled control is not "
            "a measurement of stability. Some of these may have been in this "
            "model's training data; a nearby training example is context, not "
            "independent validation."
        ),
        "summary_note": (
            "Response summaries and the drawn line cover in-domain points only. "
            "Out-of-domain and refused points are returned with their state so a "
            "consumer can show them explicitly, and are gaps in the default view."
        ),
    }


# -------------------------------------------------------------------- batch
#: Where a bare-line record stops being plausible. Longer than this on one line
#: with no ">" header is a pasted region, not a list of oligos.
MAX_LENGTH_PER_RECORD = 5000

#: The whole request, not one dimension of it. A public deployment pays for
#: sequences x conditions, so that product is what has to be bounded.
MAX_BATCH_PREDICTIONS = 2048


def parse_input(text: str, *, mode: str = "auto") -> dict:
    """Split pasted text into records, losing none of them.

    Three separate defects lived in the version this replaces, and all three had
    the same shape: an input the user supplied did not appear in the output, and
    nothing said so.

    * A FASTA header followed by no sequence produced **no record at all**,
      because the flush only fired when the buffer was non-empty. Three headers
      became two inputs, two results and an empty ``excluded`` list -- and the
      test covering this asserted the count of two, so the loss was pinned in
      place rather than caught.
    * Bare lines with no headers were **concatenated into one sequence**. The
      panel's own label promises "FASTA or one sequence per line"; ``AAAAAA`` and
      ``CCCCCC`` came back as ``AAAAAACCCCCC``.
    * Records were reconciled by identifier, so two records sharing a name were
      counted as one scored sequence while producing two rows.

    Records are therefore addressed by input **index**, which is unique by
    construction, and the identifier is carried alongside as a label. An empty
    record survives parsing precisely so that ``batch_predict`` can exclude it
    by name.

    ``mode`` is ``"fasta"``, ``"lines"``, or ``"auto"`` -- auto reads the text as
    FASTA if any line begins with ">", and as one sequence per line otherwise.
    """
    if mode not in ("auto", "fasta", "lines"):
        raise ValueError(f"unknown parse mode {mode!r}")
    lines = text.splitlines()
    if mode == "auto":
        mode = "fasta" if any(ln.lstrip().startswith(">") for ln in lines) else "lines"

    records: list[dict] = []

    if mode == "lines":
        # One record per non-blank line. Blank lines are separators, not
        # records: a trailing newline is not a sequence the user meant to send.
        for ln in lines:
            ln = ln.strip()
            if not ln:
                continue
            records.append({"index": len(records),
                            "id": f"seq_{len(records) + 1}",
                            "sequence": ln,
                            "id_source": "generated"})
        return {"mode": "lines", "records": records}

    ident: str | None = None
    buf: list[str] = []
    open_record = False

    def flush() -> None:
        nonlocal open_record
        if not open_record:
            return
        records.append({"index": len(records),
                        "id": ident or f"seq_{len(records) + 1}",
                        "sequence": "".join(buf),
                        "id_source": "fasta_header" if ident else "generated"})
        open_record = False

    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            flush()
            ident, buf = line[1:].strip() or None, []
            # Opened by the header, not by the first base. This is the line that
            # makes an empty record exist to be excluded.
            open_record = True
            continue
        if not line:
            continue
        if not open_record:
            # Sequence text before any header: a real record with a generated id
            # rather than silently discarded input.
            ident, buf, open_record = None, [], True
        buf.append(line)
    flush()
    return {"mode": "fasta", "records": records}


def _parse_fasta(text: str) -> list[dict]:
    """Backwards-compatible wrapper: records only, auto-detected mode."""
    return parse_input(text)["records"]


def batch_predict(
    pred,
    records: Sequence[dict] | str,
    conditions: Sequence[Condition] | None = None,
    *,
    heads: Sequence[str] | None = None,
    n_neighbours: int = 0,
    parse_mode: str = "auto",
) -> dict:
    """Many sequences across many conditions, with a run record.

    The run record is not release metadata. It is what makes a number in a
    figure traceable six months later: which model, which training corpus,
    which schema, which sequences were excluded and why. A batch export without
    it is a table of numbers whose provenance lives in someone's memory.

    Every supplied record is reconciled against exactly one of ``results`` or
    ``excluded``, by input index. The reconciliation is asserted before the
    response is built, so a future parsing change cannot quietly drop an input
    again.
    """
    parse_mode_used = "records"
    if isinstance(records, str):
        parsed = parse_input(records, mode=parse_mode)
        parse_mode_used = parsed["mode"]
        records = parsed["records"]
    records = [
        {**r, "index": i} for i, r in enumerate(records)
    ]
    if not records:
        raise ValueError("no sequences")
    if len(records) > MAX_BATCH_SEQUENCES:
        raise ValueError(f"{len(records)} sequences; capped at {MAX_BATCH_SEQUENCES}")
    conds = list(conditions) if conditions else [Condition()]
    if len(records) * len(conds) > MAX_BATCH_PREDICTIONS:
        raise ValueError(
            f"{len(records)} sequences x {len(conds)} conditions exceeds the "
            f"{MAX_BATCH_PREDICTIONS}-prediction budget for one batch")

    # Identifiers are the user's labels, not keys. Duplicates are reported so a
    # reader knows two rows share a name, and never merged.
    seen: dict[str, int] = {}
    duplicate_ids: list[str] = []
    for rec in records:
        ident = str(rec.get("id") or f"seq_{rec['index'] + 1}")
        seen[ident] = seen.get(ident, 0) + 1
        if seen[ident] == 2:
            duplicate_ids.append(ident)

    results: list[dict] = []
    excluded: list[dict] = []
    for rec in records:
        idx = int(rec["index"])
        raw = str(rec.get("sequence", ""))
        seq = clean(raw)
        ident = str(rec.get("id") or f"seq_{idx + 1}")
        if not seq:
            excluded.append({"record_index": idx, "id": ident,
                             "reason": "no nucleotides after cleaning",
                             "input": raw[:80]})
            continue
        if len(seq) > MAX_LENGTH_PER_RECORD:
            excluded.append({"record_index": idx, "id": ident, "reason":
                             f"length {len(seq)} exceeds {MAX_LENGTH_PER_RECORD}"})
            continue
        for cond in conds:
            out = pred.predict(seq, cond, heads=heads, n_neighbours=n_neighbours)[0]
            results.append({"id": ident, "record_index": idx, **out})

    scored_indices = {r["record_index"] for r in results}
    excluded_indices = {e["record_index"] for e in excluded}
    # Every input lands in exactly one place. This is the invariant the old
    # identifier-keyed counting could not express, let alone check.
    unaccounted = (
        {r["index"] for r in records} - scored_indices - excluded_indices
    )
    assert not unaccounted and not (scored_indices & excluded_indices), (
        f"batch reconciliation failed: {len(records)} in, "
        f"{len(scored_indices)} scored, {len(excluded_indices)} excluded")

    return {
        "schema_version": PREDICTION_SCHEMA_VERSION,
        "scan": "batch",
        "run_record": _run_record(
            pred,
            heads=sorted(heads) if heads else sorted(pred.model.heads),
            workflow="batch",
            input_mode=parse_mode_used,
            heads_requested=sorted(heads) if heads else "all",
            n_sequences_in=len(records),
            # Distinct *inputs* that produced rows, not distinct identifiers.
            # Two records sharing a name are two sequences.
            n_sequences_scored=len(scored_indices),
            n_sequences_excluded=len(excluded_indices),
            duplicate_identifiers=sorted(set(duplicate_ids)),
            n_conditions=len(conds),
            n_rows=len(results),
            conditions=[c.to_dict() for c in conds],
            condition_imputed_fields=[list(c.imputed) for c in conds],
            reconciled=True,
        ),
        "results": results,
        # Never silently dropped. A batch of 200 that returns 197 rows with no
        # explanation is a batch whose three missing sequences get noticed after
        # the figure is drawn, if at all.
        "excluded": excluded,
    }
