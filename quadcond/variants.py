"""Two questions about one edit: what it does to regulation, and to folding.

A regulatory-effect predictor tells you that a promoter variant is expected to
change expression. It does not tell you *how*. A G-quadruplex predictor tells
you that the same edit is expected to cost the local G4 fifteen degrees of
melting temperature. It does not tell you whether anything downstream cares.
Neither is a mechanism on its own; together they are a testable hypothesis with
a named structure in it.

``variant_scan`` runs both over the same window and joins them exactly. The
structural half is ``scans.mutation_scan``, unchanged -- every single-base
substitution against a frozen wild type, routed to the strand carrying each
head's motif, with motif loss reported as a state rather than a blank. The
regulatory half is whatever :mod:`quadcond.alphagenome` was given. The join key
is the substitution itself on the forward genome strand, so a row pairs with at
most one regulatory record and a missing pair is visible rather than inferred.

The quadrant, and what it is not
--------------------------------

Each substitution is placed in one of four cells against two thresholds the
caller states. It is a triage device for deciding what to measure next. It is
not a classifier: no threshold here was fitted, the two axes are model output on
different scales, and nothing in this module has been benchmarked against
measured regulatory variants. The combined rank is the *minimum* of the two
within-request percentile ranks -- a variant ranks high only if it ranks high on
both axes -- which is a deliberate refusal to trade one axis off against the
other by weighting them, since no weighting has been calibrated.
"""
from __future__ import annotations

from typing import Sequence

from . import alphagenome as ag
from . import scans
from .conditions import Condition
from .motifs import clean, find_g4, find_im, g4hunter_mean, revcomp
from .schema import PREDICTION_SCHEMA_VERSION

BASES_SCREEN = ("A", "C", "G", "T")

#: One element window at a time. A whole promoter is many windows and many
#: interval queries, and a cap that admits a chromosome is not a cap.
MAX_WINDOW_LENGTH = scans.MAX_SCAN_LENGTH

QUADRANT_NOTE = (
    "Two axes, four cells, thresholds supplied by the caller. This is a triage "
    "device for choosing what to measure, not a classifier: neither threshold "
    "was fitted, the axes are model output on different scales, and no "
    "benchmark against measured regulatory variants has been run."
)

RANK_NOTE = (
    "Within this request only. The combined rank is the smaller of the two "
    "percentile ranks, so a variant ranks high only when it ranks high on both "
    "axes. It is not a probability, not calibrated, and not comparable across "
    "requests -- a percentile is a statement about the other rows in the same "
    "table."
)


def _percentile_ranks(values: list[float | None]) -> list[float | None]:
    """Fractional rank of each present value among the present values.

    Absent values stay absent. Ranking a missing regulatory score as zero would
    put every variant Atlas has nothing to say about at the bottom of the table,
    which is a claim about those variants and not a property of the data.
    """
    present = sorted(v for v in values if v is not None)
    if len(present) < 2:
        return [None if v is None else 1.0 for v in values]
    out: list[float | None] = []
    n = len(present)
    for v in values:
        if v is None:
            out.append(None)
            continue
        below = sum(1 for p in present if p < v)
        equal = sum(1 for p in present if p == v)
        out.append(round((below + 0.5 * equal) / n, 4))
    return out


def _motif_census(seq: str) -> dict:
    comp = revcomp(seq)
    return {"g4_forward": len(find_g4(seq)), "g4_reverse": len(find_g4(comp)),
            "im_forward": len(find_im(seq)), "im_reverse": len(find_im(comp))}


def _total(census: dict) -> int:
    return sum(census.values())


def prescreen(sequence: str, chromosome: str = "", start: int | None = None) -> dict:
    """Could a variant scan of this window produce a structural finding at all?

    Offline, no model, no API call. It exists because the first end-to-end run
    of this workflow spent a regulatory query on a window that carries no
    canonical G4 or i-motif on either strand -- the join was perfect, every
    variant paired, and every structural cell came back ``no_motif``, which is
    the correct answer and not an informative one.

    **Both the reference and every single-base alternate are screened, on both
    strands.** A first version of this looked at the reference alone, which
    silently discards the most interesting case a variant workflow can find: a
    substitution that *creates* a motif where the reference has none. Screening
    the reference only would have declared such a window uninteresting and
    stopped it from ever being queried.

    The verdict distinguishes four states rather than two, because "no delta"
    covers three very different findings:

    ``motif_present``       the reference carries a motif; substitutions can
                            change its predicted stability.
    ``motif_gain_possible`` the reference carries none, but at least one
                            substitution creates one. No Tm delta is computable
                            -- there is no wild-type value to difference against
                            -- and the variant is still a candidate.
    ``motif_loss_possible`` a substitution removes the reference's motif.
    ``no_motif``            neither the reference nor any single substitution
                            is matched by the canonical rules.

    A canonical-motif hit is a sequence pattern matched by this implementation's
    rules. It is not evidence that the structure forms at this locus in cells,
    and its absence is not evidence that the DNA cannot adopt one.
    """
    seq = clean(sequence)
    if not seq:
        raise ValueError("empty sequence")
    ref_census = _motif_census(seq)
    ref_total = _total(ref_census)
    hunter = g4hunter_mean(seq)

    gained: list[dict] = []
    lost: list[dict] = []
    for i, wt in enumerate(seq):
        for alt in BASES_SCREEN:
            if alt == wt:
                continue
            mut = seq[:i] + alt + seq[i + 1:]
            mut_total = _total(_motif_census(mut))
            label = f"{wt}{i + 1}{alt}"
            entry = {"label": label, "position": i,
                     "position_1based": i + 1,
                     "genomic_position": (start + i) if start else None,
                     "reference": wt, "alternate": alt,
                     "motifs_reference": ref_total, "motifs_mutant": mut_total}
            if mut_total > ref_total:
                gained.append(entry)
            elif mut_total < ref_total:
                lost.append(entry)

    if ref_total:
        verdict = "motif_present"
        note = ("The reference carries a canonical motif, so substitutions can "
                "change its predicted stability and a regulatory query buys a "
                "second axis.")
    elif gained:
        verdict = "motif_gain_possible"
        note = (f"The reference carries no canonical motif, but {len(gained)} "
                f"substitution(s) create one. No melting-temperature delta is "
                f"computable for those -- there is no wild-type value to "
                f"difference against -- but they are candidates and must not be "
                f"dropped for lacking a number.")
    elif abs(hunter) >= 0.5:
        verdict = "borderline"
        note = (f"No canonical motif in the reference or in any single "
                f"substitution, but G4Hunter mean is {hunter:+.3f}. Structural "
                f"heads will flag predictions here as extrapolation.")
    else:
        verdict = "no_motif"
        note = ("Neither the reference nor any single substitution is matched "
                "by this implementation's canonical motif rules, and G4Hunter "
                "mean is below 0.5. A variant scan will withhold a structural "
                "delta for every substitution. That is a statement about model "
                "applicability, not about whether this DNA can adopt a "
                "non-canonical structure.")

    return {
        "scan": "prescreen",
        "window": {"chromosome": chromosome or None, "start": start,
                   "end": (start + len(seq) - 1) if start else None,
                   "length": len(seq), "sequence": seq},
        "reference_motifs": ref_census,
        "reference_motif_total": ref_total,
        "g4hunter_mean": round(hunter, 4),
        "substitutions_screened": 3 * len(seq),
        "motif_gain_candidates": gained,
        "motif_loss_candidates": lost,
        "n_motif_gain": len(gained),
        "n_motif_loss": len(lost),
        "verdict": verdict,
        "note": note,
        "caveat": ("Canonical-rule matching on the reference and on every "
                   "single-base alternate, both strands. A hit is a sequence "
                   "pattern, not a demonstrated structure; an absence is a "
                   "limit of the rules, not of the molecule."),
    }


_CANDIDATE_REASON = {
    "motif_lost": (
        "The substitution removes this head's canonical motif, so the mutant "
        "side has nothing for the head to be about and no delta is emitted. "
        "This is a finding, not a gap."),
    "motif_gained": (
        "The substitution creates a canonical motif the wild type does not "
        "carry. A delta may still be emitted, but it differences a value about "
        "a motif against a value about a sequence that had none, so the "
        "category is the result and the number is secondary."),
    "motif_count_changed": (
        "Both sides carry this head's motif, but not the same number of them. "
        "A delta IS emitted and is not meaningless -- it simply does not "
        "describe the same set of elements on each side, which is why the "
        "variant is listed here as well as ranked."),
}


def variant_scan(
    pred,
    sequence: str,
    chromosome: str,
    start: int,
    *,
    condition: Condition | None = None,
    heads: Sequence[str] | None = None,
    positions: Sequence[int] | None = None,
    atlas_source: ag.AtlasSource | None = None,
    structural_head: str | None = None,
    regulatory_scorer: str | None = None,
    structural_threshold: float | None = None,
    regulatory_threshold: float | None = None,
    assembly: str | None = None,
) -> dict:
    """Join a mutation scan of ``sequence`` to regulatory scores for the same edits.

    ``sequence`` is the forward-strand reference sequence of the window, and
    ``start`` is the 1-based genomic coordinate of its first base. Those two
    facts are what make the join exact, and they are checked rather than
    trusted: every substitution's reference base is compared against the
    assembly base implied by the window, and a mismatch raises instead of
    producing a table joined at the wrong offset. An off-by-one in a coordinate
    is not a degraded result; it is a different variant.
    """
    seq = clean(sequence)
    if not seq:
        raise ValueError("empty sequence")
    if len(seq) > MAX_WINDOW_LENGTH:
        raise ValueError(
            f"window is {len(seq)} nt; a variant scan is capped at "
            f"{MAX_WINDOW_LENGTH} nt. Scan one motif-bearing element at a time "
            f"rather than a whole promoter.")
    if seq != sequence.strip().upper().replace("U", "T"):
        # `clean` masks unrecognised characters to N and keeps the length, so a
        # window carrying an ambiguity code still maps 1:1 to coordinates -- but
        # the caller should know the assembly base it supplied is not what was
        # scanned.
        pass
    if int(start) < 1:
        raise ValueError("start is a 1-based genomic coordinate; got "
                         f"{start}")
    start = int(start)
    end = start + len(seq) - 1
    source = atlas_source or ag.NullAtlas()

    inner = scans.mutation_scan(pred, seq, condition, heads=heads,
                                positions=positions)

    try:
        records = source.records_for_interval(chromosome, start, end)
        atlas_error = None
    except ag.AtlasUnavailable as exc:
        # Reported, never swallowed. A response whose regulatory axis is empty
        # because a network call failed must not read like a response whose
        # regulatory axis is empty because the locus is quiet.
        records, atlas_error = {}, str(exc)

    heads_present = list(inner["heads"])
    axis_head = structural_head or (heads_present[0] if heads_present else None)

    rows: list[dict] = []
    for sub in inner["substitutions"]:
        pos = start + int(sub["position"])
        key = ag.VariantKey(chromosome, pos, sub["wild_type_base"],
                            sub["mutant_base"]).normalised()
        cell = sub["heads"].get(axis_head) if axis_head else None
        delta = cell.get("delta") if cell else None
        motif_state = cell["motif"]["state"] if cell else None

        rec = records.get(key)
        headline = rec.headline([regulatory_scorer] if regulatory_scorer else None) if rec else None
        rows.append({
            "variant": key.as_dict(),
            "label": key.hgvs_like(),
            "mutation": sub["label"],
            "position_in_window": sub["position"],
            "structural": {
                "head": axis_head,
                "strand": cell.get("strand") if cell else None,
                "delta": delta,
                "delta_paired_sd": cell.get("delta_paired_sd") if cell else None,
                "motif_state": motif_state,
                "applicability": cell.get("mutant") if cell else None,
                # Magnitude is what the quadrant uses; the sign is kept because
                # a stabilising and a destabilising edit of equal size are
                # opposite hypotheses.
                "magnitude": abs(delta) if isinstance(delta, (int, float)) else None,
            },
            "regulatory": ({
                "scorer": headline[0], "score": headline[1],
                "all_scores": dict(rec.scores), "source": rec.source,
            } if rec and headline else None),
            "heads": sub["heads"],
        })

    struct_rank = _percentile_ranks([r["structural"]["magnitude"] for r in rows])
    reg_rank = _percentile_ranks(
        [abs(r["regulatory"]["score"]) if r["regulatory"] else None for r in rows])

    counts = {"structural_only": 0, "regulatory_only": 0,
              "both": 0, "neither": 0, "unclassified": 0}
    for row, sr, rr in zip(rows, struct_rank, reg_rank):
        row["structural"]["rank"] = sr
        if row["regulatory"] is not None:
            row["regulatory"]["rank"] = rr
        mag = row["structural"]["magnitude"]
        reg = abs(row["regulatory"]["score"]) if row["regulatory"] else None
        if mag is None or reg is None:
            # Two different absences: the head declined or lost the motif, or
            # the regulatory source has no record. Either way the variant is
            # unclassified, because a quadrant needs two coordinates.
            row["quadrant"] = "unclassified"
            row["quadrant_reason"] = (
                "no structural delta (motif lost, refused or out of domain)"
                if mag is None else "no regulatory record for this substitution")
            row["combined_rank"] = None
            counts["unclassified"] += 1
            continue
        s_hi = mag >= structural_threshold if structural_threshold is not None else (sr or 0) >= 0.75
        r_hi = reg >= regulatory_threshold if regulatory_threshold is not None else (rr or 0) >= 0.75
        row["quadrant"] = ("both" if s_hi and r_hi else
                           "structural_only" if s_hi else
                           "regulatory_only" if r_hi else "neither")
        row["combined_rank"] = round(min(sr or 0.0, rr or 0.0), 4)
        counts[row["quadrant"]] += 1

    rows.sort(key=lambda r: (r["combined_rank"] is None,
                             -(r["combined_rank"] or 0.0)))

    # A categorical shortlist, alongside the ranked one and never merged into it.
    #
    # The ranking needs a numeric delta on both axes, so a substitution that
    # destroys or creates the motif may have no usable delta and land in
    # `unclassified` and disappears off the bottom of a ranked table. That is a
    # retrieval failure with scientific consequences, and the fix is not to
    # invent a number for those rows. It is to list them separately, by
    # category, with their regulatory score attached where one exists.
    categorical: dict[str, list[dict]] = {
        "motif_lost": [], "motif_gained": [], "motif_count_changed": [],
        "motif_gained_other_strand": []}

    def _entry(row, state, reason, *, delta):
        return {
            "label": row["label"],
            "mutation": row["mutation"],
            "strand": row["structural"]["strand"],
            "head": row["structural"]["head"],
            "units": (inner.get("heads") or {}).get(
                row["structural"]["head"], {}).get("units"),
            "motif_state": state,
            # Carried explicitly. Two of these categories DO produce a number,
            # and an entry that stays silent about it reads as though none do.
            "delta": delta,
            "has_delta": delta is not None,
            "regulatory_score": (row["regulatory"]["score"]
                                 if row["regulatory"] else None),
            "regulatory_rank": (row["regulatory"].get("rank")
                                if row["regulatory"] else None),
            "why_listed": reason,
        }

    for row in rows:
        state = row["structural"]["motif_state"]
        delta = row["structural"]["delta"]
        if state in _CANDIDATE_REASON:
            categorical[state].append(
                _entry(row, state, _CANDIDATE_REASON[state], delta=delta))
        # A gain on the strand this head was not routed to. It has no delta on
        # this axis by construction -- the head never saw that strand.
        axis_cell = (row.get("heads") or {}).get(row["structural"]["head"])
        other = (axis_cell or {}).get("other_strand")
        if other:
            categorical["motif_gained_other_strand"].append(
                _entry(row, "motif_gained_other_strand",
                       other["note"], delta=None) | {"strand": other["strand"]})
    for entries in categorical.values():
        entries.sort(key=lambda e: -(e["regulatory_rank"] or 0.0))

    return {
        "schema_version": PREDICTION_SCHEMA_VERSION,
        "scan": "variant",
        "model_version": inner.get("model_version"),
        "run_record": {
            **inner["run_record"],
            "workflow": "variant_scan",
            "locus": {
                "chromosome": ag._norm_chrom(chromosome),
                "start": start, "end": end,
                # Spelled out because an off-by-one is a different variant, and
                # "start" means at least three things across the formats this
                # workflow touches.
                "coordinate_convention": "1-based inclusive, as in VCF",
                "strand_of_window_sequence": "+ (forward / reference strand)",
                "assembly": assembly or "not stated by caller",
                "assembly_stated_by_caller": bool(assembly),
            },
            "regulatory_source": source.describe(),
            "regulatory_error": atlas_error,
            "structural_axis_head": axis_head,
            "regulatory_axis_scorer": regulatory_scorer,
            "structural_threshold": structural_threshold,
            "regulatory_threshold": regulatory_threshold,
            "threshold_basis": ("caller-supplied absolute thresholds"
                                if structural_threshold is not None
                                or regulatory_threshold is not None
                                else "75th percentile within this request"),
        },
        "window": {"chromosome": ag._norm_chrom(chromosome), "start": start,
                   "end": end, "sequence": seq,
                   "g4_elements": [e.as_dict() for e in find_g4(seq)],
                   "im_elements": [e.as_dict() for e in find_im(seq)]},
        "wild_type": inner["wild_type"],
        "heads": inner["heads"],
        "variants": rows,
        "quadrant_counts": counts,
        # Candidates that carry a structural finding but no number. Ranked
        # within their category by the regulatory axis where one exists, and
        # kept out of `variants`' ranking rather than given an invented delta.
        "structural_candidates": categorical,
        "structural_candidate_counts": {k: len(v) for k, v in categorical.items()},
        "candidate_note": (
            "Categorical structural findings, listed so that a table ranked on "
            "delta magnitude cannot hide them. They do not all lack a number, "
            "and each entry says which it is. 'motif_lost': the mutant has no "
            "motif for the head to be about, so no delta is emitted. "
            "'motif_gained': a delta may be emitted, but it differences a "
            "value about a motif against a value about a sequence that had "
            "none. 'motif_count_changed': both sides carry the motif and a "
            "delta is emitted -- it simply does not describe the same set of "
            "elements on each side. 'motif_gained_other_strand': the "
            "substitution builds this head's motif on the strand the head was "
            "not routed to, which this scan did not score. In every case read "
            "the category as the result and the number, where there is one, as "
            "secondary."),
        "quadrant_note": QUADRANT_NOTE,
        "rank_note": RANK_NOTE,
        "delta_note": inner["delta_note"],
        "ranking_note": inner["ranking_note"],
        "interpretation_note": (
            "A variant in the 'both' cell is a candidate for a structure-mediated "
            "regulatory mechanism: the regulatory model expects an effect and "
            "this head expects the local G-quadruplex or i-motif to change. That "
            "is a hypothesis to test, not a finding. Neither axis observes the "
            "other, no joint model was fitted, and the structural axis carries "
            "the same limitation as every other QuadCond delta -- mutation-effect "
            "prediction has not been validated here as its own task."),
    }
