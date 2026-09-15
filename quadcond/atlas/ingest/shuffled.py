"""Adapter: composition-matched negatives (DERIVED tier).

Generates dinucleotide-shuffled counterparts of the atlas's experimental
positives so that folding classifiers are trained against a background with
identical mono- and dinucleotide composition.  See ``quadcond.negatives`` for
why this matters.

These rows carry ``label_class="catalog"``: the negative label is an argument
from motif architecture, not a measurement of non-folding.  Any head trained on
them is a mixed-label-class head, and its model card says so.
"""
from __future__ import annotations

from ...conditions import Condition
from ...negatives import matched_negatives
from ..db import Record

SOURCE = "dinucleotide_shuffled_negatives"


def build(
    positives: list[dict],
    *,
    kind: str = "G4",
    n_per_positive: int = 1,
    reject_motif: bool = False,
    seed: int = 0,
) -> list[Record]:
    seqs = [p["sequence"] for p in positives]
    negs = matched_negatives(
        seqs, kind=kind, n_per_positive=n_per_positive,
        reject_motif=reject_motif, seed=seed,
    )
    out: list[Record] = []
    for i, neg in enumerate(negs):
        parent = positives[i // n_per_positive]
        cond = Condition.from_mapping(
            {c: parent.get(c) for c in
             ("k", "na", "li_nh4", "mg", "ph", "temperature", "crowder_pct", "strand_conc")},
            track_imputed=False,
        )
        out.append(
            Record(
                sequence=neg,
                kind=kind,
                source=SOURCE,
                evidence_tier="derived",
                label_class="catalog",
                condition=cond,
                folded=0,
                method=f"Altschul-Erikson dinucleotide shuffle of {parent.get('source_id') or 'positive'}",
                source_id=f"shuffle{i}:{parent.get('record_id')}",
                qc_flags=[
                    "derived_negative",
                    f"parent_record={parent.get('record_id')}",
                    f"reject_motif={reject_motif}",
                ],
            )
        )
    return out


def register(atlas) -> None:
    atlas.register_source(
        SOURCE,
        title="Dinucleotide-shuffled composition-matched negatives",
        evidence_tier="derived",
        notes="Generated in-house by QuadCond. Preserves exact mono- and dinucleotide "
              "counts of the paired positive; destroys four-tract register.",
    )
