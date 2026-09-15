"""What each head responds to, and which workspace it belongs in.

Two questions the interface kept answering by guessing.

**Which controls actually do anything to this result?** The app has six ion
sliders, a pH control and a temperature control. No head consumes all of them,
and some consume none. A slider that moves while a number sits still is worse
than a missing slider: it invites the reader to conclude the head is
*insensitive* to that variable, which is a scientific claim, when the truth is
that the training data never varied it. ``applicability`` already carries the
evidence -- ``n_unique`` per axis -- and nothing was reading it for this.

**Which of these twelve things is this?** Twelve heads on one screen implies
twelve comparable numbers. They are four different kinds of statement, and the
grouping is derivable from metadata the model already carries rather than from
a list typed into a component, which is what makes it stay correct when a head
is added.

Nothing here is a new claim. It is `claims.py`'s three fields and the
applicability table, read for a question the UI had been answering on its own.
"""
from __future__ import annotations

from . import claims

#: The condition axes a query can specify, in the order the interface shows them.
CONDITION_AXES = (
    "k", "na", "li_nh4", "mg", "ph", "temperature", "crowder_pct", "strand_conc",
)

#: Human labels, so the UI does not invent its own.
AXIS_LABEL = {
    "k": "K+", "na": "Na+", "li_nh4": "Li+/NH4+", "mg": "Mg2+",
    "ph": "pH", "temperature": "temperature",
    "crowder_pct": "crowding (% PEG)", "strand_conc": "strand concentration",
}

SEQUENCE_EVIDENCE = "sequence_evidence"
CONDITION_RESPONSE = "condition_response"
GENOMIC_CONTEXT = "genomic_context"
DIAGNOSTICS = "diagnostics"

WORKSPACE_LABEL = {
    SEQUENCE_EVIDENCE: "Sequence evidence",
    CONDITION_RESPONSE: "Condition response",
    GENOMIC_CONTEXT: "Genomic context",
    DIAGNOSTICS: "Advanced diagnostics",
}

WORKSPACE_NOTE = {
    SEQUENCE_EVIDENCE:
        "Folding classification, topology, melting temperature and transition "
        "pH. What the sequence itself supports, at the reference buffer each "
        "head was trained in.",
    CONDITION_RESPONSE:
        "Heads whose training data actually spans the axis being varied - so a "
        "prediction can depend on it. A head appears here because its training "
        "table varied, not because the interface has a slider. Whether the "
        "prediction does move is measured, not assumed: each series reports the "
        "observed spread across the range against the head's own error scale, "
        "and says so when the answer is 'it does not'.",
    GENOMIC_CONTEXT:
        "Antibody occupancy at a locus in chromatin. A different observable "
        "from a melting transition in a defined buffer, kept on its own axis "
        "so it cannot be read as P(folds).",
    DIAGNOSTICS:
        "Predictions inherited from another model, and labels derived from the "
        "motif definition itself. Useful for benchmarking and for exercising "
        "the pipeline; not evidence about DNA.",
}


def axis_response(head) -> dict[str, dict]:
    """Per-axis: did this head learn a response, and on what evidence?

    Three states, and collapsing any two of them loses something:

    ``varied``     the axis moved in training over ``n_unique`` values; a
                   prediction can depend on it.
    ``fixed``      every training record used one value. The head has no
                   response to apply. This is an *absence*, not an
                   extrapolation, and it is why a slider must not imply an
                   answer.
    ``unknown``    no applicability entry. Treated as ``fixed`` by the caller,
                   because assuming a response nobody recorded is the failure
                   mode this module exists to prevent.

    ``predicts`` marks the axis a Tm head outputs rather than consumes, which is
    a fourth thing again: temperature is neither varied nor fixed for `g4_tm`,
    it is the answer.
    """
    ap = head.applicability or {}
    out: dict[str, dict] = {}
    for axis in CONDITION_AXES:
        if axis == "temperature" and getattr(head, "target", "") == "tm":
            out[axis] = {
                "state": "predicts",
                "label": AXIS_LABEL[axis],
                "note": (
                    "This head predicts a melting temperature; it does not take "
                    "one. The requested temperature is the evaluation point of "
                    "the two-state model applied afterwards, not an input here."
                ),
            }
            continue
        rng = ap.get(axis)
        if not isinstance(rng, dict) or "n_unique" not in rng:
            out[axis] = {
                "state": "unknown",
                "label": AXIS_LABEL[axis],
                "note": "No training range recorded for this axis; treated as no response.",
            }
            continue
        if rng["n_unique"] <= 1:
            out[axis] = {
                "state": "fixed",
                "label": AXIS_LABEL[axis],
                "value": rng["min"],
                "note": (
                    f"Every training record used {AXIS_LABEL[axis]}="
                    f"{rng['min']:g}. This head has no response to apply - an "
                    f"absence, not an extrapolation."
                ),
            }
            continue
        out[axis] = {
            "state": "varied",
            "label": AXIS_LABEL[axis],
            "min": rng["min"],
            "max": rng["max"],
            "n_unique": rng["n_unique"],
            "note": (
                f"Training spans {rng['min']:g}-{rng['max']:g} over "
                f"{rng['n_unique']} distinct values."
            ),
        }
    return out


def responds_to(head) -> list[str]:
    """The axes a prediction from this head can actually depend on."""
    return [a for a, r in axis_response(head).items() if r["state"] == "varied"]


def workspace(head) -> str:
    """The one group this head is filed under.

    Semantics decides, and only semantics. A genomic proxy that happens to
    carry a condition column is still a genomic proxy; a distilled head that
    tracks salt beautifully is still another model's output.

    An earlier version also moved condition-responsive heads out of sequence
    evidence, which emptied that group of `g4_tm` -- the best-attested head in
    the model, and the first thing anyone looks for. Condition response is not a
    different *kind* of head, it is a second question about the same one, so it
    is an additional placement (`also_in`) rather than a reassignment.
    """
    meta = getattr(head, "training_meta", {}) or {}
    ts = claims.target_semantics(head.name, meta)
    if ts == claims.SYNTHETIC:
        return DIAGNOSTICS
    if ts in (claims.DERIVED, claims.PREDICTED):
        return DIAGNOSTICS
    if ts == claims.GENOMIC_PROXY:
        return GENOMIC_CONTEXT
    return SEQUENCE_EVIDENCE


def also_in(head) -> list[str]:
    """Extra groups this head appears in, beyond its filed one."""
    return [CONDITION_RESPONSE] if responds_to(head) else []


def head_capabilities(head) -> dict:
    """Everything the interface needs to place and gate one head."""
    axes = axis_response(head)
    ws = workspace(head)
    return {
        "workspace": ws,
        "workspace_label": WORKSPACE_LABEL[ws],
        "also_in": also_in(head),
        "axes": axes,
        "responds_to": [a for a, r in axes.items() if r["state"] == "varied"],
        # Said here so no consumer has to infer it from the field name.
        # `responds_to` is read off the training table: these axes took more
        # than one value, so a prediction *can* depend on them. That it does is
        # not implied, and is measured empirically by /scan/conditions, which
        # reports the observed spread and says "flat" when there is none. A
        # composition-matched shuffle inherits its parent row's buffer, so a
        # classification head can show wide training variation in an axis its
        # label never depended on.
        "responds_to_note": (
            "These axes varied in this head's training data, so a prediction "
            "can depend on them. That a prediction does depend on them is not "
            "implied, and is measured by /scan/conditions."
        ),
        "inert_axes": [a for a, r in axes.items() if r["state"] in ("fixed", "unknown")],
        "predicted_axes": [a for a, r in axes.items() if r["state"] == "predicts"],
        "use_conditions": bool(getattr(head, "use_conditions", False)),
    }


def workspaces(heads: dict) -> list[dict]:
    """The four groups, populated, in reading order - empty ones included.

    An empty group is served rather than dropped: a viewer that shows three
    tabs when the backend has four kinds of head is a viewer that has quietly
    hidden a category of claim.
    """
    order = [SEQUENCE_EVIDENCE, CONDITION_RESPONSE, GENOMIC_CONTEXT, DIAGNOSTICS]
    grouped: dict[str, list[str]] = {k: [] for k in order}
    for name, head in heads.items():
        grouped[workspace(head)].append(name)
        for extra in also_in(head):
            grouped[extra].append(name)
    return [
        {
            "id": key,
            "label": WORKSPACE_LABEL[key],
            "note": WORKSPACE_NOTE[key],
            "heads": sorted(grouped[key]),
        }
        for key in order
    ]
