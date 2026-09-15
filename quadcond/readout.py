"""One place that knows how to read a prediction entry.

A refused head carries no ``value``, ``probability`` or ``posterior`` key. That
is deliberate: a consumer which forgets to check gets a ``KeyError`` rather than
a number the domain checker had already rejected sitting in a figure.

The cost of that design is that every consumer has to check, and until v0.4.8
several of them did not -- because the refusal branch in ``Predictor.predict``
ran *after* the binary and multiclass branches and so never ran for a
classification head at all. Querying ``ZZZZ`` (which cleans to ``NNNN``) made
``g4_tm`` and ``im_pht`` refuse properly while ``g4_fold`` and ``im_fold``
returned calibrated probabilities carrying ``applicability.refused: true`` and
no top-level ``refused`` flag, so the batch table -- which checks the top-level
field -- rendered classifier numbers for an input the service does not consider
a sequence.

Refusal is now decided before any task-specific serialisation, for every head
type, and these helpers exist so that the consumers all check it the same way.
"""
from __future__ import annotations

from typing import Any


def is_refused(entry: dict | None) -> bool:
    """True when this head declined the query, whatever its task."""
    if not entry:
        return False
    ap = entry.get("applicability") or {}
    return bool(entry.get("refused") or ap.get("refused"))


def refusal_reason(entry: dict | None) -> str:
    if not entry:
        return ""
    ap = entry.get("applicability") or {}
    return str(entry.get("refusal_reason") or ap.get("refusal_reason") or "")


def readout(entry: dict | None) -> tuple[str, Any] | None:
    """``(kind, value)`` for a rendered entry, or ``None`` if there is nothing.

    ``kind`` is one of ``"probability"``, ``"class"`` or ``"value"``. ``None``
    means refused, absent, or a head that returned no scalar -- three states a
    caller must not paper over with a zero.
    """
    if entry is None or is_refused(entry):
        return None
    if "probability" in entry:
        return "probability", entry["probability"]
    if "posterior" in entry:
        return "class", entry.get("argmax")
    if "value" in entry:
        return "value", entry["value"]
    return None
