"""The frozen prediction record schema.

A prediction that travels -- into a figure, a spreadsheet, a viewer, someone
else's pipeline -- arrives stripped of the conversation that produced it. The
three claim-basis fields and the applicability block are what stop a CUT&Tag
proxy score from being read as P(folds) six months later, so the shape carrying
them is a contract rather than an implementation detail.

``PREDICTION_SCHEMA_VERSION`` is the version of the *document*, not of the model
or the atlas. Within major version 1 fields may be added and none may be removed,
renamed, or redefined; ``model_version`` and ``atlas_fingerprint`` inside each
record identify what produced it.
"""
from __future__ import annotations

import json
from pathlib import Path

PREDICTION_SCHEMA_VERSION = "1.1.0"
_SCHEMA_PATH = Path(__file__).with_name("prediction-v1.json")

__all__ = ["PREDICTION_SCHEMA_VERSION", "prediction_schema", "validate_prediction"]


def prediction_schema() -> dict:
    return json.loads(_SCHEMA_PATH.read_text())


def validate_prediction(record: dict, *, strict: bool = True) -> list[str]:
    """Return a list of schema violations; empty means valid.

    Falls back to a required-field walk when ``jsonschema`` is not installed, so
    the contract is still checkable in a minimal environment -- with less
    coverage, which the caller can detect from ``strict``.
    """
    try:
        import jsonschema
    except ImportError:
        if strict:
            raise
        return _shallow_check(record)
    v = jsonschema.Draft202012Validator(prediction_schema())
    return [f"{'/'.join(str(p) for p in e.path) or '<root>'}: {e.message}"
            for e in v.iter_errors(record)]


def _shallow_check(record: dict) -> list[str]:
    schema = prediction_schema()
    problems = [f"<root>: missing {k}" for k in schema["required"] if k not in record]
    req = schema["$defs"]["head_output"]["required"]
    for name, out in (record.get("predictions") or {}).items():
        problems += [f"predictions/{name}: missing {k}" for k in req if k not in out]
    return problems
