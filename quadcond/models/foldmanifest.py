"""Persist which sequence sat in which cross-validation fold.

Folds in this project are built at training time and thrown away.  That is fine
until somebody asks the question a leakage audit has to ask -- *was this
sequence held out when that number was computed?* -- at which point there is no
artifact that answers it, and the honest reply is "cannot be established".
``benchmarks/ingest_audit/gse296171_overlap.py`` fails closed on exactly that.

This module writes the missing artifact.  It does **not** change how folds are
built: it calls :func:`quadcond.models.train.build_folds`, the same function
:func:`~quadcond.models.train.train_head` calls, on the same rows.  Because that
function is deterministic in ``(rows, spec, cluster_threshold, n_folds, seed)``,
the assignment can be reconstructed for a head that was trained earlier without
retraining it -- provided the atlas has not changed underneath, which is what
the recorded ``rows_fingerprint`` is for.

    python -m quadcond.models.foldmanifest --atlas data/atlas.db

Writes the manifest into ``benchmarks/published_tools/results/``, which is
tracked -- ``data/`` is gitignored wholesale, and a manifest nobody else can
check is not an audit artifact.

The guarantee this artifact carries is narrow and worth stating precisely: it
records the split that *this code on these rows* produces.  If the atlas has
changed since a head was trained, the fingerprint will differ and the manifest
says so rather than pretending.  A manifest whose fingerprint does not match the
shipped head must not be used to clear that head of leakage.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

MANIFEST_NAME = "fold_assignments.csv"
META_NAME = "fold_assignments_meta.json"
MANIFEST_VERSION = 2

# Tracked, unlike data/, which .gitignore excludes wholesale. A manifest that
# cannot be committed cannot be checked by anyone else, which defeats it.
DEFAULT_OUT = Path("benchmarks") / "published_tools" / "results"

#: Functions whose source determines the split. If any of them changes, a
#: reconstruction performed today is not the split that ran at training time,
#: however well the row fingerprint matches.
SPLIT_CODE_UNITS = (
    ("quadcond.models.train", "build_folds"),
    ("quadcond.models.train", "_condition_groups"),
    ("quadcond.models.grouping", "cluster_sequences"),
    ("quadcond.models.grouping", "_greedy_clusters"),
    ("quadcond.models.grouping", "kmer_matrix"),
)

PROVENANCE_RECONSTRUCTED = "reconstructed"
PROVENANCE_HISTORICAL = "historical"
PROVENANCE_STALE = "stale"

#: What each state licenses. Only one of them clears a published number.
#:
#: historical     -- the training run wrote these assignments as it trained.
#:                   Proof. Clears a published number of leakage.
#: reconstructed  -- replayed by build_folds today, over rows that match the
#:                   shipped model and with the split code unchanged since.
#:                   Evidence, not proof: the seed and fold count are not
#:                   recorded anywhere either, and a match on rows plus code
#:                   does not establish them.
#: stale          -- the rows no longer match the shipped model, or the split
#:                   code has changed since the manifest was written. Describes
#:                   a partition that is not the one any shipped head used.
#:                   Licenses nothing.
#:
#: The distinction between the last two is the one that bites: a stale manifest
#: still looks complete, still names every head, and still carries a fingerprint
#: that matches the atlas. Only the code hash tells them apart.
PROVENANCE_STATES = (PROVENANCE_HISTORICAL, PROVENANCE_RECONSTRUCTED, PROVENANCE_STALE)


def classify_manifest(meta: dict, *, shipped_fingerprints=(),
                      current_code_fingerprint: str | None = None) -> dict:
    """Decide which of the three states a manifest is actually in.

    ``assignment_provenance`` in the file records what the writer claimed at
    write time. Staleness cannot be recorded then -- it is a fact about the
    world afterwards -- so it is derived here, and it overrides the claim.
    """
    claimed = meta.get("assignment_provenance")
    rows_fp = meta.get("rows_fingerprint")
    code_fp = (meta.get("split_code_fingerprint") or {}).get("combined")
    if current_code_fingerprint is None:
        current_code_fingerprint = split_code_fingerprint()["combined"]

    reasons: list[str] = []
    shipped = list(shipped_fingerprints)

    rows_ok = bool(rows_fp) and (not shipped or rows_fp in shipped)
    if not rows_fp:
        reasons.append("manifest records no rows_fingerprint")
    elif shipped and rows_fp not in shipped:
        reasons.append(
            f"rows_fingerprint {rows_fp[:16]}... matches no shipped model")

    code_ok = bool(code_fp) and code_fp == current_code_fingerprint
    if not code_fp:
        reasons.append(
            "manifest records no split_code_fingerprint, so it cannot be shown "
            "that the split code is unchanged since it was written")
    elif code_fp != current_code_fingerprint:
        reasons.append(
            f"split code has changed since the manifest was written "
            f"({code_fp[:16]}... -> {current_code_fingerprint[:16]}...)")

    if not rows_ok or not code_ok:
        state = PROVENANCE_STALE
    elif claimed == PROVENANCE_HISTORICAL:
        state = PROVENANCE_HISTORICAL
    else:
        state = PROVENANCE_RECONSTRUCTED
        reasons.append(
            "replayed rather than recorded at training time; the seed and fold "
            "count are not recorded on any shipped head, so matching rows and "
            "code is evidence of fold identity, not proof")

    return {
        "state": state,
        "claimed_provenance": claimed,
        "rows_fingerprint_matches_a_shipped_model": rows_ok,
        "split_code_unchanged": code_ok,
        "usable_to_clear_a_published_number": state == PROVENANCE_HISTORICAL,
        "reasons": reasons,
    }


def split_code_fingerprint() -> dict:
    """Hash the source of every function that decides the split.

    Matching ``rows_fingerprint`` establishes that the *data* is the same. It
    says nothing about whether the split would come out the same, because the
    grouping threshold, the clustering implementation and the fold assignment
    all live in code that can change independently of the atlas. Without this,
    a reconstruction could match the data fingerprint exactly and still
    describe folds that never existed.
    """
    import hashlib
    import importlib
    import inspect

    h = hashlib.sha256()
    per_unit = {}
    for mod_name, fn_name in SPLIT_CODE_UNITS:
        try:
            mod = importlib.import_module(mod_name)
            src = inspect.getsource(getattr(mod, fn_name))
        except (ImportError, AttributeError, OSError) as e:
            src = f"<unavailable: {e}>"
        d = hashlib.sha256(src.encode()).hexdigest()
        per_unit[f"{mod_name}.{fn_name}"] = d[:16]
        h.update(f"{mod_name}.{fn_name}={d}".encode())
    return {"combined": h.hexdigest(), "per_unit": per_unit}

FIELDS = [
    "head", "record_id", "sequence", "nucleic_acid", "group_id", "fold",
    "kind", "label_class", "evidence_tier", "source", "source_id",
]


def _get(row, key):
    """Read a column that may not exist on this row object."""
    try:
        return row[key]
    except (IndexError, KeyError):
        return None


def assignments_for_spec(rows, spec, *, cluster_threshold=0.90, n_folds=5, seed=0):
    """Fold assignment for one head, or None when the head would not train.

    Mirrors ``train_head``'s own early exit, so a head that is skipped during
    training does not acquire a phantom manifest entry.
    """
    from .train import _rows_to_arrays, build_folds

    seqs, conds, X, y, meta = _rows_to_arrays(rows, spec, True)
    if len(y) < spec.min_rows:
        return None
    groups, split, folds = build_folds(
        seqs, conds, X, y, spec,
        cluster_threshold=cluster_threshold, n_folds=n_folds, seed=seed)

    fold_of = np.full(len(y), -1, dtype=int)
    for f, (_tr, te) in enumerate(split):
        fold_of[te] = f
    # Every row must land in exactly one test fold. A -1 here would mean a row
    # that was trained on but never evaluated, which the manifest must not
    # silently record as held out.
    unassigned = int((fold_of < 0).sum())

    # _rows_to_arrays DROPS rows whose label is null or out of class, so the
    # i-th training row is not rows[i]. meta carries the atlas primary key, and
    # that is the only safe way back to the original row. Joining positionally
    # here would mislabel provenance on every head that has a null label
    # anywhere in its query, silently and plausibly.
    by_id = {_get(r, "record_id"): r for r in rows}

    out = []
    for i in range(len(y)):
        m = meta[i] if i < len(meta) else {}
        rid = m.get("id")
        src_row = by_id.get(rid)
        out.append({
            "head": spec.name,
            "record_id": rid,
            "sequence": seqs[i],
            "nucleic_acid": (_get(src_row, "nucleic_acid") or "DNA") if src_row is not None else "DNA",
            "group_id": int(groups[i]),
            "fold": int(fold_of[i]),
            "kind": spec.kind,
            "label_class": _get(src_row, "label_class") if src_row is not None else None,
            "evidence_tier": m.get("tier"),
            "source": m.get("source"),
            "source_id": _get(src_row, "source_id") if src_row is not None else None,
        })
    return {
        "rows": out,
        "n_folds": int(folds),
        "n_groups": int(len(np.unique(groups))),
        "n_rows": int(len(y)),
        "unassigned_rows": unassigned,
        "group_by": spec.group_by,
        "cluster_threshold": cluster_threshold,
        "seed": seed,
    }


def build(atlas, tasks=None, *, cluster_threshold=0.90, n_folds=5, seed=0,
          allow_predicted=False):
    """Replay fold construction for every task and collect the assignments."""
    from .base import rows_fingerprint
    from .train import DEFAULT_TASKS

    specs = list(tasks if tasks is not None else DEFAULT_TASKS)
    all_rows: list[dict] = []
    # The fingerprint must be taken over the ATLAS rows the split was built
    # from, not over the assignment records derived from them -- otherwise it
    # would not notice the atlas changing underneath, which is the one thing it
    # exists to notice.
    fingerprint_rows: list = []
    per_head: dict[str, dict] = {}
    skipped: dict[str, str] = {}

    for spec in specs:
        tiers = list(spec.tiers)
        if "predicted" in tiers and not allow_predicted:
            skipped[spec.name] = "needs the predicted tier"
            continue
        rows = atlas.query(kind=spec.kind, tiers=tiers, label=spec.target,
                           sources=spec.sources, label_classes=spec.label_classes)
        if spec.max_rows and len(rows) > spec.max_rows:
            rs = np.random.RandomState(0)
            keep = rs.choice(len(rows), spec.max_rows, replace=False)
            rows = [rows[i] for i in sorted(keep)]
        if not rows:
            skipped[spec.name] = "no rows"
            continue
        res = assignments_for_spec(
            rows, spec, cluster_threshold=cluster_threshold,
            n_folds=n_folds, seed=seed)
        if res is None:
            skipped[spec.name] = f"fewer than {spec.min_rows} labelled rows"
            continue
        all_rows.extend(res.pop("rows"))
        fingerprint_rows.extend(rows)
        per_head[spec.name] = res

    meta = {
        "manifest_version": MANIFEST_VERSION,
        "rows_fingerprint": rows_fingerprint(fingerprint_rows) if fingerprint_rows else None,
        "split_code_fingerprint": split_code_fingerprint(),
        # Replayed now, not recorded then. Nothing in a shipped head records the
        # fold assignment that actually ran, so this is the strongest claim the
        # artifact can honestly make. Promoting it to "historical" requires
        # train_head itself to write the manifest at training time.
        "assignment_provenance": PROVENANCE_RECONSTRUCTED,
        "provenance_meaning": {
            PROVENANCE_RECONSTRUCTED: (
                "Folds were replayed today by build_folds on rows read today. This "
                "reproduces the split IF the atlas, the seed, the fold count, the "
                "cluster threshold and the split code are all unchanged since "
                "training. Matching rows_fingerprint establishes only the first of "
                "those. It is evidence, not proof, of historical fold identity."),
            PROVENANCE_HISTORICAL: (
                "Folds were written by the training run itself. This is proof."),
        },
        "heads": per_head,
        "skipped_heads": skipped,
        "total_assignments": len(all_rows),
        "caveat": ("Records the split this code produces on these rows. If the atlas "
                   "has changed since a head was trained, rows_fingerprint will not "
                   "match that head's and this manifest must not be used to clear it. "
                   "Even when it does match, a reconstructed manifest cannot prove "
                   "which folds produced a published number."),
    }
    return all_rows, meta


def write(atlas, out_dir=DEFAULT_OUT, **kw):
    rows, meta = build(atlas, **kw)
    d = Path(out_dir)
    d.mkdir(parents=True, exist_ok=True)
    with (d / MANIFEST_NAME).open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    (d / META_NAME).write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return d / MANIFEST_NAME, d / META_NAME, meta


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--atlas", default="data/atlas.db")
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    ap.add_argument("--n-folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--cluster-threshold", type=float, default=0.90)
    ap.add_argument("--allow-predicted", action="store_true")
    a = ap.parse_args()

    from ..atlas import Atlas

    atlas = Atlas(a.atlas)
    try:
        csv_path, meta_path, meta = write(
            atlas, out_dir=a.out, n_folds=a.n_folds, seed=a.seed,
            cluster_threshold=a.cluster_threshold,
            allow_predicted=a.allow_predicted)
    finally:
        atlas.close()

    print(json.dumps({
        "written": [str(csv_path), str(meta_path)],
        "total_assignments": meta["total_assignments"],
        "heads": {k: {"n_rows": v["n_rows"], "n_groups": v["n_groups"],
                      "n_folds": v["n_folds"], "unassigned_rows": v["unassigned_rows"]}
                  for k, v in meta["heads"].items()},
        "skipped_heads": meta["skipped_heads"],
        "rows_fingerprint": meta["rows_fingerprint"],
        "split_code_fingerprint": meta["split_code_fingerprint"]["combined"],
        "assignment_provenance": meta["assignment_provenance"],
    }, indent=2))
    bad = {k: v["unassigned_rows"] for k, v in meta["heads"].items()
           if v["unassigned_rows"]}
    if bad:
        print(f"\nERROR: rows never held out in any fold: {bad}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
