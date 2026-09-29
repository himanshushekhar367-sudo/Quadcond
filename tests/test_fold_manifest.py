"""Tests for the persisted fold manifest.

The manifest's whole value is that it records the split *training actually
used*.  So the tests that matter are the ones that would catch it drifting:
that it comes from the same function `train_head` calls, that it is
deterministic, that every row lands in exactly one held-out fold, and that a
manifest built from the wrong atlas is rejected rather than trusted.
"""
from __future__ import annotations

import json

import numpy as np
import pytest

from quadcond.atlas import Atlas, Record
from quadcond.conditions import Condition
from quadcond.models import foldmanifest
from quadcond.models.train import DEFAULT_TASKS, TaskSpec, build_folds, _rows_to_arrays

G = "GGGTTAGGGTTAGGGTTAGGG"


@pytest.fixture
def atlas(tmp_path):
    a = Atlas(tmp_path / "t.db")
    a.register_source("src", evidence_tier="experimental")
    recs = []
    for i in range(60):
        seq = "GGG" + "TTA" * (1 + i % 5) + "GGG" + "A" * (i % 7) + "GGGTTAGGG"
        recs.append(Record(seq, "G4", "src", "experimental",
                           Condition.from_mapping({"k": 100.0, "ph": 7.0}),
                           tm=40.0 + i % 30, source_id=f"r{i}"))
    a.add(recs)
    yield a
    a.close()


def _tm_spec():
    for s in DEFAULT_TASKS:
        if s.name == "g4_tm":
            return s
    pytest.skip("no g4_tm task")


# --------------------------------------------------------------------------
def test_manifest_uses_the_same_split_training_uses(atlas):
    """The manifest must not reimplement fold construction.

    If it did, the two would drift the first time either changed, and a leakage
    audit reading a drifted manifest reports a clean bill from the wrong folds.
    """
    spec = _tm_spec()
    rows = atlas.query(kind=spec.kind, tiers=list(spec.tiers), label=spec.target)
    seqs, conds, X, y, meta = _rows_to_arrays(rows, spec, True)
    groups, split, folds = build_folds(seqs, conds, X, y, spec)

    res = foldmanifest.assignments_for_spec(rows, spec)
    assert res is not None

    expected = np.full(len(y), -1, dtype=int)
    for f, (_tr, te) in enumerate(split):
        expected[te] = f
    assert [r["fold"] for r in res["rows"]] == [int(v) for v in expected]
    assert [r["group_id"] for r in res["rows"]] == [int(g) for g in groups]


def test_every_row_is_held_out_exactly_once(atlas):
    res = foldmanifest.assignments_for_spec(rows_of(atlas), _tm_spec())
    assert res["unassigned_rows"] == 0
    assert all(r["fold"] >= 0 for r in res["rows"])


def test_manifest_is_deterministic(atlas):
    spec = _tm_spec()
    a = foldmanifest.assignments_for_spec(rows_of(atlas), spec)
    b = foldmanifest.assignments_for_spec(rows_of(atlas), spec)
    assert [r["fold"] for r in a["rows"]] == [r["fold"] for r in b["rows"]]
    assert [r["group_id"] for r in a["rows"]] == [r["group_id"] for r in b["rows"]]


def test_provenance_is_joined_by_record_id_not_position(atlas):
    """_rows_to_arrays drops unlabelled rows, so positional joins mislabel.

    A null-Tm row in the middle of the query shifts every later row's
    provenance by one if the join is positional. That is silent and plausible,
    which is what makes it dangerous.
    """
    atlas.add([Record("GGGAAAGGGAAAGGGAAAGGG", "G4", "src", "experimental",
                      Condition.from_mapping({"k": 100.0, "ph": 7.0}),
                      tm=None, folded=1, source_id="unlabelled")])
    spec = _tm_spec()
    rows = atlas.query(kind=spec.kind, tiers=list(spec.tiers), label=spec.target)
    res = foldmanifest.assignments_for_spec(rows, spec)
    by_id = {r["record_id"]: r for r in rows}
    for entry in res["rows"]:
        assert entry["record_id"] in by_id
        assert entry["sequence"] == by_id[entry["record_id"]]["sequence"]


def test_a_head_below_its_row_minimum_gets_no_phantom_entry(atlas):
    spec = TaskSpec(**{**_tm_spec().__dict__, "min_rows": 10_000})
    assert foldmanifest.assignments_for_spec(rows_of(atlas), spec) is None


def test_write_produces_both_files_and_a_fingerprint(atlas, tmp_path):
    csv_path, meta_path, meta = foldmanifest.write(atlas, out_dir=tmp_path)
    assert csv_path.exists() and meta_path.exists()
    assert meta["rows_fingerprint"]
    assert meta["manifest_version"] == foldmanifest.MANIFEST_VERSION
    written = json.loads(meta_path.read_text())
    assert written["rows_fingerprint"] == meta["rows_fingerprint"]


def test_fingerprint_is_over_atlas_rows_not_assignments(atlas, tmp_path):
    """It must change when the atlas changes, which is the only thing it is for."""
    _, _, before = foldmanifest.write(atlas, out_dir=tmp_path)
    atlas.add([Record("GGGCCCGGGCCCGGGCCCGGG", "G4", "src", "experimental",
                      Condition.from_mapping({"k": 100.0, "ph": 7.0}),
                      tm=55.0, source_id="new")])
    _, _, after = foldmanifest.write(atlas, out_dir=tmp_path)
    assert before["rows_fingerprint"] != after["rows_fingerprint"]


def rows_of(atlas):
    spec = _tm_spec()
    return atlas.query(kind=spec.kind, tiers=list(spec.tiers), label=spec.target)


# --------------------------------------------------------------------------
# provenance: reconstructed is not historical
# --------------------------------------------------------------------------
def test_manifest_labels_itself_reconstructed_not_historical(atlas, tmp_path):
    _, _, meta = foldmanifest.write(atlas, out_dir=tmp_path)
    assert meta["assignment_provenance"] == foldmanifest.PROVENANCE_RECONSTRUCTED
    assert "evidence, not proof" in meta["provenance_meaning"][
        foldmanifest.PROVENANCE_RECONSTRUCTED]


def test_split_code_fingerprint_covers_the_functions_that_decide_the_split(atlas, tmp_path):
    fp = foldmanifest.split_code_fingerprint()
    assert "quadcond.models.train.build_folds" in fp["per_unit"]
    assert "quadcond.models.grouping.cluster_sequences" in fp["per_unit"]
    for unit, digest in fp["per_unit"].items():
        assert not digest.startswith("<unavailable"), f"{unit} source not readable"


def test_split_code_fingerprint_changes_when_grouping_code_changes(monkeypatch):
    """The point of the hash is that it moves when the split logic moves."""
    import quadcond.models.grouping as grouping
    before = foldmanifest.split_code_fingerprint()["combined"]

    def replacement(sequences, threshold=0.90, k=4, *, exact_max=6000):
        raise NotImplementedError

    monkeypatch.setattr(grouping, "cluster_sequences", replacement)
    after = foldmanifest.split_code_fingerprint()["combined"]
    assert before != after


def test_default_output_is_a_tracked_directory():
    """data/ is gitignored wholesale; a manifest nobody can commit is not an audit artifact."""
    assert "published_tools" in str(foldmanifest.DEFAULT_OUT)
    assert not str(foldmanifest.DEFAULT_OUT).startswith("data")
