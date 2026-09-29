"""Tests for the GSE296171 ingest audit.

The audit's job is to refuse.  So these tests mostly check that it refuses for
the right reasons: an inverted score column, a sequence that collides with a
benchmark evaluation sequence once the molecule is erased, a mutant that cannot
be traced to its parent, and an absent fold manifest.  A check that only ever
passes would be worse than no check, because it would licence an ingest.

The real workbook is not committed (CC BY-NC, 2.4 MB), so every fixture here is
synthetic and hermetic.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import sqlite3
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
AUDIT = REPO / "benchmarks" / "ingest_audit" / "gse296171_overlap.py"

pytestmark = pytest.mark.skipif(not AUDIT.exists(), reason="audit script not present")


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("gse296171_overlap", AUDIT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# --------------------------------------------------------------------------
# sequence keys
# --------------------------------------------------------------------------
def test_rna_sequence_must_be_acgu(mod):
    assert mod.rna_seq("gggUUagggUUaggg") == "GGGUUAGGGUUAGGG"
    assert mod.rna_seq("GGGTTAGGG") is None, "T is not RNA; must not be silently accepted"
    assert mod.rna_seq("GGGNNAGGG") is None
    assert mod.rna_seq(None) is None


def test_dna_key_erases_the_molecule_on_purpose(mod):
    """The match key must collide, because a sequence-only featuriser collides."""
    assert mod.dna_key("GGGUUAGGG") == mod.dna_key("GGGTTAGGG") == "GGGTTAGGG"


def test_molecule_identity_survives_the_match(mod, tmp_path):
    """Matching on DNA must not overwrite what the row actually is."""
    wb = _workbook(tmp_path, [_row("GGGUUAGGGUUAGGGUUAGGG", -4.0, "stG4")])
    rows = mod.read_sheet(wb, mod.SHEET_SCORES)
    atlas = _atlas(tmp_path, [("GGGTTAGGGTTAGGGTTAGGG", 62.0)])
    summary, detail = mod.check_dna_rna_overlap(rows, mod.pool_column(rows), atlas)
    assert summary["distinct_sequences_overlapping"] == 1
    assert detail[0]["rna_sequence"] == "GGGUUAGGGUUAGGGUUAGGG"
    assert detail[0]["rna_molecule"] == "RNA"
    assert detail[0]["atlas_sequence"] == "GGGTTAGGGTTAGGGTTAGGG"
    assert detail[0]["atlas_molecule"] == "DNA"


# --------------------------------------------------------------------------
# statistics helpers
# --------------------------------------------------------------------------
def test_spearman_matches_a_known_value(mod):
    assert mod.spearman([1, 2, 3, 4, 5], [5, 4, 3, 2, 1]) == pytest.approx(-1.0)
    assert mod.spearman([1, 2, 3, 4, 5], [1, 2, 3, 4, 5]) == pytest.approx(1.0)
    assert mod.spearman([1, 2], [2, 1]) is None, "n<3 must not produce a correlation"


def test_spearman_handles_ties(mod):
    assert mod.spearman([1, 1, 2, 3], [1, 1, 2, 3]) == pytest.approx(1.0)


def test_correlation_is_reported_with_an_interval(mod):
    lo, hi = mod.fisher_ci(-0.6557, 55)
    assert lo < -0.6557 < hi
    assert mod.fisher_ci(0.5, 3) is None, "an interval needs n>3"


# --------------------------------------------------------------------------
# directionality -- the check that stops the label being learned backwards
# --------------------------------------------------------------------------
def _row(seq, kcl, pool, *, licl=None, pds=None, tract=3.0):
    return {
        mod_COL_SEQ: seq,
        mod_COL_KCL: kcl,
        mod_COL_LICL: licl if licl is not None else (None if kcl is None else kcl / 3.0),
        mod_COL_PDS: pds if pds is not None else (None if kcl is None else kcl - 1.5),
        "Minimum G-Tract (2 = GG, 3 = GGG, 4 = GGGG)": tract,
        _POOL: pool,
    }


def _directional_rows(sign=1.0):
    """A miniature pool with the real structure: controls flat, G4s shifted."""
    rows = []
    for i in range(12):
        rows.append(_row(f"{'A'*(10+i)}CGU", sign * 1.4, "neg_AAA", tract=2.0))
    for i in range(12):
        rows.append(_row(f"GG{'A'*i}GGAGGAGG", sign * -3.5, "nsG4_Mut", tract=3.0))
    for i, (tract, val) in enumerate([(2.0, 0.6), (3.0, -0.5), (4.0, -2.8)]):
        for j in range(6):
            rows.append(_row(f"GGG{'U'*j}AGGG{'A'*i}C", sign * val, "stG4", tract=tract))
    return rows


def test_directionality_confirms_negative_is_folded(mod):
    res = mod.check_directionality(_directional_rows(1.0), _POOL)
    assert res["status"] == "PASS"
    assert res["independent_tests_passed"] == res["independent_tests_total"] == 4


def test_directionality_rejects_an_inverted_dataset(mod):
    """If a future release flips the sign, the audit must not shrug."""
    res = mod.check_directionality(_directional_rows(-1.0), _POOL)
    assert res["status"] == "FAIL"
    assert res["failures"], "an inverted dataset must name which tests failed"


def test_directionality_does_not_consult_the_atlas(mod):
    """Direction must be established from controls alone, or it is circular.

    Enforced by signature: the check cannot read an atlas it is never handed.
    """
    import inspect
    params = list(inspect.signature(mod.check_directionality).parameters)
    assert params == ["rows", "pool_col"], (
        f"check_directionality takes {params}; giving it an atlas would make the "
        "direction test circular with the concordance result")


# --------------------------------------------------------------------------
# input integrity
# --------------------------------------------------------------------------
def test_integrity_fails_on_duplicate_sequences(mod):
    rows = [_row("GGGUUAGGG", -1.0, "stG4"), _row("GGGUUAGGG", -1.2, "stG4")]
    res = mod.check_input_integrity(rows, [], _POOL)
    assert res["status"] == "FAIL"
    assert any("unique" in f for f in res["failures"])


def test_integrity_fails_on_a_missing_score_value(mod):
    rows = [_row("GGGUUAGGG", None, "stG4")]
    res = mod.check_input_integrity(rows, [], _POOL)
    assert res["status"] == "FAIL"
    assert any("missing values" in f for f in res["failures"]), res["failures"]
    assert res["missing_values_by_score_column"][mod_COL_KCL] == 1


def test_integrity_fails_when_composition_does_not_reproduce_the_paper(mod):
    rows = [_row("GGGUUAGGG", -1.0, "stG4")]
    res = mod.check_input_integrity(rows, [], _POOL)
    assert res["status"] == "FAIL"
    assert any("composition" in f for f in res["failures"])
    assert res["pool_composition_expected"] == mod.EXPECTED_POOLS


# --------------------------------------------------------------------------
# scaffolds
# --------------------------------------------------------------------------
def test_scaffold_assignment_traces_mutants_to_their_parent(mod):
    assert mod.assign_scaffold("stG4_Pool_ns_VEGFA_1mut_100", "nsG4_Mut")[0] == "parent:VEGFA"
    assert mod.assign_scaffold("stG4_Pool_ns_U2AF2_AAAmut_7", "nsG4_Mut")[0] == "parent:U2AF2"


def test_synthetic_designs_are_not_given_a_parent(mod):
    g, basis = mod.assign_scaffold("stG4_Pool_GLV_42", "stG4")
    assert g.startswith("synthetic:")
    assert "no_parent" in basis


def test_unrecognised_identifier_is_unassigned_not_quietly_independent(mod):
    g, _ = mod.assign_scaffold("something_unexpected", "mystery")
    assert g == "UNASSIGNED"


def test_unassigned_sequences_fail_the_scaffold_gate(mod):
    rows = [_row("GGGUUAGGG", -3.0, "mystery")]
    summary, table = mod.check_scaffolds(rows, [], _POOL)
    assert summary["status"] == "FAIL"
    assert summary["unassigned_sequences"] == 1


def test_scaffold_groups_are_counted_per_parent(mod):
    muts = [{mod_COL_SEQ: f"GG{'A'*i}GGAGGAGG", "Pool identifier": f"stG4_Pool_ns_VEGFA_1mut_{i}"}
            for i in range(5)]
    rows = [_row(m[mod_COL_SEQ], -4.0, "nsG4_Mut") for m in muts]
    summary, table = mod.check_scaffolds(rows, muts, _POOL)
    assert summary["parent_scaffolds"] == 1
    assert summary["sequences_per_parent_scaffold"]["parent:VEGFA"] == 5
    assert summary["strong_folders_from_parent_scaffolds"] == 5
    assert summary["unassigned_sequences"] == 0


# --------------------------------------------------------------------------
# splits
# --------------------------------------------------------------------------
def _scaffold_table(n_per_parent=20, parents=("VEGFA", "U2AF2", "EWSR1")):
    out = []
    for p in parents:
        for i in range(n_per_parent):
            out.append({"rna_sequence": f"GG{p}{i}".replace("0", "A"),
                        "pool_class": "nsG4_Mut",
                        "pool_identifier": f"stG4_Pool_ns_{p}_1mut_{i}",
                        "group_key": f"parent:{p}", "assignment_basis": "mutant_of_parent",
                        "rt_stop_kcl": -4.0, "is_strong_folder": 1})
    return out


def test_grouped_split_removes_all_scaffold_leakage(mod):
    res = mod.check_splits(_scaffold_table(), n_splits=3)
    assert res["parent_scaffold_grouped_split"]["total_test_rows_with_scaffold_sibling_in_train"] == 0
    assert res["status"] == "PASS"


def test_random_split_leaks_and_the_audit_says_by_how_much(mod):
    res = mod.check_splits(_scaffold_table(), n_splits=3)
    leaked = res["random_sequence_level_split"]["total_test_rows_with_scaffold_sibling_in_train"]
    assert leaked > 0, "a random split over 3 scaffolds must leak"
    assert res["leakage_removed_by_grouping"] == leaked


def test_split_audit_fits_no_model(mod):
    """The comparison must be structural; an estimator here could be tuned."""
    src = AUDIT.read_text(encoding="utf-8")
    body = src.split("def check_splits")[1].split("\ndef ")[0]
    for forbidden in ("sklearn", "fit(", "predict(", "Ridge", "score("):
        assert forbidden not in body, f"check_splits must not {forbidden!r}"


# --------------------------------------------------------------------------
# leakage, fail-closed
# --------------------------------------------------------------------------
def test_missing_fold_manifest_fails_closed(mod, monkeypatch, tmp_path):
    monkeypatch.setattr(mod, "REPO", tmp_path)
    path, problems = mod.locate_fold_manifest()
    assert path is None
    assert problems and "fold assignment manifest" in problems[0]


def _manifest(tmp_path, fingerprint, *, with_meta=True, provenance="historical"):
    d = tmp_path / "data"
    d.mkdir(parents=True, exist_ok=True)
    (d / "fold_assignments.csv").write_text("head,record_id,fold\n", encoding="utf-8")
    if with_meta:
        from quadcond.models.foldmanifest import split_code_fingerprint
        (d / "fold_assignments_meta.json").write_text(
            json.dumps({"rows_fingerprint": fingerprint,
                        "assignment_provenance": provenance,
                        "split_code_fingerprint": split_code_fingerprint()}),
            encoding="utf-8")
    return d


def _model_card(tmp_path, fingerprint):
    d = tmp_path / "manuscript" / "inputs"
    d.mkdir(parents=True, exist_ok=True)
    (d / "quadcond_model.json").write_text(
        json.dumps({"dataset_fingerprint_sha256": fingerprint}), encoding="utf-8")
    return d


def test_fold_manifest_is_accepted_when_its_fingerprint_matches(mod, monkeypatch, tmp_path):
    monkeypatch.setattr(mod, "REPO", tmp_path)
    _manifest(tmp_path, "abc123")
    _model_card(tmp_path, "abc123")
    path, problems = mod.locate_fold_manifest()
    assert path is not None and not problems


def test_manifest_from_the_wrong_atlas_is_rejected(mod, monkeypatch, tmp_path):
    """A manifest that exists but describes other rows is worse than none.

    Accepting it would return a confident all-clear computed over the wrong
    folds, which is the one failure mode this gate exists to prevent.
    """
    monkeypatch.setattr(mod, "REPO", tmp_path)
    _manifest(tmp_path, "built_from_atlas_core")
    _model_card(tmp_path, "trained_on_full_atlas")
    path, problems = mod.locate_fold_manifest()
    assert path is None
    assert any("matches no shipped model" in p for p in problems)
    assert any("atlas_core" in p for p in problems)
    assert any("STALE" in p for p in problems)


def test_manifest_without_metadata_is_rejected(mod, monkeypatch, tmp_path):
    monkeypatch.setattr(mod, "REPO", tmp_path)
    _manifest(tmp_path, "abc123", with_meta=False)
    _model_card(tmp_path, "abc123")
    path, problems = mod.locate_fold_manifest()
    assert path is None
    assert any("provenance cannot be checked" in p for p in problems)


def test_manifest_is_rejected_when_no_model_card_exists(mod, monkeypatch, tmp_path):
    monkeypatch.setattr(mod, "REPO", tmp_path)
    _manifest(tmp_path, "abc123")
    path, problems = mod.locate_fold_manifest()
    assert path is None
    assert any("cannot be tied to a trained model" in p for p in problems)


def test_reconstructed_manifest_does_not_clear_the_gate(mod, monkeypatch, tmp_path):
    """Dataset identity is not fold identity.

    A manifest whose rows match the shipped model exactly still cannot prove
    which folds produced a published number, because the seed, fold count,
    cluster threshold and split code are all free to have changed since. It
    must be labelled evidence, not proof.
    """
    monkeypatch.setattr(mod, "REPO", tmp_path)
    _manifest(tmp_path, "abc123", provenance="reconstructed")
    _model_card(tmp_path, "abc123")
    path, problems = mod.locate_fold_manifest()
    assert path is None
    assert any("RECONSTRUCTED" in p for p in problems)
    assert any("evidence of fold identity, not proof" in p for p in problems)


def test_manifest_without_a_provenance_field_is_not_trusted(mod, monkeypatch, tmp_path):
    monkeypatch.setattr(mod, "REPO", tmp_path)
    d = tmp_path / "data"
    d.mkdir(parents=True, exist_ok=True)
    (d / "fold_assignments.csv").write_text("head,record_id,fold\n", encoding="utf-8")
    (d / "fold_assignments_meta.json").write_text(
        json.dumps({"rows_fingerprint": "abc123"}), encoding="utf-8")
    _model_card(tmp_path, "abc123")
    path, problems = mod.locate_fold_manifest()
    assert path is None, "an unlabelled manifest must not be assumed historical"


def test_a_matching_manifest_clears_the_fold_gate_but_not_a_collision(
        mod, monkeypatch, tmp_path):
    """Folds established is necessary, not sufficient."""
    monkeypatch.setattr(mod, "REPO", tmp_path)
    _manifest(tmp_path, "abc123")
    _model_card(tmp_path, "abc123")
    _benchmark_dir(tmp_path, ["GGGTTAGGGTTAGGGTTAGGG"])
    res = mod.check_leakage([_row("GGGUUAGGGUUAGGGUUAGGG", -4.0, "stG4")], [])
    assert res["fold_assignments_established"] is True
    assert res["status"] == "FAIL", "a real collision must still block"


def test_leakage_gate_blocks_on_a_benchmark_collision(mod, monkeypatch, tmp_path):
    _benchmark_dir(tmp_path, ["GGGTTAGGGTTAGGGTTAGGG"])
    monkeypatch.setattr(mod, "REPO", tmp_path)
    rows = [_row("GGGUUAGGGUUAGGGUUAGGG", -4.0, "stG4")]   # same string once U->T
    res = mod.check_leakage(rows, [])
    assert res["prohibited_overlap_sequences"] == 1
    assert res["status"] == "FAIL"


def test_leakage_gate_fails_closed_without_benchmark_files(mod, monkeypatch, tmp_path):
    monkeypatch.setattr(mod, "REPO", tmp_path)
    res = mod.check_leakage([_row("GGGUUAGGG", -1.0, "stG4")], [])
    assert res["status"] == "FAIL_CLOSED"
    assert not res["fold_assignments_established"]


def test_clean_dataset_still_fails_closed_without_folds(mod, monkeypatch, tmp_path):
    """No collisions is not the same as certified safe."""
    _benchmark_dir(tmp_path, ["AAAAAAAAAAAA"])
    monkeypatch.setattr(mod, "REPO", tmp_path)
    res = mod.check_leakage([_row("GGGUUAGGG", -1.0, "stG4")], [])
    assert res["prohibited_overlap_sequences"] == 0
    assert res["status"] == "FAIL_CLOSED", "absent folds must block even with zero overlap"


# --------------------------------------------------------------------------
# schema readiness
# --------------------------------------------------------------------------
def test_only_the_continuous_score_blocks_the_schema(mod):
    """Ligand fields and a new label class were over-strict; see the review note.

    An adapter may decline the +PDS arms, and `biophysical` is correct by
    schema.py's own definition for a structure measured in a defined buffer.
    What remains blocking is that there is nowhere to keep the continuous score.
    """
    res = mod.check_schema_readiness()
    assert res["missing_fields"] == ["rt_stop_score"]
    assert res["status"] == "BLOCKED"


def test_ligand_fields_are_conditional_not_required(mod):
    res = mod.check_schema_readiness()
    conditional = {c["field"]: c for c in res["conditionally_required_fields"]}
    assert set(conditional) == {"ligand", "ligand_conc"}
    for c in conditional.values():
        assert "PDS" in c["required_only_if"] or "PDS" in c["why"]


def test_condition_region_is_an_applicability_fence_not_a_label_class(mod):
    res = mod.check_schema_readiness()
    assert "biophysical is correct" in res["label_class_verdict"]
    fence = res["applicability_domain_fence_required"]
    assert fence["blocking_here"] is False
    assert "claims" in fence["enforced_by"]


def test_schema_readiness_does_not_modify_the_schema(mod):
    before = (REPO / "quadcond" / "atlas" / "schema.py").read_bytes()
    mod.check_schema_readiness()
    assert (REPO / "quadcond" / "atlas" / "schema.py").read_bytes() == before


def test_primary_conditions_are_recorded_as_reported_not_measured(mod):
    c = mod.PRIMARY_CONDITIONS
    assert c["ph"] == 8.3 and c["temperature_C"] == 40.0
    assert "not read from any data table" in c["reported_in"].lower()
    assert c["strand_conc"] is None, "a pool has no reportable strand concentration"


def test_conditions_are_in_reaction_not_stock(mod):
    """Quoting the 5X stock would overstate Mg2+ by 5.5x.

    The buffer is diluted into a 22 uL reaction. A condition-aware model given
    15 mM Mg2+ for a 2.7 mM experiment is being told something false about the
    only axis it exists to represent.
    """
    c = mod.PRIMARY_CONDITIONS
    assert c["in_reaction"]["mg_mM_total"] == pytest.approx(2.73, abs=0.01)
    assert c["in_reaction"]["monovalent_mM"] == pytest.approx(95.45, abs=0.01)
    assert c["nominal_stock"]["mg_mM"] == 15.0
    assert c["in_reaction"]["mg_mM_total"] < c["nominal_stock"]["mg_mM"]
    assert "mg_mM_free" not in c["in_reaction"], (
        "free Mg2+ is not knowable from the methods; only total may be recorded")


# --------------------------------------------------------------------------
# the audit must not be able to write to the atlas
# --------------------------------------------------------------------------
def test_atlas_is_opened_read_only(mod, tmp_path):
    db = _atlas(tmp_path, [("GGGTTAGGG", 60.0)])
    con = mod.open_atlas_readonly(db)
    with pytest.raises(sqlite3.OperationalError):
        con.execute("delete from records")
    con.close()


def test_concordance_is_labelled_association_only(mod, tmp_path):
    wb = _workbook(tmp_path, [_row("GGGUUAGGGUUAGGGUUAGGG", -4.0, "stG4")])
    rows = mod.read_sheet(wb, mod.SHEET_SCORES)
    atlas = _atlas(tmp_path, [("GGGTTAGGGTTAGGGTTAGGG", 62.0)])
    summary, _ = mod.check_dna_rna_overlap(rows, mod.pool_column(rows), atlas)
    text = summary["cross_molecule_concordance"]["interpretation"]
    assert "ASSOCIATION ONLY" in text
    assert "does NOT establish" in text


# --------------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------------
def _atlas(tmp_path, seq_tm):
    db = tmp_path / f"atlas_{abs(hash(tuple(seq_tm))) % 10**8}.db"
    con = sqlite3.connect(db)
    con.execute("""create table records (
        record_id integer primary key, sequence text, kind text, nucleic_acid text,
        label_class text, evidence_tier text, source text, method text,
        tm real, folded integer, k real, na real, ph real, temperature real)""")
    for s, tm in seq_tm:
        con.execute("insert into records (sequence,kind,nucleic_acid,label_class,"
                    "evidence_tier,source,method,tm,folded,k,na,ph,temperature) "
                    "values (?,'G4','DNA','biophysical','experimental',"
                    "'g4stab_experimental_tm','UV melting',?,1,100,0,7.0,25)", (s, tm))
    con.commit()
    con.close()
    return db


def _workbook(tmp_path, rows):
    from openpyxl import Workbook
    p = tmp_path / f"wb_{len(rows)}_{abs(hash(rows[0][mod_COL_SEQ])) % 10**6}.xlsx"
    wb = Workbook()
    ws = wb.active
    ws.title = "All Sequences Processed"
    header = list(rows[0])
    ws.append(header)
    for r in rows:
        ws.append([r.get(h) for h in header])
    wb.create_sheet("Mutation deltaWT Scores").append(["Shortened Sequences", "Pool identifier"])
    wb.save(p)
    return p


def _benchmark_dir(tmp_path, sequences):
    d = tmp_path / "benchmarks" / "published_tools" / "results"
    d.mkdir(parents=True, exist_ok=True)
    with (d / "tool_scores_all.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "sequence"])
        for i, s in enumerate(sequences):
            w.writerow([f"s{i}", s])
    with (d / "pairs_g4_tm.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["ref", "alt"])
    return d


# column names, resolved once so the row helper stays readable
mod_COL_SEQ = "Shortened Sequences"
mod_COL_KCL = "RT Stop Score in KCl (G in KCl/7dG in KCl)"
mod_COL_LICL = "RT Stop Score in LiCl (G in LiCl/7dG in LiCl)"
mod_COL_PDS = "RT Stop Score in KCl + PDS (G in KCl + PDS/7dG in KCl + PDS)"
_POOL = ("Pool Identifier (stG4 = Structural Pool, neg_AAA = Negative AAA controls, "
         "neg_ns = Negative natural sequence RBP binding site controls, nsG4_Mut = "
         "Natural Sequences and Mutational G4 Pool, pUG = pUG Oligos)")


def test_a_changed_split_code_makes_a_manifest_stale(mod, monkeypatch, tmp_path):
    """Rows can match perfectly while the partition is no longer reproducible."""
    monkeypatch.setattr(mod, "REPO", tmp_path)
    d = tmp_path / "data"
    d.mkdir(parents=True, exist_ok=True)
    (d / "fold_assignments.csv").write_text("head,record_id,fold\n", encoding="utf-8")
    (d / "fold_assignments_meta.json").write_text(json.dumps({
        "rows_fingerprint": "abc123",
        "assignment_provenance": "historical",
        "split_code_fingerprint": {"combined": "a_grouping_implementation_since_replaced",
                                   "per_unit": {}},
    }), encoding="utf-8")
    _model_card(tmp_path, "abc123")
    path, problems = mod.locate_fold_manifest()
    assert path is None, "matching rows must not clear a changed split implementation"
    assert any("STALE" in p for p in problems)
    assert any("split code has changed" in p for p in problems)
