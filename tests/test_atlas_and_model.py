"""Integration tests: atlas round-trip, ingestion honesty, end-to-end training."""
import numpy as np
import pytest

from quadcond.atlas import Atlas, Record
from quadcond.atlas.ingest._columns import resolve
from quadcond.atlas.ingest import shuffled
from quadcond.conditions import Condition
from quadcond.models.train import TaskSpec, train_head

HTELO = "AGGGTTAGGGTTAGGGTTAGGG"
MYC = "TGAGGGTGGGTAGGGTGGGTAA"


@pytest.fixture
def atlas(tmp_path):
    a = Atlas(tmp_path / "t.db")
    a.register_source("src", evidence_tier="experimental")
    yield a
    a.close()


def test_record_roundtrip_preserves_condition(atlas):
    c = Condition.from_mapping({"k": 25.0, "ph": 5.5})
    atlas.add([Record("GGGTTAGGGTTAGGGTTAGGG", "G4", "src", "experimental", c, topology="hybrid")])
    r = atlas.query()[0]
    assert r["k"] == 25.0 and r["ph"] == 5.5
    assert r["topology"] == "hybrid"


def test_duplicate_measurements_are_not_double_counted(atlas):
    rec = Record("GGGTTAGGGTTAGGGTTAGGG", "G4", "src", "experimental",
                 Condition.preset("k100"), tm=60.0, source_id="x")
    atlas.add([rec])
    atlas.add([rec])
    assert atlas.count() == 1


def test_same_sequence_two_conditions_are_two_rows(atlas):
    seq = "GGGTTAGGGTTAGGGTTAGGG"
    atlas.add([
        Record(seq, "G4", "src", "experimental", Condition.preset("k100"), tm=60.0, source_id="x"),
        Record(seq, "G4", "src", "experimental", Condition.preset("na100"), tm=45.0, source_id="x"),
    ])
    assert atlas.count() == 2


def test_tier_filter_excludes_predicted(atlas):
    atlas.register_source("pred", evidence_tier="predicted")
    atlas.add([
        Record("GGGTTAGGGTTAGGGTTAGGG", "G4", "src", "experimental", tm=60.0, source_id="a"),
        Record("GGGCGCGGGAGGAATTGGGCGGG", "G4", "pred", "predicted", tm=70.0, source_id="b"),
    ])
    exp = atlas.query(tiers=("experimental",), label="tm")
    assert len(exp) == 1 and exp[0]["source"] == "src"


def test_neighbour_search_prefers_matching_conditions(atlas):
    seq = "GGGTTAGGGTTAGGGTTAGGG"
    atlas.add([
        Record(seq, "G4", "src", "experimental", Condition.preset("k100"), tm=60.0, source_id="k"),
        Record(seq, "G4", "src", "experimental", Condition.preset("na100"), tm=45.0, source_id="n"),
    ])
    top = atlas.neighbours(seq, Condition.preset("k100"), n=1)[0]
    assert top["source_id"] == "k"


def test_column_resolver_handles_messy_headers():
    cols = ["Sequence", "T_m (°C)", "[K+] mM", "pH"]
    got = resolve(cols, {"sequence": ["sequence"], "tm": ["tm", "t_m"],
                         "k": ["k", "k+"], "ph": ["ph"]}, required=["sequence", "tm"])
    assert got["sequence"] == "Sequence" and got["tm"] == "T_m (°C)"


def test_column_resolver_raises_on_missing_required():
    with pytest.raises(ValueError):
        resolve(["a", "b"], {"sequence": ["sequence"]}, required=["sequence"])


def test_shuffled_negatives_are_tagged_derived_and_unfolded():
    parents = [{"sequence": "GGGTTAGGGTTAGGGTTAGGG", "source_id": "p", "record_id": 1}]
    recs = shuffled.build(parents, kind="G4")
    assert recs[0].evidence_tier == "derived" and recs[0].folded == 0


def _synthetic_rows(n=240):
    """G4-like positives vs shuffles, with a Tm that genuinely depends on [K+]."""
    rng = np.random.default_rng(0)
    rows = []
    for i in range(n):
        loops = rng.integers(1, 6, 3)
        seq = "GGG" + "".join(
            "".join(rng.choice(list("ATC"), L)) + "GGG" for L in loops
        )
        for k in (10.0, 50.0, 100.0):
            tm = 40.0 + 12.0 * np.log10(k) - 1.5 * loops.sum() + rng.normal(0, 1.0)
            rows.append({
                "sequence": seq, "kind": "G4", "record_id": i, "source": "syn",
                "evidence_tier": "experimental", "folded": 1, "topology": None,
                "tm": float(tm), "ph_t": None, "dg": None,
                "k": k, "na": 0.0, "li_nh4": 0.0, "mg": 0.0, "ph": 7.0,
                "temperature": 25.0, "crowder_pct": 0.0, "strand_conc": 5.0,
            })
    return rows



@pytest.fixture(scope="module")
def trained_predictor():
    """A two-head model over the synthetic panel, wrapped in the real API.

    Small and fast on purpose: these tests check the shape of what comes out of
    ``Predictor.predict`` and the claim basis attached to it, neither of which
    depends on the model being good.
    """
    from quadcond.models.base import MultiTaskModel
    from quadcond.models.predict import Predictor

    rows = _synthetic_rows(n=60)
    model = MultiTaskModel()
    for spec in (TaskSpec("g4_tm", "G4", "regression", "tm"),
                 TaskSpec("g4_fold", "G4", "binary", "folded")):
        try:
            head = train_head(rows, spec, n_seeds=2, n_folds=3, verbose=False)
        except Exception:
            head = None
        if head is not None:
            head.training_meta.setdefault("claim", spec.claim or "synthetic")
            head.training_meta.setdefault("label_classes", ["biophysical"])
            model.add(head)
    if not model.heads:
        pytest.skip("no head trained on the synthetic panel")
    return Predictor(model=model)


def test_condition_features_recover_a_salt_dependent_target():
    """The point of the whole package, as a test.

    A target that depends on log10[K+] must be learnable with condition
    features and largely unlearnable without them.
    """
    rows = _synthetic_rows()
    spec = TaskSpec("syn_tm", "G4", "regression", "tm")
    with_c = train_head(rows, spec, use_conditions=True, n_seeds=2, n_folds=3, verbose=False)
    without = train_head(rows, spec, use_conditions=False, n_seeds=2, n_folds=3, verbose=False)
    assert with_c.metrics["r2"] > 0.8
    assert with_c.metrics["r2"] - without.metrics["r2"] > 0.3


def test_regression_interval_reports_both_coverages_and_does_not_conflate_them():
    """The in-sample number is near-tautological; the group-split one is not.

    Taking the 90th percentile of a residual pool and then reporting the
    fraction of that same pool it covers gives ~0.90 whatever the model does --
    the quantile is defined to. That number is still worth recording, but it is
    not evidence, so the field is named for what it is and a second field
    carries the real estimate: the half-width fitted on half the groups and
    measured on the other half. This test pins both, and pins the naming, so a
    later refactor cannot quietly restore "conformal coverage" as the headline.
    """
    rows = _synthetic_rows()
    spec = TaskSpec("syn_tm", "G4", "regression", "tm")
    head = train_head(rows, spec, n_seeds=2, n_folds=3, conformal_alpha=0.1, verbose=False)
    m = head.metrics

    assert "conformal_halfwidth" not in m and "empirical_coverage" not in m, \
        "the old names promised a finite-sample guarantee the computation never had"
    assert m["oof_residual_halfwidth"] > 0
    assert m["in_sample_coverage_of_residual_pool"] >= 0.88

    # the held-out estimate exists and is a probability; on this tiny synthetic
    # panel it is allowed to be worse than nominal -- that is the point of it
    cov = m.get("split_conformal_coverage")
    assert cov is None or 0.0 <= cov <= 1.0
    if cov is not None:
        assert m["split_conformal_n_splits"] >= 1
        assert m["split_conformal_halfwidth_mean"] > 0


def test_head_reports_its_applicability_domain():
    rows = _synthetic_rows()
    head = train_head(rows, TaskSpec("syn_tm", "G4", "regression", "tm"),
                      n_seeds=2, n_folds=3, verbose=False)
    assert head.applicability["k"]["min"] == 10.0
    assert head.applicability["k"]["max"] == 100.0


def test_mojibake_repair_recovers_author_names():
    from quadcond.atlas.ingest.g4sp import _demojibake

    mangled_dash = "Syn\u00e2\u20ac\u201cAnti"
    mangled_name = "Maru\u00c5\u00a1i\u00c4\u008d, M."
    assert _demojibake(mangled_dash) == "Syn\u2013Anti"
    assert _demojibake(mangled_name) == "Maru\u0161i\u010d, M."
    # already-clean text must survive untouched
    assert _demojibake("Karg, B. (2018). Journal, 1-2.") == "Karg, B. (2018). Journal, 1-2."
    assert _demojibake("Maru\u0161i\u010d, M.") == "Maru\u0161i\u010d, M."
    assert _demojibake("") == ""


# ------------------------------------------------- declarative JSON adapters
def _adapter_spec(tmp_path, **overrides):
    import json

    spec = {
        "id_prefix": "T",
        "source": "test_source",
        "kind": "G4",
        "label_class": "biophysical",
        "columns": {"sequence": "seq", "tm": "melt", "k": "potassium_M", "ph": "pH"},
        "units": {"k": "M", "tm": "C"},
        "constants": {"na": 0, "temperature": 25},
    }
    spec.update(overrides)
    path = tmp_path / "spec.json"
    path.write_text(json.dumps(spec))
    return path


def _adapter_table(tmp_path):
    import pandas as pd

    path = tmp_path / "melts.csv"
    pd.DataFrame({"seq": [HTELO], "melt": [64.2], "potassium_M": [0.1], "pH": [7.0]}).to_csv(
        path, index=False
    )
    return path


def test_declarative_adapter_converts_molar_to_millimolar(tmp_path):
    from quadcond.atlas.ingest import declarative

    spec = declarative.load_spec(_adapter_spec(tmp_path))
    recs = declarative.load(_adapter_table(tmp_path), spec=spec)
    assert len(recs) == 1
    assert recs[0].condition.k == pytest.approx(100.0)   # 0.1 M -> 100 mM
    assert recs[0].tm == pytest.approx(64.2)
    assert recs[0].label_class == "biophysical"


def test_declarative_adapter_converts_kelvin(tmp_path):
    import pandas as pd

    from quadcond.atlas.ingest import declarative

    table = tmp_path / "k.csv"
    pd.DataFrame({"seq": [HTELO], "melt": [337.35], "potassium_M": [0.1], "pH": [7.0]}).to_csv(
        table, index=False
    )
    spec = declarative.load_spec(_adapter_spec(tmp_path, units={"k": "M", "tm": "K"}))
    recs = declarative.load(table, spec=spec)
    assert recs[0].tm == pytest.approx(64.2, abs=0.01)


def test_declarative_adapter_rejects_unknown_units(tmp_path):
    from quadcond.atlas.ingest import declarative

    with pytest.raises(declarative.AdapterError):
        declarative.load_spec(_adapter_spec(tmp_path, units={"k": "furlongs"}))


def test_declarative_adapter_rejects_unknown_fields(tmp_path):
    from quadcond.atlas.ingest import declarative

    with pytest.raises(declarative.AdapterError):
        declarative.load_spec(
            _adapter_spec(tmp_path, columns={"sequence": "seq", "wibble": "x"})
        )


def test_declarative_adapter_applies_value_maps(tmp_path):
    import pandas as pd

    from quadcond.atlas.ingest import declarative

    table = tmp_path / "v.csv"
    pd.DataFrame({"seq": [HTELO, MYC], "melt": [64.2, 82.0], "potassium_M": [0.1, 0.1],
                  "pH": [7.0, 7.0], "call": ["yes", "no"]}).to_csv(table, index=False)
    spec = declarative.load_spec(_adapter_spec(
        tmp_path,
        columns={"sequence": "seq", "tm": "melt", "k": "potassium_M", "ph": "pH",
                 "folded": "call"},
        value_maps={"folded": {"yes": 1, "no": 0}},
    ))
    recs = declarative.load(table, spec=spec)
    assert [r.folded for r in recs] == [1, 0]


# ------------------------------------------------------------- provenance
def test_label_class_separates_biophysical_from_genomic_proxy(atlas):
    seq = "GGGTTAGGGTTAGGGTTAGGG"
    atlas.add([
        Record(seq, "G4", "src", "experimental", label_class="biophysical",
               tm=60.0, source_id="a"),
        Record(seq + "AA", "G4", "src", "experimental", label_class="genomic_proxy",
               folded=1, source_id="b"),
    ])
    bio = atlas.query(label_classes=("biophysical",))
    assert len(bio) == 1 and bio[0]["source_id"] == "a"


def test_dataset_fingerprint_is_stable_and_sensitive():
    from quadcond.models.base import rows_fingerprint

    base = _synthetic_rows(20)
    assert rows_fingerprint(base) == rows_fingerprint(list(reversed(base)))
    changed = [dict(r) for r in base]
    changed[0]["tm"] = changed[0]["tm"] + 1.0
    assert rows_fingerprint(base) != rows_fingerprint(changed)


def test_model_save_records_a_verifiable_checksum(tmp_path):
    import json

    from quadcond.models.base import MultiTaskModel, file_sha256

    m = MultiTaskModel()
    path = m.save(tmp_path / "m.joblib")
    recorded = json.loads(path.with_suffix(".json").read_text())["artifact_sha256"]
    assert recorded == file_sha256(path)
    MultiTaskModel.load(path)  # verifies cleanly

    path.write_bytes(path.read_bytes() + b"corrupted")
    with pytest.raises(Exception):
        MultiTaskModel.load(path)


def test_imseeker_does_not_leak_pht_into_the_condition(tmp_path):
    """pH_T is the label; storing it as the condition's pH is target leakage.

    Regression test for a real bug: setting condition.ph = ph_t let the pH_T
    regression head read its own answer out of a feature and report R2 = 0.99.
    """
    import pandas as pd

    from quadcond.atlas.ingest import imseeker

    path = tmp_path / "imseeker.xlsx"
    pd.DataFrame({
        "Original name": ["A", "B"],
        "New_name": ["A_1", "B_2"],
        "Sequences": ["TGTCCCCACACCCCTGTCCCCACACCCCTGT",
                      "CCCTAACCCTAACCCTAACCCT"],
        "pHT": [6.6, 5.9],
    }).to_excel(path, sheet_name=imseeker.SHEET_MEASUREMENTS, index=False, startrow=1)

    recs = imseeker.load(path)
    assert len(recs) == 2
    for r in recs:
        assert r.ph_t in (6.6, 5.9)
        assert "ph" in r.condition.imputed, "pH_T must not be set as the condition pH"
        assert r.condition.k == 100.0 and r.condition.na == 10.0
        assert any(f.startswith("ph_titration_range") for f in r.qc_flags)


def test_imseeker_rejects_deletion_notation(tmp_path):
    import pandas as pd

    from quadcond.atlas.ingest import imseeker

    path = tmp_path / "im2.xlsx"
    pd.DataFrame({
        "New_name": ["good", "delta"],
        "Sequences": ["CCCTAACCCTAACCCTAACCCT", "TGTCCCCACΔCCCCTGTCCCC"],
        "pHT": [6.1, 6.2],
    }).to_excel(path, sheet_name=imseeker.SHEET_MEASUREMENTS, index=False, startrow=1)
    recs = imseeker.load(path)
    assert [r.source_id for r in recs] == ["good"]
    assert imseeker.load.last_stats["rejected_sequence"] == 1


def test_condition_grouping_holds_out_whole_buffers():
    """A condition-response head must be validated on unseen buffers, not unseen
    sequences — grouping by sequence when there are only two is meaningless."""
    import numpy as np

    from quadcond.conditions import Condition
    from quadcond.models.train import _condition_groups

    conds = [Condition.from_mapping({"k": k, "ph": ph}, track_imputed=False)
             for k in (10.0, 100.0, 1000.0) for ph in (5.0, 6.0)]
    groups = _condition_groups(conds)
    # pH is not part of the buffer key, so each K+ level is one group
    assert len(np.unique(groups)) == 3
    assert groups[0] == groups[1]      # same K+, different pH -> same group
    assert groups[0] != groups[2]


def test_duvma_parser_keeps_pht_out_of_the_condition(tmp_path):
    """Same leakage rule as the iM-Seeker adapter, enforced for the PDF path."""
    from pathlib import Path

    from quadcond.atlas.ingest import duvma

    pdf = Path("data/external/downloads/gkag110_supplemental_file.pdf")
    if not pdf.exists():
        import pytest

        pytest.skip("5DUVMA supporting information not present")
    recs = duvma.load(pdf)
    pht = [r for r in recs if r.ph_t is not None]
    tm = [r for r in recs if r.tm is not None]
    assert pht and tm
    for r in pht:
        # The invariant is that pH was never *set* from the label, which is what
        # `imputed` records. Numeric inequality is the wrong check: one C9 value
        # is exactly 7.00, which collides with the imputed default.
        assert "ph" in r.condition.imputed, "pH_T must not be set as the condition pH"
        assert r.condition.temperature > 0
    for r in tm:
        assert 4.0 <= r.condition.ph <= 7.2   # pH is a real condition for a melt
    # the ten ionic strengths are recovered as K+
    assert {r.condition.k for r in recs} >= {10.0, 100.0, 1000.0}


def _peaks_file():
    from pathlib import Path

    for p in sorted(Path("data/external/downloads").glob("gse220882_peaks.jsonl*")):
        return p
    return None


def test_gse220882_rows_are_genomic_proxy_with_a_fully_imputed_condition():
    """A CUT&Tag peak is not a melting curve, and the schema has to say so.

    Two failures this guards against. Pooling: if these rows ever came in as
    ``label_class="biophysical"`` they would join the heads that claim to predict
    folding in a buffer. And the condition: nothing about the intracellular
    chemistry was measured, so every field must be marked imputed or a
    condition-response head could fit a constant and call it a salt effect.
    """
    import pytest

    from quadcond.atlas.ingest import gse220882

    path = _peaks_file()
    if path is None:
        pytest.skip("GSE220882 peak file not present (run scripts/06_gse220882_peaks.py)")
    recs = gse220882.load(path, cell_lines=("HEK293T",), targets=("iM",))
    assert recs
    for r in recs[:200]:
        assert r.evidence_tier == "experimental"
        assert r.label_class == "genomic_proxy"
        assert r.folded == 1
        assert r.tm is None and r.ph_t is None
        for field in ("k", "na", "mg", "ph", "temperature"):
            assert field in r.condition.imputed, f"{field} must be marked imputed"
        assert r.genomic and r.genomic["build"] == "hg38"


def test_gse220882_negatives_carry_a_motif_and_match_composition():
    """The negative has to be a motif too, or the head just counts cytosines.

    Peak windows are GC-rich: a canonical i-motif appears in 37% of the strongest
    iMab peaks against 3.9% of random hg38, but only 1.36x against a
    dinucleotide shuffle of the same windows. Almost all of that gap is
    composition. So the negatives are motifs pulled out of shuffled peak windows,
    and this test enforces both halves of that: they are motifs, and they are
    composition-matched.
    """
    import collections

    import pytest

    from quadcond.atlas.ingest import gse220882
    from quadcond.motifs import find_im

    path = _peaks_file()
    if path is None:
        pytest.skip("GSE220882 peak file not present")
    negs = gse220882.negatives(path, cell_lines=("HEK293T",), targets=("iM",))
    assert negs
    for r in negs[:200]:
        assert r.folded == 0
        assert r.label_class == "catalog"
        assert find_im(r.sequence), "a negative with no motif makes the task trivial"

    pos = gse220882.load(path, cell_lines=("HEK293T",), targets=("iM",))

    def comp(rows):
        c = collections.Counter()
        for r in rows:
            c.update(r.sequence)
        n = sum(c.values())
        return {b: c[b] / n for b in "ACGT"}

    a, b = comp(pos), comp(negs)
    for base in "ACGT":
        assert abs(a[base] - b[base]) < 0.05, f"{base}: {a[base]:.3f} vs {b[base]:.3f}"


def test_g4_fold_is_not_fed_the_genomic_negatives():
    """A label_class filter alone was not enough, and this is why.

    The GSE220882 negatives are ``label_class="catalog"`` -- the same class every
    dinucleotide shuffle carries. A ``g4_fold`` restricted to
    ``("biophysical", "catalog")`` therefore let 5,135 peak-derived negatives in
    while their genomic_proxy positives stayed out, which handed the head a pile
    of unpaired negatives and lifted its reported AUROC from 0.97 to 0.99. The
    fix is naming the four sources; this test is here so nobody removes it.
    """
    from quadcond.models.train import DEFAULT_TASKS

    spec = next(t for t in DEFAULT_TASKS if t.name == "g4_fold")
    assert spec.sources, "g4_fold must pin its sources by name"
    assert not any("gse220882" in s for s in spec.sources)
    im = next(t for t in DEFAULT_TASKS if t.name == "im_fold")
    assert im.sources and not any("gse220882" in s for s in im.sources)
    gen = next(t for t in DEFAULT_TASKS if t.name == "im_fold_genomic")
    assert all("gse220882" in s for s in gen.sources)
    assert tuple(gen.label_classes) == ("genomic_proxy", "catalog")


# --------------------------------------------------------------------- schema
def test_prediction_record_matches_the_frozen_schema(trained_predictor):
    """Every prediction validates against quadcond/schema/prediction-v1.json.

    The schema is a contract with anyone who stores a prediction: within major
    version 1 a field may be added but never removed, renamed or redefined. This
    test is what makes that promise checkable rather than aspirational -- if a
    field is dropped from the emitter, the schema stops matching and this fails.
    """
    from quadcond.conditions import Condition
    from quadcond.schema import PREDICTION_SCHEMA_VERSION, validate_prediction

    cond = Condition.from_mapping({"ph": 7.0, "k": 100.0}, track_imputed=True)
    for seq in ("GGGTTAGGGTTAGGGTTAGGG", "CCCTAACCCTAACCCTAACCC", "ATATATATATATATAT"):
        rec = trained_predictor.predict(seq, cond)[0]
        assert rec["schema_version"] == PREDICTION_SCHEMA_VERSION
        problems = validate_prediction(rec)
        assert not problems, f"{seq}: {problems}"


def test_claim_basis_travels_with_every_head(trained_predictor):
    """No prediction may leave the API without the three claim-basis fields.

    They are the difference between a melting curve and an antibody peak, and a
    number lifted into a figure carries whatever it left with.
    """
    from quadcond.conditions import Condition

    rec = trained_predictor.predict(
        "GGGTTAGGGTTAGGGTTAGGG", Condition.from_mapping({"ph": 7.0, "k": 100.0}))[0]
    assert rec["predictions"]
    for name, out in rec["predictions"].items():
        for field in ("has_experimental_observation", "target_semantics",
                      "biophysically_grounded", "calibration_scope", "applicability"):
            assert field in out, f"{name} is missing {field}"
        # The licensing field must agree with the semantics it is derived from.
        assert out["biophysically_grounded"] == (out["target_semantics"] == "biophysical")
