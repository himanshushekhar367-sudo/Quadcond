"""Comparison-workflow contracts, checked without a trained artifact.

Every scan test that existed before this file skipped unless the canonical
27 MB model was on disk, which meant the routing, identity, refusal and
record-accounting defects those workflows shipped with were unreachable from
CI and from anyone's fresh checkout. They are not numerical defects. Which
strand a head is asked about, which mutant is scored, whether a refusal reaches
a consumer, and whether every supplied record is accounted for are all
properties of the control flow -- so they can be pinned with sentinel
estimators, and they should be, because a test that only runs where the model
happens to be present is a test the release gate cannot rely on.

The estimators here are deliberately trivial and deliberately not predictions.
They exist so that the real ``Head``, ``MultiTaskModel``, ``Predictor`` and
``scans`` code paths run end to end. Nothing in this file asserts anything about
accuracy, and no number it produces means anything.
"""
from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from quadcond import scans
from quadcond.conditions import Condition
from quadcond.models.base import Head, MultiTaskModel
from quadcond.models.predict import Predictor
from quadcond.motifs import revcomp

TEL22 = "GGGTTAGGGTTAGGGTTAGGG"
TEL22_C = revcomp(TEL22)                       # CCCTAACCCTAACCCTAACCC
#: Its own reverse complement -- the case that collapsed the strand tables.
PALINDROME = "AGGGTTAGGGTTAGGGTTAGGGCCCTAACCCTAACCCTAACCCT"

_COND_RANGE = {"min": 0.0, "max": 150.0, "n_unique": 12}


class _SentinelRegressor:
    """A deterministic function of the features. Not a prediction."""

    def __init__(self, offset: float) -> None:
        self.offset = offset

    def predict(self, X):
        X = np.asarray(X, dtype=float)
        return X.sum(axis=1) * 0.01 + self.offset


class _SentinelClassifier:
    def __init__(self, offset: float) -> None:
        self.offset = offset

    def predict_proba(self, X):
        X = np.asarray(X, dtype=float)
        p = 1.0 / (1.0 + np.exp(-(X.sum(axis=1) * 0.001 + self.offset)))
        return np.column_stack([1.0 - p, p])


def _head(name, kind, task, target, applicability, *, classes=()):
    maker = _SentinelRegressor if task == "regression" else _SentinelClassifier
    return Head(
        name=name, kind=kind, task=task, target=target,
        classes=list(classes),
        use_conditions=True,
        estimators=[maker(0.1 * i) for i in range(4)],
        conformal_q=1.0 if task == "regression" else None,
        metrics={},
        training_meta={"tiers": ["experimental"], "n_rows": 100,
                       "claim": "sentinel head for control-flow tests"},
        applicability=applicability,
    )


def _sentinel_model() -> MultiTaskModel:
    """Five heads spanning every task and both strands' motifs.

    ``im_tm_condition`` carries a ``sequence_allowlist`` holding only the C-rich
    telomeric strand, which is what the real head carries and what makes the
    strand-routing defect visible: the same duplex position is answerable from
    one strand and refused from the other.
    """
    common = {"k": dict(_COND_RANGE), "na": dict(_COND_RANGE),
              "length": {"min": 10, "max": 60}}
    heads = {
        "g4_tm": _head("g4_tm", "G4", "regression", "tm", dict(common)),
        "im_pht": _head("im_pht", "iM", "regression", "ph_t", dict(common)),
        "im_tm_condition": _head(
            "im_tm_condition", "iM", "regression", "tm",
            {**common, "sequence_allowlist": [TEL22_C, revcomp(PALINDROME)]}),
        "g4_fold": _head("g4_fold", "G4", "binary", "folded", dict(common)),
        "im_fold": _head("im_fold", "iM", "binary", "folded", dict(common)),
    }
    return MultiTaskModel(heads=heads, version="sentinel-0",
                          dataset_fingerprint="fingerprint-sentinel",
                          artifact_sha256="sha256-sentinel")


@dataclasses.dataclass
class _RecordingPredictor(Predictor):
    """A predictor that remembers which sequences each call was given."""

    calls: list = dataclasses.field(default_factory=list)

    def predict(self, sequences, condition=None, **kw):        # type: ignore[override]
        seqs = [sequences] if isinstance(sequences, str) else list(sequences)
        self.calls.append({"sequences": seqs, "heads": list(kw.get("heads") or [])})
        return super().predict(sequences, condition, **kw)


def _predictor(atlas=None) -> _RecordingPredictor:
    return _RecordingPredictor(model=_sentinel_model(), atlas=atlas)


# --------------------------------------------------------------------- F1
def test_a_condition_scan_asks_each_head_about_its_own_strand():
    """The available answer must not depend on which strand was pasted.

    ``evidence_payload`` and ``mutation_scan`` routed by motif; ``condition_scan``
    sent the typed sequence to every responsive head. So sweeping K+ over the
    G-rich telomeric strand had ``im_tm_condition`` -- a head whose allowlist
    holds the C-rich construct -- refuse at every point, while pasting the
    complement of the *same duplex position* made those identical queries
    in-domain. One duplex position, two answers, chosen by the user's clipboard.
    """
    g_typed = scans.condition_scan(_predictor(), TEL22, "k", [50.0, 100.0],
                                   n_neighbours=0)
    c_typed = scans.condition_scan(_predictor(), TEL22_C, "k", [50.0, 100.0],
                                   n_neighbours=0)

    for scan, typed in ((g_typed, TEL22), (c_typed, TEL22_C)):
        im = scan["series"]["im_tm_condition"]
        g4 = scan["series"]["g4_tm"]
        assert im["strand_sequence"] == TEL22_C, "the iM head answers on C-tracts"
        assert g4["strand_sequence"] == TEL22, "the G4 head answers on G-tracts"
        assert im["strand"] != g4["strand"], "opposite strands of one position"
        assert all(p["applicability"]["state"] != "refused" for p in im["points"]), \
            "the allowlisted construct is in-domain however the user typed it"
        assert scan["strand_sequence"] == {"+": typed, "-": revcomp(typed)}

    # And the same numbers, because it is the same molecule.
    assert ([p.get("value") for p in g_typed["series"]["im_tm_condition"]["points"]]
            == [p.get("value") for p in c_typed["series"]["im_tm_condition"]["points"]])


def test_a_condition_scan_names_the_strand_on_every_series_and_point():
    scan = scans.condition_scan(_predictor(), TEL22, "k", [50.0, 100.0],
                                n_neighbours=0)
    for name, series in scan["series"].items():
        assert series["strand"] in ("+", "-"), name
        assert series["strand_sequence"] in (TEL22, TEL22_C), name
        assert {p["strand"] for p in series["points"]} == {series["strand"]}
    assert scan["run_record"]["strand_by_head"]["g4_tm"] == "+"


def test_the_values_key_holds_the_axis_and_nothing_else():
    """The response declared ``number[]`` and shipped point dictionaries.

    A loop inside the summary reused the name that held the axis values, so the
    last head's list of points was what left the function under ``values``. The
    client's declared type said ``number[]``; a probe asking for [50, 100] got
    back two dicts.
    """
    scan = scans.condition_scan(_predictor(), TEL22, "k", [50.0, 100.0],
                                n_neighbours=0)
    assert scan["values"] == [50.0, 100.0]
    assert all(isinstance(v, float) for v in scan["values"])


# --------------------------------------------------------------------- F5
def test_a_condition_summary_covers_only_the_supported_points():
    """An unsupported tail must not decide the verdict for the supported part."""
    pred = _predictor()
    # 400 mM K+ is outside every head's recorded range, so the last two points
    # are extrapolation.
    scan = scans.condition_scan(pred, TEL22, "k", [50.0, 100.0, 300.0, 400.0],
                                n_neighbours=0)
    obs = scan["series"]["g4_tm"]["observed_response"]
    assert obs["summary_basis"] == "in_domain_points_only"
    assert obs["n_points"] == 4
    assert obs["n_supported"] == 2
    assert obs["n_out_of_domain"] == 2
    # The out-of-domain points are still returned, with their state, so a client
    # can show them explicitly rather than having them quietly dropped.
    assert len(scan["series"]["g4_tm"]["points"]) == 4
    assert scan["summary_note"]


def test_no_summary_claims_a_difference_is_unresolvable():
    """Wording check, and it is a claim check.

    Comparing a response range with a marginal residual half-width is a
    descriptive heuristic. The old verdict was the word ``flat`` and the note
    said "read this as no resolvable response", which is stronger than a range
    beside an error scale can support.
    """
    scan = scans.condition_scan(_predictor(), TEL22, "k", [50.0, 60.0],
                                n_neighbours=0)
    for series in scan["series"].values():
        note = series["observed_response"].get("note", "")
        assert "no resolvable response" not in note
        assert series["observed_response"]["verdict"] != "flat", \
            "'flat' reads as a measured insensitivity; the verdict is now " \
            "'below_error_scale' and says the comparison is descriptive"


# --------------------------------------------------------------------- F2
def test_a_self_complementary_sequence_scores_the_requested_mutant(monkeypatch):
    """The palindrome case, which collapsed the two strands into one entry.

    ``by_strand`` was keyed by the wild-type *string*. When a sequence is its own
    reverse complement the two keys are the same key, so the reverse-strand edit
    list overwrote the forward one: the cell said strand ``+`` and reported the
    mutation ``A1C``, while the sequence handed to the head was the reverse
    complement of that mutant -- a different edit at a different end of the
    molecule, with the label of the one the user asked for.
    """
    assert revcomp(PALINDROME) == PALINDROME, "the case under test"

    seen: list[list[str]] = []
    real = scans.featurize

    def spy(seqs, conds, **kw):
        seen.append(list(seqs))
        return real(seqs, conds, **kw)

    monkeypatch.setattr(scans, "featurize", spy)

    pred = _predictor()
    scan = scans.mutation_scan(pred, PALINDROME, heads=["g4_tm"], positions=[0])

    row = next(r for r in scan["substitutions"] if r["label"] == "A1C")
    cell = row["heads"]["g4_tm"]

    # The forward mutant edits position 0; the reverse mutant edits the last
    # base. They are different strings, and only one of them is the edit the
    # label names on the strand the cell reports.
    assert row["sequence"] != row["complement"]
    assert row["sequence"].startswith("C")
    assert row["complement"].endswith("G")

    # One featurize batch per head: [wild type, *mutants on that head's strand].
    assert len(seen) == 1, seen
    batch = seen[0]
    assert batch[0] == PALINDROME
    expected = [r["sequence"] if cell["strand"] == "+" else r["complement"]
                for r in scan["substitutions"]]
    assert batch[1:] == expected, (
        "the sequences actually featurized must be the mutants the cells' "
        "strand and labels name -- not their reverse complements")


def test_a_mutation_scan_keeps_both_strands_distinct():
    pred = _predictor()
    scan = scans.mutation_scan(pred, PALINDROME, heads=["g4_tm", "im_pht"],
                               positions=[0, 1])
    strands = scan["wild_type"]["strand_by_head"]
    assert set(strands.values()) <= {"+", "-"}
    assert scan["wild_type"]["strand_sequence"] == {"+": PALINDROME, "-": PALINDROME}
    for row in scan["substitutions"]:
        # Two different edits of one duplex position, never the same string.
        assert row["sequence"] != row["complement"]


# --------------------------------------------------------------------- F3
@pytest.mark.parametrize("head_name", ["g4_tm", "im_pht", "g4_fold", "im_fold"])
def test_every_head_type_refuses_a_query_that_is_not_a_sequence(head_name):
    """``ZZZZ`` cleans to ``NNNN``. No head may answer it -- including classifiers.

    The refusal branch used to sit after the binary and multiclass branches, so
    it was unreachable for a classification head: ``g4_fold`` and ``im_fold``
    returned a calibrated probability with ``applicability.refused: true`` and no
    top-level ``refused`` key, and the batch table -- which checks the top-level
    field -- rendered those numbers.
    """
    pred = _predictor()
    entry = pred.predict("ZZZZ", Condition(), n_neighbours=0)[0]["predictions"][head_name]
    assert entry["refused"] is True, "the flag every consumer keys on"
    assert entry["refusal_reason"]
    for key in ("value", "probability", "uncalibrated_probability",
                "posterior", "argmax", "confidence", "interval"):
        assert key not in entry, f"a refused head must not emit {key!r}"


def test_the_refusal_helpers_agree_with_the_serialised_entry():
    from quadcond import readout
    pred = _predictor()
    preds = pred.predict("ZZZZ", Condition(), n_neighbours=0)[0]["predictions"]
    for name, entry in preds.items():
        assert readout.is_refused(entry), name
        assert readout.readout(entry) is None, name


# --------------------------------------------------------------------- F4
def test_a_fasta_header_with_no_sequence_is_excluded_not_dropped():
    """Three headers in, three records out -- one of them excluded by name.

    The old parser only created a record when the buffer was non-empty, so an
    empty middle record produced two inputs, two results and an empty
    ``excluded`` list. The test that covered this asserted the count of two, so
    the loss was pinned in place rather than caught.
    """
    fasta = ">tel22\n" + TEL22 + "\n>empty\n\n>myc\nGGGGAGGGTGGGGAGGGTGGGG"
    parsed = scans.parse_input(fasta)
    assert parsed["mode"] == "fasta"
    assert [r["id"] for r in parsed["records"]] == ["tel22", "empty", "myc"]

    r = scans.batch_predict(_predictor(), fasta, [Condition()], heads=["g4_tm"])
    rec = r["run_record"]
    assert rec["n_sequences_in"] == 3
    assert rec["n_sequences_scored"] == 2
    assert rec["n_sequences_excluded"] == 1
    assert [e["id"] for e in r["excluded"]] == ["empty"]
    assert rec["reconciled"] is True


def test_bare_lines_are_one_sequence_each():
    """The panel promises "FASTA or one sequence per line" and meant it.

    Without headers the parser accumulated every line into one buffer, so
    ``AAAAAA`` and ``CCCCCC`` came back as the single sequence ``AAAAAACCCCCC``.
    """
    parsed = scans.parse_input("AAAAAA\nCCCCCC\n")
    assert parsed["mode"] == "lines"
    assert [r["sequence"] for r in parsed["records"]] == ["AAAAAA", "CCCCCC"]


def test_two_records_sharing_a_name_are_two_sequences():
    """Identifiers are the user's labels, not keys.

    ``n_sequences_scored`` counted distinct identifiers, so a duplicated FASTA
    header reported one scored sequence beside two rows -- a discrepancy that
    reads as a silent exclusion and is not one.
    """
    fasta = f">dup\n{TEL22}\n>dup\n{TEL22_C}"
    r = scans.batch_predict(_predictor(), fasta, [Condition()], heads=["g4_tm"])
    rec = r["run_record"]
    assert rec["n_sequences_in"] == 2
    assert rec["n_sequences_scored"] == 2
    assert rec["n_rows"] == 2
    assert rec["duplicate_identifiers"] == ["dup"]
    assert {row["record_index"] for row in r["results"]} == {0, 1}


def test_every_supplied_record_is_reconciled():
    fasta = f">a\n{TEL22}\n>b\n\n>c\nZZZZZZZZZZ\n>d\n{TEL22_C}"
    r = scans.batch_predict(_predictor(), fasta, [Condition()], heads=["g4_tm"])
    rec = r["run_record"]
    assert rec["n_sequences_in"] == 4
    assert rec["n_sequences_scored"] + rec["n_sequences_excluded"] == 4
    indices = ({row["record_index"] for row in r["results"]}
               | {e["record_index"] for e in r["excluded"]})
    assert indices == {0, 1, 2, 3}


# --------------------------------------------------------------------- F6
def test_no_interface_returns_a_combined_g4_versus_im_score():
    """The joint table and the signed preference are gone, not annotated.

    They multiplied two probabilities under an independence assumption that the
    two structures violate by construction, and subtracted two probabilities
    that were never on a common scale. A warning beside a number does not change
    what the number claims, and this one was reachable from the documented CLI
    subcommand and the Streamlit tab.
    """
    out = _predictor().competition(TEL22, Condition())
    for retired in ("joint", "preference_g4_minus_im", "interpretation"):
        assert retired not in out, f"{retired} was withdrawn in v0.4.8"
    assert out["combination_withdrawn"]["removed"]
    assert "g4" in out and "im" in out
    assert out["g_rich_strand"] == TEL22 and out["c_rich_strand"] == TEL22_C


# --------------------------------------------------------------------- F9
def test_all_three_comparisons_carry_the_same_run_record():
    """One provenance contract, not three hand-written dictionaries.

    Only ``/batch`` used to carry the model artifact hash and the atlas
    fingerprint, so a mutation or condition export could name a model version
    but not the file behind it.
    """
    pred = _predictor()
    records = [
        scans.mutation_scan(pred, TEL22, heads=["g4_tm"], positions=[0])["run_record"],
        scans.condition_scan(pred, TEL22, "k", [50.0, 100.0],
                             n_neighbours=0)["run_record"],
        scans.batch_predict(pred, [{"id": "a", "sequence": TEL22}],
                            [Condition()], heads=["g4_tm"])["run_record"],
    ]
    for rec in records:
        assert rec["quadcond_model_version"] == "sentinel-0"
        assert rec["model_artifact_sha256"] == "sha256-sentinel"
        assert rec["atlas_fingerprint"] == "fingerprint-sentinel"
        assert rec["prediction_schema_version"]
        assert rec["workflow"]


def test_the_mutation_response_states_what_the_paired_spread_is_not():
    scan = scans.mutation_scan(_predictor(), TEL22, heads=["g4_tm"], positions=[0])
    assert "not a validated interval" in scan["delta_note"]
    assert scan["run_record"]["delta_semantics"] == scan["delta_note"]


# ------------------------------------------------------------------ notes
def test_measured_neighbours_are_measurements():
    """A composition-matched shuffle is not a measurement of stability.

    ``Atlas.neighbours`` defaults to ``("experimental", "derived")``, and the
    derived tier holds distillation targets and shuffled controls. Listing one
    under "measured neighbours" is the claim-basis error the rest of this
    project's apparatus exists to prevent.
    """
    seen = {}

    class _FakeAtlas:
        def neighbours(self, sequence, condition, **kw):
            seen.update(kw)
            return []

    scans.condition_scan(_predictor(atlas=_FakeAtlas()), TEL22, "k",
                         [50.0, 100.0], n_neighbours=3)
    assert seen["tiers"] == ("experimental",)


def test_the_total_work_of_a_request_is_bounded_not_just_one_axis(monkeypatch):
    """Public serving pays for points x heads, so that product is what is capped.

    The per-axis and per-batch caps bound one dimension each of a request whose
    cost is a product. Sixty-four points is cheap against one head and expensive
    against twelve, and the old limits could not tell those apart.
    """
    pred = _predictor()
    monkeypatch.setattr(scans, "MAX_CONDITION_PREDICTIONS", 4)
    with pytest.raises(ValueError, match="budget"):
        scans.condition_scan(pred, TEL22, "k", [10.0, 20.0, 30.0, 40.0],
                             n_neighbours=0)
    monkeypatch.setattr(scans, "MAX_BATCH_PREDICTIONS", 3)
    with pytest.raises(ValueError, match="budget"):
        scans.batch_predict(pred, [{"id": str(i), "sequence": TEL22}
                                   for i in range(2)],
                            [Condition(), Condition()], heads=["g4_tm"])
