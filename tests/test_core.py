from pathlib import Path
import random
from collections import Counter

import numpy as np
import pytest

from quadcond.conditions import Condition, condition_distance
from quadcond.features import featurize, feature_names
from quadcond.models.calibration import (BinaryCalibrator, TemperatureScaler,
                                         expected_calibration_error)
from quadcond.models.grouping import cluster_sequences
from quadcond.motifs import (find_g4, find_im, find_im_graph, g4hunter_mean,
                             parse_tracts, revcomp)
from quadcond.negatives import dinucleotide_shuffle, matched_negatives
from quadcond.thermo import folded_fraction_ph, folded_fraction_thermal

HTELO = "AAAGGGTTAGGGTTAGGGTTAGGGAA"
CTELO = "AATCCCTAACCCTAACCCTAACCCTTT"
MYC = "TGAGGGTGGGTAGGGTGGGTAA"


# --------------------------------------------------------------------- motifs
def test_finds_telomeric_g4():
    els = find_g4(HTELO)
    assert len(els) == 1
    assert els[0].sequence == "GGGTTAGGGTTAGGGTTAGGG"
    assert els[0].loops == (3, 3, 3)


def test_finds_telomeric_im():
    els = find_im(CTELO)
    assert len(els) == 1
    assert els[0].loops == (3, 3, 3)


def test_g4hunter_sign_is_strand_aware():
    assert g4hunter_mean(HTELO) > 1.0
    assert g4hunter_mean(revcomp(HTELO)) < -1.0


def test_graph_search_finds_at_least_the_canonical_call():
    graph = {e.sequence for e in find_im_graph(CTELO, representative="all")}
    assert any("CCCTAACCCTAACCCTAACCC" in s for s in graph)


def test_parse_tracts_rejects_short_sequences():
    assert parse_tracts("GGGTT", "G") is None


def test_both_strand_search_reports_coordinates_on_the_forward_frame():
    seq = "AAAA" + revcomp("GGGTTAGGGTTAGGGTTAGGG") + "AAAA"
    els = find_g4(seq, both_strands=True)
    assert any(e.strand == "-" for e in els)
    for e in els:
        assert 0 <= e.start < e.end <= len(seq)


# ----------------------------------------------------------------- conditions
def test_condition_defaults_are_flagged_as_imputed():
    c = Condition.from_mapping({"k": 50.0})
    assert "k" not in c.imputed and "ph" in c.imputed


def test_ionic_strength_counts_magnesium_three_times():
    assert Condition.from_mapping({"k": 0, "mg": 10}).ionic_strength == pytest.approx(30.0)


def test_potassium_and_sodium_buffers_are_far_apart():
    assert condition_distance(Condition.preset("k100"), Condition.preset("na100")) > 2.0


def test_presets_are_not_marked_imputed():
    assert Condition.preset("physiological").imputed == ()


# ------------------------------------------------------------------- features
def test_feature_matrix_matches_declared_names():
    X = featurize([HTELO, MYC], [Condition.preset("k100")] * 2, kind="G4")
    assert X.shape == (2, len(feature_names("G4")))
    assert np.isfinite(X).all()


def test_conditions_actually_change_the_features():
    a = featurize([HTELO], [Condition.preset("k100")], kind="G4")
    b = featurize([HTELO], [Condition.preset("na100")], kind="G4")
    assert not np.allclose(a, b)


def test_sequence_only_features_ignore_conditions():
    a = featurize([HTELO], [Condition.preset("k100")], kind="G4", use_conditions=False)
    b = featurize([HTELO], [Condition.preset("na100")], kind="G4", use_conditions=False)
    assert np.allclose(a, b)


# ------------------------------------------------------------------ negatives
def test_dinucleotide_shuffle_preserves_dinucleotide_counts():
    rng = random.Random(0)
    for seq in (HTELO, MYC, CTELO):
        sh = dinucleotide_shuffle(seq, rng)
        assert Counter(seq) == Counter(sh)
        di = lambda s: Counter(a + b for a, b in zip(s, s[1:]))
        assert di(seq) == di(sh)


def test_matched_negatives_are_composition_matched():
    negs = matched_negatives([HTELO, MYC], kind="G4", n_per_positive=2, seed=1)
    assert len(negs) == 4
    assert all(len(n) in {len(HTELO), len(MYC)} for n in negs)


# ---------------------------------------------------------------- calibration
def test_isotonic_calibration_reduces_ece():
    rng = np.random.default_rng(0)
    y = rng.binomial(1, 0.3, 4000)
    p = np.clip(rng.normal(y * 0.5 + 0.25, 0.15), 0.01, 0.99)
    cal = BinaryCalibrator("isotonic").fit(p[:2000], y[:2000])
    before = expected_calibration_error(y[2000:], p[2000:])
    after = expected_calibration_error(y[2000:], cal.transform(p[2000:]))
    assert after < before / 2


def test_temperature_scaling_preserves_argmax():
    rng = np.random.default_rng(1)
    proba = rng.dirichlet([2, 1, 1], 300)
    y = rng.integers(0, 3, 300)
    t = TemperatureScaler().fit(proba, y)
    assert (t.transform(proba).argmax(1) == proba.argmax(1)).all()


# ------------------------------------------------------------------- grouping
def test_clustering_groups_near_duplicates():
    seqs = [HTELO, HTELO + "T", "T" + HTELO, "ACGTACGTACGTACGTACGT"]
    g = cluster_sequences(seqs, threshold=0.9)
    assert g[0] == g[1] == g[2]
    assert g[3] != g[0]


# --------------------------------------------------------------------- thermo
def test_folded_fraction_is_half_at_the_midpoint():
    assert folded_fraction_thermal(60.0, 60.0) == pytest.approx(0.5)
    assert folded_fraction_ph(6.0, 6.0) == pytest.approx(0.5)


def test_folded_fraction_is_monotone_in_the_right_direction():
    assert folded_fraction_thermal(20, 60) > folded_fraction_thermal(80, 60)
    assert folded_fraction_ph(5.0, 6.0) > folded_fraction_ph(7.5, 6.0)


# --------------------------------------------------------------- buffer parsing
def test_buffer_parser_handles_the_common_literature_forms():
    from quadcond.buffers import parse_buffer

    p = parse_buffer("10 mM lithium cacodylate + 100 mM KCl")
    assert p.k == 100 and p.li_nh4 == 10 and p.na == 0 and p.usable

    p = parse_buffer("20 mM Tris (pH 7.2) + 99 mM NaCl + 1 mM KCl + 5 mM MgCl2")
    assert (p.na, p.k, p.mg, p.ph) == (99, 1, 5, 7.2)

    p = parse_buffer("PBS (pH 7.4)")
    assert p.na > 150 and 4 < p.k < 5


def test_buffer_parser_does_not_mistake_sodium_for_not_applicable():
    """Regression: an `\\bn/?a\\b` pattern silently ate every bare "Na" term."""
    from quadcond.buffers import parse_buffer

    p = parse_buffer("20 mM LiCl + 100 Na")
    assert p.na == 100 and p.li_nh4 == 20


def test_buffer_parser_flags_rather_than_guesses():
    from quadcond.buffers import parse_buffer

    unknown = parse_buffer("Buffer unsure (thought to be K)")
    assert not unknown.usable and "no_cation_recovered" in unknown.flags

    partial = parse_buffer("100 mM K + unknown Li cacodylate (pH 7.2)")
    assert partial.k == 100 and partial.li_nh4 == 0
    assert any(f.startswith("unquantified_component") for f in partial.flags)

    weird = parse_buffer("50 mM Unobtainium chloride")
    assert weird.unparsed and not weird.usable


def test_buffer_parser_converts_units_and_records_assumptions():
    from quadcond.buffers import parse_buffer

    assert parse_buffer("0.1 M KCl").k == 100.0
    p = parse_buffer("20 mM Kpi + 80 mM KCl")
    assert p.k == 110.0            # 20 * 1.5 phosphate convention + 80
    assert any("phosphate" in a for a in p.assumptions)


def test_censored_and_multiphasic_melting_points_are_not_coerced():
    from quadcond.atlas.ingest.g4stab_supp import _parse_tm

    assert _parse_tm(65.5) == (65.5, [])
    for raw in (">90", "<30", "low", "51/71", "36.6/39.1"):
        tm, flags = _parse_tm(raw)
        assert tm is None and flags, f"{raw!r} must not become a number"


# ------------------------------------------------------------ claim semantics
def test_claim_basis_separates_measured_from_biophysically_measured():
    """The v0.4.0 bug, pinned so it cannot come back.

    ``grounded_in_measurements`` answered "was anything measured", and a CUT&Tag
    peak makes that true. The model card therefore printed "measurement-grounded:
    yes" for a head whose label is antibody occupancy at a locus, on the same day
    the results report called that head a proxy. Three fields replace it, and the
    only one that licenses a folding claim is the third.
    """
    from quadcond import claims

    proxy = {"sources": ["gse220882_hek293t_im_cutandtag",
                         "shuffled::gse220882_hek293t_im_cutandtag"],
             "tiers": ["derived", "experimental"],
             "label_classes": ["catalog", "genomic_proxy"]}
    s = claims.semantics("im_fold_genomic", proxy, "binary")
    assert s["has_experimental_observation"] is True
    assert s["target_semantics"] == "genomic_proxy"
    assert s["biophysically_grounded"] is False

    biophys = {"sources": ["g4stab_experimental_tm"], "tiers": ["experimental"],
               "label_classes": ["biophysical"]}
    s = claims.semantics("g4_tm", biophys, "regression")
    assert s["has_experimental_observation"] and s["biophysically_grounded"]

    distilled = {"sources": ["g4stab_webdb_predicted"], "tiers": ["predicted"],
                 "label_classes": ["biophysical"]}
    s = claims.semantics("g4_tm_distilled", distilled, "regression")
    assert s["has_experimental_observation"] is False
    assert s["target_semantics"] == "predicted"


def test_shipped_model_splits_seven_three_two():
    """The headline is counted from the artifact, not typed into a document.

    This test earned its keep: `locus_state` made the split 7/3/2 and four
    documents plus the service headline still said 7/2/2. The counted number
    disagreed with the typed one and the typed one was wrong, which is the
    entire reason the count is taken from the artifact.
    """
    import pytest

    from quadcond import assets, claims
    from quadcond.models.base import MultiTaskModel

    try:
        path = assets.resolve("model")
    except assets.AssetError:
        pytest.skip("model artifact not available in this environment")
    m = MultiTaskModel.load(path)
    t = claims.tally(m.heads)
    assert (t["biophysical"], t["genomic_proxy"], t["derived"] + t["predicted"]) == (7, 3, 2)
    assert len(m.heads) == 12
    for name in ("im_fold_genomic", "g4_fold_genomic"):
        meta = m.heads[name].training_meta
        assert meta["biophysically_grounded"] is False
        # the claim string has to disown the folding reading in words, not only
        # in a metadata flag -- the string is what ends up quoted in figures
        assert "GENOMIC-PROXY SCORE, NOT A FOLDING PROBABILITY" in meta["claim"]


def test_an_asset_with_the_right_name_and_wrong_bytes_is_refused(tmp_path):
    """Name-based trust is the failure this guards against.

    ``artifacts/quadcond_model.joblib`` is a filename several versions of this
    project have used. Loading whichever one happens to be on disk would produce
    numbers nobody could later attribute to a release, so a checksum mismatch is
    an error rather than a warning -- and it is checked before the bytes reach
    the unpickler, since a truncated joblib archive can take the process out with
    an OOM kill instead of raising.
    """
    import pytest

    from quadcond import assets

    a = next(iter(assets.assets().values()))
    impostor = tmp_path / a.filename
    impostor.write_bytes(b"not the model you are looking for")
    with pytest.raises(assets.AssetError) as e:
        assets.verify(impostor, a)
    assert "Refusing" in str(e.value)


def test_atlas_identity_survives_being_opened():
    """A read must not change what the atlas *is*.

    SQLite rewrites page metadata whenever a database is opened, so the file
    checksum of an untouched, byte-reproducible atlas differs after a single
    read. Verifying atlases by file hash would make reproducibility look broken;
    the identity is the ordered hash of the rows instead.
    """
    import pytest

    from quadcond import assets
    from quadcond.atlas import Atlas

    try:
        db = assets.resolve("atlas")
    except assets.AssetError:
        pytest.skip("atlas not available in this environment")
    before = assets.atlas_fingerprint(db)
    Atlas(db).close()          # opening writes to the meta table
    assert assets.atlas_fingerprint(db) == before


# ------------------------------------------------------------ v0.4.3 fixes
def test_an_unlisted_sequence_gets_no_value_at_all():
    """A refusal must remove the number, not annotate it.

    v0.4.2 called this "refusing rather than extrapolating" in four documents
    while still emitting `value`, still printing it from the CLI, and still
    describing the value as optionally present in its own JSON schema. The flag
    was there; the number was there next to it; every consumer was trusted to
    remember which one won.
    """
    import pytest

    from quadcond import assets
    from quadcond.conditions import Condition
    from quadcond.models.predict import Predictor

    try:
        path = assets.resolve("model")
    except assets.AssetError:
        pytest.skip("model artifact not available in this environment")
    pred = Predictor.load(path)
    if "im_tm_condition" not in pred.model.heads:
        pytest.skip("condition head not in this artifact")

    cond = Condition.from_mapping({"ph": 7.0, "k": 100.0, "temperature": 25.0})
    out = pred.predict("CCCCTCCCCTCCCCTCCCC", cond,
                       heads=["im_tm_condition"])[0]["predictions"]["im_tm_condition"]

    assert out["refused"] is True
    assert "value" not in out, "a refused head must not emit a value"
    assert "interval" not in out
    assert "folded_fraction_at_condition" not in out
    assert out["applicability"]["in_domain"] is False
    assert "no prediction to report" in out["refusal_reason"]

    # and a sequence it WAS trained on still answers
    allow = pred.model.heads["im_tm_condition"].applicability["sequence_allowlist"]
    ok = pred.predict(allow[0], cond, heads=["im_tm_condition"])[0]["predictions"]
    assert "value" in ok["im_tm_condition"]
    assert not ok["im_tm_condition"].get("refused")


def test_every_condition_field_is_range_checked():
    """crowder_pct and strand_conc were in every applicability table and checked
    by nothing, so a query at 40% PEG passed as in-domain against training data
    that never left 0%."""
    import pytest

    from quadcond import assets
    from quadcond.conditions import Condition
    from quadcond.models.predict import Predictor

    try:
        path = assets.resolve("model")
    except assets.AssetError:
        pytest.skip("model artifact not available in this environment")
    pred = Predictor.load(path)
    cond = Condition.from_mapping(
        {"ph": 7.0, "k": 100.0, "crowder_pct": 40.0, "strand_conc": 500.0})
    ap = pred.predict("AGGGTTAGGGTTAGGGTTAGGG", cond,
                      heads=["g4_tm"])[0]["predictions"]["g4_tm"]["applicability"]
    joined = " ".join(ap["warnings"])
    assert "crowder_pct" in joined, "40% crowder must not pass silently"
    assert "strand_conc" in joined, "500 uM strand must not pass silently"


def test_a_tm_head_is_not_domain_checked_on_the_query_temperature():
    """It predicts a temperature; it does not consume one.

    Every G4STAB row carries temperature=25 as a metadata default on a Tm
    measurement, so the "training range" is [25, 25] and any physiological query
    was flagged out-of-domain for asking a question this field does not answer.
    """
    import pytest

    from quadcond import assets
    from quadcond.conditions import Condition
    from quadcond.models.predict import Predictor

    try:
        path = assets.resolve("model")
    except assets.AssetError:
        pytest.skip("model artifact not available in this environment")
    pred = Predictor.load(path)
    cond = Condition.from_mapping({"ph": 7.0, "k": 100.0, "na": 10.0, "temperature": 37.0})
    ap = pred.predict("AGGGTTAGGGTTAGGGTTAGGG", cond,
                      heads=["g4_tm"])[0]["predictions"]["g4_tm"]["applicability"]
    assert ap["in_domain"] is True, ap["warnings"]


def test_each_head_is_evaluated_on_the_strand_carrying_its_motif():
    """The defect this guards is a category error, not an inaccuracy.

    The evidence endpoint detected motifs on both strands and then called every
    head on the *input* strand. For a G-rich sequence the i-motif is on the
    reverse complement, so `im_pht` -- a model of C-tract sequences -- was being
    asked about a G-tract one and returning a number for it.
    """
    import pytest

    from quadcond import assets, service

    try:
        assets.resolve("model")
    except assets.AssetError:
        pytest.skip("model artifact not available in this environment")

    out = service.evidence_payload({
        "sequence": "AGGGTTAGGGTTAGGGTTAGGG",
        "k": 100.0, "na": 10.0, "ph": 7.0, "temperature": 37.0,
    })
    by_kind = {c["kind"]: c for c in out["evidence"]}
    assert "g-quadruplex" in by_kind and "i-motif" in by_kind

    g4, im = by_kind["g-quadruplex"], by_kind["i-motif"]
    assert g4["strand"] == "+" and g4["strand_sequence"].startswith("AGGG")
    assert im["strand"] == "-", "the i-motif is on the complement of a G-rich input"
    assert im["strand_sequence"].startswith("CCC"), im["strand_sequence"]
    assert "G" not in im["strand_sequence"][:3]


def test_the_evidence_endpoint_does_not_serve_a_partition_function():
    """No `ensemble`, no per-card probability, nothing summing to one.

    exp(-dG/RT) over {G4, i-motif, unfolded} omits the duplex -- the state that
    actually competes with both at physiological conditions -- along with strand
    concentration and stoichiometry. Percentages over an incomplete state space
    are not uncertain; they are about something else.
    """
    import pytest

    from quadcond import assets, service

    try:
        assets.resolve("model")
    except assets.AssetError:
        pytest.skip("model artifact not available in this environment")

    out = service.evidence_payload({"sequence": "AGGGTTAGGGTTAGGGTTAGGG",
                                    "k": 100.0, "ph": 7.0})
    assert "ensemble" not in out
    for card in out["evidence"]:
        # Any key containing "probab", not just the literal "probability". The
        # earlier version checked the exact string and passed for two releases
        # while every card carried `probabilityOfFolding` -- a number from
        # g4_fold/im_fold sitting under g4_tm/im_pht's source, claim basis and
        # applicability. Whatever a reader took from that card, they took it
        # about the wrong model.
        stray = [k for k in card if "probab" in k.lower()]
        assert not stray, f"{card['kind']} carries {stray} under {card['source']}'s provenance"
        assert "deltaG" not in card
        assert card["foldedFraction_note"], "a folded fraction needs its caveat attached"


def test_the_discrimination_scores_carry_their_own_provenance():
    """g4_fold / im_fold answer a narrower question than the Tm and pH_T heads.

    They are reported separately, each with its own head name, strand, claim
    basis and applicability, rather than riding inside a card describing a
    different model.
    """
    out = _payload()
    by_head = {d["head"]: d for d in out["discrimination"]}
    assert set(by_head) == {"g4_fold", "im_fold"}

    for name, d in by_head.items():
        assert 0.0 <= d["probability"] <= 1.0
        assert d["claim_basis"], f"{name} has no claim basis"
        assert "applicability" in d, f"{name} has no applicability"
        assert "NOT P(this sequence folds)" in d["note"]
        assert d["strand"] in ("+", "-")

    # routed like everything else: im_fold scores the C-rich strand of a G-rich
    # input, and agrees with the structural card about which strand that is
    im_card = next(c for c in out["evidence"] if c["kind"] == "i-motif")
    assert by_head["im_fold"]["strand"] == im_card["strand"] == "-"
    assert by_head["im_fold"]["strand_sequence"] == im_card["strand_sequence"]
    assert by_head["g4_fold"]["strand"] == "+"


def test_the_peak_overlap_head_is_named_and_gated_for_what_it_is():
    """`locus_state` invited a reading the experiment does not support.

    BG4 and iMab CUT&Tag were parallel reactions on separate aliquots: a 'both'
    window is population-level peak co-localisation, not two structures on one
    molecule. It is also trained exclusively on 201-nt windows.
    """
    import pytest

    from quadcond import assets
    from quadcond.models.base import MultiTaskModel

    try:
        path = assets.resolve("model")
    except assets.AssetError:
        pytest.skip("model artifact not available in this environment")
    m = MultiTaskModel.load(path)

    assert "locus_state" not in m.heads, "the misleading name must not come back"
    head = m.heads["locus_peak_overlap_state"]
    claim = head.training_meta["claim"]
    assert "NOT CO-OCCUPANCY" in claim
    assert "PARALLEL REACTIONS" in claim
    assert "0.156" in claim, "iM-only recall belongs in the claim"
    assert head.applicability["length"]["min"] == 201.0
    assert head.applicability["length"]["max"] == 201.0


# ------------------------------------------------------------ v0.4.4 fixes
def _payload(**over):
    import pytest

    from quadcond import assets, service

    try:
        assets.resolve("model")
    except assets.AssetError:
        pytest.skip("model artifact not available in this environment")
    body = {"sequence": "AGGGTTAGGGTTAGGGTTAGGG",
            "k": 100.0, "na": 10.0, "ph": 7.0, "temperature": 37.0}
    body.update(over)
    return service.evidence_payload(body)


def test_the_genomic_proxy_heads_are_routed_by_strand_too():
    """v0.4.3 fixed the structural cards and left these two pinned to forward.

    Same category error, one layer over: `im_fold_genomic` is trained on iMab
    peak windows and was being asked about the G-rich strand of a G-rich input,
    while the i-motif card beside it correctly used the reverse complement. Two
    cards in one response describing different molecules, with nothing saying so.
    """
    out = _payload()
    by_head = {g["head"]: g for g in out["genomic_evidence"]}
    assert set(by_head) == {"g4_fold_genomic", "im_fold_genomic"}

    g4, im = by_head["g4_fold_genomic"], by_head["im_fold_genomic"]
    assert g4["strand"] == "+" and g4["strand_sequence"].startswith("AGGG")
    assert im["strand"] == "-", "iMab scores the C-rich strand"
    assert im["strand_sequence"].startswith("CCC")

    # and the proxy agrees with the structural card about which strand is which
    im_card = next(c for c in out["evidence"] if c["kind"] == "i-motif")
    assert im_card["strand"] == im["strand"]
    assert im_card["strand_sequence"] == im["strand_sequence"]


def test_a_c_rich_input_routes_every_head_the_other_way():
    """The routing must follow the motif, not the sign of a hardcoded default."""
    out = _payload(sequence="CCCTAACCCTAACCCTAACCCT")
    by_kind = {c["kind"]: c for c in out["evidence"]}
    assert by_kind["i-motif"]["strand"] == "+"
    assert by_kind["g-quadruplex"]["strand"] == "-"
    by_head = {g["head"]: g for g in out["genomic_evidence"]}
    assert by_head["im_fold_genomic"]["strand"] == "+"
    assert by_head["g4_fold_genomic"]["strand"] == "-"


def test_the_locus_head_is_asked_in_the_orientation_it_was_trained_in():
    """The ingest canonicalises every window to its G-richer strand.

    Evaluating it on whichever strand the caller typed asks it about a
    representation it never saw -- and it has no motif gate, so nothing
    downstream would have flagged the mismatch.
    """
    out = _payload(sequence="CCCTAACCCTAACCCTAACCCT")
    joint = out["locus_peak_overlap_state"]
    assert joint is not None
    assert joint["strand"] == "-", "a C-rich input must be flipped to G-rich"
    assert joint["strand_sequence"].count("G") >= joint["strand_sequence"].count("C")


def test_the_locus_head_gets_no_i_motif_motif_gate():
    """`_motif_gate` dispatched G4-or-else, and `locus` fell into the else.

    A 201-nt G-rich window -- exactly the shape the head was trained on -- was
    reported out of domain solely for lacking a four-C-tract motif.
    """
    import pytest

    from quadcond import assets
    from quadcond.conditions import Condition
    from quadcond.models.predict import Predictor

    try:
        path = assets.resolve("model")
    except assets.AssetError:
        pytest.skip("model artifact not available in this environment")
    pred = Predictor.load(path)
    head = pred.model.heads["locus_peak_overlap_state"]

    seq = (("GGGTTAGGGTTAGGGTTAGGG" + "ATCGATCGAT") * 7)[:201]
    assert len(seq) == 201

    # the nominal intracellular condition the peak rows carry
    cond = Condition.from_mapping({
        "k": 140.0, "na": 10.0, "li_nh4": 0.0, "mg": 0.5,
        "ph": 7.2, "temperature": 37.0, "crowder_pct": 0.0, "strand_conc": 5.0,
    })
    ap = pred.predict(seq, cond, heads=["locus_peak_overlap_state"])[0][
        "predictions"]["locus_peak_overlap_state"]["applicability"]
    assert ap["in_domain"] is True, ap["warnings"]
    assert not any("i-motif" in w for w in ap["warnings"])

    # length is still the real gate, and it still bites
    short = pred.predict("AGGGTTAGGGTTAGGGTTAGGG", cond,
                         heads=["locus_peak_overlap_state"])[0][
        "predictions"]["locus_peak_overlap_state"]["applicability"]
    assert short["in_domain"] is False
    assert any("length" in w for w in short["warnings"])
    assert head.applicability["length"]["min"] == 201.0


def test_the_locus_block_carries_its_applicability_to_the_client():
    """The service was sending it and the viewer's type dropped it, so the panel
    rendered a 201-nt-window posterior for a 22-nt oligo with no warning."""
    joint = _payload()["locus_peak_overlap_state"]
    assert "applicability" in joint
    assert joint["applicability"]["in_domain"] is False
    assert any("length" in w for w in joint["applicability"]["warnings"])


def test_the_evidence_response_carries_no_ensemble_calibration():
    """The ΔG calibration existed to scale a Boltzmann step this endpoint no
    longer performs. Shipping it beside strand-resolved cards would let a caller
    rebuild the partition function that was withdrawn."""
    out = _payload()
    assert "calibration" not in out
    assert "RT_kcal" not in out
    assert "ensemble" not in out


def test_the_locus_record_itself_does_not_assert_co_occupancy():
    """Inspect the RECORD, not the docstring.

    The previous version of this test read `__doc__` and passed while every
    record still carried
    ``label_is_antibody_co_occupancy_not_thermodynamic_competition`` in its
    qc_flags -- the exact persistent field the test was written to protect. A
    flag written into the atlas outlives every document that describes it, so
    the assertion has to run against what `_record()` actually emits.
    """
    from quadcond.atlas.ingest import gse220882_loci as g

    locus = {"chrom": "chr1", "start": 1000, "end": 1201, "summit": 1100,
             "window": "GGGTTAGGGTTAGGGTTAGGG" * 9 + "GGGTTAGGGTTA",
             "state": "both", "n_peaks": 2}
    seq, orientation = g.g_rich_orientation(locus["window"])
    rec = g._record(seq, "both", locus, source=g.SOURCE, source_id="chr1:1000-1201",
                    orientation=orientation, extra=[])
    flags = " ".join(rec.qc_flags)

    assert "co_occupancy" not in flags.replace("not_co_occupancy", ""), flags
    assert "peak_overlap" in flags
    assert "parallel_reactions_on_separate_aliquots" in flags
    assert "training_window_nt=201" in flags
    assert any(f.startswith("locus_peak_overlap_state=") for f in rec.qc_flags)

    # and the prose that describes the rows agrees with the rows
    doc = g.__doc__ or ""
    assert "PARALLEL REACTIONS" in doc
    assert "not co-occupancy" in doc.lower()
    assert "co-occupied" not in doc.lower()


def test_the_registered_source_notes_do_not_assert_co_occupancy():
    """The source table is written once at build time and read forever."""
    from quadcond.atlas import Atlas
    from quadcond.atlas.ingest import gse220882_loci as g

    import tempfile

    with tempfile.TemporaryDirectory() as d:
        atlas = Atlas(Path(d) / "t.db")
        try:
            g.register(atlas)
            rows = {r["source"]: dict(r) for r in
                    atlas.conn.execute("SELECT * FROM sources")}
        finally:
            atlas.close()

    for src, row in rows.items():
        blob = f"{row.get('title', '')} {row.get('notes', '')}".lower()
        assert "co-occupancy" not in blob.replace("not co-occupancy", ""), (src, blob)


def test_the_buffer_count_has_one_definition():
    """Two counts of the same thing, differing, with neither defined.

    The report's prose said 261 buffers from a live query while its atlas table
    said 255 from the model's training snapshot, because the snapshot's key
    omitted Mg2+. A buffer is one distinct (K, Na, Li/NH4, Mg, pH) combination,
    everywhere.
    """
    import sqlite3

    import pytest

    from pathlib import Path

    from quadcond.atlas import Atlas

    db = Path("data/atlas.db")
    if not db.exists():
        pytest.skip("atlas not available in this environment")

    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    try:
        direct = con.execute(
            "SELECT COUNT(DISTINCT k||'_'||na||'_'||li_nh4||'_'||mg||'_'||ph) "
            "FROM records WHERE source='g4stab_experimental_tm'").fetchone()[0]
    finally:
        con.close()

    summary = Atlas(db).summary()
    row = summary[summary["source"] == "g4stab_experimental_tm"].iloc[0]
    assert int(row["n_conditions"]) == direct, (
        "Atlas.summary() and a direct count must use the same buffer key")


def test_the_asset_manifest_head_count_is_not_prose():
    """It said "the eleven trained heads" for two releases after there were 12."""
    import json

    import pytest

    from quadcond.assets import MANIFEST_PATH

    m = json.loads(MANIFEST_PATH.read_text())
    models = [a for a in m["assets"] if a["kind"] == "model"]
    if not models:
        pytest.skip("no model asset recorded")
    for a in models:
        assert "eleven" not in a["description"].lower()
        assert "12 heads" in a["description"], a["description"]


# --------------------------------- canonical artifacts, not the presentation layer
def _sidecar(name: str = "artifacts/quadcond_model.json") -> dict:
    import json

    import pytest

    p = Path(name)
    if not p.exists():
        pytest.skip(f"{name} not present")
    return json.loads(p.read_text())


def test_the_serialized_sidecar_carries_the_corrected_buffer_count():
    """Read the ARTIFACT, not a freshly opened atlas.

    The report and the README were regenerated with 261 while both
    `quadcond_model.json` sidecars — the files the model card is built from and
    the API serves — still said 255, and the test that was meant to catch this
    opened `atlas.db` and compared it against itself. A number is corrected when
    the artifact carries it, not when the page does.
    """
    for name in ("artifacts/quadcond_model.json", "artifacts/quadcond_model_seqonly.json"):
        m = _sidecar(name)
        rows = (m.get("atlas_snapshot") or {}).get("summary") or []
        g4 = [r for r in rows if r.get("source") == "g4stab_experimental_tm"]
        if not g4:
            continue
        assert g4[0]["n_conditions"] == 261, (
            f"{name} still records {g4[0]['n_conditions']} buffers for g4stab")


def test_the_model_card_has_no_stale_buffer_count():
    import pytest

    p = Path("docs/MODEL_CARD.md")
    if not p.exists():
        pytest.skip("model card not present")
    text = p.read_text()
    assert '"n_conditions": 255' not in text
    assert "locus_peak_overlap_state" in text, "the card must describe all 12 heads"


def test_the_served_applicability_marks_tm_temperature_as_not_enforced():
    """`/info` serves the applicability table verbatim.

    `_domain_check` skips the query temperature for a Tm-target head, but the
    field stayed in the table the API publishes — so the served contract
    described a gate the code does not apply. Fixing the report's rendering did
    not fix what the API says.
    """
    import pytest

    from quadcond import assets
    from quadcond.models.base import MultiTaskModel

    try:
        path = assets.resolve("model")
    except assets.AssetError:
        pytest.skip("model artifact not available in this environment")
    m = MultiTaskModel.load(path)

    for name, head in m.heads.items():
        ap = head.applicability or {}
        if "temperature" not in ap:
            continue
        if head.target == "tm":
            assert ap.get("temperature_is_enforced") is False, (
                f"{name} publishes a temperature range with no note that it is not a gate")
            assert "not a domain gate" in ap.get("temperature_note", "").lower()
        else:
            # a pH_T head genuinely consumes temperature; it must stay enforced
            assert ap.get("temperature_is_enforced", True) is True, name


def test_info_and_domain_check_agree_about_every_gate():
    """Whatever /info publishes as a range must be what the code enforces."""
    import pytest

    from quadcond import assets, service
    from quadcond.conditions import Condition
    from quadcond.models.predict import Predictor, _domain_check

    try:
        path = assets.resolve("model")
    except assets.AssetError:
        pytest.skip("model artifact not available in this environment")
    pred = Predictor.load(path)

    # A query at the far edge of each published range must produce a warning
    # naming that field -- unless the table says the field is not enforced.
    head = pred.model.heads["g4_tm"]
    ap = head.applicability
    cond = Condition.from_mapping({
        "k": ap["k"]["max"] + 1000, "na": 0.0, "ph": 7.0, "temperature": 37.0,
    })
    out = _domain_check(head, cond, "AGGGTTAGGGTTAGGGTTAGGG")
    joined = " ".join(out["warnings"])
    assert "k" in joined, "a query past the published K+ range must be flagged"
    assert "temperature" not in joined, (
        "temperature is published as not enforced for a Tm head; it must not be flagged")

    del service  # imported to assert the module loads alongside these checks


# --------------------------------------------------------------------------- #
# v0.4.7: capabilities, scans, and the synthetic quarantine
# --------------------------------------------------------------------------- #

TEL22 = "AGGGTTAGGGTTAGGGTTAGGG"
ROOT = Path(__file__).resolve().parent.parent


def _load_model():
    """The trained artifact, or a skip. Same contract as the tests above."""
    from quadcond import assets
    from quadcond.models.base import MultiTaskModel
    try:
        return MultiTaskModel.load(assets.resolve("model"))
    except assets.AssetError:
        pytest.skip("model artifact not available in this environment")


def _load_predictor():
    from quadcond import assets
    from quadcond.models.predict import Predictor
    try:
        return Predictor.load(assets.resolve("model"))
    except assets.AssetError:
        pytest.skip("model artifact not available in this environment")
# --------------------------------------------------------------------------- #

def test_an_axis_that_never_varied_is_an_absence_not_a_range():
    model = _load_model()
    """The distinction the six sliders were papering over.

    `im_pht` is 160 measurements in one buffer. Asking it about pH is not an
    extrapolation from a range -- there is no range, and the interface must not
    imply that moving the control moves the answer.
    """
    from quadcond import capabilities as cap
    axes = cap.axis_response(model.heads["im_pht"])
    assert axes["k"]["state"] == "fixed"
    assert "no response to apply" in axes["k"]["note"]
    assert cap.responds_to(model.heads["im_pht"]) == []


def test_a_tm_head_predicts_temperature_rather_than_consuming_it():
    model = _load_model()
    from quadcond import capabilities as cap
    axes = cap.axis_response(model.heads["g4_tm"])
    assert axes["temperature"]["state"] == "predicts"
    assert "does not take one" in axes["temperature"]["note"]
    # And it is neither offered as a control nor reported as an absence.
    caps = cap.head_capabilities(model.heads["g4_tm"])
    assert "temperature" not in caps["responds_to"]
    assert "temperature" not in caps["inert_axes"]


def test_semantics_decides_the_workspace_not_the_condition_axis():
    """A genomic proxy with a condition column is still a genomic proxy.

    An earlier version filed condition-responsive heads under Condition
    response *instead of* Sequence evidence, which emptied the latter of
    `g4_tm` -- the best-attested head in the model.
    """
    model = _load_model()
    from quadcond import capabilities as cap
    groups = {w["id"]: w["heads"] for w in cap.workspaces(model.heads)}
    assert "g4_tm" in groups["sequence_evidence"]
    assert "g4_tm" in groups["condition_response"]
    assert "g4_tm_distilled" in groups["diagnostics"]
    for name in ("g4_fold_genomic", "im_fold_genomic", "locus_peak_overlap_state"):
        assert name in groups["genomic_context"]
        assert name not in groups["sequence_evidence"]


def test_a_query_of_no_nucleotides_is_refused_not_extrapolated():
    predictor = _load_predictor()
    """`clean("ZZZZ")` is `"NNNN"` -- non-empty, and it featurized happily.

    It returned a confident 25 degC melting temperature. In a batch of 200
    FASTA records that number is never eyeballed and is indistinguishable from
    a real one.
    """
    from quadcond.conditions import Condition
    e = predictor.predict("ZZZZ", Condition.preset("physiological"),
                          heads=["g4_tm"], n_neighbours=0)[0]["predictions"]["g4_tm"]
    assert e["refused"] is True
    # No value key at all: a consumer that forgets the flag gets a KeyError,
    # which is the correct failure.
    assert "value" not in e
    assert "no recognisable nucleotides" in e["refusal_reason"]


def test_partial_ambiguity_warns_but_still_answers():
    predictor = _load_predictor()
    """One N in twenty-two is a real query whose features were masked."""
    from quadcond.conditions import Condition
    e = predictor.predict("AGGGTTAGGNTTAGGGTTAGGG", Condition.preset("physiological"),
                          heads=["g4_tm"], n_neighbours=0)[0]["predictions"]["g4_tm"]
    assert not e.get("refused")
    assert "value" in e
    assert any("positions are N" in w for w in e["applicability"]["warnings"])


def test_the_mutation_scan_asks_each_head_about_its_own_strand():
    predictor = _load_predictor()
    """The category error, at a new endpoint.

    A G4 and an i-motif at one duplex position sit on opposite strands. A first
    pass at this scan mutated only the typed strand and asked three i-motif
    heads about Tel22's G-rich one.
    """
    from quadcond import scans
    from quadcond.conditions import Condition
    r = scans.mutation_scan(predictor, TEL22, Condition.preset("physiological"),
                            heads=["g4_tm", "im_pht"], positions=[0])
    assert r["wild_type"]["strand_by_head"] == {"g4_tm": "+", "im_pht": "-"}
    row = r["substitutions"][0]
    assert row["heads"]["g4_tm"]["strand"] == "+"
    assert row["heads"]["im_pht"]["strand"] == "-"
    # One duplex edit, propagated to both strands.
    assert row["complement"] == revcomp(row["sequence"])
    assert r["wild_type"]["strand_sequence"] == {"+": TEL22, "-": revcomp(TEL22)}


def test_a_mutation_that_destroys_the_motif_yields_a_state_and_no_number():
    predictor = _load_predictor()
    from quadcond import scans
    from quadcond.conditions import Condition
    r = scans.mutation_scan(predictor, TEL22, Condition.preset("physiological"),
                            heads=["g4_tm"], positions=[1])
    cells = [row["heads"]["g4_tm"] for row in r["substitutions"]]
    assert cells, "position 1 is a G in the first tract"
    for cell in cells:
        assert cell["motif"]["state"] == "motif_lost"
        # Not a zero, not a null: no key.
        assert "delta" not in cell


def test_the_delta_uncertainty_is_paired_and_says_it_is_not_validated():
    predictor = _load_predictor()
    """Not the marginal half-widths, and not sold as a validated interval.

    Wild type and mutant differ by one base, so the ensemble members' errors are
    strongly correlated and mostly cancel. Propagating the marginal intervals
    would report a difference far noisier than it is.
    """
    from quadcond import scans
    from quadcond.conditions import Condition
    r = scans.mutation_scan(predictor, TEL22, Condition.preset("physiological"),
                            heads=["g4_tm"], positions=[0])
    cell = next(row["heads"]["g4_tm"] for row in r["substitutions"]
                if "delta" in row["heads"]["g4_tm"])
    assert "delta_paired_sd" in cell and cell["n_estimators"] > 1
    # Said once for the scan, not stamped onto all 330 cells.
    assert "not a validated interval" in r["delta_note"]
    assert "delta_note" not in cell
    # The delta is the difference of the reported values, up to rounding.
    assert abs((cell["mutant_value"] - cell["wild_type_value"]) - cell["delta"]) < 0.05
    assert "No significance test is applied" in r["ranking_note"]


def test_the_scan_freezes_its_baseline_into_the_result():
    predictor = _load_predictor()
    """The panel this replaces recomputed the displayed wild type on selection
    while keeping the earlier scan's rows, so the table and its stated reference
    could disagree without saying so."""
    from quadcond import scans
    from quadcond.conditions import Condition
    cond = Condition.preset("physiological")
    r = scans.mutation_scan(predictor, TEL22, cond, heads=["g4_tm"], positions=[0])
    assert r["wild_type"]["sequence"] == TEL22
    assert r["wild_type"]["condition_detail"] == cond.to_dict()
    assert r["wild_type"]["predictions"]["g4_tm"]["value"] is not None


def test_condition_scan_names_inert_heads_rather_than_dropping_them():
    predictor = _load_predictor()
    """A head missing from a chart reads as 'not applicable'. A head listed as
    inert says which of two very different reasons applies."""
    from quadcond import scans
    from quadcond.conditions import Condition
    r = scans.condition_scan(predictor, TEL22, "k", [0, 50, 100, 150],
                             Condition.preset("physiological"), n_neighbours=0)
    inert = {h["head"] for h in r["inert_heads"]}
    # Every drawn series names the strand it was asked about, and the two
    # strand strings appear once at the top level rather than on every point.
    assert set(r["strand_sequence"]) == {"+", "-"}
    for name, series in r["series"].items():
        assert series["strand"] in ("+", "-"), name
        assert series["strand_sequence"] in r["strand_sequence"].values(), name
    assert "im_pht" in inert and "g4_topology" in inert
    assert all(h["reason"] for h in r["inert_heads"])
    assert "im_pht" not in r["series"]


def test_condition_scan_measures_whether_the_prediction_actually_moved():
    predictor = _load_predictor()
    """'The axis varied in training' licenses a dependence; it does not show one.

    A composition-matched shuffle inherits its parent row's buffer, so a
    classification head shows wide training variation in an axis its label never
    depended on. `g4_fold` is that case: its probability moves by ~1e-3 across
    the whole potassium range.
    """
    from quadcond import scans
    from quadcond.conditions import Condition
    r = scans.condition_scan(predictor, TEL22, "k", [0, 25, 50, 100, 150],
                             Condition.preset("physiological"), n_neighbours=0)
    tm = r["series"]["g4_tm"]["observed_response"]
    assert tm["verdict"] == "responds"
    assert tm["spread"] > 10, "K+ 0->150 mM is a large, measured Tm effect"

    fold = r["series"]["g4_fold"]["observed_response"]
    assert fold["quantity"] == "probability"
    assert fold["spread"] < 0.01, "this head has no real potassium response"
    # The verdict word matters. "flat" reads as a measured insensitivity to
    # potassium, which is a claim nobody made; the comparison behind it is a
    # response range set beside a marginal residual half-width, which is
    # descriptive and not a test.
    assert fold["verdict"] in ("below_error_scale", "unknown_scale")
    assert "no resolvable response" not in fold.get("note", "")
    # Summaries cover the supported part of the range only.
    assert fold["summary_basis"] == "in_domain_points_only"


def test_condition_scan_refuses_an_axis_no_head_could_answer():
    predictor = _load_predictor()
    from quadcond import scans
    with pytest.raises(ValueError, match="unknown condition axis"):
        scans.condition_scan(predictor, TEL22, "ionic_strength", [1, 2])


def test_batch_preserves_identifiers_and_reports_exclusions():
    predictor = _load_predictor()
    from quadcond import scans
    from quadcond.conditions import Condition
    fasta = ">tel22\nAGGGTTAGGGTTAGGGTTAGGG\n>empty\n\n>myc\nGGGGAGGGTGGGGAGGGTGGGG"
    r = scans.batch_predict(predictor, fasta, [Condition.preset("physiological")],
                            heads=["g4_tm"])
    ids = [row["id"] for row in r["results"]]
    assert ids == ["tel22", "myc"], "identifiers survive verbatim"
    # Three headers in, three records out. This assertion used to read
    # `n_sequences_in == 2`, which pinned the defect in place instead of
    # catching it: the parser created no record at all for a header with no
    # sequence, so an input the user supplied vanished with an empty
    # `excluded` list beside it. A batch of 200 that returns 197 rows with no
    # explanation is a batch whose three missing sequences get noticed after
    # the figure is drawn, if at all.
    assert r["run_record"]["n_sequences_in"] == 3
    assert r["run_record"]["n_sequences_scored"] == 2
    assert [e["id"] for e in r["excluded"]] == ["empty"]
    assert r["run_record"]["quadcond_model_version"]
    assert r["run_record"]["model_artifact_sha256"] is not None
    assert r["run_record"]["prediction_schema_version"]


def test_the_run_record_carries_what_a_figure_needs_six_months_later():
    predictor = _load_predictor()
    from quadcond import scans
    from quadcond.conditions import Condition
    r = scans.batch_predict(predictor, [{"id": "a", "sequence": TEL22}],
                            [Condition.preset("physiological")], heads=["g4_tm"])
    rec = r["run_record"]
    for key in ("quadcond_model_version", "prediction_schema_version",
                "heads_requested", "n_rows", "conditions",
                "condition_imputed_fields"):
        assert key in rec, key


def test_a_synthetic_tier_can_never_report_as_measurement_grounded():
    """The branch generator registered its rows as `evidence_tier="experimental"`.

    Every claim-basis field downstream would then have certified a fit to a
    generated CSV as measurement-grounded, and the model card, /info and the
    viewer's badges would all have agreed with each other.
    """
    from quadcond import claims
    meta = {"tiers": ("synthetic", "experimental"), "sources": ("synthetic_dataset",)}
    sem = claims.semantics("g4_fold", meta, "binary")
    assert sem["target_semantics"] == claims.SYNTHETIC
    assert sem["biophysically_grounded"] is False
    assert sem["has_experimental_observation"] is False


def test_the_headline_leads_with_the_synthetic_warning():
    """A reader who stops after the first clause has to have been told."""
    from quadcond import claims

    class _H:
        def __init__(self, name, tiers):
            self.name = name
            self.training_meta = {"tiers": tiers, "sources": ()}

    line = claims.headline({"a": _H("a", ("synthetic",)),
                            "b": _H("b", ("experimental",))})
    assert line.startswith("SYNTHETIC MODEL")
    assert "Not for scientific use" in line


def test_the_synthetic_builder_refuses_the_canonical_artifact_path():
    """The one path it must never produce."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "synth", ROOT / "scripts" / "18_synthetic_dev_model.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    for bad in ("artifacts/quadcond_model.joblib", "artifacts/my_model.joblib"):
        with pytest.raises(SystemExit):
            mod.checked_output(bad)
    assert mod.checked_output("artifacts/synthetic_dev_model.joblib").name \
        == "synthetic_dev_model.joblib"
