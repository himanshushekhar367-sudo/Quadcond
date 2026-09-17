"""The genomic path: G4-seq scanner, and refusal of oligo folding heads on long windows."""
from pathlib import Path

import pytest

from quadcond import assets
from quadcond.conditions import Condition

ROOT = Path(__file__).resolve().parents[1]
TEL = "GGGTTAGGGTTAGGGTTAGGG"
MYC = "TGGGGAGGGTGGGGAGGGTGGGGAAGG"


def _scanner():
    a = assets.assets().get("g4seq_scanner")
    if a is None:
        pytest.skip("no g4seq_scanner in manifest")
    try:
        path = assets.resolve_asset(a)
    except assets.AssetError:
        pytest.skip("g4seq_scanner asset not present")
    from quadcond.genome_scan import G4SeqScanner
    return G4SeqScanner.load(path)


def test_manifest_records_the_scanner_as_optional():
    a = assets.assets()["g4seq_scanner"]
    assert a.kind == "scanner" and a.required is False and a.sha256
    assert a.url and a.url.endswith(a.filename)


def test_every_published_asset_has_a_url_except_the_training_atlas():
    for a in assets.assets().values():
        if a.name == "atlas":
            assert a.url is None
        else:
            assert a.url and "/releases/download/" in a.url


def test_window_features_are_strand_symmetric():
    from quadcond.genome_scan import window_features
    from quadcond.motifs import revcomp
    s = "ACGTTGCA" * 5 + MYC + "TTACAGGA" * 10
    fw, rv = window_features([s, revcomp(s)])
    assert (abs(fw - rv) < 1e-5).all()


def test_scanner_separates_a_promoter_g4_from_background():
    sc = _scanner()
    bg = ("ATTCAGCTTAACGATTGCAATCTGATCAGTTACAGATCA" * 4)[:124]
    myc = bg[:50] + MYC + bg[50 + len(MYC):]
    p = sc.score_windows([bg, myc])
    assert p[1] > 0.5 > p[0]


def test_isolated_canonical_motif_is_rescued_and_labelled():
    sc = _scanner()
    seq = "ACGT" * 40 + "TT" + TEL + "AA" + "CATG" * 40
    regions = sc.scan(seq)
    hits = [r for r in regions if any(m["sequence"] == TEL for m in r["motifs"])]
    assert hits and hits[0]["call_basis"] in {"g4seq_scanner", "canonical_motif_rescue"}


def test_folding_head_refuses_a_genomic_window():
    from quadcond.models.predict import Predictor
    try:
        pred = Predictor.load(assets.resolve("model"), None)
    except assets.AssetError:
        pytest.skip("model asset not present")
    long = ("ACGT" * 31)
    out = pred.predict([long, TEL], Condition(k=100), heads=["g4_fold", "g4_tm"])
    g = out[0]["predictions"]["g4_fold"]
    assert g["applicability"]["refused"] and "genome-scan" in g["applicability"]["refusal_reason"]
    assert not out[1]["predictions"]["g4_fold"]["applicability"]["refused"]
    assert not out[0]["predictions"]["g4_tm"]["applicability"]["refused"]
