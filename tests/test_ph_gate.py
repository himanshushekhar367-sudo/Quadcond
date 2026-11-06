"""A transitional-pH head predicts pH; its imputed pH=7 must not gate it.

Every pH_T record in the atlas stores ph=7.0 and lists "ph" in
condition_imputed: a titration sweeps pH, so 7.0 is a placeholder. Treating it
as a fixed training condition hid the i-motif card at every pH below 7 -- where
i-motifs fold. These tests pin the exemption, keep the gate on every other
axis and head, and check the "use this head's buffer" suggestion.
"""
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from quadcond.capabilities import axis_response
from quadcond.conditions import Condition
from quadcond.models.predict import _domain_check

C_STRAND = "CCCTAACCCTAACCCTAACCCT"
FIXED = {
    "ph": {"min": 7.0, "max": 7.0, "n_unique": 1},
    "k": {"min": 100.0, "max": 100.0, "n_unique": 1},
    "na": {"min": 10.0, "max": 10.0, "n_unique": 1},
    "mg": {"min": 0.0, "max": 0.0, "n_unique": 1},
    "temperature": {"min": 25.0, "max": 25.0, "n_unique": 1},
}


def _head(target: str, applicability=None):
    return SimpleNamespace(
        name="im_pht" if target == "ph_t" else "im_tm_probe", kind="iM", target=target,
        task="regression", applicability=applicability or FIXED, training_meta={},
    )


def _cond(**kw) -> Condition:
    base = dict(k=100, na=10, li_nh4=0, mg=0, ph=7.0, temperature=25, crowder_pct=0, strand_conc=5)
    base.update(kw)
    return Condition(**base)


@pytest.mark.parametrize("ph", [7.0, 6.4, 6.0, 5.5, 4.5])
def test_ph_t_head_is_not_gated_on_its_imputed_ph(ph):
    out = _domain_check(_head("ph_t"), _cond(ph=ph), C_STRAND)
    assert out["in_domain"], out["warnings"]
    assert not any(w.startswith("ph") for w in out["warnings"])
    assert "ph" not in out["suggested_condition"]


def test_other_heads_are_still_gated_on_a_fixed_ph():
    out = _domain_check(_head("tm"), _cond(ph=5.5), C_STRAND)
    assert not out["in_domain"]
    assert any(w.startswith("ph:") for w in out["warnings"])
    assert out["suggested_condition"] == {"ph": 7.0}


def test_ph_t_head_is_still_gated_on_every_other_axis():
    out = _domain_check(_head("ph_t"), _cond(k=140, mg=1, temperature=37, ph=5.5), C_STRAND)
    assert not out["in_domain"]
    assert out["suggested_condition"] == {"k": 100.0, "mg": 0.0, "temperature": 25.0}


def test_suggestion_clamps_into_a_varied_range_and_lists_only_axes_to_move():
    ap = dict(FIXED, k={"min": 10.0, "max": 500.0, "n_unique": 9})
    out = _domain_check(_head("ph_t", ap), _cond(k=1000), C_STRAND)
    assert out["suggested_condition"] == {"k": 500.0}
    inside = _domain_check(_head("ph_t", ap), _cond(k=300), C_STRAND)
    assert inside["in_domain"] and inside["suggested_condition"] == {}


def test_capabilities_report_ph_as_the_answer_of_a_ph_t_head():
    assert axis_response(_head("ph_t"))["ph"]["state"] == "predicts"
    assert axis_response(_head("tm"))["ph"]["state"] == "fixed"
    assert axis_response(_head("tm"))["temperature"]["state"] == "predicts"


MODEL = os.environ.get("QUADCOND_MODEL")


@pytest.mark.skipif(not (MODEL and Path(MODEL).exists()), reason="release model not available")
def test_imotif_card_is_shown_and_folds_at_acidic_ph():
    from quadcond import service

    def card(ph):
        r = service.evidence_payload({"sequence": C_STRAND, "k": 100, "na": 10, "li_nh4": 0, "mg": 0,
                                      "ph": ph, "temperature": 25, "crowder_pct": 0, "strand_conc": 5})
        return next(c for c in r["evidence"] if c["kind"] == "i-motif")

    neutral, acidic, basic = card(7.0), card(5.5), card(8.0)
    assert neutral["applicability"]["in_domain"] and acidic["applicability"]["in_domain"]
    # pH_T is read at the head's own reference pH, so it cannot drift with the
    # pH the user asks about; only the folded fraction moves.
    assert neutral["predictedPhT"] == acidic["predictedPhT"] == basic["predictedPhT"]
    assert neutral["predictedPhT_interval"] == basic["predictedPhT_interval"]
    assert neutral["foldedFraction"] < 0.1 < 0.9 < acidic["foldedFraction"]
    assert "n = 2.5" in acidic["foldedFraction_note"]
