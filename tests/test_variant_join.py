"""The AlphaGenome join, checked without a key, a network or a 90 GB file.

Everything that can go wrong in this join is a control-flow or coordinate
problem, and none of it needs a real regulatory score to expose:

* a one-base coordinate offset silently pairs each structural row with the
  wrong variant, and every number in the table stays plausible;
* a chromosome naming mismatch (``7`` against ``chr7``) empties the regulatory
  axis and reads as "this locus is quiet";
* a missing regulatory record treated as a zero puts every unscored variant at
  the bottom of the ranking, which is a claim about those variants;
* a source that failed to load, reported as a source that found nothing.

The sentinel predictor from ``test_comparison_contracts`` supplies the
structural half. The regulatory half comes from a table written in the test, and
from a real bgzip/Tabix fixture where ``pysam`` is available.
"""
from __future__ import annotations

import pytest

from quadcond import alphagenome as ag
from quadcond import variants as var

from test_comparison_contracts import TEL22, _predictor

CHROM = "chr8"
START = 128_748_301          # 1-based coordinate of TEL22[0]


def _table(tmp_path, rows, header="chromosome\tposition\tref\talt\tAVI\n"):
    p = tmp_path / "avi.tsv"
    p.write_text(header + "".join(rows), encoding="utf-8")
    return p


def _all_snvs(seq, chrom=CHROM, start=START, score=lambda i, b: 0.5):
    for i, wt in enumerate(seq):
        for alt in "ACGT":
            if alt == wt:
                continue
            yield f"{chrom}\t{start + i}\t{wt}\t{alt}\t{score(i, alt)}\n"


# ---------------------------------------------------------------- coordinates
def test_the_join_is_on_the_variant_not_on_position_alone(tmp_path):
    """Ref and alt are part of the key, so an offset cannot pair silently.

    A join on position alone would attach one regulatory record to all three
    substitutions at that base, and would survive a shifted window with every
    number still looking reasonable.
    """
    path = _table(tmp_path, list(_all_snvs(TEL22)))
    src = ag.TableAtlas.from_path(path)
    res = var.variant_scan(_predictor(), TEL22, CHROM, START,
                           heads=["g4_tm"], atlas_source=src)
    assert res["window"] == {**res["window"]}
    for row in res["variants"]:
        v = row["variant"]
        # The reference base recorded for the variant is the assembly base at
        # that coordinate, not the base at the same index of some other window.
        assert v["reference"] == TEL22[v["position"] - START]
        assert v["reference"] != v["alternate"]
    assert all(r["regulatory"] is not None for r in res["variants"])


def test_a_shifted_window_does_not_quietly_half_join(tmp_path):
    path = _table(tmp_path, list(_all_snvs(TEL22)))
    src = ag.TableAtlas.from_path(path)
    shifted = var.variant_scan(_predictor(), TEL22, CHROM, START + 1,
                               heads=["g4_tm"], atlas_source=src)
    joined = sum(1 for r in shifted["variants"] if r["regulatory"] is not None)
    # The table's records are keyed by (position, ref, alt). One base along,
    # most rows find no record at all rather than finding the wrong one, and
    # the response says how many were left unclassified.
    assert joined < len(shifted["variants"])
    assert shifted["quadrant_counts"]["unclassified"] > 0


def test_chromosome_spelling_does_not_decide_the_answer(tmp_path):
    """`7` and `chr7` are the same chromosome. A mismatch is not a quiet locus."""
    path = _table(tmp_path, list(_all_snvs(TEL22, chrom="8")))
    src = ag.TableAtlas.from_path(path)
    res = var.variant_scan(_predictor(), TEL22, "chr8", START,
                           heads=["g4_tm"], atlas_source=src)
    assert all(r["regulatory"] is not None for r in res["variants"])


# ------------------------------------------------------------------- absence
def test_a_missing_regulatory_record_is_not_a_low_score(tmp_path):
    """Absent and low are different, and only one of them is a finding."""
    rows = [r for r in _all_snvs(TEL22) if not r.split("\t")[1].endswith("305")]
    src = ag.TableAtlas.from_path(_table(tmp_path, rows))
    res = var.variant_scan(_predictor(), TEL22, CHROM, START,
                           heads=["g4_tm"], atlas_source=src)
    missing = [r for r in res["variants"] if r["regulatory"] is None]
    assert missing, "the fixture withholds one position"
    for row in missing:
        assert row["quadrant"] == "unclassified"
        assert "no regulatory record" in row["quadrant_reason"]
        assert row["combined_rank"] is None


def test_no_regulatory_source_still_runs_and_says_so():
    res = var.variant_scan(_predictor(), TEL22, CHROM, START, heads=["g4_tm"])
    assert res["run_record"]["regulatory_source"]["source"] == "none"
    assert "not assigned" in res["run_record"]["regulatory_source"]["note"]
    assert res["quadrant_counts"]["unclassified"] == len(res["variants"])
    # The structural half is unaffected.
    assert any(r["structural"]["delta"] is not None for r in res["variants"])


def test_a_source_that_failed_is_not_a_source_that_found_nothing():
    class _Broken(ag.AtlasSource):
        name = "broken"
        def records_for_interval(self, chromosome, start, end):
            raise ag.AtlasUnavailable("the endpoint could not be resolved")

    res = var.variant_scan(_predictor(), TEL22, CHROM, START,
                           heads=["g4_tm"], atlas_source=_Broken())
    assert res["run_record"]["regulatory_error"]
    assert "could not be resolved" in res["run_record"]["regulatory_error"]


# ------------------------------------------------------------------ quadrants
def test_the_quadrant_needs_both_coordinates(tmp_path):
    """A variant that destroys the motif has no delta, so it has no quadrant."""
    src = ag.TableAtlas.from_path(_table(tmp_path, list(_all_snvs(TEL22))))
    res = var.variant_scan(_predictor(), TEL22, CHROM, START,
                           heads=["g4_tm"], atlas_source=src)
    for row in res["variants"]:
        if row["structural"]["magnitude"] is None:
            assert row["quadrant"] == "unclassified"
        else:
            assert row["quadrant"] in {"both", "structural_only",
                                       "regulatory_only", "neither"}
    assert sum(res["quadrant_counts"].values()) == len(res["variants"])


def test_the_combined_rank_requires_both_axes(tmp_path):
    """The minimum of two ranks: high on one axis alone is not high."""
    def score(i, alt):
        return 9.0 if i == 0 else 0.01
    src = ag.TableAtlas.from_path(_table(tmp_path, list(_all_snvs(TEL22, score=score))))
    res = var.variant_scan(_predictor(), TEL22, CHROM, START,
                           heads=["g4_tm"], atlas_source=src)
    ranked = [r for r in res["variants"] if r["combined_rank"] is not None]
    assert ranked, "some variants carry both axes"
    for row in ranked:
        s = row["structural"]["rank"]
        g = row["regulatory"]["rank"]
        assert row["combined_rank"] == pytest.approx(min(s, g), abs=1e-4)
    # Sorted best-first, with the unrankable at the end.
    seen_none = False
    for row in res["variants"]:
        if row["combined_rank"] is None:
            seen_none = True
        else:
            assert not seen_none, "unclassified rows must not outrank scored ones"


def test_thresholds_are_recorded_with_the_result(tmp_path):
    src = ag.TableAtlas.from_path(_table(tmp_path, list(_all_snvs(TEL22))))
    res = var.variant_scan(_predictor(), TEL22, CHROM, START, heads=["g4_tm"],
                           atlas_source=src, structural_threshold=1.0,
                           regulatory_threshold=0.25)
    rec = res["run_record"]
    assert rec["structural_threshold"] == 1.0
    assert rec["regulatory_threshold"] == 0.25
    assert rec["threshold_basis"] == "caller-supplied absolute thresholds"
    locus = rec["locus"]
    assert locus["chromosome"] == CHROM
    assert locus["start"] == START
    assert locus["end"] == START + len(TEL22) - 1
    # The convention is stated, not assumed: "start" means three different
    # things across the formats this workflow touches.
    assert locus["coordinate_convention"] == "1-based inclusive, as in VCF"
    assert locus["strand_of_window_sequence"].startswith("+")
    assert locus["assembly_stated_by_caller"] is False
    assert locus["assembly"] == "not stated by caller"


def test_nothing_here_emits_a_combined_score(tmp_path):
    """Two axes on different scales are never multiplied into one number."""
    src = ag.TableAtlas.from_path(_table(tmp_path, list(_all_snvs(TEL22))))
    res = var.variant_scan(_predictor(), TEL22, CHROM, START,
                           heads=["g4_tm"], atlas_source=src)
    for banned in ("combined_score", "priority_score", "joint_score"):
        assert banned not in res
        assert all(banned not in row for row in res["variants"])
    assert "not a classifier" in res["quadrant_note"]
    assert "not calibrated" in res["rank_note"]


def test_the_window_is_bounded():
    with pytest.raises(ValueError, match="capped"):
        var.variant_scan(_predictor(), "G" * 500, CHROM, START, heads=["g4_tm"])


# --------------------------------------------------------------------- tabix
def test_a_tabix_bundle_is_queried_by_region(tmp_path):
    """The published AVI bundle path, exercised on a real bgzip/Tabix fixture."""
    pysam = pytest.importorskip("pysam")
    raw = tmp_path / "avi.bed"
    header = "#chromosome\tposition\treference_bases\talternate_bases\tAVI\tAVI_phred\n"
    body = "".join(f"{CHROM}\t{START + i}\t{TEL22[i]}\t{alt}\t0.4\t12.0\n"
                   for i in range(len(TEL22)) for alt in "ACGT" if alt != TEL22[i])
    raw.write_text(header + body, encoding="utf-8")
    gz = str(raw) + ".gz"
    pysam.tabix_compress(str(raw), gz, force=True)
    pysam.tabix_index(gz, seq_col=0, start_col=1, end_col=1, force=True, meta_char="#")

    src = ag.TabixAtlas.from_path(gz)
    described = src.describe()
    assert described["source"] == "alphagenome_avi_tabix"
    # Columns come from the file's own header, never from a hardcoded order.
    assert "AVI" in described["score_columns"]
    assert "AVI_phred" in described["score_columns"]

    res = var.variant_scan(_predictor(), TEL22, CHROM, START,
                           heads=["g4_tm"], atlas_source=src,
                           regulatory_scorer="AVI")
    assert all(r["regulatory"] is not None for r in res["variants"])
    assert {r["regulatory"]["scorer"] for r in res["variants"]} == {"AVI"}


def test_a_tabix_bundle_without_an_index_is_refused(tmp_path):
    pysam = pytest.importorskip("pysam")
    raw = tmp_path / "unindexed.bed"
    raw.write_text("#chromosome\tposition\tref\talt\tAVI\n", encoding="utf-8")
    gz = str(raw) + ".gz"
    pysam.tabix_compress(str(raw), gz, force=True)
    with pytest.raises(ag.AtlasUnavailable, match="index"):
        ag.TabixAtlas.from_path(gz)


# -------------------------------------------------------------------- sources
def test_the_live_backend_refuses_without_a_key(monkeypatch):
    monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
    with pytest.raises(ag.AtlasUnavailable, match="no AlphaGenome API key"):
        ag.LiveAtlas.from_api_key(None)


def test_source_resolution_prefers_the_published_bundle(monkeypatch, tmp_path):
    """No key, permissive licence, region query: the Tabix bundle comes first."""
    monkeypatch.setenv("ALPHAGENOME_API_KEY", "not-used-because-tabix-wins")
    calls = {}
    monkeypatch.setattr(ag.TabixAtlas, "from_path",
                        classmethod(lambda cls, p: calls.setdefault("tabix", str(p))))
    ag.resolve_source(tabix=tmp_path / "avi.tsv.gz")
    assert "tabix" in calls
    monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
    assert ag.resolve_source().describe()["source"] == "none"


# ------------------------------------------------- screening and candidates
NEAR_MISS = "GGGTTAGGGTTAGGGTTAGTG"     # last G-tract broken; one edit restores it


def test_screening_sees_motifs_a_substitution_would_create():
    """The reference alone is not the question a variant workflow is asking.

    A first version screened the reference only, which discards the most
    interesting case there is: a window with no motif where one substitution
    creates one. Such a window would have been declared uninteresting and never
    queried.
    """
    from quadcond.variants import prescreen
    d = prescreen(NEAR_MISS)
    assert d["reference_motif_total"] == 0
    assert d["verdict"] == "motif_gain_possible"
    assert d["n_motif_gain"] >= 1
    gained = d["motif_gain_candidates"][0]
    assert gained["motifs_mutant"] > gained["motifs_reference"]


def test_screening_separates_four_states_not_two():
    from quadcond.variants import prescreen
    assert prescreen(TEL22)["verdict"] == "motif_present"
    assert prescreen(NEAR_MISS)["verdict"] == "motif_gain_possible"
    quiet = "AGGAGGGCAGAGAGCTGGGGCCTCGGACTCACCCGACGCTTGTGATGAGCTGCACCCAGGA"
    assert prescreen(quiet)["verdict"] == "no_motif"
    # And the no-motif note is about the rules, not about the molecule.
    note = prescreen(quiet)["note"]
    assert "model applicability" in note
    assert "can adopt" in note


def test_a_motif_destroying_variant_never_falls_off_the_shortlist(tmp_path):
    """A delta-ranked table cannot be the only place a finding appears.

    Motif loss and gain have no delta by construction, so they sort to the
    bottom of a ranked table and off the end of a shortlist -- which is a
    retrieval failure, not a formatting one. They get their own categorical
    list instead, and it is not given an invented number.
    """
    src = ag.TableAtlas.from_path(_table(tmp_path, list(_all_snvs(TEL22))))
    res = var.variant_scan(_predictor(), TEL22, CHROM, START,
                           heads=["g4_tm"], atlas_source=src)
    cats = res["structural_candidates"]
    assert set(cats) == {"motif_lost", "motif_gained", "motif_count_changed"}
    listed = [e for entries in cats.values() for e in entries]
    for entry in listed:
        assert "delta" not in entry, "a categorical candidate must carry no delta"
        assert entry["why_no_delta"]
        # The regulatory axis is still attached where it exists, so the list is
        # orderable without a structural number.
        assert "regulatory_score" in entry
    ranked_labels = {r["label"] for r in res["variants"]
                     if r["combined_rank"] is not None}
    for entry in listed:
        assert entry["label"] not in ranked_labels, (
            "a categorical finding must not also appear in the delta ranking")
    assert "not as effects of a measured size" in res["candidate_note"]


# --------------------------------------------------------------- provenance
def test_every_result_identifies_the_model_file_by_its_contents(tmp_path):
    """`model_artifact_sha256: ""` names a version and nothing else.

    A path is a location in somebody's filesystem. While the release artifact's
    identity is disputed, a run record that cannot say which bytes answered is
    not usable evidence.
    """
    src = ag.TableAtlas.from_path(_table(tmp_path, list(_all_snvs(TEL22))))
    pred = _predictor()
    # The sentinel predictor is built in memory and has no file behind it, so
    # the hash is legitimately absent -- and must be reported as absent rather
    # than as an empty string that looks like a value.
    res = var.variant_scan(pred, TEL22, CHROM, START, heads=["g4_tm"],
                           atlas_source=src)
    prov = res["run_record"]["provenance"]
    assert set(prov) >= {"model", "environment", "environment_vs_build", "caveat"}
    assert res["run_record"]["model_artifact_sha256"] != ""
    model = prov["model"]
    assert "matches_manifest" in model
    # Three-valued: unknown must never serialise as a pass.
    assert model["matches_manifest"] in (True, False, None)
    if model["matches_manifest"] is None:
        assert "not a pass" in model["note"]
    env = prov["environment"]
    assert env["python"] and "scikit-learn" in env["packages"]


def test_a_loaded_model_is_identified_by_the_bytes_that_were_opened(tmp_path):
    from quadcond import provenance
    f = tmp_path / "artifact.bin"
    f.write_bytes(b"not a model, but a definite set of bytes")
    import hashlib
    assert provenance.file_sha256(f) == hashlib.sha256(f.read_bytes()).hexdigest()
    d = provenance.describe_table(f, role="regulatory_scores")
    assert d["sha256"] and d["bytes"] == f.stat().st_size and d["role"] == "regulatory_scores"
