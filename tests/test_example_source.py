"""The workbench's chr22 example must show both axes without any configuration."""
from quadcond import alphagenome as ag
from quadcond import service


def test_table_env_is_a_configured_source(tmp_path, monkeypatch):
    table = tmp_path / "scores.tsv"
    table.write_text("chromosome\tposition\treference\talternate\tscorer\tscore\n"
                     "chr1\t100\tA\tC\tAVI_SCORE\t0.5\n")
    monkeypatch.delenv(ag.AVI_PATH_ENV, raising=False)
    monkeypatch.setenv(ag.TABLE_PATH_ENV, str(table))
    src = ag.resolve_source()
    assert src.describe()["source"] == "table"


def test_tabix_still_wins_over_a_table(tmp_path, monkeypatch):
    monkeypatch.setenv(ag.TABLE_PATH_ENV, str(tmp_path / "unused.tsv"))
    monkeypatch.setenv(ag.AVI_PATH_ENV, str(tmp_path / "missing.tsv.gz"))
    try:
        ag.resolve_source()
    except ag.AtlasUnavailable as exc:
        assert "missing.tsv.gz" in str(exc)


def test_bundled_example_covers_exactly_the_worked_example():
    src = ag.ExampleTableAtlas.from_path(ag.EXAMPLE_TABLE)
    d = src.describe()
    assert d["source"] == "example_table"
    assert d["n_variants"] == 183
    assert "only the bundled worked example" in d["note"]
    recs = src.records_for_interval("chr22", 36201667, 36201728)
    assert len(recs) == 183
    assert src.records_for_interval("chr1", 1_000_000, 1_000_100) == {}


def test_service_falls_back_to_the_example_only_when_nothing_is_configured(monkeypatch):
    for env in (ag.AVI_PATH_ENV, ag.TABLE_PATH_ENV, "ALPHAGENOME_API_KEY"):
        monkeypatch.delenv(env, raising=False)
    monkeypatch.setattr(service, "_regulatory", None)
    assert service._regulatory_source().describe()["source"] == "example_table"
    monkeypatch.setattr(service, "_regulatory", None)
