"""The atlas fingerprint must distinguish atlases that differ by molecule.

``motifs.clean`` maps U to T, so an RNA row and a DNA row carrying the same
bases are indistinguishable by sequence alone. Version 1 of the fingerprint
hashed a column list that omitted ``nucleic_acid``, which meant two atlases with
different DNA/RNA composition could produce the same identity -- and identity is
the one thing the fingerprint exists to establish.

The fix is versioned rather than in-place: the published v0.5.1 assets record v1
fingerprints, and ``atlas.db``'s recorded identity cannot be recomputed while
that file is unpublished. So v1 must keep producing exactly what it produced
before, and these tests pin that as hard as they pin the fix.
"""
from __future__ import annotations

import sqlite3

import pytest

from quadcond import assets

COLS = ("sequence", "seq_hash", "kind", "k", "na", "li_nh4", "mg", "ph",
        "temperature", "folded", "topology", "tm", "dg", "ph_t",
        "evidence_tier", "label_class", "source", "source_id",
        "condition_imputed", "nucleic_acid")


def _atlas(path, rows):
    """Minimal records table; the fingerprint only ever SELECTs from it."""
    con = sqlite3.connect(path)
    con.execute(f"CREATE TABLE records ({', '.join(c + ' TEXT' for c in COLS)})")
    con.executemany(
        f"INSERT INTO records ({', '.join(COLS)}) "
        f"VALUES ({', '.join('?' * len(COLS))})", rows)
    con.commit()
    con.close()
    return str(path)


def _row(seq, source_id, nucleic_acid):
    return (seq, f"h_{seq}", "G4", "100", "0", "0", "0", "7.0", "25",
            "1", None, "60.0", None, None, "experimental", "biophysical",
            "src", source_id, "", nucleic_acid)


# --------------------------------------------------------------------------
# the defect
# --------------------------------------------------------------------------
def test_v1_cannot_tell_dna_from_rna(tmp_path):
    """The bug, pinned. v1 must keep behaving this way; that is why v2 exists."""
    dna = _atlas(tmp_path / "dna.db", [_row("GGGTTAGGG", "r1", "DNA")])
    rna = _atlas(tmp_path / "rna.db", [_row("GGGTTAGGG", "r1", "RNA")])
    assert (assets.atlas_fingerprint(dna, version=1)
            == assets.atlas_fingerprint(rna, version=1))


def test_v2_tells_dna_from_rna(tmp_path):
    dna = _atlas(tmp_path / "dna.db", [_row("GGGTTAGGG", "r1", "DNA")])
    rna = _atlas(tmp_path / "rna.db", [_row("GGGTTAGGG", "r1", "RNA")])
    assert (assets.atlas_fingerprint(dna, version=2)
            != assets.atlas_fingerprint(rna, version=2))


def test_v2_distinguishes_composition_not_just_presence(tmp_path):
    """Same row count, same sequences, different molecule split."""
    a = _atlas(tmp_path / "a.db", [_row("GGGTTAGGG", "r1", "DNA"),
                                   _row("GGGAAAGGG", "r2", "RNA")])
    b = _atlas(tmp_path / "b.db", [_row("GGGTTAGGG", "r1", "RNA"),
                                   _row("GGGAAAGGG", "r2", "DNA")])
    assert (assets.atlas_fingerprint(a, version=2)
            != assets.atlas_fingerprint(b, version=2))


# --------------------------------------------------------------------------
# backward compatibility -- the published assets must keep verifying
# --------------------------------------------------------------------------
def test_v1_is_byte_identical_to_the_published_definition(tmp_path):
    """Recomputed here from the column list v1 shipped with.

    If this drifts, every content_fingerprint in assets_manifest.json becomes a
    mismatch and 'quadcond assets fetch' breaks for anyone already on v0.5.1.
    """
    import hashlib
    db = _atlas(tmp_path / "a.db", [_row("GGGTTAGGG", "r1", "DNA"),
                                    _row("GGGAAAGGG", "r2", "DNA")])
    con = sqlite3.connect(db)
    legacy_cols = ("sequence, seq_hash, kind, k, na, li_nh4, mg, ph, temperature, "
                   "folded, topology, tm, dg, ph_t, evidence_tier, label_class, "
                   "source, source_id, condition_imputed")
    rows = con.execute(
        f"SELECT {legacy_cols} FROM records ORDER BY source, source_id, seq_hash"
    ).fetchall()
    con.close()
    h = hashlib.sha256()
    for r in rows:
        h.update(repr(r).encode())
    assert assets.atlas_fingerprint(db, version=1) == h.hexdigest()


def test_an_asset_defaults_to_version_one(tmp_path):
    """Manifest entries written before the fix carry no version field."""
    a = assets.Asset(name="x", filename="x.db", sha256=None, bytes=1,
                     kind="atlas", description="", checksum_stable=False,
                     content_fingerprint="deadbeef")
    assert a.fingerprint_version == 1


def test_new_fingerprints_default_to_the_fixed_version():
    assert assets.FINGERPRINT_VERSION == 2


def test_verify_uses_the_version_the_asset_records(tmp_path):
    """A v1 asset must verify under v1, not be failed by the new default."""
    db = _atlas(tmp_path / "a.db", [_row("GGGTTAGGG", "r1", "DNA")])
    v1 = assets.atlas_fingerprint(db, version=1)
    asset = assets.Asset(name="atlas_core", filename="a.db", sha256=None, bytes=1,
                         kind="atlas", description="", checksum_stable=False,
                         content_fingerprint=v1, fingerprint_version=1)
    assets.verify(db, asset)                      # must not raise

    bad = assets.Asset(name="atlas_core", filename="a.db", sha256=None, bytes=1,
                       kind="atlas", description="", checksum_stable=False,
                       content_fingerprint=v1, fingerprint_version=2)
    with pytest.raises(assets.AssetError):
        assets.verify(db, bad)


# --------------------------------------------------------------------------
# details
# --------------------------------------------------------------------------
def test_null_nucleic_acid_reads_as_dna(tmp_path):
    """Rows predating the column must not get a different identity from DNA rows."""
    null = _atlas(tmp_path / "null.db", [_row("GGGTTAGGG", "r1", None)])
    dna = _atlas(tmp_path / "dna.db", [_row("GGGTTAGGG", "r1", "DNA")])
    assert (assets.atlas_fingerprint(null, version=2)
            == assets.atlas_fingerprint(dna, version=2))


def test_unknown_version_is_refused(tmp_path):
    db = _atlas(tmp_path / "a.db", [_row("GGGTTAGGG", "r1", "DNA")])
    with pytest.raises(assets.AssetError, match="unknown atlas fingerprint version"):
        assets.atlas_fingerprint(db, version=99)


def test_fingerprint_is_stable_across_reads(tmp_path):
    """The reason atlases are identified by content and not by file bytes."""
    db = _atlas(tmp_path / "a.db", [_row("GGGTTAGGG", "r1", "DNA")])
    first = assets.atlas_fingerprint(db)
    sqlite3.connect(db).execute("SELECT count(*) FROM records").fetchone()
    assert assets.atlas_fingerprint(db) == first
