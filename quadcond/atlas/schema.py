"""SQLite schema for the harmonised G4 / i-motif atlas.

Design rules the schema enforces:

* **One row = one measurement under one condition.**  The same sequence
  measured in K+ and in Na+ is two rows, not one row with two columns.  This is
  what makes the atlas condition-resolved rather than condition-collapsed.
* **Every row declares its evidence tier.**  ``experimental`` (a measurement),
  ``derived`` (constructed by us, e.g. a composition-matched shuffle), or
  ``predicted`` (output of somebody's model).  Training code must opt in to a
  tier explicitly; nothing mixes silently.
* **Every row also declares its label class**, which is a different question
  from the tier.  A BG4 ChIP peak and a UV melting curve are both
  ``experimental``, but one measures a structure in a buffer and the other
  measures antibody occupancy in a nucleus.  ``biophysical`` labels may train
  folding and stability heads; ``genomic_proxy`` labels must be trained and
  evaluated separately and never silently pooled with them; ``catalog`` entries
  (motif calls, candidate lists) are not folding evidence at all in either
  direction.
* **Every row carries provenance**: source key, DOI, and the source's own id
  (PDB code, GEO accession, table row).
* **Imputed condition fields are recorded**, so "pH 7.0" that came from a
  default is distinguishable from "pH 7.0" that was reported.
"""
from __future__ import annotations

SCHEMA_VERSION = 4

DDL = """
PRAGMA journal_mode=WAL;
PRAGMA foreign_keys=ON;

CREATE TABLE IF NOT EXISTS meta (
    key   TEXT PRIMARY KEY,
    value TEXT
);

CREATE TABLE IF NOT EXISTS sources (
    source        TEXT PRIMARY KEY,
    title         TEXT,
    doi           TEXT,
    url           TEXT,
    evidence_tier TEXT NOT NULL,
    licence       TEXT,
    retrieved_at  TEXT,
    notes         TEXT
);

CREATE TABLE IF NOT EXISTS records (
    record_id     INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence      TEXT NOT NULL,
    seq_hash      TEXT NOT NULL,
    -- 'locus' is a duplex position, not a single structure: the G-rich and
    -- C-rich strands of the same site compete for the same opened duplex, so a
    -- row asking "which structure does this LOCUS carry" cannot be filed under
    -- either one. It exists because BG4 and iMab were assayed on the same cell
    -- population -- in PARALLEL REACTIONS on separate aliquots, so what the rows
    -- carry is peak overlap, not co-occupancy: nothing in that design observes
    -- both structures on one DNA molecule.
    kind          TEXT NOT NULL CHECK (kind IN ('G4','iM','locus')),
    nucleic_acid  TEXT DEFAULT 'DNA',

    k             REAL, na REAL, li_nh4 REAL, mg REAL,
    ph            REAL, temperature REAL,
    crowder_pct   REAL, strand_conc REAL,
    condition_imputed TEXT,

    folded        INTEGER,
    topology      TEXT,
    tm            REAL,
    dg            REAL,
    ph_t          REAL,

    evidence_tier TEXT NOT NULL CHECK (evidence_tier IN ('experimental','derived','predicted')),
    label_class   TEXT NOT NULL DEFAULT 'biophysical'
                  CHECK (label_class IN ('biophysical','genomic_proxy','catalog','synthetic')),
    method        TEXT,
    source        TEXT NOT NULL REFERENCES sources(source),
    source_doi    TEXT,
    source_id     TEXT,
    organism      TEXT,
    genomic       TEXT,
    qc_flags      TEXT,
    added_at      TEXT DEFAULT CURRENT_TIMESTAMP,

    UNIQUE (seq_hash, kind, source, source_id, k, na, li_nh4, mg, ph, temperature)
);

CREATE INDEX IF NOT EXISTS idx_records_kind      ON records(kind);
CREATE INDEX IF NOT EXISTS idx_records_tier      ON records(evidence_tier);
CREATE INDEX IF NOT EXISTS idx_records_labelcls ON records(label_class);
CREATE INDEX IF NOT EXISTS idx_records_source    ON records(source);
CREATE INDEX IF NOT EXISTS idx_records_hash      ON records(seq_hash);
CREATE INDEX IF NOT EXISTS idx_records_tm        ON records(tm)   WHERE tm   IS NOT NULL;
CREATE INDEX IF NOT EXISTS idx_records_pht       ON records(ph_t) WHERE ph_t IS NOT NULL;
CREATE INDEX IF NOT EXISTS idx_records_topology  ON records(topology) WHERE topology IS NOT NULL;
"""

LABEL_COLUMNS = ("folded", "topology", "tm", "dg", "ph_t")
LABEL_CLASSES = ("biophysical", "genomic_proxy", "catalog", "synthetic")
CONDITION_COLUMNS = ("k", "na", "li_nh4", "mg", "ph", "temperature", "crowder_pct", "strand_conc")
TOPOLOGIES = ("parallel", "antiparallel", "hybrid")
