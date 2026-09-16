"""The atlas: a queryable, provenance-carrying store of G4 / i-motif observations."""
from __future__ import annotations

import hashlib
import json
import sqlite3
import threading
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

from ..conditions import Condition, condition_distance
from ..motifs import clean
from .schema import CONDITION_COLUMNS, DDL, SCHEMA_VERSION


def seq_hash(seq: str) -> str:
    return hashlib.sha1(clean(seq).encode()).hexdigest()[:16]


@dataclass
class Record:
    """One measurement of one sequence under one condition."""

    sequence: str
    kind: str
    source: str
    evidence_tier: str
    condition: Condition = field(default_factory=Condition)
    label_class: str = "biophysical"
    nucleic_acid: str = "DNA"
    folded: int | None = None
    topology: str | None = None
    tm: float | None = None
    dg: float | None = None
    ph_t: float | None = None
    method: str | None = None
    source_doi: str | None = None
    source_id: str | None = None
    organism: str | None = None
    genomic: dict | None = None
    qc_flags: list[str] = field(default_factory=list)

    def row(self) -> tuple:
        c = self.condition
        return (
            clean(self.sequence), seq_hash(self.sequence), self.kind, self.nucleic_acid,
            c.k, c.na, c.li_nh4, c.mg, c.ph, c.temperature, c.crowder_pct, c.strand_conc,
            json.dumps(list(c.imputed)),
            self.folded, self.topology, self.tm, self.dg, self.ph_t,
            self.evidence_tier, self.label_class,
            self.method, self.source, self.source_doi, self.source_id,
            self.organism,
            json.dumps(self.genomic) if self.genomic else None,
            json.dumps(self.qc_flags) if self.qc_flags else None,
        )


_INSERT = """
INSERT OR IGNORE INTO records
 (sequence, seq_hash, kind, nucleic_acid,
  k, na, li_nh4, mg, ph, temperature, crowder_pct, strand_conc, condition_imputed,
  folded, topology, tm, dg, ph_t,
  evidence_tier, label_class, method, source, source_doi, source_id, organism, genomic, qc_flags)
VALUES (?,?,?,?, ?,?,?,?,?,?,?,?,?, ?,?,?,?,?, ?,?,?,?,?,?,?,?,?)
"""


class Atlas:
    """Thin, explicit wrapper over the SQLite store.

    The connection is per-thread, not per-instance. One ``Atlas`` is built once
    and then read from every request thread of ``ThreadingHTTPServer``, and
    SQLite refuses a connection used off the thread that created it:

        ProgrammingError: SQLite objects created in a thread can only be used
        in that same thread.

    With a sentinel model and no atlas nothing ever opened a connection, so the
    whole test suite and the fixture server passed while every real ``/predict``
    after the first request returned a 500. A thread-local handle keeps the
    single-threaded callers (the CLI, the ingest scripts, ``report.py``)
    unchanged -- ``self.conn`` still reads like an attribute -- while giving
    each server thread its own.
    """

    def __init__(self, path: str | Path = "data/atlas.db"):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._local = threading.local()
        conn = self._open()
        conn.executescript(DDL)
        conn.execute(
            "INSERT OR REPLACE INTO meta(key,value) VALUES ('schema_version',?)",
            (str(SCHEMA_VERSION),),
        )
        conn.commit()
        self._vec_caches: dict = {}

    def _open(self) -> sqlite3.Connection:
        conn = sqlite3.connect(self.path)
        conn.row_factory = sqlite3.Row
        self._local.conn = conn
        return conn

    @property
    def conn(self) -> sqlite3.Connection:
        """This thread's connection, opened on first use.

        Deliberately not ``check_same_thread=False``: that silences the error
        without making concurrent use safe, and a shared cursor returning
        another thread's rows is a far worse failure than a refusal.
        """
        conn = getattr(self._local, "conn", None)
        return conn if conn is not None else self._open()

    # ------------------------------------------------------------------ write
    def register_source(
        self,
        source: str,
        *,
        title: str = "",
        doi: str = "",
        url: str = "",
        evidence_tier: str = "experimental",
        licence: str = "",
        retrieved_at: str = "",
        notes: str = "",
    ) -> None:
        self.conn.execute(
            """INSERT OR REPLACE INTO sources
               (source,title,doi,url,evidence_tier,licence,retrieved_at,notes)
               VALUES (?,?,?,?,?,?,?,?)""",
            (source, title, doi, url, evidence_tier, licence, retrieved_at, notes),
        )
        self.conn.commit()

    def add(self, records: Iterable[Record]) -> int:
        rows = [r.row() for r in records]
        if not rows:
            return 0
        before = self.count()
        self.conn.executemany(_INSERT, rows)
        self.conn.commit()
        self._vec_caches = {}
        return self.count() - before

    # ------------------------------------------------------------------- read
    def count(self, **where: Any) -> int:
        sql = "SELECT COUNT(*) FROM records"
        params: list[Any] = []
        if where:
            sql += " WHERE " + " AND ".join(f"{k}=?" for k in where)
            params = list(where.values())
        return self.conn.execute(sql, params).fetchone()[0]

    def query(
        self,
        *,
        kind: str | None = None,
        tiers: Sequence[str] | None = None,
        label_classes: Sequence[str] | None = None,
        label: str | None = None,
        sources: Sequence[str] | None = None,
        limit: int | None = None,
    ) -> list[sqlite3.Row]:
        """Fetch records, optionally restricted to those carrying a given label."""
        sql = "SELECT * FROM records WHERE 1=1"
        params: list[Any] = []
        if kind:
            sql += " AND kind=?"
            params.append(kind)
        if tiers:
            sql += f" AND evidence_tier IN ({','.join('?' * len(tiers))})"
            params += list(tiers)
        if label_classes:
            sql += f" AND label_class IN ({','.join('?' * len(label_classes))})"
            params += list(label_classes)
        if sources:
            sql += f" AND source IN ({','.join('?' * len(sources))})"
            params += list(sources)
        if label:
            sql += f" AND {label} IS NOT NULL"
        if limit:
            sql += " LIMIT ?"
            params.append(limit)
        return self.conn.execute(sql, params).fetchall()

    def composition(self) -> dict:
        """The total, and what it is made of. Never report one without the other.

        398,375 rows sounds like 398,375 measurements of DNA folding. It is not:
        the great majority are another predictor's output or generated shuffles,
        and the rows that carry a physical measurement of the structure being
        predicted are a few thousand. A headline count with no breakdown beside
        it is the single easiest way for this project to mislead someone, so the
        breakdown is computed here and every surface that prints the total pulls
        it from this method.
        """
        rows = self.conn.execute(
            "SELECT label_class, evidence_tier, COUNT(*) n FROM records "
            "GROUP BY 1, 2"
        ).fetchall()
        buckets = {"biophysical_measured": 0, "genomic_proxy": 0,
                   "predicted": 0, "derived_or_shuffle": 0}
        for lc, tier, n in rows:
            if lc == "genomic_proxy":
                buckets["genomic_proxy"] += n
            elif tier == "experimental" and lc == "biophysical":
                buckets["biophysical_measured"] += n
            elif tier == "predicted":
                buckets["predicted"] += n
            else:
                buckets["derived_or_shuffle"] += n
        total = sum(buckets.values())
        return {
            "total": total,
            **buckets,
            "labels": {
                "biophysical_measured":
                    "a physical measurement of the structure being predicted "
                    "(melting curves, pH titrations)",
                "genomic_proxy":
                    "antibody occupancy at a genomic locus (BG4 / iMab CUT&Tag) "
                    "-- an experimental observation of a different quantity",
                "predicted": "another model's output; no measurement upstream",
                "derived_or_shuffle":
                    "computed from sequence, or a composition-matched shuffle",
            },
            "one_line": (
                f"{total:,} rows: {buckets['biophysical_measured']:,} biophysical "
                f"measurements, {buckets['genomic_proxy']:,} genomic proxy, "
                f"{buckets['predicted']:,} predicted, "
                f"{buckets['derived_or_shuffle']:,} derived/shuffle"
            ),
        }

    def frame(self, **kw):
        import pandas as pd

        rows = self.query(**kw)
        return pd.DataFrame([dict(r) for r in rows])

    def summary(self):
        import pandas as pd

        q = """
        SELECT source, kind, evidence_tier, label_class, COUNT(*) AS n,
               SUM(folded IS NOT NULL)   AS n_folded,
               SUM(topology IS NOT NULL) AS n_topology,
               SUM(tm IS NOT NULL)       AS n_tm,
               SUM(ph_t IS NOT NULL)     AS n_pht,
               COUNT(DISTINCT seq_hash)  AS n_unique_seq,
               -- Mg2+ is part of the buffer. Leaving it out of this key gave
               -- 255 for g4stab where the report's own count said 261, and the
               -- page showed both numbers with neither defined. One key.
               COUNT(DISTINCT k || '_' || na || '_' || li_nh4 || '_' || mg || '_' || ph)
                 AS n_conditions
        FROM records GROUP BY source, kind, evidence_tier, label_class ORDER BY n DESC
        """
        return pd.read_sql_query(q, self.conn)

    # -------------------------------------------------------------- retrieval
    @staticmethod
    def _kmer_vector(seq: str, k: int = 3) -> np.ndarray:
        import itertools

        kmers = ["".join(p) for p in itertools.product("ACGT", repeat=k)]
        idx = {km: i for i, km in enumerate(kmers)}
        v = np.zeros(len(kmers), dtype=np.float32)
        s = clean(seq)
        for i in range(len(s) - k + 1):
            j = idx.get(s[i : i + k])
            if j is not None:
                v[j] += 1
        n = np.linalg.norm(v)
        return v / n if n else v

    def _ensure_vectors(self, kind: str | None = None, tiers=None):
        """Build (and memoise) the k-mer index used for neighbour search.

        Recomputing this per query is the difference between a millisecond and
        several seconds once the atlas holds hundreds of thousands of rows, so
        it is cached per (kind, tiers) key and invalidated on write.
        """
        key = (kind, tuple(tiers) if tiers else None)
        cache = getattr(self, "_vec_caches", None)
        if cache is None:
            cache = self._vec_caches = {}
        if key in cache:
            return cache[key]
        rows = self.query(kind=kind, tiers=tiers)
        ids = [r["record_id"] for r in rows]
        mat = (np.vstack([self._kmer_vector(r["sequence"]) for r in rows])
               if rows else np.zeros((0, 64), np.float32))
        cache[key] = (ids, mat, {r["record_id"]: r for r in rows})
        return cache[key]

    def neighbours(
        self,
        sequence: str,
        condition: Condition | None = None,
        *,
        kind: str | None = None,
        n: int = 5,
        condition_weight: float = 0.35,
        tiers: Sequence[str] | None = ("experimental", "derived"),
    ) -> list[dict]:
        """Nearest experimental neighbours in joint sequence x condition space.

        Similarity = cosine on 3-mer profiles, penalised by scaled condition
        distance.  This is what lets a user ask "what has actually been
        measured near this sequence, near these conditions?" -- the honest
        answer to a prediction they should not fully trust.
        """
        ids, mat, rows_by_id = self._ensure_vectors(kind, tiers)
        if not ids:
            return []
        q = self._kmer_vector(sequence)
        sims = mat @ q
        cond = condition or Condition()
        out: list[dict] = []
        order = np.argsort(-sims)[: max(n * 20, 50)]
        for i in order:
            r = rows_by_id[ids[i]]
            rc = Condition.from_mapping({c: r[c] for c in CONDITION_COLUMNS}, track_imputed=False)
            cd = condition_distance(cond, rc)
            score = float(sims[i]) - condition_weight * cd / (1.0 + cd)
            out.append(
                {
                    "record_id": r["record_id"],
                    "sequence": r["sequence"],
                    "kind": r["kind"],
                    "sequence_similarity": round(float(sims[i]), 4),
                    "condition_distance": round(cd, 3),
                    "joint_score": round(score, 4),
                    "condition": rc.label(),
                    "topology": r["topology"],
                    "tm": r["tm"],
                    "ph_t": r["ph_t"],
                    "folded": r["folded"],
                    "evidence_tier": r["evidence_tier"],
                    "method": r["method"],
                    "source": r["source"],
                    "source_id": r["source_id"],
                    "source_doi": r["source_doi"],
                }
            )
        out.sort(key=lambda d: -d["joint_score"])
        return out[:n]

    def close(self):
        """Close this thread's connection.

        Only this thread's: connections are per-thread and there is no registry
        of the others, deliberately. The callers that close are single-threaded
        (the CLI, the ingest scripts, test teardown), and a server thread's
        connection is released when that thread ends. A cross-thread close
        would need a lock and could pull a connection out from under a request
        mid-query, which is a worse failure than a handle that outlives its
        usefulness by a few milliseconds.
        """
        conn = getattr(self._local, "conn", None)
        if conn is not None:
            conn.close()
            self._local.conn = None
