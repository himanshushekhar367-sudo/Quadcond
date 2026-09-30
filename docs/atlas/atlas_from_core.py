#!/usr/bin/env python3
"""Assemble data/atlas.db: a rebuild as the base, plus whatever only core has.

The direction matters, and the first version of this script had it backwards.
Core cannot be the base: its `records` table carries
`CHECK (kind IN ('G4','iM'))`, while the current schema allows `'locus'` too, so
inserting the 70,628 GSE220882 locus-state rows into a copy of core fails with

    sqlite3.IntegrityError: CHECK constraint failed: kind IN ('G4','iM')

Both files report `schema_version = 4`, so the version does not distinguish them
-- the constraint was widened without bumping it. This script therefore reads the
constraint rather than trusting the version, and builds from whichever database
can actually hold the union.

    python docs/atlas/atlas_from_core.py --rebuild data/atlas_rebuild.db

Base: the rebuild (all sources buildable from the inputs on disk). Added from
core: only sources the rebuild lacks -- in practice `g4sp_topology` and its
shuffle, 2,010 rows, because `G4 Dataset.xlsx` is not on disk. The result is
written with a `meta` row recording what it is, so nothing downstream mistakes it
for the 398,375-record atlas the shipped models were trained on. See
docs/REBUILD_THE_ATLAS.md.
"""
from __future__ import annotations

import argparse
import re
import sqlite3
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def kinds_allowed(db: Path) -> set[str]:
    con = sqlite3.connect(db)
    sql = con.execute("select sql from sqlite_master where name='records'").fetchone()[0]
    con.close()
    m = re.search(r"kind\s+TEXT[^,]*CHECK\s*\(\s*kind\s+IN\s*\(([^)]*)\)", sql, re.I)
    return {v.strip().strip("'\"") for v in m.group(1).split(",")} if m else set()


def counts(db: Path) -> dict[str, int]:
    con = sqlite3.connect(db)
    out = dict(con.execute("select source, count(*) from records group by source"))
    con.close()
    return out


def kinds_used(db: Path, sources) -> set[str]:
    con = sqlite3.connect(db)
    q = "select distinct kind from records where source in ({})".format(
        ",".join("?" * len(sources)))
    out = {r[0] for r in con.execute(q, tuple(sources))}
    con.close()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rebuild", default=str(REPO / "data" / "atlas_rebuild.db"))
    ap.add_argument("--core", default=str(REPO / "data" / "atlas_core.db"))
    ap.add_argument("--out", default=str(REPO / "data" / "atlas.db"))
    a = ap.parse_args()

    rebuild, core, out = Path(a.rebuild), Path(a.core), Path(a.out)
    for p in (rebuild, core):
        if not p.exists():
            raise SystemExit(f"not found: {p}  (build the rebuild first: "
                             "python scripts/01_build_atlas.py --db data/atlas_rebuild.db)")
    if out.resolve() in (rebuild.resolve(), core.resolve()):
        raise SystemExit("--out would overwrite an input; pick another path")

    rb, cr = counts(rebuild), counts(core)
    only_core = sorted(set(cr) - set(rb))
    needed = kinds_used(core, only_core) if only_core else set()
    allowed = kinds_allowed(rebuild)
    if needed - allowed:
        raise SystemExit(f"the rebuild's schema allows kinds {sorted(allowed)} but core's extra "
                         f"sources need {sorted(needed)}; migrate the schema before merging")

    # SQLite's online backup, not a file copy. shutil.copy2 of a database that
    # has just been written copies the .db and its -wal as two separate,
    # non-atomic operations; the pair can land inconsistent, and what follows is
    # a UNIQUE-constraint error from a damaged index and then "database disk
    # image is malformed". That happened on 30 September and corrupted an output.
    src = sqlite3.connect(f"file:{rebuild}?mode=ro", uri=True)
    dst = sqlite3.connect(out)
    with dst:
        src.backup(dst)
    src.close(); dst.close()

    added = {}
    if only_core:
        con = sqlite3.connect(out)
        con.execute("attach database ? as core", (str(core),))
        cols = [r[1] for r in con.execute("PRAGMA table_info(records)")]
        keep = [c for c in cols if c != "record_id"]
        core_cols = {r[1] for r in con.execute("PRAGMA core.table_info(records)")}
        keep = [c for c in keep if c in core_cols]
        for src in only_core:
            con.execute(f"insert into records ({','.join(keep)}) "
                        f"select {','.join(keep)} from core.records where source=?", (src,))
            added[src] = cr[src]
        # Source rows travel with their records, or provenance is lost.
        src_cols = [r[1] for r in con.execute("PRAGMA table_info(sources)")]
        common = [c for c in src_cols if c in {r[1] for r in con.execute(
            "PRAGMA core.table_info(sources)")}]
        con.execute(f"insert or ignore into sources ({','.join(common)}) "
                    f"select {','.join(common)} from core.sources where source in ({','.join('?' * len(only_core))})",
                    tuple(only_core))
        con.commit()
        con.execute("detach database core")
        # Integrity check before anything downstream reads it.
        bad = con.execute("PRAGMA integrity_check").fetchone()[0]
        con.close()
        if bad != "ok":
            raise SystemExit(f"{out} failed integrity_check after the merge: {bad}")

    con = sqlite3.connect(out)
    total = con.execute("select count(*) from records").fetchone()[0]
    nsrc = con.execute("select count(distinct source) from records").fetchone()[0]
    con.execute("create table if not exists meta (key text primary key, value text)")
    note = (f"Assembled {date.today().isoformat()}: base {rebuild.name}"
            + (f", plus {len(added)} source(s) from {core.name}" if added else "")
            + f". {total} records, {nsrc} sources. NOT the 398,375-record atlas the shipped "
              "v0.5.1 models were trained on, which no longer exists. Do not quote its "
              "fingerprint as the training dataset fingerprint. See docs/REBUILD_THE_ATLAS.md.")
    con.execute("insert or replace into meta (key, value) values ('assembly_provenance', ?)",
                (note,))
    con.commit(); con.close()

    print(f"wrote {out}\n  {total:,} records, {nsrc} sources")
    for s, n in added.items():
        print(f"  + {s}: {n:,} rows from {core.name}")
    print(f"\n{note}")
    print(f"\nVerify:  python docs/atlas/compare_rebuild.py --db {out}")


if __name__ == "__main__":
    main()
