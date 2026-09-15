#!/usr/bin/env python3
"""Build the QuadCond atlas from every source reachable in this environment.

Run:  python scripts/01_build_atlas.py [--db data/atlas.db] [--g4stab-seqs 40000]
"""
from __future__ import annotations

import argparse
import hashlib
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from quadcond.atlas import Atlas, Record
from quadcond.atlas.ingest import (duvma, g4sp, g4stab_db, g4stab_supp,
                                   gse220882, gse220882_loci, imseeker,
                                   shuffled)
from quadcond.conditions import Condition
from quadcond.motifs import clean, find_im, revcomp

DOWNLOADS = Path("data/external/downloads")


def build_g4sp(atlas: Atlas) -> int:
    path = DOWNLOADS / "G4 Dataset.xlsx"
    if not path.exists():
        print(f"  ! {path} not found; skipping G4SP")
        return 0
    g4sp.register(atlas)
    recs = g4sp.load(path)
    n = atlas.add(recs)
    print(f"  + g4sp_topology: {n} records ({len(recs)} parsed)")
    return n


def build_g4stab_supp(atlas: Atlas) -> int:
    """The 2,382-measurement experimental Tm table, if it has been downloaded.

    This is the only source in the atlas that varies the buffer across the same
    sequences, so it is the one that makes a condition-aware stability head
    possible at all. Not redistributed with the repo -- see docs/GET_THE_DATA.md.
    """
    candidates = sorted(DOWNLOADS.glob("*G4STAB*Supplementary*Table*1*")) + \
                 sorted(DOWNLOADS.glob("*g4stab*supp*"))
    if not candidates:
        print("  . G4STAB Supplementary Table 1 not found; skipping "
              "(see docs/GET_THE_DATA.md). The g4_tm head cannot train without it.")
        return 0
    path = candidates[0]
    g4stab_supp.register(atlas)
    recs = g4stab_supp.load(path)
    n = atlas.add(recs)
    st = getattr(g4stab_supp.load, "last_stats", {})
    cov = getattr(g4stab_supp.load, "last_buffer_coverage", {})
    print(f"  + g4stab_experimental_tm: {n} records from {path.name}")
    print(f"      buffers parsed: {cov.get('usable', 0)}/{cov.get('rows', 0)} "
          f"({cov.get('usable_fraction', 0):.1%}); rows dropped for an unresolved "
          f"buffer: {st.get('rejected_buffer', 0)}")
    print(f"      Tm usable: {st.get('tm_usable', 0)}; censored: {st.get('tm_censored', 0)}; "
          f"multiphasic: {st.get('tm_multiphasic', 0)} (ingested as folding evidence, tm NULL)")
    return n


def build_imseeker(atlas: Atlas) -> int:
    """160 measured i-motif transitional pH values, if the file is present.

    Small, but the only i-motif folding measurements in the atlas -- and the
    difference between a pH-aware head and a placeholder.
    """
    candidates = (sorted(DOWNLOADS.glob("*Supplementary Data Set*.xlsx")) +
                  sorted(DOWNLOADS.glob("*gkae092*.xlsx")) +
                  sorted(DOWNLOADS.glob("*imseeker*")))
    if not candidates:
        print("  . iM-Seeker Supplementary Data Set not found; skipping "
              "(see docs/GET_THE_DATA.md). The im_pht head cannot train without it.")
        return 0
    path = candidates[0]
    imseeker.register(atlas)
    recs = imseeker.load(path)
    n = atlas.add(recs)
    st = getattr(imseeker.load, "last_stats", {})
    print(f"  + im_pht_biophysical: {n} records from {path.name}")
    print(f"      kept {st.get('kept', 0)}/{st.get('rows', 0)}; rejected for deletion "
          f"notation or bad bases: {st.get('rejected_sequence', 0)}; "
          f"{st.get('with_detected_motif', 0)} annotated with a detected motif")
    return n


def build_duvma(atlas: Atlas) -> int:
    """The condition axis for i-motifs: two constructs across pH x ionic strength.

    Parsed straight out of the supporting-information PDF. Small in sequences and
    large in conditions -- the exact complement of the iM-Seeker panel.
    """
    candidates = (sorted(DOWNLOADS.glob("*gkag110*.pdf")) +
                  sorted(DOWNLOADS.glob("*5duvma*")))
    if not candidates:
        print("  . 5DUVMA supporting information not found; skipping "
              "(see docs/GET_THE_DATA.md). im_pht will have no salt axis.")
        return 0
    path = candidates[0]
    duvma.register(atlas)
    recs = duvma.load(path)
    n = atlas.add(recs)
    st = getattr(duvma.load, "last_stats", {})
    n_tm = sum(1 for r in recs if r.tm is not None)
    n_ph = sum(1 for r in recs if r.ph_t is not None)
    print(f"  + duvma_im_stability_landscape: {n} records from {path.name}")
    print(f"      tables parsed: {st}")
    print(f"      {n_tm} Tm and {n_ph} pH_T values, 2 constructs, "
          f"10 ionic strengths (10-1000 mM K+)")
    return n


def build_gse220882(atlas: Atlas, cell_lines=("HEK293T",)) -> int:
    """The sequence-diversity axis for i-motifs: iMab/BG4 CUT&Tag peaks.

    Positives are the best canonical motif inside each reproducible peak window;
    negatives are the best canonical motif inside a dinucleotide shuffle of the
    same window. Both carry a four-tract motif and match on composition, so the
    head has to separate them on architecture. See the adapter's docstring.

    The input JSONL is produced by ``scripts/06_gse220882_peaks.py``, which runs
    wherever the GEO bigWigs are.
    """
    candidates = sorted(DOWNLOADS.glob("gse220882_peaks.jsonl*"))
    if not candidates:
        print("  . GSE220882 peak file not found; skipping (run "
              "scripts/06_gse220882_peaks.py first -- see docs/GET_THE_DATA.md). "
              "im_fold will have only the 160-oligo biophysical panel.")
        return 0
    path = candidates[0]
    gse220882.register(atlas, cell_lines=cell_lines)
    pos = gse220882.load(path, cell_lines=cell_lines)
    neg = gse220882.negatives(path, cell_lines=cell_lines)
    n = atlas.add(pos) + atlas.add(neg)
    st = getattr(gse220882.load, "last_stats", {})
    for key in sorted(k for k in st if k.endswith("|windows")):
        grp = key.rsplit("|", 1)[0]
        w, kept = st[key], st.get(f"{grp}|kept", 0)
        bg = gse220882.GENOMIC_BACKGROUND[grp.split("|")[1]]
        print(f"      {grp:16s} {w:>7,} reproducible peaks -> {kept:>6,} carry a "
              f"canonical motif ({kept/w:.1%}; genomic background {bg:.1%}, "
              f"{kept/w/bg:.1f}x)")
    print(f"  + gse220882 CUT&Tag: {n} records from {path.name} "
          f"({len(pos)} positives, {len(neg)} matched negatives)")
    return n


def build_gse220882_loci(atlas: Atlas) -> int:
    """Joint G4 / i-motif occupancy, from the two antibodies run in the same cells.

    This is the only source in the atlas that can say anything about the two
    structures *competing*, because it is the only one where both were assayed
    at the same positions under one condition. Everything else in the atlas
    measures one structure at a time, which is why the competition output was an
    independence proxy rather than a model.
    """
    candidates = sorted(DOWNLOADS.glob("gse220882_peaks.jsonl*"))
    if not candidates:
        print("  . GSE220882 peak file not found; skipping the locus-state source")
        return 0
    path = candidates[0]
    gse220882_loci.register(atlas)
    pos = gse220882_loci.load(path)
    neg = gse220882_loci.negatives(path)
    n = atlas.add(pos) + atlas.add(neg)
    st = getattr(gse220882_loci.load, "last_stats", {})
    total = st.get("loci", 0)
    print(f"  + gse220882 locus states: {n} records from {total:,} merged loci")
    for k in ("both", "g4_only", "im_only"):
        v = st.get(k, 0)
        print(f"      {k:<10} {v:>7,}  ({v / max(total, 1):.1%})")
    print(f"      {'neither':<10} {len(neg):>7,}  (composition-matched shuffles)")
    return n


def build_negatives(atlas: Atlas, n_per_positive: int = 1) -> int:
    """Composition-matched negatives, generated per parent source.

    Each parent source gets its own negative source, ``shuffled::<parent>``, for
    two reasons. Provenance: a negative should point at the positive it was
    derived from. And balance: 5DUVMA contributes 937 measurements of just two
    constructs, so pooling all i-motif shuffles under one name would hand a
    folding classifier 160 positives against 1,097 negatives, with the negative
    class carrying almost no extra sequence diversity for the trouble.
    """
    total = 0
    for kind in ("G4", "iM"):   # locus rows bring their own matched negatives
        parents: dict[str, list[dict]] = {}
        for r in atlas.query(kind=kind, tiers=("experimental",)):
            # genomic_proxy sources arrive with their own negatives, shuffled at
            # the peak-window level rather than at the motif level -- a much
            # tighter control. Shuffling them again here would both collide on
            # the source name and replace that control with a weaker one.
            if r["folded"] == 1 and r["label_class"] != "genomic_proxy":
                parents.setdefault(r["source"], []).append(dict(r))
        for parent_source, pos in parents.items():
            neg_source = f"shuffled::{parent_source}"
            atlas.register_source(
                neg_source,
                title=f"Dinucleotide-shuffled negatives for {parent_source}",
                evidence_tier="derived",
                notes="Altschul-Erikson shuffle of each positive in the parent source: "
                      "identical mono- and dinucleotide counts, four-tract register "
                      "destroyed. label_class=catalog -- the negative label is an "
                      "argument from motif architecture, not a measured non-folding event.",
            )
            # Not `hash()`: Python randomises string hashing per process, so
            # that seed changed on every run and the negatives -- and every
            # metric computed against them -- changed with it. A 16-sequence
            # external check swung by 0.17 between two builds of what was
            # supposed to be the same atlas before this was tracked down.
            seed = int(hashlib.sha1(parent_source.encode()).hexdigest()[:8], 16)
            recs = shuffled.build(pos, kind=kind, n_per_positive=n_per_positive,
                                  seed=seed % 10_000)
            for r in recs:
                r.source = neg_source
            n = atlas.add(recs)
            total += n
            print(f"  + {neg_source}: {n}")
    return total


def build_g4stab_predicted(atlas: Atlas, n_sequences: int, seed: int = 11) -> int:
    d = DOWNLOADS / "g4stab_db"
    files = sorted(d.glob("pred_seqs_*.csv"))
    if not files:
        print("  ! g4stab web-database shards not found; skipping")
        return 0
    g4stab_db.register(atlas)

    # Pick a deterministic random subset of sequences, then keep every ionic
    # condition for those sequences -- the condition response is the signal.
    rng = random.Random(seed)
    keep: set[str] = set()
    print("  . scanning shards for a sequence subset ...")
    first = pd.read_csv(files[0], usecols=["seq"])
    pool = first["seq"].unique().tolist()
    rng.shuffle(pool)
    keep = set(pool[:n_sequences])

    recs: list[Record] = []
    total = 0
    for f in files:
        df = pd.read_csv(f, dtype={"Chr": str})
        df = df[df["seq"].isin(keep)]
        if df.empty:
            continue
        for _, r in df.iterrows():
            cond = Condition.from_mapping(
                {"k": r["k"], "na": r["na"], "li_nh4": r["li/nh4"],
                 "ph": 7.0, "temperature": 25.0},
                track_imputed=True,
            )
            recs.append(
                Record(
                    sequence=str(r["seq"]), kind="G4",
                    source=g4stab_db.SOURCE, evidence_tier="predicted",
                    label_class="biophysical",
                    condition=cond, tm=float(r["predicted_temperature"]),
                    method="G4STAB deep-ensemble prediction",
                    source_doi=g4stab_db.DOI,
                    source_id=f"{r['Chr']}:{r['Start']}-{r['End']}{r['Strand']}",
                    organism="Homo sapiens (GRCh38)",
                    genomic={"chrom": str(r["Chr"]), "start": int(r["Start"]),
                             "end": int(r["End"]), "strand": r["Strand"],
                             "gene": r["Symbol"], "gene_type": r["Gene Type"],
                             "phastCons": float(r["phastCons"]),
                             "phyloP": float(r["phyloP"])},
                    qc_flags=["model_output_not_measurement", f"sem={r['sem']}"],
                )
            )
            if len(recs) >= 25_000:
                total += atlas.add(recs)
                recs = []
        print(f"    {f.name}: cumulative {total}")
    if recs:
        total += atlas.add(recs)
    print(f"  + g4stab_webdb_predicted: {total} records")
    return total


def build_complementary_im(atlas: Atlas, limit: int = 30_000) -> int:
    """Complementary-strand C-rich candidates from annotated genomic G4 loci.

    A genomic G4 on one strand implies a C-rich tract on the other; whether it
    folds as an i-motif is an open question, so these enter as DERIVED
    candidates (never as folding evidence) and exist mainly to give the joint
    G4/iM competition module real complementary pairs to reason about.
    """
    source = "complementary_strand_im_candidates"
    atlas.register_source(
        source,
        title="C-rich complementary-strand i-motif candidates from annotated genomic G4 loci",
        evidence_tier="derived",
        notes="Reverse complement of G4STAB web-database genomic PQS. Motif-derived "
              "candidates, NOT folding evidence. Used for joint G4/iM competition analysis.",
    )
    rows = atlas.query(kind="G4", tiers=("predicted",), limit=limit * 5)
    seen: set[str] = set()
    recs: list[Record] = []
    import json as _json

    for r in rows:
        gid = r["source_id"]
        if gid in seen:
            continue
        seen.add(gid)
        rc = revcomp(clean(r["sequence"]))
        if not find_im(rc):
            continue
        genomic = _json.loads(r["genomic"]) if r["genomic"] else None
        if genomic:
            genomic = dict(genomic)
            genomic["strand"] = "-" if genomic.get("strand") == "+" else "+"
            genomic["paired_g4"] = r["sequence"]
        recs.append(
            Record(
                sequence=rc, kind="iM", source=source, evidence_tier="derived",
                label_class="catalog",
                condition=Condition.from_mapping({"ph": 5.8}, track_imputed=True),
                folded=None,
                method="reverse complement of annotated genomic G4 locus",
                source_id=f"rc:{gid}", organism="Homo sapiens (GRCh38)",
                genomic=genomic,
                qc_flags=["motif_derived_candidate", "not_folding_evidence"],
            )
        )
        if len(recs) >= limit:
            break
    n = atlas.add(recs)
    print(f"  + complementary_strand_im_candidates: {n}")
    return n


def build_im_architecture_set(atlas: Atlas) -> int:
    """Positive/negative pairs for the i-motif architecture head.

    Positives are the complementary-strand candidates that carry a canonical
    four-C-tract architecture; negatives are their dinucleotide shuffles.  This
    trains a head that scores *architecture quality*, explicitly not folding.
    """
    cands = [dict(r) for r in atlas.query(kind="iM", sources=("complementary_strand_im_candidates",))]
    if not cands:
        return 0
    shuffled.register(atlas)
    source = "im_architecture_positives"
    atlas.register_source(
        source,
        title="i-motif architecture positives (canonical four-C-tract candidates)",
        evidence_tier="derived",
        notes="folded=1 here means 'has canonical iM architecture', NOT 'folds in solution'.",
    )
    recs = [
        Record(
            sequence=c["sequence"], kind="iM", source=source, evidence_tier="derived",
            label_class="catalog",
            condition=Condition.from_mapping({"ph": 5.8}, track_imputed=True),
            folded=1, method="canonical four-C-tract motif call",
            source_id=c["source_id"], organism=c["organism"],
            qc_flags=["architecture_label_only", "not_folding_evidence"],
        )
        for c in cands
    ]
    n = atlas.add(recs)
    neg_source = "im_architecture_negatives"
    atlas.register_source(
        neg_source,
        title="Shuffled negatives for the i-motif architecture head",
        evidence_tier="derived",
        notes="Kept under their own source so the measurement-grounded im_fold head can "
              "select the real positives and their shuffles without pulling in 12,000 "
              "motif-derived rows.",
    )
    negs = shuffled.build([dict(sequence=c["sequence"], source_id=c["source_id"],
                                record_id=c["record_id"], ph=5.8)
                           for c in cands],
                          kind="iM", n_per_positive=1, reject_motif=True, seed=13)
    for r in negs:
        r.source = neg_source
    n += atlas.add(negs)
    print(f"  + im architecture positives + matched negatives: {n}")
    return n


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default="data/atlas.db")
    ap.add_argument("--g4stab-seqs", type=int, default=40_000,
                    help="number of distinct genomic sequences to import from the "
                         "G4STAB web database (each contributes 5 ionic conditions)")
    ap.add_argument("--im-candidates", type=int, default=12_000)
    ap.add_argument("--skip-predicted", action="store_true")
    args = ap.parse_args()

    Path(args.db).unlink(missing_ok=True)
    for ext in (".db-wal", ".db-shm"):
        Path(str(args.db).replace(".db", ext)).unlink(missing_ok=True)

    atlas = Atlas(args.db)
    print("Building QuadCond atlas ->", args.db)
    build_g4sp(atlas)
    build_g4stab_supp(atlas)
    build_imseeker(atlas)
    build_duvma(atlas)
    build_gse220882(atlas)
    build_gse220882_loci(atlas)
    if not args.skip_predicted:
        build_g4stab_predicted(atlas, args.g4stab_seqs)
        build_complementary_im(atlas, args.im_candidates)
        build_im_architecture_set(atlas)
    build_negatives(atlas)

    print("\nAtlas summary")
    print(atlas.summary().to_string(index=False))
    print(f"\n{atlas.composition()['one_line']}")
    atlas.close()


if __name__ == "__main__":
    main()
