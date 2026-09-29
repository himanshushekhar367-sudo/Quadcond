#!/usr/bin/env python3
"""Audit GSE296171 (Martyr et al., NAR 2026) before any of it reaches the atlas.

The dataset is a 6,000-member synthetic RNA pool assayed by RT-stop sequencing,
published with processed per-sequence scores.  It is attractive: sequence
resolved, condition-armed (KCl / LiCl / +PDS), and large.  Three properties make
it dangerous to ingest quickly, and this script measures all three without
writing to the atlas, the model weights or the frozen benchmark results.

1. **It is RNA, and the atlas is DNA.**  Collapsing U to T -- which is what the
   existing overlap audit does, correctly, for DNA collections -- makes 64 of
   these sequences string-identical to sequences the atlas already holds with
   measured melting temperatures.  A sequence-only featuriser cannot tell the
   two molecules apart.  That is a leakage path into `g4_tm`, whose evaluation
   set contains those same strings.

2. **The row count overstates the diversity.**  1,013 of the 1,603 strong
   folders are single-nucleotide mutants of six parent sequences.  Split at
   random, a model memorises six scaffolds and reports an excellent number.
   The published benchmark already found grouped-vs-random worth 0.15 R^2 on
   one head; here the scaffold structure is far more concentrated.

3. **The readout is a proxy, on a scale nothing else in the atlas uses.**  The
   RT-stop score is a log ratio against a 7-deazaguanine control, measured in
   reverse-transcriptase buffer (pH 8.3, 15 mM Mg2+, 40 C).  The atlas schema
   has no column for a continuous folding propensity and no ligand column for
   the pyridostatin arms, so ingestion requires schema changes that this script
   only reports -- it does not make them.

Run:

    python gse296171_overlap.py /path/to/stG4_STable_2_ProcessedData.xlsx

Exit status is 0 only when every check passes.  Any FAIL leaves status 1 and the
verdict says which gate stopped it.  Nothing is merged, nothing is retrained,
and no existing file is written.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
import re
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent.parent

SHEET_SCORES = "All Sequences Processed"
SHEET_MUT = "Mutation deltaWT Scores"

# Published composition (Martyr et al. 2026).  Reproducing these is the
# integrity gate: a workbook that does not match them is not the file we audited.
EXPECTED_ROWS = 6000
EXPECTED_MUT_ROWS = 1889
EXPECTED_POOLS = {"stG4": 3213, "nsG4_Mut": 1833, "neg_ns": 709,
                  "neg_AAA": 216, "pUG": 29}

COL_SEQ = "Shortened Sequences"
COL_KCL = "RT Stop Score in KCl (G in KCl/7dG in KCl)"
COL_LICL = "RT Stop Score in LiCl (G in LiCl/7dG in LiCl)"
COL_KCL_PDS = "RT Stop Score in KCl + PDS (G in KCl + PDS/7dG in KCl + PDS)"
REQUIRED_SCORE_COLS = [COL_KCL, COL_LICL, COL_KCL_PDS]

# More negative = more folded.  Asserted here as a hypothesis and tested in
# check_directionality; never assumed by any other check.
FOLDED_IS_NEGATIVE = True
STRONG_FOLDER_CUTOFF = -2.0

# Conditions for the primary arm, from the methods and the SI legend
# ("RT Stop Scores in KCl at 40 C").  Recorded so the schema mapping is explicit
# about what is reported versus what would be imputed on ingest.
PRIMARY_CONDITIONS = {
    # NOMINAL values are the stock/anneal concentrations quoted in the methods.
    # They are NOT what the RNA experiences: the 5X buffer is diluted into a
    # 22 uL reaction, so the in-reaction figures below are the ones a
    # condition-aware model must be given. Reporting the nominal Mg2+ would
    # overstate it by 5.5x, which is not a rounding difference.
    "nominal_stock": {"monovalent_mM": 100.0, "mg_mM": 15.0},
    "in_reaction": {
        # 375 mM * 4/22 (5X buffer) + 100 mM * 6/22 (annealed RNA)
        "monovalent_mM": round(375.0 * 4 / 22 + 100.0 * 6 / 22, 2),
        # 15 mM * 4/22; TOTAL Mg2+, not free -- dNTPs chelate some of it
        "mg_mM_total": round(15.0 * 4 / 22, 2),
        "tris_mM": round(250.0 * 4 / 22, 2),
    },
    "k_arm": "monovalent is K+; the in-house 5X substitutes the salt, so the "
             "LiCl arm carries no potassium",
    "ph": 8.3,
    "ph_caveat": "set on the Tris stock at room temperature; Tris at 40 C sits "
                 "roughly 0.4 units lower",
    "temperature_C": 40.0,
    "strand_conc": None,
    "strand_conc_note": "not reported for a 6,000-member pool; must be imputed",
    "reported_in": ("derived from the preprint methods' pipetting volumes "
                    "(10.1101/2025.04.29.651304) plus the SI Figure S1 legend "
                    "for the 40 C arm. NOT read from any data table, and not "
                    "yet checked against the published methods."),
}

SEQ_RNA_OK = re.compile(r"^[ACGU]+$")
SEQ_DNA_OK = re.compile(r"^[ACGT]+$")

# stG4_Pool_ns_<PARENT>_<1mut|AAAmut>_<n>
RE_MUT = re.compile(r"^stG4_Pool_ns_(?P<parent>[A-Za-z0-9]+)_(?P<kind>1mut|AAAmut)_\d+$")
RE_PARENT_ONLY = re.compile(r"^stG4_Pool_ns_(?P<parent>[A-Za-z0-9]+)$")


# --------------------------------------------------------------------------
# sequence keys: match on DNA, but never lose the molecule
# --------------------------------------------------------------------------
def rna_seq(raw) -> str | None:
    if raw is None:
        return None
    s = re.sub(r"\s+", "", str(raw)).upper()
    return s if SEQ_RNA_OK.match(s) else None


def dna_key(raw) -> str | None:
    """The key used ONLY to find cross-molecule collisions.

    This deliberately does what a sequence-only featuriser does -- erase the
    molecule -- because the question is whether the model can tell them apart,
    not whether a biologist can.
    """
    if raw is None:
        return None
    s = re.sub(r"\s+", "", str(raw)).upper().replace("U", "T")
    return s if SEQ_DNA_OK.match(s) else None


# --------------------------------------------------------------------------
# small stats, so scipy is optional
# --------------------------------------------------------------------------
def _ranks(xs: list[float]) -> list[float]:
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    out = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            out[order[k]] = avg
        i = j + 1
    return out


def pearson(xs: list[float], ys: list[float]) -> float | None:
    n = len(xs)
    if n < 3:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = math.sqrt(sum((a - mx) ** 2 for a in xs))
    syy = math.sqrt(sum((b - my) ** 2 for b in ys))
    return None if sxx == 0 or syy == 0 else sxy / (sxx * syy)


def spearman(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 3:
        return None
    return pearson(_ranks(xs), _ranks(ys))


def fisher_ci(r: float, n: int, z: float = 1.96) -> tuple[float, float] | None:
    """Confidence interval on a correlation, so it is reported with its uncertainty."""
    if r is None or n < 4 or abs(r) >= 1.0:
        return None
    zr = 0.5 * math.log((1 + r) / (1 - r))
    se = 1.0 / math.sqrt(n - 3)
    lo, hi = zr - z * se, zr + z * se
    return (math.tanh(lo), math.tanh(hi))


def median(xs: list[float]) -> float | None:
    if not xs:
        return None
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2.0


# --------------------------------------------------------------------------
# reading
# --------------------------------------------------------------------------
def read_sheet(path: Path, sheet: str) -> list[dict]:
    try:
        from openpyxl import load_workbook
    except ImportError:
        raise SystemExit("this audit needs openpyxl; pip install openpyxl")
    wb = load_workbook(path, read_only=True, data_only=True)
    if sheet not in wb.sheetnames:
        raise SystemExit(f"sheet {sheet!r} not in {path.name}; found {wb.sheetnames}")
    ws = wb[sheet]
    it = ws.iter_rows(values_only=True)
    header = [str(c).strip() if c is not None else "" for c in next(it)]
    rows = [dict(zip(header, r)) for r in it]
    wb.close()
    return [r for r in rows if any(v is not None for v in r.values())]


def pool_column(rows: list[dict]) -> str:
    for k in rows[0]:
        if k.startswith("Pool Identifier"):
            return k
    raise SystemExit("no 'Pool Identifier' column found")


def as_float(v):
    try:
        f = float(v)
        return None if math.isnan(f) else f
    except (TypeError, ValueError):
        return None


# --------------------------------------------------------------------------
# check 1 -- input integrity
# --------------------------------------------------------------------------
def check_input_integrity(rows, mut_rows, pool_col) -> dict:
    seqs = [rna_seq(r.get(COL_SEQ)) for r in rows]
    bad = sum(1 for s in seqs if s is None)
    good = [s for s in seqs if s]
    pools = Counter(str(r.get(pool_col)) for r in rows)
    missing_cols = [c for c in REQUIRED_SCORE_COLS if c not in rows[0]]
    missing_vals = {c: sum(1 for r in rows if as_float(r.get(c)) is None)
                    for c in REQUIRED_SCORE_COLS if c in rows[0]}
    lens = [len(s) for s in good]
    res = {
        "rows": len(rows),
        "rows_expected": EXPECTED_ROWS,
        "mutation_rows": len(mut_rows),
        "mutation_rows_expected": EXPECTED_MUT_ROWS,
        "sequences_parsed_as_rna": len(good),
        "sequences_unparsable": bad,
        "distinct_sequences": len(set(good)),
        "duplicate_sequences": len(good) - len(set(good)),
        "length_min": min(lens) if lens else None,
        "length_max": max(lens) if lens else None,
        "pool_composition": dict(pools),
        "pool_composition_expected": EXPECTED_POOLS,
        "required_score_columns_missing": missing_cols,
        "missing_values_by_score_column": missing_vals,
    }
    fails = []
    if len(rows) != EXPECTED_ROWS:
        fails.append(f"row count {len(rows)} != published {EXPECTED_ROWS}")
    if len(mut_rows) != EXPECTED_MUT_ROWS:
        fails.append(f"mutation row count {len(mut_rows)} != published {EXPECTED_MUT_ROWS}")
    if bad:
        fails.append(f"{bad} sequences are not plain ACGU RNA")
    if len(good) != len(set(good)):
        fails.append("sequences are not unique")
    if missing_cols:
        fails.append(f"missing score columns: {missing_cols}")
    if any(v for v in missing_vals.values()):
        fails.append(f"missing values in score columns: {missing_vals}")
    if dict(pools) != EXPECTED_POOLS:
        fails.append("pool composition does not reproduce the published counts")
    res["status"] = "PASS" if not fails else "FAIL"
    res["failures"] = fails
    return res


# --------------------------------------------------------------------------
# check 2 -- directionality, established from controls only
# --------------------------------------------------------------------------
def check_directionality(rows, pool_col) -> dict:
    by_pool = defaultdict(list)
    for r in rows:
        v = as_float(r.get(COL_KCL))
        if v is not None:
            by_pool[str(r.get(pool_col))].append(v)

    neg = median(by_pool.get("neg_AAA", []))
    nat = median(by_pool.get("nsG4_Mut", []))

    # G-tract monotonicity, within the synthetic pool only (no natural sequences,
    # so scaffold identity cannot drive it).
    tract_col = "Minimum G-Tract (2 = GG, 3 = GGG, 4 = GGGG)"
    by_tract = defaultdict(list)
    for r in rows:
        if str(r.get(pool_col)) != "stG4":
            continue
        t, v = as_float(r.get(tract_col)), as_float(r.get(COL_KCL))
        if t is not None and v is not None:
            by_tract[int(t)].append(v)
    tract_med = {k: median(v) for k, v in sorted(by_tract.items())}

    # A G4 stabiliser should push folders further in the folded direction.
    pds_shift = median([as_float(r.get(COL_KCL_PDS)) - as_float(r.get(COL_KCL))
                        for r in rows
                        if str(r.get(pool_col)) == "nsG4_Mut"
                        and as_float(r.get(COL_KCL_PDS)) is not None
                        and as_float(r.get(COL_KCL)) is not None])
    # Li+ does not support G4; the signal should collapse relative to K+.
    licl_nat = median([as_float(r.get(COL_LICL)) for r in rows
                       if str(r.get(pool_col)) == "nsG4_Mut"
                       and as_float(r.get(COL_LICL)) is not None])

    tracts = sorted(tract_med)
    monotone = all(tract_med[a] > tract_med[b] for a, b in zip(tracts, tracts[1:])) \
        if len(tracts) > 1 else False

    ev = {
        "median_kcl_negative_control_neg_AAA": neg,
        "median_kcl_natural_g4_nsG4_Mut": nat,
        "median_kcl_by_min_g_tract_synthetic_only": tract_med,
        "g_tract_monotone_toward_negative": monotone,
        "median_pds_shift_on_natural_g4": pds_shift,
        "median_licl_natural_g4": licl_nat,
        "licl_signal_collapses_vs_kcl": (licl_nat is not None and nat is not None
                                         and abs(licl_nat) < abs(nat)),
    }
    tests = {
        "non-folding AAA control is less negative than natural G4":
            neg is not None and nat is not None and neg > nat,
        "longer minimum G-tract shifts score negative": monotone,
        "pyridostatin shifts natural G4 further negative":
            pds_shift is not None and pds_shift < 0,
        "K+ produces a stronger signal than Li+":
            bool(ev["licl_signal_collapses_vs_kcl"]),
    }
    passed = all(tests.values())
    return {
        "hypothesis": "more negative RT-stop score = more folded",
        "evidence": ev,
        "tests": tests,
        "independent_tests_passed": sum(tests.values()),
        "independent_tests_total": len(tests),
        "status": "PASS" if passed else "FAIL",
        "failures": [] if passed else
            [k for k, v in tests.items() if not v],
        "note": ("Direction is established from controls and condition shifts only. "
                 "No atlas label is consulted, so this check is independent of the "
                 "cross-molecule concordance reported separately."),
    }


# --------------------------------------------------------------------------
# check 3 -- DNA/RNA overlap against the atlas
# --------------------------------------------------------------------------
def open_atlas_readonly(path: Path) -> sqlite3.Connection:
    """Read-only connection: the audit must not be able to alter the atlas."""
    return sqlite3.connect(f"file:{path}?mode=ro", uri=True)


def check_dna_rna_overlap(rows, pool_col, atlas_path: Path) -> tuple[dict, list[dict]]:
    con = open_atlas_readonly(atlas_path)
    atlas = defaultdict(list)
    q = ("select sequence, kind, nucleic_acid, label_class, evidence_tier, source, "
         "method, tm, folded, k, na, ph, temperature from records")
    for row in con.execute(q):
        (seq, kind, na_mol, lab, tier, src, meth, tm, folded, k, na, ph, temp) = row
        key = dna_key(seq)
        if key:
            atlas[key].append({
                "atlas_sequence": seq, "atlas_molecule": na_mol or "DNA",
                "kind": kind, "label_class": lab, "evidence_tier": tier,
                "source": src, "method": meth, "tm": tm, "folded": folded,
                "k": k, "na": na, "ph": ph, "temperature": temp,
            })
    con.close()

    detail: list[dict] = []
    seq_hits: set[str] = set()
    by_source = Counter()
    tm_pairs_x: list[float] = []
    tm_pairs_y: list[float] = []
    for r in rows:
        rna = rna_seq(r.get(COL_SEQ))
        key = dna_key(r.get(COL_SEQ))
        if not key or key not in atlas:
            continue
        seq_hits.add(key)
        score = as_float(r.get(COL_KCL))
        tms = [a["tm"] for a in atlas[key] if a["tm"] is not None]
        if tms and score is not None:
            tm_pairs_x.append(median(tms))
            tm_pairs_y.append(score)
        for a in atlas[key]:
            by_source[f"{a['source']}::{a['label_class']}"] += 1
            detail.append({
                "rna_sequence": rna,
                "dna_match_key": key,
                "rna_molecule": "RNA",
                "atlas_sequence": a["atlas_sequence"],
                "atlas_molecule": a["atlas_molecule"],
                "pool_class": str(r.get(pool_col)),
                "rna_rt_stop_kcl": score,
                "rna_rt_stop_licl": as_float(r.get(COL_LICL)),
                "atlas_kind": a["kind"],
                "atlas_source": a["source"],
                "atlas_label_class": a["label_class"],
                "atlas_evidence_tier": a["evidence_tier"],
                "atlas_method": a["method"],
                "atlas_tm": a["tm"],
                "atlas_folded": a["folded"],
                "atlas_k_mM": a["k"], "atlas_na_mM": a["na"],
                "atlas_ph": a["ph"], "atlas_temperature_C": a["temperature"],
            })

    rho = spearman(tm_pairs_x, tm_pairs_y)
    r = pearson(tm_pairs_x, tm_pairs_y)
    ci = fisher_ci(rho, len(tm_pairs_x)) if rho is not None else None
    summary = {
        "atlas_file": atlas_path.name,
        "atlas_opened_readonly": True,
        "distinct_sequences_overlapping": len(seq_hits),
        "atlas_rows_touched": len(detail),
        "overlap_by_atlas_source": dict(by_source),
        "matched_tm_pairs": len(tm_pairs_x),
        "cross_molecule_concordance": {
            "spearman_rho": round(rho, 4) if rho is not None else None,
            "pearson_r": round(r, 4) if r is not None else None,
            "rho_95pct_ci": [round(c, 4) for c in ci] if ci else None,
            "n": len(tm_pairs_x),
            "interpretation": (
                "ASSOCIATION ONLY. This is a cross-molecule, cross-assay concordance "
                "between DNA melting temperature and an RNA reverse-transcriptase stop "
                "ratio measured in a different buffer. It does NOT establish that the "
                "RT-stop score measures thermodynamic stability, and must not be "
                "reported as if it did. Causal or thermodynamic interpretation belongs "
                "in the validation report, not in this audit."),
        },
        "status": "INFO",
        "note": ("Overlap is not itself a failure. It becomes one at the leakage gate "
                 "if any of these sequences is also a benchmark evaluation sequence."),
    }
    return summary, detail


# --------------------------------------------------------------------------
# check 4 -- leakage, fail-closed
# --------------------------------------------------------------------------
def load_benchmark_sets() -> tuple[dict[str, set[str]], list[str]]:
    out: dict[str, set[str]] = {}
    problems: list[str] = []
    res = REPO / "benchmarks" / "published_tools" / "results"
    tool = res / "tool_scores_all.csv"
    if tool.exists():
        s = set()
        with tool.open(encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh):
                k = dna_key(row.get("sequence", ""))
                if k:
                    s.add(k)
        out["tool_scores_all"] = s
    else:
        problems.append(f"missing benchmark evaluation set: {tool}")
    pairs = res / "pairs_g4_tm.csv"
    if pairs.exists():
        s = set()
        with pairs.open(encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh):
                for col in ("ref", "alt"):
                    k = dna_key(row.get(col, ""))
                    if k:
                        s.add(k)
        out["pairs_g4_tm"] = s
    else:
        problems.append(f"missing dTm pair set: {pairs}")
    return out, problems


SHIPPED_MODEL_CARDS = (
    Path("manuscript") / "inputs" / "quadcond_model.json",
    Path("manuscript") / "inputs" / "quadcond_model_seqonly.json",
)


def shipped_fingerprints() -> dict[str, str]:
    """dataset_fingerprint_sha256 of each shipped model, without unpickling it."""
    out: dict[str, str] = {}
    for rel in SHIPPED_MODEL_CARDS:
        card = REPO / rel
        if not card.exists():
            continue
        try:
            d = json.loads(card.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        fp = d.get("dataset_fingerprint_sha256") or d.get("dataset_fingerprint")
        if fp:
            out[rel.name] = fp
    return out


def _validate_fold_manifest(path: Path) -> tuple[Path | None, list[str]]:
    """A manifest is only usable if it describes the rows the model was trained on.

    The existence of a file called ``fold_assignments.csv`` proves nothing. One
    built from a different atlas -- the trimmed ``atlas_core.db`` rather than the
    full ``atlas.db``, say -- describes a different split over different rows, and
    using it to clear a shipped head of leakage would be worse than having no
    manifest at all: it would return a confident, wrong all-clear. So the
    fingerprint is checked against the shipped model cards, and a mismatch keeps
    the gate closed.
    """
    meta_path = path.with_name(path.stem + "_meta.json")
    if not meta_path.exists():
        return None, [
            f"Found {path.name} but no {meta_path.name} beside it. Without the "
            "metadata the manifest's provenance cannot be checked, so it cannot "
            "be used to clear any sequence."]
    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as e:
        return None, [f"{meta_path.name} is unreadable: {e}"]

    shipped = shipped_fingerprints()
    if not shipped:
        return None, [
            "No shipped model card found, so the manifest cannot be tied to a "
            "trained model. Searched: "
            + ", ".join(str(r) for r in SHIPPED_MODEL_CARDS)]

    # One definition of the three states, in quadcond.models.foldmanifest, used
    # by both the writer and this reader. A second implementation here would
    # drift, and a leakage gate reading a drifted definition is the failure this
    # whole file exists to prevent.
    sys.path.insert(0, str(REPO))
    try:
        from quadcond.models.foldmanifest import classify_manifest
    except ImportError as e:
        return None, [
            f"could not import the manifest classifier ({e}); failing closed "
            "rather than guessing which state this manifest is in"]

    verdict = classify_manifest(meta, shipped_fingerprints=list(shipped.values()))
    if not verdict["usable_to_clear_a_published_number"]:
        problems = [
            f"Fold manifest is {verdict['state'].upper()} "
            f"(claimed: {verdict['claimed_provenance']!r}).",
            *verdict["reasons"],
        ]
        if verdict["state"] == "stale":
            problems.append(
                "A stale manifest still names every head and still looks complete. "
                "Most often it was built against data/atlas_core.db instead of the "
                "full data/atlas.db, or the split code changed after it was "
                "written. Rebuild it against the atlas the model was trained on, "
                "with the current code.")
        else:
            problems.append(
                "To obtain a historical manifest, have the training run write the "
                "assignment as it trains. Until then the correct status is BLOCKED.")
        return None, problems
    return path, []


def locate_fold_manifest() -> tuple[Path | None, list[str]]:
    """Find a persisted train/validation/held-out assignment.

    There is deliberately no fallback here.  Folds in this project are built at
    training time by GroupKFold over cluster_sequences and are not written to
    disk, so 'which fold is this sequence in' cannot be answered from any
    committed artifact.  Guessing would defeat the purpose of the gate, so the
    audit fails closed and says what would have to exist.
    """
    candidates = [
        REPO / "data" / "fold_assignments.csv",
        REPO / "data" / "splits.json",
        REPO / "benchmarks" / "published_tools" / "results" / "fold_assignments.csv",
        REPO / "benchmarks" / "ingest_audit" / "results" / "fold_assignments.csv",
    ]
    # Tracked location first: data/ is gitignored, so a manifest there cannot be
    # committed and cannot be checked by a reviewer.
    candidates.sort(key=lambda c: 0 if "published_tools" in str(c) else 1)
    for c in candidates:
        if c.exists():
            return _validate_fold_manifest(c)

    def _rel(p: Path) -> str:
        try:
            return str(p.relative_to(REPO))
        except ValueError:      # a candidate outside the repo must not crash the report
            return str(p)

    return None, [
        "No persisted fold assignment manifest found. Searched: "
        + ", ".join(_rel(c) for c in candidates),
        "Folds are constructed at train time (GroupKFold over "
        "quadcond.models.grouping.cluster_sequences) and never serialised, so no "
        "committed artifact records which sequences were held out.",
    ]


def check_leakage(rows, overlap_detail) -> dict:
    cand_keys = {dna_key(r.get(COL_SEQ)) for r in rows}
    cand_keys.discard(None)
    bench, bench_problems = load_benchmark_sets()
    fold_path, fold_problems = locate_fold_manifest()

    per_set = {}
    prohibited = 0
    for name, s in bench.items():
        hit = cand_keys & s
        per_set[name] = {
            "benchmark_sequences": len(s),
            "also_in_candidate": len(hit),
            "fraction_of_benchmark_covered": round(len(hit) / len(s), 6) if s else None,
            "example_collisions": sorted(hit)[:10],
        }
        prohibited += len(hit)

    overlap_keys = {d["dna_match_key"] for d in overlap_detail}
    bench_union = set().union(*bench.values()) if bench else set()
    atlas_and_bench = overlap_keys & bench_union

    failures = list(bench_problems)
    status = "PASS"
    if fold_path is None:
        status = "FAIL_CLOSED"
        failures.extend(fold_problems)
    if bench_problems:
        status = "FAIL_CLOSED"
    if prohibited:
        status = "FAIL"

    return {
        "fold_manifest": str(fold_path) if fold_path else None,
        "fold_assignments_established": fold_path is not None,
        "benchmark_overlap": per_set,
        "prohibited_overlap_sequences": prohibited,
        "sequences_in_both_atlas_and_benchmark": len(atlas_and_bench),
        "status": status,
        "failures": failures,
        "required_action": (
            "Train/validation/held-out membership cannot be established from any "
            "committed artifact, so no sequence in this dataset may be declared "
            "training-eligible by this audit. Either serialise the fold assignment "
            "used for the frozen benchmark and rerun, or quarantine every sequence "
            "listed in gse296171_dna_rna_overlap.csv from training and rerun the "
            "benchmark with them excluded."
            if status != "PASS" else
            "No prohibited overlap detected between the candidate and the benchmark "
            "evaluation sets."),
    }


# --------------------------------------------------------------------------
# check 5 -- scaffold structure
# --------------------------------------------------------------------------
def assign_scaffold(pool_id: str, pool_class: str) -> tuple[str, str]:
    """Return (group_key, assignment_basis).

    Synthetic combinatorial designs get one group each and are flagged: they have
    no parent scaffold, so grouping them by scaffold is meaningless and they need
    sequence-similarity clustering instead.  Anything that matches no rule is
    UNASSIGNED and is counted as a failure rather than quietly given its own
    group, which would silently license it as independent.
    """
    pid = str(pool_id or "")
    m = RE_MUT.match(pid)
    if m:
        return f"parent:{m.group('parent')}", f"mutant_of_parent ({m.group('kind')})"
    m = RE_PARENT_ONLY.match(pid)
    if m:
        return f"parent:{m.group('parent')}", "parent_wildtype"
    if pool_class == "stG4" or pid.startswith("stG4_Pool_GLV"):
        return f"synthetic:{pid}", "synthetic_combinatorial_no_parent"
    if pool_class in ("neg_AAA", "neg_ns", "pUG") or pid.startswith(
            ("stG4_Pool_neg_AAA", "stG4_Pool_neg_ns", "stG4_Pool_pUG")):
        return f"control:{pool_class}:{pid}", f"control_{pool_class}"
    return "UNASSIGNED", "no rule matched"


def check_scaffolds(rows, mut_rows, pool_col) -> tuple[dict, list[dict]]:
    # The mutation sheet is the authority on which sequence is a mutant of what.
    mut_by_seq: dict[str, str] = {}
    mut_counts: Counter = Counter()
    for m in mut_rows:
        s = rna_seq(m.get(COL_SEQ))
        pid = str(m.get("Pool identifier") or "")
        if s:
            mut_by_seq[s] = pid
        g, _ = assign_scaffold(pid, "nsG4_Mut")
        mut_counts[g] += 1

    table: list[dict] = []
    groups: Counter = Counter()
    basis_counts: Counter = Counter()
    unassigned: list[str] = []
    strong_by_group: Counter = Counter()

    for r in rows:
        s = rna_seq(r.get(COL_SEQ))
        cls = str(r.get(pool_col))
        pid = mut_by_seq.get(s or "", "")
        if not pid:
            # rows outside the mutation sheet carry no explicit pool id column,
            # so fall back to the class, which the integrity check validated.
            pid = f"stG4_Pool_GLV:{s}" if cls == "stG4" else f"stG4_Pool_{cls}:{s}"
        group, basis = assign_scaffold(pid, cls)
        if group == "UNASSIGNED":
            unassigned.append(s or "<unparsable>")
        groups[group] += 1
        basis_counts[basis] += 1
        score = as_float(r.get(COL_KCL))
        strong = score is not None and score < STRONG_FOLDER_CUTOFF
        if strong:
            strong_by_group[group] += 1
        table.append({
            "rna_sequence": s, "pool_class": cls, "pool_identifier": pid,
            "group_key": group, "assignment_basis": basis,
            "rt_stop_kcl": score,
            "is_strong_folder": int(bool(strong)),
        })

    parent_groups = {g: n for g, n in groups.items() if g.startswith("parent:")}
    strong_total = sum(strong_by_group.values())
    strong_in_parents = sum(n for g, n in strong_by_group.items() if g.startswith("parent:"))

    summary = {
        "distinct_group_keys": len(groups),
        "parent_scaffolds": len(parent_groups),
        "parent_scaffold_names": sorted(parent_groups),
        "sequences_per_parent_scaffold": dict(sorted(parent_groups.items())),
        "assignment_basis_counts": dict(basis_counts),
        "unassigned_sequences": len(unassigned),
        "unassigned_examples": unassigned[:10],
        "strong_folders_total": strong_total,
        "strong_folders_from_parent_scaffolds": strong_in_parents,
        "strong_folder_concentration": (
            round(strong_in_parents / strong_total, 4) if strong_total else None),
        "effective_independent_n_strong_regime": len(
            {g for g in strong_by_group if g.startswith("parent:")}) or None,
        "status": "FAIL" if unassigned else "PASS",
        "failures": ([f"{len(unassigned)} sequences could not be assigned to a scaffold "
                      "group; they must not be treated as independent observations"]
                     if unassigned else []),
        "note": ("Synthetic combinatorial designs have no parent scaffold and are given "
                 "one group each. That is correct for provenance but NOT sufficient for "
                 "cross-validation: they are combinatorial variants of a shared design "
                 "grammar and still need sequence-similarity clustering "
                 "(quadcond.models.grouping.cluster_sequences) before use as folds."),
    }
    return summary, table


# --------------------------------------------------------------------------
# check 6 -- what the two split schemes actually do
# --------------------------------------------------------------------------
def check_splits(scaffold_table, n_splits: int = 5, seed: int = 0) -> dict:
    """Composition and leakage diagnostics for random vs grouped splits.

    No model is fitted and no label is consulted to choose anything.  The
    diagnostic is structural: under a random split, how many held-out sequences
    have a sibling from the same scaffold in the training half?  That number is
    the size of the illusion, and it needs no estimator to compute.
    """
    rows = [r for r in scaffold_table if r["rna_sequence"]]
    rnd = random.Random(seed)
    idx = list(range(len(rows)))
    rnd.shuffle(idx)

    # random split: assign each row a fold directly from the shuffled order
    random_fold = {}
    for pos, i in enumerate(idx):
        random_fold[i] = pos % n_splits

    # grouped split: whole scaffolds move together, greedily balanced
    by_group = defaultdict(list)
    for i, r in enumerate(rows):
        by_group[r["group_key"]].append(i)
    sizes = sorted(by_group.items(), key=lambda kv: -len(kv[1]))
    load = [0] * n_splits
    grouped_fold = {}
    for g, members in sizes:
        f = load.index(min(load))
        for i in members:
            grouped_fold[i] = f
        load[f] += len(members)

    def leakage(fold_map):
        """Held-out rows whose scaffold also appears in training."""
        per_fold = []
        for f in range(n_splits):
            test = [i for i in range(len(rows)) if fold_map[i] == f]
            train_groups = {rows[i]["group_key"] for i in range(len(rows))
                            if fold_map[i] != f}
            bleed = sum(1 for i in test if rows[i]["group_key"] in train_groups)
            strong_test = sum(rows[i]["is_strong_folder"] for i in test)
            strong_bleed = sum(rows[i]["is_strong_folder"] for i in test
                               if rows[i]["group_key"] in train_groups)
            per_fold.append({
                "fold": f, "test_n": len(test),
                "test_rows_with_scaffold_sibling_in_train": bleed,
                "fraction_contaminated": round(bleed / len(test), 4) if test else None,
                "strong_folders_in_test": strong_test,
                "strong_folders_contaminated": strong_bleed,
            })
        tot = sum(p["test_rows_with_scaffold_sibling_in_train"] for p in per_fold)
        strong_tot = sum(p["strong_folders_contaminated"] for p in per_fold)
        return per_fold, tot, strong_tot

    r_folds, r_tot, r_strong = leakage(random_fold)
    g_folds, g_tot, g_strong = leakage(grouped_fold)

    return {
        "n_splits": n_splits,
        "seed": seed,
        "rows_considered": len(rows),
        "random_sequence_level_split": {
            "per_fold": r_folds,
            "total_test_rows_with_scaffold_sibling_in_train": r_tot,
            "fraction_of_all_rows": round(r_tot / len(rows), 4) if rows else None,
            "strong_folders_contaminated": r_strong,
        },
        "parent_scaffold_grouped_split": {
            "per_fold": g_folds,
            "total_test_rows_with_scaffold_sibling_in_train": g_tot,
            "fraction_of_all_rows": round(g_tot / len(rows), 4) if rows else None,
            "strong_folders_contaminated": g_strong,
        },
        "leakage_removed_by_grouping": r_tot - g_tot,
        "status": "PASS" if g_tot == 0 else "FAIL",
        "failures": ([] if g_tot == 0 else
                     [f"grouped split still leaks {g_tot} rows; grouping key is wrong"]),
        "note": ("No estimator was fitted and no held-out label was read. These are "
                 "structural counts, so they cannot be inflated by model selection. "
                 "The performance consequence of the difference must be measured by "
                 "the training pipeline, not asserted here."),
    }


# --------------------------------------------------------------------------
# check 7 -- schema readiness
# --------------------------------------------------------------------------
def check_schema_readiness() -> dict:
    from importlib import util as _u
    spec = _u.spec_from_file_location("atlas_schema", REPO / "quadcond" / "atlas" / "schema.py")
    mod = _u.module_from_spec(spec)
    spec.loader.exec_module(mod)
    existing_labels = set(getattr(mod, "LABEL_CLASSES", ()))
    ddl = getattr(mod, "DDL", "")

    # Scope note (2026-09-28): this gate previously also demanded `ligand`,
    # `ligand_conc` and a new label_class. Reviewing the RNA adapter showed all
    # three were over-strict. Ligand fields are needed only if the +PDS arms are
    # ingested, and an adapter may simply decline them. And `biophysical` is
    # correct by schema.py's own definition -- a structure measured in a defined
    # buffer, as against `genomic_proxy` for in-cell assays -- so RT-stop does
    # not need a class of its own. The concern that motivated one was
    # APPLICABILITY DOMAIN, not labelling: pH 8.3 with ~2.7 mM Mg2+ is a region
    # no other atlas row occupies. That belongs to the claims module as a fence,
    # and is reported below rather than blocking the schema.
    required = [
        {"field": "rt_stop_score", "type": "REAL", "present": "rt_stop_score" in ddl,
         "blocking": True,
         "why": ("The measurement is a continuous log-ratio folding propensity. The "
                 "schema's label columns are folded (integer), tm, dg and ph_t; none "
                 "of them is this. An adapter may threshold into `folded` using the "
                 "authors' own cut-off and keep the score as a QC flag, which is "
                 "valid but discards the gradation that makes the dataset valuable.")},
    ]
    conditional = [
        {"field": "ligand", "type": "TEXT", "present": "ligand" in ddl,
         "required_only_if": "the +PDS arms are ingested",
         "why": ("PDS rows sit at the same condition as their ligand-free twin with a "
                 "different label -- the contradiction a condition-aware model must "
                 "never be shown. An adapter that excludes them needs no ligand "
                 "field.")},
        {"field": "ligand_conc", "type": "REAL", "present": "ligand_conc" in ddl,
         "required_only_if": "the +PDS arms are ingested",
         "why": "A ligand identity without its concentration is not a condition."},
    ]
    missing = [r["field"] for r in required if not r["present"]]
    return {
        "schema_version": getattr(mod, "SCHEMA_VERSION", None),
        "existing_label_classes": sorted(existing_labels),
        "nucleic_acid_column_present": "nucleic_acid" in ddl,
        "required_fields": required,
        "conditionally_required_fields": conditional,
        "missing_fields": missing,
        "label_class_verdict": (
            "biophysical is correct for an in-buffer structural measurement; RT-stop "
            "does not need a label class of its own. Carry the readout on the row "
            "(measurement=enzymatic_proxy_rt_stop) so the distinction survives."),
        "applicability_domain_fence_required": {
            "reason": ("pH 8.3 and ~2.7 mM total Mg2+ at 40 C is a condition region no "
                       "other atlas row occupies. Pooling it with CD and UV melting near "
                       "pH 7 at low Mg2+ would silently extend every head's claimed "
                       "domain."),
            "enforced_by": "quadcond.claims, not the schema",
            "blocking_here": False,
        },
        "primary_arm_conditions": PRIMARY_CONDITIONS,
        "status": "BLOCKED" if missing else "READY",
        "note": ("This check reports only. It does not modify schema.py, and the "
                 "schema must stay frozen until the overlap manifest is reviewed."),
    }


# --------------------------------------------------------------------------
def verdict(rep: dict) -> list[str]:
    out = []
    lk = rep["leakage"]
    if not lk["fold_assignments_established"]:
        out.append(
            "FAIL CLOSED: no committed artifact records train/validation/held-out "
            "membership, so this audit cannot certify any sequence as training-eligible. "
            "Serialise the fold assignment behind the frozen benchmark, or quarantine "
            "the overlap manifest from training and rerun the benchmark without it.")
    if lk["prohibited_overlap_sequences"]:
        out.append(
            f"LEAKAGE: {lk['prohibited_overlap_sequences']} candidate sequences are also "
            "benchmark evaluation sequences. Ingesting them as training data makes the "
            "published comparison partly in-sample.")
    sc = rep["scaffolds"]
    if sc.get("strong_folder_concentration"):
        out.append(
            f"{sc['strong_folder_concentration']:.1%} of strong folders come from "
            f"{sc['parent_scaffolds']} parent scaffolds. Describe this dataset by its "
            "scaffold count, not by its row count, and group folds on the scaffold.")
    sp = rep["splits"]
    out.append(
        f"A random split leaks {sp['random_sequence_level_split']['total_test_rows_with_scaffold_sibling_in_train']} "
        f"test rows to their own scaffold; grouping removes {sp['leakage_removed_by_grouping']}.")
    sr = rep["schema_readiness"]
    if sr["missing_fields"]:
        out.append(
            "SCHEMA BLOCKED: ingestion requires " + ", ".join(sr["missing_fields"]) +
            ". The schema stays frozen until the overlap manifest is reviewed.")
    out.append(
        "The DNA/RNA correlation is a cross-molecule concordance, not evidence that "
        "RT-stop measures thermodynamic stability. Keep that distinction in the "
        "validation report.")
    return out


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in fields})


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("workbook", help="stG4_STable_2_ProcessedData.xlsx")
    ap.add_argument("--atlas", default=None)
    ap.add_argument("--out", default=str(ROOT / "results"))
    ap.add_argument("--n-splits", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    path = Path(a.workbook)
    if not path.exists():
        raise SystemExit(f"not found: {path}")

    atlas_path = Path(a.atlas) if a.atlas else None
    atlas_is_subset = False
    if atlas_path is None:
        full, core = REPO / "data" / "atlas.db", REPO / "data" / "atlas_core.db"
        atlas_path = full if full.exists() else core
        atlas_is_subset = atlas_path.name == "atlas_core.db"
        if atlas_is_subset:
            print("WARNING: comparing against atlas_core.db, a subset of the full atlas. "
                  "Overlap below is a LOWER BOUND. Fetch data/atlas.db and rerun before "
                  "acting on these numbers.\n", file=sys.stderr)
    if not atlas_path.exists():
        raise SystemExit(f"no atlas at {atlas_path}")

    rows = read_sheet(path, SHEET_SCORES)
    mut_rows = read_sheet(path, SHEET_MUT)
    pool_col = pool_column(rows)

    integrity = check_input_integrity(rows, mut_rows, pool_col)
    direction = check_directionality(rows, pool_col)
    overlap, overlap_detail = check_dna_rna_overlap(rows, pool_col, atlas_path)
    leakage = check_leakage(rows, overlap_detail)
    scaffolds, scaffold_table = check_scaffolds(rows, mut_rows, pool_col)
    splits = check_splits(scaffold_table, n_splits=a.n_splits, seed=a.seed)
    schema = check_schema_readiness()

    report = {
        "dataset": "GSE296171 / Martyr et al., Nucleic Acids Research 54(16), 2026",
        "doi": "10.1093/nar/gkag821",
        "geo_accession": "GSE296171",
        "figshare": "10.6084/m9.figshare.33132299",
        "licence": "CC BY-NC",
        "workbook": str(path),
        "atlas_compared_against": atlas_path.name,
        "atlas_is_full": not atlas_is_subset,
        "overlap_is_lower_bound": atlas_is_subset,
        "input_integrity": integrity,
        "directionality": direction,
        "dna_rna_overlap": overlap,
        "leakage": leakage,
        "scaffolds": scaffolds,
        "splits": splits,
        "schema_readiness": schema,
        "nothing_was_modified": {
            "atlas": "opened read-only",
            "model_weights": "not opened",
            "frozen_benchmark_results": "read only",
            "schema": "read only, not modified",
        },
    }
    gates = {
        "input_integrity": integrity["status"],
        "directionality": direction["status"],
        "leakage": leakage["status"],
        "scaffolds": scaffolds["status"],
        "splits": splits["status"],
        "schema_readiness": schema["status"],
    }
    report["gates"] = gates
    report["audit_status"] = (
        "PASS" if all(v in ("PASS", "READY", "INFO") for v in gates.values()) else "BLOCKED")
    report["verdict"] = verdict(report)

    outdir = Path(a.out)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "gse296171_overlap_summary.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    (outdir / "gse296171_split_audit.json").write_text(
        json.dumps(splits, indent=2), encoding="utf-8")
    write_csv(outdir / "gse296171_dna_rna_overlap.csv", overlap_detail, [
        "rna_sequence", "rna_molecule", "dna_match_key", "pool_class",
        "rna_rt_stop_kcl", "rna_rt_stop_licl",
        "atlas_sequence", "atlas_molecule", "atlas_kind", "atlas_source",
        "atlas_label_class", "atlas_evidence_tier", "atlas_method",
        "atlas_tm", "atlas_folded", "atlas_k_mM", "atlas_na_mM",
        "atlas_ph", "atlas_temperature_C"])
    write_csv(outdir / "gse296171_scaffold_groups.csv", scaffold_table, [
        "rna_sequence", "pool_class", "pool_identifier", "group_key",
        "assignment_basis", "rt_stop_kcl", "is_strong_folder"])

    print(json.dumps({k: report[k] for k in
                      ("dataset", "audit_status", "gates", "verdict")}, indent=2))
    print(f"\nwritten to {outdir}:")
    for f in ("gse296171_overlap_summary.json", "gse296171_dna_rna_overlap.csv",
              "gse296171_scaffold_groups.csv", "gse296171_split_audit.json"):
        print(f"  {f}")
    print("\nNothing was merged, nothing was retrained, the schema was not modified.")
    return 0 if report["audit_status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
