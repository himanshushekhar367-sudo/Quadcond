"""Adapter: G4STAB web database (PREDICTED tier).

Liew, D., Dharmatilleke, A. D., See, E. & Yong, E. H. (2025) *G4STAB: a
multi-input deep learning model to predict G-quadruplex thermodynamic stability
based on sequence and salt concentration.* Bioinformatics 41(10):btaf545.
doi:10.1093/bioinformatics/btaf545
Files: ``data/pred_seqs_*.csv`` in github.com/donn-liew/g4stab-web-database

~1.9M rows: human genomic PQS scored by the G4STAB ensemble under four ionic
conditions, with genomic annotation (gene, phastCons, phyloP) and a per-row
standard error.

**These are model outputs, not measurements.**  They enter the atlas at the
``predicted`` tier and the training code refuses to use them as regression
targets unless you pass ``--allow-predicted``, in which case every report
labels the resulting head as distilled.  Their legitimate uses are (a) a
benchmark to compare a new model against, (b) a condition-response prior, and
(c) genomic context for candidate prioritisation.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterator

from ...conditions import Condition
from ..db import Record

SOURCE = "g4stab_webdb_predicted"
DOI = "10.1093/bioinformatics/btaf545"
URL = "https://donn-liew.github.io/g4stab-web-database/"


def iter_records(
    directory: str | Path,
    *,
    max_rows: int | None = None,
    keep_annotation: bool = True,
) -> Iterator[Record]:
    import pandas as pd

    directory = Path(directory)
    files = sorted(directory.glob("pred_seqs_*.csv"))
    if not files:
        raise FileNotFoundError(f"no pred_seqs_*.csv under {directory}")
    seen = 0
    for f in files:
        for chunk in pd.read_csv(f, chunksize=50_000):
            for _, r in chunk.iterrows():
                if max_rows is not None and seen >= max_rows:
                    return
                cond = Condition.from_mapping(
                    {
                        "k": r.get("k"),
                        "na": r.get("na"),
                        "li_nh4": r.get("li/nh4"),
                        "ph": 7.0,
                        "temperature": 25.0,
                    },
                    track_imputed=True,
                )
                genomic = None
                if keep_annotation:
                    genomic = {
                        "chrom": str(r.get("Chr")),
                        "start": int(r["Start"]) if pd.notna(r.get("Start")) else None,
                        "end": int(r["End"]) if pd.notna(r.get("End")) else None,
                        "strand": r.get("Strand"),
                        "gene": r.get("Symbol"),
                        "gene_type": r.get("Gene Type"),
                        "phastCons": float(r["phastCons"]) if pd.notna(r.get("phastCons")) else None,
                        "phyloP": float(r["phyloP"]) if pd.notna(r.get("phyloP")) else None,
                    }
                yield Record(
                    sequence=str(r["seq"]),
                    kind="G4",
                    source=SOURCE,
                    evidence_tier="predicted",
                    condition=cond,
                    tm=float(r["predicted_temperature"]),
                    method="G4STAB deep-ensemble prediction",
                    source_doi=DOI,
                    source_id=f"{r.get('Chr')}:{r.get('Start')}-{r.get('End')}{r.get('Strand')}",
                    organism="Homo sapiens (GRCh38)",
                    genomic=genomic,
                    qc_flags=["model_output_not_measurement", f"sem={r.get('sem')}"],
                )
                seen += 1


def register(atlas) -> None:
    atlas.register_source(
        SOURCE,
        title="G4STAB web database: predicted Tm for human genomic PQS under 4 ionic conditions",
        doi=DOI,
        url=URL,
        evidence_tier="predicted",
        notes="MODEL OUTPUT. Not experimental. Use for benchmarking / priors / genomic context only.",
    )
