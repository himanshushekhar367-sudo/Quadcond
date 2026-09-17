"""Shared helpers for the genome-wide QuadCond x AlphaGenome analysis."""
from __future__ import annotations

import gzip
import os
from pathlib import Path

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def read_chrom(fasta: str, chrom: str) -> str:
    """One chromosome from an indexed FASTA (.fai), upper-cased (soft-mask removed)."""
    import pysam
    with pysam.FastaFile(fasta) as fa:
        names = set(fa.references)
        name = chrom if chrom in names else chrom.replace("chr", "")
        if name not in names:
            raise KeyError(f"{chrom} not in {fasta}")
        return fa.fetch(name).upper()


def out_dir(env: str = "GW_OUT") -> Path:
    p = Path(os.environ.get(env, "gw_out")).expanduser()
    p.mkdir(parents=True, exist_ok=True)
    return p


def open_text(path: str | Path, mode: str = "rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)
