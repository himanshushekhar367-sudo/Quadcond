"""Adapter: iMab / BG4 CUT&Tag peaks in live human cells (GSE220882).

Zanin, I., Ruggiero, E., Nicoletto, G., Lago, S., Maurizio, I., Gallina, I. &
Richter, S. N. (2023) *Genome-wide mapping of i-Motifs reveals their association
with transcription regulation in live human cells.* Nucleic Acids Research
51(16):8309-8321. doi:10.1093/nar/gkad626

Input is the JSONL written by ``scripts/06_gse220882_peaks.py``, which does the
peak calling and sequence extraction where the bigWigs live. See that script for
why SEACR could not be run directly and what was done instead.

Why this source exists
----------------------
The biophysical i-motif panels are deep but narrow. iM-Seeker contributes 160
constructs that are variants of a small number of parents; 5DUVMA contributes
two sequences measured under many conditions. ``scripts/07_validate_external.py``
measured the consequence: on 14 CD-validated oligos from a different lab the
``im_fold`` head separates i-motifs from controls by **-0.020** -- it does not
transfer. These peaks are tens of thousands of independent genomic loci, which
is the missing axis.

What these rows are, and are not
--------------------------------
A CUT&Tag peak is an antibody binding site in chromatin. It is not a melting
curve. Every row here is ``evidence_tier="experimental"`` (something was
measured) and ``label_class="genomic_proxy"`` (what was measured is antibody
occupancy at a locus, not the folding of an oligonucleotide in a defined
buffer). The two axes have to stay separate: pooling a peak with a UV melt would
let a genomic-occupancy label train a head that claims to predict biophysics.
The condition attached to each row is nominal intracellular chemistry, and every
one of its fields is marked imputed, because nothing about the buffer was
measured.

Sequences: motifs, not windows
------------------------------
``06_gse220882_peaks.py`` emits a 201 nt window centred on the peak summit. A
window is not comparable to a 25 nt designed oligo, and training on windows
would trade one domain shift for another. So the ingested positive is the best
canonical motif found *inside* the window -- searched on both strands, since
CUT&Tag maps a duplex locus and the C-rich strand may be either one -- and the
window is kept in ``genomic`` as context. Windows with no canonical motif are
not ingested as positives; they are counted and reported, because that count is
itself the finding (only 37% of the strongest i-motif peaks contain one).

The negative control is matched at the motif level
--------------------------------------------------
Peak windows are heavily GC-biased, so a peak-versus-random-genome classifier
learns GC content and nothing else. Measured on 40,000 random hg38 windows, a
canonical i-motif appears in 3.9% of the genome and in 37.3% of the strongest
HEK293T iMab peaks -- 9.6x. Against a dinucleotide-preserving shuffle of the
same windows it is 1.36x. The composition explains most of the enrichment; the
arrangement explains the rest, and the rest is the only part a folding model can
usefully learn.

So :func:`negatives` builds each negative by shuffling a peak window with the
Altschul-Erikson dinucleotide shuffle and taking the best canonical motif from
*that*. Positives and negatives then have matched composition and both carry a
four-tract motif, and the head is forced to discriminate on tract and loop
architecture instead of counting cytosines.

WDLPS is included and flagged, not trained on
---------------------------------------------
The WDLPS libraries are roughly a tenth the depth of the HEK293T ones and their
peaks are barely enriched: 9.5% motif-bearing for iMab, 3.6% for BG4, against
3.9% and 2.1% genomic background. Those rows are ingested under their own source
names so they remain available and visible, and they are excluded from training
by default.
"""
from __future__ import annotations

import gzip
import json
import random
from pathlib import Path

from ...conditions import CONDITION_FIELDS, DEFAULTS, Condition
from ...motifs import clean, find_g4, find_im
from ...negatives import dinucleotide_shuffle
from ..db import Record

DOI = "10.1093/nar/gkad626"
GEO = "GSE220882"
URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220882"

# Nominal intracellular chemistry. Nothing here was measured.
NOMINAL = {"k": 140.0, "na": 10.0, "li_nh4": 0.0, "mg": 0.5,
           "ph": 7.2, "temperature": 37.0}


def nominal_condition() -> Condition:
    """The nominal condition, with *every* field marked imputed.

    ``Condition.from_mapping(..., track_imputed=True)`` marks only the fields the
    mapping left out, which is the right rule everywhere else and the wrong one
    here: these numbers are a plausible guess at intracellular chemistry, not a
    reading off an instrument, and a supplied guess is exactly as unmeasured as
    an omitted one. Marking all eight keeps that visible in every row and stops a
    condition-response head from fitting the constant and calling it a salt
    effect.
    """
    return Condition(
        **{f: NOMINAL.get(f, DEFAULTS[f]) for f in CONDITION_FIELDS},
        imputed=CONDITION_FIELDS,
    )

# Motif-bearing fraction of 40,000 random N-free 201 nt hg38 windows, measured
# with these same searchers. Reported alongside the peak rates so the enrichment
# is never quoted without its background.
GENOMIC_BACKGROUND = {"iM": 0.039, "G4": 0.021}


def source_name(cell_line: str, target: str) -> str:
    return f"gse220882_{cell_line.lower()}_{target.lower()}_cutandtag"


def _best_motif(seq: str, target: str):
    """Longest canonical motif in the window, searched on both strands."""
    finder = find_im if target == "iM" else find_g4
    els = finder(clean(seq), both_strands=True)
    if not els:
        return None
    return max(els, key=lambda e: (e.length, e.score))


def _iter_peaks(path: str | Path):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if line:
                yield json.loads(line)


def _record(peak: dict, el, *, folded: int, source: str, source_id: str,
            extra_flags: list[str]) -> Record:
    return Record(
        sequence=el.sequence,
        kind=peak["target"],
        source=source,
        evidence_tier="experimental" if folded else "derived",
        condition=nominal_condition(),
        label_class="genomic_proxy" if folded else "catalog",
        folded=folded,
        method=(f"{peak['antibody']} CUT&Tag in live {peak['cell_line']} cells; "
                f"SEACR-style AUC block calling, >= {peak['n_replicates']}/3 replicates"),
        source_doi=DOI,
        source_id=source_id,
        organism="Homo sapiens",
        genomic={
            "build": peak["genome_build"], "chrom": peak["chrom"],
            "window_start": peak["window_start"], "window_end": peak["window_end"],
            "peak_start": peak["peak_start"], "peak_end": peak["peak_end"],
            "summit": peak["summit"], "motif_offset_in_window": el.start,
            "cell_line": peak["cell_line"], "antibody": peak["antibody"],
            "geo": GEO,
        },
        qc_flags=[
            "label_is_antibody_occupancy_at_a_locus_not_a_folding_measurement",
            "condition_is_nominal_intracellular_chemistry_nothing_was_measured",
            f"cell_line={peak['cell_line']}",
            f"antibody={peak['antibody']}",
            f"tier={peak['tier']}",
            f"n_replicates={peak['n_replicates']}",
            f"peak_auc={peak['auc']}",
            f"strand={el.strand}",
            f"shares_locus_with_other_target={peak['shares_locus_with_other_target']}",
        ] + extra_flags,
    )


def load(
    path: str | Path,
    *,
    tier: str = "stringent",
    cell_lines: tuple[str, ...] = ("HEK293T",),
    targets: tuple[str, ...] = ("iM", "G4"),
    require_motif: bool = True,
    exclude_shared_loci: bool = False,
    dry_run: bool = False,
):
    """Ingest peak-derived positives.

    ``cell_lines`` defaults to HEK293T alone: see the module docstring on WDLPS
    depth. Pass ``("HEK293T", "WDLPS")`` to take both.
    """
    out: list[Record] = []
    stats: dict[str, int] = {}

    def bump(k):
        stats[k] = stats.get(k, 0) + 1

    for peak in _iter_peaks(path):
        if peak["tier"] != tier or peak["cell_line"] not in cell_lines \
                or peak["target"] not in targets:
            continue
        key = f"{peak['cell_line']}|{peak['target']}"
        bump(f"{key}|windows")
        if exclude_shared_loci and peak["shares_locus_with_other_target"]:
            bump(f"{key}|dropped_shared_locus")
            continue
        el = _best_motif(peak["sequence"], peak["target"])
        if el is None:
            bump(f"{key}|no_canonical_motif")
            if require_motif:
                continue
        bump(f"{key}|kept")
        if dry_run:
            continue
        loc = f"{peak['chrom']}:{peak['peak_start']}-{peak['peak_end']}"
        if el is None:
            continue
        out.append(_record(
            peak, el, folded=1,
            source=source_name(peak["cell_line"], peak["target"]),
            source_id=f"{loc}|{peak['antibody']}",
            extra_flags=[f"peak_window_length={len(peak['sequence'])}",
                         "sequence_is_the_best_canonical_motif_inside_the_peak_window"],
        ))

    load.last_stats = stats
    return stats if dry_run else out


def negatives(
    path: str | Path,
    *,
    tier: str = "stringent",
    cell_lines: tuple[str, ...] = ("HEK293T",),
    targets: tuple[str, ...] = ("iM", "G4"),
    n_per_positive: int = 1,
    max_attempts: int = 8,
    seed: int = 0,
):
    """Composition-matched, motif-bearing negatives.

    Each peak window is dinucleotide-shuffled and the best canonical motif is
    taken from the shuffle. A shuffle that yields no motif is retried up to
    ``max_attempts`` times and then abandoned, so the yield is below 1 per
    positive; the shortfall is reported in ``negatives.last_stats``.
    """
    rng = random.Random(seed)
    out: list[Record] = []
    stats: dict[str, int] = {}

    def bump(k):
        stats[k] = stats.get(k, 0) + 1

    for peak in _iter_peaks(path):
        if peak["tier"] != tier or peak["cell_line"] not in cell_lines \
                or peak["target"] not in targets:
            continue
        if _best_motif(peak["sequence"], peak["target"]) is None:
            continue
        key = f"{peak['cell_line']}|{peak['target']}"
        made = 0
        for _ in range(max_attempts):
            if made >= n_per_positive:
                break
            sh = dinucleotide_shuffle(clean(peak["sequence"]), rng=rng)
            el = _best_motif(sh, peak["target"])
            if el is None:
                continue
            loc = f"{peak['chrom']}:{peak['peak_start']}-{peak['peak_end']}"
            out.append(_record(
                peak, el, folded=0,
                source=f"shuffled::{source_name(peak['cell_line'], peak['target'])}",
                source_id=f"shuffle{made}:{loc}",
                extra_flags=[
                    "derived_negative",
                    "Altschul-Erikson dinucleotide shuffle of the peak window, then "
                    "the best canonical motif inside the shuffle",
                    f"parent_locus={loc}",
                ],
            ))
            made += 1
            bump(f"{key}|negatives")
        if made < n_per_positive:
            bump(f"{key}|shuffle_yielded_no_motif")

    negatives.last_stats = stats
    return out


def register(atlas, *, cell_lines=("HEK293T", "WDLPS"), targets=("iM", "G4")) -> None:
    for cell in cell_lines:
        for target in targets:
            ab = "iMab" if target == "iM" else "BG4"
            name = source_name(cell, target)
            atlas.register_source(
                name,
                title=f"{ab} CUT&Tag peaks in live {cell} cells ({GEO}), "
                      f"reduced to the best canonical {target} motif per reproducible peak",
                doi=DOI,
                url=URL,
                evidence_tier="experimental",
                notes=(
                    "label_class=genomic_proxy: the measurement is antibody occupancy at "
                    "a genomic locus, not folding of an oligonucleotide in a defined "
                    "buffer. The attached condition is nominal intracellular chemistry "
                    "and is entirely imputed. Peaks called from the deposited hg38 "
                    "bigWigs by a SEACR-style AUC block caller cut at the peak counts "
                    "the authors report in Table S1, kept only where >= 2 of 3 "
                    "replicates agree, and filtered against the ENCODE hg38 blacklist. "
                    f"Genomic background for a canonical {target} motif in a 201 nt "
                    f"window: {GENOMIC_BACKGROUND[target]:.1%}."
                ),
            )
            atlas.register_source(
                f"shuffled::{name}",
                title=f"Dinucleotide-shuffled, motif-bearing negatives for {name}",
                doi=DOI,
                url=URL,
                evidence_tier="derived",
                notes="Altschul-Erikson dinucleotide shuffle of the peak window, then "
                      "the best canonical motif inside the shuffle -- so positives and "
                      "negatives match on composition and both carry a four-tract "
                      "motif, and only architecture separates them.",
            )
