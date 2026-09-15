"""Adapter: BG4 and iMab CUT&Tag peak overlap at merged HEK293T loci.

The two structures are not independent. They occupy complementary strands of
one duplex, and forming either requires that duplex to open in a way that
presents the other strand unpaired. ``Predictor.competition`` currently
multiplies two heads' probabilities as though the events were independent,
which is the one assumption the biology most clearly violates -- and it also
compares two heads trained on different label semantics, so the product is not
on any interpretable scale.

Fixing that needs a label derived from measurement rather than from an
assumption. Zanin et al. ran BG4 and iMab CUT&Tag on the same HEK293T cell
population, so every reproducible peak from either antibody can be placed on a
common coordinate system and each window asked which antibodies had peaks there:

**This is peak overlap, not co-occupancy.** The two antibodies were used in
PARALLEL REACTIONS on separate aliquots. A window with reproducible peaks from
both shows that both structures were detectable in that cell population at that
position; no part of the design observes both antibodies bound to the same DNA
molecule, and no part of it constrains whether the two structures could coexist
there at the same instant. Earlier versions of this file said co-occupancy was
"an observation rather than an assumption", which overstates what a parallel-
reaction design can deliver: it replaces an *independence* assumption with a
*measurement of marginal detectability*, which is a real improvement and a
different claim. The head trained on it is named ``locus_peak_overlap_state``
for that reason.

=============  ===============================================================
class          meaning
=============  ===============================================================
``both``       BG4 and iMab peaks overlap this locus
``g4_only``    BG4 only
``im_only``    iMab only
``neither``    composition-matched shuffle of a locus window
=============  ===============================================================

This is still a genomic proxy -- occupancy, not folding, with the iMab caveats
that carries -- so every row is ``label_class="genomic_proxy"`` and the head
trained from it may not be read as a thermodynamic competition. What it
replaces is worse: an assumption of independence that is known to be false.

**Strand canonicalisation.** A locus has two strands and the same site reads as
G-rich one way and C-rich the other. Every window is therefore oriented to its
G-richer strand before featurisation, so a G-tract is always a G-tract and the
model is not asked to learn that a sequence and its reverse complement are the
same locus. The orientation applied is recorded per row.

**Loci, not peaks.** A site with peaks from both antibodies appears twice in
the peak file, once from each antibody's call, with slightly different
boundaries. Overlapping
peaks from both targets are merged into non-overlapping loci first; otherwise
the same site would enter the training set two or three times and inflate the
apparent sample size.
"""
from __future__ import annotations

import gzip
import json
import random
from collections import defaultdict
from pathlib import Path

from ...conditions import CONDITION_FIELDS, DEFAULTS, Condition
from ...motifs import clean, revcomp
from ...negatives import dinucleotide_shuffle
from ..db import Record

SOURCE = "gse220882_hek293t_locus_state"
NEG_SOURCE = "shuffled::gse220882_hek293t_locus_state"
DOI = "10.1093/nar/gkad626"
GEO = "GSE220882"
URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220882"

CLASSES = ("both", "g4_only", "im_only", "neither")

NOMINAL = {"k": 140.0, "na": 10.0, "li_nh4": 0.0, "mg": 0.5,
           "ph": 7.2, "temperature": 37.0}


def nominal_condition() -> Condition:
    """Nominal intracellular chemistry with every field marked imputed."""
    return Condition(
        **{f: NOMINAL.get(f, DEFAULTS[f]) for f in CONDITION_FIELDS},
        imputed=CONDITION_FIELDS,
    )


def _iter_peaks(path: str | Path, cell_line: str, tier: str):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            d = json.loads(line)
            if d["cell_line"] == cell_line and d["tier"] == tier:
                yield d


def g_rich_orientation(seq: str) -> tuple[str, str]:
    """Return the G-richer strand and which orientation that was.

    Ties go to the given strand, so the choice is deterministic.
    """
    s = clean(seq)
    rc = revcomp(s)
    return (s, "+") if s.count("G") >= rc.count("G") else (rc, "-")


def build_loci(path: str | Path, *, cell_line: str = "HEK293T",
               tier: str = "stringent") -> list[dict]:
    """Merge overlapping peaks from both antibodies into labelled loci."""
    by_chrom: dict[str, list[dict]] = defaultdict(list)
    for d in _iter_peaks(path, cell_line, tier):
        by_chrom[d["chrom"]].append(d)

    loci: list[dict] = []
    for chrom, peaks in by_chrom.items():
        peaks.sort(key=lambda d: (d["peak_start"], d["peak_end"]))
        cur: dict | None = None
        for p in peaks:
            if cur is not None and p["peak_start"] < cur["end"]:
                cur["end"] = max(cur["end"], p["peak_end"])
                cur["targets"].add(p["target"])
                cur["members"].append(p)
            else:
                if cur is not None:
                    loci.append(cur)
                cur = {"chrom": chrom, "start": p["peak_start"], "end": p["peak_end"],
                       "targets": {p["target"]}, "members": [p]}
        if cur is not None:
            loci.append(cur)

    for lo in loci:
        t = lo.pop("targets")
        lo["state"] = ("both" if t == {"G4", "iM"}
                       else "g4_only" if t == {"G4"} else "im_only")
        # represent the locus by its strongest contributing peak's window
        best = max(lo["members"], key=lambda p: p["auc"])
        lo["window"] = best["sequence"]
        lo["summit"] = best["summit"]
        lo["auc"] = best["auc"]
        lo["n_peaks"] = len(lo["members"])
        del lo["members"]
    return loci


def _record(seq: str, state: str, lo: dict, *, source: str, source_id: str,
            orientation: str, extra: list[str]) -> Record:
    return Record(
        sequence=seq,
        kind="locus",
        source=source,
        evidence_tier="experimental" if state != "neither" else "derived",
        condition=nominal_condition(),
        label_class="genomic_proxy" if state != "neither" else "catalog",
        topology=state,                      # the four-class label lives here
        folded=None,
        method="BG4 and iMab CUT&Tag in the same live HEK293T cells; "
               "reproducible peaks merged into loci",
        source_doi=DOI,
        source_id=source_id,
        organism="Homo sapiens",
        genomic={"build": "hg38", "chrom": lo["chrom"], "start": lo["start"],
                 "end": lo["end"], "summit": lo["summit"], "geo": GEO,
                 "strand_used": orientation},
        qc_flags=[
            # This flag travels into the atlas and outlives every document, so
            # it has to be right on its own. It read
            # "label_is_antibody_co_occupancy_not_thermodynamic_competition",
            # which denied the wrong thing: it ruled out competition while
            # asserting co-occupancy, and the design does not support that
            # either. Parallel reactions on separate aliquots give peak overlap
            # in a population.
            "label_is_cutandtag_peak_overlap_in_a_population",
            "not_co_occupancy_antibodies_were_parallel_reactions_on_separate_aliquots",
            "not_thermodynamic_competition_no_free_energy_was_measured",
            "training_window_nt=201",
            "condition_is_nominal_intracellular_chemistry_nothing_was_measured",
            f"locus_peak_overlap_state={state}",
            f"strand_oriented_G_rich={orientation}",
            f"peaks_merged={lo['n_peaks']}",
        ] + extra,
    )


def load(path: str | Path, *, cell_line: str = "HEK293T", tier: str = "stringent",
         dry_run: bool = False):
    loci = build_loci(path, cell_line=cell_line, tier=tier)
    stats: dict[str, int] = {"loci": len(loci)}
    for lo in loci:
        stats[lo["state"]] = stats.get(lo["state"], 0) + 1
    if dry_run:
        return stats

    out: list[Record] = []
    for lo in loci:
        seq, orient = g_rich_orientation(lo["window"])
        out.append(_record(
            seq, lo["state"], lo, source=SOURCE,
            source_id=f"{lo['chrom']}:{lo['start']}-{lo['end']}",
            orientation=orient,
            extra=[f"peak_auc={lo['auc']}"],
        ))
    load.last_stats = stats
    return out


def negatives(path: str | Path, *, cell_line: str = "HEK293T",
              tier: str = "stringent", seed: int = 0, per_class_cap: int | None = None):
    """The ``neither`` class: composition-matched shuffles of the same windows.

    Shuffling the locus windows rather than sampling random genome keeps GC and
    dinucleotide composition matched, so the head cannot separate the classes on
    base composition alone -- which, at these loci, is most of the apparent
    signal.
    """
    rng = random.Random(seed)
    loci = build_loci(path, cell_line=cell_line, tier=tier)
    if per_class_cap:
        rng.shuffle(loci)
        loci = loci[:per_class_cap]
    out: list[Record] = []
    for i, lo in enumerate(loci):
        sh = dinucleotide_shuffle(clean(lo["window"]), rng=rng)
        seq, orient = g_rich_orientation(sh)
        out.append(_record(
            seq, "neither", lo, source=NEG_SOURCE,
            source_id=f"shuffle:{lo['chrom']}:{lo['start']}-{lo['end']}",
            orientation=orient,
            extra=["derived_negative",
                   "Altschul-Erikson dinucleotide shuffle of the locus window"],
        ))
    return out


def register(atlas) -> None:
    atlas.register_source(
        SOURCE,
        title="BG4 / iMab CUT&Tag peak overlap at merged HEK293T loci "
              "(GSE220882): both / g4_only / im_only",
        doi=DOI, url=URL, evidence_tier="experimental",
        notes="PEAK OVERLAP, NOT CO-OCCUPANCY. BG4 and iMab CUT&Tag were run as "
              "parallel reactions on separate aliquots of one HEK293T population, "
              "so a 'both' window means both antibodies gave reproducible peaks "
              "over it in that population -- nothing observes two structures on the "
              "same molecule, and nothing establishes they could coexist there. "
              "Still a genomic proxy: the label is antibody occupancy at a locus, "
              "not folding, and the attached condition is nominal with every field "
              "imputed. Windows are 201 nt and oriented to the G-rich strand; "
              "overlapping peaks from both antibodies are merged into one locus "
              "before labelling.",
    )
    atlas.register_source(
        NEG_SOURCE,
        title="Composition-matched 'neither' class for the locus-state head",
        doi=DOI, url=URL, evidence_tier="derived",
        notes="Dinucleotide shuffles of the same locus windows, so the four classes "
              "match on composition and the head must separate them on arrangement.",
    )
