"""A dependency-free HTTP service exposing the calibrated heads.

Written on ``http.server`` rather than FastAPI on purpose. The NAR Web Server
rules require a server that anyone can use without logging in, registering or
leaving an email address, and the shortest path to that is a service with no
framework, no database, no session and no account model -- there is nothing to
log into because there is nothing that stores anything.

Three endpoints and no state:

``GET  /health``     liveness, model version, asset checksums
``GET  /info``       the heads, what each may claim, and the applicability tables
``POST /predict``    sequences x one condition -> calibrated outputs
``POST /evidence``   the same, shaped for the 3D viewer: strand-resolved cards,
                     each with the quantity that was measured for that structure
                     on that strand. Deliberately NOT an ensemble -- there are no
                     Boltzmann weights and the numbers do not sum to one; see
                     ``evidence_payload`` for why the partition function this
                     endpoint used to serve was withdrawn

Every prediction carries its claim basis, its calibration scope, its
applicability warnings and the condition fields the caller never supplied. That
is not decoration: the viewer draws a number next to a picture of a molecule,
which is the most persuasive place a number can appear, and a proxy score
rendered identically to a measured melting temperature is how a tool teaches
people something false.

    python -m quadcond.service --port 8765
"""
from __future__ import annotations

import argparse
import json
import math
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

from . import __version__, assets, capabilities, claims, readout, scans
from .conditions import CONDITION_FIELDS, PRESETS, Condition
from .motifs import clean, find_g4, find_im, revcomp
from .schema import PREDICTION_SCHEMA_VERSION, prediction_schema
from .thermo import DEFAULT_DH_G4, R, folded_fraction_ph, folded_fraction_thermal

MAX_SEQUENCES = 64
MAX_LENGTH = 5000

_predictor = None
_lock = threading.Lock()

#: Set by ``--allow-synthetic-model``. Off by default, and off is a refusal
#: rather than a warning.
ALLOW_SYNTHETIC_MODEL = False


class SyntheticModelRefused(assets.AssetError):
    """A model trained on generated data, loaded by a service that did not ask."""


def _check_provenance(model) -> dict | None:
    """Refuse to serve a synthetic model unless the operator opted in.

    A branch of this project carried a generator that wrote to the canonical
    artifact path and registered its rows as experimental. Had it ever run, the
    service would have served a fit to a made-up CSV with every provenance
    badge reading "measured", and nothing anywhere would have said otherwise.

    Two independent things now stop that. The evidence tier makes every head
    classify as SYNTHETIC, which no downstream formatting can undo. And this,
    which stops the artifact loading at all: a warning would be printed once
    into a log nobody reads, while the service went on answering.
    """
    snap = getattr(model, "atlas_snapshot", None) or {}
    provenance = snap.get("provenance")
    if provenance != "synthetic":
        return None
    banner = {
        "provenance": "synthetic",
        "note": snap.get("provenance_note", ""),
        "warning": "EVERY NUMBER FROM THIS SERVICE IS A SOFTWARE TEST FIXTURE.",
    }
    if not ALLOW_SYNTHETIC_MODEL:
        raise SyntheticModelRefused(
            "this model artifact was trained on generated data "
            "(atlas_snapshot.provenance == 'synthetic'). The service refuses it "
            "by default. Start with --allow-synthetic-model if you are "
            "deliberately running the software fixture; every response will be "
            "flagged.")
    print("\n*** SYNTHETIC MODEL LOADED — NOT FOR SCIENTIFIC USE ***")
    print(f"*** {banner['note']}\n")
    return banner


def predictor():
    global _predictor
    with _lock:
        if _predictor is None:
            from .models.predict import Predictor
            model = assets.resolve("model")
            try:
                db = assets.resolve("atlas")
            except assets.AssetError:
                db = None
            p = Predictor.load(model, db)
            # Before the predictor is cached, so a refusal cannot be worked
            # around by making a second request.
            p.synthetic_banner = _check_provenance(p.model)
            _predictor = p
    return _predictor


def _banner() -> dict:
    """The synthetic-model flag, on every response that carries a number.

    Spread into each payload rather than printed at startup: a caller who
    scripts against this service never sees the console, and a JSON file that
    ends up in a figure has to carry its own warning.
    """
    p = _predictor
    b = getattr(p, "synthetic_banner", None) if p is not None else None
    return {"synthetic_model": b} if b else {}


def dg_from_tm(tm_c: float, t_c: float, dh: float = DEFAULT_DH_G4) -> float:
    """Folding free energy at ``t_c`` implied by a melting temperature."""
    return -dh * (1.0 - (t_c + 273.15) / (tm_c + 273.15))


def dg_from_pht(ph_t: float, ph: float, hill: float = 2.5,
                t_c: float = 37.0) -> float:
    """Folding free energy implied by a transitional pH.

    The Hill titration gives the folded fraction; a two-state equilibrium turns
    that into a free energy. Same trick as the thermal case, in the other
    variable, which is what lets an i-motif and a quadruplex be weighted against
    each other on one scale at all.
    """
    theta = folded_fraction_ph(ph, ph_t, hill)
    theta = min(max(theta, 1e-6), 1 - 1e-6)
    return -R * (t_c + 273.15) * math.log(theta / (1 - theta))


def _atlas_composition(pred) -> dict | None:
    """What the atlas total is made of, when an atlas is attached."""
    atlas = getattr(pred, "atlas", None)
    if atlas is None:
        snap = getattr(pred.model, "atlas_snapshot", None) or {}
        return {"total": snap.get("n_records"),
                "note": "no atlas attached to this service; total is from the "
                        "model's training snapshot and carries no breakdown"} \
            if snap.get("n_records") else None
    try:
        comp = atlas.composition()
    except Exception:                                     # noqa: BLE001
        return None
    # Which atlas this is matters: atlas_core.db holds the measured rows and
    # their shuffles, atlas.db adds a quarter-million distilled ones. A
    # composition with no name on it invites the reader to assume the other one.
    comp["atlas"] = getattr(atlas, "path", None) and str(atlas.path)
    snap = getattr(pred.model, "atlas_snapshot", None) or {}
    if snap.get("n_records") and snap["n_records"] != comp["total"]:
        comp["note"] = (
            f"this is the atlas attached to the service ({comp['total']:,} rows); "
            f"the heads were trained on {snap['n_records']:,} rows. The training "
            f"corpus is identified by atlas_fingerprint on every prediction.")
    return comp


class UnsupportedPolymer(ValueError):
    """Raised for a query the atlas has no evidence for."""


def _reject_rna(sequences: list[str], polymer: str | None) -> None:
    """Refuse RNA rather than answering it from DNA measurements.

    Every measurement in the atlas is DNA. RNA G-quadruplexes are not a small
    perturbation of that: they are effectively always parallel, they are more
    stable than their DNA counterparts at the same ionic strength, and the 2'-OH
    changes the loop energetics that most of these features encode. A head
    trained on DNA melting curves will still return a confident number for an
    RNA sequence, and that number has nothing behind it.

    Refusing is the only honest option while the atlas has no RNA rows, and it
    is better than a warning: a warning next to a number in a 3D viewer is read
    as a caveat on a real answer, not as an absence of one.
    """
    if polymer and polymer.upper() == "RNA":
        raise UnsupportedPolymer(
            "RNA is not supported. Every measurement in this atlas is DNA, and "
            "rG4s differ in topology, stability and loop energetics -- a DNA-"
            "trained head would answer confidently and wrongly. Ingest rG4-seq "
            "and published rG4 melting data as their own source first."
        )
    for s in sequences:
        if "U" in clean(s.upper().replace("U", "U")) or "U" in s.upper():
            raise UnsupportedPolymer(
                "sequence contains U, so it is RNA. Not supported: the atlas is "
                "entirely DNA and no rG4 measurement stands behind any head here."
            )


def _condition(payload: dict) -> Condition:
    base = dict(PRESETS[payload["preset"]]) if payload.get("preset") in PRESETS else {}
    for f in CONDITION_FIELDS:
        if payload.get(f) is not None:
            base[f] = float(payload[f])
    return Condition.from_mapping(base, track_imputed=True)


def predict_payload(payload: dict) -> dict:
    seqs = payload.get("sequences") or ([payload["sequence"]] if payload.get("sequence") else [])
    if not seqs:
        raise ValueError("no sequences given")
    if len(seqs) > MAX_SEQUENCES:
        raise ValueError(f"at most {MAX_SEQUENCES} sequences per request")
    for s in seqs:
        if len(s) > MAX_LENGTH:
            raise ValueError(f"sequences must be at most {MAX_LENGTH} nt")
    _reject_rna(seqs, payload.get("polymer"))
    cond = _condition(payload)
    p = predictor()
    res = p.predict(seqs, cond, heads=payload.get("heads"),
                    n_neighbours=int(payload.get("neighbours", 3)))
    return {"quadcond_version": __version__,
            "model_version": p.model.version,
            **_banner(),
            "results": res}


def mutation_scan_payload(payload: dict) -> dict:
    """Every single substitution against a frozen wild type.

    The wild-type sequence and the buffer are read once, here, and travel back
    inside the response. A client that keeps a table of deltas cannot then pair
    it with a different baseline without the mismatch being visible in the file
    it kept.
    """
    seq = payload.get("sequence")
    if not seq:
        raise ValueError("no sequence given")
    _reject_rna([seq], payload.get("polymer"))
    positions = payload.get("positions")
    if positions is not None:
        if not isinstance(positions, list) or any(isinstance(p, bool) or not isinstance(p, int) for p in positions):
            raise ValueError("mutation positions must be a list of integer indices")
    res = scans.mutation_scan(
        predictor(), seq, _condition(payload),
        heads=payload.get("heads"), positions=positions)
    return {"quadcond_version": __version__, **_banner(), **res}


def condition_scan_payload(payload: dict) -> dict:
    """One sequence along one condition axis, for the heads that learned it."""
    seq = payload.get("sequence")
    if not seq:
        raise ValueError("no sequence given")
    _reject_rna([seq], payload.get("polymer"))
    axis = payload.get("axis")
    if not axis:
        raise ValueError("no axis given")
    values = payload.get("values")
    if values is None:
        # A span plus a count, so a caller does not have to enumerate points to
        # ask a simple question.
        lo, hi = payload.get("min"), payload.get("max")
        if lo is None or hi is None:
            raise ValueError("give either `values`, or `min` and `max`")
        n = max(2, min(int(payload.get("points", 12)), scans.MAX_CONDITION_POINTS))
        step = (float(hi) - float(lo)) / (n - 1)
        values = [float(lo) + step * i for i in range(n)]
    res = scans.condition_scan(
        predictor(), seq, axis, values, _condition(payload),
        heads=payload.get("heads"),
        n_neighbours=int(payload.get("neighbours", 3)))
    return {"quadcond_version": __version__, **_banner(), **res}


def variant_payload(payload: dict) -> dict:
    """One window, both axes. The regulatory source is server-configured.

    The client does not get to name a file path or hand over an API key: the
    operator configures one source when the service starts, and every request
    uses it. A prediction endpoint that will read any path a caller sends is a
    file-disclosure endpoint with a scientific interface, and a key posted in a
    request body is a key in somebody's browser history.
    """
    from . import variants as _variants

    raw = payload.get("sequence", "")
    seq = clean(raw)
    if not seq:
        raise ValueError("no sequence given")
    _reject_rna([raw], payload.get("polymer"))
    chromosome = payload.get("chromosome")
    start = payload.get("start")
    if not chromosome or start is None:
        raise ValueError(
            "a variant scan needs `chromosome` and `start` (the 1-based "
            "coordinate of the window's first base). Without them the "
            "structural rows cannot be joined to a regulatory record, and a "
            "guessed offset is a different variant.")
    res = _variants.variant_scan(
        predictor(), seq, str(chromosome), int(start),
        condition=_condition(payload),
        heads=payload.get("heads"),
        atlas_source=_regulatory_source(),
        structural_head=payload.get("structural_head"),
        regulatory_scorer=payload.get("regulatory_scorer"),
        structural_threshold=payload.get("structural_threshold"),
        regulatory_threshold=payload.get("regulatory_threshold"))
    return {"quadcond_version": __version__, **_banner(), **res}


_regulatory: object | None = None


def _regulatory_source():
    """The one regulatory source this process was started with, resolved once.

    Cached like the predictor, and for the same reason: a Tabix handle and a
    gRPC channel are both expensive to build, and a source that resolved
    differently between two requests would make two rows of one study
    incomparable without saying so.
    """
    global _regulatory
    from . import alphagenome as _ag
    with _lock:
        if _regulatory is None:
            try:
                _regulatory = _ag.resolve_source()
            except _ag.AtlasUnavailable as exc:
                # A misconfigured source degrades to "no regulatory axis", which
                # the response states, rather than taking the endpoint down.
                _regulatory = _ag.NullAtlas()
                _regulatory.startup_error = str(exc)           # type: ignore[attr-defined]
    return _regulatory


def batch_payload(payload: dict) -> dict:
    """Many sequences, many buffers, one run record.

    Accepts either structured records or pasted FASTA. Identifiers survive
    either way: a result that renumbers the caller's sequences has to be
    re-matched by hand, and hand-matching is where a row acquires the wrong
    label.
    """
    mode = str(payload.get("parse_mode", "auto"))
    records = payload.get("records")
    if records is None:
        text = payload.get("fasta") or payload.get("sequences")
        if isinstance(text, list):
            records = [{"index": i, "id": f"seq_{i + 1}", "sequence": s}
                       for i, s in enumerate(text)]
        elif isinstance(text, str):
            # `auto` reads the text as FASTA when any line starts with ">", and
            # as one sequence per line otherwise -- which is what the client's
            # input label has always promised. The previous parser did neither
            # reliably: it concatenated bare lines into a single sequence and
            # created no record at all for a header with an empty body.
            parsed = scans.parse_input(text, mode=mode)
            records = parsed["records"]
        else:
            raise ValueError("give `records`, `fasta` or `sequences`")
    _reject_rna([str(r.get("sequence", "")) for r in records], payload.get("polymer"))

    raw_conditions = payload.get("conditions")
    if raw_conditions:
        conds = [Condition.from_mapping(c) for c in raw_conditions]
    else:
        conds = [_condition(payload)]
    res = scans.batch_predict(
        predictor(), records, conds, heads=payload.get("heads"),
        n_neighbours=int(payload.get("neighbours", 0)))
    return {"quadcond_version": __version__, **_banner(), **res}


def _g_richness(seq: str) -> int:
    """G count, the tie-break the locus ingest used to canonicalise windows."""
    return seq.upper().count("G")


def _has_value(entry: dict) -> bool:
    """A refused head has no value to read. Guards every consumer of one."""
    return not entry.get("refused") and "value" in entry


LOCUS_OVERLAP_NOTE = (
    "Peak-overlap class for a 201-nt genomic window. BG4 and iMab CUT&Tag were "
    "run as parallel reactions on separate aliquots, so a window called 'both' "
    "had reproducible peaks from both antibodies in the same cell population -- "
    "it is NOT evidence that the two structures occupied the same DNA molecule, "
    "and it is not thermodynamic competition. The 'neither' class is a "
    "composition-matched shuffle, not an observed unoccupied locus. Balanced "
    "accuracy is 0.450 against a 0.250 floor and iM-only recall is 0.156. Read "
    "it as an exploratory genomic-resemblance score."
)


def evidence_payload(payload: dict) -> dict:
    """Per-strand structural evidence for one duplex position. **Not an ensemble.**

    This replaces `/ensemble`, and the change is not cosmetic. Two things were
    wrong with putting a G-quadruplex, an i-motif and an unfolded state in one
    partition function:

    **The strand.** Motifs were detected on both strands, and then every head was
    evaluated on the *input* strand. For a G-rich sequence like Tel22 the i-motif
    lives on the reverse complement, so `im_pht` -- a model of C-tract
    sequences -- was being asked about a G-tract one. The number it returned was
    a category error, not an inaccuracy. Each head is now called on the strand
    that actually carries its motif, and the response says which.

    **The partition function.** Even with the strands right, exp(-dG/RT) over
    {G4, iM, unfolded} is not the equilibrium of this system. A G-quadruplex on
    one strand and an i-motif on its complement cannot both fold without the
    duplex being open, so the competing state that dominates at physiological
    conditions -- the duplex itself -- was missing entirely, along with strand
    concentration and stoichiometry. Published work resolves G4/duplex/single-
    strand equilibria only with the complementary strand explicitly in the model,
    and reports context-dependent G4/iM coexistence and intermediates that a
    two-state-per-fold treatment cannot represent.

    So the numbers are reported as what they are: separate, strand-resolved
    predictions with their own intervals and applicability, which do **not** sum
    to one. A user who wants an occupancy has to supply the duplex model this
    tool does not have.
    """
    raw = payload.get("sequence", "")
    seq = clean(raw)
    if not seq:
        raise ValueError("no sequence given")
    _reject_rna([raw], payload.get("polymer"))
    cond = _condition(payload)
    p = predictor()
    comp = revcomp(seq)

    # Two predictions, one per strand. The G4 heads answer about whichever
    # strand carries G-tracts; the i-motif heads about whichever carries
    # C-tracts. On a duplex position those are usually different strands.
    fwd = p.predict(seq, cond, n_neighbours=0)[0]["predictions"]
    rev = p.predict(comp, cond, n_neighbours=0)[0]["predictions"]

    cards: list[dict] = []

    def _strand_for(motif_finder, head_key):
        """Pick the strand carrying this motif, and the prediction made on it."""
        f_hits = motif_finder(seq, both_strands=False)
        r_hits = motif_finder(comp, both_strands=False)
        if f_hits and (len(f_hits) >= len(r_hits)):
            return "+", seq, fwd, f_hits
        if r_hits:
            return "-", comp, rev, r_hits
        return None, None, None, []

    strand_g4, g4_seq, g4_preds, g4_hits = _strand_for(find_g4, "g4_tm")
    strand = strand_g4
    if g4_hits and "g4_tm" in g4_preds and _has_value(g4_preds["g4_tm"]):
        tm = g4_preds["g4_tm"]["value"]
        cards.append({
            "id": "g-quadruplex",
            "kind": "g-quadruplex",
            "strand": strand,
            "strand_sequence": g4_seq,
            "topology": g4_preds.get("g4_topology", {}).get("argmax", "unknown"),
            "predictedTm": round(tm, 1),
            "predictedTm_interval": g4_preds["g4_tm"].get("interval"),
            "foldedFraction": round(folded_fraction_thermal(cond.temperature, tm), 4),
            "foldedFraction_note": (
                "Fraction folded for the isolated strand under a two-state model "
                "at the predicted Tm. It is NOT occupancy at a genomic locus: the "
                "complementary strand is not in this model, and duplex formation "
                "competes with both structures."),
            "source": "head:g4_tm",
            "claim_basis": claims.semantics(
                "g4_tm", p.model.heads["g4_tm"].training_meta, "regression"),
            "applicability": g4_preds["g4_tm"]["applicability"],
            "motifs": len(g4_hits),
        })

    strand_im, im_seq, im_preds, im_hits = _strand_for(find_im, "im_pht")
    strand = strand_im
    if im_hits and "im_pht" in im_preds and _has_value(im_preds["im_pht"]):
        pht = im_preds["im_pht"]["value"]
        cards.append({
            "id": "i-motif",
            "kind": "i-motif",
            "strand": strand,
            "strand_sequence": im_seq,
            "topology": "intercalated C:C+ tetraplex",
            "predictedPhT": round(pht, 2),
            "predictedPhT_interval": im_preds["im_pht"].get("interval"),
            "foldedFraction": round(folded_fraction_ph(cond.ph, pht), 4),
            "foldedFraction_note": (
                "Fraction folded for the isolated strand from a Hill titration at "
                "the predicted transitional pH. Same caveat as the G4 card: not "
                "occupancy, and the duplex is not in the model."),
            "source": "head:im_pht",
            "claim_basis": claims.semantics(
                "im_pht", p.model.heads["im_pht"].training_meta, "regression"),
            "applicability": im_preds["im_pht"]["applicability"],
            "motifs": len(im_hits),
        })

    # `g4_fold` / `im_fold` used to ride inside the structural cards as
    # `probabilityOfFolding` -- a number from one head sitting under another
    # head's `source`, `claim_basis` and `applicability`. Whatever a reader took
    # from that card, they took it about the wrong model. These heads answer a
    # narrower question than the Tm and pH_T heads, on a constructed task, so
    # they get their own entries carrying their own basis and their own limits.
    discrimination = []
    for _name, _label, _src, _strand, _sseq in (
        ("g4_fold", "G4 vs. matched shuffles", g4_preds, strand_g4, g4_seq),
        ("im_fold", "i-motif vs. matched shuffles", im_preds, strand_im, im_seq),
    ):
        if not _src or _name not in _src or _name not in p.model.heads:
            continue
        _sem = claims.semantics(_name, p.model.heads[_name].training_meta, "binary")
        # Refusal now reaches classification heads too (it used to be
        # unreachable for them), so this can no longer index `probability`
        # blindly. A refused head is reported as refused, with no number.
        if readout.is_refused(_src[_name]):
            discrimination.append({
                "head": _name, "label": _label, "strand": _strand,
                "strand_sequence": _sseq, "refused": True,
                "refusal_reason": readout.refusal_reason(_src[_name]),
                "claim_basis": _sem,
                "applicability": _src[_name]["applicability"],
            })
            continue
        discrimination.append({
            "head": _name,
            "label": _label,
            "strand": _strand,
            "strand_sequence": _sseq,
            "probability": _src[_name]["probability"],
            "claim_basis": _sem,
            "applicability": _src[_name]["applicability"],
            "note": (
                "Discrimination between measured examples of this structure and "
                "dinucleotide-preserving shuffles of themselves, calibrated on that "
                "constructed, roughly balanced task. NOT P(this sequence folds) for "
                "any real candidate set, and not the quantity the melting "
                "temperature or transitional pH describes -- different head, "
                "different label, different applicability."),
        })

    # The locus head was trained on windows canonicalised to their G-richer
    # strand (`g_rich_orientation` in the ingest), so it must be asked in that
    # orientation. Evaluating it on whichever strand the caller typed asks it
    # about a representation it never saw -- and because it has no motif gate,
    # nothing downstream would have flagged the mismatch.
    locus_seq, locus_preds, locus_strand = (
        (seq, fwd, "+") if _g_richness(seq) >= _g_richness(comp) else (comp, rev, "-")
    )
    joint = None
    if "locus_peak_overlap_state" in locus_preds:
        e = locus_preds["locus_peak_overlap_state"]
        if readout.is_refused(e):
            joint = {
                "refused": True,
                "refusal_reason": readout.refusal_reason(e),
                "applicability": e["applicability"],
                "strand": locus_strand,
                "strand_sequence": locus_seq,
                "note": LOCUS_OVERLAP_NOTE,
            }
        else:
            joint = {
                "posterior": e["posterior"],
                "argmax": e["argmax"],
                "claim_basis": claims.semantics(
                    "locus_peak_overlap_state",
                    p.model.heads["locus_peak_overlap_state"].training_meta,
                    "multiclass"),
                "applicability": e["applicability"],
                "strand": locus_strand,
                "strand_sequence": locus_seq,
                "note": LOCUS_OVERLAP_NOTE,
            }

    # The genomic heads take the same strand routing as the structural cards,
    # and for the same reason. v0.4.3 fixed `im_pht` and left these two pinned
    # to the forward strand -- so on a G-rich input, `im_fold_genomic`, trained
    # on iMab peak windows oriented C-rich, was scoring the G-rich strand. The
    # structural card said 0.93 for the reverse strand while the proxy card
    # reported the forward-strand value beside it, and nothing in the response
    # said they were describing different molecules.
    def _genomic(name: str, label: str, finder) -> dict | None:
        if name not in fwd and name not in rev:
            return None
        f_hits, r_hits = finder(seq, both_strands=False), finder(comp, both_strands=False)
        if r_hits and len(r_hits) > len(f_hits):
            strand, src, strand_seq = "-", rev, comp
        else:
            strand, src, strand_seq = "+", fwd, seq
        if name not in src:
            return None
        if readout.is_refused(src[name]):
            return {
                "head": name, "label": label, "refused": True,
                "refusal_reason": readout.refusal_reason(src[name]),
                "strand": strand, "strand_sequence": strand_seq,
                "claim": claims.semantics(
                    name, p.model.heads[name].training_meta,
                    "binary")["calibration_scope"],
                "applicability": src[name]["applicability"],
            }
        return {
            "head": name,
            "label": label,
            "score": src[name]["probability"],
            "strand": strand,
            "strand_sequence": strand_seq,
            "claim": claims.semantics(
                name, p.model.heads[name].training_meta, "binary")["calibration_scope"],
            "applicability": src[name]["applicability"],
        }

    genomic = [
        g for g in (
            _genomic("g4_fold_genomic", "BG4 CUT&Tag occupancy", find_g4),
            _genomic("im_fold_genomic", "iMab CUT&Tag occupancy", find_im),
        ) if g is not None
    ]

    return {
        "quadcond_version": __version__,
        "model_version": p.model.version,
        "prediction_schema_version": PREDICTION_SCHEMA_VERSION,
        # Every other numerical endpoint spreads `_banner()` into its response;
        # this one did not, so a predictor carrying a synthetic-provenance
        # warning produced a flagged /predict response and an unflagged
        # /evidence response containing the same numbers in cards. "Every
        # response is flagged" was true of five endpoints out of six, and the
        # missing one is the one the viewer renders by default.
        **_banner(),
        "model_artifact_sha256": getattr(p.model, "artifact_sha256", None),
        "atlas_fingerprint": getattr(p.model, "dataset_fingerprint", None) or None,
        "sequence": seq,
        "complement": comp,
        "evidence": cards,
        "discrimination": discrimination,
        "evidence_note": (
            "Strand-resolved structural evidence, NOT a thermodynamic ensemble. "
            "These numbers do not sum to one and must not be renormalised: each "
            "is a separate prediction for one structure on one strand of an "
            "isolated oligonucleotide. The state that competes with both at "
            "physiological conditions is the duplex, which this tool does not "
            "model -- it has no strand-concentration or stoichiometry term and no "
            "duplex free energy. QuadCond v0.4.2 served these as Boltzmann "
            "weights over {G4, i-motif, unfolded}; that partition function was "
            "not the equilibrium of this system and has been withdrawn."),
        "genomic_evidence": genomic,
        "genomic_evidence_note": (
            "Antibody occupancy at genomic loci (BG4 / iMab CUT&Tag in live "
            "HEK293T cells). An experimental observation of a different quantity "
            "than folding: never P(folds), in cells or anywhere else."),
        "condition": cond.label(),
        "condition_detail": cond.to_dict(),
        "condition_imputed_fields": list(cond.imputed),
        "locus_peak_overlap_state": joint,
        # The ensemble ΔG calibration is deliberately absent. It existed to fix
        # the scale of a Boltzmann step this endpoint no longer performs, and
        # shipping it beside strand-resolved cards would invite a caller to
        # reconstruct exactly the partition function that was withdrawn.
    }


class Handler(BaseHTTPRequestHandler):
    server_version = f"quadcond/{__version__}"

    def log_message(self, fmt, *args):        # quieter than the default
        pass

    def _send(self, code: int, body: dict) -> None:
        data = json.dumps(body, default=str).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(data)))
        # The viewer is served from a different port in development. No
        # credentials are ever sent, because there is no session to send.
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.end_headers()
        self.wfile.write(data)

    def do_OPTIONS(self):
        self._send(204, {})

    # ------------------------------------------------------------ static web
    #
    # One process, one port, one URL. The prediction service and the viewer used
    # to be two servers behind a documented proxy configuration, which is a
    # deployment note nobody reads and an extra thing to get wrong on an
    # institutional host. When a built frontend is present the API serves it,
    # so `quadcond serve` is the whole web server.
    #
    # `/api` is reserved for the JSON endpoints and the static handler never
    # sees those paths; anything else that is not a file on disk falls back to
    # index.html, because the viewer is a single-page application and a deep
    # link must not 404.
    def _web_root(self):
        from . import web
        return web.root()

    def _serve_static(self, path: str) -> bool:
        import mimetypes
        root = self._web_root()
        if root is None:
            return False
        rel = path.lstrip("/") or "index.html"
        candidate = (root / rel).resolve()
        try:
            # Containment check first: a request for ../../etc/passwd is a path,
            # not a page.
            candidate.relative_to(root.resolve())
        except ValueError:
            self._send(403, {"error": "path outside the web root"})
            return True
        if candidate.is_dir():
            candidate = candidate / "index.html"
        if not candidate.is_file():
            candidate = root / "index.html"
            if not candidate.is_file():
                return False
        data = candidate.read_bytes()
        ctype = mimetypes.guess_type(str(candidate))[0] or "application/octet-stream"
        self.send_response(200)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        # Hashed asset filenames are immutable; the entry document is not.
        if "/assets/" in str(candidate).replace("\\", "/"):
            self.send_header("Cache-Control", "public, max-age=31536000, immutable")
        else:
            self.send_header("Cache-Control", "no-cache")
        self.end_headers()
        self.wfile.write(data)
        return True

    def do_GET(self):
        route = self.path.split("?", 1)[0].rstrip("/")
        # `/` is the application when one is built, and the health document when
        # nothing is. An API-only deployment keeps its old behaviour; a full
        # deployment does not make a visitor guess a path to reach the viewer.
        if route == "" and self._web_root() is not None:
            if self._serve_static("/index.html"):
                return
        if route not in ("", "/health", "/ready", "/info", "/schema/prediction") \
                and not route.startswith("/api") \
                and self._serve_static(self.path.split("?", 1)[0]):
            return
        if self.path.rstrip("/") in ("/health", "/ready", ""):
            # Liveness and readiness are separate questions, and this handler
            # used to answer both with "ok". A process that is running but whose
            # required model file is missing reported `status: "ok"` with the
            # missing asset listed in a field nothing downstream read -- so a
            # load balancer or an uptime check would have kept a service in
            # rotation that cannot answer a single prediction.
            #
            #   /health -- is the process alive? 200 whenever it can reply.
            #   /ready  -- can it serve predictions? 503 until every required
            #              asset is present and verifies.
            want_ready = self.path.rstrip("/") == "/ready"
            try:
                rows = assets.status()
            except Exception as exc:                     # noqa: BLE001
                self._send(503, {"status": "unavailable", "ready": False,
                                 "error": str(exc)})
                return
            missing = [r["name"] for r in rows
                       if r.get("required") and r.get("state") != "ok"]
            body = {
                "status": "ok",
                "ready": not missing,
                "missing_required_assets": missing,
                "quadcond_version": __version__,
                "assets": rows,
                "authentication": "none; this service has no accounts",
                "note": ("status reports liveness only; 'ready' reports whether "
                         "the required model and atlas are present and verify. "
                         "GET /ready returns 503 when they are not."),
            }
            self._send(503 if (want_ready and missing) else 200, body)
            return
        if self.path.rstrip("/") == "/info":
            try:
                p = predictor()
                self._send(200, {
                    "quadcond_version": __version__,
                    "model_version": p.model.version,
                    "model_artifact_sha256": getattr(p.model, "artifact_sha256", None),
                    "atlas_fingerprint": getattr(p.model, "dataset_fingerprint", None) or None,
                    **_banner(),
                    "headline": claims.headline(p.model.heads),
                    "heads": {
                        n: {**claims.semantics(n, h.training_meta, h.task),
                            "task": h.task, "kind": h.kind,
                            "target": h.target,
                            "units": {"tm": "degC", "ph_t": "pH units"}.get(h.target, ""),
                            "claim": h.training_meta.get("claim", ""),
                            "n_training_rows": h.training_meta.get("n_rows"),
                            "metrics": h.metrics,
                            "applicability": h.applicability,
                            # Which of the interface's controls this head can
                            # actually respond to, and which group it belongs in.
                            # Served rather than hardcoded in the client: a list
                            # of head names typed into a component is a list that
                            # goes stale the first time a head is added, and the
                            # component has no way to know it has.
                            **capabilities.head_capabilities(h)}
                        for n, h in p.model.heads.items()},
                    "workspaces": capabilities.workspaces(p.model.heads),
                    "condition_axes": [
                        {"id": a, "label": capabilities.AXIS_LABEL[a]}
                        for a in capabilities.CONDITION_AXES],
                    "presets": PRESETS,
                    # The corpus total is never served without its breakdown.
                    # "398,375 records" reads as 398,375 measurements of DNA
                    # folding; four thousand of them are.
                    "atlas_composition": _atlas_composition(p),
                    "prediction_schema_version": PREDICTION_SCHEMA_VERSION,
                    "prediction_schema_url": "/schema/prediction",
                })
            except Exception as exc:                     # noqa: BLE001
                self._send(503, {"error": str(exc)})
            return
        if self.path.rstrip("/") == "/schema/prediction":
            # Served, not just documented: a consumer can fetch the contract it
            # is being asked to rely on and validate against it in CI.
            self._send(200, prediction_schema())
            return
        self._send(404, {"error": "unknown endpoint",
                         "endpoints": ["/health", "/info", "/predict", "/evidence",
                                       "/scan/mutations", "/scan/conditions",
                                       "/scan/variant", "/batch",
                                       "/schema/prediction"]})

    def do_POST(self):
        route = self.path.rstrip("/")
        n = int(self.headers.get("Content-Length") or 0)
        if n > 1 << 20:
            self._send(413, {"error": "request too large"})
            return
        try:
            payload = json.loads(self.rfile.read(n) or b"{}")
        except json.JSONDecodeError as exc:
            self._send(400, {"error": f"invalid JSON: {exc}"})
            return
        try:
            if route == "/predict":
                self._send(200, predict_payload(payload))
            elif route == "/scan/mutations":
                self._send(200, mutation_scan_payload(payload))
            elif route == "/scan/conditions":
                self._send(200, condition_scan_payload(payload))
            elif route == "/scan/variant":
                self._send(200, variant_payload(payload))
            elif route == "/batch":
                self._send(200, batch_payload(payload))
            elif route in ("/evidence", "/ensemble"):
                # /ensemble is kept as an alias so an older viewer gets the
                # corrected payload rather than a 404 -- but the response no
                # longer contains an `ensemble` key, so a client that reads one
                # fails loudly instead of rendering strand-confused numbers.
                self._send(200, evidence_payload(payload))
            else:
                self._send(404, {"error": "unknown endpoint"})
        except ValueError as exc:
            self._send(400, {"error": str(exc)})
        except assets.AssetError as exc:
            self._send(503, {"error": str(exc)})
        except Exception as exc:                          # noqa: BLE001
            self._send(500, {"error": f"{type(exc).__name__}: {exc}"})


def serve(host: str = "127.0.0.1", port: int = 8765) -> None:
    srv = ThreadingHTTPServer((host, port), Handler)
    shown = "localhost" if host in ("127.0.0.1", "0.0.0.0") else host
    print(f"quadcond {__version__} serving on http://{shown}:{port}")
    from . import web as _web
    root = _web.root()
    if root is not None:
        print(f"  workbench  http://{shown}:{port}/   (from {root})")
    else:
        print("  workbench  not built -- API only.")
        print("             run: npm --prefix web install && npm --prefix web run build")
    print(f"  health     http://{shown}:{port}/health")
    print(f"  readiness  http://{shown}:{port}/ready    (503 until the assets verify)")
    print(f"  info       http://{shown}:{port}/info")
    print("  POST /predict  /scan/mutations  /scan/conditions  /scan/variant  /batch")
    print("  no accounts, no sessions, no stored requests")
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")


def main() -> None:
    ap = argparse.ArgumentParser(description="QuadCond prediction service")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument(
        "--allow-synthetic-model", action="store_true",
        help="serve a model trained on generated data (software fixture). "
             "Every response is flagged. Refused without this flag.")
    a = ap.parse_args()
    global ALLOW_SYNTHETIC_MODEL
    ALLOW_SYNTHETIC_MODEL = a.allow_synthetic_model
    serve(a.host, a.port)


if __name__ == "__main__":
    main()
