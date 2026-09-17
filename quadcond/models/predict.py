"""The prediction API: probabilities, intervals, applicability, and evidence.

``Predictor.predict`` never returns a bare number.  Every call returns, for
each head that applies:

* the calibrated quantity (probability / posterior / value with a conformal
  interval),
* the *claim* that head is entitled to make, in words,
* an applicability verdict for the requested condition, because a model trained
  only in 100 mM K+ has no business answering confidently at 10 mM,
* and -- when an atlas is attached -- the nearest real measurements, so the
  user can check the prediction against evidence rather than trusting it.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from .. import claims
from ..schema import PREDICTION_SCHEMA_VERSION
from ..conditions import Condition
from ..features import featurize
from ..motifs import Element, clean, find_g4, find_im, g4hunter_mean, revcomp
from ..thermo import (DEFAULT_DH_G4, DEFAULT_HILL_IM, folded_fraction_ph,
                      folded_fraction_thermal)
from .base import MultiTaskModel


#: Binary folding heads whose output on long genomic windows is refused.
GENOMIC_WINDOW_REFUSED = frozenset({"g4_fold", "g4_fold_genomic", "im_fold", "im_fold_genomic"})


def _motif_gate(head, seq: str) -> list[str]:
    """Warn when a head is being asked about a sequence class it never saw.

    A G4 head is trained on motif-containing positives and shuffles of them, so
    handing it a C-rich strand produces a number with no meaning. The i-motif
    heads are the mirror image. The number is still returned -- silently
    suppressing outputs hides bugs -- but it is flagged.

    The dispatch used to be ``if head.kind == "G4": ... else: <i-motif gate>``,
    which was correct while there were two kinds. ``locus`` made it wrong: a
    locus head answers "which antibodies had peaks over this window", a question
    that is asked of G-rich, C-rich and neither-rich windows alike, and its
    training set contains all three by construction. Falling into the else
    branch gave it a C-tract gate, so a perfectly in-domain 201-nt G-rich locus
    was reported out of domain for lacking an i-motif. The gate now dispatches
    on the kinds it knows and stays silent on the ones it does not, because a
    motif requirement invented for a head that has none is worse than no gate.
    """
    issues: list[str] = []
    if head.kind == "G4":
        if not find_g4(seq) and g4hunter_mean(seq) < 0.5:
            issues.append(
                "no canonical G4 motif and G4Hunter mean < 0.5: this head was trained "
                "only on G4-forming sequences and their shuffles, so this score is "
                "extrapolation"
            )
    elif head.kind == "iM":
        if not find_im(seq) and g4hunter_mean(seq) > -0.5:
            issues.append(
                "no canonical four-C-tract i-motif motif: this head was trained only on "
                "motif-carrying sequences and their shuffles, so this score is extrapolation"
            )
    elif head.kind == "locus":
        # No motif requirement: the training windows span G-rich, C-rich and
        # neither. What does constrain it is length, and that is enforced by the
        # applicability table -- every training window is 201 nt.
        pass
    return issues


def _domain_check(head, condition: Condition, seq: str) -> dict:
    """Is this query inside the region the head actually learned?"""
    ap = head.applicability or {}
    issues: list[str] = _motif_gate(head, seq)
    checks = {
        "k": condition.k, "na": condition.na, "li_nh4": condition.li_nh4,
        "mg": condition.mg, "ph": condition.ph, "temperature": condition.temperature,
        # Present in every Condition and in every applicability table, and until
        # now checked by neither: a query at 40% PEG or 500 uM strand went
        # through as in-domain against training data that never left 0% and 5 uM.
        "crowder_pct": condition.crowder_pct, "strand_conc": condition.strand_conc,
    }
    for field, value in checks.items():
        rng = ap.get(field)
        if not rng:
            continue
        if field == "temperature" and head.target == "tm":
            # A melting-temperature head predicts the temperature; it does not
            # take one. Every G4STAB row records temperature=25 because that is
            # a metadata default on a Tm measurement, not a condition that
            # varied, so the "training range" is [25, 25] and any physiological
            # query was being flagged out-of-domain for asking a question the
            # head does not answer with this field. The requested temperature is
            # the evaluation point of the downstream two-state model, and it is
            # range-checked there instead.
            continue
        if rng["n_unique"] == 1:
            # One distinct value is not a range. Reporting "outside [10, 10]"
            # invites the reader to think 10 was measured and 12 would be
            # interpolation; nothing here was measured as a function of this
            # field at all.
            if value != rng["min"]:
                issues.append(
                    f"{field}: every training record used {field}={rng['min']:g}, "
                    f"so this head has no {field} response to apply to your "
                    f"{field}={value:g}. Not an extrapolation -- an absence."
                )
            continue
        if value < rng["min"] or value > rng["max"]:
            issues.append(
                f"{field}={value:g} is outside the training range "
                f"[{rng['min']:g}, {rng['max']:g}]"
            )
        elif rng["n_unique"] <= 2 and value not in (rng["min"], rng["max"]):
            issues.append(
                f"{field} took only {rng['n_unique']} distinct value(s) in training "
                f"({rng['min']:g}); this head cannot resolve {field}"
            )
    # Was this a sequence at all?
    #
    # `clean` masks anything outside ACGT to N, so "ZZZZ" becomes "NNNN" -- a
    # non-empty string that passes every emptiness check downstream and
    # featurizes happily into a confident melting temperature. A batch of 200
    # FASTA records is exactly where a mangled one is never eyeballed, and the
    # number it produces is indistinguishable from a real one.
    #
    # An all-N query is refused rather than warned about: there is no sequence
    # here to be in or out of any domain. Partial N is a warning, because a real
    # sequence carrying a few ambiguity codes is a legitimate query whose
    # features were nonetheless computed over masked positions.
    cleaned = clean(seq)
    n_count = cleaned.count("N")
    all_n = bool(cleaned) and n_count == len(cleaned)
    if not all_n and n_count:
        issues.append(
            f"{n_count} of {len(cleaned)} positions are N (unrecognised "
            f"characters are masked, not dropped); features were computed over "
            f"the masked sequence")

    lr = ap.get("length")
    if lr and not (lr["min"] <= len(clean(seq)) <= lr["max"]):
        issues.append(
            f"length={len(clean(seq))} outside training range [{lr['min']:g}, {lr['max']:g}]"
        )

    # A head trained on a handful of constructs is licensed for those constructs
    # and nothing else, and the condition-range check does not catch that: ask
    # im_tm_condition about an unseen sequence at an in-range buffer and it
    # returned in_domain=True with a confident number. It also returned very
    # nearly the SAME number as for Tel21C, because a two-sequence panel gives a
    # model almost no way to depend on sequence -- so the value was not merely
    # unlicensed, it was uninformative, and nothing in the output said so.
    # Two different verdicts, and v0.4.2 conflated them. `in_domain: False` says
    # the query is outside the training range and the number beside it is an
    # extrapolation -- a caller who understands that may still want it. A
    # *refusal* says there is no number to have. The allowlist is the second
    # case: this head returned 44.82 degC for an unseen sequence and 44.73 degC
    # for a training sequence, so its output does not depend on the sequence at
    # all. Emitting a value with a flag left it to every downstream consumer to
    # remember not to plot it, and the last release's own prose called that
    # "refusing" when the CLI was still printing the number.
    refused: str | None = None
    # A refusal, not an extrapolation. `in_domain: False` means the query lies
    # outside the training range and the number beside it is an extrapolation --
    # a caller who understands that may still want it. There is no such reading
    # here: nothing was asked. Emitting a value with a flag would leave every
    # downstream consumer to remember not to plot a melting temperature for
    # "ZZZZ".
    if all_n:
        refused = (
            "the query contains no recognisable nucleotides — every position was "
            "masked to N. Unrecognised characters are masked rather than dropped, "
            "so this reached the model as a sequence of the right length and "
            "nothing else. There is no prediction to report."
        )
        issues.append(refused)
    allow = ap.get("sequence_allowlist")
    if refused is None and allow and clean(seq) not in set(allow):
        refused = (
            f"sequence is not one of the {len(allow)} constructs this head was "
            f"trained on. It learned how those specific sequences respond to "
            f"conditions, not how sequence affects the response, and it returns "
            f"essentially the same value whatever sequence it is given. There is "
            f"no prediction to report, not an uncertain one."
        )
        issues.append(refused)
    # Folding classifiers on genomic-length windows: refused, not extrapolated.
    # These heads were trained on oligos (and, for the genomic-proxy pair, on
    # short motifs cut from peaks). On 124-nt human G4-seq windows g4_fold scored
    # AUROC 0.29-0.69 -- at or below chance against motif-matching windows that
    # G4-seq did not observe (benchmarks/published_tools). A flagged number that
    # is worse than chance is not a number to hand out; the genomic question has
    # its own model trained on G4-seq.
    if (refused is None and head.name in GENOMIC_WINDOW_REFUSED
            and lr and len(cleaned) > lr["max"]):
        refused = (
            f"a {len(cleaned)}-nt window is longer than any sequence this "
            f"folding head was trained on (max {lr['max']:g} nt), and on "
            f"genomic windows it performs at or below chance. Use "
            f"`quadcond genome-scan` (trained on G4-seq, K+) to find G4 windows, "
            f"then score the motifs it reports with g4_tm / g4_topology."
        )
        issues.append(refused)
    return {
        "in_domain": not issues,
        "refused": refused is not None,
        "refusal_reason": refused,
        "warnings": issues,
        # Kept under the applicability block for callers written against
        # v0.4.0, but derived rather than read from the retired field, so a
        # genomic-proxy head reports False here as it always should have.
        "biophysically_grounded": claims.target_semantics(
            head.name, head.training_meta) == claims.BIOPHYSICAL,
    }


@dataclass
class Predictor:
    model: MultiTaskModel
    atlas: object | None = None

    # ------------------------------------------------------------------ ctors
    @classmethod
    def load(cls, model_path: str | Path, atlas_path: str | Path | None = None) -> "Predictor":
        model = MultiTaskModel.load(model_path)
        atlas = None
        if atlas_path and Path(atlas_path).exists():
            from ..atlas import Atlas

            atlas = Atlas(atlas_path)
        pred = cls(model=model, atlas=atlas)
        # The file that was actually opened, identified by its contents.
        # `MultiTaskModel.artifact_sha256` is only populated by the build script
        # that saved it, so an artifact loaded from disk usually has none -- and
        # every run record then reported `model_artifact_sha256: ""`, naming a
        # model version and nothing that identifies the bytes behind it.
        from ..provenance import file_sha256
        pred.model_path = str(model_path)
        try:
            pred.model_sha256 = file_sha256(model_path)
        except OSError:
            pred.model_sha256 = None
        pred.atlas_path = str(atlas_path) if atlas_path else None
        return pred

    # ---------------------------------------------------------------- predict
    def predict(
        self,
        sequences: str | Sequence[str],
        condition: Condition | None = None,
        *,
        heads: Sequence[str] | None = None,
        n_neighbours: int = 3,
        dh_kcal: float = DEFAULT_DH_G4,
        hill: float = DEFAULT_HILL_IM,
    ) -> list[dict]:
        if isinstance(sequences, str):
            sequences = [sequences]
        sequences = [clean(s) for s in sequences]
        cond = condition or Condition()
        wanted = set(heads) if heads else set(self.model.heads)

        results: list[dict] = [
            {
                # The document shape is a frozen contract (quadcond/schema/
                # prediction-v1.json). A record that leaves this process into a
                # figure or someone else's pipeline arrives without the
                # conversation that produced it, so it carries its own version,
                # its model identity and its claim basis.
                "schema_version": PREDICTION_SCHEMA_VERSION,
                "model_version": getattr(self.model, "version", None),
                "atlas_fingerprint": getattr(self.model, "dataset_fingerprint", None) or None,
                "sequence": s,
                "length": len(s),
                "condition": cond.label(),
                "condition_detail": cond.to_dict(),
                # Anything the caller left unset fell back to a reference default.
                # Surfacing it here stops a defaulted 100 mM K+ from being read
                # back as a measured one.
                "condition_imputed_fields": list(cond.imputed),
                "motifs": {
                    "g4_canonical": [e.as_dict() for e in find_g4(s)],
                    "im_canonical": [e.as_dict() for e in find_im(s)],
                    "g4hunter_mean": round(g4hunter_mean(s), 4),
                },
                "predictions": {},
                "evidence": [],
            }
            for s in sequences
        ]

        for name, head in self.model.heads.items():
            if name not in wanted:
                continue
            X = featurize(sequences, [cond] * len(sequences),
                          kind=head.kind, use_conditions=head.use_conditions)
            out = head.predict(X)
            for i, s in enumerate(sequences):
                sem = claims.semantics(name, head.training_meta, head.task)
                entry: dict = {
                    "head": name,
                    "kind": head.kind,
                    "task": head.task,
                    "claim": head.training_meta.get("claim", ""),
                    # The claim basis travels with every prediction, because a
                    # number lifted out of this dict into a figure or a table
                    # otherwise arrives with no way to tell a melting curve from
                    # an antibody peak. See quadcond/claims.py.
                    "has_experimental_observation": sem["has_experimental_observation"],
                    "target_semantics": sem["target_semantics"],
                    "biophysically_grounded": sem["biophysically_grounded"],
                    "calibration_scope": sem["calibration_scope"],
                    "trained_on": head.training_meta.get("tiers", []),
                    "label_classes": head.training_meta.get("label_classes") or [],
                    "n_training_rows": head.training_meta.get("n_rows"),
                    "applicability": _domain_check(head, cond, s),
                }
                if entry["applicability"].get("refused"):
                    # Refusal is decided **before** any task-specific
                    # serialisation, and for every head type. This branch used
                    # to sit third, after `binary` and `multiclass`, so it was
                    # unreachable for a classification head: querying `ZZZZ`
                    # (which cleans to `NNNN`) refused `g4_tm` and `im_pht`
                    # correctly while `g4_fold` and `im_fold` returned
                    # calibrated probabilities with the refusal buried in
                    # `applicability` and no top-level flag -- so every consumer
                    # keying on `refused`, the batch table included, rendered
                    # classifier numbers for an input the domain checker had
                    # already declared not a sequence.
                    #
                    # A refused entry carries no `probability`, `posterior`,
                    # `value`, `interval` or `folded_fraction_at_condition` key
                    # at all. A consumer that forgets to check gets a KeyError,
                    # which is the correct failure: the alternative is a number
                    # in a figure with the flag left behind in the JSON. See
                    # quadcond/readout.py for the shared accessors.
                    entry["refused"] = True
                    entry["refusal_reason"] = entry["applicability"]["refusal_reason"]
                elif head.task == "binary":
                    entry["probability"] = round(float(out["probability"][i]), 4)
                    entry["uncalibrated_probability"] = round(float(out["raw_probability"][i]), 4)
                elif head.task == "multiclass":
                    entry["posterior"] = {
                        c: round(float(out["proba"][i][j]), 4)
                        for j, c in enumerate(head.classes)
                    }
                    entry["argmax"] = out["argmax"][i]
                    entry["confidence"] = round(float(out["confidence"][i]), 4)
                else:
                    v = float(out["value"][i])
                    entry["value"] = round(v, 3)
                    entry["ensemble_std"] = round(float(out["ensemble_std"][i]), 3)
                    if "interval_low" in out:
                        entry["interval"] = [round(float(out["interval_low"][i]), 3),
                                             round(float(out["interval_high"][i]), 3)]
                        entry["interval_kind"] = out["interval_kind"]
                        entry["interval_nominal_level"] = out["interval_nominal_level"]
                        entry["interval_note"] = out["interval_note"]
                    # turn the midpoint into the quantity people actually use
                    if head.target == "tm":
                        entry["folded_fraction_at_condition"] = round(
                            folded_fraction_thermal(cond.temperature, v, dh_kcal), 4)
                        entry["transition_model"] = f"two-state van 't Hoff, dH={dh_kcal} kcal/mol"
                    elif head.target == "ph_t":
                        entry["folded_fraction_at_condition"] = round(
                            folded_fraction_ph(cond.ph, v, hill), 4)
                        entry["transition_model"] = f"Hill titration, n={hill}"
                results[i]["predictions"][name] = entry

        if self.atlas is not None and n_neighbours:
            for i, s in enumerate(sequences):
                results[i]["evidence"] = self.atlas.neighbours(
                    s, cond, n=n_neighbours,
                )
        return results

    # ------------------------------------------------------------ competition
    def competition(
        self,
        duplex_sequence: str,
        condition: Condition | None = None,
    ) -> dict:
        """Both strands of one duplex position, scored separately. **No ranking.**

        The joint table and the signed preference this method used to return are
        withdrawn, for the same reason ``/ensemble`` was withdrawn in v0.4.3 and
        with the same consequence: there is no number here that combines the two
        structures.

        What was wrong with them. ``g4_only``, ``im_only``, ``both`` and
        ``neither`` were products of two probabilities under an independence
        assumption, applied to two events that are mutually exclusive by
        construction -- forming either structure requires the duplex to open on
        the opposite strand. ``preference_g4_minus_im`` subtracted two
        probabilities that are not on a common scale: one head discriminates G4s
        from dinucleotide-preserving shuffles, the other was, whenever it was
        present, ``im_architecture`` -- a different label on a different task --
        and the difference of two numbers with different meanings does not
        acquire a meaning from being computed. The four joint entries summed to
        one, which is exactly what makes such a table read as a probability
        distribution over states.

        A warning printed beside a number does not change what the number
        claims, and this one was reachable from the documented CLI subcommand
        and the Streamlit tab, not merely left in the library as dead code. So
        the combination is gone rather than annotated. What remains is what the
        heads can actually support: each strand's predictions, on the strand
        that carries its motif, under one buffer, with their own applicability
        and claim basis. A caller who needs a ranking has to state the
        comparison they mean and validate it.
        """
        cond = condition or Condition()
        fwd = clean(duplex_sequence)
        rev = revcomp(fwd)

        g4_seq, im_seq = (fwd, rev) if fwd.count("G") >= fwd.count("C") else (rev, fwd)
        g4 = self.predict(g4_seq, cond,
                          heads=[h for h in self.model.heads if h.startswith("g4")],
                          n_neighbours=0)[0]
        im = self.predict(im_seq, cond,
                          heads=[h for h in self.model.heads if h.startswith("im")],
                          n_neighbours=0)[0]

        return {
            "status": "STRAND_RESOLVED_EVIDENCE",
            "combination_withdrawn": {
                "removed": ["joint", "preference_g4_minus_im", "interpretation"],
                "since": "0.4.8",
                "reason": (
                    "The joint table multiplied two probabilities as if forming "
                    "a G-quadruplex and forming an i-motif were independent "
                    "events at one locus, which is exactly what they are not: "
                    "the two structures occupy complementary strands and are "
                    "mutually exclusive once the duplex opens. The signed "
                    "preference subtracted two probabilities that are not on a "
                    "common scale -- different label semantics, different tasks, "
                    "different calibration sets. Neither established a ranking "
                    "or a physical comparison, and a caveat beside a number does "
                    "not change what the number claims."),
            },
            "note": (
                "Two independent, strand-resolved sets of predictions under one "
                "buffer. They are not on a common probability scale, they do not "
                "sum to anything, and nothing here ranks one structure against "
                "the other. The competing state at physiological conditions is "
                "the duplex, which this tool does not model -- it has no "
                "strand-concentration or stoichiometry term and no duplex free "
                "energy."),
            "condition": cond.label(),
            "condition_detail": cond.to_dict(),
            "g_rich_strand": g4_seq,
            "c_rich_strand": im_seq,
            "g4": g4["predictions"],
            "im": im["predictions"],
        }

    # ------------------------------------------------------------------- scan
    def scan(
        self,
        sequence: str,
        condition: Condition | None = None,
        *,
        kinds: Sequence[str] = ("G4", "iM"),
        both_strands: bool = True,
        max_elements: int = 500,
    ) -> list[dict]:
        """Find elements in a long sequence and score each one."""
        cond = condition or Condition()
        els: list[Element] = []
        if "G4" in kinds:
            els += find_g4(sequence, both_strands=both_strands)
        if "iM" in kinds:
            els += find_im(sequence, both_strands=both_strands)
        els = sorted(els, key=lambda e: (e.start, e.kind))[:max_elements]
        out: list[dict] = []
        for e in els:
            heads = [h for h in self.model.heads
                     if h.startswith("g4" if e.kind == "G4" else "im")]
            p = self.predict(e.sequence, cond, heads=heads, n_neighbours=0)[0]
            row = e.as_dict()
            row["predictions"] = p["predictions"]
            out.append(row)
        return out
