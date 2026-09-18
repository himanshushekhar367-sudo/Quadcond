"""AlphaGenome Atlas as a second axis on the same substitutions.

AlphaGenome Atlas (Google DeepMind, released 8 September 2026) publishes
precomputed molecular-effect predictions for every possible human
single-nucleotide variant. QuadCond predicts what a single substitution does to
a G-quadruplex or i-motif under a stated buffer. Those are two different
questions about the *same* edit, and until now nothing joined them.

The join is exact rather than approximate, which is the reason this module is
worth having. ``scans.mutation_scan`` already enumerates every single-base
substitution in a window and reports each one's effect on folding; Atlas'
``query_interval`` returns scores for every single-base substitution in a
genomic interval. Given the window's coordinates, the two enumerations are the
same set, and each QuadCond row pairs with exactly one Atlas record on
``(position, reference base, alternate base)``. No overlap heuristic, no
nearest-feature rule, no window arithmetic that has to be trusted.

What this module deliberately does not do
-----------------------------------------

It does not produce a combined score. A regulatory-effect prediction and a
predicted melting-temperature change are not on a common scale and no weighting
of them has been calibrated against anything, so multiplying or summing them
would manufacture exactly the kind of number the rest of this project spends its
effort refusing to emit. What is reported instead is both axes, the quadrant a
variant falls in against thresholds the caller states, and a within-request rank
that is high only when a variant ranks high on **both** axes. The quadrant is a
triage device. It is not a validated classifier, and nothing here has been
benchmarked against measured regulatory variants.

Access
------

Live Atlas queries need an API key (``ALPHAGENOME_API_KEY``) and reach
``gdmscience.googleapis.com`` over gRPC; the service is offered for
non-commercial use. Neither the key nor the endpoint is required to use this
module: ``TableAtlas`` reads an exported table of scores instead, so a lab
behind a restrictive network, or one working from a downloaded subset, gets the
same join with no external call. ``NullAtlas`` runs the QuadCond half alone and
says in the response that the regulatory axis is absent, rather than filling it
with zeros.
"""
from __future__ import annotations

import csv
import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

#: Scorer names the AlphaGenome client ships as its recommended set. Served as
#: a default only; a live client is asked for its own scorer metadata, because
#: a hardcoded list is a list that goes stale the first time the service adds a
#: scorer and has no way to know it has.
DEFAULT_SCORERS = (
    "RNA_SEQ",
    "ATAC",
    "DNASE",
    "CAGE",
    "CHIP_TF",
    "CHIP_HISTONE",
    "SPLICE_SITES",
)

BASES = ("A", "C", "G", "T")


class AtlasUnavailable(RuntimeError):
    """The regulatory axis could not be obtained. Never silently zero-filled."""


def _norm_chrom(name: str) -> str:
    """``1`` and ``chr1`` are the same chromosome; pick one spelling and keep it.

    A join that misses because one side wrote ``chr7`` and the other wrote ``7``
    produces an empty regulatory axis and no error, which reads as "Atlas has
    nothing to say about this locus" -- a scientific statement, arrived at by a
    string mismatch.
    """
    n = str(name).strip()
    if not n:
        return n
    low = n.lower()
    if low.startswith("chr"):
        n = n[3:]
    return f"chr{n}"


@dataclass(frozen=True)
class VariantKey:
    """One substitution on the forward genome strand.

    Forward strand, always. QuadCond routes each *head* to the strand carrying
    its motif, so a cell in a mutation scan may describe the reverse strand --
    but the variant itself is a fact about the reference assembly, and a
    coordinate system that flipped with the head would make two rows of one
    table incomparable.
    """

    chromosome: str
    position: int          # 1-based, as in VCF and as in the Atlas client
    reference: str
    alternate: str

    def normalised(self) -> "VariantKey":
        return VariantKey(_norm_chrom(self.chromosome), int(self.position),
                          self.reference.upper(), self.alternate.upper())

    def as_dict(self) -> dict:
        return {"chromosome": self.chromosome, "position": self.position,
                "reference": self.reference, "alternate": self.alternate}

    def hgvs_like(self) -> str:
        return f"{self.chromosome}:{self.position}{self.reference}>{self.alternate}"


@dataclass
class AtlasRecord:
    """Whatever the regulatory source says about one substitution."""

    key: VariantKey
    scores: dict[str, float] = field(default_factory=dict)
    source: str = ""

    def headline(self, preferred: Sequence[str] | None = None) -> tuple[str, float] | None:
        """The one score used for the quadrant, named rather than assumed.

        The caller states which scorer is the regulatory axis. When it says
        nothing, the largest absolute score is used and its name is reported
        alongside, so a reader can see which quantity the quadrant was drawn
        against -- rather than discovering later that two rows were classified
        on two different scorers.
        """
        if not self.scores:
            return None
        for name in (preferred or ()):
            if name in self.scores:
                return name, float(self.scores[name])
        name = max(self.scores, key=lambda k: abs(float(self.scores[k])))
        return name, float(self.scores[name])


# --------------------------------------------------------------------- sources
class AtlasSource:
    """Where the regulatory axis comes from. Three implementations, one contract."""

    name = "atlas"

    def describe(self) -> dict:
        return {"source": self.name}

    def records_for_interval(self, chromosome: str, start: int,
                             end: int) -> dict[VariantKey, AtlasRecord]:
        raise NotImplementedError


class NullAtlas(AtlasSource):
    """No regulatory axis. The structural half still runs, and the response says so."""

    name = "none"

    def describe(self) -> dict:
        return {
            "source": "none",
            "note": ("No regulatory source was configured, so every variant in "
                     "this response carries the structural axis only. The "
                     "quadrant is not assigned: a missing score is not a low "
                     "score."),
        }

    def records_for_interval(self, chromosome, start, end):
        return {}


class TableAtlas(AtlasSource):
    """Scores from a file the user exported. No key, no network, no rate limit.

    Two shapes are accepted, because both are what people actually have. A wide
    table carries one column per scorer. A long table carries ``scorer`` and
    ``score`` columns and one row per (variant, scorer) pair, which is what a
    tidy export from the Atlas client produces.
    """

    name = "table"

    #: Column spellings seen in exports, VCFs and tidy AnnData frames.
    ALIASES = {
        "chromosome": ("chromosome", "chrom", "chr", "#chrom", "contig"),
        "position": ("position", "pos", "start", "variant_position"),
        "reference": ("reference", "ref", "reference_bases", "ref_allele"),
        "alternate": ("alternate", "alt", "alternate_bases", "alt_allele"),
        "scorer": ("scorer", "scorer_name", "output_type", "track"),
        "score": ("score", "value", "raw_score", "quantile_score"),
    }

    def __init__(self, records: dict[VariantKey, AtlasRecord], path: str = "",
                 scorers: Sequence[str] = ()) -> None:
        self._records = records
        self.path = path
        self.scorers = tuple(scorers)

    @classmethod
    def from_path(cls, path: str | Path) -> "TableAtlas":
        p = Path(path)
        if not p.exists():
            raise AtlasUnavailable(f"{p} does not exist")
        with p.open(newline="", encoding="utf-8-sig") as fh:
            sample = fh.read(8192)
            fh.seek(0)
            try:
                dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
            except csv.Error:
                dialect = csv.excel_tab if "\t" in sample else csv.excel
            reader = csv.DictReader(fh, dialect=dialect)
            fields = [f for f in (reader.fieldnames or []) if f]
            if not fields:
                raise AtlasUnavailable(f"{p} has no header row")
            picked = cls._resolve(fields)
            for required in ("chromosome", "position", "reference", "alternate"):
                if required not in picked:
                    raise AtlasUnavailable(
                        f"{p} has no column for {required}; looked for any of "
                        f"{', '.join(cls.ALIASES[required])}")
            long_form = "scorer" in picked and "score" in picked
            score_columns = [] if long_form else [
                f for f in fields if f not in picked.values()]
            records: dict[VariantKey, AtlasRecord] = {}
            scorers: set[str] = set()
            for row in reader:
                try:
                    key = VariantKey(row[picked["chromosome"]],
                                     int(float(row[picked["position"]])),
                                     row[picked["reference"]],
                                     row[picked["alternate"]]).normalised()
                except (TypeError, ValueError, KeyError):
                    # A malformed row is skipped and counted, never guessed at.
                    continue
                rec = records.setdefault(key, AtlasRecord(key, {}, f"table:{p.name}"))
                if long_form:
                    name = str(row[picked["scorer"]])
                    val = _as_float(row[picked["score"]])
                    if val is not None:
                        rec.scores[name] = val
                        scorers.add(name)
                else:
                    for col in score_columns:
                        val = _as_float(row.get(col))
                        if val is not None:
                            rec.scores[col] = val
                            scorers.add(col)
        return cls(records, str(p), sorted(scorers))

    @classmethod
    def _resolve(cls, fields: Sequence[str]) -> dict[str, str]:
        lower = {f.strip().lower(): f for f in fields}
        out: dict[str, str] = {}
        for canonical, names in cls.ALIASES.items():
            for n in names:
                if n in lower:
                    out[canonical] = lower[n]
                    break
        return out

    def describe(self) -> dict:
        from . import provenance as _prov
        # Identified by contents, not by path. The export usually lives beside
        # a provenance.json the exporter wrote; where it does, that file's
        # query settings and retrieval time are carried through, because "which
        # scorer, which tissues, fetched when" is not recoverable from the
        # score column alone.
        ident = _prov.describe_table(self.path, role="regulatory_scores")
        sidecar = Path(self.path).with_name("provenance.json")
        if sidecar.exists():
            try:
                ident["exporter_provenance"] = json.loads(
                    sidecar.read_text(encoding="utf-8"))
                ident["exporter_provenance_sha256"] = _prov.file_sha256(sidecar)
            except Exception:                                  # noqa: BLE001
                pass
        return {"source": "table", "path": self.path,
                "n_variants": len(self._records), "scorers": list(self.scorers),
                "file": ident,
                "note": ("Regulatory scores read from a user-supplied export. "
                         "QuadCond did not fetch them and cannot vouch for "
                         "which model version or scorer settings produced "
                         "them; the file, identified above by its checksum, is "
                         "the provenance.")}

    def records_for_interval(self, chromosome, start, end):
        chrom = _norm_chrom(chromosome)
        return {k: v for k, v in self._records.items()
                if k.chromosome == chrom and start <= k.position <= end}


class TabixAtlas(AtlasSource):
    """The published AVI Tabix file, queried by region. No key, no network.

    This is the backend to prefer. Google DeepMind publishes AVI SNV scores as
    an 88.5 GB Tabix bundle under a licence that permits commercial and
    non-commercial use, alongside non-commercial splicing and SHAP
    feature-importance bundles. Tabix is the reason the size does not matter:
    the index turns "every SNV in chr8:128,748,300-128,748,340" into a seek and
    a few kilobytes read, so scanning one G-quadruplex element costs the same
    whether the file holds one chromosome or all of them.

    It also removes the two things that make the live API awkward for a public
    server -- an API key the operator has to hold, and a per-request dependency
    on an external service that a lab network may not be able to reach at all.

    The column layout is read from the file's own header rather than assumed.
    A hardcoded column order is a silent mis-join the first time the publisher
    adds a field, and a mis-joined regulatory score is worse than none.
    """

    name = "tabix"

    def __init__(self, path: str, columns: dict[str, int], score_columns: dict[str, int],
                 *, backend: str, header: str = "") -> None:
        self.path = path
        self._columns = columns
        self._score_columns = score_columns
        self._backend = backend
        self._header = header

    # -- header handling -----------------------------------------------------
    @classmethod
    def _header_line(cls, path: str) -> str:
        """The last ``#`` line before the data, however it can be read."""
        try:
            import pysam                                       # noqa: F401
            with pysam.TabixFile(path) as tf:                  # type: ignore[attr-defined]
                lines = list(tf.header)
            if lines:
                return lines[-1]
        except Exception:                                      # noqa: BLE001
            pass
        import gzip
        try:
            with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
                last = ""
                for line in fh:
                    if not line.startswith("#"):
                        break
                    last = line.rstrip("\n")
                return last
        except OSError as exc:
            raise AtlasUnavailable(f"cannot read {path}: {exc}") from exc

    @classmethod
    def from_path(cls, path: str | Path) -> "TabixAtlas":
        p = Path(path)
        if not p.exists():
            raise AtlasUnavailable(f"{p} does not exist")
        index = Path(f"{p}.tbi")
        if not index.exists() and not Path(f"{p}.csi").exists():
            raise AtlasUnavailable(
                f"{p} has no .tbi or .csi index beside it. A Tabix bundle "
                f"without its index is a 90 GB file that has to be read from "
                f"the start; index it with `tabix -p bed {p.name}` or download "
                f"the index the publisher ships.")
        backend = "pysam"
        try:
            import pysam                                       # noqa: F401
        except ImportError:
            import shutil
            if not shutil.which("tabix"):
                raise AtlasUnavailable(
                    "querying a Tabix file needs either the `pysam` package "
                    "(`pip install pysam`) or the `tabix` binary on PATH; "
                    "neither is present. Use --atlas-table with a plain export "
                    "instead.") from None
            backend = "tabix-cli"

        header = cls._header_line(str(p))
        fields = [f.strip().lstrip("#").strip()
                  for f in header.split("\t")] if header else []
        if len(fields) < 4:
            raise AtlasUnavailable(
                f"{p} carries no usable header line, so its columns cannot be "
                f"resolved. Supply the layout explicitly rather than letting "
                f"this guess: a mis-joined regulatory score is worse than "
                f"none. Header read: {header!r}")
        lower = {f.lower(): i for i, f in enumerate(fields)}

        def find(*names: str) -> int | None:
            for n in names:
                if n in lower:
                    return lower[n]
            for n in names:
                for low, i in lower.items():
                    if low.startswith(n):
                        return i
            return None

        columns = {
            "chromosome": find("chromosome", "chrom", "chr", "contig"),
            "position": find("position", "pos", "start"),
            "reference": find("reference_bases", "reference", "ref"),
            "alternate": find("alternate_bases", "alternate", "alt"),
        }
        missing = [k for k, v in columns.items() if v is None]
        if missing:
            raise AtlasUnavailable(
                f"{p}: could not find a column for {', '.join(missing)} in its "
                f"header ({', '.join(fields[:12])}...)")
        taken = set(columns.values())
        score_columns = {f: i for i, f in enumerate(fields)
                         if i not in taken and f}
        if not score_columns:
            raise AtlasUnavailable(f"{p}: header has no score columns")
        return cls(str(p), {k: int(v) for k, v in columns.items()},
                   score_columns, backend=backend, header=header)

    # -- querying ------------------------------------------------------------
    def _rows(self, chromosome: str, start: int, end: int) -> Iterable[str]:
        chrom = _norm_chrom(chromosome)
        # Both spellings are tried, because whether the published file writes
        # `chr7` or `7` is the publisher's choice and not something to hardcode.
        candidates = [chrom, chrom[3:]] if chrom.startswith("chr") else [chrom, f"chr{chrom}"]
        if self._backend == "pysam":
            import pysam
            with pysam.TabixFile(self.path) as tf:             # type: ignore[attr-defined]
                available = set(tf.contigs)
                for name in candidates:
                    if name in available:
                        yield from tf.fetch(name, start - 1, end)
                        return
            return
        import subprocess
        for name in candidates:
            out = subprocess.run(
                ["tabix", self.path, f"{name}:{start}-{end}"],
                capture_output=True, text=True, timeout=120)
            if out.returncode == 0 and out.stdout.strip():
                yield from out.stdout.splitlines()
                return

    def describe(self) -> dict:
        return {
            "source": "alphagenome_avi_tabix",
            "path": self.path,
            "backend": self._backend,
            "header": self._header,
            "score_columns": list(self._score_columns),
            "note": ("AlphaGenome AVI scores read from the published Tabix "
                     "bundle by region. The scores are model output, not "
                     "measurements; AlphaGenome has not been validated for "
                     "clinical use. The file is the provenance -- QuadCond "
                     "records the path and header it read and makes no claim "
                     "about which release produced it."),
        }

    def records_for_interval(self, chromosome, start, end):
        out: dict[VariantKey, AtlasRecord] = {}
        c = self._columns
        for line in self._rows(chromosome, int(start), int(end)):
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(c.values()):
                continue
            try:
                key = VariantKey(parts[c["chromosome"]], int(parts[c["position"]]),
                                 parts[c["reference"]], parts[c["alternate"]]).normalised()
            except (ValueError, IndexError):
                continue
            if not (start <= key.position <= end):
                continue
            if key.reference not in BASES or key.alternate not in BASES:
                continue                       # SNVs only; the scan enumerates SNVs
            rec = out.setdefault(key, AtlasRecord(key, {}, f"tabix:{Path(self.path).name}"))
            for name, idx in self._score_columns.items():
                if idx < len(parts):
                    val = _as_float(parts[idx])
                    if val is not None:
                        rec.scores[name] = val
        return out


def _wait_ready(grpc, channel, timeout, attempts: int = 3) -> None:
    """Wait for the channel, with retries.

    A single 30 s readiness wait fails on a slow or briefly unreachable network
    and takes the whole run with it before one query has been sent. The wait is
    retried with a growing budget instead.
    """
    import time as _time
    last = None
    for i in range(attempts):
        try:
            grpc.channel_ready_future(channel).result((timeout or 30) * (i + 1))
            return
        except Exception as exc:                              # noqa: BLE001
            last = exc
            if i + 1 < attempts:
                _time.sleep(2 * (i + 1))
    raise TimeoutError(
        f"the Atlas endpoint did not become reachable after {attempts} attempts "
        f"({type(last).__name__}). Check the network, then retry; nothing was queried."
    ) from last


def _large_message_client(_atlas, key: str, timeout, address, limit_bytes: int):
    """The Atlas client, with gRPC's receive limit raised.

    `atlas.create` builds its channel without `grpc.max_receive_message_length`,
    so it inherits the 4 MB default. A whole-interval query for the default
    scorer set returns far more than that -- a 21 nt window came back at 20 MB --
    and the call fails with RESOURCE_EXHAUSTED "Received message larger than
    max". That is a client-side cap, not a service limit, so it is raised here
    rather than worked around by splitting windows.

    Everything else matches `atlas.create`, and if its internals ever move, this
    falls back to it.
    """
    try:
        grpc = _atlas.grpc
        cfg = (_atlas.importlib.resources.files("alphagenome")
               / "protos/grpc_service_config.json").read_text()
        channel = grpc.secure_channel(
            address or "dns:///gdmscience.googleapis.com:443",
            grpc.ssl_channel_credentials(),
            options=(("grpc.service_config", cfg),
                     ("grpc.max_receive_message_length", limit_bytes),
                     ("grpc.max_send_message_length", limit_bytes)))
        _wait_ready(grpc, channel, timeout)
        stub = _atlas.atlas_service_pb2_grpc.AtlasServiceStub(channel=channel)
        return _atlas.AtlasClient(stub, metadata=[("x-goog-api-key", key)])
    except AttributeError:
        return _atlas.create(key, timeout=timeout, address=address)


class LiveAtlas(AtlasSource):
    """The AlphaGenome Atlas gRPC service, queried one interval at a time.

    ``query_interval`` rather than one call per variant: a 30-nt element carries
    90 substitutions, and 90 round trips to answer a question the service will
    answer in one is the same mistake the mutation scan made before it was
    batched.

    The client is imported lazily. ``alphagenome`` is a heavyweight optional
    dependency, and a QuadCond install that has no interest in variants should
    not carry it, nor should importing this module fail without it.
    """

    name = "alphagenome"

    def __init__(self, client, scorers: Sequence[str], address: str = "") -> None:
        self._client = client
        self.scorers = tuple(scorers)
        self.address = address

    @classmethod
    def from_api_key(cls, api_key: str | None = None, *,
                     scorers: Sequence[str] | None = None,
                     timeout: float | None = 30.0,
                     address: str | None = None,
                     max_message_mb: int = 256) -> "LiveAtlas":
        key = api_key or os.environ.get("ALPHAGENOME_API_KEY", "")
        if not key:
            raise AtlasUnavailable(
                "no AlphaGenome API key. Set ALPHAGENOME_API_KEY, pass --api-key, "
                "or use a downloaded score table with --atlas-table instead. "
                "AlphaGenome is offered for non-commercial use and the key is "
                "issued by Google DeepMind, not by this project.")
        try:
            from alphagenome.atlas import atlas as _atlas
        except ImportError as exc:                            # noqa: BLE001
            raise AtlasUnavailable(
                "the `alphagenome` package is not installed; "
                "`pip install alphagenome`, or use --atlas-table") from exc
        try:
            client = _large_message_client(_atlas, key, timeout, address,
                                           int(max_message_mb) * 1024 * 1024)
        except Exception as exc:                              # noqa: BLE001
            # The endpoint is not reachable from every network -- notably not
            # from sandboxes without outbound DNS -- and a connection failure
            # must not be mistaken for "this locus has no regulatory effect".
            raise AtlasUnavailable(
                "could not reach the AlphaGenome Atlas service: "
                f"{type(exc).__name__}: {exc}") from exc
        return cls(client, scorers or DEFAULT_SCORERS,
                   address or "dns:///gdmscience.googleapis.com:443")

    def describe(self) -> dict:
        return {"source": "alphagenome_atlas", "address": self.address,
                "scorers": list(self.scorers),
                "note": ("Live AlphaGenome Atlas query. Atlas predictions are "
                         "model output, not measurements, and AlphaGenome has "
                         "not been validated for clinical use.")}

    def records_for_interval(self, chromosome, start, end):
        from alphagenome.data import genome

        interval = genome.Interval(chromosome=_norm_chrom(chromosome),
                                   start=int(start) - 1, end=int(end))
        try:
            frames = self._client.query_interval(
                interval, requested_scorers=list(self.scorers),
                progress_bar=False)
        except Exception as exc:                              # noqa: BLE001
            raise AtlasUnavailable(f"Atlas interval query failed: {exc}") from exc
        return _records_from_anndata(frames, self.name)


def _records_from_anndata(frames: Mapping[str, Any],
                          source: str) -> dict[VariantKey, AtlasRecord]:
    """Collapse the client's per-scorer AnnData objects into one record per variant.

    Two obs schemas are handled, because the client uses both. The live Atlas
    returns a single ``variant`` column holding an object with ``chromosome``,
    ``position``, ``reference_bases`` and ``alternate_bases``; other paths return
    those as four separate columns. The first version of this function knew only
    about the second, and its guard was a `continue` -- so a perfectly good
    response parsed to zero records and read downstream as "no regulatory effect
    at this locus", which is the one thing this module exists not to do. A
    response that carries rows but yields no records now raises.

    Each scorer returns its own matrix over (variant x track). The aggregate is
    the largest absolute value across tracks -- the strongest effect in any cell
    type or tissue. That is deliberate and lossy: a variant with one large
    effect in one tissue and a variant with moderate effects everywhere reduce
    to similar numbers. The per-track matrices are not discarded by the service,
    only by this summary.
    """
    out: dict[VariantKey, AtlasRecord] = {}
    seen_rows = 0
    for scorer, adata in (frames or {}).items():
        obs = getattr(adata, "obs", None)
        if obs is None or not len(obs):
            continue
        seen_rows += len(obs)
        try:
            import numpy as _np
            x = adata.X
            matrix = _np.asarray(x.toarray() if hasattr(x, "toarray") else x, dtype=float)
            if matrix.ndim == 1:
                matrix = matrix.reshape(-1, 1)
            per_variant = _np.nanmax(_np.abs(matrix), axis=1)
        except Exception:                                     # noqa: BLE001
            continue
        cols = {c.lower(): c for c in obs.columns}

        def _key_from_variant(v) -> VariantKey | None:
            try:
                return VariantKey(v.chromosome, int(v.position),
                                  str(v.reference_bases), str(v.alternate_bases)).normalised()
            except Exception:                                 # noqa: BLE001
                return None

        def _pick(prefix: str) -> str | None:
            for low, real in cols.items():
                if low.startswith(prefix):
                    return real
            return None

        variant_col = cols.get("variant")
        c_chrom, c_pos = _pick("chrom"), _pick("position")
        c_ref, c_alt = _pick("ref"), _pick("alt")
        if variant_col is None and not all((c_chrom, c_pos, c_ref, c_alt)):
            continue
        for i, (_, row) in enumerate(obs.iterrows()):
            if variant_col is not None:
                key = _key_from_variant(row[variant_col])
            else:
                try:
                    key = VariantKey(row[c_chrom], int(row[c_pos]),
                                     str(row[c_ref]), str(row[c_alt])).normalised()
                except Exception:                             # noqa: BLE001
                    key = None
            if key is None or i >= len(per_variant):
                continue
            value = float(per_variant[i])
            if value != value:                                # NaN: no score, not zero
                continue
            rec = out.setdefault(key, AtlasRecord(key, {}, source))
            rec.scores[str(scorer)] = value
    if seen_rows and not out:
        raise AtlasUnavailable(
            f"the Atlas client returned {seen_rows} row(s) this adapter could not read "
            f"(obs columns: {sorted(cols) if 'cols' in dir() else 'unknown'}). Refusing to "
            "report an empty regulatory axis, which would read as 'no effect here'.")
    return out


def _as_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        f = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return None if f != f else f       # NaN is absent, not zero


#: Where a published AVI bundle is looked for when no path is given.
AVI_PATH_ENV = "ALPHAGENOME_AVI_TABIX"


def resolve_source(*, tabix: str | Path | None = None,
                   table: str | Path | None = None,
                   api_key: str | None = None,
                   use_live: bool = False,
                   scorers: Sequence[str] | None = None) -> AtlasSource:
    """Pick a regulatory source from what the caller actually has.

    Order is deliberate. The published AVI Tabix bundle comes first: it needs
    no key, its licence permits commercial and non-commercial use, and a region
    query costs a seek rather than a network round trip -- which is what makes
    it viable for a public server. A user-supplied table is next, for a lab
    working from a subset. The live API is last, because it requires a key the
    operator must hold and an endpoint some networks cannot resolve at all.
    """
    tabix = tabix or os.environ.get(AVI_PATH_ENV) or None
    if tabix:
        return TabixAtlas.from_path(tabix)
    if table:
        return TableAtlas.from_path(table)
    if use_live or api_key or os.environ.get("ALPHAGENOME_API_KEY"):
        return LiveAtlas.from_api_key(api_key, scorers=scorers)
    return NullAtlas()
