"""quadcond -- command line interface."""
from __future__ import annotations

import argparse
import json
import sys

from . import readout
from .conditions import PRESETS, Condition
from .motifs import iter_fasta

# Not literal paths any more: an installed quadcond is rarely run from the
# repository root, and silently defaulting to "artifacts/quadcond_model.joblib"
# turned "file not found" into the first thing a new user saw. None means
# "resolve through quadcond.assets", which searches the working copy, then the
# cache, and verifies whatever it finds before opening it.
DEFAULT_DB = None
DEFAULT_MODEL = None


def _resolve(a, kind: str):
    """Turn --model / --db into a verified path, or exit with the reason."""
    from . import assets

    explicit = getattr(a, "model" if kind == "model" else "db", None)
    try:
        return assets.resolve(kind, explicit)
    except assets.AssetError as exc:
        sys.exit(f"{exc}")


# --------------------------------------------------------------------------- #
def _condition_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("experimental condition")
    g.add_argument("--preset", choices=sorted(PRESETS),
                   help="named buffer preset (overridden by explicit flags)")
    g.add_argument("--k", type=float, help="K+ (mM)")
    g.add_argument("--na", type=float, help="Na+ (mM)")
    g.add_argument("--li-nh4", type=float, help="Li+/NH4+ (mM)")
    g.add_argument("--mg", type=float, help="Mg2+ (mM)")
    g.add_argument("--ph", type=float)
    g.add_argument("--temperature", type=float, help="degrees C")
    g.add_argument("--crowder", type=float, dest="crowder_pct", help="%% w/v crowder")
    g.add_argument("--strand-conc", type=float, help="oligo concentration (uM)")


def _condition(a: argparse.Namespace) -> Condition:
    """Build the query condition, tracking what the caller never supplied.

    ``track_imputed=True`` matters here and not only on ingestion. A user who
    types ``--k 100`` and nothing else gets pH 7.0 and 25 C from the reference
    defaults, and every prediction then depends on two numbers they did not
    choose. The fields they left unset are listed back on the result.
    """
    base = dict(PRESETS[a.preset]) if getattr(a, "preset", None) else {}
    for f in ("k", "na", "li_nh4", "mg", "ph", "temperature", "crowder_pct", "strand_conc"):
        v = getattr(a, f, None)
        if v is not None:
            base[f] = v
    return Condition.from_mapping(base, track_imputed=True)


def _emit(obj, as_json: bool, printer=None) -> None:
    if as_json or printer is None:
        print(json.dumps(obj, indent=2, default=str))
    else:
        printer(obj)


# --------------------------------------------------------------------------- #
def _load_predictor(a):
    """Load the model, and the atlas too when one is available.

    A missing atlas is not fatal -- it only removes nearest-measurement
    retrieval -- so it is resolved separately and its absence is reported once
    rather than aborting the run.
    """
    from . import assets
    from .models.predict import Predictor

    assets.warn_on_drift()
    model = _resolve(a, "model")
    try:
        db = assets.resolve("atlas", getattr(a, "db", None))
    except assets.AssetError:
        db = None
        if getattr(a, "db", None):
            sys.exit(f"--db {a.db} could not be used; run 'quadcond assets status'")
    return Predictor.load(model, db)


def cmd_predict(a) -> None:

    seqs = list(a.sequences)
    if a.fasta:
        seqs += [s for _, s in iter_fasta(a.fasta)]
    if not seqs:
        sys.exit("no sequences given (positional args or --fasta)")
    pred = _load_predictor(a)
    res = pred.predict(seqs, _condition(a), n_neighbours=a.neighbours)
    _emit(res, a.json, _print_predictions)


def _print_predictions(res) -> None:
    for r in res:
        print(f"\n{r['sequence']}  ({r['length']} nt)   @ {r['condition']}")
        if r.get("condition_imputed_fields"):
            print("  not supplied, defaulted: "
                  + ", ".join(r["condition_imputed_fields"]))
        print(f"  G4Hunter mean {r['motifs']['g4hunter_mean']:+.3f}   "
              f"canonical G4 motifs: {len(r['motifs']['g4_canonical'])}   "
              f"canonical iM motifs: {len(r['motifs']['im_canonical'])}")
        for name, p in r["predictions"].items():
            line = f"  {name:<18}"
            # Refused first, for every task. A refused entry carries no scalar
            # key of any kind, so the historic probability/posterior/value chain
            # fell through to `p["value"]` and raised.
            if readout.is_refused(p):
                print(line + f"refused -- {readout.refusal_reason(p)}")
                continue
            if "probability" in p:
                line += f"p = {p['probability']:.3f}"
            elif "posterior" in p:
                best = p["argmax"]
                line += f"{best} ({p['confidence']:.3f})  " + \
                        "  ".join(f"{k}={v:.3f}" for k, v in p["posterior"].items())
            else:
                line += f"{p['value']:.2f}"
                if "interval" in p:
                    line += (f"  [{p['interval'][0]:.2f}, {p['interval'][1]:.2f}] "
                             f"{p['interval_nominal_level']:.0%} "
                             f"{p.get('interval_kind', '')}")
                if "folded_fraction_at_condition" in p:
                    line += f"   folded fraction {p['folded_fraction_at_condition']:.3f}"
            if not p["applicability"]["in_domain"]:
                line += "   [OUT OF DOMAIN]"
            ts = p.get("target_semantics")
            if ts == "genomic_proxy":
                line += "   [GENOMIC PROXY -- not a folding probability]"
            elif not p["applicability"]["biophysically_grounded"]:
                line += f"   [{ts} -- not biophysically grounded]"
            print(line)
            for w in p["applicability"]["warnings"]:
                print(f"      ! {w}")
        if r["evidence"]:
            print("  nearest measured neighbours:")
            for e in r["evidence"]:
                lab = []
                if e["topology"]:
                    lab.append(f"topology={e['topology']}")
                if e["tm"] is not None:
                    lab.append(f"Tm={e['tm']}")
                if e["ph_t"] is not None:
                    lab.append(f"pH_T={e['ph_t']}")
                print(f"      {e['sequence'][:40]:<40} sim={e['sequence_similarity']:.3f} "
                      f"{e['condition']}  {' '.join(lab)}  [{e['evidence_tier']}] "
                      f"{e['source_id'] or ''}")


def cmd_scan(a) -> None:
    from .models.predict import Predictor

    pred = Predictor.load(_resolve(a, 'model'), None)
    cond = _condition(a)
    records = []
    if a.fasta:
        items = list(iter_fasta(a.fasta))
    else:
        items = [("query", a.sequence)]
    for name, seq in items:
        for row in pred.scan(seq, cond, both_strands=not a.forward_only):
            row["seq_name"] = name
            records.append(row)
    if a.out:
        import pandas as pd

        flat = []
        for r in records:
            d = {k: v for k, v in r.items() if k != "predictions"}
            for hn, hp in r["predictions"].items():
                if readout.is_refused(hp):
                    d[f"{hn}_refused"] = readout.refusal_reason(hp)
                    continue
                if "probability" in hp:
                    d[f"{hn}_p"] = hp["probability"]
                elif "posterior" in hp:
                    d[f"{hn}_call"] = hp["argmax"]
                    for c, v in hp["posterior"].items():
                        d[f"{hn}_{c}"] = v
                else:
                    d[f"{hn}_value"] = hp["value"]
            flat.append(d)
        pd.DataFrame(flat).to_csv(a.out, index=False)
        print(f"wrote {len(flat)} elements -> {a.out}")
    else:
        _emit(records, True)


def cmd_competition(a) -> None:
    """Both strands of one duplex position. No combined score, by design.

    The joint table and the signed G4-minus-iM preference this subcommand used
    to print were withdrawn in v0.4.8. They combined two probabilities that are
    not on a common scale under an independence assumption that the two
    structures violate by construction. The subcommand is kept -- the
    strand-resolved evidence is the useful half and scripts call it -- and the
    output now says what was removed and why, so a pipeline reading the old keys
    fails visibly rather than silently reading a different quantity.
    """
    print("Strand-resolved evidence for one duplex position. The joint table "
          "and the\nsigned preference were withdrawn in v0.4.8: nothing in this "
          "output ranks one\nstructure against the other. See docs/CLAIMS.md.",
          file=sys.stderr)
    from .models.predict import Predictor

    pred = Predictor.load(_resolve(a, 'model'), None)
    _emit(pred.competition(a.sequence, _condition(a)), True)


def cmd_variant(a) -> None:
    """One genomic window, two axes: regulatory effect and structural effect.

    The window is the unit rather than a single variant because that is what
    both halves are efficient at. A mutation scan already enumerates every
    substitution in the window in one batched pass, and a Tabix region query
    returns every AVI record in the same span for the cost of a seek. Asking
    variant-by-variant would be slower on both sides and would answer a
    narrower question -- the useful output is the ranking *within* the element.
    """
    from . import alphagenome as ag
    from . import variants as var

    sequence = a.sequence
    if a.screen:
        # Before the model, before the key. A window with no canonical motif on
        # either strand returns `no_motif` for every substitution -- correct,
        # and not worth 3n regulatory queries to discover.
        seq = a.sequence or (_window_from_fasta(a.fasta, a.chromosome, a.start, a.end)
                             if a.fasta else "")
        if not seq:
            raise SystemExit("--screen needs a SEQUENCE or --fasta with coordinates")
        _emit(var.prescreen(seq, a.chromosome, a.start), True)
        return
    try:
        source = ag.resolve_source(tabix=a.avi_tabix, table=a.atlas_table,
                                   api_key=a.api_key, use_live=a.live)
    except ag.AtlasUnavailable as exc:
        print(f"regulatory source unavailable: {exc}", file=sys.stderr)
        print("continuing with the structural axis only.", file=sys.stderr)
        source = ag.NullAtlas()

    if a.fasta:
        # The window can come from an indexed assembly, so a caller with
        # coordinates does not have to paste bases and risk an off-by-one
        # between what they pasted and what they said the coordinates were.
        sequence = _window_from_fasta(a.fasta, a.chromosome, a.start, a.end)
    if not sequence:
        raise SystemExit("give a SEQUENCE, or --fasta with --chromosome/--start/--end")

    from .models.predict import Predictor
    pred = Predictor.load(_resolve(a, 'model'), None)
    res = var.variant_scan(
        pred, sequence, a.chromosome, a.start, condition=_condition(a),
        heads=a.heads.split(",") if a.heads else None,
        atlas_source=source,
        structural_head=a.structural_head,
        regulatory_scorer=a.regulatory_scorer,
        structural_threshold=a.structural_threshold,
        regulatory_threshold=a.regulatory_threshold)
    _emit(res, True)


def _window_from_fasta(path: str, chromosome: str, start: int, end: int | None) -> str:
    """Pull one window out of an indexed FASTA, 1-based inclusive."""
    if end is None:
        raise SystemExit("--fasta needs --end as well as --start")
    try:
        import pysam
    except ImportError:
        raise SystemExit(
            "--fasta needs `pysam` (`pip install pysam`) and a .fai index "
            "beside the assembly. Paste the window sequence instead if that is "
            "easier.") from None
    with pysam.FastaFile(path) as fh:                          # type: ignore[attr-defined]
        names = set(fh.references)
        chrom = chromosome if chromosome in names else (
            chromosome[3:] if chromosome.startswith("chr") and chromosome[3:] in names
            else f"chr{chromosome}")
        if chrom not in names:
            raise SystemExit(f"{chromosome} is not in {path}")
        return fh.fetch(chrom, int(start) - 1, int(end)).upper()


def cmd_sweep(a) -> None:
    """Vary one condition variable and report the response -- the core use case."""
    import numpy as np

    from .models.predict import Predictor

    pred = Predictor.load(_resolve(a, 'model'), None)
    values = np.linspace(a.start, a.stop, a.steps)
    rows = []
    for v in values:
        cond = _condition(a).replace(**{a.vary: float(v)})
        res = pred.predict(a.sequence, cond, n_neighbours=0)[0]
        row = {a.vary: float(v)}
        for name, p in res["predictions"].items():
            if readout.is_refused(p):
                row[f"{name}:refused"] = readout.refusal_reason(p)
                continue
            if "probability" in p:
                row[name] = p["probability"]
            elif "posterior" in p:
                row.update({f"{name}:{k}": val for k, val in p["posterior"].items()})
            else:
                row[name] = p["value"]
                if "folded_fraction_at_condition" in p:
                    row[f"{name}:folded_fraction"] = p["folded_fraction_at_condition"]
        row["in_domain"] = all(
            p["applicability"].get("in_domain", False)
            for p in res["predictions"].values()
            if not readout.is_refused(p))
        rows.append(row)
    if a.out:
        import pandas as pd

        pd.DataFrame(rows).to_csv(a.out, index=False)
        print(f"wrote {len(rows)} rows -> {a.out}")
    else:
        import pandas as pd

        print(pd.DataFrame(rows).to_string(index=False))


def _print_composition(atlas) -> None:
    """The total never appears without what it is made of."""
    c = atlas.composition()
    print(f"\ntotal records: {c['total']:,}")
    for key in ("biophysical_measured", "genomic_proxy", "predicted",
                "derived_or_shuffle"):
        print(f"  {c[key]:>9,}  {key:<21} {c['labels'][key]}")
    print("\n  Only the first bucket is a measurement of the structure being "
          "predicted.\n  See docs/CLAIMS.md for which head may claim what.")


def cmd_atlas(a) -> None:
    from .atlas import Atlas
    from .atlas.ingest import ADAPTERS

    atlas = Atlas(a.db)
    if a.atlas_cmd == "summary":
        df = atlas.summary()
        print(df.to_string(index=False))
        _print_composition(atlas)
    elif a.atlas_cmd == "neighbours":
        _emit(atlas.neighbours(a.sequence, _condition(a), kind=a.kind, n=a.n), True)
    elif a.atlas_cmd == "sources":
        for r in atlas.conn.execute("SELECT * FROM sources ORDER BY evidence_tier, source"):
            print(f"[{r['evidence_tier']:<12}] {r['source']}")
            if r["title"]:
                print(f"               {r['title']}")
            if r["doi"]:
                print(f"               doi:{r['doi']}")
            if r["notes"]:
                print(f"               {r['notes']}")
    elif a.atlas_cmd == "ingest-json":
        from .atlas.ingest import declarative

        spec = declarative.load_spec(a.spec)
        if a.dry_run:
            _emit(declarative.preview(spec, a.path), True)
            return
        declarative.register(atlas, spec)
        recs = declarative.load(a.path, spec=spec)
        n = atlas.add(recs)
        print(f"ingested {n} new records from {a.path} using {a.spec} "
              f"({len(recs)} parsed, source '{spec['source']}')")
    elif a.atlas_cmd == "ingest":
        mod = ADAPTERS[a.adapter]
        kw = {}
        if a.adapter == "lab":
            if not a.source:
                sys.exit("--source is required for the lab adapter")
            kw["source"] = a.source
        if a.dry_run:
            print(json.dumps(mod.load(a.path, dry_run=True, **kw), indent=2))
            return
        if a.adapter == "lab":
            mod.register(atlas, a.source, doi=a.doi or "")
        else:
            mod.register(atlas)
        recs = mod.load(a.path, **kw)
        n = atlas.add(recs)
        print(f"ingested {n} new records from {a.path} via '{a.adapter}' "
              f"({len(recs)} parsed)")
    atlas.close()


def cmd_train(a) -> None:
    from .atlas import Atlas
    from .models.train import DEFAULT_TASKS, train_all

    tasks = DEFAULT_TASKS
    if a.tasks:
        tasks = [t for t in DEFAULT_TASKS if t.name in set(a.tasks)]
    atlas = Atlas(a.db)
    model = train_all(atlas, tasks=tasks, use_conditions=not a.no_conditions,
                      allow_predicted=not a.no_predicted,
                      n_seeds=a.seeds, n_folds=a.folds)
    print("saved ->", model.save(a.out))
    atlas.close()


def cmd_report(a) -> None:
    from .report import build_report

    path = build_report(_resolve(a, 'model'), _resolve(a, 'atlas'), a.out,
                       ablation=a.ablation)
    print("wrote", path)


def cmd_assets(a) -> None:
    from . import assets

    if a.asset_cmd == "status":
        rows = assets.status()
        print(f"cache: {assets.cache_dir()}")
        print(f"{'asset':<16}{'kind':<8}{'state':<20}{'size':>9}  path")
        for r in rows:
            print(f"{r['name']:<16}{r['kind']:<8}{r['state']:<20}"
                  f"{r['size_mb']:>7.1f} MB  {r['path'] or '-'}")
        missing = [r for r in rows if r["required"] and r["state"] != "ok"]
        if missing:
            sys.exit("\nrequired assets missing or unverified: "
                     + ", ".join(r["name"] for r in missing)
                     + "\nrun: quadcond assets fetch")
        return

    if a.asset_cmd == "path":
        for kind in ("model", "atlas"):
            try:
                print(f"{kind}: {assets.resolve(kind)}")
            except assets.AssetError as exc:
                print(f"{kind}: unresolved -- {str(exc).splitlines()[0]}")
        return

    try:
        got = assets.fetch(getattr(a, "name", None),
                           dest=getattr(a, "dest", None),
                           force=getattr(a, "force", False))
    except assets.AssetError as exc:
        sys.exit(str(exc))
    if not got:
        sys.exit("nothing fetched. If the manifest has no URLs yet, copy the files "
                 f"into {assets.cache_dir()} or pass --model/--db explicitly.")


def cmd_serve(a) -> None:
    from . import service
    service.ALLOW_SYNTHETIC_MODEL = a.allow_synthetic_model
    service.serve(a.host, a.port)


def cmd_info(a) -> None:
    from . import __version__, assets, claims
    from .models.base import MultiTaskModel

    if a.json:
        m = MultiTaskModel.load(_resolve(a, "model"))
        out = m.summary()
        out["claim_basis"] = {
            n: claims.semantics(n, h.training_meta, h.task)
            for n, h in m.heads.items()
        }
        out["assets"] = assets.status()
        out["built_with"] = assets.manifest().get("built_with", {})
        out["environment_drift"] = assets.environment_drift()
        _emit(out, True)
        return

    print(f"quadcond {__version__}")
    rows = assets.status()
    drift = assets.environment_drift()
    print("\nassets")
    for r in rows:
        req = "" if r["required"] else "  (optional)"
        print(f"  {r['name']:<14}{r['state']:<20}{r['size_mb']:>7.1f} MB{req}")
        if r["path"]:
            print(f"                 {r['path']}")
    if any(r["state"] != "ok" and r["required"] for r in rows):
        print("\n  required assets are missing or do not verify. Run:  quadcond assets fetch")
        return

    model_path = _resolve(a, "model")
    m = MultiTaskModel.load(model_path)
    # Computed here, not read off the object: the checksum of an artifact cannot
    # be stored inside that artifact -- pickling happens before the bytes exist
    # to hash. It lives in the sidecar .json and in the asset manifest, and the
    # authoritative answer is the file on disk.
    print(f"\nmodel v{m.version}  sha256 {assets.sha256(model_path)[:16]}...")
    print(f"training-row fingerprint {m.dataset_fingerprint[:16]}...")
    print(f"\n{claims.headline(m.heads)}\n")
    print(f"  {'head':<20}{'claim basis':<16}{'metric':>12}   grouped by")
    key = {"binary": "auroc", "multiclass": "balanced_accuracy", "regression": "r2"}
    for n, h in m.heads.items():
        sem = claims.semantics(n, h.training_meta, h.task)
        v = h.metrics.get(key[h.task])
        flag = "" if sem["biophysically_grounded"] else "   <- not a biophysical claim"
        print(f"  {n:<20}{sem['target_semantics']:<16}{v:>12.3f}   "
              f"{h.training_meta.get('group_by', '?')}{flag}")
    if drift:
        print("\nlibrary drift from the artifact build")
        for d in drift:
            mark = "  (affects results)" if d["critical"] else ""
            print(f"  {d['package']:<14}installed {d['installed'] or '-':<10}"
                  f"built with {d['built_with']}{mark}")
    print("\nWhat each head may and may not be read as: docs/CLAIMS.md")


# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="quadcond",
        description="Condition-aware, calibrated G-quadruplex / i-motif prediction "
                    "backed by a provenance-carrying atlas.",
    )
    p.add_argument("--model", default=DEFAULT_MODEL)
    p.add_argument("--db", default=DEFAULT_DB)
    sub = p.add_subparsers(dest="cmd", required=True)

    q = sub.add_parser("predict", help="score one or more sequences")
    q.add_argument("sequences", nargs="*")
    q.add_argument("--fasta")
    q.add_argument("--neighbours", type=int, default=3)
    q.add_argument("--json", action="store_true")
    _condition_args(q)
    q.set_defaults(func=cmd_predict)

    q = sub.add_parser("scan", help="find and score every element in a long sequence")
    q.add_argument("sequence", nargs="?", default="")
    q.add_argument("--fasta")
    q.add_argument("--out")
    q.add_argument("--forward-only", action="store_true")
    _condition_args(q)
    q.set_defaults(func=cmd_scan)

    q = sub.add_parser(
        "competition",
        help="predict both strands of one locus separately; no combined score")
    q.add_argument("sequence")
    _condition_args(q)
    q.set_defaults(func=cmd_competition)

    q = sub.add_parser(
        "variant",
        help="join a mutation scan to AlphaGenome regulatory scores for the same edits")
    q.add_argument("sequence", nargs="?", default="",
                   help="forward-strand window sequence; omit when using --fasta")
    q.add_argument("--chromosome", required=True)
    q.add_argument("--start", type=int, required=True,
                   help="1-based genomic coordinate of the window's first base")
    q.add_argument("--end", type=int, help="only needed with --fasta")
    q.add_argument("--fasta", help="indexed assembly to pull the window from")
    q.add_argument("--avi-tabix", dest="avi_tabix",
                   help="published AlphaGenome AVI Tabix bundle (no API key needed)")
    q.add_argument("--atlas-table", dest="atlas_table",
                   help="a CSV/TSV export of regulatory scores")
    q.add_argument("--api-key", dest="api_key",
                   help="AlphaGenome API key for a live Atlas query")
    q.add_argument("--screen", action="store_true",
                   help="motif check only: no model, no API call, no key. Says "
                        "whether a variant scan of this window could produce "
                        "any structural delta before a query is spent on it")
    q.add_argument("--live", action="store_true",
                   help="query the live Atlas service (needs a key)")
    q.add_argument("--heads", help="comma-separated head names")
    q.add_argument("--structural-head", dest="structural_head",
                   help="which head supplies the structural axis (default: the first)")
    q.add_argument("--regulatory-scorer", dest="regulatory_scorer",
                   help="which scorer supplies the regulatory axis")
    q.add_argument("--structural-threshold", dest="structural_threshold", type=float,
                   help="absolute |delta| above which a variant is structurally high")
    q.add_argument("--regulatory-threshold", dest="regulatory_threshold", type=float,
                   help="absolute score above which a variant is regulatorily high")
    _condition_args(q)
    q.set_defaults(func=cmd_variant)

    q = sub.add_parser("sweep", help="response of every head to one condition variable")
    q.add_argument("sequence")
    q.add_argument("--vary", default="k",
                   choices=["k", "na", "li_nh4", "mg", "ph", "temperature", "crowder_pct"])
    q.add_argument("--start", type=float, default=0.0)
    q.add_argument("--stop", type=float, default=150.0)
    q.add_argument("--steps", type=int, default=16)
    q.add_argument("--out")
    _condition_args(q)
    q.set_defaults(func=cmd_sweep)

    q = sub.add_parser("atlas", help="inspect or extend the atlas")
    asub = q.add_subparsers(dest="atlas_cmd", required=True)
    asub.add_parser("summary").set_defaults(func=cmd_atlas)
    asub.add_parser("sources").set_defaults(func=cmd_atlas)
    n = asub.add_parser("neighbours")
    n.add_argument("sequence")
    n.add_argument("--kind", choices=["G4", "iM"])
    n.add_argument("-n", type=int, default=5)
    _condition_args(n)
    n.set_defaults(func=cmd_atlas)
    i = asub.add_parser("ingest")
    i.add_argument("adapter", choices=["g4sp", "g4stab-supp", "imseeker", "lab"])
    i.add_argument("path")
    i.add_argument("--source", help="dataset key (required for the lab adapter)")
    i.add_argument("--doi", default="")
    i.add_argument("--dry-run", action="store_true",
                   help="show which columns were matched and stop")
    i.set_defaults(func=cmd_atlas)
    j = asub.add_parser("ingest-json",
                        help="ingest any table using a declarative JSON adapter "
                             "(units, constants and value maps declared in the spec)")
    j.add_argument("spec", help="path to the adapter JSON")
    j.add_argument("path", help="path to the data table (csv/tsv/xlsx)")
    j.add_argument("--dry-run", action="store_true",
                   help="show the mapping and unit conversions, then stop")
    j.set_defaults(func=cmd_atlas)

    q = sub.add_parser("train", help="(re)train heads from the atlas")
    q.add_argument("--out", default=DEFAULT_MODEL)
    q.add_argument("--tasks", nargs="*")
    q.add_argument("--seeds", type=int, default=5)
    q.add_argument("--folds", type=int, default=5)
    q.add_argument("--no-conditions", action="store_true",
                   help="ablation: drop all condition features")
    q.add_argument("--no-predicted", action="store_true",
                   help="refuse to train any head on the predicted tier")
    q.set_defaults(func=cmd_train)

    q = sub.add_parser("report", help="build the HTML dashboard")
    q.add_argument("--out", default="artifacts/quadcond_report.html")
    q.add_argument("--ablation", default="artifacts/quadcond_model_seqonly.joblib")
    q.set_defaults(func=cmd_report)

    q = sub.add_parser(
        "serve",
        help="run the HTTP service, and the browser workbench if one is built")
    q.add_argument("--host", default="127.0.0.1",
                   help="0.0.0.0 to accept connections from other machines")
    q.add_argument("--port", type=int, default=8765)
    q.add_argument(
        "--allow-synthetic-model", action="store_true",
        help="serve a model trained on generated data (software fixture). "
             "Every response is flagged. Refused without this flag.")
    q.set_defaults(func=cmd_serve)

    q = sub.add_parser("info", help="version, assets, heads and what each may claim")
    q.add_argument("--json", action="store_true")
    q.set_defaults(func=cmd_info)

    q = sub.add_parser("assets",
                       help="locate, verify or download the models and the atlas")
    asub2 = q.add_subparsers(dest="asset_cmd", required=True)
    asub2.add_parser("status", help="where each asset is and whether it verifies") \
        .set_defaults(func=cmd_assets)
    fz = asub2.add_parser("fetch", help="download the recorded assets and verify them")
    fz.add_argument("--name", help="fetch only this asset")
    fz.add_argument("--dest", help="directory to fetch into (default: the user cache)")
    fz.add_argument("--force", action="store_true", help="re-download even if present")
    fz.set_defaults(func=cmd_assets)
    asub2.add_parser("path", help="print the resolved model and atlas paths") \
        .set_defaults(func=cmd_assets)
    return p


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
