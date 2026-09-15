"""Locating the trained models and the atlas, and refusing to load the wrong ones.

The code is 250 KB and the artifacts are 300 MB. Shipping them together makes a
PyPI release absurd and a `pip install` slow for everyone who only wants the
motif finders; shipping them apart means a fresh install has to find them, and
"find them" is where a scientific tool quietly starts producing numbers from
whatever file happened to be on disk.

So resolution is explicit and ordered, and every candidate is checksummed against
:data:`MANIFEST` before it is opened:

1. ``--model`` / ``--db``, or the ``QUADCOND_MODEL`` / ``QUADCOND_DB``
   environment variables -- an explicit choice, used as given.
2. ``./artifacts/`` and ``./data/`` -- a working copy of the repository.
3. The user cache (``$QUADCOND_CACHE``, else the platform cache directory),
   where ``quadcond assets fetch`` puts downloads.

A checksum mismatch is an error, never a warning. The reason is not tamper
paranoia: it is that ``artifacts/quadcond_model.joblib`` is a name several
versions of this project have used, and a v0.4.0 artifact loaded by v0.5 code
would answer every question with numbers whose provenance nobody could later
reconstruct. Better to stop.

The verification also happens **before** the file is unpickled, which is not
merely tidy. A truncated joblib archive does not fail politely -- it can hand the
allocator a garbage length prefix and take the process out with an OOM kill,
which looks like a memory bug rather than a corrupt download.
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import sys
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

MANIFEST_PATH = Path(__file__).resolve().parent / "assets_manifest.json"

CACHE_ENV = "QUADCOND_CACHE"
MODEL_ENV = "QUADCOND_MODEL"
DB_ENV = "QUADCOND_DB"


@dataclass(frozen=True)
class Asset:
    name: str
    filename: str
    sha256: str
    bytes: int
    kind: str            # "model" | "atlas"
    description: str
    url: str | None = None
    required: bool = True
    # SQLite rewrites page metadata whenever a database is opened, so an atlas
    # file's checksum changes after a mere read. Checking it on every load would
    # make a perfectly reproducible atlas fail verification the second time it is
    # used. Atlas identity is the ordered hash of its rows instead -- the same
    # fingerprint scripts/08_results_report.py quotes -- and the file hash is
    # used only to validate a fresh download.
    checksum_stable: bool = True
    content_fingerprint: str | None = None

    @property
    def size_mb(self) -> float:
        return self.bytes / 1e6


def _load_manifest() -> dict:
    if not MANIFEST_PATH.exists():
        return {"release": "unknown", "assets": []}
    return json.loads(MANIFEST_PATH.read_text())


def manifest() -> dict:
    return _load_manifest()


def assets() -> dict[str, Asset]:
    m = _load_manifest()
    out = {}
    for a in m.get("assets", []):
        out[a["name"]] = Asset(**a)
    return out


def cache_dir() -> Path:
    """Where fetched assets live. Honours ``QUADCOND_CACHE`` first."""
    env = os.environ.get(CACHE_ENV)
    if env:
        return Path(env).expanduser()
    if sys.platform == "darwin":
        base = Path.home() / "Library" / "Caches"
    elif os.name == "nt":
        base = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
    else:
        base = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
    return base / "quadcond"


def sha256(path: str | Path, *, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


class AssetError(RuntimeError):
    """Raised when an asset is missing, or present with the wrong contents."""


def search_paths(asset: Asset) -> list[Path]:
    """Every place this asset could legitimately be, in resolution order."""
    sub = "artifacts" if asset.kind == "model" else "data"
    return [
        Path.cwd() / sub / asset.filename,
        cache_dir() / asset.filename,
    ]


def atlas_fingerprint(path: str | Path) -> str:
    """The stable identity of an atlas: SHA-256 over its ordered rows."""
    import sqlite3

    con = sqlite3.connect(f"file:{Path(path)}?mode=ro", uri=True)
    try:
        rows = con.execute(
            "SELECT sequence, seq_hash, kind, k, na, li_nh4, mg, ph, temperature, "
            "folded, topology, tm, dg, ph_t, evidence_tier, label_class, source, "
            "source_id, condition_imputed FROM records "
            "ORDER BY source, source_id, seq_hash"
        ).fetchall()
    finally:
        con.close()
    h = hashlib.sha256()
    for r in rows:
        h.update(repr(r).encode())
    return h.hexdigest()


def verify(path: str | Path, asset: Asset) -> None:
    """Raise unless ``path`` holds exactly the asset the manifest records.

    Models are identified by their file bytes. Atlases are identified by their
    contents, for the SQLite reason recorded on :class:`Asset`.
    """
    if not asset.checksum_stable:
        if not asset.content_fingerprint:
            return
        got = atlas_fingerprint(path)
        if got != asset.content_fingerprint:
            raise AssetError(
                f"{path} is not the {asset.name} this build expects.\n"
                f"  expected content fingerprint {asset.content_fingerprint}\n"
                f"  found    content fingerprint {got}\n"
                "Refusing to use it. Rebuild with scripts/01_build_atlas.py, or "
                "pass the atlas you mean with --db."
            )
        return
    got = sha256(path)
    if got != asset.sha256:
        raise AssetError(
            f"{path} is not the {asset.name} this build expects.\n"
            f"  expected sha256 {asset.sha256}\n"
            f"  found    sha256 {got}\n"
            "Refusing to load it. A file with the right name and the wrong contents "
            "is how a release ends up reporting another version's numbers. Run "
            "'quadcond assets fetch' for the matching artifact, or pass the one you "
            "mean explicitly."
        )


def resolve(kind: str, explicit: str | None = None, *, check: bool = True) -> Path:
    """Find one asset, verify it, and return its path.

    ``explicit`` (a ``--model`` / ``--db`` value or an environment variable) is
    taken as given and still checksummed, so pointing at a stale artifact fails
    loudly rather than silently.
    """
    known = [a for a in assets().values() if a.kind == kind]
    if not known:
        raise AssetError(f"no {kind} asset is recorded in {MANIFEST_PATH}")
    asset = known[0]

    return resolve_asset(asset, explicit, check=check)


def configured_paths(asset: Asset, explicit: str | None = None) -> list[Path]:
    primary = next((a for a in assets().values() if a.kind == asset.kind), None)
    env = os.environ.get(MODEL_ENV if asset.kind == "model" else DB_ENV) if primary == asset else None
    candidate = explicit or env
    return [Path(candidate).expanduser()] if candidate else search_paths(asset)


def resolve_asset(asset: Asset, explicit: str | None = None, *, check: bool = True) -> Path:
    """Resolve exactly the paths also inspected by readiness."""
    paths = configured_paths(asset, explicit)
    for p in paths:
        if p.exists():
            if check:
                verify(p, asset)
            return p
    raise AssetError(f"could not find {asset.name}; looked in: " + ", ".join(map(str, paths)))


def status() -> list[dict]:
    """One row per asset: where it is, whether it verifies."""
    rows = []
    for a in assets().values():
        found, state, where = None, "missing", None
        for p in configured_paths(a):
            if p.exists():
                found = p
                break
        if found is not None:
            where = str(found)
            try:
                verify(found, a)
                state = "ok"
            except (AssetError, OSError, ValueError) as exc:
                state = "verification failed"
        rows.append({
            "name": a.name, "kind": a.kind, "filename": a.filename,
            "size_mb": round(a.size_mb, 1), "required": a.required,
            "state": state, "path": where, "sha256": a.sha256,
        })
    return rows


def fetch(name: str | None = None, *, dest: Path | None = None,
          force: bool = False, progress=print) -> list[Path]:
    """Download the recorded assets into the cache and verify each one.

    Downloads to a temporary name and only moves it into place once the checksum
    matches, so an interrupted transfer can never leave a half-file sitting where
    a later run would find it and trust it.
    """
    dest = Path(dest or cache_dir())
    dest.mkdir(parents=True, exist_ok=True)
    out: list[Path] = []

    for a in assets().values():
        if name and a.name != name:
            continue
        target = dest / a.filename
        if target.exists() and not force:
            try:
                verify(target, a)
                progress(f"  {a.name}: already present and verified ({target})")
                out.append(target)
                continue
            except AssetError:
                progress(f"  {a.name}: present but does not verify; re-fetching")
        if not a.url:
            progress(
                f"  {a.name}: no download URL in the manifest. This build was cut "
                "before the assets were published; copy the file into "
                f"{dest} or pass its path explicitly."
            )
            continue
        tmp = target.with_suffix(target.suffix + ".part")
        progress(f"  {a.name}: fetching {a.size_mb:.0f} MB from {a.url}")
        try:
            with urllib.request.urlopen(a.url) as r, open(tmp, "wb") as fh:
                shutil.copyfileobj(r, fh)
        except (urllib.error.URLError, OSError) as exc:
            tmp.unlink(missing_ok=True)
            raise AssetError(f"{a.name}: download failed: {exc}") from exc
        # A freshly downloaded file has not been opened by SQLite yet, so its
        # byte checksum is meaningful here even for an atlas -- this is the one
        # point where the transfer itself can be validated.
        #
        # But only against a recorded checksum. Atlas entries carry
        # `sha256: null` by design (their identity is `content_fingerprint`, an
        # ordered hash of their rows, because SQLite rewrites page metadata on
        # open), and comparing a real digest with None made every atlas download
        # fail as a mismatch. Where there is no byte checksum to compare, the
        # transfer is validated against the recorded content fingerprint
        # instead, which is the identity that build actually pinned.
        got = sha256(tmp)
        if a.sha256:
            if got != a.sha256:
                tmp.unlink(missing_ok=True)
                raise AssetError(
                    f"{a.name}: download does not match the manifest\n"
                    f"  expected {a.sha256}\n  got      {got}"
                )
        elif a.content_fingerprint:
            try:
                seen = atlas_fingerprint(tmp)
            except Exception as exc:                       # noqa: BLE001
                tmp.unlink(missing_ok=True)
                raise AssetError(
                    f"{a.name}: downloaded file could not be fingerprinted: {exc}"
                ) from exc
            if seen != a.content_fingerprint:
                tmp.unlink(missing_ok=True)
                raise AssetError(
                    f"{a.name}: download does not match the manifest\n"
                    f"  expected content fingerprint {a.content_fingerprint}\n"
                    f"  found    content fingerprint {seen}"
                )
        else:
            tmp.unlink(missing_ok=True)
            raise AssetError(
                f"{a.name}: the manifest records neither a sha256 nor a content "
                f"fingerprint, so a download cannot be verified. Refusing to "
                f"install an unverifiable asset."
            )
        tmp.replace(target)
        progress(f"  {a.name}: verified -> {target}")
        out.append(target)
    return out


# --------------------------------------------------------------------------- #
# Environment drift
# --------------------------------------------------------------------------- #
# A pickled estimator is not a specification. scikit-learn will happily unpickle
# a model built by an earlier version and then behave differently, and it says so
# through InconsistentVersionWarning -- one line inside a warnings stack trace,
# which is exactly the kind of thing that gets scrolled past. The versions the
# artifacts were built with are recorded in the manifest so the mismatch can be
# reported as a sentence instead.
#
# This is a warning rather than an error on purpose: a newer scikit-learn is the
# normal state of a fresh install and usually harmless. Reproducing the published
# metrics to the last decimal is what requires matching versions.
CRITICAL_PACKAGES = ("scikit-learn", "xgboost")


def environment_drift() -> list[dict]:
    """Packages whose installed version differs from the artifact build."""
    import importlib.metadata as md

    built = _load_manifest().get("built_with") or {}
    out = []
    for pkg, want in built.items():
        try:
            got = md.version(pkg)
        except md.PackageNotFoundError:
            got = None
        if got != want:
            out.append({"package": pkg, "built_with": want, "installed": got,
                        "critical": pkg in CRITICAL_PACKAGES})
    return out


def warn_on_drift(printer=None) -> list[dict]:
    """Report version drift once, in one readable line per package."""
    import sys as _sys

    drift = environment_drift()
    critical = [d for d in drift if d["critical"]]
    if not critical:
        return drift
    emit = printer or (lambda m: print(m, file=_sys.stderr))
    names = ", ".join(
        f"{d['package']} {d['installed'] or 'missing'} (artifact built with {d['built_with']})"
        for d in critical
    )
    emit(
        f"note: {names}. The artifacts load and are usable, but exact "
        "reproduction of the published metrics needs the recorded versions -- "
        "'quadcond info --json' lists them under built_with."
    )
    return drift
