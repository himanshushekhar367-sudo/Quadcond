"""What produced a number, recorded at the moment it is produced.

The run record used to read ``model_artifact_sha256: ""``. The field existed,
the exporters carried it, and it was empty -- because it was read from
``MultiTaskModel.artifact_sha256``, which a model only carries if the build
script that saved it wrote one in. The artifact loaded from disk had none, so
every export named a model *version* and nothing that identified the file.

That is a material omission at the best of times and an indefensible one while
the release artifact's identity is disputed: a path names a location in
somebody else's filesystem, not a set of bytes. This module computes the hash of
the file actually opened, and records the rest of what a second person would
need to reproduce the run -- library versions included, because a model pickled
under one scikit-learn and unpickled under another is not guaranteed to be the
same estimator, and the absence of a warning is not evidence that it is.
"""
from __future__ import annotations

import hashlib
import platform
import sys
from importlib import metadata
from pathlib import Path

#: Packages whose version can change a prediction. Recorded, not asserted.
RELEVANT_PACKAGES = ("scikit-learn", "xgboost", "joblib", "numpy", "scipy")


def file_sha256(path: str | Path, *, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def environment() -> dict:
    """Interpreter and library versions, as installed right now."""
    versions: dict[str, str | None] = {}
    for name in RELEVANT_PACKAGES:
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "packages": versions,
    }


def build_environment() -> dict:
    """The versions the manifest says the artifact was built under."""
    try:
        from . import assets
        return dict(assets.manifest().get("built_with", {}))
    except Exception:                                          # noqa: BLE001
        return {}


def environment_matches_build() -> dict:
    """Per-package comparison of the running environment against the build one.

    Reported as a per-package verdict rather than one boolean, because "the
    environment differs" is not actionable and "scikit-learn differs, xgboost
    does not" is. A mismatch on the calibrator library matters for the binary
    and multiclass heads; the regression heads used for Tm and pH_T are XGBoost
    and are not affected by it. Neither statement is a guarantee about
    numerical equivalence -- cross-version loading of pickled estimators is
    unsupported, and a silent load is not proof of a correct one.
    """
    running = environment()["packages"]
    built = build_environment()
    out = {}
    for name, want in built.items():
        got = running.get(name)
        out[name] = {"built_with": want, "running": got,
                     "same": (got is not None and str(got) == str(want))}
    return out


def model_identity(pred, *, manifest_asset: str = "model") -> dict:
    """Which model file answered, and whether it is the one the release expects.

    ``matches_manifest`` is three-valued on purpose. ``True`` and ``False`` mean
    the hash was compared; ``None`` means it could not be -- no path was
    recorded, or the manifest carries no hash for this asset -- and an unknown
    must not be serialised as a pass.
    """
    model = pred.model
    path = getattr(pred, "model_path", None)
    sha = getattr(pred, "model_sha256", None)
    if sha is None and path:
        try:
            sha = file_sha256(path)
        except OSError:
            sha = None

    expected_sha = expected_version = None
    try:
        from . import assets
        asset = assets.assets().get(manifest_asset)
        if asset is not None:
            expected_sha = asset.sha256
        expected_version = assets.manifest().get("release")
    except Exception:                                          # noqa: BLE001
        pass

    matches = None if (sha is None or not expected_sha) else (sha == expected_sha)
    return {
        "model_version_reported_by_artifact": getattr(model, "version", None),
        "model_file": str(path) if path else None,
        "model_file_sha256": sha,
        "dataset_fingerprint": getattr(model, "dataset_fingerprint", None) or None,
        "manifest_release": expected_version,
        "manifest_expected_model_sha256": expected_sha,
        "matches_manifest": matches,
        "note": (
            "A path names a location, not a set of bytes; the hash above is of "
            "the file that was actually opened. `matches_manifest: false` means "
            "this run used an artifact the release does not describe -- the "
            "numbers may still be the ones intended, but that has not been "
            "established here and metadata agreement is not estimator identity."
            if matches is False else
            "`matches_manifest: null` means the comparison could not be made, "
            "which is not a pass."
            if matches is None else
            "The loaded file matches the hash the release manifest expects."),
    }


def run_provenance(pred, *, sources: dict | None = None, **extra) -> dict:
    """One block, attached to any result that leaves this process."""
    rec = {
        "model": model_identity(pred),
        "environment": environment(),
        "environment_vs_build": environment_matches_build(),
        "caveat": (
            "Recorded, not verified. Matching library versions make a run "
            "reproducible; they do not by themselves establish that a model "
            "loaded under them computes what it computed when it was fitted."),
    }
    if sources:
        rec["inputs"] = sources
    rec.update(extra)
    return rec


def describe_table(path: str | Path, **extra) -> dict:
    """Identify an input file by its contents, with its size and modification time."""
    p = Path(path)
    try:
        st = p.stat()
        return {"path": str(p), "sha256": file_sha256(p), "bytes": st.st_size,
                "mtime_utc": __import__("datetime").datetime.utcfromtimestamp(
                    st.st_mtime).isoformat() + "Z", **extra}
    except OSError as exc:
        return {"path": str(p), "error": str(exc), **extra}
