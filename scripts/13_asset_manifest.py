#!/usr/bin/env python3
"""Regenerate quadcond/assets_manifest.json from the files on disk.

The manifest is what makes ``quadcond assets verify`` mean anything: it records
the SHA-256 the build expects, and a file with the right name and the wrong
contents is refused rather than loaded. That guarantee is only as good as the
manifest being *current*, and until now the manifest was maintained by hand --
so the v0.4.2 integrity pass rewrote both model artifacts, the checksums went
stale, and the service began refusing to start with its own model. A hand-
maintained checksum is a checksum that will eventually be wrong.

Two identity schemes, deliberately:

- **models** are identified by file bytes. They are written once and read
  many times, so the bytes are stable.
- **atlases** are identified by an ordered hash of their rows. SQLite rewrites
  page metadata when a database is opened -- even read-only, even by a query
  that changes nothing -- so the file hash of an atlas changes after a mere
  read, and a byte checksum would report tampering every time someone looked
  at it.

    python scripts/13_asset_manifest.py            # rewrite the manifest
    python scripts/13_asset_manifest.py --check    # exit 1 if it is stale
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from quadcond import __version__
from quadcond.assets import MANIFEST_PATH, atlas_fingerprint


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


SEARCH = ("artifacts", "data", ".")


def locate(filename: str, root: Path) -> Path | None:
    """Assets are not all in one place: models are built into ``artifacts/``,
    atlases into ``data/``. Look where each actually lands rather than making
    the caller pass two roots."""
    for d in (root, *(Path(x) for x in SEARCH)):
        p = d / filename
        if p.exists():
            return p
    return None


def head_count(path: Path) -> int | None:
    """Read the head count from the artifact's own sidecar summary.

    The description said "the eleven trained heads" for two releases after there
    were twelve, because it was prose in a JSON file that nothing regenerated.
    """
    side = path.with_suffix(".json")
    if not side.exists():
        return None
    try:
        heads = json.loads(side.read_text()).get("heads")
    except json.JSONDecodeError:
        return None
    return len(heads) if heads else None


def entry(spec: dict, root: Path) -> dict | None:
    path = locate(spec["filename"], root)
    if path is None:
        return None
    out = dict(spec)
    n = head_count(path) if spec["kind"] == "model" else None
    if n:
        base = "The sequence-only ablation counterpart" if "seqonly" in spec["filename"] \
            else "The trained heads (sequence + conditions)"
        suffix = " Needed only to reproduce the ablation table." \
            if "seqonly" in spec["filename"] else ""
        out["description"] = f"{base}: {n} heads.{suffix}"
    out["bytes"] = path.stat().st_size
    if spec["kind"] == "atlas":
        out["sha256"] = None
        out["checksum_stable"] = False
        out["content_fingerprint"] = atlas_fingerprint(path)
    else:
        out["sha256"] = sha256(path)
        out["checksum_stable"] = True
        out["content_fingerprint"] = None
    return out


def build(root: Path, previous: dict) -> dict:
    assets, missing = [], []
    for spec in previous["assets"]:
        e = entry(spec, root)
        if e is None:
            # Keep the recorded entry. A machine that has not downloaded the
            # optional atlas must not, by regenerating the manifest, publish one
            # that no longer knows the atlas exists -- the first run of this
            # script did exactly that and silently dropped both atlases.
            missing.append(spec["filename"])
            assets.append(dict(spec))
        else:
            assets.append(e)
    built_with = previous.get("built_with") or {}
    try:
        import sklearn

        built_with = {**built_with, "scikit-learn": sklearn.__version__}
    except ImportError:
        pass
    return {
        "release": __version__,
        "note": previous.get("note", ""),
        "assets": assets,
        "built_with": built_with,
        "_missing": missing,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--artifacts", default="artifacts")
    ap.add_argument("--check", action="store_true",
                    help="do not write; exit 1 if the manifest is out of date")
    a = ap.parse_args()

    previous = json.loads(MANIFEST_PATH.read_text())
    fresh = build(Path(a.artifacts), previous)
    missing = fresh.pop("_missing")

    old = {x["filename"]: x for x in previous["assets"]}
    changed = []
    for x in fresh["assets"]:
        o = old.get(x["filename"], {})
        key = "content_fingerprint" if x["kind"] == "atlas" else "sha256"
        if o.get(key) != x[key]:
            changed.append((x["filename"], o.get(key), x[key]))

    for f in missing:
        print(f"  [absent]  {f}  (kept in the manifest; not present in {a.artifacts}/)")
    for name, was, now in changed:
        print(f"  [changed] {name}\n      was {was}\n      now {now}")
    if not changed:
        print("  every asset present matches the manifest")

    if a.check:
        if changed:
            print("\nmanifest is stale; run without --check to rewrite it")
        sys.exit(1 if changed else 0)

    if previous.get("release") != fresh["release"]:
        print(f"  [release] {previous.get('release')} -> {fresh['release']}")
    MANIFEST_PATH.write_text(json.dumps(fresh, indent=2) + "\n")
    print(f"\nwrote {MANIFEST_PATH}")


if __name__ == "__main__":
    main()
