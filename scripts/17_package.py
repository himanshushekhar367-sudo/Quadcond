#!/usr/bin/env python3
"""Build the source archive.

Written down rather than typed at a shell because the archive has been wrong
twice in ways a one-off command makes easy: directories at mode 0555, so an
extracted tree could not be written into, and a silently missing data
directory that a fresh checkout's test suite reads. Both were invisible
until someone else extracted the file.

Every entry is normalised on the way in -- 0755 for directories, 0644 for
regular files (0755 if it was executable), owner root:root, a deterministic
order -- so two runs of this script on the same tree produce byte-identical
archives, and an extraction is writable without ``chmod -R``.

    python3 scripts/17_package.py [--out DIR]
"""
from __future__ import annotations

import argparse
import fnmatch
import gzip
import re
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
NAME = "quadcond"

# Anything matching these is not source.
#
# Two groups, for two different reasons. The first is build residue and local
# state. The second is data this repository may read but must not redistribute:
# `data/external/downloads/` holds publisher supplementary files fetched under
# their own terms, and the trained `.joblib` estimators are 27 MB apiece with
# six dated backups beside them. Both are reproduced by documented steps --
# `docs/GET_THE_DATA.md` for the data, `scripts/02_train.py` for the models --
# and the JSON sidecars that carry every reported metric *are* shipped, so the
# archive stays checkable without shipping either.
EXCLUDE = [
    "__pycache__", "*.pyc", "*.pyo", ".git", ".pytest_cache", ".mypy_cache",
    "node_modules",
    "*.db", "*.db-journal", "*.db-wal", "*.db-shm", "*.sqlite", "*.sqlite3",
    "*.tar.gz", "*.tgz", "*.zip",
    "dist", "build", "*.egg-info", ".venv", "venv",
    "data/external/downloads", "data/external/downloads/*",
    "*.joblib", "*.joblib.*",
    ".DS_Store", "*.swp", "*~",
]

# A source archive that grows past this is shipping something it should not.
# The check exists because "why is the tarball 263 MB" is a question that only
# gets asked once the file is already on someone else's disk.
MAX_ARCHIVE_BYTES = 8 * 1024 * 1024


def excluded(rel: str) -> bool:
    """True if the path, or any directory on the way to it, is excluded.

    Matching every segment rather than just the full path and the basename is
    the whole correctness of this function: `*/build/*` does not match
    `build/lib/x.py`, because `build` is at the start, and `*.egg-info` does not
    match `quadcond.egg-info/PKG-INFO`, because the basename is `PKG-INFO`. Both
    slipped into an archive that looked right in the listing.
    """
    parts = rel.split("/")
    for pat in EXCLUDE:
        if fnmatch.fnmatch(rel, pat):
            return True
        if "/" not in pat and any(fnmatch.fnmatch(part, pat) for part in parts):
            return True
    return False


def entries(root: Path) -> list[Path]:
    """Every shipped path, parents included, in a stable order."""
    found: set[Path] = set()
    for path in sorted(root.rglob("*")):
        rel = path.relative_to(root).as_posix()
        if excluded(rel):
            continue
        if path.is_symlink() or not (path.is_file() or path.is_dir()):
            continue
        found.add(path)
        for parent in path.relative_to(root).parents:
            if parent.as_posix() not in (".", ""):
                found.add(root / parent)
    return sorted(found, key=lambda p: p.relative_to(root).as_posix())


def normalise(info: tarfile.TarInfo) -> tarfile.TarInfo:
    # 0555 directories are the reason this function exists: an extracted tree
    # that cannot be written into is not a source release.
    info.uid = info.gid = 0
    info.uname = info.gname = "root"
    info.mtime = 0
    info.pax_headers = {}
    if info.isdir():
        info.mode = 0o755
    else:
        info.mode = 0o755 if info.mode & 0o100 else 0o644
    return info


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT.parent), help="directory to write the archive into")
    args = ap.parse_args()

    # Read the version out of the source rather than importing the package: this
    # script has to work in a bare checkout, before any install.
    text = (ROOT / NAME / "__init__.py").read_text(encoding="utf-8")
    match = re.search(r'^__version__ = "([^"]+)"', text, re.M)
    if match is None:
        raise SystemExit(f"no __version__ in {NAME}/__init__.py")
    out = Path(args.out) / f"{NAME}-v{match.group(1)}-source.tar.gz"
    out.parent.mkdir(parents=True, exist_ok=True)

    paths = entries(ROOT)
    with out.open('wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as compressed:
        with tarfile.open(fileobj=compressed, mode='w') as tar:
            for path in paths:
                rel = path.relative_to(ROOT).as_posix()
                tar.add(path, arcname=f"{NAME}/{rel}", recursive=False, filter=normalise)

    size = out.stat().st_size
    files = sum(1 for p in paths if p.is_file())
    if size > MAX_ARCHIVE_BYTES:
        biggest = sorted(
            ((p.stat().st_size, p.relative_to(ROOT).as_posix()) for p in paths if p.is_file()),
            reverse=True,
        )[:10]
        out.unlink()
        listing = "\n".join(f"  {n / 1e6:8.2f} MB  {r}" for n, r in biggest)
        raise SystemExit(
            f"archive would be {size:,} bytes, over the {MAX_ARCHIVE_BYTES:,} limit.\n"
            f"Largest entries:\n{listing}\n"
            "Add what does not belong to EXCLUDE, or raise the limit deliberately."
        )
    print(f"wrote {out} ({size:,} bytes; {files} files, {len(paths) - files} dirs)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
