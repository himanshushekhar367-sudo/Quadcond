#!/usr/bin/env python3
"""Run the suite and record the outcome where the report can read it.

The report needs to state a test result. Three earlier versions of that line
were wrong, each in the same direction:

1. a literal ("52 tests passing") that nobody updated;
2. a count of ``def test_`` in the tree, which reports a healthy number through
   an entirely red suite, because counting functions establishes that tests
   exist and nothing else;
3. a JUnit file this script wrote but never *removed* -- so when pytest was
   unavailable the command correctly failed, the previous run's XML survived
   untouched, and the report went on printing "73 passing" from it.

The third is the subtlest and the worst: a stale green artifact is
indistinguishable from a fresh one, and it is produced by exactly the situation
where you most want to know something is wrong.

So: the old file is deleted before pytest runs, the new one is written to a
temporary path and moved into place only on a successful *run* (green or red --
a red suite is a real result and must be published), and the report treats a
missing file as "not recorded" rather than as anything.

    python scripts/15_run_tests.py
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

OUT = Path("artifacts/pytest-results.xml")


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)

    # Delete first. If this script cannot produce a fresh result, the report
    # must say "not recorded" -- never repeat the last good one.
    stale = OUT.exists()
    OUT.unlink(missing_ok=True)

    fd, tmp_name = tempfile.mkstemp(dir=OUT.parent, suffix=".xml")
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        cmd = [sys.executable, "-m", "pytest", "tests/", "-q",
               f"--junitxml={tmp}"] + sys.argv[1:]
        try:
            rc = subprocess.call(cmd)
        except OSError as exc:
            print(f"could not run pytest: {exc}")
            if stale:
                print(f"removed the previous {OUT} rather than leave a stale "
                      f"result for the report to read")
            sys.exit(2)

        if not tmp.exists() or tmp.stat().st_size == 0:
            print("pytest produced no JUnit output; nothing recorded.")
            if stale:
                print(f"removed the previous {OUT}.")
            sys.exit(rc or 2)

        # Atomic: the report never sees a half-written file.
        tmp.replace(OUT)
        print(f"\nwrote {OUT}")
        if rc != 0:
            print("suite is RED -- the report will say so on its front page.")
        sys.exit(rc)
    finally:
        tmp.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
