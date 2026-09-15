"""Where the built viewer lives, if it was built.

The web assets are not part of the Python distribution and are not required to
use it: the CLI, the Python API and the JSON endpoints all work with no
frontend present. When a build *is* present, the service serves it from the same
port as the API, so a deployment is one process and one URL rather than two
servers and a proxy rule.

Resolution order, first hit wins:

1. ``QUADCOND_WEB_ROOT`` -- an explicit path, for a packaged deployment.
2. ``quadcond/_web`` inside the installed package -- where a wheel build copies
   ``web/dist``.
3. ``web/dist`` beside the repository root -- where ``npm run build`` leaves it
   in a source checkout.
"""
from __future__ import annotations

import os
from pathlib import Path

WEB_ROOT_ENV = "QUADCOND_WEB_ROOT"


def candidates() -> list[Path]:
    here = Path(__file__).resolve().parent
    found = []
    env = os.environ.get(WEB_ROOT_ENV)
    if env:
        found.append(Path(env).expanduser())
    found.append(here / "_web")
    found.append(here.parent / "web" / "dist")
    return found


def root() -> Path | None:
    """The directory holding ``index.html``, or ``None`` if nothing is built."""
    for path in candidates():
        if (path / "index.html").is_file():
            return path
    return None


def describe() -> dict:
    r = root()
    return {
        "web_root": str(r) if r else None,
        "served": r is not None,
        "searched": [str(p) for p in candidates()],
        "note": ("The viewer is optional. With no build present the service "
                 "answers its JSON endpoints and serves no pages; run "
                 "`npm --prefix web run build` to produce one."),
    }
