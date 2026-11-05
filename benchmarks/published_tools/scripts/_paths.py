"""Where the benchmark scripts find things.

Every location can be overridden with an environment variable. The defaults
assume this repository's layout, with the third-party tools cloned into
benchmarks/published_tools/tools/.

    QC_TOOLS        cloned tool repositories        default <published_tools>/tools
    QC_WORK         intermediate tables and scores  default <published_tools>/work
    QUADCOND_MODEL  QuadCond model file             default <repo>/artifacts/quadcond_model.joblib
    QUADCOND_DB     QuadCond atlas                  default <repo>/data/atlas_core.db
"""
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent            # benchmarks/published_tools/scripts
REPO = HERE.parents[2]                            # repository root
TOOLS = Path(os.environ.get("QC_TOOLS", HERE.parent / "tools"))
WORK = Path(os.environ.get("QC_WORK", HERE.parent / "work"))
MODEL = Path(os.environ.get("QUADCOND_MODEL", REPO / "artifacts" / "quadcond_model.joblib"))
ATLAS = Path(os.environ.get("QUADCOND_DB", REPO / "data" / "atlas_core.db"))
