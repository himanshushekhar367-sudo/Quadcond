"""Ingestion adapters.

Each adapter is responsible for one source and for being honest about it:
what tier the data belongs to, which condition fields were reported vs
imputed, and what QC flags apply.
"""
from . import g4sp, g4stab_db, g4stab_supp, imseeker, lab, shuffled

ADAPTERS = {
    "g4sp": g4sp,
    "g4stab-db": g4stab_db,
    "g4stab-supp": g4stab_supp,
    "imseeker": imseeker,
    "lab": lab,
    "shuffled": shuffled,
}

__all__ = ["ADAPTERS", "g4sp", "g4stab_db", "g4stab_supp", "imseeker", "lab", "shuffled"]
