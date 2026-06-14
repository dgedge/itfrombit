"""checks/multipart.py — multi-part completeness (HARD).

A DRIFT entry that retires K>=2 distinct things must register at least K signature_strings,
so the sync check sweeps every leg. Fires only on a HUMAN-set multipart_count (the extractor's
guess lives in auto_hint, not here) — so this never nags about un-curated entries; it catches
the M9-style failure where one leg of a multi-part retirement is left unswept.
"""
from __future__ import annotations

from ..findings import Finding, HARD

CHECK_ID = "multipart"


def run(corpus, registry) -> list:
    out = []
    for c in registry:
        if c.kind != "drift-entry" or c.multipart_count < 2:
            continue
        if len(c.signature_strings) < c.multipart_count:
            out.append(Finding(CHECK_ID, HARD,
                c.location.get("file", "DRIFT.md"), c.location.get("line", 0),
                f"{c.id} declares {c.multipart_count} retired legs but only "
                f"{len(c.signature_strings)} signature(s) registered → leg(s) unswept",
                "multi-part under-propagation risk (the M9 failure mode)"))
    return out
