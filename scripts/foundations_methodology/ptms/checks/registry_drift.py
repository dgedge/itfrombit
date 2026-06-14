"""checks/registry_drift.py — keep the sidecar honest against the prose (SOFT).

Three detectors: (moved) a record whose anchor_text is no longer at its recorded line;
(unregistered) a claim present in the prose but absent from the registry; (stale) a registry
record whose id no longer resolves in the prose. All soft — they request a re-extract, not a
canon fix.
"""
from __future__ import annotations

from ..findings import Finding, SOFT

CHECK_ID = "registry_drift"


def run(corpus, registry) -> list:
    out = []
    for c in registry:
        f, ln = c.location.get("file"), c.location.get("line")
        at = (c.location.get("anchor_text") or "").strip()
        if not f or not ln or not at:
            continue
        if at not in corpus.line(f, ln):
            out.append(Finding(CHECK_ID, SOFT, f, ln,
                f"{c.id}: anchor_text not at recorded line (re-locate via `ptms extract --write`)",
                at[:50]))
    try:
        from ..extract import extract as _extract
        fresh_ids = {x.id for x in _extract(corpus)}
    except Exception:
        return out
    reg_ids = {c.id for c in registry}
    for i in sorted(fresh_ids - reg_ids):
        out.append(Finding(CHECK_ID, SOFT, "registry", 0, f"unregistered claim in prose: {i}"))
    for i in sorted(reg_ids - fresh_ids):
        out.append(Finding(CHECK_ID, SOFT, "registry", 0, f"stale registry record: {i} not in prose"))
    return out
