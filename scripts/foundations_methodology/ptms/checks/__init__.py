"""checks/ — read-only check modules for `ptms check`. Each exposes CHECK_ID and
run(corpus, registry) -> list[Finding].

discover() returns ONLY the modules listed in ORDER (the read-only suite), in leverage order.
A check must be registered in ORDER to run under `check` — this deliberately excludes
`evidence_regression` (Phase 2), whose run() signature takes a manifest and which is the
EXECUTING command (`ptms verify` imports it directly), not part of the read-only suite.
"""
from __future__ import annotations

import importlib

# leverage order, and the allow-list of read-only checks (the report re-sorts by severity anyway)
ORDER = [
    "sync_propagation",     # HARD — killer feature
    "multipart",            # HARD
    "tier_contradiction",   # HARD
    "citation_integrity",   # mixed
    "numeric_tieouts",      # HARD
    "tier_coverage",        # SOFT (collapsed)
    "registry_drift",       # SOFT
    "contamination",        # SOFT — live claims on retracted ground (graph-based; surfaced, not blocking)
]


def discover() -> list:
    out = []
    for name in ORDER:
        m = importlib.import_module(f"{__name__}.{name}")
        if hasattr(m, "run") and hasattr(m, "CHECK_ID"):
            out.append(m)
    return out
