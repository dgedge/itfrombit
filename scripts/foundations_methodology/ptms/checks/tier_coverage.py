"""checks/tier_coverage.py — §14 rows lacking a tier tag (SOFT, collapsed in the report).

§16.1 mandates a per-row tier in the §14 falsifiable-signatures table. This is the audit's
documented backlog (~63 of 69 rows untagged today) — soft, and collapsed by the renderer so it
is one line, not 63 alarms. A row intentionally category-covered can be silenced by setting
`tier` on its registry record (or a future `tier_waived` flag).
"""
from __future__ import annotations

from ..findings import Finding, SOFT

CHECK_ID = "tier_coverage"


def run(corpus, registry) -> list:
    out = []
    for c in registry:
        if c.kind == "§14-row" and not c.tier:
            out.append(Finding(CHECK_ID, SOFT,
                c.location.get("file", "ANCHOR.md"), c.location.get("line", 0),
                "§14 row carries no [Tier:] tag (§16.1 mandates a per-row tier)",
                c.location.get("anchor_text", "")))
    return out
