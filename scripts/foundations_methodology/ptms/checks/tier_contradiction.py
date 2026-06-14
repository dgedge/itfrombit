"""checks/tier_contradiction.py — a DRIFT-downgraded claim still labelled Locked/derived (HARD).

For each DRIFT entry that retired/reduced a target (status retired/superseded/withdrawn/
proposition/open) with explicit `supersedes[]` links, check whether the target ANCHOR claim's
block still carries a Locked/derived/prediction label without a co-located retirement marker.
Requires human-confirmed `supersedes[]`, so it cannot fire on coincidental keyword reuse.
"""
from __future__ import annotations

from ..findings import Finding, HARD
from ..registry import index_by_id

CHECK_ID = "tier_contradiction"
LOCK_WORDS = ("[TIER: LOCKED]", "DERIVED", "PREDICTION", "LOCKED")
RETIRE_MARKERS = ("RETIRED", "WITHDRAWN", "SUPERSEDED", "FLAGGED", "REMOVED FROM LOCKED")
DOWNGRADE = ("retired", "superseded", "withdrawn", "proposition", "open")


def run(corpus, registry) -> list:
    idx = index_by_id(registry)
    out = []
    for c in registry:
        if c.kind != "drift-entry" or c.status not in DOWNGRADE:
            continue
        code = c.id.split(":", 1)[1].upper()
        for tgt_id in c.supersedes:
            tgt = idx.get(tgt_id)
            if not tgt:
                continue
            f, ln = tgt.location.get("file", "ANCHOR.md"), tgt.location.get("line", 0)
            if not ln:
                continue
            s, e = corpus.enclosing_block(f, ln)
            block = "\n".join(corpus.lines(f)[s - 1:e]).upper()
            if any(w in block for w in LOCK_WORDS) and not (
                    any(m in block for m in RETIRE_MARKERS) or code in block):
                out.append(Finding(CHECK_ID, HARD, f, ln,
                    f"{tgt_id} still labelled Locked/derived/prediction but {c.id} reduced it "
                    f"({c.status})", f"block {s}-{e}"))
    return out
