"""checks/contamination.py — retraction-contamination surfaced in the read-only suite (SOFT).

A LIVE ANCHOR claim that cites a retired / superseded / withdrawn / flagged node is "built on
retracted ground" — e.g. a claim marked *closed* while a mechanism it depends on was later
ruled a no-go (the item-106-vs-item-87 class). The dependency-graph logic lives in `graph.py`
(tested by `tests/test_graph.py`); this thin wrapper promotes it into `ptms check` so it can no
longer hide in `graph`-only output that nobody runs.

SOFT by design — never blocks a push. Status is auto-extracted, so a *partially*-retired target
(one withdrawn sub-leg of an otherwise-live item) can over-flag, and a citer that references a
retired node precisely to mark it superseded is legitimate. So this SURFACES candidates for
review. Elevate to HARD only once extraction distinguishes full from partial retirement.
"""
from __future__ import annotations

from ..findings import Finding, SOFT
from .. import graph as G

CHECK_ID = "contamination"


def run(corpus, registry) -> list:
    g = G.build(registry)
    contam = G.contamination(g)
    out = []
    for tgt in sorted(contam, key=lambda t: (-len(contam[t]), t)):
        node = g.nodes[tgt]
        f = node.location.get("file", "ANCHOR.md")
        ln = node.location.get("line", 0)
        citers = contam[tgt]
        shown = ", ".join(citers[:6]) + (f" (+{len(citers) - 6})" if len(citers) > 6 else "")
        out.append(Finding(
            CHECK_ID, SOFT, f, ln,
            f"{tgt} ({node.status}) cited by {len(citers)} LIVE claim(s): {shown} "
            f"— built on retracted ground; confirm the cited part is still valid",
            "graph: retraction-contamination"))
    return out
