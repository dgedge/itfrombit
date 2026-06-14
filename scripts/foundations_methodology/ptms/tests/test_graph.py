#!/usr/bin/env python3
"""Phase-3 self-test — the claim-graph analyses, on a fixture registry with known structure.

Asserts the exact contamination set (recall + precision: the live ANCHOR citer is flagged; the
retired citer, the DRIFT-entry citer, and the self-ref are excluded; citing a LIVE node is not
contamination) and the frontier ranking. exit 0 iff the graph logic holds.
"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))  # .../ai_methodology

from ptms.registry import Claim
from ptms import graph as G

fails = []
def check(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    if not c:
        fails.append(m)


reg = [
    # item 90 RETIRED, and self-references itself (must be dropped from edges)
    Claim(id="anchor:§15:90", kind="§15-item", status="retired",
          cited_refs=["§15 item 90"], location={"file": "ANCHOR.md", "line": 10}),
    Claim(id="anchor:§15:91", kind="§15-item", status="live",
          location={"file": "ANCHOR.md", "line": 11}),
    Claim(id="anchor:§15:92", kind="§15-item", status="live",          # LIVE ANCHOR citer of 90 -> contamination
          cited_refs=["§15 item 90"], location={"file": "ANCHOR.md", "line": 12}),
    Claim(id="anchor:§15:93", kind="§15-item", status="superseded",    # retired citer of 90 -> excluded
          cited_refs=["§15 item 90"], location={"file": "ANCHOR.md", "line": 13}),
    Claim(id="drift:Z1", kind="drift-entry", status="superseded",      # DRIFT-entry citer of 90 -> excluded
          cited_refs=["§15 item 90"], location={"file": "DRIFT.md", "line": 1}),
    Claim(id="anchor:§15:94", kind="§15-item", status="live",          # cites LIVE 91 -> not contamination
          cited_refs=["§15 item 91"], location={"file": "ANCHOR.md", "line": 14}),
]
g = G.build(reg)
contam = G.contamination(g)

check(set(contam) == {"anchor:§15:90"}, "only the retired item 90 is a contamination target")
check(contam.get("anchor:§15:90") == ["anchor:§15:92"],
      "only the LIVE ANCHOR citer (92) flagged; retired-citer 93, DRIFT Z1, and self-ref excluded")
check("anchor:§15:90" not in g.out["anchor:§15:90"], "self-reference dropped from edges")
check(all(t != "anchor:§15:91" for t in contam), "citing a LIVE node (91) is not contamination")

fr = G.frontier(g, top=5)
check(fr[0] == ("anchor:§15:90", 3), f"frontier top = item 90, in-degree 3 (got {fr[0]})")

check(G.resolve_ref("§15 item 999", g.idx) is None, "unknown §15 item -> unresolved")
check(G.resolve_ref("DRIFT Z1", g.idx) == "drift:Z1", "DRIFT ref resolves")
check(not g.unresolved, "no dangling refs in the fixture")

# DOT smoke: valid digraph wrapper, the contamination edge present
dot = G.to_dot(g, scope="contamination")
check(dot.startswith("digraph ptms {") and dot.rstrip().endswith("}"), "DOT is a well-formed digraph")
check('"anchor:§15:92" -> "anchor:§15:90"' in dot, "contamination edge 92->90 present in DOT")

print()
if fails:
    print(f"SELFTEST(graph): {len(fails)} FAILED")
    sys.exit(1)
print("SELFTEST(graph): exit 0 — graph logic verified (contamination recall+precision, frontier, DOT).")
