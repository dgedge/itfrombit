"""graph.py — the claim dependency graph + analyses + Graphviz DOT (Phase 3). Read-only, stdlib.

Nodes = registry claims; edges = `cited_refs` resolved to registry ids (A references B). Two
analyses: retraction-contamination (LIVE ANCHOR claims citing a retired/superseded claim — built
on retracted ground) and the open frontier (the §15 targets ranked by how much depends on them).
Exploratory: no pass/fail, no exit code (the `graph` command always exits 0).
"""
from __future__ import annotations

import re
from collections import Counter

from .registry import index_by_id

ITEM = re.compile(r"§15 item (\d+)")
DREF = re.compile(r"DRIFT ([A-Z]\d+)")
SEC = re.compile(r"§(\d+(?:\.\d+)?)")
RETIRED = ("retired", "superseded", "withdrawn", "flagged")
ANCHOR_KINDS = ("§15-item", "§14-row", "anchor-section")


def resolve_ref(ref: str, idx: dict):
    """Map a cross-ref string to a registry id (or None). §x.y falls back to its parent §x."""
    m = ITEM.match(ref)
    if m:
        cid = f"anchor:§15:{m.group(1)}"
        return cid if cid in idx else None
    m = DREF.match(ref)
    if m:
        cid = f"drift:{m.group(1)}"
        return cid if cid in idx else None
    m = SEC.match(ref)
    if m:
        full = f"anchor:sec:{m.group(1)}"
        if full in idx:
            return full
        major = f"anchor:sec:{m.group(1).split('.')[0]}"
        return major if major in idx else None
    return None


class Graph:
    def __init__(self, registry):
        self.idx = index_by_id(registry)
        self.nodes = {c.id: c for c in registry}
        self.out = {c.id: set() for c in registry}   # A -> {ids A references}
        self.inn = {c.id: set() for c in registry}   # B <- {ids referencing B}
        self.unresolved = []                          # (citer_id, raw_ref)
        for c in registry:
            for ref in c.cited_refs:
                t = resolve_ref(ref, self.idx)
                if t is None:
                    self.unresolved.append((c.id, ref))
                    continue
                if t == c.id:
                    continue                          # self-reference
                self.out[c.id].add(t)
                self.inn[t].add(c.id)

    def edge_count(self):
        return sum(len(v) for v in self.out.values())


def build(registry) -> Graph:
    return Graph(registry)


def _sources(g: Graph):
    """contamination sources: nodes that are retired-ish in status, plus any supersedes-target."""
    src = {nid for nid, c in g.nodes.items() if c.status in RETIRED}
    src |= {t for c in g.nodes.values() for t in c.supersedes if t in g.nodes}
    return src


def contamination(g: Graph) -> dict:
    """{retired_target_id: sorted [live ANCHOR claim ids citing it]} — excludes retired citers,
    DRIFT-entry citers (their job is to discuss retirements), and self-refs."""
    out = {}
    for tgt in _sources(g):
        citers = []
        for a in g.inn.get(tgt, ()):
            ca = g.nodes[a]
            if ca.status in RETIRED or ca.kind not in ANCHOR_KINDS:
                continue
            if tgt in ca.supersedes:                  # defensive (DRIFT entries already excluded)
                continue
            citers.append(a)
        if citers:
            out[tgt] = sorted(citers)
    return out


def frontier(g: Graph, top: int = 15):
    """§15 open-target items ranked by in-degree (most depended-on => most blocking)."""
    items = [(nid, len(g.inn[nid])) for nid, c in g.nodes.items() if c.kind == "§15-item"]
    items.sort(key=lambda x: (-x[1], x[0]))
    return items[:top]


def render_md(g: Graph, top: int = 15) -> str:
    contam = contamination(g)
    total = sum(len(v) for v in contam.values())
    L = ["# PTMS claim graph — analysis", "",
         f"nodes: {len(g.nodes)}   resolved ref-edges: {g.edge_count()}   "
         f"unresolved refs: {len(g.unresolved)}", ""]

    L.append(f"## Retraction-contamination — live ANCHOR claims citing a retired claim "
             f"({total} citations across {len(contam)} retired targets)")
    L.append("_Candidates for review: a citer leans on a node whose status is retired/superseded/"
             "withdrawn/flagged. Status is auto-extracted, so a partially-retired item (one withdrawn "
             "sub-leg in an otherwise-live item) may over-flag — confirm against the cited part._")
    if not contam:
        L.append("- none")
    for tgt in sorted(contam, key=lambda t: (-len(contam[t]), t)):
        c = g.nodes[tgt]
        cites = contam[tgt]
        shown = ", ".join(cites[:8]) + (f" … (+{len(cites) - 8})" if len(cites) > 8 else "")
        L.append(f"- **{tgt}** ({c.status}, {c.location.get('file', '')}:{c.location.get('line', 0)}) "
                 f"— cited by {len(cites)}: {shown}")
    L.append("")

    L.append(f"## Open frontier / blockers — §15 targets by dependents (top {top})")
    for nid, deg in frontier(g, top):
        c = g.nodes[nid]
        L.append(f"- {nid} (in-degree {deg}) — {c.location.get('anchor_text', '')[:60]}")
    L.append("")

    L.append("## Summary")
    L.append(f"- by kind: {dict(Counter(c.kind for c in g.nodes.values()))}")
    L.append(f"- by status: {dict(Counter(c.status for c in g.nodes.values()))}")
    if g.unresolved:
        sample = "; ".join(f"{a} → {r}" for a, r in g.unresolved[:5])
        L.append(f"- dangling refs ({len(g.unresolved)}): {sample}")
    return "\n".join(L)


_STATUS_COLOR = {"live": "white", "retired": "#ffb3b3", "withdrawn": "#ffb3b3",
                 "superseded": "#ffd9b3", "flagged": "#fff2b3", "open": "#b3d9ff",
                 "proposition": "#e6e6e6"}
_KIND_SHAPE = {"§15-item": "ellipse", "§14-row": "box", "drift-entry": "diamond",
               "anchor-section": "folder"}


def _short(nid: str) -> str:
    """Compact node label for the DOT figure: anchor:§15:86 -> §15:86, drift:M9 -> M9."""
    return nid.replace("anchor:", "").replace("drift:", "")


def to_dot(g: Graph, scope: str = "contamination", top: int = 15) -> str:
    contam = contamination(g)
    contam_edges = {(a, tgt) for tgt, v in contam.items() for a in v}
    if scope == "contamination":
        # top-N most-cited retired targets + their live citers; draw ONLY the contamination edges
        # (the incidental gray dependencies among these nodes are dropped for legibility).
        targets = [t for t, _ in sorted(((t, len(v)) for t, v in contam.items()),
                                         key=lambda x: (-x[1], x[0]))[:top]]
        keep = set(targets) | {a for t in targets for a in contam[t]}
        edges = [(a, t, True) for t in targets for a in contam[t]]
    elif scope == "frontier":
        fr = [nid for nid, _ in frontier(g, top)]
        keep = set(fr) | {a for nid in fr for a in g.inn[nid]}
        edges = [(a, b, (a, b) in contam_edges) for a in sorted(keep) for b in sorted(g.out[a]) if b in keep]
    else:
        keep = set(g.nodes)
        edges = [(a, b, (a, b) in contam_edges) for a in sorted(keep) for b in sorted(g.out[a]) if b in keep]
    graph_attrs = ("rankdir=TB; pack=true; nodesep=0.3; ranksep=0.6;"   # bipartite clusters, packed landscape
                   if scope == "contamination" else "rankdir=LR; nodesep=0.25; ranksep=0.7;")
    L = ["digraph ptms {", f"  {graph_attrs}",
         '  node [style=filled, fontsize=10, fontname="Helvetica"];']
    for nid in sorted(keep):
        c = g.nodes[nid]
        L.append(f'  "{nid}" [shape={_KIND_SHAPE.get(c.kind, "ellipse")}, '
                 f'fillcolor="{_STATUS_COLOR.get(c.status, "white")}", label="{_short(nid)}"];')
    for a, b, is_contam in edges:
        attr = ' [color="red", penwidth=2]' if is_contam else ' [color="gray75"]'
        L.append(f'  "{a}" -> "{b}"{attr};')
    L.append("}")
    return "\n".join(L)
