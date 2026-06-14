#!/usr/bin/env python3
r"""Item 99 deep-dive: the two untested joints of the 4.8.8 uniqueness proof chain.

E2 (and its 2026-06-10 finish) re-verified the flatness harness and re-credited the
selector to Lemma 3.1's C4. This script audits the two joints NEITHER pass tested:

  JOINT 1 — the R4-edge reading. Lemma 3.1 builds the internal graph on {LQ, I3, chi, W}
  from "3 explicit operational couplings": (LQ,I3) via CNOT, (chi,W) via R2, and
  (I3,chi) via R4. But R4 = NOT[(LQ=0) AND (I3=0) AND (chi=1)] is a THREE-bit
  constraint — reading it as the single edge (I3,chi) is a choice. Does the C4 (and
  hence 4.8.8) survive the other readings?

  JOINT 2 — the degree-3 budget. Lemma 3.1's 2-regularity comes from "degree-3 lattice
  budget minus 1 external port". Where does degree 3 live in the axiom ledger, and what
  happens without it?

Plus the Lemma-3.3 commutativity premise (no constraint couples a G-bit to a C-bit),
verified from the constraint supports. Self-asserting; exit 0 = verified."""
import itertools as it

# ---------------- constraint supports (from sec 2.2 verbatim) ----------------
SUPPORT = {
    "R1":   {"G0", "G1"},
    "R2":   {"chi", "W"},
    "R3":   {"LQ", "C0", "C1"},
    "R4":   {"LQ", "I3", "chi"},
    "CNOT": {"LQ", "I3"},
}
G_bits, C_bits = {"G0", "G1"}, {"C0", "C1"}
assert all(not (s & G_bits and s & C_bits) for s in SUPPORT.values())
print("[0] Lemma 3.3 premise VERIFIED from the constraint supports: no constraint couples")
print("    a generation bit to a colour bit (R1 is G-only; R3 is LQ+colour; R4 is")
print("    LQ+I3+chi) — the commutativity step's input holds exactly.")

# ---------------- JOINT 1: the R4-edge reading ----------------
# (a) does R4's own dependency structure single out (I3,chi)? R4 forbids exactly one
# cell of the 8 on (LQ, I3, chi); a single forbidden cell makes EVERY pair conditionally
# dependent — no pair is singled out by the constraint's correlations alone.
cells = [c for c in it.product((0, 1), repeat=3) if not (c[0] == 0 and c[1] == 0 and c[2] == 1)]
assert len(cells) == 7
def pair_conditionally_dependent(i, j):
    k = 3 - i - j
    for ck in (0, 1):                                  # condition on the third bit
        sub = [(c[i], c[j]) for c in cells if c[k] == ck]
        for a in (0, 1):
            for b in (0, 1):
                pa = sum(1 for s in sub if s[0] == a) / len(sub)
                pb = sum(1 for s in sub if s[1] == b) / len(sub)
                pab = sum(1 for s in sub if s == (a, b)) / len(sub)
                if abs(pab - pa * pb) > 1e-12:
                    return True
    return False
dep = {(i, j): pair_conditionally_dependent(i, j) for i, j in [(0, 1), (0, 2), (1, 2)]}
assert all(dep.values())
print("\n[1a] R4 forbids exactly one cell of 8 on (LQ, I3, chi); ALL three pairs are")
print("     conditionally dependent — the (I3,chi) edge is NOT singled out by R4's own")
print("     correlation structure. The edge assignment is a READING, so test all readings:")

# (b) enumerate the readings: mandatory edges (LQ,I3) [CNOT], (chi,W) [R2]; R4 -> one of
# its three internal pairs; then complete to a 2-regular SIMPLE graph if possible.
V = ["LQ", "I3", "chi", "W"]
mandatory = [("LQ", "I3"), ("chi", "W")]
def complete_2regular(edges):
    deg = {v: sum(v in e for e in edges) for v in V}
    if any(d > 2 for d in deg.values()) or len(set(map(frozenset, edges))) < len(edges):
        return None                                    # multi-edge or over-degree
    need = [v for v in V if deg[v] == 1]
    if len(need) == 2 and frozenset(need) not in set(map(frozenset, edges)):
        return edges + [tuple(need)]
    return None
results = {}
for r4pair in [("LQ", "I3"), ("LQ", "chi"), ("I3", "chi")]:
    g = complete_2regular(mandatory + [r4pair])
    if g is None:
        results[r4pair] = ("EXCLUDED (multi-edge / over-degree)", None, None)
        continue
    # classify the 2-regular graph: C4 iff connected
    adj = {v: [u for e in g for u in e if v in e and u != v] for v in V}
    seen, stack = {V[0]}, [V[0]]
    while stack:
        for u in adj[stack.pop()]:
            if u not in seen:
                seen.add(u); stack.append(u)
    cyc = ["LQ"]
    while len(cyc) < 4:
        cyc.append(next(u for u in adj[cyc[-1]] if u != (cyc[-2] if len(cyc) > 1 else None)))
    antip = {frozenset((cyc[0], cyc[2])), frozenset((cyc[1], cyc[3]))}
    results[r4pair] = ("C4" if len(seen) == 4 else "2xC2", g, antip)
for r4pair, (verdict, g, antip) in results.items():
    extra = f"  antipodal pairs {sorted(sorted(p) for p in antip)}" if antip else ""
    print(f"     R4 -> edge {r4pair}: {verdict}{extra}")
assert results[("LQ", "I3")][0].startswith("EXCLUDED")
assert results[("LQ", "chi")][0] == "C4" and results[("I3", "chi")][0] == "C4"
a1, a2 = results[("I3", "chi")][2], results[("LQ", "chi")][2]
assert a1 == {frozenset(("LQ", "chi")), frozenset(("I3", "W"))}
assert a2 == {frozenset(("LQ", "W")), frozenset(("I3", "chi"))}
print("     -> BOTH viable readings give C4: Lemma 3.1 and the {4,f,f}->f=8 uniqueness")
print("        are ROBUST to the R4-edge ambiguity (a freedom no prior pass tested).")
print("        BUT the antipodal pairs differ: canon's axis anchor (horizontal {chi,LQ},")
print("        vertical {I3,W}) holds ONLY under the (I3,chi) reading — the axis-")
print("        assignment layer is convention-conditioned on 'R4 contributes the edge")
print("        between the two bits not already coupled to LQ'. Named, presentation-level.")

# ---------------- JOINT 2: the degree-3 budget ----------------
# 3-regular simple graphs on 4 vertices: K4 is the unique one (each vertex degree 3).
edges_all = [tuple(sorted(e)) for e in it.combinations(V, 2)]
threereg = []
for r in range(len(edges_all) + 1):
    for sub in it.combinations(edges_all, r):
        deg = {v: sum(v in e for e in sub) for v in V}
        if all(d == 3 for d in deg.values()):
            threereg.append(sub)
assert len(threereg) == 1 and len(threereg[0]) == 6   # K4, girth 3
print("\n[2] Degree-3 budget (d_int = 3-1 = 2): with a degree-4 budget instead (d_int = 3)")
print("    the internal graph is forced to K4 (unique 3-regular on 4 vertices, girth 3) —")
print("    no C4 plaquette, the 4.8.8 chain collapses. Trivalence is LOAD-BEARING, and its")
print("    status is: DISCLOSED PRE-GEOMETRIC INPUT (sec 0 table: '3-regular qubit network';")
print("    sec 1: 'crystallises spontaneously from the 3-regular qubit network') — but it is")
print("    ABSENT from item 99's own postulate list (Postulates 0-3). The input ledger of")
print("    the uniqueness theorem is incomplete by one named, already-disclosed premise.")

# ---------------- flatness harness (E2, re-asserted) ----------------
from fractions import Fraction as Fr
assert [f for f in range(3, 100) if Fr(1) - Fr(3, 2) + Fr(1, 4) + Fr(2, f) == 0] == [8]

print("""
VERDICT (item 99 after the deep-dive):
  * GREEN STANDS and is STRENGTHENED: the C4 — the load-bearing selector per E2's
    re-credit — is robust to the R4-edge ambiguity (both viable readings yield C4;
    the third is structurally excluded), and the Lemma-3.3 commutativity premise is
    verified exactly from the constraint supports.
  * Two findings to record: (i) the axis-anchoring layer (which antipodal pair is
    horizontal/vertical) is conditioned on the (I3,chi) reading of R4 — a named
    presentation-level convention, not load-bearing for uniqueness; (ii) the theorem's
    input list should include TRIVALENCE — disclosed at sec 0/sec 1 but missing from
    the item's postulate enumeration; without it (degree-4 budget) the internal graph
    is K4 and the 4.8.8 chain collapses.
  * Net: 4.8.8 uniqueness = derivation from FIVE disclosed inputs (spatial extent,
    8-bit partition, R1-R4, weak CNOT, trivalence) + one presentation convention.
    No reverse-engineering; the input ledger is now complete.""")
print("ALL ASSERTIONS PASSED — every number above is verified. exit 0")
