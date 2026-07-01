#!/usr/bin/env python3
"""
z3_tetra_audit.py
=================
Reproducibility artifact for the §12.1 / item-110 reconciliation (2026-06-05):
the "no Z^3 tetrahedron" theorem (§12.1) and the He-4 K4 tetrahedron (item 110)
are NOT in contradiction -- they concern different adjacency relations and object
scales. Self-asserting.

(1) The Z^3 single-cell FACE-adjacency graph is bipartite -> no 4 single cells are
    mutually face-adjacent (no single-cell K4). This is exactly §12.1.
(2) A regular K4 tetrahedron DOES embed in Z^3 via FACE-DIAGONALS: the 4 alternating
    cube corners are mutually at distance sqrt(2). So "no tetrahedron in Z^3" holds
    only for face-adjacency, not geometrically.
(3) The gauge web L(Z^3) is NOT bipartite (the 6 edges at a vertex form a K6 with
    triangles) -- so "if L(Z^3) is bipartite" is a false premise; §12.1 is about the
    Z^3 cell graph, not L(Z^3).

Item 110's K4 is 4 COMPOSITE baryons (12 cells), to which the single-cell bipartite
argument does not apply (6 pairwise baryon-contacts among 12 cells need no single-
cell K4). numpy only.
"""
import numpy as np
import itertools

fails = []
def check(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    if not c: fails.append(m)

# (1) Z^3 cell face-adjacency: bipartite + no single-cell K4
def face_adj(a, b): return sum(abs(x-y) for x, y in zip(a, b)) == 1
nbhd = list(itertools.product(range(-1, 2), repeat=3))
bip = all((sum(a) % 2) != (sum(b) % 2) for a in nbhd for b in nbhd if face_adj(a, b))
check(bip, "Z^3 cell face-adjacency graph is bipartite (every face-edge crosses parity classes)")
cells = list(itertools.product(range(3), repeat=3))
has_k4 = any(all(face_adj(q[i], q[j]) for i in range(4) for j in range(i+1, 4))
             for q in itertools.combinations(cells, 4))
check(not has_k4, "no 4 single cells are mutually face-adjacent (no single-cell K4) -- this is §12.1")

# (2) K4 via face-diagonals: 4 alternating cube corners, all pairwise distance sqrt(2)
tet = [(0, 0, 0), (1, 1, 0), (1, 0, 1), (0, 1, 1)]
d = sorted({round(float(np.linalg.norm(np.subtract(a, b))), 4) for a, b in itertools.combinations(tet, 2)})
check(d == [round(2**0.5, 4)], f"4 alternating cube corners are mutually distance sqrt2 = regular K4 in Z^3 (face-diagonal)  {d}")
check(all(sum(t) % 2 == 0 for t in tet), "those 4 corners share parity (independent set in the face graph, K4 only via face-diagonals)")

# (3) L(Z^3) is NOT bipartite: 6 edges at a vertex are mutually adjacent in the line graph (K6 -> triangles)
star = [frozenset([(0, 0, 0), dd]) for dd in
        [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]]
allshare = all(len(e1 & e2) >= 1 for e1, e2 in itertools.combinations(star, 2))
check(allshare, "the 6 Z^3 edges at a vertex are mutually adjacent in L(Z^3) (a K6) -> L(Z^3) NOT bipartite")

print("\n  => §12.1 (single-cell face-adjacency, bipartite, C4) and item 110 (4 composite baryons, K4)")
print("     are about different objects/adjacencies; no contradiction. The open piece is the explicit")
print("     Z^3-embeddable 12-cell packing realizing item 110's 6 inter-baryon interfaces (asserted, not shown).")

print("\n" + "=" * 64)
import sys
if fails:
    print(f"RESULT: {len(fails)} FAILED"); [print("  -", f) for f in fails]; sys.exit(1)
print("RESULT: exit 0 -- §12.1 and item 110 reconciled (different adjacency/scale);")
print("the 'irreconcilable contradiction' and the 'L(Z^3) bipartite' premise both fail.")
print("=" * 64)
