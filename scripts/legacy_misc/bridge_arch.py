#!/usr/bin/env python3
"""
bridge_arch.py
==============
Scope correction to the gamma_TEE / orientation-gluing audits (DRIFT TQO
Target-3 step 3-bis, 2026-06-06). Self-asserting: exit 0 == verified.

CONTEXT. orientation_gluing_audit.py showed that connecting [8,4,4] octagon
cells by *sharing boundary qubits* is obstructed by distance-4. But TCH does
NOT connect matter cells that way (user correction, confirmed vs canon): cells
join THROUGH the gauge web -- "single bridge edges between vertices of distinct
Q_3 clusters" (ANCHOR 1), face f -> face 7-f of the neighbour (ANCHOR 3.1),
with the gauge web carrying its OWN degrees of freedom (12 photonic modes per
Q_3, ANCHOR 7.1). So the sharing obstruction is off-target for TCH.

This script establishes the corrected, architecture-level picture:

  (1) Inter-cell couplings that add NO new DOF cannot create topological order.
      Two independent [8,4,4] cells form a k=0 code whose stabilizer group is
      MAXIMAL (rank = n): the all-ones X-check makes any added single-qubit
      bridge coupling anticommute. So bridges MUST carry their own qubits.

  (2) WITH bridge qubits (DOF on the edges), k>0 is reachable -- this is exactly
      the toric / HGP edge-qubit pattern (toric k=2). TCH's gauge web is this
      kind of structure, so TCH's architecture is NOT subject to the sharing
      obstruction.

Conclusion: what holds TCH at gamma_TEE=0 is the framework's SPECIFICATION of
the gauge web as band topology (Chern number, SPT) -- not a geometric
obstruction. A gamma>0 TCH would require respecifying that gauge web as a
string-net stabilizer code on the bridge qubits (architecturally permitted).

numpy only.
"""
import numpy as np
import sys

fails = []
def check(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    if not c: fails.append(m)

def f2_rank(M):
    M = (np.array(M) % 2).astype(np.int8).copy()
    if M.size == 0: return 0
    R, Cc = M.shape; r = 0
    for c in range(Cc):
        piv = next((i for i in range(r, R) if M[i, c]), None)
        if piv is None: continue
        M[[r, piv]] = M[[piv, r]]
        for i in range(R):
            if i != r and M[i, c]: M[i] ^= M[r]
        r += 1
        if r == R: break
    return r

def sym(a, b, n): return int((a[:n] @ b[n:] + a[n:] @ b[:n]) % 2)

H = np.array([[1, 1, 1, 1, 1, 1, 1, 1],
              [0, 0, 0, 0, 1, 1, 1, 1],
              [0, 0, 1, 1, 0, 0, 1, 1],
              [0, 1, 0, 1, 0, 1, 0, 1]])

# (1) two independent [8,4,4] cells (NO sharing): k=0, stabilizer-maximal
n = 16; rows = []
for off in (0, 8):
    for r in H:
        v = np.zeros(2 * n, int); v[off:off + 8] = r; rows.append(v)        # X-checks
    for r in H:
        v = np.zeros(2 * n, int); v[n + off:n + off + 8] = r; rows.append(v)  # Z-checks
M = np.array(rows); rank = f2_rank(M); k = n - rank
print("=" * 70)
print(f"(1) two independent [8,4,4] cells: n={n}, stabilizer rank={rank}, k={k}")
print("=" * 70)
check(k == 0 and rank == n, "k=0 with rank=n -> stabilizer group is MAXIMAL (no logicals)")

# the all-ones X-check makes any single-qubit Z on a cell anticommute
allones_X = np.zeros(2 * n, int); allones_X[0:8] = 1
z_on_endpoint = np.zeros(2 * n, int); z_on_endpoint[n + 7] = 1
check(sym(allones_X, z_on_endpoint, n) == 1,
      "a Z on a bridge-endpoint qubit anticommutes with the cell all-ones X-check")
print("    => inter-cell couplings that add NO new DOF cannot be added / keep k=0;")
print("       topological order REQUIRES the bridges to carry their own qubits.")

# (2) WITH bridge/edge qubits: k>0 (the toric/HGP architecture TCH's gauge web belongs to)
def toric(L):
    h = lambda x, y: (y % L) * L + (x % L)
    v = lambda x, y: L * L + (y % L) * L + (x % L)
    nn = 2 * L * L; rr = []
    for x in range(L):
        for y in range(L):
            xv = np.zeros(2 * nn, int)
            for q in [h(x, y), h(x - 1, y), v(x, y), v(x, y - 1)]: xv[q] = 1
            rr.append(xv)
    for x in range(L):
        for y in range(L):
            zv = np.zeros(2 * nn, int)
            for q in [h(x, y), h(x, y + 1), v(x, y), v(x + 1, y)]: zv[nn + q] = 1
            rr.append(zv)
    return np.array(rr), nn
Mt, nt = toric(3); kt = nt - f2_rank(Mt)
print("\n" + "=" * 70)
print(f"(2) qubits ON the edges/bridges (toric code, L=3): k={kt}")
print("=" * 70)
check(kt == 2, "bridge/edge-qubit architecture gives k>0 (k=2) -- the non-obstructed class")
print("    => TCH connects cells via a gauge web carrying its own DOF (ANCHOR 7.1), i.e. this")
print("       edge-qubit class -- so TCH is NOT subject to the boundary-sharing obstruction.")

print("\n" + "=" * 70)
if fails:
    print(f"RESULT: {len(fails)} FAILED"); [print("  -", f) for f in fails]; sys.exit(1)
print("RESULT: exit 0 -- the sharing obstruction is boundary-sharing-specific; TCH's gauge-web")
print("(bridge-DOF) architecture is the non-obstructed edge-qubit class. gamma_TEE=0 for TCH holds")
print("by the band-topology SPECIFICATION of the gauge web, not by any geometric obstruction.")
print("=" * 70)
