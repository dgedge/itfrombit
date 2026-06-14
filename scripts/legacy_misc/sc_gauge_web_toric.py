#!/usr/bin/env python3
"""
sc_gauge_web_toric.py
=====================
The topological ("string-net") alternative for the TCH gauge sector, computed
on the framework's OWN gauge lattice (ANCHOR 7.3: the dual of the truncated-cube
gauge web is the Simple Cubic Bravais lattice). Self-asserting: exit 0 == verified.

Construction: the Z_2 (toric) gauge code on an LxLxL periodic SC lattice.
  qubits   : on EDGES (3 per vertex: +x,+y,+z), N = 3 L^3
  stars A_v: X on the 6 edges at vertex v        (X-type, 1 per vertex)
  plaq B_p : Z on the 4 edges of each unit square (Z-type, 3 per vertex)

RESULT (verified):
  - k = 3 logical qubits, CONSTANT in L (L=2,3,4) -> intrinsic topological order
    (size-independent ground-state degeneracy = the hallmark gamma>0 signature).
    k = 3 = b_1(T^3), the first Betti number of the 3-torus -- a pure topological
    fact of the periodic lattice. (It is NOT "3 generations"; reading it that way
    would repeat the K=16 numerology error -- see DRIFT, 2026-06-06.)
  - Logical Z = a 1D non-contractible FLUX LOOP (weight L) -- Wilson-loop-like,
    matching the framework's native flux-string language.
  - Logical X = a 2D MEMBRANE (weight L^2), the e/m-dual partner.

HONEST PHYSICS BOUND (NOT a free swap for the framework -- [proposition], not computed here):
  The Z_2 toric phase is GAPPED -- it has NO massless photon. The framework's
  emergent photon is the massless mode of a U(1) Chern/Wilson gauge sector
  (ANCHOR 7.1/7.3), whose deconfined Coulomb phase has gamma_TEE = 0. So massless
  photon (U(1) Coulomb, gamma=0) and toric topological order (Z_2, gamma>0, gapped)
  are MUTUALLY EXCLUSIVE in this gauge sector. gamma_TEE=0 for TCH is therefore not
  a free specification choice to flip: it is tied to having a massless photon at all.

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

def build(L):
    e = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    def vid(x, y, z): return ((x % L) * L + (y % L)) * L + (z % L)
    def eid(x, y, z, d): return 3 * vid(x, y, z) + d
    N = 3 * L ** 3
    stars, plaqs = [], []
    for x in range(L):
        for y in range(L):
            for z in range(L):
                s = np.zeros(N, int)
                s[eid(x, y, z, 0)] = 1; s[eid(x - 1, y, z, 0)] = 1
                s[eid(x, y, z, 1)] = 1; s[eid(x, y - 1, z, 1)] = 1
                s[eid(x, y, z, 2)] = 1; s[eid(x, y, z - 1, 2)] = 1
                stars.append(s)
                for (d1, d2) in [(0, 1), (1, 2), (2, 0)]:
                    p = np.zeros(N, int)
                    p[eid(x, y, z, d1)] = 1
                    p[eid(x + e[d1][0], y + e[d1][1], z + e[d1][2], d2)] = 1
                    p[eid(x + e[d2][0], y + e[d2][1], z + e[d2][2], d1)] = 1
                    p[eid(x, y, z, d2)] = 1
                    plaqs.append(p)
    return np.array(stars), np.array(plaqs), N, eid

# --- k = 3 constant in L => topological order ---
print("=" * 72)
print("Z_2 toric code on the SC gauge web (ANCHOR 7.3): logical-qubit count")
print("=" * 72)
for L in (2, 3, 4):
    A, B, N, eid = build(L)
    comm = (A @ B.T % 2 == 0).all()
    k = N - f2_rank(A) - f2_rank(B)
    print(f"  L={L}: N={N} edge-qubits, {A.shape[0]} stars, {B.shape[0]} plaquettes; CSS-commute={comm}, k={k}")
    check(comm, f"L={L}: stars (X) and plaquettes (Z) commute -> valid CSS code")
    check(k == 3, f"L={L}: k=3 (constant in L -> intrinsic topological order, gamma>0)")

# --- logical operators on L=3 ---
L = 3; A, B, N, eid = build(L)
zloop = np.zeros(N, int)
for x in range(L): zloop[eid(x, 0, 0, 0)] = 1                 # 1D Z flux loop wrapping x
xmem = np.zeros(N, int)
for y in range(L):
    for z in range(L): xmem[eid(0, y, z, 0)] = 1             # 2D X membrane at x=0

print("\n" + "=" * 72)
print("Logical operator geometry (L=3)")
print("=" * 72)
zl_comm_stars = sum(int(A[i] @ zloop % 2) for i in range(A.shape[0])) == 0
zl_not_stab = f2_rank(np.vstack([B, zloop])) != f2_rank(B)
print(f"  Z-loop (1D flux string, weight {int(zloop.sum())}): commutes-with-stars={zl_comm_stars}, non-contractible={zl_not_stab}")
check(zl_comm_stars and zl_not_stab, "logical Z = non-contractible 1D FLUX LOOP (Wilson-loop-like)")

xm_comm_plaq = sum(int(xmem @ B[i] % 2) for i in range(B.shape[0])) == 0
xm_not_stab = f2_rank(np.vstack([A, xmem])) != f2_rank(A)
xm_anti_zl = int(xmem @ zloop % 2) == 1
print(f"  X-membrane (2D, weight {int(xmem.sum())}): commutes-with-plaqs={xm_comm_plaq}, non-contractible={xm_not_stab}, anticommutes-with-its-Z-loop={xm_anti_zl}")
check(xm_comm_plaq and xm_not_stab and xm_anti_zl, "logical X = 2D MEMBRANE, e/m-dual to the Z flux loop")

print("\n" + "=" * 72)
if fails:
    print(f"RESULT: {len(fails)} FAILED"); [print("  -", f) for f in fails]; sys.exit(1)
print("RESULT: exit 0 -- a Z_2 string-net on the SC gauge web is topologically ordered")
print("(k=3, 1D flux loops + 2D membranes). This is the concrete gamma>0 route on the")
print("framework's own lattice. PHYSICS BOUND (proposition): this phase is GAPPED -- no")
print("massless photon; the framework's U(1) photon needs the Coulomb phase (gamma=0). So")
print("topological order and a massless photon are mutually exclusive in this gauge sector.")
print("=" * 72)
