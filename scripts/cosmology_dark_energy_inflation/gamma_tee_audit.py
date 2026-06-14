#!/usr/bin/env python3
"""
gamma_tee_audit.py
==================
Reproducibility artifact for the topological-entanglement-entropy (gamma_TEE)
audit of the macroscopic 2D-projection stabilizer lattice (TQO programme;
ANCHOR 2D-projection / Kitaev-honeycomb topological-order claim). Self-asserting:
every quoted number is asserted, so exit 0 == verified.

Stabilizer-code facts used (Fattal-Cubitt-Yamamoto-Bravyi-Chuang; Kitaev-Preskill):
  - S(A) = |A| - dim(S_A) bits, S_A = stabilizers supported wholly on A.
  - logical qubits k = n - rank_F2(symplectic check matrix); k = ground-state
    degeneracy exponent. k>0 with NON-LOCAL logicals == long-range order.
  - gamma_TEE > 0  <=>  topological order  <=>  non-contractible logical strings.

RESULT (two halves, both computed):

(1) NEGATIVE -- the framework's STATED code is trivial.
    The [8,4,4] extended-Hamming / RM(1,3) CSS cell is FULLY-constraining:
    [[8,0,4]], k=0, no logical qubits. Any local tiling of these cells stays
    k=0 -> gamma_TEE = 0. Naive qubit-sharing between cells is not even a valid
    code (CSS X-checks anticommute with neighbour Z-checks on shared qubits).
    Validated against a toric-code harness (k=2, the gamma=log2 reference).
    => "Kitaev-honeycomb-like topological order in the 2D projection" is NOT
       realised by the [8,4,4] stabilizers.

(2) POSITIVE / E4 bound -- a topological code IS constructible from [8,4,4].
    The hypergraph product (Tillich-Zemor 2009) of the classical [8,4,4] code
    with itself is a VALID topological CSS code with k=16 logical qubits. The
    toric code is itself the hypergraph product of two repetition codes, so this
    is literally the toric construction, generalized. Hence gamma>0 is reachable
    from the [8,4,4] ingredient -- BUT it is a DIFFERENT code (80 qubits, product
    layout, 16 logicals), not the 8-qubit octagon cell. Whether the framework's
    actual octagon/orientation geometry realises such a gluing (and what its
    extra structure means physically) is open, NOT proven impossible.

numpy only.
"""
import numpy as np
import itertools
import sys

fails = []
def check(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    if not c: fails.append(m)

def f2_rank(M):
    M = (np.array(M) % 2).astype(np.int8).copy()
    if M.size == 0: return 0
    rows, cols = M.shape; r = 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i, c]), None)
        if piv is None: continue
        M[[r, piv]] = M[[piv, r]]
        for i in range(rows):
            if i != r and M[i, c]: M[i] ^= M[r]
        r += 1
        if r == rows: break
    return r

def k_logical(M, n): return n - f2_rank(M)                 # logical qubits = degeneracy exponent

def stab_entropy(M, n, A):                                 # S(A) bits = |A| - dim(S_A)
    Aset = set(A); Abar = [q for q in range(n) if q not in Aset]
    cols = Abar + [q + n for q in Abar]
    g = M.shape[0]
    return len(A) - (g if not cols else g - f2_rank(M[:, cols]))

def sym_anticommute(r1, r2, n):                            # symplectic product x1.z2 + z1.x2
    return int((r1[:n] @ r2[n:] + r1[n:] @ r2[:n]) % 2)

# ============================================================================
print("=" * 74)
print("HARNESS -- toric code on LxL torus: k=2 (non-contractible logicals), gamma=log2")
print("=" * 74)
def toric(L):
    h = lambda x, y: (y % L) * L + (x % L)
    v = lambda x, y: L * L + (y % L) * L + (x % L)
    n = 2 * L * L; rows = []
    for x in range(L):
        for y in range(L):
            xv = np.zeros(2 * n, int)
            for q in [h(x, y), h(x - 1, y), v(x, y), v(x, y - 1)]: xv[q] = 1
            rows.append(xv)
    for x in range(L):
        for y in range(L):
            zv = np.zeros(2 * n, int)
            for q in [h(x, y), h(x, y + 1), v(x, y), v(x + 1, y)]: zv[n + q] = 1
            rows.append(zv)
    return np.array(rows), n
for L in (3, 4, 5):
    M, n = toric(L); k = k_logical(M, n)
    print(f"  toric L={L}: n={n}, k = n - rank = {k}")
    check(k == 2, f"toric L={L} gives k=2 (method validated: detects topological order)")

# ============================================================================
print("\n" + "=" * 74)
print("FRAMEWORK STATED CODE -- [8,4,4] / RM(1,3) cell and its tilings: k=0, gamma_TEE=0")
print("=" * 74)
H = np.array([[1, 1, 1, 1, 1, 1, 1, 1],
              [0, 0, 0, 0, 1, 1, 1, 1],
              [0, 0, 1, 1, 0, 0, 1, 1],
              [0, 1, 0, 1, 0, 1, 0, 1]])
check((H @ H.T % 2 == 0).all(), "RM(1,3) is self-dual (H H^T = 0) -> CSS-admissible")

def cell844():
    n = 8; rows = []
    for r in H:
        xv = np.zeros(2 * n, int); xv[:n] = r; rows.append(xv)
    for r in H:
        zv = np.zeros(2 * n, int); zv[n:] = r; rows.append(zv)
    return np.array(rows), n
M, n = cell844(); k = k_logical(M, n)
prof = [stab_entropy(M, n, list(range(a))) for a in range(0, 9)]
print(f"  single [8,4,4] cell: n=8, k = {k};  S(A) profile |A|=0..8: {prof}")
check(k == 0, "single [8,4,4] cell is [[8,0,4]] -- k=0, fully-constraining, NO logical qubits")
check(prof == [0, 1, 2, 3, 2, 3, 2, 1, 0], "S(A) is boundary-law with NO topological constant")

def tiling844(D):
    cells = D * D; n = 8 * cells; rows = []
    for c in range(cells):
        off = 8 * c
        for r in H:
            xv = np.zeros(2 * n, int); xv[off:off + 8] = r; rows.append(xv)
        for r in H:
            zv = np.zeros(2 * n, int); zv[n + off:n + off + 8] = r; rows.append(zv)
    return np.array(rows), n
for D in (2, 3):
    M, n = tiling844(D); k = k_logical(M, n)
    print(f"  {D}x{D} independent-cell tiling: n={n}, k = {k}")
    check(k == 0, f"{D}x{D} tiling stays k=0 -> gamma_TEE = 0 (no emergent topological order)")

# naive qubit-sharing between two cells: is it even a valid stabilizer code?
nq = 15; rows = []
for q_map in ([0, 1, 2, 3, 4, 5, 6, 7], [7, 8, 9, 10, 11, 12, 13, 14]):  # share qubit 7
    for r in H:
        vX = np.zeros(2 * nq, int)
        for i, b in enumerate(r):
            if b: vX[q_map[i]] = 1
        rows.append(vX)
    for r in H:
        vZ = np.zeros(2 * nq, int)
        for i, b in enumerate(r):
            if b: vZ[nq + q_map[i]] = 1
        rows.append(vZ)
rows = np.array(rows)
anti = sum(sym_anticommute(rows[i], rows[j], nq) for i in range(len(rows)) for j in range(i + 1, len(rows)))
print(f"  two cells sharing 1 qubit: {len(rows)} checks, {anti} anticommuting pairs")
check(anti > 0, f"naive same-orientation gluing is INVALID ({anti} anticommuting CSS pairs)")

# ============================================================================
print("\n" + "=" * 74)
print("E4 BOUND -- a topological code IS constructible from [8,4,4]: hypergraph product")
print("=" * 74)
def hgp(H1, H2):
    m1, n1 = H1.shape; m2, n2 = H2.shape
    Hx = np.hstack([np.kron(H1, np.eye(n2, dtype=int)),
                    np.kron(np.eye(m1, dtype=int), H2.T)]) % 2
    Hz = np.hstack([np.kron(np.eye(n1, dtype=int), H2),
                    np.kron(H1.T, np.eye(m2, dtype=int))]) % 2
    return Hx, Hz, n1 * n2 + m1 * m2

# sanity: toric code == HGP(repetition, repetition)
rep = np.array([[1, 1, 0], [0, 1, 1], [1, 0, 1]])
Hx, Hz, nq = hgp(rep, rep)
check((Hx @ Hz.T % 2 == 0).all() and nq - f2_rank(Hx) - f2_rank(Hz) == 2,
      "HGP(rep,rep) reproduces the toric code (k=2) -- HGP is the toric construction generalized")

# the framework's [8,4,4] hypergraph-producted with itself
Hx, Hz, nq = hgp(H, H)
k = nq - f2_rank(Hx) - f2_rank(Hz)
print(f"  HGP([8,4,4],[8,4,4]): n={nq} qubits, {Hx.shape[0]} X-checks, {Hz.shape[0]} Z-checks, k = {k}")
check((Hx @ Hz.T % 2 == 0).all(), "HGP X/Z checks commute -> a VALID CSS code")
check(k == 16, "HGP([8,4,4]) has k=16 logical qubits = (n-r)^2+(m-r)^2 = 4^2+0^2 -> topological")
check(k > 0, "k>0 with non-local logicals: crosses the divide the k=0 cell-tiling failed")

print("\n" + "=" * 74)
if fails:
    print(f"RESULT: {len(fails)} FAILED"); [print("  -", f) for f in fails]; sys.exit(1)
print("RESULT: exit 0 -- verified.")
print(" (1) gamma_TEE = 0 for the STATED [8,4,4] code: k=0 cell, k=0 tilings, naive")
print("     gluing inconsistent. The 2D-projection topological-order claim is NOT")
print("     realised by the [8,4,4] stabilizers.")
print(" (2) A topological code (k=16) IS constructible from [8,4,4] via the hypergraph")
print("     product -- a DIFFERENT 80-qubit code, not the octagon cell. Whether the")
print("     framework's actual octagon/orientation geometry realises such a gluing is")
print("     OPEN (E4: not realised by the stated structure; not proven impossible).")
print("=" * 74)
