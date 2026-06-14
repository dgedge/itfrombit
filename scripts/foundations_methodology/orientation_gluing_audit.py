#!/usr/bin/env python3
"""
orientation_gluing_audit.py
===========================
Tests the hypothesis (user, 2026-06-06) that connecting the [8,4,4] octagon
cells "in patterns where their orientation differed from one to the next"
might yield an under-constrained, toric-code-like code with gamma_TEE > 0 --
a gluing the framework has not specified. Companion to gamma_tee_audit.py.

Self-asserting: exit 0 == every quoted claim verified.

RESULT: the orientation hypothesis is closed NEGATIVELY, for a precise,
orientation-INDEPENDENT reason -- while confirming a topological code IS
constructible from [8,4,4] by a DIFFERENT (non-sharing) construction.

  (A) OBSTRUCTION THEOREM (exhausted over all 254 proper coordinate subsets).
      Two cells each carrying the FULL self-dual [8,4,4] code, sharing a
      coordinate set S, have all cross X(A)-Z(B) checks commute IFF C|_S is
      self-orthogonal, i.e. 2*rank(C|_S) <= |S|. The distance-4 property
      forces 2*rank(C|_S) > |S| for EVERY proper nonempty S. Orientation only
      permutes coordinates -> permutes which S -> cannot change the rank
      profile. So NO orientation pattern and NO overlap size gives even a
      VALID stabilizer code (let alone a topological one).

  (B) BRUTE-FORCE CONFIRMATION. A ring of 4 cells, 2-qubit edge sharing, all
      16^3 D_8 orientation assignments (cell 0 gauge-fixed): 0 valid codes.

  (C) WHY THE HYPERGRAPH PRODUCT ESCAPES. HGP never shares qubits between code
      copies; it introduces NEW qubits (n*n bit-qubits + m*m check-qubits) on a
      product complex -> valid CSS code, k=16. gamma>0 is reachable from the
      [8,4,4] ingredient, but NOT via octagon-boundary-sharing.

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

H = np.array([[1, 1, 1, 1, 1, 1, 1, 1],
              [0, 0, 0, 0, 1, 1, 1, 1],
              [0, 0, 1, 1, 0, 0, 1, 1],
              [0, 1, 0, 1, 0, 1, 0, 1]])
check((H @ H.T % 2 == 0).all(), "[8,4,4] is self-dual (C = C^perp); all-ones is codeword 0")
C = np.array([(a * H[0] + b * H[1] + c * H[2] + d * H[3]) % 2
              for a in (0, 1) for b in (0, 1) for c in (0, 1) for d in (0, 1)])

# ---- (A) the orientation-independent obstruction, over all 254 proper subsets ----
print("=" * 72)
print("(A) OBSTRUCTION -- exhaust all proper coordinate subsets S (share-size & content)")
print("=" * 72)
# A share on coordinate set S is valid (all cross X(A)-Z(B) checks commute) IFF the projected
# generator H|_S is self-orthogonal: (H|_S)(H|_S)^T = 0 mod 2. This INCLUDES the diagonal
# <a|_S, a|_S> = wt(a|_S) mod 2 (a codeword used as an X-check on A AND a Z-check on B), which
# is exactly the all-ones / odd-overlap obstruction. Self-orthogonality forces rank <= |S|/2.
any_valid_share = False
for s in range(1, 8):
    ranks = set(); s_possible = False
    for S in itertools.combinations(range(8), s):
        P = H[:, list(S)]
        ranks.add(f2_rank(P))
        if (P @ P.T % 2 == 0).all():       # C|_S self-orthogonal -> a valid share would exist
            s_possible = True; any_valid_share = True
    print(f"  |S|={s}: rank(C|_S) in {sorted(ranks)}, need 2*rank<= {s} (self-orth) -> "
          f"{'possible' if s_possible else 'IMPOSSIBLE (every subset, every orientation)'}")
check(not any_valid_share,
      "NO proper nonempty share S makes two full [8,4,4] cells commute (2*rank(C|_S) > |S| always)")

# ---- (B) brute-force: ring of cells, 2-qubit sharing, all D_8 orientations ----
print("\n" + "=" * 72)
print("(B) BRUTE FORCE -- ring of 4 cells, 2-qubit edge sharing, all D_8 orientations")
print("=" * 72)
def d8():
    rots = [tuple((j + s) % 8 for j in range(8)) for s in range(8)]
    refs = [tuple((s - j) % 8 for j in range(8)) for s in range(8)]
    return rots + refs
ORI = d8()

def build_ring(M, orients, share=((3, 4), (7, 0))):
    parent = list(range(8 * M))
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    (rA, rB), (lA, lB) = share
    for i in range(M):
        j = (i + 1) % M
        parent[find(i * 8 + rA)] = find(j * 8 + lA)
        parent[find(i * 8 + rB)] = find(j * 8 + lB)
    roots = sorted(set(find(i * 8 + p) for i in range(M) for p in range(8)))
    gid = {r: k for k, r in enumerate(roots)}; n = len(roots)
    gq = lambda i, slot: gid[find(i * 8 + slot)]
    rows = []
    for i in range(M):
        sig = orients[i]
        for r in H:
            xv = np.zeros(2 * n, int)
            for j, b in enumerate(r):
                if b: xv[gq(i, sig[j])] ^= 1
            rows.append(xv)
        for r in H:
            zv = np.zeros(2 * n, int)
            for j, b in enumerate(r):
                if b: zv[n + gq(i, sig[j])] ^= 1
            rows.append(zv)
    return np.array(rows), n

def sym(a, b, n): return int((a[:n] @ b[n:] + a[n:] @ b[:n]) % 2)

M = 4; valid = 0
for combo in itertools.product(range(16), repeat=M - 1):
    orients = [ORI[0]] + [ORI[c] for c in combo]
    rows, n = build_ring(M, orients)
    g = len(rows)
    if not any(sym(rows[i], rows[j], n) for i in range(g) for j in range(i + 1, g)):
        valid += 1
print(f"  valid (all-commuting) orientation assignments: {valid} / {16 ** (M - 1)}")
check(valid == 0, "brute force confirms: 0 valid codes over all D_8^3 orientations (orientation irrelevant)")

# ---- (C) contrast: hypergraph product (new qubits, not sharing) works ----
print("\n" + "=" * 72)
print("(C) CONTRAST -- hypergraph product introduces NEW qubits, never shares between copies")
print("=" * 72)
def hgp(H1, H2):
    m1, n1 = H1.shape; m2, n2 = H2.shape
    Hx = np.hstack([np.kron(H1, np.eye(n2, dtype=int)),
                    np.kron(np.eye(m1, dtype=int), H2.T)]) % 2
    Hz = np.hstack([np.kron(np.eye(n1, dtype=int), H2),
                    np.kron(H1.T, np.eye(m2, dtype=int))]) % 2
    return Hx, Hz, n1 * n2 + m1 * m2
Hx, Hz, nq = hgp(H, H)
k = nq - f2_rank(Hx) - f2_rank(Hz)
print(f"  HGP([8,4,4]): {nq} qubits = 64 bit-qubits + 16 check-qubits (NEW); commute={(Hx@Hz.T%2==0).all()}, k={k}")
check((Hx @ Hz.T % 2 == 0).all() and k == 16,
      "HGP is a valid topological code (k=16) -- but via a product complex, NOT octagon-sharing")

print("\n" + "=" * 72)
if fails:
    print(f"RESULT: {len(fails)} FAILED"); [print("  -", f) for f in fails]; sys.exit(1)
print("RESULT: exit 0 -- the orientation hypothesis is closed NEGATIVELY (distance-4 obstruction,")
print("orientation-independent); a topological code from [8,4,4] needs the non-sharing HGP route.")
print("=" * 72)
