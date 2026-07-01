#!/usr/bin/env python3
"""
TQO Target 3 (topological entanglement entropy) -- STEP 1: the [8,4,4] CSS matter cell.

Exact stabiliser computation, no hardcoding. Establishes:
  (1) the [8,4,4] CSS cell (X- and Z-stabilisers both from the self-dual RM(1,3)
      generator matrix GEN, cf. smg_code_projection.py / ps_su4_fock_smg.py) is a
      [[8,0,4]] code: k = 8 - rank(Hx) - rank(Hz) = 0 -- a UNIQUE pure stabiliser
      ground state (|S| = 2^8). CSS validity: RM(1,3) self-dual => all stabs commute.
  (2) its EXACT entanglement entropy via the Fattal et al. (2004) stabiliser formula
      S_A = |A| - dim_F2(Cx cap V_A) - dim_F2(Cz cap V_A), with Cx=Cz=RM(1,3),
      so S_A = |A| - 2*dim(RM(1,3) supported entirely in A).
  (3) the Kitaev-Preskill tripartite combination
      gamma = S_A+S_B+S_C - S_AB - S_BC - S_CA + S_ABC  (finite-size, on 8 qubits).

CONCLUSION (printed): k=0 => trivial gapped stabiliser state => NO intrinsic
topological order at the CELL level (TO needs ground-space degeneracy growing with
system topology). So gamma_TEE != 0 (the TQO foundation) must come from the
macroscopic inter-cell gauge-web lattice on Z^3 (x) Q3 -- and is GATED on an
explicit macroscopic stabiliser model, which canon describes only qualitatively
(L(TCH), gauge web), not as a diagonalisable operator set.

numpy only.
"""
import numpy as np
from itertools import combinations
from collections import Counter

N = 8
# RM(1,3) generator matrix == the framework's R-constraint supports (same GEN as the
# committed scripts smg_code_projection.py, su4_smg_real_test.py, ps_su4_fock_smg.py).
GEN = np.array([
    [1, 1, 1, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 1, 1, 1, 1],
    [0, 0, 1, 1, 0, 0, 1, 1],
    [0, 1, 0, 1, 0, 1, 0, 1],
], dtype=int) % 2


def rank_f2(M):
    M = (M.copy() % 2).astype(int)
    r = 0
    rows, cols = M.shape
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i, c]), None)
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        r += 1
        if r == rows:
            break
    return r


def codewords(G):
    k = G.shape[0]
    out = []
    for m in range(1 << k):
        v = np.zeros(N, dtype=int)
        for b in range(k):
            if (m >> b) & 1:
                v ^= G[b]
        out.append(tuple(v % 2))
    return set(out)


# --- CSS validity + code parameters ---------------------------------------
# commutation: every pair of X-support and Z-support must overlap in even parity.
overlaps = [(int(np.dot(GEN[i], GEN[j]) % 2)) for i in range(4) for j in range(4)]
assert max(overlaps) == 0, "stabilisers fail to commute (odd overlap)"
rx, rz = rank_f2(GEN), rank_f2(GEN)
k = N - rx - rz
print(f"[8,4,4] CSS cell: rank(Hx)={rx}, rank(Hz)={rz}, n={N}  ->  k = {k}")
print(f"  => {'UNIQUE pure stabiliser ground state' if k == 0 else f'{2**k}-fold degenerate ground space'}")
print(f"  |S| = 2^(rx+rz) = {2**(rx+rz)}   (a pure state on {N} qubits requires |S| = 2^{N} = {2**N})")
assert k == 0 and rx + rz == N, "expected [[8,0,4]]"

Cx = codewords(GEN)
Cz = codewords(GEN)
we = Counter(sum(c) for c in Cx)
print(f"  Cx = Cz = RM(1,3): {len(Cx)} codewords; weight enumerator = {dict(sorted(we.items()))}")


# --- exact entanglement entropy -------------------------------------------
def dim_supported_in(C, A):
    Aset = set(A)
    sub = [c for c in C if all(c[i] == 0 for i in range(N) if i not in Aset)]
    return int(round(np.log2(len(sub)))) if sub else 0


def S(A):
    return len(A) - dim_supported_in(Cx, A) - dim_supported_in(Cz, A)


print("\nExact entanglement entropy S_A by region size (bits):")
for size in range(0, N + 1):
    vals = [S(list(A)) for A in combinations(range(N), size)]
    print(f"  |A|={size}: S_A in [{min(vals)}, {max(vals)}],  mean {np.mean(vals):.3f}")
print("  (S_A = |A| for |A|<=3: no codeword of distance 4 fits -> maximally entangled small")
print("   regions = good-code entanglement. Drops appear at |A|=4 supports of weight-4 codewords.)")


# --- Kitaev-Preskill tripartite combination -------------------------------
def KP(A, B, C):
    return (S(A) + S(B) + S(C)
            - S(A + B) - S(B + C) - S(C + A)
            + S(A + B + C))


print("\nKitaev-Preskill tripartite gamma (finite-size combination on 8 qubits --")
print("NOT the thermodynamic invariant; included only to show the cell carries no")
print("size-stable topological term):")
for A, B, C in [([0, 1], [2, 3], [4, 5]),
                ([0, 1, 2], [3, 4], [5, 6]),
                ([0], [1, 2, 3], [4, 5, 6]),
                ([0, 1], [2, 3, 4], [5, 6, 7])]:
    print(f"  A={A} B={B} C={C}:  gamma = {KP(A, B, C)}")

print("""
CONCLUSION
==========
 * The [8,4,4] CSS matter cell is a [[8,0,4]] code (k=0): a UNIQUE gapped stabiliser
   state. k=0 => NO topological ground-space degeneracy => a TRIVIAL gapped phase,
   i.e. NO intrinsic topological order AT THE CELL LEVEL. The high S_A for small A is
   good-code entanglement (distance 4), not long-range topological order.
 * Hence the TQO reframe's foundation -- gamma_TEE != 0, emergent gauge bosons/anyons --
   CANNOT originate in the cell. It must come from the macroscopic inter-cell
   gauge-web lattice on Z^3 (x) Q3.
 * That macroscopic stabiliser Hamiltonian is described only qualitatively in canon
   (line graph L(TCH); 'gauge web'; Wilson action on the SC dual). It is NOT yet
   specified as an explicit, diagonalisable operator set; the existing strip scripts
   (css_2d_strip_chargeblock_krylov_audit.py) are charge-block Krylov audits, not
   TEE-ready stabiliser lattices.
 * THEREFORE: a genuine macroscopic gamma_TEE is GATED on first writing down the
   explicit macroscopic stabiliser model. Until that exists, 'the TCH vacuum is
   topologically ordered' is NOT earned by computation. Target 3 is well-posed and
   exact, but blocked at step 2 (the lattice model), not step 1 (the cell).
""")
