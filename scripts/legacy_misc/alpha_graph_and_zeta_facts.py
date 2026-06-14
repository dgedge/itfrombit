#!/usr/bin/env python3
"""
Two independent verifications feeding the §5.4 / Part 12 alpha story.

PART 1 — GRAPH INVARIANTS of the 16-node 2-octagon interaction graph.
  Reconstructs the graph from the prose (two 8-cycles joined by two "rung"
  edges that, with one boundary edge from each octagon, close a 4-node bridge
  square), respecting the stated bipartition A={evens}, B={odds}. Checks:
  V=16, E=18, beta_1=3, bipartite, Kirchhoff spanning-tree count, the dual
  3x3 Laplacian det, the chiral +/-lambda spectral pairing, and the
  "connected 4-point trace" Tr[A^4]-[Tr A^2]^2.
  NOTE: the trace identity depends on the exact edge list, which the framework
  does not commit explicitly; this is ONE faithful reconstruction. Flagged.

PART 2 — WHERE zeta(3) LIVES.  Numerically confirms the Sommerfield-Petermann
  2-loop electron g-2 coefficient c2 = 197/144 + pi^2/12 - (pi^2/2)ln2
  + (3/4)zeta(3) equals the textbook -0.328478965..., and that DELETING the
  zeta(3) term breaks it. This is the VERTEX (a_e) coefficient. Contrasted with
  the QED running coefficients (vacuum polarization / beta function), which are
  RATIONAL through 3 loops (zeta(3) enters Adler at 4 loops, beta at 5 loops;
  arXiv:1005.2058, five-loop QED beta).

Self-asserting. exit 0 == all asserts pass. Needs numpy + mpmath.
"""
import numpy as np, mpmath as mp
mp.mp.dps = 30

print("="*70); print("PART 1 — 16-node 2-octagon graph invariants"); print("="*70)
# Build adjacency. Oct1 = 0..7 cycle; Oct2 = 8..15 cycle; rungs 1-8 and 9-0.
A = np.zeros((16,16), int)
def edge(i,j): A[i,j]=A[j,i]=1
for k in range(8):                      # octagon 1 : 0-1-...-7-0
    edge(k,(k+1)%8)
for k in range(8):                      # octagon 2 : 8-9-...-15-8
    edge(8+k, 8+(k+1)%8)
edge(1,8); edge(9,0)                    # two rungs closing the 4-cycle 0-1-8-9-0
V = 16; E = int(A.sum()//2)
beta1 = E - V + 1
print(f"  V={V}  E={E}  beta_1=E-V+1={beta1}")
# bipartite with stated A={evens},B={odds}?
bip = all(A[i,j]==0 for i in range(16) for j in range(16) if (i%2)==(j%2))
print(f"  bipartite under A=evens/B=odds: {bip}")
# bridge nodes of the 4-cycle 0-1-8-9
print(f"  bridge square cycle 0-1-8-9-0 present:",
      all([A[0,1],A[1,8],A[8,9],A[9,0]]), " (n_b=4)")
# Kirchhoff: spanning trees = any cofactor of the graph Laplacian
L = np.diag(A.sum(1)) - A
tau = round(float(np.linalg.det(L[1:,1:])))
print(f"  Kirchhoff spanning-tree count tau = {tau}")
# dual 3x3 Laplacian quoted in Part 12 §5 / ANCHOR §7.8
Ldual = np.array([[8,0,-1],[0,8,-1],[-1,-1,4]])
print(f"  dual 3x3 Laplacian det = {round(float(np.linalg.det(Ldual)))} (quoted 240)")
# chiral pairing: bipartite adjacency has symmetric +/- spectrum
ev = np.sort(np.linalg.eigvalsh(A.astype(float)))
sym = np.allclose(ev, -ev[::-1], atol=1e-9)
print(f"  adjacency spectrum symmetric (+/-lambda chiral pairing): {sym}")
# connected 4-point trace
trA2 = np.trace(np.linalg.matrix_power(A,2))
trA4 = np.trace(np.linalg.matrix_power(A,4))
conn = trA4 - trA2**2
print(f"  Tr[A^2]={trA2}  Tr[A^4]={trA4}  connected Tr[A^4]-[Tr A^2]^2 = {conn}")
print(f"  (framework claims this equals -240; reconstruction-dependent — see note)")

assert (V,E,beta1)==(16,18,3), "graph V/E/beta1 mismatch"
assert bip, "not bipartite under stated coloring"
assert round(float(np.linalg.det(Ldual)))==240, "dual Laplacian det != 240"
assert sym, "spectrum not chiral-symmetric"
# tau and conn are reported (reconstruction-dependent); not hard-asserted to 240.

print("\n  REPORTED (reconstruction-dependent): this faithful edge list gives")
print(f"  tau={tau}, connected-trace={conn}. If these differ from 240, the 240")
print("  claims need the framework's explicit edge list to adjudicate. The dual")
print("  Laplacian det=240 (from the quoted face structure) is solid regardless.")

print("\n" + "="*70); print("PART 2 — where zeta(3) lives in QED"); print("="*70)
pi, ln2, z3 = mp.pi, mp.log(2), mp.zeta(3)
c2_full = mp.mpf(197)/144 + pi**2/12 - (pi**2/2)*ln2 + mp.mpf(3)/4*z3
c2_noz3 = mp.mpf(197)/144 + pi**2/12 - (pi**2/2)*ln2            # drop the zeta(3) term
textbook = mp.mpf('-0.328478965579')   # Petermann 1957 / Sommerfield 1957, a_e 2-loop
print(f"  Sommerfield-Petermann c2 (with 3/4 zeta(3)) = {mp.nstr(c2_full,12)}")
print(f"  textbook value                              = {mp.nstr(textbook,12)}")
print(f"  match: {mp.almosteq(c2_full, textbook, abs_eps=mp.mpf('1e-9'))}")
print(f"  same coefficient WITHOUT the zeta(3) term   = {mp.nstr(c2_noz3,12)}")
print(f"  -> dropping zeta(3) shifts c2 by {mp.nstr(abs(c2_full-c2_noz3),4)} (=(3/4)zeta3): essential.")
print( "  c2 is the electron VERTEX F2(0) (a_e). The QED RUNNING (vacuum")
print( "  polarization/Adler/beta) is RATIONAL through 3 loops; zeta(3) enters")
print( "  Adler at 4 loops, beta at 5 loops (arXiv:1005.2058). So restoring an")
print( "  integral in the alpha^-1 (vacuum-pol) channel cannot and need not")
print( "  produce zeta(3): wrong Green's function (ANCHOR item 138 obstruction a).")

assert mp.almosteq(c2_full, textbook, abs_eps=mp.mpf('1e-9')), "c2 != textbook"
assert abs(c2_full - c2_noz3) > mp.mpf('0.9'), "zeta(3) term not essential"
print("\nexit 0 == graph invariants (V,E,beta1,bipartite,dual-det,chiral) and the")
print("zeta(3)-in-vertex facts all verified.")
