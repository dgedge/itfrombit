#!/usr/bin/env python3
"""
Chiral-lattice-gauge frontier: checks for the Lüscher-overlap / Pati-Salam synthesis
(DRIFT K4 consolidation). Two operator-exact results, asserted:

  (1) ADMISSIBILITY: the walk Wilson kernel D_W = V.alpha + W_body*beta VANISHES at all 4
      EVEN Brillouin-zone corners (V=0 and W_body=0 there) -> sgn(gamma5 D_W) is undefined
      -> the walk-sourced overlap is NON-ADMISSIBLE; the 4 even-corner modes are unlifted
      momentum-space TASTES (distinct BZ momenta), NOT an internal SU(4) gauge-colour 4-plet.
      A single chiral species requires an IMPORTED (face-direction) Wilson term.

  (2) COMPLEXITY: the chiral fermion content is COMPLEX (R != R-bar). The one-generation
      hypercharge multiset {Y} is NOT symmetric under Y -> -Y -> complex fermion determinant
      -> Lüscher's non-Abelian phase obstruction is PRESENT (Pf^4 = det^2 would need a
      pseudo-real *full* rep; the joint (4,2,1)+(4bar,1,2) is complex).

Net: the frontier does NOT close; bare-substrate anomaly cancellations stand; continuum closure
is conditional on the unsolved non-Abelian chiral measure + imported exp-local overlap. numpy only.
"""
import numpy as np, itertools
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex); I2 = np.eye(2, dtype=complex)
kron = np.kron
al = [kron(sx, sx), kron(sx, sy), kron(sx, sz)]; beta = kron(sz, I2)


def V(k):
    kx, ky, kz = k
    return 8 * np.array([np.sin(kx) * np.cos(ky) * np.cos(kz), np.cos(kx) * np.sin(ky) * np.cos(kz), np.cos(kx) * np.cos(ky) * np.sin(kz)])
def Wbody(k):
    kx, ky, kz = k; return 8 * (1 - np.cos(kx) * np.cos(ky) * np.cos(kz))
def DW(k): return sum(V(k)[i] * al[i] for i in range(3)) + Wbody(k) * beta


print("(1) walk Wilson kernel ||D_W|| at the 8 BZ corners:")
even_norms, odd_norms = [], []
for c in itertools.product([0.0, np.pi], repeat=3):
    parity = int(round(sum(np.array(c) / np.pi))) % 2
    n = np.linalg.norm(DW(c)); (even_norms if parity == 0 else odd_norms).append(n)
    print(f"   k/pi={tuple(int(round(x/np.pi)) for x in c)} [{'EVEN' if parity==0 else 'odd '}]: ||D_W||={n:.3f}")
assert max(even_norms) < 1e-9 and min(odd_norms) > 1.0
print("   -> D_W=0 at all 4 EVEN corners -> sgn(0) undefined -> NON-ADMISSIBLE overlap; 4 unlifted")
print("      tastes (distinct BZ momenta), not a gauge colour 4-plet.")

print("\n(2) one-generation hypercharge multiset (Y=Q-T3) symmetry under Y->-Y:")
G0, G1, LQ, C0, C1, I3, CHI, W = range(8)
def Zf(c): return 1 if c[I3] == 0 else -1
def sZc(c): return -3 if (c[C0], c[C1]) == (0, 0) else -1
def Q(c): return 0.5 * Zf(c) + (1 / 3) * sZc(c) + 0.5
def T3(c): return (0.5 if c[I3] == 0 else -0.5) if c[CHI] == 0 else 0.0
Ys = []
for LQv, C0v, C1v, I3v in itertools.product([0, 1], repeat=4):
    if LQv == 0 and (C0v, C1v) == (0, 0) and I3v == 0: continue          # R4: no nu_R
    if LQv == 1 and (C0v, C1v) == (0, 0): continue                       # coloured quark
    if LQv == 0 and (C0v, C1v) != (0, 0): continue                       # colourless lepton
    Ys.append(round(Q([0, 0, LQv, C0v, C1v, I3v, 0, 0]) - T3([0, 0, LQv, C0v, C1v, I3v, 0, 0]), 4))
from collections import Counter
cnt = Counter(Ys); neg = Counter(round(-y, 4) for y in Ys)
print(f"   {{Y}}        = {sorted(cnt.items())}")
print(f"   {{-Y}}       = {sorted(neg.items())}")
print(f"   symmetric? {cnt == neg}  -> {'REAL/vector-like' if cnt==neg else 'COMPLEX/CHIRAL'}")
assert cnt != neg     # complex/chiral -> complex determinant -> Lüscher obstruction present
print("   -> COMPLEX content -> complex determinant -> non-Abelian Lüscher phase obstruction PRESENT.")
print("\nALL ASSERTS PASSED. Frontier does not close; see DRIFT K4 consolidation.")
