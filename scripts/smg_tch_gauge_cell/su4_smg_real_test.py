#!/usr/bin/env python3
"""
Approach B (Pati-Salam SU(4)_C) SMG mirror gap — GENUINE operator + local-ED test.

Replaces the retracted hardcoded `ps_su4` toy script (which did `gap = 2.0 if t>=2 else t`,
no Hamiltonian). Two real results, both asserted:

  (1) OBSTRUCTION-1 (colour-centre): the CSS X-stabiliser does NOT commute with the colour
      centre for SU(4) (Z4 cycle, Z2xZ2 generators) -- *identically* to SU(3) (Z3). Reason is
      Pauli algebra (X-stab is Pauli-X on the colour bits; the centre is Pauli-Z; X anticommutes
      with Z), group-independent. So SU(4) being "register-native" does NOT make the SMG
      gauge-covariant (= the gauge-dressing tautology again, check_gauge_dressing.py).

  (2) OBSTRUCTION-2 (kinetic collapse): the local SMG gap collapses under genuine fermion
      hopping (2.0 -> ~0.13 at t=4), reproducing the n/q=6 strip value, for SU(4) as for SU(3).

Conclusion: the SU(4) route does NOT close the chiral-lattice-gauge frontier. numpy only.
"""
import numpy as np
G0, G1, LQ, C0, C1, I3, CHI, W = range(8)
N = 8; DIM = 256
ALL = [tuple((s >> (N - 1 - b)) & 1 for b in range(N)) for s in range(DIM)]
GEN = [[1, 1, 1, 1, 1, 1, 1, 1], [0, 0, 0, 0, 1, 1, 1, 1], [0, 0, 1, 1, 0, 0, 1, 1], [0, 1, 0, 1, 0, 1, 0, 1]]


def idx(bits): return sum(bits[b] << (N - 1 - b) for b in range(N))
def jw(bit):
    M = np.zeros((DIM, DIM), complex)
    for s, b in enumerate(ALL):
        t = list(b); t[bit] ^= 1; M[idx(tuple(t)), s] = (-1) ** sum(b[:bit])
    return M
def par(bit): return np.diag([1 if b[bit] == 0 else -1 for b in ALL]).astype(complex)
def prod(ops):
    M = np.eye(DIM, dtype=complex)
    for o in ops: M = M @ o
    return M
def gap(H):
    ev = np.linalg.eigvalsh(H); lv = []
    for v in ev:
        if not lv or v - lv[-1] > 1e-7: lv.append(float(v))
    return lv[1] - lv[0]


g = [jw(b) for b in range(N)]; p = [par(b) for b in range(N)]
xs = [prod([g[b] for b in range(N) if row[b]]) for row in GEN]
zs = [prod([p[b] for b in range(N) if row[b]]) for row in GEN]

H0 = -sum(zs + xs)
print(f"GUARDRAIL local CSS gap = {gap(H0):.4f} (expect 2.0)")
assert abs(gap(H0) - 2.0) < 1e-6


def color_perm(cmap):
    P = np.zeros((DIM, DIM), complex)
    for s, b in enumerate(ALL):
        t = list(b); t[C0], t[C1] = cmap[(b[C0], b[C1])]; P[idx(tuple(t)), s] = 1
    return P


Z3 = color_perm({(0, 1): (1, 0), (1, 0): (1, 1), (1, 1): (0, 1), (0, 0): (0, 0)})   # SU(3) 3-cycle
Z4 = color_perm({(0, 0): (0, 1), (0, 1): (1, 0), (1, 0): (1, 1), (1, 1): (0, 0)})   # SU(4) 4-cycle
Zc0, Zc1 = par(C0), par(C1)


def maxcomm(ops, gen): return max(np.linalg.norm(o @ gen - gen @ o) for o in ops)


print("\n(1) X-stabiliser vs colour-centre commutators (0 would mean gauge-covariant):")
cZ3, cZ4, cc0, cc1 = maxcomm(xs, Z3), maxcomm(xs, Z4), maxcomm(xs, Zc0), maxcomm(xs, Zc1)
print(f"   SU(3) Z3 cycle : {cZ3:.3f}     SU(4) Z4 cycle : {cZ4:.3f}")
print(f"   Z2(C0) gen     : {cc0:.3f}    Z2(C1) gen     : {cc1:.3f}")
print("   -> nonzero for BOTH SU(3) and SU(4): X anticommutes with the Z-type centre (Pauli),")
print("      so SU(4) is NOT 'native'; the SMG is not gauge-covariant for either group.")
assert cZ3 > 1.0 and cZ4 > 1.0 and cc0 > 1.0 and cc1 > 1.0      # SU(4) fails just like SU(3)
assert abs(cZ3 - cZ4) < 1e-6                                    # identical failure

print("\n(2) SMG gap vs genuine kinetic hopping (Majorana bilinears):")
def hop(pairs):
    H = sum(1j * g[a] @ g[b] for a, b in pairs); return (H + H.conj().T) / 2
Hk = hop([(C0, C1), (I3, LQ), (CHI, G0), (W, G1)])
print(f"   ||[H_hop, X-stab]|| = {max(np.linalg.norm(Hk @ x - x @ Hk) for x in xs):.3f} (nonzero -> competes)")
print(f"   {'t':>5} {'gap':>8}")
gaps = {}
for t in [0.0, 1.0, 2.0, 4.0, 8.0]:
    gaps[t] = gap(H0 + t * Hk); print(f"   {t:>5.1f} {gaps[t]:>8.4f}")
assert gaps[0.0] > 1.99 and gaps[4.0] < 0.5        # collapses under hopping

print("\nRESULT: SU(4) fails obstruction-1 (operator-exact, = SU(3)) and obstruction-2 (gap")
print("collapses under hopping). The Pati-Salam route does NOT close the frontier.")
print("ALL ASSERTS PASSED.")
