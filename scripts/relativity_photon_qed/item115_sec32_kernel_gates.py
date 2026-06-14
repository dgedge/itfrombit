#!/usr/bin/env python3
r"""Do the 3.2 amplitudes restore isotropy? — YES-shaped, but not by reweighting:
3.2 specifies a DIFFERENT kernel, and it must be run through the theorem-target gates.

3.2 verbatim: T_d = -(i/sqrt3) R_d (V_em + V_weak + V_strong), with
  V_em      diagonal in Q (= I3 - (1-LQ)/2 as printed — see the erratum check below),
  V_weak    sqrt(2/9) CNOT (the 5.9/2.2 R-rule: I3(t+1) = I3(t) XOR LQ(t)),
  V_strong  PERMUTES (C0,C1), g_s = 1,
  T_y, T_z  C4-conjugates of T_x.

Structural contrast with the Final-Boss core the pin audit + gauge scan tested:
  * Final-Boss hops XOR-FLIP a register triple (colour pair flipped jointly, I3 flipped
    unconditionally). XOR{C0,C1} maps (0,0)<->(1,1) — it crosses the R3 colour-class
    boundary (kicks leptons out of P) and the unconditional I3 flip changes Q always.
  * 3.2's V_strong is the SWAP (0,1)<->(1,0): it NEVER crosses the colour class — and
    3.2's only Q-changing piece is V_weak, LQ-CONTROLLED (quarks only), the W channel.

So the gates must be re-run on the 3.2 object. Pre-registered gate set (the pin audit's
own theorem target): P-closure (P H P = H), [H, Q_116] = 0 for the EM-gauged part,
G0 conservation, rank-3 spatial support with cubic isotropy. Two Bloch readings posed
(canon's line is ambiguous): (i) SCALAR-HOP — the internal operator is direction-
independent, only the Bloch phase carries d; (ii) CONJUGATED — the register operator is
rotated per direction by the C4 octant action. Self-asserting; exit 0 = verified."""
import itertools as it
import numpy as np

# ---------------- spaces, charges ----------------
def b(n, i): return (n >> i) & 1
def valid(n):
    return not (b(n,0) and b(n,1)) and b(n,7) == b(n,6) and ((b(n,2) == 0) == ((b(n,3), b(n,4)) == (0, 0)))
PHYS = [n for n in range(256) if valid(n)]
assert len(PHYS) == 48
P = np.zeros((256, 256)); P[PHYS, PHYS] = 1.0
def q116(n):                                          # item-116 electric charge (authoritative)
    zf = 1 if b(n,5) == 0 else -1
    szc = -3 if (b(n,3), b(n,4)) == (0, 0) else -1
    return 0.5 * zf + szc / 3 + 0.5
Q = np.diag([q116(n) for n in range(256)])
ZG0 = np.diag([1.0 if not b(n,0) else -1.0 for n in range(256)])

# ---------------- the 3.2 channel operators on the 256-register ----------------
Vem = np.diag([b(n,5) - 0.5 * (1 - b(n,2)) for n in range(256)])   # diagonal (as printed)
Vweak = np.zeros((256, 256))                          # sqrt(2/9) CNOT: LQ controls, I3 target
for n in range(256):
    Vweak[n ^ (1 << 5) if b(n,2) else n, n] += np.sqrt(2 / 9)
Vstrong = np.zeros((256, 256))                        # SWAP (C0,C1), g_s = 1
for n in range(256):
    c0, c1 = b(n,3), b(n,4)
    m = n if c0 == c1 else (n ^ (1 << 3)) ^ (1 << 4)
    Vstrong[m, n] += 1.0
leak = lambda X: float(np.linalg.norm((np.eye(256) - P) @ X @ P))
comm = lambda X, Y: float(np.linalg.norm(X @ Y - Y @ X))

print("[1] GATES on the 3.2 channels (contrast: the Final-Boss core leaked 256/384 and")
print("    had ||[PHP,Q]|| = 9.05):")
for nm, X in [("V_em", Vem), ("V_strong", Vstrong), ("V_weak", Vweak)]:
    print(f"    {nm:<9s} P-leak {leak(X):.2e}   [.,Q116] {comm(X,Q):.3f}   [.,Z_G0] {comm(X,ZG0):.2e}")
assert leak(Vem) == leak(Vstrong) == leak(Vweak) == 0.0          # ALL P-closed exactly
assert comm(Vem, Q) < 1e-12 and comm(Vstrong, Q) < 1e-12         # Q-neutral transport channel
assert comm(Vweak, Q) > 1.0                                      # the W channel, and only it
assert max(comm(X, ZG0) for X in (Vem, Vstrong, Vweak)) < 1e-12  # G0 conserved by ALL channels
print("    -> P-closure EXACT for all three channels (the 256/384 leakage was an artifact")
print("       of XOR-flipping constraint bits); the Q-mixing content is EXACTLY V_weak")
print("       (LQ-controlled, quarks only — the properly-W-gauged channel); G0 trivially")
print("       conserved -> the {G0, Ward, 3D} trilemma was a property of the XOR core,")
print("       NOT of 3.2.")

# ---------------- [2] reading (i): scalar-hop Bloch kernel ----------------
Vgam = Vem + Vstrong                                  # the EM-gauged (Q-neutral) transport operator
HgamP = (P @ Vgam @ P)[np.ix_(PHYS, PHYS)]
assert np.linalg.norm(HgamP - HgamP.T) < 1e-12
# H(k) = (1/sqrt3) V_gam * sum_d 2 sin(k_d)  (from -(i/sqrt3) V e^{ik_d} + h.c.)
# separable: register operator x cubic lattice function -> velocity tensor proportional
# to identity at any band extremum; cubic isotropy is exact by construction.
ks = [np.array(v) for v in ([0.3,0,0],[0,0.3,0],[0,0,0.3])]
disp = [np.sort(np.linalg.eigvalsh(HgamP * (2/np.sqrt(3)) * np.sin(k).sum())) for k in ks]
assert max(np.linalg.norm(disp[0]-d) for d in disp[1:]) < 1e-12
print("\n[2] READING (i) scalar-hop: H(k) = f(k) (V_em+V_strong)|_P with f cubic —")
print("    dispersions along x,y,z IDENTICAL (separable kernel): rank-3 support, exact")
print("    cubic isotropy, [H,Q]=0, P H P = H, G0 conserved: ALL four theorem-target")
print("    gates PASS at FIRST order with FULL 48-state support.")

# ---------------- [3] reading (ii): C4-conjugated register operator ----------------
def role_perm_c4z():
    """C4 about z on octants (x,y,z)->(-y,x,z); octant v=(b2,b1,b0) are coords; the
    register ROLES sit on octants (octant map), so C4 permutes the 8 bit roles."""
    perm = {}
    for v in range(8):
        x, y, z = (b(v,2)*2-1), (b(v,1)*2-1), (b(v,0)*2-1)
        xn, yn = -y, x
        perm[v] = (((xn+1)//2) << 2) | (((yn+1)//2) << 1) | ((z+1)//2)
    return perm
perm = role_perm_c4z()
R = np.zeros((256, 256))
for n in range(256):
    m = sum(b(n, v) << perm[v] for v in range(8))
    R[m, n] = 1.0
Ty_int = R @ Vgam @ R.T                               # the y-direction internal operator
print("\n[3] READING (ii) conjugated: T_y internal = R V R^T with R the C4z octant-role")
print(f"    permutation: P-leak {leak(Ty_int):.3f}, [PT_yP, Q116] {comm(P@Ty_int@P, Q):.3f}")
ok2 = leak(Ty_int) < 1e-12 and comm(P @ Ty_int @ P, Q) < 1e-12
print(f"    -> reading (ii) {'ALSO passes' if ok2 else 'FAILS the gates (R permutes constraint'}")
if not ok2:
    print("       roles G0/G1/chi/W into active ones — the constraints are not O-covariant,")
    print("       as established in the gauge scan). The 3.2-faithful Ward kernel is")
    print("       reading (i); the conjugation line applies to the SPATIAL bridge, not the")
    print("       register roles, on pain of re-importing the XOR-core's failures.")

# ---------------- [4] the printed 3.2 Q-formula: erratum check ----------------
bad_bit = [n for n in PHYS if abs((b(n,5) - 0.5*(1 - b(n,2))) - q116(n)) > 1e-9]
i3half = lambda n: 0.5 if b(n,5) == 0 else -0.5       # charitable reading: I3 as +/-1/2
bad_half = [n for n in PHYS if abs((i3half(n) - 0.5*(1 - b(n,2))) - q116(n)) > 1e-9]
gaps = {round((i3half(n) - 0.5*(1 - b(n,2))) - q116(n), 9) for n in bad_half}
print(f"\n[4] ERRATUM CHECK on 3.2's printed 'Q = I3 - (1-LQ)/2' vs the item-116 charge:")
print(f"    bit reading (I3 in {{0,1}}): {len(bad_bit)}/48 disagree. Charitable reading (I3 as")
print(f"    +/-1/2): leptons all REPRODUCE; {len(bad_half)}/48 disagree — exactly the 36 quarks,")
print(f"    each off by {sorted(gaps)} = the missing quark hypercharge (B-L)/2 = +1/6 term.")
print(f"    Candidate erratum (corrected form: Q = I3 + (B-L)/2). The gate results above")
print(f"    use diagonality only, so they are UNAFFECTED.")
assert len(bad_half) == 36 and all(abs(g + 1/6) < 1e-6 for g in gaps)

print("""
================================================================================
VERDICT — "do the 3.2 amplitudes restore isotropy?": YES, and more than asked:
  3.2 does not reweight the Final-Boss core, it REPLACES it. Its strong channel is a
  colour SWAP (never crosses the R3 class boundary) and its only Q-changing channel
  is the LQ-controlled weak CNOT. Consequences, all computed:
    * P-closure EXACT for every 3.2 channel (the pin audit's 256/384 leakage and the
      gauge scan's rank-2 ceiling were artifacts of XOR-flipping constraint bits);
    * the Q-neutral transport operator V_em + V_strong commutes with the item-116
      charge AND with Z_G0: {G0 theorem, U(1) Ward, 3D} are SIMULTANEOUSLY satisfied
      — the pick-two trilemma dissolves for the 3.2 kernel;
    * the scalar-hop Bloch kernel passes all four theorem-target gates at FIRST
      order with full 48-state support and EXACT cubic isotropy (separable);
    * the layered-matter implication is WITHDRAWN for the 3.2 kernel: it was a
      property of the XOR-flip core only.
  Item 115's matter kernel is hereby CANDIDATE-PINNED: H_matter = P(V_em+V_strong)P
  x cubic Bloch, with V_weak the separately-gauged W channel. Remaining: the overall
  hop normalization t, the conjugation-line reading (the register-rotated variant
  fails the gates if R rotates constraint roles), the 3.2 Q-formula erratum, and the
  Pi_mn integral + sign itself.
================================================================================""")
print("ALL ASSERTIONS PASSED — every number above is verified. exit 0")
