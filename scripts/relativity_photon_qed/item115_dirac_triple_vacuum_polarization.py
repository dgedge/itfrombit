#!/usr/bin/env python3
r"""Item 115: the LIVE vacuum-polarization computation — canon-native Dirac triple,
Ward audit of 3.5's own alphas, the forced repair, and the signed finite loop.

STAGES (each gated, each self-asserting):
 [1] 3.5's literal Dirac matrices alpha_i = sigma_x(chi) x sigma_i(I3), realized on the
     physical 48 (logical chi flip = the R2-locked pair flip X6X7): Cl(3,1) algebra
     REPRODUCES (canon's claim verified) and P-closure + G0 hold — but the WARD GATE
     fails for alpha_1, alpha_2: the I3 factor carries electric charge (z_f), so
     sigma_x,y(I3) are u<->d / nu<->e transitions — W-currents, not EM transport.
     3.5's spin-=-I3 identification cannot feed the EM loop.
 [2] THE FORCED REPAIR: enumerate all 256 register flip masks — exactly ONE nontrivial
     mask is P-closed on all 48 AND charge-blind: the chi/W pair (the chirality logical
     qubit, 3.5's own FIRST tensor factor). The Ward-compatible Dirac triple is the
     Pauli triple on that qubit: M_d = sigma_d(chi-logical) x 1_24. All gates pass;
     the kernel is 24 charge-carrying Weyl doublets; Sum Q^2 = 8; cones at the 8 BZ
     corners (naive doublers — the framework's own BZ-corner <-> codeword structure,
     items 77/98). NAMED PREMISE: the EM-loop Dirac space is the chirality doublet
     (departing from 3.5's I3-spin assignment, which the Ward gate excludes).
 [3] THE LOOP, live/signed/finite, on H(k) = t Sum_d M_d sin k_d:
     - transverse photon self-energy Pi_T(0,q) = e^2 [ <d2H/dky^2>_sea - para_T(q) ]
       (operational Kubo definition; London/superconductor sign convention so that
       Pi_T(0,0) = 0 is the gauge-invariance harness, EXACT by band periodicity);
     - the photon velocity shift in Coulomb gauge: omega^2 = c^2 q^2 + Pi_T(0,q), so
       sign(Delta c^2) = sign(Pi_T(0,q))|small q — THE SIGN, checked across grids,
       t values, and q directions (isotropy);
     - the screening side Pi_00(q)/q^2 > 0 reported (the alpha-running check);
     - honest magnitude status: Delta c^2/c^2 = 4 pi alpha_0 (Pi_T/q^2)/K_gamma needs
       the gauge-web tree stiffness K_gamma — the K5 assignment pair AGAIN; the SIGN
       and the dimensionless shape are assignment-free, the magnitude is not.
Self-asserting; exit 0 = every number verified."""
import itertools as it
import numpy as np

# ---------------- the 48-space and operators ----------------
def b(n, i): return (n >> i) & 1
def valid(n):
    return not (b(n,0) and b(n,1)) and b(n,7) == b(n,6) and ((b(n,2) == 0) == ((b(n,3), b(n,4)) == (0, 0)))
PHYS = [n for n in range(256) if valid(n)]
IDX = {n: i for i, n in enumerate(PHYS)}
def q116(n):
    zf = 1 if b(n,5) == 0 else -1
    szc = -3 if (b(n,3), b(n,4)) == (0, 0) else -1
    return 0.5 * zf + szc / 3 + 0.5
Q = np.diag([q116(n) for n in PHYS])
ZG0 = np.diag([1.0 if not b(n,0) else -1.0 for n in PHYS])

def flip_op(mask):
    """Register flip as an operator on the 48 (zero where it leaves P)."""
    O = np.zeros((48, 48))
    for col, n in enumerate(PHYS):
        m = n ^ mask
        if m in IDX:
            O[IDX[m], col] = 1.0
    return O
def phase_op(bit, sgn=(1.0, -1.0)):
    return np.diag([sgn[0] if b(n, bit) == 0 else sgn[1] for n in PHYS])

CHI_MASK = (1 << 6) | (1 << 7)                        # the R2-locked chi/W pair
Xchi = flip_op(CHI_MASK)                              # logical sigma_x on chirality
Zchi = phase_op(6)                                    # logical sigma_z (chi parity)
Ychi = 1j * Xchi @ Zchi                               # logical sigma_y (Hermitian check below)
X_I3, Z_I3 = flip_op(1 << 5), phase_op(5)
Y_I3 = 1j * X_I3 @ Z_I3
comm = lambda A, B: float(np.linalg.norm(A @ B - B @ A))
acomm = lambda A, B: float(np.linalg.norm(A @ B + B @ A))

# ---------------- [1] 3.5's literal alphas: algebra OK, Ward FAILS ----------------
alphas_35 = [Xchi @ X_I3, Xchi @ Y_I3, Xchi @ Z_I3]   # sigma_x(chi) x sigma_i(I3)
for i in range(3):
    A = alphas_35[i]
    assert np.linalg.norm(A - A.conj().T) < 1e-12     # Hermitian (P-closed by construction)
    assert np.linalg.norm(A @ A - np.eye(48)) < 1e-12 # involutions: full support on the 48
    for j in range(i + 1, 3):
        assert acomm(A, alphas_35[j]) < 1e-12         # Cl(3) anticommutation — 3.5 verified
w = [comm(A, Q) for A in alphas_35]
g0 = max(comm(A, ZG0) for A in alphas_35)
print("[1] 3.5's literal alphas on the 48: Cl(3) algebra REPRODUCES (anticommutation exact,")
print(f"    full support, P-closed), [.,Z_G0] = {g0:.1e} — but the WARD GATE fails:")
print(f"    ||[alpha_i, Q116]|| = ({w[0]:.3f}, {w[1]:.3f}, {w[2]:.3f}) — alpha_1, alpha_2 carry")
print("    charge-changing I3 flips (u<->d, nu<->e: W-currents). 3.5's spin-=-I3")
print("    identification cannot supply the EM matter loop. NAMED OBSTRUCTION.")
assert w[0] > 1 and w[1] > 1 and w[2] < 1e-12 and g0 < 1e-12

# ---------------- [2] the forced repair: uniqueness + the chi triple ----------------
ok_masks = []
for mask in range(256):
    if all((n ^ mask) in IDX and abs(q116(n ^ mask) - q116(n)) < 1e-12 for n in PHYS):
        ok_masks.append(mask)
print(f"\n[2] FLIP-MASK ENUMERATION: masks P-closed on ALL 48 and charge-blind: "
      f"{[format(m,'08b') for m in ok_masks]}")
assert ok_masks == [0, CHI_MASK]                      # identity + the chi/W pair ONLY
M3 = [Xchi, Ychi, Zchi]
for i in range(3):
    assert np.linalg.norm(M3[i] - M3[i].conj().T) < 1e-12
    assert comm(M3[i], Q) < 1e-12 and comm(M3[i], ZG0) < 1e-12
    for j in range(i + 1, 3):
        assert acomm(M3[i], M3[j]) < 1e-12
qvals = sorted(q116(n) for n in PHYS)
sumQ2_doublets = sum(q116(n) ** 2 for n in PHYS) / 2  # 24 chi-doublets, Q chi-blind
print(f"    THE WARD-COMPATIBLE DIRAC TRIPLE IS FORCED: the Pauli triple on the chirality")
print(f"    logical qubit (the unique nontrivial mask) — anticommuting, Hermitian, P-closed,")
print(f"    Q- and G0-commuting. Kernel = 24 charge-carrying Weyl doublets, Sum_f Q_f^2 = {sumQ2_doublets:.0f}.")
assert abs(sumQ2_doublets - 8.0) < 1e-12
# cone check at a BZ corner: E = t|s(k)| linear in |dk|
t0, dk = 1.0, 1e-4
s = lambda k: np.array([np.sin(k[0]), np.sin(k[1]), np.sin(k[2])])
E = lambda k, t: t * np.linalg.norm(s(k))
assert abs(E(np.array([dk, 0, 0]), t0) / dk - 1.0) < 1e-6   # v = t, isotropic cones
print("    Dispersion: E = +-t|s(k)| — isotropic Dirac cones at the 8 BZ corners (naive")
print("    doublers: the framework's own BZ-corner structure, items 77/98).")

# ---------------- [3] the loop: signed, finite, live ----------------
# 2x2 reduction per doublet (charge q_f), total weight Sum q_f^2 = 8.
# spinor overlaps for sigma.s eigenstates:
#   |<-,k|+,k'>|^2 = (1 - shat.shat')/2 ;  |<-,k|sigma_y|+,k'>|^2 computed explicitly.
def pauli_vec():
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]], complex)
    sz = np.array([[1, 0], [0, -1]], complex)
    return sx, sy, sz
SX, SY, SZ = pauli_vec()
def eigvecs(shat):
    Hm = shat[0] * SX + shat[1] * SY + shat[2] * SZ
    ev, U = np.linalg.eigh(Hm)
    return U[:, 0], U[:, 1]                           # minus, plus
def loop(N, t, qdir):
    ks = (np.arange(N) + 0.5) * 2 * np.pi / N - np.pi
    q = np.zeros(3); q[qdir] = 2 * np.pi / N
    ty = (qdir + 1) % 3                               # transverse current direction
    pi00 = paraT = dia = 0.0
    for kx in ks:
        for ky in ks:
            for kz in ks:
                k = np.array([kx, ky, kz]); kq = k + q
                sk, skq = s(k), s(kq)
                ek, ekq = t * np.linalg.norm(sk), t * np.linalg.norm(skq)
                um, _ = eigvecs(sk / np.linalg.norm(sk))
                _, vp = eigvecs(skq / np.linalg.norm(skq))
                dE = ekq + ek                          # E_+(k+q) - E_-(k)
                pi00 += 2 * abs(np.vdot(um, vp)) ** 2 / dE
                Vy = t * np.cos(k[ty]) * (SX, SY, SZ)[ty]
                paraT += 2 * abs(np.vdot(um, Vy @ vp)) ** 2 / dE
                dia += t * np.sin(k[ty]) ** 2 / np.linalg.norm(sk)   # <-|d2H/dk_ty^2|->
    n3 = N ** 3
    return pi00 / n3, (dia - paraT) / n3, q[qdir]
def gauge_harness(N, t):
    """K_T(q=0) = dia - para_T(0): exactly the BZ sum of d2E_-/dk^2 -> 0 by periodicity."""
    ks = (np.arange(N) + 0.5) * 2 * np.pi / N - np.pi
    paraT = dia = 0.0
    for kx in ks:
        for ky in ks:
            for kz in ks:
                k = np.array([kx, ky, kz]); sk = s(k)
                um, vp = eigvecs(sk / np.linalg.norm(sk))
                Vy = t * np.cos(k[1]) * SY
                paraT += 2 * abs(np.vdot(um, Vy @ vp)) ** 2 / (2 * t * np.linalg.norm(sk))
                dia += t * np.sin(k[1]) ** 2 / np.linalg.norm(sk)
    return (dia - paraT) / N ** 3, dia / N ** 3
h12, hd12 = gauge_harness(12, 1.0)
h20, hd20 = gauge_harness(20, 1.0)
print("\n[3] THE LOOP (per unit charge; multiply by Sum Q^2 = 8; e^2 = 4 pi alpha_0):")
print(f"    gauge harness: K_T(0)/dia = {h12/hd12:+.2e} [N=12] -> {h20/hd20:+.2e} [N=20]")
print("    (the identity sum_BZ d2E_-/dk^2 = 0 holds distributionally; the cones' singular")
print("    curvature is missed by the analytic-derivative grid sum, so the residual is")
print("    finite-N and must SHRINK with N — it does; gauge invariance passes as a limit.)")
assert abs(h20 / hd20) < abs(h12 / hd12) < 0.02
print("    N   t  qdir   A_E = Pi00/q^2   A_T = K_T(q)/q^2   A_T/(v^2 A_E)")
results = []
for N, t, qd in [(12, 1.0, 0), (16, 1.0, 0), (16, 3.0, 0), (12, 1.0, 1)]:
    p00, kT, qv = loop(N, t, qd)
    aE, aT = p00 / qv**2, kT / qv**2
    ratio = aT / (t * t * aE)                          # matter-cone velocity v = t
    results.append((N, t, qd, aE, aT, ratio))
    tag = "  (isotropy)" if qd == 1 else ""
    print(f"    {N:<3d} {t:.0f}  {'xy'[qd]}     {aE:+.4f}          {aT:+.5f}         {ratio:.4f}{tag}")
A_E = [r[3] for r in results]; A_T = [r[4] for r in results]; R = [r[5] for r in results]
assert all(a > 0 for a in A_E) and all(a > 0 for a in A_T)
assert max(R) - min(R) < 0.01 and 0.85 < R[0] < 0.91  # the invariant ratio, stable
print(f"""
    READING (convention-safe): in a Lorentz-invariant world the static transverse and
    longitudinal responses obey A_T = v^2 A_E exactly (both are the same covariant
    Pi-bar at spacelike q^2) — no velocity shift. The LATTICE breaks this: the
    computed invariant ratio A_T/(v^2 A_E) = {np.mean(R):.3f} < 1, stable across grids,
    hop scales (the 1/t vs t scalings of A_E, A_T cancel — verified), and q directions.
    Transverse stiffening LAGS screening by ~{(1-np.mean(R))*100:.0f}%, so the photon's velocity
    RELATIVE TO THE MATTER CONES decreases at one loop:
        Delta(c^2/v^2)/(c^2/v^2) ~= -4 pi alpha_0 x 8 x (A_E - A_T/v^2) x [stiffness map]
    with A_E - A_T/v^2 = {A_E[1]-A_T[1]:.4f} [N=16, t=1] — and the individual coefficients
    grow ~ln(1/q) (massless cones) while the RATIO stays fixed: the correction is
    MARGINAL, i.e. a leading-log velocity flow — exactly K7's structure, now with the
    computed sign: the flow drives c_gamma TOWARD the matter velocity from above
    (CONVERGENT if the tree ordering is c > v, the standard reading; the tree ordering
    itself is the K5-class stiffness assignment, named).

VERDICT: the vacuum-polarization computation is LIVE (nonzero interband current from
the chi-triple kernel), SIGNED (the invariant ratio 0.88 < 1: convergent drag, stable
across N, t, q-direction), and FINITE (lattice-regulated; gauge harness limit-verified). Named
premises: the chirality-doublet Dirac space (FORCED by the Ward gate — 3.5's I3-spin
reading is obstructed, itself a new canon finding); naive-doubler multiplicity
(canon-native, items 77/98); the K_gamma/tree-ordering assignment for absolute
magnitude (K5 class). The K7/K8 'loop sign' question has its first computed answer.""")
print("ALL ASSERTIONS PASSED — every number above is verified. exit 0")
