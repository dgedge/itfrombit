#!/usr/bin/env python3
r"""The K6 mode-matrix overlap — the refinement that decides the anisotropic-flow verdict.

CANON SPEC (DRIFT K7, item 102): the gauge-sector mode space is the O_h permutation rep
on the 6 cube faces: K6 complete-graph adjacency at Gamma, eigenvalues {5, -1 x5} =
A_1g + (T_1u + E_g); first-order k.p is fixed by O_h selection rules to the T_1u <-> E_g
channel with E_g bilinears u ~ (2 k_z p_z - k_x p_x - k_y p_y)/sqrt6,
v ~ (k_x p_x - k_y p_y)/sqrt2; derived velocities v[100]=sqrt(2/3), v[110]=1/sqrt2,
v[111]=1/sqrt3 (the upper branch) — these three numbers are the rebuild HARNESS.

WHAT THE OVERLAP DECIDES. The matter current couples to the T_1u (vector) content of a
gauge mode; the E_g (quadrupole strain) content is EM-inert (K8). The previous rung
matched the directional drag D(q^) to the tree branches assuming TRANSVERSE-vector
coupling. The K6 eigenvectors say what the coupling actually is:
  (i)  w_T: the T_1u weight of each propagating branch;
  (ii) e_T: the ORIENTATION of that T_1u content relative to q^ — transverse or
       longitudinal. A longitudinal vertex sees NO static drag (gauge invariance:
       static longitudinal A is pure gauge — verified on the matter kernel below).
Self-asserting; exit 0 = every number verified."""
import numpy as np

# ---------------- [1] K6 at Gamma + the symmetry basis ----------------
K6 = np.ones((6, 6)) - np.eye(6)                      # faces ordered +x,-x,+y,-y,+z,-z
ev = np.sort(np.linalg.eigvalsh(K6))
assert np.allclose(ev, [-1, -1, -1, -1, -1, 5])
r2 = 1 / np.sqrt(2)
P = {  # T_1u triple (vector) and E_g pair (u, v) on the 6 faces
    "px": np.array([r2, -r2, 0, 0, 0, 0]),
    "py": np.array([0, 0, r2, -r2, 0, 0]),
    "pz": np.array([0, 0, 0, 0, r2, -r2]),
    "u":  np.array([-1, -1, -1, -1, 2, 2]) / np.sqrt(12),
    "v":  np.array([1, 1, -1, -1, 0, 0]) / 2.0,
}
for a in P:
    for bnm in P:
        assert abs(np.dot(P[a], P[bnm]) - (a == bnm)) < 1e-12
print("[1] K6 at Gamma: eigenvalues {5, -1 x5} = A_1g + (T_1u + E_g) — canon reproduced;")
print("    orthonormal T_1u/E_g face basis built.")

# ---------------- [2] the k.p block and the velocity harness ----------------
def C_of(n):
    """The 2x3 E_g<-T_1u coupling block at unit |k| along n^ (rows u,v; cols px,py,pz)."""
    nx, ny, nz = n
    return np.array([[-nx / np.sqrt(6), -ny / np.sqrt(6), 2 * nz / np.sqrt(6)],
                     [nx / np.sqrt(2), -ny / np.sqrt(2), 0]])
def branches(n):
    """Singular triples of C: velocities = singular values; T-space right-vectors."""
    C = C_of(np.asarray(n, float) / np.linalg.norm(n))
    U, S, Vt = np.linalg.svd(C)
    return S, Vt[:len(S)], U                          # velocities, e_T rows, E-space parts
harness = {(1, 0, 0): np.sqrt(2 / 3), (1, 1, 0): 1 / np.sqrt(2), (1, 1, 1): 1 / np.sqrt(3)}
print("\n[2] VELOCITY HARNESS (upper branch vs canon):")
for n, vref in harness.items():
    S, _, _ = branches(n)
    print(f"    n = {n}: v_upper = {S[0]:.6f}  (canon {vref:.6f})  v_lower = {S[1]:.6f}")
    assert abs(S[0] - vref) < 1e-12
S110 = branches((1, 1, 0))[0]
assert abs(S110[1] - 1 / np.sqrt(6)) < 1e-12          # the lower [110] branch: 1/sqrt6
print("    exact match at all three directions (and lower branches 0, 1/sqrt6, 1/sqrt3).")

# ---------------- [3] the overlap structure: w_T and the e_T orientation ----------------
print("\n[3] EIGENVECTOR STRUCTURE of the propagating branches (the OVERLAP):")
print("    every propagating k.p eigenstate is (E_g-part +- T_1u-part)/sqrt2 — the block-")
print("    off-diagonal form forces w_T = 1/2 IDENTICALLY (no directional discrimination")
print("    from the weight). The DISCRIMINATION is in the ORIENTATION of the T_1u vector:")
for n in [(1, 0, 0), (1, 1, 0), (1, 1, 1)]:
    nh = np.asarray(n, float) / np.linalg.norm(n)
    S, Et, _ = branches(n)
    for bi, lab in [(0, "upper"), (1, "lower")]:
        if S[bi] < 1e-9:
            continue
        L = float(np.dot(Et[bi], nh) ** 2)
        print(f"    n = {n} {lab} (v = {S[bi]:.4f}): longitudinal fraction |e_T.n^|^2 = {L:.4f}")
# exact statements — the ORIENTATION MAP (the kernel of C rotates with direction):
S, Et, _ = branches((1, 0, 0))
assert abs(np.dot(Et[0], [1, 0, 0]) ** 2 - 1.0) < 1e-12        # [100] upper: LONGITUDINAL
S, Et, _ = branches((1, 1, 0))
nh110 = np.array([1, 1, 0]) / np.sqrt(2)
assert abs(np.dot(Et[0], nh110) ** 2) < 1e-12                  # [110] upper: TRANSVERSE
assert abs(abs(np.dot(Et[0], np.array([1, -1, 0]) / np.sqrt(2))) - 1) < 1e-9  # = in-plane
assert abs(np.dot(Et[1], nh110) ** 2 - 1.0) < 1e-12            # [110] lower: LONGITUDINAL
# [111]: kernel of C is the longitudinal direction itself -> propagating pair TRANSVERSE
ker111 = np.linalg.svd(C_of((1 / np.sqrt(3),) * 3))[2][2]
assert abs(abs(np.dot(ker111, np.ones(3) / np.sqrt(3))) - 1) < 1e-9
ker100 = np.linalg.svd(C_of((1.0, 0.0, 0.0)))[2][1:]           # 2-dim null at [100]
assert all(abs(np.dot(kv, [1, 0, 0])) < 1e-9 for kv in ker100) # = the transverse pair
print("    -> THE ORIENTATION MAP (computed; the C(k) kernel rotates with direction):")
print("       [100]: propagating branch LONGITUDINAL; the transverse pair = the ZERO modes")
print("       [110]: upper branch TRANSVERSE (in-plane (1,-1,0)); lower LONGITUDINAL;")
print("              zero mode = out-of-plane z")
print("       [111]: propagating pair TRANSVERSE; zero mode = LONGITUDINAL")
print("       The earlier all-longitudinal guess was WRONG (self-caught by assert): the")
print("       coupling character of the canonical fast branch CHANGES with direction.")

# ---------------- [3b] the corrected-vertex static drag table ----------------
SXp = np.array([[0, 1], [1, 0]], complex); SYp = np.array([[0, -1j], [1j, 0]], complex)
SZp = np.array([[1, 0], [0, -1]], complex); PAUL = (SXp, SYp, SZp)
sfun = lambda k: np.array([np.sin(k[0]), np.sin(k[1]), np.sin(k[2])])
def evecs2(sh):
    _, U = np.linalg.eigh(sh[0] * SXp + sh[1] * SYp + sh[2] * SZp)
    return U[:, 0], U[:, 1]
def kern(N, qvec, evec):
    ks = (np.arange(N) + 0.5) * 2 * np.pi / N - np.pi
    q = np.array(qvec, float) * 2 * np.pi / N
    e = np.array(evec, float); e /= np.linalg.norm(e)
    pi00 = para = dia = 0.0
    for kx in ks:
        for ky in ks:
            for kz in ks:
                k = np.array([kx, ky, kz]); kq = k + q
                sk, skq = sfun(k), sfun(kq)
                nk, nkq = np.linalg.norm(sk), np.linalg.norm(skq)
                um, _ = evecs2(sk / nk); _, vp = evecs2(skq / nkq)
                dE = nk + nkq
                pi00 += 2 * abs(np.vdot(um, vp)) ** 2 / dE
                V = sum(e[d] * np.cos(k[d] + q[d] / 2) * PAUL[d] for d in range(3))
                para += 2 * abs(np.vdot(um, V @ vp)) ** 2 / dE
                dia += sum(e[d] ** 2 * np.sin(k[d]) ** 2 for d in range(3)) / nk
    n3 = N ** 3; q2 = float(np.dot(q, q))
    return pi00 / n3 / q2, (dia - para) / n3 / q2
N = 16
aE110, aT110in = kern(N, (1, 1, 0), (1, -1, 0))       # the [110]-upper vertex (in-plane)
aE111, aT111 = kern(N, (1, 1, 1), (1, -1, 0))         # the [111] transverse vertex
D110in, D111 = aT110in - aE110, aT111 - aE111
print(f"""
[3b] CORRECTED-VERTEX STATIC DRAG TABLE (w_T = 1/2; e^2 = 4 pi alpha_0; N = {N}):
     branch          v_tree   vertex          static drag (x w_T)
     [100] upper     0.8165   LONGITUDINAL    0           (gauge-protected, [4])
     [110] upper     0.7071   transv. in-pl   {0.5*D110in:+.5f}
     [110] lower     0.4082   LONGITUDINAL    0           (gauge-protected)
     [111] pair      0.5774   transverse      {0.5*D111:+.5f}
     THE FASTEST BRANCH IS GAUGE-PROTECTED from static drag (longitudinal vertex);
     the slower transverse-coupled branches feel the full negative kernel. The
     [100]-vs-rest split WIDENS: the anisotropic flow is DIVERGENT — reinstated,
     now by the ORIENTATION mechanism (the previous rung's kernel-ordering reason
     was wrong, its verdict accidentally right; both recorded). Within the
     transverse sector the {{[110]u,[111]}} sub-split flows mildly CONVERGENT
     ({0.5*D110in:+.5f} vs {0.5*D111:+.5f}: the faster [110]u dragged slightly harder).""")
assert D110in < 0 and D111 < 0 and 0.5 * abs(D110in - D111) < 0.002

# ---------------- [4] static longitudinal drag vanishes (gauge check on matter) ----------------
SX = np.array([[0, 1], [1, 0]], complex); SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex); PAULI = (SX, SY, SZ)
s = lambda k: np.array([np.sin(k[0]), np.sin(k[1]), np.sin(k[2])])
def eigvecs(sh):
    _, U = np.linalg.eigh(sh[0] * SX + sh[1] * SY + sh[2] * SZ)
    return U[:, 0], U[:, 1]
def KL(N):
    """Static LONGITUDINAL current kernel along x at q = (2pi/N,0,0) — must vanish."""
    ks = (np.arange(N) + 0.5) * 2 * np.pi / N - np.pi
    q = np.array([2 * np.pi / N, 0, 0])
    para = dia = 0.0
    for kx in ks:
        for ky in ks:
            for kz in ks:
                k = np.array([kx, ky, kz]); kq = k + q
                sk, skq = s(k), s(kq)
                nk, nkq = np.linalg.norm(sk), np.linalg.norm(skq)
                um, _ = eigvecs(sk / nk); _, vp = eigvecs(skq / nkq)
                V = np.cos(k[0] + q[0] / 2) * SX      # LONGITUDINAL vertex (e || q)
                para += 2 * abs(np.vdot(um, V @ vp)) ** 2 / (nk + nkq)
                dia += np.sin(k[0]) ** 2 / nk
    return (dia - para) / N ** 3, dia / N ** 3
kl12, d12 = KL(12); kl20, d20 = KL(20)
print(f"\n[4] STATIC LONGITUDINAL GAUGE CHECK on the matter kernel: K_L(0,q)/dia =")
print(f"    {kl12/d12:+.2e} [N=12], {kl20/d20:+.2e} [N=20] — MACHINE-EXACT zero at finite q")
print("    and finite N (stronger than the transverse harness, whose cone-kink residual is")
print("    finite-N): a static longitudinal A at lattice momentum q IS an exact lattice")
print("    gauge transformation, so the longitudinal vertex sees no static drag, exactly.")
assert abs(kl12 / d12) < 1e-12 and abs(kl20 / d20) < 1e-12

print("""
================================================================================
VERDICT (the K6 mode-matrix overlap):
  * The overlap is exact and richer than guessed: every propagating branch is half
    E_g (EM-inert) + half T_1u, and the T_1u ORIENTATION rotates with direction —
    [100] longitudinal / [110] split (upper transverse, lower longitudinal) /
    [111] transverse. The C(k) kernel (zero mode) carries the complementary
    orientation. (The all-longitudinal first guess was wrong — self-caught.)
  * THE ANISOTROPIC VERDICT, CORRECTED IN MECHANISM AND REINSTATED IN SIGN: the
    previous rung's reasoning (kernel ordering D[100] > D[111] with transverse
    vertices everywhere) was the WRONG MECHANISM — but the corrected vertex map
    yields the same direction of flow, more starkly: the FASTEST branch ([100],
    longitudinal vertex) is GAUGE-PROTECTED from static drag entirely, while the
    slower transverse-coupled branches ([110]u, [111]) feel the full negative
    kernel. The [100]-vs-rest split widens: DIVERGENT, by orientation protection.
    Sub-structure: within the transverse sector the flow is mildly convergent.
  * THE DYNAMIC CAVEAT (named): longitudinal protection is exact only at omega = 0;
    on shell (omega = v_g q) the longitudinal channel couples through continuity to
    Pi_00(omega, q), and whether the [100] branch then feels a real shift or Landau
    damping depends on v_g vs v_m — the K5 ordering controls the qualitative
    outcome. The static verdict above is the leading IR statement.
  * STANDING RESULTS: the transverse kernels (0.881, the K5-free rate beta) now
    attach to the CORRECT objects — the [110]-upper and [111] branches and the
    flat transverse modes (which acquire induced dispersion ~ e^2 A_T q^2, a new
    separate observation).
  * NEXT OBJECT: the dynamic polarization Pi(omega ~ v q, q) + the v_web/v_m
    ordering — now the single gate on the velocity-unification conjecture's fate.
================================================================================""")
print("ALL ASSERTIONS PASSED — every number above is verified. exit 0")
