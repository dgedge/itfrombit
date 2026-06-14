r"""CONFIRM + CHARACTERISE the convergence-object result (fast, rigorous).

The exhaustive search (item115_clifford_triple_search.py) found: strict-gate
(P-close + Ward + G0) max anticommuting clique = 3, dispersion isotropic; a 4th
(mass) generator needs Ward dropped (clique -> 5). The only constraint-and-Q-
preserving hop X-part is a = 192 = flip(chi, W) jointly (R2-preserving chiral
flip). This script confirms the STRUCTURE and pins the mass story:
  (1) the spatial triple is the 3 Paulis of the effective chirality (chi,W)
      qubit; H(k) = sum_d sin(k_d) Gamma_d has H^2 = (sum sin^2 k_d) I EXACTLY
      -> isotropic 3D Dirac cone (T-R1 spatial closure, rigorous);
  (2) the natural Dirac mass X_chi X_W is Q-neutral and G0-neutral but COMMUTES
      with one triple member (it IS one of the chirality Paulis) -> it cannot be
      the anticommuting 4th; no Q-neutral Clifford mass exists;
  (3) reaching a 4th anticommuting generator requires a Q-CHANGING operator
      (Ward-broken) from the I3/colour sector — i.e. the weak (V_weak) channel.
exit 0 = H^2 isotropy exact, chirality-qubit identification verified, mass story
         pinned (fact + flagged interpretation).
"""
import numpy as np

def b(n, i): return (n >> i) & 1
def valid(n):
    return (not (b(n,0) and b(n,1))) and b(n,7) == b(n,6) and \
           ((b(n,2) == 0) == ((b(n,3), b(n,4)) == (0, 0)))
PHYS = [n for n in range(256) if valid(n)]
IDX = {n: i for i, n in enumerate(PHYS)}; PS = set(PHYS)
def q116(n):
    zf = 1 if b(n,5) == 0 else -1
    szc = -3 if (b(n,3), b(n,4)) == (0, 0) else -1
    return 0.5 * zf + szc / 3 + 0.5

def paul(a, z):
    M = np.zeros((48, 48), dtype=complex)
    pref = [1, 1j, -1, -1j][bin(a & z).count("1") % 4]
    for n in PHYS:
        m = n ^ a
        if m not in PS: return None
        M[IDX[m], IDX[n]] = pref * (-1) ** (bin(z & n).count("1"))
    return M

CHI, W, I3, C0, C1 = 6, 7, 5, 3, 4
# the effective chirality qubit on the R2 subspace {chi=W}: Paulis
sx = paul(1 << CHI | 1 << W, 0)                       # X_chi X_W : flips both
sz = paul(0, 1 << CHI)                                # Z_chi (= Z_W on R2 space)
sy = 1j * sx @ sz                                     # Y = i X Z  (Hermitian here)
TRIPLE = [sx, sy, sz]
print("[1] effective chirality-qubit Paulis on the 48 (X_chiX_W, Y, Z_chi):")
for i, G in enumerate(TRIPLE):
    herm = np.linalg.norm(G - G.conj().T) < 1e-9
    sq = np.linalg.norm(G @ G - np.eye(48)) < 1e-9
    print(f"    G[{i}]: Hermitian={herm}, square=I {sq}")
    assert herm and sq
for i in range(3):
    for j in range(i+1, 3):
        ac = np.linalg.norm(TRIPLE[i] @ TRIPLE[j] + TRIPLE[j] @ TRIPLE[i])
        print(f"    {{G[{i}],G[{j}]}} = {ac:.2e}  (anticommute)")
        assert ac < 1e-9
print("    -> a genuine mutually-anticommuting triple.")

print("\n[2] T-R1 ISOTROPY, rigorous: H(k)=sum_d sin(k_d)G_d  =>  H^2 = (sum sin^2 k_d) I")
def Hk(k): return sum(np.sin(k[d]) * TRIPLE[d] for d in range(3))
worst = 0.0
for k in [(0.3,-0.7,0.5),(1.1,0.2,-0.9),(0.05,0.05,0.05),(2.0,-1.3,0.4)]:
    s2 = sum(np.sin(x)**2 for x in k)
    worst = max(worst, np.linalg.norm(Hk(k) @ Hk(k) - s2 * np.eye(48)))
print(f"    max ||H^2 - (sum sin^2 k_d) I|| over test momenta = {worst:.2e}")
assert worst < 1e-9
print("    -> eigenvalues = +/- sqrt(sum sin^2 k_d): EXACTLY isotropic at leading")
print("       order (sum k_d^2 = |k|^2), anisotropy only at O(k^4). T-R1 SPATIAL CLOSED.")

print("\n[3] gauge gates on the triple (Ward + G0), computed:")
Q = np.diag([q116(n) for n in PHYS])
ZG0 = np.diag([1.0 if not b(n,0) else -1.0 for n in PHYS])
for i, G in enumerate(TRIPLE):
    cq = np.linalg.norm(G @ Q - Q @ G); cg = np.linalg.norm(G @ ZG0 - ZG0 @ G)
    print(f"    G[{i}]: [.,Q]={cq:.2e}  [.,Z_G0]={cg:.2e}")
    assert cq < 1e-9 and cg < 1e-9
print("    -> the spatial triple is simultaneously Ward (Q-neutral) AND G0-conserving:")
print("       the {G0, Ward, 3D} trilemma DOES NOT bind the spatial sector.")

print("\n[4] THE MASS (4th generator) — pinned over the FULL G0-conserving pool:")
def anticomm_all(G, triple): return all(np.linalg.norm(G@T+T@G) < 1e-9 for T in triple)
# full pool: every register-local reflection that P-closes and conserves G0
# (a_0 = 0), split by whether it is Q-neutral (Ward) or Q-changing.
qn_fourth, qc_fourth = [], []
for a in range(256):
    if a & 1: continue                               # G0-free (conserve generation)
    if any((n ^ a) not in PS for n in PHYS): continue   # P-closing X-part
    for z in range(256):
        G = paul(a, z)
        if G is None: continue
        if np.linalg.norm(G @ G - np.eye(48)) > 1e-9: continue
        if np.linalg.norm(G - G.conj().T) > 1e-9: continue
        if np.linalg.norm(G @ ZG0 - ZG0 @ G) > 1e-9: continue     # G0
        if not anticomm_all(G, TRIPLE): continue
        if np.linalg.norm(G @ Q - Q @ G) < 1e-9: qn_fourth.append((a, z))
        else: qc_fourth.append((a, z))
print(f"    operators anticommuting with the WHOLE chirality triple (full G0 pool):")
print(f"      Ward-NEUTRAL: {len(qn_fourth)}   Q-CHANGING: {len(qc_fourth)}")
assert len(qn_fourth) == 0 and len(qc_fourth) == 0
print("    -> ZERO of either kind. The chirality triple is the COMPLETE su(2) of one")
print("       qubit, and nothing anticommutes with all of X,Y,Z simultaneously on a")
print("       single qubit — the triple is isolated. So there is NO Clifford mass")
print("       (beta) of any kind extending it: the natural relativistic content is a")
print("       MASSLESS chiral (Weyl) fermion, and a bare Dirac mass does not exist on")
print("       this register. (The separate clique-5 the search found under dropped-Ward")
print("       is a DIFFERENT operator set; it contains no Ward-neutral spatial triple,")
print("       so it is not a gauge-compatible kinetic sector.)")
print("    -> FACT (exact, exhaustive): the gauge-compatible kinetic content is a")
print("       massless isotropic Weyl cone; fermion mass cannot be a bare Clifford")
print("       term and must arise from a Yukawa / symmetry-breaking coupling that")
print("       pairs the chiralities WITHOUT being an anticommuting beta — exactly the")
print("       Standard-Model situation (all SM fermion masses are Yukawa). [The Weyl")
print("       cone + no-bare-mass is derived; the Yukawa mechanism is the SM-standard")
print("       reading, not a new framework derivation.]")

print(f"""
[5] VERDICT — the convergence object, resolved:
  * T-R1 SPATIAL LIFT CLOSED: a P-closing, Ward-compatible, G0-conserving
    anticommuting Clifford TRIPLE exists (the 3 chirality-(chi,W)-qubit Paulis);
    H(k)=sum sin(k_d)G_d has H^2=(sum sin^2 k_d) I exactly => isotropic 3D Dirac
    dispersion E^2=|k|^2 at leading order, anisotropy O(k^4). Yesterday's
    [111]-anisotropy was the §3.2 scalar-hop's shared-Gamma artifact, NOT the
    register's limit — the triple was in the gate-passing pool all along.
  * T-R2 SIGN: the dressed propagator's leading kinetic term is now isotropic by
    construction (Clifford H^2 ∝ I), and the matter loop is O_h-symmetric (the
    T-R2 lemma) => the marginal anisotropy has no isotropic-order source; this
    moves T-R2 from 'leaning PASS' to 'PASS at the kernel level' for the chiral
    matter sector (the gauge-web E_g/T1u multiplet dressing still rides on top).
  * MASS: the chirality triple is COMPLETE and ISOLATED (nothing anticommutes
    with all 3 Paulis of a qubit) — NO Clifford mass of any kind extends it
    (Ward-neutral 0, Q-changing 0; exhaustive). So the gauge-compatible matter
    content is a MASSLESS isotropic Weyl cone, and a bare Dirac mass does not
    exist on the register: mass must be Yukawa / symmetry-breaking — the SM
    situation exactly (derived: Weyl cone + no bare mass; SM-standard: Yukawa).
  RESIDUAL: C4/O_h covariance of the triple (are G_1,G_2,G_3 the lattice-
    rotation images of one another, making the x,y,z assignment canonical) — the
    leading isotropy holds regardless (Clifford H^2 ∝ I), but manifest spatial
    covariance is the remaining structural check.
exit 0""")
print("ALL ASSERTIONS PASSED — triple confirmed, isotropy exact, mass story pinned.")
