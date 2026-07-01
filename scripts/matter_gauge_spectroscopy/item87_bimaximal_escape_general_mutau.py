#!/usr/bin/env python3
r"""ITEM 87/106 — chasing the "1-3-zero / general mu-tau" escape from the interlacing no-go
(item87_mutau_gim_bimaximal_attempt.py). Does a GENERAL mu-tau matrix (m23 != 0, i.e. WITHOUT the
topological-GIM texture zero) reach the bimaximal base, and is it substrate-motivated?

REFINEMENT of the no-go. The interlacing no-go was SPECIFIC to the GIM-zero case (m23 = e = 0): there
the mu-tau-odd eigenvalue c = d is a diagonal entry of the 2x2 block, so it is forced to the MIDDLE
mass. With e != 0 the odd eigenvalue is d - e (the 2x2-block diagonal is d + e) -- DECOUPLED from
interlacing -- so the odd state CAN be the heaviest (nu_3) -> theta_13 = 0. So the bimaximal base is
NOT a universal no-go; it is reachable by a general mu-tau matrix. We construct it and ask the two
honest questions: (1) is theta_12 then PREDICTED or FREE? (2) is the required m23 substrate-motivated?

Construction. mu-tau matrix m_nu = [[A, b, b],[b, d, e],[b, e, d]]. Odd state (0,1,-1)/sqrt2 has
eigenvalue d - e (set = m3, the heaviest -> theta_13 = 0, theta_23 = 45). The even 2x2 block in
basis {(1,0,0),(0,1,1)/sqrt2} is [[A, sqrt2 b],[sqrt2 b, d+e]] with eigenvalues {m1, m2} and a FREE
mixing angle theta_12. Parametrise the block by theta_12: [[A,F],[F,d+e]] = R(theta_12) diag(m1,m2) R^T,
then b = F/sqrt2, d = ((d+e)+m3)/2, e = ((d+e)-m3)/2.

Self-asserting; numpy only. exit 0 = the escape is mapped: base achievable, theta_12 FREE, m23 LARGE
(GIM-disfavoured); the no-go is confined to the GIM-zero sub-case.
"""
import sys
import numpy as np

_ok = True
def check(name, cond):
    global _ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    _ok = _ok and bool(cond)

def angles_deg(U):
    s13sq = abs(U[0, 2])**2
    s12sq = abs(U[0, 1])**2 / max(1 - s13sq, 1e-15)
    s23sq = abs(U[1, 2])**2 / max(1 - s13sq, 1e-15)
    return (np.degrees(np.arcsin(np.sqrt(np.clip(s12sq, 0, 1)))),
            np.degrees(np.arcsin(np.sqrt(np.clip(s23sq, 0, 1)))),
            np.degrees(np.arcsin(np.sqrt(np.clip(s13sq, 0, 1)))))

def build_mutau(theta12_deg, m1, m2, m3):
    """mu-tau matrix with odd eigenvalue = m3 (=> theta_13=0,theta_23=45) and 2x2 mixing = theta12."""
    t = np.radians(theta12_deg)
    R = np.array([[np.cos(t), np.sin(t)], [-np.sin(t), np.cos(t)]])
    B = R @ np.diag([m1, m2]) @ R.T          # 2x2 even block [[A,F],[F,d+e]]
    A, F, Dpe = B[0, 0], B[0, 1], B[1, 1]    # d+e = Dpe
    b = F / np.sqrt(2)
    d = (Dpe + m3) / 2.0
    e = (Dpe - m3) / 2.0                      # = m23
    return np.array([[A, b, b], [b, d, e], [b, e, d]]), e

def diag_angles(m):
    w, U = np.linalg.eigh(m)
    o = np.argsort(np.abs(w))
    return angles_deg(U[:, o]), np.abs(w[o])

m1, m2, m3 = 0.79e-3, 8.72e-3, 50.2e-3       # eV, Koide spectrum

# ---------------------------------------------------------------------------
# [1] the bimaximal base (theta_12 = 45) IS reachable by a general mu-tau matrix
# ---------------------------------------------------------------------------
print("[1] general mu-tau matrix targeting the bimaximal base (theta_12=45, theta_23=45, theta_13=0):")
m_bi, e_bi = build_mutau(45.0, m1, m2, m3)
(t12, t23, t13), w = diag_angles(m_bi)
print(f"    realised angles = ({t12:.2f}, {t23:.2f}, {t13:.2f}) deg; masses(meV)={np.round(w*1e3,2)}")
print(f"    required m23 = e = {e_bi*1e3:+.2f} meV   (entries all real: {np.all(np.isreal(m_bi))})")
check("a general mu-tau matrix REACHES the bimaximal base with the Koide masses (NOT a no-go)",
      abs(t12-45) < 0.5 and abs(t23-45) < 0.5 and t13 < 1e-6
      and np.allclose(np.sort(w), np.sort([m1, m2, m3])))

# ---------------------------------------------------------------------------
# [2] but theta_12 is a FREE parameter (1 DOF after the 3 eigenvalues): scan it
# ---------------------------------------------------------------------------
print("\n[2] theta_12 scan -- is it predicted or free? (all keep Koide masses, theta_13=0, theta_23=45)")
free = True
for tin in [10, 20, 30, 33.44, 45]:
    m, _ = build_mutau(tin, m1, m2, m3)
    (a12, a23, a13), ww = diag_angles(m)
    massok = np.allclose(np.sort(ww), np.sort([m1, m2, m3]))
    print(f"    theta_12 in = {tin:5.1f} -> out ({a12:.1f},{a23:.1f},{a13:.1f}); Koide masses kept: {massok}")
    free = free and abs(a12 - tin) < 0.5 and a13 < 1e-6 and abs(a23-45) < 0.5 and massok
check("theta_12 is a FREE parameter (any value realisable with the SAME masses) -> underdetermined", free)

# ---------------------------------------------------------------------------
# [3] the required m23 is LARGE -- exactly what topological GIM (DeltaW=2) SUPPRESSES
# ---------------------------------------------------------------------------
print("\n[3] is the required m23 substrate-motivated?")
print(f"    bimaximal needs |m23| = {abs(e_bi)*1e3:.1f} meV, comparable to the HEAVIEST mass "
      f"({m3*1e3:.0f} meV) -- a LARGE Gen2-Gen3 (DeltaW=2) coupling.")
print(f"    but topological GIM (the Hamming theorem, quark |V_13|=0) SUPPRESSES DeltaW=2 couplings")
print(f"    -> m23 -> 0, which is exactly the interlacing-FAIL case (item87_mutau_gim_..., theta_12->0).")
m_gim, _ = build_mutau(45.0, m1, m2, m3);
# force the GIM limit e->0 on the same even-block, show it collapses to the middle-mass/no-go
A = m_bi[0,0]; b = m_bi[0,1]; Dpe = m_bi[1,1] + m_bi[1,2]   # d+e
m_e0 = np.array([[A, b, b],[b, Dpe/2, 0],[b, 0, Dpe/2]])    # same even block but e=0 (GIM)
(g12, g23, g13), _ = diag_angles(m_e0)
print(f"    GIM limit e=0 -> angles ({g12:.1f},{g23:.1f},{g13:.1f}) -- theta_12 collapses, base lost.")
check("the bimaximal needs a LARGE m23 (>= m2) that GIM suppresses; the GIM limit e=0 loses the base",
      abs(e_bi) > m2 and g12 < 1.0)

# ---------------------------------------------------------------------------
print(f"""
--- VERDICT (1-3-zero / general mu-tau escape: MAPPED) ---
The escape EXISTS, refining the earlier no-go: a GENERAL mu-tau matrix (m23 != 0) DOES reach the
bimaximal base (theta_13=0, theta_23=45, theta_12=45) with the Koide masses and all-real entries --
so the bimaximal base is NOT a universal no-go; the interlacing no-go is confined to the GIM-zero
sub-case (m23=0).
BUT the escape is NOT a derivation, for two independent reasons:
  (1) theta_12 is a FREE parameter -- the 2x2 even block has 1 leftover DOF after the three masses,
      so ANY theta_12 (10..45 deg) is realisable with the same spectrum. theta_12=45 (bimaximal) is a
      CHOICE, not forced. (This is the achievability-script underdetermination, now at the matrix level.)
  (2) the required m23 is LARGE (|m23| ~ {abs(e_bi)*1e3:.0f} meV, ~ the heaviest mass) -- a big Gen2-Gen3
      (DeltaW=2) coupling, which the framework's OWN topological-GIM theorem SUPPRESSES (DeltaW=2 -> 0).
      The substrate doesn't merely fail to pick the bimaximal texture; its GIM structure pushes m23->0,
      i.e. toward the interlacing-FAIL case. The substrate actively DISFAVOURS the bimaximal base.
NET: the bimaximal twist theta_13=delta/sqrt2 stays MATCHED but its base is achievable-only-by-tuning
(theta_12 free) AND GIM-disfavoured (needs a large DeltaW=2 coupling GIM suppresses). The "1-3-zero
escape" does not yield a substrate derivation -- it converts the no-go into "achievable but
underdetermined + GIM-disfavoured", which is the honest final status of the PMNS-angle derivation.
TIER: escape mapped; bimaximal base achievable-not-derived (theta_12 free; large m23 GIM-suppressed).
""")
print("\n" + ("ALL CHECKS PASSED -- escape mapped: achievable, theta_12 free, m23 large (GIM-disfavoured)"
              if _ok else "SOME CHECKS FAILED"))
sys.exit(0 if _ok else 1)
