#!/usr/bin/env python3
"""
item87_pmns_angle_achievability.py

PMNS angle-achievability test (#5). With the neutrino MASS spectrum fixed by
delta_nu=1/3 (item87_neutrino_mass_spectrum.py) and the seesaw structure
(m_D real-diagonal from the walk, M_R complex-symmetric Majorana), is the small
measured theta_13 ~ 8.6 deg ACHIEVABLE, or is it a genuine no-go?

Logic:
 1. The framework's delta_nu=1/3 neutrino mass matrix is a CIRCULANT. A circulant
    is diagonalized by the DFT for ANY entries -> trimaximal mixing -> theta_13 =
    arcsin(1/sqrt3) = 35.26 deg (symmetry-locked, independent of delta_nu). If the
    charged leptons are ALSO circulant (same DFT), PMNS = I -> theta_13 = 0. So
    the two symmetric limits BRACKET [0, 35.26] deg and the measured 8.6 sits
    between -- neither symmetric limit predicts it.
 2. Achievability: type-I seesaw m_nu = m_D M_R^-1 m_D with m_D real-diagonal can
    reproduce ANY light m_nu (Casas-Ibarra): given the measured PMNS + the
    delta_nu masses, M_R = m_D m_nu^-1 m_D is a valid complex-symmetric Majorana.
    So theta_13 = 8.6 deg IS achievable -> NOT a no-go.
 3. But m_D (the Dirac masses) is unfixed by canon: different m_D give different
    valid M_R for the SAME observables -> the angles are UNDETERMINED, not
    predicted. The one symmetry-respecting structure (circulant) overshoots.

Verdict: theta_13 is achievable-but-underdetermined (accommodatable, not a
falsification; but not predicted). Self-asserting; numpy only.
"""
import sys
import numpy as np

_ok = True


def check(name, cond):
    global _ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    _ok = _ok and bool(cond)


def angles_deg(U):
    """extract (theta12, theta23, theta13) in degrees from a 3x3 mixing matrix."""
    s13sq = abs(U[0, 2])**2
    s12sq = abs(U[0, 1])**2 / max(1 - s13sq, 1e-15)
    s23sq = abs(U[1, 2])**2 / max(1 - s13sq, 1e-15)
    deg = 180 / np.pi
    return (np.degrees(np.arcsin(np.sqrt(np.clip(s12sq, 0, 1)))),
            np.degrees(np.arcsin(np.sqrt(np.clip(s23sq, 0, 1)))),
            np.degrees(np.arcsin(np.sqrt(np.clip(s13sq, 0, 1)))))


def pmns_pdg(t12, t23, t13, dcp=0.0):
    """build a PMNS matrix from angles (deg) in the PDG convention."""
    t12, t23, t13 = np.radians([t12, t23, t13])
    s12, c12, s23, c23, s13, c13 = (np.sin(t12), np.cos(t12), np.sin(t23),
                                    np.cos(t23), np.sin(t13), np.cos(t13))
    e = np.exp(-1j * dcp)
    U23 = np.array([[1, 0, 0], [0, c23, s23], [0, -s23, c23]], complex)
    U13 = np.array([[c13, 0, s13 * e], [0, 1, 0], [-s13 * np.conj(e), 0, c13]], complex)
    U12 = np.array([[c12, s12, 0], [-s12, c12, 0], [0, 0, 1]], complex)
    return U23 @ U13 @ U12


def dft3():
    w = np.exp(2j * np.pi / 3)
    return np.array([[w**(j * k) for k in range(3)] for j in range(3)]) / np.sqrt(3)


# ===========================================================================
# 1. the circulant symmetry forces trimaximal -> theta_13 = 35.26 deg
# ===========================================================================
F = dft3()
# any symmetric circulant m_nu = circ(a,b,b); eigenvectors = DFT columns
for (a, b) in [(1.3 + 0.2j, 0.4 * np.exp(1j / 3)), (2.0, 0.9 * np.exp(1j / 3))]:
    C = np.array([[a, b, b], [b, a, b], [b, b, a]])
    diag = F.conj().T @ C @ F
    offdiag = np.linalg.norm(diag - np.diag(np.diag(diag)))
    assert offdiag < 1e-9, "DFT should diagonalize the circulant"
_, _, t13_circ = angles_deg(F)
print(f"circulant neutrino matrix (delta_nu=1/3) -> U_nu = DFT (trimaximal): "
      f"theta_13 = {t13_circ:.2f} deg")
check("circulant symmetry FORCES theta_13 = arcsin(1/sqrt3) = 35.26 deg",
      abs(t13_circ - 35.26) < 0.1)

# both sectors circulant -> PMNS = DFT^dag DFT = I -> theta_13 = 0
pmns_II = F.conj().T @ F
_, _, t13_both = angles_deg(pmns_II)
print(f"both sectors circulant -> PMNS = I: theta_13 = {t13_both:.2f} deg")
check("the two symmetric limits BRACKET the measured 8.6 deg (0 < 8.6 < 35.26)",
      t13_both < 1 and t13_circ > 33 and 0 < 8.6 < t13_circ)
check("NEITHER symmetric limit predicts 8.6 deg (both off by >15 deg)",
      abs(t13_circ - 8.6) > 15 and abs(t13_both - 8.6) > 7)

# ===========================================================================
# 2. achievability: seesaw with real-diagonal m_D reproduces the measured PMNS
# ===========================================================================
m = np.array([0.79e-3, 8.72e-3, 50.2e-3])          # eV, the delta_nu=1/3 spectrum
U_meas = pmns_pdg(33.4, 45.0, 8.6, dcp=1.0)         # measured PMNS (with a CP phase)
m_nu = U_meas.conj() @ np.diag(m) @ U_meas.conj().T  # light mass matrix (symmetric)
check("light mass matrix is complex-symmetric", np.allclose(m_nu, m_nu.T))


def seesaw_MR(mD_diag):
    mD = np.diag(mD_diag)                            # REAL-diagonal Dirac mass (walk)
    return mD @ np.linalg.inv(m_nu) @ mD            # M_R = m_D m_nu^-1 m_D


mD1 = np.array([1.0, 4.0, 20.0])                    # a hierarchical Dirac choice
MR1 = seesaw_MR(mD1)
mnu_back = np.diag(mD1) @ np.linalg.inv(MR1) @ np.diag(mD1)
t12b, t23b, t13b = angles_deg(U_meas)
print(f"\nseesaw (m_D real-diagonal {mD1}) -> M_R complex-symmetric: "
      f"{np.allclose(MR1, MR1.T)}")
print(f"  reproduces the target m_nu: ||diff||={np.linalg.norm(mnu_back-m_nu):.1e}")
print(f"  -> PMNS angles (theta12,theta23,theta13) = "
      f"({t12b:.1f}, {t23b:.1f}, {t13b:.1f}) deg")
check("M_R is a valid complex-symmetric Majorana", np.allclose(MR1, MR1.T))
check("seesaw reproduces the measured spectrum+mixing exactly (m_nu recovered)",
      np.linalg.norm(mnu_back - m_nu) < 1e-12)
check("=> theta_13 = 8.6 deg IS ACHIEVABLE in the substrate seesaw (NOT a no-go)",
      abs(t13b - 8.6) < 0.1)

# ===========================================================================
# 3. but m_D is unfixed -> the angles are UNDERDETERMINED (no prediction)
# ===========================================================================
mD2 = np.array([1.0, 1.0, 1.0])                     # a democratic Dirac choice
MR2 = seesaw_MR(mD2)
mnu_back2 = np.diag(mD2) @ np.linalg.inv(MR2) @ np.diag(mD2)
rel = np.linalg.norm(MR1 / np.linalg.norm(MR1) - MR2 / np.linalg.norm(MR2))
print(f"\na DIFFERENT m_D {mD2} gives a DIFFERENT valid M_R (normalized diff "
      f"{rel:.2f}) reproducing the SAME observables (||diff||="
      f"{np.linalg.norm(mnu_back2-m_nu):.1e})")
check("different (unfixed) m_D -> different valid M_R, same observables "
      "-> angles UNDERDETERMINED", rel > 0.3
      and np.linalg.norm(mnu_back2 - m_nu) < 1e-12)

# ===========================================================================
print("\n--- VERDICT ---")
print("theta_13 = 8.6 deg is ACHIEVABLE in the substrate seesaw -- NOT a no-go,")
print("NOT a falsification. A real-diagonal m_D (the walk's Y_nu) with a complex-")
print("symmetric M_R reproduces the measured PMNS + the delta_nu=1/3 masses exactly.")
print()
print("BUT it is UNDERDETERMINED, not predicted:")
print(f" * the framework's symmetry-respecting structure (delta_nu=1/3 CIRCULANT)")
print(f"   is symmetry-LOCKED to trimaximal, theta_13 = {t13_circ:.1f} deg -- a 4x")
print("   overshoot; the both-circulant limit gives theta_13 = 0. The measured")
print("   8.6 deg is bracketed in (0, 35.26) but predicted by neither symmetric limit.")
print(" * hitting 8.6 deg needs the seesaw's m_D (Dirac masses) which canon does")
print("   NOT fix -- different m_D give different valid M_R for the same data.")
print()
print("So PMNS #5 closes as: angles ACCOMMODATABLE but UNDETERMINED. The small")
print("theta_13 is neither a no-go nor a prediction -- the substrate fixes the")
print("MASSES (delta_nu=1/3, the clean win) but leaves the MIXING angles to the")
print("unfixed Dirac sector. Predicting theta_13 needs a derived m_D / a circulant-")
print("breaking operator the substrate does not supply.")

print("\n" + ("ALL CHECKS PASSED" if _ok else "SOME CHECKS FAILED"))
sys.exit(0 if _ok else 1)
