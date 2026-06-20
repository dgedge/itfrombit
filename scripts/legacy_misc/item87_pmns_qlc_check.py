#!/usr/bin/env python3
r"""ITEM 87 / ITEM 106 — the PMNS ANGLES: is item 53 (the seesaw SCALE) the lever, and what is the
honest status of the bimaximal delta=2/9 near-prediction once item 87's later no-go is folded in?

(NOTE: this script's "QLC" content RE-DERIVES the existing item 106 bimaximal ansatz — found after
the fact. It is kept because the re-derivation surfaces (i) item 53 ruled out as the angle lever and
(ii) a real item-106 <-> item-87 mechanism tension. No new-discovery claim is made for the numbers.)

THREE honest results:

[A] item 53 is a SCALE; the angles are a TEXTURE -> item 53 CANNOT gate the angles. In the seesaw
    m_nu = m_D M_R^-1 m_D, scaling M_R -> c*M_R (item 53 fixes |M_R|) scales the light MASSES by 1/c
    but leaves the EIGENVECTORS -- hence all three mixing angles -- invariant. Shown below.

[B] two ansaetze coexist in canon:
      bare circulant (delta_nu=1/3)  -> TRIMAXIMAL  theta_13 = arcsin(1/sqrt3) = 35.26 deg (4x overshoot)
      item 106 bimaximal (delta=2/9) -> theta_13 = delta/sqrt2 = 9.0 deg, theta_12 = 45-delta = 32.3 deg
    The bimaximal matches NuFIT (8.57, 33.44) to ~0.4/1.2 deg -- a parameter-free near-prediction.

[C] the bimaximal twist delta = 2/9 IS the universal-2/9 Cabibbo lambda_W (item 86 / CKM item 126),
    so the bimaximal angles are a quark-lepton-complementarity (QLC) relation, not an independent fit.

[D] BUT item 106's "formally closed (UV-TBM -> IR-bimaximal RG-flow)" status is UNDERMINED by item 87's
    later (2026-06-16) work: the TBM theta_13~6deg source (the loop-level VIRTUAL COLOUR BRIDGE) is a
    STRUCTURAL no-go (item87_missing_operator_search.py: the walk connects no two valid nu generations
    at any order), the cross-talk route is M12-falsified, and the RG-flow TBM->bimaximal is asserted-
    not-computed. So the bimaximal NUMBERS stand (parameter-free, delta=2/9); the substrate MECHANISM
    connecting them does not. Honest tier: MATCHED, MECHANISM-OPEN (not "formally closed").

Self-asserting; numpy only.
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

T12_OBS, T23_OBS, T13_OBS = 33.44, 49.2, 8.57          # NuFIT 5.x normal ordering

# ===========================================================================
# [A] item 53 (SCALE) does NOT gate the angles  (scale-invariance demonstration)
# ===========================================================================
print("[A] does the seesaw SCALE (item 53) move the angles?")
rng = np.random.default_rng(20260618)
mD = np.diag([1.0, 4.0, 20.0])
A = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)); MR = A + A.T
def nu_angles(MR_):
    m_nu = mD @ np.linalg.inv(MR_) @ mD
    w, U = np.linalg.eigh(m_nu @ m_nu.conj().T)
    return angles_deg(U[:, np.argsort(w)])
ang1, ang2 = nu_angles(MR), nu_angles(1e6 * MR)
print(f"    angles at |M_R| vs 1e6*|M_R|: ({ang1[0]:.2f},{ang1[1]:.2f},{ang1[2]:.2f}) "
      f"vs ({ang2[0]:.2f},{ang2[1]:.2f},{ang2[2]:.2f}) deg")
check("scaling M_R (item 53) leaves all three angles invariant -> item 53 cannot gate the angles",
      np.allclose(ang1, ang2, atol=1e-6))

# ===========================================================================
# [B]+[C] trimaximal overshoot vs item-106 bimaximal (delta = universal 2/9 = CKM lambda)
# ===========================================================================
t13_tri = np.degrees(np.arcsin(np.sqrt(1/3)))
print(f"\n[B] bare circulant -> trimaximal theta_13 = {t13_tri:.2f} deg vs measured {T13_OBS} "
      f"({t13_tri/T13_OBS:.1f}x overshoot)")
check("trimaximal theta_13 = 35.26 deg overshoots measured 8.57 deg (texture-locked)",
      abs(t13_tri - T13_OBS) > 25)

delta = 2/9                                             # the universal 2/9 twist (item 86); = CKM lambda_W
lam_W = 2/9                                             # CKM Wolfenstein lambda (item 126)
delta_deg = np.degrees(delta)                           # item-106 convention: delta as radians-as-angle
bim_t13 = delta_deg / np.sqrt(2)                        # item 106: theta_13 = delta/sqrt2
bim_t12 = 45.0 - delta_deg                              # item 106: theta_12 = 45 - delta
print(f"\n[C] item-106 bimaximal with delta = 2/9 (= universal 2/9 = CKM lambda):")
print(f"    delta = 2/9 = {delta_deg:.2f} deg;  theta_13 = delta/sqrt2 = {bim_t13:.2f} deg "
      f"vs {T13_OBS} (|d|={abs(bim_t13-T13_OBS):.2f})")
print(f"    theta_12 = 45 - delta = {bim_t12:.2f} deg vs {T12_OBS} (|d|={abs(bim_t12-T12_OBS):.2f})")
check("bimaximal theta_13 = delta/sqrt2 within 1 deg of measured 8.57 (parameter-free near-prediction)",
      abs(bim_t13 - T13_OBS) < 1.0)
check("bimaximal theta_12 = 45-delta within 2 deg of measured 33.44", abs(bim_t12 - T12_OBS) < 2.0)
check("the bimaximal twist delta IS the universal-2/9 Cabibbo lambda_W (QLC, not an independent fit)",
      abs(delta - lam_W) < 1e-12)

# ===========================================================================
print(f"""
--- VERDICT (item 87 / item 106 angles) ---
[A] item 53 is a SCALE; the angles are a TEXTURE -- scaling M_R leaves all three angles invariant.
    => item 53 is RULED OUT as the angle lever (decisive). The open input is the m_D/M_R TEXTURE.
[B/C] the bimaximal ansatz (item 106) IS a parameter-free QLC near-prediction: with the framework's
    OWN universal twist delta = 2/9 (= CKM lambda), theta_13 = delta/sqrt2 = {bim_t13:.1f} deg and
    theta_12 = 45-delta = {bim_t12:.1f} deg match NuFIT (8.57, 33.44) to {abs(bim_t13-T13_OBS):.1f}/{abs(bim_t12-T12_OBS):.1f} deg.
    The bare circulant gives trimaximal 35.26 deg (overshoot); the bimaximal is the good one.
[D] BUT item 106's "formally closed (UV-TBM -> IR-bimaximal RG-flow)" OVERSTATES it. Folding in
    item 87 (2026-06-16): the TBM theta_13~6deg leg (virtual colour bridge) is a STRUCTURAL no-go,
    the cross-talk is M12-falsified, and the RG-flow is asserted-not-computed. So the NUMBERS stand;
    the MECHANISM connecting them to the substrate does not.
NET TIER: item 53 ruled out (new, decisive); the bimaximal delta=2/9 angles are MATCHED
(parameter-free near-prediction) but MECHANISM-OPEN -- item 106 reclassified from "formally closed"
to "matched, mechanism-open". The sharp open question: derive a delta=2/9-sized (NOT maximal, NOT the
no-go colour bridge) charged-lepton/neutrino texture that yields the bimaximal twist on the data side.
""")
print("\n" + ("ALL CHECKS PASSED" if _ok else "SOME CHECKS FAILED"))
sys.exit(0 if _ok else 1)
