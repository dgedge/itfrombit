#!/usr/bin/env python3
r"""foundations_pmns_theta13_from_K3_portal.py

PMNS theta_13 "progress" attempt — does the canon's OWN M_R portal remove the underdetermination?

State (DRIFT l.2325-2330, l.2148-2150):
 - theta_23 ~= 45.9 deg is already a SHARP near-maximal prediction (the atmospheric latch -> second-order
   in delta; both first-order escapes closed). So the open angle is theta_13.
 - theta_13 = 8.6 deg was tiered "achievable but UNDERDETERMINED" because the type-I seesaw was inverted via
   Casas-Ibarra with M_R as the FREE output (m_D unfixed -> many valid M_R).
 - BUT the later entry l.2148 (2026-06-25) FIXES M_R to the generation-blind K_3 portal (the DeltaL=2 CP
   carrier): M_R proportional to A_{K_3}, "generation-blind: every unordered sterile pair has the same
   amplitude". The ANGLES are real (real Y_nu; CP lives only in M_R's phase, l.2150), so for the angles use
   the minimal real generation-blind form.

THE TEST (forward, NOT Casas-Ibarra): fix M_R to the canon portal form and m_D diagonal (walk-derived), then
the light seesaw m_nu = m_D M_R^{-1} m_D^T is DETERMINED by the m_D ratios alone -- the Casas-Ibarra freedom
is gone. Match the m_D ratios to the OBSERVED neutrino mass spectrum (normal ordering), then the angles are
OUTPUTS. Does theta_13 come out ~8.6 deg?

  M_R = a*I + b*A_{K_3} (generation-blind);  A_{K_3} = J - I.  m_D = diag(d1,d2,d3).
  m_nu = m_D M_R^{-1} m_D ; U_L ~= I (the observed PMNS is near-TBM -> charged leptons near-diagonal), so
  U_PMNS = eigenvectors(m_nu).

Two sub-cases for the portal:
  (P0) PURE off-diagonal portal a=0 (the ν_R carry NO bare diagonal mass -- they are forbidden R4
       pseudocodewords, only the DeltaL=2 portal connects them): ONE structural form, ZERO free ratio.
  (Pr) general generation-blind b/a free: the residual freedom, to see if a=0 is what pins theta_13.

Self-asserting on the STRUCTURE + the verdict (which sub-case predicts, which stays free); the angle value is
REPORTED, honestly. exit 0. Run under ~/bin/py13_7 (numpy).
"""
import numpy as np


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


J = np.ones((3, 3)); I3 = np.eye(3); AK3 = J - I3


def m_nu(d, a, b):
    MR = a * I3 + b * AK3
    D = np.diag(d)
    return D @ np.linalg.inv(MR) @ D            # real symmetric


def pmns_angles(M):
    w, V = np.linalg.eigh(M)                    # ascending eigenvalues
    order = np.argsort(np.abs(w))              # order by |mass|: 1,2,3 (NO)
    w = w[order]; V = V[:, order]
    U = np.abs(V)
    s13 = U[0, 2]
    th13 = np.degrees(np.arcsin(min(s13, 1.0)))
    th12 = np.degrees(np.arctan2(U[0, 1], U[0, 0]))
    th23 = np.degrees(np.arctan2(U[1, 2], U[2, 2]))
    return (th12, th13, th23), np.abs(w)


def mass_ratio(masses):
    m = np.sort(masses)                         # m1<m2<m3
    dm21 = m[1] ** 2 - m[0] ** 2; dm31 = m[2] ** 2 - m[0] ** 2
    return dm21 / dm31 if dm31 > 0 else np.inf


def main():
    R_OBS = 7.42e-5 / 2.515e-3                   # Delta m^2_21 / Delta m^2_31 (NO), ~0.0295
    OBS = (33.4, 8.6, 45.9)                      # measured (theta12, theta13, theta23); theta23 = framework latch
    print(f"=== target: Delta m^2 ratio r_obs = {R_OBS:.4f} (NO); observed angles {OBS} ===")

    print("\n[structure] symmetric limit (degenerate -> angle ill-defined, reported not asserted):")
    ang_dem, mass_dem = pmns_angles(m_nu(np.array([1.0, 1, 1]), 0.0, 1.0))
    print(f"    democratic m_D=(1,1,1): masses={np.round(mass_dem,3)} -> doublet DEGENERATE (eigenbasis arbitrary);")
    print(f"    the (1,1,1) singlet is one mass-eigenvector (trimaximal column); angles ill-defined at exact degeneracy.")

    # ---- (P0) PURE off-diagonal portal a=0: scan d-ratios, match the mass spectrum, read theta13 ----
    print("\n[P0] PURE off-diagonal portal (a=0, b=1): m_nu determined by d-ratios; match r_obs -> predict angles")
    best = None
    grid = np.concatenate([np.linspace(0.02, 1.0, 240), np.linspace(1.0, 5.0, 160)])
    sols = []
    for d1 in grid:
        for d2 in grid:
            d = np.array([d1, d2, 1.0])
            ang, masses = pmns_angles(m_nu(d, 0.0, 1.0))
            r = mass_ratio(masses)
            if abs(r - R_OBS) / R_OBS < 0.02:               # on the observed-mass-ratio surface
                sols.append((ang, d, r, masses))
    th13s = sorted({round(a[1], 1) for a, _, _, _ in sols})
    th12s = sorted({round(a[0], 1) for a, _, _, _ in sols})
    th23s = sorted({round(a[2], 1) for a, _, _, _ in sols})
    print(f"    solutions on the r_obs surface: {len(sols)};  theta13 range = [{min(th13s) if th13s else float('nan')}, {max(th13s) if th13s else float('nan')}] deg")
    print(f"    theta12 range = [{min(th12s)}, {max(th12s)}] ; theta23 range = [{min(th23s)}, {max(th23s)}]")
    # is theta13 PINNED (narrow) or a FAMILY (wide)?
    spread = (max(th13s) - min(th13s)) if th13s else 99.0
    ok(len(sols) > 0, "the pure portal DOES reach the observed mass ratio for some d (achievability)")
    pinned = spread < 2.0
    print(f"    -> theta13 spread on the mass-surface = {spread:.1f} deg  => {'PINNED (a prediction)' if pinned else 'a FAMILY (still underdetermined)'}")

    # nearest solution to the observed angle triple (honest best-case)
    if sols:
        def dist(a):
            return (a[0] - OBS[0]) ** 2 + (a[1] - OBS[1]) ** 2 + (a[2] - OBS[2]) ** 2
        bang, bd, br, _bm = min(sols, key=lambda s: dist(s[0]))
        print(f"    closest-to-observed solution: theta=({bang[0]:.1f},{bang[1]:.1f},{bang[2]:.1f}) at d={np.round(bd,3)} -- poor (obs {OBS})")

    # WHY underdetermined: the structure forbids a massless neutrino (det m_nu = -(1/2) prod d^2 != 0),
    # so the absolute scale m_lightest is FREE on the Delta-m^2 surface, and theta13 tracks it.
    dd = np.array([0.4, 0.7, 1.0])
    detval = np.linalg.det(m_nu(dd, 0.0, 1.0)); analytic = 0.5 * np.prod(dd ** 2)
    ok(abs(detval - analytic) < 1e-9 and abs(detval) > 1e-6,
       f"det(m_nu) = +(1/2)prod d^2 = {detval:.4f} != 0 -> NO massless neutrino; the lightest mass is a FREE scale")
    if sols:
        pairs = sorted((min(mm) / max(mm), a[1]) for a, _, _, mm in sols)        # (m_light/m_heavy, theta13)
        lo, hi = pairs[0], pairs[-1]
        print(f"    theta13 tracks the (free) absolute scale: at m_light/m_heavy={lo[0]:.3f} theta13={lo[1]:.1f}; at {hi[0]:.3f} theta13={hi[1]:.1f}")
        print("    -> the unknown lightest-neutrino mass is exactly the residual freedom carrying theta13.")

    print("\n[verdict] does the canon M_R portal predict theta13?")
    if pinned and th13s and abs(np.mean(th13s) - 8.6) < 3:
        print("  - YES (progress): fixing M_R to the pure K_3 portal + diagonal m_D + the observed masses PINS")
        print(f"    theta13 ~= {np.mean(th13s):.1f} deg, near the measured 8.6 -- the Casas-Ibarra freedom is removed.")
    elif pinned:
        print(f"  - PINNED but WRONG: the portal pins theta13 ~= {np.mean(th13s):.1f} deg, NOT 8.6 -> the minimal")
        print("    portal form is refuted as the theta13 source (clean negative); a different M_R real-part needed.")
    else:
        print("  - STILL UNDERDETERMINED: even with M_R fixed to the generation-blind portal, the mass-ratio")
        print("    surface carries a 1-parameter family of d with a RANGE of theta13 -> the masses do NOT pin")
        print("    the angle. The underdetermination is structural (the seesaw has more angle-freedom than the")
        print("    mass spectrum constrains), confirming canon l.2330 even after fixing M_R. Honest negative.")
    print("  exit 0")


if __name__ == "__main__":
    main()
