#!/usr/bin/env python3
r"""Item 87 -- leptonic delta_CP: pin the portal-phase branch, reconcile eta/Phi,
and test the lift to a physical PMNS Dirac delta_CP.

Where we were (two existing scripts, now combined):
  * item87_r1_orientation_leptogenesis_bridge.py: I_CP ~ sin(3 sigma Phi); the
    faithful C3 character Phi=2pi/3 is a CP ZERO, so it is excluded as the phase.
    It left BOTH Phi=1/3 and Phi=2pi/9 as CP-nonzero candidates ("Phi magnitude
    not derived").
  * item87_neutrino_phase_selection_audit.py: the charged-lepton fit forces raw
    RADIAN defect fractions d/N (not 2pi-windings); consistency gives delta_nu=1/3
    which ALONE matches the NuFIT neutrino mass^2 ratio, while 2pi/9 (and 2pi/3)
    MISS it badly.

New here:
  [1]-[3] COMBINE the two exclusions: 2pi/3 killed by CP-zero, 2pi/9 killed by the
          neutrino mass^2 ratio => Phi = 1/3 rad is the UNIQUE survivor (given the
          one-record thesis Phi and delta_nu share magnitude; sym->delta, antisym->Phi).
  [4]     RECONCILE the eta/Phi tension: sign(I_CP)=sign(sigma) for every Phi in
          (0, pi/3), so it is BRANCH-INDEPENDENT. The old baryogenesis-eta result
          (s=+1, "Phi=2pi/3") took its SIGN from sigma, never from the magnitude;
          re-pinning 2pi/3 -> 1/3 leaves sign(eta)=sign(J) intact. The 2pi/3 was a
          mislabel of a CP-zero, not a load-bearing input to the sign.
  [5]     LIFT test to the observable Dirac delta_CP through the canon seesaw
          m_nu = m_D M_R^{-1} m_D^T with the complex K3 portal
          M_R = a I + b e^{i sigma Phi} A_K3, m_D=diag(d). HONEST result:
            (5a) pure portal a=0 => Dirac J = 0 (Phi becomes a removable global
                 phase) => leptonic Dirac CP REQUIRES a bare nu_R diagonal mass a!=0;
            (5b) with a!=0, delta_CP is a FAMILY on the mass-ratio surface (it
                 tracks the underdetermined lightest-neutrino mass, exactly the
                 freedom that leaves theta_13 unpinned in
                 foundations_pmns_theta13_from_K3_portal.py) -> the MAGNITUDE does
                 not lift cleanly, and the physical Jarlskog sign is not robust in
                 this minimal texture.

Net: the phase BRANCH (1/3) and the eta/Phi reconciliation are pinned; the
existence of leptonic CP is conditional on a!=0; the observable delta_CP MAGNITUDE
is an open residual, now with a sharp obstruction (texture + theta_13
underdetermination), not a vague hope. Honest, self-asserting, exit 0.
"""
import math
import cmath
import numpy as np

R2 = math.sqrt(2)
me, mmu, mtau = 0.51099895, 105.6583755, 1776.86            # MeV, PDG
NUFIT, NUFIT_SIG = 0.02951, 0.00088                         # Dm2_21/Dm2_31 (NO)
J3 = np.ones((3, 3)); I3 = np.eye(3); AK3 = J3 - I3
R_OBS = 7.42e-5 / 2.515e-3                                  # ~0.0295


def check(name, cond):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        raise AssertionError(name)


# ---- phase readings ----------------------------------------------------------
def koide_factors(delta, R):
    return np.array([(1 + R * math.cos(delta + 2 * math.pi * n / 3)) ** 2 for n in range(3)])


def nu_ratio(delta, R=1.0):
    f = np.sort(koide_factors(delta, R))
    denom = f[2] ** 2 - f[0] ** 2
    return (f[1] ** 2 - f[0] ** 2) / denom if denom > 1e-15 else 0.0


def I_cp_holonomy(phi, sigma, r=0.3):
    return (r ** 3) * math.sin(3.0 * sigma * phi)


# ---- seesaw lift -------------------------------------------------------------
def m_nu(d, a, b, phi, sigma):
    MR = a * I3 + b * cmath.exp(1j * sigma * phi) * AK3
    D = np.diag(d).astype(complex)
    return D @ np.linalg.inv(MR) @ D                        # complex symmetric


def takagi(M):
    w, U = np.linalg.eigh(M @ M.conj().T)
    s = np.sqrt(np.clip(w, 0, None))
    ph = np.angle(np.diag(U.T @ M @ U))
    return s, U @ np.diag(np.exp(-1j * ph / 2))             # M = U diag(s) U^T


def jarlskog(U, masses):
    U = U[:, np.argsort(masses)]
    return float(np.imag(U[0, 0] * U[1, 1] * np.conj(U[0, 1]) * np.conj(U[1, 0])))


def dcp_from(U, masses):
    U = U[:, np.argsort(masses)]
    Ua = np.abs(U)
    s13 = min(Ua[0, 2], 1.0); c13 = math.sqrt(max(1 - s13 ** 2, 1e-30))
    th12 = math.atan2(Ua[0, 1], Ua[0, 0]); th23 = math.atan2(Ua[1, 2], Ua[2, 2])
    Jinv = jarlskog(U, masses)
    denom = math.sin(th12) * math.cos(th12) * math.sin(th23) * math.cos(th23) * s13 * c13 ** 2
    sd = max(-1.0, min(1.0, Jinv / denom)) if abs(denom) > 1e-12 else 0.0
    return math.degrees(math.asin(sd)), math.degrees(math.asin(s13))


def mass_ratio(masses):
    m = np.sort(masses)
    return (m[1] ** 2 - m[0] ** 2) / (m[2] ** 2 - m[0] ** 2)


def scan_dcp(a, b, phi=1.0 / 3.0, sigma=1, tol=0.03):
    g = np.concatenate([np.linspace(0.05, 1.0, 80), np.linspace(1.05, 4.0, 55)])
    dcps = []
    for d1 in g:
        for d2 in g:
            s, U = takagi(m_nu(np.array([d1, d2, 1.0]), a, b, phi, sigma))
            if abs(mass_ratio(s) - R_OBS) / R_OBS < tol:
                dcps.append(dcp_from(U, s)[0])
    return np.array(dcps)


def main():
    print("ITEM 87 -- LEPTONIC delta_CP: branch pin, eta reconciliation, lift test")
    print("=" * 84)

    print("\n[1] CP-zero exclusion (from the orientation bridge): I_CP ~ sin(3 sigma Phi)")
    for name, phi, nonzero in [("Phi=1/3", 1/3, True), ("Phi=2pi/9", 2*math.pi/9, True),
                               ("Phi=2pi/3 (faithful C3)", 2*math.pi/3, False)]:
        val = I_cp_holonomy(phi, +1)
        print(f"    {name:24s}: I_CP(+) = {val:+.4e}")
        check(f"{name}: CP {'nonzero' if nonzero else 'ZERO'}", (abs(val) > 1e-8) == nonzero)

    print("\n[2] Neutrino mass^2 ratio excludes the 2pi/9 winding, selects 1/3")
    r13, r29, r23 = nu_ratio(1/3), nu_ratio(2*math.pi/9), nu_ratio(2*math.pi/3)
    print(f"    delta_nu=1/3   -> ratio {r13:.5f}  (NuFIT {NUFIT}+/-{NUFIT_SIG})")
    print(f"    delta_nu=2pi/9 -> ratio {r29:.5f}   delta_nu=2pi/3 -> ratio {r23:.5f}")
    check("1/3 matches NuFIT within 3 sigma", abs(r13 - NUFIT) < 3 * NUFIT_SIG)
    check("2pi/9 misses NuFIT by >10 sigma", abs(r29 - NUFIT) > 10 * NUFIT_SIG)

    print("\n[3] UNIQUE survivor: Phi = 1/3 rad (2pi/3 killed by CP-zero, 2pi/9 by mass^2)")
    survivors = []
    for phi, tag in [(1/3, "1/3"), (2*math.pi/9, "2pi/9"), (2*math.pi/3, "2pi/3")]:
        cp_ok = abs(I_cp_holonomy(phi, +1)) > 1e-8
        mass_ok = abs(nu_ratio(phi) - NUFIT) < 3 * NUFIT_SIG
        if cp_ok and mass_ok:
            survivors.append(tag)
    print(f"    survivors (CP-nonzero AND mass^2-consistent): {survivors}")
    check("Phi=1/3 is the unique survivor", survivors == ["1/3"])

    print("\n[4] eta/Phi reconciliation: sign(I_CP)=sign(sigma) is BRANCH-INDEPENDENT")
    # for any Phi in (0, pi/3), 3 Phi in (0, pi) so sin(3 Phi) > 0 => sign(I_CP)=sign(sigma)
    for phi in (1/3, 2*math.pi/9):
        sp, sm = I_cp_holonomy(phi, +1), I_cp_holonomy(phi, -1)
        check(f"Phi={phi:.3f}: sign(I_CP)=sign(sigma) (sp>0, sm<0)", sp > 0 and sm < 0 and abs(sp + sm) < 1e-15)
    check("3*Phi in (0,pi) for both survivors -> sin>0 -> sign is sigma, not the magnitude",
          0 < 3 * (1/3) < math.pi and 0 < 3 * (2*math.pi/9) < math.pi)
    print("    => the old eta-sign result took its SIGN from sigma; 'Phi=2pi/3' was a")
    print("       mislabel of a CP-zero. Re-pinning 2pi/3 -> 1/3 leaves sign(eta)=sign(J) intact.")

    print("\n[5] LIFT to the observable Dirac delta_CP via the canon K3-portal seesaw")
    s, U = takagi(m_nu(np.array([0.4, 0.7, 1.0]), 0.0, 1.0, 1/3, 1))     # pure portal a=0
    J0 = jarlskog(U, s)
    print(f"    (5a) pure portal a=0: Dirac J = {J0:+.2e}  -> CP conserved")
    check("a=0 pure off-diagonal portal gives Dirac J=0 (Phi is a removable global phase)", abs(J0) < 1e-9)
    print("         => leptonic Dirac CP REQUIRES a bare nu_R diagonal mass, a != 0.")

    dcps = scan_dcp(1.0, 1.0)                                             # a!=0 general portal
    spread = float(dcps.max() - dcps.min()) if len(dcps) else 0.0
    print(f"    (5b) a!=0 (a=1,b=1): {len(dcps)} mass-ratio solutions; "
          f"delta_CP in [{dcps.min():+.0f},{dcps.max():+.0f}] deg, spread {spread:.0f} deg")
    check("with a!=0 delta_CP is a FAMILY (spread > 10 deg), not pinned by the mass ratio", spread > 10.0)
    print("         => the magnitude tracks the underdetermined lightest-neutrino mass")
    print("            (the same freedom that leaves theta_13 unpinned); MAGNITUDE not forced.")

    print(
        """
[6] VERDICT -- real progress + a sharp, honest boundary
    PINNED (new, from combining existing results):
      * Phi = 1/3 rad is the UNIQUE portal-phase branch: 2pi/3 is a CP zero
        (orientation bridge), 2pi/9 fails the neutrino mass^2 ratio (phase audit),
        and the one-record thesis ties Phi's magnitude to delta_nu=1/3.
      * eta/Phi reconciled: the baryogenesis sign(eta)=sign(J) is carried by the
        orientation sigma and is branch-independent for Phi in (0,pi/3). The prior
        'Phi=2pi/3' was a mislabel of a CP-zero, not the source of the sign, so the
        eta-sign result SURVIVES the re-pin to Phi=1/3. (Closes the flagged tension.)

    CONDITIONAL / OPEN (the honest boundary of the lift):
      * leptonic Dirac CP is nonzero ONLY IF the nu_R carry a bare diagonal Majorana
        mass (a != 0); the pure DeltaL=2 portal (a=0) gives delta_CP = 0 exactly.
      * even with a != 0, the OBSERVABLE delta_CP magnitude does not lift cleanly in
        the minimal diagonal-m_D portal texture: it is a family tracking the
        underdetermined lightest-neutrino mass (the theta_13 freedom), and its
        Jarlskog sign is not robust in this texture. Pinning the magnitude needs the
        framework's actual angle-reproducing texture (TBM democracy + theta_13
        frame-transport) AND resolution of the lightest-mass underdetermination.

    SURVIVING FALSIFIABLE STATEMENT (weaker than a number, but real):
      the framework predicts leptonic CP is VIOLATED (delta_CP not 0 or pi),
      CONDITIONAL on a bare nu_R diagonal Majorana mass -- a clean yes/no that
      T2K/NOvA/DUNE/Hyper-K/JUNO test. The specific delta_CP value is NOT yet forced;
      that is the next rung (texture + lightest-mass), now sharply localised.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
