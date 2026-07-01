#!/usr/bin/env python3
r"""Item 87 -- neutrino-phase selection audit: which reading of the Koide/portal
phase survives, and what the irreducible residual actually is.

The open gate (ANCHOR sec 5.8 / item 87) is "derive Phi, do not insert it."  The
canon carries THREE candidate phase values and tests two of them as equals:

    (A) raw defect-to-plaquette ratio in radians:  delta = d/N   (lepton 2/9, nu 3/9=1/3)
    (B) Berry / holonomy "fraction of a full turn": delta = 2*pi*(d/N)
    (C) single-valued generation-triangle holonomy: Phi = 2*pi/9  (3*Phi = 2*pi/3 = omega)

These are NOT interchangeable: (A) is a raw radian, (B)/(C) are windings (fractions
of 2*pi).  The point of this script: the framework's OWN verified charged-lepton
fit adjudicates between them, because the lepton mass spectrum is reproduced to
0.006-0.007% ONLY for one reading.  That collapses the phase ambiguity and
RE-LOCATES the residual.

Mechanism split (confirmed from item87_deltaL2_holonomy_coupling.py):
  * mass hierarchy  <- Dirac sector (seesaw); the Majorana portal does not set it;
  * CP invariant    <- the K3 complex-symmetric Majorana portal M_H = I + r e^{i s Phi} A_K3
                       (uniform edge phase -> degenerate Majorana spectrum, so it is
                        a CP carrier, not a mass-hierarchy source).
So "the neutrino phase" wears two hats (mass-circulant delta_nu and CP-portal Phi);
the audit tests whether one reading serves both.
"""
import math
import numpy as np

R2 = math.sqrt(2)
me, mmu, mtau = 0.51099895, 105.6583755, 1776.86          # MeV, PDG
NUFIT, NUFIT_SIG = 0.02951, 0.00088                        # Dm2_21/Dm2_31


def check(cond, msg):
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def koide_factors(delta, R):
    return np.array([(1 + R * math.cos(delta + 2 * math.pi * n / 3)) ** 2 for n in range(3)])


def lepton_errors(delta, R=R2):
    f = koide_factors(delta, R)
    order = np.argsort(f)
    pred = (mtau / f.max()) * f
    return abs(pred[order[0]] - me) / me, abs(pred[order[1]] - mmu) / mmu


def nu_ratio(delta, R=1.0):
    f = np.sort(koide_factors(delta, R))
    denom = f[2] ** 2 - f[0] ** 2
    return (f[1] ** 2 - f[0] ** 2) / denom if denom > 1e-15 else 0.0


# --- CP portal (verbatim from item87_deltaL2_holonomy_coupling.py) ---------------
H_DIRAC = np.diag(np.array([0.5, 0.8, 1.3]) ** 2)


def a_k3():
    return np.ones((3, 3), dtype=complex) - np.eye(3, dtype=complex)


def m_holonomy(phi, sigma=1, r=0.5):
    return np.eye(3, dtype=complex) + r * np.exp(1j * sigma * phi) * a_k3()


def cp_invariant(M):
    Md = M.conj().T
    return float(np.imag(np.trace(H_DIRAC @ Md @ M @ Md @ H_DIRAC @ M)))


def main():
    print("ITEM 87 -- NEUTRINO PHASE SELECTION AUDIT")
    print("=" * 84)

    print("\n[1] The verified lepton fit selects RAW RADIANS (d/N), not a winding")
    e_raw, mu_raw = lepton_errors(2 / 9)
    e_turn, mu_turn = lepton_errors(2 * math.pi * 2 / 9)
    print(f"    delta_l = 2/9 rad        -> m_e err {e_raw*100:.3f}%, m_mu err {mu_raw*100:.3f}%")
    print(f"    delta_l = 2*pi*(2/9) rad -> m_e err {e_turn*100:.1f}%, m_mu err {mu_turn*100:.1f}%")
    check(e_raw < 1e-4 and mu_raw < 1e-4, "reading (A) d/N RADIANS fits leptons to <0.01%")
    check(e_turn > 0.1, "reading (B) 2*pi*(d/N) winding FAILS the leptons (>10%)")

    print("\n[2] Same reading carries the neutrino mass ratio; the windings do not")
    r_raw = nu_ratio(1 / 3)
    r_turn = nu_ratio(2 * math.pi / 3)
    r_berry = nu_ratio(2 * math.pi / 9)
    print(f"    delta_nu = 1/3 rad   -> ratio {r_raw:.5f}  (NuFIT {NUFIT}+/-{NUFIT_SIG})")
    print(f"    delta_nu = 2*pi/3    -> ratio {r_turn:.5f}")
    print(f"    delta_nu = 2*pi/9    -> ratio {r_berry:.5f}")
    check(abs(r_raw - NUFIT) < 3 * NUFIT_SIG, "(A) delta_nu=1/3 rad fits NuFIT ratio within 3 sigma")
    check(abs(r_turn - NUFIT) > 10 * NUFIT_SIG, "(B) delta_nu=2*pi/3 winding FAILS the ratio")
    check(abs(r_berry - NUFIT) > 10 * NUFIT_SIG, "(C) delta_nu=2*pi/9 single-valued holonomy FAILS the ratio")

    print("\n[3] At the mass-selected phase, the K3 portal still carries CP (existence)")
    Ip = cp_invariant(m_holonomy(1 / 3, +1))
    Im = cp_invariant(m_holonomy(1 / 3, -1))
    print(f"    I_CP(+sigma) = {Ip:+.4e}   I_CP(-sigma) = {Im:+.4e}")
    check(abs(Ip) > 1e-3, "K3 portal at Phi=1/3 gives nonzero CP invariant")
    check(abs(Ip + Im) < 1e-12, "CP sign reverses with orientation sigma")

    print(
        """
[4] VERDICT
    The phase ambiguity collapses on the framework's OWN verified evidence:
    the charged-lepton spectrum (0.006-0.007%) forces delta = d/N read as RAW
    RADIANS (2/9), NOT as a fraction of a full turn.  Consistency then fixes the
    neutrino reading to delta_nu = 3/9 = 1/3 rad, which alone matches the NuFIT
    mass-squared ratio; the two genuine "winding" readings (2*pi/3 and the
    single-valued-holonomy 2*pi/9) BOTH miss the data badly.  So:

      * EXISTENCE is settled: the K3 Majorana portal supplies a nonzero,
        orientation-odd CP carrier (verified above) -- the "no inter-generation
        B_nu" obstruction is resolved by construction.
      * The PHASE READING is settled by the lepton fit: d/N in radians, which
        EXCLUDES the Berry/holonomy-as-2*pi-fraction route (it would quantize to
        2*pi*k/9 and miss every mass observable).  Given the shared nu_R-defect
        origin, Phi = delta_nu = 1/3 selects the baryon CP sign up to sigma.

    IRREDUCIBLE RESIDUAL (sharpened, honest):
      (a) WHY a dimensionless defect fraction d/N equals a phase in RADIANS.
          This is NOT a winding/holonomy (those are excluded here), so the
          standard Berry-phase derivation cannot supply it; a non-winding
          "defect-fraction = radian" mechanism is required and is not derived.
      (b) the global orientation sigma (R15 sign pointer; conditionally tied to
          the CKM sign, ANCHOR sec 5.8 R15 addendum).
    This does NOT close item 87.  It removes two of the three candidate phases on
    evidence and re-points the open problem from "derive Phi" to "derive the
    d/N-as-radians identification + sigma".  Phi=1/3 is still selected by fit
    consistency, not derived from boot mechanics.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
