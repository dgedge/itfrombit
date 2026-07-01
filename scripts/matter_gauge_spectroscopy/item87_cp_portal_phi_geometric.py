#!/usr/bin/env python3
r"""Item 87 -- the genuine CP phase Phi: a geometric (holonomy-quantized) phase of
the Z3 generation bundle, NOT the CP-trivial mass-shape delta.

Pivot.  delta (mass circulant) is established CP-TRIVIAL (sets real masses, raw
d/N radians, winding EXCLUDED by the lepton fit).  Phi (the K3 Majorana portal,
M_H = I + r e^{i sigma Phi} A_K3) is the CP-VIOLATING baryogenesis phase -- a
different object.  The canon conflates them by testing Phi = delta_nu = 1/3, and
floats three values (1/3, 2pi/9, 2pi/3).  This script applies the one principle
that distinguishes Phi from delta:

  Phi lives on the Z3 GENERATION BUNDLE.  For the generation wavefunction to be
  single-valued under the Z3 generation symmetry, transport around the K3
  generation triangle must return a Z3 element:  total holonomy = e^{3 i Phi} in
  {1, omega, omega^2}.  Hence  3 Phi = 2 pi k / 3  =>  Phi = 2 pi k / 9
  (the Z9 ladder).

This is the OPPOSITE of delta: for delta the winding/holonomy reading was excluded
(the fit needs raw d/N); for Phi the holonomy quantization is REQUIRED (it is a
bundle CP phase).  The ladder EXCLUDES Phi = 1/3 rad (not a Z3 holonomy) -- so the
Phi = delta_nu conflation fails on principle -- and keeps the geometric candidates
2pi/9 (holonomy omega) and 2pi/3 (holonomy 1, the C3 recovery character, R15).
"""
import math
import numpy as np

H_DIRAC = np.diag(np.array([0.5, 0.8, 1.3]) ** 2)  # weak-basis probe (item87)
OMEGA = np.exp(2j * math.pi / 3)
Z3 = np.array([1.0, OMEGA, OMEGA ** 2])


def check(cond, msg):
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def a_k3():
    return np.ones((3, 3), dtype=complex) - np.eye(3, dtype=complex)


def m_portal(phi, sigma=1, r=0.5):
    return np.eye(3, dtype=complex) + r * np.exp(1j * sigma * phi) * a_k3()


def cp_invariant(M):
    Md = M.conj().T
    return float(np.imag(np.trace(H_DIRAC @ Md @ M @ Md @ H_DIRAC @ M)))


def triangle_holonomy(M):
    return M[0, 1] * M[1, 2] * M[2, 0]


def in_Z3(z):
    return float(np.min(np.abs(z / abs(z) - Z3)))  # distance of the unit phase to nearest Z3 element


def main():
    print("ITEM 87 -- Phi as a Z3-bundle geometric CP phase")
    print("=" * 80)

    cands = {
        "0 (CP-conserving)": 0.0,
        "delta_nu = 1/3 rad": 1.0 / 3.0,
        "2pi/9 (k=1)": 2 * math.pi / 9,
        "pi/3": math.pi / 3,
        "2pi/3 (k=3, C3 char)": 2 * math.pi / 3,
        "pi (CP-conserving)": math.pi,
    }

    print("\n[1] CP content: Phi in {0, pi} conserve CP; generic Phi violates it")
    for name, phi in cands.items():
        I = cp_invariant(m_portal(phi))
        viol = abs(I) > 1e-6
        print(f"    Phi={name:<22s} I_CP = {I:+.4e}   CP-violating={viol}")
        if phi in (0.0, math.pi):
            check(abs(I) < 1e-9, f"Phi={name}: CP-conserving (I_CP=0)")
        else:
            check(viol, f"Phi={name}: CP-violating")

    print("\n[2] Z3-bundle holonomy quantization: e^{3 i Phi} must be a Z3 element")
    for name, phi in cands.items():
        h = triangle_holonomy(m_portal(phi))
        dist = in_Z3(h)
        on = dist < 1e-9
        print(f"    Phi={name:<22s} holonomy e^(3iPhi) Z3-distance = {dist:.4f}  on-ladder={on}")
    # the discriminating checks:
    check(in_Z3(triangle_holonomy(m_portal(1.0 / 3.0))) > 1e-2,
          "Phi=1/3 rad is NOT on the Z3 ladder -> EXCLUDED as a bundle CP phase")
    check(in_Z3(triangle_holonomy(m_portal(2 * math.pi / 9))) < 1e-9,
          "Phi=2pi/9 IS on the ladder (holonomy = omega)")
    check(in_Z3(triangle_holonomy(m_portal(2 * math.pi / 3))) < 1e-9,
          "Phi=2pi/3 IS on the ladder (holonomy = 1, the C3 character)")

    print("\n[3] Opposite quantization of the two phases (they are different objects)")
    # delta_mass: raw d/N radians, winding EXCLUDED by the lepton fit (prior scripts).
    # Phi_CP:     holonomy-quantized to the Z9 ladder, winding REQUIRED.
    check(in_Z3(triangle_holonomy(m_portal(1.0 / 3.0))) > 1e-2
          and in_Z3(triangle_holonomy(m_portal(2 * math.pi / 9))) < 1e-9,
          "delta's value (1/3) is off-ladder; Phi's geometric values are on-ladder")

    print("\n[4] Orientation: the CP sign flips with sigma (baryogenesis carrier)")
    for phi in (2 * math.pi / 9, 2 * math.pi / 3):
        Ip, Im = cp_invariant(m_portal(phi, +1)), cp_invariant(m_portal(phi, -1))
        check(abs(Ip + Im) < 1e-9 and abs(Ip) > 1e-6, f"Phi={phi:.4f}: I_CP(+s) = -I_CP(-s) != 0")

    print(
        """
[5] VERDICT (a principled narrowing of Phi -- not a unique closure)
    Phi is the CP-VIOLATING Majorana phase of the Z3 generation bundle, a
    different object from the CP-TRIVIAL mass-shape delta.  Single-valuedness of
    the Z3 generation wavefunction forces the K3 triangle holonomy e^{3 i Phi}
    into Z3, i.e. Phi on the Z9 ladder Phi = 2 pi k / 9.  Consequences:
      * Phi = delta_nu = 1/3 rad is EXCLUDED: 1/3 is off-ladder (its holonomy
        e^{i} is not a Z3 element).  The canon's Phi = delta_nu identification is a
        category error -- it borrows the raw, off-ladder mass-shape value for a
        bundle CP phase.  This RESOLVES the three-way ambiguity by deleting the
        non-geometric candidate.
      * The surviving geometric candidates are Phi = 2pi/9 (holonomy omega,
        nontrivial) and Phi = 2pi/3 (holonomy 1, trivial; the R15 C3 recovery
        character).  Both are CP-violating and orientation-odd (good baryon-sign
        carriers).
      * delta and Phi are quantized by OPPOSITE principles -- delta raw (winding
        excluded by the fit), Phi holonomy-quantized (winding required) -- which is
        itself the reason they are distinct phases and should not be equated.

    MECHANISM SYNTHESIS (grounded in canon lines 351 / 743 / 3460): the two phases
    arise from OPPOSITE structures, which is why one is pi-FREE and the other pi-FUL.
      * delta (mass-shape, pi-free rational) is a CHIRAL ANOMALY coefficient.
        ANCHOR line 3460 already FALSIFIED the 2pi Aharonov-Bohm holonomy route
        (it gives delta=4pi/27, carrying pi) and hypothesised delta=2/9 from the
        Fujikawa anomaly -- whose pi^2 denominator natively cancels geometric pi to
        leave rationals -- with the chiral Q^3 anomaly already EXACTLY -2/9 (line
        351). My winding-excluded / k=1 / raw-radian result independently confirms
        the pi-free, non-holonomic character. (Anomaly mechanism: canon hypothesis-
        stage, but the right KIND, and a strictly better candidate than my earlier
        arc-subtension, which is superseded by it.)
      * Phi (CP, pi-ful holonomy) is the Z3 recovery-character phase 2pi/3
        (line 743, R15 C3 character) -- a genuine winding on the generation bundle.
      So: delta = anomaly (pi-free), Phi = holonomy (pi-ful). Opposite mechanisms,
      opposite pi-content; this is the deep reason the Phi = delta_nu conflation in
      the portal iff-script fails, and it reconciles the canon's own split (the
      portal script floats 1/3 & 2pi/9; R15 + the baryogenesis sign use 2pi/3).

    REMAINING (open): pin 2pi/9 vs 2pi/3 -- the nontrivial-holonomy (Z9 faithful)
    vs trivial-holonomy (C3 character) reading -- from boot mechanics, and derive
    the global orientation sigma (R15, conditionally CKM-tied).  This narrows item
    87's CP gate from "derive a free phase" to "choose one of two Z3-ladder rungs
    + an orientation"; it does not close it.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
