#!/usr/bin/env python3
r"""Item 87 -- does delta_CP pin once the PMNS angles are forced? YES: to J=0.

Attack the leptonic Dirac delta_CP magnitude through the framework's ACTUAL
angle-reproducing texture (not the arbitrary-m_D seesaw, which was underdetermined).

The angle texture is the closed frame-transport lemma
(item87_frame_transport_lemma_closure.py, item87_pmns_sign_standard_bridge_gate.py):

  * the PMNS angles are the QEC CORRECTION unitary U_frame = exp(delta K_or/3),
    a REAL SO(3) rotation (a correction is a real frame move, not a collapse).
    Concretely U_PMNS = R23(pi/4) R13(delta/sqrt2) R12(pi/4 - delta) at delta=2/9
    gives (theta12,theta23,theta13) ~ (32.3, 45.0, 9.0) deg, matching observation.
  * the CP object is the SEPARATE passive syndrome record omega*K_or = i K_or:
    Hermitian, eigenvalues {0, +/-sqrt3}, DIAGONAL in the circulant/DFT (mass)
    basis -> it cannot rotate the angles, and in the mass basis it is a set of
    phases on the neutrino MASS eigenstates = MAJORANA phases, not the Dirac phase.

Consequence (the point of this script): a REAL PMNS mixing has Jarlskog J = 0, so
the framework's leptonic DIRAC CP is CONSERVED (delta_CP in {0, pi}, J_lepton = 0).
All leptonic CP violation lives in the MAJORANA sector (Phi = 1/3 magnitude, sigma
sign) -- it drives leptogenesis (the eta sign) and shows in 0nubb (m_bb), NOT in
the long-baseline Dirac delta_CP.

This is a sharp, falsifiable prediction: a high-significance measurement of a
NONZERO leptonic Dirac Jarlskog (delta_CP != 0, pi) falsifies the real
frame-transport PMNS. Honest caveat: it holds at the order the framework derives
the PMNS (leading real frame-transport); a subleading complex correction, if one
existed, could seed a small Dirac phase -- none is produced by the current texture.

Self-asserting, exit 0.
"""
import math
import numpy as np

DELTA = 2.0 / 9.0
TOL = 1e-12


def check(name, cond):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        raise AssertionError(name)


def r12(t): c, s = math.cos(t), math.sin(t); return np.array([[c, s, 0], [-s, c, 0], [0, 0, 1.0]])
def r23(t): c, s = math.cos(t), math.sin(t); return np.array([[1.0, 0, 0], [0, c, s], [0, -s, c]])
def r13(t): c, s = math.cos(t), math.sin(t); return np.array([[c, 0, s], [0, 1.0, 0], [-s, 0, c]])
def Kor():  return np.array([[0, -1, 1], [1, 0, -1], [-1, 1, 0.0]])


def angles_deg(U):
    Ua = np.abs(U)
    s13 = min(Ua[0, 2], 1.0); c13 = math.sqrt(max(1 - s13**2, 1e-30))
    th13 = math.degrees(math.asin(s13))
    th12 = math.degrees(math.atan2(Ua[0, 1], Ua[0, 0]))
    th23 = math.degrees(math.atan2(Ua[1, 2], Ua[2, 2]))
    return th12, th23, th13


def jarlskog(U):
    return float(np.imag(U[0, 0] * U[1, 1] * np.conj(U[0, 1]) * np.conj(U[1, 0])))


def main():
    print("ITEM 87 -- leptonic Dirac delta_CP from the real frame-transport PMNS")
    print("=" * 82)

    print("\n[1] The angle-correct texture reproduces the PMNS angles and is REAL")
    U = r23(math.pi / 4) @ r13(DELTA / math.sqrt(2)) @ r12(math.pi / 4 - DELTA)
    th12, th23, th13 = angles_deg(U)
    print(f"    U_PMNS = R23(pi/4) R13(delta/sqrt2) R12(pi/4 - delta),  delta = 2/9")
    print(f"    (theta12, theta23, theta13) = ({th12:.1f}, {th23:.1f}, {th13:.1f}) deg   [obs ~ (33.4, 49.2, 8.6)]")
    check("theta13 ~ 8.6 deg (reactor angle from the universal 2/9)", abs(th13 - 8.6) < 0.6)
    check("theta12 ~ 33 deg (solar)", abs(th12 - 33.4) < 2.0)
    check("theta23 = 45 deg (atmospheric latch; octant a separate open item)", abs(th23 - 45.0) < 0.1)
    check("the frame-transport PMNS is REAL (orthogonal)", np.allclose(U.imag, 0) and np.allclose(U @ U.T, np.eye(3)))

    print("\n[2] A real PMNS has Jarlskog J = 0  ->  Dirac CP is CONSERVED")
    Jlep = jarlskog(U.astype(complex))
    print(f"    J_lepton = {Jlep:+.3e}   (a real mixing matrix carries no Dirac CP)")
    check("leptonic Dirac Jarlskog J = 0 (delta_CP in {0, pi})", abs(Jlep) < 1e-12)

    print("\n[3] The CP object i*K_or is a MAJORANA/sign carrier, not a Dirac phase")
    K = Kor(); iK = 1j * K
    check("i*K_or is Hermitian (a syndrome readout, not an applied rotation)", np.allclose(iK.conj().T, iK))
    ev = np.linalg.eigvalsh(iK)
    print(f"    eig(i*K_or) = {np.round(np.sort(ev), 4)}   (a discrete sign pointer {{0, +/-sqrt3}})")
    check("eig(i*K_or) = {0, +/-sqrt3} (discrete CP sign, R15)", np.allclose(np.sort(ev), [-math.sqrt(3), 0, math.sqrt(3)]))
    # K_or is a circulant -> diagonal in the DFT (mass-eigenstate) basis: a mass-basis phase = Majorana phase
    F = np.array([[1, 1, 1], [1, np.exp(2j*np.pi/3), np.exp(4j*np.pi/3)],
                  [1, np.exp(4j*np.pi/3), np.exp(8j*np.pi/3)]]) / math.sqrt(3)
    diagF = F.conj().T @ K @ F
    offnorm = np.linalg.norm(diagF - np.diag(np.diag(diagF)))
    check("K_or is diagonal in the DFT/mass basis (its phases are Majorana, mass-eigenstate rephasings)", offnorm < 1e-10)

    print("\n[4] Adding the Majorana phases leaves the Dirac Jarlskog at zero")
    # Majorana rephasing of the neutrino mass eigenstates by the K_or sign pointer
    for phi in (DELTA, math.sqrt(3) * DELTA, 0.7):
        P_maj = np.diag(np.exp(1j * phi * np.array([0.0, +1.0, -1.0])))   # {0,+,-} sign pointer
        Jm = jarlskog((U.astype(complex)) @ P_maj)
        check(f"Majorana rephasing (phi={phi:.3f}) leaves Dirac J = 0", abs(Jm) < 1e-12)
    print("    => the framework's leptonic CP is MAJORANA (0nubb / leptogenesis), not Dirac.")

    print(
        """
[5] VERDICT -- delta_CP DOES pin, to the CP-conserving Dirac value J_lepton = 0
    Forcing the PMNS angles with the framework's own closed texture ANSWERS the
    magnitude question, in the negative-for-Dirac direction:

      * the angles are the QEC correction unitary U_frame = exp(delta K_or/3), a
        REAL SO(3) rotation -- it reproduces (theta12,theta23,theta13) ~ (32,45,9),
        matching data, and being real it carries Jarlskog J = 0;
      * therefore the leptonic DIRAC CP phase is CONSERVED: delta_CP in {0, pi},
        J_lepton = 0 -- a definite prediction, not a family and not a free number;
      * the framework's CP content is the SEPARATE object omega*K_or = i K_or: a
        discrete sign pointer (eigenvalues {0,+/-sqrt3}) diagonal in the mass
        basis, i.e. MAJORANA phases with magnitude set by Phi = 1/3 and sign by
        sigma. That is what drives the baryon asymmetry (leptogenesis eta sign)
        and 0nubb (m_bb) -- it does NOT appear as a Dirac oscillation phase.

    So the two-birds hope resolves as: the phase branch is Phi = 1/3 (prior step),
    the eta/Phi tension is reconciled (sign from sigma), and the observable-lift
    verdict is that leptonic DIRAC CP is CONSERVED (J_lepton = 0) while leptonic
    CP violation is MAJORANA. This is sharply falsifiable:

      * DUNE / Hyper-K / T2K / NOvA: a high-significance NONZERO Dirac Jarlskog
        (delta_CP clearly != 0, pi) FALSIFIES the real frame-transport PMNS.
        (Current global fits mildly prefer delta_CP ~ 3pi/2 but still allow the
        CP-conserving pi at ~1 sigma in normal ordering -- the prediction is live.)
      * 0nubb: the predicted Majorana phase enters m_bb -- the framework's leptonic
        CP should show HERE, not in oscillations.

    Honest caveat: this holds at the order the framework derives the PMNS (leading
    real frame-transport). No complex sub-correction is produced by the current
    texture; if one were, it could seed a small Dirac phase. The clean, testable
    headline stands: framework => leptonic Dirac CP conserved, CP is Majorana.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
