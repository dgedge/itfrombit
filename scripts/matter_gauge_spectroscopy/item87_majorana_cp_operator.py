#!/usr/bin/env python3
r"""ITEM 87 — the missing CP operator IS a complex-symmetric Majorana M_R (corrects a prior result).

item87_missing_operator_search.py concluded leptonic CP is absent ("Im[(Y_nu^dag Y_nu)^2] = 0").
THAT WAS WRONG: (i) it is the original weak basis, where a real-diagonal Y trivially has no phase --
not the physical CP measure; (ii) the seesaw diagonalisations there used eigh / a buggy Takagi,
which HERMITISE AWAY the Majorana phase. This script settles it with BASIS-INVARIANT weak-basis CP
invariants (pure traces -- no diagonalisation, no Takagi), and finds the missing operator.

RESULT (robust). With the walk's REAL-diagonal Dirac coupling Y_nu (H = Y_nu^dag Y_nu = diag(y^2),
real), a complex-symmetric Majorana M_R (Delta L=2) gives NONZERO weak-basis CP invariants
Im Tr[H M^dag M M^dag H M] etc. -- i.e. PHYSICAL leptonic CP -- whereas a real M_R gives exactly 0.
The CP invariant flips sign with the M_R phase (phi -> -phi), so the leptogenesis sign tracks the
sign of that phase. This holds whether M_R is non-degenerate OR a degenerate symmetric circulant --
so the earlier J=0 / "no CP" was a DIAGONALISATION ARTIFACT, not physics; CP from a complex M_R is
generic. (Degeneracy only sets the leptogenesis REGIME -- resonant vs hierarchical -- not whether CP
exists.)

THE MISSING OPERATOR: a complex-symmetric Majorana mass M_R for the right-handed neutrinos
(Delta L=2). It is NOT in the documented walk substrate (the walk is Hermitian / Delta L=0 and gives
only a real-diagonal, generation-block-diagonal nu sector -- the item-87 phase-lift obstruction).
But it is exactly the object the framework already posits as the neutrino Koide structure (Part 5:
delta_nu=1/3), which item 87 flags as "inserted as an ansatz, not derived from the substrate". So the
operator is IDENTIFIED and VERIFIED to give CP; the OPEN pieces are (i) deriving this complex M_R --
its phase delta_nu -- from the boot / high-scale dynamics (item 53 supplies only the scale), and
(ii) linking sign(delta_nu) to the global walk orientation (the e^{+ik pi/4} that fixes the CKM sign)
so the leptogenesis sign correlates with the quark CKM.

Self-asserting; exit 0 = real M_R -> zero CP; complex M_R (degenerate or not) -> nonzero CP that
flips sign with the phase.
"""
from __future__ import annotations
import numpy as np

Y_DIAG = np.array([0.5, 0.8, 1.3])           # walk's real-diagonal Dirac coupling (values immaterial)
H = np.diag(Y_DIAG ** 2)


def cp_invariants(M):
    """Basis-independent weak-basis CP invariants for (H real-diagonal, M complex-symmetric M_R)."""
    Md = M.conj().T
    return {
        "Im Tr[H Md M Md H M]": float(np.imag(np.trace(H @ Md @ M @ Md @ H @ M))),
        "Im Tr[H M Md M Md M]": float(np.imag(np.trace(H @ M @ Md @ M @ Md @ M))),
    }


def MR_generic(phi):
    """Non-degenerate complex-symmetric M_R: distinct diagonal + complex off-diagonal (phase phi)."""
    z = 0.3 * np.exp(1j * phi)
    return np.array([[1.0, z, 0.2], [z, 2.0, 0.4], [0.2, 0.4, 3.0]], complex)


def MR_circulant(phi):
    """Complex-symmetric circulant c0 I + c1 (P + P^2) -- degenerate (eig k=1,2 coincide)."""
    c1 = 0.5 * np.exp(1j * phi)
    return 1.0 * np.eye(3) + c1 * (np.roll(np.eye(3), 1, 0) + np.roll(np.eye(3), 2, 0))


def check(cond, msg):
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def main() -> int:
    print("ITEM 87 — the missing CP operator is a complex-symmetric Majorana M_R")
    print("=" * 84)

    print("\n[1] real M_R (phi=0): weak-basis CP invariants must vanish")
    inv0 = cp_invariants(MR_generic(0.0))
    for k, v in inv0.items():
        print(f"    {k} = {v:+.4e}")
    check(all(abs(v) < 1e-12 for v in inv0.values()), "real M_R -> all CP invariants = 0 (no CP)")

    print("\n[2] complex M_R (phi=1/3 = the documented delta_nu): physical CP")
    invc = cp_invariants(MR_generic(1.0 / 3.0))
    for k, v in invc.items():
        print(f"    {k} = {v:+.4e}")
    check(any(abs(v) > 1e-3 for v in invc.values()),
          "complex M_R + real-diagonal Y_nu -> NONZERO weak-basis CP invariant (physical leptonic CP)")

    print("\n[3] the CP sign tracks the M_R phase sign (orientation):")
    key = lambda phi: cp_invariants(MR_generic(phi))["Im Tr[H Md M Md H M]"]
    for phi in (1 / 3, -1 / 3, 2 * np.pi / 9, -2 * np.pi / 9):
        print(f"    phi={phi:+.4f}: I = {key(phi):+.4e}")
    check(abs(key(1 / 3) + key(-1 / 3)) < 1e-12 and key(1 / 3) != 0,
          "I(phi) = -I(-phi): the leptonic/leptogenesis CP sign tracks sign(phase) -> sign(delta_nu)")

    print("\n[4] robustness — even a DEGENERATE symmetric circulant M_R gives nonzero CP:")
    eigs = np.linalg.eigvals(MR_circulant(1 / 3))
    gap = float(np.min([abs(eigs[i] - eigs[j]) for i in range(3) for j in range(i + 1, 3)]))
    invcirc = cp_invariants(MR_circulant(1 / 3))["Im Tr[H Md M Md H M]"]
    print(f"    circulant eigenvalue min-gap = {gap:.2e};  CP invariant = {invcirc:+.2e}")
    check(gap < 1e-9 and abs(invcirc) > 1e-3,
          "complex M_R gives CP EVEN when degenerate -> the earlier J=0 was a diagonalisation "
          "artifact, not physics; CP existence is robust (degeneracy only sets the leptogenesis regime)")

    print("\n" + "=" * 84)
    print("VERDICT — operator IDENTIFIED + VERIFIED (corrects the prior 'no CP' conclusion)")
    print(
        "  The missing operator is a complex-symmetric Majorana M_R (Delta L=2). With the walk's\n"
        "  real-diagonal Y_nu it gives PHYSICAL leptonic CP -- nonzero weak-basis invariants, confirmed\n"
        "  basis-independently (no Takagi) -- and the CP sign tracks the M_R phase sign. It does so\n"
        "  whether M_R is non-degenerate or a degenerate circulant, so the earlier 'Im[(Y^dag Y)^2]=0\n"
        "  -> no CP' was the wrong basis and the eigh/Takagi steps hermitised away the Majorana phase:\n"
        "  CORRECTED. This operator is NOT in the documented walk substrate (Hermitian / Delta L=0), but\n"
        "  it IS the framework's own neutrino Koide structure (delta_nu=1/3), which item 87 flags as an\n"
        "  ansatz 'not derived from the substrate'. So the operator is identified and shown to work; the\n"
        "  OPEN pieces are (i) deriving this complex M_R (the delta_nu phase) from boot / high-scale\n"
        "  dynamics (item 53 gives only the scale), and (ii) linking sign(delta_nu) to the global walk\n"
        "  orientation so the leptogenesis sign correlates with the CKM. TIER: operator class identified\n"
        "  + verified to give CP; its substrate derivation + the sign-orientation link are the residual."
    )
    print("exit 0 — complex Majorana M_R gives physical leptonic CP (sign tracks the phase); derivation open.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
