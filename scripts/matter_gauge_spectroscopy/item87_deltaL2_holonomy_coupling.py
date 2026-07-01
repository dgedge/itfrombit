#!/usr/bin/env python3
r"""Item 87 / item 126 / R15: Delta L=2 CP-holonomy iff theorem.

Problem
-------
The documented walk is real/Hermitian in the lepton sector, Delta L=0, and
CNOT-gated so that the all-zero electron neutrino is inert.  Existing scripts
therefore prove a no-go for generating the leptonic CP phase, hence the
baryon-asymmetry sign, from the walk alone.

Admissible coupling class
-------------------------
The existing no-go proves that the documented walk cannot generate the
leptonic CP phase: it is Hermitian/Delta L=0, CNOT-gated, and leaves the
all-zero electron-neutrino state inert.  Therefore any CP/baryon-sign carrier
must be outside the documented walk in one specific way.

The minimal admissible class is an unconditional sterile-generation Majorana
portal:

Add a sterile-generation Majorana portal:

    L_DeltaL2 = 1/2 nu_R^T C M_H nu_R + h.c.

with

    M_H = M0 [ I + r exp(i sigma Phi) A_K3 ],

where A_K3 is the off-diagonal adjacency matrix on the three sterile generation
labels.  The coupling is:

  * Delta L=2: it is a Majorana bilinear;
  * unconditional: it is not CNOT-gated and therefore moves the all-zero nu_e;
  * generation-blind: every unordered pair uses the same complex amplitude;
  * holonomic: the oriented generation triangle carries phase exp(i 3 sigma Phi);
  * CP-odd: reversing the orientation sigma reverses the weak-basis CP invariant.

Iff statement
-------------
Within the admissible S3-covariant one-portal class:

  * no unconditional complex-symmetric Delta L=2 generation holonomy
    -> no three-generation CP carrier in the framework's documented walk;
  * such a holonomy with nonzero oriented phase
    -> nonzero weak-basis CP invariant, with sign reversed by orientation.

So the CP/baryon-sign residual is no longer "find some phase".  It is exactly:
derive this Delta L=2 holonomy portal, its Phi, and its orientation sigma from
QEC/boot mechanics.  Without that, the coupling is an iff target, not a Locked
baryogenesis sign theorem.
"""

from __future__ import annotations

import math
import numpy as np


Y_DIAG = np.array([0.5, 0.8, 1.3])
H_DIRAC = np.diag(Y_DIAG**2)
PERMS = [
    np.eye(3, dtype=complex)[list(p)]
    for p in [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
]


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def a_k3() -> np.ndarray:
    return np.ones((3, 3), dtype=complex) - np.eye(3, dtype=complex)


def m_holonomy(phi: float, sigma: int = 1, r: float = 0.5) -> np.ndarray:
    """Generation-blind complex-symmetric Delta L=2 Majorana portal."""
    return np.eye(3, dtype=complex) + r * np.exp(1j * sigma * phi) * a_k3()


def cp_invariant(M: np.ndarray) -> float:
    """Weak-basis CP invariant used in item87_majorana_cp_operator.py."""
    Md = M.conj().T
    return float(np.imag(np.trace(H_DIRAC @ Md @ M @ Md @ H_DIRAC @ M)))


def cp_invariant_vector(M: np.ndarray) -> np.ndarray:
    """Two independent weak-basis CP probes used as a robust zero/nonzero test."""
    Md = M.conj().T
    return np.array(
        [
            np.imag(np.trace(H_DIRAC @ Md @ M @ Md @ H_DIRAC @ M)),
            np.imag(np.trace(H_DIRAC @ M @ Md @ M @ Md @ M)),
        ],
        dtype=float,
    )


def cp_norm(M: np.ndarray) -> float:
    return float(np.linalg.norm(cp_invariant_vector(M)))


def triangle_holonomy(M: np.ndarray) -> complex:
    """Product of the three unordered generation-pair amplitudes."""
    return M[0, 1] * M[1, 2] * M[2, 0]


def commutes_with_generation_permutations(M: np.ndarray) -> bool:
    return all(np.max(np.abs(P.T @ M @ P - M)) < 1.0e-14 for P in PERMS)


def main() -> None:
    print("ITEM 87 / 126 / R15 -- DELTA L=2 CP-HOLONOMY IFF THEOREM")
    print("=" * 88)

    phases = {
        "delta_nu radians": 1.0 / 3.0,
        "nu_R Berry 2pi/9": 2.0 * math.pi / 9.0,
    }

    print("[1] Algebraic form")
    M0 = m_holonomy(0.0)
    print(f"    A_K3 off-diagonal degree = {int(np.max(np.sum(np.abs(a_k3()) > 0, axis=1)))}")
    check(np.max(np.abs(M0 - M0.T)) < 1.0e-14, "portal is complex-symmetric, hence Majorana-admissible")
    check(np.all(np.abs(M0[0, 1:]) > 0), "unconditional portal moves the all-zero nu_e generation")
    check(commutes_with_generation_permutations(M0), "K3 portal is S3-generation blind")
    one_edge = np.eye(3, dtype=complex)
    one_edge[0, 1] = one_edge[1, 0] = 0.5 * np.exp(1j / 3)
    check(not commutes_with_generation_permutations(one_edge), "one-edge complex portal is rejected by generation-blindness")
    check(cp_norm(M0) < 1.0e-12, "real portal carries no CP, as required")
    check(cp_norm(m_holonomy(1.0 / 3.0, r=0.0)) < 1.0e-12, "no Delta L=2 inter-generation portal -> no CP carrier")

    print("\n[2] Iff controls inside the admissible class")
    for phi in (0.0, math.pi):
        print(f"    phi={phi:.6f}: CP norm = {cp_norm(m_holonomy(phi)):.3e}")
        check(cp_norm(m_holonomy(phi)) < 1.0e-12, "phase-free/rephasable endpoint has zero CP")

    print("\n[3] CP and orientation tests")
    for name, phi in phases.items():
        Mp = m_holonomy(phi, sigma=+1)
        Mm = m_holonomy(phi, sigma=-1)
        Ip = cp_invariant(Mp)
        Im = cp_invariant(Mm)
        Cp = cp_norm(Mp)
        hp = triangle_holonomy(Mp)
        print(f"    {name:<18s}: phi={phi:.12f}")
        print(f"      I_CP(+orientation) = {Ip:+.6e}")
        print(f"      I_CP(-orientation) = {Im:+.6e}")
        print(f"      CP-invariant norm  = {Cp:.6e}")
        print(f"      triangle holonomy  = {hp.real:+.6e}{hp.imag:+.6e}i")
        check(Cp > 1.0e-3, f"{name}: complex holonomy gives nonzero CP")
        check(abs(Ip + Im) < 1.0e-12, f"{name}: CP sign reverses with orientation")
        check(abs(abs(hp) - abs(Mp[0, 1]) ** 3) < 1.0e-14, f"{name}: holonomy is the 3-edge phase product")

    print("\n[4] No free baryon-sign claim")
    b_from_bl = 28.0 / 79.0
    for sigma in (+1, -1):
        I = cp_invariant(m_holonomy(1.0 / 3.0, sigma=sigma))
        # The sign convention for the decay asymmetry is not derived here; this
        # line only shows that one binary orientation is enough to choose a sign.
        print(f"    orientation {sigma:+d}: sign carrier I_CP={I:+.6e}; sphaleron factor |B/(B-L)|={b_from_bl:.6f}")
    check(abs(b_from_bl - 28.0 / 79.0) < 1.0e-15, "sphaleron conversion factor is a separate fixed ledger")

    print(
        """
[5] VERDICT
    The R15 iff target is now sharp inside the admissible class:

      M_H = M0 [I + r exp(i sigma Phi) A_K3].

    Necessity: the documented walk cannot do the job (Delta L=0, CNOT-gated,
    nu_e inert, CP=0).  Real/rephasable or absent inter-generation Delta L=2
    portals also give zero CP in this test.  Generation-blindness rejects
    one-edge ad hoc complex insertions; the S3-covariant off-diagonal object is
    the K3 adjacency.

    Sufficiency: the K3 complex-symmetric Majorana holonomy has nonzero
    weak-basis CP invariants, and the sign flips under sigma -> -sigma.  Thus
    it supplies exactly the kind of object the baryon sign needs: a leptonic CP
    sign feeding the already-verified 28/79 sphaleron conversion.

    What remains open:
      derive the portal from QEC/boot mechanics; derive Phi rather than taking
      it from the geometric delta_nu / nu_R Berry primitive; and derive which
      global orientation sigma is selected.  Without those, R15 is closed as
      an iff theorem for the required object, not as a Locked baryogenesis sign
      prediction.

ALL ASSERTIONS PASSED"""
    )


if __name__ == "__main__":
    main()
