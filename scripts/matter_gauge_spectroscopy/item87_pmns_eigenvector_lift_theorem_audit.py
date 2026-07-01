#!/usr/bin/env python3
r"""ITEM 87 -- PMNS eigenvector-lift theorem audit.

Question
--------
Can the current framework prove that lepton-sector recovery bills K_or as a
real so(3) PMNS eigenvector rotation with mean-cycle normalization K_or/3,
rather than merely as the Hermitian i K_or CP-orientation pointer?

Answer
------
Not yet from the current canon.  The possible proof route is now precise.

The existing R15 work proves a sign readout:

    omega K_or

is an S3-scalar Stinespring record when omega is the global orientation line.
That closes the CP-orientation pointer conditionally.  It does not, by itself,
prove that the same record is billed as a coherent real frame update.

There are three inequivalent maps:

  1. CP/mass insertion: H = i K_or is Hermitian and circulant.  It is
     diagonalized by the DFT and does not produce the PMNS QLC eigenvectors.

  2. Reset/jump readout: a directed-cycle Kraus/jump is either a conditioned
     120-degree generation permutation or, unconditioned, a dissipative CPTP
     channel.  It is not the small coherent connection exp(delta K_or/3).

  3. Frame transport: a coherent correction/polar unitary

         U_frame(delta) = exp(delta K_or/3)

     is a real SO(3) rotation on the generation frame.  Combined with the
     standard shear and the atmospheric latch, it gives the QLC tangent.

Thus the missing theorem is not "K_or is readable" (already conditional) and
not "K_or exists" (already found).  It is:

    The lepton recovery instrument has a coherent frame-transport polar factor,
    applied before reset, whose sign-sector generator is the mean boundary
    current K_or/3; the low-energy PMNS matrix bills this frame holonomy.

Representation theory cannot fix the scalar 1/3: q K_or transforms correctly
for every real q.  The factor 1/3 needs a service-normalization lemma, such as
"one completed oriented boundary winding is averaged over the three R1
generation corners".
"""

from __future__ import annotations

from itertools import permutations
from pathlib import Path
import math

import numpy as np
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parents[1]
TOL = 1.0e-10
DELTA = 2.0 / 9.0


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def require_text(path: str, needles: tuple[str, ...]) -> None:
    text = (ROOT / path).read_text(encoding="utf-8")
    missing = [needle for needle in needles if needle not in text]
    if missing:
        raise AssertionError(f"{path} missing required text: {missing}")


def perms() -> list[tuple[tuple[int, ...], np.ndarray, int]]:
    out: list[tuple[tuple[int, ...], np.ndarray, int]] = []
    eye = np.eye(3)
    for perm in permutations(range(3)):
        inv = sum(1 for i in range(3) for j in range(i + 1, 3) if perm[i] > perm[j])
        out.append((perm, eye[list(perm)], -1 if inv % 2 else 1))
    return out


PERMS = perms()


def skew(axis: np.ndarray) -> np.ndarray:
    x, y, z = axis
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])


def k_or() -> np.ndarray:
    return skew(np.ones(3))


def cycle() -> np.ndarray:
    return np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])


def dft3() -> np.ndarray:
    w = np.exp(2j * np.pi / 3.0)
    return np.array([[w ** (j * k) for k in range(3)] for j in range(3)], dtype=complex) / np.sqrt(3.0)


def angles_deg(u: np.ndarray) -> np.ndarray:
    s13sq = abs(u[0, 2]) ** 2
    s12sq = abs(u[0, 1]) ** 2 / max(1.0 - s13sq, 1.0e-15)
    s23sq = abs(u[1, 2]) ** 2 / max(1.0 - s13sq, 1.0e-15)
    return np.array(
        [np.degrees(np.arcsin(np.sqrt(np.clip(x, 0.0, 1.0)))) for x in (s12sq, s23sq, s13sq)]
    )


def fmt(values: np.ndarray) -> tuple[float, float, float]:
    return tuple(float(round(v, 3)) for v in values)


def r12(theta: float) -> np.ndarray:
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])


def r23(theta: float) -> np.ndarray:
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]])


def r13(theta: float) -> np.ndarray:
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


def base_matrix() -> np.ndarray:
    return r23(np.pi / 4.0) @ r12(np.pi / 4.0)


def qlc_matrix() -> np.ndarray:
    return r23(np.pi / 4.0) @ r13(DELTA / math.sqrt(2.0)) @ r12(np.pi / 4.0 - DELTA)


def derivative_coefficients(axis: np.ndarray) -> np.ndarray:
    x, y, z = axis
    return np.array([-z, -(x + y) / math.sqrt(2.0), (y - x) / math.sqrt(2.0)])


def standard_axis_for_q(q: float) -> np.ndarray:
    return np.array([-q - 0.5, 0.5 - q, 2.0 * q])


def is_sign_operator(mat: np.ndarray) -> bool:
    return all(np.linalg.norm(p.T @ mat @ p - parity * mat) < TOL for _perm, p, parity in PERMS)


def section_imported_canon() -> None:
    print("\n[1] Imported canon clauses")
    require_text(
        "python_code/r15_global_orientation_sign_pointer_theorem.py",
        (
            "global substrate-orientation line",
            "omega*K_R1 is an S3-scalar environment record",
            "source remains one coherent generation-singlet port",
        ),
    )
    require_text(
        "python_code/item87_pmns_sign_standard_bridge_gate.py",
        (
            "normalized as K_or/3",
            "eigenvector-lift/normalization",
        ),
    )
    require_text(
        "python_code/item123_sterile_generation_singlet_source.py",
        (
            "the source is one coherent port",
            "A generation-blind environment has one Stinespring label",
        ),
    )
    require_text(
        "python_code/item87_lepton_cp_obstruction.py",
        (
            "Leptonic CP vanishes in every physical lepton sector",
            "existence gap",
        ),
    )
    check(True, "readout, source-port, and PMNS bridge audits are present")


def section_cp_pointer_not_lift() -> None:
    print("\n[2] Hermitian CP pointer is not the PMNS eigenvector lift")
    k = k_or()
    h_cp = 1j * k / 3.0
    f = dft3()
    diag = f.conj().T @ h_cp @ f
    offdiag = np.linalg.norm(diag - np.diag(np.diag(diag)))
    print(f"  ||F^dag (iK_or/3) F - diagonal|| = {offdiag:.3e}")
    print(f"  DFT angles                         = {fmt(angles_deg(f))}")
    check(np.linalg.norm(h_cp - h_cp.conj().T) < TOL, "i K_or/3 is Hermitian")
    check(offdiag < TOL, "i K_or/3 stays in the circulant/DFT mass-operator algebra")
    check(abs(angles_deg(f)[2] - 35.26438968) < 1.0e-6, "DFT eigenvectors give trimaximal theta13, not PMNS QLC")


def section_reset_not_lift() -> None:
    print("\n[3] Directed reset/jump readout is not the small coherent frame connection")
    p = cycle()
    k = k_or()
    theta_120 = 2.0 * math.pi / (3.0 * math.sqrt(3.0))
    p_from_k = expm(theta_120 * k)
    p_from_minus_k = expm(-theta_120 * k)
    close = min(np.linalg.norm(p_from_k - p), np.linalg.norm(p_from_minus_k - p))
    print(f"  nearest ||exp(+-2pi K/(3sqrt3)) - P_cycle|| = {close:.3e}")
    check(close < TOL, "a conditioned directed-cycle correction is a 120-degree permutation")

    ket0 = np.array([1.0, 0.0, 0.0])
    rho0 = np.outer(ket0, ket0)
    eps = 0.01
    rho_uncond = (1.0 - eps) * rho0 + eps * (p @ rho0 @ p.T)
    purity = float(np.trace(rho_uncond @ rho_uncond))
    rho_frame = expm(eps * k / 3.0) @ rho0 @ expm(eps * k / 3.0).T
    frame_purity = float(np.trace(rho_frame @ rho_frame))
    print(f"  unconditioned jump purity at eps=0.01 = {purity:.6f}")
    print(f"  coherent frame purity at eps=0.01     = {frame_purity:.6f}")
    check(purity < 1.0 - 1.0e-4, "unconditioned reset/jump channel is dissipative")
    check(abs(frame_purity - 1.0) < TOL, "coherent frame transport preserves purity")
    check(np.linalg.norm(p - expm(eps * k / 3.0)) > 1.0, "a completed jump is not the small K_or/3 connection")


def section_normalization_not_forced() -> None:
    print("\n[4] Representation theory does not fix the mean-cycle coefficient")
    k = k_or()
    candidates = {
        "raw edge q=1": 1.0,
        "unit-axis q=1/sqrt3": 1.0 / math.sqrt(3.0),
        "mean-cycle q=1/3": 1.0 / 3.0,
    }
    for name, q in candidates.items():
        sign_ok = is_sign_operator(q * k)
        axis = standard_axis_for_q(q) + q * np.ones(3)
        deriv = derivative_coefficients(axis)
        print(f"  {name:24s}: sign-covariant={sign_ok}, d12={deriv[0]:+.6f}, d23={deriv[1]:+.6f}, d13={deriv[2]:+.6f}")
        check(sign_ok, f"{name} has the correct sign-representation covariance")
    check(abs(derivative_coefficients(standard_axis_for_q(1.0 / 3.0) + np.ones(3) / 3.0)[0] + 1.0) < TOL, "only the mean-cycle value gives the QLC solar derivative")
    check(abs(derivative_coefficients(standard_axis_for_q(1.0) + np.ones(3))[0] + 3.0) < TOL, "raw sign current gives a different but equally equivariant normalization")


def section_conditional_lift() -> None:
    print("\n[5] Conditional frame-transport theorem")
    q = 1.0 / 3.0
    total_axis = standard_axis_for_q(q) + q * np.ones(3)
    deriv = derivative_coefficients(total_axis)
    u_frame = base_matrix() @ expm(DELTA * skew(total_axis))
    qlc = qlc_matrix()
    print(f"  total axis                  = {total_axis}")
    print(f"  derivatives/radian          = ({deriv[0]:+.6f},{deriv[1]:+.6f},{deriv[2]:+.6f})")
    print(f"  frame-lift finite angles    = {fmt(angles_deg(u_frame))}")
    print(f"  exact ordered QLC angles    = {fmt(angles_deg(qlc))}")
    check(np.linalg.norm(total_axis - np.array([-0.5, 0.5, 1.0])) < TOL, "mean-cycle sign plus standard shear gives the QLC tangent")
    check(abs(deriv[0] + 1.0) < TOL and abs(deriv[1]) < TOL and abs(deriv[2] - 1.0 / math.sqrt(2.0)) < TOL, "first-order QLC derivatives follow")
    check(np.linalg.norm(angles_deg(u_frame) - angles_deg(qlc)) < 1.1, "single-exponential lift remains close to finite ordered QLC")


def main() -> int:
    print("ITEM 87 -- PMNS EIGENVECTOR-LIFT THEOREM AUDIT")
    print("=" * 96)
    section_imported_canon()
    section_cp_pointer_not_lift()
    section_reset_not_lift()
    section_normalization_not_forced()
    section_conditional_lift()
    print(
        """
VERDICT
  Current canon does not yet prove the PMNS eigenvector-lift theorem.

  It proves that the oriented boundary can be read as a sign record, and the
  global orientation line can supply that sign pointer conditionally.  But a
  record readout is not automatically a real SO(3) frame update.  The Hermitian
  iK_or pointer stays in the DFT/circulant mass algebra, and a reset/jump
  channel is dissipative or a completed 120-degree permutation, not the small
  coherent connection exp(delta K_or/3).

  The live route is narrow and plausible:
    prove a coherent-recovery polar-factor lemma saying that lepton-sector
    recovery applies a pre-reset frame transport U_frame=exp(delta K_or/3),
    with one completed oriented boundary winding averaged over the three R1
    corners.  Under that lemma, the sign+standard bridge gives the QLC tangent.

exit 0 -- eigenvector-lift not closed; reduced to a coherent frame-transport plus mean-cycle normalization lemma.
"""
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
