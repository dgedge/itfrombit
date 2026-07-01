#!/usr/bin/env python3
r"""ITEM 87 -- can the PMNS coherent pre-reset polar factor be shown?

Question
--------
Can we show that lepton recovery has a coherent pre-reset polar factor

    U_frame = exp(delta K_or / 3)

rather than only a Hermitian CP pointer or a dissipative reset?

Answer
------
We can show a conditional uniqueness theorem, not an unconditional current-canon
proof.

Assume:

  A1. Before Landauer reset, lepton recovery has a coherent one-port
      Stinespring branch on the three R1 generation corners.
  A2. The only oriented generation datum available to that branch is the
      global-orientation-contracted R1 boundary cochain K_or.
  A3. The pre-reset branch is frame transport, so its polar unitary is generated
      by a real antisymmetric generator on the generation qutrit.
  A4. One completed oriented boundary winding is billed as one circulation
      shared over the three R1 corners.

Then:

  * S3 covariance forces the real antisymmetric generator to be q K_or.
  * A4 fixes q = 1/3.
  * The universal service primitive supplies delta = 2/9.
  * Hence the polar factor is uniquely exp(delta K_or / 3), up to the global
    orientation sign.

What remains open is A1/A3 as a physical lepton-recovery theorem.  The current
canon proves a sign record is readable; it does not prove that this record is
implemented as coherent frame transport before reset.
"""

from __future__ import annotations

from itertools import permutations
from pathlib import Path
import math

import numpy as np
from scipy.linalg import expm, polar


ROOT = Path(__file__).resolve().parents[1]
DELTA = 2.0 / 9.0
TOL = 1.0e-10


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


def sign_covariant_skew_basis() -> list[np.ndarray]:
    """Solve P^T G P = sgn(P) G for skew-symmetric real 3x3 matrices."""

    # Parametrize skew matrix by axis vector a=(x,y,z), i.e. G=skew(a).
    rows: list[np.ndarray] = []
    basis = [skew(np.eye(3)[i]) for i in range(3)]
    for _perm, p, parity in PERMS:
        moved = [p.T @ b @ p - parity * b for b in basis]
        for i in range(3):
            for j in range(3):
                rows.append([m[i, j] for m in moved])
    mat = np.array(rows, dtype=float)
    _u, s, vh = np.linalg.svd(mat)
    null = vh[np.sum(s > TOL) :]
    return [sum(coeff * b for coeff, b in zip(vec, basis)) for vec in null]


def circulation_per_corner(generator: np.ndarray) -> float:
    """Mean absolute oriented boundary current leaving each generation corner."""

    offdiag_abs = np.abs(generator.copy())
    np.fill_diagonal(offdiag_abs, 0.0)
    # For K_or, each corner has two signed incident entries of magnitude 1, but
    # one oriented unit of circulation passes the corner per full directed cycle.
    # The normalization below counts total directed circulation / three corners.
    total_directed_edges = float(np.sum(np.triu(offdiag_abs)))
    return total_directed_edges / 3.0


def angles_deg(u: np.ndarray) -> tuple[float, float, float]:
    s13sq = abs(u[0, 2]) ** 2
    s12sq = abs(u[0, 1]) ** 2 / max(1.0 - s13sq, 1.0e-15)
    s23sq = abs(u[1, 2]) ** 2 / max(1.0 - s13sq, 1.0e-15)
    return tuple(float(np.degrees(np.arcsin(np.sqrt(np.clip(x, 0.0, 1.0))))) for x in (s12sq, s23sq, s13sq))


def r12(theta: float) -> np.ndarray:
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])


def r23(theta: float) -> np.ndarray:
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]])


def base_matrix() -> np.ndarray:
    return r23(np.pi / 4.0) @ r12(np.pi / 4.0)


def derivative_coefficients(axis: np.ndarray) -> np.ndarray:
    x, y, z = axis
    return np.array([-z, -(x + y) / math.sqrt(2.0), (y - x) / math.sqrt(2.0)])


def standard_axis_for_q(q: float) -> np.ndarray:
    return np.array([-q - 0.5, 0.5 - q, 2.0 * q])


def section_current_canon_boundary() -> None:
    print("\n[1] Current-canon boundary")
    require_text(
        "python_code/r15_global_orientation_sign_pointer_theorem.py",
        (
            "omega*K_R1 is an S3-scalar environment record",
            "fixed orientation line adds no stochastic port count",
        ),
    )
    require_text(
        "python_code/item87_pmns_eigenvector_lift_theorem_audit.py",
        (
            "record readout is not automatically a real SO(3) frame update",
            "coherent-recovery polar-factor lemma",
        ),
    )
    require_text(
        "python_code/item87_pmns_eigenvector_lift_status.py",
        (
            "SOLE REMAINING PRIMITIVE: the eigenvector-LIFT theorem",
            "as a real SO(3) rotation",
        ),
    )
    check(True, "current canon has sign readout, but marks frame transport as the remaining primitive")


def section_covariant_generator_uniqueness() -> np.ndarray:
    print("\n[2] S3 covariance forces the sign generator")
    basis = sign_covariant_skew_basis()
    k = k_or()
    print(f"  dim sign-covariant skew generators = {len(basis)}")
    check(len(basis) == 1, "there is a unique sign-covariant real antisymmetric generator")
    b = basis[0]
    scale = np.vdot(k, b).real / np.vdot(k, k).real
    residual = np.linalg.norm(b - scale * k)
    print(f"  basis generator scale relative to K_or = {scale:.6f}; residual={residual:.3e}")
    check(residual < TOL, "the generator is proportional to K_or")
    return k


def section_mean_cycle_normalization(k: np.ndarray) -> np.ndarray:
    print("\n[3] Mean-cycle normalization fixes K_or/3 only after a service lemma")
    for label, q in (("raw edge", 1.0), ("unit-axis", 1.0 / math.sqrt(3.0)), ("mean-cycle", 1.0 / 3.0)):
        gen = q * k
        print(f"  {label:10s} q={q:.6f}: mean circulation/corner={circulation_per_corner(gen):.6f}")
    q = 1.0 / 3.0
    check(abs(circulation_per_corner(q * k) - 1.0 / 3.0) < TOL, "K_or/3 bills one total oriented circulation over three corners")
    check(abs(circulation_per_corner(k) - 1.0) < TOL, "raw K_or bills one circulation at each corner, three times the mean-cycle PMNS value")
    return q * k


def section_polar_factor(k_mean: np.ndarray) -> None:
    print("\n[4] Conditional coherent pre-reset polar factor")
    u_target = expm(DELTA * k_mean)
    positive_record = np.diag([0.97, 1.01, 1.03])
    kraus = u_target @ positive_record
    u_polar, p_polar = polar(kraus)
    print(f"  ||U_polar - exp(delta K_or/3)|| = {np.linalg.norm(u_polar - u_target):.3e}")
    print(f"  ||P_polar - positive_record||   = {np.linalg.norm(p_polar - positive_record):.3e}")
    check(np.linalg.norm(u_polar - u_target) < TOL, "if the pre-reset branch has this coherent factor, polar decomposition recovers it exactly")
    check(np.linalg.norm(u_polar.T @ u_polar - np.eye(3)) < TOL, "the polar factor is real orthogonal")
    check(abs(np.linalg.det(u_polar) - 1.0) < TOL, "the polar factor is in SO(3)")


def section_angle_consequence(k_mean: np.ndarray) -> None:
    print("\n[5] Consequence for the PMNS tangent")
    q = 1.0 / 3.0
    total_axis = standard_axis_for_q(q) + q * np.ones(3)
    d = derivative_coefficients(total_axis)
    u = base_matrix() @ expm(DELTA * skew(total_axis))
    print(f"  total QLC axis from mean-cycle sign + standard shear = {total_axis}")
    print(f"  derivatives/radian (12,23,13) = ({d[0]:+.6f},{d[1]:+.6f},{d[2]:+.6f})")
    print(f"  finite single-exponential angles = {tuple(round(x, 3) for x in angles_deg(u))}")
    check(np.linalg.norm(k_mean - k_or() / 3.0) < TOL, "mean sign current is K_or/3")
    check(abs(d[0] + 1.0) < TOL and abs(d[1]) < TOL and abs(d[2] - 1.0 / math.sqrt(2.0)) < TOL, "QLC tangent follows under the coherent-lift lemma")


def main() -> int:
    print("ITEM 87 -- COHERENT PRE-RESET POLAR FACTOR ATTEMPT")
    print("=" * 92)
    section_current_canon_boundary()
    k = section_covariant_generator_uniqueness()
    k_mean = section_mean_cycle_normalization(k)
    section_polar_factor(k_mean)
    section_angle_consequence(k_mean)
    print(
        """
VERDICT
  We can show the theorem conditionally and uniquely:
    if lepton recovery has a coherent one-port pre-reset frame-transport branch,
    S3 covariance and the global orientation line force its real antisymmetric
    generator to be proportional to K_or; the mean-cycle service normalization
    fixes the proportionality to K_or/3; delta=2/9 then gives
        U_frame = exp(delta K_or / 3).

  We cannot honestly mark this as an unconditional current-canon proof.  The
  current canon proves sign readout and one coherent source port, but it does
  not yet prove that the sign record is implemented as a coherent pre-reset
  frame-transport polar factor.  That remains the exact physical lemma to add.

exit 0 -- conditional uniqueness shown; existence of the coherent lepton-recovery frame branch remains open.
"""
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
