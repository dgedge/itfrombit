#!/usr/bin/env python3
r"""ITEM 87 -- PMNS sign+standard eigenvector-rotation bridge gate.

Question
--------
Can the K_or/sign-pointer progress be driven into the missing PMNS
eigenvector-rotation primitive?

Result
------
Partly, but only with one clearly named new bridge.

The representation theorem already says that K_or is the oriented R1 boundary
cochain.  As a real antisymmetric matrix it is the sign-direction generator

    K_or = A(1,1,1).

It also supplies the Hodge orientation on the S3 standard plane:

    K_or^2 = -3 P_standard.

Thus an oriented R1 boundary can, in principle, convert a standard-generation
contrast into a real SO(3) eigenvector connection.  That is exactly the missing
kind of object: not a scalar support and not a Hermitian CP phase label.

The coefficient is the load-bearing part.  If a completed oriented boundary
service is billed as one circulation shared by the three generation corners,
the real sign generator is K_or/3, i.e. axis q(1,1,1) with q=1/3.  Combine that
with:

  * a standard-plane shear,
  * the universal delta=2/9 reactor primitive d theta_13 / d delta = 1/sqrt(2),
  * the atmospheric latch d theta_23 / d delta = 0.

Then the standard component is forced:

    s = (-5/6, 1/6, 2/3),

and the total tangent is

    s + (1/3,1/3,1/3) = (-1/2, 1/2, 1),

which is exactly the QLC tangent found in the previous audit.  The solar shift
d theta_12 / d delta = -1 is then a consequence, not an additional fit.

What is not yet proved by the older canon is the eigenvector-lift clause:
the physical sterile/charged-lepton readout must bill the oriented R1 boundary
as a real SO(3) rotation generator with mean-cycle normalization K_or/3.  A
Hermitian i K_or CP pointer does not do this.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.linalg import expm


DELTA = 2.0 / 9.0
TOL = 1.0e-10


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def skew(axis: np.ndarray) -> np.ndarray:
    x, y, z = axis
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])


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


def qlc_angles() -> np.ndarray:
    return np.array([45.0 - math.degrees(DELTA), 45.0, math.degrees(DELTA) / math.sqrt(2.0)])


def angles_deg(u: np.ndarray) -> np.ndarray:
    s13sq = abs(u[0, 2]) ** 2
    s12sq = abs(u[0, 1]) ** 2 / max(1.0 - s13sq, 1.0e-15)
    s23sq = abs(u[1, 2]) ** 2 / max(1.0 - s13sq, 1.0e-15)
    return np.array(
        [np.degrees(np.arcsin(np.sqrt(np.clip(x, 0.0, 1.0)))) for x in (s12sq, s23sq, s13sq)]
    )


def fmt(values: np.ndarray) -> tuple[float, float, float]:
    return tuple(float(round(v, 3)) for v in values)


def derivative_coefficients(axis: np.ndarray) -> np.ndarray:
    """First-order PMNS angle derivatives per connection radian at U0.

    The positive reactor branch is chosen.
    """

    x, y, z = axis
    return np.array([-z, -(x + y) / math.sqrt(2.0), (y - x) / math.sqrt(2.0)])


def standard_axis_for_q(q: float) -> np.ndarray:
    """Standard shear forced by reactor primitive and atmospheric latch.

    Conditions:
      x + y + z = 0,
      y - x = 1,
      total d theta_23 = 0 after adding q(1,1,1).
    """

    x = -q - 0.5
    y = 0.5 - q
    z = 2.0 * q
    return np.array([x, y, z])


def section_orientation_hodge() -> None:
    print("\n[1] Oriented R1 boundary as real Hodge operator")
    ones = np.ones(3)
    k_or = skew(ones)
    p_standard = np.eye(3) - np.outer(ones, ones) / 3.0
    print(f"  K_or =\n{k_or}")
    print(f"  ||K_or^2 + 3 P_standard|| = {np.linalg.norm(k_or @ k_or + 3.0 * p_standard):.3e}")
    check(np.linalg.norm(k_or @ ones) < TOL, "K_or kills the generation singlet")
    check(np.linalg.norm(k_or @ k_or + 3.0 * p_standard) < TOL, "K_or/sqrt(3) is a complex structure on the standard plane")
    check(np.linalg.matrix_rank(k_or) == 2, "K_or acts only on the two-dimensional standard plane")


def section_mean_cycle_bridge() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    print("\n[2] Mean-cycle sign bridge plus standard shear")
    q = 1.0 / 3.0
    sign_axis = q * np.ones(3)
    standard_axis = standard_axis_for_q(q)
    total_axis = sign_axis + standard_axis
    deriv = derivative_coefficients(total_axis)
    print(f"  sign axis K_or/3             = {sign_axis}")
    print(f"  forced standard axis          = {standard_axis}")
    print(f"  total bridge axis             = {total_axis}")
    print(f"  derivatives/radian (12,23,13) = ({deriv[0]:.6f},{deriv[1]:.6f},{deriv[2]:.6f})")
    check(abs(standard_axis.sum()) < TOL, "standard component is in x+y+z=0")
    check(abs(deriv[2] - 1.0 / math.sqrt(2.0)) < TOL, "universal 2/9 reactor primitive is retained")
    check(abs(deriv[1]) < TOL, "atmospheric latch is satisfied at first order")
    check(abs(deriv[0] + 1.0) < TOL, "solar QLC shift follows from q=1/3 plus the latch")
    check(np.linalg.norm(total_axis - np.array([-0.5, 0.5, 1.0])) < TOL, "bridge axis is exactly the QLC tangent axis")
    return sign_axis, standard_axis, total_axis


def section_normalization_audit() -> None:
    print("\n[3] Why the sign coefficient is the load-bearing primitive")
    candidates = {
        "raw K_or edge coefficient q=1": 1.0,
        "unit-axis coefficient q=1/sqrt(3)": 1.0 / math.sqrt(3.0),
        "mean-cycle coefficient q=1/3": 1.0 / 3.0,
    }
    for name, q in candidates.items():
        axis = standard_axis_for_q(q) + q * np.ones(3)
        deriv = derivative_coefficients(axis)
        print(f"  {name:36s}: d12={deriv[0]:+.6f}, d23={deriv[1]:+.6f}, d13={deriv[2]:+.6f}")
    check(abs(derivative_coefficients(standard_axis_for_q(1.0) + np.ones(3))[0] + 3.0) < TOL, "raw K_or gives a three-times solar shift")
    unit_q = 1.0 / math.sqrt(3.0)
    unit_axis = standard_axis_for_q(unit_q) + unit_q * np.ones(3)
    check(abs(derivative_coefficients(unit_axis)[0] + math.sqrt(3.0)) < TOL, "unit-norm sign gives a sqrt(3) solar shift")
    mean_axis = standard_axis_for_q(1.0 / 3.0) + np.ones(3) / 3.0
    check(abs(derivative_coefficients(mean_axis)[0] + 1.0) < TOL, "mean-cycle K_or/3 is the normalization that predicts QLC solar shift")


def section_finite_angles(total_axis: np.ndarray) -> None:
    print("\n[4] Finite delta=2/9 prediction")
    u_exp = base_matrix() @ expm(DELTA * skew(total_axis))
    qlc = qlc_matrix()
    print(f"  bridge exponential angles = {fmt(angles_deg(u_exp))}")
    print(f"  exact QLC product angles  = {fmt(angles_deg(qlc))}")
    print(f"  QLC target formula        = {fmt(qlc_angles())}")
    print(f"  angle mismatch norm       = {np.linalg.norm(angles_deg(u_exp) - qlc_angles()):.3f} deg")
    check(np.linalg.norm(angles_deg(u_exp) - qlc_angles()) < 1.1, "single-exponential bridge remains near the exact QLC product")
    check(np.linalg.norm(angles_deg(qlc) - qlc_angles()) < TOL, "exact QLC product is still a finite ordered-rotation construction")


def main() -> int:
    print("ITEM 87 -- PMNS SIGN+STANDARD BRIDGE GATE")
    print("=" * 92)
    section_orientation_hodge()
    _sign_axis, _standard_axis, total_axis = section_mean_cycle_bridge()
    section_normalization_audit()
    section_finite_angles(total_axis)
    print(
        """
VERDICT
  The sign+standard primitive can be driven to a clean conditional theorem.

  If the oriented R1 boundary K_or is lifted from a CP-orientation label to a
  real eigenvector-rotation current, and if one completed sign winding is
  normalized as K_or/3 across the three generation corners, then the universal
  2/9 standard shear plus the atmospheric latch uniquely gives the QLC tangent.
  The solar shift is then predicted rather than separately assigned.

  The remaining unproved step is exactly the eigenvector-lift/normalization
  theorem: current canon proves the sign pointer can read K_or, but it does not
  yet prove that this readout is billed as a real SO(3) PMNS rotation with
  mean-cycle normalization.

exit 0 -- sign+standard PMNS bridge conditionally driven; eigenvector-lift remains the live primitive.
"""
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
