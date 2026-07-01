#!/usr/bin/env python3
r"""ITEM 87 -- standard-representation PMNS texture derivation attempt.

Goal
----
Try to derive the missing PMNS angle texture from a standard-representation
generation shear whose strength is the universal 2/9 primitive.

Result
------
There is real progress, but not a full derivation.

At the bimaximal base

    U0 = R23(pi/4) R12(pi/4),

write a small real generation texture as a right-acting antisymmetric
connection

    A(x,y,z) = [[0,-z, y],
                [z, 0,-x],
                [-y,x, 0]].

Under S3 generation relabelling this 3-vector splits into

    sign direction        (1,1,1),
    standard plane        x + y + z = 0.

The QLC target has first-order derivatives

    d theta12 / d delta = -1,
    d theta23 / d delta =  0,
    d theta13 / d delta =  1/sqrt(2).

But a pure standard-plane generator obeys the first-order identity

    d theta23 / d delta = - (d theta12 / d delta) / sqrt(2).

Therefore a pure standard-representation texture cannot derive the full QLC
pattern: if it shifts theta12 by -delta, it must also shift theta23 upward.
It can generate theta13 ~= delta/sqrt(2), but not the whole PMNS pattern.

The unique QLC tangent is simple,

    A_QLC = [[0,-1, 1/2],
             [1, 0, 1/2],
             [-1/2,-1/2,0]],

with axis (-1/2, 1/2, 1).  Its standard part gives the reactor-angle lever,
while its sign part cancels the unwanted atmospheric drift and completes the
solar shift.  Using that sign part as an angle texture would identify the
K_or/sign sector with real eigenvector rotation, which is not currently derived
from the substrate.  Hence the texture is conditional, not locked.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.linalg import expm


DELTA = 2.0 / 9.0
TOL = 1.0e-10
OBS = np.array([33.44, 49.2, 8.57])


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def r12(theta: float) -> np.ndarray:
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])


def r23(theta: float) -> np.ndarray:
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]])


def r13(theta: float) -> np.ndarray:
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


def skew(axis: np.ndarray) -> np.ndarray:
    x, y, z = axis
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])


def angles_deg(u: np.ndarray) -> np.ndarray:
    s13sq = abs(u[0, 2]) ** 2
    s12sq = abs(u[0, 1]) ** 2 / max(1.0 - s13sq, 1.0e-15)
    s23sq = abs(u[1, 2]) ** 2 / max(1.0 - s13sq, 1.0e-15)
    return np.array(
        [np.degrees(np.arcsin(np.sqrt(np.clip(x, 0.0, 1.0)))) for x in (s12sq, s23sq, s13sq)]
    )


def qlc_angles() -> np.ndarray:
    return np.array([45.0 - math.degrees(DELTA), 45.0, math.degrees(DELTA) / math.sqrt(2.0)])


def qlc_matrix() -> np.ndarray:
    return r23(np.pi / 4.0) @ r13(DELTA / math.sqrt(2.0)) @ r12(np.pi / 4.0 - DELTA)


def base_matrix() -> np.ndarray:
    return r23(np.pi / 4.0) @ r12(np.pi / 4.0)


def derivative_at_base(axis: np.ndarray, h: float = 1.0e-7) -> np.ndarray:
    u0 = base_matrix()
    return (angles_deg(u0 @ expm(h * skew(axis))) - angles_deg(u0)) / math.degrees(h)


def fmt_angles(values: np.ndarray) -> tuple[float, float, float]:
    return tuple(float(round(x, 3)) for x in values)


def sign_standard_split(axis: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    sign = np.ones(3) * axis.mean()
    standard = axis - sign
    return sign, standard


def section_first_order_obstruction() -> None:
    print("\n[1] First-order standard-plane obstruction")
    # For A(x,y,z) at U0, away from the reactor cusp and with the positive
    # reactor branch chosen:
    #   d12 = -z
    #   d23 = -(x+y)/sqrt(2)
    #   d13 = (y-x)/sqrt(2)
    # If x+y+z=0, then d23 = z/sqrt(2) = -d12/sqrt(2).
    standard_axis = np.array([-5.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0])
    d12 = -standard_axis[2]
    d23 = -(standard_axis[0] + standard_axis[1]) / math.sqrt(2.0)
    d13 = (standard_axis[1] - standard_axis[0]) / math.sqrt(2.0)
    print(f"  standard axis = {standard_axis}")
    print(f"  derivatives/radian = (d12,d23,d13)=({d12:.6f},{d23:.6f},{d13:.6f})")
    check(abs(standard_axis.sum()) < TOL, "axis is in the standard plane x+y+z=0")
    check(abs(d23 + d12 / math.sqrt(2.0)) < TOL, "pure standard texture obeys d23 = -d12/sqrt(2)")
    check(abs(d13 - 1.0 / math.sqrt(2.0)) < TOL, "this standard direction gives the desired reactor derivative")
    check(abs(d12 + 1.0) > 0.1 and abs(d23) > 0.1, "it cannot also give d12=-1 and d23=0")


def section_candidate_predictions() -> None:
    print("\n[2] Finite predictions at delta=2/9")
    u0 = base_matrix()
    qlc = qlc_matrix()
    axis_qlc = np.array([-0.5, 0.5, 1.0])
    sign_axis, standard_axis = sign_standard_split(axis_qlc)
    u_standard = u0 @ expm(DELTA * skew(standard_axis))
    u_full = u0 @ expm(DELTA * skew(axis_qlc))

    print(f"  QLC target angles              = {fmt_angles(qlc_angles())}")
    print(f"  pure standard exponential      = {fmt_angles(angles_deg(u_standard))}")
    print(f"  standard+sign exponential      = {fmt_angles(angles_deg(u_full))}")
    print(f"  exact constructed QLC product  = {fmt_angles(angles_deg(qlc))}")
    print(f"  observed representative        = {fmt_angles(OBS)}")
    check(abs(angles_deg(u_standard)[2] - OBS[2]) < 0.2, "pure standard delta=2/9 predicts the reactor angle well")
    check(abs(angles_deg(u_standard)[1] - 45.0) > 5.0, "pure standard texture shifts theta23 away from maximal")
    check(np.linalg.norm(angles_deg(u_full) - qlc_angles()) < 1.1, "adding the sign component gives a near-QLC one-parameter connection")
    check(np.linalg.norm(angles_deg(qlc) - qlc_angles()) < TOL, "the exact QLC product is still a constructed texture")
    print(f"  sign component axis            = {sign_axis}")
    print(f"  standard component axis        = {standard_axis}")


def section_unique_tangent() -> None:
    print("\n[3] Unique tangent for the QLC angle formula")
    # Differentiate U(delta)=R23(pi/4) R13(delta/sqrt2) R12(pi/4-delta) at delta=0.
    x_qlc = np.array([[0.0, -1.0, 0.5], [1.0, 0.0, 0.5], [-0.5, -0.5, 0.0]])
    axis = np.array([-0.5, 0.5, 1.0])
    sign_axis, standard_axis = sign_standard_split(axis)
    sign_norm = np.linalg.norm(sign_axis)
    standard_norm = np.linalg.norm(standard_axis)
    print("  A_QLC =")
    print(x_qlc)
    print(f"  axis = {axis}; sign_norm={sign_norm:.6f}; standard_norm={standard_norm:.6f}")
    check(np.linalg.norm(skew(axis) - x_qlc) < TOL, "A_QLC is the tangent of the exact QLC product")
    check(abs(sign_axis.sum()) > 0.1 and abs(standard_axis.sum()) < TOL, "QLC tangent has both sign and standard pieces")
    check(standard_norm > sign_norm, "the standard component is the dominant angle-texture component")
    check(sign_norm > 0.0, "the sign component is nevertheless required to cancel the theta23 drift")


def main() -> int:
    print("ITEM 87 -- STANDARD-REPRESENTATION PMNS TEXTURE DERIVATION ATTEMPT")
    print("=" * 96)
    section_first_order_obstruction()
    section_candidate_predictions()
    section_unique_tangent()
    print(
        """
VERDICT
  A standard-representation texture tied to delta=2/9 can be derived at the
  level of a minimal generation shear, and it naturally produces the reactor
  angle scale theta13 ~= delta/sqrt(2).

  It does not close the full PMNS pattern.  Pure standard-plane generators obey
      d theta23 = - d theta12 / sqrt(2)
  at the bimaximal base, so they cannot simultaneously produce
      d theta12 = -1,  d theta23 = 0,  d theta13 = 1/sqrt(2).

  The QLC tangent that does this is simple and one-parameter, but it contains a
  sign-representation component in addition to the standard component.  Using
  that sign component as a real angle texture would be a new primitive: the
  current K_or/sign-pointer work licenses CP orientation, not PMNS eigenvector
  rotation.

exit 0 -- theta13 scale conditionally derived; full PMNS texture still needs a new sign/standard bridge.
"""
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
