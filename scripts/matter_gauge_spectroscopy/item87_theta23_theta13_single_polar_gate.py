#!/usr/bin/env python3
r"""Item 87 -- theta13 and theta23 octant from the single-polar PMNS transport.

Question
--------
Can the framework turn the old "theta23 second-octant, conditional on
frame-transport closure" bracket into a clean prediction?

Result
------
Yes, at the current item-87 tier, but only under the physically specific
single-polar reading of the QEC correction.

The current canon has three ingredients:

  1. Dynamic R1 / Hasse-edge recovery supplies the oriented generation-boundary
     record and the sector contact magnitude delta = 2/9 for the charged-lepton
     / PMNS correction scale.
  2. The frame-transport lemma says the PMNS correction is the coherent
     pre-reset QEC correction unitary, i.e. one real SO(3) polar factor, not a
     measurement collapse or a post-hoc Euler-angle report.
  3. The mean-cycle sign component is K_or/3.  Combining that with the standard
     plane shear and the atmospheric latch fixes the bridge tangent

         a = (-1/2, 1/2, 1).

Therefore the physical leading PMNS correction is

    U = U0 exp[delta A(a)],     U0 = R23(pi/4) R12(pi/4),     delta = 2/9.

This gives theta13 = 8.93 deg and theta23 = 46.01 deg.  If the finite
correction scale is instead normalized to the observed reactor angle theta13 =
8.6 deg, the same one-parameter branch gives theta23 = 45.93 deg.  The octant
is second, but near-maximal.

Important distinction
---------------------
The often-useful Euler shorthand

    R23(pi/4) R13(delta/sqrt(2)) R12(pi/4-delta)

has theta23 = 45 deg exactly.  It is a real angle-reporting representative and
is good enough for the Jarlskog=0 / Dirac-CP statement, but it is not a single
QEC correction polar factor.  Once the frame-transport lemma is read as a
single coherent recovery apply, the single exponential is the physical branch
and the second-octant lift is forced by the noncommutativity of the correction.

Exit 0 means the leading-order PMNS prediction is:

    theta13 ~= 8.9 deg, theta23 ~= 46.0 deg (second octant, near maximal),
    J_lepton = 0.

It does not force theta23 ~= 49 deg; reaching that value requires breaking the
atmospheric latch or adding an unfixed charged-lepton/Dirac-sector correction.
"""

from __future__ import annotations

import math

import numpy as np


DELTA = 2.0 / 9.0
AXIS = np.array([-0.5, 0.5, 1.0])
TOL = 1.0e-10


def check(message: str, condition: bool) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def skew(axis: np.ndarray) -> np.ndarray:
    x, y, z = axis
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])


def expm_so3(matrix: np.ndarray) -> np.ndarray:
    """Rodrigues exponential for a real antisymmetric 3x3 matrix."""

    axis = np.array([matrix[2, 1], matrix[0, 2], matrix[1, 0]])
    theta = float(np.linalg.norm(axis))
    if theta < 1.0e-15:
        return np.eye(3)
    k = matrix / theta
    return np.eye(3) + math.sin(theta) * k + (1.0 - math.cos(theta)) * (k @ k)


def r12(theta: float) -> np.ndarray:
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])


def r23(theta: float) -> np.ndarray:
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]])


def r13(theta: float) -> np.ndarray:
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


def pmns_angles_deg(matrix: np.ndarray) -> tuple[float, float, float]:
    u = np.abs(matrix)
    theta13 = math.degrees(math.asin(min(1.0, float(u[0, 2]))))
    theta12 = math.degrees(math.atan2(float(u[0, 1]), float(u[0, 0])))
    theta23 = math.degrees(math.atan2(float(u[1, 2]), float(u[2, 2])))
    return theta12, theta23, theta13


def base_pmns() -> np.ndarray:
    return r23(math.pi / 4.0) @ r12(math.pi / 4.0)


def single_polar(delta: float = DELTA, axis: np.ndarray = AXIS) -> np.ndarray:
    return base_pmns() @ expm_so3(delta * skew(axis))


def euler_shorthand(delta: float = DELTA) -> np.ndarray:
    return r23(math.pi / 4.0) @ r13(delta / math.sqrt(2.0)) @ r12(math.pi / 4.0 - delta)


def jarlskog(matrix: np.ndarray) -> float:
    u = matrix.astype(complex)
    return float(np.imag(u[0, 0] * u[1, 1] * np.conj(u[0, 1]) * np.conj(u[1, 0])))


def fit_delta_to_theta13(target_deg: float) -> tuple[float, tuple[float, float, float]]:
    best: tuple[float, float, tuple[float, float, float]] | None = None
    for delta in np.linspace(0.0, 0.35, 35001):
        angles = pmns_angles_deg(single_polar(float(delta)))
        error = abs(angles[2] - target_deg)
        if best is None or error < best[0]:
            best = (error, float(delta), angles)
    assert best is not None
    return best[1], best[2]


def derivative_coefficients(axis: np.ndarray) -> np.ndarray:
    x, y, z = axis
    return np.array([-z, -(x + y) / math.sqrt(2.0), (y - x) / math.sqrt(2.0)])


def section_forced_axis() -> None:
    print("\n[1] The sign+standard bridge fixes the tangent")
    sign_component = float(np.mean(AXIS))
    derivative = derivative_coefficients(AXIS)
    print(f"    axis a = {AXIS}")
    print(f"    mean/sign component = {sign_component:.12f}")
    print(
        "    first-order angle derivatives per radian "
        f"(theta12, theta23, theta13) = ({derivative[0]:+.6f}, {derivative[1]:+.6f}, {derivative[2]:+.6f})"
    )
    check("mean/sign component is K_or/3 = 1/3", abs(sign_component - 1.0 / 3.0) < TOL)
    check("atmospheric latch: first-order d theta23 = 0", abs(derivative[1]) < TOL)
    check("reactor primitive: d theta13 = 1/sqrt(2)", abs(derivative[2] - 1.0 / math.sqrt(2.0)) < TOL)
    check("solar QLC derivative: d theta12 = -1", abs(derivative[0] + 1.0) < TOL)


def section_single_polar_prediction() -> None:
    print("\n[2] Single coherent QEC correction: U0 exp(delta A(a))")
    angles = pmns_angles_deg(single_polar())
    print(
        "    delta = 2/9 gives "
        f"(theta12, theta23, theta13) = ({angles[0]:.3f}, {angles[1]:.3f}, {angles[2]:.3f}) deg"
    )
    check("theta13 lands in the reactor-angle window", 8.3 < angles[2] < 9.2)
    check("theta23 is in the second octant", angles[1] > 45.0)
    check("theta23 is near-maximal, not a 49-degree branch", angles[1] < 46.5)
    check("single-polar PMNS is real and orthogonal", np.allclose(single_polar() @ single_polar().T, np.eye(3)))
    check("there is no leptonic Dirac CP in the leading real transport", abs(jarlskog(single_polar())) < TOL)

    fitted_delta, fitted_angles = fit_delta_to_theta13(8.6)
    print(
        "    if normalized to theta13=8.6 deg: "
        f"delta={fitted_delta:.5f}, theta23={fitted_angles[1]:.3f} deg"
    )
    check("reactor-normalized branch remains second-octant", fitted_angles[1] > 45.0)
    check("reactor-normalized branch remains near-maximal", fitted_angles[1] < 46.2)


def section_controls() -> None:
    print("\n[3] Controls: product shorthand and latch-breaking routes")
    shorthand_angles = pmns_angles_deg(euler_shorthand())
    print(
        "    Euler shorthand gives "
        f"(theta12, theta23, theta13) = ({shorthand_angles[0]:.3f}, {shorthand_angles[1]:.3f}, {shorthand_angles[2]:.3f}) deg"
    )
    check("Euler shorthand is real, so it is sufficient for J=0", abs(jarlskog(euler_shorthand())) < TOL)
    check("Euler shorthand does not force an octant", abs(shorthand_angles[1] - 45.0) < TOL)

    negative_angles = pmns_angles_deg(single_polar(axis=-AXIS))
    print(
        "    reversed correction axis gives "
        f"(theta12, theta23, theta13) = ({negative_angles[0]:.3f}, {negative_angles[1]:.3f}, {negative_angles[2]:.3f}) deg"
    )
    check("reversed-axis branch is excluded by the solar angle", negative_angles[0] > 50.0)

    # Reproduce the latch-breaking cost: hold reactor and solar levers but demand theta23 ~= 49 deg.
    delta_fit, _ = fit_delta_to_theta13(8.6)
    needed_d23 = math.radians(49.0 - 45.0) / delta_fit
    xy_sum = -needed_d23 * math.sqrt(2.0)
    y = (xy_sum + 1.0) / 2.0
    x = y - 1.0
    latch_broken = np.array([x, y, 1.0])
    broken_sign = float(np.mean(latch_broken))
    print(
        "    theta23=49 while keeping reactor+solar levers needs "
        f"a'=({x:.3f}, {y:.3f}, 1.000), mean={broken_sign:.3f}"
    )
    check("49-degree branch breaks K_or/3 mean-cycle normalization", abs(broken_sign - 1.0 / 3.0) > 0.1)


def section_branch_verdict() -> None:
    print("\n[4] Status")
    print(
        """
    The octant becomes forced only after choosing the physical frame-transport
    object.  A QEC correction is one coherent pre-reset apply, so the
    single-polar branch U0 exp(delta A(a)) is the natural reading of the closed
    frame-transport lemma.  On that branch theta13 and the second-octant,
    near-maximal theta23 follow from the same fixed tangent and delta=2/9.

    If one rejects the single-polar reading and treats the PMNS texture only as
    an Euler-angle shorthand, theta13 remains reproduced but theta23 stays
    exactly maximal and the octant is not predicted.  Thus the remaining
    caveat is not a free angle: it is the physical identification of the QEC
    correction with a single coherent polar factor, which current canon already
    favours.
"""
    )


def main() -> int:
    print("ITEM 87 -- theta13 + theta23 octant from single-polar PMNS transport")
    print("=" * 88)
    section_forced_axis()
    section_single_polar_prediction()
    section_controls()
    section_branch_verdict()
    print("ALL ASSERTIONS PASSED")
    print(
        "\nVERDICT: conditionally forced prediction -- theta13 ~= 8.9 deg and "
        "theta23 ~= 46.0 deg, second octant, near maximal.  No in-framework route "
        "to theta23 ~= 49 deg without breaking the atmospheric latch."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
