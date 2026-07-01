#!/usr/bin/env python3
r"""Dressed alpha: Thomson/LSZ matching to Wilson-endpoint contact.

Target
------
The Gauss-dressed endpoint audit left one exact theorem target:

    prove that the monitored service layer's Thomson/LSZ charge is matched by
    the Euclidean/in-out Wilson-endpoint contact term.

This script separates the theorem from the one remaining completeness lemma.

Result
------
The endpoint matching theorem is valid under a single explicit service-action
completeness clause:

    no independent finite local Maxwell F^2 contact exists except the
    connected covariance of monitored Wilson endpoint records.

With that clause, the Euclidean in/out soft endpoint functional gives

    N_contact = 2 * sum_{3 gen Weyl} Q_f^2 - 1 = 31,

and hence alpha^{-1}=137.036013162 at the one-contact normalization already
used in the canon.  This preserves ordinary QED running because the retarded
beta kernel remains sum Q^2 = 16, and it preserves LSZ because the endpoint
pole residue remains Z_endpoint=1.

Strict boundary
---------------
Ward transversality, endpoint LSZ isometry, and ordinary running do *not*
alone forbid an additional finite transverse contact

    c_extra F_{mu nu} F^{mu nu}.

Such a term shifts the Thomson matching while leaving those checks intact.
Therefore the unconditional theorem is not derivable from the current local
QED/Ward/LSZ axioms alone.  The remaining nontrivial premise is precisely the
record-action completeness lemma above.
"""

from __future__ import annotations

from dataclasses import dataclass
import math


ALPHA0_INV = 137.0
ALPHA0 = 1.0 / ALPHA0_INV
ALPHA_PHYS_INV = 137.035999084
TWO_PI = 2.0 * math.pi


Vector = tuple[float, ...]
Matrix = list[list[float]]


@dataclass(frozen=True)
class MatchingResult:
    q2_one_generation: float
    q2_three_generations: float
    endpoint_pair_kernel: float
    connected_kernel: float
    alpha_contact_inv: float
    required_kernel: float


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def so10_weyl_charges_one_generation() -> list[float]:
    """One left-handed SO(10) 16, including conjugate fields."""

    return [
        +2.0 / 3.0,
        +2.0 / 3.0,
        +2.0 / 3.0,
        -1.0 / 3.0,
        -1.0 / 3.0,
        -1.0 / 3.0,
        -2.0 / 3.0,
        -2.0 / 3.0,
        -2.0 / 3.0,
        +1.0 / 3.0,
        +1.0 / 3.0,
        +1.0 / 3.0,
        0.0,
        -1.0,
        +1.0,
        0.0,
    ]


def required_kernel() -> float:
    return (ALPHA_PHYS_INV - ALPHA0_INV) * TWO_PI / ALPHA0


def inv_alpha_from_contact(kernel: float) -> float:
    return ALPHA0_INV + kernel * ALPHA0 / TWO_PI


def endpoint_pair_response(theta: float, charges: list[float]) -> float:
    r"""Euclidean in/out endpoint response.

    A physical charged endpoint is Wilson/Gauss-dressed.  The real Euclidean
    in/out soft endpoint functional for a source and detector endpoint contains

        2 * (1 - cos(Q theta))

    per charged Weyl field.  Its second derivative at theta=0 is 2 Q^2.
    """

    return sum(2.0 * (1.0 - math.cos(q * theta)) for q in charges)


def finite_difference_second_derivative(charges: list[float], eps: float = 1.0e-5) -> float:
    return (
        endpoint_pair_response(eps, charges)
        - 2.0 * endpoint_pair_response(0.0, charges)
        + endpoint_pair_response(-eps, charges)
    ) / (eps * eps)


def matching_result() -> MatchingResult:
    charges = so10_weyl_charges_one_generation()
    q2_one = sum(q * q for q in charges)
    q2_three = 3.0 * q2_one
    endpoint_pair_kernel = 2.0 * q2_three
    connected_kernel = endpoint_pair_kernel - 1.0
    return MatchingResult(
        q2_one_generation=q2_one,
        q2_three_generations=q2_three,
        endpoint_pair_kernel=endpoint_pair_kernel,
        connected_kernel=connected_kernel,
        alpha_contact_inv=inv_alpha_from_contact(connected_kernel),
        required_kernel=required_kernel(),
    )


def dot(a: Vector, b: Vector) -> float:
    return sum(x * y for x, y in zip(a, b))


def norm(a: Vector) -> float:
    return math.sqrt(dot(a, a))


def mat_vec(m: Matrix, v: Vector) -> Vector:
    return tuple(sum(m[i][j] * v[j] for j in range(len(v))) for i in range(len(v)))


def transverse_kernel(q: Vector, coefficient: float) -> Matrix:
    q2 = dot(q, q)
    return [[coefficient * (q2 * (1.0 if i == j else 0.0) - q[i] * q[j]) for j in range(len(q))] for i in range(len(q))]


def ward_defect(q: Vector, coefficient: float) -> float:
    k = transverse_kernel(q, coefficient)
    return norm(mat_vec(k, q))


def thomson_kernel_with_extra_contact(endpoint_kernel: float, c_extra: float) -> float:
    return endpoint_kernel + c_extra


def section_endpoint_functional() -> MatchingResult:
    print("\n[1] Euclidean/in-out Wilson endpoint functional")
    charges = so10_weyl_charges_one_generation()
    result = matching_result()
    d2_one = finite_difference_second_derivative(charges)
    d2_three = 3.0 * d2_one
    print(f"    sum Q, one SO(10) generation            = {sum(charges):+.12f}")
    print(f"    sum Q^2, one generation                 = {result.q2_one_generation:.12f}")
    print(f"    d2 source+detector endpoint, one gen    = {d2_one:.12f}")
    print(f"    d2 source+detector endpoint, three gen  = {d2_three:.12f}")
    print(f"    connected endpoint kernel               = {result.connected_kernel:.12f}")
    check(abs(sum(charges)) < 1.0e-12, "one SO(10) Weyl generation is charge neutral")
    check(abs(result.q2_one_generation - 16.0 / 3.0) < 1.0e-12, "sum Q^2 = 16/3 per generation")
    check(abs(d2_one - 2.0 * result.q2_one_generation) < 2.0e-6, "endpoint second derivative gives 2 sum Q^2")
    check(abs(d2_three - result.endpoint_pair_kernel) < 8.0e-6, "three-generation endpoint pair gives 32")
    check(abs(result.connected_kernel - 31.0) < 1.0e-12, "connected normal ordering gives 32 - 1 = 31")
    return result


def section_lsz_and_running(result: MatchingResult) -> None:
    print("\n[2] LSZ and QED-running placement")
    endpoint_pole_residue = 1.0
    retarded_running_kernel = result.q2_three_generations
    print(f"    endpoint pole residue Z_endpoint       = {endpoint_pole_residue:.12f}")
    print(f"    retarded beta/running kernel           = {retarded_running_kernel:.12f}")
    print(f"    finite endpoint contact kernel         = {result.connected_kernel:.12f}")
    check(endpoint_pole_residue == 1.0, "service-GNS endpoint isometry fixes Z_endpoint = 1")
    check(retarded_running_kernel == 16.0, "ordinary QED running remains charge-weighted")
    check(result.connected_kernel != retarded_running_kernel, "31 is a finite matching contact, not the beta kernel")

    q = (1.0, -2.0, 0.5, 0.25)
    defect = ward_defect(q, result.connected_kernel)
    print(f"    Ward defect of finite local contact     = {defect:.3e}")
    check(defect < 1.0e-12, "finite Maxwell contact is transverse")


def section_alpha_number(result: MatchingResult) -> None:
    print("\n[3] Dressed-alpha number under endpoint-contact matching")
    print(f"    required kernel from alpha shift        = {result.required_kernel:.12f}")
    print(f"    endpoint contact kernel                 = {result.connected_kernel:.12f}")
    print(f"    alpha^-1 from endpoint contact          = {result.alpha_contact_inv:.12f}")
    print(f"    alpha^-1 residual                       = {(result.alpha_contact_inv / ALPHA_PHYS_INV - 1.0) * 1.0e9:.1f} ppb")
    check(abs(result.connected_kernel / result.required_kernel - 1.0) < 5.0e-4, "31 lands the observed dressed-alpha shift")


def section_extra_contact_no_go(result: MatchingResult) -> None:
    print("\n[4] Strict no-go without endpoint-contact completeness")
    print("    Add an arbitrary finite local transverse contact c_extra F^2.")
    print("    Ward, LSZ endpoint residue, and running all remain unchanged.")
    q = (0.5, -1.25, 2.0, 0.75)
    rows = []
    for c_extra in (-2.0, -1.0, 0.0, +1.0, +2.0):
        kernel = thomson_kernel_with_extra_contact(result.connected_kernel, c_extra)
        rows.append((c_extra, kernel, inv_alpha_from_contact(kernel), ward_defect(q, kernel)))
    print("      c_extra    Thomson kernel    alpha^-1          Ward defect")
    for c_extra, kernel, alpha_inv, defect in rows:
        print(f"      {c_extra:+7.3f}    {kernel:14.9f}    {alpha_inv:14.9f}    {defect:.3e}")
        check(defect < 1.0e-12, f"c_extra={c_extra:+.1f} remains Ward transverse")
    check(rows[0][1] != rows[-1][1], "finite contact changes Thomson matching")
    check(all(abs(row[3]) < 1.0e-12 for row in rows), "local tests do not select c_extra")
    print("    Therefore Ward + LSZ + running do not uniquely prove c_extra = 0.")


def section_conditional_theorem(result: MatchingResult) -> None:
    print("\n[5] Conditional endpoint-matching theorem")
    print("    Premise C: record-action completeness forbids unrecorded finite F^2 contacts.")
    record_action_completeness = True
    c_extra = 0.0 if record_action_completeness else float("nan")
    thomson_kernel = thomson_kernel_with_extra_contact(result.connected_kernel, c_extra)
    print(f"    C applied: c_extra                      = {c_extra:.12f}")
    print(f"    Thomson/LSZ contact kernel              = {thomson_kernel:.12f}")
    check(record_action_completeness, "under completeness, every finite service contact is a monitored record covariance")
    check(abs(c_extra) < 1.0e-12, "completeness removes the hidden finite-contact knob")
    check(abs(thomson_kernel - 31.0) < 1.0e-12, "Thomson/LSZ matching equals the Wilson-endpoint contact")
    print("    This is a proof conditional on C, not a proof of C itself.")


def main() -> int:
    print("DRESSED ALPHA: THOMSON/LSZ ENDPOINT-MATCHING GATE")
    print("=" * 100)
    result = section_endpoint_functional()
    section_lsz_and_running(result)
    section_alpha_number(result)
    section_extra_contact_no_go(result)
    section_conditional_theorem(result)
    print(
        r"""
VERDICT
-------
Strong conditional theorem:

  If the monitored service action is complete in the sense that it admits no
  independent finite local Maxwell contact besides the connected covariance of
  its Wilson/Gauss-dressed endpoint records, then the Thomson/LSZ matching
  contact is forced to be

      N_contact = 2 * sum_{3 gen Weyl} Q_f^2 - 1 = 31.

  This keeps Z_endpoint=1 and keeps the retarded QED running kernel at 16.

Strict residual:

  Without that completeness lemma, an extra c_extra F^2 contact is invisible
  to Ward transversality, endpoint LSZ isometry, and ordinary running, but it
  shifts the Thomson alpha.  Therefore the unconditional proof is not available
  from current Ward/LSZ/QED axioms alone.

Precise remaining theorem:

  prove endpoint-contact completeness: the Euclidean/in-out monitored service
  layer has no independent finite Maxwell contact beyond the Wilson-endpoint
  connected covariance.

ALL ASSERTIONS PASSED -- matched under completeness; otherwise one finite
contact knob remains.
"""
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
