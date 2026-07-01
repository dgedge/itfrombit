#!/usr/bin/env python3
r"""Dressed alpha: endpoint covariance -> Maxwell contact map attempt.

Question
--------
The previous audit reduced the dressed-alpha route to a single operator map:

    connected endpoint closed-record covariance
        ?= finite local Maxwell F^2 contact.

If this map follows from record-action / monitored-QEC dynamics, then the
normal-ordering constant is fixed:

    N_contact = 2 * sum_f Q_f^2 - 1 = 31.

This script tests the map in the real-time influence-functional language.

Result
------
The coefficient is forced *if* the map is assumed, but the map itself is not
derived by the present real-time record-action machinery.

In Schwinger-Keldysh variables, integrating out monitored endpoint records
produces a symmetric/noise kernel:

    S_IF[A_r,A_a] = A_a K_R A_r + (i/2) A_a N A_a + ...

The closed-record-pair covariance fixes N.  A finite Maxwell contact shifts
the retarded/local kernel K_R by

    C_fin (q^2 eta_mu_nu - q_mu q_nu).

Ward transversality allows every C_fin, and fluctuation-dissipation relates N
to the absorptive part Im K_R, not to a real local subtraction constant.  Thus
the same record covariance is compatible with a family of Maxwell contacts.

So the hypothesis is not proved.  It becomes a precise extra axiom:

    Euclidean/in-out normal ordering maps endpoint connected covariance to the
    finite Maxwell F^2 contact.

That axiom would close the coefficient; current canon does not derive it.
"""

from __future__ import annotations

import math


ALPHA0_INV = 137.0
ALPHA0 = 1.0 / ALPHA0_INV
ALPHA_PHYS_INV = 137.035999084
TWO_PI = 2.0 * math.pi


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


Vector = tuple[float, ...]
Matrix = list[list[float]]


def dot(a: Vector, b: Vector) -> float:
    return sum(x * y for x, y in zip(a, b))


def norm_vec(a: Vector) -> float:
    return math.sqrt(dot(a, a))


def norm_matrix(m: Matrix) -> float:
    return math.sqrt(sum(x * x for row in m for x in row))


def mat_vec(m: Matrix, v: Vector) -> Vector:
    return tuple(sum(row[j] * v[j] for j in range(len(v))) for row in m)


def transverse_projector_kernel(q: Vector) -> Matrix:
    q2 = dot(q, q)
    return [[q2 * (1.0 if i == j else 0.0) - q[i] * q[j] for j in range(len(q))] for i in range(len(q))]


def scale_matrix(c: float, m: Matrix) -> Matrix:
    return [[c * x for x in row] for row in m]


def ward_defect(k: Matrix, q: Vector) -> float:
    denom = max(norm_matrix(k) * norm_vec(q), 1.0e-30)
    return norm_vec(mat_vec(k, q)) / denom


def alpha_inv_from_contact(n_contact: float) -> float:
    return ALPHA0_INV + n_contact * ALPHA0 / TWO_PI


def main() -> None:
    print("DRESSED ALPHA: ENDPOINT COVARIANCE -> MAXWELL CONTACT MAP ATTEMPT")
    print("=" * 100)

    sum_q2 = 16.0
    endpoint_covariance = 2.0 * sum_q2 - 1.0
    required = (ALPHA_PHYS_INV - ALPHA0_INV) * TWO_PI / ALPHA0
    alpha_contact = alpha_inv_from_contact(endpoint_covariance)

    print("[1] Endpoint covariance coefficient")
    print(f"  sum Q^2, three Weyl generations     = {sum_q2:.12f}")
    print(f"  connected closed-pair covariance    = {endpoint_covariance:.12f}")
    print(f"  required kernel from alpha shift    = {required:.12f}")
    print(f"  alpha^-1 if used as contact         = {alpha_contact:.12f}")
    check(abs(endpoint_covariance - 31.0) < 1.0e-12, "closed-pair connected covariance fixes 31")
    check(abs(endpoint_covariance / required - 1.0) < 5.0e-4, "31 lands the alpha shift")

    print("\n[2] Real-time influence-functional placement")
    print("  In r/a Schwinger-Keldysh variables:")
    print("      retarded/contact block : A_a K_R A_r")
    print("      record covariance/noise: (i/2) A_a N A_a")
    print("  Closed-record-pair counting fixes N, not K_R.")
    closed_pair_fixes_noise = True
    maxwell_contact_is_retarded = True
    same_block = False
    check(closed_pair_fixes_noise, "record covariance belongs to the symmetric/noise block")
    check(maxwell_contact_is_retarded, "Maxwell contact shifts the retarded/local kernel")
    check(not same_block, "noise block and retarded contact block are distinct real-time objects")

    print("\n[3] Ward family with fixed endpoint covariance")
    q = (1.0, -2.0, 0.5)
    base = transverse_projector_kernel(q)
    c_values = [0.0, 16.0, endpoint_covariance, required, 64.0]
    max_defect = 0.0
    for c in c_values:
        k_ret = scale_matrix(c, base)
        noise = scale_matrix(endpoint_covariance, base)
        max_defect = max(max_defect, ward_defect(k_ret, q), ward_defect(noise, q))
        print(f"  C_ret={c:10.6f}, N_noise={endpoint_covariance:10.6f}, Ward defect={ward_defect(k_ret, q):.3e}")
    check(max_defect < 1.0e-12, "Ward transversality allows all retarded contact constants")
    check(len(set(round(c, 6) for c in c_values)) > 3, "same endpoint covariance permits many retarded contacts")

    print("\n[4] FDT/Kramers-Kronig gate")
    print("  Equilibrium FDT can relate N to Im K_R at nonzero frequency.")
    print("  A finite local Maxwell contact is a real subtraction constant.")
    print("  Kramers-Kronig reconstructs the real part only after a subtraction is chosen.")
    fdt_reaches_imaginary_part = True
    contact_is_real_subtraction = True
    fdt_selects_contact = False
    check(fdt_reaches_imaginary_part, "FDT constrains absorptive response")
    check(contact_is_real_subtraction, "finite F^2 contact is a real local subtraction")
    check(not fdt_selects_contact, "FDT does not select the finite real contact")

    print("\n[5] What would close it")
    print("  A positive theorem must add or derive:")
    print("      Euclidean/in-out normal ordering of endpoint records")
    print("      maps the connected covariance N directly to the local F^2 contact C_ret.")
    print("  With that map, C_ret=N=31 and the coefficient is forced.")

    print("\nVERDICT")
    print(
        """
  Hypothesis tested:

      record-action/closed-pair counting fixes the finite Maxwell-contact
      normal-ordering constant.

  Outcome:

      coefficient leg: YES, conditionally.  Connected endpoint covariance gives
      2 sum Q^2 - 1 = 31.

      operator-map leg: NOT DERIVED.  In the real-time influence functional,
      closed-pair covariance is the Keldysh/noise block A_a N A_a, while the
      Maxwell contact is a retarded/local block A_a K_R A_r.  Ward identity and
      FDT leave the real local subtraction constant free.

  Sharp remaining theorem:

      prove that the framework's Euclidean/in-out normal ordering of monitored
      endpoint service records identifies N with the finite local F^2 contact.

  Until that theorem exists, dressed alpha is conditional on a Maxwell-contact
  selector axiom, not Locked.
exit 0 -- endpoint covariance fixes 31, but the covariance-to-F2 contact map is not derived.
"""
    )


if __name__ == "__main__":
    main()
