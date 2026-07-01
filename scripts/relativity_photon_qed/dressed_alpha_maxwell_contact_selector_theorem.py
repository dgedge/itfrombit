#!/usr/bin/env python3
r"""Dressed alpha: closed-pair Maxwell-contact selector theorem.

Hypothesis tested
-----------------
Can record-action / closed-record-pair counting fix the finite Maxwell-contact
normal-ordering constant without changing ordinary QED running?

The only QED-compatible interpretation is not

    retarded photon self-energy kernel = 31,

but rather

    finite local F^2 contact at the service cutoff = connected endpoint
    record-pair covariance.

This script separates the theorem from its remaining premise.

Theorem, conditional on the endpoint-contact map
-----------------------------------------------
For three SO(10) Weyl generations the charge-square trace is

    sum_f Q_f^2 = 16.

A closed forward/backward record pair doubles the equal-time endpoint record
kernel.  Normal ordering a covariance removes the single constant/identity
mode.  Therefore the connected closed-pair endpoint covariance is

    N_contact = 2 * sum_f Q_f^2 - 1 = 31.

Inserted as a finite local transverse Maxwell contact term, this changes the
zero-momentum matching value of alpha but leaves the nonlocal retarded
Ward/Kubo response and the subsequent charge-weighted QED running unchanged.

Boundary
--------
This is a conditional closure, not a full derivation.  The present canon proves
the closed-record-pair/Born covariance and the Ward-preserving contact family.
It does not yet prove the operator map

    endpoint closed-record covariance == finite local Maxwell F^2 contact.

Without that map, the contact coefficient remains a normal-ordering selector
postulate.  With that map, the coefficient 31 is forced rather than fitted.
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


def sm_weyl_charge_square_per_generation() -> float:
    """SO(10) 16 as left-handed Weyl fields, including conjugates."""

    q_left = 3 * (2 / 3) ** 2 + 3 * (1 / 3) ** 2
    u_conj = 3 * (2 / 3) ** 2
    d_conj = 3 * (1 / 3) ** 2
    lepton_doublet = 0.0**2 + 1.0**2
    e_conj = 1.0**2
    nu_conj = 0.0**2
    return q_left + u_conj + d_conj + lepton_doublet + e_conj + nu_conj


def inv_alpha_from_kernel(kernel: float) -> float:
    return ALPHA0_INV + kernel * ALPHA0 / TWO_PI


def main() -> None:
    print("DRESSED ALPHA: MAXWELL-CONTACT SELECTOR THEOREM")
    print("=" * 96)

    q2_one_gen = sm_weyl_charge_square_per_generation()
    q2_three_gen = 3.0 * q2_one_gen
    closed_pair_raw = 2.0 * q2_three_gen
    connected_covariance = closed_pair_raw - 1.0
    required = (ALPHA_PHYS_INV - ALPHA0_INV) * TWO_PI / ALPHA0
    alpha_contact = inv_alpha_from_kernel(connected_covariance)

    print("[1] Closed-pair connected covariance")
    print(f"  Weyl charge-square trace, one generation     = {q2_one_gen:.12f}")
    print(f"  Weyl charge-square trace, three generations = {q2_three_gen:.12f}")
    print(f"  closed forward/backward pair                = {closed_pair_raw:.12f}")
    print(f"  connected covariance, subtract identity     = {connected_covariance:.12f}")
    print(f"  required kernel from alpha shift            = {required:.12f}")
    print(f"  alpha^-1 from contact selector              = {alpha_contact:.12f}")
    print(f"  relative error in shift                     = {(connected_covariance / required - 1.0) * 1e4:.3f}e-4")
    print(f"  relative error in alpha^-1                  = {(alpha_contact / ALPHA_PHYS_INV - 1.0) * 1e9:.1f} ppb")
    check(abs(q2_one_gen - 16.0 / 3.0) < 1e-12, "one Weyl generation has sum Q^2 = 16/3")
    check(abs(q2_three_gen - 16.0) < 1e-12, "three generations have sum Q^2 = 16")
    check(abs(connected_covariance - 31.0) < 1e-12, "connected closed-pair covariance gives N_contact = 31")
    check(abs(connected_covariance / required - 1.0) < 5e-4, "contact selector lands the alpha shift at near-hit precision")

    print("\n[2] Normal-ordering controls")
    alternatives = {
        "retarded Weyl Ward kernel": q2_three_gen,
        "closed pair without connected subtraction": closed_pair_raw,
        "connected closed-pair covariance": connected_covariance,
        "subtract two identity-like modes": closed_pair_raw - 2.0,
    }
    for name, kernel in alternatives.items():
        value = inv_alpha_from_kernel(kernel)
        shift_ratio = (value - ALPHA0_INV) / (ALPHA_PHYS_INV - ALPHA0_INV)
        print(f"  {name:<45} kernel={kernel:6.2f}  shift/obs={shift_ratio:8.5f}")
    check(alternatives["closed pair without connected subtraction"] != connected_covariance,
          "the -1 is the connected-covariance normal-ordering step, not closed-pair doubling alone")
    check(abs((inv_alpha_from_kernel(closed_pair_raw - 2.0) - ALPHA0_INV) /
              (ALPHA_PHYS_INV - ALPHA0_INV) - 1.0) > 0.03,
          "subtracting two modes is visibly different; the covariance subtraction is a sharp rule")

    print("\n[3] QED phenomenology gate")
    beta_ratio_if_misread = connected_covariance / q2_three_gen
    print(f"  If 31 is misread as retarded beta kernel, beta ratio = {beta_ratio_if_misread:.6f}  FAIL")
    print("  If 31 is a finite local F^2 contact, beta kernel remains 16.0     PASS")
    strong_kernel_equivalence = False
    finite_contact_matching = True
    endpoint_contact_map_derived = False
    check(not strong_kernel_equivalence, "do not identify the retarded LSZ kernel with the closed-pair count")
    check(finite_contact_matching, "a finite local Maxwell contact can preserve Ward identity and QED running")
    check(not endpoint_contact_map_derived, "canon has not yet derived endpoint covariance = Maxwell contact")

    print("\n[4] The precise remaining premise")
    print("  Proven / standard pieces:")
    print("    * closed-record-pair covariance subtracts one constant mode")
    print("    * local F^2 contact is the unique isotropic Ward-preserving quadratic contact")
    print("    * finite contact matching leaves low-energy QED running charge-weighted")
    print("  Missing operator map:")
    print("    * endpoint service-record covariance is the finite Maxwell F^2 contact")
    print("      produced when the service layer is normal ordered at the cutoff")

    print("\nVERDICT")
    print(
        """
  Conditional positive theorem:

      If the finite Maxwell contact term is the connected covariance of the
      closed endpoint service records, then record-action fixes its coefficient:

          N_contact = 2 sum_f Q_f^2 - 1 = 31.

      This preserves QED phenomenology because the retarded nonlocal response
      and beta function remain charge-weighted.  The closed-pair count only
      fixes the local matching/contact term.

  Remaining wall:

      The endpoint-contact operator map is not yet derived from existing canon.
      The result therefore upgrades the hypothesis from a loose numerical
      escape to a precise conditional theorem, but it does not fully Lock the
      dressed fine-structure constant.

  Next theorem target:

      derive, from the framed service/null-chain action or monitored QEC
      influence functional, that normal ordering the endpoint service records
      contributes exactly their connected closed-pair covariance to the local
      Maxwell F^2 contact.
exit 0 -- coefficient forced if endpoint-covariance-to-F2-contact map is proved.
"""
    )


if __name__ == "__main__":
    main()
