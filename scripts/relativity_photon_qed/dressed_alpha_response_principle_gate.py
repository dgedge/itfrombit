#!/usr/bin/env python3
r"""Dressed alpha: can a closed-pair response principle preserve QED?

Question
--------
The closed-record-pair arithmetic is compelling:

    three Weyl generations:  sum_f Q_f^2 = 16,
    closed pair minus one:   2 * 16 - 1 = 31,

and the kernel 31 gives the observed 137 -> 137.036 shift at near-hit
precision.  The remaining long-shot proposal is therefore:

    prove that the retarded Thomson/LSZ charge residue bills the closed-pair
    kernel 31, while preserving ordinary QED phenomenology.

Result
------
The strong principle cannot hold.  If the retarded photon self-energy kernel is
31 rather than the charge-weighted Weyl value 16, the QED vacuum-polarisation
coefficient and running are multiplied by 31/16.  That is not a convention; it
changes the response observable.

There is one mathematically consistent escape, but it is weaker than the
requested theorem: the closed-pair kernel can be installed as a finite local
Maxwell contact term / subtraction constant at the service cutoff.  A local
transverse F^2 counterterm preserves Ward identities and leaves the subsequent
charge-weighted QED beta function intact.  But this is a boundary matching
axiom, not a theorem from the present record-action machinery.  Keldysh/FDT
relates noise to the absorptive part of retarded response, not to the finite
subtraction constant.  The subtraction constant is exactly where renormalisation
conditions live.

So the positive route is:

    not "retarded kernel = closed-pair kernel",
    but "record closure fixes the finite Maxwell contact term".

The current canon does not prove that contact-term selector.  This script
records the strongest compatible form and its remaining missing clause.
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


def inv_alpha_from_kernel(kernel: float) -> float:
    return ALPHA0_INV + kernel * ALPHA0 / TWO_PI


def main() -> None:
    print("DRESSED ALPHA: CLOSED-PAIR RESPONSE PRINCIPLE GATE")
    print("=" * 96)

    q2_weyl = 16.0
    closed_pair = 2.0 * q2_weyl - 1.0
    required = (ALPHA_PHYS_INV - ALPHA0_INV) * TWO_PI / ALPHA0
    alpha_closed = inv_alpha_from_kernel(closed_pair)

    print("[1] Arithmetic target")
    print(f"  required one-loop-style kernel        = {required:.9f}")
    print(f"  Weyl charge-weighted kernel           = {q2_weyl:.9f}")
    print(f"  closed-pair kernel 2*16-1             = {closed_pair:.9f}")
    print(f"  alpha^-1 from closed-pair contact     = {alpha_closed:.9f}")
    print(f"  residual                              = {(alpha_closed / ALPHA_PHYS_INV - 1) * 1e9:.1f} ppb")
    check(abs(closed_pair - 31.0) < 1e-12, "closed-pair arithmetic reproduces N=31")
    check(abs(alpha_closed / ALPHA_PHYS_INV - 1) < 2e-7, "N=31 is numerically near the observed shift")

    print("\n[2] Strong retarded-kernel principle")
    beta_ratio = closed_pair / q2_weyl
    print("  Hypothesis: replace the retarded Ward/Kubo self-energy kernel")
    print("  by the closed-pair kernel itself.")
    print(f"  beta / Weyl beta coefficient ratio    = {beta_ratio:.6f}")
    print("  This would nearly double vacuum polarisation and QED running.")
    check(beta_ratio > 1.9, "strong kernel replacement changes QED phenomenology")
    print("  verdict: FAIL as a physical QED response principle.")

    print("\n[3] Finite-contact escape")
    print("  Decompose the transverse inverse propagator schematically as")
    print("      D_T^{-1}(q^2) = Z_contact q^2 + Pi_ret(q^2),")
    print("  with q_mu Pi^{mu nu}=0.  A finite local F^2 contact term shifts")
    print("  Z_contact but does not alter the charge-weighted beta function below")
    print("  the matching scale.")
    preserves_ward = True
    leaves_running = True
    selected_by_current_axioms = False
    print(f"  contact kernel chosen from closed pair = {closed_pair:.1f}")
    print("  beta kernel remains charge-weighted    = 16.0")
    check(preserves_ward, "local F^2 contact term preserves Ward transversality")
    check(leaves_running, "finite contact matching leaves ordinary low-energy QED running intact")
    check(not selected_by_current_axioms, "current axioms do not select this finite subtraction")

    print("\n[4] Why FDT/Keldysh does not prove the selector")
    print("  Closed record pairs are symmetric/noise data.  FDT relates noise to")
    print("  the absorptive spectral density of the retarded response.  Kramers-Kronig")
    print("  reconstructs the real retarded part only up to a subtraction constant.")
    print("  The proposed N=31 contact is precisely such a subtraction constant.")
    fdt_fixes_absorptive_part = True
    kk_needs_subtraction = True
    subtraction_selected = False
    check(fdt_fixes_absorptive_part, "equilibrium noise can constrain Im retarded response")
    check(kk_needs_subtraction, "the local Maxwell contact term is a subtraction constant")
    check(not subtraction_selected, "closed-record-pair theorem does not currently select the subtraction")

    print("\nVERDICT")
    print(
        """
  No positive proof under current canon.

  The requested strong theorem,

      retarded Thomson/LSZ residue = closed-pair kernel 31,

  is incompatible with standard QED phenomenology if it is read as the
  self-energy / beta-function kernel.  It replaces the charge-weighted
  retarded current response by a symmetric record-noise count and multiplies
  the Weyl vacuum-polarisation coefficient by 31/16.

  The only compatible form is weaker and very sharp:

      the closed-pair kernel fixes a finite, local, Ward-preserving Maxwell
      contact term at the service cutoff, while the nonlocal retarded response
      and running remain the usual charge-weighted QED objects.

  That would preserve QED phenomenology, but it is not derived.  It is a new
  normal-ordering / matching-selector axiom unless the record-action machinery
  can prove why closed record-pair counting chooses the finite F^2 subtraction
  constant.  The remaining target is therefore not LSZ-kernel equivalence; it is
  a Maxwell-contact selector theorem.
exit 0 -- strong response principle fails; finite-contact escape is compatible but unselected.
"""
    )


if __name__ == "__main__":
    main()
