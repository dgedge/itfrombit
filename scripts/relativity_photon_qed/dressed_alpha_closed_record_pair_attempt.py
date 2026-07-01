#!/usr/bin/env python3
r"""Dressed alpha: closed-record-pair rescue attempt.

This is the one new lever supplied by the recent reconstruction work.  The
Born/record-action result says that probabilities are closed forward/backward
record pairs.  Numerically that is exactly the sort of factor the dressed-alpha
problem wants:

    physical shift needs  N_eff ~= 31,
    charge-weighted Weyl kernel over three generations is 16,
    2 * 16 - 1 = 31.

So the question is sharp:

    Can the physical low-energy QED coupling bill a closed-record-pair kernel
    2*sum(Q^2) - 1, rather than the ordinary retarded Ward/Kubo kernel sum(Q^2)?

Verdict
-------
No under the current axioms.  The closed-record pair is the right object for
Born probabilities/noise (amplitude times conjugate amplitude).  The measured
fine-structure constant is a retarded electromagnetic response: the Peierls
current-current Hessian / photon self-energy.  Retarded response uses the
commutator/Kubo moment, not the closed Keldysh norm.  Doubling the Weyl kernel
therefore changes the observable class.  The additional "-1" subtraction is
also unselected by the Ward identity or record-pair theorem.

The attempt is useful because it explains why the old N1=31 fit is so
seductive: it is the closed-record-pair doubling of the real charge-weighted
Weyl kernel, with one extra subtraction.  But that is not QED alpha unless a
new theorem says that the Thomson/LSZ charge is measured by a closed-pair
noise kernel rather than by the retarded Ward kernel.
"""

from __future__ import annotations

import math


ALPHA0_INV = 137.0
ALPHA0 = 1.0 / ALPHA0_INV
ALPHA_PHYS_INV = 137.035999084
DELTA = ALPHA_PHYS_INV - ALPHA0_INV
TWO_PI = 2.0 * math.pi


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def inv_alpha_from_kernel(kernel: float) -> float:
    return ALPHA0_INV + kernel * ALPHA0 / TWO_PI


def main() -> None:
    print("DRESSED ALPHA: CLOSED-RECORD-PAIR RESCUE ATTEMPT")
    print("=" * 96)

    n_required = DELTA * TWO_PI / ALPHA0
    sm_dirac = 8.0
    sm_weyl = 16.0
    closed_pair = 2.0 * sm_weyl
    closed_pair_minus_identity = closed_pair - 1.0

    print("[1] Numerical attraction")
    print(f"  target delta(alpha^-1)                  = {DELTA:.9f}")
    print(f"  required one-loop-style kernel          = {n_required:.6f}")
    print(f"  charge-weighted Dirac kernel            = {sm_dirac:.1f}")
    print(f"  charge-weighted Weyl/SO(10) kernel      = {sm_weyl:.1f}")
    print(f"  closed-record-pair doubled Weyl kernel  = {closed_pair:.1f}")
    print(f"  doubled Weyl minus one identity channel = {closed_pair_minus_identity:.1f}")
    ia_31 = inv_alpha_from_kernel(closed_pair_minus_identity)
    print(f"  resulting alpha^-1                      = {ia_31:.9f}")
    print(f"  residual vs observed                    = {(ia_31 / ALPHA_PHYS_INV - 1) * 1e9:.1f} ppb")
    check(30.8 < n_required < 31.1, "observed shift asks for N_eff ~= 31")
    check(closed_pair_minus_identity == 31.0, "2*16-1 gives the historical N1=31")
    check(abs(ia_31 / ALPHA_PHYS_INV - 1) < 2e-7, "2*16-1 reproduces alpha^-1 at near-hit precision")

    print("\n[2] Observable-class test")
    print("  Closed-record pair:")
    print("    probability/noise object, schematically  A A*  or Keldysh symmetric correlator.")
    print("  QED alpha:")
    print("    retarded response object, schematically  i theta(t)<[J(t),J(0)]>  plus diamagnetic term.")
    print("  Ward/Kubo identities fix the latter as the Peierls-current Hessian.")
    print("  A closed-pair doubling is therefore not automatically a response coefficient.")

    # These booleans encode the current theorem status.  They are intentionally
    # explicit rather than hidden in a fit: the reconstruction scripts close Born
    # probabilities, while the dressed-alpha Ward scripts close the response
    # identity negatively.
    born_pair_closes_probabilities = True
    ward_kubo_alpha_is_retarded = True
    closed_pair_equals_retarded_response = False
    check(born_pair_closes_probabilities, "closed-record-pair theorem applies to Born probabilities")
    check(ward_kubo_alpha_is_retarded, "physical alpha is the retarded Ward/Kubo Peierls response")
    check(not closed_pair_equals_retarded_response, "no current theorem identifies closed-pair noise with retarded alpha")

    print("\n[3] The -1 subtraction is not selected")
    for subtract in (0, 1, 2):
        kernel = closed_pair - subtract
        ia = inv_alpha_from_kernel(kernel)
        print(
            f"  2*16-{subtract:<1d} -> kernel={kernel:4.1f}, "
            f"alpha^-1={ia:.9f}, shift/target={(ia - ALPHA0_INV) / DELTA:.6f}"
        )
    print("  The -1 choice is numerically excellent, but it is not fixed by the")
    print("  closed-record-pair theorem.  Record closure supplies the factor of two")
    print("  for probabilities; it does not supply a unique photon/identity subtraction")
    print("  for a retarded electromagnetic coupling.")
    check(abs((inv_alpha_from_kernel(32.0) - ALPHA0_INV) / DELTA - 1.0) > 0.02, "2*16 without subtraction is not the observed value")

    print("\nVERDICT")
    print(
        """
  Progress, but not closure:

    * The recent closed-record-pair result explains the numerical form of the
      old N1=31 hit: it is 2*(charge-weighted Weyl kernel 16) minus one.

    * That makes N1=31 less mysterious, but it does not make it physical QED
      alpha.  The record-pair theorem licenses Born probabilities / symmetric
      noise.  Low-energy alpha is a retarded Ward/Kubo response coefficient.
      Those are different real-time components of the Schwinger-Keldysh
      structure.

    * To promote this route, the missing theorem would have to be explicit:
      the Thomson/LSZ charge must bill the closed-record symmetric kernel
      2*sum(Q^2)-1 rather than the retarded Peierls Hessian.  That is a new
      EM response principle, not a corollary of Born/record closure.

  Therefore the new foundations result narrows the dressed-alpha mystery but
  does not solve it.  The residual is now: retarded Ward response versus
  closed-record noise.  The bare alpha0=1/137 remains derived; the exact
  137.035999... dressing remains bounded and unpromoted.
exit 0 -- closed-record-pair rescue is a named near-hit, rejected as current QED alpha.
"""
    )


if __name__ == "__main__":
    main()
