#!/usr/bin/env python3
r"""Dressed alpha: why the second-order endpoint term is not form-independent.

Question tested
---------------
Can the remaining dressed-alpha residual be derived *form-independently* from the
second connected endpoint covariance in the monitored Wilson-endpoint layer?

The leading endpoint-contact theorem fixes

    c1 = 2 sum Q^2 - 1 = 31,

so the one-contact prediction is

    alpha^{-1} = 137 + 31 alpha0/(2 pi) = 137.036013162,

102.7 ppb high.  The natural second-order charge structure is Q^4:

    c2 = -K2 sum Q^4,  sum Q^4 = 88/9,

where K2 is the scalar endpoint two-loop/vertex kernel.  The preceding scripts
showed that K2 = 1 (bare fourth moment) improves the residual to ~6 ppb, while
simple connected-cumulant weights miss.  This script isolates the reason:

    the endpoint charge algebra fixes sum Q^4, but it does not fix K2.

Equivalently, there is a continuous family of endpoint influence functionals with
the same leading c1 = 31 and the same Q^4 charge tensor, but different K2.  The
observed value asks for K2 ~= 1.06, a modest loop/normal-ordering weight.  That is
a real second-order endpoint loop integral, not a form-independent consequence of
the record covariance alone.

So this is a no-go for the *form-independent* closure, plus a precise target for
the genuine remaining theorem.
"""
from __future__ import annotations

import math
from fractions import Fraction

TWO_PI = 2.0 * math.pi
ALPHA0_INV = 137.0
ALPHA0 = 1.0 / ALPHA0_INV
X = ALPHA0 / TWO_PI

# CODATA 2018 is the convention used by the existing dressed-alpha scripts; the
# 2022 value shows that the conclusion is insensitive to the latest convention.
CODATA_2018 = 137.035999084
CODATA_2022 = 137.035999177

ok = True


def check(label: str, condition: bool) -> None:
    global ok
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    ok = ok and bool(condition)


def ppb(value: float, ref: float = CODATA_2018) -> float:
    return (value / ref - 1.0) * 1.0e9


def charges_one_generation() -> list[Fraction]:
    # One SO(10) 16, written as Weyl charge labels.  Triplication is done below.
    return (
        [Fraction(2, 3)] * 3
        + [Fraction(-1, 3)] * 3
        + [Fraction(-2, 3)] * 3
        + [Fraction(1, 3)] * 3
        + [Fraction(0), Fraction(-1), Fraction(1), Fraction(0)]
    )


def alpha_inv_with_kernel(k2: float) -> float:
    """Plain alpha0-power expansion with c2 = -k2 * sumQ4."""
    sum_q2 = 16.0
    sum_q4 = 88.0 / 9.0
    c1 = 2.0 * sum_q2 - 1.0
    c2 = -k2 * sum_q4
    return ALPHA0_INV + c1 * X + c2 * X * X


def kernel_required(ref: float) -> float:
    sum_q2 = 16.0
    sum_q4 = 88.0 / 9.0
    c1 = 2.0 * sum_q2 - 1.0
    leading = ALPHA0_INV + c1 * X
    c2_required = (ref - leading) / (X * X)
    return abs(c2_required) / sum_q4


def main() -> int:
    print("DRESSED ALPHA: second-order form-independence no-go / loop target")
    print("=" * 82)

    q = charges_one_generation() * 3
    sum_q = sum(q)
    sum_q2 = sum(x * x for x in q)
    sum_q4 = sum(x ** 4 for x in q)
    print("\n[1] Endpoint charge tensors")
    print(f"    records = {len(q)}; sum Q = {sum_q}; sum Q^2 = {sum_q2}; sum Q^4 = {sum_q4}")
    check("charge neutral", sum_q == 0)
    check("leading tensor fixed: sum Q^2 = 16", sum_q2 == 16)
    check("second charge tensor fixed: sum Q^4 = 88/9", sum_q4 == Fraction(88, 9))

    leading = alpha_inv_with_kernel(0.0)
    bare_q4 = alpha_inv_with_kernel(1.0)
    k2_2018 = kernel_required(CODATA_2018)
    k2_2022 = kernel_required(CODATA_2022)
    exact_2018 = alpha_inv_with_kernel(k2_2018)
    exact_2022 = alpha_inv_with_kernel(k2_2022)

    print("\n[2] Same leading contact, variable second-order endpoint kernel")
    print(f"    K2 = 0      (one-contact only)      alpha^-1 = {leading:.12f}  {ppb(leading):+7.1f} ppb")
    print(f"    K2 = 1      (bare Q^4 moment)       alpha^-1 = {bare_q4:.12f}  {ppb(bare_q4):+7.1f} ppb")
    print(f"    K2 = {k2_2018:.6f} (CODATA-2018 target)   alpha^-1 = {exact_2018:.12f}  {ppb(exact_2018):+7.1f} ppb")
    print(f"    K2 = {k2_2022:.6f} (CODATA-2022 target)   alpha^-1 = {exact_2022:.12f}  {ppb(exact_2022, CODATA_2022):+7.1f} ppb")
    check("all K2 choices preserve the leading c1 = 31 endpoint contact", abs((leading - ALPHA0_INV) / X - 31.0) < 1e-12)
    check("bare Q^4 reaches the few-ppb near-miss, not closure", 3.0 < abs(ppb(bare_q4)) < 10.0)
    check("data require a modest 5-8 percent loop weight over the bare Q^4 moment", 1.05 < k2_2022 < k2_2018 < 1.08)

    print("\n[3] Consequence for the proposed form-independent derivation")
    print("    Charge algebra gives the tensor:")
    print("        c2 = -K2 * sum Q^4.")
    print("    It does not give the scalar K2.  K2=1 is the bare fourth moment;")
    print("    K2~=1.06 is the observed endpoint-loop target; both have the same")
    print("    monitored endpoint records and the same leading 31.")

    # A continuous family with identical leading endpoint algebra but different
    # second-order predictions proves non-identifiability from the covariance form
    # alone.  The separation between K2=1 and K2=Ktarget is the exact missing theorem.
    family = [0.8, 1.0, k2_2022, k2_2018, 1.2]
    ppbs = [ppb(alpha_inv_with_kernel(k)) for k in family]
    print("\n    K2 family residuals (same c1 and same Q^4 charge tensor):")
    for k, residual in zip(family, ppbs):
        print(f"        K2={k:8.6f} -> residual {residual:+8.1f} ppb")
    check("varying K2 changes the residual by tens of ppb while leaving the charge tensors fixed",
          max(ppbs) - min(ppbs) > 30.0)

    delta_c2 = (k2_2018 - 1.0) * float(sum_q4)
    print(
        f"""
[4] VERDICT
    The pasted attempt is a real advance, but it does not close dressed alpha.

    What is derived form-independently:
      * the leading endpoint-contact coefficient 31;
      * the second-order charge structure Q^4, with sum Q^4 = 88/9;
      * the sign and scale of the correction, reducing 102.7 ppb to a few ppb.

    What is not derived:
      * the scalar endpoint two-loop kernel K2.  CODATA asks for
        K2 = {k2_2018:.6f} (2018 convention) or {k2_2022:.6f} (2022 convention),
        i.e. an extra delta c2 = {delta_c2:.3f} over the bare Q^4 coefficient.

    Therefore a purely form-independent connected-covariance/counting theorem is
    impossible at second order: the covariance fixes the Q^4 tensor but leaves the
    loop kernel free.  The honest remaining theorem is narrower and more physical:

        compute the monitored Wilson-endpoint internal-photon/vertex kernel and
        prove K2 ~= 1.06 from that dynamics, with no fitted C2.

    The old 24/7 term is not rescued; it is a fitted stand-in for this endpoint
    loop kernel.  The new target is real, small, and well-defined, but it is a
    loop calculation rather than another combinatorial cumulant.
exit 0"""
    )
    print("ALL CHECKS PASSED" if ok else "CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
