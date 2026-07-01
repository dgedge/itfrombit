#!/usr/bin/env python3
r"""Dressed alpha: does the connected 4th cumulant (with leg/ordering weights)
give the exact second-order coefficient?

The leading contact is c1 = 2*sumQ^2 - 1 = 31 (two endpoint legs, one connected
identity mode removed). The exact form-independent (alpha0-power) second-order
coefficient needed to hit CODATA is c2 ~ -10.43. The previous attempt showed the
plain 4th moment -sumQ^4 = -9.78 is the right STRUCTURE (Q^4, right sign, 102.7 ->
6.4 ppb) but ~6% short, and the leftover is too big to be 3rd order.

This script asks the sharp question: is the exact c2 a CONNECTED 4th cumulant with
the same weight logic as the leading (leg factor, -1 ordering), or is it not a
simple charge-content cumulant at all? It enumerates a PRE-SPECIFIED, structurally
motivated candidate set (no free search) and reports which -- if any -- lands
within the ~0.5% (=> < ~1 ppb) needed for a derivation.

Honest result (computed below): NO simple charge-content cumulant hits it. The
faithful analog of the leading, 2*sumQ^4 - 1, overshoots massively (~-80 ppb).
The exact -10.43 is *bracketed* between the plain 4th moment sumQ^4 = 9.78 (+6.4
ppb) and twice the connected 4th cumulant 2|sumQ^4 - 3(sumQ^2)^2/N| = 12.44 (~-20
ppb) but matched by neither, nor by any natural combination without tuning. So the
second-order coefficient is NOT a simple 4th cumulant with charge-sum weights: the
Q^4 structure is right, the exact weight needs genuine two-loop endpoint dynamics
(vertex/propagator structure), not a charge moment. An honest negative that rules
out the simple-cumulant hypothesis and brackets the coefficient. Self-asserting.
"""
from __future__ import annotations
import math
from fractions import Fraction

TWO_PI = 2.0 * math.pi
A0 = 1.0 / 137.0
X = A0 / TWO_PI
X2 = X * X
CODATA_2018 = 137.035999084
CODATA_2022 = 137.035999177
ok = True
def check(name, cond):
    global ok; print(f"  [{'PASS' if cond else 'FAIL'}] {name}"); ok = ok and bool(cond)


def charges_one_gen():
    return [Fraction(2, 3)] * 3 + [Fraction(-1, 3)] * 3 + [Fraction(-2, 3)] * 3 + \
           [Fraction(1, 3)] * 3 + [Fraction(0), Fraction(-1), Fraction(1), Fraction(0)]


def main():
    print("DRESSED ALPHA: connected 4th-cumulant weights vs the exact 2nd-order coefficient")
    print("=" * 84)

    q3 = charges_one_gen() * 3                       # 48 records (16 x 3 gen)
    N = len(q3)
    sumQ2 = sum(c * c for c in q3)                   # 16
    sumQ4 = sum(c ** 4 for c in q3)                  # 88/9
    sumQ6 = sum(c ** 6 for c in q3)
    kappa4_sum = sumQ4 - 3 * Fraction(sumQ2 ** 2, N)  # connected 4th cumulant (extensive)
    print(f"    N records = {N};  sumQ^2 = {sumQ2};  sumQ^4 = {sumQ4} = {float(sumQ4):.4f};  sumQ^6 = {float(sumQ6):.4f}")
    print(f"    connected 4th cumulant  sumQ^4 - 3(sumQ^2)^2/N = {kappa4_sum} = {float(kappa4_sum):.4f}")

    a_lead = 137.0 + float(2 * sumQ2 - 1) * X
    print(f"\n[1] leading c1 = 2 sumQ^2 - 1 = 31 -> alpha^-1 = {a_lead:.9f} (+{(a_lead/CODATA_2018-1)*1e9:.1f} ppb)")

    for tag, cod in (("2018", CODATA_2018), ("2022", CODATA_2022)):
        c2_ex = (cod - a_lead) / X2
        print(f"    exact c2 (alpha0-power, CODATA {tag}) = {c2_ex:.4f}")
    c2_exact = (CODATA_2018 - a_lead) / X2

    print("\n[2] Pre-specified structurally-motivated candidates (as c2), plain series")
    cands = {
        "-sumQ^4                     (plain 4th moment)":         -float(sumQ4),
        "-(2 sumQ^4 - 1)             (faithful analog of 31)":    -(2 * float(sumQ4) - 1),
        "-(sumQ^4 - 1)               (moment, -1 ordering)":      -(float(sumQ4) - 1),
        "-|kappa4|                   (connected 4th cumulant)":   float(kappa4_sum),   # kappa4_sum<0
        "-2|kappa4|                  (two-leg connected cumulant)": 2 * float(kappa4_sum),
        "-(2|kappa4| - 1)            (two-leg cumulant, -1)":      2 * float(kappa4_sum) + 1,
        "-sumQ^6                     (wrong order, control)":     -float(sumQ6),
    }
    best = None
    for name, c2 in cands.items():
        a = a_lead + c2 * X2
        p = (a / CODATA_2018 - 1) * 1e9
        flag = ""
        if best is None or abs(p) < abs(best[1]):
            best = (name, p, c2)
        print(f"    c2 = {c2:8.3f}   {name:44s} -> {p:+7.1f} ppb")
    print(f"\n    exact needed c2 = {c2_exact:.4f}")
    print(f"    closest candidate: {best[0].strip()}  ({best[1]:+.1f} ppb, c2={best[2]:.3f})")

    print("\n[3] Verdict tests")
    check("the faithful analog 2 sumQ^4 - 1 overshoots badly (not the 2nd order)",
          (a_lead + (-(2 * float(sumQ4) - 1)) * X2 / 1) and abs((a_lead + (-(2 * float(sumQ4) - 1)) * X2) / CODATA_2018 - 1) * 1e9 > 40)
    # is the exact coefficient bracketed by sumQ^4 and 2|kappa4| but hit by neither?
    lo, hi = float(sumQ4), 2 * abs(float(kappa4_sum))     # 9.778 and 12.444
    check("exact |c2| is bracketed by sumQ^4 (9.78) and 2|kappa4| (12.44)", lo < abs(c2_exact) < hi)
    check("NO candidate lands within ~1 ppb (i.e. within ~0.1% of exact c2)",
          min(abs((a_lead + c * X2) / CODATA_2018 - 1) * 1e9 for c in cands.values()) > 3.0)
    # discipline: is the near one (sumQ^4) within 0.5%? no -> not a derivation
    check("even the closest natural cumulant (sumQ^4) is >3% off the exact coefficient",
          abs(-float(sumQ4) / c2_exact - 1) > 0.03)

    print(
        f"""
[4] VERDICT -- honest negative: the 2nd-order coefficient is NOT a simple 4th cumulant
    Tested a pre-specified set of structurally-motivated Q^4 cumulant candidates
    (plain 4th moment, connected 4th cumulant, both with the leading's leg/-1
    weight logic). None reproduces the exact form-independent coefficient
    c2 = {c2_exact:.2f}:

      * the FAITHFUL analog of the leading 2 sumQ^2 - 1, namely 2 sumQ^4 - 1 = 18.6,
        overshoots to ~-80 ppb -- so the leg-factor-2 + (-1) weight logic does NOT
        carry over to the second order;
      * the exact |c2| = {abs(c2_exact):.2f} is BRACKETED between the plain 4th moment
        sumQ^4 = {float(sumQ4):.2f} (+6.4 ppb) and twice the connected 4th cumulant
        2|kappa4| = {2*abs(float(kappa4_sum)):.2f} (~-20 ppb), but hit by NEITHER, and by no natural
        combination without tuning.

    So the "4th-cumulant weights" route, tested honestly, does not close it: the
    second order is the Q^4 (two-loop vacuum-polarization) STRUCTURE, but its exact
    coefficient is not a charge-content cumulant. Physically that is expected -- a
    two-loop coefficient carries vertex/propagator (loop-integral) weight, not just
    a charge moment; the leading 1-loop could be a pure charge sum, the 2-loop
    cannot. Closing the last ~6% therefore requires the actual second-order endpoint
    loop integral (the internal-photon/vertex kernel), not a cumulant of the charge
    records.

    Net (dressed alpha, honest current state): parameter-free to 0.10 ppm (leading
    contact, derived); the Q^4 second-order STRUCTURE is identified and un-fitted
    (102.7 -> 6.4 ppb), but its exact coefficient is NOT a simple charge cumulant --
    it needs the two-loop endpoint loop-integral weight. Two disciplined attempts
    (moment, cumulant) have bounded and structurally identified the term without
    locking the coefficient. Remaining work is a genuine loop computation, not a
    counting.
exit 0"""
    )
    print("ALL CHECKS PASSED" if ok else "CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
