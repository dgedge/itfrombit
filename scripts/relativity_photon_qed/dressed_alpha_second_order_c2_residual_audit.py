#!/usr/bin/env python3
r"""Dressed alpha: the 102.7 ppb one-contact residual and the second-order C2.

Context (2026-07-01): the endpoint-contact completeness theorem
(`dressed_alpha_endpoint_contact_completeness_theorem.py`) closed the LEADING
Maxwell contact inside the monitored Wilson-endpoint layer: its coefficient is
forced to N1 = 2*sum_{3 gen Weyl} Q^2 - 1 = 31, with no in-layer knob. The
one-contact value alpha^-1 = 137 + N1*alpha0/(2pi) = 137.036013 sits 102.7 ppb
above CODATA. The canon's older two-term form
    alpha^-1(alpha^-1 - 137) = N1/(2pi) - |C2|/(2pi)^2 * 1/alpha^-1,  |C2| = 24/7,
lands at 3 ppb. This script does the next honest step on the residual: it
separates what is DERIVED from what is FIT, sizes the residual, and audits C2.

Findings (all computed below):
  * the LINEAR one-contact alpha^-1 = 137 + 31*alpha0/(2pi) = 137.036013 is
    parameter-free and form-independent (137=T(16)+1, 31=2*sumQ^2-1, 2pi all
    derived) -> 102.7 ppb = 0.10 ppm from CODATA. This is the genuine result.
  * the residual is O(alpha/2pi): the second-order shift is ~1/3 of alpha/2pi of
    the leading shift -- the natural size of a next-order QED term, not a knob.
  * of the 102.7 ppb, ~65 ppb is absorbed by merely CHOOSING the self-consistent
    quadratic form (the delta^2 term), which the canon itself flags as an ansatz;
    C2=24/7 then fits the last ~38 ppb.
  * C2 is a FIT, not razor-selected: many low-complexity rationals land within
    the datum's ~1% constraint (16.3 formula-freedom). So the 3 ppb match is not
    parameter-free; the parameter-free floor is the 102.7 ppb leading contact.

Net: dressed alpha is parameter-free to 0.10 ppm (leading contact, derived); the
open theorem is a FORM-INDEPENDENT second-order term from the 2nd-order endpoint
covariance in the monitored layer -- not a fitted C2. Self-asserting. exit 0.
"""
from __future__ import annotations
import math
from fractions import Fraction

TWO_PI = 2.0 * math.pi
A0 = 1.0 / 137.0
N1 = 31                      # 2*sum_{3 gen Weyl}Q^2 - 1, endpoint-contact theorem
CODATA = 137.035999084      # alpha^-1 (CODATA 2018; framework's value)
C2_CANON = Fraction(24, 7)

ok = True
def check(name, cond):
    global ok; print(f"  [{'PASS' if cond else 'FAIL'}] {name}"); ok = ok and bool(cond)

def ppb(x):
    return (x / CODATA - 1.0) * 1e9


def solve_recursive(c2):
    """Self-consistent alpha^-1(alpha^-1 - 137) = N1/2pi - c2/(2pi)^2 /alpha^-1."""
    x = 137.036
    for _ in range(200):
        rhs = N1 / TWO_PI - (c2 / TWO_PI**2) / x
        x = 137.0 + rhs / x           # from x(x-137)=rhs
    return x


def main():
    print("DRESSED ALPHA: 102.7 ppb one-contact residual + second-order C2 audit")
    print("=" * 82)

    print("\n[1] Leading one-contact (LINEAR, parameter-free) = the genuine result")
    a_lin = 137.0 + N1 * A0 / TWO_PI
    print(f"    alpha^-1 = 137 + 31*alpha0/(2pi) = {a_lin:.9f}   ({ppb(a_lin):+.1f} ppb = {ppb(a_lin)/1000:+.3f} ppm)")
    check("N1=31 = 2*sumQ^2-1 (endpoint-contact theorem)", N1 == 2 * 16 - 1)
    check("linear one-contact matches the canon 137.036013 value", abs(a_lin - 137.036013162) < 2e-6)
    check("leading contact is parameter-free at ~0.10 ppm (102-103 ppb)", 100 < abs(ppb(a_lin)) < 105)

    print("\n[2] The residual is O(alpha/2pi): the natural next-order size")
    lead_shift = a_lin - 137.0
    resid = CODATA - a_lin                     # what the 2nd order must supply
    print(f"    leading shift (137 -> )      = {lead_shift:.9f}")
    print(f"    residual to CODATA           = {resid:.3e}   ({ppb(a_lin):+.1f} ppb)")
    ratio = abs(resid) / lead_shift
    print(f"    residual / leading shift     = {ratio:.3e}   (alpha/2pi = {A0/TWO_PI:.3e}; ratio/(alpha/2pi) = {ratio/(A0/TWO_PI):.3f})")
    check("residual is O(alpha/2pi) (0.1--1 x alpha/2pi of the leading shift)", 0.1 < ratio / (A0 / TWO_PI) < 1.0)

    print("\n[3] Decompose the 102.7 ppb: form-choice vs C2 fit")
    a_form = solve_recursive(0.0)              # self-consistent quadratic, C2=0
    a_c2 = solve_recursive(float(C2_CANON))    # + C2=24/7
    print(f"    linear (C2=0, no self-consistency)      = {a_lin:.9f}  ({ppb(a_lin):+.1f} ppb)")
    print(f"    self-consistent form, C2=0              = {a_form:.9f}  ({ppb(a_form):+.1f} ppb)  <- form choice absorbs ~{ppb(a_lin)-ppb(a_form):.0f} ppb")
    print(f"    self-consistent form, C2=24/7           = {a_c2:.9f}  ({ppb(a_c2):+.1f} ppb)  <- C2 fits the last ~{ppb(a_form)-ppb(a_c2):.0f} ppb")
    check("the self-consistent FORM (an ansatz) already absorbs most of the 102.7 ppb", ppb(a_lin) - ppb(a_form) > 40)
    check("adding C2=24/7 reaches the few-ppb regime", abs(ppb(a_c2)) < 6)

    print("\n[4] Is C2=24/7 razor-selected or a fit? (16.3 search-space audit)")
    # exact C2 that lands on CODATA in the recursive form
    lo, hi = 0.0, 20.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if solve_recursive(mid) > CODATA else (lo, mid)
    c2_exact = 0.5 * (lo + hi)
    print(f"    exact C2 for CODATA = {c2_exact:.5f}   (24/7 = {float(C2_CANON):.5f}, off by {100*(float(C2_CANON)/c2_exact-1):+.2f}%)")
    # datum constraint: |C2| pinned only to ~1% (ANCHOR L630). Count low-complexity rationals within 1%.
    band = 0.01
    competitors = sorted({Fraction(p, q) for q in range(1, 13) for p in range(1, 60)
                          if abs(float(Fraction(p, q)) / c2_exact - 1.0) < band}, key=float)
    print(f"    low-complexity rationals (q<=12) within {band*100:.0f}% of C2_exact: {len(competitors)}")
    print(f"      e.g. {[str(c) for c in competitors[:8]]}")
    check("24/7 is within the datum band", C2_CANON in competitors)
    check("C2 is a FIT: many low-complexity rationals fit the 1% datum band (not razor-unique)", len(competitors) >= 4)

    print(
        f"""
[5] VERDICT -- the residual sharpened: derived floor + fitted refinement
    * DERIVED, parameter-free, form-independent: the leading one-contact
      alpha^-1 = 137 + 31*alpha0/(2pi) = {a_lin:.6f}, i.e. CODATA to {abs(ppb(a_lin)):.0f} ppb
      (0.10 ppm), from 137=T(16)+1, N1=31=2*sumQ^2-1 (endpoint-contact theorem),
      and the loop factor 2pi. No fit. This is the genuine dressed-alpha result.
    * The remaining 102.7 ppb is a true O(alpha/2pi) second-order term (verified
      size). Of it, ~{ppb(a_lin)-ppb(a_form):.0f} ppb is absorbed merely by choosing the
      self-consistent quadratic FORM (an ansatz, per canon L630), and C2=24/7
      fits the last ~{ppb(a_form)-ppb(a_c2):.0f} ppb -- but C2 is NOT razor-selected ({len(competitors)}
      low-complexity rationals fit the ~1% datum band). So the 3 ppb match is a
      FIT, not a derivation.
    * OPEN THEOREM (sharpened): derive the second-order term FORM-INDEPENDENTLY
      from the 2nd-order connected endpoint covariance in the monitored Wilson-
      endpoint layer (the natural next cumulant after the 31), reproducing the
      ~100 ppb without the ansatz form or a fitted C2. Only that closes the
      dressed value below 0.10 ppm as a derivation.

    Bottom line: dressed alpha is now honestly PARAMETER-FREE to 0.10 ppm
    (leading contact, derived); the sub-0.1-ppm refinement to 3 ppb rests on an
    ansatz form + a fitted C2, and the closure target is the form-independent
    second-order endpoint covariance.
exit 0"""
    )
    print("ALL CHECKS PASSED" if ok else "CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
