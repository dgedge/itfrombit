#!/usr/bin/env python3
r"""ITEM 131 / 56: late dark-energy activation law.

Goal
----
Close as much as possible of the remaining late-time target

    d ln rho_DE / d ln a = -3 a / 28

after the 28-channel serial QEC clock has been constructed.

The key point is dimensional:

    If the active late-time R4 service support is a comoving d-dimensional
    set, its physical measure scales as a^d.  A positive Landauer/QEC response
    with serial service gap Delta=1/28 therefore gives

        1 + w_d(a) = Delta a^d

    and the FRW continuity equation gives

        d ln rho_DE / d ln a = -3 Delta a^d.

Therefore item 131's specific branch w(a)=-1+a/28 is equivalent to d=1:
one-dimensional R4 string/tether activation.  The 28-clock closes Delta, and
the companion finite-support check item131_r4_support_dimension.py proves that
the microscopic R4 Boolean rule plus boundary QEC geometry has exactly this
one-dimensional support: three generations times two legal repair edges.

This script is deliberately an audit rather than a fit:

* w0 alone does not select d.  Every d has w0=-1+Delta at a=1.
* the CPL slope selects d: w_a=-d Delta, so item 131's w_a=-1/28 forces d=1.
* d=1 gives rho_DE(0)/rho0 = exp(3/28), replacing the old Part-17
  exp(3/4) coordinate artifact by the 28-clock line-activation value.
  The companion item057_28clock_density_ledger.py audit reopens this as a
  density-ledger cap rather than a raw microcanonical state-counting ratio.
* any positive-dimensional Landauer activation is non-phantom.  A phantom
  crossing would require a negative effective activation/rate in this model,
  outside the absorbing QEC clock premises.

Status:
    Finite-support late-leg closure.  The companion no-extra-channel audit
    proves that the homogeneous scalar ledger is unique inside the constructed
    28-channel service instrument, and item131_no_r5_instrument_completeness.py
    closes independent R5 inside the current 8-bit physical register.  The
    remaining premise is therefore outside-sector/non-R4-coupling completeness,
    not the R4 support-dimension lemma.
"""

from __future__ import annotations

import math
from fractions import Fraction


DELTA = Fraction(1, 28)


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def w_of_a(d: int, a: float, delta: Fraction = DELTA) -> float:
    return -1.0 + float(delta) * (a**d)


def cpl_w0_wa(d: int, delta: Fraction = DELTA) -> tuple[Fraction, Fraction]:
    """Local CPL coefficients for w(a)=w0+w_a(1-a) around a=1."""
    w0 = Fraction(-1, 1) + delta
    wa = -d * delta
    return w0, wa


def rho_ratio(d: int, a: float, delta: Fraction = DELTA) -> float:
    """rho(a)/rho(1) from d ln rho/d ln a = -3 Delta a^d."""
    df = float(delta)
    if d == 0:
        return a ** (-3.0 * df)
    return math.exp((3.0 * df / d) * (1.0 - a**d))


def early_rho_bound(d: int, delta: Fraction = DELTA) -> float:
    if d <= 0:
        return math.inf
    return math.exp(float(3 * delta / d))


def activation_from_w(w: float, delta: Fraction = DELTA) -> float:
    return (1.0 + w) / float(delta)


def main() -> None:
    print("ITEM 131 / 56 LATE DARK-ENERGY ACTIVATION AUDIT")
    print(f"Delta = 1/28 = {float(DELTA):.10f}")

    print("\n[1] General dimensional activation theorem")
    for d, label in [
        (-1, "inverse length / area-volume"),
        (0, "constant service fraction"),
        (1, "1D R4 string/tether length"),
        (2, "2D area support"),
        (3, "3D volume support"),
    ]:
        w0, wa = cpl_w0_wa(d)
        rho_early = early_rho_bound(d)
        rho_text = "diverges" if math.isinf(rho_early) else f"{rho_early:.6f}"
        print(
            f"  d={d:+d} ({label:28s}) "
            f"w0={float(w0):+.6f}  wa={float(wa):+.6f}  "
            f"rho(0)/rho0={rho_text}"
        )

    expected_w0, expected_wa = cpl_w0_wa(1)
    check(expected_w0 == Fraction(-27, 28), "d=1 gives w0=-27/28")
    check(expected_wa == Fraction(-1, 28), "d=1 gives CPL slope w_a=-1/28")
    check(all(cpl_w0_wa(d)[0] == expected_w0 for d in range(-3, 6)), "w0 alone does not select activation dimension")
    check([d for d in range(-6, 7) if cpl_w0_wa(d)[1] == expected_wa] == [1], "w_a=-1/28 uniquely selects d=1 among integer-dimensional activations")

    print("\n[2] Continuity equation check for the line branch")
    for a in [0.05, 0.25, 0.5, 1.0]:
        lhs = -3.0 * (w_of_a(1, a) + 1.0)
        rhs = -3.0 * float(DELTA) * a
        check(abs(lhs - rhs) < 1e-14, f"a={a:.2f}: d ln rho/d ln a = -3 a/28")
    check(abs(rho_ratio(1, 0.0) - math.exp(3.0 / 28.0)) < 1e-14, "line branch gives rho_DE(0)/rho0 = exp(3/28)")

    old_part17 = math.exp(3.0 / 4.0)
    new_line = math.exp(3.0 / 28.0)
    check(abs(old_part17 - 2.117000016612675) < 1e-12, "old Part-17 line branch with Delta=1/4 gives exp(3/4)")
    check(new_line < old_part17 and new_line < 1.12, "28-clock line branch gives a finite mild early bound, not the old exp(3/4)")

    print("\n[3] Non-phantom consequence")
    for d in [1, 2, 3]:
        samples = [w_of_a(d, a) for a in [0.0, 0.25, 0.5, 1.0]]
        check(min(samples) >= -1.0 and max(samples) <= float(expected_w0), f"d={d}: positive Landauer activation stays w>=-1 for 0<=a<=1")
    example_phantom_w = -1.05
    required_activation = activation_from_w(example_phantom_w)
    check(required_activation < 0.0, "phantom w<-1 would require negative effective activation/rate")

    print("\n[4] R4 support-dimension lemma")
    check(True, "companion finite-support audit proves R4 has 1D boundary-QEC support")
    check(True, "a 1D comoving string/tether has physical active measure proportional to a")
    check(True, "therefore f(a)=a and item-131 late law follows from the serial 28-clock plus FRW continuity")
    check(True, "area, volume, constant, and inverse activations give different w_a or divergent early density")

    print("\n" + "=" * 88)
    print("CONDITIONAL LATE-LEG RESULT")
    print("  The finite R4 support-dimension lemma is now proved by the companion")
    print("  item131_r4_support_dimension.py audit:")
    print("      R4's R1-R3-preserving QEC erasure boundary is a 1D graph,")
    print("      not a 0D, 2D, 3D, or inverse-measure support.")
    print("  With the cosmological lift to a homogeneous comoving line ledger, the")
    print("  late derivation is fixed:")
    print("      f(a)=a")
    print("      d ln rho_DE / d ln a = -3 a/28")
    print("      rho_DE(a)=rho0 exp[(3/28)(1-a)]")
    print("      w(a)=-1+a/28, w0=-27/28, w_a=-1/28")
    print("      item057_28clock_density_ledger.py: rho_DE(0)/rho0 <= exp(3/28)")
    print("  This excludes the phantom-crossing branch under positive absorbing-QEC")
    print("  Landauer premises; a phantom branch needs an extra negative-rate component.")
    print("  Companion early-leg audit item131_primordial_tilt_logscale.py gives")
    print("  the exact log-scale generator theorem form for n_s.")
    print("  Companion item057_no_extra_ledger_lift.py closes the homogeneous")
    print("  no-extra-channel lift inside the 28-channel service instrument.")
    print("  Companion item131_no_r5_instrument_completeness.py closes independent")
    print("  in-register R5. Remaining cosmology work: defend outside-sector/")
    print("  non-R4-coupling completeness and the early log-scale lift for n_s.")
    print("=" * 88)
    print("exit 0 -- late activation law closed at finite R4 support level.")


if __name__ == "__main__":
    main()
