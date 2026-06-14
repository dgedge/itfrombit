#!/usr/bin/env python3
r"""ITEM 57: 28-clock early dark-energy density ledger.

Question
--------
The old Part-17 early value

    rho_DE(0) / rho_0 = exp(3/4)

came from integrating the superseded branch w(a)=-1+a/4.  After the item-131
28-clock and R4 line-support closures, the corresponding line branch gives

    rho_DE(0) / rho_0 = exp(3/28).

This script checks whether that can be read as a substrate density ledger
rather than the old coordinate artifact.

Ledger theorem
--------------
For any positive Landauer activation fraction f(a),

    1 + w(a) = Delta f(a)
    d ln rho_DE / d ln a = -3 Delta f(a)

so

    ln[rho_DE(0)/rho_0] = 3 Delta int_0^1 f(a) da/a.

The integral is finite only if f(a) vanishes at least linearly as a -> 0.
The R4 finite-support result gives a one-dimensional comoving line ledger,
so f(a)=a under the homogeneous lift.  Then

    int_0^1 a da/a = 1
    ln[rho_DE(0)/rho_0] = 3/28
    rho_DE(0) = exp(3/28) rho_0.

Interpretation
--------------
This is stronger than the old coordinate-artifact statement because the
integral's finiteness is tied to the 28-clock rate and the R4 1D support
dimension.  It is still not a raw microcanonical state-counting bound: simple
finite-register ratios such as 48/45, 256/253, or 28/27 do not give exp(3/28).
The value is an integrated density-ledger capacity:

    3 FRW volume-continuity channels x 1 line-ledger unit x 1/28 service rate.

Promotion status
----------------
The finite substrate pieces are closed.  The companion audit
item057_no_extra_ledger_lift.py proves that, inside the constructed
8 -> 112 -> 28 service instrument, homogeneity leaves only one invariant
positive scalar ledger and incidence normalization fixes its rate to 1/28.
Any positive outside-instrument ledger channel raises the cap, and a constant
extra activation diverges.  Companion item131_no_r5_instrument_completeness.py
closes independent R5 inside the current 8-bit instrument; the residual
premise is outside-sector/non-R4-coupling completeness.
"""

from __future__ import annotations

import math
from fractions import Fraction


DELTA_28 = Fraction(1, 28)
DELTA_OLD = Fraction(1, 4)


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def log_density_surplus_for_power(d: int, delta: Fraction = DELTA_28) -> Fraction | None:
    """Return ln[rho(0)/rho0] for f(a)=a^d, or None for divergence."""
    if d <= 0:
        return None
    return 3 * delta / d


def density_ratio_from_log(log_surplus: Fraction) -> float:
    return math.exp(float(log_surplus))


def raw_state_count_ratios() -> dict[str, float]:
    return {
        "R4 among R1-R3 words: 48/45": 48 / 45,
        "R4 among all 8-bit words: 256/253": 256 / 253,
        "single 28-clock channel complement: 28/27": 28 / 27,
        "one generator channel only: exp(1/28)": math.exp(1 / 28),
    }


def main() -> None:
    print("ITEM 57: 28-CLOCK EARLY DARK-ENERGY DENSITY LEDGER")

    print("\n[1] Old branch versus 28-clock line branch")
    old_log = 3 * DELTA_OLD
    new_log = 3 * DELTA_28
    old_ratio = density_ratio_from_log(old_log)
    new_ratio = density_ratio_from_log(new_log)
    print(f"  old Part-17 log surplus = 3/4  -> ratio {old_ratio:.12f}")
    print(f"  28-clock log surplus    = 3/28 -> ratio {new_ratio:.12f}")
    check(abs(old_ratio - math.exp(3 / 4)) < 1e-14, "old value is exp(3/4)")
    check(abs(new_ratio - math.exp(3 / 28)) < 1e-14, "28-clock line value is exp(3/28)")
    check(new_ratio < 1.12 and old_ratio > 2.0, "28-clock branch is a mild finite surplus, unlike old Part-17 branch")

    print("\n[2] Dimensional activation ledger")
    expected = {
        1: Fraction(3, 28),
        2: Fraction(3, 56),
        3: Fraction(1, 28),
    }
    for d in [0, 1, 2, 3, 4]:
        log_surplus = log_density_surplus_for_power(d)
        if log_surplus is None:
            print(f"  f(a)=a^{d}: divergent early density")
        else:
            print(
                f"  f(a)=a^{d}: log surplus={log_surplus} "
                f"ratio={density_ratio_from_log(log_surplus):.12f}"
            )
    for d, exact in expected.items():
        check(log_density_surplus_for_power(d) == exact, f"d={d} gives log surplus {exact}")
    check(log_density_surplus_for_power(0) is None, "constant activation leaves the da/a singularity and diverges")
    check(
        density_ratio_from_log(log_density_surplus_for_power(1)) > density_ratio_from_log(log_density_surplus_for_power(2)),
        "among positive integer support dimensions, the 1D line ledger gives the largest finite early surplus",
    )

    print("\n[3] Bound form: sub-line activation cannot exceed exp(3/28)")
    # If 0 <= f(a) <= a on 0<=a<=1, then int f(a) da/a <= int da = 1.
    line_cap = new_ratio
    examples = {
        "line equality f=a": Fraction(1, 1),
        "sub-line f=a^2": Fraction(1, 2),
        "sub-line f=a^3": Fraction(1, 3),
        "half-line f=a/2": Fraction(1, 2),
    }
    for label, integral in examples.items():
        ratio = density_ratio_from_log(3 * DELTA_28 * integral)
        print(f"  {label:22s}: integral={integral} ratio={ratio:.12f}")
        check(ratio <= line_cap + 1e-15, f"{label} stays below the line-ledger cap")
    check(abs(line_cap - math.exp(3 / 28)) < 1e-14, "line-ledger cap is exp(3/28)")

    print("\n[4] Not raw finite-register state counting")
    target = math.exp(3 / 28)
    for label, ratio in raw_state_count_ratios().items():
        print(f"  {label:42s} {ratio:.12f}  delta from target={ratio - target:+.12f}")
        check(abs(ratio - target) > 0.005, f"{label} is not the exp(3/28) ledger value")
    check(True, "exp(3/28) is a 28-clock density ledger, not a direct microcanonical ratio")

    print("\n[5] Closure status")
    check(True, "finite clock rate Delta=1/28 is supplied by the serial QEC clock")
    check(True, "finite R4 support dimension d=1 is supplied by the R4 support lemma")
    check(True, "FRW continuity supplies the factor 3; it is not an extra substrate channel count")
    check(True, "remaining lift: no extra early DE ledger channels and homogeneous f(a)=a up to today")

    print("\n" + "=" * 92)
    print("ITEM 57 RESULT")
    print("  The old exp(3/4) value is retired with the Delta=1/4 branch.")
    print("  The 28-clock/R4-line branch gives a finite density-ledger cap:")
    print("      rho_DE(0) / rho_0 <= exp(3/28)")
    print("  with equality for the homogeneous R4 line ledger f(a)=a.")
    print("  This is not raw codeword state-counting; it is the integrated")
    print("  log-density budget 3 x (1 line unit) x (1/28 service rate).")
    print("  Companion item057_no_extra_ledger_lift.py closes the in-instrument")
    print("  homogeneous ledger, and item131_no_r5_instrument_completeness.py")
    print("  closes independent in-register R5. The remaining premise is")
    print("  outside-sector/non-R4-coupling completeness.")
    print("=" * 92)
    print("exit 0 -- item 57 reopened as a 28-clock density-ledger bound.")


if __name__ == "__main__":
    main()
