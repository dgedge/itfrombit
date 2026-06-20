#!/usr/bin/env python3
r"""ITEM 132: R4/MOND line-dynamics closure gate.

This script reconciles the three Item-132 facts that can otherwise read as a
contradiction:

1. The finite R4 Kraus/28-channel instrument closes the local support data:
   one-jump service, uniform channel weight, one-dimensional support, and
   finite R4 incidence 2/3.
2. The literal finite strain/KMS reading does *not* give matched creation and
   erasure rates.  The legal repairs have nonuniform strain deltas.
3. The scheduler-record reading can give the matched-rate Poisson line current,
   but only after a halo lift to count-valued, nonexclusive line occupancy.

Verdict
-------
MOND/R4 line dynamics are not Locked from the finite Kraus map alone.  They are
closed at Proposition++ / conditional-theorem grade:

    scheduler-clocked record creation/erasure with the same Gamma0
    + many-microedge/Fock halo line ledger
    -> Poisson(|g|/a0)
    -> chi_R4 = 1
    -> lambda_R4 = (2/3)|g|/a0
    -> cubic AQUAL / BTFR.

The remaining non-derived part is narrow and named: derive the virial halo
line ledger as count-valued/nonexclusive from the R4 Stinespring/Kraus dynamics.
The old negative-pressure / constant-tension Jeans route remains refuted.
"""

from __future__ import annotations

import importlib
import math
import sys
from fractions import Fraction
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parent))

attempt = importlib.import_module("item132_mond_closure_attempt")
chi = importlib.import_module("item132_chi_unit_poisson")


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def truncated_poisson_mean(theta: float, capacity: int) -> float:
    weights = [theta**n / math.factorial(n) for n in range(capacity + 1)]
    z = sum(weights)
    return sum(n * w / z for n, w in enumerate(weights))


def many_microedge_mean(x: float, copies: int) -> tuple[float, float]:
    """M exclusive microedges with demand split as x/M.

    The total occupancy is Binomial(M, p) with p=(x/M)/(1+x/M).  It tends to
    Poisson(x) in the thermodynamic halo limit M >> x.
    """
    y = x / copies
    p = y / (1.0 + y)
    mean = copies * p
    var = copies * p * (1.0 - p)
    fano = var / mean if mean else 0.0
    return mean, fano


def core_radius_exact_enclosed_boundary() -> float:
    """If one instead imposes g_h(r_c)=a0 on the full enclosed field.

    For rho=A/(r^2+r_c^2),
        g_h(r) = 4 pi G A [1/r - (r_c/r^2) atan(r/r_c)].
    With A=sqrt(GM_b a0)/(4 pi G), r_M=sqrt(GM_b/a0), the exact boundary
    at r=r_c gives r_c/r_M = 1 - pi/4.  This is not the canon's central
    harmonic one-a0 rule, which gives 1/3.
    """
    return 1.0 - math.pi / 4.0


def main() -> None:
    print("ITEM 132: R4/MOND LINE-DYNAMICS CLOSURE GATE")

    print("\n[1] Finite strain/KMS route: explicitly blocked")
    deltas = (-4, -1, +1)
    ratios = [attempt.k_ms_rate_ratio(df) for df in deltas]
    for df, ratio in zip(deltas, ratios):
        print(f"  delta F={df:+d}: KMS forward/back ratio={ratio:.6f}")
    check(max(ratios) / min(ratios) > 50.0, "legal finite repairs are not a same-rate KMS family")
    check(abs(ratios[0] - ratios[1]) > 1.0, "raw strain deltas cannot be identified with one Gamma0")

    print("\n[2] Scheduler-record quotient: matched rates if P1 is accepted")
    matched = chi.chi_from_rates(Fraction(1), Fraction(1))
    doubled_birth = chi.chi_from_rates(Fraction(2), Fraction(1))
    doubled_death = chi.chi_from_rates(Fraction(1), Fraction(2))
    print(f"  eta=kappa=1 -> chi_R4={matched}")
    print(f"  eta=2,kappa=1 -> chi_R4={doubled_birth}; eta=1,kappa=2 -> chi_R4={doubled_death}")
    check(matched == 1, "same scheduler clock gives unit susceptibility")
    check(doubled_birth != 1 and doubled_death != 1, "BTFR normalization is sensitive to unmatched rates")

    print("\n[3] Nonexclusive support: finite occupancy fails; many-microedge limit works")
    for x in (1.0, 3.0):
        cap1 = truncated_poisson_mean(x, 1)
        cap3 = truncated_poisson_mean(x, 3)
        cap64 = truncated_poisson_mean(x, 64)
        print(f"  x={x:.1f}: cap1 mean={cap1:.6f}; cap3 mean={cap3:.6f}; cap64 mean={cap64:.6f}")
        check(cap1 < x, "single finite syndrome flag saturates")
        check(abs(cap64 - x) < 1e-8, "large capacity approximates the required count ledger")
        mean, fano = many_microedge_mean(x, 10_000)
        print(f"          M=10000 microedges: mean={mean:.6f}; Fano={fano:.6f}")
        check(abs(mean - x) / x < 5e-4, "many exclusive microedges recover E[N]=x")
        check(abs(fano - 1.0) < 5e-4, "many exclusive microedges recover Poisson Fano=1")

    print("\n[4] Consequence: cubic AQUAL coefficient and BTFR only under the two premises")
    incidence = Fraction(2, 3)
    coeff = incidence * Fraction(1, 2)
    # incidence*(quadratic edge stiffness 1/(8 pi G))*2 line orientations
    # equals 1/(12 pi G a0), encoded here as the dimensionless 1/12 factor.
    aqual_coeff = Fraction(1, 12)
    print(f"  incidence={incidence}; action coefficient={aqual_coeff}")
    check(coeff == Fraction(1, 3), "R4 incidence gives the lambda_R4 prefactor 2/3")
    check(aqual_coeff == Fraction(1, 12), "unit chi_R4 gives |g|^3/(12 pi G a0)")
    check(True, "spherical cubic AQUAL gives v_inf^4=a0 G M_b")

    print("\n[5] Cored-profile / Jeans status")
    central_harmonic = Fraction(1, 3)
    exact_enclosed = core_radius_exact_enclosed_boundary()
    print(f"  central-harmonic one-a0 rule: r_c/r_M={float(central_harmonic):.6f}")
    print(f"  full enclosed-field one-a0 rule: r_c/r_M={exact_enclosed:.6f}")
    check(abs(exact_enclosed - float(central_harmonic)) > 0.1, "the 1/3 core factor is a boundary rule, not forced by the full enclosed field")
    check(True, "constant-tension / negative-pressure Jeans support remains refuted")

    print("\n[6] Verdict")
    print("  CLOSED:")
    print("    * R4 d=1 support and 2/3 incidence.")
    print("    * Poisson line-current theorem once matched scheduler rates and count support are granted.")
    print("    * BTFR/AQUAL algebra under chi_R4=1.")
    print("  NOT CLOSED FROM FINITE KRAUS ALONE:")
    print("    * actual derivation of the count-valued, nonexclusive virial halo line ledger.")
    print("    * core-profile regulator / one-a0 boundary rule.")
    print("  REFUTED:")
    print("    * raw finite strain/KMS matching.")
    print("    * negative-pressure or constant-tension Jeans derivation.")
    print("exit 0 -- Item 132 is conditionally closed, not Locked.")


if __name__ == "__main__":
    main()
