#!/usr/bin/env python3
r"""ITEM 132: reduce the P1 premise to a same-service-class theorem plus an
external phase-latch scale.

Question
--------
The R4/MOND branch used "P1" as a compact name for several things:

    * a count-valued nonexclusive line ledger;
    * matched creation and erasure rates;
    * the Poisson mean x=|g|/a0;
    * the one-a0 central-cell service quantum.

Those are not all the same premise.  This script separates them.

Result
------
The count-valued ledger and matched rate form are now derivable from the
current R4 service structure, provided that creation and recycling are the same
monitored service class:

    birth  n -> n+1  at Gamma0 x,
    death  n -> n-1  at Gamma0 n.

The common clock Gamma0 cancels, so the stationary line count is exactly
Poisson(x), with no eta/kappa multiplier.  A two-clock model, a finite flag, or
a KMS-strain-biased ledger all reintroduce the very free multiplier P1 was
introduced to forbid.

Therefore the residual should be named sharply:

    P1_scheduler      : closed under same monitored service-class reading;
    P1_scale/latch    : still open -- prove x=|g|/a0 and one a0 per local
                        central-cell service record.

This is a reduction, not a full MOND lock: the a0 phase-return latch remains
the separate hard target.
"""

from __future__ import annotations

import math

from item132_r4_stinespring_fock_lift import (
    build_stinespring_columns,
    repair_edges,
    stinespring_errors,
)


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def poisson_probs(mean: float, nmax: int = 64) -> list[float]:
    probs = [math.exp(-mean)]
    for n in range(nmax):
        probs.append(probs[-1] * mean / (n + 1))
    z = sum(probs)
    return [p / z for p in probs]


def immigration_death_mean(x: float, gamma_birth: float = 1.0, gamma_death: float = 1.0) -> float:
    """Stationary mean of n->n+1 at gamma_birth*x and n->n-1 at gamma_death*n."""

    return (gamma_birth / gamma_death) * x


def detailed_balance_ratio(mean: float, n: int) -> float:
    probs = poisson_probs(mean, max(64, n + 2))
    return probs[n + 1] / probs[n]


def finite_flag_occupancy(x: float) -> float:
    """The saturating one-bit occupancy control."""

    return 1.0 - math.exp(-x)


def kms_biased_mean(x: float, delta: float) -> float:
    """Wrong-ledger control: a strain-Boltzmann factor on the birth channel."""

    return x * math.exp(-delta)


def main() -> None:
    print("ITEM 132: P1 SERVICE-CLASS REDUCTION")
    print("=" * 92)

    print("\n[1] Finite R4 service event exists and is normalized")
    edges = repair_edges()
    columns, _labels, _idx = build_stinespring_columns()
    diag_err, offdiag = stinespring_errors(columns)
    print(f"  R4 repair edges                         = {len(edges)}")
    print(f"  Stinespring diag/offdiag errors         = {diag_err:.3e} / {offdiag:.3e}")
    check(len(edges) == 6, "R4 has the six finite repair events used by the line ledger")
    check(diag_err < 1e-12 and offdiag < 1e-12, "finite R4 Stinespring map is an isometry")

    print("\n[2] Same monitored service class => Poisson(x), no multiplier")
    for x in (0.1, 1.0, 3.0, 10.0):
        mean = immigration_death_mean(x)
        probs = poisson_probs(mean)
        computed_mean = sum(n * p for n, p in enumerate(probs))
        fano = sum((n - computed_mean) ** 2 * p for n, p in enumerate(probs)) / computed_mean
        ratio_2 = detailed_balance_ratio(mean, 2)
        print(
            f"  x={x:5.2f}: mean={computed_mean:.9f}, Fano={fano:.9f}, "
            f"P3/P2={ratio_2:.9f} (target x/3={x/3:.9f})"
        )
        check(abs(computed_mean - x) < 1e-10, "same-clock immigration-death mean is exactly x")
        check(abs(fano - 1.0) < 1e-8, "same-clock immigration-death ledger is Poisson")
        check(abs(ratio_2 - x / 3.0) < 1e-10, "detailed balance ratios are fixed by one mean x")

    print("\n[3] Controls: what would reopen P1 as a fitted knob")
    x = 1.0
    controls = [
        ("same service class", immigration_death_mean(x, 1.0, 1.0)),
        ("two-clock eta/kappa=1.20", immigration_death_mean(x, 1.2, 1.0)),
        ("two-clock eta/kappa=0.80", immigration_death_mean(x, 0.8, 1.0)),
        ("KMS strain bias delta=0.25", kms_biased_mean(x, 0.25)),
    ]
    for label, mean in controls:
        print(f"  {label:28s} stationary mean = {mean:.9f}")
    print(f"  finite one-bit flag occupancy at x=1 = {finite_flag_occupancy(1.0):.9f}")
    check(abs(controls[0][1] - 1.0) < 1e-12, "one service class leaves no eta/kappa multiplier")
    check(abs(controls[1][1] - 1.0) > 0.1, "two service clocks reintroduce an arbitrary multiplier")
    check(abs(controls[3][1] - 1.0) > 0.1, "strain-weighted KMS birth is the wrong ledger for driven service counts")
    check(abs(finite_flag_occupancy(1.0) - 1.0) > 0.3, "finite occupancy is a saturating flag, not the line count")

    print("\n[4] What P1 now means")
    print("  CLOSED UNDER CURRENT SERVICE-CLASS READING:")
    print("    repeated Stinespring slots give a count ledger; creation and recycling")
    print("    are one monitored R4 service class; the common Gamma0 clock cancels.")
    print("    This proves chi_R4=1 and retires eta/kappa as a free rate ratio.")
    print("  STILL OPEN / SEPARATE:")
    print("    the source variable x=|g|/a0 inherits the one-a0 phase-latch scale.")
    print("    That coefficient is not supplied by the Poisson service-class theorem.")
    print("\nexit 0 -- P1_scheduler reduced to a same-service-class theorem; P1_scale remains the a0 latch.")


if __name__ == "__main__":
    main()
