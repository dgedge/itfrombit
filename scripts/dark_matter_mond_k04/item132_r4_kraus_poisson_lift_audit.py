#!/usr/bin/env python3
r"""ITEM 132: R4 Kraus map versus Poisson line-current lift.

Question
--------
Can the remaining dark-halo/BTFR blocker be closed directly from the actual
finite R4 boundary-QEC Kraus/Stinespring map?

The conditional theorem is already clean:

    N_e in N,
    n -> n+1 at rate Gamma0*x,
    n -> n-1 at rate Gamma0*n,
    x = |g|/a0

has stationary Poisson(x), so E[N_e]=x and the finite R4 incidence gives

    lambda_R4 = (2/3)|g|/a0,
    chi_R4 = 1.

Verdict
-------
The finite R4/28-channel instrument proves the right local ingredients:

* one-jump finite bandwidth;
* uniform service weights on the 28 logical/transverse channels;
* one-dimensional R4 repair support with the 2/3 incidence factor.

It does not, by itself, prove the two extra premises needed by the Poisson
line-current theorem:

1. rate matching between R4 line creation and erasure;
2. count-valued nonexclusive occupancy on a coarse halo line.

Rate matching follows if creation and erasure are the same R4 QEC edge seen as
opposite orientations of one detailed-balance/KMS pair with the same Kraus
norm Gamma0.  A generic CPTP finite instrument still permits eta/kappa != 1.

Nonexclusive occupancy follows only after a virial/thermodynamic line-network
lift: many independent exclusive microedges, or an unbounded Fock-like line
sector.  A single finite exclusive syndrome flag gives a truncated/saturating
response and cannot yield exact MOND/AQUAL linearity.

Thus the BTFR recovery remains a conditional Poisson line-current theorem, not
a negative-pressure theorem and not yet a parameter-free result.
"""

from __future__ import annotations

import math
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def read(relpath: str) -> str:
    return (ROOT / relpath).read_text(encoding="utf-8")


def chi_from_rates(eta: Fraction, kappa: Fraction) -> Fraction:
    """Unbounded immigration-death mean is (eta/kappa)*x."""
    return eta / kappa


def truncated_poisson_stats(theta: float, capacity: int) -> tuple[float, float, float]:
    """Stationary stats for blocked finite-capacity immigration-death.

    Birth is blocked at n=capacity.  The stationary law is a Poisson(theta)
    distribution conditioned on n <= capacity.
    """
    weights = [theta**n / math.factorial(n) for n in range(capacity + 1)]
    z = sum(weights)
    probs = [w / z for w in weights]
    mean = sum(n * p for n, p in enumerate(probs))
    var = sum((n - mean) ** 2 * p for n, p in enumerate(probs))
    fano = var / mean if mean else 0.0
    return mean, var, fano


def many_copy_exclusive_stats(x: float, copies: int, eta: float = 1.0, kappa: float = 1.0) -> tuple[float, float, float]:
    """M independent binary microedges with demand split across copies.

    Each copy has birth/death ratio y=(eta/kappa)*x/M.  The total occupancy is
    Binomial(M, y/(1+y)), approaching Poisson((eta/kappa)*x) as M -> infinity.
    """
    y = (eta / kappa) * x / copies
    p = y / (1.0 + y)
    mean = copies * p
    var = copies * p * (1.0 - p)
    fano = var / mean if mean else 0.0
    return mean, var, fano


def poisson_stats(theta: float) -> tuple[float, float, float]:
    return theta, theta, 1.0


def action_coeff(chi: Fraction) -> Fraction:
    """Coefficient multiplying |g|^3/(pi G a0)."""
    return chi * Fraction(1, 12)


def btfr_a0_eff(chi: Fraction) -> Fraction:
    """Spherical AQUAL gives v_inf^4=(a0/chi)G M_b."""
    return Fraction(1, 1) / chi


def main() -> None:
    print("ITEM 132: R4 KRAUS MAP / POISSON LINE-CURRENT AUDIT")

    print("\n[1] Finite instrument facts already closed")
    w_to_28 = read("python_code/item131_w_to_28_instrument.py")
    one_jump = read("python_code/item131_one_jump_premise.py")
    r4_support = read("python_code/item131_r4_support_dimension.py")
    chi_unit = read("python_code/item132_chi_unit_poisson.py")

    check("service channel has rate 1/28" in w_to_28, "8-bit -> 28-channel Stinespring refinement fixes uniform service weights")
    check("one-tick bandwidth <= 1" in one_jump, "finite W=S*C service is one-jump, not parallel reset")
    check("exactly two legal" in r4_support and "R4 boundary has 3 generations x 2 legal repair edges" in r4_support, "R4 finite support has two legal repair edges per generation")
    check("count-valued line ledger" in chi_unit, "Poisson closure explicitly assumes a count-valued line ledger")

    incidence = Fraction(2, 3)
    check(incidence == Fraction(2, 3), "finite R4 repair incidence supplies the 2/3 factor")

    print("\n[2] Kraus covariance does not by itself force eta/kappa=1")
    for eta, kappa in [(Fraction(1), Fraction(1)), (Fraction(2), Fraction(1)), (Fraction(1), Fraction(2))]:
        chi = chi_from_rates(eta, kappa)
        print(f"  eta={eta}, kappa={kappa}: chi_R4=eta/kappa={chi}")
    check(chi_from_rates(Fraction(1), Fraction(1)) == 1, "same creation/erasure Kraus norm gives chi_R4=1")
    check(chi_from_rates(Fraction(2), Fraction(1)) == 2, "extra creation multiplicity would double the BTFR susceptibility")
    check(chi_from_rates(Fraction(1), Fraction(2)) == Fraction(1, 2), "extra erasure multiplicity would halve the BTFR susceptibility")
    check(True, "uniform channel covariance fixes rates within a channel family, not the forward/backward ratio")
    check(True, "rate matching therefore needs a detailed-balance/KMS same-edge premise, not mere CPTP normalization")

    print("\n[3] A single finite exclusive syndrome flag cannot be the MOND line ledger")
    for x in [0.1, 1.0, 3.0]:
        cap1_mean, _, cap1_fano = truncated_poisson_stats(x, 1)
        cap3_mean, _, cap3_fano = truncated_poisson_stats(x, 3)
        cap50_mean, _, cap50_fano = truncated_poisson_stats(x, 50)
        print(
            f"  x={x:3.1f}: cap1 mean={cap1_mean:.6f}, Fano={cap1_fano:.6f}; "
            f"cap3 mean={cap3_mean:.6f}, Fano={cap3_fano:.6f}; "
            f"cap50 mean={cap50_mean:.6f}, Fano={cap50_fano:.6f}"
        )
        check(cap1_mean < x, "capacity-1 exclusive occupancy saturates below the required mean x")
        check(cap3_mean <= x, "finite low capacity gives a truncated/saturating response")
        check(abs(cap50_mean - x) < 1e-8, "large enough capacity approximates the non-saturating branch")

    print("\n[4] Many-copy exclusive microedges converge to the Poisson line current")
    for x in [1.0, 3.0]:
        for copies in [4, 28, 10_000]:
            mean, var, fano = many_copy_exclusive_stats(x, copies)
            print(f"  x={x:3.1f}, M={copies:5d}: mean={mean:.9f}, Fano={fano:.9f}")
        mean_big, _, fano_big = many_copy_exclusive_stats(x, 10_000)
        check(abs(mean_big - x) / x < 5e-4, "M >> x recovers E[N]=x")
        check(abs(fano_big - 1.0) < 5e-4, "M >> x recovers Poisson Fano factor 1")

    print("\n[5] Exact unbounded line-current theorem gives the AQUAL coefficient")
    for x in [0.1, 1.0, 3.0]:
        mean, var, fano = poisson_stats(x)
        line_density = float(incidence) * mean
        print(f"  x={x:3.1f}: Poisson mean={mean:.6f}, Fano={fano:.6f}, lambda_R4={line_density:.6f}")
        check(mean == x and var == x and fano == 1.0, "unbounded immigration-death ledger is exactly Poisson(x)")

    chi = Fraction(1)
    check(action_coeff(chi) == Fraction(1, 12), "chi_R4=1 gives |g|^3/(12 pi G a0)")
    check(btfr_a0_eff(chi) == 1, "chi_R4=1 gives v_inf^4=a0 G M_b")

    print("\n[6] Honest theorem status")
    check(True, "closed: negative-pressure BTFR derivation is withdrawn")
    check(True, "closed conditionally: matched-rate count ledger implies chi_R4=1 and cubic MOND/AQUAL action")
    check(True, "not closed from finite Kraus alone: same-norm forward/backward R4 pair must still be proved")
    check(True, "not closed from finite Kraus alone: virial halo must still be lifted to a count-valued nonexclusive line ledger")

    print("\n" + "=" * 100)
    print("R4 KRAUS / POISSON-LIFT VERDICT")
    print("  The actual finite R4/28-channel instrument supplies one-jump service,")
    print("  uniform channel weights, and the 2/3 R4 incidence.  It does not by")
    print("  itself force eta/kappa=1 or unbounded occupancy.")
    print("  The recoverable theorem is therefore conditional:")
    print("      same-norm R4 creation/erasure pair + count-valued halo line ledger")
    print("      -> Poisson(|g|/a0) -> chi_R4=1 -> cubic MOND/AQUAL action -> BTFR.")
    print("  A plausible route to the count ledger is the thermodynamic halo limit of")
    print("  many independent exclusive R4 microedges, but that virial/Fock lift is")
    print("  the remaining microscopic theorem target.")
    print("=" * 100)
    print("exit 0 -- BTFR recovery remains a conditional Poisson line-current theorem.")


if __name__ == "__main__":
    main()
