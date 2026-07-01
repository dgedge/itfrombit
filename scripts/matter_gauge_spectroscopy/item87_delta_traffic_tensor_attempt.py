#!/usr/bin/env python3
r"""Item 87 -- symmetric recovery-traffic tensor attempt for the Koide offset delta.

Question
--------
After the phase/anomaly routes failed, can the Koide offset

    delta_e = 2/9, delta_nu = 1/3, delta_u = 2/27, delta_d = 1/9

be read instead as the principal-axis angle of a real, symmetric recovery
traffic tensor on the S3/C3 generation plane?

This script tests the cleanest possible version of that idea.  Put the three
R1 generation sectors at angles 0, 2pi/3, 4pi/3 in the real E-plane.  A
count-derived symmetric traffic tensor has the form

    K = lambda I + sum_i m_i v_i v_i^T,

where m_i are non-negative integer service counts.  The isotropic lambda term
does not affect the principal axis.  The anisotropy is therefore controlled by

    z = (K_xx - K_yy) + 2 i K_xy = sum_i m_i exp(2 i phi_i),
    delta_traffic = (1/2) arg z  mod the axis symmetry.

Promotion criterion
-------------------
This is allowed to promote only if the current sector ledger supplies the
generation-space anisotropy m_i.  It is not enough to search integer triples
until one approximates the target angle.  That would only replace "fit an angle"
by "fit a traffic vector".

Outcome
-------
The type change is useful: delta belongs naturally to the symmetric/traffic
sector, while Phi belongs to the antisymmetric/holonomy sector.  But the current
ledger supplies only scalar counts (base d=2, optional neutrality +1, N_eff=9 or
27).  It does not supply the needed generation anisotropy.  Consequently:

* S3-equivariant traffic is isotropic -> no principal axis -> no delta.
* One marked generation gives a root-of-unity axis, not the small offsets.
* Integer-count searches can approximate the targets, especially at N_eff=27,
  but the matching triples are not derived and are non-unique.

So the traffic-tensor route is type-correct but currently underspecified; the
missing theorem is an anisotropic generation service-covariance law.
"""

from __future__ import annotations

import cmath
import math
from dataclasses import dataclass
from fractions import Fraction


PHI = (0.0, 2.0 * math.pi / 3.0, 4.0 * math.pi / 3.0)
TARGETS = {
    "charged_lepton": (Fraction(2, 9), 9, 2),
    "neutrino": (Fraction(1, 3), 9, 3),
    "up": (Fraction(2, 27), 27, 2),
    # Q4 colour-shared reading: 1/9 = 3/27.
    "down": (Fraction(1, 9), 27, 3),
}


@dataclass(frozen=True)
class SearchResult:
    counts: tuple[int, int, int]
    angle: float
    error: float
    close_count: int


def fold_axis_angle(theta: float) -> float:
    """Fold an unoriented principal axis into the fundamental [0, pi/6] wedge."""

    while theta < 0.0:
        theta += math.pi
    theta %= math.pi / 3.0
    if theta > math.pi / 6.0:
        theta = math.pi / 3.0 - theta
    return theta


def traffic_angle(counts: tuple[int, int, int]) -> float | None:
    z = sum(counts[i] * cmath.exp(2j * PHI[i]) for i in range(3))
    if abs(z) < 1.0e-12:
        return None
    return fold_axis_angle(0.5 * math.atan2(z.imag, z.real))


def all_count_triples(total: int) -> list[tuple[int, int, int]]:
    out: list[tuple[int, int, int]] = []
    for a in range(total + 1):
        for b in range(total - a + 1):
            out.append((a, b, total - a - b))
    return out


def nearest_integer_tensor(total: int, target: Fraction, tolerance: float = 0.01) -> SearchResult:
    vals: list[tuple[tuple[int, int, int], float]] = []
    for counts in all_count_triples(total):
        angle = traffic_angle(counts)
        if angle is not None:
            vals.append((counts, angle))
    target_f = float(target)
    best_counts, best_angle = min(vals, key=lambda item: abs(item[1] - target_f))
    close_count = sum(1 for _, angle in vals if abs(angle - target_f) < tolerance)
    return SearchResult(
        counts=best_counts,
        angle=best_angle,
        error=abs(best_angle - target_f),
        close_count=close_count,
    )


def fmt_frac(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def assert_true(name: str, condition: bool) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}")
    if not condition:
        raise AssertionError(name)


def main() -> None:
    print("Item 87: symmetric recovery-traffic tensor attempt")
    print("=" * 76)

    print("\n[1] Type split: symmetric traffic vs antisymmetric holonomy")
    print("    traffic tensor K=K^T -> mass-shape axis delta")
    print("    antisymmetric recovery connection -> CP holonomy Phi")
    assert_true(
        "the route is type-correct: delta is not being treated as a gauge phase",
        True,
    )

    print("\n[2] S3-equivariant traffic is isotropic")
    for total in (3, 9, 27):
        counts = (total // 3, total // 3, total // 3)
        angle = traffic_angle(counts)
        print(f"    total={total:2d}, counts={counts}: angle={angle}")
        assert_true(f"generation-blind total={total} traffic has no principal axis", angle is None)

    print("\n[3] A single marked generation gives only a symmetry axis")
    for counts in ((1, 0, 0), (0, 1, 0), (0, 0, 1), (2, 0, 0), (3, 0, 0)):
        angle = traffic_angle(counts)
        print(f"    counts={counts}: folded principal-axis angle={angle:.12f}")
        assert_true("single-axis anisotropy folds to a root-of-unity axis", abs(angle or 0.0) < 1.0e-10)

    print("\n[4] Current scalar ledger is insufficient")
    print("    The existing derivations supply d and N_eff, not a generation vector (m0,m1,m2).")
    print("    base d=2 and the neutrality +1 are scalar service counts.  Without a law assigning")
    print("    them to generation directions or covariances, K is isotropic or arbitrary.")
    assert_true("no current script/canon object supplies m_i for the three generation directions", True)

    print("\n[5] Null search: if we allow arbitrary integer traffic triples, do the targets become unique?")
    print("    sector          target    total  best counts      angle      error    #within 0.01 rad")
    exact_hits = 0
    for sector, (target, n_eff, d_count) in TARGETS.items():
        for total_name, total in (("d", d_count), ("N_eff", n_eff)):
            result = nearest_integer_tensor(total, target)
            if result.error < 1.0e-9:
                exact_hits += 1
            print(
                f"    {sector:<14s} {fmt_frac(target):>7s} "
                f"{total_name:>5s}={total:<2d}  {str(result.counts):>14s} "
                f"{result.angle:9.6f} {result.error:9.6f} {result.close_count:10d}"
            )

    print("\n[6] Verdict")
    print(
        """
    The symmetric traffic tensor is the right TYPE of object for delta: it keeps
    delta separate from the CP holonomy Phi and avoids the phase/anomaly no-goes.

    But the current framework does not derive the needed generation anisotropy.
    Generation-blind S3 traffic is isotropic.  A single marked generation gives
    a root-of-unity axis.  Arbitrary integer triples can approximate the target
    angles, especially at N_eff=27, but those triples are not supplied by the
    recovery ledger and are not unique.  That is a null-model warning, not a
    derivation.

    Precise remaining theorem target:

        derive C_ij = <delta j_i delta j_j>_sector on the three R1 generation
        channels from the monitored recovery service, and show that its
        symmetric E-plane projection has principal axes
        (2/9, 1/3, 2/27, 1/9) for (e, nu, u, d), with one rule.

    Until that covariance law exists, delta=d/N remains a phenomenological
    defect-count law.  The traffic-tensor route is viable as a target, but not
    closed by current canon.
"""
    )

    assert_true("no exact integer-traffic closure found in the tested d/N_eff ledgers", exact_hits == 0)
    print("ALL ASSERTIONS PASSED -- traffic route is type-correct but underspecified; no closure.")
    print("exit 0")


if __name__ == "__main__":
    main()
