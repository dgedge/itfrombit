#!/usr/bin/env python3
"""
Neutrino branch-selection verifier.

Question:
  Can the framework promote "branch-conditional NH" to an NH theorem by
  proving that the substrate selects the local delta_nu ~= 1/3 branch and
  excludes the remote IH-compatible branches?

Result:
  Conditional theorem verified; unconditional theorem not proved by current
  canon.

What is verified:
  1. The Koide/circulant phase has D3 generation-ring redundancy:
       delta -> delta + 2pi/3  (cyclic generation relabel)
       delta -> -delta         (reflection)
     so the physical phase can be folded to the principal domain [0, pi/3].

  2. In that principal domain, the observed oscillation ratio has:
       local NH branch near delta ~= 0.33059, close to 1/3;
       remote IH branch near delta ~= 1.041, close to pi/3.

  3. If the substrate defect-density rule is promoted to the principal
     branch selector
       delta_nu lies in the Voronoi cell of d_eff/N = 3/9
       i.e. |delta_nu - 1/3| <= 1/(2*9),
     then all ratio-compatible points in the selected cell are NH and no
     IH-compatible branch survives.

What is not proved:
  Current canon gives d_eff/N = 3/9 as a rational defect density and, after
  the L-W audit, only a transcendental phase near that value. It does not yet
  prove the phase-lift/principal-Voronoi selector. That selector is exactly
  the missing lemma needed to turn branch-conditional NH into a full NH theorem.

Run:
  python3 verify_neutrino_branch_selection.py
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass


R_NU = 1.0
DEFECT_CENTER = 1.0 / 3.0
N_PLAQUETTE = 9
VORONOI_HALF_WIDTH = 1.0 / (2.0 * N_PLAQUETTE)
DOMAIN_MAX = math.pi / 3.0

DM21_EXP = 7.41e-5
DM31_NO_EXP = 2.511e-3
RATIO_EXP = DM21_EXP / DM31_NO_EXP
RATIO_SIGMA = 0.00088


@dataclass(frozen=True)
class Hit:
    delta: float
    canonical: float
    perm: tuple[int, int, int]
    ratio: float
    ordering: str


def canonical_delta(delta: float) -> float:
    """Fold delta by D3 generation-ring symmetry to [0, pi/3]."""

    period = 2.0 * math.pi / 3.0
    x = delta % period
    if x > DOMAIN_MAX:
        x = period - x
    # Avoid -0.0 noise.
    return abs(x)


def mass_factors(delta: float, R: float = R_NU) -> list[float]:
    return [(1.0 + R * math.cos(delta + 2.0 * math.pi * n / 3.0)) ** 2 for n in range(3)]


def hits_at(delta: float, sigma_window: float = 1.0) -> list[Hit]:
    """Ratio-compatible assignments at a given delta."""

    factors = mass_factors(delta)
    out: list[Hit] = []
    for perm in itertools.permutations(range(3)):
        m1, m2, m3 = [factors[i] for i in perm]
        dm21 = m2 * m2 - m1 * m1
        dm31 = m3 * m3 - m1 * m1
        if dm21 <= 0.0 or abs(dm31) < 1e-15:
            continue
        ratio = dm21 / abs(dm31)
        if abs(ratio - RATIO_EXP) <= sigma_window * RATIO_SIGMA:
            out.append(
                Hit(
                    delta=delta,
                    canonical=canonical_delta(delta),
                    perm=perm,
                    ratio=ratio,
                    ordering="NH" if dm31 > 0.0 else "IH",
                )
            )
    return out


def scan_delta(
    lo: float,
    hi: float,
    steps: int,
    sigma_window: float = 1.0,
) -> list[Hit]:
    hits: list[Hit] = []
    for i in range(steps + 1):
        delta = lo + (hi - lo) * i / steps
        hits.extend(hits_at(delta, sigma_window=sigma_window))
    return hits


def summarize(label: str, hits: list[Hit]) -> None:
    print(label)
    if not hits:
        print("  no ratio-compatible hits")
        return
    nh = sum(h.ordering == "NH" for h in hits)
    ih = len(hits) - nh
    print("  total hits = %d, NH = %d, IH = %d" % (len(hits), nh, ih))
    print("  canonical delta range = %.12f .. %.12f" % (
        min(h.canonical for h in hits),
        max(h.canonical for h in hits),
    ))
    print("  permutations = %s" % sorted({h.perm for h in hits}))


def main() -> None:
    print("=" * 78)
    print("Neutrino branch-selection verifier")
    print("=" * 78)
    print("target ratio Delta m21^2 / |Delta m31^2| = %.8f +/- %.8f" % (
        RATIO_EXP,
        RATIO_SIGMA,
    ))
    print("D3 principal domain: [0, pi/3] = [0, %.12f]" % DOMAIN_MAX)
    print("defect-density center d_eff/N = 3/9 = %.12f" % DEFECT_CENTER)
    print("principal Voronoi cell: [%.12f, %.12f]" % (
        DEFECT_CENTER - VORONOI_HALF_WIDTH,
        DEFECT_CENTER + VORONOI_HALF_WIDTH,
    ))
    print()

    local_hits = scan_delta(
        DEFECT_CENTER - 0.08,
        DEFECT_CENTER + 0.08,
        steps=20000,
        sigma_window=1.0,
    )
    summarize("A. Local TCH neighborhood, delta=1/3 +/- 0.08", local_hits)
    print()

    voronoi_hits = scan_delta(
        DEFECT_CENTER - VORONOI_HALF_WIDTH,
        DEFECT_CENTER + VORONOI_HALF_WIDTH,
        steps=20000,
        sigma_window=1.0,
    )
    summarize("B. Principal defect-density Voronoi cell", voronoi_hits)
    assert all(h.ordering == "NH" for h in voronoi_hits)
    print("  assertion: all ratio-compatible points in selected cell are NH")
    print()

    global_hits = scan_delta(0.0, 2.0 * math.pi, steps=240000, sigma_window=1.0)
    summarize("C. Full phase circle before substrate branch selection", global_hits)
    ih_hits = [h for h in global_hits if h.ordering == "IH"]
    assert ih_hits, "expected remote IH branches in the unconstrained phase circle"
    nearest_ih = min(ih_hits, key=lambda h: abs(h.canonical - DEFECT_CENTER))
    print("  nearest IH canonical delta = %.12f" % nearest_ih.canonical)
    print("  distance from 3/9 center   = %.12f" % abs(nearest_ih.canonical - DEFECT_CENTER))
    print("  distance in Voronoi half-widths = %.3f" % (
        abs(nearest_ih.canonical - DEFECT_CENTER) / VORONOI_HALF_WIDTH
    ))
    print()

    print("=" * 78)
    print("VERDICT")
    print("=" * 78)
    print("  Conditional theorem verified:")
    print("    If the substrate selects the principal defect-density branch")
    print("    around d_eff/N = 3/9, then the remote IH branches are excluded")
    print("    and the neutrino ordering is NH.")
    print()
    print("  Not yet an unconditional theorem:")
    print("    The current framework proves/counts d_eff/N = 3/9 as the rational")
    print("    defect density, but after the L-W audit it does not yet prove the")
    print("    phase-lift rule that the actual transcendental Koide phase must lie")
    print("    in that Voronoi cell. That missing phase-lift/principal-branch lemma")
    print("    is the remaining load-bearing step.")


if __name__ == "__main__":
    main()
