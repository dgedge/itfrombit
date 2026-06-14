#!/usr/bin/env python3
"""
Audit the TCH Part 05 neutrino mass-ordering claim.

Question:
  Does the Part 05 neutrino Koide/circulant block structurally force normal
  ordering, and does that survive the DRIFT M9 / Lindemann-Weierstrass
  recharacterisation of delta_nu = 1/3 from exact rational phase to
  transcendental-near-1/3?

This script deliberately separates four statements:

  A. Dimensionless spectrum from the documented formula.
  B. Permutation audit: which assignment of raw nodes to (nu1,nu2,nu3)
     satisfies the solar convention Delta m21^2 > 0 and the observed small
     ratio Delta m21^2 / |Delta m31^2|?
  C. Local robustness around the TCH branch delta_nu ~= 1/3.
  D. Absolute scale: what is predicted only after importing/anchoring
     |Delta m31^2|, and what remains open (mu_nu dimensionalisation).

No external packages are required.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass


# Part 05 / NuFIT 5.3 central values used in the repo text.
DM21_EXP = 7.41e-5
DM31_NO_EXP = 2.511e-3
RATIO_EXP = DM21_EXP / DM31_NO_EXP
RATIO_SIGMA = 0.00088  # quoted in Part 05 for Delta m21^2 / Delta m31^2

R_TCH = 1.0
DELTA_TCH = 1.0 / 3.0


@dataclass(frozen=True)
class Assignment:
    """Raw-node assignment to physical labels (nu1, nu2, nu3)."""

    perm: tuple[int, int, int]
    dm21: float
    dm31: float
    ratio: float

    @property
    def ordering(self) -> str:
        return "NH" if self.dm31 > 0 else "IH"

    @property
    def rel_ratio_error(self) -> float:
        return abs(self.ratio - RATIO_EXP) / RATIO_EXP


def raw_amplitudes(R: float, delta: float) -> list[float]:
    """Koide amplitude a_n = 1 + R cos(delta + 2 pi n/3)."""

    return [1.0 + R * math.cos(delta + 2.0 * math.pi * n / 3.0) for n in range(3)]


def mass_factors(R: float, delta: float) -> list[float]:
    """Dimensionless mass factors m_n / mu = a_n^2."""

    return [a * a for a in raw_amplitudes(R, delta)]


def assignments(R: float, delta: float) -> list[Assignment]:
    """
    Enumerate physical label assignments.

    The neutrino convention fixes Delta m21^2 > 0.  The sign of Delta m31^2
    then distinguishes normal (positive) from inverted (negative) ordering.
    """

    mf = mass_factors(R, delta)
    out: list[Assignment] = []
    for perm in itertools.permutations(range(3)):
        m1, m2, m3 = [mf[i] for i in perm]
        dm21 = m2 * m2 - m1 * m1
        dm31 = m3 * m3 - m1 * m1
        if dm21 <= 0.0 or abs(dm31) < 1e-15:
            continue
        out.append(Assignment(perm, dm21, dm31, dm21 / abs(dm31)))
    return out


def fit_delta_bisection(
    R: float,
    lo: float = 0.25,
    hi: float = 0.40,
    target: float = RATIO_EXP,
    steps: int = 80,
) -> float:
    """
    Fit the local NH TCH branch for the ratio.

    On the branch near delta = 1/3, the matching permutation is (1,2,0):
      raw node 1 -> nu1, raw node 2 -> nu2, raw node 0 -> nu3.
    """

    perm = (1, 2, 0)

    def branch_ratio(delta: float) -> float:
        mf = mass_factors(R, delta)
        m1, m2, m3 = [mf[i] for i in perm]
        return (m2 * m2 - m1 * m1) / abs(m3 * m3 - m1 * m1)

    flo = branch_ratio(lo) - target
    fhi = branch_ratio(hi) - target
    if flo * fhi > 0:
        raise RuntimeError("bisection bracket does not straddle the target")

    for _ in range(steps):
        mid = 0.5 * (lo + hi)
        fmid = branch_ratio(mid) - target
        if flo * fmid <= 0:
            hi = mid
            fhi = fmid
        else:
            lo = mid
            flo = fmid
    return 0.5 * (lo + hi)


def scaled_masses(R: float, delta: float, dm31_abs: float = DM31_NO_EXP) -> tuple[float, list[float]]:
    """Fix mu from |Delta m31^2| on the local NH assignment and return masses."""

    perm = (1, 2, 0)
    mf = mass_factors(R, delta)
    ordered = [mf[i] for i in perm]
    denom = ordered[2] * ordered[2] - ordered[0] * ordered[0]
    mu = math.sqrt(dm31_abs / denom)
    return mu, [mu * x for x in ordered]


def local_scan(
    R: float,
    center: float = DELTA_TCH,
    half_width: float = 0.08,
    grid: int = 16001,
    sigma_window: float = 1.0,
) -> list[tuple[float, Assignment]]:
    """Find assignments within a ratio sigma-window near the TCH phase branch."""

    lo = center - half_width
    hi = center + half_width
    hits: list[tuple[float, Assignment]] = []
    for j in range(grid):
        delta = lo + (hi - lo) * j / (grid - 1)
        for a in assignments(R, delta):
            if abs(a.ratio - RATIO_EXP) <= sigma_window * RATIO_SIGMA:
                hits.append((delta, a))
    return hits


def global_scan(
    R: float,
    grid: int = 120001,
    rel_tol: float = 0.03,
) -> list[tuple[float, Assignment]]:
    """
    Coarse global scan over [0, 2pi).

    This is not the TCH claim; it checks whether the raw circulant has remote
    phase branches that can mimic the ratio with other orderings.
    """

    hits: list[tuple[float, Assignment]] = []
    for j in range(grid):
        delta = 2.0 * math.pi * j / grid
        for a in assignments(R, delta):
            if a.rel_ratio_error <= rel_tol:
                hits.append((delta, a))
    return hits


def local_r_delta_scan(
    r_lo: float = 0.90,
    r_hi: float = 1.10,
    delta_center: float = DELTA_TCH,
    delta_half_width: float = 0.08,
    r_grid: int = 161,
    delta_grid: int = 321,
    sigma_window: float = 1.0,
) -> list[tuple[float, float, Assignment]]:
    """Small 2D robustness scan around the TCH point (R, delta) = (1, 1/3)."""

    hits: list[tuple[float, float, Assignment]] = []
    d_lo = delta_center - delta_half_width
    d_hi = delta_center + delta_half_width
    for i in range(r_grid):
        R = r_lo + (r_hi - r_lo) * i / (r_grid - 1)
        for j in range(delta_grid):
            delta = d_lo + (d_hi - d_lo) * j / (delta_grid - 1)
            for a in assignments(R, delta):
                if abs(a.ratio - RATIO_EXP) <= sigma_window * RATIO_SIGMA:
                    hits.append((R, delta, a))
    return hits


def print_assignment_table(R: float, delta: float) -> None:
    print("B. Permutation audit at R=1, delta=1/3")
    print("   perm gives raw-node labels (nu1, nu2, nu3)")
    print("   perm        ordering   ratio       rel.error")
    for a in sorted(assignments(R, delta), key=lambda x: (x.rel_ratio_error, x.perm)):
        print(
            "   %-11s %-8s %.8f   %.2f%%"
            % (str(a.perm), a.ordering, a.ratio, 100.0 * a.rel_ratio_error)
        )


def main() -> None:
    print("=" * 78)
    print("TCH neutrino ordering audit")
    print("=" * 78)
    print("Inputs from repo Part 05:")
    print("  R_nu = 1")
    print("  delta_nu branch = 1/3 (read after DRIFT M9 as transcendental-near-1/3)")
    print("  target ratio Delta m21^2 / |Delta m31^2| = %.8f +/- %.8f" % (RATIO_EXP, RATIO_SIGMA))
    print()

    amps = raw_amplitudes(R_TCH, DELTA_TCH)
    mf = mass_factors(R_TCH, DELTA_TCH)
    print("A. Raw Part 05 spectrum")
    print("   n   amplitude a_n      mass factor a_n^2")
    for n, (a, m) in enumerate(zip(amps, mf)):
        print("   %d   %.12f      %.12f" % (n, a, m))
    print("   raw order: n=1 light, n=2 middle, n=0 heavy")
    print()

    print_assignment_table(R_TCH, DELTA_TCH)
    print()

    delta_fit = fit_delta_bisection(R_TCH)
    a_fit = [a for a in assignments(R_TCH, delta_fit) if a.perm == (1, 2, 0)][0]
    print("C. Phase-retraction robustness")
    print("   exact ratio fit on the local TCH branch:")
    print("     delta_fit = %.12f" % delta_fit)
    print("     delta_fit - 1/3 = %.12e" % (delta_fit - DELTA_TCH))
    print("     ordering = %s, perm = %s, ratio = %.8f" % (a_fit.ordering, a_fit.perm, a_fit.ratio))
    hits_1s = local_scan(R_TCH, sigma_window=1.0)
    nh_1s = sum(1 for _, a in hits_1s if a.ordering == "NH")
    ih_1s = len(hits_1s) - nh_1s
    if hits_1s:
        lo = min(d for d, _ in hits_1s)
        hi = max(d for d, _ in hits_1s)
        perms = sorted({a.perm for _, a in hits_1s})
        print("   local +/-0.08 rad, 1-sigma ratio hits:")
        print("     delta range = %.12f .. %.12f" % (lo, hi))
        print("     NH hits = %d, IH hits = %d, perms = %s" % (nh_1s, ih_1s, perms))
    print("   interpretation: exact rational 1/3 is not needed; the whole local")
    print("   ratio-compatible TCH branch is normal ordering.")
    hits_2d = local_r_delta_scan()
    nh_2d = sum(1 for _, _, a in hits_2d if a.ordering == "NH")
    ih_2d = len(hits_2d) - nh_2d
    if hits_2d:
        rlo = min(R for R, _, _ in hits_2d)
        rhi = max(R for R, _, _ in hits_2d)
        dlo = min(delta for _, delta, _ in hits_2d)
        dhi = max(delta for _, delta, _ in hits_2d)
        perms = sorted({a.perm for _, _, a in hits_2d})
        print("   local two-parameter guardrail, R=0.9..1.1 and delta=1/3+/-0.08:")
        print("     R range among 1-sigma hits = %.6f .. %.6f" % (rlo, rhi))
        print("     delta range among 1-sigma hits = %.12f .. %.12f" % (dlo, dhi))
        print("     NH hits = %d, IH hits = %d, perms = %s" % (nh_2d, ih_2d, perms))
    print()

    hits_global = global_scan(R_TCH)
    nh_global = sum(1 for _, a in hits_global if a.ordering == "NH")
    ih_global = len(hits_global) - nh_global
    print("D. Global branch caveat")
    print("   scanning delta in [0,2pi) at 3% relative-ratio tolerance:")
    print("     NH hits = %d, IH hits = %d" % (nh_global, ih_global))
    print("   interpretation: NH is not a theorem of the abstract R=1 circulant.")
    print("   It is a theorem only after selecting the TCH local phase branch near 1/3.")
    print()

    mu, masses = scaled_masses(R_TCH, DELTA_TCH)
    mu_fit, masses_fit = scaled_masses(R_TCH, delta_fit)
    print("E. Absolute scale audit")
    print("   using delta=1/3 and fixing mu from |Delta m31^2|:")
    print("     mu = %.12e eV" % mu)
    print("     masses = [%.6f, %.6f, %.6f] meV" % tuple(1000.0 * m for m in masses))
    print("     sum = %.6f meV" % (1000.0 * sum(masses)))
    print("   using delta_fit and fixing mu from |Delta m31^2|:")
    print("     mu = %.12e eV" % mu_fit)
    print("     masses = [%.6f, %.6f, %.6f] meV" % tuple(1000.0 * m for m in masses_fit))
    print("     sum = %.6f meV" % (1000.0 * sum(masses_fit)))
    print("   interpretation: the dimensionless ordering/ratio is tested here; the")
    print("   absolute mass scale is still imported through |Delta m31^2| unless")
    print("   item-53 / seesaw-dimensionalisation supplies mu_nu.")
    print()

    print("=" * 78)
    print("VERDICT")
    print("=" * 78)
    print("  Progress: the Part 05 local TCH branch does make a clean NH prediction.")
    print("  The sign survives DRIFT M9: replacing exact delta_nu=1/3 by a nearby")
    print("  transcendental phase that fits the ratio keeps the same permutation")
    print("  (raw n=1,n=2,n=0) -> (nu1,nu2,nu3), hence Delta m31^2 > 0.")
    print()
    print("  Limit: this is not yet a full ordering theorem from bit-budget alone.")
    print("  It still assumes the Part 05 neutrino mass block, the TCH local phase")
    print("  branch near 1/3, and the mass-eigenstate labeling/reconciliation with")
    print("  PMNS eigenvectors. The absolute scale remains open.")


if __name__ == "__main__":
    main()
