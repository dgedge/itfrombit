#!/usr/bin/env python3
r"""Topological obstruction to horizon fast scrambling.

Question
--------
Can the existing local horizon-QEC service graph be "slightly enriched" into a
fast scrambler, or does fast scrambling require a genuinely nonlocal service
primitive?

Answer
------
A bounded-degree graph embedded on a genus-g horizon obeys a separator bound.
An O(1)-gap expander, by contrast, has linear bisection width.  Combining the
two gives

    g + 1 >= gamma^2 N / 128,

up to the standard separator constant, for a normalized Markov gap gamma.  Thus
an O(1)-gap bounded-degree horizon-cell graph requires extensive genus or
non-spatial address links.  A genus-0 Schwarzschild horizon cannot supply it.

This script does not prove a new expander.  It hardens the negative/conditional
BH fast-scrambling status:

  * any bounded local Moore/all-contact stencil still fails;
  * any finite-range geometric stencil fails unless the range grows like the
    horizon diameter, at which point the degree grows with N;
  * a true O(1)-gap service graph is necessarily a new nonlocal address-space
    primitive, not an unexploited consequence of V_cell/V_Sch or local KMS.
"""

from __future__ import annotations

import math


LAMBDA_QCD_GEV = 0.332
M_PLANCK_GEV = 1.220890e19
GEV_PER_KG = 5.60959e26
M_SUN_GEV = 1.98892e30 * GEV_PER_KG


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def horizon_radius(m_gev: float) -> float:
    # Natural units: G = 1/M_P^2, R_s = 2GM.
    return 2.0 * m_gev / M_PLANCK_GEV**2


def horizon_cells(m_solar: float) -> float:
    m_gev = m_solar * M_SUN_GEV
    rs = horizon_radius(m_gev)
    # Canonical cell count used by the existing BH scrambling audits.
    return 16.0 * math.pi * rs * rs * LAMBDA_QCD_GEV**2


def expander_cut_lower_bound(n: float, degree: int, gamma: float) -> float:
    r"""Minimum bisection edge cut from normalized gap gamma.

    Cheeger: Phi >= gamma/2.  For |S| = N/2 in a d-regular graph, the boundary
    edge count is at least Phi*d*|S| = gamma*d*N/4.
    """

    return gamma * degree * n / 4.0


def separator_cut_upper_bound(n: float, degree: int, genus: float) -> float:
    r"""Genus-g separator converted to a bounded-degree edge cut.

    Gilbert-Hutchinson-Tarjan form: O(sqrt((g+1)N)); we keep the common
    2*sqrt(2(g+1)N) constant because only scaling and order are load-bearing.
    """

    return degree * 2.0 * math.sqrt(2.0 * (genus + 1.0) * n)


def minimum_genus_for_gap(n: float, gamma: float) -> float:
    # From gamma*d*N/4 <= d*2*sqrt(2(g+1)N).
    return max(0.0, gamma * gamma * n / 128.0 - 1.0)


def torus_range_gap_bound(n: float, radius: int) -> float:
    r"""Small-k upper bound for a finite-range geometric stencil on an LxL torus.

    For any bounded stencil with Chebyshev radius R, the k=(1,0) Fourier mode
    changes by O((R/L)^2).  The constant is deliberately generous; if this
    upper bound still tends to zero, no finite local stencil can be an expander.
    """

    l_side = math.sqrt(n)
    return min(1.0, (math.pi * radius / l_side) ** 2)


def degree_for_chebyshev_radius(radius: int) -> int:
    return (2 * radius + 1) ** 2 - 1


def main() -> None:
    print("BH FAST-SCRAMBLING TOPOLOGICAL OBSTRUCTION")
    print("=" * 86)

    gamma = 0.30
    degree = 8
    print(f"\n[1] Expander bisection versus horizon separator (gamma={gamma}, d={degree})")
    for n in (10**4, 10**6, 10**8):
        cut_low = expander_cut_lower_bound(n, degree, gamma)
        cut_sphere = separator_cut_upper_bound(n, degree, genus=0)
        g_min = minimum_genus_for_gap(n, gamma)
        print(
            f"    N={n:>9.3e}: expander cut >= {cut_low:>10.3e}; "
            f"genus-0 separator <= {cut_sphere:>10.3e}; "
            f"required genus >= {g_min:>10.3e}"
        )
    check(
        expander_cut_lower_bound(10**8, degree, gamma)
        > 1.0e2 * separator_cut_upper_bound(10**8, degree, genus=0),
        "genus-0 horizon cannot carry an O(1)-gap bounded-degree expander at large N",
    )

    print("\n[2] Astrophysical horizon scale: extensive nonlocality required")
    for m_solar in (3.0, 30.0, 4.3e6):
        n = horizon_cells(m_solar)
        g_min = minimum_genus_for_gap(n, gamma)
        cut_low = expander_cut_lower_bound(n, degree, gamma)
        print(
            f"    M={m_solar:g} Msun: N_H={n:.3e}; "
            f"required genus >= {g_min:.3e}; "
            f"linear bisection contacts >= {cut_low:.3e}"
        )
        check(g_min / n > 6.0e-4, "required genus is extensive in horizon cell count")

    print("\n[3] Bounded geometric stencil check")
    n = 1.0e8
    for radius in (1, 2, 3, 5, 10):
        deg = degree_for_chebyshev_radius(radius)
        gap_bound = torus_range_gap_bound(n, radius)
        print(f"    R={radius:2d}: degree={deg:4d}, gap upper-bound ~ {gap_bound:.3e}")
    check(
        torus_range_gap_bound(n, radius=10) < 1.0e-5,
        "finite-range all-contact/Moore-like stencils still have collapsing gap",
    )

    print("\n[4] What would be needed to get gamma=0.3 geometrically?")
    l_side = math.sqrt(n)
    needed_radius = math.sqrt(gamma) * l_side / math.pi
    needed_degree = degree_for_chebyshev_radius(math.ceil(needed_radius))
    print(
        f"    On N={n:.1e} cells, a geometric stencil would need R~{needed_radius:.1f}, "
        f"degree~{needed_degree:.3e}."
    )
    check(
        needed_degree / n > 0.03,
        "geometric fast scrambling needs horizon-scale fanout, not a finite service alphabet",
    )

    print(
        r"""
[5] VERDICT
    The fast-scrambling frontier is now sharper.

    Current canon gives a genus-0, bounded-degree, local horizon service graph:
    direct-sum V_Sch shell labels plus one-bit/local KMS service.  By the
    separator obstruction, no such graph can have an O(1) spectral gap.  Enlarging
    the local Moore/all-contact stencil changes only the constant, not the
    scaling.

    Therefore fast scrambling is not merely "not yet derived" from the local
    channel.  It requires one of two new structures:

      A. an extensive nonlocal contact topology on the horizon (genus/long links
         scaling with N_H), or
      B. a non-spatial address-space expander service layer whose links are not
         geometric horizon adjacencies.

    Neither is currently in canon.  Without A or B, the framework predicts
    no substrate-native fast scrambling.  That is a falsifiable black-hole
    signature, not a bookkeeping inconvenience.

ALL CHECKS PASSED"""
    )


if __name__ == "__main__":
    main()
