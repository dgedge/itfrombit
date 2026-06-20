#!/usr/bin/env python3
r"""BLACK-HOLE HORIZON SURFACE SCRAMBLING GAP.

This script does not claim to derive the full localized-mass QEC steady state.
It closes a narrower fork:

    Can a purely local horizon-surface service graph produce black-hole fast
    scrambling, or must the QEC service layer be expander/active-address-like?

Result:

  * A nearest-neighbour 2D horizon surface has normalized Markov gap

        gap_local(L) = (1 - cos(2 pi / L)) / 2 ~ pi^2 / N

    for N=L^2 cells.  Its mixing time is O(N log N), not O(log N).  Purely
    local surface diffusion is therefore too slow by O(N).

  * A parallel active-address/expander service graph has O(1) gap and hence
    O(log N) mixing in service ticks, the required fast-scrambling scaling.

  * A serial global demux that services only one cell per whole horizon tick
    also fails: its gap is the expander gap divided by N.

The remaining theorem is thus sharply stated: the localized horizon QEC
scheduler must realize an O(1)-gap service graph per cell-tick, not a local
torus walk and not a globally serial address queue.
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
    return 2.0 * m_gev / M_PLANCK_GEV**2


def horizon_cells(m_solar: float) -> float:
    m_gev = m_solar * M_SUN_GEV
    rs = horizon_radius(m_gev)
    return 16.0 * math.pi * rs * rs * LAMBDA_QCD_GEV**2


def local_torus_gap_from_n(n: float) -> float:
    """Normalized random-walk gap on an LxL torus, analytically continued."""
    L = math.sqrt(n)
    # Stable form of (1 - cos(2*pi/L))/2 for astronomical L.
    return math.sin(math.pi / L) ** 2


def ramanujan_gap(degree: int) -> float:
    """Best-case normalized adjacency gap scale for a d-regular expander."""
    return 1.0 - 2.0 * math.sqrt(degree - 1.0) / degree


def mixing_time_bound(n: float, gap: float, eps: float = 1.0e-2) -> float:
    """Standard reversible-chain scale t ~ gap^{-1} log(N/eps)."""
    return math.log(n / eps) / gap


def main() -> None:
    print("[1] Markov gaps for candidate horizon service graphs")
    exp_gap_8 = ramanujan_gap(8)
    complete_gap = 1.0
    print(f"    8-regular expander O(1) gap >= {exp_gap_8:.6f}")
    print(f"    complete active-address gap    = {complete_gap:.6f}")
    check(0.3 < exp_gap_8 < 0.4, "8-regular expander has an O(1) normalized gap")

    print("\n[2] Scaling table")
    header = (
        "    N_cells        local gap      local mix       expander mix"
        "    serial-expander mix"
    )
    print(header)
    for n in (16**2, 32**2, 64**2, 128**2, 1.0e8):
        local_gap = local_torus_gap_from_n(n)
        local_mix = mixing_time_bound(n, local_gap)
        exp_mix = mixing_time_bound(n, exp_gap_8)
        serial_mix = mixing_time_bound(n, exp_gap_8 / n)
        print(
            f"    {n:11.3e}  {local_gap:12.3e}  {local_mix:12.3e}"
            f"  {exp_mix:12.3e}  {serial_mix:12.3e}"
        )
    check(
        mixing_time_bound(1.0e8, local_torus_gap_from_n(1.0e8))
        > 1.0e6 * mixing_time_bound(1.0e8, exp_gap_8),
        "local 2D surface mixing is parametrically slower than expander service",
    )

    print("\n[3] Astrophysical horizon sizes")
    for m_solar in (3.0, 30.0, 4.3e6):
        n = horizon_cells(m_solar)
        local_gap = local_torus_gap_from_n(n)
        local_mix = mixing_time_bound(n, local_gap)
        exp_mix = mixing_time_bound(n, exp_gap_8)
        serial_mix = mixing_time_bound(n, exp_gap_8 / n)
        print(f"    M={m_solar:g} Msun: N_H={n:.3e}, log(N_H)={math.log(n):.2f}")
        print(f"      local surface mix       ~ {local_mix:.3e} cell-service ticks")
        print(f"      parallel expander mix   ~ {exp_mix:.3e} cell-service ticks")
        print(f"      serial global demux mix ~ {serial_mix:.3e} horizon ticks")
        check(local_mix / exp_mix > 1.0e35, "astrophysical local diffusion is not fast scrambling")
        check(serial_mix / exp_mix > 1.0e35, "serial demux is also not fast scrambling")

    print(
        """
[4] VERDICT
    The black-hole QEC scrambling target is no longer vague.  A horizon
    scheduler built only from local nearest-neighbour surface updates gives
    t_mix ~ N_H log N_H and fails the fast-scrambling scaling by O(N_H).
    A globally serial demux fails for the same reason.

    The only finite-service route left by this audit is an O(1)-gap service
    graph: active-address / expander-like mixing in parallel cell-service time.
    Therefore the remaining theorem is exact:

        derive an O(1)-gap horizon service-current graph from the QEC scheduler,
        or concede that the present local surface dynamics is not a fast
        scrambler.

exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
