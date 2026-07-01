#!/usr/bin/env python3
r"""ITEM 123 / 132: zero-mode mass as R4 polarisation, not MOND double-counting.

Question
--------
The CMB branch wants the R4 zero-mode reservoir to be pressureless mobile
dark matter.  The galaxy branch wants an active R4/MOND response.  Adding both
as independent halo components double-counts.  A more promising reading is:

    N_zero   = the conserved rest-mass reservoir,
    N_active = a baryon-coupled response / polarisation of that same reservoir.

Then active R4 is not an extra mass component.  It is a constitutive relation
for how the zero-mode reservoir is arranged around baryons.

This script tests the accounting and the simplest consequence.  If the R4
polarisation profile is the cored MOND/phantom profile already derived in
item132_halo_closure.py,

    rho_z(r) = A / (r^2 + r_c^2),
    A = sqrt(G M_b a0) / (4 pi G),
    r_c = r_M/3,      r_M = sqrt(G M_b/a0),

then the finite zero-mode budget fixes a truncation radius rather than a
depletion factor:

    M_z(R_t) / M_b = (Omega_zero/Omega_b) / f_ret,

where f_ret is the fraction of the fair cosmic baryon share retained in the
visible baryonic galaxy.  The f_ret=1 row is the compact "all baryons retained"
limit; smaller f_ret gives a larger halo edge, as in ordinary abundance
matching.

Result
------
The polarisation branch is algebraically consistent and removes the >95%
depletion problem, but it is not a full galaxy-formation theorem.  It replaces
the ordinary NFW/fair-sample halo by an R4-polarised zero-mode halo.  The local
shape theorem is supplied by the companion scripts:

    geodesic zero-mode dust + baryonic R4 susceptibility
      -> rho_z(r) = A/(r^2+r_c^2) up to the finite reservoir edge.

item123_zero_mode_r4_constitutive_law.py reduces the cored profile to the
minimal shell-capture latch, and
item123_zero_mode_shell_capture_latch_theorem.py derives that latch from P1.
Thus MOND/RAR is not an extra force added to CDM if P1 is accepted; it is the
effective field of the same zero-mode mass distribution.  If P1 fails, the
adopted conservative branch remains ordinary zero-mode CDM and active R4/MOND
must stay subleading.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from item123_cmb_cold_matter_closure_status import OMEGA_B_H2, density_for_power
from item132_halo_closure import A0, G, MSUN, amplitude, core_radius, halo_mass, r_mond


KPC = 3.0856775814913673e19


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def exchange_energy(n_zero: float, n_active: float, epsilon: float = 1.0) -> float:
    """Rest-count Hamiltonian for a closed zero/active reservoir."""

    return epsilon * (n_zero + n_active)


def active_fraction(x: float, n_tot: float) -> float:
    """Same-cell exchange fraction from x/(N_tot+x)."""

    return x / (n_tot + x)


def mass_ratio_cored(y: float) -> float:
    """M_halo(y r_M) / M_b for r_c=r_M/3 and A=sqrt(G M_b a0)/(4 pi G)."""

    return y - (1.0 / 3.0) * math.atan(3.0 * y)


def solve_y_for_ratio(target: float, lo: float = 1.0e-9, hi: float = 100.0) -> float:
    while mass_ratio_cored(hi) < target:
        hi *= 2.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if mass_ratio_cored(mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def halo_ratio_direct(m_b: float, y: float) -> float:
    r_m = r_mond(m_b)
    amp = amplitude(m_b)
    rc = core_radius(m_b)
    return halo_mass(np.asarray([y * r_m]), amp, rc)[0] / m_b


def deep_mond_edge_quality(y: float) -> float:
    """Return g_h^2/(a0 g_b) for the cored halo at radius y r_M.

    The ideal asymptotic MOND value is 1.  The cored finite-reservoir edge is
    not required to be exactly asymptotic; this reports the deviation.
    """

    ratio = mass_ratio_cored(y)
    # g_h/g_b = M_h/M_b for a spherical halo at the same radius.
    # g_b/a0 = 1/y^2, so g_h^2/(a0 g_b) = (ratio/y)^2.
    return (ratio / y) ** 2


@dataclass(frozen=True)
class GalaxyRow:
    m_b_msun: float
    f_ret: float
    r_m_kpc: float
    r_c_kpc: float
    r_t_zero_kpc: float
    r_t_dark_kpc: float
    y_zero: float
    y_dark: float


def galaxy_row(m_b_msun: float, zero_to_baryon: float, dark_to_baryon: float, f_ret: float) -> GalaxyRow:
    m_b = m_b_msun * MSUN
    r_m = r_mond(m_b)
    y_zero = solve_y_for_ratio(zero_to_baryon / f_ret)
    y_dark = solve_y_for_ratio(dark_to_baryon / f_ret)
    return GalaxyRow(
        m_b_msun=m_b_msun,
        f_ret=f_ret,
        r_m_kpc=r_m / KPC,
        r_c_kpc=core_radius(m_b) / KPC,
        r_t_zero_kpc=y_zero * r_m / KPC,
        r_t_dark_kpc=y_dark * r_m / KPC,
        y_zero=y_zero,
        y_dark=y_dark,
    )


def main() -> None:
    print("ITEM 123 / 132: ZERO-MODE AS R4 POLARISATION BRANCH")
    print("=" * 100)

    omega_nur, omega_zero, omega_dark, zeq = density_for_power(1)
    zero_to_baryon = omega_zero / OMEGA_B_H2
    dark_to_baryon = omega_dark / OMEGA_B_H2
    y_zero = solve_y_for_ratio(zero_to_baryon)
    y_dark = solve_y_for_ratio(dark_to_baryon)

    print("\n[1] Rest-count exchange: active response is not extra rest mass")
    e0 = exchange_energy(1000.0, 0.0)
    e1 = exchange_energy(990.0, 10.0)
    e2 = exchange_energy(1000.0 - math.pi, math.pi)
    print(f"  E(N0=1000, Na=0)        = {e0:.12f}")
    print(f"  E(N0=990,  Na=10)       = {e1:.12f}")
    print(f"  E(N0=1000-pi, Na=pi)    = {e2:.12f}")
    check(abs(e0 - e1) < 1.0e-12, "zero->active exchange preserves rest-count energy")
    check(abs(e0 - e2) < 1.0e-12, "continuous response variable carries no extra T00 in the rest ledger")
    print("  Therefore N_active can be a response coordinate only if the stress-energy")
    print("  remains T00 = epsilon (N_zero+N_active), not epsilon N_zero plus a")
    print("  second independent MOND mass ledger.")

    print("\n[2] Finite reservoir edge from the cored R4 polarisation profile")
    print(f"  omega_b h^2             = {OMEGA_B_H2:.6f}")
    print(f"  omega_nuR h^2           = {omega_nur:.6f}")
    print(f"  omega_zero h^2          = {omega_zero:.6f}")
    print(f"  omega_dark h^2          = {omega_dark:.6f}")
    print(f"  z_eq                    = {zeq:.1f}")
    print(f"  zero-mode/baryon ratio  = {zero_to_baryon:.6f}")
    print(f"  total-dark/baryon ratio = {dark_to_baryon:.6f}")
    print("  cored R4 profile mass ratio:")
    print("      M_z(y r_M)/M_b = y - (1/3) atan(3y)")
    print(f"  y_t for zero-mode only  = {y_zero:.6f} r_M")
    print(f"  y_t for all dark matter = {y_dark:.6f} r_M")
    check(abs(mass_ratio_cored(y_zero) / zero_to_baryon - 1.0) < 1.0e-12, "zero-mode budget fixes a unique truncation radius")
    check(abs(mass_ratio_cored(y_dark) / dark_to_baryon - 1.0) < 1.0e-12, "total dark budget likewise fixes a unique outer edge")
    check(4.0 < y_zero < 6.0, "zero-mode edge is a few MOND radii, not a 96% depletion")
    check(5.0 < y_dark < 7.0, "including nu_R moves the edge outward by an O(1) amount")

    print("\n[3] Galaxy-scale numbers and the baryon-retention dependence")
    print("  The f_ret=1 row is the all-cosmic-baryons-retained lower edge.  Real")
    print("  galaxies generally retain less than their fair baryon share, so the")
    print("  zero-mode reservoir per visible baryon is larger by 1/f_ret.")
    print("  M_b [Msun]      f_ret  r_M[kpc]  y_t,z  R_t,z[kpc]  y_t,d  R_t,dark[kpc]")
    for m_b_msun in (1.0e9, 1.0e10, 6.0e10, 1.0e11):
        for f_ret in (1.0, 0.5, 0.25, 0.10):
            row = galaxy_row(m_b_msun, zero_to_baryon, dark_to_baryon, f_ret)
            print(
                f"  {row.m_b_msun:10.3e}  {row.f_ret:5.2f} "
                f"{row.r_m_kpc:8.3f} {row.y_zero:6.2f} {row.r_t_zero_kpc:11.3f}"
                f" {row.y_dark:6.2f} {row.r_t_dark_kpc:14.3f}"
            )
    mw_low_edge = galaxy_row(6.0e10, zero_to_baryon, dark_to_baryon, 1.0)
    mw_retained = galaxy_row(6.0e10, zero_to_baryon, dark_to_baryon, 0.25)
    check(35.0 < mw_low_edge.r_t_zero_kpc < 60.0, "all-baryons-retained Milky-Way-like row is a compact lower-edge limit")
    check(120.0 < mw_retained.r_t_zero_kpc < 220.0, "25% baryon retention gives a virial-scale Milky-Way-like edge")

    print("\n[4] Does the edge still look MOND-like?")
    for y in (1.0, 3.0, y_zero, y_dark, 10.0, 30.0):
        print(
            f"  y={y:7.3f}: M_h/M_b={mass_ratio_cored(y):8.4f}, "
            f"g_h^2/(a0 g_b)={deep_mond_edge_quality(y):7.4f}, "
            f"x=|g|/a0~1/y={1.0/y:7.4f}"
        )
    check(deep_mond_edge_quality(y_zero) > 0.75, "zero-mode edge remains within O(20%) of the asymptotic MOND relation")
    check(deep_mond_edge_quality(30.0) > 0.96, "far outside the core the cored profile tends to the MOND asymptote")

    print("\n[5] Active fraction is a response, not depletion")
    n_tot = 1.0e6
    for y in (1.0, y_zero, y_dark):
        x = 1.0 / y
        f = active_fraction(x, n_tot)
        print(f"  y={y:7.3f}: x=|g|/a0={x:.6f}, active fraction x/(Ntot+x)={f:.3e}")
    check(active_fraction(1.0 / y_zero, n_tot) < 1.0e-6, "polarisation needs only a tiny active fraction of a large reservoir")
    print("  This is why the polarisation branch avoids the depletion contradiction:")
    print("  the zero-mode mass remains present, but the active coordinate is not")
    print("  another mass budget.  It is the response/readout arranging the reservoir.")

    print("\n[6] Branch verdict")
    print("  WORKS AS ACCOUNTING:")
    print("    N_active may be a baryon-coupled response/readout of N_zero without")
    print("    adding a second rest-mass ledger.  The finite zero-mode abundance then")
    print("    fixes a halo edge R_t ~= %.3f r_M for the cored R4 profile." % y_zero)
    print("  LOCAL SHAPE THEOREM STATUS:")
    print("    the companion constitutive-law and shell-capture-latch scripts derive")
    print("    rho=A/(r^2+r_c^2) from R4 p=3 plus the P1 two-state service latch.")
    print("    The polarised branch is therefore closed at local P1-theorem grade.")
    print("  KILL CONDITION:")
    print("    if P1 fails, or if full galaxy formation forces ordinary fair-sample")
    print("    CDM/NFW halos plus independent R4 response, the old double-counting")
    print("    objection returns.  Otherwise the remaining work is phenomenology.")
    print("exit 0 -- polarised-zero-mode branch closes locally conditional on P1.")


if __name__ == "__main__":
    main()
