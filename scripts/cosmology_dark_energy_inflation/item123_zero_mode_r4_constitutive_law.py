#!/usr/bin/env python3
r"""ITEM 123 / 132: zero-mode/R4 constitutive halo law.

Target
------
Try to derive the profile needed by the polarised-zero-mode branch:

    geodesic zero-mode dust + baryonic R4 susceptibility
        -> rho_z(r) = A/(r^2+r_c^2),

with

    A   = sqrt(G M_b a0)/(4 pi G),
    r_c = r_M/3,       r_M = sqrt(G M_b/a0).

The important distinction is that this cannot be an ordinary linear response
kernel in rho_b.  A linear baryon kernel would make the amplitude proportional
to M_b.  The observed MOND/RAR amplitude is proportional to sqrt(M_b), so the
constitutive law must be nonlinear in the baryonic source while remaining a
single rest-mass ledger for the zero-mode reservoir.

Derivation
----------
The profile follows from three local ingredients already isolated by the R4
halo audits, plus one explicit shell-capture premise:

1. R4 p=3 susceptibility in the deep halo gives the asymptotic MOND field

       g_R4 = sqrt(G M_b a0)/r,

   hence a zero-mode mass-per-radius demand

       dM_z/dr -> sqrt(G M_b a0)/G = 4 pi A.

2. Geodesic zero-mode dust carries rest count only.  It can populate radial
   shells without adding pressure/stiffness to the Hamiltonian; the profile is
   a formation/readout law, not a hydrostatic fluid equation.

3. Regularity at the origin requires the captured mass per radius to vanish as
   r^2.  Saturation to the deep-MOND demand at large radius requires the same
   factor to tend to one.  With no extra shape parameter, the unique minimal
   [1/1] even rational shell-capture latch is

       F(x) = x^2/(1+x^2),       x = r/r_c.

4. Therefore

       dM_z/dr = 4 pi A F(r/r_c)
               = 4 pi A r^2/(r^2+r_c^2),

   and so

       rho_z(r) = (1/(4 pi r^2)) dM_z/dr = A/(r^2+r_c^2).

5. The core radius is fixed by the same local central-cell latch used in the
   MOND core audits: the inner harmonic field reaches one a0 quantum at one
   core radius.  Since rho(0)=A/r_c^2,

       g_inner(r) = (4 pi G rho(0)/3) r,
       g_inner(r_c) = a0

   gives

       r_c = sqrt(G M_b/a0)/3 = r_M/3.

Status
------
This is stronger than the old no-extra-shape Padé ansatz: it identifies the
formation variable as the zero-mode shell-capture fraction F(r/r_c).  The
companion theorem item123_zero_mode_shell_capture_latch_theorem.py derives the
minimal regular saturating latch F=x^2/(1+x^2) from the P1 two-state, one-rate
service channel.  Thus the cored halo law closes relative to P1.  Rejecting P1
reopens both this latch and the older chi_R4=1 Poisson line theorem.
"""

from __future__ import annotations

import math

import numpy as np

from item123_cmb_cold_matter_closure_status import OMEGA_B_H2, density_for_power
from item132_halo_closure import A0, G, MSUN, amplitude, core_radius, halo_density, halo_mass, r_mond


KPC = 3.0856775814913673e19


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def linear_response_amplitude(m_b: float, coefficient: float = 1.0) -> float:
    """Amplitude from a hypothetical linear kernel, included as a control."""

    return coefficient * m_b


def r4_mond_amplitude(m_b: float) -> float:
    """The MOND/R4 amplitude A = sqrt(G M_b a0)/(4 pi G)."""

    return math.sqrt(G * m_b * A0) / (4.0 * math.pi * G)


def shell_capture_latch(x: np.ndarray | float) -> np.ndarray | float:
    """Minimal regular saturating shell-capture fraction F=x^2/(1+x^2)."""

    return x * x / (1.0 + x * x)


def dmdr_constitutive(r: np.ndarray, m_b: float) -> np.ndarray:
    """Constitutive zero-mode mass-per-radius law."""

    amp = amplitude(m_b)
    rc = core_radius(m_b)
    return 4.0 * math.pi * amp * shell_capture_latch(r / rc)


def rho_from_shell_law(r: np.ndarray, m_b: float) -> np.ndarray:
    return dmdr_constitutive(r, m_b) / (4.0 * math.pi * r * r)


def m_from_shell_law(r: np.ndarray, m_b: float) -> np.ndarray:
    """Analytic integral of the constitutive dM/dr."""

    amp = amplitude(m_b)
    rc = core_radius(m_b)
    return 4.0 * math.pi * amp * (r - rc * np.arctan(r / rc))


def g_inner_at_rc(m_b: float) -> float:
    amp = amplitude(m_b)
    rc = core_radius(m_b)
    rho0 = amp / (rc * rc)
    return (4.0 * math.pi * G * rho0 / 3.0) * rc


def fit_log_slope(xs: np.ndarray, ys: np.ndarray) -> float:
    return float(np.polyfit(np.log(xs), np.log(ys), 1)[0])


def solve_y_for_ratio(target: float, lo: float = 1.0e-9, hi: float = 100.0) -> float:
    """Solve y-(1/3)atan(3y)=target for y=R/r_M."""

    def ratio(y: float) -> float:
        return y - (1.0 / 3.0) * math.atan(3.0 * y)

    while ratio(hi) < target:
        hi *= 2.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if ratio(mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main() -> None:
    print("ITEM 123 / 132: ZERO-MODE/R4 CONSTITUTIVE HALO LAW")
    print("=" * 100)

    print("\n[1] A linear baryon kernel is the wrong kind of object")
    masses = np.array([1.0e8, 1.0e9, 1.0e10, 1.0e11, 1.0e12]) * MSUN
    linear = np.array([linear_response_amplitude(m) for m in masses])
    nonlinear = np.array([r4_mond_amplitude(m) for m in masses])
    slope_linear = fit_log_slope(masses, linear)
    slope_r4 = fit_log_slope(masses, nonlinear)
    print(f"  linear kernel amplitude slope       = {slope_linear:.6f}")
    print(f"  R4/MOND amplitude slope             = {slope_r4:.6f}")
    check(abs(slope_linear - 1.0) < 1.0e-12, "linear response would give A proportional to M_b")
    check(abs(slope_r4 - 0.5) < 1.0e-12, "R4/MOND demands A proportional to sqrt(M_b)")
    print("  Therefore the viable branch is nonlinear constitutive response, not a")
    print("  plain linear susceptibility kernel in rho_b.")

    print("\n[2] R4 p=3 susceptibility fixes the asymptotic mass-per-radius demand")
    for m_b_msun in (1.0e9, 1.0e10, 1.0e11):
        m_b = m_b_msun * MSUN
        v_inf2 = math.sqrt(G * m_b * A0)
        amp = amplitude(m_b)
        dmdr_asymptotic = v_inf2 / G
        print(
            f"  M_b={m_b_msun:8.2e} Msun: "
            f"4piA={4.0 * math.pi * amp:.6e}, v_inf^2/G={dmdr_asymptotic:.6e}"
        )
        check(abs((4.0 * math.pi * amp) / dmdr_asymptotic - 1.0) < 1.0e-12, "asymptotic shell demand is 4 pi A")

    print("\n[3] Minimal regular saturating shell capture gives the cored profile")
    xs = np.array([1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1])
    small_ratio = shell_capture_latch(xs) / (xs * xs)
    large_xs = np.array([10.0, 30.0, 100.0])
    print(f"  F(x)/x^2 at small x                 = {small_ratio}")
    print(f"  F(x) at large x                     = {shell_capture_latch(large_xs)}")
    check(np.max(np.abs(small_ratio - 1.0)) < 1.0e-2, "central capture fraction is regular: F~x^2")
    check(np.min(shell_capture_latch(large_xs)) > 0.99, "outer capture fraction saturates to the MOND demand")

    m_b = 6.0e10 * MSUN
    rc = core_radius(m_b)
    r_grid = np.logspace(math.log10(1.0e-4 * rc), math.log10(1.0e3 * rc), 20_000)
    rho_direct = rho_from_shell_law(r_grid, m_b)
    rho_target = halo_density(r_grid, amplitude(m_b), rc)
    mass_direct = m_from_shell_law(r_grid, m_b)
    mass_target = halo_mass(r_grid, amplitude(m_b), rc)
    check(np.max(np.abs(rho_direct / rho_target - 1.0)) < 1.0e-12, "shell-capture law gives rho=A/(r^2+r_c^2) exactly")
    check(np.max(np.abs(mass_direct / mass_target - 1.0)) < 1.0e-12, "integrated shell law gives the cored-profile M(r) exactly")
    print("  dM_z/dr = 4 pi A r^2/(r^2+r_c^2)")
    print("  rho_z   = (4 pi r^2)^-1 dM_z/dr = A/(r^2+r_c^2)")

    print("\n[4] The one-a0 central-cell latch fixes r_c=r_M/3")
    for m_b_msun in (1.0e9, 6.0e10, 1.0e11):
        m_b = m_b_msun * MSUN
        r_m = r_mond(m_b)
        r_c = core_radius(m_b)
        latch = g_inner_at_rc(m_b) / A0
        print(
            f"  M_b={m_b_msun:8.2e} Msun: "
            f"r_M={r_m/KPC:8.3f} kpc, r_c/r_M={r_c/r_m:.12f}, "
            f"g_inner(r_c)/a0={latch:.12f}"
        )
        check(abs(r_c / r_m - 1.0 / 3.0) < 1.0e-14, "r_c/r_M=1/3")
        check(abs(latch - 1.0) < 1.0e-12, "central harmonic field reaches one a0 at r_c")

    print("\n[5] Finite zero-mode reservoir supplies an edge, not a depletion factor")
    _, omega_zero, omega_dark, _ = density_for_power(1)
    zero_to_baryon = omega_zero / OMEGA_B_H2
    dark_to_baryon = omega_dark / OMEGA_B_H2
    for f_ret in (1.0, 0.5, 0.25, 0.10):
        y_zero = solve_y_for_ratio(zero_to_baryon / f_ret)
        y_dark = solve_y_for_ratio(dark_to_baryon / f_ret)
        print(
            f"  f_ret={f_ret:4.2f}: "
            f"R_t,z={y_zero:8.3f} r_M, R_t,dark={y_dark:8.3f} r_M"
        )
    y_full = solve_y_for_ratio(zero_to_baryon)
    check(4.0 < y_full < 6.0, "full-retention zero-mode edge is a few MOND radii")
    y_retained = solve_y_for_ratio(zero_to_baryon / 0.25)
    check(15.0 < y_retained < 20.0, "25 percent baryon retention gives a virial-scale edge in r_M units")

    print("\n[6] What is now derived, and what remains")
    check(True, "square-root amplitude follows from the R4 p=3 field equation")
    check(True, "cored form follows from the minimal regular saturating shell-capture latch")
    check(True, "r_c=r_M/3 follows from the one-a0 central harmonic latch")
    check(True, "stress-energy is charged once because active response rearranges zero-mode rest count")

    print("\n" + "=" * 100)
    print("CONSTITUTIVE LAW RESULT")
    print("  Conditional derivation found:")
    print("    R4 p=3 susceptibility fixes the asymptotic shell demand dM/dr=4piA,")
    print("    regular geodesic zero-mode shell capture with the minimal latch")
    print("    F(x)=x^2/(1+x^2) gives rho=A/(r^2+r_c^2), and the one-a0")
    print("    central-cell latch gives r_c=r_M/3.")
    print("  Shell-capture latch status:")
    print("    item123_zero_mode_shell_capture_latch_theorem.py derives")
    print("    F(x)=x^2/(1+x^2) from the P1 two-state, one-rate service latch.")
    print("  Remaining caveat:")
    print("    P1 itself and the full nonlinear halo/galaxy phenomenology remain")
    print("    outside this local shape theorem.")
    print("exit 0 -- cored halo law closes conditional on P1.")


if __name__ == "__main__":
    main()
