#!/usr/bin/env python3
r"""ITEM 123 / 132: nonlinear zero-mode/R4 galaxy-formation bracket.

Question
--------
Does the polarised-zero-mode branch make a concrete galaxy-scale object, or
only a local MOND analogy?

Setup
-----
Use the P1-conditional cored law

    rho_z(r) = A/(r^2+r_c^2),
    A = sqrt(G M_b a0)/(4 pi G),
    r_c = r_M/3,   r_M = sqrt(G M_b/a0),

and let the finite zero-mode reservoir set an outer edge

    M_z(R_t) / M_b = (Omega_zero/Omega_b) / f_ret,

where f_ret is the visible-baryon retention fraction.  This is a nonlinear
formation bracket, not a fit to SPARC: it predicts how the same zero-mode
rest-count reservoir can arrange itself into an R4/MOND-like halo without
double-counting an independent active-R4 mass ledger.

Result
------
For ordinary baryon-retention fractions, the edge lies at virial-like radii,
the halo acceleration is already close to the MOND asymptote over the observed
outer disk, and the residual curve is a finite-edge correction rather than a
new force.  The branch is internally viable at analytic bracket grade.  Full
phenomenology still requires disk geometry, abundance matching, nonlinear
structure, and real galaxy likelihoods.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from item123_cmb_cold_matter_closure_status import OMEGA_B_H2, density_for_power


KPC = 3.0856775814913673e19
G = 6.67430e-11
C = 299_792_458.0
MPC = 3.0856775814913673e22
MSUN = 1.98847e30
H0 = 67.4 * 1000.0 / MPC
A0 = C * H0 / (2.0 * math.pi)


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def m_ratio(y: float) -> float:
    """M_z(y r_M)/M_b for r_c=r_M/3."""

    return y - (1.0 / 3.0) * math.atan(3.0 * y)


def r_mond(m_b: float) -> float:
    return math.sqrt(G * m_b / A0)


def amplitude(m_b: float) -> float:
    return math.sqrt(G * m_b * A0) / (4.0 * math.pi * G)


def core_radius(m_b: float) -> float:
    return r_mond(m_b) / 3.0


def halo_mass(r: float, amp: float, rc: float) -> float:
    return 4.0 * math.pi * amp * (r - rc * math.atan(r / rc))


def solve_y(target: float) -> float:
    lo, hi = 1e-9, 10.0
    while m_ratio(hi) < target:
        hi *= 2.0
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        if m_ratio(mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def g_baryon_spherical(y: float) -> float:
    """Dimensionless baryon acceleration g_b/a0 at radius y r_M."""

    return 1.0 / (y * y)


def g_halo_dimensionless(y: float, y_edge: float) -> float:
    """Dimensionless halo acceleration g_z/a0, truncated at y_edge."""

    y_eff = min(y, y_edge)
    enclosed_ratio = m_ratio(y_eff)
    return enclosed_ratio / (y * y)


def mond_asymptote_dimensionless(y: float) -> float:
    """Deep-MOND halo acceleration sqrt(a0 g_b)/a0."""

    return 1.0 / y


def linspace(lo: float, hi: float, n: int) -> list[float]:
    if n == 1:
        return [lo]
    step = (hi - lo) / (n - 1)
    return [lo + i * step for i in range(n)]


def rms_log_residual(y_values: list[float], y_edge: float) -> float:
    residuals = [
        math.log10(g_halo_dimensionless(y, y_edge)) - math.log10(mond_asymptote_dimensionless(y))
        for y in y_values
    ]
    return math.sqrt(sum(r * r for r in residuals) / len(residuals))


def max_fractional_residual(y_values: list[float], y_edge: float) -> float:
    return max(abs(g_halo_dimensionless(y, y_edge) / mond_asymptote_dimensionless(y) - 1.0) for y in y_values)


def central_density(m_b: float) -> float:
    amp = amplitude(m_b)
    rc = core_radius(m_b)
    return amp / (rc * rc)


def log_slope(xs: list[float], ys: list[float]) -> float:
    lx = [math.log(x) for x in xs]
    ly = [math.log(y) for y in ys]
    mx = sum(lx) / len(lx)
    my = sum(ly) / len(ly)
    return sum((x - mx) * (y - my) for x, y in zip(lx, ly)) / sum((x - mx) ** 2 for x in lx)


@dataclass(frozen=True)
class BracketRow:
    m_b_msun: float
    f_ret: float
    y_edge: float
    r_m_kpc: float
    r_c_kpc: float
    r_edge_kpc: float
    rho0_msun_pc3: float
    rms_outer_dex: float
    max_outer_frac: float


def row(m_b_msun: float, f_ret: float, zero_to_baryon: float) -> BracketRow:
    m_b = m_b_msun * MSUN
    y_edge = solve_y(zero_to_baryon / f_ret)
    r_m = r_mond(m_b)
    rc = core_radius(m_b)
    y_outer = linspace(2.0, min(0.8 * y_edge, 12.0), 200)
    rho0_si = central_density(m_b)
    # kg/m^3 -> Msun/pc^3
    pc = KPC / 1000.0
    rho0_msun_pc3 = rho0_si * pc**3 / MSUN
    return BracketRow(
        m_b_msun=m_b_msun,
        f_ret=f_ret,
        y_edge=y_edge,
        r_m_kpc=r_m / KPC,
        r_c_kpc=rc / KPC,
        r_edge_kpc=y_edge * r_m / KPC,
        rho0_msun_pc3=rho0_msun_pc3,
        rms_outer_dex=rms_log_residual(y_outer, y_edge),
        max_outer_frac=max_fractional_residual(y_outer, y_edge),
    )


def direct_mass_check(m_b_msun: float, y_edge: float, target_ratio: float) -> float:
    m_b = m_b_msun * MSUN
    r_m = r_mond(m_b)
    mass = halo_mass(y_edge * r_m, amplitude(m_b), core_radius(m_b))
    return mass / m_b / target_ratio


def main() -> None:
    print("ITEM 123 / 132: NONLINEAR ZERO-MODE GALAXY BRACKET")
    print("=" * 100)

    omega_nur, omega_zero, omega_dark, zeq = density_for_power(1)
    zero_to_baryon = omega_zero / OMEGA_B_H2
    dark_to_baryon = omega_dark / OMEGA_B_H2
    print("\n[1] Finite reservoir ratios")
    print(f"  omega_b h^2       = {OMEGA_B_H2:.6f}")
    print(f"  omega_zero h^2    = {omega_zero:.6f}")
    print(f"  omega_nuR h^2     = {omega_nur:.6f}")
    print(f"  omega_dark h^2    = {omega_dark:.6f}")
    print(f"  z_eq              = {zeq:.1f}")
    print(f"  zero/baryon       = {zero_to_baryon:.6f}")
    print(f"  dark/baryon       = {dark_to_baryon:.6f}")
    check(4.0 < zero_to_baryon < 4.5, "zero-mode reservoir is the dominant mobile halo share")

    print("\n[2] Edge radii and outer-disk MOND residuals")
    print("  M_b[Msun]    f_ret   r_M[kpc]  r_c[kpc]  y_edge  R_edge[kpc] rho0[Msun/pc3] RMS[dex] max|frac|")
    rows: list[BracketRow] = []
    for m_b_msun in (1e8, 1e9, 1e10, 6e10, 1e11):
        for f_ret in (1.0, 0.5, 0.25, 0.1):
            br = row(m_b_msun, f_ret, zero_to_baryon)
            rows.append(br)
            print(
                f"  {br.m_b_msun:9.2e}  {br.f_ret:5.2f} "
                f"{br.r_m_kpc:8.3f} {br.r_c_kpc:8.3f} {br.y_edge:7.3f} "
                f"{br.r_edge_kpc:11.3f} {br.rho0_msun_pc3:13.4e} "
                f"{br.rms_outer_dex:8.4f} {br.max_outer_frac:9.4f}"
            )
    mw_quarter = next(br for br in rows if br.m_b_msun == 6e10 and br.f_ret == 0.25)
    check(120.0 < mw_quarter.r_edge_kpc < 220.0, "Milky-Way-like 25% retention edge is virial-scale")
    check(mw_quarter.rms_outer_dex < 0.08, "outer-disk acceleration is close to the MOND asymptote before the edge")

    print("\n[3] Nonlinear square-root scaling is preserved across the bracket")
    masses = [1e8 * MSUN, 1e9 * MSUN, 1e10 * MSUN, 1e11 * MSUN]
    amps = [amplitude(m) for m in masses]
    rcs = [core_radius(m) for m in masses]
    rho0s = [central_density(m) for m in masses]
    slope_amp = log_slope(masses, amps)
    slope_rc = log_slope(masses, rcs)
    slope_rho0 = log_slope(masses, rho0s)
    print(f"  A slope versus M_b        = {slope_amp:.6f} (target 1/2)")
    print(f"  r_c slope versus M_b      = {slope_rc:.6f} (target 1/2)")
    print(f"  rho0 slope versus M_b     = {slope_rho0:.6f} (target -1/2)")
    check(abs(slope_amp - 0.5) < 1e-12, "R4 amplitude has the MOND sqrt(M_b) scaling")
    check(abs(slope_rc - 0.5) < 1e-12, "core radius scales as sqrt(M_b)")
    check(abs(slope_rho0 + 0.5) < 1e-12, "central density falls as M_b^-1/2")

    print("\n[4] Mass accounting is exact: edge consumes the finite reservoir")
    for f_ret in (1.0, 0.25, 0.1):
        y = solve_y(zero_to_baryon / f_ret)
        ratio = direct_mass_check(6e10, y, zero_to_baryon / f_ret)
        print(f"  f_ret={f_ret:4.2f}: y_edge={y:8.4f}, direct/target mass ratio={ratio:.12f}")
        check(abs(ratio - 1.0) < 1e-11, "finite edge exactly consumes the requested zero-mode budget")

    print("\n[5] What this does and does not close")
    print("  CLOSED AT BRACKET GRADE:")
    print("    A nonlinear zero-mode/R4 constitutive profile is internally consistent,")
    print("    consumes the finite zero-mode budget at a unique edge, preserves the")
    print("    sqrt(M_b) RAR amplitude, and gives virial-scale edges for plausible")
    print("    baryon-retention fractions.")
    print("  NOT CLOSED:")
    print("    disk geometry, baryon-feedback formation history, scatter, clusters,")
    print("    satellite/subhalo statistics, and full likelihoods.  This is the")
    print("    analytic formation bracket, not a hydrodynamical galaxy model.")
    print("\nexit 0 -- nonlinear polarised-zero-mode galaxy branch is viable at analytic bracket grade.")


if __name__ == "__main__":
    main()
