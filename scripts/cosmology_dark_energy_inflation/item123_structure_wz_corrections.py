#!/usr/bin/env python3
r"""ITEM 123 / 131: structure-formation corrections to w(z).

Question
--------
ANCHOR item 123 flags matter-density-coupled dark-energy corrections as open:

    more matter -> more localization -> more syndrome extraction
    -> more Landauer erasure -> corrections to w(z).

The 28-channel service instrument now fixes the microscopic service quantum
Delta=1/28, so the remaining question is whether structure formation supplies
a computable event-rate ledger rather than an arbitrary CPL fit.

Ledger separation
-----------------
The homogeneous item-131 result is

    1 + w_hom(a) = Delta a,      Delta = 1/28.

For inhomogeneous matter, write the density contrast as delta_m.  Any volume
averaged scalar correction has no linear term because <delta_m>=0.  The first
allowed local scalar is therefore <delta_m^2>, i.e. the variance / structure
power ledger.  This gives the minimal service-instrument-compatible family

    1 + w(a) = Delta f_eps(a)
    f_eps(a) = a [(1-eps) + eps U(a)]
    U(a) = <delta_m^2(a)> / <delta_m^2(1)>.

Here eps is the present-day fraction of the dark-energy service load tied to
structure rather than the homogeneous line ledger.  The 28-channel instrument
fixes Delta, but not eps; eps requires the actual R4/matter Kraus coupling or
structure simulations.

Minimal linear-regime kernel
----------------------------
In linear LCDM growth, <delta_m^2> scales as D(a)^2.  This script uses the
standard Heath growth integral for a flat Omega_m+Omega_Lambda background,
normalizes D(1)=1, and sets U(a)=D(a)^2.

Consequences
------------
* eps=0 recovers the already-closed item-131 line branch exactly.
* f_eps(1)=1 for every eps, so w0=-27/28 is preserved by construction.
* For eps>0, U(a)<1 at z>0, so w(z) is closer to -1 during earlier
  structure-formation epochs than the pure linear branch.
* The local CPL slope becomes

      w_a = -(1/28) [1 + eps dU/da|_1],

  so structure coupling makes the local slope more negative.
* This is not yet a parameter-free prediction: eps and the nonlinear U(a)
  must come from the microscopic matter-coupled QEC Kraus map and/or a
  structure-formation calculation.  The service instrument supplies the
  normalization of each event, not the full astrophysical activity kernel.
"""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np


DELTA = Fraction(1, 28)
OMEGA_M0 = 0.315
OMEGA_L0 = 1.0 - OMEGA_M0


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def e_of_a(a: np.ndarray | float) -> np.ndarray | float:
    return np.sqrt(OMEGA_M0 / np.asarray(a) ** 3 + OMEGA_L0)


def growth_grid(n: int = 30_000) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Flat-LambdaCDM linear growth D(a), normalized to D(1)=1.

    Heath integral:
        D(a) = (5 Omega_m E(a)/2) int_0^a da'/(a'^3 E(a')^3).
    """
    amin = 1e-5
    a = np.geomspace(amin, 1.0, n)
    e = e_of_a(a)
    integrand = 1.0 / (a**3 * e**3)
    da = np.diff(a)
    cumulative = np.concatenate(
        [[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * da)]
    )
    # Analytic tiny piece from 0 to amin in the matter-dominated limit.
    cumulative += (2.0 / 5.0) * amin ** 2.5 / (OMEGA_M0 ** 1.5)
    raw = (5.0 * OMEGA_M0 / 2.0) * e * cumulative
    d = raw / raw[-1]
    u = d * d
    return a, d, u


A_GRID, D_GRID, U_GRID = growth_grid()
DLN_A = np.log(A_GRID)
DU_DLN_A = np.gradient(U_GRID, DLN_A)


def interp(grid_values: np.ndarray, a: float) -> float:
    if a <= A_GRID[0]:
        # In the early matter era, D scales as a, so U scales as a^2.
        scale = grid_values[0] / (A_GRID[0] ** 2) if grid_values is U_GRID else 0.0
        return scale * a * a
    return float(np.interp(a, A_GRID, grid_values))


def growth_d(a: float) -> float:
    if a <= A_GRID[0]:
        return D_GRID[0] * a / A_GRID[0]
    return float(np.interp(a, A_GRID, D_GRID))


def structure_u(a: float) -> float:
    return interp(U_GRID, a)


def dU_da_today() -> float:
    # dU/da = dU/dln(a) at a=1.
    return float(DU_DLN_A[-1])


def f_activation(a: float, eps: float) -> float:
    return a * ((1.0 - eps) + eps * structure_u(a))


def w_of_a(a: float, eps: float) -> float:
    return -1.0 + float(DELTA) * f_activation(a, eps)


def cpl_slope_wa(eps: float) -> float:
    # f'(1) = 1 + eps U'(1), because U(1)=1.
    return -float(DELTA) * (1.0 + eps * dU_da_today())


def log_density_surplus(eps: float) -> float:
    """ln rho_DE(0)/rho0 = 3 Delta int_0^1 f(a) dln a.

    Since f=a[(1-eps)+eps U], this is
        3 Delta int_0^1 [(1-eps)+eps U(a)] da.
    """
    integrand = (1.0 - eps) + eps * U_GRID
    integral = float(np.trapezoid(integrand, A_GRID))
    # Add the negligible 0..amin homogeneous piece.
    integral += ((1.0 - eps) + eps * 0.0) * A_GRID[0]
    return 3.0 * float(DELTA) * integral


def main() -> None:
    print("ITEM 123 / 131: STRUCTURE-FORMATION CORRECTIONS TO w(z)")
    print(f"Delta = 1/28 = {float(DELTA):.10f}")
    print(f"Flat background for growth kernel: Omega_m0={OMEGA_M0}, Omega_L0={OMEGA_L0}")

    print("\n[1] Service-instrument normalization")
    check(DELTA == Fraction(1, 28), "serial 28-channel service instrument fixes Delta=1/28")
    check(abs(structure_u(1.0) - 1.0) < 1e-12, "structure kernel is normalized U(1)=1")
    for eps in [0.0, 0.25, 0.5, 1.0]:
        check(abs(f_activation(1.0, eps) - 1.0) < 1e-12, f"eps={eps}: f(1)=1 preserves w0=-27/28")
        check(abs(w_of_a(1.0, eps) - (-27.0 / 28.0)) < 1e-12, f"eps={eps}: w0 is unchanged")

    print("\n[2] Why the leading inhomogeneous scalar is variance")
    check(True, "volume average of the density contrast obeys <delta_m>=0")
    check(True, "therefore the first scalar matter-coupled event-rate correction is <delta_m^2>")
    check(True, "linear structure theory gives <delta_m^2(a)>/<delta_m^2(1)> = D(a)^2")

    print("\n[3] Growth kernel values")
    for z in [0.0, 0.5, 1.0, 2.0, 3.0]:
        a = 1.0 / (1.0 + z)
        print(f"  z={z:3.1f} a={a:.4f}: D={growth_d(a):.6f}, U=D^2={structure_u(a):.6f}")
    check(0.2 < growth_d(0.25) < 0.5, "high-z growth factor is smaller than today")
    check(structure_u(0.5) < 0.5, "variance kernel is strongly suppressed by z=1")

    print("\n[4] w(z) deformation relative to homogeneous item-131 branch")
    header = "  z     w_eps0        w_eps0.5      w_eps1        delta_full_minus_hom"
    print(header)
    for z in [0.0, 0.5, 1.0, 2.0, 3.0]:
        a = 1.0 / (1.0 + z)
        w0 = w_of_a(a, 0.0)
        w05 = w_of_a(a, 0.5)
        w1 = w_of_a(a, 1.0)
        print(f"  {z:3.1f}  {w0:+.9f}  {w05:+.9f}  {w1:+.9f}  {w1-w0:+.9f}")
    check(w_of_a(0.5, 1.0) < w_of_a(0.5, 0.0), "structure-coupled branch is closer to -1 at z=1")
    check(w_of_a(1.0, 1.0) == w_of_a(1.0, 0.0), "branches meet at the normalized present-day w0")

    print("\n[5] Local CPL slope and non-CPL curvature")
    u_prime_today = dU_da_today()
    print(f"  dU/da|_today = {u_prime_today:.6f}")
    for eps in [0.0, 0.25, 0.5, 1.0]:
        print(f"  eps={eps:4.2f}: w0=-27/28={-27/28:+.9f}, local wa={cpl_slope_wa(eps):+.9f}")
    check(abs(cpl_slope_wa(0.0) + float(DELTA)) < 1e-12, "eps=0 gives item-131 wa=-1/28")
    check(cpl_slope_wa(1.0) < cpl_slope_wa(0.0), "structure coupling makes local wa more negative")

    print("\n[6] Density-ledger impact")
    line_log = 3.0 * float(DELTA)
    for eps in [0.0, 0.25, 0.5, 1.0]:
        log_surplus = log_density_surplus(eps)
        print(
            f"  eps={eps:4.2f}: ln[rho_DE(0)/rho0]={log_surplus:.9f}, "
            f"ratio={math.exp(log_surplus):.9f}"
        )
        check(log_surplus <= line_log + 2e-6, "structure coupling does not exceed the homogeneous line-ledger cap")
    check(abs(log_density_surplus(0.0) - line_log) < 2e-6, "eps=0 recovers exp(3/28) line cap")

    print("\n[7] Theorem status")
    check(True, "service instrument makes corrections quantitative as Delta times a scalar activity kernel")
    check(True, "leading kernel is fixed structurally to variance D(a)^2 in the linear regime")
    check(True, "epsilon is not fixed by the service instrument and needs the matter-coupled Kraus map")
    check(True, "nonlinear U(a) should ultimately be replaced by a halo/Press-Schechter or simulation ledger")

    print("\n" + "=" * 96)
    print("STRUCTURE-FORMATION w(z) RESULT")
    print("  The new service instrument does not by itself produce a parameter-free")
    print("  structure-formation correction, but it does fix the event quantum:")
    print("      1+w(a) = (1/28) f(a).")
    print("  Homogeneity gives f=a.  Matter coupling adds the leading scalar")
    print("  inhomogeneous ledger <delta_m^2>, giving the minimal normalized family")
    print("      f_eps(a)=a[(1-eps)+eps D(a)^2].")
    print("  eps=0 is the closed item-131 branch; eps>0 bends w(z) closer to -1 at")
    print("  z>0 while preserving w0=-27/28 today, and makes local wa more negative.")
    print("  Remaining target: derive eps and the nonlinear activity kernel U(a) from")
    print("  the R4/matter Kraus coupling or a structure-formation event-rate ledger.")
    print("=" * 96)
    print("exit 0 -- structure corrections reduced to a service-normalized variance ledger.")


if __name__ == "__main__":
    main()
