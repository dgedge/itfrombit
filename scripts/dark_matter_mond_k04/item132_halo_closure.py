#!/usr/bin/env python3
r"""ITEM 132 constructive closure harness: R4/string halo profile.

This is deliberately narrower than item132_halo_audit.py.  It asks:

1. Can the square-root BTFR amplitude be obtained without inserting
   v^4 = G M_b a0 by hand?
2. Can the r_c = r_M/3 core-radius factor be tied to a substrate boundary
   condition rather than fitted?
3. Does a literal constant-tension Jeans equation derive the cored profile?

Verdict encoded by the checks below:

* Square-root amplitude closes conditionally if the R4 network has the
  critical 3D string action

      E[g] = int |g|^3 / (12 pi G a0) d^3x,

  whose Euler-Lagrange equation is

      div((|g|/a0) g) = 4 pi G rho_b.

  For spherical baryonic mass this gives g = sqrt(G M_b a0)/r, hence
  A = sqrt(G M_b a0)/(4 pi G) and the BTFR.  This is the MOND/AQUAL
  p=3 equation, but the p=3 exponent is now the load-bearing substrate
  target: critical one-dimensional string flux embedded in three space
  dimensions / conformal scale invariance.  The companion audit
  item132_p3_r4_action.py shows that, once this local scale-covariant R4
  action premise is granted, p=3 is forced by D=3 scale covariance and by the
  flat-curve/BTFR scaling.  Companion item132_r4_local_action_lift.py reduces
  the microscopic R4-to-action lift to a local line-density law:
  lambda_R4=(2/3)|g|/a0.  Companion item132_r4_line_density_dynamics.py shows
  generic stable dynamics only gets lambda_R4=(2/3) chi_R4 |g|/a0; chi_R4=1 is
  the remaining normalization target.  Companion item132_chi_unit_poisson.py
  gives chi_R4=1 conditionally from a W=S*C one-jump, count-valued
  boundary-QEC immigration-death ledger with matched creation/erasure rate
  Gamma0.

* The factor 1/3 closes conditionally if the core boundary is where the
  central harmonic R4 acceleration reaches one horizon-walk acceleration
  quantum a0.  With rho = A/(r^2 + r_c^2), the central field is

      g_h(r) = (4 pi G A / (3 r_c^2)) r.

  Setting g_h(r_c) = a0 gives r_c = sqrt(G M_b/a0)/3 = r_M/3.

* A literal Jeans/static-stress equation with constant radial string tension
  does not derive the cored profile.  The cored profile requires a running
  effective radial tension tau(r) ~ v_inf^2 log(r/r_c) at large radius.
  Therefore the profile should be treated as a geometric flux/action result,
  not as a constant-EOS hydrostatic fluid.
"""

from __future__ import annotations

import math
import numpy as np


G = 6.67430e-11
c = 299_792_458.0
MPC = 3.0856775814913673e22
MSUN = 1.98847e30
H0 = 67.4 * 1000.0 / MPC
A0 = c * H0 / (2.0 * math.pi)


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def r_mond(m_b: float) -> float:
    return math.sqrt(G * m_b / A0)


def amplitude(m_b: float) -> float:
    """Halo rho=A/(r^2+r_c^2) amplitude implied by the p=3 R4 flux law."""
    return math.sqrt(G * m_b * A0) / (4.0 * math.pi * G)


def core_radius(m_b: float) -> float:
    """Core radius implied by the one-a0 central harmonic boundary condition."""
    return r_mond(m_b) / 3.0


def halo_mass(r: np.ndarray, amp: float, rc: float) -> np.ndarray:
    return 4.0 * math.pi * amp * (r - rc * np.arctan(r / rc))


def halo_v2(r: np.ndarray, amp: float, rc: float) -> np.ndarray:
    return G * halo_mass(r, amp, rc) / r


def halo_density(r: np.ndarray, amp: float, rc: float) -> np.ndarray:
    return amp / (r * r + rc * rc)


def required_tau_over_v2(x: np.ndarray) -> np.ndarray:
    """Required tau(r)/v_inf^2 for radial stress p_r=-tau rho, p_t=0.

    For rho=A/(r^2+r_c^2), Jeans/static-stress equilibrium gives

      tau' + [2 r_c^2/(r(r^2+r_c^2))] tau = g_h(r).

    In dimensionless x=r/r_c and g_h = v_inf^2/r * [1-atan(x)/x],

      tau/v_inf^2 = (1+x^2)/x^2 *
          integral_0^x u/(1+u^2) * [1-atan(u)/u] du.

    A constant tau would make this flat.  It is not flat; asymptotically it
    grows logarithmically.
    """
    out = []
    grid = np.logspace(-8, math.log10(max(float(np.max(x)), 1e-7)), 200_000)
    integrand = (grid / (1.0 + grid * grid)) * (1.0 - np.arctan(grid) / grid)
    # cumulative trapezoid without scipy dependency
    dx = np.diff(grid)
    cumulative = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * dx)])
    for val in x:
        integ = np.interp(val, grid, cumulative)
        out.append(((1.0 + val * val) / (val * val)) * integ)
    return np.asarray(out)


def main() -> None:
    print("ITEM 132 HALO CLOSURE HARNESS")
    print(f"a0 = c H0 / 2pi = {A0:.6e} m/s^2 at H0=67.4 km/s/Mpc")

    print("\n[1] Square-root amplitude from the critical R4/string p=3 flux law")
    masses = np.array([1e8, 1e9, 1e10, 1e11, 1e12]) * MSUN
    amps = np.array([amplitude(m) for m in masses])
    vinf4 = (4.0 * math.pi * G * amps) ** 2
    slope_a = np.polyfit(np.log(masses), np.log(amps), 1)[0]
    slope_btfr = np.polyfit(np.log(masses), np.log(vinf4), 1)[0]
    check(abs(slope_a - 0.5) < 1e-12, f"A scales as M_b^{slope_a:.3f}; this is the needed square-root amplitude")
    check(abs(slope_btfr - 1.0) < 1e-12, f"v_inf^4 scales as M_b^{slope_btfr:.3f}; BTFR slope is exact")
    for mb in [1e9 * MSUN, 1e11 * MSUN]:
        r = np.logspace(18, 22, 64)
        g_b = G * mb / (r * r)
        g_r4 = math.sqrt(G * mb * A0) / r
        check(np.allclose(g_r4 * g_r4, A0 * g_b), "spherical p=3 equation gives g_R4^2 = a0 g_b exactly")

    print("\n[2] Cored profile, explicit M(r), v_c(r), and the 1/3 core factor")
    mb = 1e11 * MSUN
    amp = amplitude(mb)
    rc = core_radius(mb)
    rM = r_mond(mb)
    vflat2 = math.sqrt(G * mb * A0)
    print(f"  M_b = 1e11 Msun: r_M = {rM / (1000 * 3.0856775814913673e16):.3f} kpc, r_c = {rc / (1000 * 3.0856775814913673e16):.3f} kpc")
    check(abs(rc / rM - 1.0 / 3.0) < 1e-15, "r_c/r_M = 1/3 from g_inner(r_c)=a0")
    g_inner_at_rc = (4.0 * math.pi * G * amp / (3.0 * rc * rc)) * rc
    check(abs(g_inner_at_rc / A0 - 1.0) < 1e-12, "central harmonic acceleration at r_c is exactly one a0 quantum")
    rs = np.array([0.01, 0.1, 1.0, 10.0, 100.0]) * rc
    v2 = halo_v2(rs, amp, rc)
    check(abs(v2[-1] / vflat2 - 1.0) < 0.02, "v_h^2(r) approaches sqrt(G M_b a0) on the flat plateau")
    # Finite-difference Poisson consistency: dM/dr = 4 pi r^2 rho.
    rgrid = np.logspace(math.log10(0.01 * rc), math.log10(100.0 * rc), 10_000)
    dmdr = np.gradient(halo_mass(rgrid, amp, rc), rgrid)
    shell = 4.0 * math.pi * rgrid * rgrid * halo_density(rgrid, amp, rc)
    check(np.max(np.abs(dmdr / shell - 1.0)) < 2e-3, "M(r) differentiates back to rho(r)=A/(r^2+r_c^2)")
    print("  rho_h(r) = A/(r^2+r_c^2)")
    print("  M_h(r)   = 4 pi A [r - r_c atan(r/r_c)]")
    print("  v_h^2(r) = 4 pi G A [1 - (r_c/r) atan(r/r_c)]")

    print("\n[3] Jeans/static-stress check: constant radial string tension still fails")
    xs = np.array([0.1, 1.0, 3.0, 10.0, 100.0])
    tau = required_tau_over_v2(xs)
    for x, t in zip(xs, tau):
        print(f"  required tau({x:5.1f} r_c) / v_inf^2 = {t:.4f}")
    check(tau[-1] > 5.0 * tau[1], "required effective tension runs with radius; it is not a constant string EOS")
    growth = tau[-1] - tau[-2]
    check(abs(growth - math.log(xs[-1] / xs[-2])) < 0.25, "large-radius tau/v_inf^2 has logarithmic growth")

    print("\n" + "=" * 92)
    print("CONDITIONAL CLOSURE:")
    print("  * Square-root amplitude closes if R4 is governed by the critical p=3 string/action")
    print("    div((|g|/a0) g)=4 pi G rho_b.  That gives A=sqrt(G M_b a0)/(4 pi G).")
    print("  * The 1/3 core factor closes if the core boundary is the one-a0 central harmonic")
    print("    acceleration condition.  The 1/3 is the spherical Poisson 4pi/3 coefficient.")
    print("  * Companion item132_p3_r4_action.py shows p=3 is forced once R4")
    print("    has a local scale-covariant action in D=3. Companion")
    print("    item132_r4_local_action_lift.py reduces that action premise to")
    print("    local line density, and item132_r4_line_density_dynamics.py shows")
    print("    generic dynamics leaves lambda_R4=(2/3) chi_R4 |g|/a0.")
    print("    item132_chi_unit_poisson.py gives chi_R4=1 if the W=S*C boundary-QEC")
    print("    line ledger has matched creation/erasure rates and Poisson count support.")
    print("    What remains open for full closure is deriving those rate-matching/count")
    print("    lemmas intrinsically, plus the one-a0 core boundary rule.")
    print("    A literal constant-tension Jeans equation does not do it.")
    print("=" * 92)


if __name__ == "__main__":
    main()
