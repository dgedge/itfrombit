#!/usr/bin/env python3
r"""Post-final-tier audit of the remaining CC generation-vertex residual.

The final-tier gate (`cc_clause24_proof.py`) removed the old "half + 28-clock"
theorem target.  The one-insertion Rayleigh-Schrodinger operator map already
gives

    Delta(-ln gamma) = alpha0 * C_loop,
    C_loop = (1/3) <1/|sin k|> = 0.3035627046...

The stationary active-demux branch requires

    C_target = 0.3047502315...

so the remaining miss is a real coefficient residual:

    Delta C = 0.0011875269  (0.389672% of C_target),
    rho_Lambda(C_loop) = 1.001872 rho_obs.

This script tests non-tunable next-order mechanisms that are still admissible
after the final-tier proof:

  1. exact finite-coupling evaluation of the same one-insertion Pauli operator;
  2. alpha-scheme substitution / dressed-alpha controls;
  3. the old finite source-shape and determinant/resummation diagnostics.

Exit 0 means the residual is reproduced and those closure routes fail as
derivations.  It does not mean the CC is Locked.
"""

from __future__ import annotations

import math

import numpy as np

from cc_generation_vertex_item115_loop import extrapolate_c, rho_for_c
from cc_generation_vertex_next_order_audit import weighted_generation_c
from register_handoff_form_selection import (
    ALPHA0,
    BASE_GAMMA,
    solve_gamma_for_stationary_queue,
    solve_target,
)


PHYSICAL_ALPHA_INV_CONTROL = 137.035999084


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def finite_coupling_c(n_grid: int, g: float, chunk: int = 40) -> float:
    """Exact symmetrized finite-g coefficient for H = s.sigma + g sigma_x.

    For each k, the lower-band vacuum contribution is -|s +/- g e_x|.  The
    positive barrier coefficient is the negative of the symmetrized vacuum
    shift:

        C_g(k) = (|s+g e_x| + |s-g e_x| - 2|s|) / (2 g^2).

    In the g -> 0 limit this becomes

        (1 - s_x^2/|s|^2) / (2 |s|),

    equal by cubic isotropy to C_loop = (1/3)<1/|s|>.
    """

    ks = (np.arange(n_grid) + 0.5) * 2.0 * math.pi / n_grid - math.pi
    sins = np.sin(ks)
    sin2 = sins * sins
    total = 0.0
    for i in range(0, n_grid, chunk):
        sx = sins[i : i + chunk][:, None, None]
        s2_rest = sin2[None, :, None] + sin2[None, None, :]
        e0 = np.sqrt(sx * sx + s2_rest)
        ep = np.sqrt((sx + g) * (sx + g) + s2_rest)
        em = np.sqrt((sx - g) * (sx - g) + s2_rest)
        total += float(np.sum((ep + em - 2.0 * e0) / (2.0 * g * g)))
    return total / n_grid**3


def continuum_c_loop() -> tuple[float, float]:
    c_lin, c_quad, _ = extrapolate_c([64, 96, 128, 160, 192, 256])
    return 0.5 * (c_lin + c_quad), 0.5 * abs(c_lin - c_quad)


def main() -> None:
    print("CC NEXT-ORDER RESIDUAL AUDIT")

    _, q1_target = solve_target()
    gamma_star = solve_gamma_for_stationary_queue(q1_target)
    delta_req = -math.log(gamma_star / BASE_GAMMA)
    c_target = delta_req / ALPHA0
    c_loop, c_err = continuum_c_loop()
    delta_c = c_target - c_loop
    rho_loop = rho_for_c(c_loop, q1_target)

    print("\n[1] The residual after the final-tier operator map")
    print(f"  Delta_req                         = {delta_req:.15f}")
    print(f"  C_target                          = {c_target:.15f}")
    print(f"  C_loop                            = {c_loop:.15f} +/- {c_err:.2e}")
    print(f"  Delta C = C_target - C_loop        = {delta_c:.15f}")
    print(f"  C_loop/C_target - 1               = {(c_loop / c_target - 1.0) * 100:+.6f}%")
    print(f"  rho(C_loop)/rho_obs               = {rho_loop:.9f}")
    check(c_err < 1e-7, "loop extrapolation error is far below the coefficient residual")
    check(0.001 < delta_c < 0.0013, "the residual is the real -0.39% C-level miss")
    check(abs(math.log(rho_loop)) < 0.002, "the final-tier loop remains a 0.2%-in-rho landing")

    print("\n[2] Exact finite-coupling evaluation of the same operator")
    n_fc = 160
    c0_grid = weighted_generation_c(n_fc, "point")
    finite_rows = [
        ("g = alpha0", ALPHA0),
        ("g = sqrt(alpha0)", math.sqrt(ALPHA0)),
        ("g = 2 sqrt(alpha0)", 2.0 * math.sqrt(ALPHA0)),
    ]
    print(f"  reference point-vertex C(N={n_fc}) = {c0_grid:.12f}")
    print(f"  {'case':<20s} {'C_finite':>14s} {'gain':>13s} {'gain/resid':>12s}")
    max_gain_fraction = 0.0
    for label, g in finite_rows:
        c_finite = finite_coupling_c(n_fc, g)
        gain = c_finite - c0_grid
        gain_fraction = gain / delta_c
        max_gain_fraction = max(max_gain_fraction, gain_fraction)
        print(f"  {label:<20s} {c_finite:14.12f} {gain:13.6e} {gain_fraction:11.4%}")
    check(max_gain_fraction < 0.05, "finite-coupling nonlinearity supplies under 5% of the required residual even in the stress case")
    check(finite_coupling_c(n_fc, math.sqrt(ALPHA0)) < c_target, "finite coupling does not close the C target")

    print("\n[3] Alpha-scheme controls")
    alpha_physical = 1.0 / PHYSICAL_ALPHA_INV_CONTROL
    c_target_physical_alpha = delta_req / alpha_physical
    alpha_required = delta_req / c_loop
    print(f"  canonical alpha0^-1              = {1.0 / ALPHA0:.9f}")
    print(f"  physical/dressed alpha^-1 ctrl   = {PHYSICAL_ALPHA_INV_CONTROL:.9f}")
    print(f"  alpha_eff/alpha0 required         = {alpha_required / ALPHA0:.12f}")
    print(f"  required alpha_eff^-1             = {1.0 / alpha_required:.9f}")
    print(f"  C_target with physical alpha ctrl = {c_target_physical_alpha:.12f}")
    print(f"  miss with physical alpha ctrl     = {(c_loop / c_target_physical_alpha - 1.0) * 100:+.6f}%")
    check(alpha_required / ALPHA0 > 1.003, "closing by alpha alone would need a +0.39% alpha enhancement")
    check(alpha_physical / ALPHA0 < 1.0, "the physical/dressed-alpha control shifts in the wrong direction")
    check(abs((alpha_physical / ALPHA0) - 1.0) < 0.001, "physical-alpha shift is too small even before its sign")

    print("\n[4] Source-shape and resummation diagnostics, retained only as controls")
    point_96 = weighted_generation_c(96, "point")
    edge3_96 = weighted_generation_c(96, "coherent_3_edge")
    edge6_96 = weighted_generation_c(96, "coherent_6_edge")
    print(f"  point vertex N=96                = {point_96:.12f}")
    print(f"  coherent 3-edge star N=96        = {edge3_96:.12f}")
    print(f"  coherent 6-edge star N=96        = {edge6_96:.12f}")
    check(edge3_96 < point_96 and edge6_96 < point_96, "coherent source geometry suppresses rather than raises C")

    resummations = [
        ("linear", c_loop),
        ("-ln(1-alpha C)/alpha", -math.log(1.0 - ALPHA0 * c_loop) / ALPHA0),
        ("-1/2 ln(1-2 alpha C)/alpha", -0.5 * math.log(1.0 - 2.0 * ALPHA0 * c_loop) / ALPHA0),
        ("C/(1-alpha C)", c_loop / (1.0 - ALPHA0 * c_loop)),
        ("C + alpha C^2", c_loop + ALPHA0 * c_loop * c_loop),
    ]
    best = min(resummations, key=lambda row: abs(row[1] - c_target))
    for label, c_eff in resummations:
        print(f"  {label:<32s} C={c_eff:.12f} miss={(c_eff / c_target - 1.0) * 100:+.6f}%")
    print(f"  best ordinary control             = {best[0]}")
    check(abs(best[1] / c_target - 1.0) > 0.001, "ordinary resummations still miss by more than 0.1% in C")

    print(
        """
VERDICT
  The post-final-tier open piece is genuinely narrower than the old theorem
  target.  Billing, half-bubble, and k-uniform vertex normalization are no
  longer the live problem; the remaining object is:

      Delta C = 0.001187526879799
      C_target/C_loop - 1 = 0.3911966%.

  The obvious canon-admissible next-order closures fail:
    * exact finite-coupling evaluation of the same Pauli insertion supplies
      only percent-level pieces of the residual;
    * the dressed/physical-alpha control is too small and shifts in the wrong
      direction;
    * coherent finite source shapes suppress C;
    * ordinary determinant/resummation controls are not selected and still miss.

  Live theorem target:
      derive a genuine next-order generation-vertex renormalization
      delta C = +0.0011875269
  from the monitored QEC/register algebra, or concede the non-horizon CC route
  as a 0.2%-grade near-closure rather than a Locked prediction.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
