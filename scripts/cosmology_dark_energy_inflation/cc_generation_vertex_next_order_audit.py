#!/usr/bin/env python3
r"""Next-order audit for the CC generation-vertex residual.

Prior state:
    The item-115 chirality/Weyl loop gives

        C_loop = (1/3)<1/|sin k|> = 0.3035627046

    for the generation-vertex barrier

        Delta(-ln gamma) = alpha0 * C.

    The stationary active-demux CC branch requires

        C_target = 0.3047502315,

    so the evaluated one-loop coefficient is low by 0.3897%, while already
    giving rho_Lambda = 1.00187 rho_obs.

This script asks what the remaining 0.39% can and cannot be:

  1. It is not numerical convergence.
  2. It is not supplied by a finite edge-star source shape; coherent edge
     smearing suppresses the loop instead of raising it, while incoherent
     edge-label billing collapses back to the point vertex by |e^{ik/2}|^2=1.
  3. Generic determinant/resummation forms improve the number but do not select
     a unique correction.
  4. Before the final-tier proof, the residual was naturally parameterized as an
     external-leg / vertex normalization,

        C_eff = C_loop * (1 + alpha0 * z),

     with z_req = 0.535939.  The canon-shaped candidate

        z = 1/2 + 1/28 = 15/28

     lands at rho = 1.0000008 rho_obs.

     This is now diagnostic history.  The later final-tier proof
     (`cc_clause24_proof.py`) shows that the one-insertion RS shift already is
     C_loop: the half is not an independent normalization, and the 28-clock
     factor is not a live scheduler theorem target.

exit 0 means the narrowing is reproduced.  It does not mean the cosmological
constant is Locked.
"""

from __future__ import annotations

import math

import numpy as np

from register_handoff_form_selection import (
    ALPHA0,
    BASE_GAMMA,
    solve_gamma_for_stationary_queue,
    solve_target,
)
from cc_generation_vertex_item115_loop import extrapolate_c, rho_for_c


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def weighted_generation_c(n_grid: int, mode: str, chunk: int = 36) -> float:
    """Single-axis Pauli generation bubble with optional finite source form.

    The local point vertex has form factor 1.  A coherent 3-edge star uses the
    positive incident-edge amplitude (e^{ikx/2}+e^{iky/2}+e^{ikz/2})/3.  A
    coherent 6-edge star uses the signed-edge amplitude
    (cos(kx/2)+cos(ky/2)+cos(kz/2))/3.  Both are normalized to one at k=0.

    Incoherent edge billing is not implemented as a separate weight because each
    labelled edge has |e^{ik_i/2}|^2=1, so the average over labels is exactly the
    point vertex.
    """

    ks = (np.arange(n_grid) + 0.5) * 2.0 * math.pi / n_grid - math.pi
    sin2 = np.sin(ks) ** 2
    phases = np.exp(0.5j * ks)
    cos_half = np.cos(0.5 * ks)

    total = 0.0
    for i in range(0, n_grid, chunk):
        sx = sin2[i : i + chunk][:, None, None]
        s2tot = sx + sin2[None, :, None] + sin2[None, None, :]
        base = (1.0 - sx / s2tot) / np.sqrt(s2tot) / 2.0

        if mode == "point":
            weight = 1.0
        elif mode == "coherent_3_edge":
            amp = (
                phases[i : i + chunk][:, None, None]
                + phases[None, :, None]
                + phases[None, None, :]
            ) / 3.0
            weight = np.abs(amp) ** 2
        elif mode == "coherent_6_edge":
            amp = (
                cos_half[i : i + chunk][:, None, None]
                + cos_half[None, :, None]
                + cos_half[None, None, :]
            ) / 3.0
            weight = amp**2
        else:
            raise ValueError(mode)

        total += float(np.sum(base * weight))

    return total / n_grid**3


def main() -> None:
    print("CC GENERATION-VERTEX NEXT-ORDER AUDIT")

    _, q1_target = solve_target()
    gamma_star = solve_gamma_for_stationary_queue(q1_target)
    c_target = -math.log(gamma_star / BASE_GAMMA) / ALPHA0
    c_lin, c_quad, rows = extrapolate_c([64, 96, 128, 160, 192, 256])
    c_loop = 0.5 * (c_lin + c_quad)
    c_err = abs(c_lin - c_quad) / 2.0
    rho_loop = rho_for_c(c_loop, q1_target)

    print("\n[1] Target and verified one-loop residual")
    print(f"  C_target                      = {c_target:.12f}")
    print(f"  C_loop extrapolated            = {c_loop:.12f} +/- {c_err:.2e}")
    print(f"  C_loop/C_target - 1            = {(c_loop / c_target - 1.0) * 100:+.6f}%")
    print(f"  rho(C_loop)/rho_obs            = {rho_loop:.9f}")
    print("  convergence rows:")
    for n, c in rows:
        print(f"    N={n:<3d} C={c:.12f} rho={rho_for_c(c, q1_target):.9f}")
    check(c_err < 1e-7, "extrapolation uncertainty is far below the residual")
    check(abs(c_loop / c_target - 1.0) > 0.003, "the 0.39% coefficient residual is real")
    check(abs(math.log(rho_loop)) < 0.002, "the existing loop already lands at 0.2%-in-rho grade")

    print("\n[2] Finite source-shape tests")
    point_96 = weighted_generation_c(96, "point")
    edge3_96 = weighted_generation_c(96, "coherent_3_edge")
    edge6_96 = weighted_generation_c(96, "coherent_6_edge")
    print(f"  point vertex, N=96             = {point_96:.12f}")
    print(f"  coherent 3-edge star, N=96     = {edge3_96:.12f}")
    print(f"  coherent 6-edge star, N=96     = {edge6_96:.12f}")
    print("  incoherent labelled-edge billing = point vertex exactly (|e^{ik/2}|^2=1)")
    check(abs(point_96 - dict(rows)[96]) < 1e-12, "weighted point vertex matches the prior loop")
    check(edge3_96 < 0.65 * point_96, "coherent 3-edge smearing suppresses the coefficient")
    check(edge6_96 < 0.50 * point_96, "coherent 6-edge smearing suppresses the coefficient further")
    check(edge3_96 < c_target and edge6_96 < c_target, "finite coherent source shape cannot supply the positive residual")

    print("\n[3] Resummation families")
    resummations = [
        ("linear one-loop", c_loop),
        ("-ln(1-alpha C)/alpha", -math.log(1.0 - ALPHA0 * c_loop) / ALPHA0),
        ("ln(1+alpha C)/alpha", math.log(1.0 + ALPHA0 * c_loop) / ALPHA0),
        ("-1/2 ln(1-2 alpha C)/alpha", -0.5 * math.log(1.0 - 2.0 * ALPHA0 * c_loop) / ALPHA0),
        ("C/(1-alpha C)", c_loop / (1.0 - ALPHA0 * c_loop)),
        ("C + alpha C^2/2", c_loop + ALPHA0 * c_loop * c_loop / 2.0),
        ("C + alpha C^2", c_loop + ALPHA0 * c_loop * c_loop),
    ]
    print(f"  {'form':<34s} {'C_eff':>14s} {'miss_C':>12s} {'rho':>12s}")
    best = None
    for name, c_eff in resummations:
        rho = rho_for_c(c_eff, q1_target)
        miss = c_eff / c_target - 1.0
        if best is None or abs(math.log(rho)) < abs(math.log(best[2])):
            best = (name, c_eff, rho)
        print(f"  {name:<34s} {c_eff:14.12f} {miss * 100:+11.6f}% {rho:12.9f}")
    assert best is not None
    check(best[0] in {"-1/2 ln(1-2 alpha C)/alpha", "C/(1-alpha C)", "C + alpha C^2"}, "ordinary resummations improve but do not identify the correction")
    check(abs(best[1] / c_target - 1.0) > 0.001, "ordinary resummations still miss by more than 0.1% in coefficient")

    print("\n[4] External-leg / vertex-normalization parameterization")
    z_req_linear = (c_target / c_loop - 1.0) / ALPHA0
    z_req_exp = math.log(c_target / c_loop) / ALPHA0
    print(f"  z_req in C(1+alpha z)           = {z_req_linear:.12f}")
    print(f"  z_req in C exp(alpha z)         = {z_req_exp:.12f}")
    print(f"  z_req - 1/2                     = {z_req_linear - 0.5:.12f}")
    print(f"  1/28                            = {1.0 / 28.0:.12f}")
    print(f"  z_28 = 1/2 + 1/28 = 15/28       = {15.0 / 28.0:.12f}")

    z_candidates = [
        ("required z", z_req_linear),
        ("half external leg", 0.5),
        ("half + 28-clock", 0.5 + 1.0 / 28.0),
        ("4/7 straddle", 4.0 / 7.0),
        ("golden-ratio control", (math.sqrt(5.0) - 1.0) / 2.0),
    ]
    print(f"\n  {'z candidate':<24s} {'z':>14s} {'C_lin':>14s} {'rho_lin':>12s} {'rho_exp':>12s}")
    for name, z in z_candidates:
        c_linear = c_loop * (1.0 + ALPHA0 * z)
        c_exp = c_loop * math.exp(ALPHA0 * z)
        print(
            f"  {name:<24s} {z:14.12f} {c_linear:14.12f} "
            f"{rho_for_c(c_linear, q1_target):12.9f} {rho_for_c(c_exp, q1_target):12.9f}"
        )

    c_half28 = c_loop * (1.0 + ALPHA0 * (0.5 + 1.0 / 28.0))
    rho_half28 = rho_for_c(c_half28, q1_target)
    check(abs(math.log(rho_half28)) < 1e-5, "half-plus-28 external normalization lands at ppm-in-rho grade")
    check(abs((z_req_linear - 0.5) - 1.0 / 28.0) < 3e-4, "the post-half residual is numerically a 28-clock unit")

    print(
        """
VERDICT
  The next-order audit does not promote the cosmological constant.  Its
  half-plus-28 parameterization is retained only as diagnostic history after
  the final-tier proof.

  Excluded or disfavoured:
    * numerical error: the loop extrapolation error is orders below the residual;
    * finite coherent edge-star geometry: it suppresses the loop instead of
      raising it;
    * generic determinant/resummation forms: they improve the number but do not
      select a unique canon object.

  Historical diagnostic:
      C_eff = C_loop * (1 + alpha0 * z), with z_req = 0.535939.

  The canon-shaped candidate z = 1/2 + 1/28 = 15/28 lands at
      rho_Lambda = 1.0000008 rho_obs
  under the linear external-leg reading.  But `cc_clause24_proof.py` supersedes
  this as a theorem target: there is no independent half external leg to prove.

  Live target after the final-tier gate:
      derive the real next-order residual
      delta C = C_target - C_loop = 0.0011875269...

  from the monitored QEC/register algebra, or leave the non-horizon CC route as
  a 0.2%-grade near-closure rather than a Locked cosmological constant.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
