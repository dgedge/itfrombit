#!/usr/bin/env python3
r"""Audit the fixed-connected-object routes for the remaining CC residual.

The final-tier generation-vertex map gives

    C_loop   = 0.3035627046...
    C_target = 0.3047502315...

so the live residual is

    delta C = +0.0011875269...

If this is a genuine next-order correction to

    Delta(-ln gamma) = alpha0 * C,

then the order-alpha0 connected object must supply

    C2_required = delta C / alpha0 = 0.16269118...

This script checks the three obvious "new fixed object" families now named in
canon:

  1. a connected two-loop gauge-web vertex;
  2. a non-Gaussian bath cumulant of the already-forced chi/W Pauli loop;
  3. a QEC/Kraus normal-ordering term.

Exit 0 means: none of these is presently a canon-derived closure.  It is a
negative/decision audit, not a cosmological-constant promotion.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import mpmath as mp
import numpy as np

from cc_generation_vertex_item115_loop import extrapolate_c, rho_for_c
from register_handoff_form_selection import (
    ALPHA0,
    BASE_GAMMA,
    solve_gamma_for_stationary_queue,
    solve_target,
)


mp.mp.dps = 40


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def target_numbers() -> tuple[float, float, float, float]:
    _, q1_target = solve_target()
    gamma_star = solve_gamma_for_stationary_queue(q1_target)
    c_target = -math.log(gamma_star / BASE_GAMMA) / ALPHA0
    c_lin, c_quad, _ = extrapolate_c([64, 96, 128, 160, 192, 256])
    c_loop = 0.5 * (c_lin + c_quad)
    delta_c = c_target - c_loop
    return q1_target, c_target, c_loop, delta_c


def graph_invariants() -> tuple[int, int, int, int, int, int]:
    """Return V, E, beta1, tau, dual det, bare connected trace."""

    a = np.zeros((16, 16), dtype=int)

    def edge(i: int, j: int) -> None:
        a[i, j] = a[j, i] = 1

    for k in range(8):
        edge(k, (k + 1) % 8)
        edge(8 + k, 8 + ((k + 1) % 8))
    edge(1, 8)
    edge(9, 0)

    v = 16
    e = int(a.sum() // 2)
    beta1 = e - v + 1
    lap = np.diag(a.sum(axis=1)) - a
    tau = round(float(np.linalg.det(lap[1:, 1:])))
    dual = np.array([[8, 0, -1], [0, 8, -1], [-1, -1, 4]], dtype=float)
    dual_det = round(float(np.linalg.det(dual)))
    tr2 = int(np.trace(np.linalg.matrix_power(a, 2)))
    tr4 = int(np.trace(np.linalg.matrix_power(a, 4)))
    connected_trace = tr4 - tr2 * tr2
    return v, e, beta1, tau, dual_det, connected_trace


@dataclass(frozen=True)
class BathStats:
    n_grid: int
    mean: float
    raw2: float
    var: float
    kappa3: float
    kappa4: float


def bath_stats(n_grid: int, chunk: int = 48) -> BathStats:
    """Cumulants of the leading one-point Pauli loop density.

    X(k) = (1 - sin^2(k_x)/|sin k|^2) / (2 |sin k|)

    The mean is C_loop on one cubic axis.  If the same bath variable supplied a
    non-Gaussian cumulant correction to log <exp(gX)>, the next coefficient
    would be kappa_2/2.
    """

    ks = (np.arange(n_grid) + 0.5) * 2.0 * math.pi / n_grid - math.pi
    sin2 = np.sin(ks) ** 2
    raw = [0.0, 0.0, 0.0, 0.0]
    count = 0
    for i in range(0, n_grid, chunk):
        sx = sin2[i : i + chunk][:, None, None]
        s2 = sx + sin2[None, :, None] + sin2[None, None, :]
        x = (1.0 - sx / s2) / (2.0 * np.sqrt(s2))
        count += x.size
        xp = x
        for p in range(4):
            raw[p] += float(np.sum(xp))
            xp = xp * x
    m1, m2, m3, m4 = (r / count for r in raw)
    var = m2 - m1 * m1
    k3 = m3 - 3.0 * m1 * m2 + 2.0 * m1**3
    central4 = m4 - 4.0 * m1 * m3 + 6.0 * m1 * m1 * m2 - 3.0 * m1**4
    k4 = central4 - 3.0 * var * var
    return BathStats(n_grid, m1, m2, var, k3, k4)


def qed_two_loop_vertex_half() -> float:
    """Half the magnitude of the textbook two-loop electron g-2 coefficient."""

    pi = mp.pi
    c2 = (
        mp.mpf(197) / 144
        + pi**2 / 12
        - (pi**2 / 2) * mp.log(2)
        + mp.mpf(3) * mp.zeta(3) / 4
    )
    return float(abs(c2) / 2)


def main() -> None:
    print("CC FIXED-CONNECTED-OBJECT AUDIT")

    q1_target, c_target, c_loop, delta_c = target_numbers()
    c2_req = delta_c / ALPHA0
    print("\n[0] Target")
    print(f"  C_target                         = {c_target:.15f}")
    print(f"  C_loop                           = {c_loop:.15f}")
    print(f"  delta C required                  = {delta_c:.15f}")
    print(f"  C2_required = deltaC/alpha0       = {c2_req:.12f}")
    print(f"  rho(C_loop)/rho_obs               = {rho_for_c(c_loop, q1_target):.9f}")
    check(0.162 < c2_req < 0.164, "required next-order coefficient is the known 0.1627 target")

    print("\n[1] Connected two-loop gauge-web vertex")
    v, e, beta1, tau, dual_det, conn = graph_invariants()
    print(f"  16-node graph: V={v}, E={e}, beta1={beta1}, tau={tau}, dual_det={dual_det}")
    print(f"  bare connected trace Tr[A^4]-Tr[A^2]^2 = {conn}")
    print("  solid invariant: tau=240; unpinned object: the Green function in Tr[(AG)^4].")
    graph_controls = [
        ("1/tau", 1.0 / tau),
        ("beta1/tau", beta1 / tau),
        ("V/tau", v / tau),
        ("E/tau", e / tau),
        ("sqrt(beta1/tau)", math.sqrt(beta1 / tau)),
    ]
    print(f"  {'graph-only control':<22s} {'value':>12s} {'value/C2_req':>14s}")
    for label, val in graph_controls:
        print(f"  {label:<22s} {val:12.9f} {val / c2_req:13.3f}")
    half_qed = qed_two_loop_vertex_half()
    print(f"  |Sommerfield-Petermann c2|/2      = {half_qed:.12f}")
    print(f"  relative to target                = {(half_qed / c2_req - 1.0) * 100:+.3f}%")
    check(tau == 240 and dual_det == 240, "the graph's solid 240 invariant is reproduced")
    check(conn == -1156, "bare-adjacency connected trace is not the claimed -240 object")
    check(abs(half_qed / c2_req - 1.0) > 0.005, "QED two-loop half is only a near number, not an exact hit")
    print("  Verdict: no closure.  The near QED vertex number is the wrong Green's function;")
    print("  the canon graph supplies counts, not a fixed generation-vertex tensor.")

    print("\n[2] Non-Gaussian bath cumulant of the forced Pauli loop")
    rows = [bath_stats(n) for n in (64, 96, 128, 160, 192)]
    print(f"  {'N':>4s} {'mean':>13s} {'var/2':>13s} {'raw2/2':>13s} {'kappa3/6':>13s}")
    for row in rows:
        print(
            f"  {row.n_grid:4d} {row.mean:13.9f} {row.var / 2:13.9f} "
            f"{row.raw2 / 2:13.9f} {row.kappa3 / 6:13.9f}"
        )
    finest = rows[-1]
    print(f"  required C2                       = {c2_req:.9f}")
    print(f"  kappa2/2 at N={finest.n_grid:<3d}              = {finest.var / 2:.9f}")
    print(f"  raw second moment/2 at N={finest.n_grid:<3d}    = {finest.raw2 / 2:.9f}")
    check(finest.var / 2 < 0.13 * c2_req, "the genuine second cumulant is an order too small")
    check(finest.raw2 / 2 < 0.45 * c2_req, "even the uncentered second moment is far too small")
    check(rows[-1].kappa3 > 1.7 * rows[0].kappa3, "higher cumulants are cutoff-sensitive, not fixed constants")
    print("  Verdict: no closure.  The fixed second cumulant is too small; higher")
    print("  cumulants are UV/cutoff sensitive for this Weyl-point integrand.")

    print("\n[3] QEC/Kraus normal-ordering term")
    r_walk = 2.0 / 3.0
    r_uniform8 = 5.0 / 8.0
    kraus_controls = [
        ("exact Lindblad exponential", 0.0),
        ("walk-active amplitude Kraus r^2/4", r_walk * r_walk / 4.0),
        ("walk-active no-jump log r^2/2", r_walk * r_walk / 2.0),
        ("8-bit amplitude Kraus r^2/4", r_uniform8 * r_uniform8 / 4.0),
        ("8-bit no-jump log r^2/2", r_uniform8 * r_uniform8 / 2.0),
    ]
    print(f"  {'normal-ordering control':<36s} {'C2':>12s} {'miss':>12s}")
    best = None
    for label, val in kraus_controls:
        miss = val / c2_req - 1.0
        print(f"  {label:<36s} {val:12.9f} {miss:11.2%}")
        if best is None or abs(miss) < abs(best[2]):
            best = (label, val, miss)
    assert best is not None
    check(abs(best[2]) > 0.15, "standard Kraus normal-ordering controls miss by far more than the closure tolerance")
    print("  Verdict: no closure.  The channel has O(tau^2) discretization choices,")
    print("  not a unique dimensionless scalar.  Picking a value between r^2/4 and")
    print("  r^2/2 would be a fitted normal-ordering convention.")

    print(
        """
[4] Consolidated verdict
  I do not find a real Locked closure in the currently pinned canon objects.

  The remaining target is still:
      delta C = +0.001187526879799
      C2_required = 0.162691182532  if the object is order-alpha0 beyond C_loop.

  The three obvious connected-object routes fail for different reasons:
    * gauge-web two-loop: graph counts exist, but the actual two-loop vertex and
      Green function are not specified; the QED two-loop near-number is the
      wrong observable/Green's function;
    * non-Gaussian bath: the fixed second cumulant is too small, and higher
      cumulants are cutoff-sensitive rather than continuum constants;
    * QEC/Kraus normal-ordering: exact Lindblad gives zero while first-order
      Kraus conventions give r^2/4 or r^2/2; the target lies between them with
      no canon-selected normal-ordering map.

  Therefore a future positive theorem must first derive a new operator/tensor,
  not just a scalar prefactor.  Until then the non-horizon CC route remains
  rho_Lambda = 1.001872 rho_obs near-closure, not a Locked prediction.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
