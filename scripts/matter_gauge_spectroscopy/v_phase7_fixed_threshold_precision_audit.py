#!/usr/bin/env python3
r"""v_phase7_fixed_threshold_precision_audit.py

Phase 7B: fixed-threshold precision audit for the EW VEV route.

This script is deliberately conservative.  It does *not* claim to implement
the full Buttazzo--Degrassi--Giardino--Giudice--Sala--Salvio--Strumia
NNLO/3-loop analysis.  Instead it makes the next honest step beyond the older
M_Z pole-proxy scripts:

  * use a fixed MSbar top-threshold convention at mu = M_t, with central
    NNLO-style boundary values from the near-criticality literature;
  * run the same critical-boundary calculation in that fixed convention;
  * apply the Phase-7A broken-triad normal-ordering scale kappa = 1/sqrt(3);
  * report how much residual remains.

The result is a guardrail: the broken-triad scale remains close, but the exact
landing moves at the percent level when the threshold convention changes.  So
the substrate scale is now derived locally, while the final numerical VEV still
requires a genuine fixed-scheme multi-loop SM threshold/RGE calculation.

Reference convention:
    arXiv:1307.3536 extracts SM MSbar parameters with full 2-loop NNLO
    threshold precision and extrapolates with full 3-loop NNLO RGEs.  The
    central values below are the standard Mt-scale values used as a lightweight
    fixed-threshold proxy here.
"""

from __future__ import annotations

import math
import sys


LOOP = 16.0 * math.pi * math.pi
ALPHA0 = 1.0 / 137.036
M_P = 1.2209e19
V_OBS = 246.22

# Fixed MSbar top-threshold proxy.  These are not fitted by this script.
M_T_MATCH = 173.34
G_Y_MT = 0.35830
G_2_MT = 0.64779
G_3_MT = 1.16660
Y_T_MT = 0.93690

N_STEPS = 9000
ok = True


def check(name: str, cond: bool) -> None:
    global ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    ok = ok and bool(cond)


def betas_1loop(g_y: float, g_2: float, g_3: float, y_t: float, lam: float) -> tuple[float, ...]:
    """One-loop SM beta functions in non-GUT hypercharge normalization.

    This script intentionally keeps the RGE order visible.  The threshold
    convention is fixed; the missing item is the full multi-loop RGE.
    """
    dg_y = (41.0 / 6.0) * g_y**3 / LOOP
    dg_2 = (-19.0 / 6.0) * g_2**3 / LOOP
    dg_3 = -7.0 * g_3**3 / LOOP
    dy_t = y_t * (
        4.5 * y_t**2
        - (17.0 / 12.0) * g_y**2
        - 2.25 * g_2**2
        - 8.0 * g_3**2
    ) / LOOP
    dlam = (
        24.0 * lam**2
        - 6.0 * y_t**4
        + 0.375 * (2.0 * g_2**4 + (g_2**2 + g_y**2) ** 2)
        + 12.0 * lam * y_t**2
        - 3.0 * lam * (3.0 * g_2**2 + g_y**2)
    ) / LOOP
    return dg_y, dg_2, dg_3, dy_t, dlam


def rk4_step(state: tuple[float, ...], dt: float) -> tuple[float, ...]:
    def add(s: tuple[float, ...], k: tuple[float, ...], scale: float) -> tuple[float, ...]:
        return tuple(si + scale * ki for si, ki in zip(s, k))

    k1 = betas_1loop(*state)
    k2 = betas_1loop(*add(state, k1, 0.5 * dt))
    k3 = betas_1loop(*add(state, k2, 0.5 * dt))
    k4 = betas_1loop(*add(state, k3, dt))
    return tuple(si + dt * (a + 2.0 * b + 2.0 * c + d) / 6.0 for si, a, b, c, d in zip(state, k1, k2, k3, k4))


def run_to_planck(lambda_mt: float) -> float:
    state = (G_Y_MT, G_2_MT, G_3_MT, Y_T_MT, lambda_mt)
    dt = math.log(M_P / M_T_MATCH) / N_STEPS
    for _ in range(N_STEPS):
        state = rk4_step(state, dt)
    return state[4]


def critical_lambda_mt() -> float:
    lo, hi = 0.04, 0.28
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if run_to_planck(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def lambda_required() -> float:
    return (ALPHA0**8 / (V_OBS / M_P)) ** 2


def species(lambda_mt: float) -> list[tuple[str, float, float, float]]:
    gz2 = G_2_MT * G_2_MT + G_Y_MT * G_Y_MT
    return [
        ("W", 6.0, G_2_MT * G_2_MT / 4.0, 5.0 / 6.0),
        ("Z", 3.0, gz2 / 4.0, 5.0 / 6.0),
        ("top", -12.0, Y_T_MT * Y_T_MT / 2.0, 3.0 / 2.0),
        ("H", 1.0, 3.0 * lambda_mt, 3.0 / 2.0),
        ("Goldstone", 3.0, lambda_mt, 3.0 / 2.0),
    ]


def delta_lambda_cw(kappa: float, lambda_mt: float) -> float:
    return sum(n * a * a * (math.log(a / (kappa * kappa)) - c) / LOOP for _, n, a, c in species(lambda_mt))


def vev_ratio(lambda_eff: float) -> float:
    return (ALPHA0**8 / math.sqrt(lambda_eff)) / (V_OBS / M_P)


def solve_kappa(lambda_mt: float, target_delta: float) -> float:
    lo, hi = 0.05, 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if delta_lambda_cw(mid, lambda_mt) < target_delta:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main() -> int:
    print("=" * 96)
    print("PHASE 7B: FIXED-THRESHOLD PRECISION AUDIT")
    print("=" * 96)

    lam_crit = critical_lambda_mt()
    lam_req = lambda_required()
    k_triad = 1.0 / math.sqrt(3.0)
    delta_triad = delta_lambda_cw(k_triad, lam_crit)
    lam_eff = lam_crit + delta_triad
    ratio_no_cw = vev_ratio(lam_crit)
    ratio_triad = vev_ratio(lam_eff)
    k_req = solve_kappa(lam_crit, lam_req - lam_crit)

    print("\n[0] Fixed top-threshold inputs")
    print(f"  mu_match = M_t = {M_T_MATCH:.2f} GeV")
    print(f"  gY={G_Y_MT:.5f}, g2={G_2_MT:.5f}, g3={G_3_MT:.5f}, y_t={Y_T_MT:.5f}")
    print("  RGE used here: one-loop only.  Threshold convention: fixed MSbar Mt proxy.")
    check("threshold inputs are in the expected MSbar ranges",
          0.35 < G_Y_MT < 0.37 and 0.63 < G_2_MT < 0.67 and 1.1 < G_3_MT < 1.25 and 0.90 < Y_T_MT < 0.97)

    print("\n[1] Critical-boundary run in the fixed threshold convention")
    print(f"  lambda_crit(M_t), one-loop RGE to lambda(M_P)=0 = {lam_crit:.9f}")
    print(f"  lambda_required by exact alpha_0^8 VEV       = {lam_req:.9f}")
    print(f"  pred/obs without finite CW                   = {ratio_no_cw:.9f}")
    check("fixed-threshold critical lambda is close to, but distinct from, the M_Z proxy value",
          0.145 < lam_crit < 0.152)

    print("\n[2] Broken-triad finite CW normal ordering")
    print(f"  kappa_triad = 1/sqrt(3)                       = {k_triad:.9f}")
    print(f"  Delta lambda_CW(kappa_triad)                 = {delta_triad:+.9f}")
    print(f"  lambda_eff                                  = {lam_eff:.9f}")
    print(f"  pred/obs with triad finite CW                = {ratio_triad:.9f}")
    print(f"  exact kappa in this mixed-order proxy         = {k_req:.9f}")
    print(f"  kappa_triad / kappa_exact - 1                = {100.0 * (k_triad / k_req - 1.0):+.3f}%")
    check("triad normal ordering remains percent-close in a fixed threshold convention",
          abs(ratio_triad - 1.0) < 0.02)
    check("fixed-threshold proxy is not an exact lock", abs(ratio_triad - 1.0) > 0.005)

    print("\nVERDICT")
    print("  The fixed-threshold audit confirms the right scale but prevents overclaim.")
    print("  Moving from the older M_Z pole proxy to an Mt MSbar threshold convention")
    print("  shifts the exact kappa target at the percent level.  Therefore Phase 7A")
    print("  derives the local substrate normal-ordering scale, but the final Higgs VEV")
    print("  number is not locked by this mixed-order calculation.  A real closure needs")
    print("  a full fixed-scheme run: two-loop threshold matching plus three-loop SM")
    print("  RGEs, or an equivalent public precision-SM implementation.")

    if ok:
        print("\nALL CHECKS PASSED")
        return 0
    print("\nCHECKS FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
