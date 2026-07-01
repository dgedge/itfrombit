#!/usr/bin/env python3
r"""v_phase5_precision_normalisation_target.py

Precision normalisation target for the EW VEV route.

This does not replace a proper multi-loop Standard Model matching calculation.
It isolates the one-loop target in the current conventions:

    alpha_0^8 / sqrt(lambda_crit(y_t)) = v/M_P.

Here lambda_crit(y_t) is the weak-scale quartic that runs to lambda(M_P)=0.
The script solves for the top-Yukawa proxy y_t* that would make the one-loop
critical-boundary prediction land exactly, and reports the equivalent finite
Coleman-Weinberg amplitude required at the central proxy.

Exit 0 means the residual is quantified, not closed.
"""

from __future__ import annotations

import math
import sys


LOOP = 16.0 * math.pi * math.pi
ALPHA0 = 1.0 / 137.036
M_P = 1.2209e19
M_Z = 91.1876
V_OBS = 246.22
M_H = 125.25
LAMBDA_PHYS = M_H * M_H / (2.0 * V_OBS * V_OBS)

G_Y0 = 0.3573
G_20 = 0.6517
G_30 = 1.2200
Y_T0 = 0.9500
N_STEPS = 6000
ok = True


def check(name: str, cond: bool) -> None:
    global ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    ok = ok and bool(cond)


def betas(g_y: float, g_2: float, g_3: float, y_t: float, lam: float) -> tuple[float, ...]:
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

    k1 = betas(*state)
    k2 = betas(*add(state, k1, 0.5 * dt))
    k3 = betas(*add(state, k2, 0.5 * dt))
    k4 = betas(*add(state, k3, dt))
    return tuple(si + dt * (a + 2.0 * b + 2.0 * c + d) / 6.0 for si, a, b, c, d in zip(state, k1, k2, k3, k4))


def run_to_planck(lambda_ew: float, y_t0: float) -> float:
    state = (G_Y0, G_20, G_30, y_t0, lambda_ew)
    dt = math.log(M_P / M_Z) / N_STEPS
    for _ in range(N_STEPS):
        state = rk4_step(state, dt)
    return state[4]


def critical_lambda(y_t0: float) -> float:
    lo, hi = 0.04, 0.28
    for _ in range(44):
        mid = 0.5 * (lo + hi)
        if run_to_planck(mid, y_t0) < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def v_ratio_from_lambda(lam: float) -> float:
    return (ALPHA0**8 / math.sqrt(lam)) / (V_OBS / M_P)


def solve_y_star() -> float:
    lo, hi = 0.90, 1.02
    for _ in range(44):
        mid = 0.5 * (lo + hi)
        ratio = v_ratio_from_lambda(critical_lambda(mid))
        # ratio decreases as y_t increases in this one-loop convention.
        if ratio > 1.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main() -> int:
    print("=" * 96)
    print("PHASE 5B: PRECISION RG / CW NORMALISATION TARGET")
    print("=" * 96)

    lambda_req = (ALPHA0**8 / (V_OBS / M_P)) ** 2
    lambda_central = critical_lambda(Y_T0)
    ratio_central = v_ratio_from_lambda(lambda_central)
    cw_amp_central = 1.0 / ratio_central
    y_star = solve_y_star()
    lambda_star = critical_lambda(y_star)
    ratio_star = v_ratio_from_lambda(lambda_star)
    dlog_dy = (math.log(v_ratio_from_lambda(critical_lambda(Y_T0 + 0.005)))
               - math.log(v_ratio_from_lambda(critical_lambda(Y_T0 - 0.005)))) / 0.010

    print("\n[0] Exact one-loop target in current conventions")
    print(f"  lambda required by exact alpha_0^8 VEV = {lambda_req:.9f}")
    print(f"  central proxy y_t0                     = {Y_T0:.6f}")
    print(f"  lambda_crit(y_t0)                      = {lambda_central:.9f}")
    print(f"  pred/obs at central proxy              = {ratio_central:.9f}")
    print(f"  finite CW amplitude needed centrally   = {cw_amp_central:.9f}")
    print(f"  delta ln(C_CW)                         = {math.log(cw_amp_central):+.9f}")

    print("\n[1] Top-threshold equivalent")
    print(f"  y_t* such that pred/obs=1              = {y_star:.9f}")
    print(f"  delta y_t* from central proxy          = {y_star - Y_T0:+.9f}")
    print(f"  lambda_crit(y_t*)                      = {lambda_star:.9f}")
    print(f"  pred/obs at y_t*                       = {ratio_star:.9f}")
    print(f"  d ln(pred/obs) / d y_t near central    = {dlog_dy:+.6f}")
    check("one-loop target y_t* exists in a plausible threshold range", 0.95 < y_star < 0.99)
    check("y_t* makes the one-loop critical-boundary VEV land", abs(ratio_star - 1.0) < 1.0e-5)
    check("central finite CW amplitude is modest rather than an OOM repair", 0.93 < cw_amp_central < 0.98)

    print("\n[2] Interpretation")
    print("  The same residual can be expressed as either:")
    print("    (a) a modest finite CW/field-normalisation amplitude C_CW ~= 0.952")
    print("        at the central one-loop proxy, or")
    print("    (b) a top-threshold/matching shift y_t0 -> y_t* ~= 0.969 in this")
    print("        one-loop convention.")
    print("  Because y_t0 and lambda_crit are scheme-dependent at this level, neither")
    print("  number is a canon constant.  They define the precision calculation that")
    print("  would be needed for a real closure.")

    print("\nVERDICT")
    print("  The normalisation residual is now a concrete precision target.  It is not")
    print("  evidence for a new free prefactor, but it is also not closed by the present")
    print("  one-loop calculation.  Closure requires either a proper multi-loop SM")
    print("  threshold run in a fixed scheme, or a substrate theorem for the finite")
    print("  Coleman-Weinberg normal-ordering constant.")

    if ok:
        print("\nALL CHECKS PASSED")
        return 0
    print("\nCHECKS FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
