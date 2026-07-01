#!/usr/bin/env python3
r"""v_phase4_critical_prefactor_audit.py

Audit the remaining RG/Coleman-Weinberg prefactor in the EW VEV route.

The earlier phase-3 formula used the measured weak-scale quartic,

    v/M_P = alpha_0^8 / sqrt(lambda_EW),

and landed about 11% high.  This script asks whether the residual is better
understood as a precision RG boundary problem rather than as a free O(1)
Coleman-Weinberg multiplier.

We use the same one-loop SM beta functions as v_phase3_quartic.py.  The critical
boundary lambda(M_P)=0 selects a weak-scale lambda_crit.  In this one-loop
scheme lambda_crit ~= 0.144, which reduces the VEV mismatch from 11.1% to 5.1%.
Moreover the result is strongly top-threshold sensitive: a +0.02 shift in the
EW top-Yukawa proxy brings the ratio to 0.997.  That is not a derivation, but it
shows the residual sits in the known precision RG/top-threshold slot rather than
requiring a new fitted prefactor.

Exit 0 means: sharpened, not closed.  A closure still needs a precise multi-loop
SM matching calculation, or a genuine substrate theorem fixing the finite
Coleman-Weinberg normalisation.
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

# Same one-loop proxy inputs as v_phase3_quartic.py.
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


def run_to_planck(lambda_ew: float, y_t0: float = Y_T0) -> tuple[float, ...]:
    state = (G_Y0, G_20, G_30, y_t0, lambda_ew)
    dt = math.log(M_P / M_Z) / N_STEPS
    for _ in range(N_STEPS):
        state = rk4_step(state, dt)
    return state


def critical_lambda(y_t0: float = Y_T0) -> float:
    lo, hi = 0.04, 0.28
    for _ in range(44):
        mid = 0.5 * (lo + hi)
        if run_to_planck(mid, y_t0)[4] < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def v_ratio_from_lambda(lam: float) -> float:
    return (ALPHA0**8 / math.sqrt(lam)) / (V_OBS / M_P)


def main() -> int:
    print("=" * 96)
    print("PHASE 4B: RG / COLEMAN-WEINBERG PREFACTOR AUDIT")
    print("=" * 96)

    lambda_req = (ALPHA0**8 / (V_OBS / M_P)) ** 2
    lambda_mp_phys = run_to_planck(LAMBDA_PHYS)[4]
    lambda_crit = critical_lambda()
    ratio_phys = v_ratio_from_lambda(LAMBDA_PHYS)
    ratio_crit = v_ratio_from_lambda(lambda_crit)

    print("\n[0] Baseline")
    print(f"  lambda_phys(EW)              = {LAMBDA_PHYS:.9f}")
    print(f"  lambda(M_P) from lambda_phys = {lambda_mp_phys:+.9f}")
    print(f"  lambda_required by exact alpha_0^8 VEV = {lambda_req:.9f}")
    print(f"  alpha_0^8/sqrt(lambda_phys) / observed = {ratio_phys:.6f}")

    print("\n[1] Critical-boundary value in the same one-loop scheme")
    print(f"  lambda_crit(EW) such that lambda(M_P)=0 = {lambda_crit:.9f}")
    print(f"  alpha_0^8/sqrt(lambda_crit) / observed = {ratio_crit:.6f}")
    print(f"  lambda_crit / lambda_required          = {lambda_crit / lambda_req:.6f}")
    check("critical boundary improves the VEV mismatch", abs(ratio_crit - 1.0) < abs(ratio_phys - 1.0))
    check("central one-loop critical boundary does not close the prefactor exactly",
          0.02 < abs(ratio_crit - 1.0) < 0.08)

    print("\n[2] Top-threshold sensitivity")
    print(f"  {'delta y_t':>9s} {'lambda_crit':>14s} {'pred/obs':>10s}")
    scan = []
    for delta in [-0.03, -0.02, -0.01, 0.0, 0.01, 0.02, 0.03]:
        lc = critical_lambda(Y_T0 + delta)
        rr = v_ratio_from_lambda(lc)
        scan.append((delta, lc, rr))
        print(f"  {delta:+9.3f} {lc:14.9f} {rr:10.6f}")
    best = min(scan, key=lambda x: abs(x[2] - 1.0))
    print(f"  closest sampled threshold: delta y_t={best[0]:+.3f}, pred/obs={best[2]:.6f}")
    check("a small top-threshold shift can absorb the remaining percent-level residual",
          abs(best[0]) <= 0.03 and abs(best[2] - 1.0) < 0.01)

    print("\nVERDICT")
    print("  The prefactor is sharpened, not derived.  The central one-loop critical")
    print("  quartic reduces the VEV error from 11.1% to about 5.1%, and ordinary")
    print("  top-threshold sensitivity can move the result through the observed value.")
    print("  Therefore the residual should be treated as a precision RG / matching")
    print("  problem, not as an arbitrary Coleman-Weinberg multiplier.  A real closure")
    print("  still needs a multi-loop threshold calculation or a substrate theorem for")
    print("  the finite CW normalisation.")

    if ok:
        print("\nALL CHECKS PASSED")
        return 0
    print("\nCHECKS FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
