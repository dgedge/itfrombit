#!/usr/bin/env python3
r"""v_program_higgs_quartic_boundary_candidate.py

Candidate theorem for the Higgs physical quartic / pole-mass residual.

The previous v-program audit separated two different objects:

    * the VEV/zero-momentum normal-ordering object, where the broken EW
      service quotient has rank 3 and kappa = 1/sqrt(3);
    * the physical Higgs pole quartic, which still needed a boundary or
      pole-matching theorem.

This script tests the sharpest theorem-shaped candidate:

    lambda(M_P) = - N_EW * alpha_0, with N_EW = 4.

Interpretation:

    At the UV boundary the scalar quartic sees the full monitored
    SU(2)_L x U(1)_Y electroweak service alphabet: W1,W2,W3,B, hence four
    service directions.  The negative sign is the normal-ordered
    fermion/record-action instability direction of the Higgs potential.  After
    EWSB, the VEV normal-ordering scale is different: Q = T3 + Y is the
    unbroken null direction, leaving the rank-3 broken quotient and
    kappa = 1/sqrt(3).  Thus the 4-count and the 3-count are not competing
    readings; they apply to different sides of the problem.

This is not a public multi-loop SM calculation.  It uses the same fixed
M_t-threshold one-loop proxy as v_phase7_fixed_threshold_precision_audit.py.
The point is to decide whether there is a plausible finite theorem target,
not to claim precision closure.

Exit 0 means:

    * C = 4 gives the observed Higgs pole mass at sub-percent grade in the
      fixed-threshold proxy;
    * adjacent integer counts C = 3 and C = 5 miss at percent grade;
    * using the v-program predicted VEV carries the result back to the
      existing ~1% VEV residual, so the candidate is a quartic theorem, not by
      itself a full absolute-mass lock.
"""

from __future__ import annotations

import math
import sys


LOOP = 16.0 * math.pi * math.pi
ALPHA0 = 1.0 / 137.036
M_P = 1.2209e19
V_OBS = 246.22
M_H_OBS = 125.25
V_RATIO_PHASE7 = 0.988719312

# Fixed MSbar top-threshold proxy, as in phase 7B.
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


def relerr(pred: float, obs: float) -> float:
    return pred / obs - 1.0


def percent(pred: float, obs: float) -> float:
    return 100.0 * relerr(pred, obs)


def betas_1loop(g_y: float, g_2: float, g_3: float, y_t: float, lam: float) -> tuple[float, ...]:
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


def solve_lambda_for_boundary(c_count: float) -> float:
    target = -c_count * ALPHA0
    lo, hi = 0.04, 0.28
    for _ in range(70):
        mid = 0.5 * (lo + hi)
        if run_to_planck(mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def higgs_mass(lambda_low: float, v: float) -> float:
    return math.sqrt(2.0 * lambda_low) * v


def main() -> int:
    print("=" * 96)
    print("HIGGS QUARTIC BOUNDARY CANDIDATE: lambda(M_P) = -4 alpha_0")
    print("=" * 96)

    lambda_obs = M_H_OBS * M_H_OBS / (2.0 * V_OBS * V_OBS)
    lambda_mp_obs = run_to_planck(lambda_obs)
    c_obs = -lambda_mp_obs / ALPHA0

    print("\n[0] Fixed-threshold proxy readout of the observed Higgs")
    print(f"  lambda_obs(low)                    = {lambda_obs:.9f}")
    print(f"  lambda(M_P) from observed Higgs     = {lambda_mp_obs:+.9f}")
    print(f"  -lambda(M_P)/alpha_0                = {c_obs:.6f}")
    check("observed Higgs corresponds to an integer-near EW service count", abs(c_obs - 4.0) < 0.06)

    print("\n[1] Integer service-count boundary scan")
    print(f"  {'C':>3s} {'lambda_low':>14s} {'mH(v_obs)':>12s} {'err':>10s} {'mH(v_prog)':>12s} {'err':>10s}")
    rows = []
    v_prog = V_RATIO_PHASE7 * V_OBS
    for c_count in (3.0, 4.0, 5.0):
        lam = solve_lambda_for_boundary(c_count)
        mh_obs_v = higgs_mass(lam, V_OBS)
        mh_prog_v = higgs_mass(lam, v_prog)
        rows.append((c_count, lam, mh_obs_v, mh_prog_v))
        print(
            f"  {c_count:3.0f} {lam:14.9f} {mh_obs_v:12.6f} "
            f"{percent(mh_obs_v, M_H_OBS):+9.3f}% {mh_prog_v:12.6f} "
            f"{percent(mh_prog_v, M_H_OBS):+9.3f}%"
        )

    c4 = rows[1]
    c3 = rows[0]
    c5 = rows[2]
    check("C=4 gives sub-percent Higgs mass with observed v in this proxy",
          abs(relerr(c4[2], M_H_OBS)) < 0.005)
    check("adjacent integer C=3 is distinguishably worse", abs(relerr(c3[2], M_H_OBS)) > 0.015)
    check("adjacent integer C=5 is distinguishably worse", abs(relerr(c5[2], M_H_OBS)) > 0.015)
    check("using v-program v leaves the known percent-scale VEV residual",
          0.005 < abs(relerr(c4[3], M_H_OBS)) < 0.02)

    print("\n[2] Interpretation")
    print("  Candidate theorem:")
    print("      lambda(M_P) = -4 alpha_0")
    print("  where 4 is the full pre-EWSB EW service alphabet W1,W2,W3,B.")
    print("  This is distinct from the Phase-7 kappa=1/sqrt(3), where the unbroken")
    print("  electromagnetic direction Q=T3+Y has been quotiented out, leaving three")
    print("  broken directions for the VEV normal-ordering scale.")
    print("\n  Status:")
    print("    * promising theorem target: integer-near, sector-native, and adjacent")
    print("      integer controls are worse in the fixed-threshold proxy;")
    print("    * not a closure: the result must be repeated in a public fixed-scheme")
    print("      multi-loop threshold/RGE calculation, and the negative boundary sign")
    print("      still needs an operator proof from the record-action scalar sector.")

    if ok:
        print("\nALL CHECKS PASSED")
        return 0
    print("\nCHECKS FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
