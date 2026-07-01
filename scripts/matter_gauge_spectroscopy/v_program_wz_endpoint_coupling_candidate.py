#!/usr/bin/env python3
r"""v_program_wz_endpoint_coupling_candidate.py

Can the vector masses m_W and m_Z be tightened without importing arbitrary
fixed-threshold gauge couplings?

The absolute-mass audit showed that a control run using only

    sin^2(theta_W)(M_P) = 3/8,  alpha_0 = 1/137.036

misses W/Z by O(10%).  That does not by itself prove the vector sector is free:
W and Z are pole masses, so the physical formula uses the EW-scale electric
coupling and the on-shell weak-angle readout.

This script tests the sharpest non-fitting candidate:

    sin^2(theta_W)_pole = 2/9,
    alpha_EW = alpha(M_Z)  (now reduced in ew_alpha_mz_from_framework_dressed_alpha.py),
    v = v_phase7.

Important: this does NOT revive the retired M9 claim that 2/9 is a bare/UV
Weinberg angle.  The UV charge trace remains 3/8.  The candidate is instead an
IR/pole endpoint theorem: after EWSB and ordinary gauge running, the on-shell
service quotient read by the W/Z pole ledger is 2/9.  If the endpoint theorem is
not proven, the branch is only a standard-SM-input postdiction.

The audit also checks a tempting finite-count shortcut

    alpha_EW^{-1} ?= alpha_0^{-1} - N

and rejects promoting it: N=7 and N=8 are both close, and the required value is
non-integer in this proxy.  Thus the EM running/subtraction remains a separate
public-SM or dressed-alpha theorem, not an integer byte win.
"""

from __future__ import annotations

import math
import sys


ALPHA0_INV = 137.036
ALPHA0 = 1.0 / ALPHA0_INV
ALPHA_MZ_INV_PUBLIC = 128.95
ALPHA_MZ_PUBLIC = 1.0 / ALPHA_MZ_INV_PUBLIC
M_P = 1.2209e19
V_OBS = 246.22
V_RATIO_PHASE7 = 0.988719312
V_PHASE7 = V_RATIO_PHASE7 * V_OBS
M_W_OBS = 80.377
M_Z_OBS = 91.1876

ok = True


def check(name: str, cond: bool) -> None:
    global ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    ok = ok and bool(cond)


def relerr(pred: float, obs: float) -> float:
    return pred / obs - 1.0


def pct(pred: float, obs: float) -> float:
    return 100.0 * relerr(pred, obs)


def masses_from_alpha_s2(alpha: float, s2: float, v: float = V_PHASE7) -> tuple[float, float, float, float]:
    e = math.sqrt(4.0 * math.pi * alpha)
    s = math.sqrt(s2)
    c = math.sqrt(1.0 - s2)
    g2 = e / s
    gy = e / c
    gz = e / (s * c)
    return 0.5 * g2 * v, 0.5 * gz * v, g2, gy


def planck_boundary_gauge_run(mu: float, alpha: float = ALPHA0) -> tuple[float, float, float]:
    """One-loop toy: sin^2(theta_W)=3/8 at M_P plus low-energy alpha_0.

    Kept as the control from v_program_absolute_mass_closure_audit.py.
    """
    e2 = 4.0 * math.pi * alpha
    b_y = 41.0 / 6.0
    b_2 = -19.0 / 6.0
    log_mp_mu = math.log(M_P / mu)
    a_y = 2.0 * b_y * log_mp_mu / (16.0 * math.pi * math.pi)
    a_2 = 2.0 * b_2 * log_mp_mu / (16.0 * math.pi * math.pi)
    inv_g2_high = (1.0 / e2 - a_y - a_2) / (1.0 + 5.0 / 3.0)
    g2_high_sq = 1.0 / inv_g2_high
    g2_low = 1.0 / math.sqrt(1.0 / g2_high_sq + a_2)
    gy_low = 1.0 / math.sqrt(5.0 / (3.0 * g2_high_sq) + a_y)
    sin2 = gy_low * gy_low / (g2_low * g2_low + gy_low * gy_low)
    return gy_low, g2_low, sin2


def alpha_inv_required(s2: float, v: float, mass: float, kind: str) -> float:
    s = math.sqrt(s2)
    c = math.sqrt(1.0 - s2)
    if kind == "W":
        e = 2.0 * mass * s / v
    elif kind == "Z":
        e = 2.0 * mass * s * c / v
    else:
        raise ValueError(kind)
    return 4.0 * math.pi / (e * e)


def main() -> int:
    print("=" * 96)
    print("W/Z ENDPOINT COUPLING CANDIDATE")
    print("=" * 96)
    print(f"Phase-7 v = {V_PHASE7:.6f} GeV")

    print("\n[1] Control: UV 3/8 boundary plus alpha0 still fails")
    gy_run, g2_run, sin2_run = planck_boundary_gauge_run(V_OBS)
    mw_run = 0.5 * g2_run * V_OBS
    mz_run = 0.5 * math.sqrt(g2_run * g2_run + gy_run * gy_run) * V_OBS
    print(f"  one-loop run: gY={gy_run:.6f}, g2={g2_run:.6f}, sin^2={sin2_run:.6f}")
    print(f"  m_W={mw_run:.6f} GeV ({pct(mw_run, M_W_OBS):+.3f}%), "
          f"m_Z={mz_run:.6f} GeV ({pct(mz_run, M_Z_OBS):+.3f}%)")
    check("3/8+alpha0 Planck-boundary control misses W/Z by O(10%)",
          abs(relerr(mw_run, M_W_OBS)) > 0.05 and abs(relerr(mz_run, M_Z_OBS)) > 0.05)

    print("\n[2] Bare alpha with an IR 2/9 endpoint is not enough")
    s2_endpoint = 2.0 / 9.0
    mw_bare, mz_bare, g2_bare, gy_bare = masses_from_alpha_s2(ALPHA0, s2_endpoint)
    print(f"  sin^2_endpoint=2/9, alpha0^-1={ALPHA0_INV:.3f}")
    print(f"  m_W={mw_bare:.6f} GeV ({pct(mw_bare, M_W_OBS):+.3f}%), "
          f"m_Z={mz_bare:.6f} GeV ({pct(mz_bare, M_Z_OBS):+.3f}%)")
    check("bare alpha endpoint branch still misses at percent level",
          abs(relerr(mw_bare, M_W_OBS)) > 0.02 and abs(relerr(mz_bare, M_Z_OBS)) > 0.02)

    print("\n[3] EW-scale alpha plus on-shell 2/9 endpoint")
    mw_ep, mz_ep, g2_ep, gy_ep = masses_from_alpha_s2(ALPHA_MZ_PUBLIC, s2_endpoint)
    ratio_pred = math.sqrt(1.0 - s2_endpoint)
    ratio_obs = M_W_OBS / M_Z_OBS
    print(f"  alpha(M_Z)^-1={ALPHA_MZ_INV_PUBLIC:.3f}, sin^2_endpoint=2/9")
    print(f"  g2={g2_ep:.6f}, gY={gy_ep:.6f}")
    print(f"  m_W={mw_ep:.6f} GeV ({pct(mw_ep, M_W_OBS):+.3f}%), "
          f"m_Z={mz_ep:.6f} GeV ({pct(mz_ep, M_Z_OBS):+.3f}%)")
    print(f"  M_W/M_Z predicted={ratio_pred:.9f}, observed={ratio_obs:.9f}, "
          f"err={100.0 * (ratio_pred / ratio_obs - 1.0):+.3f}%")
    check("EW-alpha + 2/9 endpoint gives W/Z at sub-percent grade",
          abs(relerr(mw_ep, M_W_OBS)) < 0.005 and abs(relerr(mz_ep, M_Z_OBS)) < 0.005)
    check("2/9 endpoint gives the W/Z ratio at per-mille grade",
          abs(ratio_pred / ratio_obs - 1.0) < 0.001)

    print("\n[4] Required alpha_EW under the 2/9 endpoint")
    ainv_w = alpha_inv_required(s2_endpoint, V_PHASE7, M_W_OBS, "W")
    ainv_z = alpha_inv_required(s2_endpoint, V_PHASE7, M_Z_OBS, "Z")
    ainv_req = 0.5 * (ainv_w + ainv_z)
    print(f"  alpha^-1 required by W = {ainv_w:.6f}")
    print(f"  alpha^-1 required by Z = {ainv_z:.6f}")
    print(f"  mean required          = {ainv_req:.6f}")
    print(f"  public alpha(M_Z)^-1   = {ALPHA_MZ_INV_PUBLIC:.6f} "
          f"({100.0 * (ALPHA_MZ_INV_PUBLIC / ainv_req - 1.0):+.3f}% vs required)")
    check("required alpha_EW is close to the standard EW-scale electromagnetic coupling",
          abs(ALPHA_MZ_INV_PUBLIC / ainv_req - 1.0) < 0.01)

    print("\n[5] Discrete alpha0^-1 - N subtraction is not unique")
    rows = []
    print(f"  {'N':>3s} {'alpha_inv':>12s} {'mW err':>10s} {'mZ err':>10s}")
    for n in range(5, 11):
        alpha_inv = ALPHA0_INV - n
        mw, mz, _g2, _gy = masses_from_alpha_s2(1.0 / alpha_inv, s2_endpoint)
        score = max(abs(relerr(mw, M_W_OBS)), abs(relerr(mz, M_Z_OBS)))
        rows.append((score, n, alpha_inv, mw, mz))
        print(f"  {n:3d} {alpha_inv:12.6f} {pct(mw, M_W_OBS):+9.3f}% {pct(mz, M_Z_OBS):+9.3f}%")
    best = min(rows)
    print(f"  best integer subtraction in this scan: N={best[1]}, max error={100*best[0]:.3f}%")
    check("integer byte-subtraction candidates are not uniquely selected",
          sum(1 for score, *_ in rows if score < 0.01) >= 2)
    check("exact required subtraction is non-integer",
          abs((ALPHA0_INV - ainv_req) - round(ALPHA0_INV - ainv_req)) > 0.2)

    print("\n[6] Theorem status")
    print("  Candidate endpoint theorem:")
    print("      W/Z pole masses use alpha_EW = alpha(M_Z) and sin^2(theta_W)_on-shell = 2/9.")
    print("  This preserves M9: 2/9 is NOT a UV/bare Weinberg angle; the UV charge trace")
    print("  remains 3/8 with standard running.  What is new is the possible endpoint")
    print("  readout of the pole ledger after EWSB.")
    print("\n  Remaining work:")
    print("    * derive the post-EWSB on-shell 2/9 service quotient, or keep it retired;")
    print("    * use ew_alpha_mz_from_framework_dressed_alpha.py for alpha_EW;")
    print("      the remaining alpha-side caveat is first-principles hadronic VP,")
    print("      a strong-sector precision problem rather than an EW anchor;")
    print("    * repeat with fixed-scheme pole-mass matching alongside the Higgs quartic.")

    if ok:
        print("\nALL CHECKS PASSED")
        return 0
    print("\nCHECKS FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
