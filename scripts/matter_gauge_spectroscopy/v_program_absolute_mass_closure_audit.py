#!/usr/bin/env python3
r"""v_program_absolute_mass_closure_audit.py

Can the v-program be closed to absolute m_H, m_W and m_Z?

This script tests the strongest honest reading of the current canon:

1.  Complete-cell billing is finite-theorem closed under the scalar protected
    pointer identification (Phase 6A).

2.  The broken-triad normal-ordering scale kappa=1/sqrt(3) is locally derived
    from the EW service algebra (Phase 7A).

3.  A fixed Mt-threshold one-loop proxy (Phase 7B) gives

        v_pred / v_obs = 0.9887,

    so a proper public multi-loop SM-RGE calculation plausibly has only a
    percent-level job left for the VEV.

Question: does that already imply absolute m_H, m_W, m_Z to few-percent grade?

Verdict
-------
Not as a full intrinsic closure.

* If the standard fixed-threshold low-energy gauge couplings are accepted as
  the dimensionless SM readout, then m_W and m_Z inherit the v-program error and
  land at ~1--2%.  This is a conditional postdiction, not a new gauge-coupling
  derivation.

* m_H is still the quartic/pole-matching bottleneck.  Using the one-loop
  critical quartic gives a ~6% high Higgs mass; using the finite-CW lambda_eff
  gives the identity m_H=sqrt(2) alpha_0^8 M_P and is ~11% high; using the
  observed quartic makes m_H few-percent only tautologically.

* A Planck-boundary 3/8 + alpha_0 one-loop gauge run does not close W/Z either:
  it gives m_W,m_Z high by ~8--12%.  Thus W/Z need the standard EW running/
  threshold input (or a new gauge-coupling theorem), not just the existing
  v-program.

The v-program is therefore narrowed but not locked: complete-cell billing is
closed, the local normal-ordering scale is derived, and v/W/Z are few-percent
conditional on fixed-scheme SM dimensionless couplings.  The remaining true
closure target is a public fixed-scheme multi-loop computation plus an EW
coupling/quartic boundary theorem, not only a finite cell theorem.
"""

from __future__ import annotations

import math


ALPHA0 = 1.0 / 137.036
M_P = 1.2209e19
V_OBS = 246.22
M_H_OBS = 125.25
M_W_OBS = 80.377
M_Z_OBS = 91.1876

# Phase-7B fixed Mt-threshold proxy values.
G_Y_MT = 0.35830
G_2_MT = 0.64779
LAMBDA_CRIT_MT = 0.148641944
DELTA_LAMBDA_TRIAD = 0.013992574
LAMBDA_EFF_MT = LAMBDA_CRIT_MT + DELTA_LAMBDA_TRIAD
V_RATIO_TRIAD = 0.988719312


def check(name: str, cond: bool) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        raise AssertionError(name)


def relerr(pred: float, obs: float) -> float:
    return pred / obs - 1.0


def percent(pred: float, obs: float) -> float:
    return 100.0 * relerr(pred, obs)


def planck_boundary_gauge_run(mu: float, alpha: float = ALPHA0) -> tuple[float, float, float]:
    """One-loop toy run from sin^2(theta_W)=3/8 at M_P plus low-energy alpha.

    This is a stringent control: it asks whether the existing charge-forced
    GQW boundary and alpha_0 are enough to derive g_2,g_Y at the EW scale.
    They are not; the output misses W/Z by O(10%).
    """
    e2 = 4.0 * math.pi * alpha
    b_y = 41.0 / 6.0
    b_2 = -19.0 / 6.0
    log_mp_mu = math.log(M_P / mu)
    a_y = 2.0 * b_y * log_mp_mu / (16.0 * math.pi * math.pi)
    a_2 = 2.0 * b_2 * log_mp_mu / (16.0 * math.pi * math.pi)
    # At the high boundary sin^2=3/8, so gY^2=(3/5)g2^2.
    # At the low scale, 1/e^2=1/gY^2+1/g2^2.
    inv_g2_high = (1.0 / e2 - a_y - a_2) / (1.0 + 5.0 / 3.0)
    g2_high_sq = 1.0 / inv_g2_high
    g2_low = 1.0 / math.sqrt(1.0 / g2_high_sq + a_2)
    gy_low = 1.0 / math.sqrt(5.0 / (3.0 * g2_high_sq) + a_y)
    sin2 = gy_low * gy_low / (g2_low * g2_low + gy_low * gy_low)
    return gy_low, g2_low, sin2


def main() -> int:
    print("V-PROGRAM ABSOLUTE MASS CLOSURE AUDIT")
    print("=" * 86)

    print("\n[1] Phase-7 v readout")
    v_pred = V_RATIO_TRIAD * V_OBS
    print(f"    v_pred = {v_pred:.6f} GeV  ({percent(v_pred, V_OBS):+.3f}% vs obs)")
    check("Phase-7 fixed-threshold proxy gives v at percent grade", abs(relerr(v_pred, V_OBS)) < 0.02)

    print("\n[2] W/Z if standard fixed-threshold gauge couplings are accepted")
    m_w = 0.5 * G_2_MT * v_pred
    m_z = 0.5 * math.sqrt(G_2_MT * G_2_MT + G_Y_MT * G_Y_MT) * v_pred
    print(f"    m_W = {m_w:.6f} GeV  ({percent(m_w, M_W_OBS):+.3f}% vs obs)")
    print(f"    m_Z = {m_z:.6f} GeV  ({percent(m_z, M_Z_OBS):+.3f}% vs obs)")
    check("conditional W mass is few-percent once g2 is supplied", abs(relerr(m_w, M_W_OBS)) < 0.03)
    check("conditional Z mass is few-percent once gY,g2 are supplied", abs(relerr(m_z, M_Z_OBS)) < 0.03)

    print("\n[3] Higgs mass alternatives")
    lambda_obs = M_H_OBS * M_H_OBS / (2.0 * V_OBS * V_OBS)
    m_h_from_obs_lambda = math.sqrt(2.0 * lambda_obs) * v_pred
    m_h_from_crit = math.sqrt(2.0 * LAMBDA_CRIT_MT) * v_pred
    m_h_from_eff = math.sqrt(2.0 * LAMBDA_EFF_MT) * v_pred
    m_h_byte_identity = math.sqrt(2.0) * M_P * ALPHA0**8
    print(f"    using observed lambda(v):       {m_h_from_obs_lambda:.6f} GeV  ({percent(m_h_from_obs_lambda, M_H_OBS):+.3f}%)")
    print(f"    using lambda_crit(M_t):         {m_h_from_crit:.6f} GeV  ({percent(m_h_from_crit, M_H_OBS):+.3f}%)")
    print(f"    using lambda_eff=crit+CW:       {m_h_from_eff:.6f} GeV  ({percent(m_h_from_eff, M_H_OBS):+.3f}%)")
    print(f"    identity sqrt(2) alpha_0^8 M_P: {m_h_byte_identity:.6f} GeV  ({percent(m_h_byte_identity, M_H_OBS):+.3f}%)")
    check("Higgs is few-percent only if the observed low-energy quartic is supplied",
          abs(relerr(m_h_from_obs_lambda, M_H_OBS)) < 0.02)
    check("one-loop critical quartic does not yet close m_H to few-percent",
          abs(relerr(m_h_from_crit, M_H_OBS)) > 0.03)
    check("finite-CW lambda_eff is not the physical Higgs quartic",
          abs(relerr(m_h_from_eff, M_H_OBS)) > abs(relerr(m_h_from_crit, M_H_OBS)))

    print("\n[4] Control: Planck-boundary 3/8 plus alpha_0 does not derive W/Z")
    gy_run, g2_run, sin2_run = planck_boundary_gauge_run(V_OBS)
    mw_run = 0.5 * g2_run * V_OBS
    mz_run = 0.5 * math.sqrt(g2_run * g2_run + gy_run * gy_run) * V_OBS
    print(f"    gY={gy_run:.6f}, g2={g2_run:.6f}, sin^2={sin2_run:.6f}")
    print(f"    m_W={mw_run:.6f} GeV ({percent(mw_run, M_W_OBS):+.3f}%), "
          f"m_Z={mz_run:.6f} GeV ({percent(mz_run, M_Z_OBS):+.3f}%)")
    check("simple Planck-boundary gauge run misses W/Z at O(10%)",
          abs(relerr(mw_run, M_W_OBS)) > 0.05 and abs(relerr(mz_run, M_Z_OBS)) > 0.05)

    print(
        """
[5] VERDICT
    The v-program does not fully close to absolute m_H, m_W, m_Z from the
    existing canon alone.

    What is closed/tight:
      * complete-cell billing is finite-theorem closed under the scalar pointer
        endpoint identification;
      * kappa=1/sqrt(3) is derived from the broken EW service rank;
      * the fixed-threshold proxy gives v at 1.1%;
      * W/Z are few-percent conditional on accepting standard fixed-threshold
        low-energy gauge couplings.

    What remains genuinely open:
      * m_H needs the physical quartic/pole matching.  The one-loop critical
        quartic is still ~6% high, while lambda_eff is a VEV-normalisation
        object, not the physical Higgs quartic.
      * W/Z need a real EW gauge-coupling boundary or the standard public SM
        running/matching inputs.  The existing 3/8 + alpha_0 Planck-boundary
        toy misses by O(10%).

    Therefore the correct status is: v-program substantially strengthens the
    second-anchor story, but a true "one-anchor absolute EW masses" lock still
    requires a fixed-scheme public multi-loop SM calculation and/or a new
    dimensionless EW coupling/quartic theorem.  The remaining gap is bounded and
    technical, but not gone.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
