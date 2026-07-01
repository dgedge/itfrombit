#!/usr/bin/env python3
r"""Consolidate the EW/nuclear second-anchor frontier.

The current v-program has several strong pieces:

  * complete-cell billing / filled-byte selection;
  * the broken EW service scale kappa = 1/sqrt(3);
  * a scalar-quartic candidate lambda(M_P) = -4 alpha0;
  * a vector endpoint candidate: framework-run alpha(M_Z) plus
    sin^2(theta_W)_pole = 2/9, now supplied by the pole-exposure projector.

This audit asks what those pieces do and do not close when combined.  The
important discipline is to keep the three boundary maps separate:

  V-map:  complete-cell billing + CW/RG threshold matching -> v,
  H-map:  scalar quartic/pole boundary -> m_H,
  Z-map:  EW pole coupling/endpoint readout -> m_W,m_Z.

Exit 0 means the frontier is tightened, not closed:

  * the old "irreducible EW/top second anchor" wording is too pessimistic;
  * the absolute EW masses are not locked by the current canon, because the
    V/H maps still need fixed-scheme precision and scalar boundary/pole
    theorems, while the Z-map quotient leg is closed only at tree/pole-projector
    grade and still needs precision pole matching;
  * nuclear MeV-scale coefficients remain QCD/contact many-body residuals,
    not a third absolute anchor.
"""

from __future__ import annotations

import importlib
import math
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "python_code"))

abs_mass = importlib.import_module("v_program_absolute_mass_closure_audit")
higgs = importlib.import_module("v_program_higgs_quartic_boundary_candidate")
wz = importlib.import_module("v_program_wz_endpoint_coupling_candidate")
ew_alpha = importlib.import_module("ew_alpha_mz_from_framework_dressed_alpha")
pole = importlib.import_module("v_program_wz_pole_exposure_operator")


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def pct(pred: float, obs: float) -> float:
    return 100.0 * (pred / obs - 1.0)


def main() -> None:
    print("EW / NUCLEAR SECOND-ANCHOR PRECISION FRONTIER AUDIT")
    print("=" * 96)

    print("[1] Owner-script anchors")
    anchors = [
        (
            "complete-cell scalar pointer",
            ROOT / "python_code/v_phase6_complete_cell_billing_scalar_pointer.py",
            "only vacuum and filled byte are O_h-fixed pointer codewords",
        ),
        (
            "broken-triad service scale",
            ROOT / "python_code/v_phase7_broken_triad_service_scale.py",
            "EW service Gram matrix has rank 3",
        ),
        (
            "absolute mass closure not locked",
            ROOT / "python_code/v_program_absolute_mass_closure_audit.py",
            "does not fully close to absolute m_H, m_W, m_Z",
        ),
        (
            "Higgs quartic boundary target",
            ROOT / "python_code/v_program_higgs_quartic_boundary_candidate.py",
            "lambda(M_P) = - N_EW * alpha_0",
        ),
        (
            "W/Z endpoint-coupling target",
            ROOT / "python_code/v_program_wz_endpoint_coupling_candidate.py",
            "sin^2(theta_W)_pole = 2/9",
        ),
        (
            "alpha(M_Z) from framework alpha(0)",
            ROOT / "python_code/ew_alpha_mz_from_framework_dressed_alpha.py",
            "alpha(M_Z) is the framework's alpha(0), not an EW import",
        ),
        (
            "Z-map lock audit",
            ROOT / "python_code/v_program_wz_zmap_lock_audit.py",
            "endpoint quotient projector closed",
        ),
        (
            "W/Z pole endpoint exposure operator",
            ROOT / "python_code/v_program_wz_pole_exposure_operator.py",
            "Tr E_B = 2,  Tr E_W = 7",
        ),
        (
            "nuclear is contact residual, not third anchor",
            ROOT / "python_code/item113_t1_t2_local_map_theorems.py",
            "closed-record-pair energy unit eps_sat=2 alpha0 Lambda",
        ),
    ]
    for label, path, needle in anchors:
        text = path.read_text()
        ok = needle in text
        print(f"    [{'PASS' if ok else 'FAIL'}] {label}")
        if not ok:
            raise AssertionError(label)

    print("\n[2] V-map: complete-cell billing plus current fixed-threshold proxy")
    v_pred = abs_mass.V_RATIO_TRIAD * abs_mass.V_OBS
    print(f"    v_pred = {v_pred:.6f} GeV ({pct(v_pred, abs_mass.V_OBS):+.3f}% vs v_obs)")
    check(abs(v_pred / abs_mass.V_OBS - 1.0) < 0.02, "current v-program lands at percent grade")
    check(abs(v_pred / abs_mass.V_OBS - 1.0) > 0.005, "current v-program is not yet a precision EW input")

    print("\n[3] H-map: quartic boundary and its coupling to the VEV")
    lambda_c4 = higgs.solve_lambda_for_boundary(4.0)
    lambda_c3 = higgs.solve_lambda_for_boundary(3.0)
    lambda_c5 = higgs.solve_lambda_for_boundary(5.0)
    mh_c4_vobs = higgs.higgs_mass(lambda_c4, higgs.V_OBS)
    mh_c4_vprog = higgs.higgs_mass(lambda_c4, v_pred)
    lambda_needed_vprog = higgs.M_H_OBS**2 / (2.0 * v_pred**2)
    c_needed_vprog = -higgs.run_to_planck(lambda_needed_vprog) / higgs.ALPHA0
    c_needed_vobs = -higgs.run_to_planck(
        higgs.M_H_OBS**2 / (2.0 * higgs.V_OBS**2)
    ) / higgs.ALPHA0
    print(f"    C=3 -> m_H(v_obs) = {higgs.higgs_mass(lambda_c3, higgs.V_OBS):.6f} GeV")
    print(f"    C=4 -> m_H(v_obs) = {mh_c4_vobs:.6f} GeV ({pct(mh_c4_vobs, higgs.M_H_OBS):+.3f}%)")
    print(f"    C=5 -> m_H(v_obs) = {higgs.higgs_mass(lambda_c5, higgs.V_OBS):.6f} GeV")
    print(f"    C=4 -> m_H(v_prog)= {mh_c4_vprog:.6f} GeV ({pct(mh_c4_vprog, higgs.M_H_OBS):+.3f}%)")
    print(f"    observed v reads C = {c_needed_vobs:.6f}")
    print(f"    v-program v reads C = {c_needed_vprog:.6f}")
    check(abs(mh_c4_vobs / higgs.M_H_OBS - 1.0) < 0.005, "C=4 is a sharp Higgs target if observed v is supplied")
    check(abs(mh_c4_vprog / higgs.M_H_OBS - 1.0) > 0.005, "C=4 does not by itself absorb the current VEV residual")
    check(3.4 < c_needed_vprog < 3.7, "with v-program v, exact Higgs would require non-integer effective C")

    print("\n[4] Z-map: vector pole endpoint and EW coupling")
    s2_endpoint = 2.0 / 9.0
    alpha_mz_inv_framework = ew_alpha.alpha_mz_inv(ew_alpha.DELTA_LEP_FULL, ew_alpha.DELTA_HAD5)
    alpha_mz_framework = 1.0 / alpha_mz_inv_framework
    mw_ep, mz_ep, g2_ep, gy_ep = wz.masses_from_alpha_s2(alpha_mz_framework, s2_endpoint)
    s2_exact_ratio = 1.0 - (wz.M_W_OBS / wz.M_Z_OBS) ** 2
    alpha_inv_req_w = wz.alpha_inv_required(s2_endpoint, v_pred, wz.M_W_OBS, "W")
    alpha_inv_req_z = wz.alpha_inv_required(s2_endpoint, v_pred, wz.M_Z_OBS, "Z")
    alpha_inv_req_mean = 0.5 * (alpha_inv_req_w + alpha_inv_req_z)
    alpha_inv_req_exact = wz.alpha_inv_required(s2_exact_ratio, v_pred, wz.M_W_OBS, "W")
    print(f"    framework alpha(M_Z)^-1 = {alpha_mz_inv_framework:.6f} from alpha(0)+SM running")
    print(f"    endpoint 2/9 + framework alpha(M_Z): m_W={mw_ep:.6f} ({pct(mw_ep, wz.M_W_OBS):+.3f}%), m_Z={mz_ep:.6f} ({pct(mz_ep, wz.M_Z_OBS):+.3f}%)")
    print(f"    exact pole-ratio sin^2 = {s2_exact_ratio:.9f}; 2/9 = {s2_endpoint:.9f} ({100*(s2_exact_ratio/s2_endpoint-1):+.3f}% in sin^2)")
    print(f"    alpha^-1 required at 2/9: W={alpha_inv_req_w:.6f}, Z={alpha_inv_req_z:.6f}, mean={alpha_inv_req_mean:.6f}")
    print(f"    alpha^-1 required at exact ratio = {alpha_inv_req_exact:.6f}; public alpha(M_Z)^-1={wz.ALPHA_MZ_INV_PUBLIC:.6f}")
    check(abs(alpha_mz_inv_framework / wz.ALPHA_MZ_INV_PUBLIC - 1.0) < 5e-4, "alpha(M_Z) is reduced to framework alpha(0) + SM running")
    check(pole.uv_charge_trace_sin2().numerator == 3 and pole.uv_charge_trace_sin2().denominator == 8,
          "UV trace remains 3/8 while pole projector supplies the endpoint quotient")
    check(abs(mw_ep / wz.M_W_OBS - 1.0) < 0.005 and abs(mz_ep / wz.M_Z_OBS - 1.0) < 0.005, "endpoint candidate gives sub-percent W/Z with framework-run alpha(M_Z)")
    check(abs(s2_exact_ratio / s2_endpoint - 1.0) < 0.005, "2/9 is close to the pole ratio")
    check(abs(wz.ALPHA_MZ_INV_PUBLIC / alpha_inv_req_mean - 1.0) > 0.003, "the remaining Z-map mismatch is precision pole matching, not alpha(M_Z) import")

    print("\n[5] Failure controls")
    gy_run, g2_run, _sin2_run = abs_mass.planck_boundary_gauge_run(abs_mass.V_OBS)
    mw_uv = 0.5 * g2_run * abs_mass.V_OBS
    mz_uv = 0.5 * math.sqrt(g2_run * g2_run + gy_run * gy_run) * abs_mass.V_OBS
    mw_bare, mz_bare, *_ = wz.masses_from_alpha_s2(wz.ALPHA0, s2_endpoint)
    print(f"    UV 3/8 + alpha0:       m_W={mw_uv:.3f} ({pct(mw_uv, wz.M_W_OBS):+.2f}%), m_Z={mz_uv:.3f} ({pct(mz_uv, wz.M_Z_OBS):+.2f}%)")
    print(f"    IR 2/9 + bare alpha0:  m_W={mw_bare:.3f} ({pct(mw_bare, wz.M_W_OBS):+.2f}%), m_Z={mz_bare:.3f} ({pct(mz_bare, wz.M_Z_OBS):+.2f}%)")
    check(abs(mw_uv / wz.M_W_OBS - 1.0) > 0.05, "UV 3/8 + alpha0 remains dead for W")
    check(abs(mw_bare / wz.M_W_OBS - 1.0) > 0.02, "IR endpoint still needs EW-scale alpha, not bare alpha0")

    print("\n[6] Consolidated status")
    print("    The EW second-anchor problem has split into three finite precision gates:")
    print("      V-map: complete-cell billing + kappa=1/sqrt(3) + fixed-scheme SM thresholds -> v;")
    print("      H-map: scalar record-action quartic boundary / pole matching -> m_H;")
    print("      Z-map: alpha(M_Z) reduced to dressed-alpha/SM running; pole-projector 2/9 -> m_W,m_Z.")
    print("    Nuclear MeV coefficients do not form a third absolute anchor: item 113")
    print("    reduces them to alpha0*Lambda contact energy plus local TCH contact maps,")
    print("    with many-body/shell physics still open.")

    print(
        """
================================================================================
VERDICT
  Tightened, not locked.  The old statement "EW/top second anchor" is too
  pessimistic because the VEV scale is now a conditional radiative/complete-cell
  prediction.  But absolute m_H,m_W,m_Z are not yet one-anchor outputs.

  What would lock the sector is now precise:
    (1) a public fixed-scheme multi-loop V-map calculation for v;
    (2) an operator proof of the scalar boundary lambda(M_P)=-4 alpha0, with
        the pole-mass matching done in the same scheme;
    (3) precision W/Z pole matching around the finite endpoint projector;
        alpha(M_Z) is reduced to dressed-alpha plus standard SM running, and
        the 2/9 endpoint quotient is supplied by the LSZ pole-exposure operator.

  Until those land, the EW sector is a bounded dimensionless-boundary problem,
  not an unconstrained dimensionful anchor and not a locked prediction.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
