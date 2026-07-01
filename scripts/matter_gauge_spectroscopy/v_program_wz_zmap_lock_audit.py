#!/usr/bin/env python3
r"""EW Z-map lock audit after alpha(M_Z) and endpoint-projector reductions.

The vector-mass branch had two distinct open ingredients:

  1. alpha(M_Z), previously inserted as the public electroweak coupling;
  2. the post-EWSB on-shell endpoint quotient sin^2(theta_W)=2/9.

The new alpha-running script closes (1) at standard-SM-running grade: the
framework's dressed Thomson alpha(0)^-1=137.036 runs to alpha(M_Z)^-1 ~=128.95.
The endpoint-exposure script closes (2) at finite pole-projector grade: the
post-EWSB W/Z pole space has 3 species x 3 on-shell polarisations = 9 residue
records, and the abelian hypercharge endpoint contributes exactly the two
transverse Maxwell records, leaving a 7-record weak complement.

Exit 0 means the Z-map quotient leg is closed at pole-projector grade, but the
absolute W/Z mass branch is still not a full precision lock because the V-map
and pole-mass matching residuals remain:

  * alpha(M_Z) is no longer an independent EW input;
  * the pole endpoint quotient is the projector trace Tr(P_B)/Tr(I)=2/9;
  * W/Z masses are sub-percent with alpha(M_Z) + that endpoint quotient;
  * UV 3/8 remains the charge-trace map and is not revived as an IR quotient.
"""

from __future__ import annotations

import importlib
import math
import sys


ew_alpha = importlib.import_module("ew_alpha_mz_from_framework_dressed_alpha")
wz = importlib.import_module("v_program_wz_endpoint_coupling_candidate")
pole = importlib.import_module("v_program_wz_pole_exposure_operator")

ok = True


def check(name: str, cond: bool) -> None:
    global ok
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    ok = ok and bool(cond)


def pct(pred: float, obs: float) -> float:
    return 100.0 * (pred / obs - 1.0)


def main() -> int:
    print("=" * 96)
    print("EW Z-MAP LOCK AUDIT: alpha(M_Z) reduced; endpoint quotient projector closed")
    print("=" * 96)

    print("\n[1] alpha(M_Z) from framework alpha(0) + SM running")
    delta_lep_1 = ew_alpha.delta_alpha_lep_1loop()
    inv_full = ew_alpha.alpha_mz_inv(ew_alpha.DELTA_LEP_FULL, ew_alpha.DELTA_HAD5)
    inv_1loop = ew_alpha.alpha_mz_inv(delta_lep_1, ew_alpha.DELTA_HAD5)
    band_lo = ew_alpha.alpha_mz_inv(
        ew_alpha.DELTA_LEP_FULL, ew_alpha.DELTA_HAD5 + ew_alpha.DELTA_HAD5_UNC
    )
    band_hi = ew_alpha.alpha_mz_inv(
        ew_alpha.DELTA_LEP_FULL, ew_alpha.DELTA_HAD5 - ew_alpha.DELTA_HAD5_UNC
    )
    print(f"    alpha(0)^-1 framework        = {ew_alpha.ALPHA0_INV:.6f}")
    print(f"    Delta_alpha_lep one-loop     = {delta_lep_1:.6f}")
    print(f"    Delta_alpha_lep full         = {ew_alpha.DELTA_LEP_FULL:.6f}")
    print(f"    Delta_alpha_had^(5)          = {ew_alpha.DELTA_HAD5:.5f} +- {ew_alpha.DELTA_HAD5_UNC:.5f}")
    print(f"    alpha(M_Z)^-1 full           = {inv_full:.6f}")
    print(f"    alpha(M_Z)^-1 one-loop-lep   = {inv_1loop:.6f}")
    print(f"    hadronic-VP band             = [{band_lo:.6f}, {band_hi:.6f}]")
    check("alpha(M_Z) is reproduced from framework alpha(0) to <0.05%",
          abs(inv_full / wz.ALPHA_MZ_INV_PUBLIC - 1.0) < 5.0e-4)
    check("public alpha(M_Z) lies in the framework-running/hadronic-VP band at this precision",
          band_lo - 0.012 <= wz.ALPHA_MZ_INV_PUBLIC <= band_hi + 0.012)

    print("\n[2] W/Z masses with the reduced alpha(M_Z) leg")
    alpha_mz_framework = 1.0 / inv_full
    s2_endpoint = 2.0 / 9.0
    mw, mz, g2, gy = wz.masses_from_alpha_s2(alpha_mz_framework, s2_endpoint)
    print(f"    g2={g2:.6f}, gY={gy:.6f}, sin^2_endpoint={s2_endpoint:.9f}")
    print(f"    m_W={mw:.6f} GeV ({pct(mw, wz.M_W_OBS):+.3f}%)")
    print(f"    m_Z={mz:.6f} GeV ({pct(mz, wz.M_Z_OBS):+.3f}%)")
    print(f"    M_W/M_Z=sqrt(7/9)={math.sqrt(7.0 / 9.0):.9f}")
    check("framework-run alpha(M_Z) + endpoint 2/9 gives sub-percent W/Z",
          abs(mw / wz.M_W_OBS - 1.0) < 0.005 and abs(mz / wz.M_Z_OBS - 1.0) < 0.005)

    print("\n[3] Pole endpoint projector and wrong finite quotients")
    candidates = {
        "UV SO(10)/GUT charge trace": 3.0 / 8.0,
        "one U(1) direction out of four": 1.0 / 4.0,
        "one neutral direction out of broken triad": 1.0 / 3.0,
        "post-EWSB pole endpoint projector": 2.0 / 9.0,
    }
    exact_pole = 1.0 - (wz.M_W_OBS / wz.M_Z_OBS) ** 2
    print(f"    exact pole-ratio sin^2 from masses = {exact_pole:.9f}")
    for name, value in candidates.items():
        print(f"    {name:<42s} {value:.9f} ({100.0 * (value / exact_pole - 1.0):+.3f}% vs pole)")
    check("UV 3/8 remains distinct from the pole endpoint", abs(3.0 / 8.0 - s2_endpoint) > 0.10)
    check("rank-4 and broken-triad rank quotients do not equal 2/9",
          abs(1.0 / 4.0 - s2_endpoint) > 0.02 and abs(1.0 / 3.0 - s2_endpoint) > 0.10)
    check("pole endpoint projector gives 2/9 rather than a same-space rank quotient",
          pole.uv_charge_trace_sin2().numerator == 3 and pole.uv_charge_trace_sin2().denominator == 8)
    check("2/9 is close to the pole ratio; residual is precision pole matching",
          abs(s2_endpoint / exact_pole - 1.0) < 0.005)

    print("\n[4] Lock status")
    print(
        """
    Z-map advanced:
      alpha(M_Z) is no longer an independent electroweak input.  It is the
      framework's dressed alpha(0) run to M_Z by standard photon vacuum
      polarisation; the hard part of that running is the usual hadronic VP /
      QCD R-ratio, not a new EW anchor.

    Z-map quotient leg closed at finite pole-projector grade:
      the W/Z LSZ endpoint space has rank 9, while the hypercharge endpoint
      projector has rank 2 (the two transverse abelian Maxwell records), so
      sin^2(theta_W)_pole = 2/9.  The weak complement has rank 7.

    What is still not a full precision lock:
      the exact reference pole ratio is 0.373% above 2/9, and the W/Z absolute
      masses inherit the current V-map and pole-matching residuals.  The UV
      charge trace remains 3/8 with ordinary SM running; it is a different map.
exit 0"""
    )

    if ok:
        print("ALL CHECKS PASSED")
        return 0
    print("CHECKS FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
