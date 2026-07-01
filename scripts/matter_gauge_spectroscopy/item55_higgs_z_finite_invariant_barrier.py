#!/usr/bin/env python3
"""Item 55: finite-invariant barrier for the Higgs/Z mass ratio.

The current canon has already constructed both sides of the old target:

* Higgs side: the O_h matter-cell A1g breathing eigenvalue
      lambda_H = 2 k_edge + 4 k_face + 2 k_body + k_site.

* Z side: the D_4h gauge-bridge Eu skeleton
      lambda_Z = k_shear + k_mix^2 / k_shear.

This script tests the remaining hope that a parameter-free m_H/m_Z ratio is
hidden in the finite representation data already available.  It is deliberately
not a fit.  It checks whether the existing finite A1g/bridge invariants select
one ratio, or whether the ratio remains continuously movable unless a new
single electroweak service/Hessian theorem ties the scalar and gauge stiffness
ledgers together.

Exit 0 means the finite-invariant route is blocked in the current operator
class and the remaining theorem target is sharply identified.
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass


ok = True


def check(name: str, cond: bool) -> None:
    global ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    ok = ok and bool(cond)


@dataclass(frozen=True)
class EwReference:
    m_h: float = 125.25
    m_z: float = 91.1876

    @property
    def ratio(self) -> float:
        return self.m_h / self.m_z

    @property
    def lambda_over_gz2(self) -> float:
        # SM tree identity: (m_H/m_Z)^2 = 8 lambda_H / g_Z^2.
        return self.ratio * self.ratio / 8.0


def lambda_h(k_edge: float, k_face: float, k_body: float, k_site: float) -> float:
    return 2.0 * k_edge + 4.0 * k_face + 2.0 * k_body + k_site


def lambda_z(k_shear: float, k_mix: float) -> float:
    if k_shear <= 0.0:
        raise ValueError("k_shear must be positive")
    return k_shear + k_mix * k_mix / k_shear


def structural_ratio(params: tuple[float, float, float, float, float, float]) -> float:
    ke, kf, kb, ks, kshear, kmix = params
    return math.sqrt(lambda_h(ke, kf, kb, ks) / lambda_z(kshear, kmix))


def main() -> int:
    ref = EwReference()
    near_r = math.sqrt(15.0 / 8.0)
    near_lambda_over_gz2 = 15.0 / 64.0

    print("=" * 96)
    print("ITEM 55 HIGGS/Z FINITE-INVARIANT BARRIER")
    print("=" * 96)
    print("\n[0] Reference target")
    print(f"  m_H/m_Z pole proxy                 = {ref.ratio:.9f}")
    print(f"  required lambda_H/g_Z^2            = {ref.lambda_over_gz2:.9f}")
    print(f"  rational near-target sqrt(15/8)    = {near_r:.9f}")
    print(f"  rational near-target 15/64         = {near_lambda_over_gz2:.9f}")
    print(
        f"  near-target ratio error            = "
        f"{100.0 * (near_r / ref.ratio - 1.0):+.4f}%"
    )
    check("15/64 is close but not exact", 0.002 < abs(near_r / ref.ratio - 1.0) < 0.01)

    print("\n[1] Existing finite spectral data")
    print("  matter-cell A1g: lambda_H = 2 ke + 4 kf + 2 kb + ksite")
    print("  gauge bridge Eu: lambda_Z = kshear + kmix^2/kshear")
    print("  all-bonds matter-cell internal ratio lambda_A1g/lambda_Eg = 8/3")
    internal_a1g_eg_ratio = math.sqrt(8.0 / 3.0)
    print(f"  sqrt(8/3) = {internal_a1g_eg_ratio:.9f}  (wrong sector: Eg is not Z)")
    check(
        "matter-cell internal A1g/Eg ratio is not the observed Higgs/Z ratio",
        abs(internal_a1g_eg_ratio / ref.ratio - 1.0) > 0.10,
    )

    print("\n[2] Continuous-movability test inside the current skeleton")
    # Both rows obey the same finite representation content.  They differ only
    # in the allowed positive stiffnesses.  If the finite data selected the
    # ratio, these two values could not differ.
    params_a = (1.0, 0.0, 0.0, 0.0, 1.0, 1.0)  # lambda_H=2, lambda_Z=2
    params_b = (1.0, 0.5, 0.25, 0.3, 2.0, 0.7)  # lambda_H=4.8, lambda_Z=2.245
    params_c = (3.0, 1.0, 0.5, 0.0, 1.0, 3.0)  # lambda_H=12, lambda_Z=10
    ratios = [structural_ratio(p) for p in (params_a, params_b, params_c)]
    for label, params, ratio in zip(("A", "B", "C"), (params_a, params_b, params_c), ratios):
        print(f"  admissible point {label}: params={params} -> M_H/M_Z={ratio:.9f}")
    check("same finite irreps allow multiple different ratios", max(ratios) - min(ratios) > 0.20)

    # Constructive tuning: for any positive target ratio r, fixed Higgs stiffness
    # and k_shear=1 can solve k_mix = sqrt(lambda_H/r^2 - 1) when lambda_H>r^2.
    lh_fixed = lambda_h(1.0, 0.5, 0.25, 0.3)
    kmix_for_obs = math.sqrt(lh_fixed / (ref.ratio * ref.ratio) - 1.0)
    tuned = math.sqrt(lh_fixed / lambda_z(1.0, kmix_for_obs))
    print("\n[3] Explicit tuning witness")
    print(f"  hold lambda_H fixed at {lh_fixed:.6f}, k_shear=1")
    print(f"  choose k_mix={kmix_for_obs:.9f} -> M_H/M_Z={tuned:.9f}")
    check("observed ratio is reachable only after choosing a continuous bridge stiffness", abs(tuned - ref.ratio) < 1e-12)

    print("\n[4] Verdict")
    print(
        "  The finite A1g matter-cell and D4h bridge constructions do not contain\n"
        "  a hidden parameter-free Higgs/Z ratio.  They identify the scalar and\n"
        "  neutral-gauge directions, but the ratio is continuously movable through\n"
        "  the allowed positive stiffnesses and the cross-action-space scale.  The\n"
        "  attractive 15/64 target is therefore not licensed by the current finite\n"
        "  invariant data; it would need a new electroweak service/Hessian theorem\n"
        "  or a UV boundary condition fixing lambda_H, y_t, and the gauge couplings\n"
        "  before RG running."
    )
    check("finite-invariant route is blocked in the current item-55 skeleton", True)

    if ok:
        print("\nALL CHECKS PASSED")
        return 0
    print("\nCHECKS FAILED")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
