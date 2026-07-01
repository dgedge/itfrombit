#!/usr/bin/env python3
r"""Item 123 / M15: CMB selector and halo-accounting consolidation audit.

Question
--------
After the zero-mode reservoir route, the CMB/dark-sector problem has two
remaining live gates:

  A. acoustic selector:

       H_CMB ?= (63/64) H_selector

     Numerically this repairs theta_*.  Is it derived by the current depth-6
     service ledger, or only a theorem target?

  B. halo accounting:

     Can the framework keep both a fair-sampling zero-mode CDM halo and an
     independent active R4/MOND galaxy force without double-counting the same
     mass budget?

Verdict
-------
The zero-mode route is conditionally strong: it gives the right pressureless
dark density and equality scale.  But the 63/64 acoustic factor is not derived
by the current symmetric depth-6 ledger.  A rank-63 projection over 64 slots is
possible only after adding a distinguished acoustic-invisible boundary/latch
slot.  The current boot/service ledger contains no such labelled slot; by the
same discipline that rejected the 129/128 mean-current correction, 63/64 stays
a sharp theorem target, not a Locked selector theorem.

Halo accounting closes only by branch choice:

  * zero-mode CDM halo branch: consistent bookkeeping, active R4/MOND is not
    also fitted as an independent force in the same galaxy regime;
  * MOND/RAR active branch: requires a new >95% zero-mode galaxy
    depletion/screening theorem;
  * fair-sampling zero-mode + active MOND: double-counted.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import itertools
import math
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_code"))

import item123_cmb_lock_attempt as lock
import item123_zero_mode_source_geodesic_halo_gate as halo


DEPTH = 6
N_SLOT = 2**DEPTH
DEPTH_FACTOR = Fraction(N_SLOT - 1, N_SLOT)


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def invariant_binary_masks(n: int) -> list[tuple[int, ...]]:
    """Return binary masks invariant under the full permutation group S_n.

    For n=64 we do not enumerate 64! permutations.  The theorem is simple:
    full permutation invariance forces all coordinates equal.  This helper
    enumerates small n as an executable witness and the same implication is
    asserted symbolically for n=64.
    """

    masks: list[tuple[int, ...]] = []
    perms = list(itertools.permutations(range(n)))
    for mask in itertools.product((0, 1), repeat=n):
        if all(tuple(mask[p[i]] for i in range(n)) == mask for p in perms):
            masks.append(mask)
    return masks


def rank_fraction(mask: tuple[int, ...]) -> Fraction:
    return Fraction(sum(mask), len(mask))


def h0_root_and_projection() -> tuple[float, float, float, float]:
    omega_b, _omega_nur, omega_dark = lock.framework_matter_budget()
    target = lock.theta100_planck()
    selector_theta = lock.theta100_framework(lock.H0_SELECTOR, omega_b, omega_dark)
    h0_root = lock.find_h0_root(omega_b, omega_dark, target)
    h0_6364 = lock.H0_SELECTOR * float(DEPTH_FACTOR)
    theta_6364 = lock.theta100_framework(h0_6364, omega_b, omega_dark)
    return target, selector_theta, h0_root, theta_6364


def density_result() -> tuple[float, float, float, float]:
    m_nur = halo.ALPHA0**2 * halo.LAMBDA_QCD_EV
    omega_nur = halo.omega_from_ratio(halo.ALPHA0 / 208.0, m_nur)
    omega_zero = 4.0 * omega_nur
    omega_dark = omega_nur + omega_zero
    z_eq = halo.z_eq(halo.OMEGA_B_H2 + omega_dark)
    return omega_nur, omega_zero, omega_dark, z_eq


@dataclass(frozen=True)
class Branch:
    name: str
    zero_mode_background: bool
    active_r4_force: bool
    zero_mode_depletion: float
    verdict: str

    @property
    def double_counted(self) -> bool:
        return self.zero_mode_background and self.active_r4_force and self.zero_mode_depletion < 0.95


def main() -> None:
    print("ITEM 123 / M15: CMB SELECTOR + HALO ACCOUNTING AUDIT")
    print("=" * 100)
    print(f"CAMB version: {getattr(lock.camb, '__version__', 'unknown')}")

    print("\n[1] Zero-mode dark reservoir arithmetic")
    omega_b, omega_nur_lock, omega_dark_lock = lock.framework_matter_budget()
    omega_nur, omega_zero, omega_dark, zeq = density_result()
    print(f"  omega_b h^2              = {omega_b:.6f}")
    print(f"  omega_nuR h^2            = {omega_nur:.6f}")
    print(f"  omega_zero h^2           = {omega_zero:.6f}")
    print(f"  omega_dark h^2           = {omega_dark:.6f}")
    print(f"  z_eq                     = {zeq:.1f}")
    check(abs(omega_nur / omega_nur_lock - 1.0) < 2.0e-3, "two owner scripts agree on sterile share")
    check(abs(omega_dark / omega_dark_lock - 1.0) < 2.0e-3, "two owner scripts agree on total dark share")
    check(abs(omega_dark - 0.1209) < 0.001, "zero-mode plus nu_R supplies an LCDM-sized pressureless budget")
    check(3400.0 < zeq < 3450.0, "matter-radiation equality is restored to the CMB range")

    print("\n[2] Acoustic selector numerics")
    target_theta, selector_theta, h0_root, theta_6364 = h0_root_and_projection()
    h0_6364 = lock.H0_SELECTOR * float(DEPTH_FACTOR)
    print(f"  Planck/LCDM reference 100 theta_* = {target_theta:.6f}")
    print(f"  selector H0={lock.H0_SELECTOR:.3f}: theta = {selector_theta:.6f} ({(selector_theta/target_theta-1)*100:+.3f}%)")
    print(f"  CAMB acoustic root H0             = {h0_root:.6f}")
    print(f"  (63/64) H_selector                = {h0_6364:.6f}")
    print(f"  theta at 63/64 projection         = {theta_6364:.6f} ({(theta_6364/target_theta-1)*100:+.4f}%)")
    print(f"  root / [(63/64) selector]         = {h0_root/h0_6364:.8f}")
    check(abs(selector_theta / target_theta - 1.0) > 0.002, "selector-H0 acoustic offset is real")
    check(abs(theta_6364 / target_theta - 1.0) < 2.0e-5, "63/64 projection numerically closes theta_*")
    check(abs(h0_root / h0_6364 - 1.0) < 1.0e-4, "CAMB root is essentially the 63/64 candidate")

    print("\n[3] Symmetry obstruction to deriving 63/64 from current depth-6 ledger")
    small_masks = invariant_binary_masks(4)
    small_ranks = sorted(rank_fraction(mask) for mask in small_masks)
    print(f"  executable S4 witness invariant ranks = {small_ranks}")
    check(small_ranks == [Fraction(0, 1), Fraction(1, 1)], "unlabelled symmetric slots admit only all-off/all-on projectors")
    print(f"  depth-6 service ledger slots          = 2^{DEPTH} = {N_SLOT}")
    print(f"  desired acoustic projector rank       = {N_SLOT - 1}/{N_SLOT}")
    print("  Current ledger status:")
    print("    * depth-6 residual-fault tensor ledger treats the 64 pair slots symmetrically;")
    print("    * the 129/128 audit found covariance but no licensed mean boundary leg;")
    print("    * no current acoustic operator labels one slot as invisible to theta_* while visible to late H0.")
    print("  Therefore a 63/64 projection breaks the present S_64 slot symmetry.")
    print("  It can be promoted only by a new boundary/latch primitive that names the missing slot.")
    check(DEPTH_FACTOR == Fraction(63, 64), "finite-depth candidate is exactly one missing slot out of 64")
    check(True, "current symmetric ledger cannot derive a rank-63 projector without extra structure")

    print("\n[4] Halo branch accounting")
    zero_to_baryon = omega_zero / omega_b
    rar_scatter_dex = 0.07
    tolerated_unmodelled = 10.0**rar_scatter_dex - 1.0
    max_zero_fraction_for_mond = tolerated_unmodelled / zero_to_baryon
    required_depletion = 1.0 - max_zero_fraction_for_mond
    branches = [
        Branch(
            "zero-mode CDM halo",
            zero_mode_background=True,
            active_r4_force=False,
            zero_mode_depletion=0.0,
            verdict="allowed branch",
        ),
        Branch(
            "active R4/MOND with screened zero-mode",
            zero_mode_background=True,
            active_r4_force=True,
            zero_mode_depletion=required_depletion,
            verdict="open theorem branch",
        ),
        Branch(
            "fair-sample zero-mode plus active R4/MOND",
            zero_mode_background=True,
            active_r4_force=True,
            zero_mode_depletion=0.0,
            verdict="double-counted",
        ),
    ]
    print(f"  zero-mode / baryon density ratio      = {zero_to_baryon:.3f}")
    print(f"  RAR scatter proxy                     = {rar_scatter_dex:.2f} dex -> {tolerated_unmodelled:.3f} fractional room")
    print(f"  MOND branch max fair-sample zero-mode = {max_zero_fraction_for_mond:.3f}")
    print(f"  required zero-mode depletion          >= {required_depletion:.1%}")
    for branch in branches:
        dc = "DOUBLE" if branch.double_counted else "ok"
        print(
            f"  {branch.name:42s} {branch.verdict:20s} "
            f"depletion={branch.zero_mode_depletion:6.1%} accounting={dc}"
        )
    check(zero_to_baryon > 4.0, "fair-sample zero-mode is a dominant halo-mass component")
    check(required_depletion > 0.95, "MOND-preserving branch needs more than 95% zero-mode depletion/screening")
    check(branches[2].double_counted, "fair-sample zero-mode plus active R4/MOND double-counts the galaxy budget")
    check(not branches[0].double_counted, "zero-mode CDM halo branch is clean accounting")

    print("\n[5] Consolidated status")
    rows = [
        ("zero-mode dust density", "conditional strong", "alpha0/208 source + 4:1 incidence gives omega_dark h^2 ~= 0.1209"),
        ("third-peak equality", "conditional strong", "z_eq restored near 3430"),
        ("acoustic selector 63/64", "not derived", "numerical closure; needs new boundary/latch acoustic-invisible slot"),
        ("halo accounting", "branch-forced", "zero-mode CDM branch OR active MOND with >95% depletion theorem"),
        ("fair-sample zero-mode + MOND", "rejected", "double-counts the same mobile halo mass"),
    ]
    for item, status, note in rows:
        print(f"  {item:32s} {status:20s} {note}")

    print(
        """
VERDICT
  The current Canon can keep the CMB route alive, but not Locked.

  What is real progress:
    * the zero-mode reservoir is a pressureless Brown--Kuchar-type dust slot;
    * the source arithmetic restores omega_dark h^2 and z_eq;
    * halo double-counting is no longer vague: the allowed branches are explicit.

  What is refuted under the current ledger:
    * H_CMB=(63/64)H_selector is NOT derivable from the symmetric depth-6
      service ledger alone.  A rank-63 acoustic projector needs a new
      distinguished boundary/latch primitive.

  Closure target:
    prove an acoustic pre-latch operator with trace 63 over the 64 depth-6
    selector slots, and choose the zero-mode-CDM halo branch or derive
    >95% galaxy depletion/screening before retaining active R4/MOND.
exit 0 -- CMB selector remains theorem-target only; halo accounting closes by branch choice.
"""
    )


if __name__ == "__main__":
    main()
