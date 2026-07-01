#!/usr/bin/env python3
r"""ITEM 123 / 132: halo non-double-counting residual theorem.

The CMB/structure branch now uses the R4 zero-mode reservoir as CDM-like dust:

    omega_zero ~= 0.0967,   omega_nuR ~= 0.0242,   omega_dark ~= 0.1209.

Item 132 separately gives a conditional R4 Poisson/Fock line ledger that can
produce MOND/AQUAL behaviour.  The danger is double-counting: using the
zero-mode as the mobile halo mass and also fitting active R4 as an independent
baryonic MOND force in the same galaxy.

This script makes the accounting rule explicit.

Result
------
In current canon, the non-double-counted branch is:

    N_zero + nu_R = the mobile halo mass used by CMB/linear/nonlinear structure.
    N_active/R4   = exchange/polarisation on that reservoir, not an independent
                    second halo or independent MOND acceleration fit.

If one wants an independent active-R4 MOND/RAR branch instead, the fair-sample
zero-mode halo must be depleted/screened by >95%.  Current operator inventory
does not supply that depletion.  Therefore the live prediction is not
"LCDM halo + MOND"; it is zero-mode-CDM halo with at most a subleading active
R4 polarisation residual.

This is an accounting theorem, not a galaxy-formation simulation.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import sys

ROOT = __import__("pathlib").Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_code"))

from item123_cmb_cold_matter_closure_status import OMEGA_B_H2, RAR_SCATTER_DEX, density_for_power
from item123_halo_branch_fork_resolution import required_depletion_fraction


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


@dataclass(frozen=True)
class HaloLedger:
    name: str
    omega_h2: float
    gravitational_mass: bool
    independent_mond_force: bool
    allowed_with_zero_mode_halo: bool
    meaning: str


def ledgers() -> list[HaloLedger]:
    omega_nur, omega_zero, _omega_dark, _zeq = density_for_power(1)
    return [
        HaloLedger(
            "nu_R sterile relic",
            omega_nur,
            gravitational_mass=True,
            independent_mond_force=False,
            allowed_with_zero_mode_halo=True,
            meaning="subdominant mobile dust/WDM-like particle share",
        ),
        HaloLedger(
            "N_zero reservoir",
            omega_zero,
            gravitational_mass=True,
            independent_mond_force=False,
            allowed_with_zero_mode_halo=True,
            meaning="dominant Brown--Kuchar CDM-like mobile halo mass",
        ),
        HaloLedger(
            "active R4 line count",
            0.0,
            gravitational_mass=False,
            independent_mond_force=True,
            allowed_with_zero_mode_halo=False,
            meaning="late exchanged/polarisation line response; mean-zero as background",
        ),
        HaloLedger(
            "K04 fossil debris",
            0.0,
            gravitational_mass=False,
            independent_mond_force=False,
            allowed_with_zero_mode_halo=True,
            meaning="pinned/subdominant substrate fossil, not mobile halo mass",
        ),
    ]


def rar_room_from_scatter(scatter_dex: float) -> float:
    """Fractional acceleration room corresponding to a dex scatter budget."""

    return 10.0**scatter_dex - 1.0


def double_count_excess(zero_to_baryon: float) -> float:
    """If a MOND force already supplies the baryonic discrepancy, adding a
    fair-sample zero-mode halo adds this many baryon accelerations/masses."""

    return zero_to_baryon


def depletion_needed(excess: float, room: float) -> float:
    return 1.0 - room / excess


def active_r4_exchange_fraction(x_over_ntot: float) -> float:
    """Only admitted same-cell exchange: active fraction x/(N_tot+x)."""

    return x_over_ntot / (1.0 + x_over_ntot)


def x_over_ntot_for_fraction(fraction: float) -> float:
    return fraction / (1.0 - fraction)


def main() -> None:
    print("ITEM 123 / 132: HALO NON-DOUBLE-COUNTING RESIDUAL THEOREM")
    print("=" * 100)

    omega_nur, omega_zero, omega_dark, zeq = density_for_power(1)
    zero_to_baryon, scatter_room_old, max_fraction_old, depletion_old = required_depletion_fraction()
    room = rar_room_from_scatter(RAR_SCATTER_DEX)
    excess = double_count_excess(zero_to_baryon)
    depletion = depletion_needed(excess, room)
    max_fraction = 1.0 - depletion

    print("\n[1] Cosmic halo budget")
    print(f"  omega_b h^2       = {OMEGA_B_H2:.6f}")
    print(f"  omega_nuR h^2     = {omega_nur:.6f}")
    print(f"  omega_zero h^2    = {omega_zero:.6f}")
    print(f"  omega_dark h^2    = {omega_dark:.6f}")
    print(f"  z_eq              = {zeq:.1f}")
    print(f"  zero/baryon       = {zero_to_baryon:.3f}")
    check(abs(zero_to_baryon - omega_zero / OMEGA_B_H2) < 1e-12, "zero/baryon ratio computed from the same source ledger")
    check(omega_zero > 4.0 * OMEGA_B_H2, "zero-mode is a dominant halo component, not a perturbation")

    print("\n[2] Ledger roles")
    total_mass = 0.0
    for ledger in ledgers():
        if ledger.gravitational_mass:
            total_mass += ledger.omega_h2
        print(
            f"  {ledger.name:22s} omega={ledger.omega_h2:9.6f} "
            f"mass={str(ledger.gravitational_mass):5s} "
            f"independent_MOND={str(ledger.independent_mond_force):5s}"
        )
        print(f"    {ledger.meaning}")
    check(abs(total_mass - omega_dark) < 1e-12, "only nu_R plus N_zero carry the mobile dark mass in the adopted branch")

    print("\n[3] Non-double-counting inequality")
    print(f"  RAR scatter proxy                         = {RAR_SCATTER_DEX:.2f} dex -> room {room:.3f}")
    print(f"  fair-sample zero-mode excess if MOND is already fitted = {excess:.3f} baryon masses")
    print(f"  max allowed fair-sample fraction in MOND branch        = {max_fraction:.3f}")
    print(f"  required zero-mode depletion/screening                 = {depletion:.1%}")
    check(abs(room - scatter_room_old) < 1e-12, "scatter convention matches item123_halo_branch_fork_resolution.py")
    check(abs(max_fraction - max_fraction_old) < 1e-12, "max fraction convention matches prior fork audit")
    check(abs(depletion - depletion_old) < 1e-12, "depletion convention matches prior fork audit")
    check(depletion > 0.95, "an independent active-R4 MOND fit needs >95% zero-mode depletion")

    print("\n[4] Same-cell exchange cannot be the depletion theorem")
    need = x_over_ntot_for_fraction(depletion)
    print("  admitted exchange form: f_active = x/(N_tot+x)")
    print(f"  f={depletion:.3f} requires x/N_tot={need:.2f}")
    for ratio in (1.0e-3, 1.0e-2, 1.0e-1, 1.0, need):
        print(f"    x/N_tot={ratio:8.3g} -> f_active={active_r4_exchange_fraction(ratio):.5f}")
    check(need > 20.0, "depletion would require x much larger than the reservoir")
    check(active_r4_exchange_fraction(1.0e-2) < 0.011, "Poisson/MOND reservoir limit N_tot>>x gives negligible depletion")
    print("  The Poisson/Fock MOND line ledger needs the reservoir limit N_tot >> x.")
    print("  The depletion branch needs x/N_tot >> 1. These are incompatible readings")
    print("  of the same admitted exchange operator.")

    print("\n[5] Branch verdict")
    branches = {
        "zero-mode-CDM": "ADOPTED: CMB + structure + halo mass; active R4 only residual/polarisation",
        "independent MOND": "BLOCKED: needs >95% depletion/screening operator not in current canon",
        "zero-mode + independent MOND": "REJECTED: double-counts the same galaxy discrepancy",
    }
    for name, verdict in branches.items():
        print(f"  {name:24s} {verdict}")

    print(
        """
  Consequence:
    The live galaxy/halo programme is a standard zero-mode-CDM halo programme
    with possible subleading R4 polarisation, not an independent MOND force
    added on top.  To reopen the MOND branch, canon must derive a local
    depletion/screening operator that removes >95% of N_zero from galaxies
    without destroying the R4 Poisson line ledger.
exit 0 -- halo non-double-counting is resolved by residual/polarisation accounting.
"""
    )


if __name__ == "__main__":
    main()
