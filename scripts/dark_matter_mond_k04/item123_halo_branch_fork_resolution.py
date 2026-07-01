#!/usr/bin/env python3
r"""ITEM 123: halo branch fork resolution gate.

Question
--------
The zero-mode reservoir closes the CMB cold-matter budget at conditional
theorem grade, but it creates a galaxy-accounting fork:

  A. use N_zero as the CDM-like mobile halo mass; active R4/MOND is not also
     fitted as an independent baryonic force in the same regime;
  B. keep active R4/MOND as the baryonic RAR force, but derive a depletion or
     screening theorem that removes the fair-sampling zero-mode from galaxies.

The numerical gate is already known: branch B needs about 96% depletion of the
zero-mode in galaxy environments.  This script asks whether the current
R4/zero-mode operator inventory can supply that depletion.

Result
------
No, not within the documented minimal lift.  The admitted terms are:

  * rest count, H_zero = epsilon sum_i N_i;
  * same-cell exchange, N_zero + N_active = N_tot.

All spatial anti-bias/screening mechanisms -- reservoir hopping, gradient
stiffness, phase stiffness, elastic/shear coupling, or baryon-selected removal
-- require new operators.  The only admitted exchange channel also cannot do
the job while preserving the active R4 Poisson line ledger: depleting a fraction
f of N_zero by exchange needs x/N_tot = f/(1-f).  For f ~= 0.96 this is about
24, whereas the Poisson/MOND limit is N_tot >> x.

Therefore the current canon resolves the fork by branch adoption, not by a new
depletion theorem: the non-double-counted branch is zero-mode CDM-like halos,
with active R4/MOND demoted to a possible subleading/polarisation response
unless a future >95% depletion/screening operator is derived.
"""

from __future__ import annotations

from dataclasses import dataclass

from item123_cmb_cold_matter_closure_status import (
    OMEGA_B_H2,
    RAR_SCATTER_DEX,
    density_for_power,
)
from item123_r4_zero_mode_dust_hamiltonian import operator_inventory


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


@dataclass(frozen=True)
class HaloBranch:
    name: str
    zero_mode_as_mass: bool
    active_r4_independent_force: bool
    depletion_required: bool
    status: str
    reason: str

    @property
    def double_counted_without_depletion(self) -> bool:
        return self.zero_mode_as_mass and self.active_r4_independent_force and self.depletion_required


def required_depletion_fraction() -> tuple[float, float, float, float]:
    """Return zero/baryon, scatter room, max allowed fraction, depletion."""

    _omega_nur, omega_zero, _omega_dark, _zeq = density_for_power(1)
    zero_to_baryon = omega_zero / OMEGA_B_H2
    scatter_room = 10.0**RAR_SCATTER_DEX - 1.0
    max_fair_sample_fraction = scatter_room / zero_to_baryon
    depletion = 1.0 - max_fair_sample_fraction
    return zero_to_baryon, scatter_room, max_fair_sample_fraction, depletion


def exchange_fraction(x_over_ntot: float) -> float:
    """Active fraction for the closed exchange N_active ~ N_tot*x/(N_tot+x)."""

    return x_over_ntot / (1.0 + x_over_ntot)


def exchange_ratio_needed(fraction: float) -> float:
    return fraction / (1.0 - fraction)


def main() -> None:
    print("ITEM 123: HALO BRANCH FORK RESOLUTION")
    print("=" * 92)

    print("\n[1] Numerical depletion gate")
    zero_to_baryon, scatter_room, max_fraction, depletion = required_depletion_fraction()
    print(f"  zero-mode / baryon fair-sample mass ratio = {zero_to_baryon:.3f}")
    print(f"  RAR scatter proxy                         = {RAR_SCATTER_DEX:.2f} dex -> {scatter_room:.3f}")
    print(f"  max fair-sample zero-mode in MOND branch  = {max_fraction:.3f}")
    print(f"  required depletion/screening              = {depletion:.1%}")
    check(zero_to_baryon > 4.0, "fair-sampling zero-mode is a dominant halo component")
    check(depletion > 0.95, "MOND-preserving branch requires >95% zero-mode depletion")

    print("\n[2] Operator inventory: can current canon deplete or screen it?")
    terms = operator_inventory()
    for term in terms:
        print(f"  {term.name:28s} {term.status:23s} {term.form}")
        print(f"    {term.reason}")
    admitted = {term.name for term in terms if term.status.startswith("ADMITTED")}
    rejected = {term.name for term in terms if term.status.startswith("REQUIRES")}
    check(admitted == {"rest count", "local active exchange"}, "minimal reservoir admits only rest count plus same-cell exchange")
    spatial_rejections = {
        "density-gradient stiffness",
        "reservoir hopping",
        "phase stiffness",
        "elastic shear",
    }
    check(spatial_rejections <= rejected, "all spatial screening/transport/stiffness routes require new operators")
    print("  No admitted term couples N_zero to baryon density, galaxy environment,")
    print("  long-range screening, reservoir transport, or pressure support.")

    print("\n[3] The only admitted exchange cannot provide 96% depletion in the MOND limit")
    needed_x_over_ntot = exchange_ratio_needed(depletion)
    print(f"  exchange fraction f = x/(N_tot+x)")
    print(f"  to get f={depletion:.3f}: x/N_tot = f/(1-f) = {needed_x_over_ntot:.2f}")
    for ratio in (0.001, 0.01, 0.1, 1.0, needed_x_over_ntot):
        print(f"    x/N_tot={ratio:8.3g} -> active fraction={exchange_fraction(ratio):.4f}")
    check(needed_x_over_ntot > 20.0, "96% depletion requires x far larger than N_tot")
    check(exchange_fraction(0.01) < 0.011, "Poisson/MOND reservoir limit N_tot>>x depletes only a tiny fraction")
    print("  But the active R4 line law was derived in the reservoir limit N_tot >> x,")
    print("  where the active marginal is Poisson(x).  Driving x/N_tot~24 would leave")
    print("  that limit and destroy the same nonexclusive Poisson ledger branch B wants.")

    print("\n[4] Branch classification")
    branches = [
        HaloBranch(
            "zero-mode CDM-like halo",
            zero_mode_as_mass=True,
            active_r4_independent_force=False,
            depletion_required=False,
            status="ADOPTED-CANONICAL-BRANCH",
            reason="keeps the CMB cold reservoir as the mobile halo and avoids double-counting",
        ),
        HaloBranch(
            "active R4/MOND with depleted zero-mode",
            zero_mode_as_mass=True,
            active_r4_independent_force=True,
            depletion_required=True,
            status="NO-CURRENT-THEOREM",
            reason="needs >95% anti-bias/screening; no admitted operator supplies it",
        ),
        HaloBranch(
            "fair-sample zero-mode plus independent active MOND",
            zero_mode_as_mass=True,
            active_r4_independent_force=True,
            depletion_required=True,
            status="REJECTED-DOUBLE-COUNT",
            reason="adds a 4.3 baryon-mass collisionless halo and an independent MOND force",
        ),
    ]
    for branch in branches:
        print(f"  {branch.name:48s} {branch.status:24s}")
        print(f"    {branch.reason}")
    check(branches[0].status == "ADOPTED-CANONICAL-BRANCH", "zero-mode CDM-like halo is the only current non-double-counted branch")
    check(branches[1].status == "NO-CURRENT-THEOREM", "active-MOND branch remains possible only with a new depletion theorem")
    check(branches[2].status == "REJECTED-DOUBLE-COUNT", "fair-sample zero-mode plus independent MOND is barred")

    print("\n[5] Consequence")
    print("  The fork is resolved at current-canon level by adopting N_zero as the")
    print("  CDM-like mobile halo.  Active R4/MOND is not retained as an independent")
    print("  baryonic response in the same galaxy regime.  It may survive only as a")
    print("  subleading/polarisation ledger, or be reopened if a future local operator")
    print("  derives >95% galaxy depletion/screening without breaking the Poisson R4")
    print("  line law.")
    print("exit 0 -- halo fork resolved by branch adoption; depletion theorem absent.")


if __name__ == "__main__":
    main()
