#!/usr/bin/env python3
r"""ITEM 123: search for the missing pressureless CMB component.

The CMB third-peak gate is now quantitative.  With baryons plus the cold
sterile-nu_R relic, the framework has

    omega_m h^2 = 0.0224 + 0.024 = 0.0464,

so matter-radiation equality lands at recombination.  Restoring the LCDM
equality redshift requires an extra component

    omega_x h^2 = 0.096,
    |w_x| < 0.0045,
    c_s^2 < O(10^-3)

on third-peak scales.  This script asks which canon-native branches can meet
those gates.

Result
------
No existing branch passes.  The only phenomenologically viable form is a
constrained pressureless shift-charge, i.e. an AeST/Brown-Kuchar-style dust
variable.  A framework-native version would have to be a new R4 service-phase
or Stueckelberg charge with three derived clauses:

  1. exact shift symmetry, so d_mu J^mu = 0 and rho_x scales as a^-3;
  2. constrained dust Hamiltonian, not a normal scalar, so w = c_s^2 = 0;
  3. abundance fixed to omega_x h^2 = 0.096 by the already-declared 80% R4
     dark share, not by a fitted initial condition.

That is a sharp theorem target, not a closure.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from r4_eos_cmb_resolution import (
    OMEGA_B_H2,
    OMEGA_C_H2_NEEDED,
    OMEGA_NUR_H2,
    OMEGA_R4_H2,
)


T_CMB = 2.7255
N_EFF = 3.044
OMEGA_GAMMA_H2 = 2.469e-5 * (T_CMB / 2.7255) ** 4
NEUTRINO_FACTOR = 1.0 + (7.0 / 8.0) * (4.0 / 11.0) ** (4.0 / 3.0) * N_EFF
OMEGA_R_H2 = OMEGA_GAMMA_H2 * NEUTRINO_FACTOR
Z_REC = 1090.0
H0_OVER_C_MPC = 0.674 / 2997.92458


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def z_eq(omega_m_h2: float) -> float:
    return omega_m_h2 / OMEGA_R_H2 - 1.0


def w_tolerance(error_fraction: float = 0.10) -> float:
    return math.log1p(error_fraction) / (3.0 * math.log1p(Z_REC))


def aH_over_c_at_rec(omega_m_h2: float, h: float = 0.674) -> float:
    omega_m = omega_m_h2 / h**2
    omega_r = OMEGA_R_H2 / h**2
    e_z = math.sqrt(omega_m * (1.0 + Z_REC) ** 3 + omega_r * (1.0 + Z_REC) ** 4)
    return H0_OVER_C_MPC * e_z / (1.0 + Z_REC)


def sound_speed_bound(k_mpc: float, omega_m_h2: float) -> float:
    return aH_over_c_at_rec(omega_m_h2) / k_mpc


@dataclass(frozen=True)
class Candidate:
    name: str
    omega_h2: float | None
    background: bool
    conserved: bool
    pressureless: bool
    low_sound_speed: bool
    abundance_fixed: bool
    canon_derived: bool
    note: str

    @property
    def abundance_ok(self) -> bool:
        return self.omega_h2 is not None and abs(self.omega_h2 - OMEGA_R4_H2) < 1.0e-6

    @property
    def phenomenology_ok(self) -> bool:
        return (
            self.background
            and self.conserved
            and self.pressureless
            and self.low_sound_speed
            and self.abundance_ok
        )

    @property
    def canon_ok(self) -> bool:
        return self.phenomenology_ok and self.abundance_fixed and self.canon_derived


def frozen_network_w(p_dimension: int) -> float:
    return -p_dimension / 3.0


def scalar_k_essence_sound_speed(power_n: float) -> float:
    """For P(X) proportional X^n, w = c_s^2 = 1/(2n-1)."""
    return 1.0 / (2.0 * power_n - 1.0)


def main() -> None:
    print("ITEM 123: SEARCH FOR THE MISSING PRESSURELESS CMB SLOT")
    print("=" * 92)

    omega_lcdm = OMEGA_B_H2 + OMEGA_C_H2_NEEDED
    omega_fw = OMEGA_B_H2 + OMEGA_NUR_H2
    omega_needed = omega_lcdm - omega_fw
    zeq_fw = z_eq(omega_fw)
    zeq_lcdm = z_eq(omega_lcdm)
    zeq_with_x = z_eq(omega_fw + omega_needed)
    wtol = w_tolerance(0.10)
    cs2_010 = sound_speed_bound(0.10, omega_lcdm) ** 2
    cs2_020 = sound_speed_bound(0.20, omega_lcdm) ** 2

    print("\n[1] Target, recomputed")
    print(f"  omega_b h^2       = {OMEGA_B_H2:.4f}")
    print(f"  omega_nuR h^2     = {OMEGA_NUR_H2:.4f}")
    print(f"  missing omega_x   = {omega_needed:.4f}")
    print(f"  R4 declared share = {OMEGA_R4_H2:.4f}")
    print(f"  current z_eq      = {zeq_fw:.0f}")
    print(f"  LCDM z_eq         = {zeq_lcdm:.0f}")
    print(f"  with omega_x      = {zeq_with_x:.0f}")
    print(f"  |w_x| bound       = {wtol:.4f} for <10% density error by recombination")
    print(f"  c_s^2 bounds      = {cs2_010:.3e} at k=0.10/Mpc, {cs2_020:.3e} at k=0.20/Mpc")
    check(abs(omega_needed - OMEGA_R4_H2) < 1.0e-12, "the missing slot is exactly the declared 80% R4 dark share")
    check(abs(zeq_with_x - zeq_lcdm) < 1.0, "adding omega_x restores z_eq to LCDM")
    check(wtol < 0.005 and cs2_020 < 1.0e-3, "the component must be dust to high accuracy")

    print("\n[2] Existing canon branches")
    candidates = [
        Candidate(
            "sterile nu_R relic",
            OMEGA_NUR_H2,
            True,
            True,
            True,
            True,
            True,
            True,
            "valid CMB cold anchor, but only 0.024; fourfold short for the extra slot",
        ),
        Candidate(
            "active R4 Stinespring halo line",
            OMEGA_R4_H2,
            False,
            False,
            False,
            False,
            True,
            True,
            "Poisson(|g|/a0) halo response; homogeneous recombination has |g|=0",
        ),
        Candidate(
            "cumulative R4 service archive",
            None,
            True,
            False,
            False,
            False,
            False,
            True,
            "count history is time-integrated, not local conserved rho~a^-3",
        ),
        Candidate(
            "frozen R4 line network",
            OMEGA_R4_H2,
            True,
            True,
            frozen_network_w(1) == 0.0,
            False,
            True,
            True,
            "p=1 frozen network gives w=-1/3, the retired string branch",
        ),
        Candidate(
            "frozen K04 wall network",
            OMEGA_R4_H2,
            True,
            True,
            frozen_network_w(2) == 0.0,
            False,
            True,
            True,
            "p=2 frozen network gives w=-2/3",
        ),
        Candidate(
            "dilute 0D point-defect gas",
            OMEGA_R4_H2,
            True,
            True,
            True,
            True,
            False,
            False,
            "phenomenologically dust, but no canon derivation of species or abundance",
        ),
        Candidate(
            "external AeST pressureless shift-charge",
            OMEGA_R4_H2,
            True,
            True,
            True,
            True,
            False,
            False,
            "would solve the CMB, but is an import until derived from substrate mechanics",
        ),
        Candidate(
            "R4 service-phase / Stueckelberg dust",
            OMEGA_R4_H2,
            True,
            True,
            True,
            True,
            False,
            False,
            "best live target: derive shift symmetry, dust constraint, and abundance from R4 service algebra",
        ),
    ]

    for c in candidates:
        omega = "unset" if c.omega_h2 is None else f"{c.omega_h2:.4f}"
        print(
            f"  {c.name:38s} omega={omega:>6s} "
            f"bg={str(c.background):5s} cons={str(c.conserved):5s} "
            f"dust={str(c.pressureless):5s} cs={str(c.low_sound_speed):5s} "
            f"fixed={str(c.abundance_fixed):5s} derived={str(c.canon_derived):5s} "
            f"phen={str(c.phenomenology_ok):5s} canon={str(c.canon_ok):5s}"
        )
        print(f"    - {c.note}")

    check(not any(c.canon_ok for c in candidates), "no existing canon branch supplies the pressureless slot")
    live = [c for c in candidates if c.phenomenology_ok and not c.canon_ok]
    check(any(c.name == "R4 service-phase / Stueckelberg dust" for c in live), "the live route is a new constrained shift-charge/dust theorem")

    print("\n[3] Why a normal scalar is not enough")
    for n in (1, 2, 10, 100):
        cs2 = scalar_k_essence_sound_speed(float(n))
        print(f"  P(X) ~ X^{n:<3d}: w = c_s^2 = {cs2:.4f}")
    n_required = (1.0 + 1.0 / cs2_020) / 2.0
    print(f"  To get c_s^2 < {cs2_020:.3e} needs n > {n_required:.0f}, an extreme singular limit.")
    print("  Therefore the viable object is not a generic scalar field. It must be constrained")
    print("  dust: a shift variable plus a Lagrange multiplier / service-number density.")
    check(n_required > 800, "normal k-essence needs a singular high-power limit to mimic dust")

    print("\n[4] Minimal completion theorem target")
    clauses = [
        ("shift symmetry", "derive a service phase theta with theta -> theta + const and d_mu J^mu=0"),
        ("dust Hamiltonian", "derive a constrained action with p=0 and c_s^2=0, not P(X) radiation/stiffness"),
        ("abundance", "fix the conserved charge to omega_x h^2=0.096 from the 80% R4 ledger share"),
        ("separation", "avoid double-counting the late MOND line-current and homogeneous R4 dark-energy ledger"),
        ("Boltzmann", "insert the component in a Boltzmann code and recover the third peak"),
    ]
    for i, (name, clause) in enumerate(clauses, 1):
        print(f"  {i}. {name:16s}: {clause}")

    print("\nVERDICT")
    print("  The search finds no already-derived pressureless component. The only branch that")
    print("  has the right CMB phenomenology is a constrained shift-charge/dust sector carrying")
    print("  the existing 0.096 R4 share. If canon can derive it as an R4 service-phase")
    print("  Stueckelberg/Brown-Kuchar dust variable, z_eq is restored and the CMB route becomes")
    print("  live. Without those clauses it is an AeST-class import, not a framework result.")
    print("exit 0 -- pressureless-slot search complete; live target named, not closed.")


if __name__ == "__main__":
    main()
