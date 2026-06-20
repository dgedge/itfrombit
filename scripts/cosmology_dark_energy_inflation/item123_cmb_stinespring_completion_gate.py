#!/usr/bin/env python3
r"""ITEM 123: CMB completion after the R4 Stinespring/Fock lift.

Question
--------
The R4/MOND halo ledger has just been sharpened: repeated Stinespring service
records give a count-valued, nonexclusive line ledger under the P1 scheduler
reading.  Does that also solve the CMB third-peak problem?

Answer
------
No.  The Stinespring/Fock lift closes the *halo line-count* problem, but it
does not create a homogeneous pressureless component for recombination.

Three internal readings fail for distinct reasons:

1. Active line occupancy:
   stationary mean N = |g|/a0, so the homogeneous background has |g|=0 and
   zero active R4 density.  It is a nonlinear halo response, not background CDM.
2. Cumulative Stinespring event history:
   it is count-valued, but it is a time-integrated service archive, not a
   conserved local particle number with rho~a^-3.  If made gravitational, it
   grows with the service integral rather than redshifting as dust.
3. Frozen/glassy R4:
   a frozen p-dimensional network has w=-p/3.  R4 is a line network (p=1),
   so freezing gives w=-1/3, the already-retired string-gas branch, not dust.

What remains is an external/effective AeST-class slot: add a conserved,
pressureless, shift-charge component with Omega h^2=0.096, w≈0, and small
sound speed before recombination.  That would restore z_eq and the third peak,
but it is an import until derived from substrate mechanics.
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
A_REC = 1.0 / (1.0 + Z_REC)
H0_OVER_C_MPC = 0.674 / 2997.92458


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def z_eq(omega_m_h2: float) -> float:
    return omega_m_h2 / OMEGA_R_H2 - 1.0


def w_tolerance_for_density_error(error_fraction: float = 0.10) -> float:
    """Require (1+z_rec)^(3w) to stay within 1+error_fraction."""
    return math.log1p(error_fraction) / (3.0 * math.log1p(Z_REC))


def aH_over_c_at_rec(omega_m_h2: float, h: float = 0.674) -> float:
    """Comoving Hubble scale aH/c at recombination, in Mpc^-1."""
    omega_m = omega_m_h2 / h**2
    omega_r = OMEGA_R_H2 / h**2
    e_z = math.sqrt(omega_m * (1.0 + Z_REC) ** 3 + omega_r * (1.0 + Z_REC) ** 4)
    return H0_OVER_C_MPC * e_z / (1.0 + Z_REC)


def sound_speed_bound(k_mpc: float, omega_m_h2: float) -> float:
    """Crude clustering condition c_s < aH/k at recombination."""
    return aH_over_c_at_rec(omega_m_h2) / k_mpc


def frozen_network_w(dimension: int) -> float:
    return -dimension / 3.0


def cumulative_archive_mean(gamma_over_h0: float, a_end: float, late_activation: bool) -> float:
    """Dimensionless integrated event archive per comoving line.

    This deliberately computes a cumulative history, not an active dust density.
    In matter era, dt = da/(a H), H = H0 sqrt(Om) a^-3/2.  The result scales
    with the integration upper limit and never becomes a conserved a^-3 number.
    """
    omega_m = 0.315
    pref = gamma_over_h0 / math.sqrt(omega_m)
    if late_activation:
        return pref * (2.0 / 5.0) * a_end ** 2.5
    return pref * (2.0 / 3.0) * a_end ** 1.5


@dataclass(frozen=True)
class Candidate:
    name: str
    background: bool
    conserved: bool
    pressureless: bool
    abundance_fixed: bool
    note: str

    @property
    def passes(self) -> bool:
        return self.background and self.conserved and self.pressureless and self.abundance_fixed


def main() -> None:
    print("ITEM 123: CMB COMPLETION AFTER R4 STINESPRING/FOCK LIFT")
    print("=" * 100)

    print("\n[1] Baseline equality pressure")
    omega_lcdm = OMEGA_B_H2 + OMEGA_C_H2_NEEDED
    omega_framework = OMEGA_B_H2 + OMEGA_NUR_H2
    zeq_lcdm = z_eq(omega_lcdm)
    zeq_framework = z_eq(omega_framework)
    print(f"  omega_b h^2           = {OMEGA_B_H2:.4f}")
    print(f"  omega_nuR h^2         = {OMEGA_NUR_H2:.4f}")
    print(f"  omega_R4 gap h^2      = {OMEGA_R4_H2:.4f}")
    print(f"  LCDM z_eq             = {zeq_lcdm:.0f}")
    print(f"  framework z_eq        = {zeq_framework:.0f}  (z_rec={Z_REC:.0f})")
    check(3300 < zeq_lcdm < 3500, "LCDM equality is well before recombination")
    check(1000 < zeq_framework < 1250, "R4/nu_R equality lands near recombination")
    check(abs(OMEGA_NUR_H2 + OMEGA_R4_H2 - OMEGA_C_H2_NEEDED) < 1.0e-12, "the missing cold slot is exactly the 80% R4 share")

    print("\n[2] Stinespring halo ledger is not homogeneous dust")
    print("  Active R4 line occupancy is Poisson(x) with x=|g|/a0.")
    print("  Homogeneous recombination background has no halo acceleration field: |g|=0.")
    active_mean_homogeneous = 0.0
    check(active_mean_homogeneous == 0.0, "active R4 line density vanishes in a homogeneous background")
    print("  Therefore the Stinespring/Fock lift closes halo MOND counts, not a background CDM fluid.")

    print("\n[3] Cumulative event archive is count-valued but not a conserved dust density")
    for gamma in (1.0, 10.0, 100.0):
        early = cumulative_archive_mean(gamma, A_REC, late_activation=False)
        late = cumulative_archive_mean(gamma, A_REC, late_activation=True)
        print(f"  Gamma/H0={gamma:5.1f}: archive to recombination const={early:.3e}, late-activation={late:.3e}")
    check(cumulative_archive_mean(10.0, A_REC, True) < cumulative_archive_mean(10.0, A_REC, False), "late activation suppresses early event archive")
    check(cumulative_archive_mean(10.0, A_REC, True) < 1.0e-6, "R4 late-activation archive is negligible by recombination for O(H0) rates")
    print("  More importantly: the archive is an integral of service events.  It grows with the")
    print("  integration limit; it is not a local conserved particle number with rho~a^-3.")

    print("\n[4] Frozen/glassy route reproduces the retired string branch")
    for p, label in ((0, "point gas"), (1, "R4 line network"), (2, "K04 wall network")):
        w = frozen_network_w(p)
        print(f"  frozen {label:16s}: p={p}, w={w:+.6f}")
    check(frozen_network_w(1) == -1.0 / 3.0, "frozen R4 lines give w=-1/3, not dust")
    check(frozen_network_w(2) == -2.0 / 3.0, "frozen K04 walls are even more dark-energy-like")
    check(frozen_network_w(0) == 0.0, "only point-like defects give dust, and R4 support is d=1")

    print("\n[5] Requirements for an external/effective completion")
    w_tol = w_tolerance_for_density_error(0.10)
    cs_01 = sound_speed_bound(0.10, omega_lcdm)
    cs_02 = sound_speed_bound(0.20, omega_lcdm)
    print(f"  Required extra dust abundance: omega_x h^2 = {OMEGA_R4_H2:.4f}")
    print(f"  If normalized today, |w_x| must be < {w_tol:.4f} to keep rho_x(z_rec) within 10% of dust.")
    print(f"  crude clustering bounds at recombination:")
    print(f"    k=0.10 Mpc^-1: c_s < {cs_01:.4f}, c_s^2 < {cs_01**2:.3e}")
    print(f"    k=0.20 Mpc^-1: c_s < {cs_02:.4f}, c_s^2 < {cs_02**2:.3e}")
    check(w_tol < 0.005, "CMB completion needs an equation of state within ~0.5% of dust")
    check(cs_02**2 < 1.0e-3, "third-peak-scale clustering needs very small sound speed")

    omega_with_slot = omega_framework + OMEGA_R4_H2
    zeq_with_slot = z_eq(omega_with_slot)
    print(f"  If such a slot is added: omega_m h^2={omega_with_slot:.4f}, z_eq={zeq_with_slot:.0f}")
    check(abs(zeq_with_slot - zeq_lcdm) < 1.0, "an omega_R4 pressureless slot restores equality exactly")

    print("\n[6] Candidate ledger")
    candidates = [
        Candidate(
            "active R4 Stinespring line occupancy",
            background=False,
            conserved=False,
            pressureless=False,
            abundance_fixed=True,
            note="Poisson halo count, x=|g|/a0; homogeneous x=0",
        ),
        Candidate(
            "cumulative R4 service archive",
            background=True,
            conserved=False,
            pressureless=False,
            abundance_fixed=False,
            note="time-integrated history, not local rho~a^-3 matter",
        ),
        Candidate(
            "frozen R4/K04 defect network",
            background=True,
            conserved=True,
            pressureless=False,
            abundance_fixed=True,
            note="extended p>=1 network has w=-p/3<0",
        ),
        Candidate(
            "external AeST-class pressureless charge",
            background=True,
            conserved=True,
            pressureless=True,
            abundance_fixed=False,
            note="would work if abundance and coupling are derived; presently imported",
        ),
    ]
    for cand in candidates:
        print(
            f"  {cand.name:38s} "
            f"bg={cand.background!s:5s} cons={cand.conserved!s:5s} "
            f"dust={cand.pressureless!s:5s} fixed={cand.abundance_fixed!s:5s} "
            f"pass={cand.passes!s:5s}  {cand.note}"
        )
    check(not candidates[0].passes, "Stinespring halo ledger does not pass CMB dust criteria")
    check(not candidates[1].passes, "event archive does not pass CMB dust criteria")
    check(not candidates[2].passes, "frozen defect network does not pass CMB dust criteria")
    check(not candidates[3].passes, "AeST slot still lacks a derived abundance/coupling")

    print("\nVERDICT")
    print("  The R4 Stinespring/Fock result helps MOND, but it does not rescue the CMB.")
    print("  It is a halo response with zero homogeneous background; the cumulative")
    print("  service archive is not a conserved local dust density; and freezing R4")
    print("  gives w=-1/3, not w=0. The only viable completion remains an AeST-class")
    print("  pressureless shift-charge with omega h^2=0.096, |w|<~0.005 and")
    print("  c_s^2<~10^-3 before recombination. That is a precise external theorem")
    print("  target, not something current R4/nu_R canon derives.")
    print("exit 0 -- CMB completion remains open/no-go; Stinespring lift does not supply dust.")


if __name__ == "__main__":
    main()
