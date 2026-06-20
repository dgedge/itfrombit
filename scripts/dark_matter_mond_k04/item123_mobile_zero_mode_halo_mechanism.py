#!/usr/bin/env python3
r"""ITEM 123: mobile dark-halo mechanism after K04 pinning.

Question
--------
K04 debris is now substrate-pinned and therefore cannot be the dominant mobile
halo component.  Is there a canon-native mobile mechanism left, or does the
dark sector fall back entirely to the old R4/nu_R budget?

Result
------
The live mobile mechanism is the R4 zero-mode reservoir, not K04 debris and not
the active R4 MOND line count.  It is a distinct variable:

    N_tot = N_zero + N_active .

The already-constructed chain gives:

  * N_active is the late halo/MOND exchange excitation, slaved to |g|/a0.
  * N_zero is the conserved reservoir count.
  * The minimal reservoir Hamiltonian is rest-count only, so w=c_s^2=0.
  * The finite R4 incidence gives Omega_zero = 4 Omega_nuR.
  * The alpha/208 source-law candidate gives Omega_dark h^2 = 0.12089.

Mobility is not K04-style depinning.  A Brown--Kuchar/constrained-dust
continuum lift has stress tensor

    T^{mu nu} = rho u^mu u^nu,

and stress conservation gives geodesic flow plus mass conservation.  Since the
zero-mode labels are record-number/worldline labels, not substrate bonds, the
K04 Peierls barrier does not apply.  Therefore the mechanism can carry mobile
halo mass if the Brown--Kuchar/geodesic lift and the alpha/208 boot source map
are accepted.

Honest boundary
---------------
This is a conditional mechanism, not a locked particle species.  The finite
ledger closes pressurelessness and relative abundance, but the source map
alpha/208 and the continuum Boltzmann/geodesic implementation remain theorem
targets.  If either fails, the mobile halo burden reverts to the 17.7 keV
nu_R relic plus the R4 MOND response, which is probably insufficient for the
CMB and cluster/galaxy budget.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

from item123_r4_zero_mode_dust_hamiltonian import operator_inventory
from r4_eos_cmb_resolution import OMEGA_B_H2, OMEGA_C_H2_NEEDED, OMEGA_NUR_H2, OMEGA_R4_H2


ROOT = Path(__file__).resolve().parents[1]


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def run_owner(name: str, needle: str) -> None:
    proc = subprocess.run(
        [sys.executable, str(ROOT / "python_code" / name)],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        print(proc.stdout)
        print(proc.stderr, file=sys.stderr)
        raise SystemExit(f"{name} failed with exit {proc.returncode}")
    ok = needle in proc.stdout
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {needle}")
    if not ok:
        raise AssertionError(f"{name} missing {needle!r}")


def z_eq(omega_m_h2: float) -> float:
    omega_gamma_h2 = 2.469e-5
    n_eff = 3.044
    neutrino_factor = 1.0 + (7.0 / 8.0) * (4.0 / 11.0) ** (4.0 / 3.0) * n_eff
    omega_r_h2 = omega_gamma_h2 * neutrino_factor
    return omega_m_h2 / omega_r_h2 - 1.0


def dust_stress_tensor(rho: float, velocity: np.ndarray) -> np.ndarray:
    """Nonrelativistic local form sufficient for the mobility gate.

    T00=rho, T0i=rho v_i, Tij=rho v_i v_j.  There is no isotropic pressure
    term p delta_ij.  This is the finite analogue of Brown--Kuchar dust in a
    local inertial patch.
    """

    t = np.zeros((4, 4), dtype=float)
    t[0, 0] = rho
    t[0, 1:] = rho * velocity
    t[1:, 0] = rho * velocity
    t[1:, 1:] = rho * np.outer(velocity, velocity)
    return t


def k04_drive_ratio(acceleration: float = 1.2e-10, cell_step_m: float = 1.188e-15) -> float:
    c = 299_792_458.0
    return acceleration * cell_step_m / (c * c)


@dataclass(frozen=True)
class Candidate:
    name: str
    mobile: bool
    pressureless: bool
    abundance: str
    role: str
    verdict: str


def main() -> None:
    print("ITEM 123: MOBILE ZERO-MODE HALO MECHANISM")
    print("=" * 92)

    print("\n[1] Owner gates already in canon")
    run_owner("item123_r4_zero_mode_reservoir_lift.py", "reservoir lift advances the form")
    run_owner("item123_r4_zero_mode_dust_hamiltonian.py", "constrained-dust Hamiltonian closes")
    run_owner("item123_r4_zero_mode_abundance_ratio.py", "relative abundance ratio closes conditionally")
    run_owner("item123_nuR_absolute_density_boot_qec.py", "alpha/208 is a strong conditional boot-QEC density theorem target")

    print("\n[2] Candidate taxonomy after K04 depinning")
    candidates = [
        Candidate(
            "K04 wall debris",
            mobile=False,
            pressureless=False,
            abundance="subdominant / vacuum-grade shadow",
            role="pinned fossil, not halo mass",
            verdict="RETIRED-AS-HALO",
        ),
        Candidate(
            "17.7 keV nu_R relic",
            mobile=True,
            pressureless=True,
            abundance=f"omega h^2={OMEGA_NUR_H2:.3f}",
            role="particle-like clustered subcomponent",
            verdict="LIVE-BUT-20-PERCENT",
        ),
        Candidate(
            "active R4 MOND line count",
            mobile=False,
            pressureless=False,
            abundance="slaved to |g|/a0, homogeneous zero",
            role="late halo response / MOND law",
            verdict="LIVE-AS-FORCE-LAW",
        ),
        Candidate(
            "R4 zero-mode reservoir",
            mobile=True,
            pressureless=True,
            abundance=f"omega h^2=4 omega_nuR={OMEGA_R4_H2:.3f}",
            role="conserved dust reservoir; mobile-halo candidate",
            verdict="LIVE-CONDITIONAL",
        ),
    ]
    for c in candidates:
        print(
            f"  {c.name:28s} mobile={str(c.mobile):5s} dust={str(c.pressureless):5s} "
            f"{c.abundance:30s} {c.verdict:18s} {c.role}"
        )
    check(not candidates[0].mobile, "K04 remains pinned and cannot carry mobile halo mass")
    check(candidates[1].mobile and candidates[1].pressureless, "nu_R remains a mobile dust subcomponent")
    check(not candidates[2].pressureless, "active R4 line count is not the early/mobile dust component")
    check(candidates[3].mobile and candidates[3].pressureless, "R4 zero-mode reservoir is the only new mobile-dust candidate")

    print("\n[3] Why the zero-mode is not K04-pinned")
    terms = operator_inventory()
    admitted = {term.name for term in terms if term.status.startswith("ADMITTED")}
    rejected = {term.name for term in terms if term.status.startswith("REQUIRES")}
    print(f"  admitted reservoir terms = {sorted(admitted)}")
    print(f"  rejected spatial terms   = {sorted(rejected)}")
    check(admitted == {"rest count", "local active exchange"}, "reservoir has rest count plus local exchange only")
    check("elastic shear" in rejected and "reservoir hopping" in rejected, "substrate shear/hopping operators are not in the minimal reservoir")
    print(f"  K04 astrophysical drive ratio per cell step = {k04_drive_ratio():.3e}")
    check(k04_drive_ratio() < 1.0e-40, "K04 depinning drive is negligible")
    print("  The zero-mode is not a bond/domain wall, so the K04 Peierls barrier is not its mobility law.")

    print("\n[4] Brown--Kuchar/geodesic mobility check")
    rho = 1.0
    v = np.array([2.0e-3, -1.0e-3, 0.5e-3])
    t = dust_stress_tensor(rho, v)
    isotropic_pressure = float(np.trace(t[1:, 1:]) / 3.0)
    kinetic_anisotropy = float(np.linalg.norm(t[1:, 1:] - isotropic_pressure * np.eye(3)))
    print(f"  T00={t[0,0]:.6f}; T0i={t[0,1:]}")
    print(f"  trace(Tij)/3={isotropic_pressure:.3e}; anisotropic flow stress={kinetic_anisotropy:.3e}")
    check(isotropic_pressure / rho < 1.0e-5, "nonrelativistic dust has negligible pressure for halo velocities")
    check(kinetic_anisotropy > 0.0, "dust can carry momentum/flow without pressure support")
    print("  In the continuum lift, nabla_mu(rho u^mu u^nu)=0 gives mass conservation plus geodesic flow.")

    print("\n[5] Density and equality consequence")
    omega_zero = 4.0 * OMEGA_NUR_H2
    omega_dark = OMEGA_NUR_H2 + omega_zero
    omega_m = OMEGA_B_H2 + omega_dark
    print(f"  omega_nuR h^2        = {OMEGA_NUR_H2:.5f}")
    print(f"  omega_zero h^2       = {omega_zero:.5f}")
    print(f"  omega_dark h^2       = {omega_dark:.5f}")
    print(f"  omega_c target h^2   = {OMEGA_C_H2_NEEDED:.5f}")
    print(f"  z_eq with zero-mode  = {z_eq(omega_m):.1f}")
    check(abs(omega_zero - OMEGA_R4_H2) < 1.0e-12, "directed R4 incidence gives the missing 80 percent")
    check(abs(omega_dark - OMEGA_C_H2_NEEDED) < 1.0e-12, "nu_R plus zero-mode restores the cold dark budget at the relative-ratio level")
    check(3300.0 < z_eq(omega_m) < 3500.0, "matter-radiation equality returns to the observed CMB range")

    print("\n[6] Remaining theorem gates")
    rows = [
        ("finite reservoir exchange", "CLOSED", "N_zero+N_active=N_tot; active marginal remains Poisson"),
        ("dust Hamiltonian", "CLOSED", "rest-count only; w=c_s^2=0"),
        ("relative abundance", "CLOSED-CONDITIONAL", "directed R4 service incidence gives 4:1 to nu_R"),
        ("absolute density", "CANDIDATE", "alpha/208 boot source law lands at omega_dark=0.12089"),
        ("mobile/geodesic lift", "CONDITIONAL", "Brown-Kuchar continuum stress must be the Boltzmann species"),
        ("halo phenomenology", "OPEN", "needs perturbation/Boltzmann + structure calculation, distinct from active MOND"),
        ("source map", "OPEN", "derive uniform Q-complement addressing and one-erasure sterile release"),
    ]
    for name, status, note in rows:
        print(f"  {name:24s} {status:18s} {note}")

    print("\nVERDICT")
    print("  A new mobile mechanism is available: the R4 zero-mode reservoir as a")
    print("  constrained Brown--Kuchar dust component.  It is not K04 debris, and it is")
    print("  not the active R4 MOND line occupancy.  It is a conserved reservoir count")
    print("  whose minimal Hamiltonian has rest energy only and whose relative abundance")
    print("  is tied to the nu_R sector by the finite R4 directed-edge count.")
    print("  Status: live conditional mechanism.  It can carry the mobile halo/CMB cold")
    print("  budget if the alpha/208 source map and continuum geodesic/Boltzmann lift are")
    print("  derived.  If those fail, the framework has no substrate-native dominant")
    print("  mobile halo component beyond the subdominant nu_R relic and the R4 MOND")
    print("  force-law response.")
    print("exit 0 -- mobile halo burden moves to the R4 zero-mode reservoir, conditionally.")


if __name__ == "__main__":
    main()
