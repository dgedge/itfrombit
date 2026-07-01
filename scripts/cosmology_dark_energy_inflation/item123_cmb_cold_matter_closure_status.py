#!/usr/bin/env python3
r"""ITEM 123: cold-matter / CMB third-peak closure status audit.

This is a consolidation gate, not another fit.

The framework's CMB problem used to be a hard factor-five deficit:

    sterile nu_R only:       omega_c h^2 ~= 0.024
    observed cold slot:      omega_c h^2 ~= 0.120

The live repair is

    omega_nuR      from n_nuR/n_gamma = alpha0/208 and m_nuR=alpha0^2 Lambda_QCD,
    omega_zero     = 4 omega_nuR from directed R4 service-edge incidence,
    omega_dark     = 5 omega_nuR.

This script answers the narrower status question:

  * What is actually closed by the current substrate algebra?
  * What is only conditional on R14/item-79 alpha0 billing and item-118 mass?
  * Does the zero-mode double-count active R4/MOND?

Verdict
-------
The third-peak cold-matter budget is conditionally closed:

  * |Q|=208, uniform Q addressing, one S3-singlet sterile release port, the
    Brown--Kuchar dust form, and the 4:1 zero-mode/nu_R split are derived by
    companion scripts.
  * The alpha0 power is not selected by the CMB number alone.  It is licensed
    only if the R14/item-79 monitored-service billing theorem is admitted:
    one non-unitary sterile release bills one bare alpha0 firing.
  * The absolute landing also inherits the item-118 sterile mass
    m_nuR=alpha0^2 Lambda_QCD.
  * Halo double-counting is avoided only by branch discipline: either use the
    zero-mode as the CDM-like mobile halo mass, or retain active R4/MOND only
    after deriving strong galaxy depletion/screening of the zero-mode.
"""

from __future__ import annotations

from dataclasses import dataclass
import math


ALPHA0 = 1.0 / 137.035999
LAMBDA_QCD_EV = 0.332e9
T_CMB = 2.7255
ZETA3 = 1.2020569031595943
KB_EV_PER_K = 8.617333262e-5
HBARC_EV_CM = 1.973269804e-5
RHO_CRIT_H2_EV_CM3 = 1.05371e4
N_EFF = 3.044

OMEGA_B_H2 = 0.0224
OMEGA_C_OBS_H2 = 0.1200
OMEGA_NUR_TARGET_H2 = 0.0240
RAR_SCATTER_DEX = 0.07


@dataclass(frozen=True)
class Gate:
    name: str
    status: str
    reason: str


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def photon_density_cm3(t_cmb: float = T_CMB) -> float:
    return 2.0 * ZETA3 / math.pi**2 * (KB_EV_PER_K * t_cmb / HBARC_EV_CM) ** 3


def omega_from_ratio(n_over_ngamma: float, mass_ev: float) -> float:
    return n_over_ngamma * photon_density_cm3() * mass_ev / RHO_CRIT_H2_EV_CM3


def omega_r_h2() -> float:
    omega_gamma = 2.469e-5 * (T_CMB / 2.7255) ** 4
    neutrino_factor = 1.0 + (7.0 / 8.0) * (4.0 / 11.0) ** (4.0 / 3.0) * N_EFF
    return omega_gamma * neutrino_factor


def z_eq(omega_m_h2: float) -> float:
    return omega_m_h2 / omega_r_h2() - 1.0


def density_for_power(power: int, denominator: int = 208) -> tuple[float, float, float, float]:
    """Return omega_nuR, omega_zero, omega_dark, z_eq for alpha0^power/denominator."""

    mass = ALPHA0**2 * LAMBDA_QCD_EV
    omega_nur = omega_from_ratio(ALPHA0**power / denominator, mass)
    omega_zero = 4.0 * omega_nur
    omega_dark = omega_nur + omega_zero
    return omega_nur, omega_zero, omega_dark, z_eq(OMEGA_B_H2 + omega_dark)


def denominator_hits(target: float, tolerance: float = 0.01) -> list[int]:
    mass = ALPHA0**2 * LAMBDA_QCD_EV
    hits: list[int] = []
    for denominator in range(1, 513):
        omega_nur = omega_from_ratio(ALPHA0 / denominator, mass)
        if abs(omega_nur / target - 1.0) < tolerance:
            hits.append(denominator)
    return hits


def required_denominator(target: float) -> float:
    mass = ALPHA0**2 * LAMBDA_QCD_EV
    return ALPHA0 / (target * RHO_CRIT_H2_EV_CM3 / (mass * photon_density_cm3()))


def main() -> None:
    print("ITEM 123: CMB COLD-MATTER / THIRD-PEAK CLOSURE STATUS")
    print("=" * 96)

    print("\n[1] Cold-matter budget")
    omega_nur, omega_zero, omega_dark, zeq_full = density_for_power(1)
    zeq_nur_only = z_eq(OMEGA_B_H2 + omega_nur)
    print(f"  alpha0                 = {ALPHA0:.9e}")
    print(f"  n_gamma(T0)            = {photon_density_cm3():.3f} cm^-3")
    print(f"  m_nuR=alpha0^2 Lambda  = {ALPHA0**2 * LAMBDA_QCD_EV / 1e3:.4f} keV")
    print(f"  omega_nuR h^2          = {omega_nur:.6f}")
    print(f"  omega_zero h^2         = {omega_zero:.6f}")
    print(f"  omega_dark h^2         = {omega_dark:.6f}")
    print(f"  z_eq with nu_R only    = {zeq_nur_only:.1f}")
    print(f"  z_eq with zero-mode    = {zeq_full:.1f}")
    check(abs(omega_dark / OMEGA_C_OBS_H2 - 1.0) < 0.01, "nu_R plus 4:1 zero-mode lands within 1% of the cold dark slot")
    check(1000.0 < zeq_nur_only < 1250.0, "sterile-only equality is near recombination: the old third-peak failure")
    check(3400.0 < zeq_full < 3450.0, "zero-mode budget restores LCDM-like equality")

    print("\n[2] Alpha-power audit")
    rows = []
    for power in (0, 1, 2):
        row = (power,) + density_for_power(power)
        rows.append(row)
        print(
            f"  p={power}: n_nuR/n_gamma=alpha0^{power}/208, "
            f"omega_nuR={row[1]:.6f}, omega_dark={row[3]:.6f}, "
            f"dark/obs={row[3] / OMEGA_C_OBS_H2:.3f}"
        )
    p0, p1, p2 = rows
    check(p0[3] / p1[3] > 130.0, "p=0 overproduces p=1 by one bare-alpha inverse")
    check(p1[3] / p2[3] > 130.0, "p=2 underproduces p=1 by one bare-alpha inverse")
    check(abs(p1[3] / OMEGA_C_OBS_H2 - 1.0) < 0.01, "p=1 is the only p in {0,1,2} compatible with the CMB cold slot")
    print("  Interpretation: this discriminates the power once the candidate family is")
    print("  given.  It is not an independent substrate proof of p=1.  The proof comes")
    print("  only from R14/item-79: one non-unitary sterile release bills one bare-alpha")
    print("  monitored-service firing.")

    print("\n[3] Denominator audit")
    hits = denominator_hits(OMEGA_NUR_TARGET_H2)
    req = required_denominator(OMEGA_NUR_TARGET_H2)
    print(f"  inverse denominator for omega_nuR={OMEGA_NUR_TARGET_H2:.3f}: {req:.3f}")
    print(f"  integer denominators <=512 within 1%: {hits}")
    check(208 in hits, "|Q|=208 is in the observational 1% band")
    check(len(hits) > 1, "the number alone is not unique; |Q|=208 must be structurally derived")
    print("  Current structural status: |Q|=256-48=208, connected monitored Q-service")
    print("  gives uniform addressing, and the sterile source is one S3 singlet port.")

    print("\n[4] Halo double-counting branch audit")
    zero_to_baryon = omega_zero / OMEGA_B_H2
    dark_to_baryon = omega_dark / OMEGA_B_H2
    scatter_room = 10.0**RAR_SCATTER_DEX - 1.0
    max_zero_fraction_with_active_mond = scatter_room / zero_to_baryon
    required_depletion = 1.0 - max_zero_fraction_with_active_mond
    print(f"  zero-mode / baryon mass ratio       = {zero_to_baryon:.3f}")
    print(f"  total dark / baryon mass ratio      = {dark_to_baryon:.3f}")
    print(f"  RAR scatter proxy                   = {RAR_SCATTER_DEX:.2f} dex -> {scatter_room:.3f} fractional room")
    print(f"  if active MOND is retained, fair-sample zero-mode fraction <= {max_zero_fraction_with_active_mond:.3f}")
    print(f"  required galaxy depletion/screening >= {required_depletion:.1%}")
    check(zero_to_baryon > 4.0, "fair-sample zero-mode is a dominant mobile halo component")
    check(required_depletion > 0.95, "active-MOND branch needs >95% zero-mode galaxy depletion/screening")
    print("  Allowed accounting branch A: zero-mode is the CDM-like mobile halo mass;")
    print("  active R4/MOND is not also fitted as an independent galaxy force.")
    print("  Allowed accounting branch B: keep active R4/MOND only after deriving")
    print("  strong galaxy depletion/screening of the zero-mode.")
    print("  Forbidden branch: fair-sample zero-mode plus independent active MOND.")

    print("\n[5] Closure ledger")
    gates = [
        Gate(
            "third-peak cold budget",
            "CONDITIONALLY CLOSED",
            "p=1, m_nuR=alpha0^2 Lambda, and 4:1 incidence give omega_dark h^2=0.1209 and z_eq~=3430",
        ),
        Gate(
            "Q denominator",
            "DERIVED",
            "|Q|=256-48=208 and monitored connected Q-service gives uniform addressing",
        ),
        Gate(
            "sterile port count",
            "DERIVED",
            "generation-blind source algebra forces one S3 bright port, not three generation ports",
        ),
        Gate(
            "alpha0 power",
            "CONDITIONAL ON R14/ITEM-79",
            "CMB selects p=1 inside the candidate family, but the substrate theorem is one non-unitary firing bills one bare alpha0",
        ),
        Gate(
            "zero-mode dust form",
            "DERIVED AT EFFECTIVE-ACTION GRADE",
            "minimal reservoir inventory admits rest count and local exchange only, giving Brown-Kuchar w=cs^2=0 dust",
        ),
        Gate(
            "halo double-counting",
            "BRANCH-CLOSED, PHENOMENOLOGY OPEN",
            "accounting is clean only after choosing zero-mode CDM halos or deriving >95% galaxy depletion for active MOND",
        ),
        Gate(
            "full CMB likelihood / theta star",
            "NOT CLOSED HERE",
            "third-peak equality is supplied; acoustic-scale selector and full perturbation implementation remain separate gates",
        ),
    ]
    for gate in gates:
        print(f"  {gate.name:32s} {gate.status:36s} {gate.reason}")

    print("\nVERDICT")
    print("  The cold-matter / third-peak issue can be promoted from sharp no-go to")
    print("  conditional closure.  The old factor-five deficit is gone if the current")
    print("  R14/item-79 bare-alpha billing map and item-118 sterile mass are admitted.")
    print("  What remains is not a free density fit but three explicit residuals:")
    print("    1. the result inherits the R14/item-79 reconstruction floor for one")
    print("       bare-alpha non-unitary sterile firing;")
    print("    2. the absolute density inherits m_nuR=alpha0^2 Lambda_QCD;")
    print("    3. galaxy phenomenology must choose zero-mode CDM halos, or prove")
    print("       zero-mode depletion before retaining active R4/MOND.")
    print("  Therefore the book's risk is narrowed but not retired: the CMB third-peak")
    print("  budget is conditionally closed, while the full cosmological claim still")
    print("  lives or dies on the source/mass anchors and the halo branch.")
    print("exit 0 -- CMB cold-matter closure status: conditional closure, not locked.")


if __name__ == "__main__":
    main()
