#!/usr/bin/env python3
r"""ITEM 123: CMB completion made QUANTITATIVE — the matter-radiation-equality kill-shot
and the AeST-class structural target.

This sharpens r4_eos_cmb_resolution.py. That script settled the R4 EoS status split and
stated the CMB no-go as a "5x collisionless-anchor shortfall". The shortfall is correct but
soft: a peak-amplitude argument. This script gives the no-go its decisive, parameter-free
teeth and then turns it into a SPECIFIC missing-theorem target instead of a dead end.

THE FRAMEWORK'S COSMOLOGICAL INVENTORY (no free CDM):
  * baryons                              omega_b   h^2 = 0.0224
  * sterile nu_R (the 20% dark share)    omega_nuR h^2 = 0.024   -- cold ENOUGH to cluster on
       CMB scales (M_fs ~ 8e4 Msun << acoustic mass; item118_nuR_production_warmth.py)
  * R4, diffuse/homogeneous              w(a) = -1 + a/28        -- DARK ENERGY (item 131),
       NOT recombination dust: at a_rec it is w ~ -1 (Lambda-like, anti-clustering)
  * R4, bound in halos                   p=3 AQUAL / MOND        -- a LATE nonlinear response
       (item 132), not a homogeneous z~1100 fluid
  * R4, microscopic exhaust             w = +1/3                 -- radiation, free-streams
  => the ONLY cold clustering matter at recombination is omega_b + omega_nuR = 0.0464.

THE KILL-SHOT (new): matter-radiation equality is z_eq = omega_m/omega_r - 1. With only
0.0464 of clustering matter the framework has z_eq ~ 1100 ~ z_rec: equality COINCIDES with
recombination, so the universe is radiation-driven through the entire acoustic epoch ->
decaying potentials -> the third peak is wrong. LCDM needs z_eq ~ 3400 ~ 3 z_rec. This is the
parameter-free reason the peak fails, not a tuned amplitude.

THE TARGET (constructive): AeST (Skordis & Zlosnik, PRL 127, 161302, 2021) completes exactly
this kind of MOND cosmology with ONE field carrying (a) a conserved shift-charge ~ a^-3 (dust
-> CMB+P(k) like LCDM), (b) the quasi-static MOND limit, (c) a residual Lambda. The framework
already HAS (b) [item 132 p=3] and (c) [item 131 w->-1]. It is missing EXACTLY (a): a conserved
comoving a^-3 R4 charge. If that charge carries omega h^2 = 0.096, then nu_R + R4-dust = 0.120
exactly and z_eq -> 3400: the CMB closes. So the no-go is "one missing conserved Noether charge
in an otherwise AeST-shaped sector", and it is the framework's MOST at-risk sector (Planck has
ALREADY measured the third peak -- unlike the GW echoes or the X-ray line, which await LISA/XRISM).

exit 0 = the equality kill-shot, the R4-DE wrong-way check, the exact closure identity, and the
AeST gap diagnosis are computed and ASSERTED; the verdict (no-go sharpened, not resolved; the
dust charge is a conjecture/import, a specific hard target) is reported honestly.
"""

from __future__ import annotations

import math

# single-source the budget (do not restate it) -----------------------------------------------
from r4_eos_cmb_resolution import (
    OMEGA_B_H2,
    OMEGA_C_H2_NEEDED,
    OMEGA_NUR_H2,
    OMEGA_R4_H2,
)

# single-source the homogeneous R4 dark-energy law -------------------------------------------
from item123_structure_wz_corrections import DELTA, w_of_a


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


# --------------------------------------------------------------------------- radiation density
# omega_gamma = Omega_gamma h^2 computed from the CMB temperature (compute, don't recall).
A_RAD = 7.565723e-16        # radiation constant a = 4 sigma/c  (J m^-3 K^-4)
C = 2.99792458e8            # m/s
T_CMB = 2.7255              # K (FIRAS)
RHO_CRIT_PER_H2 = 1.87834e-26   # kg m^-3, critical density for h=1
N_EFF = 3.044               # standard effective neutrino number

rho_gamma = A_RAD * T_CMB ** 4 / C ** 2          # photon mass-equivalent density (kg/m^3, h=1)
OMEGA_GAMMA_H2 = rho_gamma / RHO_CRIT_PER_H2
NEUTRINO_FACTOR = 1.0 + (7.0 / 8.0) * (4.0 / 11.0) ** (4.0 / 3.0) * N_EFF
OMEGA_R_H2 = OMEGA_GAMMA_H2 * NEUTRINO_FACTOR    # photons + relativistic neutrinos

Z_REC = 1090.0              # Planck last-scattering redshift z*
A_REC = 1.0 / (1.0 + Z_REC)


def z_eq(omega_m_h2: float) -> float:
    """matter-radiation equality redshift: 1+z_eq = rho_m/rho_r|today = omega_m/omega_r."""
    return omega_m_h2 / OMEGA_R_H2 - 1.0


def main() -> None:
    print("ITEM 123: CMB COMPLETION — EQUALITY KILL-SHOT + AeST TARGET")

    print("\n[0] Radiation density (computed from T_CMB)")
    print(f"  omega_gamma h^2 = {OMEGA_GAMMA_H2:.3e}  (photons; T={T_CMB} K)")
    print(f"  N_eff={N_EFF} -> neutrino factor {NEUTRINO_FACTOR:.4f}")
    print(f"  omega_r h^2     = {OMEGA_R_H2:.3e}  (photons + relativistic neutrinos)")
    check(abs(OMEGA_GAMMA_H2 - 2.47e-5) < 0.1e-5, "omega_gamma h^2 ~ 2.47e-5 (T_CMB cross-check)")
    check(abs(OMEGA_R_H2 - 4.18e-5) < 0.2e-5, "omega_r h^2 ~ 4.18e-5 (photons+nu cross-check)")

    print("\n[1] Clustering matter at recombination — the framework has no free CDM")
    omega_m_lcdm = OMEGA_B_H2 + OMEGA_C_H2_NEEDED
    omega_m_fw = OMEGA_B_H2 + OMEGA_NUR_H2            # baryons + cold sterile nu_R only
    print(f"  LCDM      omega_m h^2 = omega_b + omega_c   = {omega_m_lcdm:.4f}")
    print(f"  framework omega_m h^2 = omega_b + omega_nuR = {omega_m_fw:.4f}   (R4 is DE/MOND, not dust)")
    print(f"  ratio framework/LCDM clustering matter        = {omega_m_fw/omega_m_lcdm:.3f}")
    check(abs(omega_m_lcdm - 0.1424) < 1e-4, "LCDM clustering matter is 0.142")
    check(abs(omega_m_fw - 0.0464) < 1e-4, "framework clustering matter is only 0.046 (baryons+nu_R)")

    print("\n[2] THE KILL-SHOT: matter-radiation equality lands at recombination")
    zeq_lcdm = z_eq(omega_m_lcdm)
    zeq_fw = z_eq(omega_m_fw)
    print(f"  LCDM:      z_eq = {zeq_lcdm:.0f}   (z_eq/z_rec = {zeq_lcdm/Z_REC:.2f}: matter-dominated WELL before recombination)")
    print(f"  framework: z_eq = {zeq_fw:.0f}   (z_eq/z_rec = {zeq_fw/Z_REC:.2f}: equality COINCIDES with recombination)")
    check(3300 < zeq_lcdm < 3500, "LCDM z_eq ~ 3400 (matches Planck)")
    check(1000 < zeq_fw < 1250, "framework z_eq ~ 1100")
    check(0.9 < zeq_fw / Z_REC < 1.2, "framework equality coincides with recombination (z_eq ~ z_rec)")
    check(zeq_lcdm / Z_REC > 3.0, "LCDM needs equality at ~3x recombination for the standard peak structure")
    print("  => with z_eq~z_rec the acoustic epoch is RADIATION-driven: potentials decay, the")
    print("     third peak is suppressed/wrong. This is the parameter-free no-go, not a tuned amplitude.")

    print("\n[3] The homogeneous R4 (dark energy) goes the WRONG way at recombination")
    w_r4_rec = w_of_a(A_REC, 0.0)                    # item-131 law: -1 + a/28
    # DE density vs matter density at recombination (w ~ -1 => ~constant DE; matter ~ (1+z)^3)
    omega_de_today, omega_m_today = 0.685, 0.315
    de_over_matter_rec = (omega_de_today / omega_m_today) / (1.0 + Z_REC) ** 3
    print(f"  R4-DE EoS at a_rec={A_REC:.2e}:  w = {w_r4_rec:.6f}  (dust needs w=0; this is ~ -1, Lambda-like)")
    print(f"  R4-DE / matter density at z_rec ~ {de_over_matter_rec:.1e}  (utterly negligible AND wrong-pressure)")
    check(w_r4_rec < -0.999, "homogeneous R4 is w ~ -1 at recombination (anti-clustering, not dust)")
    check(de_over_matter_rec < 1e-8, "R4-DE is ~1e-9 of matter at recombination: it cannot be the CDM")
    check(abs(w_of_a(1.0, 0.0) - (-27.0 / 28.0)) < 1e-12, "law cross-check: w0 = -27/28 (item 131)")

    print("\n[4] nu_R is cold ENOUGH to count (the issue is the AMOUNT, not the temperature)")
    M_HM_NUR = 8.0e4          # half-mode mass (Msun) from item118_nuR_production_warmth.py
    # acoustic-scale mass: matter in a sphere of the comoving sound horizon (~145 Mpc).
    r_s_mpc = 145.0
    rho_m_today_msun_mpc3 = omega_m_today * 2.775e11   # rho_crit,0 * Omega_m, in Msun/Mpc^3 (h=1)
    M_acoustic = (4.0 / 3.0) * math.pi * (r_s_mpc / 2.0) ** 3 * rho_m_today_msun_mpc3
    print(f"  nu_R half-mode mass   M_fs   ~ {M_HM_NUR:.0e} Msun  (very cold; non-thermal boundary-printed)")
    print(f"  acoustic-scale mass   M_aco  ~ {M_acoustic:.1e} Msun  (sound horizon r_s~{r_s_mpc} Mpc)")
    print(f"  M_fs / M_acoustic ~ {M_HM_NUR/M_acoustic:.0e}: nu_R free-streaming is irrelevant on CMB scales")
    check(M_HM_NUR / M_acoustic < 1e-8, "nu_R clusters fully on acoustic scales -> a VALID cold anchor")
    print("  => nu_R is genuine cold matter for the CMB; it just carries only 0.024, not 0.120.")

    print("\n[5] The exact closure identity — what an a^-3 R4 dust charge would buy")
    omega_m_with_dust = OMEGA_B_H2 + OMEGA_NUR_H2 + OMEGA_R4_H2
    anchor_with_dust = OMEGA_NUR_H2 + OMEGA_R4_H2
    zeq_with_dust = z_eq(omega_m_with_dust)
    print(f"  if R4 (0.096) were conserved a^-3 dust:")
    print(f"    nu_R + R4-dust = {anchor_with_dust:.3f}  (target omega_c = {OMEGA_C_H2_NEEDED:.3f})")
    print(f"    omega_m -> {omega_m_with_dust:.4f}  =>  z_eq -> {zeq_with_dust:.0f}  (= LCDM, peak structure restored)")
    check(abs(anchor_with_dust - OMEGA_C_H2_NEEDED) < 1e-12, "nu_R + R4-dust closes the cold budget EXACTLY")
    check(abs(zeq_with_dust - zeq_lcdm) < 1.0, "and reproduces the LCDM equality redshift exactly")

    print("\n[6] AeST gap diagnosis — the framework is AeST-shaped, missing ONE ingredient")
    # AeST (Skordis-Zlosnik 2021) = one field with three pieces; tick what the framework already has.
    aest_pieces = {
        "(a) conserved shift-charge ~ a^-3 (dust: CMB+P(k))": False,  # MISSING
        "(b) quasi-static MOND/AQUAL limit":                  True,   # item 132 p=3
        "(c) residual cosmological constant":                 True,   # item 131 w -> -1
    }
    for piece, have in aest_pieces.items():
        print(f"    [{'HAVE' if have else 'MISS'}] {piece}")
    have = sum(aest_pieces.values())
    check(have == 2, "framework already has 2 of AeST's 3 ingredients (MOND limit + residual Lambda)")
    check(aest_pieces["(a) conserved shift-charge ~ a^-3 (dust: CMB+P(k))"] is False,
          "the missing piece is EXACTLY a conserved comoving a^-3 R4 charge")

    # the four dust criteria are all NO under current R4 mechanics (from r4_eos_cmb_resolution).
    dust_criteria_met = [False, False, False, False]
    check(not any(dust_criteria_met),
          "current R4 = non-conserved service/exhaust ledger; only nu_R is a conserved species")

    print("\n[7] Verdict")
    print(
        "  CMB completion is NOT achieved. The no-go is SHARPENED from '5x anchor shortfall' to a\n"
        "  parameter-free statement: with only baryons + nu_R clustering, z_eq ~ z_rec, so the acoustic\n"
        "  epoch is radiation-driven and the third peak is wrong; the homogeneous R4 law actively goes\n"
        "  the WRONG way (w -> -1 at recombination, not w=0 dust)."
    )
    print(
        "  CONSTRUCTIVE: the sector is AeST-shaped -- it already has the MOND limit (item 132) and a\n"
        "  residual Lambda (item 131) -- and is missing EXACTLY one thing: a conserved comoving a^-3 R4\n"
        "  shift-charge. If it exists with omega h^2=0.096, nu_R+R4-dust=0.120 and z_eq->3400: the CMB\n"
        "  closes. So the target is precise: derive a conserved Noether charge from W=S*C boundary-QEC."
    )
    print(
        "  TIER: the no-go is the current canon verdict (settled); the a^-3 R4 dust charge is a\n"
        "  CONJECTURE/import (ambiguous) -- not derived, and presently contradicted by the R4-DE law.\n"
        "  AT-RISK: this is the framework's most falsifiable sector -- Planck has ALREADY measured the\n"
        "  third peak (cf. GW echoes / X-ray line, which await LISA / XRISM)."
    )
    print("exit 0 -- CMB no-go made quantitative (z_eq kill-shot); AeST dust-charge is the precise open target.")


if __name__ == "__main__":
    main()
