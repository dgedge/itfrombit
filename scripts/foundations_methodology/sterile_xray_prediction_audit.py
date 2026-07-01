#!/usr/bin/env python3
r"""Sterile-nu_R X-ray target: line position, abundance, and flux ledger.

This script combines the item-118 line-position/mixing audit with the item-123
absolute-density chain.  It is deliberately tiered:

  * line energy: sharp conditional target from m_s = alpha0^2 Lambda_QCD;
  * abundance: conditional boot-QEC target from n_s/n_gamma = alpha0/208 and
    the R4 zero-mode 4:1 incidence, giving f_s = 1/5 of the dark density;
  * flux: a conditional theorem in the one-denominator R4-enforcement branch.
    The companion script item118_nuR_mixing_theorem.py proves the finite
    neutral-edge coefficient is exactly one and tightens the Schur sector to a
    single bright spectral moment K_B.  The remaining named premise is the
    service-normalised value K_B=1/v_R4.

The output is a pre-registration-friendly observable ledger, not a reduction of
X-ray data.
"""

from __future__ import annotations

import math


ALPHA0 = 1.0 / 137.035999
LAMBDA_QCD_GEV = 0.332
LAMBDA_QCD_EV = LAMBDA_QCD_GEV * 1e9
V_R4_EV = 246.0e9
T_CMB = 2.7255
ZETA3 = 1.2020569031595943
KB_EV_PER_K = 8.617333262e-5
HBARC_EV_CM = 1.973269804e-5
RHO_CRIT_H2_EV_CM3 = 1.05371e4
MSUN_G = 1.98847e33
EV_G = 1.78266192e-33
MPC_CM = 3.0856775814913673e24
KEV_ERG = 1.602176634e-9
T_HUBBLE_S = 4.35e17


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def photon_density_cm3() -> float:
    return 2.0 * ZETA3 / math.pi**2 * (KB_EV_PER_K * T_CMB / HBARC_EV_CM) ** 3


def decay_rate_sinv(sin2_2theta: float, mass_kev: float) -> float:
    # Pal-Wolfenstein convention used in the existing item-118 audit.
    return 1.38e-32 * (sin2_2theta / 1e-7) * mass_kev**5


def main() -> None:
    print("STERILE-nu_R X-RAY TARGET AUDIT")
    print("=" * 88)

    mass_ev = ALPHA0**2 * LAMBDA_QCD_EV
    mass_kev = mass_ev / 1e3
    line_kev = mass_kev / 2.0
    n_gamma = photon_density_cm3()
    n_over_ngamma = ALPHA0 / 208.0
    n_sterile = n_over_ngamma * n_gamma
    omega_sterile = n_sterile * mass_ev / RHO_CRIT_H2_EV_CM3
    omega_dark = 5.0 * omega_sterile
    f_sterile_dark = omega_sterile / omega_dark

    print("\n[1] Hard line-position and density targets")
    print(f"  m_s = alpha0^2 Lambda_QCD = {mass_kev:.4f} keV")
    print(f"  E_gamma = m_s/2           = {line_kev:.4f} keV")
    print(f"  n_s/n_gamma = alpha0/208  = {n_over_ngamma:.9e}")
    print(f"  n_s                       = {n_sterile:.6f} cm^-3")
    print(f"  Omega_s h^2               = {omega_sterile:.6f}")
    print(f"  Omega_dark h^2 = 5 Omega_s= {omega_dark:.6f}")
    print(f"  f_s = Omega_s/Omega_dark  = {f_sterile_dark:.3f}")
    check(17.5 < mass_kev < 17.9, "sterile mass is 17.7 keV-class")
    check(8.75 < line_kev < 8.95, "radiative line is fixed near 8.85 keV")
    check(abs(omega_sterile - 0.02418) / 0.02418 < 0.01, "sterile density matches the alpha0/208 chain")
    check(abs(f_sterile_dark - 0.2) < 1e-12, "sterile branch is exactly 20% of the paired dark ledger")

    print("\n[2] Mixing and decay-rate tiers")
    kappa_r4 = 1.0
    theta_r4 = mass_ev / V_R4_EV
    sin2_r4 = 4.0 * theta_r4**2
    gamma_r4 = decay_rate_sinv(sin2_r4, mass_kev)
    tau_r4 = 1.0 / gamma_r4
    print(f"  R4 Schur moment: kappa = v_R4 K_B = {kappa_r4:.1f}")
    print(f"  R4-scale theorem branch: theta = kappa m_s/v_R4 = {theta_r4:.3e}")
    print(f"  sin^2(2 theta) = {sin2_r4:.3e}")
    print(f"  Gamma = {gamma_r4:.3e} s^-1")
    print(f"  tau   = {tau_r4:.3e} s = {tau_r4/T_HUBBLE_S:.2e} Hubble times")
    check(sin2_r4 < 1e-10, "R4-scale mixing is below representative current X-ray bounds")
    check(tau_r4 / T_HUBBLE_S > 1e10, "sterile relic is cosmologically stable at the predicted mixing")

    print("\n[3] Observer-facing flux normalization")
    mass_g = mass_ev * EV_G
    particles_per_msun = MSUN_G / mass_g
    photons_per_s_per_msun_sterile = gamma_r4 * particles_per_msun
    photons_per_s_per_msun_dark = f_sterile_dark * photons_per_s_per_msun_sterile
    energy_erg = line_kev * KEV_ERG
    energy_lum_per_msun_sterile = photons_per_s_per_msun_sterile * energy_erg
    d_ref = 100.0 * MPC_CM
    m_ref_dark = 1.0e14
    flux_ref = photons_per_s_per_msun_dark * m_ref_dark / (4.0 * math.pi * d_ref**2)
    print(f"  photons/s per sterile Msun = {photons_per_s_per_msun_sterile:.3e}")
    print(f"  photons/s per total-dark Msun (f_s=0.2) = {photons_per_s_per_msun_dark:.3e}")
    print(f"  energy luminosity per sterile Msun = {energy_lum_per_msun_sterile:.3e} erg/s")
    print("  aperture-integrated local-target scaling:")
    print(f"    F_gamma = {flux_ref:.3e} ph cm^-2 s^-1")
    print("              * (M_dark / 1e14 Msun) * (100 Mpc / D)^2")
    print("              * (sin^2(2theta) / 2.07e-14)")
    print("              * (f_s / 0.2)")
    check(1e-12 < flux_ref < 1e-10, "reference cluster-scale flux has the expected order of magnitude")

    print("\n[4] Falsification logic")
    print("  Energy branch:")
    print("    A robust sterile-like line from this branch must be at 8.84 keV, redshifted")
    print("    by the source. A different sterile line energy does not support this branch.")
    print("  Flux branch:")
    print("    In the service-normalised R4 branch, K_B=1/v_R4 and theta=m_s/v_R4;")
    print("    non-detection below the flux normalization above falsifies that branch")
    print("    for the target. If the kappa=1 premise is rejected, non-detection")
    print("    only tightens the mixing-angle upper bound.")
    print("  Bright-line failure:")
    print("    A claimed 8.84 keV dark line implying sin^2(2theta) >> 2e-14 would refute")
    print("    the kappa=1 R4 mixing branch even if it preserved the mass identification.")

    print("\n[5] Verdict")
    print("  The line energy and abundance are strong conditional targets.  The flux is")
    print("  also a single-number target inside the one-denominator R4-enforcement")
    print("  branch: the finite register fixes the neutral edge and coefficient, and")
    print("  the Schur sector has only one remaining scalar, kappa=v_R4 K_B.")
    print("exit 0 -- sterile X-ray prediction ledger complete")


if __name__ == "__main__":
    main()
