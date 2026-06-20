#!/usr/bin/env python3
r"""Trans-Lambda_QCD null-chain precision residual map.

K31/K32 close the support/endpoint question at framed-causal-set smooth-field
EFT grade:

  * a TeV photon is not a Bloch mode of the QCD-spaced IR crystal;
  * a finer oscillator support would reopen the cosmological-constant ledger;
  * the viable high-energy object is one normalized null-chain external leg
    carrying total P^mu, with endpoint LSZ residue Z=1.

This script addresses the narrower residuals that remain after that closure.
It does not claim to be a gamma-ray propagation code or a full QED event
generator.  It gives the durable bookkeeping:

  1. precision QED phenomenology is standard normalized-leg EFT work;
  2. finite-density effects are plasma / opacity self-energies, not support
     failures;
  3. Lorentz-violation tests become null bounds on residual dispersion
     coefficients of the null-chain leg.  Any O(1) per-service-link dispersion
     is excluded; the canon default is exact null propagation, i.e. coefficient
     zero.
"""

from __future__ import annotations

import math


ALPHA0_INV = 137.035999084
ALPHA0 = 1.0 / ALPHA0_INV
LAMBDA_GEV = 0.3317
TEV_GEV = 1000.0
HBARC_GEV_FM = 0.1973269804
ME_EV = 510_998.95
CM_PER_GPC = 3.0856775814913673e27
C_CM_S = 2.99792458e10


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def support_numbers(e_gev: float = TEV_GEV) -> dict[str, float]:
    e_bz = math.pi * LAMBDA_GEV
    lambda_gamma = e_gev / math.pi
    return {
        "e_bz": e_bz,
        "e_over_bz": e_gev / e_bz,
        "lambda_gamma": lambda_gamma,
        "lambda_gamma_over_lambda": lambda_gamma / LAMBDA_GEV,
        "rho_inflation": (lambda_gamma / LAMBDA_GEV) ** 4,
        "service_units": math.ceil(e_gev / LAMBDA_GEV),
    }


def running_alpha_one_loop(q_gev: float) -> tuple[float, int]:
    """One-loop threshold orientation, not a precision electroweak fit."""

    fermions = [
        (0.00051099895, -1.0, 1),
        (0.1056583755, -1.0, 1),
        (1.77686, -1.0, 1),
        (0.00216, 2.0 / 3.0, 3),
        (0.00467, -1.0 / 3.0, 3),
        (0.093, -1.0 / 3.0, 3),
        (1.27, 2.0 / 3.0, 3),
        (4.18, -1.0 / 3.0, 3),
        (172.76, 2.0 / 3.0, 3),
    ]
    sigma = 0.0
    active = 0
    for mass, charge, nc in fermions:
        if q_gev > mass:
            sigma += nc * charge * charge * math.log(q_gev / mass)
            active += 1
    alpha_inv_q = ALPHA0_INV - (2.0 / (3.0 * math.pi)) * sigma
    return alpha_inv_q, active


def finite_density(e_gev: float = TEV_GEV, n_e_cm3: float = 1.0e-7) -> dict[str, float]:
    e_ev = e_gev * 1.0e9
    omega_p_ev = 3.713e-11 * math.sqrt(n_e_cm3)
    plasma_dv = -0.5 * (omega_p_ev / e_ev) ** 2
    epsilon_thr_ev = ME_EV * ME_EV / e_ev
    cmb_epsilon_ev = 2.701 * 8.617333262e-5 * 2.7255
    e_thr_cmb_tev = (ME_EV * ME_EV / cmb_epsilon_ev) / 1.0e12
    return {
        "omega_p_ev": omega_p_ev,
        "plasma_dv": plasma_dv,
        "epsilon_thr_ev": epsilon_thr_ev,
        "cmb_epsilon_ev": cmb_epsilon_ev,
        "e_thr_cmb_tev": e_thr_cmb_tev,
    }


def tof_bound_zeta(
    e_gev: float,
    distance_gpc: float,
    delta_t_s: float,
    power: int,
    scale_gev: float = LAMBDA_GEV,
) -> float:
    r"""Bound zeta_n in |Delta v|/c = |zeta_n| (E/scale)^n.

    This is a benchmark sensitivity calculator, not a replacement for a current
    SME/global likelihood.  It is useful because it shows how severe the null
    prediction must be if the scale is Lambda_QCD rather than M_Pl.
    """

    t_prop = distance_gpc * CM_PER_GPC / C_CM_S
    return (delta_t_s / t_prop) / ((e_gev / scale_gev) ** power)


def ir_photon_anisotropy(e_gev: float) -> float:
    """SC Maxwell envelope artifact Delta v/v ~ (E/Lambda)^2/12."""

    return (e_gev / LAMBDA_GEV) ** 2 / 12.0


def main() -> None:
    print("TRANS-LAMBDA_QCD PRECISION RESIDUAL MAP")
    print("=" * 86)

    print("\n[1] Support/LSZ status recap")
    s = support_numbers()
    print(f"  QCD-crystal Brillouin ceiling pi Lambda = {s['e_bz']:.3f} GeV")
    print(f"  1 TeV / ceiling                         = {s['e_over_bz']:.1f}")
    print(f"  finer oscillator scale for 1 TeV support = {s['lambda_gamma']:.1f} GeV")
    print(f"  vacuum-density inflation if oscillatorized = {s['rho_inflation']:.3e}")
    print(f"  null-chain service units at 1 TeV          = {s['service_units']}")
    check(s["e_over_bz"] > 900.0, "TeV photons are not QCD-crystal Bloch modes")
    check(s["rho_inflation"] > 1.0e11, "fine oscillator support would reopen the CC ledger")

    print("\n[2] Precision QED layer")
    for q in (91.1876, 1000.0, 100_000.0):
        alpha_inv, active = running_alpha_one_loop(q)
        print(f"  orientation: alpha^-1({q:g} GeV) ~ {alpha_inv:.3f} with {active} active thresholds")
    print("  meaning: process-level predictions are standard normalized-leg EFT tasks")
    print("           (higher loops, detector acceptances, Compton/pair kernels), not")
    print("           new microscopic support or LSZ-normalization gates.")

    print("\n[3] Finite-density / propagation layer")
    fd = finite_density()
    print(f"  IGM hbar omega_p (n_e=1e-7 cm^-3) = {fd['omega_p_ev']:.3e} eV")
    print(f"  1 TeV plasma Delta v/v             = {fd['plasma_dv']:.3e}")
    print(f"  pair target threshold at 1 TeV      = {fd['epsilon_thr_ev']:.3f} eV")
    print(f"  CMB pair threshold                  = {fd['e_thr_cmb_tev']:.1f} TeV")
    check(abs(fd["plasma_dv"]) < 1.0e-45, "IGM plasma speed shift is negligible for TeV timing")
    check(0.1 < fd["epsilon_thr_ev"] < 10.0, "TeV opacity is an EBL/IR transfer problem")

    print("\n[4] Lorentz-violation null bounds as coefficient gates")
    for e_gev, label in ((1.0e3, "1 TeV"), (1.0e5, "100 TeV")):
        b1 = tof_bound_zeta(e_gev, distance_gpc=1.0, delta_t_s=1.0, power=1)
        b2 = tof_bound_zeta(e_gev, distance_gpc=1.0, delta_t_s=1.0, power=2)
        print(f"  benchmark {label}, 1 Gpc, 1 s:")
        print(f"    |zeta_1| in Delta v/c=zeta_1 E/Lambda  < {b1:.3e}")
        print(f"    |zeta_2| in Delta v/c=zeta_2 (E/Lambda)^2 < {b2:.3e}")
        check(b1 < 1.0e-18, f"{label}: O(1) linear per-service LV is excluded")
        check(b2 < 1.0e-20, f"{label}: O(1) quadratic per-service LV is excluded")

    optical = ir_photon_anisotropy(2.0e-9)
    tev_if_misread = ir_photon_anisotropy(1.0e3)
    print("\n[5] Category check against the IR photon lattice artifact")
    print(f"  optical SC-envelope anisotropy  = {optical:.3e}")
    print(f"  TeV SC-envelope misread would be = {tev_if_misread:.3e} (nonsense: outside support)")
    check(optical < 1.0e-17, "IR photon envelope remains below the optical SME-scale line")
    check(tev_if_misread > 1.0e5, "applying the IR lattice artifact at TeV is a category error")

    print(
        r"""
[6] VERDICT
  The trans-Lambda_QCD residual is no longer "can TeV photons exist?" or
  "what is the LSZ normalization?"  At current canon grade those are closed by
  the framed causal-set/null-chain endpoint.

  What remains is a phenomenology programme:
    P1. precision QED processes with one normalized external null-chain leg;
    P2. astrophysical transfer functions (EBL opacity, plasma, detector response);
    P3. Lorentz-violation null tests, parameterized as residual zeta_n.  Because
        the relevant unit is Lambda_QCD, even tiny measured delays would force
        zeta_n to be absurdly small.  Therefore the canon prediction is not an
        O(1) Planck-suppressed dispersion; it is exact null-chain propagation
        plus ordinary medium effects.  Any nonzero intrinsic per-service
        dispersion must be derived and then confronted with the coefficient
        gates above.

  Thus: support/LSZ closed; precision phenomenology and LV bounds open as
  ordinary external calculations around the closed endpoint representation.
exit 0 -- residual map complete."""
    )


if __name__ == "__main__":
    main()
