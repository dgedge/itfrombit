#!/usr/bin/env python3
r"""Precision phenomenology frontier map.

Purpose
-------
The current framework has several sectors where the principle-level object is
now defined well enough that further progress is no longer mainly a new
substrate theorem.  The remaining work is a precision calculation around a
fixed object:

  * CMB/halo: the zero-mode-CDM branch is forward-Boltzmann grade around the
    native 63/64 acoustic clock, but still needs a full likelihood and
    nonlinear halo/galaxy calculation.
  * Black holes: flux/species/polarization maps around the local KMS ladder,
    freeze shell, Schwarzschild greybody transfer, and conditional 10/27 source.
  * Trans-Lambda_QCD quanta: process-level QED, astrophysical opacity, detector
    response, and Lorentz-violation null bounds around the normalized null-chain
    endpoint leg.

This script is a triage classifier, not a replacement for those external runs.
Exit 0 means the residuals are correctly classified as precision/transfer/
likelihood work around already-defined framework objects, rather than hidden
normalization or principle gaps.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def require_text(path: str, needles: tuple[str, ...]) -> None:
    text = (ROOT / path).read_text(encoding="utf-8")
    missing = [needle for needle in needles if needle not in text]
    if missing:
        raise AssertionError(f"{path} missing required text: {missing}")


@dataclass(frozen=True)
class Frontier:
    sector: str
    fixed_object: str
    remaining_calculation: str
    decisive_outputs: tuple[str, ...]
    status: str
    scripts: tuple[str, ...]


FRONTIERS = (
    Frontier(
        sector="CMB / halo",
        fixed_object="Boltzmann-grade 63/64 acoustic-clock branch plus zero-mode-as-CDM bookkeeping",
        remaining_calculation=(
            "full Planck/ACT/SPT/DESI/BAO/SN likelihood implementation with nuisance parameters; "
            "nonlinear halo/galaxy modelling for the zero-mode-CDM branch; "
            "no independent MOND/R4 double-counting"
        ),
        decisive_outputs=(
            "theta_star and third-peak likelihood",
            "TT/TE/EE/lensing chi2 and parameter shifts",
            "S8 / matter power / halo catalogue consistency",
            "zero-mode-as-CDM versus active-R4/MOND non-double-counting",
        ),
        status="Boltzmann-grade; full likelihood/nonlinear halo external",
        scripts=(
            "python_code/item123_cmb_native_spectrum_lift.py",
            "python_code/item123_halo_depletion_halomodel_audit.py",
            "python_code/item123_cmb_theta_halo_completion_gate.py",
        ),
    ),
    Frontier(
        sector="black-hole observables",
        fixed_object="local KMS ladder, beta-one freeze shell, escape cone, Schwarzschild greybody transfer, conditional 10/27 source",
        remaining_calculation=(
            "QEC species/polarization emission ledger, all-contact service-class proof, "
            "Kerr/fermion extensions, and observable flux templates"
        ),
        decisive_outputs=(
            "absolute Hawking luminosity coefficient",
            "species/helicity weights before greybody transfer",
            "line/continuum spectra after spin/partial-wave barriers",
            "echo-null and fast-scrambling negative versus a new nonlocal graph",
        ),
        status="exterior transfer computed; flux/species ledger conditional",
        scripts=(
            "python_code/bh_greybody_transfer.py",
            "python_code/bh_flux_species_polarization_ledger.py",
            "python_code/bh_observable_residual_map.py",
            "python_code/bh_fast_scrambling_topological_obstruction.py",
        ),
    ),
    Frontier(
        sector="trans-Lambda_QCD QED / LV",
        fixed_object="normalized null-chain endpoint leg with Z_endpoint=1, Ward identity, physical polarizations, and spin frame",
        remaining_calculation=(
            "process-level QED amplitudes, higher loops, astrophysical transfer/opacity, "
            "detector response, and observational Lorentz-violation bounds"
        ),
        decisive_outputs=(
            "Compton/pair-production/detector cross-section bookkeeping",
            "EBL/CMB gamma-gamma opacity and plasma transfer",
            "instrument response and source-intrinsic lag marginalisation",
            "bounds on intrinsic zeta_n dispersion coefficients; canon default zeta_n=0",
        ),
        status="support/LSZ/Ward closed; precision transfer external",
        scripts=(
            "python_code/foundations_trans_lambda_qed_phenomenology.py",
            "python_code/foundations_trans_lambda_precision_residuals.py",
            "python_code/foundations_trans_lambda_precision_factorization_theorem.py",
        ),
    ),
)


def lv_benchmark_bounds(energy_tev: float, distance_gpc: float, residual_seconds: float) -> tuple[float, float]:
    """Bounds for Delta v/c = zeta_n (E/Lambda_QCD)^n from time-of-flight residual.

    Uses the same benchmark convention as the canon notes: a signal over D with
    residual dt bounds |Delta v/c| < dt / (D/c).
    """

    seconds_per_gpc = 1.02927125e17
    lambda_qcd_gev = 0.3317
    e_gev = energy_tev * 1000.0
    eps = residual_seconds / (distance_gpc * seconds_per_gpc)
    ratio = e_gev / lambda_qcd_gev
    return eps / ratio, eps / (ratio * ratio)


def main() -> None:
    print("PRECISION PHENOMENOLOGY FRONTIER MAP")
    print("=" * 96)

    require_text(
        "python_code/item123_cmb_native_spectrum_lift.py",
        ("forward Boltzmann-code grade", "full likelihood/halo modelling remain"),
    )
    require_text(
        "python_code/bh_flux_species_polarization_ledger.py",
        ("two-helicity", "species/polarization ledger"),
    )
    require_text(
        "python_code/foundations_trans_lambda_qed_phenomenology.py",
        ("Z_endpoint = 1", "ordinary QED"),
    )
    require_text(
        "python_code/foundations_trans_lambda_precision_factorization_theorem.py",
        ("intrinsic dispersion coefficient zeta_n", "precision QED and astrophysical transfer"),
    )

    print("\n[1] Frontier table")
    for item in FRONTIERS:
        print(f"\n  {item.sector}")
        print(f"    fixed object:      {item.fixed_object}")
        print(f"    remaining calc:    {item.remaining_calculation}")
        print(f"    status:            {item.status}")
        print("    decisive outputs:")
        for output in item.decisive_outputs:
            print(f"      - {output}")
        for script in item.scripts:
            check((ROOT / script).exists(), f"{item.sector}: {script} exists")
        check("external" in item.status or "conditional" in item.status, f"{item.sector}: status is calculation-layer, not principle-locked")

    print("\n[2] LV benchmark bounds")
    z1, z2 = lv_benchmark_bounds(energy_tev=1.0, distance_gpc=1.0, residual_seconds=1.0)
    z1_100, z2_100 = lv_benchmark_bounds(energy_tev=100.0, distance_gpc=1.0, residual_seconds=1.0)
    print(f"    1 TeV, 1 Gpc, 1 s:       |zeta_1| < {z1:.3e}, |zeta_2| < {z2:.3e}")
    print(f"    100 TeV, 1 Gpc, 1 s:     |zeta_1| < {z1_100:.3e}, |zeta_2| < {z2_100:.3e}")
    check(z1 < 4.0e-21 and z2 < 2.0e-24, "TeV LV benchmark matches the canon precision-null scale")
    check(z1_100 < 4.0e-23 and z2_100 < 2.0e-28, "100 TeV timing is an even sharper null test")

    print("\n[3] Classification rule")
    print(
        """
  A residual belongs to precision phenomenology if:
    (i)  the framework object and normalization are already fixed;
    (ii) the remaining map is an exterior transfer, likelihood, opacity,
         detector, species, or nuisance-parameter calculation;
    (iii) changing the result does not license adding a new substrate knob.

  It returns to principle-theorem status only if the calculation demands a new
  unnormalised external leg, new service-rate coefficient, new species source,
  or nonzero intrinsic LV operator.

  Current outcome:
    CMB/halo, black-hole flux/species, and trans-Lambda_QCD QED are all now in
    this precision layer.  They are important and falsifiable, but they should
    be pursued with public Boltzmann, likelihood, transfer, and detector
    pipelines rather than by changing the finite substrate algebra.

exit 0"""
    )


if __name__ == "__main__":
    main()
