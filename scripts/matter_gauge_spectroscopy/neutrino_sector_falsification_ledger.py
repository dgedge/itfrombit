#!/usr/bin/env python3
r"""Neutrino-sector falsification ledger.

This script gathers the currently canonical neutrino-sector readouts into one
small numerical table:

* active-neutrino branch: normal ordering, light-Majorana 0nu beta beta band;
* sterile branch: 17.7 keV nu_R line position and gauge-forbidden mixing scale;
* direct-lab branch: KATRIN/TRISTAN is far above the predicted mixing;
* baryogenesis branch: alpha^4 scale plus record-alphabet correction.

It is intentionally not a new derivation.  The purpose is to expose the compact
experimental cluster and to remove one ambiguous wording trap: canon is not
"Dirac-like rather than Majorana."  The active neutrino is Majorana in the code
identity statement; the sterile nu_R is a separate pseudocodeword branch with
gauge-forbidden Dirac mixing.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


ALPHA = 1.0 / 137.035999
LAMBDA_QCD_EV = 0.332e9
V_R4_EV = 246.0e9

# Canonical item-87/Part-5 normal-ordering masses after the Delta m31 anchor.
M1_EV = 0.00079
M2_EV = 0.00871
M3_EV = 0.05012

# Representative current PMNS first-row weights.  The exact phases are not
# specified by canon, so m_bb is quoted as the phase envelope from these three
# positive terms.
SIN2_THETA12 = 0.304
SIN2_THETA13 = 0.0222

ETA_OBS = 6.120e-10
S_PER_PHOTON = 7.04
S1_RECORD = math.log(8 * 137)


@dataclass(frozen=True)
class Gate:
    observable: str
    canonical: str
    test: str
    status: str


def active_mass_readouts() -> dict[str, float]:
    c12 = 1.0 - SIN2_THETA12
    c13 = 1.0 - SIN2_THETA13

    terms = (
        M1_EV * c12 * c13,
        M2_EV * SIN2_THETA12 * c13,
        M3_EV * SIN2_THETA13,
    )
    mbb_max = sum(terms)
    mbb_min = max(0.0, max(terms) - (sum(terms) - max(terms)))
    mbeta = math.sqrt(
        M1_EV * M1_EV * c12 * c13
        + M2_EV * M2_EV * SIN2_THETA12 * c13
        + M3_EV * M3_EV * SIN2_THETA13
    )
    dm21 = M2_EV * M2_EV - M1_EV * M1_EV
    dm31 = M3_EV * M3_EV - M1_EV * M1_EV
    return {
        "sum_m": M1_EV + M2_EV + M3_EV,
        "m_beta": mbeta,
        "mbb_min": mbb_min,
        "mbb_max": mbb_max,
        "r_dm": dm21 / dm31,
    }


def sterile_readouts() -> dict[str, float]:
    m_s_ev = ALPHA * ALPHA * LAMBDA_QCD_EV
    theta = m_s_ev / V_R4_EV
    return {
        "m_s_ev": m_s_ev,
        "line_kev": m_s_ev / 2.0e3,
        "theta": theta,
        "ue4_sq": theta * theta,
        "sin2_2theta": 4.0 * theta * theta,
        "m_dirac_ev": theta * m_s_ev,
    }


def baryogenesis_readouts() -> dict[str, float]:
    eta_leading = (3.0 / 14.0) * ALPHA**4
    eta_record = eta_leading * S_PER_PHOTON / S1_RECORD
    return {
        "eta_alpha4": eta_leading,
        "eta_record": eta_record,
        "eta_record_over_obs": eta_record / ETA_OBS,
        "s1_record": S1_RECORD,
    }


def print_gate_table(gates: list[Gate]) -> None:
    print("\n[4] COMPACT FALSIFICATION TABLE")
    print("-" * 96)
    print(f"{'observable':<26} {'canonical readout':<33} {'test':<25} status")
    print("-" * 96)
    for g in gates:
        print(f"{g.observable:<26} {g.canonical:<33} {g.test:<25} {g.status}")


def main() -> None:
    active = active_mass_readouts()
    sterile = sterile_readouts()
    baryo = baryogenesis_readouts()

    print("NEUTRINO-SECTOR FALSIFICATION LEDGER")
    print("=" * 96)
    print("[1] active neutrino branch: normal ordering + Majorana identity")
    print(f"    m1,m2,m3        = {M1_EV*1e3:.2f}, {M2_EV*1e3:.2f}, {M3_EV*1e3:.2f} meV")
    print(f"    Sigma m_nu      = {active['sum_m']*1e3:.2f} meV")
    print(f"    Delta m21^2/m31^2 = {active['r_dm']:.5f}")
    print(f"    m_beta          = {active['m_beta']*1e3:.2f} meV")
    print(
        "    m_bb envelope   = "
        f"{active['mbb_min']*1e3:.2f}..{active['mbb_max']*1e3:.2f} meV "
        "(light-Majorana exchange, phases free)"
    )
    print("    interpretation  = active neutrino is Majorana in canon; sterile nu_R is a separate branch.")

    print("\n[2] sterile nu_R branch: line position and mixing")
    print(f"    m_s             = {sterile['m_s_ev']/1e3:.2f} keV")
    print(f"    E_gamma         = {sterile['line_kev']:.2f} keV")
    print(f"    theta=m_s/v_R4  = {sterile['theta']:.3e}")
    print(f"    |U_e4|^2        = {sterile['ue4_sq']:.3e}")
    print(f"    sin^2(2theta)   = {sterile['sin2_2theta']:.3e}")
    print(f"    m_D             = {sterile['m_dirac_ev']*1e3:.3f} meV")
    print("    interpretation  = non-detection expected; a bright 8.84 keV line would falsify this mixing reading.")

    print("\n[3] baryogenesis branch: alpha^4 scale and record-alphabet correction")
    print(f"    (3/14) alpha^4          = {baryo['eta_alpha4']:.4e}")
    print(f"    s1=ln(8*137)            = {baryo['s1_record']:.6f}")
    print(f"    eta_record              = {baryo['eta_record']:.4e}")
    print(f"    eta_record / eta_obs    = {baryo['eta_record_over_obs']:.5f}")
    print("    interpretation          = numerically sharp, but still conditional on the event-current theorem.")

    gates = [
        Gate(
            "ordering",
            f"NH, Sigma={active['sum_m']*1e3:.1f} meV",
            "oscillation/cosmology",
            "sharp branch",
        ),
        Gate(
            "0nu beta beta",
            f"m_bb={active['mbb_min']*1e3:.1f}-{active['mbb_max']*1e3:.1f} meV",
            "LEGEND/nEXO/CUPID",
            "allowed; likely below near-term reach",
        ),
        Gate(
            "sterile X-ray",
            f"E_gamma={sterile['line_kev']:.2f} keV",
            "XRISM/NuSTAR/Athena",
            "line position sharp; flux tiny",
        ),
        Gate(
            "sterile mixing",
            f"sin2 2th={sterile['sin2_2theta']:.1e}",
            "X-ray/direct beta",
            "bright detection >1e-12 fails",
        ),
        Gate(
            "direct keV sterile",
            f"|Ue4|^2={sterile['ue4_sq']:.1e}",
            "KATRIN/TRISTAN",
            "well below 1e-6 reach",
        ),
        Gate(
            "baryogenesis",
            f"eta={baryo['eta_record']:.3e}",
            "event-ledger proof",
            "data matched; theorem still gated",
        ),
        Gate(
            "galaxy warmth",
            "17.7 keV nu_R is cold",
            "RAR/galaxy halos",
            "tension if sizable in galaxies",
        ),
    ]
    print_gate_table(gates)

    assert abs(sterile["m_s_ev"] / 1e3 - 17.7) < 0.2
    assert 1.0e-15 < sterile["sin2_2theta"] < 1.0e-13
    assert active["sum_m"] < 0.061
    assert active["mbb_max"] < 0.005
    assert 0.99 < baryo["eta_record_over_obs"] < 1.01
    print("\nexit 0 -- neutrino-sector ledger is internally consistent.")


if __name__ == "__main__":
    main()
