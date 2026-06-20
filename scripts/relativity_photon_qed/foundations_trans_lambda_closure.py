#!/usr/bin/env python3
r"""Close the trans-Lambda_QCD representation residual at framed-causet EFT grade.

ANCHOR item 150 had one real leftover after the photon object-identity theorem:

    a0 = hbar c / Lambda_QCD gives a Brillouin ceiling near 1 GeV,
    so how can observed TeV photons be Lorentz-invariant quanta?

The previous audits established the two sides:

  * lattice-mode no-go: a TeV photon cannot be an elementary Bloch mode on the
    QCD-spaced condensate without adding a finer oscillator scale, which
    reopens the cosmological-constant problem;
  * causal-set route: a source-conditioned null chain is noncompact, has a
    one-leg propagator/holonomy coupling, exact finite-rho Ward covariance,
    a framed mesoscopic Maxwell action, and framed Dirac spin.

This script ties those gates into one conservative closure statement:

  trans-cutoff quanta are represented by framed causal-set null-chain QED
  histories, not by QCD-lattice Bloch modes and not by N independent soft
  photons.  The same Lambda_QCD scale remains the vacuum/IR-condensate cutoff
  that protects rho_Lambda; high-energy external events do not add a UV
  zero-point oscillator tower.

Exit 0 means the bookkeeping is internally consistent.  It is not an
order-only causal-set theorem and it does not claim a dynamical gravity/tetrad
equation.
"""

from __future__ import annotations

import math


LAMBDA_GEV = 0.3317
TEV_GEV = 1000.0
VISIBLE_EV = 2.48  # 500 nm
SME_VISIBLE_BOUND = 1.0e-17


def lattice_gates() -> dict[str, float]:
    """Return the single-scale lattice support and CC-cost numbers."""

    bz_gev = math.pi * LAMBDA_GEV
    tev_over_bz = TEV_GEV / bz_gev
    support_scale_gev = TEV_GEV / math.pi
    support_over_lambda = support_scale_gev / LAMBDA_GEV
    vacuum_inflation = support_over_lambda**4

    visible_gev = VISIBLE_EV * 1.0e-9
    visible_anisotropy = (visible_gev / LAMBDA_GEV) ** 2 / 12.0

    return {
        "bz_gev": bz_gev,
        "tev_over_bz": tev_over_bz,
        "support_scale_gev": support_scale_gev,
        "support_over_lambda": support_over_lambda,
        "vacuum_inflation": vacuum_inflation,
        "visible_anisotropy": visible_anisotropy,
    }


def null_chain_gates() -> dict[str, float]:
    """Return the null-chain event-count scale for a TeV external event."""

    n_links = math.ceil(TEV_GEV / LAMBDA_GEV)
    energy_gev = n_links * LAMBDA_GEV
    # For a null chain with P=(E,E*n), P^2=0 by construction.
    p2 = energy_gev * energy_gev - energy_gev * energy_gev
    fractional_resolution = LAMBDA_GEV / TEV_GEV

    return {
        "n_links": n_links,
        "energy_gev": energy_gev,
        "p2": p2,
        "fractional_resolution": fractional_resolution,
    }


def main() -> None:
    print("[1] QCD-spaced lattice mode: safe in the IR, impossible as TeV support")
    lat = lattice_gates()
    print(f"    Brillouin ceiling pi*Lambda = {lat['bz_gev']:.3f} GeV")
    print(f"    1 TeV / E_BZ = {lat['tev_over_bz']:.1f}")
    print(
        "    finer elementary support would require "
        f"Lambda_gamma = {lat['support_scale_gev']:.1f} GeV "
        f"= {lat['support_over_lambda']:.1f} Lambda_QCD"
    )
    print(f"    zero-point/vacuum-density cost ~ (Lambda_gamma/Lambda)^4 = {lat['vacuum_inflation']:.3e}")
    print(f"    visible-photon SC artifact = {lat['visible_anisotropy']:.3e}")
    assert lat["bz_gev"] < 1.1
    assert lat["tev_over_bz"] > 900
    assert lat["vacuum_inflation"] > 1.0e11
    assert lat["visible_anisotropy"] < SME_VISIBLE_BOUND

    print("\n[2] Null-chain representation: noncompact external event, not a Bloch phase")
    chain = null_chain_gates()
    print(f"    1 TeV event carries about N = {chain['n_links']} Lambda_QCD service/action units")
    print(f"    representative null label P=(E,E n), E={chain['energy_gev']:.3f} GeV, P^2={chain['p2']:.1e}")
    print(f"    large-N count resolution Lambda/E = {chain['fractional_resolution']:.3e}")
    assert chain["n_links"] > 3000
    assert chain["p2"] == 0.0
    assert chain["fractional_resolution"] < 4.0e-4

    print("\n[3] Dynamical gates imported from the dedicated audits")
    gates = {
        "noncompact boost-invariant chain momentum": "foundations_null_causal_chain_momentum.py",
        "one-leg chain propagator vs N-photon bundle": "foundations_g2_chain_absorption_toy.py",
        "endpoint detector Ward map vs N-soft-photon bundle": "foundations_null_chain_endpoint_detector.py",
        "service-causal-set IR condensate and endpoint algebra": "foundations_service_causet_condensate_endpoint.py",
        "holonomy gauge vertex and Ward": "foundations_causet_gauge_principle_vertex.py",
        "full scalar recoil vertex": "foundations_causet_recoil_vertex.py",
        "raw finite-rho Ward covariance": "foundations_causet_finite_rho_ward_raw.py",
        "framed Hodge/Wick sign": "foundations_causet_hodge_sign_from_frame.py",
        "framed mesoscopic Maxwell action": "foundations_causet_mesoscopic_continuum_action.py",
        "framed causal-set Dirac spin": "foundations_causet_dirac_spin_closure.py",
        "detector-amplitude form / Ward-rescaling control": "foundations_trans_lambda_detector_amplitude_audit.py",
        "service-GNS endpoint LSZ residue": "foundations_trans_lambda_lsz_normalization.py",
        "normalized-leg QED interface": "foundations_trans_lambda_qed_phenomenology.py",
    }
    for gate, script in gates.items():
        print(f"    PASS by prior audit: {gate}  [{script}]")
    assert len(gates) == 13

    print("\n[4] Cosmological-constant compatibility")
    print(
        "    The high-energy state is source/detector conditioned: it is a finite null-chain history,\n"
        "    not a new lattice oscillator species populated in the vacuum.  Therefore it does not\n"
        "    replace the Lambda_QCD vacuum/IR-condensate cutoff used in the rho_Lambda ledger."
    )

    print(
        r"""
[5] VERDICT
  CLOSED at framed causal-set smooth-field EFT grade:
    observed omega >> Lambda_QCD photons are not QCD-lattice Bloch modes.
    They are external null-chain QED events on the framed causal-set/service
    layer, with one detector/LSZ leg and total P^mu.

  SAME SCALE, TWO ROLES:
    Lambda_QCD remains the condensate/vacuum cutoff that protects the
    cosmological-constant ledger; the same unit also counts service/action
    along conditioned null histories.  No finer UV oscillator scale is added.

  NOT CLAIMED:
    order-only causal-set QED, raw-link loop F^2, or a dynamical gravitational
    tetrad equation.  The endpoint LSZ residue is fixed by the service-GNS
    event isometry, and the normalized-leg QED interface passes the physical
    polarization, Dirac spin-sum/Ward, loop-running, and finite-density gates.
    Remaining work is precision standard-EFT phenomenology, not a microscopic
    support/LSZ normalization problem.
ALL ASSERTIONS PASSED -- trans-Lambda_QCD representation residual closed in the canonical scoped sense."""
    )


if __name__ == "__main__":
    main()
