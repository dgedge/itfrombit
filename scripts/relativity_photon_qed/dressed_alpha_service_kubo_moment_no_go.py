#!/usr/bin/env python3
r"""Dressed alpha: service-occupation -> Ward/Kubo moment-map audit.

Question
--------
The service-occupation route gives the right sign/form and lands at
0.9956 of the observed dressed-alpha shift, but it is not licensed as the
physical low-energy QED coupling.  The only credible way to license it without
adding a fitted coefficient is a fluctuation-dissipation/Kramers-Kronig map:

    monitored service occupation  ->  current spectral density  ->  Re Pi(0).

This script tests that map on the actual bridge/Wilson-web operators.

Result
------
Negative closure under the current axioms.

The monitored service occupation is a diagonal service-label Born weight in
the quasi-stationary bridge state.  The Ward/Kubo electromagnetic residue is an
inverse-energy spectral moment of the Peierls current:

    Re Pi(0) ~ sum_g |<g|J>|^2 S_eta(omega_g) / omega_g.

The zeroth spectral moment, inverse-energy moment, coherent-current branch,
and service-label occupation are all different objects on the actual
bridge-web spectrum.  They would coincide only after inserting a new
normal-ordering/energy-denominator rule.  The required multiplier is the same
q_eff^2 ~= 2.54 found by the sector-billing audit, so the service occupation
cannot be promoted by a hidden FDT identity already present in canon.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

import bridge_web_lindblad_keldysh_poles as poles
import dressed_alpha_bridge_web_open_system as bw
import dressed_alpha_monitor_web_continuum_dos as cont
import dressed_alpha_rate_readout_theorem_audit as rr
import dressed_alpha_ward_kubo_peierls_observable_audit as wk


ROOT = Path(__file__).resolve().parents[1]
TARGET = bw.DELTA_TARGET
M_CONTINUUM = 12_000


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def require_text(path: str, needles: tuple[str, ...]) -> None:
    text = (ROOT / path).read_text(encoding="utf-8")
    missing = [needle for needle in needles if needle not in text]
    if missing:
        raise AssertionError(f"{path} missing required text: {missing}")


def current_spectral_moments() -> dict[str, float]:
    """Peierls-current spectral moments on the continuum Wilson-web DOS."""

    h2, _pairs, idx, _bas = bw.build_pair_system()
    evals_unit, evecs = np.linalg.eigh(h2)
    evals = cont.aud.T_M * evals_unit
    e_ref = float(evals.min())
    omega = evals - e_ref
    current = bw.current_portal(idx, bw.ETA_PIN)
    amps = evecs.T @ current
    groups = poles.eigen_groups(evals)
    s_eta = cont.make_smooth_s_eta(M_CONTINUUM)

    rows: list[tuple[float, float, float]] = []
    ground_weight = 0.0
    for group in groups:
        weight = float(np.sum(amps[group] ** 2))
        w = float(np.mean(omega[group]))
        if w < 1e-10:
            ground_weight += weight
            continue
        spectral = weight * s_eta(w, bw.ETA_PIN)
        if spectral > 0.0:
            rows.append((w, weight, spectral))

    omegas = np.array([r[0] for r in rows])
    weights = np.array([r[2] for r in rows])
    m_minus1 = float(np.sum(weights / omegas))
    m0 = float(np.sum(weights))
    m_plus1 = float(np.sum(weights * omegas))
    mean = float(np.sum(weights * omegas) / m0)
    variance = float(np.sum(weights * (omegas - mean) ** 2) / m0)
    eff_denominator = m0 / m_minus1
    return {
        "n_lines": float(len(rows)),
        "j_norm": float(np.sum(amps**2)),
        "ground_weight": ground_weight,
        "m_minus1": m_minus1,
        "m0": m0,
        "m_plus1": m_plus1,
        "mean_omega": mean,
        "std_omega": math.sqrt(variance),
        "cv_omega": math.sqrt(variance) / mean,
        "eff_denominator": eff_denominator,
        "omega_min": float(np.min(omegas)),
        "omega_max": float(np.max(omegas)),
    }


def service_ratio() -> tuple[float, float]:
    gamma_web, _gamma_escape = rr.continuum_width_operator()
    _rho, obs, residual = rr.qss_state(1.0, gamma_web)
    return float(obs["delta"] / TARGET), residual


def main() -> None:
    require_text(
        "python_code/dressed_alpha_sector_billing_subtraction_theorem.py",
        (
            "The service occupation is a diagonal monitored pointer observable",
            "The EM vertex residue is the Ward/Kubo Peierls-current Hessian",
        ),
    )
    require_text(
        "holographic_circlette_book/part_3_standard_model/ch13_fine_structure_constant.tex",
        (
            "service-occupation route",
            "Ward/Kubo moment-map audit rules out licensing it",
        ),
    )

    service, residual = service_ratio()
    ward = wk.ward_kubo_dispersion()
    moments = current_spectral_moments()
    ward_web = float(ward["ratio_web"])
    multiplier_service_to_ward = service / ward_web
    q_eff = math.sqrt(1.0 / ward_web)

    print("DRESSED ALPHA: SERVICE-OCCUPATION -> WARD/KUBO MOMENT-MAP AUDIT")
    print("=" * 96)
    print("[1] The two candidate readings")
    print(f"  monitored service occupation / observed shift = {service:.6f}")
    print(f"  QSS residual                                  = {residual:.3e}")
    print(f"  Ward/Kubo web self-energy / observed shift    = {ward_web:.6f}")
    print(f"  service / Ward-web multiplier                 = {multiplier_service_to_ward:.6f}")
    print(f"  q_eff needed if billed as Peierls charge       = {q_eff:.6f}")
    check(abs(service - 1.0) < 0.01, "service occupation is the near-hit")
    check(ward_web < 0.5, "physical Ward/Kubo slot undershoots the observed shift")
    check(multiplier_service_to_ward > 2.0, "service-to-Ward map would need an O(1) new multiplier")

    print("\n[2] Spectral moment hierarchy")
    print(f"  nonzero current spectral lines          = {moments['n_lines']:.0f}")
    print(f"  sum_g |<g|J>|^2                         = {moments['j_norm']:.6f}")
    print(f"  omega~0 current weight                  = {moments['ground_weight']:.3e}")
    print(f"  M_-1 = sum weight*S(omega)/omega        = {moments['m_minus1']:.6f}")
    print(f"  M_0  = sum weight*S(omega)              = {moments['m0']:.6f}")
    print(f"  M_+1 = sum weight*S(omega)*omega        = {moments['m_plus1']:.6f}")
    print(f"  effective denominator M_0/M_-1          = {moments['eff_denominator']:.6f} Lambda")
    print(f"  weighted omega mean +/- std             = {moments['mean_omega']:.6f} +/- {moments['std_omega']:.6f}")
    print(f"  coefficient of variation                = {moments['cv_omega']:.3f}")
    print(f"  omega support                           = {moments['omega_min']:.6f} .. {moments['omega_max']:.6f}")
    print("  A hidden FDT identification would have to collapse this hierarchy to a")
    print("  selected energy denominator. The actual spectrum is broad, so the")
    print("  inverse-energy Kubo moment is not a disguised occupation count.")
    check(moments["cv_omega"] > 0.10, "current spectral support is not monoenergetic")
    check(moments["omega_max"] / moments["omega_min"] > 2.0, "Kubo denominator varies across the support")

    print("\n[3] Negative closure")
    print(
        """
  The requested positive theorem would need:

      service occupation = Ward/Kubo Peierls residue.

  The actual algebra gives:

      service occupation = diagonal monitored pointer Born weight,
      Ward/Kubo residue   = inverse-energy current spectral moment.

  They have the same sign and both use the Peierls bridge/web object, which is
  why the service route is physically suggestive. But the map between them is
  not forced by Ward, Kubo, FDT, or the current service algebra. The missing
  factor is q_eff^2 ~= 2.54, equivalent to a new charge normalisation or a new
  finite normal-ordering subtraction.

  Verdict: the service-occupation route is closed negatively under current
  axioms. It remains the best internal near-hit, but it is not physical QED
  alpha unless canon adds a new EM sector-billing / normal-ordering theorem.
exit 0 -- no hidden Ward/Kubo moment map licenses the service occupation.
"""
    )


if __name__ == "__main__":
    main()
