#!/usr/bin/env python3
"""Rate/readout theorem audit for the dressed-alpha CPTP route.

The continuum CPTP calculation is numerically sharp:

    native pair-basis monitor, rate = 1
    + Peierls/Wilson web width operator
    + two current-link service-label occupation readout
    -> delta/target = 0.9956.

This script audits the remaining theorem half:

  1. Is monitor rate = 1 licensed independently of alpha?
  2. Is the two Peierls current-link occupation the selected dressed-alpha
     observable, or just the observable that lands near the number?

The answer is deliberately conditional.  The internal service-form theorem can
be made tight: rate=1 is the R14/item-79 native service interrogation clock, and
the two current-link projector occupation is the unique readout satisfying
simultaneously:

  * alpha is a monitored service-label Born weight, not a flux hazard;
  * the Peierls/Ward vertex has support only on bridge links (0,8), (1,9);
  * the native monitor pointer basis is the pair-label basis.

What this does not prove is the external sector-billing statement that the
measured low-energy QED coupling must bill this monitored service-label
observable rather than a different self-energy/flux observable.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

import dressed_alpha_bridge_web_open_system as bw
import dressed_alpha_monitor_web_cptp_audit as aud
import dressed_alpha_monitor_web_continuum_dos as cont


ROOT = Path(__file__).resolve().parents[1]
TARGET = bw.DELTA_TARGET
M_CONTINUUM = 12_000


@dataclass(frozen=True)
class ObservableRow:
    name: str
    scale: float
    delta: float
    ratio: float
    status: str
    reason: str


def require_text(path: str, needles: tuple[str, ...]) -> None:
    text = (ROOT / path).read_text(encoding="utf-8")
    missing = [needle for needle in needles if needle not in text]
    if missing:
        raise AssertionError(f"{path} missing required text: {missing}")


def continuum_width_operator():
    s_eta = cont.make_smooth_s_eta(M_CONTINUUM)
    gamma_web, gamma_escape, bright, dark = aud.grouped_web_width_operator(cont.H2, cont.IDX, s_eta)
    assert bright == 32
    assert dark == 104
    return gamma_web, gamma_escape


def qss_state(rate: float, gamma_web: np.ndarray):
    cand = aud.MonitorCandidate("native pair-basis monitor", "pair_site", rate, "audit", "")
    matvec = aud.make_matvec(cont.H, gamma_web, cand, cont.LABELS, cont.CMINUS, cont.CPLUS)
    eigenvalue, eigenvector, residual, solver = aud.scipy_dominant(matvec, cont.D * cont.D)
    rho = eigenvector.reshape((cont.D, cont.D))
    rho = 0.5 * (rho + rho.conj().T)
    trace = float(np.trace(rho).real)
    if trace < 0:
        rho = -rho
        trace = -trace
    rho /= trace
    obs = aud.qss_observables(eigenvalue, eigenvector, cont.IDX, cont.D)
    mode = aud.row_mode(obs, residual)
    assert solver == "scipy-eigs"
    assert mode == "qss"
    return rho, obs, residual


def alpha_delta(scale: float) -> tuple[float, float]:
    delta = 137.0 / scale - 137.0
    return delta, delta / TARGET


def observable_rows(rho: np.ndarray, gamma_web: np.ndarray) -> list[ObservableRow]:
    diag = np.diag(rho).real
    uniform = 1.0 / cont.D
    i08 = cont.IDX[(0, 8)]
    i19 = cont.IDX[(1, 9)]
    i09 = cont.IDX[(0, 9)]
    i18 = cont.IDX[(1, 8)]
    p_minus = np.outer(cont.CMINUS, cont.CMINUS)
    p_plus = np.outer(cont.CPLUS, cont.CPLUS)

    candidates = [
        (
            "two Peierls link projectors",
            0.5 * (diag[i08] + diag[i19]) / uniform,
            "selected",
            "service-label Born weight + Peierls on-link support",
        ),
        (
            "individual link (0,8)",
            diag[i08] / uniform,
            "selected control",
            "one of the two symmetry-related current links",
        ),
        (
            "individual link (1,9)",
            diag[i19] / uniform,
            "selected control",
            "one of the two symmetry-related current links",
        ),
        (
            "off-link cross pair (0,9)",
            diag[i09] / uniform,
            "rejected",
            "Peierls current has no single-link element on this pair",
        ),
        (
            "off-link cross pair (1,8)",
            diag[i18] / uniform,
            "rejected",
            "Peierls current has no single-link element on this pair",
        ),
        (
            "coherent J_- portal projector",
            float(np.trace(p_minus @ rho).real) / uniform,
            "rejected",
            "not a native pair-label service projector",
        ),
        (
            "coherent J_+ portal projector",
            float(np.trace(p_plus @ rho).real) / uniform,
            "rejected",
            "orthogonal flux branch, not the eta=-1 service label",
        ),
        (
            "total web escape hazard",
            float(np.trace(gamma_web @ rho).real) / (float(np.trace(gamma_web).real) / cont.D),
            "rejected",
            "flux/self-energy observable, not service-label alpha weight",
        ),
        (
            "diagonal web hazard",
            float(np.sum(np.diag(gamma_web).real * diag)) / (float(np.trace(gamma_web).real) / cont.D),
            "rejected",
            "hazard-weighted flux, not a one-hot service projector",
        ),
    ]

    rows = []
    for name, scale, status, reason in candidates:
        delta, ratio = alpha_delta(scale)
        rows.append(ObservableRow(name, scale, delta, ratio, status, reason))
    return rows


def rate_ratio(rate: float, gamma_web: np.ndarray) -> float:
    _rho, obs, _res = qss_state(rate, gamma_web)
    return obs["delta"] / TARGET


def find_upper_crossing(gamma_web: np.ndarray) -> float:
    low, high = 1.0, 1.1
    y_low = rate_ratio(low, gamma_web) - 1.0
    y_high = rate_ratio(high, gamma_web) - 1.0
    assert y_low * y_high <= 0
    for _ in range(10):
        mid = 0.5 * (low + high)
        y_mid = rate_ratio(mid, gamma_web) - 1.0
        if y_low * y_mid <= 0:
            high = mid
            y_high = y_mid
        else:
            low = mid
            y_low = y_mid
    return 0.5 * (low + high)


def provenance_checks() -> None:
    require_text(
        "python_code/alpha0_count_rate_theorem.py",
        (
            "one-hot monitored service observable",
            "one label per interrogation",
            "one lattice interrogation per tick",
            "Gamma = alpha0 Lambda",
        ),
    )
    require_text(
        "python_code/item79_unital_channel.py",
        (
            "site-basis dephasing jumps",
            "unique fixed point",
            "P(emission)=1/137",
        ),
    )
    require_text(
        "python_code/omega_em_portal_vertex.py",
        (
            "Peierls current lives ON LINKS",
            "(0,8) and (1,9)",
            "off-link cross pairs (0,9),(1,8)",
        ),
    )
    require_text(
        "python_code/vacpol_action.py",
        (
            "Photon-fermion vertex -- FORCED to Peierls minimal coupling",
            "the minimal (Peierls) vertex is the UNIQUE one satisfying current",
        ),
    )


def main() -> None:
    provenance_checks()
    gamma_web, gamma_escape = continuum_width_operator()
    rho_native, obs_native, residual_native = qss_state(1.0, gamma_web)
    ratio_native = obs_native["delta"] / TARGET
    ratio_matter = rate_ratio(aud.T_M, gamma_web)
    upper_crossing = find_upper_crossing(gamma_web)
    rows = observable_rows(rho_native, gamma_web)

    print("DRESSED-ALPHA RATE + READOUT THEOREM AUDIT")
    print("=" * 104)
    print("[1] Rate provenance")
    print("  Existing R14/item-79 theorem gives one service-label interrogation per")
    print("  substrate tick.  In dimensionless Lambda units this is monitor rate = 1.")
    print("  The Grover matter scale t_m/t_web = 1/3 belongs to the coherent Hamiltonian,")
    print("  not to the syndrome/service monitor.")
    print()
    print(f"  continuum Gamma_web                 = {gamma_escape:.6e} Lambda")
    print(f"  native monitor rate                 = 1.000000 -> delta/target {ratio_native:.6f}")
    print(f"  matter hopping scale as monitor     = {aud.T_M:.6f} -> delta/target {ratio_matter:.6f}")
    print(f"  fitted upper crossing control       = {upper_crossing:.6f}")
    print(f"  upper crossing - native rate        = {upper_crossing - 1.0:+.6f}")
    print("  Reading: rate=1 is licensed only if the dressed-alpha bridge uses the")
    print("  same native service clock as the bare alpha0 monitor.  Choosing the exact")
    print("  crossing rate would be a fit unless independently derived.")

    print("\n[2] Readout observable controls at native rate = 1")
    print(f"  {'observable':<32} {'scale':>11} {'delta':>12} {'/target':>10} {'status':<16} reason")
    for row in rows:
        print(
            f"  {row.name:<32} {row.scale:>11.6f} {row.delta:>12.6e}"
            f" {row.ratio:>10.6f} {row.status:<16} {row.reason}"
        )

    selected = next(row for row in rows if row.name == "two Peierls link projectors")
    link_08 = next(row for row in rows if row.name == "individual link (0,8)")
    link_19 = next(row for row in rows if row.name == "individual link (1,9)")
    hazard = next(row for row in rows if row.name == "total web escape hazard")
    jminus = next(row for row in rows if row.name == "coherent J_- portal projector")

    print("\n[3] Conditional theorem status")
    print("  Internal service-form theorem: PASS, conditionally.")
    print("    If dressed alpha is a monitored service-label Born weight, the native")
    print("    pair-label projector basis selects occupations, not flux hazards.")
    print("    Peierls/Ward support then restricts the readout to the two on-link")
    print("    bridge pair labels (0,8) and (1,9).")
    print("  External sector-billing theorem: OPEN.")
    print("    The audit does not prove that measured low-energy QED alpha must bill")
    print("    this service-label occupation rather than a self-energy or flux observable.")

    print("\nVERDICT")
    print("  The remaining theorem half can be narrowed but not fully discharged here.")
    print("  Rate=1 is not a numerical fit inside this route: it is the existing native")
    print("  service-clock theorem, while t_m=1/3 is a Hamiltonian scale and gives the")
    print("  wrong monitor object.  The two-current-link occupation is the unique")
    print("  readout satisfying the service-label, Peierls-support, and pointer-basis")
    print("  constraints.  But the decisive external statement remains open: physical")
    print("  dressed QED must be shown to bill that monitored occupation observable.")
    print("  Without that sector-billing theorem, the result stays conditional, not")
    print("  a final derivation of 137.036.")

    assert abs(ratio_native - 0.995613) < 5e-5
    assert ratio_matter > 1.25
    assert 1.0 < upper_crossing < 1.05
    assert abs(selected.ratio - ratio_native) < 1e-12
    assert abs(link_08.ratio - link_19.ratio) < 1e-10
    assert abs(hazard.ratio - selected.ratio) > 0.5
    assert abs(jminus.ratio - selected.ratio) > 0.2
    assert residual_native < 1e-12
    print("ALL ASSERTS PASSED")


if __name__ == "__main__":
    main()
