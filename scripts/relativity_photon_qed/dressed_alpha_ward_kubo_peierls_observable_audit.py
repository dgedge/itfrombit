#!/usr/bin/env python3
r"""Dressed alpha: Ward/Kubo/Peierls observable audit.

Question
--------
Can the near-perfect dressed-alpha bridge/web result be promoted to a physical
low-energy QED observable?

The monitored CPTP route gives, at the native service tick,

    two Peierls-link service-label occupation -> delta/target ~= 0.9956.

But measured low-energy alpha is not normally an arbitrary pointer occupation.
Ward/Kubo says the charge is billed through the Peierls current-current /
photon self-energy slot.  This script puts the candidate observables on the
same continuum bridge/web substrate and asks whether those slots coincide.

Verdict
-------
They do not.  The service occupation is the internally selected monitored
pointer observable, but the Ward/Kubo/Peierls self-energy observable is a
different operator and a different number.  The web-dressed Kramers-Kronig
self-energy gives the standard one-loop-size undershoot, not the 0.9956
near-hit.  Therefore dressed alpha is not closed by the current Canon; a new
sector-billing theorem would have to identify the EM vertex residue with the
monitored pointer occupation, or supply a pinned Ward subtraction/normalisation
that changes the self-energy magnitude without fitting the observed shift.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

import bridge_web_lindblad_keldysh_poles as poles
import dressed_alpha_bridge_web_open_system as bw
import dressed_alpha_monitor_web_continuum_dos as cont
import dressed_alpha_rate_readout_theorem_audit as rr


ROOT = Path(__file__).resolve().parents[1]
TARGET = bw.DELTA_TARGET
M_CONTINUUM = 12_000


@dataclass(frozen=True)
class RatioRow:
    name: str
    ratio: float
    slot: str
    note: str


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def require_text(path: str, needles: tuple[str, ...]) -> None:
    text = (ROOT / path).read_text(encoding="utf-8")
    missing = [needle for needle in needles if needle not in text]
    if missing:
        raise AssertionError(f"{path} missing required text: {missing}")


def centered(mat: np.ndarray) -> np.ndarray:
    return mat - np.trace(mat) / mat.shape[0] * np.eye(mat.shape[0])


def affine_residual(target: np.ndarray, basis: np.ndarray) -> float:
    """Best Hilbert-Schmidt residual for target ~= a*basis + b*I."""

    d = target.shape[0]
    ident = np.eye(d)
    columns = np.stack([basis.reshape(-1), ident.reshape(-1)], axis=1)
    coeff, *_ = np.linalg.lstsq(columns, target.reshape(-1), rcond=None)
    fit = coeff[0] * basis + coeff[1] * ident
    denom = np.linalg.norm(centered(target))
    return float(np.linalg.norm(target - fit) / denom)


def alpha_delta_ratio(scale: float) -> float:
    delta = 137.0 / scale - 137.0
    return delta / TARGET


def native_qss_objects():
    gamma_web, gamma_escape = rr.continuum_width_operator()
    rho, obs, residual = rr.qss_state(1.0, gamma_web)
    d = cont.D
    i08 = cont.IDX[(0, 8)]
    i19 = cont.IDX[(1, 9)]
    p_link = np.zeros((d, d))
    p_link[i08, i08] = 0.5
    p_link[i19, i19] = 0.5
    p_minus = np.outer(cont.CMINUS, cont.CMINUS)
    p_plus = np.outer(cont.CPLUS, cont.CPLUS)
    return rho, obs, residual, gamma_web, gamma_escape, p_link, p_minus, p_plus


def scale_from_operator(rho: np.ndarray, operator: np.ndarray) -> float:
    d = rho.shape[0]
    return float(np.trace(operator @ rho).real) / (float(np.trace(operator).real) / d)


def qss_observable_rows(
    rho: np.ndarray,
    gamma_web: np.ndarray,
    p_link: np.ndarray,
    p_minus: np.ndarray,
    p_plus: np.ndarray,
) -> list[RatioRow]:
    d = rho.shape[0]
    hazard_scale = float(np.trace(gamma_web @ rho).real) / (float(np.trace(gamma_web).real) / d)
    return [
        RatioRow(
            "two Peierls-link service occupation",
            alpha_delta_ratio(scale_from_operator(rho, p_link)),
            "monitored pointer",
            "the near-hit route",
        ),
        RatioRow(
            "coherent eta=-1 current branch",
            alpha_delta_ratio(scale_from_operator(rho, p_minus)),
            "current branch",
            "Peierls current projector, not pointer label",
        ),
        RatioRow(
            "coherent eta=+1 current branch",
            alpha_delta_ratio(scale_from_operator(rho, p_plus)),
            "current branch",
            "orthogonal flux branch",
        ),
        RatioRow(
            "web width / escape hazard",
            alpha_delta_ratio(hazard_scale),
            "self-energy hazard",
            "absorptive slot, not service occupation",
        ),
    ]


def ward_kubo_dispersion():
    """Return Ward/Kubo dispersive self-energy ratios on the continuum DOS.

    Im Pi is the Peierls-current web emission spectrum.  Re Pi(0) is its
    Kramers-Kronig/dispersive sum.  This is the physical self-energy slot,
    distinct from the monitored pointer occupation.
    """

    h2, _pairs, idx, _bas = bw.build_pair_system()
    evals_unit, evecs = np.linalg.eigh(h2)
    evals = cont.aud.T_M * evals_unit
    e_ref = float(evals.min())
    omega = evals - e_ref
    current = bw.current_portal(idx, bw.ETA_PIN)
    amps = evecs.T @ current
    groups = poles.eigen_groups(evals)
    s_eta = cont.make_smooth_s_eta(M_CONTINUUM)

    bare = 0.0
    web = 0.0
    emission = 0.0
    ground_weight = 0.0
    for group in groups:
        weight = float(np.sum(amps[group] ** 2))
        w = float(np.mean(omega[group]))
        if w < 1e-10:
            ground_weight += weight
            continue
        spectral = s_eta(w, bw.ETA_PIN)
        bare += weight / w
        web += weight * spectral / w
        emission += weight * spectral

    gv2 = bw.ALPHA0 * bw.G_PORTAL * bw.G_PORTAL
    delta_bare = 137.0 * gv2 * bare
    delta_web = 137.0 * gv2 * web
    gamma_from_im_pi = 2.0 * np.pi * gv2 * emission
    return {
        "j_norm": float(np.sum(amps**2)),
        "j_h_j": float(current @ h2 @ current) * cont.aud.T_M,
        "ground_weight": ground_weight,
        "bare_sum": bare,
        "web_sum": web,
        "gamma_from_im_pi": gamma_from_im_pi,
        "delta_bare": delta_bare,
        "delta_web": delta_web,
        "ratio_bare": delta_bare / TARGET,
        "ratio_web": delta_web / TARGET,
    }


def main() -> None:
    print("DRESSED ALPHA: WARD/KUBO/PEIERLS OBSERVABLE AUDIT")
    print("=" * 100)

    require_text(
        "python_code/vacpol_action.py",
        (
            "Photon-fermion vertex -- FORCED to Peierls minimal coupling",
            "the minimal (Peierls) vertex is the UNIQUE one satisfying current",
            "Ward identity Pi.1 = 0",
        ),
    )
    require_text(
        "recent_papers/foundations/foundations.tex",
        (
            "Ward/Kubo billed through the Peierls current or photon self-energy slot",
            "service-label occupation",
        ),
    )

    rho, obs, residual, gamma_web, gamma_escape, p_link, p_minus, p_plus = native_qss_objects()
    rows = qss_observable_rows(rho, gamma_web, p_link, p_minus, p_plus)
    ward = ward_kubo_dispersion()

    print("\n[1] Native monitor+web QSS observables")
    print(f"  continuum Gamma_web               = {gamma_escape:.6e} Lambda")
    print(f"  native monitor pole residual       = {residual:.3e}")
    print(f"  {'observable':<38} {'delta/target':>13} {'slot':<20} note")
    for row in rows:
        print(f"  {row.name:<38} {row.ratio:>13.6f} {row.slot:<20} {row.note}")

    service = next(row for row in rows if row.name.startswith("two Peierls"))
    jminus = next(row for row in rows if row.name.startswith("coherent eta=-1"))
    hazard = next(row for row in rows if row.name.startswith("web width"))

    print("\n[2] Operator identity checks")
    print(f"  ||P_link - 1/2(P_J- + P_J+)||       = {np.linalg.norm(p_link - 0.5 * (p_minus + p_plus)):.3e}")
    print(f"  affine residual Gamma_web -> P_link = {affine_residual(gamma_web, p_link):.6f}")
    print(f"  affine residual P_J- -> P_link      = {affine_residual(p_minus, p_link):.6f}")
    print(f"  affine residual P_link -> P_J-      = {affine_residual(p_link, p_minus):.6f}")
    print(f"  commutator ||[Gamma_web,P_link]||   = {np.linalg.norm(gamma_web @ p_link - p_link @ gamma_web):.6e}")
    print("  Reading: service occupation is an incoherent pointer average of two current")
    print("  branches.  It is not the eta=-1 current branch and not the web self-energy.")

    print("\n[3] Ward/Kubo/Peierls dispersive self-energy")
    print(f"  sum_g |<E_g|J>|^2                  = {ward['j_norm']:.6f}")
    print(f"  <J|H|J>                            = {ward['j_h_j']:.3e}")
    print(f"  omega~0 current weight             = {ward['ground_weight']:.3e}")
    print(f"  bare KK sum   sum weight/omega      = {ward['bare_sum']:.6f}")
    print(f"  web KK sum    sum weight*S/omega    = {ward['web_sum']:.6f}")
    print(f"  Im-Pi cross-check Gamma             = {ward['gamma_from_im_pi']:.6e} Lambda")
    print(f"  Ward/Kubo bare loop delta/target    = {ward['ratio_bare']:.6f}")
    print(f"  Ward/Kubo web-dressed delta/target  = {ward['ratio_web']:.6f}")
    print(f"  normalisation needed for web slot   = {1.0 / ward['ratio_web']:.6f} x")

    print("\n[4] Decision")
    print(
        """
  The physical Ward/Kubo/Peierls slot is suitable in form but it does not
  rescue the 137.036 near-hit:

    * the near-hit is the monitored service-label occupation;
    * the coherent Peierls current branch gives a different response;
    * the absorptive web hazard gives a different response;
    * the dispersive Ward/Kubo self-energy gives the standard one-loop-sized
      undershoot, not the occupation value.

  Therefore the current Canon cannot identify measured low-energy QED alpha
  with the service occupation by calculation.  To close dressed alpha, one of
  two new things is needed:

    A. a sector-billing theorem mapping the EM vertex residue to the monitored
       pointer occupation despite the Ward/Kubo slot distinction; or
    B. a pinned Ward/Kubo subtraction/normalisation plus higher-loop kernel
       that moves the physical self-energy from the one-loop undershoot to the
       observed shift without using the observed shift as input.

  Until then, the service-occupation route is an excellent internal alpha-escape
  diagnostic, but not a derivation of the physical dressed electromagnetic
  coupling.
exit 0 -- Ward/Kubo/Peierls observable refutes current service-occupation closure; dressed alpha remains open.
"""
    )

    check(residual < 1e-12, "native monitor QSS is resolved")
    check(abs(service.ratio - 0.995613) < 5.0e-5, "service occupation reproduces the known near-hit")
    check(abs(jminus.ratio - service.ratio) > 0.2, "Peierls current branch is not the service occupation")
    check(abs(hazard.ratio - service.ratio) > 0.5, "web hazard/self-energy is not the service occupation")
    check(affine_residual(gamma_web, p_link) > 0.9, "web width operator is not affine-equivalent to service projector")
    check(abs(ward["gamma_from_im_pi"] / gamma_escape - 1.0) < 1.0e-12, "Im Pi reproduces the continuum web escape rate")
    check(0.30 < ward["ratio_web"] < 0.50, "web-dressed Ward/Kubo self-energy is a one-loop-sized undershoot")
    check(abs(ward["ratio_web"] - service.ratio) > 0.5, "Ward/Kubo self-energy does not equal the service near-hit")
    print("ALL ASSERTS PASSED")


if __name__ == "__main__":
    main()
