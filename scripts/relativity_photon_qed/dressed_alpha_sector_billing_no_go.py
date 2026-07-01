#!/usr/bin/env python3
"""Sector-billing no-go audit for the dressed-alpha CPTP route.

Target theorem:

    measured low-energy QED alpha bills the monitored two Peierls-link
    service-label occupation, rather than a flux/self-energy/hazard observable.

This script tests the theorem in the strongest form currently available in the
repo.  It uses three independent constraints:

  1. Ward/Kubo: the electromagnetic coupling is read from the Peierls current
     vertex / photon self-energy slot, not from an arbitrary diagonal pointer
     projector.
  2. Operator identity: the service projector would need to be the same
     low-energy billing operator as the current branch or web width, up to a
     fixed affine normalization.
  3. Item-79 emission billing: the licensed bare-alpha readout is an explicit
     emission-mode probability; the confined two-link occupation must track it
     with a fixed factor in a dressed emission model.

The result is a no-go for the theorem as stated.  The two-link occupation is a
well-selected internal service observable, but the present formalism does not
derive that measured QED bills it.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

import dressed_alpha_bridge_web_open_system as bw
import dressed_alpha_emission_readout_theorem as em
import dressed_alpha_monitor_web_continuum_dos as cont
import dressed_alpha_monitor_web_cptp_audit as aud
import dressed_alpha_rate_readout_theorem_audit as rr


ROOT = Path(__file__).resolve().parents[1]


def require_text(path: str, needles: tuple[str, ...]) -> None:
    text = (ROOT / path).read_text(encoding="utf-8")
    missing = [needle for needle in needles if needle not in text]
    if missing:
        raise AssertionError(f"{path} missing required text: {missing}")


def centered(mat: np.ndarray) -> np.ndarray:
    d = mat.shape[0]
    return mat - np.trace(mat) / d * np.eye(d)


def affine_residual(target: np.ndarray, basis: np.ndarray) -> tuple[float, float, float]:
    """Best Hilbert-Schmidt fit target ~= a*basis + b*I."""
    d = target.shape[0]
    ident = np.eye(d)
    columns = np.stack([basis.reshape(-1), ident.reshape(-1)], axis=1)
    coeff, *_ = np.linalg.lstsq(columns, target.reshape(-1), rcond=None)
    fit = coeff[0] * basis + coeff[1] * ident
    denom = np.linalg.norm(centered(target))
    resid = np.linalg.norm(target - fit) / denom if denom else float("nan")
    return float(resid), float(coeff[0].real), float(coeff[1].real)


def native_state_and_operators():
    gamma_web, _gamma_escape = rr.continuum_width_operator()
    rho, obs, residual = rr.qss_state(1.0, gamma_web)
    d = cont.D
    i08 = cont.IDX[(0, 8)]
    i19 = cont.IDX[(1, 9)]
    p_link = np.zeros((d, d))
    p_link[i08, i08] = 0.5
    p_link[i19, i19] = 0.5
    p_minus = np.outer(cont.CMINUS, cont.CMINUS)
    p_plus = np.outer(cont.CPLUS, cont.CPLUS)
    return rho, obs, residual, gamma_web, p_link, p_minus, p_plus


def delta_ratio_from_scale(scale: float) -> float:
    delta = 137.0 / scale - 137.0
    return delta / bw.DELTA_TARGET


def native_observable_ratios(
    rho: np.ndarray,
    gamma_web: np.ndarray,
    p_link: np.ndarray,
    p_minus: np.ndarray,
    p_plus: np.ndarray,
) -> dict[str, float]:
    d = rho.shape[0]
    uniform = 1.0 / d
    scale_service = float(np.trace(p_link @ rho).real) / uniform
    scale_jminus = float(np.trace(p_minus @ rho).real) / uniform
    scale_jplus = float(np.trace(p_plus @ rho).real) / uniform
    scale_hazard = float(np.trace(gamma_web @ rho).real) / (float(np.trace(gamma_web).real) / d)
    return {
        "service two-link occupation": delta_ratio_from_scale(scale_service),
        "Ward eta=-1 current branch": delta_ratio_from_scale(scale_jminus),
        "orthogonal eta=+1 branch": delta_ratio_from_scale(scale_jplus),
        "web self-energy / hazard": delta_ratio_from_scale(scale_hazard),
    }


def explicit_emission_ratios() -> tuple[float, float, float, float]:
    rows = []
    for w in (0.3, 0.5, 1.0):
        for kappa in (1e-3, 3e-3, 1e-2):
            delta_em, delta_link, _p_em, _p_link, residual = em.dressed_explicit(w, kappa, 1.0)
            rows.append(delta_link / delta_em)
            assert residual < 1e-8
    arr = np.array(rows)
    return float(arr.mean()), float(arr.std()), float(arr.min()), float(arr.max())


def main() -> None:
    require_text(
        "python_code/vacpol_action.py",
        (
            "Photon-fermion vertex -- FORCED to Peierls minimal coupling",
            "the minimal (Peierls) vertex is the UNIQUE one satisfying current",
            "The residual 'freedom' is precisely the framework's probability(1/137) <->",
        ),
    )
    require_text(
        "python_code/alpha0_count_rate_theorem.py",
        (
            "one-hot monitored service observable",
            "sector-specific proof that a given process bills this observable",
        ),
    )

    rho, obs, residual, gamma_web, p_link, p_minus, p_plus = native_state_and_operators()

    print("DRESSED-ALPHA SECTOR-BILLING NO-GO AUDIT")
    print("=" * 96)
    print("[1] Ward/Kubo billing slot")
    print("  The repo's current QED action fixes the photon vertex by Peierls minimal")
    print("  coupling.  Therefore the measured low-energy EM coupling lives in the")
    print("  current-current / photon self-energy slot unless a new sector-billing map")
    print("  identifies that slot with the monitored service projector.")

    print("\n[2] Operator identity tests on the 136-state bridge")
    resid_gamma_to_service, a_gs, b_gs = affine_residual(gamma_web, p_link)
    resid_jminus_to_service, a_ms, b_ms = affine_residual(p_minus, p_link)
    resid_service_to_jminus, a_sm, b_sm = affine_residual(p_link, p_minus)
    incoherent_identity = np.linalg.norm(p_link - 0.5 * (p_minus + p_plus))
    comm_gamma_service = np.linalg.norm(gamma_web @ p_link - p_link @ gamma_web)

    print(f"  Gamma_web ~= a P_link + b I     residual {resid_gamma_to_service:.6f}  a={a_gs:.3e} b={b_gs:.3e}")
    print(f"  P(J_-)    ~= a P_link + b I     residual {resid_jminus_to_service:.6f}  a={a_ms:.3e} b={b_ms:.3e}")
    print(f"  P_link    ~= a P(J_-) + b I     residual {resid_service_to_jminus:.6f}  a={a_sm:.3e} b={b_sm:.3e}")
    print(f"  ||P_link - 1/2(P(J_-)+P(J_+))|| = {incoherent_identity:.2e}")
    print(f"  ||[Gamma_web, P_link]||          = {comm_gamma_service:.6e}")
    print("  Reading: the selected service readout is the incoherent average of the two")
    print("  current branches.  The eta=-1 Ward branch and the web self-energy operator")
    print("  are not affine re-normalisations of that service projector.")

    print("\n[3] Same native QSS, different low-energy billing observables")
    ratios = native_observable_ratios(rho, gamma_web, p_link, p_minus, p_plus)
    for name, ratio in ratios.items():
        print(f"  {name:<34} delta/target = {ratio:.6f}")

    print("\n[4] Explicit-emission item-79 billing test")
    mean, std, lo, hi = explicit_emission_ratios()
    print(f"  confined-current delta / licensed-emission delta = {mean:.3f} +/- {std:.3f}")
    print(f"  range across (portal strength, escape rate)       = {lo:.3f} .. {hi:.3f}")
    print("  A licensed equivalence would require equality, or at least a fixed")
    print("  coupling-independent proportionality.  It fails this test.")

    print("\nVERDICT")
    print("  The requested theorem is refuted in the current formalism.")
    print("  The two Peierls-link occupation is internally selected as a service-label")
    print("  observable, but measured low-energy QED alpha is Ward/Kubo billed through")
    print("  the Peierls current/self-energy slot.  The service projector is not the")
    print("  current branch, not the web hazard operator, and not the item-79 emission")
    print("  probability in a dressed explicit-emission model.")
    print("  To make the route a derivation, one must add and prove a new sector-billing")
    print("  axiom mapping the EM vertex residue to this monitored pointer occupation.")
    print("  Without that new map, using the two-link occupation for physical dressed")
    print("  alpha remains an assignment, even though it is a well-motivated one.")

    assert residual < 1e-12
    assert abs(obs["delta"] / bw.DELTA_TARGET - ratios["service two-link occupation"]) < 1e-12
    assert resid_gamma_to_service > 0.5
    assert resid_jminus_to_service > 0.5
    assert resid_service_to_jminus > 0.5
    assert incoherent_identity < 1e-12
    assert abs(ratios["web self-energy / hazard"] - ratios["service two-link occupation"]) > 0.5
    assert std / mean > 0.5
    print("ALL ASSERTS PASSED")


if __name__ == "__main__":
    main()
