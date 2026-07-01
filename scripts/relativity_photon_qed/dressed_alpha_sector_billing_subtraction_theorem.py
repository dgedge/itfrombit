#!/usr/bin/env python3
r"""Dressed alpha: sector-billing / Ward-subtraction theorem audit.

Target
------
Try to close the remaining dressed-alpha gap by one of two routes:

  A. Sector-billing theorem:
     identify the physical low-energy EM vertex residue with the monitored
     two-Peierls-link service occupation.

  B. Pinned Ward/Kubo normalisation/subtraction:
     keep the physical Peierls current-current slot, but derive the missing
     finite factor/subtraction that moves the web-dressed self-energy to the
     observed shift.

Result
------
Under the current canon axioms both routes fail as closure theorems.

The service occupation is a diagonal monitored pointer observable.  The
physical EM coupling is the second derivative of the action with respect to
the Peierls gauge phase, i.e. the Ward/Kubo current-current plus diamagnetic
slot.  These are different operators, not affine re-normalisations of one
another.  A finite Ward/Kubo multiplication by the needed factor is equivalent
to assigning an effective charge q_eff != 1 to the Peierls phase; that is a new
sector charge, not the derived alpha0=1/137 service charge.

So this is a theorem-shaped negative: dressed alpha cannot be promoted by a
sector-billing map or finite Ward/Kubo normalisation already present in the
current framework.  Closure requires a genuinely new bridge principle.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

import dressed_alpha_bridge_web_open_system as bw
import dressed_alpha_emission_readout_theorem as em
import dressed_alpha_monitor_web_continuum_dos as cont
import dressed_alpha_rate_readout_theorem_audit as rr
import dressed_alpha_ward_kubo_peierls_observable_audit as wk


ROOT = Path(__file__).resolve().parents[1]
TARGET = bw.DELTA_TARGET


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


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
    cols = np.stack([basis.reshape(-1), ident.reshape(-1)], axis=1)
    coeff, *_ = np.linalg.lstsq(cols, target.reshape(-1), rcond=None)
    fit = coeff[0] * basis + coeff[1] * ident
    denom = np.linalg.norm(centered(target))
    return float(np.linalg.norm(target - fit) / denom)


def native_operators() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    gamma_web, _gamma_escape = rr.continuum_width_operator()
    rho, _obs, residual = rr.qss_state(1.0, gamma_web)
    assert residual < 1e-12

    d = cont.D
    p_link = np.zeros((d, d))
    p_link[cont.IDX[(0, 8)], cont.IDX[(0, 8)]] = 0.5
    p_link[cont.IDX[(1, 9)], cont.IDX[(1, 9)]] = 0.5
    p_minus = np.outer(cont.CMINUS, cont.CMINUS)
    p_plus = np.outer(cont.CPLUS, cont.CPLUS)
    return rho, gamma_web, p_link, p_minus, p_plus


def explicit_emission_ratio_spread() -> tuple[float, float, float, float]:
    ratios = []
    for w in (0.3, 0.5, 1.0):
        for kappa in (1e-3, 3e-3, 1e-2):
            delta_em, delta_link, _p_em, _p_link, residual = em.dressed_explicit(w, kappa, 1.0)
            assert residual < 1e-8
            ratios.append(delta_link / delta_em)
    arr = np.array(ratios)
    return float(arr.mean()), float(arr.std()), float(arr.min()), float(arr.max())


def ward_scaling_numbers() -> dict[str, float]:
    ward = wk.ward_kubo_dispersion()
    ratio_web = float(ward["ratio_web"])
    needed_multiplier = 1.0 / ratio_web
    q_eff = math.sqrt(needed_multiplier)
    counterterm = TARGET * (1.0 - ratio_web)
    return {
        "ratio_web": ratio_web,
        "needed_multiplier": needed_multiplier,
        "q_eff": q_eff,
        "counterterm": counterterm,
        "counterterm_over_target": counterterm / TARGET,
    }


def main() -> None:
    require_text(
        "python_code/alpha0_record_pair_symmetry_theorem.py",
        (
            "RECORD-pairs (symmetric, 136), not MATTER-pairs",
            "Sym^2(16)+1",
        ),
    )
    require_text(
        "python_code/vacpol_action.py",
        (
            "Photon-fermion vertex -- FORCED to Peierls minimal coupling",
            "the minimal (Peierls) vertex is the UNIQUE one satisfying current",
            "Ward identity Pi.1 = 0",
        ),
    )
    require_text(
        "python_code/dressed_alpha_status_consolidation.py",
        (
            "Bare 137 is derived.  The dressed shift is not currently derived.",
            "requires a new sector-billing theorem",
        ),
    )

    rho, gamma_web, p_link, p_minus, p_plus = native_operators()
    ward = ward_scaling_numbers()

    print("DRESSED ALPHA: SECTOR-BILLING / WARD-SUBTRACTION THEOREM AUDIT")
    print("=" * 100)
    print("[1] Sector-billing operator identity")
    print(f"  ||P_link - 1/2(P_J- + P_J+)||       = {np.linalg.norm(p_link - 0.5 * (p_minus + p_plus)):.3e}")
    print(f"  affine residual Gamma_web -> P_link = {affine_residual(gamma_web, p_link):.6f}")
    print(f"  affine residual P_J- -> P_link      = {affine_residual(p_minus, p_link):.6f}")
    print(f"  affine residual P_link -> P_J-      = {affine_residual(p_link, p_minus):.6f}")
    print("  Interpretation: the monitor observable is the incoherent pointer average")
    print("  of two current branches.  The EM vertex residue is a Peierls-current")
    print("  second derivative of the action, not this pointer projector.")
    check(np.linalg.norm(p_link - 0.5 * (p_minus + p_plus)) < 1e-12, "service projector is the branch average")
    check(affine_residual(gamma_web, p_link) > 0.9, "web self-energy is not an affine service projector")
    check(affine_residual(p_minus, p_link) > 0.6, "current branch is not an affine service projector")

    print("\n[2] Explicit emission billing")
    mean, std, lo, hi = explicit_emission_ratio_spread()
    print(f"  confined-link delta / licensed-emission delta = {mean:.3f} +/- {std:.3f}")
    print(f"  range across emission models                  = {lo:.3f} .. {hi:.3f}")
    print("  A sector-billing theorem would require equality, or at least one fixed")
    print("  proportionality.  The ratio is model-dependent.")
    check(std / mean > 0.5, "confined-link/emission ratio is not fixed")

    print("\n[3] Pinned Ward/Kubo normalisation test")
    print(f"  web-dressed Ward/Kubo delta/target        = {ward['ratio_web']:.6f}")
    print(f"  multiplier needed to hit observed shift   = {ward['needed_multiplier']:.6f}")
    print(f"  equivalent Peierls charge q_eff           = {ward['q_eff']:.6f}")
    print(f"  finite counterterm needed in delta(alpha^-1) = {ward['counterterm']:.9f}")
    print(f"  counterterm / observed shift              = {ward['counterterm_over_target']:.6f}")
    print("  A multiplicative rescue is just q_eff^2.  Since the Peierls phase is the")
    print("  unit-charge gauge derivative and alpha0=1/137 already fixes the service")
    print("  charge, q_eff != 1 is a new charge-normalisation postulate.  A finite")
    print("  additive rescue is a local F^2 counterterm; gauge invariance permits it")
    print("  but supplies no value for it.")
    check(abs(ward["q_eff"] - 1.0) > 0.5, "needed Ward/Kubo rescaling is not the unit Peierls charge")
    check(0.5 < ward["counterterm_over_target"] < 0.7, "subtraction is a large fitted finite term")

    print("\n[4] Candidate near-constants are not a theorem")
    candidates = {
        "sqrt(2*pi)": math.sqrt(2.0 * math.pi),
        "8/pi": 8.0 / math.pi,
        "13/5": 13.0 / 5.0,
        "33/13": 33.0 / 13.0,
        "137/54": 137.0 / 54.0,
    }
    for name, value in candidates.items():
        ratio = value * ward["ratio_web"]
        print(f"  {name:<10} factor={value:.6f} -> delta/target={ratio:.6f}")
    print("  Some simple-looking constants sit numerically nearby.  None is selected by")
    print("  the Ward identity, Peierls continuity equation, or service projector")
    print("  algebra.  Promoting one would be exactly the N1=31 failure mode again.")

    print("\nTHEOREM STATUS")
    print(
        """
  Negative closure under current axioms:

    1. The EM vertex residue is the Ward/Kubo Peierls-current Hessian of the
       action.  Its normalisation is fixed by the unit Peierls phase and the
       bare service charge alpha0.

    2. The monitored two-link occupation is a valid internal service-label
       pointer observable, but it is an incoherent branch average.  It is not
       the current-current/self-energy operator and is not the licensed
       emission probability in explicit emission models.

    3. The finite Ward/Kubo multiplier needed to hit the physical shift is
       q_eff^2 = 2.537609.  That is equivalent to inserting a new effective
       charge q_eff = 1.592 instead of deriving it.  The additive alternative
       is a finite F^2 counterterm of 0.606 of the observed shift, likewise
       unselected by the present algebra.

  Therefore neither requested positive theorem is derivable from the current
  canon.  The route can close only by adding a new, explicit EM sector-billing
  principle or a new normal-ordering/renormalisation rule that fixes the
  finite Ward/Kubo term independently of the observed alpha.
exit 0 -- current framework gives a no-go, not a dressed-alpha closure.
"""
    )


if __name__ == "__main__":
    main()
