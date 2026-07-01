#!/usr/bin/env python3
r"""Minimal U(1) plaquette action toy.

The previous plaquette examples established that a closed Wilson loop is a
detector-readable record.  This script adds a separate idea:

    a closed flux can also carry a dynamical weight.

For one U(1) plaquette, write the loop as

    W = exp(i Phi),

and assign the Wilson plaquette action

    S(Phi) = beta (1 - cos Phi).

This has the small-flux expansion

    S(Phi) = beta Phi^2/2 + O(Phi^4),

the finite toy analogue of an F^2 action.

Record-grammar distinction
--------------------------

    Readout:   detector probabilities depend on cos(Phi).
    Dynamics:  flux histories/configurations are weighted by exp[-S(Phi)].

Those are different grammar roles.  "Closed flux is readable" is not yet
"closed flux has Maxwell-like dynamics"; an action is an additional rule.

This script is not a derivation of Maxwell theory.  It is the smallest exact
calculation that separates gauge-invariant readout from action weighting.
"""

from __future__ import annotations

import math

import numpy as np


def wrap_angle(phi: float) -> float:
    """Map phi to (-pi, pi]."""

    return float((phi + math.pi) % (2.0 * math.pi) - math.pi)


def wilson_action(phi: float, beta: float) -> float:
    return beta * (1.0 - math.cos(phi))


def bright_dark(phi: float) -> tuple[float, float]:
    """Aharonov-Bohm bright/dark readout for a single U(1) loop flux."""

    return (1.0 + math.cos(phi)) / 2.0, (1.0 - math.cos(phi)) / 2.0


def discrete_flux_grid(n: int) -> np.ndarray:
    """Midpoint grid on [-pi, pi), avoiding duplicated endpoints."""

    return np.array([-math.pi + (k + 0.5) * 2.0 * math.pi / n for k in range(n)], dtype=float)


def action_weighted_distribution(beta: float, n: int = 4096) -> tuple[np.ndarray, np.ndarray]:
    phis = discrete_flux_grid(n)
    actions = np.array([wilson_action(float(phi), beta) for phi in phis])
    # Shift by the minimum action for numerical stability.
    weights = np.exp(-(actions - float(np.min(actions))))
    probs = weights / float(np.sum(weights))
    return phis, probs


def expectation(phis: np.ndarray, probs: np.ndarray, fn) -> float:
    return float(np.sum(probs * np.array([fn(float(phi)) for phi in phis])))


def circular_variance_from_cos(mean_cos: float) -> float:
    return 1.0 - mean_cos


def modified_bessel_i1(x: float, tol: float = 1e-16) -> float:
    """Modified Bessel I_1(x) by its rapidly convergent power series.

    I_1(x) = sum_m (x/2)^(2m+1) / (m! (m+1)!).

    NumPy provides i0 in this environment but not i1, and this avoids a SciPy
    dependency for a one-line analytic check.
    """

    term = x / 2.0
    total = term
    m = 0
    while True:
        # term_{m+1} / term_m = (x/2)^2 / ((m+1)(m+2)).
        term *= (x * x / 4.0) / ((m + 1) * (m + 2))
        total += term
        m += 1
        if abs(term) < tol * max(1.0, abs(total)):
            return total
        if m > 200:
            raise RuntimeError("I1 series did not converge")


def bessel_ratio_i1_i0(beta: float) -> float:
    """Exact continuum <cos Phi> for weight exp[beta cos Phi]."""

    return modified_bessel_i1(beta) / float(np.i0(beta))


def assert_close(name: str, value: float, target: float, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<72s} value={value:.12g} target={target:.12g} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def assert_less(name: str, value: float, bound: float) -> None:
    print(f"  {name:<72s} value={value:.12g} bound={bound:.12g}")
    if not value < bound:
        raise AssertionError(name)


def main() -> None:
    print("Minimal U(1) plaquette action toy")
    print("=" * 88)

    print("\n[1] Readout and action are different functions of the same closed flux")
    phi = 0.73
    beta = 2.4
    p_bright, p_dark = bright_dark(phi)
    action = wilson_action(phi, beta)
    assert_close("P_bright", p_bright, (1.0 + math.cos(phi)) / 2.0)
    assert_close("P_dark", p_dark, (1.0 - math.cos(phi)) / 2.0)
    assert_close("S(Phi)", action, beta * (1.0 - math.cos(phi)))
    print("  -> readout gives detector probabilities; action gives configuration weight.")

    print("\n[2] Small-flux action becomes quadratic")
    for small_phi in [1.0e-3, -2.0e-3, 5.0e-3]:
        exact = wilson_action(small_phi, beta)
        quadratic = beta * small_phi * small_phi / 2.0
        relative_error = abs(exact - quadratic) / quadratic
        assert_less(f"relative error at Phi={small_phi:+.1e}", relative_error, 3.0e-6)
    print("  -> this is the one-plaquette toy analogue of an F^2 term.")

    print("\n[3] Gauge/periodic invariance belongs to closed flux")
    shifted = phi + 2.0 * math.pi
    assert_close("S(Phi+2pi)=S(Phi)", wilson_action(shifted, beta), action)
    assert_close("P_bright(Phi+2pi)=P_bright(Phi)", bright_dark(shifted)[0], p_bright)
    assert_close("wrapped flux", wrap_angle(shifted), wrap_angle(phi))
    print("  -> the action weights a closed compact loop variable, not an open link phase.")

    print("\n[4] Action weighting concentrates flux as beta grows")
    beta_values = [0.0, 0.5, 2.0, 8.0]
    last_mean_cos = -1.0
    for b in beta_values:
        phis, probs = action_weighted_distribution(b)
        mean_cos = expectation(phis, probs, math.cos)
        mean_action = expectation(phis, probs, lambda x: wilson_action(x, b))
        exact_mean_cos = bessel_ratio_i1_i0(b)
        assert_close(f"<cos Phi> beta={b}", mean_cos, exact_mean_cos, tol=2.0e-7)
        print(
            f"  beta={b:<4.1f}  <cos Phi>={mean_cos:.9f}  "
            f"circ_var={circular_variance_from_cos(mean_cos):.9f}  <S>={mean_action:.9f}"
        )
        if mean_cos + 1e-12 < last_mean_cos:
            raise AssertionError("mean cos decreased as beta increased")
        last_mean_cos = mean_cos
    print("  -> high beta selects small flux; low beta leaves broad flux fluctuations.")

    print("\n[5] The detector readout averaged over dynamics is a derived ensemble quantity")
    for b in [0.0, 1.0, 4.0]:
        phis, probs = action_weighted_distribution(b)
        mean_bright = expectation(phis, probs, lambda x: bright_dark(x)[0])
        exact_mean_bright = (1.0 + bessel_ratio_i1_i0(b)) / 2.0
        assert_close(f"<P_bright> beta={b}", mean_bright, exact_mean_bright, tol=2.0e-7)
    print("  -> observable readout plus action weighting gives an ensemble prediction.")

    print(
        """
Verdict
-------
The plaquette-action toy adds one grammar distinction:

  closed Wilson loop = detector-readable flux record;
  Wilson plaquette action = dynamical weight over flux records.

The action S(Phi)=beta(1-cos Phi) is compact, gauge-invariant, and quadratic at
small flux.  Larger beta suppresses flux fluctuations, and the action-weighted
mean detector signal follows the exact Bessel ratio I1(beta)/I0(beta).

Boundary
--------
This is not a Maxwell derivation.  It is the smallest exact calculation that
separates "the loop can be read" from "the loop is dynamically weighted by an
F^2-like action."  A real Maxwell/Yang-Mills derivation would still need a
principled measure, continuum limit, and field normalization.
"""
    )


if __name__ == "__main__":
    main()
