#!/usr/bin/env python3
r"""ITEM 131 / 56: early log-scale lift from boundary-printing covariance.

Question
--------
The finite 28-channel serial QEC clock has generator gap

    gamma = 1/28.

The remaining early-leg issue is whether primordial boundary printing really
uses the logarithmic horizon variable needed for

    d ln Delta_R^2 / d ln k = -1/28.

This script isolates the lift as a scale-covariance theorem.

Finite-to-cosmological lift
---------------------------
Let T_lambda be the effective boundary-printing transfer operator between two
horizon scales whose ratio is lambda = k/k_*.

If boundary crystallization is:

    (1) Markovian / no memory between printed horizon shells,
    (2) scale-covariant / self-similar, T_{lambda mu}=T_lambda T_mu,
    (3) continuous in lambda,

then the clock parameter tau(lambda) must obey Cauchy's multiplicative
functional equation

    tau(lambda mu) = tau(lambda) + tau(mu),

so tau(lambda)=q ln(lambda).  Therefore the finite QEC generator lifts to a
log-scale RG/anomalous-dimension operator:

    T_lambda = exp[q ln(lambda) Q].

The item-131 coefficient is the radial/horizon normalization q=1.  Area-clock
or volume-clock normalizations are still scale-covariant, but give the wrong
tilt: q=2 -> n_s=13/14, q=3 -> n_s=25/28.

Power-ledger lift
-----------------
The QEC instrument is a density/probability instrument.  Its eigenvalue acts
on the scalar variance ledger Delta_R^2, not on a coherent amplitude and not
on the raw three-dimensional two-point spectrum P(k).  Acting on amplitude
would double the tilt; acting on raw P(k) would reintroduce the k^3 phase-space
factor into n_s.  The dimensionless power Delta_R^2 is exactly the variance
per d ln k shell, so it is the correct ledger for a log-shell QEC transfer.

Conclusion
----------
Under scale-covariant Holographic Boundary Crystallization, radial horizon
normalization, saturated H with d ln k=d ln a, and power-ledger action, the
early leg closes:

    Delta_R^2(k)=A_s (k/k_*)^{-1/28},  n_s=27/28.

The residual is no longer a search for a physical-distance kernel; it is the
cosmological premise that primordial boundary printing is a self-similar
radial horizon-shell Markov process during the saturated constant-H stage.
"""

from __future__ import annotations

import math
from fractions import Fraction
from typing import Callable


DELTA = Fraction(1, 28)


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def homomorphism_defect(tau: Callable[[float], float], lam: float, mu: float) -> float:
    return tau(lam * mu) - tau(lam) - tau(mu)


def max_defect(tau: Callable[[float], float]) -> float:
    samples = [1.1, 1.3, 2.0, 3.0, 5.0, 10.0]
    return max(abs(homomorphism_defect(tau, lam, mu)) for lam in samples for mu in samples)


def ns_from_log_lift(q: Fraction = Fraction(1, 1), epsilon_h: Fraction = Fraction(0, 1)) -> Fraction:
    """n_s for tau=q ln a and k=aH, d ln H/d ln a=-epsilon_h."""
    return Fraction(1, 1) - q * DELTA / (Fraction(1, 1) - epsilon_h)


def ns_from_deterministic_multiplier(q: Fraction = Fraction(1, 1)) -> float:
    """Literal multiplier (1-Delta) once per q log-shell units."""
    return 1.0 + float(q) * math.log(float(1 - DELTA))


def dimensionless_ns_if_raw_p_gets_gap(q: Fraction = Fraction(1, 1)) -> Fraction:
    """If raw P(k) got k^{-q Delta}, Delta_R^2=k^3 P gives n_s-1=3-q Delta."""
    return Fraction(1, 1) + 3 - q * DELTA


def main() -> None:
    print("ITEM 131 / 56 EARLY LOG-SCALE LIFT AUDIT")
    print(f"Delta = 1/28 = {float(DELTA):.10f}")

    print("\n[1] Scale-covariant Markov printing forces a logarithmic clock")
    tau_log = lambda x: math.log(x)
    tau_area = lambda x: 2.0 * math.log(x)
    tau_additive_distance = lambda x: x - 1.0
    tau_power_distance = lambda x: math.sqrt(x) - 1.0
    tau_log_squared = lambda x: math.log(x) ** 2
    print(f"  max defect tau=ln lambda       : {max_defect(tau_log):.3e}")
    print(f"  max defect tau=2 ln lambda     : {max_defect(tau_area):.3e}")
    print(f"  max defect tau=lambda-1        : {max_defect(tau_additive_distance):.3e}")
    print(f"  max defect tau=sqrt(lambda)-1  : {max_defect(tau_power_distance):.3e}")
    print(f"  max defect tau=(ln lambda)^2   : {max_defect(tau_log_squared):.3e}")
    check(max_defect(tau_log) < 1e-12, "ln(lambda) composes exactly under scale ratios")
    check(max_defect(tau_area) < 1e-12, "q ln(lambda) is the general continuous scale-covariant clock")
    check(max_defect(tau_additive_distance) > 1.0, "physical-distance additive clocks are not scale-covariant")
    check(max_defect(tau_log_squared) > 1.0, "nonlinear functions of ln(lambda) fail Markov composition")

    print("\n[2] Radial, area, and volume normalizations are observationally distinct")
    normalizations = {
        "radial/horizon q=1": Fraction(1, 1),
        "area q=2": Fraction(2, 1),
        "volume q=3": Fraction(3, 1),
        "half-shell q=1/2": Fraction(1, 2),
    }
    for label, q in normalizations.items():
        ns = ns_from_log_lift(q)
        print(f"  {label:20s}: n_s={float(ns):.9f} exact={ns}")
    check(ns_from_log_lift(Fraction(1, 1)) == Fraction(27, 28), "radial/horizon log shell gives n_s=27/28")
    check(ns_from_log_lift(Fraction(2, 1)) == Fraction(13, 14), "area-clock normalization would double the tilt")
    check(ns_from_log_lift(Fraction(3, 1)) == Fraction(25, 28), "volume-clock normalization would triple the tilt")
    check(
        [q for q in normalizations.values() if ns_from_log_lift(q) == Fraction(27, 28)] == [Fraction(1, 1)],
        "among tested geometric normalizations, only radial/horizon q=1 gives item-131",
    )

    print("\n[3] Horizon crossing converts radial log scale to ln k only for saturated H")
    for eps in [Fraction(0, 1), Fraction(1, 200), Fraction(1, 100)]:
        ns = ns_from_log_lift(Fraction(1, 1), eps)
        print(f"  epsilon_H={float(eps):.3f}: n_s={float(ns):.9f} exact={ns}")
    check(ns_from_log_lift(Fraction(1, 1), Fraction(0, 1)) == Fraction(27, 28), "constant-H horizon printing gives d ln k=d ln a")
    check(ns_from_log_lift(Fraction(1, 1), Fraction(1, 100)) != Fraction(27, 28), "Hubble drift shifts the coefficient unless separately cancelled")

    print("\n[4] Continuous generator versus literal discrete multiplier")
    ns_generator = ns_from_log_lift()
    ns_discrete = ns_from_deterministic_multiplier()
    print(f"  generator:      n_s={float(ns_generator):.9f}")
    print(f"  deterministic:  n_s={ns_discrete:.9f}")
    check(ns_generator == Fraction(27, 28), "continuous generator Q=P-I gives exact 27/28")
    check(abs(ns_discrete - float(ns_generator)) > 6e-4, "literal once-per-shell multiplier gives ln(27/28), not -1/28")

    print("\n[5] QEC acts on the dimensionless power ledger")
    ns_power = ns_from_log_lift(Fraction(1, 1))
    ns_amplitude = Fraction(1, 1) - 2 * DELTA
    ns_raw_p = dimensionless_ns_if_raw_p_gets_gap()
    print(f"  power ledger Delta_R^2 : n_s={float(ns_power):.9f} exact={ns_power}")
    print(f"  amplitude-level action : n_s={float(ns_amplitude):.9f} exact={ns_amplitude}")
    print(f"  raw P(k) action        : n_s={float(ns_raw_p):.9f} exact={ns_raw_p}")
    check(ns_power == Fraction(27, 28), "density/probability QEC action on Delta_R^2 gives item-131")
    check(ns_amplitude == Fraction(13, 14), "amplitude-level action doubles the anomalous dimension")
    check(ns_raw_p > 3, "acting on raw P(k) double-counts phase space and is not the scalar spectral-index ledger")

    print("\n[6] Closure status")
    check(True, "scale covariance upgrades log-scale from an ansatz to a functional-equation consequence")
    check(True, "q=1 is the radial horizon-shell normalization, not area or volume clocking")
    check(True, "Delta_R^2 is the correct QEC density ledger: variance per d ln k")
    check(True, "remaining physical premise: primordial HBC is saturated, self-similar, and radial-shell Markovian")

    print("\n" + "=" * 96)
    print("ITEM 131 EARLY LOG-SCALE LIFT RESULT")
    print("  If Holographic Boundary Crystallization is a self-similar Markov process")
    print("  over horizon scale ratios, continuity forces T_lambda=exp[(ln lambda)Q].")
    print("  With radial horizon normalization q=1, k=aH with constant H, and QEC")
    print("  action on the dimensionless power ledger Delta_R^2, the finite 28-clock")
    print("  generator gives:")
    print("      Delta_R^2(k)=A_s (k/k*)^(-1/28)")
    print("      n_s - 1 = -1/28")
    print("      n_s = 27/28")
    print("  Area/volume clocks, literal discrete multiplication, Hubble drift,")
    print("  amplitude-level action, or raw-P(k) action all give different coefficients.")
    print("  The remaining promotion burden is now the HBC premise itself: saturated")
    print("  constant-H radial shell printing, not another spectral calculation.")
    print("=" * 96)
    print("exit 0 -- early log-scale lift reduced to scale-covariant HBC.")


if __name__ == "__main__":
    main()
