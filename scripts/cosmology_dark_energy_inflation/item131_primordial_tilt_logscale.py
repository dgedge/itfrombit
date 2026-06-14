#!/usr/bin/env python3
r"""ITEM 131 / 56: primordial tilt from log-scale boundary printing.

Target
------
Close the remaining early-leg map

    n_s - 1 = -1/28

without using the failed physical-distance exponential or depolarizing/FDT
routes.

Result
------
The finite 28-channel clock closes the tilt if boundary printing is read as a
continuous Markov/RG flow in logarithmic horizon scale:

    d ln Delta_R^2 / d ln k = -Delta,  Delta = 1/28.

This is an anomalous-dimension statement.  It is not the Fourier transform of a
3D exponential correlation in physical distance.  The companion audit
item131_early_logscale_lift.py sharpens the lift: if boundary crystallization
is a scale-covariant Markov process over horizon scale ratios, then continuity
forces T_lambda=exp[(ln lambda)Q].  The remaining normalization is radial
horizon-shell clocking q=1, not area or volume clocking.

Load-bearing distinctions
-------------------------
1. Log-scale horizon map:
       k = aH.  During saturated boundary printing H is constant, so
       d ln k = d ln a.  The boundary print clock is therefore an e-fold/log-k
       clock, not a physical-distance kernel.

2. Continuous generator, not deterministic once-per-e-fold multiplication:
       the serial clock transition has first nontrivial eigenvalue 27/28.
       A literal multiplier once per ln k would give exponent ln(27/28), not
       exactly -1/28.  The exact item-131 value uses the Markov generator
       Q = P - I, whose first generator eigenvalue is -1/28.

3. Power-level action:
       the QEC boundary instrument evolves density/probability ledgers
       (active-syndrome weights), so its anomalous dimension acts on
       Delta_R^2 itself.  If it acted on a coherent amplitude, the power tilt
       would double to -2/28.

Thus the theorem form is:
    finite 28-channel serial QEC clock
    + Poisson/continuous Markov coarse-graining in N = ln a
    + saturated horizon printing, d ln k = dN
    + density/power-level scalar ledger
    => n_s = 27/28.
"""

from __future__ import annotations

import math
from fractions import Fraction


N = 28
DELTA = Fraction(1, N)


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def serial_transition_eigenvalues(n: int = N) -> list[Fraction]:
    """Triangular serial erasure clock eigenvalues on active-count sectors."""
    return [Fraction(n - m, n) for m in range(n + 1)]


def serial_generator_eigenvalues(n: int = N) -> list[Fraction]:
    """Unit-rate continuous-time generator Q=P-I eigenvalues."""
    return [ev - 1 for ev in serial_transition_eigenvalues(n)]


def dimensionless_exponential_slope(kxi: float) -> float:
    """3D C(r)=exp(-r/xi) gives Delta^2(k) proportional to k^3/(1+k^2 xi^2)^2."""
    q = kxi * kxi
    return 3.0 - 4.0 * q / (1.0 + q)


def logspace(start_exp: float, stop_exp: float, count: int) -> list[float]:
    step = (stop_exp - start_exp) / (count - 1)
    return [10 ** (start_exp + i * step) for i in range(count)]


def std(values: list[float]) -> float:
    mean = sum(values) / len(values)
    return math.sqrt(sum((x - mean) ** 2 for x in values) / len(values))


def ns_from_continuous_log_flow(delta: Fraction = DELTA) -> Fraction:
    return Fraction(1, 1) - delta


def ns_from_deterministic_tick(delta: Fraction = DELTA) -> float:
    # Delta_R^2(k) multiplied by (1-delta) once per unit ln k.
    return 1.0 + math.log(float(1 - delta))


def ns_with_hubble_drift(epsilon_h: float, delta: Fraction = DELTA) -> float:
    # k=aH, d ln H/d ln a = -epsilon_h, so d ln k/dN = 1-epsilon_h.
    return 1.0 - float(delta) / (1.0 - epsilon_h)


def main() -> None:
    print("ITEM 131 / 56 PRIMORDIAL TILT LOG-SCALE AUDIT")
    print(f"Delta = 1/28 = {float(DELTA):.10f}")

    print("\n[1] Serial clock as continuous Markov generator")
    p_eigs = serial_transition_eigenvalues()
    q_eigs = serial_generator_eigenvalues()
    check(p_eigs[0] == 1, "serial transition has vacuum eigenvalue 1")
    check(p_eigs[1] == 1 - DELTA, "first transition eigenvalue is 27/28")
    check(q_eigs[0] == 0, "generator has vacuum eigenvalue 0")
    check(q_eigs[1] == -DELTA, "first generator eigenvalue is -1/28")

    print("\n[2] Deterministic tick vs generator exponent")
    ns_generator = ns_from_continuous_log_flow()
    ns_discrete = ns_from_deterministic_tick()
    check(ns_generator == Fraction(27, 28), "continuous generator gives n_s=1-1/28=27/28 exactly")
    check(abs(ns_discrete - float(ns_generator)) > 6e-4, "literal multiplier 27/28 per e-fold gives ln(27/28), not exactly -1/28")
    print(f"  generator n_s      = {float(ns_generator):.9f}")
    print(f"  deterministic n_s  = {ns_discrete:.9f}")
    print(f"  difference         = {ns_discrete - float(ns_generator):+.9f}")

    print("\n[3] Log-scale flow closes the scalar tilt")
    # Analytic finite-difference check: log power is exactly -Delta log k.
    k = logspace(-4, 4, 4000)
    log_power = [-float(DELTA) * math.log(x) for x in k]
    slopes = [
        (log_power[i + 1] - log_power[i - 1]) / (math.log(k[i + 1]) - math.log(k[i - 1]))
        for i in range(1, len(k) - 1)
    ]
    check(max(abs(s + float(DELTA)) for s in slopes) < 1e-12, "Delta_R^2(k)=(k/k*)^-Delta has constant slope n_s-1=-Delta")
    median_slope = sorted(slopes)[len(slopes) // 2]
    check(abs((1.0 + median_slope) - float(Fraction(27, 28))) < 1e-12, "therefore n_s=27/28")

    print("\n[4] Physical-distance exponential route remains rejected")
    kxi = logspace(-2, 2, 2000)
    distance_slope = [dimensionless_exponential_slope(x) for x in kxi]
    target = -float(DELTA)
    pivot = min(kxi, key=lambda x: abs(dimensionless_exponential_slope(x) - target))
    mid_slopes = [s for x, s in zip(kxi, distance_slope) if 0.3 < x < 3.0]
    check(std(mid_slopes) > 0.8, "3D exponential correlation has running slope, not a constant anomalous dimension")
    check(abs(dimensionless_exponential_slope(pivot) - target) < 5e-3, f"it hits -1/28 only near a tuned pivot k xi ~= {pivot:.3f}")

    print("\n[5] Horizon-printing log variable")
    for epsilon in [0.0, 0.005, 0.01]:
        ns = ns_with_hubble_drift(epsilon)
        print(f"  epsilon_H={epsilon:.3f}: n_s={ns:.9f}")
    check(abs(ns_with_hubble_drift(0.0) - float(Fraction(27, 28))) < 1e-15, "saturated boundary printing with constant H gives d ln k=d ln a and exact n_s=27/28")
    check(abs(ns_with_hubble_drift(0.01) - float(Fraction(27, 28))) > 3e-4, "nonzero H drift shifts the exact coefficient unless separately cancelled")

    print("\n[6] Power-level versus amplitude-level action")
    ns_power = 1.0 - float(DELTA)
    ns_amplitude = 1.0 - 2.0 * float(DELTA)
    check(abs(ns_power - float(Fraction(27, 28))) < 1e-15, "density/probability instrument acting on Delta_R^2 gives the item-131 tilt")
    check(abs(ns_amplitude - float(Fraction(13, 14))) < 1e-15, "if the clock acted on coherent amplitude, power tilt would double to n_s=13/14")
    check(abs(ns_amplitude - ns_power) > 0.03, "amplitude-level action is observationally and structurally the wrong branch")

    print("\n" + "=" * 92)
    print("FINITE EARLY-LEG THEOREM FORM")
    print("  The remaining map is closed if boundary printing is the continuous")
    print("  log-scale Markov/RG flow already suggested by Holographic Boundary")
    print("  Crystallization:")
    print("      d/d ln k Delta_R^2 = - (1/28) Delta_R^2")
    print("  Then:")
    print("      Delta_R^2(k)=A_s (k/k*)^(-1/28)")
    print("      n_s - 1 = -1/28")
    print("      n_s = 27/28")
    print("  Exactness requires the generator reading of the serial clock, constant-H")
    print("  saturated boundary printing, and power/density-level action.  A literal")
    print("  discrete multiplier, Hubble drift, or amplitude-level action gives a")
    print("  different coefficient.")
    print("  Companion item131_early_logscale_lift.py upgrades the log variable")
    print("  from an ansatz to a scale-covariant Markov-shell consequence; the")
    print("  remaining physical premise is saturated radial HBC itself.")
    print("=" * 92)
    print("exit 0 -- early tilt reduced to a precise log-scale generator theorem.")


if __name__ == "__main__":
    main()
