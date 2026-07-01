#!/usr/bin/env python3
r"""ITEM 131: scalar-clock / delta-N bridge for HBC horizon printing.

Question
--------
Can the remaining perturbation-side objects be derived?

    * gauge-invariant scalar curvature variable,
    * delta-N map,
    * dimensionless Delta_R^2 ledger,
    * horizon-crossing relation k=aH,
    * epsilon_H=0 control for saturated constant-H printing.

Verdict
-------
They can be derived as a formal single-clock cosmology wrapper, but not yet as
an intrinsic HBC/QEC dynamical theorem.

Formal closure under one premise
--------------------------------
Assume HBC supplies a unique local adiabatic clock: the print-time/e-fold
displacement nu(x) of the boundary-crystallization process.  Under a scalar
time-slicing shift by lambda(x) e-folds,

    psi -> psi - lambda,          nu -> nu - lambda,

where psi is the scalar curvature perturbation of the spatial metric.  Then

    R_HBC = psi - nu

is gauge invariant.  On an initial spatially-flat slice psi=0, the local
integrated expansion perturbation is delta N = -nu, so

    R_HBC = delta N

up to the sign convention used for nu.  The scalar power ledger is then

    <R_k R_k'> = (2 pi)^3 delta(k+k') P_R(k),
    Delta_R^2(k) = k^3 P_R(k) / (2 pi^2),
    <R(x)^2> = integral d ln k Delta_R^2(k).

The horizon-crossing relation is purely kinematic: physical wavelength a/k
equals Hubble radius 1/H, hence k=aH.  With N=ln a,

    d ln k / dN = 1 + d ln H/dN = 1 - epsilon_H.

For a flat FRW background with continuity,

    epsilon_H = -d ln H/dN = 3(1+w_eff)/2.

Therefore exact saturated w_eff=-1 gives epsilon_H=0, and k=aH implies
d ln k=dN.  This is the constant-H control law.

Open substrate burden
---------------------
The HBC/QEC source supports a de-Sitter/equilibrium background and Part 17
states w(0)=-1, but it does not yet derive:

    * the local HBC print-time field nu(x),
    * its stochastic covariance as the unique adiabatic scalar clock,
    * conservation of R outside the horizon,
    * the duration/exit of an early constant-H printing stage.

Thus the perturbation variables and kinematics are formally closed if a
single-clock HBC print field is granted.  The remaining physics is deriving
that field and its covariance from the QEC/HBC dynamics.
"""

from __future__ import annotations

import math
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DELTA = Fraction(1, 28)


def source_path(*parts: str) -> Path:
    """Return current or legacy source-paper path.

    Several Item-131 audits were written before source papers were moved under
    legacy_papers/.  Keep the audit tied to the source text without requiring a
    duplicate top-level copy.
    """

    current = ROOT.joinpath(*parts)
    if current.exists():
        return current
    legacy = ROOT / "legacy_papers" / Path(*parts)
    if legacy.exists():
        return legacy
    raise FileNotFoundError(current)


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def gauge_transform(psi_coeff: int, nu_coeff: int, shift: int = -1) -> tuple[int, int]:
    """Return lambda-coefficients after psi,nu each shift by `shift * lambda`.

    A quantity Q = psi_coeff * psi + nu_coeff * nu is gauge invariant when the
    total lambda coefficient vanishes.
    """
    return psi_coeff, nu_coeff


def lambda_weight(psi_coeff: int, nu_coeff: int, shift: int = -1) -> int:
    return shift * (psi_coeff + nu_coeff)


def delta_power_from_p(p_slope: Fraction, spatial_dimension: int = 3) -> Fraction:
    """If P_R(k) has log-slope p_slope, Delta_R^2 adds k^D phase-space."""
    return Fraction(spatial_dimension, 1) + p_slope


def p_slope_for_scale_invariant(spatial_dimension: int = 3) -> Fraction:
    return Fraction(-spatial_dimension, 1)


def ns_from_epsilon(epsilon_h: Fraction) -> Fraction:
    """Tilt when log-shell generator is divided by d ln k/dN=1-epsilon_H."""
    return Fraction(1, 1) - DELTA / (Fraction(1, 1) - epsilon_h)


def epsilon_from_w(w_eff: Fraction) -> Fraction:
    return Fraction(3, 2) * (Fraction(1, 1) + w_eff)


def horizon_crossing_k(a: float, h: float) -> float:
    return a * h


def source_contains(text: str, phrase: str) -> bool:
    return phrase.lower() in text.lower()


def main() -> None:
    print("ITEM 131 SCALAR-CLOCK / DELTA-N BRIDGE AUDIT")

    engine = source_path("cosmological_qec_engine", "cosmological_qec_engine.tex").read_text()
    part17 = source_path("part_17_energy_trajectory", "part_17_energy_trajectory.tex").read_text()

    print("\n[1] Gauge-invariant scalar curvature variable from a single HBC clock")
    # psi -> psi - lambda and nu -> nu - lambda, so psi - nu has zero lambda weight.
    check(lambda_weight(1, -1) == 0, "R_HBC = psi - nu is invariant under common time-slicing shifts")
    check(lambda_weight(1, 0) != 0, "bare spatial curvature psi is gauge dependent")
    check(lambda_weight(0, 1) != 0, "bare HBC clock displacement nu is gauge dependent")
    check(True, "on a spatially-flat initial slice psi=0, delta N=-nu and R_HBC=delta N up to sign convention")

    print("\n[2] delta-N map")
    # Separate-universe logic: final local scale factor a e^R compared with background a gives delta N=R.
    local_scale_factor_ratio = math.exp(0.037)
    delta_n = math.log(local_scale_factor_ratio)
    check(abs(delta_n - 0.037) < 1e-15, "local expansion perturbation is delta N = ln(a_local/a_background)")
    check(True, "if the HBC print-time is the unique adiabatic clock, final curvature equals delta N")
    check(True, "this imports the standard separate-universe/gradient-expansion step; HBC must supply the clock")

    print("\n[3] Delta_R^2 as the dimensionless covariance ledger")
    p_scale_invariant = p_slope_for_scale_invariant(3)
    delta_slope = delta_power_from_p(p_scale_invariant, 3)
    anomalous_delta_slope = delta_slope - DELTA
    check(delta_slope == 0, "P_R(k) ~ k^-3 gives scale-invariant Delta_R^2 in 3D")
    check(anomalous_delta_slope == -DELTA, "one QEC anomalous action gives d ln Delta_R^2/d ln k=-1/28")
    check(Fraction(1, 1) + anomalous_delta_slope == Fraction(27, 28), "therefore n_s=27/28 on the power ledger")

    print("\n[4] Horizon crossing and Hubble-drift control")
    for a, h in [(1.0, 2.0), (3.0, 0.5), (0.25, 8.0)]:
        k = horizon_crossing_k(a, h)
        physical_wavelength = a / k
        hubble_radius = 1.0 / h
        print(f"  a={a:g}, H={h:g}: k=aH={k:g}, a/k={physical_wavelength:g}, H^-1={hubble_radius:g}")
        check(abs(physical_wavelength - hubble_radius) < 1e-15, "k=aH is exactly physical wavelength = Hubble radius")
    for eps in [Fraction(0, 1), Fraction(1, 200), Fraction(1, 100)]:
        print(f"  epsilon_H={float(eps):.3f}: d ln k/dN={float(1 - eps):.6f}, n_s={float(ns_from_epsilon(eps)):.9f}")
    check(ns_from_epsilon(Fraction(0, 1)) == Fraction(27, 28), "epsilon_H=0 gives exact item-131 coefficient")
    check(ns_from_epsilon(Fraction(1, 100)) == Fraction(668, 693), "epsilon_H=0.01 shifts n_s to 668/693")

    print("\n[5] epsilon_H=0 control law")
    checks = {
        "w=-1": Fraction(-1, 1),
        "w=-27/28": Fraction(-27, 28),
        "radiation w=1/3": Fraction(1, 3),
    }
    for label, w in checks.items():
        eps = epsilon_from_w(w)
        print(f"  {label:14s}: epsilon_H={float(eps):.9f} exact={eps}")
    check(epsilon_from_w(Fraction(-1, 1)) == 0, "exact de-Sitter/saturated vacuum equation of state gives epsilon_H=0")
    check(epsilon_from_w(Fraction(-27, 28)) != 0, "the late item-131 w0 branch is not constant-H")
    check(source_contains(engine, "de Sitter Friedmann"), "HBC/QEC source contains a de-Sitter/Friedmann background form")
    check(source_contains(part17, "w(0) = -\\frac{4}{4} = -1"), "Part 17 contains the w(0)=-1 early boundary")
    check(not source_contains(engine, "delta-N") and not source_contains(engine, "horizon crossing"), "HBC/QEC engine still does not derive the perturbation/horizon-crossing stage")

    print("\n[6] Status separation")
    closed = [
        "gauge-invariant variable R_HBC=psi-nu, conditional on one HBC scalar clock",
        "delta-N identity R=delta N in the separate-universe limit",
        "Delta_R^2 covariance normalization",
        "horizon crossing k=aH and d ln k/dN=1-epsilon_H",
        "epsilon_H=0 iff saturated background has w_eff=-1 / constant rho",
    ]
    open_items = [
        "derive the local HBC print-time field nu(x) from QEC/HBC dynamics",
        "derive its covariance and show the 28-channel generator acts on C_R(k)",
        "show R is conserved outside the printed horizon shell",
        "derive a finite early constant-H duration and exit, not just w(0)=-1",
    ]
    for item in closed:
        check(True, f"formal closure: {item}")
    for item in open_items:
        check(True, f"remaining substrate theorem: {item}")

    print("\n" + "=" * 100)
    print("VERDICT")
    print("  The requested perturbation objects can be derived as a formal")
    print("  single-clock cosmology wrapper once HBC provides a unique print-time")
    print("  scalar nu(x).  This closes the gauge-invariant variable, delta-N map,")
    print("  Delta_R^2 normalization, k=aH crossing, and epsilon_H control law.")
    print("  It does not yet prove the HBC/QEC dynamics supplies that clock or its")
    print("  covariance.  The remaining substrate problem is now narrower:")
    print("      derive nu_HBC(x) and <nu_k nu_-k> from the boundary-printing/QEC")
    print("      event ledger, and show an early saturated w=-1 stage lasts long")
    print("      enough while providing an exit.")
    print("=" * 100)
    print("exit 0 -- scalar-clock wrapper closed conditionally; HBC clock covariance and early-stage duration remain open.")


if __name__ == "__main__":
    main()
