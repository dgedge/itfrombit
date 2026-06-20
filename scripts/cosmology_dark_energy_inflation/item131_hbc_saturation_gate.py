#!/usr/bin/env python3
r"""ITEM 131: HBC saturation gate for tilt and amplitude.

This script tests the strongest version of the current Item-131 question:

    Can the existing HBC/QEC printer dynamics derive both

        N_shell alpha_0^4 = C_F = 4/3
        S_j(k=aH) = 1

    without adding a new maximal-throughput / no-horizon-covariance axiom?

Result
------
Not yet.  The existing canon proves a useful conditional theorem but not a
Locked theorem.

What is forced:

  * The finite HBC service instrument gives a serial 28-clock, so the tilt
    generator is fixed once HBC is a saturated radial log-shell flow.
  * The post-decoder topology-changing current is weight 4, hence its event
    probability is alpha_0^4.
  * The colour-restoring load selected by the admitted channels is C_F=4/3.
  * HBC entropy/area saturation fixes the total printed entropy capacity.

What is not forced:

  * The topology-changing colour-load ledger can run below capacity while the
    entropy/area ledger remains saturated, unless one proves that every printed
    entropy unit is a colour-restoring topological commit or proves a backlog
    latch that drives the colour current to marginality.
  * The nonzero scalar Fourier mode has S_j=1 under local/exchangeable service,
    but a connected intensity covariance at the horizon wave number changes the
    amplitude while leaving the 28 tilt and the homogeneous entropy ledger
    untouched.

Thus saturation has a precise theorem target:

    HBC/QEC saturation theorem =
      (A) all scalar-relevant printed entropy is carried by the
          post-decoder colour-restoring weight-4 current;
      (B) coherent constant-H printing is a finite critical latch, so
          lambda_shell=N_shell alpha_0^4 reaches C_F and exits for overload;
      (C) the service-current lift is product/exchangeable after the
          homogeneous shell clock is removed, so no connected covariance exists
          at k=aH.

Under A-C, Item 131 closes:

    A_nu = (3/4) alpha_0^4,   n_s = 27/28.

Without A-C, the present axioms admit explicit countermodels.
"""

from __future__ import annotations

import cmath
import math
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ALPHA0 = 1.0 / 137.0
C_F = 4.0 / 3.0
DELTA = Fraction(1, 28)


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def contains(path: str, phrase: str) -> bool:
    return phrase in (ROOT / path).read_text(encoding="utf-8")


def amplitude_from_colour_utilization(utilization: float, s_j: float = 1.0, f_eff: float = 1.0) -> float:
    """A_nu = F_eff S_j / N_shell with N_shell alpha^4 = utilization*C_F."""
    n_shell = utilization * C_F / ALPHA0**4
    return f_eff * s_j / n_shell


def ns_from_28_clock() -> Fraction:
    return Fraction(1, 1) - DELTA


@dataclass(frozen=True)
class Countermodel:
    name: str
    entropy_saturated: bool
    colour_utilization: float
    horizon_covariance_sj: float
    why_allowed: str

    @property
    def amplitude(self) -> float:
        return amplitude_from_colour_utilization(self.colour_utilization, self.horizon_covariance_sj)


Vector = tuple[int, int, int]


def torus_points(n: int) -> list[Vector]:
    return list(product(range(n), repeat=3))


def dot(a: Vector, b: Vector) -> int:
    return sum(x * y for x, y in zip(a, b))


def phase(n: int, k: Vector, x: Vector) -> complex:
    return cmath.exp(-2j * math.pi * dot(k, x) / n)


def mode_structure_factor_with_horizon_covariance(
    n: int,
    k: Vector,
    mean_per_site: float,
    mode_covariance_strength: float,
) -> float:
    """S_j for Poisson service plus a coherent stochastic rate component at k."""
    pts = torus_points(n)
    m = len(pts)

    def profile(x: Vector) -> float:
        return phase(n, k, x).real

    total = 0.0 + 0.0j
    for x in pts:
        for y in pts:
            poisson = mean_per_site if x == y else 0.0
            coherent = mode_covariance_strength * mean_per_site**2 * profile(x) * profile(y)
            cov = poisson + coherent
            total += phase(n, k, x) * phase(n, k, y).conjugate() * cov
    variance = float((total / (m**2 * mean_per_site**2)).real)
    n_shell = m * mean_per_site
    return variance * n_shell


def main() -> None:
    print("ITEM 131 HBC SATURATION GATE")

    print("\n[1] Source gates already in canon")
    source_checks = [
        (
            "28-clock finite instrument",
            "python_code/item131_w_to_28_instrument.py",
            "service channel has rate 1/28",
        ),
        (
            "mode-local radial crossing",
            "python_code/item131_mode_local_radial_crossing.py",
            "n_s = 27/28",
        ),
        (
            "C_F colour-restoring load",
            "python_code/item131_cf_stop_rule_closure_attempt.py",
            "2*(8/12)=4/3",
        ),
        (
            "current T5b structure factor",
            "python_code/item131_t5b_whiteness_lemma.py",
            "S_j(k=aH)=1",
        ),
        (
            "saturated-printer residual",
            "python_code/item131_saturated_printer_dynamics_attempt.py",
            "capacity inequality lambda_shell <= C_F",
        ),
        (
            "HBC entropy/area saturation",
            "cosmological_qec_engine/cosmological_qec_engine.tex",
            "Bekenstein--Hawking entropy bound saturation against the incoming Landauer flux",
        ),
    ]
    for label, path, phrase in source_checks:
        check(contains(path, phrase), label)

    print("\n[2] Conditional closure if the saturation theorem is supplied")
    a_sat = amplitude_from_colour_utilization(1.0)
    print(f"  n_s                  = {ns_from_28_clock()} = {float(ns_from_28_clock()):.12f}")
    print(f"  N_shell              = C_F alpha0^-4 = {C_F / ALPHA0**4:.6e}")
    print(f"  A_nu                 = (3/4) alpha0^4 = {a_sat:.12e}")
    check(abs(a_sat - 0.75 * ALPHA0**4) < 1e-24, "saturated colour load gives A_nu=(3/4)alpha0^4")

    print("\n[3] Countermodels admitted by the current axioms")
    countermodels = [
        Countermodel(
            name="entropy-saturated, half colour-load",
            entropy_saturated=True,
            colour_utilization=0.5,
            horizon_covariance_sj=1.0,
            why_allowed="entropy capacity can be filled by non-topological print capacity unless the ledger bridge is proved",
        ),
        Countermodel(
            name="entropy-saturated, near colour-load",
            entropy_saturated=True,
            colour_utilization=0.9,
            horizon_covariance_sj=1.0,
            why_allowed="ordinary stable queue utilisation below one is not excluded by HBC entropy saturation alone",
        ),
        Countermodel(
            name="colour-saturated, horizon covariance",
            entropy_saturated=True,
            colour_utilization=1.0,
            horizon_covariance_sj=1.2,
            why_allowed="a k=aH service-rate modulation preserves homogeneous entropy saturation but changes Pi_k variance",
        ),
    ]
    for model in countermodels:
        ratio = model.amplitude / a_sat
        print(f"  {model.name:38s} A/A_sat={ratio:7.3f}  {model.why_allowed}")
        check(model.entropy_saturated, f"{model.name}: HBC entropy saturation is still satisfied by construction")
        check(abs(model.amplitude - a_sat) > 1.0e-12, f"{model.name}: amplitude differs from the target")
    print("  These are not proposed physics.  They are independence witnesses:")
    print("  the current axioms do not logically force the Item-131 amplitude.")

    print("\n[4] Explicit horizon-covariance danger")
    n = 9
    k = (1, 0, 0)
    mean = 17.0
    s_white = mode_structure_factor_with_horizon_covariance(n, k, mean, 0.0)
    s_modulated = mode_structure_factor_with_horizon_covariance(n, k, mean, 0.02)
    print(f"  product-local current:             S_j={s_white:.12f}")
    print(f"  2% stochastic horizon-mode current: S_j={s_modulated:.12f}")
    check(abs(s_white - 1.0) < 1.0e-10, "local Poisson service gives S_j=1")
    check(s_modulated > 2.0, "a small coherent horizon-rate component materially changes S_j")
    print("  Therefore 'no horizon-scale covariance' is not cosmetic; it is part of")
    print("  the amplitude theorem.")

    print("\n[5] The finite latch theorem that would close Item 131")
    clauses = [
        (
            "ledger identity",
            "all scalar-relevant printed entropy is post-decoder colour-restoring weight-4 current",
            "kills subcritical colour-utilisation countermodels",
        ),
        (
            "critical latch",
            "finite saturated constant-H printing runs at marginal capacity lambda_shell=C_F",
            "turns the capacity inequality into equality without an adjustable utilization",
        ),
        (
            "exit law",
            "lambda_shell>C_F destroys coherent printing or triggers reheating/handoff",
            "makes equality a finite latch rather than an eternal unstable M/M/1 queue",
        ),
        (
            "spatial whiteness",
            "after removing the homogeneous shell clock, service events are local/exchangeable",
            "sets S_j(k=aH)=1 and excludes the horizon-covariance countermodel",
        ),
    ]
    for name, theorem, consequence in clauses:
        print(f"  {name:18s}: {theorem}")
        print(f"      -> {consequence}")

    print("\n[6] Promotion table")
    rows = [
        ("tilt finite clock", "CLOSED", "28-channel serial service gives the finite generator"),
        ("tilt cosmological lift", "CONDITIONAL", "requires saturated radial log-shell Markov HBC and constant H"),
        ("C_F coefficient", "SUBSTANTIALLY CLOSED", "post-filtered colour-silent channels give 2*(8/12)=4/3"),
        ("N_shell equality", "OPEN", "needs ledger identity + critical latch"),
        ("S_j(k=aH)", "OPEN", "needs no connected horizon-scale covariance"),
        ("A_nu", "CONDITIONAL", "lands at (3/4)alpha0^4 only if the two open clauses close"),
    ]
    for item, status, note in rows:
        print(f"  {item:24s} {status:22s} {note}")

    print("\nVERDICT")
    print("  The programme is now very narrow.  Existing HBC/QEC mechanics do not yet")
    print("  derive the full saturation premise.  They derive the finite 28-clock, the")
    print("  alpha^4 topology-current power, and the C_F capacity/load.  What remains")
    print("  is one finite latch theorem plus one spatial-whiteness theorem.  If both")
    print("  are proved, tilt and amplitude close together; if either fails, Item 131")
    print("  stays a conditional recovery candidate.")
    print("exit 0 -- HBC saturation reduced to ledger-identity/critical-latch/no-horizon-covariance gates.")


if __name__ == "__main__":
    main()
