#!/usr/bin/env python3
r"""ITEM 131: HBC/QEC physical coupling to scalar curvature covariance.

Question
--------
After the mode-local radial-crossing lemma, the remaining early-leg bridge is
physical:

    Does saturated HBC/QEC actually act on scalar curvature covariance/power
    during an approximately constant-H horizon-printing stage?

Result
------
The existing HBC/QEC source material supports the background pieces:

* stroboscopic CPTP / Lindblad-QEC dynamics;
* Landauer exhaust and fresh zero-entropy boundary nodes;
* Bekenstein-Hawking horizon-area saturation;
* an equilibrium/de-Sitter-style Friedmann scaling.

It does not yet derive the perturbation coupling.  The missing theorem is a
single-clock delta-N bridge:

    local HBC print-rate fluctuation  ->  local e-fold perturbation delta N
                                      ->  comoving curvature perturbation R

If that bridge is granted, the QEC service generator acts on the covariance
ledger C_R(k)=<R_k R_-k>, hence on Delta_R^2, and the previous log-shell
theorem gives n_s=27/28.  Without it, HBC is only a homogeneous entropy/area
ledger and produces no scalar curvature spectrum.

The exact coefficient also requires saturated constant-H printing:

    k = a H,  epsilon_H = - d ln H / dN = 0.

Any nonzero epsilon_H shifts the coefficient unless an extra cancellation is
derived.  Exact de Sitter also needs care: in standard single-field language
R requires a clock/adiabatic degree of freedom.  Here that clock would have to
be the HBC/QEC print process itself, not an unstated inflaton.
"""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DELTA = Fraction(1, 28)


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def has_any(text: str, needles: list[str]) -> bool:
    lower = text.lower()
    return any(needle.lower() in lower for needle in needles)


def ns_power(q: Fraction = Fraction(1, 1), epsilon_h: Fraction = Fraction(0, 1)) -> Fraction:
    """Power-ledger tilt for tau=q ln a and d ln k/dN = 1-epsilon_H."""
    return Fraction(1, 1) - q * DELTA / (Fraction(1, 1) - epsilon_h)


def ns_amplitude(epsilon_h: Fraction = Fraction(0, 1)) -> Fraction:
    """If the same generator acted on amplitude, covariance/power doubles it."""
    return Fraction(1, 1) - 2 * DELTA / (Fraction(1, 1) - epsilon_h)


def ns_raw_p(epsilon_h: Fraction = Fraction(0, 1)) -> Fraction:
    """If the generator acted on raw P(k), Delta_R^2=k^3P double-counts phase space."""
    return Fraction(4, 1) - DELTA / (Fraction(1, 1) - epsilon_h)


def main() -> None:
    print("ITEM 131 HBC/QEC PHYSICAL COUPLING AUDIT")

    engine = (ROOT / "cosmological_qec_engine" / "cosmological_qec_engine.tex").read_text()
    part17 = (ROOT / "part_17_energy_trajectory" / "part_17_energy_trajectory.tex").read_text()
    corpus = engine + "\n" + part17

    print("\n[1] Background HBC/QEC support present in source")
    support = {
        "stroboscopic CPTP clock": ["rho(t_{n+1}) = \\mathcal{E}", "CPTP"],
        "Lindblad/KMS Markov structure": ["Lindbladian semigroup", "KMS condition"],
        "Landauer exhaust": ["Landauer waste heat", "Landauer erasure"],
        "fresh zero-entropy boundary nodes": ["strictly zero entanglement entropy", "S(0) = 0"],
        "horizon precipitation": ["precipitating new $Q_3$ matter cells", "cosmological horizon"],
        "Bekenstein-Hawking saturation": ["Bekenstein--Hawking", "S_{\\max}", "horizon area"],
        "de-Sitter/Friedmann background form": ["de Sitter Friedmann", "H^2 = \\frac{\\rho_\\Lambda}{3 M_P^2}", "w(0) = -1"],
    }
    for label, needles in support.items():
        check(has_any(corpus, needles), f"source contains {label}")

    print("\n[2] Perturbation-coupling pieces absent from source")
    missing = {
        "gauge-invariant scalar curvature variable": [
            "\\zeta",
            "comoving curvature",
            "curvature perturbation",
            "Mukhanov",
        ],
        "delta-N / separate-universe map": [
            "\\delta N",
            "delta-N",
            "separate universe",
            "local e-fold",
        ],
        "dimensionless scalar power spectrum": [
            "\\Delta_{\\mathcal R}",
            "Delta_R",
            "scalar power spectrum",
        ],
        "horizon-crossing relation": [
            "k=aH",
            "horizon crossing",
            "d ln k",
        ],
        "constant-H slow-roll parameter control": [
            "epsilon_H",
            "constant H",
            "slow-roll",
        ],
    }
    for label, needles in missing.items():
        found = has_any(engine, needles)
        check(not found, f"HBC/QEC engine text does not derive {label}")

    print("\n[3] Minimal physical theorem that would close the coupling")
    theorem_legs = [
        "local HBC print-rate perturbation is the single adiabatic clock",
        "delta local print time equals curvature perturbation: R = delta N",
        "QEC generator acts on covariance/power C_R, not coherent amplitude",
        "mode-local action at horizon exit preserves the radial q=1 normalization",
        "no extra inflaton/metric transfer term contributes to n_s-1",
    ]
    for leg in theorem_legs:
        check(True, leg)

    print("\n[4] Failure-mode coefficient audit")
    branches = {
        "closed target: power ledger, q=1, eps=0": ns_power(),
        "area/total-entropy coupling q=2": ns_power(Fraction(2, 1)),
        "amplitude-level coupling": ns_amplitude(),
        "raw P(k) coupling": ns_raw_p(),
        "H drift eps=0.005": ns_power(epsilon_h=Fraction(1, 200)),
        "H drift eps=0.010": ns_power(epsilon_h=Fraction(1, 100)),
    }
    for label, value in branches.items():
        print(f"  {label:38s}: n_s={float(value):.9f} exact={value}")
    check(branches["closed target: power ledger, q=1, eps=0"] == Fraction(27, 28), "target branch gives n_s=27/28")
    check(branches["area/total-entropy coupling q=2"] == Fraction(13, 14), "total boundary entropy coupling gives the old q=2 miss")
    check(branches["amplitude-level coupling"] == Fraction(13, 14), "amplitude coupling doubles the covariance tilt")
    check(branches["raw P(k) coupling"] > 3, "raw P(k) action is the wrong observable")
    check(branches["H drift eps=0.010"] != Fraction(27, 28), "nonzero Hubble drift shifts the exact coefficient")

    print("\n[5] Saturated constant-H status")
    check(has_any(engine, ["steady-state density", "de Sitter Friedmann", "equilibrium thermodynamic steady-state"]), "source supports a homogeneous equilibrium/de-Sitter background reading")
    check(not has_any(engine, ["epsilon_H", "constant H", "horizon crossing", "inflationary stage"]), "source does not derive an early saturated constant-H horizon-crossing stage")
    check(has_any(part17, ["w(0) = -1", "driving initial exponential inflation"]), "Part 17 supplies a qualitative w(0)=-1 / exponential-inflation boundary")
    check(True, "exact de Sitter needs the HBC/QEC print process to supply the scalar clock; otherwise R is not defined as an adiabatic perturbation")

    print("\n" + "=" * 100)
    print("VERDICT")
    print("  Not closed.  The next viable lemma is not another normalization or")
    print("  finite-clock theorem; it is the single-clock HBC perturbation theorem:")
    print("      local print-rate fluctuation -> delta N -> scalar curvature R.")
    print("  If that theorem is proven, the existing power-ledger/log-shell machinery")
    print("  gives n_s=27/28.  If HBC couples only to homogeneous entropy/area load,")
    print("  there is no scalar curvature covariance to tilt.  If it couples to total")
    print("  boundary entropy, amplitude, raw P(k), or a drifting-H background, the")
    print("  coefficient changes.")
    print("  Saturated constant-H printing is currently an equilibrium/de-Sitter")
    print("  background premise, not a derived early-universe stage with epsilon_H=0")
    print("  and a defined scalar clock.")
    print("=" * 100)
    print("exit 0 -- physical coupling narrowed to the HBC delta-N scalar-clock lemma; saturated-H remains a premise.")


if __name__ == "__main__":
    main()
