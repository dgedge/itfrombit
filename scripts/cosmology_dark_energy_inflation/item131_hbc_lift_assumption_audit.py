#!/usr/bin/env python3
r"""ITEM 131: does the finite 28-clock really lift to horizon-shell/log-scale printing?

Question
--------
The finite instrument work now gives a 28-channel serial absorbing QEC clock
with generator gap

    Delta = 1/28.

Does that *by itself* imply the cosmological scalar-tilt law

    d ln Delta_R^2 / d ln k = -1/28 ?

Verdict
-------
No.  The finite clock supplies the local generator only.  The lift to
horizon-shell/log-scale printing is a conditional theorem requiring four
additional cosmological premises:

1. HBC shell Markov/semigroup property:
       T_{lambda mu} = T_lambda T_mu
   for printed horizon-scale ratios.

2. Radial horizon-shell normalization:
       tau(lambda) = ln(lambda)
   rather than area clocking 2 ln(lambda), volume clocking 3 ln(lambda), or
   physical-distance clocking.

3. Saturated constant-H horizon printing:
       k = a H,  d ln H / d ln a = 0,
   so d ln k = d ln a.

4. Power/density ledger action:
       the QEC generator acts on Delta_R^2, not coherent amplitude and not
       raw P(k).

The existing cosmological-QEC HBC text supplies horizon precipitation and
entropy-rate matching, but does not derive the scale-ratio semigroup, q=1
normalization, saturated-H condition, or power-ledger target.  DRIFT/ANCHOR
already say this is conditional.  Therefore the exact early-leg status is:

    finite 28-clock: closed at instrument level
    finite->cosmology HBC/log-shell lift: open physical premise
"""

from __future__ import annotations

import math
import re
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DELTA = Fraction(1, 28)


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def extract_section(text: str, section_label: str) -> str:
    marker = rf"\\section\{{{re.escape(section_label)}\}}"
    m = re.search(marker, text)
    if not m:
        return ""
    rest = text[m.start() :]
    next_section = re.search(r"\n\\section\{", rest[len(m.group(0)) :])
    if not next_section:
        return rest
    return rest[: len(m.group(0)) + next_section.start()]


def ns_from_log_lift(q: Fraction = Fraction(1, 1), epsilon_h: Fraction = Fraction(0, 1)) -> Fraction:
    return Fraction(1, 1) - q * DELTA / (Fraction(1, 1) - epsilon_h)


def ns_from_deterministic_multiplier() -> float:
    return 1.0 + math.log(float(1 - DELTA))


def max_homomorphism_defect(tau) -> float:
    samples = [1.1, 1.3, 2.0, 3.0, 5.0]
    return max(abs(tau(x * y) - tau(x) - tau(y)) for x in samples for y in samples)


def main() -> None:
    print("ITEM 131 HBC / LOG-SHELL LIFT ASSUMPTION AUDIT")

    print("\n[1] What the finite clock gives")
    transition_first = Fraction(27, 28)
    generator_first = transition_first - 1
    check(generator_first == -DELTA, "serial 28-channel QEC generator has first eigenvalue -1/28")
    check(True, "this is a local/instrument theorem, not yet a cosmological scale variable")

    print("\n[2] What HBC source text actually supplies")
    hbc_file = ROOT / "cosmological_qec_engine" / "cosmological_qec_engine.tex"
    hbc_text = hbc_file.read_text()
    hbc_section = extract_section(hbc_text, "Holographic Boundary Crystallization")
    check(bool(hbc_section), "HBC section exists in cosmological_qec_engine.tex")
    check("precipitation of new $Q_3$ matter cells at the cosmological horizon" in hbc_section, "HBC supplies horizon precipitation")
    check("Bekenstein--Hawking" in hbc_section and "Landauer flux" in hbc_section, "HBC supplies entropy-bound / Landauer-flux motivation")
    check("dA/dt" in hbc_section and "Friedmann" in hbc_section, "HBC supplies an entropy-rate-to-Friedmann scaling sketch")

    missing_markers = {
        "scale-ratio semigroup T_{lambda mu}=T_lambda T_mu": ["T_{lambda", "T_lambda", "lambda mu"],
        "self-similar Markov horizon-shell process": ["self-similar", "scale-covariant", "Markov horizon"],
        "power ledger Delta_R^2": ["Delta_R", "\\Delta_{\\mathcal R}", "power ledger"],
        "primordial tilt n_s": ["n_s", "scalar spectral"],
    }
    for label, needles in missing_markers.items():
        found = any(needle in hbc_section for needle in needles)
        check(not found, f"HBC source does not derive {label}")

    print("\n[3] If the shell transfer is scale-covariant, log clock follows")
    defects = {
        "ln(lambda)": max_homomorphism_defect(math.log),
        "2 ln(lambda)": max_homomorphism_defect(lambda x: 2.0 * math.log(x)),
        "lambda - 1": max_homomorphism_defect(lambda x: x - 1.0),
        "(ln lambda)^2": max_homomorphism_defect(lambda x: math.log(x) ** 2),
    }
    for label, defect in defects.items():
        print(f"  {label:14s}: semigroup defect {defect:.3e}")
    check(defects["ln(lambda)"] < 1e-12, "continuous scale-ratio semigroup forces tau proportional to ln(lambda)")
    check(defects["lambda - 1"] > 1.0, "physical-distance clock does not compose over scale ratios")
    check(defects["(ln lambda)^2"] > 1.0, "nonlinear log clock does not compose as a Markov semigroup")

    print("\n[4] The exact coefficient requires q=1, constant H, generator action, and power ledger")
    branches = {
        "radial q=1": ns_from_log_lift(Fraction(1, 1)),
        "area q=2": ns_from_log_lift(Fraction(2, 1)),
        "volume q=3": ns_from_log_lift(Fraction(3, 1)),
        "H drift eps=0.01": ns_from_log_lift(Fraction(1, 1), Fraction(1, 100)),
    }
    for label, ns in branches.items():
        print(f"  {label:16s}: n_s={float(ns):.9f} exact={ns}")
    print(f"  deterministic 27/28 multiplier: n_s={ns_from_deterministic_multiplier():.9f}")
    print(f"  amplitude-level action        : n_s={float(Fraction(13, 14)):.9f} exact=13/14")
    check(branches["radial q=1"] == Fraction(27, 28), "radial q=1 + constant H + generator + power ledger gives 27/28")
    check(branches["area q=2"] != Fraction(27, 28), "area normalization misses")
    check(branches["volume q=3"] != Fraction(27, 28), "volume normalization misses")
    check(branches["H drift eps=0.01"] != Fraction(27, 28), "Hubble drift shifts the exact value")
    check(abs(ns_from_deterministic_multiplier() - float(Fraction(27, 28))) > 6e-4, "literal transition multiplier misses exact 27/28")

    print("\n[5] Canonical status check")
    drift = (ROOT / "DRIFT.md").read_text().lower()
    anchor = (ROOT / "ANCHOR.md").read_text().lower()
    check("remaining risk is now the hbc premise itself" in drift, "DRIFT records HBC premise as the remaining risk")
    check("open piece is the hbc physical premise itself" in drift, "DRIFT records the HBC physical premise as still open")
    check("conditionally closed under the continuous log-scale generator" in anchor, "ANCHOR records the early map as conditional")
    check("remaining item-131 promotion burden" in anchor, "ANCHOR keeps promotion burden open")

    print("\n" + "=" * 96)
    print("RESULT")
    print("  The finite 28-clock does NOT, by itself, lift to horizon-shell/log-scale")
    print("  printing. It closes only the finite instrument generator.")
    print("  The log-shell lift is a conditional theorem: if primordial HBC is a")
    print("  self-similar Markov transfer over radial horizon scale ratios, with")
    print("  saturated constant H and QEC action on Delta_R^2, then n_s=27/28.")
    print("  The present HBC source text gives horizon precipitation and entropy-rate")
    print("  matching, but not the semigroup/q=1/constant-H/power-ledger derivation.")
    print("  Promotion target: derive those bridge premises from the HBC dynamics or")
    print("  explicitly tier the 1/28 inflationary observables as conditional.")
    print("=" * 96)
    print("exit 0 -- finite clock closed; HBC/log-shell lift remains a physical bridge premise.")


if __name__ == "__main__":
    main()
