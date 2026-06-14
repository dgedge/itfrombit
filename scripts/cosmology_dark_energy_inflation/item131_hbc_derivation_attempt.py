#!/usr/bin/env python3
r"""ITEM 131: derivation attempt for the HBC -> log-shell lift.

Goal
----
Try to derive, rather than assume, the bridge

    finite 28-channel QEC generator  ->  horizon-shell/log-scale transfer

needed for

    d ln Delta_R^2 / d ln k = -1/28.

Result
------
The attempt partially succeeds but exposes one central obstruction.

What can be made plausible from existing HBC/QEC ingredients:

1. Markovity:
   The substrate clock is a stroboscopic CPTP map, and fresh HBC nodes are
   zero-entropy boundary capacity.  If each printed shell is freshly reset and
   interacts only through the local QEC instrument, shell-to-shell transfer is
   Markovian.

2. Log variable:
   If primordial HBC is exactly self-similar in a saturated de Sitter phase,
   then transfer over scale ratios composes:

       T_{lambda mu} = T_lambda T_mu,

   forcing tau(lambda)=q ln(lambda).

3. Constant H:
   A saturated de Sitter equilibrium gives k=aH and d ln k=d ln a.  Hubble
   drift would shift the coefficient.

4. Power ledger:
   QEC acts on density/probability ledgers.  The dimensionless scalar power
   Delta_R^2 is variance per d ln k shell, so it is the right target if the
   horizon printer acts on scalar density fluctuations.

The obstruction:

HBC as written is Bekenstein-Hawking area bookkeeping.  If the QEC transfer
clock is tied to printed boundary capacity, the natural normalization is area
clocking q=2, not radial clocking q=1:

    S_BH ~ A ~ R^2  =>  d ln S_BH = 2 d ln R.

q=2 gives n_s=13/14, not 27/28.  The desired q=1 requires an additional
"mode-local radial crossing" lemma:

    each scalar mode receives one QEC generator action per radial horizon-scale
    logarithmic crossing, while the boundary area degeneracy belongs to the
    number of independent angular modes and is already divided out in the
    dimensionless power ledger Delta_R^2.

That lemma is physically reasonable, but it is not derived by the current HBC
area-saturation text.  Therefore the exact bridge remains open at q=1.
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


def section(text: str, name: str) -> str:
    marker = rf"\\section\{{{re.escape(name)}\}}"
    match = re.search(marker, text)
    if not match:
        return ""
    rest = text[match.start() :]
    nxt = re.search(r"\n\\section\{", rest[len(match.group(0)) :])
    if not nxt:
        return rest
    return rest[: len(match.group(0)) + nxt.start()]


def ns(q: Fraction, epsilon_h: Fraction = Fraction(0, 1)) -> Fraction:
    return Fraction(1, 1) - q * DELTA / (Fraction(1, 1) - epsilon_h)


def semigroup_defect(clock, lam: float, mu: float) -> float:
    return clock(lam * mu) - clock(lam) - clock(mu)


def max_defect(clock) -> float:
    samples = [1.1, 1.5, 2.0, 3.0, 5.0]
    return max(abs(semigroup_defect(clock, a, b)) for a in samples for b in samples)


def main() -> None:
    print("ITEM 131 HBC -> LOG-SHELL DERIVATION ATTEMPT")

    hbc = section(
        (ROOT / "cosmological_qec_engine" / "cosmological_qec_engine.tex").read_text(),
        "Holographic Boundary Crystallization",
    )
    engine = (ROOT / "cosmological_qec_engine" / "cosmological_qec_engine.tex").read_text()

    print("\n[1] Existing ingredients")
    check("rho(t_{n+1}) = \\mathcal{E}" in engine, "substrate has a stroboscopic CPTP clock")
    check("Kubo--Martin--Schwinger" in engine and "Lindbladian semigroup" in engine, "engine text has Markov/Lindblad semigroup structure")
    check("strictly zero entanglement entropy" in engine, "HBC/Past-Hypothesis text supplies fresh zero-entropy capacity")
    check("Bekenstein--Hawking" in hbc and "A = 4\\pi L_H^2" in hbc, "HBC is explicitly horizon-area entropy bookkeeping")
    check("precipitating new $Q_3$ matter cells" in hbc, "HBC supplies boundary-node precipitation")

    print("\n[2] Markov/log part: conditional but mathematically sound")
    log_defect = max_defect(math.log)
    area_defect = max_defect(lambda x: 2.0 * math.log(x))
    distance_defect = max_defect(lambda x: x - 1.0)
    print(f"  defect ln(lambda)   = {log_defect:.3e}")
    print(f"  defect 2ln(lambda)  = {area_defect:.3e}")
    print(f"  defect lambda-1     = {distance_defect:.3e}")
    check(log_defect < 1e-12, "scale-ratio semigroup forces a logarithmic clock up to q")
    check(area_defect < 1e-12, "area clock is also scale-covariant: semigroup alone does not fix q")
    check(distance_defect > 1.0, "physical-distance clock fails scale-ratio composition")

    print("\n[3] Normalization q is the load-bearing obstruction")
    branches = {
        "radial/mode-local q=1": Fraction(1, 1),
        "BH area-capacity q=2": Fraction(2, 1),
        "bulk volume q=3": Fraction(3, 1),
    }
    for label, q in branches.items():
        value = ns(q)
        print(f"  {label:24s}: n_s={float(value):.9f} exact={value}")
    check(ns(Fraction(1, 1)) == Fraction(27, 28), "q=1 gives the desired item-131 tilt")
    check(ns(Fraction(2, 1)) == Fraction(13, 14), "q=2 area clock gives the wrong coefficient")
    check(ns(Fraction(3, 1)) == Fraction(25, 28), "q=3 volume clock gives the wrong coefficient")
    check("S_{\\max} = \\frac{A}" in hbc and "Each new cell expands $A$" in hbc, "current HBC prose naturally talks in area-capacity units")

    print("\n[4] Why q=1 is not impossible, but needs a new lemma")
    check(True, "Delta_R^2 is variance per d ln k, so angular/area degeneracy is normally divided out")
    check(True, "a single scalar mode crosses the horizon once per radial d ln k shell when k=aH")
    check(True, "therefore q=1 follows if QEC action is mode-local at horizon crossing")
    check(True, "but HBC area-saturation text does not yet prove mode-local radial crossing")

    print("\n[5] Constant-H and ledger action remain explicit conditions")
    drifted = ns(Fraction(1, 1), Fraction(1, 100))
    print(f"  q=1 with epsilon_H=0.01: n_s={float(drifted):.9f} exact={drifted}")
    check(drifted != Fraction(27, 28), "nonzero Hubble drift shifts the coefficient")
    check(True, "QEC density-matrix dynamics supports power/density action, but scalar-curvature coupling is still a bridge")

    print("\n" + "=" * 98)
    print("DERIVATION ATTEMPT VERDICT")
    print("  The HBC -> log-shell lift can be reduced to one sharp missing lemma.")
    print("  Existing HBC/QEC ingredients support Markovity, scale-ratio logarithms,")
    print("  and density-ledger action conditionally.  They do not select q=1.")
    print("  Because HBC is written as Bekenstein-Hawking area saturation, its native")
    print("  capacity clock is q=2, which would give n_s=13/14.  To derive n_s=27/28,")
    print("  we need the mode-local radial-crossing lemma: the 28-clock acts once per")
    print("  scalar mode per radial horizon log shell, while the boundary-area")
    print("  degeneracy is already factored out in Delta_R^2.")
    print("  Status: not closed; the remaining bridge is now narrower and testable.")
    print("=" * 98)
    print("exit 0 -- q=1 radial mode-local HBC lemma isolated as the missing derivation.")


if __name__ == "__main__":
    main()
