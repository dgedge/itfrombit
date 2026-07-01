#!/usr/bin/env python3
r"""Correlated-null falsification ledger.

The framework makes several absence predictions that should be read together,
not as isolated "nothing happens" claims:

* frozen constants: no secular drift of G or the code-level fine-structure
  constant alpha_0;
* strong-CP closure: neutron EDM only at the tiny weak-leakage floor;
* boundary-printing inflation: no observable primordial tensor background;
* no high-scale de Sitter graviton vacuum: any B-mode detection at ordinary
  high-scale-inflation levels breaks the one-bit scalar-printer branch.

The script keeps the numbers deliberately simple.  It is a public-facing
attack-surface ledger, not a re-analysis of the experiments.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class NullPrediction:
    name: str
    mechanism: str
    prediction: float
    comparison: float
    unit: str
    positive_signal: str
    damage_scope: str

    @property
    def margin(self) -> float:
        if self.prediction == 0.0:
            return float("inf")
        return self.comparison / abs(self.prediction)


def fmt(x: float) -> str:
    if x == 0.0:
        return "0"
    if x == float("inf"):
        return "inf"
    return f"{x:.2e}"


def main() -> None:
    # Conservative, rounded external comparison scales:
    # * |Gdot/G|: current LLR-scale nulls are O(1e-14/yr) in newer analyses,
    #   depending on ephemeris, solar-mass-loss modelling, and data selection.
    # * |alphadot/alpha|: Rosenband-class clock limit O(1e-17/yr); newer clocks
    #   improve frequency ratios but the long-baseline alpha-dot null remains the
    #   clean public benchmark.
    # * neutron EDM: current world limit |d_n| < 1.8e-26 e cm; near-future searches
    #   target 1e-28..1e-30 e cm.
    # * tensors: current CMB limit O(1e-2); the framework-killing practical
    #   threshold is a robust primordial r >= 1e-3.
    rows = [
        NullPrediction(
            name="G drift",
            mechanism="frozen gravitational ledger after crystallisation lock-in",
            prediction=0.0,
            comparison=1.0e-14,
            unit="yr^-1",
            positive_signal="reproducible secular Gdot/G",
            damage_scope="global constant-ledger failure",
        ),
        NullPrediction(
            name="alpha drift",
            mechanism="frozen code alpha_0; only ordinary energy-scale running",
            prediction=0.0,
            comparison=2.3e-17,
            unit="yr^-1",
            positive_signal="clock/quasar drift requiring a changing code coupling",
            damage_scope="global code-ledger failure",
        ),
        NullPrediction(
            name="neutron EDM",
            mechanism="bare theta_bar=0; residual only weak CKM leakage",
            prediction=1.0e-31,
            comparison=1.8e-26,
            unit="e cm",
            positive_signal="d_n >> 1e-30 e cm without an added CP source",
            damage_scope="strong-CP closure failure",
        ),
        NullPrediction(
            name="primordial tensors",
            mechanism="scalar boundary-printer; no squeezed graviton vacuum",
            prediction=2.1e-9,
            comparison=1.0e-3,
            unit="r",
            positive_signal="primordial B-modes at r >= 1e-3",
            damage_scope="boundary-printing inflation failure",
        ),
        NullPrediction(
            name="high-scale tensor background",
            mechanism="one-bit event unit selects no-squeezing branch",
            prediction=2.1e-9,
            comparison=3.6e-2,
            unit="r",
            positive_signal="standard high-scale inflationary tensor background",
            damage_scope="one-bit printer / no-squeezing premise failure",
        ),
    ]

    print("CORRELATED NULL PREDICTIONS")
    print("=" * 110)
    print("[1] ledger")
    print(f"{'null':<28} {'prediction':>12} {'comparison':>12} {'unit':<8} {'margin':>12} mechanism")
    for row in rows:
        margin = "infinite" if row.margin == float("inf") else f"{row.margin:.2e}x"
        print(
            f"{row.name:<28} {fmt(row.prediction):>12} {fmt(row.comparison):>12} "
            f"{row.unit:<8} {margin:>12} {row.mechanism}"
        )

    print("\n[2] if a null turns positive")
    for row in rows:
        print(f"    {row.name:<28} -> {row.damage_scope}: {row.positive_signal}")

    print("\n[3] correlation structure")
    print("    A. Frozen-ledger nulls:")
    print("       G drift and alpha drift are both forbidden because they are crystallisation/code readouts,")
    print("       not live cosmological service rates. A positive drift in either is global damage.")
    print("    B. Phase-null:")
    print("       neutron EDM is suppressed because the bare strong-CP phase is geometrically absent;")
    print("       a large EDM is not a small correction, it means the phase-null leaked or never held.")
    print("    C. Tensor-null:")
    print("       no observable primordial B-modes follows from scalar boundary printing plus no")
    print("       high-scale squeezed graviton vacuum. A robust r >= 1e-3 is the sharpest kill switch.")
    print("    D. What makes the signature recognizable:")
    print("       LCDM+generic-BSM expects at least some of these nulls to be independent accidents;")
    print("       here they are linked by frozen ledgers and absence of live high-scale UV oscillators.")

    # Guardrails: keep the public claims in their intended regimes.
    assert rows[2].comparison / rows[2].prediction > 1.0e5
    assert rows[3].comparison / rows[3].prediction > 1.0e5
    assert rows[4].comparison / rows[4].prediction > 1.0e7
    assert rows[0].prediction == 0.0 and rows[1].prediction == 0.0
    print("\nexit 0 -- correlated-null ledger is internally consistent.")


if __name__ == "__main__":
    main()
