#!/usr/bin/env python3
r"""R13 status audit: service-history measure as record action.

This is a status consolidation, not a new derivation.  It checks the current
shape of the R13 record-action proposal after the 2026-06-19 closures:

    P[gamma] ∝ exp(-A_rec[gamma])
    A_rec = A_traffic - 1/2 Delta S_rec + i Phi_rec + Sum_a lambda_a Q_a.

Earlier wording treated Born, bare alpha0, the traffic clock, and holonomy as
one open bundle.  Current canon separates them:

  * entropy/Crooks/KMS forced legs: closed by record_action_principle.py;
  * Born square: conditionally closed by Gleason + closed-record-pair/Naimark;
  * bare alpha0 = 1/137: closed at foundationally grounded grade by R14;
  * native traffic coefficients: exhausted/free -> conditional/forced/dead;
  * traffic-clock scope: closed for the native cosmological/QEC selector clock;
    the stronger one-clock universalisation is rejected by the two-anchor
    ledger;
  * recovery holonomy/CP: sharpened to K_or plus faithful C3 phase, but still
    open until a physical sign-pointer bridge is derived.

Exit 0 means the current R13 status is internally consistent with that split:
the measure skeleton is conditionally closed, while the full principle remains
open through second-anchor sector matching and the holonomy sign-pointer sector.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Leg:
    name: str
    status: str
    evidence: str
    residual: str


LEGS = [
    Leg(
        "record unit s1 = ln(8x137)",
        "closed-forced",
        "record_action_principle.py; record_alphabet_derivation.py; R12/R14",
        "none at the R13-measure level",
    ),
    Leg(
        "Crooks split / antisymmetric entropy",
        "closed-forced",
        "record_action_principle.py monitored CTMC: traffic cancels in P/Pbar",
        "none; this fixes only ratios, not absolute rates",
    ),
    Leg(
        "KMS equilibrium limit",
        "closed-forced",
        "record_action_principle.py + DRIFT B1 half-Boltzmann scheduler",
        "none for equilibrium form",
    ),
    Leg(
        "Born square / record probabilities",
        "closed-conditional",
        "substrate_born_rule.py + born_closed_record_pair.py",
        "conditional on complex QEC Hilbert substrate / noncontextual record tests",
    ),
    Leg(
        "bare alpha0 service rate",
        "closed-foundational",
        "alpha0_count_rate_theorem.py + alpha0_record_pair_symmetry_theorem.py",
        "dressed alpha and sector use are separate audits, not the bare R13 seed",
    ),
    Leg(
        "native traffic relative coefficients",
        "closed/conditional ledger",
        "traffic audit corrections: K04, CC loop, HBC A_nu, BH 10/27, item-126 T3",
        "conditional premises remain sector-local; no ordinary native fitted coefficient",
    ),
    Leg(
        "traffic clock",
        "native-closed/universal-rejected",
        "A1/N_lock closes the native clock; ew_two_saturation_anchors.py gives an independent top/EW anchor",
        "sector matching remains, but a single universal traffic multiplier is no longer the right target",
    ),
    Leg(
        "holonomy Phi_rec / CP sign",
        "open-sign-pointer",
        "item87 audits plus r12_r13_r15_promotion_gate_attempt.py identify K_or and Phi=2*pi/3 conditionally",
        "derive the sign-representation recovery pointer and absolute orientation sigma from QEC/boot mechanics",
    ),
    Leg(
        "constraint charges Q_a",
        "separate-binary",
        "traffic_multiplier_audit.py separates M15 dust existence from K_rel",
        "binary existence and halo double-counting, not a continuous traffic multiplier",
    ),
]


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def main() -> None:
    print("R13 -- SERVICE-HISTORY MEASURE AS RECORD ACTION: CURRENT STATUS")
    print("=" * 88)
    for leg in LEGS:
        print(f"{leg.status:24s}  {leg.name}")
        print(f"    evidence : {leg.evidence}")
        print(f"    residual : {leg.residual}")

    statuses = {leg.name: leg.status for leg in LEGS}
    print("\n[checks]")
    check(statuses["Born square / record probabilities"].startswith("closed"), "Born is no longer an R13 open leg")
    check(statuses["bare alpha0 service rate"].startswith("closed"), "bare alpha0 is no longer an R13 open leg")
    check(statuses["traffic clock"] == "native-closed/universal-rejected", "native clock is narrowed and one-clock universalisation is rejected")
    check(statuses["holonomy Phi_rec / CP sign"] == "open-sign-pointer", "holonomy/CP remains the genuine R13 sign-pointer frontier")
    check(
        all(leg.status != "free" for leg in LEGS if "native traffic" in leg.name),
        "native record-action coefficients are not ordinary fitted free multipliers",
    )
    check(
        statuses["constraint charges Q_a"] == "separate-binary",
        "conserved-charge existence is not counted as a traffic multiplier",
    )

    print(
        """
[verdict]
  R13 should now be stated as:

      Measure skeleton conditionally closed.
      Full record-action principle now has one rejected overextension and one
      live phase frontier:
        (1) no single universal traffic clock across native and EW/top anchors;
        (2) recovery holonomy / CP still needs the sign-pointer bridge.

  The older bundle "traffic clock + holonomy + Born + alpha0" is stale.  Born
  and bare alpha0 have moved to the reconstructed record calculus.  Native
  traffic coefficients have either converted to sector-native invariants,
  stayed conditional under named sector premises, or been closed negatively.
  What remains is not another count: it is the second-anchor sector-matching
  problem, plus the time-complex holonomy layer.
"""
    )
    print("exit 0 -- R13 status: native clock closed, universal clock rejected; holonomy sign-pointer remains.")


if __name__ == "__main__":
    main()
