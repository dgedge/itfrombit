#!/usr/bin/env python3
r"""R12/R13 theorem-boundary audit: record alphabet and record action.

Purpose
-------
Separate three things that are easy to conflate:

  1. derived subclauses of the R12 per-event record alphabet;
  2. the maximum-caliber/exponential form of a history measure once constraints
     are named;
  3. a genuine substrate theorem that those constraints are complete and that
     this action is the universal dynamical principle.

The current canon supports (1) strongly and (2) as standard inference
mathematics, but not (3).  The 2026-06-21 promotion-gate audit closes the
one-record/factorization gate at finite service-instrument grade, while leaving
hidden pointer algebras outside the admitted event instrument.  Therefore this
script exits 0 only if the audit keeps R12/R13 at the correct tier: a
conditional skeleton with derived pieces and one closed finite-instrument gate,
not a universal least-action theorem.
"""

from __future__ import annotations

from dataclasses import dataclass
import math


ADDRESS = 8
CHANNEL = 137
RECORD_ALPHABET = ADDRESS * CHANNEL
S1 = math.log(RECORD_ALPHABET)


@dataclass(frozen=True)
class Clause:
    name: str
    status: str
    evidence: str
    residual: str


R12_CLAUSES = [
    Clause(
        "address alphabet",
        "derived",
        "record_content_from_syndrome.py: QND syndrome increments are the 8 distinct vertex triples",
        "none for the alphabet size 8",
    ),
    Clause(
        "address equidistribution",
        "conditional-theorem",
        "r12_address_equidistribution_theorem.py: RM(1,3) translation symmetry -> convolution -> uniform",
        "requires no service dynamics that spontaneously breaks the cell translation symmetry",
    ),
    Clause(
        "service channel 137",
        "derived-on-reconstruction-floor",
        "alpha0_count_rate_theorem.py + alpha0_record_pair_symmetry_theorem.py",
        "inherits the record-pair/reconstruction floor; dressed-alpha is separate",
    ),
    Clause(
        "address x channel factorization",
        "closed-finite-instrument",
        "r12_r13_r15_promotion_gate_attempt.py: C^8 tensor C^137 event instrument is bijective",
        "hidden timing/environment labels would be new pointer algebras, not part of this event instrument",
    ),
    Clause(
        "one committed record per event",
        "closed-finite-instrument",
        "r12_r13_r15_promotion_gate_attempt.py: joint one-hot projectors resolve the identity once",
        "applies to minimal irreversible commits of the admitted pointer inventory",
    ),
]


R13_CLAUSES = [
    Clause(
        "maximum-caliber exponential form",
        "mathematical-theorem-given-constraints",
        "maximum entropy over histories with additive constraints gives P[gamma] proportional to exp(-A_rec)",
        "does not derive which constraints the substrate must include",
    ),
    Clause(
        "record unit s1",
        "closed-finite-instrument",
        "R12 record alphabet gives s1=ln(8x137)",
        "inherits the finite-instrument event convention; external pointer labels would define a different instrument",
    ),
    Clause(
        "Crooks/KMS entropy split",
        "forced",
        "record_action_principle.py: traffic cancels in forward/reverse ratios; KMS is the equilibrium limit",
        "fixes ratios and equilibrium form, not the absolute traffic clock",
    ),
    Clause(
        "Born square",
        "closed-conditional",
        "born_closed_record_pair.py and substrate_born_rule.py",
        "conditional on the complex QEC Hilbert substrate and noncontextual record tests",
    ),
    Clause(
        "bare alpha0 service rate",
        "closed-on-reconstruction-floor",
        "R14 record-pair/equipartition theorem",
        "sector billing maps and dressed alpha remain separate",
    ),
    Clause(
        "traffic clock universalisation",
        "rejected-as-universal",
        "r12_r13_r15_promotion_gate_attempt.py + ew_two_saturation_anchors.py",
        "native selector clock survives, but EW/top and QCD give two irreducible saturated anchors",
    ),
    Clause(
        "recovery holonomy Phi_rec",
        "open-sign-pointer",
        "item87/R13 audits identify K_or and the faithful C3 phase as the sharp target",
        "physical sign-representation pointer and absolute orientation are not derived from QEC/boot mechanics",
    ),
]


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def normalized_exponential(beta: float, observables: list[float]) -> list[float]:
    weights = [math.exp(-beta * obs) for obs in observables]
    total = sum(weights)
    return [weight / total for weight in weights]


def entropy(probs: list[float]) -> float:
    return -sum(p * math.log(p) for p in probs if p > 0.0)


def expectation(probs: list[float], obs: list[float]) -> float:
    return sum(p * x for p, x in zip(probs, obs))


def print_clauses(title: str, clauses: list[Clause]) -> None:
    print(f"\n{title}")
    for clause in clauses:
        print(f"  {clause.name:34s} {clause.status:36s}")
        print(f"    evidence : {clause.evidence}")
        print(f"    residual : {clause.residual}")


def main() -> None:
    print("R12/R13 RECORD-ACTION THEOREM-BOUNDARY AUDIT")
    print("=" * 96)

    print("\n[1] R12 record alphabet arithmetic")
    print(f"  address alphabet       = {ADDRESS}")
    print(f"  service channel        = {CHANNEL}")
    print(f"  record alphabet M      = {RECORD_ALPHABET}")
    print(f"  s1 = ln(M)             = {S1:.6f}")
    print(f"  exp(-s1)               = {math.exp(-S1):.9f}")
    check(RECORD_ALPHABET == 1096, "M = 8 x 137 = 1096")
    check(abs(S1 - 6.999422) < 1.0e-6, "s1 is the canonical ln(8x137) record unit")
    check(abs(math.exp(-S1) - 1.0 / RECORD_ALPHABET) < 1.0e-15, "record weight is 1/M if the record-action reading is adopted")

    print_clauses("[2] R12 clause ledger", R12_CLAUSES)
    r12_closed_grade = all(
        clause.status in {"derived", "forced", "conditional-theorem", "derived-on-reconstruction-floor", "closed-finite-instrument"}
        for clause in R12_CLAUSES
    )
    check(r12_closed_grade, "R12 finite-instrument record clauses are all derived or conditionally derived")
    check(
        any(clause.name == "address equidistribution" and clause.status == "conditional-theorem" for clause in R12_CLAUSES),
        "address uniformity is stronger than a fit but still carries the no-symmetry-breaking premise",
    )

    print("\n[3] Maximum-caliber form: theorem given constraints, not substrate dynamics")
    observables = [0.0, 1.0, 2.0, 4.0]
    beta = 0.73
    probs = normalized_exponential(beta, observables)
    print("  toy histories with additive observable O =", observables)
    print("  beta =", beta)
    print("  P =", [round(p, 6) for p in probs])
    print(f"  entropy H(P)           = {entropy(probs):.6f}")
    print(f"  expected O             = {expectation(probs, observables):.6f}")
    for i in range(len(observables) - 1):
        lhs = math.log(probs[i] / probs[i + 1])
        rhs = beta * (observables[i + 1] - observables[i])
        check(abs(lhs - rhs) < 1.0e-12, "exponential-family log ratios match beta times action difference")
    print("  This proves only the inference form: once the additive constraints are named,")
    print("  the least-biased history measure is exponential.  It does not prove that the")
    print("  named R13 action is the complete physical action.")

    print_clauses("[4] R13 clause ledger", R13_CLAUSES)
    r13_locked = all(
        clause.status in {"derived", "forced", "closed-conditional", "closed-finite-instrument"}
        for clause in R13_CLAUSES
    )
    check(not r13_locked, "R13 is not locked: universal clock is rejected and recovery holonomy needs a sign pointer")
    check(
        any(clause.name == "maximum-caliber exponential form" for clause in R13_CLAUSES),
        "the exponential measure is an inference theorem only after constraints are supplied",
    )
    check(
        any(clause.status.startswith("open") for clause in R13_CLAUSES),
        "at least one live holonomy clause remains in the least-action record principle",
    )

    print("\n[5] Tier verdict")
    print("  R12 status: derived core plus finite-instrument event-binding/factorization.")
    print("  R13 status: maximum-caliber skeleton with forced thermodynamic legs, but")
    print("  not a universal substrate least-action theorem.")
    print("  Promotion-gate status:")
    print("    A. one-record-per-event and address x channel factorization close at")
    print("       finite service-instrument grade;")
    print("    B. one universal traffic clock is rejected by the two-anchor ledger;")
    print("    C. recovery holonomy reduces to a missing sign-representation pointer.")
    print("  Therefore the record action is a sharpened organizing principle with")
    print("  a closed event-instrument gate, not a locked universal dynamical law.")
    print("exit 0 -- finite-instrument gate promoted; R13 universal action still not theorem.")


if __name__ == "__main__":
    main()
