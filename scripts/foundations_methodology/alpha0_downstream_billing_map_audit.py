#!/usr/bin/env python3
r"""R14 downstream audit: what remains after bare alpha_0 = 1/137 closes.

The R14 closure removes the common bare-rate question:

    alpha_0 = 1/137

is now the Born weight of one firing projector in the monitored 137-label
service register, with the 137 count itself selected by clonable record-pairs.
This script classifies the remaining alpha_0-dependent downstream
conditionalities into two and only two live buckets:

  1. sector billing maps:
       does a given sector actually bill THIS monitored service observable, and
       with which event multiplicity?
  2. dressed-alpha:
       how does the bare 1/137 service probability become the physical
       1/137.035999... QED coupling?

It is a bookkeeping theorem, not a new numerical fit.  Exit 0 means every named
sector is assigned to one of the buckets and the residuals are not conflated.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


ALPHA0_INV = 137
ALPHA_PHYS_INV = 137.035999084


@dataclass(frozen=True)
class Sector:
    name: str
    use: str
    status: str
    evidence: str
    residual: str


SECTORS = [
    Sector(
        "R14 bare service observable",
        "one firing projector in the one-hot 137-label service monitor",
        "closed",
        "alpha0_count_rate_theorem.py + alpha0_record_pair_symmetry_theorem.py",
        "none at bare-alpha level; only reconstruction-floor assumptions",
    ),
    Sector(
        "CMB / sterile zero-mode source alpha0/208",
        "one S3-singlet sterile release port, one non-unitary firing, uniform Q-address 1/208",
        "billing_map_closed_conditional",
        "sterile_release_billing.py",
        "absolute density still uses m_nuR = alpha0^2 Lambda and Boltzmann/halo non-double-counting",
    ),
    Sector(
        "mass-sector bit weight m_d - m_u = alpha Lambda",
        "one active I3 syndrome/gauge-bridge erasure per confined tick",
        "billing_map_closed_conditional",
        "ANCHOR 5.9 + item79_loop_closure.py",
        "physical comparison depends on confinement-scale running and dressed-alpha convention",
    ),
    Sector(
        "gravity / M_P alpha dressing",
        "one P->Q monitored jump in the virtual graviton loop; coherent return unbilled",
        "billing_map_closed_conditional",
        "item119_jump_channel.py",
        "Planck hierarchy still needs horizon/selector relations and the bare-vs-dressed convention",
    ),
    Sector(
        "CC generation vertex",
        "one monitored register touch on the chi/W generation vertex",
        "billing_map_closed_conditional",
        "cc_monitored_billing_operator_algebra.md + cc_billing_exponentiation_selection.py",
        "next-order C residual is observationally sub-floor / operator-class exhausted, not a billing-map gap",
    ),
    Sector(
        "HBC scalar amplitude",
        "weight-4 topology commits billed by alpha0^4 inside the local scalar-printer ledger",
        "billing_map_conditional_on_scalar_printer",
        "item131_hbc_stop_rule_proof.py",
        "scalar-printer identification and no-new-nonlocal-scalar-source premise",
    ),
    Sector(
        "baryogenesis eta",
        "portal-licensed colourless weight-4 logical commits, alpha0^4 per exhaust photon",
        "billing_map_closed_conditional",
        "scheduler_alpha_composition.py + record_alphabet_derivation.py",
        "B-sign/CP and exhaust-photon bookkeeping, not the alpha0 power",
    ),
    Sector(
        "dressed alpha",
        "escape / relaxation dressing of the monitored bridge into the physical QED coupling",
        "open",
        "dressed_alpha_nonunital.py",
        "derive Gamma_escape / bridge_gap from the gauge-web density of states; current script gives sign/form only",
    ),
]


def main() -> None:
    print("R14 DOWNSTREAM ALPHA0 BILLING-MAP AUDIT")
    print("=" * 88)
    print(f"bare alpha0^-1     = {ALPHA0_INV}")
    print(f"physical alpha^-1  = {ALPHA_PHYS_INV:.9f}")
    print(f"dressed shift      = {ALPHA_PHYS_INV - ALPHA0_INV:.9f}")
    print()

    for sector in SECTORS:
        print(f"- {sector.name}")
        print(f"    use      : {sector.use}")
        print(f"    status   : {sector.status}")
        print(f"    evidence : {sector.evidence}")
        print(f"    residual : {sector.residual}")

    statuses = {s.status for s in SECTORS}
    assert "closed" in statuses
    assert "open" in statuses
    assert sum(1 for s in SECTORS if s.name == "dressed alpha") == 1
    assert all("free bare" not in s.residual.lower() for s in SECTORS)

    # Sector map sanity checks: the named residuals should not be mistaken for
    # a reopened bare alpha0 rate.
    closed_or_conditional = [s for s in SECTORS if s.name != "dressed alpha"]
    assert all(s.status in {
        "closed",
        "billing_map_closed_conditional",
        "billing_map_conditional_on_scalar_printer",
    } for s in closed_or_conditional)

    # The dressed shift is small enough to be convention-level in many current
    # predictions, but it is not zero and not derived by the bare count theorem.
    shift = ALPHA_PHYS_INV - ALPHA0_INV
    assert Fraction(ALPHA0_INV, 1) == 137
    assert 0.03 < shift < 0.04

    print("\n" + "=" * 88)
    print("VERDICT")
    print("  R14 closed the shared bare-rate problem.  The remaining alpha0")
    print("  conditionality is not one fuzzy foundation issue: it is a checklist.")
    print("  Most named downstream sectors already have a sector billing map at")
    print("  conditional/model grade.  Their live caveats are sector premises")
    print("  (halo/CMB, scalar-printer identity, horizon/selector input, CP sign),")
    print("  not a free bare alpha0 rate.  The single genuinely separate alpha")
    print("  problem is the dressed shift 1/137 -> 1/137.035999..., whose sign/form")
    print("  is sharpened but whose magnitude still requires the photon escape /")
    print("  bridge-relaxation calculation.")
    print("exit 0")


if __name__ == "__main__":
    main()
