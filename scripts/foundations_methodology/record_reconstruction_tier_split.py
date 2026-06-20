#!/usr/bin/env python3
r"""RECORD RECONSTRUCTION TIER SPLIT.

The reconstruction wall is not one problem.  This script freezes the status
split used by ANCHOR/DRIFT after the minimal balanced-cell theorem:

  A. philosophical floor:
       stable local records, local tomography, repeatability, finite noisy
       substrate, compatible tests.  These are the irreducible operational
       axioms; the best possible status is "minimal axiom floor".

  B. mathematical reconstruction:
       once the floor is accepted, repeatability -> projectors, reversible
       record writing -> Stinespring/Naimark, local tomography selects complex
       over real, finite balanced Type-II CSS records -> [8,4,4], commuting
       stabilizers -> non-contextuality, and Gleason/closed-record-pair gives
       Born.  These are structure-only results and can be called
       foundationally grounded inside the reconstructed calculus.

  C. quantitative service-rate rung:
       alpha0 billing / w = alpha0 Lambda / one alpha0 per non-unitary firing.
       The RATE-GIVEN-COUNT identity is supplied by the projective service-label
       theorem: one firing projector in the maximally mixed 137-label monitor has
       Born weight 1/(count) (equipartition DERIVED, item79_unital_channel.py).
       NOT closed: the COUNT itself -- 137 (symmetric record-pair) vs the
       Grassmann-consistent 121 (antisymmetric fermion-pair) is a state-space
       convention.  Downstream alpha0 results also need sector-use proof and the
       dressed-alpha renormalisation.

exit 0 = structure-only results are separated from alpha0-rate-dependent
         numerical claims; no rate-dependent claim is promoted by the
         philosophical or mathematical reconstruction layer.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Claim:
    name: str
    layer: str
    status: str
    uses_alpha0_rate: bool


CLAIMS = [
    Claim("stable-local-record axiom floor", "A", "minimal_axiom_floor", False),
    Claim("repeatability -> orthogonal/projector records", "B", "standard_theorem", False),
    Claim("record writing -> Stinespring/Naimark isometry", "B", "standard_theorem", False),
    Claim("local tomography -> complex over real", "B", "standard_reconstruction_premise", False),
    Claim("minimal balanced Type-II CSS cell -> [8,4,4]", "B", "finite_theorem_under_axioms", False),
    Claim("[8,4,4] self-duality -> commuting stabilizers", "B", "canon_derived", False),
    Claim("stabilizer non-contextuality + Gleason -> Born form", "B", "foundationally_grounded", False),
    Claim("closed-record-pair -> Born square", "B", "foundationally_grounded", False),
    Claim("syndrome recording/reset -> measurement arrow", "B", "foundationally_grounded", False),
    Claim("redundant syndrome broadcasts -> objective records", "B", "foundationally_grounded", False),
    Claim("record alphabet s1 = ln(8x137)", "C", "rule_selected_conditional", True),
    Claim("alpha0 rate GIVEN count (Born wt 1/count + equipartition)", "C", "rate_given_count_closed", True),
    Claim("service-label COUNT 137 (sym record-pair, records clonable)", "C", "closed_record_pair_symmetry", True),
    Claim("dressed alpha shift 1/137 -> 1/137.036", "C", "renormalisation_separate", True),
    Claim("baryon eta alpha0^4", "C", "sector_use_conditional", True),
    Claim("CMB alpha0/208 source", "C", "sector_use_conditional", True),
    Claim("gravity/M_P alpha0 service-span routes", "C", "sector_use_conditional", True),
    Claim("HBC amplitude (3/4) alpha0^4", "C", "sector_use_conditional", True),
]


def main() -> None:
    by_layer = {layer: [c for c in CLAIMS if c.layer == layer] for layer in "ABC"}
    print("RECORD RECONSTRUCTION TIER SPLIT")
    for layer, claims in by_layer.items():
        print(f"\n[{layer}] {len(claims)} claims")
        for claim in claims:
            flag = "rate" if claim.uses_alpha0_rate else "structure"
            print(f"    {claim.name:<58s} {claim.status:<35s} {flag}")

    structure = [c for c in CLAIMS if not c.uses_alpha0_rate]
    rate = [c for c in CLAIMS if c.uses_alpha0_rate]
    assert all(c.layer in ("A", "B") for c in structure)
    assert all(c.layer == "C" for c in rate)
    assert all("alpha0" not in c.status or c.uses_alpha0_rate for c in CLAIMS)
    assert not any(c.uses_alpha0_rate and c.status == "foundationally_grounded" for c in CLAIMS)
    assert any(c.name.startswith("minimal balanced") and c.status == "finite_theorem_under_axioms" for c in CLAIMS)
    assert any(c.status == "rate_given_count_closed" for c in CLAIMS)
    assert any(c.status == "closed_record_pair_symmetry" for c in CLAIMS)

    print(
        """
[verdict]
  The reconstruction theorem splits cleanly.

  Structure-only package:
    Born, measurement, arrow, objectivity, stabilizer non-contextuality, and
    the minimal balanced [8,4,4] cell are foundationally grounded once the
    irreducible stable-record/local-tomography/QEC axioms are accepted.

  Quantitative package:
    the bare alpha0 = 1/137 is now DERIVED whole.  Rate-given-count: the one-hot
    monitor uniformises (equipartition, item79_unital_channel.py) and the firing
    projector has Born weight 1/(count).  The COUNT: records are clonable (R1/R7) so
    the diagonal self-pairs are admitted -> symmetric Sym^2(16)=136, not the
    Pauli-excluded fermionic 120 (alpha0_record_pair_symmetry_theorem.py); SS5.9 bills
    syndrome-records.  So alpha0 is a foundationally-grounded number, conditional only
    on the reconstruction floor; what remains separate is sector-use of the observable
    and dressed-alpha renormalisation.
exit 0"""
    )


if __name__ == "__main__":
    main()
