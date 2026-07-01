#!/usr/bin/env python3
"""Minimal local-record architecture audit for Target A.

This script tests the precise statement needed after
target_a_classification.py:

  * A 7-bit compressed encoding of the SM active content exists if the
    lepton/quark predicate is merged into the colour plane.
  * That encoding is not comparable to the framework's local-record
    architecture, because it has no primitive one-bit LQ observable and its
    sterile exclusion must inspect the colour/lepton plane.
  * Once the comparable-locality axioms are imposed -- independent generation
    pair, colour plane, LQ bit, I3, chirality, and repeatable W readout -- the
    bit-count lower bound is 8, and the existing 8-bit classification isolates
    the SM active rule up to relabelling.

Exit 0 means these finite checks reproduce the printed theorem.
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations, product

import target_a_classification as target_a


def primitive_record_lower_bound() -> dict[str, int]:
    """Return the primitive role widths for comparable local architectures."""

    return {
        "generation_pair": 2,  # punctured 2-bit rule gives three generations
        "colour_plane": 2,     # 00 singlet + three nonzero colour labels
        "lq_bit": 1,           # independent lepton/quark observable
        "weak_i3": 1,
        "chirality": 1,
        "weak_readout_w": 1,   # repeatable W record, separate from chi
    }


def count_disjoint_role_assignments(nbits: int) -> int:
    """Count labelled local-record role assignments with unordered pairs.

    The generation pair and colour plane are unordered pairs; LQ, I3, chi, W
    are labelled one-bit records.  The count is intentionally simple: the
    lower-bound claim is just disjoint additivity of primitive record storage.
    """

    if nbits < sum(primitive_record_lower_bound().values()):
        return 0

    bits = range(nbits)
    total = 0
    for gen in combinations(bits, 2):
        rest1 = [b for b in bits if b not in gen]
        for colour in combinations(rest1, 2):
            rest2 = [b for b in rest1 if b not in colour]
            # Ordered choices for LQ, I3, chi, W.
            for lq in rest2:
                rest3 = [b for b in rest2 if b != lq]
                for i3 in rest3:
                    rest4 = [b for b in rest3 if b != i3]
                    for chi in rest4:
                        rest5 = [b for b in rest4 if b != chi]
                        for _w in rest5:
                            total += 1
    return total


# Seven-bit compressed model:
# (G0,G1,CL0,CL1,I3,chi,W), where CL=00 is lepton and CL!=00 are colours.
CW7 = list(product([0, 1], repeat=7))


def r1_7(c: tuple[int, ...]) -> bool:
    return not (c[0] == 1 and c[1] == 1)


def r2_7(c: tuple[int, ...]) -> bool:
    return c[6] == c[5]


def cl_pair(c: tuple[int, ...]) -> tuple[int, int]:
    return (c[2], c[3])


def r4_7(c: tuple[int, ...]) -> bool:
    return not (cl_pair(c) == (0, 0) and c[4] == 0 and c[5] == 1)


def zf_7(c: tuple[int, ...]) -> int:
    return 1 if c[4] == 0 else -1


def sum_zc_7(c: tuple[int, ...]) -> int:
    return -3 if cl_pair(c) == (0, 0) else -1


def q_7(c: tuple[int, ...]) -> F:
    return F(1, 2) * zf_7(c) + F(1, 3) * sum_zc_7(c) + F(1, 2)


def t3_7(c: tuple[int, ...]) -> F:
    return F(1, 2) * zf_7(c) if c[5] == 0 else F(0)


def y_eff_7(c: tuple[int, ...]) -> F:
    yy = q_7(c) - t3_7(c)
    return yy if c[5] == 0 else -yy


def triplet_7(c: tuple[int, ...]) -> bool:
    return cl_pair(c) != (0, 0)


def doublet_7(c: tuple[int, ...]) -> bool:
    return c[5] == 0


def content_7() -> list[tuple[int, ...]]:
    return [c for c in CW7 if r1_7(c) and r2_7(c) and r4_7(c)]


def anomalies_7(states: list[tuple[int, ...]]) -> tuple[F | int, ...]:
    a1 = sum(y_eff_7(c) for c in states if triplet_7(c))
    a2 = F(1, 2) * sum(y_eff_7(c) for c in states if doublet_7(c))
    a3 = sum(y_eff_7(c) ** 3 for c in states)
    a4 = sum(y_eff_7(c) for c in states)
    a5 = sum(1 for c in states if triplet_7(c) and c[5] == 0) - sum(
        1 for c in states if triplet_7(c) and c[5] == 1
    )
    a6 = (sum(1 for c in states if doublet_7(c)) // 2) % 2
    return (a1, a2, a3, a4, a5, a6)


def vectorlike_7(states: list[tuple[int, ...]]) -> bool:
    s = set(states)
    for c in states:
        partner = list(c)
        partner[5] = 1 - partner[5]
        partner[6] = partner[5]
        if tuple(partner) not in s:
            return False
    return True


def su2_complete_7(states: list[tuple[int, ...]]) -> bool:
    by_base: dict[tuple[int, int, tuple[int, int], int], set[int]] = {}
    for c in states:
        if c[5] != 0:
            continue
        key = (c[0], c[1], cl_pair(c), c[5])
        by_base.setdefault(key, set()).add(c[4])
    return all(vals == {0, 1} for vals in by_base.values())


def has_one_bit_lq_separator_7() -> bool:
    """Can any single bit separate CL=00 from the three colour states?"""

    sector_states = list(product([0, 1], repeat=2))
    lepton = {(0, 0)}
    for bit in (0, 1):
        for value in (0, 1):
            selected = {cl for cl in sector_states if cl[bit] == value}
            if selected == lepton:
                return True
    return False


def fixed_r4_classification() -> tuple[int, int]:
    masks = target_a.all_sector_masks()
    cases = [(mask, target_a.canonical_r4(mask)) for mask in masks]
    counts = target_a.classify_cases(cases)
    return counts["anomaly_free"], counts["chiral_anomaly_free"]


def wide_r4_survivors() -> list[tuple[frozenset[target_a.Triple], target_a.Exclusion]]:
    masks = target_a.all_sector_masks()
    exclusions = [("lq", lq, i3, 1) for lq in (0, 1) for i3 in (0, 1)]
    survivors = []
    for mask in masks:
        for exclusion in exclusions:
            states = target_a.content(mask, exclusion)
            if (
                target_a.anomaly_free(states)
                and not target_a.vectorlike(states)
                and target_a.su2_complete(states)
            ):
                survivors.append((mask, exclusion))
    return survivors


def role_map_check() -> tuple[int, int, int]:
    role_maps = target_a.role_maps()
    anomaly_free = 0
    sm_isomorphic = 0
    for rm in role_maps:
        states = target_a.content(target_a.SM_SECTOR, target_a.SM_EXCLUSION, rm)
        if target_a.anomaly_free(states, rm):
            anomaly_free += 1
        if target_a.sm_isomorphic(states, rm):
            sm_isomorphic += 1
    return len(role_maps), anomaly_free, sm_isomorphic


def main() -> int:
    widths = primitive_record_lower_bound()
    lower_bound = sum(widths.values())

    print("Target-A minimal local-record architecture audit")
    print("\n[1] primitive local-record lower bound")
    for name, width in widths.items():
        print(f"    {name:18s}: {width}")
    print(f"    lower bound       : {lower_bound} bits")
    for nbits in range(1, 9):
        print(f"    assignments at {nbits} bits: {count_disjoint_role_assignments(nbits)}")

    states7 = content_7()
    print("\n[2] explicit 7-bit compressed encoding")
    print(f"    active states     : {len(states7)}")
    print(f"    anomalies         : {anomalies_7(states7)}")
    print(f"    vectorlike        : {vectorlike_7(states7)}")
    print(f"    SU(2)-complete    : {su2_complete_7(states7)}")
    print(f"    one-bit LQ split  : {has_one_bit_lq_separator_7()}")
    print("    verdict           : valid compact encoding, not comparable local-record architecture")

    fixed_af, fixed_chiral_af = fixed_r4_classification()
    survivors = wide_r4_survivors()
    role_total, role_af, role_iso = role_map_check()

    print("\n[3] 8-bit comparable local-record architecture")
    print(f"    fixed-R4 anomaly-free rules        : {fixed_af}")
    print(f"    fixed-R4 chiral anomaly-free rules : {fixed_chiral_af}")
    print(f"    wide-R4 chiral SU2 AF survivors    : {len(survivors)}")
    for mask, exclusion in survivors:
        print(f"      survivor: sector={tuple(sorted(mask))} exclusion={exclusion}")
    print(f"    role maps tested                   : {role_total}")
    print(f"    role-map anomaly-free              : {role_af}")
    print(f"    role-map SM-isomorphic             : {role_iso}")

    assert lower_bound == 8
    assert count_disjoint_role_assignments(7) == 0
    assert count_disjoint_role_assignments(8) > 0

    assert len(states7) == 45
    assert anomalies_7(states7) == (0, 0, 0, 0, 0, 0)
    assert not vectorlike_7(states7)
    assert su2_complete_7(states7)
    assert not has_one_bit_lq_separator_7()

    assert fixed_af == 2
    assert fixed_chiral_af == 1
    assert len(survivors) == 2
    assert all(target_a.sm_rule_like(mask, exclusion) for mask, exclusion in survivors)
    assert role_total == 3360
    assert role_af == 3360
    assert role_iso == 3360

    print("\nVERDICT:")
    print("  A 7-bit SM encoding exists, so bare bit-count minimality is false.")
    print("  Under comparable local-record axioms, the primitive-role lower bound is")
    print("  8 bits, and the 8-bit Target-A classification isolates the SM active")
    print("  rule up to LQ-name and colour/name relabelling.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
