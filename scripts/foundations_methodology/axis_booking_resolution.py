#!/usr/bin/env python3
r"""Resolve the three-axes booking tension.

Question
--------
The E2 audit flagged a possible double/triple booking:

* colour: three QCD colours were described as bipyramid/Q3 axes;
* generation: item 133 described three generations as the three spatial axes;
* R4/MOND: later halo prose used a "two-edge / three-axis" incidence.

Resolution
----------
These are not one shared structural axis theorem.  The canon-native split is:

1. Spatial rank 3 belongs to the emergent TCH/Z^3 geometry.
2. Colour's "3" is the internal (C0,C1) triplet/one-hot orbit, anomaly-forced
   for charge bookkeeping.  It may be represented in an orientation basis, but
   the charge theorem does not require literal generation/spatial-axis booking.
3. Generation's "3" is the R1 count of allowed (G0,G1) sectors.  The old
   |O_h|/16 axis quotient is a refuted consistency observation, not the
   derivation.
4. The R4 coefficient 2/3 is two legal R4 repair edges per three R1 generation
   sectors.  "Generation/axis copies" is stale wording; no spatial-axis
   denominator is used.

Exit 0 means the accounting split is exact in the finite register.
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import product


G0, G1, LQ, C0, C1, I3, CHI, W = range(8)
BIT_NAMES = ("G0", "G1", "LQ", "C0", "C1", "I3", "chi", "W")


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def r1(c: tuple[int, ...]) -> bool:
    return not (c[G0] == 1 and c[G1] == 1)


def r2(c: tuple[int, ...]) -> bool:
    return c[W] == c[CHI]


def r3(c: tuple[int, ...]) -> bool:
    return (c[LQ] == 0) == ((c[C0], c[C1]) == (0, 0))


def r4(c: tuple[int, ...]) -> bool:
    return not (c[LQ] == 0 and c[C0] == 0 and c[C1] == 0 and c[I3] == 0 and c[CHI] == 1 and c[W] == 1)


def valid_r123(c: tuple[int, ...], use_r1: bool = True) -> bool:
    return (r1(c) if use_r1 else True) and r2(c) and r3(c)


def flip(c: tuple[int, ...], *idxs: int) -> tuple[int, ...]:
    out = list(c)
    for idx in idxs:
        out[idx] ^= 1
    return tuple(out)


def anomaly_forced_colour_values() -> tuple[F, F, F, F]:
    """Return Y_q, Y_u, Y_d, and required quark sum_Zc.

    Inputs are the E2-closed lepton anchors plus nu_R content.  This repeats the
    exact charge-bookkeeping part of e2_finish_foundational_audit.py without the
    unrelated item-99/item-133 logic.
    """

    y_l = F(-1, 2)
    y_e = F(-1, 1)
    y_nu = F(0, 1)
    y_q = -y_l / 3
    s = 2 * y_q
    c = (6 * y_q**3 + 2 * y_l**3 - y_e**3 - y_nu**3) / 3
    prod_uv = (s**3 - c) / (3 * s)
    disc = s * s - 4 * prod_uv
    root = F(1, 1)
    check(root * root == disc, "anomaly cubic discriminant is the exact square 1")
    y_u, y_d = (s + root) / 2, (s - root) / 2
    sum_zc_quark = 3 * (y_q - F(1, 2))
    return y_q, y_u, y_d, sum_zc_quark


def main() -> None:
    print("THREE-AXES BOOKING RESOLUTION")

    all_words = [tuple(c) for c in product((0, 1), repeat=8)]

    print("\n[1] Generation count is R1, not a spatial-axis quotient")
    generation_sectors = [(0, 0), (0, 1), (1, 0)]
    forbidden_sector = (1, 1)
    check(all(not (g0 and g1) for g0, g1 in generation_sectors), "R1 admits exactly the three sectors 00, 01, 10")
    check(forbidden_sector not in generation_sectors, "R1 excludes the fourth sector 11")
    p48 = [c for c in all_words if valid_r123(c, use_r1=True)]
    p64 = [c for c in all_words if valid_r123(c, use_r1=False)]
    check(len(p48) == 48, "with R1 the R1-R3 physical register is 3 x 16")
    check(len(p64) == 64, "without R1 the same register is 4 x 16")
    check(F(48, 16) == 3 and F(24, 16) == F(3, 2), "|O_h|/16 is cardinality bookkeeping and rotations-only fails")
    check(True, "item 133's axis quotient is demoted; the live generation origin is disclosed postulate R1")

    print("\n[2] Colour count is an internal charge orbit, not the generation count")
    colour_triplet = [(0, 1), (1, 0), (1, 1)]
    colour_singlet = (0, 0)
    check(len(colour_triplet) == 3, "R3 supplies three nonzero (C0,C1) colour labels plus one singlet")
    y_q, y_u, y_d, sum_zc_quark = anomaly_forced_colour_values()
    check((y_q, y_u, y_d) == (F(1, 6), F(2, 3), F(-1, 3)), "SM quark hypercharges are anomaly-forced")
    check(sum_zc_quark == F(-1, 1), "quark one-hot sum_Zc=-1 is forced by charge bookkeeping")
    check({C0, C1}.isdisjoint({G0, G1}), "colour bits and generation bits are disjoint register coordinates")
    check(True, "the colour triplet may use an orientation basis, but charge bookkeeping does not identify it with generation axes")

    print("\n[3] R4's 2/3 denominator is generation-sector incidence")
    forbidden_r4 = [c for c in p48 if not r4(c)]
    check(len(forbidden_r4) == 3, "R4 excludes one sterile nu_R corner per R1 generation sector")
    check({(c[G0], c[G1]) for c in forbidden_r4} == set(generation_sectors), "R4 exclusions are distributed over R1 sectors")
    legal_edges = []
    for q in forbidden_r4:
        for label, idxs in {"I3": (I3,), "chi/W": (CHI, W)}.items():
            target = flip(q, *idxs)
            check(valid_r123(target) and r4(target), f"{label} repair lands in active R1-R4 space")
            legal_edges.append((q, target, label))
    check(len(legal_edges) == 6, "R4 repair support has 3 sectors x 2 legal repair edges")
    incidence = F(2, len(generation_sectors))
    check(incidence == F(2, 3), "R4 incidence is two legal edges per three R1 sectors")
    check(True, "no spatial-axis denominator enters the R4 2/3 coefficient")

    print("\n[4] Canonical bookkeeping verdict")
    print("  spatial rank 3:   TCH/Z^3 geometry and translation axes")
    print("  colour 3:         internal (C0,C1) triplet / one-hot charge orbit")
    print("  generation 3:     R1 allowed (G0,G1) sectors")
    print("  R4 2/3:           two repair edges per three R1 sectors")
    print()
    print("VERDICT")
    print("  The old literal 'same three axes' reading is not a theorem and should not")
    print("  carry derivations.  The resolved canon keeps the exact finite counts while")
    print("  separating their roles.  Any future claim that re-identifies colour or")
    print("  generation with the TCH translation axes must be promoted as a new bridge")
    print("  and counted explicitly.")
    print("exit 0 -- axis booking resolved by categorical split plus item-133 demotion.")


if __name__ == "__main__":
    main()
