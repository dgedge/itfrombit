#!/usr/bin/env python3
r"""v_phase5_ewsb_billing_map_gate.py

Physical billing-map gate for the electroweak VEV route.

Phase 4 proved a finite theorem:

    if the EWSB observable is the complete filled-cell record, k=8 is forced.

This script tests the stronger question: does R4 enforcement itself force that
billing map?  It does not.  The R4 syndrome is a low-bit condition inside the
R1--R3 register, and the R4 repair channel has local low-Hamming repairs.  So
"remove the sterile nu_R corner" and "create the complete filled matter cell"
are different billing maps:

    R4 syndrome / repair       -> low-order local record
    filled-cell formation      -> all-eight product, alpha_0^8

Exit 0 means the gate is sharpened.  The byte-power route survives only with
the explicit physical identification that the Higgs VEV bills complete-cell
formation, not merely local R4 syndrome detection/repair.
"""

from __future__ import annotations

from itertools import combinations, product
import sys


G0, G1, LQ, C0, C1, I3, CHI, W = range(8)
BIT_NAMES = ("G0", "G1", "LQ", "C0", "C1", "I3", "chi", "W")
ok = True


def check(name: str, cond: bool) -> None:
    global ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    ok = ok and bool(cond)


def r1(c: tuple[int, ...]) -> bool:
    return not (c[G0] == 1 and c[G1] == 1)


def r2(c: tuple[int, ...]) -> bool:
    return c[W] == c[CHI]


def r3(c: tuple[int, ...]) -> bool:
    return (c[LQ] == 0) == ((c[C0], c[C1]) == (0, 0))


def r4(c: tuple[int, ...]) -> bool:
    return not (c[LQ] == 0 and c[C0] == 0 and c[C1] == 0 and c[I3] == 0 and c[CHI] == 1 and c[W] == 1)


def valid_r123(c: tuple[int, ...]) -> bool:
    return r1(c) and r2(c) and r3(c)


def valid_active(c: tuple[int, ...]) -> bool:
    return valid_r123(c) and r4(c)


def flip(c: tuple[int, ...], *idxs: int) -> tuple[int, ...]:
    out = list(c)
    for idx in idxs:
        out[idx] ^= 1
    return tuple(out)


def bitstr(c: tuple[int, ...]) -> str:
    return "".join(str(x) for x in c)


def minimal_r4_syndrome_conjunctions(words: list[tuple[int, ...]], active: list[tuple[int, ...]]) -> list[tuple[int, ...]]:
    """Find smallest common fixed-bit subsets that match all R4-forbidden words and no active word."""
    common = [i for i in range(8) if len({w[i] for w in words}) == 1]
    for size in range(1, len(common) + 1):
        winners = []
        for subset in combinations(common, size):
            values = {i: words[0][i] for i in subset}
            if all(all(w[i] == v for i, v in values.items()) for w in words) and not any(
                all(a[i] == v for i, v in values.items()) for a in active
            ):
                winners.append(subset)
        if winners:
            return winners
    return []


def main() -> int:
    print("=" * 96)
    print("PHASE 5A: PHYSICAL EWSB BILLING-MAP GATE")
    print("=" * 96)

    all_words = [tuple(c) for c in product((0, 1), repeat=8)]
    r123 = [c for c in all_words if valid_r123(c)]
    active = [c for c in r123 if r4(c)]
    forbidden = [c for c in r123 if not r4(c)]
    print("\n[0] Register census")
    print(f"  R1--R3 words: {len(r123)}")
    print(f"  active R1--R4 words: {len(active)}")
    print(f"  R4-forbidden sterile corners: {len(forbidden)}")
    for c in forbidden:
        print(f"    {bitstr(c)}  gen={(c[G0], c[G1])}")
    check("R4 excludes exactly one sterile corner per generation", len(forbidden) == 3)

    print("\n[1] R4 syndrome is low-bit, not all-eight")
    winners = minimal_r4_syndrome_conjunctions(forbidden, active)
    print(f"  minimal conjunction size excluding active words while selecting all R4 corners: {len(winners[0])}")
    for subset in winners[:8]:
        label = ", ".join(f"{BIT_NAMES[i]}={forbidden[0][i]}" for i in subset)
        print(f"    {label}")
    check("R4 syndrome can be detected by a strict subset of the byte", len(winners[0]) < 8)
    check("minimal common R4 syndrome has size three in the R1--R3 register", len(winners[0]) == 3)

    print("\n[2] Local R4 repair edges exist")
    repair_specs = {"I3": (I3,), "locked chi/W": (CHI, W)}
    repair_lengths = []
    for q in forbidden:
        for label, bits in repair_specs.items():
            target = flip(q, *bits)
            dist = sum(a != b for a, b in zip(q, target))
            repair_lengths.append(dist)
            print(f"  {bitstr(q)} --{label:12s}--> {bitstr(target)}  distance={dist}  active={valid_active(target)}")
            check(f"{label} repair lands in active space", valid_active(target))
    check("R4 repair has lower-than-byte Hamming support", min(repair_lengths) == 1 and max(repair_lengths) == 2)

    print("\n[3] Billing-map implication")
    print("  If the action bills only R4 syndrome detection, the natural support is 3 bits.")
    print("  If the action bills local R4 repair, the support is 1 or 2 bit flips.")
    print("  Neither route gives alpha_0^8.  The alpha_0^8 route is instead the")
    print("  complete-cell formation bill: all eight record modes must be created.")
    check("R4 enforcement alone cannot force the alpha_0^8 bill", len(winners[0]) != 8 and max(repair_lengths) < 8)

    print("\nVERDICT")
    print("  The physical billing map is not derivable from finite R4 enforcement alone.")
    print("  R4 supplies the reason the filled-cell channel matters, but its own local")
    print("  syndrome/repair algebra is too cheap.  Therefore the remaining theorem is")
    print("  precisely the physical identification: the Higgs VEV must bill complete")
    print("  matter-cell formation, not just R4 syndrome detection or repair.")

    if ok:
        print("\nALL CHECKS PASSED")
        return 0
    print("\nCHECKS FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
