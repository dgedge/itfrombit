#!/usr/bin/env python3
"""Item 86 / DRIFT M9: split the retired universal-2/9 claim.

The old claim joined three things under one "universal 2/9" banner:

  1. chiral QED Q^3 trace per generation;
  2. Weinberg sin^2(theta_W)=2/9;
  3. Koide / Berry phase delta=2/9.

The recoverable result is only (1).  This script rebuilds the actual R1-R3
single-generation 16-state register and verifies, with exact Fraction
arithmetic, that the left- and right-chiral charge-cube traces are both -2/9.
The two old comparison legs are explicitly kept out of the locked result:
Weinberg is charge-forced to 3/8 at the GUT normalisation, and Koide delta is a
phase/defect-ratio claim already separated by DRIFT M9.
"""

from __future__ import annotations

import itertools
from fractions import Fraction as F


G0, G1, LQ, C0, C1, I3, CHI, W = range(8)
BIT_NAMES = ("G0", "G1", "LQ", "C0", "C1", "I3", "chi", "W")


def r1(c: tuple[int, ...]) -> bool:
    return not (c[G0] and c[G1])


def r2(c: tuple[int, ...]) -> bool:
    return c[W] == c[CHI]


def r3(c: tuple[int, ...]) -> bool:
    return (c[LQ] == 0) == ((c[C0], c[C1]) == (0, 0))


def r4(c: tuple[int, ...]) -> bool:
    return not (c[LQ] == 0 and c[I3] == 0 and c[CHI] == 1)


def zf(c: tuple[int, ...]) -> int:
    return 1 if c[I3] == 0 else -1


def sum_zc(c: tuple[int, ...]) -> int:
    return -3 if (c[C0], c[C1]) == (0, 0) else -1


def charge(c: tuple[int, ...]) -> F:
    # ANCHOR section 2.8 charge map: Q = 1/2 Z_f + 1/3 sum_i Z_ci + 1/2 Z_p, with Z_p=+1 for matter.
    return F(1, 2) * zf(c) + F(1, 3) * sum_zc(c) + F(1, 2)


def species(c: tuple[int, ...]) -> str:
    if c[LQ] == 0:
        return "nu" if c[I3] == 0 else "e"
    return "u" if c[I3] == 0 else "d"


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def fmt_state(c: tuple[int, ...]) -> str:
    return "".join(str(x) for x in c)


def main() -> None:
    print("ITEM 86 UNIVERSAL-2/9 RECOVERY SPLIT")

    all_words = list(itertools.product((0, 1), repeat=8))
    gen0 = [c for c in all_words if r1(c) and r2(c) and r3(c) and (c[G0], c[G1]) == (0, 0)]
    active = [c for c in gen0 if r4(c)]
    print("\n[1] Canonical one-generation register")
    print(f"  R1-R3 states in generation 00: {len(gen0)}")
    print(f"  R1-R4 active SM states       : {len(active)}")
    check(len(gen0) == 16, "R1-R3 give one SO(10)-sized 16-state generation including nu_R")
    check(len(active) == 15, "R4 removes exactly one sterile nu_R state from the active basis")
    nu_r = [c for c in gen0 if not r4(c)]
    check(len(nu_r) == 1 and species(nu_r[0]) == "nu" and nu_r[0][CHI] == 1, "excluded state is the sterile nu_R corner")
    print(f"  excluded state: {fmt_state(nu_r[0])} ({', '.join(BIT_NAMES[i] for i, b in enumerate(nu_r[0]) if b)})")

    print("\n[2] Exact chiral Q^3 arithmetic")
    rows = []
    for chi in (0, 1):
        side = "L" if chi == 0 else "R"
        states = [c for c in gen0 if c[CHI] == chi]
        trace = sum(charge(c) ** 3 for c in states)
        rows.append((side, trace))
        display = ", ".join(f"{species(c)}:{charge(c)}" for c in states)
        print(f"  chi={chi} ({side}) states: {display}")
        print(f"  sum Q_{side}^3 = {trace}")
        check(trace == F(-2, 9), f"sum Q_{side}^3 is exactly -2/9")
    check(rows[0][1] - rows[1][1] == 0, "vector-like QED chiral cancellation sum Q_L^3 - sum Q_R^3 = 0")

    print("\n[3] Why this leg is locked but not universal")
    lhs = F(-1) ** 3 + 3 * F(2, 3) ** 3 + 3 * F(-1, 3) ** 3
    print(f"  expanded identity: (-1)^3 + 3(2/3)^3 + 3(-1/3)^3 = {lhs}")
    check(lhs == F(-2, 9), "expanded lepton+colour arithmetic equals -2/9")

    tr_t3_sq = F(2, 1)
    tr_y2 = F(10, 3)
    gut_sin2 = tr_t3_sq / (tr_t3_sq + tr_y2)
    print(f"  Weinberg charge trace: Tr(T3^2)/(Tr(T3^2)+Tr((Y/2)^2)) = {gut_sin2}")
    check(gut_sin2 == F(3, 8), "charge-forced GUT Weinberg value is 3/8, not 2/9")
    check(gut_sin2 != F(2, 9), "Weinberg 2/9 leg is retired")

    old_legs = {
        "chiral_Q3_trace": "LOCKED exact arithmetic identity",
        "Weinberg_sin2": "RETIRED; charge-forced 3/8 plus RG running",
        "Koide_delta": "NOT LOCKED as exact rational; phase/defect-ratio status per DRIFT M9",
    }
    for name, status in old_legs.items():
        print(f"  {name:18s}: {status}")
    locked = [name for name, status in old_legs.items() if status.startswith("LOCKED")]
    check(locked == ["chiral_Q3_trace"], "only the chiral Q^3 trace survives as a locked 2/9 result")

    print(
        "\nVERDICT\n"
        "  The universal-2/9 banner is not recoverable as a cross-sector invariant.\n"
        "  The recoverable theorem is narrower and cleaner: the canonical R1-R3\n"
        "  charge register forces sum Q_L^3 = sum Q_R^3 = -2/9 exactly per\n"
        "  generation, with vector-like QED cancellation.  Weinberg and Koide\n"
        "  are separate retired/conditional legs and must not be cited as support\n"
        "  for a universal invariant.\n"
        "exit 0 -- chiral Q^3 leg recovered; universal claim remains retracted."
    )


if __name__ == "__main__":
    main()
