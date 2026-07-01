#!/usr/bin/env python3
r"""v_phase6_complete_cell_billing_scalar_pointer.py

Phase 6A: complete-cell billing as a scalar protected-pointer theorem.

Phase 4 proved:

    if EWSB bills the complete filled cell, the leading channel is alpha_0^8.

Phase 5 proved:

    finite R4 syndrome/repair does not itself force that bill.

This script tests the strongest finite replacement: if the electroweak order
parameter is a *single stable local pointer record* in the minimal balanced
record cell, and it is a scalar under the cube/octahedral cell symmetry, which
record can it be?

The cell code is RM(1,3) = [8,4,4], represented as affine Boolean functions on
the vertices of Q_3.  The physical O_h action on the cube vertices is generated
by coordinate permutations and coordinate flips, acting on codewords by
permuting their eight positions.

Result:

    O_h-fixed RM(1,3) pointer records are exactly 00000000 and 11111111.

Thus a non-vacuum scalar protected endpoint is uniquely the filled byte.  Its
formation from the empty endpoint requires all eight records to be written, so
the bill is alpha_0^8.

This is a theorem conditional on a precise physical identification:

    Higgs VEV = scalar protected pointer endpoint of the local matter cell.

It does not follow from R4 syndrome repair alone, and it does not derive the
SM Higgs Lagrangian.  It is the finite record-theoretic content of the
"complete-cell billing" premise.
"""

from __future__ import annotations

from itertools import permutations, product
import sys


ok = True
N = 8


def check(name: str, cond: bool) -> None:
    global ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    ok = ok and bool(cond)


VERTS = list(product((0, 1), repeat=3))


def bitstr(word: tuple[int, ...]) -> str:
    return "".join(str(x) for x in word)


def rm13_words() -> list[tuple[int, ...]]:
    """RM(1,3): affine functions a0 + a.x on F_2^3."""
    words = []
    for a0 in (0, 1):
        for a in product((0, 1), repeat=3):
            words.append(tuple((a0 ^ (a[0] & x[0]) ^ (a[1] & x[1]) ^ (a[2] & x[2])) for x in VERTS))
    return words


def cube_automorphisms() -> list[tuple[int, ...]]:
    """The 48 cube vertex permutations: x -> P x + b over F_2."""
    out = []
    for perm in permutations(range(3)):
        for b in product((0, 1), repeat=3):
            p = []
            for x in VERTS:
                y = tuple(x[perm[i]] ^ b[i] for i in range(3))
                p.append(VERTS.index(y))
            out.append(tuple(p))
    return sorted(set(out))


def act(word: tuple[int, ...], perm: tuple[int, ...]) -> tuple[int, ...]:
    """Pullback action on vertex-labelled records."""
    out = [0] * len(word)
    for old, new in enumerate(perm):
        out[new] = word[old]
    return tuple(out)


def orbit(word: tuple[int, ...], group: list[tuple[int, ...]]) -> set[tuple[int, ...]]:
    return {act(word, g) for g in group}


def main() -> int:
    print("=" * 96)
    print("PHASE 6A: COMPLETE-CELL BILLING AS SCALAR PROTECTED-POINTER THEOREM")
    print("=" * 96)

    words = rm13_words()
    group = cube_automorphisms()
    print("\n[0] Cell and symmetry")
    print(f"  RM(1,3) words: {len(words)}")
    print(f"  cube/O_h vertex permutations: {len(group)}")
    weights = sorted({sum(w) for w in words})
    print(f"  codeword weights: {weights}")
    check("RM(1,3) has 16 words", len(words) == 16)
    check("cube symmetry group has 48 vertex permutations", len(group) == 48)
    check("weight enumerator is 1 + 14 x^4 + x^8",
          [sum(1 for w in words if sum(w) == k) for k in range(N + 1)] == [1, 0, 0, 0, 14, 0, 0, 0, 1])

    print("\n[1] Scalar pointer records")
    fixed = [w for w in words if all(act(w, g) == w for g in group)]
    for w in fixed:
        print(f"  fixed: {bitstr(w)}  weight={sum(w)}")
    check("only vacuum and filled byte are O_h-fixed pointer codewords",
          sorted(bitstr(w) for w in fixed) == ["00000000", "11111111"])
    nonvac = [w for w in fixed if any(w)]
    check("unique non-vacuum scalar pointer endpoint is the filled byte",
          len(nonvac) == 1 and sum(nonvac[0]) == 8)

    print("\n[2] Why weight-4 records do not bill a scalar VEV endpoint")
    weight4 = [w for w in words if sum(w) == 4]
    orbits = {}
    for w in weight4:
        key = tuple(sorted(bitstr(x) for x in orbit(w, group)))
        orbits.setdefault(key, []).append(w)
    print(f"  weight-4 protected records: {len(weight4)}")
    orbit_sizes = sorted(len(k) for k in orbits)
    print(f"  O_h orbit count among weight-4 records: {len(orbits)}")
    print(f"  orbit sizes: {orbit_sizes}")
    check("weight-4 records are protected but orientation-bearing, not scalar pointer endpoints",
          len(weight4) == 14 and len(orbits) == 3 and min(orbit_sizes) > 1 and all(w not in fixed for w in weight4))

    print("\n[3] Billing implication")
    print("  Empty scalar endpoint  = 00000000.")
    print("  Filled scalar endpoint = 11111111.")
    print("  The transition changes all eight record bits, so a per-record alpha_0")
    print("  monitored-service bill gives alpha_0^8.  Lower-k records are real")
    print("  protected records, but they carry an orientation/hyperplane label and")
    print("  cannot be the scalar cell endpoint.")
    check("scalar protected-pointer transition has Hamming support 8",
          sum(a != b for a, b in zip((0,) * N, nonvac[0])) == 8)

    print("\nVERDICT")
    print("  Complete-cell billing is finite-theorem closed under the scalar protected")
    print("  pointer identification of the Higgs order parameter.  The theorem does")
    print("  not come from R4 repair; it comes from the RM(1,3) record cell plus")
    print("  O_h scalarity.  The remaining physics premise is now sharply named:")
    print("  the Higgs VEV must be this scalar protected pointer endpoint, rather")
    print("  than a lower orientation-bearing record or a local R4 repair current.")

    if ok:
        print("\nALL CHECKS PASSED")
        return 0
    print("\nCHECKS FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
