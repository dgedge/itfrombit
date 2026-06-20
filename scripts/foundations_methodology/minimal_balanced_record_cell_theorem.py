#!/usr/bin/env python3
r"""MINIMAL BALANCED RECORD CELL THEOREM.

Question:
    Among finite binary local record systems with repeatable readout, commuting
    read/write stabilizers, local-tomography compatibility, and distance-4
    erasure protection, is the unique minimal balanced cell the self-dual
    doubly-even [8,4,4] code?

Answer proved here, with the hypotheses made explicit:

    If "balanced record cell" means a binary Type-II CSS record cell:

      * one binary local register per site;
      * a read/write CSS pair X(C), Z(C) built from the same classical code C;
      * all read/write stabilizers commute, so C is self-orthogonal;
      * the cell is balanced/complete, so dim C = n/2 and C = C^\perp;
      * reset/readout is even under closed record-pair reversal, so C is
        doubly even;
      * erasure protection requires d(C) >= 4;

    then the minimal block length is n = 8, and the unique cell at n = 8
    up to coordinate relabelling is the extended Hamming / RM(1,3)
    self-dual doubly-even [8,4,4] code.

What this does NOT prove:
    It does not derive binary records, complex Hilbert space, or local
    tomography from nothing.  Those are the upstream Record Reconstruction
    premises.  It closes the finite coding-theory rung once the balanced
    Type-II CSS record-cell setting is accepted.

exit 0 = no smaller Type-II distance-4 self-dual binary code exists; every
         length-8 Type-II distance-4 code is permutation-equivalent to RM(1,3);
         RM(1,3) has weight enumerator 1 + 14 y^4 + y^8; and X(C),Z(C) form
         a commuting full-rank [[8,0,4]] read/write stabilizer cell.
"""

from __future__ import annotations

import itertools
from collections import Counter


def wt(v: int) -> int:
    return v.bit_count()


def dot(a: int, b: int) -> int:
    return wt(a & b) & 1


def gf2_rank(vectors: list[int], n: int) -> int:
    """Rank over F_2 for integer bit-vectors."""
    basis = [0] * n
    rank = 0
    for x in vectors:
        y = x
        while y:
            p = y.bit_length() - 1
            if basis[p]:
                y ^= basis[p]
            else:
                basis[p] = y
                rank += 1
                break
    return rank


def span(gens: list[int]) -> frozenset[int]:
    out = {0}
    for g in gens:
        out |= {x ^ g for x in list(out)}
    return frozenset(out)


def min_distance(code: frozenset[int]) -> int:
    return min(wt(c) for c in code if c)


def self_orthogonal(code: frozenset[int]) -> bool:
    words = list(code)
    return all(dot(a, b) == 0 for a in words for b in words)


def doubly_even(code: frozenset[int]) -> bool:
    return all(wt(c) % 4 == 0 for c in code)


def is_type_ii_distance4(code: frozenset[int], n: int) -> bool:
    return (
        len(code) == 2 ** (n // 2)
        and self_orthogonal(code)
        and doubly_even(code)
        and min_distance(code) >= 4
    )


def permute_bits(v: int, perm: tuple[int, ...]) -> int:
    """Return v with old bit i sent to new bit perm[i]."""
    out = 0
    for i, j in enumerate(perm):
        if (v >> i) & 1:
            out |= 1 << j
    return out


def canonical_under_s8(code: frozenset[int]) -> tuple[int, ...]:
    best: tuple[int, ...] | None = None
    for perm in itertools.permutations(range(8)):
        image = tuple(sorted(permute_bits(c, perm) for c in code))
        if best is None or image < best:
            best = image
    assert best is not None
    return best


def exhaustive_small_lengths() -> None:
    print("[1] no smaller balanced Type-II distance-4 cell:")
    for n in (2, 4, 6):
        k = n // 2
        seen: set[frozenset[int]] = set()
        hits: list[frozenset[int]] = []
        vectors = list(range(1, 1 << n))
        for gens in itertools.combinations(vectors, k):
            if gf2_rank(list(gens), n) != k:
                continue
            code = span(list(gens))
            if code in seen:
                continue
            seen.add(code)
            if is_type_ii_distance4(code, n):
                hits.append(code)
        print(f"    n={n}: checked {len(seen)} k={k} subspaces; Type-II d>=4 hits = {len(hits)}")
        assert not hits
    print("    -> minimal Type-II self-dual doubly-even distance-4 length is at least 8.")


def length8_uniqueness() -> frozenset[int]:
    print("\n[2] length-8 uniqueness up to relabelling:")
    ones = (1 << 8) - 1
    weight4 = [v for v in range(1, 1 << 8) if wt(v) == 4]
    codes: set[frozenset[int]] = set()

    # Every Type-II self-dual code is even, so 1^n is in C^\perp = C.
    # At n=8, all nonzero words besides 1^8 have weight 4.  Three weight-4
    # generators plus 1^8 are therefore enough to enumerate the class.
    for triple in itertools.combinations(weight4, 3):
        gens = [ones, *triple]
        if gf2_rank(gens, 8) != 4:
            continue
        code = span(gens)
        if is_type_ii_distance4(code, 8):
            codes.add(code)

    classes = {canonical_under_s8(code) for code in codes}
    print(f"    labelled Type-II [8,4,4] codes found = {len(codes)}")
    print(f"    coordinate-permutation equivalence classes = {len(classes)}")
    assert len(codes) == 30
    assert len(classes) == 1
    print("    -> unique length-8 balanced cell: extended Hamming / RM(1,3).")
    return frozenset(next(iter(codes)))


def rm13_code() -> tuple[list[int], frozenset[int]]:
    # positions 0..7 are cube vertices; rows are {1, b2, b1, b0}.
    gens = [
        0b11111111,  # all ones
        0b11110000,  # b2
        0b11001100,  # b1
        0b10101010,  # b0
    ]
    return gens, span(gens)


def verify_rm13() -> None:
    print("\n[3] canonical RM(1,3) representative:")
    gens, code = rm13_code()
    gram_zero = all(dot(a, b) == 0 for a in gens for b in gens)
    rank = gf2_rank(gens, 8)
    enumerator = Counter(wt(c) for c in code)
    print(f"    rank = {rank}; |C| = {len(code)}; C subset C^perp = {gram_zero}")
    print(f"    weight enumerator = {dict(sorted(enumerator.items()))}")
    print(f"    minimum distance = {min_distance(code)}")
    assert rank == 4
    assert len(code) == 16
    assert gram_zero
    assert len(code) == 2 ** (8 - rank)  # self-orthogonal + dim n/2 => self-dual
    assert enumerator == Counter({0: 1, 4: 14, 8: 1})
    assert min_distance(code) == 4
    assert doubly_even(code)


def verify_css_read_write() -> None:
    print("\n[4] read/write CSS stabilizer compatibility:")
    gens, _ = rm13_code()
    # Symplectic vectors are (x | z).  X(g)=(g|0), Z(g)=(0|g).
    stabs = [(g, 0) for g in gens] + [(0, g) for g in gens]
    noncommuting = []
    for i, (x1, z1) in enumerate(stabs):
        for j, (x2, z2) in enumerate(stabs):
            if i >= j:
                continue
            if dot(x1, z2) ^ dot(z1, x2):
                noncommuting.append((i, j))
    symplectic_as_ints = [x | (z << 8) for x, z in stabs]
    stab_rank = gf2_rank(symplectic_as_ints, 16)
    print(f"    stabilizer generators = {len(stabs)}; symplectic rank = {stab_rank}")
    print(f"    non-commuting read/write generator pairs = {len(noncommuting)}")
    assert not noncommuting
    assert stab_rank == 8
    print("    -> X(C), Z(C) define a commuting full-rank [[8,0,4]] record cell.")


def bridge_readout() -> None:
    print("\n[5] record-alphabet bridge bookkeeping (not part of uniqueness proof):")
    records = 16
    symmetric_pairs = records * (records + 1) // 2
    print(f"    classical record words |C| = {records}")
    print(f"    unordered-with-diagonal record pairs Sym^2(16) = {symmetric_pairs}")
    print(f"    + idle channel => {symmetric_pairs + 1} = 137")
    assert symmetric_pairs == 136
    assert symmetric_pairs + 1 == 137
    print("    -> this recovers the alpha0 channel count only after the monitored-pair/idle convention is accepted.")


def main() -> None:
    print("MINIMAL BALANCED RECORD CELL THEOREM")
    exhaustive_small_lengths()
    length8_uniqueness()
    verify_rm13()
    verify_css_read_write()
    bridge_readout()
    print(
        """
[verdict]
  Conditional theorem proved:
    In the binary balanced Type-II CSS record-cell class, distance-4 erasure
    protection forces minimal length n=8, and the unique minimal cell up to
    relabelling is the self-dual doubly-even [8,4,4] extended Hamming/RM(1,3)
    code.  Its read/write stabilizers commute and form the full [[8,0,4]]
    syndrome-record cell.

  Boundary of the theorem:
    The result does not derive the upstream reconstruction premises: binary
    local records, complex/local-tomographic Hilbert composition, or the
    monitored service-pair convention that turns 16 records into 136+1=137.
    It closes the finite balanced-cell rung once those premises are admitted.
exit 0"""
    )


if __name__ == "__main__":
    main()
