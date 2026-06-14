#!/usr/bin/env python3
r"""Finite-cell construction of the section-11.5 horizon isometry V.

Question attacked:
  Canon left item 11 open as "explicitly constructing V from the [8,4,4]
  stabilizer generators against a Schwarzschild background".

What this script closes:
  * the finite [8,4,4] / Q3 stabilizer-ledger isometry;
  * the exact role of H_vacuum: it is the one global-complement/intercept
    degree that the 12-edge strain syndrome cannot see;
  * compatibility with the Schwarzschild/KMS strain weighting used by item 10,
    because the strain F is a function of the emitted syndrome record.

What remains open:
  * the radial Schwarzschild steady state / localized-mass scheduler that says
    which cells enter this finite map, with which rates, across the horizon
    shell.  The static algebraic V is constructed; the full V_Sch(M) channel
    over a horizon geometry is still a QEC-steady-state problem.

Mathematical form:
  The ideal [8,4,4] code is RM(1,3), the affine functions

      c(x) = a0 + a1 x1 + a2 x2 + a3 x3,  x in F2^3.

  The Q3 edge syndrome is the coboundary

      delta_e(c) = c(x) + c(x + e_k) = a_k.

  Thus the 12-edge horizon ledger records the three gradient bits a1,a2,a3
  redundantly over the four parallel edges in each direction.  It is blind to
  the affine intercept a0, i.e. the global complement.  Therefore:

      V_code |a1,a2,a3>_B |a0>_vac = |delta(a)>_syndrome |a0>_latch

  is an isometry, while the syndrome-only map is not.

  The same coboundary construction extends to the full 8-bit Q3 cell:

      V_cell |s modulo ALL>_B |gamma>_vac
          = |delta(s)>_syndrome |gamma>_latch,

  with 128 complement classes times one vacuum/complement qubit.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import product


ROOT = __file__.rsplit("/python_code/", 1)[0]
ANCHOR = open(ROOT + "/ANCHOR.md", encoding="utf-8").read()

N_VERTS = 8
ALL = (1 << N_VERTS) - 1
EDGES = [(a, a ^ (1 << k), k) for a in range(8) for k in range(3) if a < (a ^ (1 << k))]


def check(cond: bool, msg: str) -> None:
    if not cond:
        raise AssertionError(msg)


def bit(n: int, i: int) -> int:
    return (n >> i) & 1


def syndrome_bits(s: int) -> tuple[int, ...]:
    return tuple(bit(s, i) ^ bit(s, j) for i, j, _axis in EDGES)


def syndrome_word(s: int) -> int:
    out = 0
    for n, b in enumerate(syndrome_bits(s)):
        out |= b << n
    return out


def strain(s: int) -> int:
    return sum(syndrome_bits(s))


def affine_codeword(a0: int, a1: int, a2: int, a3: int) -> int:
    """RM(1,3) codeword, with vertex labels x=(x1,x2,x3)=bits of address."""
    out = 0
    for x in range(8):
        val = a0 ^ (a1 & bit(x, 0)) ^ (a2 & bit(x, 1)) ^ (a3 & bit(x, 2))
        out |= val << x
    return out


def code_basis_index(a1: int, a2: int, a3: int, a0: int) -> int:
    """Domain column index for H_B(gradient) tensor H_vac(intercept)."""
    return (a1 << 0) | (a2 << 1) | (a3 << 2) | (a0 << 3)


def record_row(syn: int, latch: int, extra: int = 0, extra_bits: int = 0) -> int:
    """Radiation basis row: 12 syndrome bits + latch + optional carried bits."""
    return syn | (latch << 12) | (extra << 13)


def gf2_rank(words: list[int], nbits: int) -> int:
    rows = words[:]
    rank = 0
    for b in range(nbits - 1, -1, -1):
        pivot = next((i for i, row in enumerate(rows) if (row >> b) & 1), None)
        if pivot is None:
            continue
        row = rows.pop(pivot)
        rows = [r ^ row if (r >> b) & 1 else r for r in rows]
        rank += 1
    return rank


def complement_representatives() -> list[int]:
    return [s for s in range(256) if s < (s ^ ALL)]


def verify_unique_columns(rows: list[int], label: str) -> None:
    counts = Counter(rows)
    collisions = [r for r, c in counts.items() if c > 1]
    print(f"    {label}: columns={len(rows)}, image={len(counts)}, collisions={len(collisions)}")
    check(not collisions, f"{label} is not injective")


def main() -> None:
    print("[0] Canon anchors for the section-11.5 target")
    anchors = [
        ("11.5 record channel", "12-edge $Q_3$ stabilizer operators"),
        ("11.5 QND extraction", "QND) syndrome extraction"),
        ("11.5 target V", r"\mathcal{H}_B \otimes \mathcal{H}_{\mathrm{vacuum}} \to \mathcal{H}_R"),
        ("step-0 kernel theorem", r"\ker\delta=\{\emptyset,{\rm ALL}\}"),
        ("item-10 KMS ladder", r"g_{\mathcal Q}(F)e^{-\beta F}"),
    ]
    for label, text in anchors:
        ok = text in ANCHOR
        print(f"    [{'PASS' if ok else 'FAIL'}] {label}")
        check(ok, label)

    print("\n[1] Q3 coboundary / strain map")
    print(f"    vertices={N_VERTS}, edges={len(EDGES)}")
    incidence_rows = [(1 << i) | (1 << j) for i, j, _axis in EDGES]
    rank = gf2_rank(incidence_rows, N_VERTS)
    kernel = [s for s in range(256) if syndrome_word(s) == 0]
    image = {syndrome_word(s) for s in range(256)}
    print(f"    rank(delta)={rank}; |ker(delta)|={len(kernel)}; |im(delta)|={len(image)}")
    print(f"    kernel={list(map(lambda x: format(x, '08b'), kernel))}")
    check(rank == 7, "connected Q3 incidence rank should be 7")
    check(kernel == [0, ALL], "kernel should be exactly {0, ALL}")
    check(len(image) == 128, "image should have 128 cut syndromes")

    print("\n[2] Ideal [8,4,4] code from affine F2^3 generators")
    code = [affine_codeword(*coeffs) for coeffs in product([0, 1], repeat=4)]
    weights = Counter(c.bit_count() for c in code)
    print(f"    code size={len(set(code))}; weight enumerator={dict(sorted(weights.items()))}")
    check(len(set(code)) == 16, "RM(1,3) should have 16 words")
    check(dict(sorted(weights.items())) == {0: 1, 4: 14, 8: 1}, "wrong [8,4,4] weight enumerator")
    nonzero_dists = sorted((x ^ y).bit_count() for x in code for y in code if x != y)
    check(nonzero_dists[0] == 4, "minimum distance should be 4")

    gradients_to_syndrome: dict[tuple[int, int, int], int] = {}
    for a1, a2, a3 in product([0, 1], repeat=3):
        syns = set()
        for a0 in [0, 1]:
            s = affine_codeword(a0, a1, a2, a3)
            bits = syndrome_bits(s)
            by_axis = defaultdict(list)
            for b, (_i, _j, axis) in zip(bits, EDGES):
                by_axis[axis].append(b)
            check(all(len(set(vals)) == 1 for vals in by_axis.values()),
                  "affine-code syndrome should be constant on parallel edges")
            check([by_axis[k][0] for k in range(3)] == [a1, a2, a3],
                  "edge syndrome should equal affine gradient")
            syns.add(syndrome_word(s))
        check(len(syns) == 1, "intercept/global complement must be syndrome-blind")
        gradients_to_syndrome[(a1, a2, a3)] = next(iter(syns))
    print("    delta(c) records exactly the affine gradient bits a1,a2,a3;")
    print("    the affine intercept a0 is exactly the vacuum/global-complement bit.")

    print("\n[3] Explicit V_code matrix as a partial permutation")
    # V_code maps 16 basis states |gradient>|intercept> to 13-bit radiation
    # rows |12-edge syndrome>|intercept latch>.
    v_code_rows: list[int] = []
    for a1, a2, a3, a0 in product([0, 1], repeat=4):
        s = affine_codeword(a0, a1, a2, a3)
        col = code_basis_index(a1, a2, a3, a0)
        row = record_row(syndrome_word(s), a0)
        v_code_rows.append(row)
        # The formula itself is the explicit matrix: V[row, col] = 1.
        check(col < 16 and row < (1 << 13), "basis index out of range")
    verify_unique_columns(v_code_rows, "V_code: H_gradient(8) x H_vac(2) -> H_R(2^13)")
    no_vac_rows = [syndrome_word(affine_codeword(a0, a1, a2, a3))
                   for a1, a2, a3, a0 in product([0, 1], repeat=4)]
    print(f"    syndrome-only control: columns=16, image={len(set(no_vac_rows))}")
    check(len(set(no_vac_rows)) == 8, "without vacuum latch complements collapse")

    print("\n[4] Full 8-bit cell extension")
    reps = complement_representatives()
    check(len(reps) == 128, "there should be 128 complement classes")
    v_cell_rows: list[int] = []
    for class_id, rep in enumerate(reps):
        for gamma in [0, 1]:
            s = rep ^ (ALL if gamma else 0)
            row = record_row(syndrome_word(s), gamma)
            v_cell_rows.append(row)
    verify_unique_columns(v_cell_rows, "V_cell: H_(Q3/ALL)(128) x H_vac(2) -> H_R(2^13)")
    full_no_vac = [syndrome_word(s) for s in range(256)]
    duplicates = sum(c - 1 for c in Counter(full_no_vac).values() if c > 1)
    print(f"    full syndrome-only control: columns=256, image={len(set(full_no_vac))}, duplicates={duplicates}")
    check(len(set(full_no_vac)) == 128 and duplicates == 128,
          "full syndrome-only map should identify every complement pair")

    print("\n[5] Schwarzschild/KMS compatibility gate")
    # The item-10 thermal factor uses F.  Since F is the Hamming weight of the
    # emitted syndrome word, the background Boltzmann weight is readable on R.
    for s in range(256):
        check(strain(s) == syndrome_word(s).bit_count(),
              "strain should be the syndrome Hamming weight")
    code_strain_table = Counter()
    for a0, a1, a2, a3 in product([0, 1], repeat=4):
        s = affine_codeword(a0, a1, a2, a3)
        code_strain_table[strain(s)] += 1
    print(f"    ideal-code strain table from V_code records: {dict(sorted(code_strain_table.items()))}")
    print("    F(s)=|delta(s)| for every full-cell state, so any Schwarzschild/KMS")
    print("    rate depending on F is diagonal in the radiation syndrome record.")

    print("""
[VERDICT]
    The finite algebraic part of section 11.5 is constructed explicitly.
    V_code is a 16-column partial-permutation isometry from the ideal [8,4,4]
    affine content into the 12-edge horizon syndrome ledger plus the one
    vacuum/intercept latch.  V_cell extends the same construction to the full
    8-bit Q3 register as 128 complement classes times the same vacuum latch.

    The no-vacuum controls fail exactly as expected: the 12-edge syndrome
    collapses global-complement pairs.  Thus H_vacuum is not optional; it is
    the single bit required by ker(delta)={0,ALL}.

    The Schwarzschild/KMS strain weighting is compatible with this V because
    the Hawking strain F is a function of the emitted syndrome record.  What
    remains open is not the finite stabilizer isometry, but the geometric
    steady-state problem: deriving the localized-mass horizon scheduler/radial
    shell channel V_Sch(M) that composes these cell maps across an actual
    Schwarzschild background.
exit 0""")


if __name__ == "__main__":
    main()
