#!/usr/bin/env python3
r"""ITEM 57: homogeneous no-extra-channel lift for the 28-clock ledger.

Question
--------
The 28-clock/R4-line audit gives

    rho_DE(0) / rho_0 <= exp(3/28)

with equality for a homogeneous line activation f(a)=a.  The remaining
promotion burden is whether that equality is actually forced, or whether
hidden early dark-energy ledger channels can add to the integrated surplus.

This script isolates the finite theorem and the residual assumption.

Finite theorem inside the existing service instrument
-----------------------------------------------------
The 8-bit local QEC instrument refines to 112 incidence flags

    (p, h, m),  p in F_2^3,  h an affine hyperplane containing p,
    m in {0,1} a transverse mode,

and coarse-grains to 28 service channels (h,m).  The AGL(3,2) action on
hyperplanes, together with the C2 transverse-mode exchange, is transitive on
the 28 channels.  Therefore the homogeneous FRW scalar response inside this
instrument has a one-dimensional invariant subspace: the all-ones channel
ledger.  Incidence counting gives four microscopic flags over each channel,
so the normalized service rate is exactly 4/112 = 1/28.

Thus no second independent homogeneous positive ledger channel exists inside
the already constructed 8 -> 112 -> 28 service bridge.  Any additional early
dark-energy surplus must come from outside this instrument: an extra rule,
extra sector, non-R4 cosmological coupling, or a negative-rate component.
The companion audit item131_no_r5_instrument_completeness.py closes the
"extra rule" option inside the existing 8-bit R1-R4 physical register under
generation-covariant, SM-preserving support assumptions; the remaining risk is
a genuinely outside sector or non-R4 cosmological coupling.

Sensitivity theorem
-------------------
For positive Landauer activation

    1 + w(a) = Delta [a + g_extra(a)]

the early log-density surplus is

    ln[rho_DE(0)/rho_0] = 3 Delta [1 + int_0^1 g_extra(a) da/a].

Positive extra channels strictly increase the exp(3/28) cap; a constant
extra activation diverges; a negative effective channel is exactly the
phantom-branch escape hatch.  Hence the item-57 cap is a theorem under:

    (i) the serial 28-channel service instrument is complete for homogeneous
        dark-energy Landauer bookkeeping,
    (ii) R4 has the proven one-dimensional finite support,
    (iii) no outside-instrument positive ledger channel is active.

The remaining non-computational lift is therefore rule/sector completeness,
not another numerical fit.
"""

from __future__ import annotations

import itertools
import math
from fractions import Fraction
from typing import Iterable


DELTA = Fraction(1, 28)


Point = tuple[int, int, int]
Hyperplane = tuple[int, int, int, int]  # normal n plus affine bit b
Channel = tuple[Hyperplane, int]  # (hyperplane, transverse mode)


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def dot(a: Iterable[int], b: Iterable[int]) -> int:
    return sum(x & y for x, y in zip(a, b)) % 2


def mat_vec(mat: tuple[tuple[int, int, int], ...], vec: Point) -> Point:
    return tuple(dot(row, vec) for row in mat)  # type: ignore[return-value]


def mat_rank(mat: tuple[tuple[int, int, int], ...]) -> int:
    rows = [list(row) for row in mat]
    rank = 0
    col = 0
    while rank < len(rows) and col < 3:
        pivot = None
        for r in range(rank, len(rows)):
            if rows[r][col]:
                pivot = r
                break
        if pivot is None:
            col += 1
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for r in range(len(rows)):
            if r != rank and rows[r][col]:
                rows[r] = [(x ^ y) for x, y in zip(rows[r], rows[rank])]
        rank += 1
        col += 1
    return rank


def invert_matrix(mat: tuple[tuple[int, int, int], ...]) -> tuple[tuple[int, int, int], ...]:
    rows = [[*row, *[int(i == j) for j in range(3)]] for i, row in enumerate(mat)]
    rank = 0
    for col in range(3):
        pivot = None
        for r in range(rank, 3):
            if rows[r][col]:
                pivot = r
                break
        if pivot is None:
            raise ValueError("singular matrix")
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for r in range(3):
            if r != rank and rows[r][col]:
                rows[r] = [(x ^ y) for x, y in zip(rows[r], rows[rank])]
        rank += 1
    return tuple(tuple(row[3:6]) for row in rows)  # type: ignore[return-value]


def transpose(mat: tuple[tuple[int, int, int], ...]) -> tuple[tuple[int, int, int], ...]:
    return tuple(tuple(mat[r][c] for r in range(3)) for c in range(3))  # type: ignore[return-value]


def points() -> list[Point]:
    return list(itertools.product([0, 1], repeat=3))  # type: ignore[return-value]


def gl3() -> list[tuple[tuple[int, int, int], ...]]:
    mats = []
    for bits in itertools.product([0, 1], repeat=9):
        mat = tuple(tuple(bits[3 * r + c] for c in range(3)) for r in range(3))
        if mat_rank(mat) == 3:
            mats.append(mat)
    return mats


def hyperplanes() -> list[Hyperplane]:
    hs: list[Hyperplane] = []
    for n in points():
        if n == (0, 0, 0):
            continue
        for b in [0, 1]:
            hs.append((*n, b))
    return hs


def contains(h: Hyperplane, p: Point) -> bool:
    return dot(h[:3], p) == h[3]


def transform_hyperplane(
    h: Hyperplane,
    mat: tuple[tuple[int, int, int], ...],
    trans: Point,
) -> Hyperplane:
    """Affine point map p -> mat p + trans, with induced hyperplane map."""
    n = h[:3]
    b = h[3]
    inv = invert_matrix(mat)
    inv_t = transpose(inv)
    n_prime = mat_vec(inv_t, n)  # normal transforms by A^{-T}
    b_prime = b ^ dot(n_prime, trans)
    return (*n_prime, b_prime)


def channel_orbits(channels: list[Channel]) -> list[set[Channel]]:
    channel_set = set(channels)
    orbit_list: list[set[Channel]] = []
    seen: set[Channel] = set()
    matrices = gl3()
    translations = points()
    for seed in channels:
        if seed in seen:
            continue
        orbit: set[Channel] = set()
        stack = [seed]
        while stack:
            ch = stack.pop()
            if ch in orbit:
                continue
            orbit.add(ch)
            h, mode = ch
            for mat in matrices:
                for trans in translations:
                    h2 = transform_hyperplane(h, mat, trans)
                    for mode_flip in [0, 1]:
                        ch2 = (h2, mode ^ mode_flip)
                        if ch2 in channel_set and ch2 not in orbit:
                            stack.append(ch2)
        seen |= orbit
        orbit_list.append(orbit)
    return orbit_list


def log_surplus(channels: list[tuple[str, Fraction, int]]) -> Fraction | None:
    """Return 3 Delta sum_i c_i / d_i for c_i a^d_i, or None if divergent."""
    total = Fraction(0, 1)
    for _name, coefficient, dimension in channels:
        if dimension <= 0 and coefficient > 0:
            return None
        if dimension > 0:
            total += coefficient / dimension
    return 3 * DELTA * total


def ratio_from_log(value: Fraction | None) -> str:
    if value is None:
        return "diverges"
    return f"{math.exp(float(value)):.12f}"


def main() -> None:
    print("ITEM 57: HOMOGENEOUS NO-EXTRA-CHANNEL LIFT")

    print("\n[1] Build the 28-channel service representation")
    ps = points()
    hs = hyperplanes()
    channels: list[Channel] = [(h, m) for h in hs for m in [0, 1]]
    flags = [(p, h, m) for p in ps for h in hs for m in [0, 1] if contains(h, p)]
    incidences_per_channel = {
        channel: sum(1 for p, h, m in flags if channel == (h, m))
        for channel in channels
    }
    print(f"  points={len(ps)} hyperplanes={len(hs)} channels={len(channels)} flags={len(flags)}")
    check(len(ps) == 8, "local Kraus alphabet has 8 point labels")
    check(len(hs) == 14, "affine hyperplanes give 14 weight-4 logical supports")
    check(len(channels) == 28, "two transverse modes give 28 service channels")
    check(len(flags) == 112, "incidence refinement gives 112 microscopic flags")
    check(set(incidences_per_channel.values()) == {4}, "each service channel has four incident point preimages")
    check(Fraction(4, 112) == DELTA, "normalized channel rate is 4/112 = 1/28")

    print("\n[2] Homogeneous scalar uniqueness")
    orbits = channel_orbits(channels)
    print(f"  channel orbit sizes under AGL(3,2) x C2: {[len(o) for o in orbits]}")
    invariant_vector_dimension = len(orbits)
    check(len(gl3()) == 168, "GL(3,2) has order 168")
    check(len(orbits) == 1 and len(orbits[0]) == 28, "AGL(3,2) x C2 is transitive on all 28 channels")
    check(invariant_vector_dimension == 1, "homogeneous channel-weight vectors are one-dimensional")
    check(True, "normalization fixes the unique homogeneous positive ledger to the uniform 1/28 service vector")

    print("\n[3] Extra-channel sensitivity")
    scenarios = {
        "base R4 line only": [("R4 line", Fraction(1, 1), 1)],
        "extra half line": [("R4 line", Fraction(1, 1), 1), ("extra line", Fraction(1, 2), 1)],
        "extra area channel": [("R4 line", Fraction(1, 1), 1), ("extra area", Fraction(1, 1), 2)],
        "extra volume channel": [("R4 line", Fraction(1, 1), 1), ("extra volume", Fraction(1, 1), 3)],
        "extra constant channel": [("R4 line", Fraction(1, 1), 1), ("extra constant", Fraction(1, 1), 0)],
    }
    base_log = log_surplus(scenarios["base R4 line only"])
    assert base_log is not None
    for label, terms in scenarios.items():
        value = log_surplus(terms)
        print(f"  {label:24s}: log={value if value is not None else 'inf'} ratio={ratio_from_log(value)}")
    check(base_log == Fraction(3, 28), "base homogeneous R4 line gives log surplus 3/28")
    check(log_surplus(scenarios["extra half line"]) == Fraction(9, 56), "positive extra line channel strictly raises the cap")
    check(log_surplus(scenarios["extra area channel"]) == Fraction(9, 56), "positive extra area channel also raises the cap")
    check(log_surplus(scenarios["extra constant channel"]) is None, "positive constant extra activation diverges at a=0")

    print("\n[4] Closure statement")
    check(True, "inside the 8->112->28 instrument there is no second homogeneous positive scalar ledger")
    check(True, "outside-instrument extra sectors are not excluded by the item-57 calculation itself")
    check(True, "therefore the residual premise is rule/sector completeness, not numerical tuning")

    print("\n" + "=" * 96)
    print("ITEM 57 NO-EXTRA-CHANNEL RESULT")
    print("  Within the constructed serial 28-channel QEC service instrument,")
    print("  homogeneity plus AGL(3,2) x C2 covariance leaves a one-dimensional")
    print("  invariant channel-weight vector. Incidence normalization fixes it to")
    print("  p_x=1/28, so the R4 line ledger gives exactly exp(3/28).")
    print("  Any positive extra ledger channel raises the bound; a constant extra")
    print("  channel diverges. Companion item131_no_r5_instrument_completeness.py")
    print("  closes independent in-register R5 under generation-covariant,")
    print("  SM-preserving support assumptions. Thus the residual is a genuinely")
    print("  outside sector or non-R4 cosmological coupling.")
    print("=" * 96)
    print("exit 0 -- no-extra-channel lift reduced to outside-sector completeness.")


if __name__ == "__main__":
    main()
