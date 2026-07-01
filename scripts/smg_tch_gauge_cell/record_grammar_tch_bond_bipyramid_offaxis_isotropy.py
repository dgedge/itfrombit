#!/usr/bin/env python3
r"""Primitive off-axis Wilson-loop isotropy gate on the bond-bipyramid sheet.

Purpose
-------
``record_grammar_tch_bond_bipyramid_isotropy.py`` tested cubic-axis rotations
and in-sheet axis swaps.  That is not yet the off-axis question: the licensed
flat sheet is triangular, with three primitive macro directions

    e0 = (1, 0),    e1 = (0, 1),    e2 = (1, -1),

plus their reverses.  A substrate-scale Lorentz/rotational fingerprint would
first show up as a direction-dependent Wilson area coefficient on rhombi built
from these primitive directions.

This script constructs primitive rhombus loop families from the three unordered
pairs

    (e0, e1), (e0, e2), (e1, e2),

places them inside the finite sheet so every macro face is complete, and checks

    A_ref(R,T) = 4 R T,
    P(R,T)     = 4(R+T)-2,
    chi_area   = 4 sigma,     sigma = -log(beta/18),

with zero spread across the primitive directions.

Boundary
--------
This is a finite, leading-strong-coupling primitive-direction gate.  It does
not prove full continuum SO(3), nonprimitive diagonal-path isotropy, exact
roughening entropy, weak-coupling scaling, or Lorentz restoration.  Its result
is narrower but useful: the corrected bond-bipyramid flat sheet has no forced
bare string-tension anisotropy among its primitive off-axis directions.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from dataclasses import dataclass

from record_grammar_tch_bond_bipyramid_bulk import (
    BondComplex,
    build_bond_complex,
    bonds_box,
    canonical_edge,
    xor_masks,
)
from record_grammar_tch_bond_bipyramid_edge_loop_creutz import creutz
from record_grammar_tch_bond_bipyramid_isoperimetry import largest_connected_plane
from record_grammar_tch_bond_bipyramid_static_potential import (
    BETA,
    boundary_stats,
    face_sum,
    u_fundamental,
)


SLAB_DIMS = (10, 10, 3)
MAX_RECT = 4
MAX_CREUTZ = 3
EXPECTED_PRIMITIVES = {(1, 0), (0, 1), (1, -1)}


Vec2 = tuple[int, int]


@dataclass(frozen=True)
class Family:
    name: str
    u: Vec2
    v: Vec2


@dataclass(frozen=True)
class OffAxisRow:
    family: str
    r: int
    t: int
    area: int
    perimeter: int
    y_area: float
    y_with_perimeter: float


def assert_equal(name: str, value: int | set[Vec2], target: int | set[Vec2]) -> None:
    print(f"  {name:<78s} value={value} target={target}")
    if value != target:
        raise AssertionError(name)


def assert_true(name: str, value: bool) -> None:
    print(f"  {name:<78s} value={value}")
    if not value:
        raise AssertionError(name)


def assert_close(name: str, value: float, target: float, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<78s} value={value:.12g} target={target:.12g} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def face_edges(face: tuple[tuple[int, int, int], ...]) -> set[tuple[tuple[int, int, int], tuple[int, int, int]]]:
    return {
        canonical_edge(face[left], face[right])
        for left, right in ((0, 1), (1, 2), (0, 2))
    }


def macro_sheet(
    complex_: BondComplex,
) -> tuple[dict[Vec2, tuple[int, ...]], set[Vec2], set[Vec2]]:
    """Return complete macro faces and primitive neighbour offsets on one flat sheet."""

    graph = largest_connected_plane(complex_)
    in_axes = tuple(axis for axis in (0, 1, 2) if axis != graph.axis)
    by_macro: dict[Vec2, list[int]] = defaultdict(list)
    face_macro: dict[int, Vec2] = {}
    for face_i in graph.nodes:
        sums = face_sum(complex_.faces[face_i])
        macro = (sums[in_axes[0]] // 6, sums[in_axes[1]] // 6)
        by_macro[macro].append(face_i)
        face_macro[face_i] = macro

    complete = {macro for macro, faces in by_macro.items() if len(faces) == 4}
    complete_by_macro = {
        macro: tuple(sorted(faces))
        for macro, faces in by_macro.items()
        if macro in complete
    }

    offsets: Counter[Vec2] = Counter()
    for face_i, neighbours in graph.adj.items():
        for face_j in neighbours:
            if face_i >= face_j:
                continue
            left = face_macro[face_i]
            right = face_macro[face_j]
            delta = (right[0] - left[0], right[1] - left[1])
            if delta == (0, 0):
                continue
            # Canonicalise to the positive representative used by the existing
            # sheet coordinates.
            if delta[0] < 0 or (delta[0] == 0 and delta[1] < 0):
                delta = (-delta[0], -delta[1])
            offsets[delta] += 1

    primitives = set(offsets)
    return complete_by_macro, complete, primitives


def parallelogram_points(u: Vec2, v: Vec2, r: int, t: int) -> tuple[Vec2, ...]:
    return tuple(
        (i * u[0] + j * v[0], i * u[1] + j * v[1])
        for i in range(r)
        for j in range(t)
    )


def choose_offset(points: tuple[Vec2, ...], complete: set[Vec2]) -> Vec2:
    min_x = min(x for x, _ in complete)
    max_x = max(x for x, _ in complete)
    min_y = min(y for _, y in complete)
    max_y = max(y for _, y in complete)
    candidates: list[tuple[int, int, Vec2]] = []
    for ox in range(min_x - min(x for x, _ in points), max_x - max(x for x, _ in points) + 1):
        for oy in range(min_y - min(y for _, y in points), max_y - max(y for _, y in points) + 1):
            shifted = tuple((x + ox, y + oy) for x, y in points)
            if not all(item in complete for item in shifted):
                continue
            margin = min(
                min(x - min_x, max_x - x, y - min_y, max_y - y)
                for x, y in shifted
            )
            # Secondary tie-breaker keeps the patch near the centre.
            centre_penalty = abs((min_x + max_x) - 2 * ox) + abs((min_y + max_y) - 2 * oy)
            candidates.append((margin, -centre_penalty, (ox, oy)))
    if not candidates:
        raise AssertionError(f"no complete interior offset for points={points}")
    return max(candidates)[2]


def patch_faces(
    by_macro: dict[Vec2, tuple[int, ...]],
    complete: set[Vec2],
    family: Family,
    r: int,
    t: int,
) -> tuple[int, ...]:
    points = parallelogram_points(family.u, family.v, r, t)
    offset = choose_offset(points, complete)
    faces: list[int] = []
    for x, y in points:
        faces.extend(by_macro[(x + offset[0], y + offset[1])])
    return tuple(sorted(faces))


def determinant(left: Vec2, right: Vec2) -> int:
    return left[0] * right[1] - left[1] * right[0]


def rows_for_family(complex_: BondComplex, by_macro: dict[Vec2, tuple[int, ...]], complete: set[Vec2], family: Family) -> list[OffAxisRow]:
    sigma = -math.log(u_fundamental(BETA))
    mu_control = 0.37  # arbitrary perimeter term used only to verify Creutz cancellation.
    rows: list[OffAxisRow] = []
    det = abs(determinant(family.u, family.v))
    assert_equal(f"{family.name} primitive rhombus determinant", det, 1)
    for r in range(1, MAX_RECT + 1):
        for t in range(1, MAX_RECT + 1):
            patch = patch_faces(by_macro, complete, family, r, t)
            boundary = xor_masks(complex_.face_masks, patch)
            components, even, max_degree, perimeter = boundary_stats(boundary, complex_.edges)
            expected_area = 4 * r * t
            expected_perimeter = 4 * (r + t) - 2
            assert_equal(f"{family.name} R={r},T={t} area", len(patch), expected_area)
            assert_equal(f"{family.name} R={r},T={t} perimeter", perimeter, expected_perimeter)
            assert_true(f"{family.name} R={r},T={t} boundary is even", even)
            assert_equal(f"{family.name} R={r},T={t} boundary components", components, 1)
            assert_equal(f"{family.name} R={r},T={t} boundary max degree", max_degree, 2)
            y_area = sigma * expected_area
            rows.append(
                OffAxisRow(
                    family=family.name,
                    r=r,
                    t=t,
                    area=expected_area,
                    perimeter=perimeter,
                    y_area=y_area,
                    y_with_perimeter=y_area + mu_control * perimeter + 1.23,
                )
            )
    return rows


def main() -> None:
    print("Primitive off-axis Wilson-loop isotropy gate")
    print("=" * 100)
    sigma = -math.log(u_fundamental(BETA))
    expected_chi = 4.0 * sigma
    print(f"  slab dimensions={SLAB_DIMS}; beta={BETA:.3f}; sigma=-ln(beta/18)={sigma:.9f}")
    print(f"  primitive-rhombus table R,T<= {MAX_RECT}; Creutz uses R,T<= {MAX_CREUTZ}")

    complex_ = build_bond_complex("bond-offaxis-isotropy", bonds_box(*SLAB_DIMS))
    by_macro, complete, primitives = macro_sheet(complex_)
    assert_equal("sheet primitive macro-neighbour offsets", primitives, EXPECTED_PRIMITIVES)

    families = (
        Family("e0_e1_axis_rhombus", (1, 0), (0, 1)),
        Family("e0_e2_offaxis_rhombus", (1, 0), (1, -1)),
        Family("e1_e2_offaxis_rhombus", (0, 1), (1, -1)),
    )

    print("\n[1] Primitive rhombus Creutz readings")
    print("  family                  chi_area      chi_with_perimeter")
    chi_values: list[float] = []
    chi_perim_values: list[float] = []
    for family in families:
        rows = rows_for_family(complex_, by_macro, complete, family)
        table = {(row.r, row.t): row for row in rows}
        for r in range(1, MAX_CREUTZ + 1):
            for t in range(1, MAX_CREUTZ + 1):
                chi_area = creutz(table, "y_area", r, t)
                chi_perim = creutz(table, "y_with_perimeter", r, t)
                assert_close(f"{family.name} R={r},T={t} bare Creutz", chi_area, expected_chi)
                assert_close(f"{family.name} R={r},T={t} perimeter-cancelled Creutz", chi_perim, expected_chi)
                chi_values.append(chi_area)
                chi_perim_values.append(chi_perim)
        print(f"  {family.name:<24s} {expected_chi:>12.9f} {expected_chi:>20.9f}")

    print("\n[2] Off-axis spread")
    area_spread = max(chi_values) - min(chi_values)
    perim_spread = max(chi_perim_values) - min(chi_perim_values)
    assert_close("primitive-direction bare string-tension spread", area_spread, 0.0, tol=1e-12)
    assert_close("primitive-direction perimeter-cancelled spread", perim_spread, 0.0, tol=1e-12)

    print(
        """
VERDICT:
  PASS.  The licensed bond-bipyramid sheet is a triangular macro-sheet with
  primitive neighbour directions (1,0), (0,1), and (1,-1).  Wilson rhombi built
  from all three primitive direction pairs have the same finite identities:

      A_ref = 4RT,    P = 4(R+T)-2,    chi_area = 4 sigma.

  An arbitrary additive perimeter term cancels in the same Creutz combination,
  and the residual primitive-direction spread is zero to numerical precision.
  Therefore the leading strong-coupling bare Wilson string tension has no
  forced anisotropy among primitive off-axis directions on the licensed sheet.

  Scope: this is not yet continuum SO(3).  Nonprimitive graph directions,
  exact roughening entropy off axis, weak-coupling scaling, and physical
  Lorentz-violation bounds remain separate tests.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
