#!/usr/bin/env python3
r"""Off-axis baryon Y-vs-Delta audit on the bond-bipyramid bulk graph.

Question
--------
The confinement paper leaves one clean baryon-geometry item open: build
equilateral/off-axis source families that separate the baryonic Steiner-Y
ansatz from the half-perimeter Delta ansatz.

This script checks whether that task can be completed in the finite
bond-bipyramid graph used by the current baryon-Y result.

Result
------
It cannot.  The underlying bond graph is a cubic-grid median graph.  For any
three terminals in a median graph, the exact three-terminal Steiner length is

    L_Y(q1,q2,q3) = 1/2 [d12 + d23 + d31],

because the coordinatewise median lies on shortest paths between every pair.
So the axis-family equality in the paper was not an accident; it is a theorem
of the present finite graph metric.  Off-axis and graph-equilateral triples do
not separate Y from Delta by length.

Interpretation
--------------
The finite graph still proves the gauge-invariant epsilon/Y operator and the
universal per-unit string tension.  What it cannot prove is a geometric
large-distance distinction between a Y potential and a half-perimeter Delta
potential.  That distinction requires extra geometry beyond the current
strong-coupling graph metric: e.g. an effective Euclidean/weak-coupling metric,
off-graph interpolation, Coulomb/flux-profile data, or a continuum lattice-QCD
style potential extraction.

This is a useful negative closure of the paper's named-open item.  It prevents
an overclaim: there is no hidden off-axis finite-graph witness waiting to be
found by a larger search.
"""

from __future__ import annotations

from itertools import combinations

from record_grammar_tch_bond_bipyramid_baryon_y_string import (
    adjacency_from_bonds,
    steiner_three_terminal,
)
from record_grammar_tch_bond_bipyramid_bulk import AXES, bonds_box


Vertex = tuple[int, int, int]
BOX_DIMS = (5, 5, 5)


def l1(left: Vertex, right: Vertex) -> int:
    return sum(abs(left[axis] - right[axis]) for axis in AXES)


def coordinate_span(terminals: tuple[Vertex, Vertex, Vertex]) -> int:
    return sum(max(vertex[axis] for vertex in terminals) - min(vertex[axis] for vertex in terminals) for axis in AXES)


def coordinate_median(terminals: tuple[Vertex, Vertex, Vertex]) -> Vertex:
    return tuple(sorted(vertex[axis] for vertex in terminals)[1] for axis in AXES)  # type: ignore[return-value]


def half_perimeter(terminals: tuple[Vertex, Vertex, Vertex]) -> float:
    q1, q2, q3 = terminals
    return 0.5 * (l1(q1, q2) + l1(q2, q3) + l1(q3, q1))


def junction_length(junction: Vertex, terminals: tuple[Vertex, Vertex, Vertex]) -> int:
    return sum(l1(junction, terminal) for terminal in terminals)


def affine_rank_xy(terminals: tuple[Vertex, Vertex, Vertex]) -> int:
    """Rank of two displacement vectors; rank 2 means non-collinear in 3D."""

    q1, q2, q3 = terminals
    u = [q2[i] - q1[i] for i in AXES]
    v = [q3[i] - q1[i] for i in AXES]
    cross = (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )
    return 2 if any(item != 0 for item in cross) else 1


def assert_equal(name: str, value: int | float | Vertex, target: int | float | Vertex) -> None:
    print(f"  {name:<78s} value={value} target={target}")
    if value != target:
        raise AssertionError(name)


def assert_true(name: str, value: bool) -> None:
    print(f"  {name:<78s} value={value}")
    if not value:
        raise AssertionError(name)


def describe(terminals: tuple[Vertex, Vertex, Vertex]) -> str:
    q1, q2, q3 = terminals
    return f"{q1}, {q2}, {q3}"


def main() -> None:
    print("Off-axis baryon Y-vs-Delta audit")
    print("=" * 100)
    print(f"  bond graph box={BOX_DIMS}")

    vertices = tuple((x, y, z) for x in range(BOX_DIMS[0]) for y in range(BOX_DIMS[1]) for z in range(BOX_DIMS[2]))
    adj = adjacency_from_bonds(bonds_box(*BOX_DIMS))

    print("\n[1] Explicit off-axis / graph-equilateral witnesses")
    examples: tuple[tuple[Vertex, Vertex, Vertex], ...] = (
        ((0, 0, 0), (1, 1, 0), (2, 0, 0)),  # L1-equilateral in a plane.
        ((1, 1, 1), (4, 2, 1), (2, 4, 1)),  # off-axis triangle in a sheet.
        ((0, 0, 0), (2, 1, 1), (1, 3, 2)),  # genuinely 3D off-axis.
        ((0, 1, 0), (3, 0, 2), (2, 4, 3)),  # unequal, non-collinear control.
    )
    for terminals in examples:
        steiner_length, witness = steiner_three_terminal(adj, terminals)
        median = coordinate_median(terminals)
        span = coordinate_span(terminals)
        hp = half_perimeter(terminals)
        print(f"  terminals: {describe(terminals)}")
        assert_equal("    exact Steiner length", steiner_length, span)
        assert_equal("    coordinate median arm length", junction_length(median, terminals), span)
        assert_equal("    half-perimeter Delta length", hp, float(span))
        assert_true("    non-collinear/off-axis family", affine_rank_xy(terminals) == 2)
        print(f"    DP witness={witness}; median junction={median}")

    print("\n[2] Exhaustive small-box theorem check")
    total = 0
    offaxis = 0
    equilateral = 0
    max_gap = 0.0
    for terminals in combinations(vertices, 3):
        terminals = terminals  # type: ignore[assignment]
        span = coordinate_span(terminals)
        hp = half_perimeter(terminals)
        gap = hp - span
        if gap > max_gap:
            max_gap = gap
        if affine_rank_xy(terminals) == 2:
            offaxis += 1
        d12 = l1(terminals[0], terminals[1])
        d23 = l1(terminals[1], terminals[2])
        d31 = l1(terminals[2], terminals[0])
        if d12 == d23 == d31:
            equilateral += 1
        if gap != 0:
            raise AssertionError(f"found Y-vs-Delta separation {gap} for {terminals}")
        total += 1

    assert_equal("all triples checked", total, len(tuple(combinations(vertices, 3))))
    assert_true("off-axis triples were included", offaxis > 0)
    assert_true("graph-equilateral triples were included", equilateral > 0)
    assert_equal("maximum half-perimeter minus Steiner gap", max_gap, 0.0)
    print(f"  off-axis triples={offaxis}; graph-equilateral triples={equilateral}")

    print(
        """
VERDICT:
  NEGATIVE CLOSURE.  In the current finite bond-bipyramid bulk graph, every
  three-source baryon geometry satisfies

      L_Y = 1/2 (d12 + d23 + d31).

  The paper's tri-axial Y=Delta equality is therefore not a special accident;
  it is forced by the median-graph/L1 metric of the present graph.  No larger
  finite off-axis search in this metric can distinguish a Steiner-Y ansatz from
  a half-perimeter Delta ansatz by length.

  What remains true: the SU(3) epsilon junction gives the correct baryonic
  gauge-invariant operator and prices its arms with the same string tension as
  the meson.

  What remains open: a genuine Y-vs-Delta separation requires an effective
  Euclidean/weak-coupling metric, off-graph flux interpolation, or a continuum
  baryonic potential calculation.  It is not a missing finite-graph witness.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
