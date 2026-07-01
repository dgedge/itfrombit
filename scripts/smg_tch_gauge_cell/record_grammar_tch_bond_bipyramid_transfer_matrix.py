#!/usr/bin/env python3
r"""Exact fixed-width transfer count for bond-bipyramid Wilson loops.

Purpose
-------
``record_grammar_tch_bond_bipyramid_edge_loop_creutz.py`` found a sector lower
bound

    N_sector(R,T) = 2^((R-1)(T-1))

for simple edge loops on the licensed flat offset sheet.  That lower bound was
already enough to refute an O(1) bond-dimension reading, but it left a
sector-grade caveat: the exact minimum-surface count ``N_min`` had not been
computed.

This script computes the exact minimum count in the fixed-width transfer strip
spanned by the loop.  It writes the problem as a binary record-action energy:
each bond-bipyramid cell in the strip is a 0/1 flip variable, and each face
contributes one unit when the reference patch parity disagrees with the cell
flip parity.  Contracting that local tensor network by min-sum/count
elimination gives the exact minimum area and exact degeneracy for the strip.

The result is

    N_min^strip(R,T) = 2^((R-1)(T-1)),

so the lower-bound sector is not hiding a smaller fixed-strip theorem: it is
the exact fixed-width transfer count.  Equivalently, at fixed width R the
transfer multiplier per added time column is 2^(R-1).  The ledger is therefore
2D-tensor-network tractable, not fixed-bond-dimension: it is polynomial along a
fixed-width transfer direction and exponential in loop width.

Boundary
--------
The script is exact for the tight fixed-width strip defined by the loop's
coordinate rectangle in the finite 8 x 8 x 3 bond-bipyramid slab.  It does not
count wider-margin or full-slab roughening surfaces.  Those belong to the
transverse-thickness / roughening-entropy frontier.  The point here is narrower:
remove the sector-grade caveat for the fixed-width transfer observable and
measure the exact width law.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass

from record_grammar_tch_bond_bipyramid_bulk import (
    BondComplex,
    build_bond_complex,
    bonds_box,
)
from record_grammar_tch_bond_bipyramid_edge_loop_creutz import (
    SLAB_DIMS,
    plane_face_coordinates,
    rectangle_patch,
)
from record_grammar_tch_bond_bipyramid_static_potential import (
    face_vertex_bbox,
    cell_vertex_bbox,
)


TRANSFER_MAX_RECT = 4


@dataclass(frozen=True)
class Factor:
    scope: tuple[int, ...]
    table: dict[int, tuple[int, int]]


@dataclass(frozen=True)
class TransferRow:
    r: int
    t: int
    area: int
    exact_count: int
    log_count: float
    promotion_bits: int
    transfer_bond_dim: int
    variables: int
    terms: int
    max_bag: int


def assert_equal(name: str, value: int, target: int) -> None:
    print(f"  {name:<78s} value={value} target={target}")
    if value != target:
        raise AssertionError(name)


def assert_close(name: str, value: float, target: float, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<78s} value={value:.12g} target={target:.12g} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def face_to_cells(complex_: BondComplex) -> dict[int, tuple[int, ...]]:
    out: dict[int, list[int]] = defaultdict(list)
    for cell_i, mask in enumerate(complex_.cell_masks):
        bit = 0
        work = mask
        while work:
            if work & 1:
                out[bit].append(cell_i)
            work >>= 1
            bit += 1
    return {face_i: tuple(cells) for face_i, cells in out.items()}


def tight_strip_cells(complex_: BondComplex, patch: tuple[int, ...]) -> tuple[int, ...]:
    """Cells whose bounding boxes intersect the loop rectangle bounding box."""

    patch_lo, patch_hi = face_vertex_bbox(complex_, patch)
    selected: list[int] = []
    for cell_i, cell_mask in enumerate(complex_.cell_masks):
        cell_lo, cell_hi = cell_vertex_bbox(complex_, cell_mask)
        if all(cell_hi[axis] >= patch_lo[axis] and cell_lo[axis] <= patch_hi[axis] for axis in range(3)):
            selected.append(cell_i)
    return tuple(selected)


def parity_factor(scope: tuple[int, ...], start_bit: int) -> Factor:
    table: dict[int, tuple[int, int]] = {}
    for assignment in range(1 << len(scope)):
        table[assignment] = (start_bit ^ (assignment.bit_count() & 1), 1)
    return Factor(scope=tuple(sorted(scope)), table=table)


def combine_factors(factors: list[Factor]) -> Factor:
    if not factors:
        return Factor(scope=(), table={0: (0, 1)})
    union_scope = tuple(sorted(set().union(*(factor.scope for factor in factors))))
    union_positions = {var: i for i, var in enumerate(union_scope)}
    projections = [
        tuple(union_positions[var] for var in factor.scope)
        for factor in factors
    ]
    table: dict[int, tuple[int, int]] = {}
    for assignment in range(1 << len(union_scope)):
        energy = 0
        count = 1
        for factor, positions in zip(factors, projections):
            factor_assignment = 0
            for j, pos in enumerate(positions):
                if assignment & (1 << pos):
                    factor_assignment |= 1 << j
            factor_energy, factor_count = factor.table[factor_assignment]
            energy += factor_energy
            count *= factor_count
        table[assignment] = (energy, count)
    return Factor(scope=union_scope, table=table)


def min_fill_order(n_variables: int, terms: list[tuple[tuple[int, ...], int]]) -> tuple[tuple[int, ...], int]:
    graph: dict[int, set[int]] = {var: set() for var in range(n_variables)}
    for scope, _ in terms:
        for left in scope:
            for right in scope:
                if left < right:
                    graph[left].add(right)
                    graph[right].add(left)

    remaining = set(range(n_variables))
    order: list[int] = []
    max_bag = 0
    while remaining:
        def score(var: int) -> tuple[int, int, int]:
            neighbours = [item for item in graph[var] if item in remaining]
            fill = 0
            for i, left in enumerate(neighbours):
                for right in neighbours[i + 1 :]:
                    if right not in graph[left]:
                        fill += 1
            return fill, len(neighbours), var

        var = min(remaining, key=score)
        neighbours = [item for item in graph[var] if item in remaining]
        max_bag = max(max_bag, len(neighbours) + 1)
        for i, left in enumerate(neighbours):
            for right in neighbours[i + 1 :]:
                graph[left].add(right)
                graph[right].add(left)
        remaining.remove(var)
        order.append(var)
    return tuple(order), max_bag


def exact_min_count(n_variables: int, terms: list[tuple[tuple[int, ...], int]]) -> tuple[int, int, int]:
    """Min-sum/count contraction of the fixed-width strip tensor network."""

    factors = [parity_factor(scope, start_bit) for scope, start_bit in terms]
    order, max_bag = min_fill_order(n_variables, terms)

    for var in order:
        involved = [factor for factor in factors if var in factor.scope]
        if not involved:
            continue
        factors = [factor for factor in factors if var not in factor.scope]
        combined = combine_factors(involved)
        rest_scope = tuple(item for item in combined.scope if item != var)
        rest_positions = [combined.scope.index(item) for item in rest_scope]
        new_table: dict[int, tuple[int, int]] = {}
        for assignment, (energy, count) in combined.table.items():
            rest_assignment = 0
            for j, pos in enumerate(rest_positions):
                if assignment & (1 << pos):
                    rest_assignment |= 1 << j
            current = new_table.get(rest_assignment)
            if current is None or energy < current[0]:
                new_table[rest_assignment] = (energy, count)
            elif energy == current[0]:
                new_table[rest_assignment] = (energy, current[1] + count)
        factors.append(Factor(scope=rest_scope, table=new_table))

    combined = combine_factors(factors)
    if combined.scope:
        raise AssertionError(f"uneliminated variables remain: {combined.scope}")
    min_energy, degeneracy = combined.table[0]
    return min_energy, degeneracy, max_bag


def transfer_terms(
    complex_: BondComplex,
    incident: dict[int, tuple[int, ...]],
    patch: tuple[int, ...],
) -> tuple[int, list[tuple[tuple[int, ...], int]], int]:
    selected_cells = tight_strip_cells(complex_, patch)
    cell_map = {cell_i: dense_i for dense_i, cell_i in enumerate(selected_cells)}

    active_faces = set(patch)
    for cell_i in selected_cells:
        cell_mask = complex_.cell_masks[cell_i]
        bit = 0
        work = cell_mask
        while work:
            if work & 1:
                active_faces.add(bit)
            work >>= 1
            bit += 1

    patch_set = set(patch)
    terms: list[tuple[tuple[int, ...], int]] = []
    constant = 0
    for face_i in active_faces:
        scope = tuple(cell_map[cell_i] for cell_i in incident[face_i] if cell_i in cell_map)
        start_bit = 1 if face_i in patch_set else 0
        if scope:
            terms.append((tuple(sorted(scope)), start_bit))
        else:
            constant += start_bit
    return len(selected_cells), terms, constant


def measure_rows() -> list[TransferRow]:
    complex_ = build_bond_complex("bond-slab-8x8x3", bonds_box(*SLAB_DIMS))
    _, _, coords = plane_face_coordinates(complex_)
    incident = face_to_cells(complex_)
    rows: list[TransferRow] = []
    for r in range(1, TRANSFER_MAX_RECT + 1):
        for t in range(1, TRANSFER_MAX_RECT + 1):
            patch = rectangle_patch(coords, r, t)
            n_variables, terms, constant = transfer_terms(complex_, incident, patch)
            min_energy, exact_count, max_bag = exact_min_count(n_variables, terms)
            min_energy += constant
            expected_area = 4 * r * t
            expected_count = 2 ** ((r - 1) * (t - 1))
            assert_equal(f"R={r},T={t} exact strip min area", min_energy, expected_area)
            assert_equal(f"R={r},T={t} exact strip N_min", exact_count, expected_count)
            rows.append(
                TransferRow(
                    r=r,
                    t=t,
                    area=expected_area,
                    exact_count=exact_count,
                    log_count=math.log(exact_count),
                    promotion_bits=(r - 1) * (t - 1),
                    transfer_bond_dim=2 ** max(0, r - 1),
                    variables=n_variables,
                    terms=len(terms),
                    max_bag=max_bag,
                )
            )
    return rows


def main() -> None:
    print("Exact fixed-width transfer count for bond-bipyramid Wilson loops")
    print("=" * 100)
    print(f"  slab dimensions={SLAB_DIMS}; exact rectangle table R,T<= {TRANSFER_MAX_RECT}")
    rows = measure_rows()

    print("\n[1] Exact fixed-strip transfer table")
    print("  R T  A_ref  N_min_exact  logN      promo_bits  D_width  vars terms max_bag")
    for row in rows:
        print(
            f"  {row.r:>1d} {row.t:>1d}"
            f" {row.area:>6d}"
            f" {row.exact_count:>12d}"
            f" {row.log_count:>8.3f}"
            f" {row.promotion_bits:>11d}"
            f" {row.transfer_bond_dim:>8d}"
            f" {row.variables:>5d}"
            f" {row.terms:>5d}"
            f" {row.max_bag:>7d}"
        )

    print("\n[2] Width-law gates")
    by_rt = {(row.r, row.t): row for row in rows}
    for r in range(1, TRANSFER_MAX_RECT + 1):
        multiplier = 2 ** max(0, r - 1)
        for t in range(1, TRANSFER_MAX_RECT):
            ratio = by_rt[(r, t + 1)].exact_count // by_rt[(r, t)].exact_count
            assert_equal(f"fixed R={r}, transfer multiplier T={t}->{t+1}", ratio, multiplier)

    chi_entropy = (
        by_rt[(1, 1)].log_count
        + by_rt[(2, 2)].log_count
        - by_rt[(2, 1)].log_count
        - by_rt[(1, 2)].log_count
    )
    assert_close("exact entropy Creutz increment", chi_entropy, math.log(2.0))

    max_bag = max(row.max_bag for row in rows)
    print(f"\n  max elimination bag on tested fixed strips: {max_bag}")
    assert_equal("tested fixed-strip max bag", max_bag, 18)

    print(
        """
VERDICT:
  PASS.  The fixed-width transfer contraction gives the exact minimum-surface
  count in the tight strip:

      N_min^strip(R,T) = 2^((R-1)(T-1)).

  Thus the sector count in the Creutz script is exact for the fixed-strip
  transfer observable, not merely a sampled lower bound.  Its entropy density
  contributes Delta_R Delta_T log N_min = log 2, and the transfer law at fixed
  width is multiplication by 2^(R-1) per added column.

  The ledger class is therefore pinned: 2D tensor-network tractable at fixed
  width, exponential in loop width, and not O(1)-bond-dimension.  The remaining
  frontier is wider-margin/full-slab roughening entropy, not the exact
  fixed-strip N_min theorem.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
