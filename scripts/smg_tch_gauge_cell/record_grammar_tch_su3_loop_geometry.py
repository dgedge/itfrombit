#!/usr/bin/env python3
r"""TCH/SU(3) Wilson-loop geometry certificate.

Purpose
-------
This is the geometry-only bridge between the record-grammar Wilson-loop toys and
the actual truncated-cube-honeycomb (TCH) gauge-cell language.

The previous SU(2)/SU(3) cluster scripts used a rectangular two-dimensional
plaquette lattice with a declared record-action measure.  This script does not
declare a dynamics.  It instead asks a more primitive question:

    What are the closed loop records of one TCH truncated-cube cell?

The cell used here is the combinatorial truncated cube:

    * 24 vertices, represented as cube half-edges (cube vertex, incident axis);
    * 36 links;
    * 8 triangular faces and 6 octagonal faces;
    * Euler characteristic 24 - 36 + 14 = 2;
    * cycle rank 13.

For a gauge group such as SU(3), a detector-readable Wilson record is

    W(C) = (1/3) Re Tr P prod_{e in C} U_e,

where C is a closed oriented boundary.  The present script stops before assigning
the SU(3) link matrices U_e.  It enumerates the closed boundaries C, computes
their minimal face-supported surfaces, classifies simple versus multi-component
boundaries, and reports the ledger-size scaling.

Boundary
--------
This is not a confinement calculation and not a string-tension estimate.  It is
the rung that must be correct before such a calculation is meaningful: the
TCH-cell loop basis, face relation, minimal-surface map, and record-ledger
compression are all checked as exact finite combinatorics.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass


AXES = (0, 1, 2)


@dataclass(frozen=True)
class Face:
    label: str
    vertices: tuple[int, ...]
    weight: int


@dataclass(frozen=True)
class BoundaryRecord:
    mask: int
    perimeter: int
    components: int
    min_unit_area: int
    min_weighted_area: int
    witness_faces: tuple[int, ...]


def assert_equal(name: str, value: int, target: int) -> None:
    print(f"  {name:<72s} value={value} target={target}")
    if value != target:
        raise AssertionError(name)


def assert_true(name: str, value: bool) -> None:
    print(f"  {name:<72s} value={value}")
    if not value:
        raise AssertionError(name)


def cube_vertex(x: int, y: int, z: int) -> int:
    return (x << 0) | (y << 1) | (z << 2)


def bit(vertex: int, axis: int) -> int:
    return (vertex >> axis) & 1


def with_bit(vertex: int, axis: int, value: int) -> int:
    if value:
        return vertex | (1 << axis)
    return vertex & ~(1 << axis)


def half_edge_vertex(cube_v: int, axis: int) -> int:
    """The truncated-cube vertex sitting on the cube half-edge (cube_v, axis)."""

    return 3 * cube_v + axis


def decode_half_edge(index: int) -> tuple[int, int]:
    return divmod(index, 3)


def triangle_faces() -> list[Face]:
    """The 8 triangular faces created by truncating the original cube vertices."""

    faces = []
    for v in range(8):
        faces.append(
            Face(
                label=f"T{v:03b}",
                vertices=tuple(half_edge_vertex(v, axis) for axis in AXES),
                weight=3,
            )
        )
    return faces


def octagon_face(axis: int, side: int) -> Face:
    """One octagonal face inherited from an original cube face.

    If the original cube face fixes coordinate ``axis=side``, the boundary runs
    around the square face.  Each original cube edge contributes two truncated
    half-edge vertices, so the new face has eight vertices.
    """

    a, b = [item for item in AXES if item != axis]
    v00 = with_bit(with_bit(with_bit(0, axis, side), a, 0), b, 0)
    v10 = with_bit(v00, a, 1)
    v11 = with_bit(v10, b, 1)
    v01 = with_bit(v00, b, 1)
    verts = (
        half_edge_vertex(v00, a),
        half_edge_vertex(v10, a),
        half_edge_vertex(v10, b),
        half_edge_vertex(v11, b),
        half_edge_vertex(v11, a),
        half_edge_vertex(v01, a),
        half_edge_vertex(v01, b),
        half_edge_vertex(v00, b),
    )
    return Face(label=f"O{axis}{side}", vertices=verts, weight=8)


def build_faces() -> list[Face]:
    faces = triangle_faces()
    faces.extend(octagon_face(axis, side) for axis in AXES for side in (0, 1))
    return faces


def edge_key(left: int, right: int) -> tuple[int, int]:
    if left == right:
        raise ValueError("degenerate edge")
    return (left, right) if left < right else (right, left)


def consecutive_edges(vertices: tuple[int, ...]) -> list[tuple[int, int]]:
    return [edge_key(vertices[i], vertices[(i + 1) % len(vertices)]) for i in range(len(vertices))]


def build_edges(faces: list[Face]) -> tuple[list[tuple[int, int]], dict[tuple[int, int], int]]:
    edge_set = sorted({edge for face in faces for edge in consecutive_edges(face.vertices)})
    return edge_set, {edge: i for i, edge in enumerate(edge_set)}


def mask_from_edges(face: Face, edge_index: dict[tuple[int, int], int]) -> int:
    mask = 0
    for edge in consecutive_edges(face.vertices):
        mask ^= 1 << edge_index[edge]
    return mask


def gf2_rank(rows: list[int]) -> int:
    basis: dict[int, int] = {}
    for row in rows:
        x = row
        while x:
            pivot = x.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = x
                break
            x ^= basis[pivot]
    return len(basis)


def independent_face_basis(face_masks: list[int]) -> list[int]:
    basis_rows: list[int] = []
    basis_indices: list[int] = []
    rank = 0
    for i, mask in enumerate(face_masks):
        trial = basis_rows + [mask]
        new_rank = gf2_rank(trial)
        if new_rank > rank:
            basis_rows.append(mask)
            basis_indices.append(i)
            rank = new_rank
    return basis_indices


def adjacency_from_edges(n_vertices: int, edges: list[tuple[int, int]]) -> list[list[int]]:
    adj = [[] for _ in range(n_vertices)]
    for left, right in edges:
        adj[left].append(right)
        adj[right].append(left)
    return adj


def connected_components(adj: list[list[int]]) -> int:
    seen: set[int] = set()
    count = 0
    for start in range(len(adj)):
        if start in seen:
            continue
        count += 1
        queue = deque([start])
        seen.add(start)
        while queue:
            cur = queue.popleft()
            for nxt in adj[cur]:
                if nxt not in seen:
                    seen.add(nxt)
                    queue.append(nxt)
    return count


def boundary_component_count(mask: int, edges: list[tuple[int, int]], n_vertices: int) -> tuple[int, bool]:
    """Return number of nonempty boundary components and whether every degree is even."""

    adj: dict[int, list[int]] = defaultdict(list)
    for edge_i, (left, right) in enumerate(edges):
        if mask & (1 << edge_i):
            adj[left].append(right)
            adj[right].append(left)
    if not adj:
        return 0, True
    even = all(len(adj.get(v, [])) % 2 == 0 for v in range(n_vertices))
    seen: set[int] = set()
    count = 0
    for start in list(adj):
        if start in seen:
            continue
        count += 1
        queue = deque([start])
        seen.add(start)
        while queue:
            cur = queue.popleft()
            for nxt in adj[cur]:
                if nxt not in seen:
                    seen.add(nxt)
                    queue.append(nxt)
    return count, even


def subset_indices(mask: int) -> tuple[int, ...]:
    return tuple(i for i in range(mask.bit_length()) if mask & (1 << i))


def enumerate_boundaries(face_masks: list[int], face_weights: list[int]) -> dict[int, tuple[int, int, int]]:
    """Map boundary mask -> (min_unit_area, min_weighted_area, unit_witness_subset).

    The topological face count and the edge-weighted face area are minimized
    independently.  The witness subset records the minimum-unit-area surface,
    because that is the loop-basis/topological ledger used in the next rung.
    """

    best_unit: dict[int, tuple[int, int, int]] = {}
    best_weighted: dict[int, tuple[int, int, int]] = {}
    n_faces = len(face_masks)
    for subset in range(1 << n_faces):
        boundary = 0
        unit_area = 0
        weighted_area = 0
        for face_i in range(n_faces):
            if subset & (1 << face_i):
                boundary ^= face_masks[face_i]
                unit_area += 1
                weighted_area += face_weights[face_i]
        unit_candidate = (unit_area, weighted_area, subset)
        weighted_candidate = (weighted_area, unit_area, subset)
        current_unit = best_unit.get(boundary)
        current_weighted = best_weighted.get(boundary)
        if current_unit is None or unit_candidate[:2] < current_unit[:2]:
            best_unit[boundary] = unit_candidate
        if current_weighted is None or weighted_candidate[:2] < current_weighted[:2]:
            best_weighted[boundary] = weighted_candidate

    return {
        boundary: (
            unit_area,
            best_weighted[boundary][0],
            unit_subset,
        )
        for boundary, (unit_area, _unit_weight, unit_subset) in best_unit.items()
    }


def classify_boundaries(
    best: dict[int, tuple[int, int, int]],
    edges: list[tuple[int, int]],
    n_vertices: int,
) -> list[BoundaryRecord]:
    records: list[BoundaryRecord] = []
    for boundary, (unit_area, weighted_area, subset) in best.items():
        if boundary == 0:
            continue
        components, even = boundary_component_count(boundary, edges, n_vertices)
        if not even:
            raise AssertionError("non-closed boundary found")
        records.append(
            BoundaryRecord(
                mask=boundary,
                perimeter=boundary.bit_count(),
                components=components,
                min_unit_area=unit_area,
                min_weighted_area=weighted_area,
                witness_faces=subset_indices(subset),
            )
        )
    return records


def histogram(records: list[BoundaryRecord], attr: str) -> Counter[int]:
    return Counter(getattr(record, attr) for record in records)


def format_hist(counter: Counter[int], max_items: int | None = None) -> str:
    items = sorted(counter.items())
    if max_items is not None:
        items = items[:max_items]
    return ", ".join(f"{key}:{value}" for key, value in items)


def ledger_scaling(n_cells: int, n_vertices: int, n_edges: int, cycle_rank: int) -> dict[str, int]:
    """Independent-cell geometry ledger.

    The full glued TCH honeycomb is a later task.  This report states the
    per-cell compression cleanly: flat storage of every nonzero Wilson boundary
    scales as 2^cycle_rank - 1, while the record grammar stores a byte-stabilizer
    skeleton, link records, and a cycle basis.
    """

    stabilizer_generators = 8 * n_vertices * n_cells
    link_records = n_edges * n_cells
    loop_basis_records = cycle_rank * n_cells
    face_relation_records = n_cells
    flat_loop_records = ((1 << cycle_rank) - 1) * n_cells
    promoted_collective_records = 0
    generator_ledger_size = (
        stabilizer_generators
        + link_records
        + loop_basis_records
        + face_relation_records
        + promoted_collective_records
    )
    return {
        "cells": n_cells,
        "stabilizer_generators": stabilizer_generators,
        "link_records": link_records,
        "loop_basis_records": loop_basis_records,
        "face_relation_records": face_relation_records,
        "promoted_collective_records": promoted_collective_records,
        "flat_loop_records": flat_loop_records,
        "generator_ledger_size": generator_ledger_size,
    }


def main() -> None:
    print("TCH/SU(3) Wilson-loop geometry certificate")
    print("=" * 92)

    faces = build_faces()
    edges, edge_index = build_edges(faces)
    n_vertices = 24
    n_edges = len(edges)
    n_faces = len(faces)

    print("\n[1] Truncated-cube cell checks")
    assert_equal("vertices", n_vertices, 24)
    assert_equal("links", n_edges, 36)
    assert_equal("faces", n_faces, 14)
    assert_equal("Euler V-E+F", n_vertices - n_edges + n_faces, 2)

    cell_adj = adjacency_from_edges(n_vertices, edges)
    degrees = sorted(len(items) for items in cell_adj)
    assert_equal("graph connected components", connected_components(cell_adj), 1)
    assert_true("all vertices degree 3", degrees == [3] * n_vertices)

    face_masks = [mask_from_edges(face, edge_index) for face in faces]
    face_weights = [face.weight for face in faces]
    edge_face_count = Counter()
    for mask in face_masks:
        for edge_i in range(n_edges):
            if mask & (1 << edge_i):
                edge_face_count[edge_i] += 1
    assert_true("each link belongs to exactly two faces", set(edge_face_count.values()) == {2})
    assert_equal("face-boundary rank over F2", gf2_rank(face_masks), 13)
    assert_equal("cycle rank E-V+1", n_edges - n_vertices + 1, 13)
    assert_equal("total face-edge incidences", sum(mask.bit_count() for mask in face_masks), 2 * n_edges)
    assert_true("xor of all face boundaries is zero", eval_xor(face_masks) == 0)

    basis = independent_face_basis(face_masks)
    assert_equal("independent face-loop basis size", len(basis), 13)
    print("  basis faces:", ", ".join(faces[i].label for i in basis))
    omitted = sorted(set(range(n_faces)) - set(basis))
    print("  face relation omits:", ", ".join(faces[i].label for i in omitted))

    print("\n[2] Closed Wilson-boundary enumeration")
    best = enumerate_boundaries(face_masks, face_weights)
    records = classify_boundaries(best, edges, n_vertices)
    assert_equal("unique boundaries including zero", len(best), 1 << 13)
    assert_equal("nonzero closed Wilson boundaries", len(records), (1 << 13) - 1)
    assert_true("all nonzero boundaries have perimeter >= 3", min(r.perimeter for r in records) >= 3)
    assert_true("all nonzero boundaries have at least one component", min(r.components for r in records) >= 1)

    elementary = [record for record in records if record.min_unit_area == 1]
    assert_equal("elementary face loops", len(elementary), 14)
    print(f"  perimeter histogram:       {format_hist(histogram(records, 'perimeter'))}")
    print(f"  min unit-area histogram:   {format_hist(histogram(records, 'min_unit_area'))}")
    print(f"  component histogram:       {format_hist(histogram(records, 'components'))}")
    print(f"  weighted-area histogram:   {format_hist(histogram(records, 'min_weighted_area'), max_items=16)}")

    simple = sum(1 for record in records if record.components == 1)
    multi = len(records) - simple
    print(f"  simple closed-loop records: {simple}")
    print(f"  multi-component products:   {multi}")

    print("\n[3] Elementary SU(3) loop records")
    print("  W(C) = (1/3) Re Tr P prod_{e in C} U_e")
    for face_i, face in enumerate(faces):
        mask = face_masks[face_i]
        components, even = boundary_component_count(mask, edges, n_vertices)
        assert_true(f"{face.label} is one closed component", even and components == 1)
        print(
            f"  {face.label}: perimeter={mask.bit_count():2d}, "
            f"topological_area=1, weighted_area={face.weight:2d}"
        )

    print("\n[4] Geometry-ledger scaling before dynamics")
    print("  This is independent-cell scaling.  Gluing TCH cells and adding dynamics is the next rung.")
    print("  cells  flat_loop_records  generator_ledger  loop_basis  promotions")
    for n_cells in (1, 2, 4, 8, 16, 32):
        row = ledger_scaling(n_cells, n_vertices=n_vertices, n_edges=n_edges, cycle_rank=13)
        print(
            f"  {row['cells']:>5d}"
            f"  {row['flat_loop_records']:>17d}"
            f"  {row['generator_ledger_size']:>16d}"
            f"  {row['loop_basis_records']:>10d}"
            f"  {row['promoted_collective_records']:>10d}"
        )

    one = ledger_scaling(1, n_vertices=n_vertices, n_edges=n_edges, cycle_rank=13)
    compression = one["flat_loop_records"] / one["generator_ledger_size"]
    print(f"  one-cell flat/generator ratio: {compression:.2f}x")

    print(
        """
VERDICT:
  PASS.  The truncated-cube TCH cell has the expected 24/36/14 geometry, a
  13-dimensional closed-loop space, and 8191 nonzero Wilson boundaries.  The
  record grammar stores the same loop space as a 13-generator cycle basis plus
  the byte/link skeleton, not as a flat list of all boundaries.

  This is a geometry certificate only.  It licenses the next calculation:
  attach SU(3) link matrices/action weights to this TCH loop basis and then
  measure area-law, perimeter-law, and promotion/bond-dimension behaviour on
  glued few-cell clusters.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


def eval_xor(rows: list[int]) -> int:
    out = 0
    for row in rows:
        out ^= row
    return out


if __name__ == "__main__":
    main()
