#!/usr/bin/env python3
r"""Bond-bipyramid TCH isoperimetry gate for honest Wilson loops.

Purpose
-------
``record_grammar_tch_bond_bipyramid_static_potential.py`` deliberately retired
one tempting Wilson-loop estimator: its rectangular interior face family has
``P_edges=2A``.  Area and perimeter are therefore collinear regressors, so no
Creutz ratio or least-squares estimator can distinguish a string-tension term
from a perimeter/screening term on that family.

This script asks the geometric question before building a new loop observable:

    Does the corrected bond-bipyramid slab contain flat 2D face sheets at all?

It tests the largest connected face plane in L x L x 3 slabs.  A genuine flat
sheet should have

    lambda_2(normalized Laplacian) ~ 1/L^2,
    boundary/area -> 0 for large patches,
    internal-edge density I/A -> 3/2.

An expander-like plane would instead keep a spectral gap bounded below and
would force perimeter proportional to area.  In that case a clean dynamical
string-tension observable would not be available on this geometry.

Boundary
--------
This is a geometry gate, not a Yang-Mills calculation.  It does not count
minimal spanning surfaces and does not produce a physical string tension.  It
decides whether the next observable, simple edge-loop Wilson families on a
flat sheet, is licensed by the bond-bipyramid geometry.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict, deque
from dataclasses import dataclass

import numpy as np

from record_grammar_tch_bond_bipyramid_bulk import (
    AXES,
    BondComplex,
    build_bond_complex,
    bonds_box,
    canonical_edge,
    xor_masks,
)
from record_grammar_tch_bond_bipyramid_static_potential import (
    boundary_stats,
    face_sum,
    rectangular_patch,
)


SIZES = (4, 6, 8, 10, 12)
THICKNESS = 3


@dataclass(frozen=True)
class PlaneGraph:
    axis: int
    key: int
    nodes: tuple[int, ...]
    adj: dict[int, set[int]]
    complex_: BondComplex


@dataclass(frozen=True)
class PlaneRow:
    size: int
    axis: int
    key: int
    faces: int
    lambda2: float
    lambda2_n: float
    lambda2_l2: float
    full_boundary: int
    full_p_over_a: float
    full_p_over_sqrt_a: float
    full_internal_over_a: float
    best_ball_area: int
    best_ball_boundary: int
    best_ball_p_over_sqrt_a: float
    bad_rect_p_over_a: float


def assert_true(name: str, value: bool) -> None:
    print(f"  {name:<78s} value={value}")
    if not value:
        raise AssertionError(name)


def assert_less(name: str, value: float, bound: float) -> None:
    print(f"  {name:<78s} value={value:.12g} bound={bound:.12g}")
    if not value < bound:
        raise AssertionError(name)


def assert_greater(name: str, value: float, bound: float) -> None:
    print(f"  {name:<78s} value={value:.12g} bound={bound:.12g}")
    if not value > bound:
        raise AssertionError(name)


def build_plane_adjacency(complex_: BondComplex, axis: int, key: int) -> dict[int, set[int]]:
    nodes = [
        face_i
        for face_i, face in enumerate(complex_.faces)
        if complex_.face_counts[face] == 2 and face_sum(face)[axis] == key
    ]
    edge_faces: dict[tuple[tuple[int, int, int], tuple[int, int, int]], list[int]] = defaultdict(list)
    for face_i in nodes:
        face = complex_.faces[face_i]
        for left, right in ((0, 1), (1, 2), (0, 2)):
            edge_faces[canonical_edge(face[left], face[right])].append(face_i)

    adj = {face_i: set() for face_i in nodes}
    for incident in edge_faces.values():
        for left in incident:
            if left not in adj:
                continue
            for right in incident:
                if right != left and right in adj:
                    adj[left].add(right)
    return adj


def components(adj: dict[int, set[int]]) -> list[tuple[int, ...]]:
    out: list[tuple[int, ...]] = []
    seen: set[int] = set()
    for start in adj:
        if start in seen:
            continue
        queue = deque([start])
        seen.add(start)
        comp: list[int] = []
        while queue:
            cur = queue.popleft()
            comp.append(cur)
            for nxt in adj[cur]:
                if nxt not in seen:
                    seen.add(nxt)
                    queue.append(nxt)
        out.append(tuple(sorted(comp)))
    return out


def largest_connected_plane(complex_: BondComplex) -> PlaneGraph:
    best: PlaneGraph | None = None
    for axis in AXES:
        keys = sorted({face_sum(face)[axis] for face in complex_.faces if complex_.face_counts[face] == 2})
        for key in keys:
            adj = build_plane_adjacency(complex_, axis, key)
            if not adj:
                continue
            comp = max(components(adj), key=len)
            comp_adj = {node: {nxt for nxt in adj[node] if nxt in comp} for node in comp}
            candidate = PlaneGraph(axis=axis, key=key, nodes=comp, adj=comp_adj, complex_=complex_)
            if best is None or len(candidate.nodes) > len(best.nodes):
                best = candidate
    if best is None:
        raise AssertionError("no connected face plane found")
    return best


def normalized_laplacian_gap(graph: PlaneGraph) -> float:
    index = {node: i for i, node in enumerate(graph.nodes)}
    n = len(graph.nodes)
    adjacency = np.zeros((n, n), dtype=float)
    for node in graph.nodes:
        for nbr in graph.adj[node]:
            adjacency[index[node], index[nbr]] = 1.0
    degree = adjacency.sum(axis=1)
    dinv = np.diag([1.0 / math.sqrt(item) if item > 0.0 else 0.0 for item in degree])
    laplacian = np.eye(n) - dinv @ adjacency @ dinv
    values = np.linalg.eigvalsh(laplacian)
    return float(values[1])


def boundary_edges(complex_: BondComplex, face_indices: set[int] | tuple[int, ...]) -> int:
    return xor_masks(complex_.face_masks, tuple(sorted(face_indices))).bit_count()


def internal_edges_per_face(area: int, boundary: int) -> float:
    # Triangular patch identity: P = 3A - 2I.
    return (3.0 * area - boundary) / (2.0 * area)


def bfs_ball_records(graph: PlaneGraph) -> list[tuple[int, int, float]]:
    records: list[tuple[int, int, float]] = []
    for start in graph.nodes:
        seen = {start}
        frontier = {start}
        while frontier:
            area = len(seen)
            boundary = boundary_edges(graph.complex_, seen)
            records.append((area, boundary, boundary / math.sqrt(area)))
            nxt: set[int] = set()
            for node in frontier:
                for nbr in graph.adj[node]:
                    if nbr not in seen:
                        nxt.add(nbr)
            seen |= nxt
            frontier = nxt
    return records


def bad_rectangular_ratio(complex_: BondComplex) -> float:
    patch = rectangular_patch(complex_, 2, 2)
    boundary = xor_masks(complex_.face_masks, patch)
    _, even, _, perimeter = boundary_stats(boundary, complex_.edges)
    if not even:
        raise AssertionError("bad rectangular control was not an even boundary")
    return perimeter / len(patch)


def measure_size(size: int) -> PlaneRow:
    complex_ = build_bond_complex(f"bond-slab-{size}x{size}x{THICKNESS}", bonds_box(size, size, THICKNESS))
    graph = largest_connected_plane(complex_)
    gap = normalized_laplacian_gap(graph)
    area = len(graph.nodes)
    full_boundary = boundary_edges(complex_, graph.nodes)
    balls = bfs_ball_records(graph)
    candidates = [item for item in balls if 0.20 * area <= item[0] <= 0.80 * area]
    if not candidates:
        candidates = balls
    best_area, best_boundary, best_p_sqrt = min(candidates, key=lambda item: item[2])
    return PlaneRow(
        size=size,
        axis=graph.axis,
        key=graph.key,
        faces=area,
        lambda2=gap,
        lambda2_n=gap * area,
        lambda2_l2=gap * size * size,
        full_boundary=full_boundary,
        full_p_over_a=full_boundary / area,
        full_p_over_sqrt_a=full_boundary / math.sqrt(area),
        full_internal_over_a=internal_edges_per_face(area, full_boundary),
        best_ball_area=best_area,
        best_ball_boundary=best_boundary,
        best_ball_p_over_sqrt_a=best_p_sqrt,
        bad_rect_p_over_a=bad_rectangular_ratio(complex_),
    )


def main() -> None:
    print("Bond-bipyramid TCH isoperimetry gate for Wilson loops")
    print("=" * 96)

    print("\n[1] Largest connected face-plane scaling")
    rows = [measure_size(size) for size in SIZES]
    print("  L  axis key faces  lambda2     lambda2*n  lambda2*L^2  P/A(full) I/A(full) bestA bestP P/sqrtA  badP/A")
    for row in rows:
        print(
            f"  {row.size:>2d}"
            f"  {row.axis:>4d}"
            f" {row.key:>3d}"
            f" {row.faces:>5d}"
            f"  {row.lambda2:>9.6f}"
            f"  {row.lambda2_n:>9.6f}"
            f"  {row.lambda2_l2:>11.6f}"
            f"  {row.full_p_over_a:>9.6f}"
            f"  {row.full_internal_over_a:>9.6f}"
            f" {row.best_ball_area:>5d}"
            f" {row.best_ball_boundary:>5d}"
            f" {row.best_ball_p_over_sqrt_a:>7.3f}"
            f" {row.bad_rect_p_over_a:>7.3f}"
        )

    gaps = [row.lambda2 for row in rows]
    full_p_over_a = [row.full_p_over_a for row in rows]
    full_internal = [row.full_internal_over_a for row in rows]
    bad_ratios = [row.bad_rect_p_over_a for row in rows]

    print("\n[2] Gates")
    assert_true("normalized spectral gap decreases with L", all(a > b for a, b in zip(gaps, gaps[1:])))
    assert_less("largest-plane lambda2*n stays bounded", max(row.lambda2_n for row in rows), 2.2)
    assert_less("largest-plane lambda2*L^2 stays bounded", max(row.lambda2_l2 for row in rows), 0.75)
    assert_true("full-plane boundary/area decreases with L", all(a > b for a, b in zip(full_p_over_a, full_p_over_a[1:])))
    assert_true("full-plane internal-edge density increases with L", all(a < b for a, b in zip(full_internal, full_internal[1:])))
    assert_greater("largest tested internal-edge density moves toward flat 3/2 limit", full_internal[-1], 1.40)
    assert_true("bad rectangular family remains collinear P/A=2", all(abs(item - 2.0) < 1e-12 for item in bad_ratios))

    print(
        """
VERDICT:
  PASS.  The bad rectangular interior-face family is genuinely unusable for an
  area/perimeter separation: it keeps P/A=2 exactly.  But that is not a global
  obstruction of the bond-bipyramid slab.  The largest connected offset face
  planes have a normalized spectral gap that falls like 1/L^2, boundary/area
  that decreases with L, and internal-edge density moving toward the triangular
  flat-sheet limit I/A=3/2.

  Therefore an honest edge-loop Wilson observable is geometrically licensed.
  The current static-potential patch family should remain a lower-bound
  entropy probe only.  The next dynamical rung is a simple edge-loop family on
  these connected planes, with Creutz ratios for the area/perimeter split and a
  fixed-thickness transfer matrix for exact N_min.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
