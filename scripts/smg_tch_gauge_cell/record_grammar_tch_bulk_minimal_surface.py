#!/usr/bin/env python3
r"""Minimal 3D TCH bulk cluster: Wilson surfaces and ledger growth.

Purpose
-------
The shell certificates are exact but topologically special: the boundary is an
S^2 shell, so every loop has only the two shell caps.  This script keeps the
same truncated-cube cell geometry but stops discarding the glued octagonal
faces.  Those shared faces become interior plaquettes of a finite 3D bulk
cluster.

The goal is modest and explicit:

    build the smallest 3D TCH bulk family;
    verify the chain complex d2 d3 = 0;
    solve minimal Wilson surfaces using boundary and interior faces;
    report the surface degeneracy and generator-ledger scaling.

This is the first stress test of the bounded-ledger claim away from the
pure-boundary S^2 shell cap structure.  In this first linear bulk family,
small boundary loops do not yet shortcut through the interior; that negative is
recorded explicitly.  Interior Wilson records do exist once loops are allowed
inside the bulk.

What is not claimed
-------------------
This is not continuum Yang-Mills, not a finite-temperature Polyakov loop, and
not a physical QCD string tension.  It produces a lattice-unit strong-coupling
surface law on a small 3D TCH ball.  A physical string tension needs a lattice
spacing and a continuum/beta scaling argument.
"""

from __future__ import annotations

import math
from collections import defaultdict, deque
from dataclasses import dataclass
from itertools import combinations

import numpy as np

from record_grammar_tch_glued_surface_selector import (
    GLUE_LEFT_FACE,
    GLUE_RIGHT_FACE,
    UnionFind,
    global_vertex,
    local_faces_by_label,
)
from record_grammar_tch_glued_wilson_action_average import BETA, u_fundamental
from record_grammar_tch_su3_loop_geometry import Face, build_faces, consecutive_edges, gf2_rank


MAX_SURFACE_AREA = 4


@dataclass(frozen=True)
class BulkComplex:
    cells: int
    vertices: int
    edges: list[tuple[int, int]]
    faces: list[Face]
    face_masks: list[int]
    edge_index: dict[tuple[int, int], int]
    cell_masks: list[int]
    boundary_face_indices: tuple[int, ...]
    interior_face_indices: tuple[int, ...]


@dataclass(frozen=True)
class SurfaceEntry:
    area: int
    count: int
    witness: tuple[int, ...]


@dataclass(frozen=True)
class BulkRow:
    cells: int
    boundary_patch_area: int
    perimeter: int
    shell_area: int
    shell_count: int
    bulk_area: int
    bulk_count: int
    uses_interior: bool
    minus_log_leading: float


@dataclass(frozen=True)
class InteriorRow:
    cells: int
    patch_area: int
    perimeter: int
    bulk_area: int
    bulk_count: int
    uses_interior: bool


def assert_equal(name: str, value: int, target: int) -> None:
    print(f"  {name:<78s} value={value} target={target}")
    if value != target:
        raise AssertionError(name)


def assert_true(name: str, value: bool) -> None:
    print(f"  {name:<78s} value={value}")
    if not value:
        raise AssertionError(name)


def assert_less(name: str, value: float, bound: float) -> None:
    print(f"  {name:<78s} value={value:.12g} bound={bound:.12g}")
    if not value < bound:
        raise AssertionError(name)


def assert_close(name: str, value: float, target: float, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<78s} value={value:.12g} target={target:.12g} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def canonical_edge(left: int, right: int) -> tuple[int, int]:
    if left == right:
        raise ValueError("degenerate edge")
    return (left, right) if left < right else (right, left)


def build_bulk_complex(n_cells: int) -> BulkComplex:
    if n_cells < 1:
        raise ValueError("n_cells must be positive")

    local_faces = build_faces()
    by_label = local_faces_by_label()
    left_face = by_label[GLUE_LEFT_FACE]
    right_face = by_label[GLUE_RIGHT_FACE]

    uf = UnionFind.build(24 * n_cells)
    for cell in range(n_cells - 1):
        for left_v, right_v in zip(left_face.vertices, reversed(right_face.vertices), strict=True):
            uf.union(global_vertex(cell, left_v), global_vertex(cell + 1, right_v))

    root_to_new: dict[int, int] = {}

    def canon_vertex(old: int) -> int:
        root = uf.find(old)
        if root not in root_to_new:
            root_to_new[root] = len(root_to_new)
        return root_to_new[root]

    faces: list[Face] = []
    cell_face_index: dict[tuple[int, str], int] = {}
    interior: list[int] = []
    boundary: list[int] = []

    for cell in range(n_cells - 1):
        verts = tuple(canon_vertex(global_vertex(cell, vertex)) for vertex in left_face.vertices)
        face_i = len(faces)
        faces.append(Face(label=f"I{cell}-{cell + 1}:{GLUE_LEFT_FACE}/{GLUE_RIGHT_FACE}", vertices=verts, weight=8))
        interior.append(face_i)
        cell_face_index[(cell, GLUE_LEFT_FACE)] = face_i
        cell_face_index[(cell + 1, GLUE_RIGHT_FACE)] = face_i

    for cell in range(n_cells):
        for face in local_faces:
            if (cell, face.label) in cell_face_index:
                continue
            verts = tuple(canon_vertex(global_vertex(cell, vertex)) for vertex in face.vertices)
            face_i = len(faces)
            faces.append(Face(label=f"C{cell}:{face.label}", vertices=verts, weight=face.weight))
            boundary.append(face_i)
            cell_face_index[(cell, face.label)] = face_i

    edge_set = sorted({canonical_edge(a, b) for face in faces for a, b in consecutive_edges(face.vertices)})
    edge_index = {edge: i for i, edge in enumerate(edge_set)}
    face_masks: list[int] = []
    for face in faces:
        mask = 0
        for edge in consecutive_edges(face.vertices):
            mask ^= 1 << edge_index[edge]
        face_masks.append(mask)

    cell_masks: list[int] = []
    for cell in range(n_cells):
        mask = 0
        for face in local_faces:
            mask ^= 1 << cell_face_index[(cell, face.label)]
        cell_masks.append(mask)

    return BulkComplex(
        cells=n_cells,
        vertices=len(root_to_new),
        edges=edge_set,
        faces=faces,
        face_masks=face_masks,
        edge_index=edge_index,
        cell_masks=cell_masks,
        boundary_face_indices=tuple(boundary),
        interior_face_indices=tuple(interior),
    )


def xor_masks(masks: list[int], indices: tuple[int, ...]) -> int:
    out = 0
    for index in indices:
        out ^= masks[index]
    return out


def subset_perimeter(mask: int) -> int:
    return mask.bit_count()


def face_adjacency(face_masks: list[int], allowed: tuple[int, ...]) -> dict[int, set[int]]:
    edge_faces: dict[int, list[int]] = defaultdict(list)
    for face_i in allowed:
        mask = face_masks[face_i]
        edge = 0
        while mask:
            if mask & 1:
                edge_faces[edge].append(face_i)
            mask >>= 1
            edge += 1

    adj = {face_i: set() for face_i in allowed}
    for incident in edge_faces.values():
        for a, b in combinations(incident, 2):
            adj[a].add(b)
            adj[b].add(a)
    return adj


def connected_patches(face_masks: list[int], allowed: tuple[int, ...], max_size: int) -> list[tuple[int, ...]]:
    adj = face_adjacency(face_masks, allowed)
    patches: set[tuple[int, ...]] = set()
    for start in allowed:
        queue: deque[tuple[int, ...]] = deque([(start,)])
        while queue:
            patch = queue.popleft()
            patches.add(tuple(sorted(patch)))
            if len(patch) >= max_size:
                continue
            frontier = sorted({n for face in patch for n in adj[face]} - set(patch))
            for nxt in frontier:
                new_patch = tuple(sorted((*patch, nxt)))
                if new_patch not in patches:
                    queue.append(new_patch)
    return sorted(patches, key=lambda item: (len(item), item))


def enumerate_min_surfaces(face_masks: list[int], face_indices: tuple[int, ...], max_area: int) -> dict[int, SurfaceEntry]:
    table: dict[int, SurfaceEntry] = {0: SurfaceEntry(area=0, count=1, witness=())}
    for area in range(1, max_area + 1):
        for combo in combinations(face_indices, area):
            boundary = xor_masks(face_masks, combo)
            current = table.get(boundary)
            if current is None:
                table[boundary] = SurfaceEntry(area=area, count=1, witness=combo)
            elif current.area == area:
                table[boundary] = SurfaceEntry(area=area, count=current.count + 1, witness=current.witness)
    return table


def topology_checks(complex_: BulkComplex) -> None:
    rank_d2 = gf2_rank(complex_.face_masks)
    rank_d3 = gf2_rank(complex_.cell_masks)
    euler_ball = complex_.vertices - len(complex_.edges) + len(complex_.faces) - complex_.cells
    assert_equal(f"N={complex_.cells} Euler characteristic of 3-ball", euler_ball, 1)
    assert_equal(f"N={complex_.cells} rank d3", rank_d3, complex_.cells)
    assert_equal(f"N={complex_.cells} rank d2 = F - C", rank_d2, len(complex_.faces) - complex_.cells)
    for cell_i, cell_mask in enumerate(complex_.cell_masks):
        face_indices = tuple(i for i in range(len(complex_.faces)) if cell_mask & (1 << i))
        assert_equal(f"N={complex_.cells} cell {cell_i} has 14 faces", len(face_indices), 14)
        assert_equal(f"N={complex_.cells} d2 d3 cell {cell_i}", xor_masks(complex_.face_masks, face_indices), 0)


def bulk_rows(cells: int, max_patch_size: int = MAX_SURFACE_AREA) -> list[BulkRow]:
    complex_ = build_bulk_complex(cells)
    shell_table = enumerate_min_surfaces(complex_.face_masks, complex_.boundary_face_indices, max_patch_size)
    bulk_table = enumerate_min_surfaces(
        complex_.face_masks,
        tuple(range(len(complex_.faces))),
        max_patch_size,
    )
    u = u_fundamental(BETA)
    rows: list[BulkRow] = []
    seen: set[int] = set()
    for patch in connected_patches(complex_.face_masks, complex_.boundary_face_indices, max_patch_size):
        boundary = xor_masks(complex_.face_masks, patch)
        if boundary == 0 or boundary in seen:
            continue
        seen.add(boundary)
        shell = shell_table[boundary]
        bulk = bulk_table[boundary]
        uses_interior = any(face_i in complex_.interior_face_indices for face_i in bulk.witness)
        rows.append(
            BulkRow(
                cells=cells,
                boundary_patch_area=len(patch),
                perimeter=subset_perimeter(boundary),
                shell_area=shell.area,
                shell_count=shell.count,
                bulk_area=bulk.area,
                bulk_count=bulk.count,
                uses_interior=uses_interior,
                minus_log_leading=-(math.log(bulk.count) + bulk.area * math.log(u)),
            )
        )
    return rows


def interior_rows(cells: int, max_patch_size: int = MAX_SURFACE_AREA) -> list[InteriorRow]:
    complex_ = build_bulk_complex(cells)
    bulk_table = enumerate_min_surfaces(
        complex_.face_masks,
        tuple(range(len(complex_.faces))),
        max_patch_size,
    )
    rows: list[InteriorRow] = []
    seen: set[int] = set()
    for patch in connected_patches(complex_.face_masks, tuple(range(len(complex_.faces))), max_patch_size):
        if not any(face_i in complex_.interior_face_indices for face_i in patch):
            continue
        boundary = xor_masks(complex_.face_masks, patch)
        if boundary == 0 or boundary in seen:
            continue
        seen.add(boundary)
        bulk = bulk_table[boundary]
        rows.append(
            InteriorRow(
                cells=cells,
                patch_area=len(patch),
                perimeter=subset_perimeter(boundary),
                bulk_area=bulk.area,
                bulk_count=bulk.count,
                uses_interior=any(face_i in complex_.interior_face_indices for face_i in bulk.witness),
            )
        )
    return rows


def fit_bulk_area(rows: list[BulkRow]) -> tuple[float, float, float, float]:
    x = np.array([[row.bulk_area, row.perimeter, 1.0] for row in rows], dtype=float)
    y = np.array([row.minus_log_leading for row in rows], dtype=float)
    coeff, *_ = np.linalg.lstsq(x, y, rcond=None)
    resid = y - x @ coeff
    return float(coeff[0]), float(coeff[1]), float(coeff[2]), float(np.max(np.abs(resid)))


def ledger_row(complex_: BulkComplex, max_degeneracy: int) -> dict[str, int]:
    rank_d2 = gf2_rank(complex_.face_masks)
    rank_d3 = gf2_rank(complex_.cell_masks)
    promotions = math.ceil(math.log2(max(1, max_degeneracy)))
    generator = 8 * complex_.vertices + len(complex_.edges) + rank_d2 + rank_d3 + 1 + promotions
    return {
        "rank_d2": rank_d2,
        "rank_d3": rank_d3,
        "generator": generator,
        "promotions": promotions,
        "bond_dimension_proxy": 2**promotions,
    }


def main() -> None:
    print("Minimal 3D TCH bulk cluster: Wilson surfaces and ledger growth")
    print("=" * 104)
    u = u_fundamental(BETA)
    sigma = -math.log(u)
    print(f"  beta={BETA:.3f}; u_F=beta/18={u:.9f}; sigma_lattice=-ln(u_F)={sigma:.9f}")
    print("  surfaces are enumerated up to area 4, enough for the boundary patches tested here.")

    print("\n[1] Bulk chain-complex checks")
    for cells in range(1, 7):
        complex_ = build_bulk_complex(cells)
        topology_checks(complex_)
        assert_equal(f"N={cells} boundary face count", len(complex_.boundary_face_indices), 12 * cells + 2)
        assert_equal(f"N={cells} interior face count", len(complex_.interior_face_indices), max(0, cells - 1))

    print("\n[2] Boundary-loop survey: does the linear bulk beat the S2 shell cap?")
    all_rows: list[BulkRow] = []
    max_degeneracy_by_cells: dict[int, int] = {}
    shortcut_by_cells: dict[int, int] = {}
    interior_by_cells: dict[int, int] = {}
    for cells in range(1, 7):
        rows = bulk_rows(cells)
        all_rows.extend(rows)
        max_degeneracy = max(row.bulk_count for row in rows)
        shortcuts = sum(row.bulk_area < row.shell_area for row in rows)
        interior = sum(row.uses_interior for row in rows)
        max_degeneracy_by_cells[cells] = max_degeneracy
        shortcut_by_cells[cells] = shortcuts
        interior_by_cells[cells] = interior
        fit_sigma, fit_mu, fit_const, fit_resid = fit_bulk_area(rows)
        print(
            f"  N={cells}: boundary_loops={len(rows):5d}, interior_min={interior:4d}, "
            f"shortcuts={shortcuts:4d}, max_deg={max_degeneracy:3d}, "
            f"sigma_fit={fit_sigma:.9f}, mu={fit_mu:.2e}, max_resid={fit_resid:.3e}"
        )
        assert_close(f"N={cells} leading bulk sigma", fit_sigma, sigma, tol=1e-10)
        assert_less(f"N={cells} perimeter term", abs(fit_mu), 1e-10)
        assert_less(f"N={cells} fit residual from degeneracy bookkeeping", fit_resid, 2.0)

    assert_equal("linear-chain boundary loops using interior minima", sum(interior_by_cells.values()), 0)
    assert_equal("linear-chain boundary shortcuts at area <= 4", sum(shortcut_by_cells.values()), 0)
    assert_less("max observed minimum-surface degeneracy", float(max(max_degeneracy_by_cells.values())), 16.0)

    print("\n[3] Interior Wilson records: loops that the shell cannot represent")
    max_interior_degeneracy = 1
    interior_record_count = 0
    for cells in range(2, 7):
        rows = interior_rows(cells)
        interior_record_count += len(rows)
        interior_min = sum(row.uses_interior for row in rows)
        max_deg = max(row.bulk_count for row in rows) if rows else 1
        max_interior_degeneracy = max(max_interior_degeneracy, max_deg)
        print(
            f"  N={cells}: interior_loop_records={len(rows):4d}, "
            f"interior_min={interior_min:4d}, max_deg={max_deg:3d}"
        )
        assert_true(f"N={cells} has interior Wilson records", len(rows) > 0)
        assert_true(f"N={cells} some interior records minimize on interior plaquettes", interior_min > 0)

    assert_true("interior Wilson records exist in the bulk family", interior_record_count > 0)
    assert_less("max interior minimum-surface degeneracy", float(max_interior_degeneracy), 16.0)

    print("\n[4] Ledger scaling in the first 3D bulk family")
    print("  cells  V    E    F   C3  rank_d2 rank_d3  max_deg promotions bondD generator")
    ledger_points: list[tuple[int, int]] = []
    for cells in (1, 2, 3, 4, 5, 6, 8):
        complex_ = build_bulk_complex(cells)
        if cells <= 6:
            max_deg = max(max_degeneracy_by_cells[cells], max_interior_degeneracy)
        else:
            max_deg = max(max(max_degeneracy_by_cells.values()), max_interior_degeneracy)
        led = ledger_row(complex_, max_deg)
        ledger_points.append((cells, led["generator"]))
        print(
            f"  {cells:>5d}"
            f"  {complex_.vertices:>3d}"
            f"  {len(complex_.edges):>4d}"
            f"  {len(complex_.faces):>4d}"
            f"  {complex_.cells:>3d}"
            f"  {led['rank_d2']:>7d}"
            f"  {led['rank_d3']:>7d}"
            f"  {max_deg:>7d}"
            f"  {led['promotions']:>10d}"
            f"  {led['bond_dimension_proxy']:>5d}"
            f"  {led['generator']:>9d}"
        )
        assert_less(f"N={cells} bulk generator ledger grows linearly", led["generator"], 290 * cells)
        assert_less(f"N={cells} bulk promotion count bounded in tested family", led["promotions"], 4)

    slope, intercept = np.polyfit(
        np.array([point[0] for point in ledger_points], dtype=float),
        np.array([point[1] for point in ledger_points], dtype=float),
        deg=1,
    )
    print(f"  fitted generator ledger: {slope:.6f} N + {intercept:.6f}")
    assert_close("bulk generator ledger slope", float(slope), 169.0, tol=1e-10)
    assert_close("bulk generator ledger intercept", float(intercept), 74.0, tol=1e-10)

    print(
        """
VERDICT:
  PASS.  Retaining the glued octagons as interior plaquettes gives a genuine
  finite 3D TCH bulk ball: chi=1, rank(d3)=C, rank(d2)=F-C, and d2 d3=0 on
  every cell.  Boundary loops of area <=4 on this linear chain do not shortcut
  through the interior; their minima remain shell-cap minima.  But interior
  Wilson records exist, minimize on interior plaquettes, and are invisible to
  the boundary-only shell model.  The leading strong-coupling weight still gives
  sigma_lattice=-ln(beta/18), with degeneracy recorded separately.

  In this first bulk family the generator ledger remains linear,
  L_gen=169N+74, and the observed minimum-surface degeneracy needs only a
  bounded promotion register.  This is a bulk stress test of the record-grammar
  tractability claim, not yet a continuum-QCD or physical-string-tension result.
  The next bulk rung should use a branched/periodic 3+1D complex to test true
  boundary shortcuts and genuine Polyakov loops.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
