#!/usr/bin/env python3
r"""Simple-cubic SU(3) static-potential baseline.

Purpose
-------
The naive fat TCH block fails the 3-ball topology gate.  Before debugging the
TCH bulk-completion rule, this script checks the Wilson-loop extractor on a
licensed 3D complex: a finite L x L x L simple-cubic block.

This matters for two reasons.

1. It isolates method from geometry.  On a valid cubical ball, the leading
   strong-coupling Wilson rule must return

       -log <W(R,T)> = sigma R T,     V(R) = sigma R,

   for rectangular loops whose unique minimal surface is the flat rectangle.

2. It probes the "which lattice carries colour?" fork.  The current canon's
   geometry correction treats the simple-cubic web as the macroscopic gauge
   dual more cleanly than a naive two-cell truncated-cube/octahedron tiling.
   This script therefore supplies a dual-web baseline while the primal TCH
   bulk-completion theorem remains open.

Scope
-----
This is still leading strong coupling in lattice units:

    u_F = beta / 18,       sigma_lattice = -log u_F.

It is a validator for the static-potential extractor and the bounded-ledger
watch, not a continuum QCD calculation.
"""

from __future__ import annotations

import math
from collections import defaultdict, deque
from dataclasses import dataclass
from itertools import combinations

from record_grammar_tch_glued_wilson_action_average import BETA, u_fundamental
from record_grammar_tch_su3_loop_geometry import gf2_rank


AXES = (0, 1, 2)
MAX_RECT = 2
MAX_SURFACE_AREA = MAX_RECT * MAX_RECT


@dataclass(frozen=True)
class CubicalComplex:
    size: int
    cells: int
    vertices: int
    edges: list[tuple[int, int]]
    faces: list[tuple[int, int, int, int]]
    face_masks: list[int]
    cell_masks: list[int]
    face_index: dict[tuple[int, int, int, int], int]


@dataclass(frozen=True)
class SurfaceEntry:
    area: int
    count: int
    witness: tuple[int, ...]


@dataclass(frozen=True)
class RectRow:
    size: int
    normal_axis: int
    plane_position: int
    r_axis: int
    t_axis: int
    r: int
    t: int
    perimeter: int
    min_area: int
    degeneracy: int
    minus_log: float


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


def other_axes(axis: int) -> tuple[int, int]:
    out = tuple(item for item in AXES if item != axis)
    assert len(out) == 2
    return out  # type: ignore[return-value]


def vertex_index(coord: tuple[int, int, int], size: int) -> int:
    x, y, z = coord
    return x + (size + 1) * (y + (size + 1) * z)


def face_key(normal_axis: int, plane: int, a_coord: int, b_coord: int) -> tuple[int, int, int, int]:
    return normal_axis, plane, a_coord, b_coord


def face_vertices(key: tuple[int, int, int, int], size: int) -> tuple[int, int, int, int]:
    normal_axis, plane, a_coord, b_coord = key
    a_axis, b_axis = other_axes(normal_axis)

    def coord(da: int, db: int) -> tuple[int, int, int]:
        out = [0, 0, 0]
        out[normal_axis] = plane
        out[a_axis] = a_coord + da
        out[b_axis] = b_coord + db
        return tuple(out)  # type: ignore[return-value]

    return (
        vertex_index(coord(0, 0), size),
        vertex_index(coord(1, 0), size),
        vertex_index(coord(1, 1), size),
        vertex_index(coord(0, 1), size),
    )


def canonical_edge(left: int, right: int) -> tuple[int, int]:
    if left == right:
        raise ValueError("degenerate edge")
    return (left, right) if left < right else (right, left)


def consecutive_edges(vertices: tuple[int, ...]) -> list[tuple[int, int]]:
    return [canonical_edge(vertices[i], vertices[(i + 1) % len(vertices)]) for i in range(len(vertices))]


def build_cubical_complex(size: int) -> CubicalComplex:
    if size < 1:
        raise ValueError("size must be positive")

    faces: list[tuple[int, int, int, int]] = []
    face_index: dict[tuple[int, int, int, int], int] = {}
    for normal_axis in AXES:
        for plane in range(size + 1):
            for a_coord in range(size):
                for b_coord in range(size):
                    key = face_key(normal_axis, plane, a_coord, b_coord)
                    face_index[key] = len(faces)
                    faces.append(key)

    edge_set = sorted(
        {
            edge
            for key in faces
            for edge in consecutive_edges(face_vertices(key, size))
        }
    )
    edge_index = {edge: i for i, edge in enumerate(edge_set)}
    face_masks: list[int] = []
    for key in faces:
        mask = 0
        for edge in consecutive_edges(face_vertices(key, size)):
            mask ^= 1 << edge_index[edge]
        face_masks.append(mask)

    cell_masks: list[int] = []
    for x in range(size):
        for y in range(size):
            for z in range(size):
                base = (x, y, z)
                mask = 0
                for axis in AXES:
                    a_axis, b_axis = other_axes(axis)
                    for side in (0, 1):
                        plane = base[axis] + side
                        a_coord = base[a_axis]
                        b_coord = base[b_axis]
                        mask ^= 1 << face_index[face_key(axis, plane, a_coord, b_coord)]
                cell_masks.append(mask)

    return CubicalComplex(
        size=size,
        cells=size**3,
        vertices=(size + 1) ** 3,
        edges=edge_set,
        faces=faces,
        face_masks=face_masks,
        cell_masks=cell_masks,
        face_index=face_index,
    )


def xor_masks(masks: list[int], indices: tuple[int, ...]) -> int:
    out = 0
    for index in indices:
        out ^= masks[index]
    return out


def boundary_components(mask: int, edges: list[tuple[int, int]]) -> tuple[int, bool]:
    adj: dict[int, list[int]] = defaultdict(list)
    for edge_i, (left, right) in enumerate(edges):
        if mask & (1 << edge_i):
            adj[left].append(right)
            adj[right].append(left)
    if not adj:
        return 0, True
    even = all(len(items) % 2 == 0 for items in adj.values())
    seen: set[int] = set()
    components = 0
    for start in adj:
        if start in seen:
            continue
        components += 1
        queue: deque[int] = deque([start])
        seen.add(start)
        while queue:
            cur = queue.popleft()
            for nxt in adj[cur]:
                if nxt not in seen:
                    seen.add(nxt)
                    queue.append(nxt)
    return components, even


def face_adjacency(complex_: CubicalComplex) -> dict[int, set[int]]:
    edge_faces: dict[int, list[int]] = defaultdict(list)
    for face_i, mask in enumerate(complex_.face_masks):
        edge_i = 0
        while mask:
            if mask & 1:
                edge_faces[edge_i].append(face_i)
            mask >>= 1
            edge_i += 1

    adj = {face_i: set() for face_i in range(len(complex_.faces))}
    for incident in edge_faces.values():
        for a, b in combinations(incident, 2):
            adj[a].add(b)
            adj[b].add(a)
    return adj


def connected_patches(complex_: CubicalComplex, max_area: int) -> list[tuple[int, ...]]:
    adj = face_adjacency(complex_)
    patches: set[tuple[int, ...]] = set()
    for start in range(len(complex_.faces)):
        queue: deque[tuple[int, ...]] = deque([(start,)])
        while queue:
            patch = queue.popleft()
            patches.add(patch)
            if len(patch) >= max_area:
                continue
            frontier = sorted({n for face in patch for n in adj[face]} - set(patch))
            for nxt in frontier:
                new_patch = tuple(sorted((*patch, nxt)))
                if new_patch not in patches:
                    queue.append(new_patch)
    return sorted(patches, key=lambda item: (len(item), item))


def min_surface_table(complex_: CubicalComplex, max_area: int) -> dict[int, SurfaceEntry]:
    table: dict[int, SurfaceEntry] = {0: SurfaceEntry(area=0, count=1, witness=())}
    for patch in connected_patches(complex_, max_area):
        boundary = xor_masks(complex_.face_masks, patch)
        area = len(patch)
        current = table.get(boundary)
        if current is None:
            table[boundary] = SurfaceEntry(area=area, count=1, witness=patch)
        elif current.area == area:
            table[boundary] = SurfaceEntry(area=area, count=current.count + 1, witness=current.witness)
    return table


def rectangle_patch(
    complex_: CubicalComplex,
    normal_axis: int,
    plane: int,
    r_axis: int,
    t_axis: int,
    r: int,
    t: int,
) -> tuple[int, ...]:
    if {normal_axis, r_axis, t_axis} != set(AXES):
        raise ValueError("normal, R, and T axes must be distinct")
    span = other_axes(normal_axis)
    patch = []
    for rr in range(r):
        for tt in range(t):
            coords = {r_axis: rr, t_axis: tt}
            patch.append(
                complex_.face_index[
                    face_key(normal_axis, plane, coords[span[0]], coords[span[1]])
                ]
            )
    return tuple(sorted(patch))


def topology_checks(complex_: CubicalComplex) -> None:
    rank_d2 = gf2_rank(complex_.face_masks)
    rank_d3 = gf2_rank(complex_.cell_masks)
    euler = complex_.vertices - len(complex_.edges) + len(complex_.faces) - complex_.cells
    assert_equal(f"L={complex_.size} Euler characteristic", euler, 1)
    assert_equal(f"L={complex_.size} rank d3", rank_d3, complex_.cells)
    assert_equal(f"L={complex_.size} rank d2 = F-C", rank_d2, len(complex_.faces) - complex_.cells)
    for cell_i, cell_mask in enumerate(complex_.cell_masks):
        face_indices = tuple(i for i in range(len(complex_.faces)) if cell_mask & (1 << i))
        assert_equal(f"L={complex_.size} cell {cell_i} face count", len(face_indices), 6)
        assert_equal(f"L={complex_.size} d2 d3 cell {cell_i}", xor_masks(complex_.face_masks, face_indices), 0)


def rectangle_rows(size: int) -> list[RectRow]:
    complex_ = build_cubical_complex(size)
    table = min_surface_table(complex_, MAX_SURFACE_AREA)
    sigma = -math.log(u_fundamental(BETA))
    rows: list[RectRow] = []
    plane = size // 2
    for normal_axis in AXES:
        span_axes = other_axes(normal_axis)
        for r_axis, t_axis in (span_axes, tuple(reversed(span_axes))):
            for r in range(1, min(MAX_RECT, size) + 1):
                for t in range(1, min(MAX_RECT, size) + 1):
                    patch = rectangle_patch(complex_, normal_axis, plane, r_axis, t_axis, r, t)
                    boundary = xor_masks(complex_.face_masks, patch)
                    components, even = boundary_components(boundary, complex_.edges)
                    if not even or components != 1:
                        continue
                    entry = table.get(boundary)
                    if entry is None:
                        continue
                    rows.append(
                        RectRow(
                            size=size,
                            normal_axis=normal_axis,
                            plane_position=plane,
                            r_axis=r_axis,
                            t_axis=t_axis,
                            r=r,
                            t=t,
                            perimeter=boundary.bit_count(),
                            min_area=entry.area,
                            degeneracy=entry.count,
                            minus_log=sigma * entry.area - math.log(entry.count),
                        )
                    )
    return rows


def ledger_row(complex_: CubicalComplex, max_degeneracy: int) -> dict[str, int]:
    rank_d2 = gf2_rank(complex_.face_masks)
    rank_d3 = gf2_rank(complex_.cell_masks)
    promotions = math.ceil(math.log2(max(1, max_degeneracy)))
    # Same conservative endpoint-ledger convention used in the TCH scripts:
    # 8 record generators per vertex, link records, independent plaquette
    # relations, cell closures, one string-tension record, and promotions.
    generator = 8 * complex_.vertices + len(complex_.edges) + rank_d2 + rank_d3 + 1 + promotions
    return {
        "rank_d2": rank_d2,
        "rank_d3": rank_d3,
        "promotions": promotions,
        "bond_dimension_proxy": 2**promotions,
        "generator": generator,
    }


def main() -> None:
    print("Simple-cubic SU(3) static-potential baseline")
    print("=" * 96)
    u = u_fundamental(BETA)
    sigma = -math.log(u)
    print(f"  beta={BETA:.3f}; u_F=beta/18={u:.9f}; sigma_lattice=-ln(u_F)={sigma:.9f}")
    print(f"  rectangular loops use R,T <= {MAX_RECT}; connected-surface table uses area <= {MAX_SURFACE_AREA}.")

    print("\n[1] Cubical 3-ball topology checks")
    for size in (1, 2, 3):
        topology_checks(build_cubical_complex(size))

    print("\n[2] Rectangular Wilson loops and static potential")
    all_sigmas: list[float] = []
    max_degeneracy_by_size: dict[int, int] = {}
    for size in (2, 3):
        rows = rectangle_rows(size)
        print(f"  L={size}: accepted rectangular loops={len(rows)}")
        assert_true(f"L={size} has rectangular Wilson loops", len(rows) > 0)
        max_deg = max(row.degeneracy for row in rows)
        max_degeneracy_by_size[size] = max_deg
        assert_equal(f"L={size} max minimal-surface degeneracy", max_deg, 1)
        by_key: dict[tuple[int, int, int, int], dict[tuple[int, int], RectRow]] = defaultdict(dict)
        for row in rows:
            key = (row.normal_axis, row.plane_position, row.r_axis, row.t_axis)
            by_key[key][(row.r, row.t)] = row
        for key, table in sorted(by_key.items()):
            for r in (1, 2):
                one = table.get((r, 1))
                two = table.get((r, 2))
                if one is None or two is None:
                    continue
                v_r = two.minus_log - one.minus_log
                print(
                    f"    L={size} normal={key[0]} plane={key[1]} Raxis={key[2]} "
                    f"Taxis={key[3]} R={r}: V={v_r:.9f}, sigma*R={sigma * r:.9f}"
                )
                assert_close("static-potential T-decay", v_r, sigma * r, tol=1e-12)
                all_sigmas.append(v_r / r)
    assert_less("orientation sigma spread", max(all_sigmas) - min(all_sigmas), 1e-12)

    print("\n[3] Dual-web ledger scaling")
    print("  L  cells  V    E    F    rank_d2 rank_d3 max_deg promotions bondD generator")
    for size in (1, 2, 3, 4):
        complex_ = build_cubical_complex(size)
        max_deg = max_degeneracy_by_size.get(size, 1)
        row = ledger_row(complex_, max_deg)
        print(
            f"  {size:>1d}"
            f"  {complex_.cells:>5d}"
            f"  {complex_.vertices:>3d}"
            f"  {len(complex_.edges):>4d}"
            f"  {len(complex_.faces):>4d}"
            f"  {row['rank_d2']:>7d}"
            f"  {row['rank_d3']:>7d}"
            f"  {max_deg:>7d}"
            f"  {row['promotions']:>10d}"
            f"  {row['bond_dimension_proxy']:>5d}"
            f"  {row['generator']:>9d}"
        )
        assert_equal(f"L={size} promotion count", row["promotions"], 0)
        assert_equal(f"L={size} bond-dimension proxy", row["bond_dimension_proxy"], 1)
        assert_less(f"L={size} generator ledger / vertex volume", row["generator"] / ((size + 1) ** 3), 40.0)

    print(
        """
VERDICT:
  PASS.  On a licensed simple-cubic 3-ball, the same Wilson-loop extractor gives
  the expected leading strong-coupling static potential V(R)=sigma_lattice R,
  with sigma_lattice=-ln(beta/18), zero observed minimal-surface degeneracy for
  the tested R,T<=2 loops, and no promotions.

  This does not solve the primal TCH fat-block completion problem.  It proves
  the static-potential machinery works on a valid 3D gauge web and keeps open
  the framework fork: colour dynamics may live naturally on the dual
  simple-cubic gauge web, while the closed-cell Gauss/triality statements live
  on the primal TCH shell.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
