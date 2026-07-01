#!/usr/bin/env python3
r"""Bond-centred bipyramid TCH fattening: sheets, slabs, and Wilson surfaces.

Purpose
-------
The naive fat block made from corner-stacked truncated cubes fails the 3-ball
topology gate.  That was the wrong fattening.  The current canon's corrected
matter geometry is bond-centred:

    cube centre -- square face -- adjacent cube centre.

Each bond carries an oblate square bipyramid.  Combinatorially this is an
octahedron: 6 vertices, 12 edges, 8 triangular faces, and Q3 face adjacency.
The honest finite complexes therefore follow the bond graph, not a
corner-sharing stack of truncated cubes.

This script constructs finite bond-bipyramid complexes and checks:

1. the local cell is the Q3/[8,4,4] matter-cell alphabet;
2. lines, planar bond sheets, and slabs are valid finite 3-balls;
3. sheets and slabs contain shared interior triangular faces, hence genuine
   Wilson-surface records beyond the boundary shell;
4. the small-surface record ledger has bounded degeneracy/promotions on the
   tested sheets/slabs.

It deliberately stops before claiming a physical QCD string tension.  The
output licenses the corrected "fat" geometry and says what the next
static-potential run must use.
"""

from __future__ import annotations

import math
from collections import defaultdict, deque
from dataclasses import dataclass
from itertools import combinations, product

from record_grammar_tch_su3_loop_geometry import gf2_rank


AXES = (0, 1, 2)
MAX_SURFACE_AREA = 4


Vertex = tuple[int, int, int]
Face = tuple[Vertex, Vertex, Vertex]
Bond = tuple[tuple[int, int, int], int]


@dataclass(frozen=True)
class BondComplex:
    name: str
    bonds: tuple[Bond, ...]
    vertices: tuple[Vertex, ...]
    edges: tuple[tuple[Vertex, Vertex], ...]
    faces: tuple[Face, ...]
    face_masks: tuple[int, ...]
    cell_masks: tuple[int, ...]
    face_counts: dict[Face, int]


@dataclass(frozen=True)
class SurfaceEntry:
    area: int
    count: int
    witness: tuple[int, ...]


def assert_equal(name: str, value: int, target: int) -> None:
    print(f"  {name:<78s} value={value} target={target}")
    if value != target:
        raise AssertionError(name)


def assert_true(name: str, value: bool) -> None:
    print(f"  {name:<78s} value={value}")
    if not value:
        raise AssertionError(name)


def assert_less_equal(name: str, value: int, bound: int) -> None:
    print(f"  {name:<78s} value={value} bound={bound}")
    if value > bound:
        raise AssertionError(name)


def add(left: Vertex, right: Vertex) -> Vertex:
    return tuple(left[i] + right[i] for i in AXES)  # type: ignore[return-value]


def scale2(coord: tuple[int, int, int]) -> Vertex:
    return tuple(2 * coord[i] for i in AXES)  # type: ignore[return-value]


def unit(axis: int, scale: int = 1) -> Vertex:
    out = [0, 0, 0]
    out[axis] = scale
    return tuple(out)  # type: ignore[return-value]


def canonical_edge(left: Vertex, right: Vertex) -> tuple[Vertex, Vertex]:
    if left == right:
        raise ValueError("degenerate edge")
    return (left, right) if left < right else (right, left)


def bond_cell_faces(coord: tuple[int, int, int], axis: int) -> tuple[Face, ...]:
    """Eight triangular faces of the bond-centred oblate square bipyramid.

    Coordinates are doubled so the two apexes sit on integer cube centres and
    the equator vertices sit on the intervening cube face.
    """

    other = [item for item in AXES if item != axis]
    a_axis, b_axis = other
    midpoint = add(scale2(coord), unit(axis))
    equator: list[Vertex] = []
    for a_sign, b_sign in ((1, 1), (-1, 1), (-1, -1), (1, -1)):
        vertex = list(midpoint)
        vertex[a_axis] += a_sign
        vertex[b_axis] += b_sign
        equator.append(tuple(vertex))  # type: ignore[arg-type]

    apexes = (scale2(coord), add(scale2(coord), unit(axis, 2)))
    faces: list[Face] = []
    for apex in apexes:
        for i in range(4):
            faces.append(tuple(sorted((apex, equator[i], equator[(i + 1) % 4]))))  # type: ignore[arg-type]
    return tuple(faces)


def build_bond_complex(name: str, bonds: list[Bond]) -> BondComplex:
    face_index: dict[Face, int] = {}
    face_counts: dict[Face, int] = defaultdict(int)
    faces: list[Face] = []
    cell_masks: list[int] = []

    for coord, axis in bonds:
        mask = 0
        for face in bond_cell_faces(coord, axis):
            face_counts[face] += 1
            if face_counts[face] > 2:
                raise AssertionError(f"non-manifold face shared by >2 cells: {face}")
            if face not in face_index:
                face_index[face] = len(faces)
                faces.append(face)
            mask ^= 1 << face_index[face]
        cell_masks.append(mask)

    vertices = tuple(sorted({vertex for face in faces for vertex in face}))
    edges = tuple(sorted({canonical_edge(face[i], face[j]) for face in faces for i, j in ((0, 1), (1, 2), (0, 2))}))
    edge_index = {edge: i for i, edge in enumerate(edges)}

    face_masks: list[int] = []
    for face in faces:
        mask = 0
        for i, j in ((0, 1), (1, 2), (0, 2)):
            mask ^= 1 << edge_index[canonical_edge(face[i], face[j])]
        face_masks.append(mask)

    return BondComplex(
        name=name,
        bonds=tuple(bonds),
        vertices=vertices,
        edges=edges,
        faces=tuple(faces),
        face_masks=tuple(face_masks),
        cell_masks=tuple(cell_masks),
        face_counts=dict(face_counts),
    )


def bonds_box(nx: int, ny: int, nz: int, axes: tuple[int, ...] = AXES) -> list[Bond]:
    dims = (nx, ny, nz)
    bonds: list[Bond] = []
    for coord in product(*(range(dim) for dim in dims)):
        for axis in axes:
            if coord[axis] + 1 < dims[axis]:
                bonds.append((coord, axis))  # type: ignore[arg-type]
    return bonds


def xor_masks(masks: tuple[int, ...], indices: tuple[int, ...]) -> int:
    out = 0
    for index in indices:
        out ^= masks[index]
    return out


def topology_checks(complex_: BondComplex) -> tuple[int, int]:
    rank_d2 = gf2_rank(list(complex_.face_masks))
    rank_d3 = gf2_rank(list(complex_.cell_masks))
    euler = len(complex_.vertices) - len(complex_.edges) + len(complex_.faces) - len(complex_.cell_masks)
    assert_equal(f"{complex_.name} Euler characteristic", euler, 1)
    assert_equal(f"{complex_.name} rank d3", rank_d3, len(complex_.cell_masks))
    assert_equal(f"{complex_.name} rank d2 = F-C", rank_d2, len(complex_.faces) - len(complex_.cell_masks))
    face_counts_ok = True
    closures_ok = True
    for cell_mask in complex_.cell_masks:
        face_indices = tuple(i for i in range(len(complex_.faces)) if cell_mask & (1 << i))
        face_counts_ok = face_counts_ok and len(face_indices) == 8
        closures_ok = closures_ok and xor_masks(complex_.face_masks, face_indices) == 0
    assert_true(f"{complex_.name} all cells have 8 faces", face_counts_ok)
    assert_true(f"{complex_.name} all cells satisfy d2 d3 = 0", closures_ok)
    return rank_d2, rank_d3


def local_face_adjacency_is_q3() -> None:
    faces = bond_cell_faces((0, 0, 0), 0)
    adj = {i: set() for i in range(len(faces))}
    for i, j in combinations(range(len(faces)), 2):
        if len(set(faces[i]) & set(faces[j])) == 2:
            adj[i].add(j)
            adj[j].add(i)

    assert_equal("local bond-cell face count", len(faces), 8)
    assert_equal("local face-adjacency edge count", sum(len(items) for items in adj.values()) // 2, 12)
    assert_true("local face-adjacency is 3-regular", all(len(items) == 3 for items in adj.values()))

    colour: dict[int, int] = {}
    bipartite = True
    for start in adj:
        if start in colour:
            continue
        colour[start] = 0
        queue: deque[int] = deque([start])
        while queue:
            cur = queue.popleft()
            for nxt in adj[cur]:
                if nxt in colour:
                    bipartite = bipartite and colour[nxt] != colour[cur]
                else:
                    colour[nxt] = 1 - colour[cur]
                    queue.append(nxt)
    assert_true("local face-adjacency bipartite", bipartite)
    assert_equal("local face-adjacency connected vertices", len(colour), 8)

    distances = []
    for start in adj:
        seen = {start: 0}
        queue = deque([start])
        while queue:
            cur = queue.popleft()
            for nxt in adj[cur]:
                if nxt not in seen:
                    seen[nxt] = seen[cur] + 1
                    queue.append(nxt)
        distances.append(max(seen.values()))
    assert_equal("local face-adjacency diameter", max(distances), 3)


def face_adjacency(complex_: BondComplex) -> dict[int, set[int]]:
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
        for left, right in combinations(incident, 2):
            adj[left].add(right)
            adj[right].add(left)
    return adj


def connected_patches(complex_: BondComplex, max_area: int) -> set[tuple[int, ...]]:
    adj = face_adjacency(complex_)
    patches: set[tuple[int, ...]] = set()
    for start in range(len(complex_.faces)):
        queue: deque[tuple[int, ...]] = deque([(start,)])
        while queue:
            patch = queue.popleft()
            patches.add(patch)
            if len(patch) >= max_area:
                continue
            frontier = sorted({nxt for face in patch for nxt in adj[face]} - set(patch))
            for nxt in frontier:
                new_patch = tuple(sorted((*patch, nxt)))
                if new_patch not in patches:
                    queue.append(new_patch)
    return patches


def min_surface_table(complex_: BondComplex, max_area: int) -> dict[int, SurfaceEntry]:
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


def ledger_row(complex_: BondComplex, max_degeneracy: int) -> dict[str, int]:
    rank_d2 = gf2_rank(list(complex_.face_masks))
    rank_d3 = gf2_rank(list(complex_.cell_masks))
    promotions = math.ceil(math.log2(max(1, max_degeneracy)))
    generator = 8 * len(complex_.vertices) + len(complex_.edges) + rank_d2 + rank_d3 + 1 + promotions
    return {
        "rank_d2": rank_d2,
        "rank_d3": rank_d3,
        "interior_faces": sum(1 for face in complex_.faces if complex_.face_counts[face] == 2),
        "boundary_faces": sum(1 for face in complex_.faces if complex_.face_counts[face] == 1),
        "max_degeneracy": max_degeneracy,
        "promotions": promotions,
        "bond_dimension_proxy": 2**promotions,
        "generator": generator,
    }


def surface_checks(complex_: BondComplex) -> int:
    table = min_surface_table(complex_, MAX_SURFACE_AREA)
    max_degeneracy = max(entry.count for entry in table.values())
    interior_faces = [i for i, face in enumerate(complex_.faces) if complex_.face_counts[face] == 2]
    assert_true(f"{complex_.name} has shared interior triangular faces", len(interior_faces) > 0)
    for face_i in interior_faces[: min(5, len(interior_faces))]:
        boundary = complex_.face_masks[face_i]
        entry = table[boundary]
        assert_equal(f"{complex_.name} interior single-face surface area", entry.area, 1)
    assert_less_equal(f"{complex_.name} max small-surface degeneracy", max_degeneracy, 2)
    print(
        f"  {complex_.name}: small Wilson boundaries={len(table)}, "
        f"interior faces={len(interior_faces)}, max degeneracy={max_degeneracy}"
    )
    return max_degeneracy


def main() -> None:
    print("Bond-centred bipyramid TCH fattening")
    print("=" * 96)

    print("\n[1] Local matter-cell alphabet")
    local_face_adjacency_is_q3()
    print("  interpretation: the bond-centred filler cell is the 8-face Q3/[8,4,4] matter-cell alphabet.")

    shapes = [
        ("single-bond-cell", bonds_box(2, 1, 1, axes=(0,))),
        ("bond-line-3", bonds_box(4, 1, 1, axes=(0,))),
        ("bond-sheet-2x2", bonds_box(2, 2, 1, axes=(0, 1))),
        ("bond-sheet-3x3", bonds_box(3, 3, 1, axes=(0, 1))),
        ("bond-slab-2x2x2", bonds_box(2, 2, 2)),
        ("bond-slab-3x2x2", bonds_box(3, 2, 2)),
        ("bond-slab-3x3x2", bonds_box(3, 3, 2)),
    ]

    complexes = [build_bond_complex(name, bonds) for name, bonds in shapes]

    print("\n[2] Bond-following sheets and slabs close as finite 3-balls")
    ranks: dict[str, tuple[int, int]] = {}
    for complex_ in complexes:
        ranks[complex_.name] = topology_checks(complex_)

    print("\n[3] Wilson-surface records on sheets and slabs")
    degeneracy_by_name: dict[str, int] = {}
    for name in ("bond-sheet-3x3", "bond-slab-2x2x2", "bond-slab-3x2x2"):
        complex_ = next(item for item in complexes if item.name == name)
        degeneracy_by_name[name] = surface_checks(complex_)

    print("\n[4] Ledger scaling under the corrected bond geometry")
    print("  complex             cells  V    E    F   intF  rank_d2 rank_d3 max_deg promotions bondD generator")
    for complex_ in complexes:
        max_deg = degeneracy_by_name.get(complex_.name, 1)
        row = ledger_row(complex_, max_deg)
        print(
            f"  {complex_.name:<18s}"
            f" {len(complex_.cell_masks):>5d}"
            f" {len(complex_.vertices):>3d}"
            f" {len(complex_.edges):>4d}"
            f" {len(complex_.faces):>4d}"
            f" {row['interior_faces']:>5d}"
            f" {row['rank_d2']:>8d}"
            f" {row['rank_d3']:>7d}"
            f" {row['max_degeneracy']:>7d}"
            f" {row['promotions']:>10d}"
            f" {row['bond_dimension_proxy']:>5d}"
            f" {row['generator']:>9d}"
        )
        assert_less_equal(f"{complex_.name} promotion count", row["promotions"], 1)
        assert_less_equal(f"{complex_.name} bond-dimension proxy", row["bond_dimension_proxy"], 2)
        assert_less_equal(f"{complex_.name} generator / cell bound", math.ceil(row["generator"] / len(complex_.cell_masks)), 80)

    print(
        """
VERDICT:
  PASS.  The geometrically honest fattening follows the bond-centred
  bipyramid graph, not a corner-sharing stack of truncated cubes.  The local
  cell has the 8 triangular faces and Q3 face-adjacency of the [8,4,4] matter
  alphabet.  Finite bond lines, planar bond sheets, and slabs all pass the
  3-ball topology gate.

  The sheets and slabs contain shared interior triangular faces, so they carry
  genuine Wilson-surface records beyond a boundary shell.  On the tested
  small-surface ledger, degeneracy is at most 2, requiring at most one
  collective promotion and bond-dimension proxy 2.

  This licenses the corrected bond-following geometry for the next Wilson
  calculation.  It does not yet compute a continuum string tension or replace
  the separate dual-cubic static-potential baseline.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
