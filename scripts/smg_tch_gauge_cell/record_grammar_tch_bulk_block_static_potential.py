#!/usr/bin/env python3
r"""Topology gate for a fat LxLxL TCH bulk block.

Purpose
-------
The linear bulk-chain script is a valid 3D TCH ball, but it is too thin to
show a static potential: small boundary loops remain boundary dominated.  The
obvious next move is to fatten the geometry into an L x L x L block, then
measure rectangular Wilson loops and extract V(R).

This script checks the prerequisite before doing that physics.  It asks whether
the naive construction

    "glue every plus octagon to the neighbouring minus octagon"

already defines a valid finite 3-ball for L=2.  The test searches the full
axis-wise cyclic-shift/reversal class of octagon identifications:

    reverse flags: 2^3,
    cyclic shifts: 8^3,
    total conventions: 4096.

Verdict
-------
No convention in this class gives the required cellular ball:

    d2 d3 = 0,  chi = 1,  rank d3 = C,  rank d2 = F-C.

So the fat-block static-potential run is not yet licensed.  The missing object
is a real TCH bulk-completion rule: a multi-axis orientation/recoupling rule,
or additional filler cells, that makes a thick block a valid 3-complex.  This
negative is the useful result.  It prevents a topologically invalid geometry
from being used to claim confinement physics.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import product

from record_grammar_tch_glued_surface_selector import UnionFind, global_vertex, local_faces_by_label
from record_grammar_tch_su3_loop_geometry import build_faces, consecutive_edges, gf2_rank


AXES = (0, 1, 2)
OCTAGON_VERTICES = 8


@dataclass(frozen=True)
class CandidateStats:
    size: int
    euler: int
    rank_d2: int
    rank_d3: int
    d2_d3_zero: bool
    vertices: int
    edges: int
    faces: int
    cells: int

    @property
    def is_ball(self) -> bool:
        return (
            self.euler == 1
            and self.d2_d3_zero
            and self.rank_d3 == self.cells
            and self.rank_d2 == self.faces - self.cells
        )


def assert_equal(name: str, value: int, target: int) -> None:
    print(f"  {name:<78s} value={value} target={target}")
    if value != target:
        raise AssertionError(name)


def assert_true(name: str, value: bool) -> None:
    print(f"  {name:<78s} value={value}")
    if not value:
        raise AssertionError(name)


def assert_false(name: str, value: bool) -> None:
    print(f"  {name:<78s} value={value}")
    if value:
        raise AssertionError(name)


def face_label(axis: int, side: int) -> str:
    return f"O{axis}{side}"


def coord_to_index(coord: tuple[int, int, int], size: int) -> int:
    x, y, z = coord
    return x + size * (y + size * z)


def index_to_coord(index: int, size: int) -> tuple[int, int, int]:
    z, rem = divmod(index, size * size)
    y, x = divmod(rem, size)
    return x, y, z


def neighbor(coord: tuple[int, int, int], axis: int) -> tuple[int, int, int]:
    out = list(coord)
    out[axis] += 1
    return tuple(out)  # type: ignore[return-value]


def canonical_edge(left: int, right: int) -> tuple[int, int]:
    if left == right:
        raise ValueError("degenerate edge")
    return (left, right) if left < right else (right, left)


def xor_selected(face_masks: list[int], cell_mask: int) -> int:
    out = 0
    for face_i, face_mask in enumerate(face_masks):
        if cell_mask & (1 << face_i):
            out ^= face_mask
    return out


def candidate_stats(
    size: int,
    shifts: tuple[int, int, int],
    reverse_flags: tuple[int, int, int],
) -> CandidateStats | None:
    local_faces = build_faces()
    by_label = local_faces_by_label()
    n_cells = size**3

    uf = UnionFind.build(24 * n_cells)
    for cell_index in range(n_cells):
        coord = index_to_coord(cell_index, size)
        for axis in AXES:
            if coord[axis] + 1 >= size:
                continue
            nbr_index = coord_to_index(neighbor(coord, axis), size)
            plus = list(by_label[face_label(axis, 1)].vertices)
            minus = list(by_label[face_label(axis, 0)].vertices)
            if reverse_flags[axis]:
                minus = list(reversed(minus))
            shift = shifts[axis]
            minus = minus[shift:] + minus[:shift]
            for plus_v, minus_v in zip(plus, minus, strict=True):
                uf.union(global_vertex(cell_index, plus_v), global_vertex(nbr_index, minus_v))

    root_to_new: dict[int, int] = {}

    def canon_vertex(old: int) -> int:
        root = uf.find(old)
        if root not in root_to_new:
            root_to_new[root] = len(root_to_new)
        return root_to_new[root]

    faces: list[tuple[str, tuple[int, ...]]] = []
    cell_face_index: dict[tuple[int, str], int] = {}

    # Shared octagons are added once and assigned to both incident cells.
    for cell_index in range(n_cells):
        coord = index_to_coord(cell_index, size)
        for axis in AXES:
            if coord[axis] + 1 >= size:
                continue
            nbr_index = coord_to_index(neighbor(coord, axis), size)
            plus_label = face_label(axis, 1)
            minus_label = face_label(axis, 0)
            verts = tuple(canon_vertex(global_vertex(cell_index, v)) for v in by_label[plus_label].vertices)
            if len(set(verts)) != len(verts):
                return None
            face_i = len(faces)
            faces.append((f"I{axis}", verts))
            cell_face_index[(cell_index, plus_label)] = face_i
            cell_face_index[(nbr_index, minus_label)] = face_i

    for cell_index in range(n_cells):
        for face in local_faces:
            if (cell_index, face.label) in cell_face_index:
                continue
            verts = tuple(canon_vertex(global_vertex(cell_index, v)) for v in face.vertices)
            if len(set(verts)) != len(verts):
                return None
            face_i = len(faces)
            faces.append((face.label, verts))
            cell_face_index[(cell_index, face.label)] = face_i

    try:
        edges = sorted(
            {
                canonical_edge(left, right)
                for _, vertices in faces
                for left, right in consecutive_edges(vertices)
            }
        )
    except ValueError:
        return None

    edge_index = {edge: i for i, edge in enumerate(edges)}
    face_masks: list[int] = []
    for _, vertices in faces:
        mask = 0
        try:
            for left, right in consecutive_edges(vertices):
                mask ^= 1 << edge_index[canonical_edge(left, right)]
        except ValueError:
            return None
        face_masks.append(mask)

    cell_masks: list[int] = []
    for cell_index in range(n_cells):
        mask = 0
        for face in local_faces:
            mask ^= 1 << cell_face_index[(cell_index, face.label)]
        cell_masks.append(mask)

    euler = len(root_to_new) - len(edges) + len(faces) - n_cells
    return CandidateStats(
        size=size,
        euler=euler,
        rank_d2=gf2_rank(face_masks),
        rank_d3=gf2_rank(cell_masks),
        d2_d3_zero=all(xor_selected(face_masks, cell_mask) == 0 for cell_mask in cell_masks),
        vertices=len(root_to_new),
        edges=len(edges),
        faces=len(faces),
        cells=n_cells,
    )


def search_l2_blocks() -> tuple[int, int, Counter[tuple[int, bool, bool, bool]]]:
    tested = 0
    valid = 0
    histogram: Counter[tuple[int, bool, bool, bool]] = Counter()
    for reverse_flags in product((0, 1), repeat=3):
        for shifts in product(range(OCTAGON_VERTICES), repeat=3):
            tested += 1
            stats = candidate_stats(size=2, shifts=shifts, reverse_flags=reverse_flags)  # type: ignore[arg-type]
            if stats is None:
                continue
            valid += 1
            histogram[
                (
                    stats.euler,
                    stats.d2_d3_zero,
                    stats.rank_d2 == stats.faces - stats.cells,
                    stats.rank_d3 == stats.cells,
                )
            ] += 1
            if stats.is_ball:
                raise AssertionError(
                    f"unexpected valid 2x2x2 ball for reverse={reverse_flags}, shifts={shifts}: {stats}"
                )
    return tested, valid, histogram


def main() -> None:
    print("Fat TCH block topology gate")
    print("=" * 88)
    print("  searched gluing class: axis-wise octagon reversal x cyclic shift")
    print("  reverse flags: 2^3; cyclic shifts: 8^3; total conventions: 4096")

    print("\n[1] Single-cell reference ball")
    single = candidate_stats(size=1, shifts=(0, 0, 0), reverse_flags=(1, 1, 1))
    assert_true("single-cell stats exists", single is not None)
    assert single is not None
    assert_equal("single-cell Euler characteristic", single.euler, 1)
    assert_equal("single-cell rank d3", single.rank_d3, single.cells)
    assert_equal("single-cell rank d2 = F-C", single.rank_d2, single.faces - single.cells)
    assert_true("single-cell d2 d3 = 0", single.d2_d3_zero)
    assert_true("single cell is a finite 3-ball", single.is_ball)

    print("\n[2] Naive 2x2x2 all-octagon block search")
    default = candidate_stats(size=2, shifts=(0, 0, 0), reverse_flags=(1, 1, 1))
    assert_true("default 2x2x2 stats exists", default is not None)
    assert default is not None
    print(
        "  default convention: "
        f"chi={default.euler}, V={default.vertices}, E={default.edges}, "
        f"F={default.faces}, C={default.cells}, "
        f"rank_d2={default.rank_d2}, rank_d3={default.rank_d3}, d2d3={default.d2_d3_zero}"
    )
    assert_false("default convention is a valid 3-ball", default.is_ball)

    tested, valid, histogram = search_l2_blocks()
    assert_equal("2x2x2 conventions tested", tested, 2**3 * OCTAGON_VERTICES**3)
    assert_true("some nondegenerate conventions exist", valid > 0)
    print(f"  nondegenerate conventions: {valid}/{tested}")
    print("  leading topology signatures (chi, d2d3, rank_d2=F-C, rank_d3=C):")
    for key, count in histogram.most_common(10):
        print(f"    {key}: {count}")

    assert_true("no searched 2x2x2 convention was a valid 3-ball", True)

    print(
        """
VERDICT:
  PASS AS AN OBSTRUCTION.  The single truncated-cube cell is a valid finite
  3-ball, and the earlier linear glued chain is valid, but a fat 2x2x2 block
  cannot be obtained by naively gluing all neighbouring octagonal faces under
  any axis-wise reversal/shift convention in the 4096-rule class searched here.

  Therefore the rectangular Wilson-loop static-potential calculation is
  deliberately not run on this geometry.  The next theorem is a TCH bulk
  completion rule: either a missing multi-axis orientation/recoupling rule or
  additional filler cells that make a thick block a valid cellular 3-complex.
  Only after that gate passes should one extract V(R)=sigma R or test
  promotion/bond-dimension scaling in a genuine fat 3D block.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
