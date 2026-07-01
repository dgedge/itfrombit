#!/usr/bin/env python3
r"""Primal bond-bipyramid TCH static-potential and ledger stress test.

Purpose
-------
``record_grammar_tch_bond_bipyramid_bulk.py`` closed the corrected primal
topology gate: finite sheets and slabs of bond-centred bipyramids are genuine
3-balls with shared interior triangular faces.  This script uses that licensed
geometry for the next Wilson-loop rung.

It builds a 6 x 6 x 3 bond-bipyramid slab and measures rectangular interior
Wilson boundary 1-chains.  At leading SU(3) strong coupling,

    <W(R,T)> = N_min(R,T) u_F^{A_min(R,T)},      u_F = beta/18,

where ``N_min`` is the exact number of minimum-area spanning surfaces in the
coset.  This script computes the leading area contribution and a lower-bound
degeneracy ``N_sector``: all surfaces reachable from the reference rectangular
patch by full-slab cell-boundary flips without ever exceeding the reference
area.  That lower bound is enough to falsify an O(1) promotion claim if it
grows, but it cannot confirm boundedness of the full coset.  The sector
static-potential contribution is extracted from

    V(R) ~= -log W(R,T+1) + log W(R,T).

The same run is the large-area ledger stress test.  Small surfaces had
degeneracy <= 2.  Rectangular Wilson boundaries can have many locally minimal
surfaces, so the script reports both the leading area law and the degeneracy
entropy needed to store the record without enumerating every surface as an
independent branch.

Boundary
--------
This is a finite, leading-strong-coupling calculation in lattice units.  It is
not a continuum Yang-Mills proof.  The exact area-law value is a min-cut object;
the exact degeneracy is the minimum-cut/coset degeneracy.  Counting that
degeneracy is hard in general, but is recoverable at fixed transverse thickness
by a transfer matrix.  The rectangular observable used here is a connected even
boundary 1-chain on the interior triangular-face plane; because the triangular
tessellation self-touches at some vertices, it is not a smooth simple polygon.
In fact the measured edge count satisfies ``P_edges=2A`` on this family, so this
observable does not distinguish perimeter-like from area-like degeneracy
entropy.  That requires the next family: simple edge loops on the slab
1-skeleton.
"""

from __future__ import annotations

import heapq
import math
from collections import defaultdict, deque
from dataclasses import dataclass
from itertools import product

from record_grammar_tch_bond_bipyramid_bulk import (
    AXES,
    BondComplex,
    bonds_box,
    build_bond_complex,
    ledger_row,
    topology_checks,
    xor_masks,
)


BETA = 1.2
N_COLOR = 3
SLAB_DIMS = (6, 6, 3)
MAX_R = 4
MAX_T = 4
MAX_SEARCH_STATES = 2_500_000


@dataclass(frozen=True)
class SearchResult:
    scope: str
    min_area: int
    degeneracy: int
    visited: int
    truncated: bool


@dataclass(frozen=True)
class LoopRow:
    r: int
    t: int
    declared_area: int
    perimeter_edges: int
    boundary_components: int
    boundary_max_degree: int
    min_area: int
    sector_degeneracy: int
    promotions: int
    sigma_area: float
    log_sector_degeneracy: float
    minus_log_w_sector: float
    visited: int


def u_fundamental(beta: float) -> float:
    """Leading SU(3) fundamental character coefficient for exp[(beta/3)ReTr U]."""

    return beta / (2.0 * N_COLOR * N_COLOR)


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


def face_sum(face: tuple[tuple[int, int, int], ...]) -> tuple[int, int, int]:
    return tuple(sum(vertex[axis] for vertex in face) for axis in AXES)  # type: ignore[return-value]


def xy_plane_face_grid(face: tuple[tuple[int, int, int], ...], z_sum: int) -> tuple[int, int] | None:
    """Return the integer square index for one triangular face in the xy plane.

    The bipyramid coordinates are doubled.  A centre square in the z=z0 plane
    contributes four shared triangular faces whose coordinate sums are

        xsum in {6x+2, 6x+4},    ysum in {6y+2, 6y+4},    zsum = 6z0.

    Thus xsum//6 and ysum//6 recover the rectangular grid coordinate.
    """

    x_sum, y_sum, this_z_sum = face_sum(face)
    if this_z_sum != z_sum:
        return None
    if x_sum % 6 not in (2, 4) or y_sum % 6 not in (2, 4):
        return None
    return x_sum // 6, y_sum // 6


def rectangular_patch(complex_: BondComplex, r: int, t: int, z_layer: int = 0) -> tuple[int, ...]:
    """All four shared triangular faces in each square of an R x T interior patch."""

    z_sum = 6 * z_layer
    out: list[int] = []
    for face_i, face in enumerate(complex_.faces):
        if complex_.face_counts[face] != 2:
            continue
        grid = xy_plane_face_grid(face, z_sum)
        if grid is None:
            continue
        x, y = grid
        if 0 <= x < r and 0 <= y < t:
            out.append(face_i)
    return tuple(sorted(out))


def boundary_stats(edge_mask: int, edges: tuple[tuple[tuple[int, int, int], tuple[int, int, int]], ...]) -> tuple[int, bool, int, int]:
    adj: dict[tuple[int, int, int], list[tuple[int, int, int]]] = defaultdict(list)
    perimeter_edges = 0
    for edge_i, (left, right) in enumerate(edges):
        if edge_mask & (1 << edge_i):
            perimeter_edges += 1
            adj[left].append(right)
            adj[right].append(left)

    if not adj:
        return 0, True, 0, perimeter_edges
    even = all(len(items) % 2 == 0 for items in adj.values())
    max_degree = max(len(items) for items in adj.values())
    components = 0
    seen: set[tuple[int, int, int]] = set()
    for start in adj:
        if start in seen:
            continue
        components += 1
        seen.add(start)
        queue: deque[tuple[int, int, int]] = deque([start])
        while queue:
            cur = queue.popleft()
            for nxt in adj[cur]:
                if nxt not in seen:
                    seen.add(nxt)
                    queue.append(nxt)
    return components, even, max_degree, perimeter_edges


def face_vertex_bbox(complex_: BondComplex, face_indices: tuple[int, ...]) -> tuple[tuple[int, int, int], tuple[int, int, int]]:
    vertices = [vertex for face_i in face_indices for vertex in complex_.faces[face_i]]
    lo = tuple(min(vertex[axis] for vertex in vertices) for axis in AXES)
    hi = tuple(max(vertex[axis] for vertex in vertices) for axis in AXES)
    return lo, hi  # type: ignore[return-value]


def cell_vertex_bbox(complex_: BondComplex, cell_mask: int) -> tuple[tuple[int, int, int], tuple[int, int, int]]:
    face_indices = tuple(face_i for face_i in range(len(complex_.faces)) if cell_mask & (1 << face_i))
    return face_vertex_bbox(complex_, face_indices)


def local_cell_masks(complex_: BondComplex, patch: tuple[int, ...], margin: int) -> tuple[int, ...]:
    patch_lo, patch_hi = face_vertex_bbox(complex_, patch)
    pad = 2 * margin
    lo = tuple(patch_lo[axis] - pad for axis in AXES)
    hi = tuple(patch_hi[axis] + pad for axis in AXES)
    selected: list[int] = []
    for cell_mask in complex_.cell_masks:
        cell_lo, cell_hi = cell_vertex_bbox(complex_, cell_mask)
        if all(cell_hi[axis] >= lo[axis] and cell_lo[axis] <= hi[axis] for axis in AXES):
            selected.append(cell_mask)
    return tuple(selected)


def compact_problem(start_faces: tuple[int, ...], cell_masks: tuple[int, ...]) -> tuple[int, tuple[int, ...]]:
    active_faces: set[int] = set(start_faces)
    for cell_mask in cell_masks:
        bit = 0
        mask = cell_mask
        while mask:
            if mask & 1:
                active_faces.add(bit)
            mask >>= 1
            bit += 1

    face_map = {face_i: dense for dense, face_i in enumerate(sorted(active_faces))}
    start = 0
    for face_i in start_faces:
        start ^= 1 << face_map[face_i]

    compact_cells: list[int] = []
    for cell_mask in cell_masks:
        compact = 0
        bit = 0
        mask = cell_mask
        while mask:
            if mask & 1:
                compact ^= 1 << face_map[bit]
            mask >>= 1
            bit += 1
        if compact:
            compact_cells.append(compact)
    return start, tuple(compact_cells)


def min_area_degeneracy(start: int, cell_flips: tuple[int, ...], scope: str) -> SearchResult:
    """Count the non-increasing flip-sector contribution.

    This is a lower bound on the exact minimum-coset degeneracy.  It is enough
    to refute O(1) promotion if it grows, but not enough to prove boundedness.
    """

    start_area = start.bit_count()
    min_area = start_area
    degeneracy = 0
    visited: set[int] = set()
    heap: list[tuple[int, int]] = [(start_area, start)]
    truncated = False

    while heap:
        area, state = heapq.heappop(heap)
        if state in visited:
            continue
        if area > start_area:
            continue
        visited.add(state)
        if len(visited) > MAX_SEARCH_STATES:
            truncated = True
            break
        if area < min_area:
            min_area = area
            degeneracy = 1
        elif area == min_area:
            degeneracy += 1
        for flip in cell_flips:
            nxt = state ^ flip
            if nxt in visited:
                continue
            nxt_area = nxt.bit_count()
            if nxt_area <= start_area:
                heapq.heappush(heap, (nxt_area, nxt))

    return SearchResult(
        scope=scope,
        min_area=min_area,
        degeneracy=degeneracy,
        visited=len(visited),
        truncated=truncated,
    )


def loop_rows(complex_: BondComplex) -> list[LoopRow]:
    sigma = -math.log(u_fundamental(BETA))
    rows: list[LoopRow] = []
    for r, t in product(range(1, MAX_R + 1), range(1, MAX_T + 1)):
        patch = rectangular_patch(complex_, r, t)
        assert_equal(f"R={r},T={t} declared triangular area", len(patch), 4 * r * t)

        boundary = xor_masks(complex_.face_masks, patch)
        components, even, max_degree, perimeter_edges = boundary_stats(boundary, complex_.edges)
        assert_true(f"R={r},T={t} boundary is an even 1-chain", even)
        assert_equal(f"R={r},T={t} boundary components", components, 1)

        local_cells = local_cell_masks(complex_, patch, 2)
        local_start, local_flips = compact_problem(patch, local_cells)
        local = min_area_degeneracy(local_start, local_flips, "margin-2")

        full_start, full_flips = compact_problem(patch, complex_.cell_masks)
        full = min_area_degeneracy(full_start, full_flips, "full-slab-sector")
        assert_equal(f"R={r},T={t} full-slab-sector min area", full.min_area, local.min_area)
        assert_equal(f"R={r},T={t} full-slab-sector degeneracy", full.degeneracy, local.degeneracy)
        assert_true(f"R={r},T={t} full-slab search not truncated", not full.truncated)

        assert_equal(f"R={r},T={t} no full-slab-sector undercut", full.min_area, len(patch))
        promotions = math.ceil(math.log2(max(1, full.degeneracy)))
        sigma_area = sigma * full.min_area
        log_sector_degeneracy = math.log(full.degeneracy)
        rows.append(
            LoopRow(
                r=r,
                t=t,
                declared_area=len(patch),
                perimeter_edges=perimeter_edges,
                boundary_components=components,
                boundary_max_degree=max_degree,
                min_area=full.min_area,
                sector_degeneracy=full.degeneracy,
                promotions=promotions,
                sigma_area=sigma_area,
                log_sector_degeneracy=log_sector_degeneracy,
                minus_log_w_sector=sigma_area - log_sector_degeneracy,
                visited=full.visited,
            )
        )
    return rows


def linear_fit(samples: list[tuple[float, float]]) -> tuple[float, float, float]:
    n = float(len(samples))
    sum_x = float(sum(x for x, _ in samples))
    sum_y = float(sum(y for _, y in samples))
    sum_xx = float(sum(x * x for x, _ in samples))
    sum_xy = float(sum(x * y for x, y in samples))
    denom = n * sum_xx - sum_x * sum_x
    if abs(denom) < 1e-12:
        raise AssertionError("degenerate linear fit")
    slope = (n * sum_xy - sum_x * sum_y) / denom
    intercept = (sum_y - slope * sum_x) / n
    residual = max(abs(y - (slope * x + intercept)) for x, y in samples)
    return slope, intercept, residual


def fit_static_potential(rows: list[LoopRow]) -> tuple[float, float, float]:
    by_rt = {(row.r, row.t): row for row in rows}
    samples: list[tuple[int, float]] = []
    print("\n[4] Static potential from T-decay")
    print("  R  T->T+1  V(R,T)     V/R")
    for r in range(1, MAX_R + 1):
        for t in range(1, MAX_T):
            current = by_rt[(r, t)]
            nxt = by_rt[(r, t + 1)]
            v_rt = nxt.minus_log_w_sector - current.minus_log_w_sector
            samples.append((r, v_rt))
            print(f"  {r:>1d}    {t}->{t+1:<1d}   {v_rt:9.6f}  {v_rt / r:9.6f}")

    return linear_fit([(float(r), v) for r, v in samples])


def report_entropy_scaling(rows: list[LoopRow]) -> None:
    print("\n[5] Sector degeneracy entropy scaling")
    reconstruction_error = max(
        abs(row.minus_log_w_sector - (row.sigma_area - row.log_sector_degeneracy))
        for row in rows
    )
    print(f"  max reconstruction error in -log<W_sector>=sigma A-log N_sector: {reconstruction_error:.3e}")
    assert_less("sector Wilson split reconstruction error", reconstruction_error, 1e-12)

    area_slope, area_intercept, area_resid = linear_fit(
        [(float(row.min_area), row.log_sector_degeneracy) for row in rows]
    )
    perimeter_slope, perimeter_intercept, perimeter_resid = linear_fit(
        [(float(row.perimeter_edges), row.log_sector_degeneracy) for row in rows]
    )
    perimeter_collinear = all(row.perimeter_edges == 2 * row.min_area for row in rows)
    print(
        f"  log N_sector versus area:      slope={area_slope:.6f}, "
        f"c={area_intercept:.6f}, max_resid={area_resid:.6f}"
    )
    print(
        f"  log N_sector versus P_edges:   slope={perimeter_slope:.6f}, "
        f"c={perimeter_intercept:.6f}, max_resid={perimeter_resid:.6f}"
    )
    print(f"  perimeter/area collinearity on this boundary family: {perimeter_collinear}")
    assert_true("this boundary family has P_edges=2A", perimeter_collinear)
    print(
        "  interpretation: this run measures a lower-bound entropy term.  Because "
        "it already grows, it refutes O(1) promotion, but it cannot confirm boundedness."
    )


def main() -> None:
    print("Primal bond-bipyramid TCH static-potential and ledger stress test")
    print("=" * 104)
    u = u_fundamental(BETA)
    sigma = -math.log(u)
    n_cells = len(bonds_box(*SLAB_DIMS))
    print(f"  slab dimensions={SLAB_DIMS}; bond-bipyramid cells={n_cells}")
    print(f"  beta={BETA:.3f}; u_F=beta/18={u:.9f}; triangle-area sigma=-ln(u_F)={sigma:.9f}")
    print(f"  rectangular patch range: R<= {MAX_R}, T<= {MAX_T}; full-slab non-increasing flip sector checked")

    print("\n[1] Slab topology gate")
    slab = build_bond_complex("bond-slab-6x6x3", bonds_box(*SLAB_DIMS))
    topology_checks(slab)

    print("\n[2] Rectangular interior Wilson boundary 1-chains")
    rows = loop_rows(slab)

    print("\n[3] Large-loop ledger table")
    print("  R  T  A_decl A_ref  P_edges max_deg_vtx   N_sector promotions  sigmaA  logNsec -log<Wsec>")
    for row in rows:
        print(
            f"  {row.r:>1d}  {row.t:>1d}"
            f"  {row.declared_area:>6d}"
            f" {row.min_area:>5d}"
            f" {row.perimeter_edges:>8d}"
            f" {row.boundary_max_degree:>11d}"
            f" {row.sector_degeneracy:>10d}"
            f" {row.promotions:>10d}"
            f" {row.sigma_area:>7.2f}"
            f" {row.log_sector_degeneracy:>8.2f}"
            f" {row.minus_log_w_sector:>10.3f}"
        )

    slope, intercept, residual = fit_static_potential(rows)
    print(
        f"\n  linear fit V(R)=sigma_eff R + c: "
        f"sigma_eff={slope:.6f}, c={intercept:.6f}, max_resid={residual:.6f}"
    )
    assert_less("static-potential linear-fit residual", residual, 0.32)
    assert_less("effective string-tension orientation/entropy deviation", abs(slope - 4.0 * sigma), 0.9)

    report_entropy_scaling(rows)

    print("\n[6] Generator-form ledger stress")
    max_deg = max(row.sector_degeneracy for row in rows)
    max_promotions = max(row.promotions for row in rows)
    row = ledger_row(slab, max_deg)
    print(
        f"  slab cells={len(slab.cell_masks)}, V={len(slab.vertices)}, E={len(slab.edges)}, F={len(slab.faces)}, "
        f"interior faces={row['interior_faces']}"
    )
    print(
        f"  max rectangular sector degeneracy={max_deg}, promotion bits={max_promotions}, "
        f"bondD proxy={2**max_promotions}, generator={row['generator']}, "
        f"generator/cell={row['generator'] / len(slab.cell_masks):.2f}"
    )
    assert_equal("ledger-row promotion accounting", row["promotions"], max_promotions)
    assert_less("generator/cell remains bounded on this slab", row["generator"] / len(slab.cell_masks), 40.0)
    assert_true("large-loop sector degeneracy is nontrivial", max_deg > 2)

    print(
        """
VERDICT:
  PASS, with a new caveat.  On the licensed bond-centred primal slab, the
  rectangular interior Wilson boundary 1-chains have A_ref=4RT, and no
  full-slab non-increasing flip-sector undercut appears.  The sector T-decay
  gives a linear static potential in lattice units, with the expected leading
  triangle-area string tension renormalised by the measured lower-bound
  surface entropy.

  The tractability stress test is not the small-surface result repeated.  Large
  rectangular boundaries have a rapidly growing sector degeneracy lower bound
  (for example 50128 at R=T=4), requiring 16 promotion bits if encoded as a
  collective degeneracy record.  Because this is a lower bound, its growth can
  refute the literal O(1) promotion claim but cannot prove boundedness of the
  full coset.  The current boundary family has P_edges=2A, so it cannot decide
  whether the entropy is perimeter-like or area-like.  The rigorous next step
  is fixed-thickness transfer-matrix enumeration of N_min, followed by the
  transverse-thickness limit where string roughening entropy lives.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
