#!/usr/bin/env python3
r"""Finite isotropy gate for bond-bipyramid edge-loop Wilson records.

Purpose
-------
The fixed-strip transfer theorem gives an exact edge-loop count on one licensed
flat offset sheet,

    A_min^strip(R,T) = 4RT,
    N_min^strip(R,T) = 2^((R-1)(T-1)).

This script asks whether that finite Wilson observable is direction-dependent
inside the corrected bond-bipyramid bulk.  It builds a symmetric 6 x 6 x 6
bond-bipyramid slab, chooses one central connected offset face plane for each
normal direction, and evaluates the same coordinate-rectangle edge loops with
both in-plane coordinate orderings.

For each of the six readings it verifies:

    A_ref(R,T) = 4RT,
    P(R,T)     = 4(R+T)-2,
    N_min^strip(R,T) = 2^((R-1)(T-1)).

It then compares the Creutz values

    chi_area    = 4 sigma,
    chi_entropy = log 2,
    chi_sector  = 4 sigma - log 2

across the rotated readings.

Boundary
--------
This is a finite, leading-strong-coupling, tight-strip isotropy gate.  It tests
the cubic-axis / in-sheet-rotation symmetry of the licensed bond-bipyramid
edge-loop observable.  It does not prove full continuum SO(3) restoration,
diagonal-loop isotropy, weak-coupling scaling, or the wider-margin roughening
entropy.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from record_grammar_tch_bond_bipyramid_bulk import (
    AXES,
    BondComplex,
    build_bond_complex,
    bonds_box,
    xor_masks,
)
from record_grammar_tch_bond_bipyramid_isoperimetry import (
    PlaneGraph,
    build_plane_adjacency,
    components,
)
from record_grammar_tch_bond_bipyramid_static_potential import (
    BETA,
    boundary_stats,
    face_sum,
    u_fundamental,
)
from record_grammar_tch_bond_bipyramid_transfer_matrix import (
    exact_min_count,
    face_to_cells,
    transfer_terms,
)


SLAB_DIMS = (6, 6, 6)
MAX_RECT = 3
MAX_CREUTZ = 2


@dataclass(frozen=True)
class OrientedPlane:
    normal_axis: int
    key: int
    in_axes: tuple[int, int]
    swap: bool
    coords: tuple[tuple[int, int, int], ...]


@dataclass(frozen=True)
class IsoRow:
    normal_axis: int
    key: int
    in_axes: tuple[int, int]
    swap: bool
    r: int
    t: int
    area: int
    perimeter: int
    exact_count: int
    log_count: float
    y_area: float
    y_sector: float


def assert_equal(name: str, value: int, target: int) -> None:
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


def central_connected_plane(complex_: BondComplex, normal_axis: int) -> PlaneGraph:
    keys = sorted(
        {face_sum(face)[normal_axis] for face in complex_.faces if complex_.face_counts[face] == 2}
    )
    mid = 0.5 * (keys[0] + keys[-1])
    candidates: list[tuple[int, float, int, tuple[int, ...], dict[int, set[int]]]] = []
    for key in keys:
        adj = build_plane_adjacency(complex_, normal_axis, key)
        if not adj:
            continue
        comp = max(components(adj), key=len)
        comp_adj = {node: {nxt for nxt in adj[node] if nxt in comp} for node in comp}
        candidates.append((len(comp), -abs(key - mid), key, comp, comp_adj))
    if not candidates:
        raise AssertionError(f"no plane for normal axis {normal_axis}")
    _, _, key, comp, comp_adj = max(candidates)
    return PlaneGraph(axis=normal_axis, key=key, nodes=comp, adj=comp_adj, complex_=complex_)


def plane_coords(complex_: BondComplex, graph: PlaneGraph, swap: bool) -> OrientedPlane:
    in_axes = tuple(axis for axis in AXES if axis != graph.axis)
    if swap:
        in_axes = (in_axes[1], in_axes[0])
    raw: list[tuple[int, int, int]] = []
    for face_i in graph.nodes:
        sums = face_sum(complex_.faces[face_i])
        raw.append((sums[in_axes[0]] // 6, sums[in_axes[1]] // 6, face_i))
    min_a = min(a for a, _, _ in raw)
    min_b = min(b for _, b, _ in raw)
    shifted = tuple(sorted((a - min_a, b - min_b, face_i) for a, b, face_i in raw))
    return OrientedPlane(
        normal_axis=graph.axis,
        key=graph.key,
        in_axes=in_axes,  # type: ignore[arg-type]
        swap=swap,
        coords=shifted,
    )


def rectangle_patch(plane: OrientedPlane, r: int, t: int) -> tuple[int, ...]:
    return tuple(sorted(face_i for a, b, face_i in plane.coords if 0 <= a < r and 0 <= b < t))


def oriented_rows(complex_: BondComplex, plane: OrientedPlane) -> list[IsoRow]:
    incident = face_to_cells(complex_)
    sigma = -math.log(u_fundamental(BETA))
    rows: list[IsoRow] = []
    label = f"normal={plane.normal_axis},key={plane.key},axes={plane.in_axes},swap={plane.swap}"
    for r in range(1, MAX_RECT + 1):
        for t in range(1, MAX_RECT + 1):
            patch = rectangle_patch(plane, r, t)
            boundary = xor_masks(complex_.face_masks, patch)
            components_count, even, max_degree, perimeter = boundary_stats(boundary, complex_.edges)
            expected_area = 4 * r * t
            expected_perimeter = 4 * (r + t) - 2
            expected_count = 2 ** ((r - 1) * (t - 1))
            assert_equal(f"{label} R={r},T={t} area", len(patch), expected_area)
            assert_equal(f"{label} R={r},T={t} perimeter", perimeter, expected_perimeter)
            assert_true(f"{label} R={r},T={t} boundary is even", even)
            assert_equal(f"{label} R={r},T={t} boundary components", components_count, 1)
            assert_equal(f"{label} R={r},T={t} boundary max degree", max_degree, 2)

            n_variables, terms, constant = transfer_terms(complex_, incident, patch)
            min_energy, exact_count, _ = exact_min_count(n_variables, terms)
            min_energy += constant
            assert_equal(f"{label} R={r},T={t} exact strip area", min_energy, expected_area)
            assert_equal(f"{label} R={r},T={t} exact strip count", exact_count, expected_count)
            log_count = math.log(exact_count)
            y_area = sigma * expected_area
            rows.append(
                IsoRow(
                    normal_axis=plane.normal_axis,
                    key=plane.key,
                    in_axes=plane.in_axes,
                    swap=plane.swap,
                    r=r,
                    t=t,
                    area=expected_area,
                    perimeter=perimeter,
                    exact_count=exact_count,
                    log_count=log_count,
                    y_area=y_area,
                    y_sector=y_area - log_count,
                )
            )
    return rows


def creutz(table: dict[tuple[int, int], IsoRow], attr: str, r: int, t: int) -> float:
    return (
        getattr(table[(r, t)], attr)
        + getattr(table[(r + 1, t + 1)], attr)
        - getattr(table[(r + 1, t)], attr)
        - getattr(table[(r, t + 1)], attr)
    )


def main() -> None:
    print("Bond-bipyramid edge-loop isotropy gate")
    print("=" * 100)
    sigma = -math.log(u_fundamental(BETA))
    print(f"  slab dimensions={SLAB_DIMS}; beta={BETA:.3f}; sigma=-ln(beta/18)={sigma:.9f}")
    print(f"  exact rectangle table R,T<= {MAX_RECT}; Creutz ratios use R,T<= {MAX_CREUTZ}")

    complex_ = build_bond_complex("bond-slab-6x6x6", bonds_box(*SLAB_DIMS))
    oriented_planes: list[OrientedPlane] = []
    for normal_axis in AXES:
        graph = central_connected_plane(complex_, normal_axis)
        for swap in (False, True):
            oriented_planes.append(plane_coords(complex_, graph, swap))

    print("\n[1] Rotated exact edge-loop readings")
    print("  normal key in_axes swap  chi_area      chi_logN      chi_sector")
    chi_area_values: list[float] = []
    chi_entropy_values: list[float] = []
    chi_sector_values: list[float] = []
    expected_area_chi = 4.0 * sigma
    expected_entropy_chi = math.log(2.0)
    expected_sector_chi = expected_area_chi - expected_entropy_chi
    for plane in oriented_planes:
        rows = oriented_rows(complex_, plane)
        table = {(row.r, row.t): row for row in rows}
        for r in range(1, MAX_CREUTZ + 1):
            for t in range(1, MAX_CREUTZ + 1):
                chi_area = creutz(table, "y_area", r, t)
                chi_entropy = creutz(table, "log_count", r, t)
                chi_sector = creutz(table, "y_sector", r, t)
                assert_close("rotated bare Creutz value", chi_area, expected_area_chi)
                assert_close("rotated entropy Creutz value", chi_entropy, expected_entropy_chi)
                assert_close("rotated sector Creutz value", chi_sector, expected_sector_chi)
                chi_area_values.append(chi_area)
                chi_entropy_values.append(chi_entropy)
                chi_sector_values.append(chi_sector)
        print(
            f"  {plane.normal_axis:>6d}"
            f" {plane.key:>3d}"
            f" {plane.in_axes!s:>7s}"
            f" {str(plane.swap):>4s}"
            f"  {expected_area_chi:>11.9f}"
            f"  {expected_entropy_chi:>11.9f}"
            f"  {expected_sector_chi:>11.9f}"
        )

    print("\n[2] Isotropy spreads")
    area_spread = max(chi_area_values) - min(chi_area_values)
    entropy_spread = max(chi_entropy_values) - min(chi_entropy_values)
    sector_spread = max(chi_sector_values) - min(chi_sector_values)
    assert_close("bare string-tension Creutz spread", area_spread, 0.0, tol=1e-12)
    assert_close("entropy Creutz spread", entropy_spread, 0.0, tol=1e-12)
    assert_close("sector Wilson Creutz spread", sector_spread, 0.0, tol=1e-12)

    print(
        """
VERDICT:
  PASS.  On the symmetric bond-bipyramid slab, the licensed edge-loop Wilson
  observable is isotropic under the tested cubic-axis rotations and in-plane
  axis swaps.  Each rotated reading has A_ref=4RT, P=4(R+T)-2, and exact
  fixed-strip N_min=2^((R-1)(T-1)).

  The bare leading strong-coupling Creutz value, the entropy correction, and
  the entropy-dressed sector Creutz value all have zero spread to numerical
  precision:

      chi_area    = 4 sigma,
      chi_entropy = log 2,
      chi_sector  = 4 sigma - log 2.

  This is a finite cubic-axis isotropy gate.  It does not prove continuum
  SO(3), diagonal-loop isotropy, weak-coupling scaling, or thick-string
  roughening isotropy.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
