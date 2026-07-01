#!/usr/bin/env python3
r"""Edge-loop Creutz gate on the licensed bond-bipyramid flat sheet.

Purpose
-------
``record_grammar_tch_bond_bipyramid_isoperimetry.py`` showed that the
bond-bipyramid slab contains connected offset face planes with flat-sheet
isoperimetry.  The earlier rectangular interior-face family had ``P=2A`` and
therefore could not separate a string-tension area term from a perimeter term.

This script builds the first honest edge-loop Wilson family on the licensed
offset plane.  For coordinate rectangles in an 8 x 8 x 3 slab it verifies

    A_ref(R,T) = 4 R T,
    P(R,T)     = 4(R+T)-2,

so area and perimeter are independent regressors.  It then computes a
full-slab non-increasing flip-sector degeneracy lower bound, ``N_sector``,
and checks the split

    -log W_sector = sigma A_ref - log N_sector.

Finally it applies the Creutz mixed difference.  The bare area term gives
``chi = 4 sigma`` and a pure perimeter term cancels.  The measured sector
entropy contributes its own mixed difference, so the sector Wilson value gives

    chi_sector = 4 sigma - Delta_R Delta_T log N_sector.

Boundary
--------
This is still a finite leading-strong-coupling sector calculation.  It does not
claim the exact global ``N_min`` or a continuum QCD string tension.  Its job is
to prove that the observable family is now well posed: perimeter and area
separate, Creutz cancellation works, and the entropy term is visible as a
measured ledger correction rather than being hidden inside a forced fit.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from record_grammar_tch_bond_bipyramid_bulk import (
    BondComplex,
    build_bond_complex,
    bonds_box,
    xor_masks,
)
from record_grammar_tch_bond_bipyramid_isoperimetry import largest_connected_plane
from record_grammar_tch_bond_bipyramid_static_potential import (
    BETA,
    boundary_stats,
    compact_problem,
    face_sum,
    min_area_degeneracy,
    u_fundamental,
)


SLAB_DIMS = (8, 8, 3)
MAX_RECT = 5
MAX_CREUTZ = 4


@dataclass(frozen=True)
class LoopRow:
    r: int
    t: int
    area_ref: int
    perimeter: int
    sector_degeneracy: int
    log_sector_degeneracy: float
    y_area: float
    y_sector: float


def assert_close(name: str, value: float, target: float, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<78s} value={value:.12g} target={target:.12g} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def assert_equal(name: str, value: int, target: int) -> None:
    print(f"  {name:<78s} value={value} target={target}")
    if value != target:
        raise AssertionError(name)


def assert_true(name: str, value: bool) -> None:
    print(f"  {name:<78s} value={value}")
    if not value:
        raise AssertionError(name)


def plane_face_coordinates(complex_: BondComplex) -> tuple[int, int, list[tuple[int, int, int]]]:
    graph = largest_connected_plane(complex_)
    axes = tuple(axis for axis in (0, 1, 2) if axis != graph.axis)
    coords: list[tuple[int, int, int]] = []
    buckets_a = []
    buckets_b = []
    for face_i in graph.nodes:
        sums = face_sum(complex_.faces[face_i])
        a_bucket = sums[axes[0]] // 6
        b_bucket = sums[axes[1]] // 6
        buckets_a.append(a_bucket)
        buckets_b.append(b_bucket)
        coords.append((a_bucket, b_bucket, face_i))
    min_a = min(buckets_a)
    min_b = min(buckets_b)
    shifted = [(a - min_a, b - min_b, face_i) for a, b, face_i in coords]
    return axes[0], axes[1], shifted


def rectangle_patch(coords: list[tuple[int, int, int]], r: int, t: int) -> tuple[int, ...]:
    return tuple(sorted(face_i for a, b, face_i in coords if 0 <= a < r and 0 <= b < t))


def loop_rows(complex_: BondComplex) -> list[LoopRow]:
    _, _, coords = plane_face_coordinates(complex_)
    sigma = -math.log(u_fundamental(BETA))
    rows: list[LoopRow] = []

    for r in range(1, MAX_RECT + 1):
        for t in range(1, MAX_RECT + 1):
            patch = rectangle_patch(coords, r, t)
            boundary = xor_masks(complex_.face_masks, patch)
            components, even, max_degree, perimeter = boundary_stats(boundary, complex_.edges)
            assert_equal(f"R={r},T={t} area_ref", len(patch), 4 * r * t)
            assert_equal(f"R={r},T={t} additive perimeter", perimeter, 4 * (r + t) - 2)
            assert_true(f"R={r},T={t} boundary is even", even)
            assert_equal(f"R={r},T={t} boundary components", components, 1)
            assert_equal(f"R={r},T={t} boundary max degree", max_degree, 2)

            start, flips = compact_problem(patch, complex_.cell_masks)
            sector = min_area_degeneracy(start, flips, "full-slab-sector")
            assert_equal(f"R={r},T={t} no sector undercut", sector.min_area, len(patch))
            assert_true(f"R={r},T={t} sector search not truncated", not sector.truncated)

            expected_sector = 2 ** ((r - 1) * (t - 1))
            assert_equal(f"R={r},T={t} sector degeneracy pattern", sector.degeneracy, expected_sector)
            log_n = math.log(sector.degeneracy)
            y_area = sigma * len(patch)
            rows.append(
                LoopRow(
                    r=r,
                    t=t,
                    area_ref=len(patch),
                    perimeter=perimeter,
                    sector_degeneracy=sector.degeneracy,
                    log_sector_degeneracy=log_n,
                    y_area=y_area,
                    y_sector=y_area - log_n,
                )
            )
    return rows


def fit_area_perimeter(rows: list[LoopRow], attr: str) -> tuple[float, float, float, float]:
    matrix = np.array([[row.area_ref, row.perimeter, 1.0] for row in rows], dtype=float)
    vector = np.array([getattr(row, attr) for row in rows], dtype=float)
    coeff, *_ = np.linalg.lstsq(matrix, vector, rcond=None)
    residual = vector - matrix @ coeff
    return float(coeff[0]), float(coeff[1]), float(coeff[2]), float(np.max(np.abs(residual)))


def creutz(table: dict[tuple[int, int], LoopRow], attr: str, r: int, t: int) -> float:
    return (
        getattr(table[(r, t)], attr)
        + getattr(table[(r + 1, t + 1)], attr)
        - getattr(table[(r + 1, t)], attr)
        - getattr(table[(r, t + 1)], attr)
    )


def main() -> None:
    print("Edge-loop Creutz gate on the licensed bond-bipyramid flat sheet")
    print("=" * 104)
    sigma = -math.log(u_fundamental(BETA))
    print(f"  slab dimensions={SLAB_DIMS}; beta={BETA:.3f}; sigma=-ln(beta/18)={sigma:.9f}")
    print(f"  loop table: R,T<= {MAX_RECT}; Creutz ratios use R,T<= {MAX_CREUTZ}")

    complex_ = build_bond_complex("bond-slab-8x8x3", bonds_box(*SLAB_DIMS))
    rows = loop_rows(complex_)

    print("\n[1] Edge-loop table")
    print("  R T  A_ref  P  N_sector  logNsec  y_area   y_sector")
    for row in rows:
        print(
            f"  {row.r:>1d} {row.t:>1d}"
            f" {row.area_ref:>6d}"
            f" {row.perimeter:>3d}"
            f" {row.sector_degeneracy:>9d}"
            f" {row.log_sector_degeneracy:>8.3f}"
            f" {row.y_area:>8.3f}"
            f" {row.y_sector:>10.3f}"
        )

    print("\n[2] Area/perimeter regression")
    area_coeff, perim_coeff, const, resid = fit_area_perimeter(rows, "y_area")
    assert_close("bare area fit coefficient", area_coeff, sigma, tol=1e-12)
    assert_close("bare perimeter fit coefficient", perim_coeff, 0.0, tol=1e-12)
    assert_close("bare fit constant", const, 0.0, tol=1e-12)
    assert_close("bare fit max residual", resid, 0.0, tol=1e-10)

    sec_area, sec_perim, sec_const, sec_resid = fit_area_perimeter(rows, "y_sector")
    expected_area = sigma - math.log(2.0) / 4.0
    expected_perim = math.log(2.0) / 4.0
    expected_const = -math.log(2.0) / 2.0
    assert_close("sector fit area coefficient", sec_area, expected_area, tol=1e-12)
    assert_close("sector fit perimeter coefficient", sec_perim, expected_perim, tol=1e-12)
    assert_close("sector fit constant", sec_const, expected_const, tol=1e-12)
    assert_close("sector fit max residual", sec_resid, 0.0, tol=1e-10)

    print("\n[3] Creutz cancellation")
    table = {(row.r, row.t): row for row in rows}
    chi_area_target = 4.0 * sigma
    chi_entropy_target = math.log(2.0)
    chi_sector_target = chi_area_target - chi_entropy_target
    for r in range(1, MAX_CREUTZ + 1):
        for t in range(1, MAX_CREUTZ + 1):
            chi_area = creutz(table, "y_area", r, t)
            chi_entropy = creutz(table, "log_sector_degeneracy", r, t)
            chi_sector = creutz(table, "y_sector", r, t)
            print(
                f"  R={r},T={t}: "
                f"chi_area={chi_area:.9f}, "
                f"chi_logN={chi_entropy:.9f}, "
                f"chi_sector={chi_sector:.9f}"
            )
            assert_close("Creutz bare area term", chi_area, chi_area_target, tol=1e-12)
            assert_close("Creutz sector entropy term", chi_entropy, chi_entropy_target, tol=1e-12)
            assert_close("Creutz sector Wilson term", chi_sector, chi_sector_target, tol=1e-12)

    print(
        """
VERDICT:
  PASS.  The licensed offset-plane edge loops have independently varying area
  and perimeter: A_ref=4RT and P=4(R+T)-2.  The Creutz mixed difference
  therefore becomes a meaningful estimator on this family.  It returns the
  bare leading strong-coupling area value 4 sigma, cancels a perimeter term,
  and separately exposes the measured sector entropy contribution.

  The sector degeneracy follows N_sector=2^((R-1)(T-1)) on the tested table, so
  the sector entropy is area-law: Delta_R Delta_T log N_sector = ln 2.  Since
  N_sector is a lower bound on the exact minimum-surface count N_min, the exact
  count cannot be bounded or perimeter-only.  The literal O(1) promotion /
  bond-dimension reading is refuted for large loops.

  The exact fixed-thickness transfer matrix remains the right next computation,
  but its job is now quantitative: measure the exact entropy density and
  width-dependence of the 2D tensor-network ledger, not rescue boundedness.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
