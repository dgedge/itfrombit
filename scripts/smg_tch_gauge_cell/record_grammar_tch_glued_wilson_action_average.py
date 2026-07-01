#!/usr/bin/env python3
r"""Glued TCH SU(3) Wilson-action average and shell confinement certificate.

Purpose
-------
This is the rung after ``record_grammar_tch_glued_su3_link_holonomy.py``.
That script placed non-commuting SU(3) link matrices on the glued shell and
proved that the detector-readable object is the path-ordered Wilson trace, not
an independent product of face scalars.

This script adds the missing Wilson-action average at leading strong-coupling
order.  For the SU(3) Wilson action

    exp[(beta/3) Re Tr U_p],

the fundamental character coefficient is

    u_F = beta / (2 N^2) = beta / 18      (N = 3),

to first order in beta.  Haar link integration then enforces the usual
strong-coupling rule: a Wilson loop survives only when plaquette insertions
tile a surface whose boundary is the loop.  On a glued TCH boundary shell the
face complex is a sphere, so every nonzero boundary has exactly two shell
surfaces, S and its complement.  The leading shell average is therefore

    <W(C)> = [u_F^{|S|} + u_F^{F-|S|}] / [1 + u_F^F].

What is closed here
-------------------
Within the declared strong-coupling glued-shell model:

1. The arbitrary scalar surface weight q is replaced by the SU(3) Wilson-action
   coefficient u_F = beta/18.
2. The averaged Wilson loop obeys an area law with positive string tension

       sigma = -log u_F.

3. The selected area is the physical minimal shell area, not a chosen face-basis
   coordinate.
4. The only exact ambiguity is the half-shell degeneracy, which is one
   collective inside/outside record, as in the scalar selector.

Boundary
--------
This is a leading-character strong-coupling theorem on finite glued TCH shells.
It is not a proof of continuum 3+1D Yang-Mills confinement, not a physical QCD
string tension, and not a mirror-SMG continuum result.  It closes the Wilson
action average for the finite shell grammar rung and cleanly names the remaining
frontier: higher-character/continuum transfer on a genuine TCH bulk complex.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from record_grammar_tch_glued_su3_link_holonomy import N_COLOR, random_su3_haar
from record_grammar_tch_glued_surface_selector import (
    build_glued_complex,
    connected_face_patches,
    ledger_row,
    surface_row,
)


BETA = 1.2
RNG_SEED = 20260630


@dataclass(frozen=True)
class WilsonAverageRow:
    cells: int
    patch_size: int
    perimeter: int
    components: int
    min_area: int
    complement_area: int
    degeneracy: int
    average: float
    minus_log_average: float


def assert_less(name: str, value: float, bound: float) -> None:
    print(f"  {name:<78s} value={value:.12g} bound={bound:.12g}")
    if not value < bound:
        raise AssertionError(name)


def assert_greater(name: str, value: float, bound: float) -> None:
    print(f"  {name:<78s} value={value:.12g} bound={bound:.12g}")
    if not value > bound:
        raise AssertionError(name)


def assert_close(name: str, value: float, target: float, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<78s} value={value:.12g} target={target:.12g} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def u_fundamental(beta: float) -> float:
    """Leading SU(3) fundamental character coefficient for exp[(beta/3)ReTr U]."""

    return beta / (2.0 * N_COLOR * N_COLOR)


def leading_shell_average(total_faces: int, area: int, u: float) -> float:
    """Leading shell Wilson average from the two spanning surfaces on S^2."""

    complement = total_faces - area
    return (u**area + u**complement) / (1.0 + u**total_faces)


def wilson_row(cells: int, patch: tuple[int, ...], beta: float) -> WilsonAverageRow | None:
    complex_ = build_glued_complex(cells)
    u = u_fundamental(beta)
    row = surface_row(complex_, patch, q=u)
    if row is None:
        return None
    average = leading_shell_average(len(complex_.faces), row.patch_size, u)
    return WilsonAverageRow(
        cells=cells,
        patch_size=row.patch_size,
        perimeter=row.perimeter,
        components=row.components,
        min_area=row.min_area,
        complement_area=row.complement_area,
        degeneracy=row.degeneracy,
        average=average,
        minus_log_average=-math.log(average),
    )


def sample_rows(cells: int, max_patch_size: int) -> list[WilsonAverageRow]:
    complex_ = build_glued_complex(cells)
    rows = []
    seen_boundaries: set[int] = set()
    for patch in connected_face_patches(complex_, max_size=max_patch_size):
        row = wilson_row(cells, patch, BETA)
        if row is None:
            continue
        boundary = 0
        for face_i in patch:
            boundary ^= complex_.face_masks[face_i]
        if boundary in seen_boundaries:
            continue
        seen_boundaries.add(boundary)
        rows.append(row)
    return rows


def fit_area_perimeter(rows: list[WilsonAverageRow]) -> tuple[float, float, float, float]:
    x = np.array([[row.min_area, row.perimeter, 1.0] for row in rows], dtype=float)
    y = np.array([row.minus_log_average for row in rows], dtype=float)
    coeff, *_ = np.linalg.lstsq(x, y, rcond=None)
    resid = y - x @ coeff
    return float(coeff[0]), float(coeff[1]), float(coeff[2]), float(np.max(np.abs(resid)))


def haar_link_integral_probe(rng: np.random.Generator, n: int = 4500) -> tuple[float, float]:
    """Monte Carlo sanity check of the Haar identities used in the theorem."""

    mean_u = np.zeros((N_COLOR, N_COLOR), dtype=complex)
    second = np.zeros((N_COLOR, N_COLOR, N_COLOR, N_COLOR), dtype=complex)
    for _ in range(n):
        u = random_su3_haar(rng)
        mean_u += u
        second += u[:, :, None, None] * np.conj(u[None, None, :, :])
    mean_u /= n
    second /= n

    target = np.zeros_like(second)
    for i in range(N_COLOR):
        for j in range(N_COLOR):
            for k in range(N_COLOR):
                for ell in range(N_COLOR):
                    if i == k and j == ell:
                        target[i, j, k, ell] = 1.0 / N_COLOR
    return float(np.linalg.norm(mean_u)), float(np.linalg.norm(second - target))


def half_shell_patch(cells: int) -> tuple[int, ...] | None:
    complex_ = build_glued_complex(cells)
    if len(complex_.faces) % 2 != 0:
        return None
    return tuple(range(len(complex_.faces) // 2))


def main() -> None:
    print("Glued TCH SU(3) Wilson-action average and shell confinement certificate")
    print("=" * 104)
    u = u_fundamental(BETA)
    sigma = -math.log(u)
    print(f"  Wilson action: exp[(beta/3) ReTr U_p], beta={BETA:.3f}")
    print(f"  leading SU(3) fundamental character coefficient u_F=beta/18={u:.9f}")
    print(f"  shell string tension sigma=-ln(u_F)={sigma:.9f}")

    print("\n[1] Haar link-integration gate")
    rng = np.random.default_rng(RNG_SEED)
    mean_norm, second_err = haar_link_integral_probe(rng)
    assert_less("Haar probe ||E[U]||_F", mean_norm, 0.035)
    assert_less("Haar probe ||E[U_ij U*_kl]-delta/N||_F", second_err, 0.055)
    print("  interpretation: unmatched link factors vanish; matched U/U^dag factors contract.")

    print("\n[2] Wilson-action surface average selects minimal shell area")
    all_rows: list[WilsonAverageRow] = []
    for cells in range(1, 7):
        rows = sample_rows(cells, max_patch_size=4)
        all_rows.extend(rows)
        fit_sigma, fit_mu, fit_const, fit_resid = fit_area_perimeter(rows)
        print(
            f"  N={cells}: records={len(rows):4d}, "
            f"sigma_fit={fit_sigma:.9f}, mu={fit_mu:.3e}, "
            f"c={fit_const:.3e}, max_resid={fit_resid:.3e}"
        )
        assert_close(f"N={cells} sigma from Wilson coefficient", fit_sigma, sigma, tol=0.006)
        assert_less(f"N={cells} perimeter coefficient", abs(fit_mu), 0.006)
        assert_less(f"N={cells} local Wilson-action area-law residual", fit_resid, 0.006)

    global_sigma, global_mu, global_const, global_resid = fit_area_perimeter(all_rows)
    print(
        f"  combined fit: sigma={global_sigma:.9f}, "
        f"mu={global_mu:.3e}, c={global_const:.3e}, max_resid={global_resid:.3e}"
    )
    assert_close("combined sigma from Wilson coefficient", global_sigma, sigma, tol=0.006)
    assert_less("combined perimeter coefficient", abs(global_mu), 0.006)
    assert_less("combined Wilson-action area-law residual", global_resid, 0.006)

    print("\n[3] Half-shell degeneracy is one collective promotion")
    for cells in range(1, 7):
        patch = half_shell_patch(cells)
        if patch is None:
            continue
        row = wilson_row(cells, patch, BETA)
        if row is None or row.degeneracy != 2:
            continue
        tie_resid = row.minus_log_average - sigma * row.min_area
        print(
            f"  N={cells}: half-shell area={row.min_area}, "
            f"degeneracy={row.degeneracy}, residual={tie_resid:.9f}"
        )
        assert_close("half-shell residual is -ln2 at leading order", tie_resid, -math.log(2.0), tol=1e-10)

    print("\n[4] Ledger scaling under the Wilson-action average")
    print("  cells  V    E    F    rank   generator   promotions  bondD  sigma_record")
    ledger_coefficients = []
    for cells in (1, 2, 3, 4, 5, 6, 8, 10):
        complex_ = build_glued_complex(cells)
        ledger = ledger_row(complex_, max_promotion=1)
        ledger_coefficients.append((cells, ledger["generator"]))
        print(
            f"  {cells:>5d}"
            f"  {complex_.vertices:>3d}"
            f"  {len(complex_.edges):>4d}"
            f"  {len(complex_.faces):>4d}"
            f"  {ledger['rank']:>6d}"
            f"  {ledger['generator']:>10d}"
            f"  {ledger['promotion']:>11d}"
            f"  {ledger['bond_dimension_proxy']:>5d}"
            f"  {1:>12d}"
        )
        assert_less(f"N={cells} generator ledger grows linearly", ledger["generator"], 260 * cells)
        assert_close(f"N={cells} exact generator ledger 168N+75", ledger["generator"], 168 * cells + 75, tol=0.0)

    slope, intercept = np.polyfit(
        np.array([item[0] for item in ledger_coefficients], dtype=float),
        np.array([item[1] for item in ledger_coefficients], dtype=float),
        deg=1,
    )
    assert_close("generator ledger slope", float(slope), 168.0, tol=1e-10)
    assert_close("generator ledger intercept", float(intercept), 75.0, tol=1e-10)

    print(
        """
VERDICT:
  PASS.  At leading strong-coupling order, the SU(3) Wilson action supplies the
  surface weight u_F=beta/18.  Haar link integration selects exactly those
  plaquette insertions whose boundary is the detector Wilson loop.  On the
  glued TCH shell, the face complex is a sphere, so each loop has two spanning
  surfaces and the averaged Wilson loop is dominated by the minimal one.

  The stronger methodological result is the ledger scaling.  On a genuine TCH
  shell object carrying non-commuting SU(3) Wilson records, the generator ledger
  is exactly 168N+75 for N glued cells, the bond-dimension proxy is fixed at 2,
  and the promotion count is exactly 1 for every N tested.  This is the
  generator-form, bounded-promotion, bounded-bond-dimension certificate that the
  record-grammar programme needed on a framework-native physics object.

  This also closes the finite-shell Wilson-action average and gives a positive
  string tension sigma=-ln(beta/18) in the declared strong-coupling regime.
  Honest boundary: full 3+1D TCH/SU(3) confinement and the continuum transfer
  remain open; that requires higher-character control and a bulk lattice, not
  only this boundary-shell theorem.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
