#!/usr/bin/env python3
r"""Non-primitive / tilted Wilson-loop isotropy and roughening entropy.

Purpose
-------
``record_grammar_tch_bond_bipyramid_offaxis_isotropy.py`` closed the *primitive*
off-axis gate: rhombi built from the three nearest-neighbour macro directions
(1,0), (0,1), (1,-1) (all determinant 1) share one bare string tension.  This
script goes beyond primitive rhombi -- to Wilson loops whose spanning
parallelogram is either non-primitive (|det|>1) or tilted to a
non-nearest-neighbour angle -- and asks the two questions left open by that gate:
does the string tension stay controlled (isotropic), and does the worldsheet
roughening entropy stay controlled (area-law) off the primitive axes?

Families tested (u, v spanning vectors; det = |det(u,v)|):

    axis       (1,0)/(0,1)    det 1   reference; N_sector = 2^((R-1)(T-1))
    tilt-unit  (1,1)/(0,1)    det 1   det-1 but tilted to a non-primitive angle
    diamond    (1,1)/(1,-1)   det 2   non-primitive, 45-degree rotated
    tilt3      (2,1)/(-1,1)   det 3   non-primitive, a third angle

For each loop the licensed engine returns the reference area, the exact minimum
surface area (undercut test), and the non-increasing flip-sector degeneracy
N_sector (a lower bound on the worldsheet count, the same quantity used by the
edge-loop Creutz gate).

Findings asserted here
----------------------
1.  No undercut: the exact minimum area equals the flat reference area
    4*det*R*T for every tilted/non-primitive loop -- no diagonal shortcut.
2.  String tension is isotropic: because the area is exactly 4*det*R*T, the
    Creutz mixed difference gives chi_area = 4*det*sigma, so
    chi_area/(4*det) = sigma = -log(beta/18) for every family.
3.  Roughening entropy is controlled and follows an exact closed form:
        axis      N_sector = 2^((R-1)(T-1))
        tilt-unit N_sector = 2^(max(0,(R-1)(T-2)))
        diamond   N_sector = 2^((2R-1)(T-1))
        tilt3     N_sector = 2^((3R-1)(T-1))
    (validated to L=4 in the investigation; here through R,T<=3).  Every form has
    leading term det*L^2 = (cell count), so the entropy DENSITY
    log2(N_sector)/cells -> 1 bit per cell for all orientations and
    non-primitivities: the bulk roughening entropy is isotropic.  The
    orientation dependence sits entirely in the subleading (perimeter) term.

Boundary
--------
Leading strong coupling; finite slab; R,T<=3; the entropy is the computable
flip-sector lower bound, not the exact global N_min; "isotropic" means across the
tested rational directions, not a continuum SO(2) theorem.  What is closed is the
paper's named-open item: beyond primitive rhombi the bare tension and the bulk
roughening-entropy density both stay controlled and isotropic.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from record_grammar_tch_bond_bipyramid_bulk import build_bond_complex, bonds_box, xor_masks
from record_grammar_tch_bond_bipyramid_static_potential import (
    BETA, u_fundamental, boundary_stats, compact_problem, min_area_degeneracy, local_cell_masks)
from record_grammar_tch_bond_bipyramid_offaxis_isotropy import macro_sheet, choose_offset

SLAB_DIMS = (10, 10, 3)
MAX_RT = 3
MARGIN = 3
SIGMA = -math.log(u_fundamental(BETA))


@dataclass(frozen=True)
class Family:
    name: str
    u: tuple[int, int]
    v: tuple[int, int]
    det: int
    f: object  # callable (R,T) -> log2 sector degeneracy


FAMILIES = (
    Family("axis      (1,0)/(0,1)", (1, 0), (0, 1), 1, lambda R, T: (R - 1) * (T - 1)),
    Family("tilt-unit (1,1)/(0,1)", (1, 1), (0, 1), 1, lambda R, T: max(0, (R - 1) * (T - 2))),
    Family("diamond   (1,1)/(1,-1)", (1, 1), (1, -1), 2, lambda R, T: (2 * R - 1) * (T - 1)),
    Family("tilt3     (2,1)/(-1,1)", (2, 1), (-1, 1), 3, lambda R, T: (3 * R - 1) * (T - 1)),
)


def assert_equal(name, value, target):
    print(f"  {name:<74s} value={value} target={target}")
    if value != target:
        raise AssertionError(name)


def assert_close(name, value, target, tol=1e-12):
    err = abs(value - target)
    print(f"  {name:<74s} value={value:.12g} target={target:.12g} err={err:.2e}")
    if err > tol:
        raise AssertionError(name)


def assert_true(name, value):
    print(f"  {name:<74s} value={value}")
    if not value:
        raise AssertionError(name)


def enclosed_cells(u, v, R, T):
    """Integer macro-cells in the half-open parallelogram spanned by R*u and T*v."""
    Ux, Uy, Vx, Vy = R * u[0], R * u[1], T * v[0], T * v[1]
    detM = Ux * Vy - Uy * Vx
    cx = (0, Ux, Vx, Ux + Vx)
    cy = (0, Uy, Vy, Uy + Vy)
    cells = []
    for x in range(min(cx), max(cx) + 1):
        for y in range(min(cy), max(cy) + 1):
            sN, tN = Vy * x - Vx * y, -Uy * x + Ux * y
            inside = (0 <= sN < detM and 0 <= tN < detM) if detM > 0 else (detM < sN <= 0 and detM < tN <= 0)
            if inside:
                cells.append((x, y))
    return cells, abs(detM)


def measure(complex_, by_macro, complete, fam, R, T):
    cells, detM = enclosed_cells(fam.u, fam.v, R, T)
    assert detM == fam.det * R * T, f"{fam.name} det {detM} != {fam.det*R*T}"
    off = choose_offset(tuple(cells), complete)
    faces = []
    for (x, y) in cells:
        faces.extend(by_macro[(x + off[0], y + off[1])])
    patch = tuple(sorted(faces))
    comps, even, _, perim = boundary_stats(xor_masks(complex_.face_masks, patch), complex_.edges)
    start, flips = compact_problem(patch, local_cell_masks(complex_, patch, MARGIN))
    res = min_area_degeneracy(start, flips, fam.name)
    return dict(ncells=len(cells), area=len(patch), perim=perim, comps=comps, even=even,
                min_area=res.min_area, Nsec=res.degeneracy, trunc=res.truncated)


def main():
    print("Non-primitive / tilted Wilson-loop isotropy + roughening entropy")
    print("=" * 100)
    print(f"  slab={SLAB_DIMS}; beta={BETA}; sigma=-ln(beta/18)={SIGMA:.9f}; R,T<= {MAX_RT}")
    complex_ = build_bond_complex("nonprimitive", bonds_box(*SLAB_DIMS))
    by_macro, complete, _ = macro_sheet(complex_)

    chi_over_det = []
    for fam in FAMILIES:
        print(f"\n[{fam.name}]  det={fam.det}")
        tab = {}
        for R in range(1, MAX_RT + 1):
            for T in range(1, MAX_RT + 1):
                m = measure(complex_, by_macro, complete, fam, R, T)
                tab[(R, T)] = m
                assert_equal(f"  R={R},T={T} area = 4*det*RT", m["area"], 4 * fam.det * R * T)
                assert_true(f"  R={R},T={T} no undercut (min_area==area)", m["min_area"] == m["area"])
                assert_true(f"  R={R},T={T} not truncated", not m["trunc"])
                assert_equal(f"  R={R},T={T} N_sector = 2^f", m["Nsec"], 2 ** fam.f(R, T))
        # connectivity: every loop is a single closed boundary except the thin det-1
        # tilt T=1 rows (diagonal-disconnected); report and assert the rest.
        connected = [(R, T) for (R, T), m in tab.items() if m["comps"] == 1]
        assert_true("  >=7 of 9 loops are single closed boundaries", len(connected) >= 7)
        # tension via Creutz on the (exact) area; area=4*det*RT makes this exact.
        y = {k: SIGMA * m["area"] for k, m in tab.items()}
        chi = y[(1, 1)] + y[(2, 2)] - y[(2, 1)] - y[(1, 2)]
        assert_close("  Creutz chi_area/(4*det) = sigma", chi / (4 * fam.det), SIGMA)
        chi_over_det.append(chi / (4 * fam.det))

    print("\n[tension isotropy across families]")
    assert_close("  spread of chi_area/(4*det)", max(chi_over_det) - min(chi_over_det), 0.0)

    print("\n[bulk roughening-entropy density -> 1 bit/cell (from the closed forms)]")
    for fam in FAMILIES:
        for L in (10, 200):
            dens = fam.f(L, L) / (fam.det * L * L)
            print(f"  {fam.name}  L={L:>3d}  log2N/cells = {dens:.4f}")
        assert_true(f"  {fam.name} density(L=200) > 0.98", fam.f(200, 200) / (fam.det * 200 * 200) > 0.98)
        assert_true(f"  {fam.name} density increasing 10->200",
                    fam.f(200, 200) / (fam.det * 200 ** 2) > fam.f(10, 10) / (fam.det * 100))

    print(
        """
VERDICT:
  PASS.  Beyond primitive rhombi -- for non-primitive loops (det 2, 3) and a
  det-1 loop tilted to a non-primitive angle -- the licensed sheet shows:
    (1) no undercut: the flat surface stays minimal (min area = 4*det*RT);
    (2) an exactly isotropic bare string tension, chi_area/(4*det) = sigma;
    (3) a controlled roughening entropy with exact closed forms whose leading
        term is det*L^2 = cell count, so the entropy density is 1 bit per cell
        for every orientation -- the bulk roughening entropy is isotropic, with
        the orientation dependence confined to the subleading perimeter term.
  This closes the nonprimitive-direction / off-axis roughening-entropy frontier
  at leading strong coupling on the finite sheet.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
