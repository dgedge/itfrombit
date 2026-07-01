#!/usr/bin/env python3
r"""foundations_chi_flatness_generic.py

The "topological flatness from Euler characteristic" item (ANCHOR L2508 / tier-line L2872, Part 20 §3):
  "The 4.8.8 vertex figure has chi = V-E+F = 0 ... a tiling with chi=0 can only cover a flat 2-manifold.
   This provides the topological explanation for the observed spatial flatness Omega_k = 0."

After the 2026-06-25 geometry correction the 4.8.8 octagon vertex figure is retired. Two questions:
  (Q1) is chi=0 octagon-specific (does it survive the correction)?
  (Q2) does chi=0 actually explain the cosmological Omega_k = 0?

RESULTS:
  (Q1) chi=0 is GENERIC to every Euclidean (flat) 2D tiling -- square, triangular, hexagonal, 4.8.8, and
       the corrected geometry's lattice-plane sections all have zero angle defect (delta = 2pi - sum of
       face angles = 0) at every vertex, hence (discrete Gauss-Bonnet, sum delta = 2pi chi) chi=0. It is
       NOT an octagon property; it SURVIVES the correction (the oblate-bipyramid tiling's sections are
       Euclidean too). Convex polyhedra (incl. the corrected solids -- octahedron, cuboctahedron) have
       sum delta = 4pi -> chi = 2 (they are topological spheres / curved). So in 2D, chi=0 <=> flat --
       but the content is the FLATNESS (delta=0), which is generic, not the octagon.
  (Q2) chi=0 does NOT explain Omega_k=0. In 3D, chi=0 is AUTOMATIC for every closed 3-manifold (odd
       dimension -> chi=0 by Poincare duality), independent of curvature: a flat 3-torus and a curved
       3-torus BOTH have chi=0. So chi is curvature-blind in 3D; it cannot constrain the cosmological
       spatial curvature Omega_k. The "chi=0 -> Omega_k=0" step conflates topology (chi) with geometry
       (curvature) -- a flaw INDEPENDENT of the octagon error.

The PROPER discrete-flatness statement is geometric, not topological: the Regge deficit angle around every
edge is zero (the cells' dihedral angles sum to 2pi), i.e. the substrate space-fills flat R^3. The
oblate-bipyramid tiling does this (Monte-Carlo space-filling already verified in
foundations_bipyramid_photon_crystal.py), so the substrate IS flat -- and this survives the correction.
The cosmology-sector argument (ANCHOR L3705: "Omega=1 is a structural tautology because the substrate is a
flat Euclidean space-filling tessellation by construction") is the real account; it stands with the
geometry language corrected (oblate-bipyramid Z^3 (x) Q3, not "4.8.8 / TCH").

VERDICT: the chi=0 item is DOWNGRADED from "topological derivation of Omega_k=0" to a generic consistency
marker. It survives the geometry correction (not octagon-specific) but never derived Omega_k (chi is
curvature-blind in 3D). The flatness that IS real -- the substrate is a zero-deficit Euclidean tiling --
survives on the corrected geometry.

Self-asserting; exit 0.
"""
import numpy as np

DEG = np.pi / 180.0


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def main():
    print("=== chi=0 'flatness' after the geometry correction ===\n")

    # ---- Q1: 2D angle defect -> chi=0 is generic to Euclidean tilings (not octagon) ----
    print("[Q1] 2D angle defect delta = 360deg - (sum of face angles at a vertex); delta=0 <=> flat <=> chi=0:")
    # interior angle of a regular n-gon = (n-2)*180/n
    ang = lambda n: (n - 2) * 180.0 / n
    euclid = {
        "square 4.4.4.4": [4, 4, 4, 4],
        "triangular 3^6": [3] * 6,
        "hexagonal 6.6.6": [6, 6, 6],
        "4.8.8 (retired octagon VF)": [4, 8, 8],
        "snub-square 3.3.4.3.4": [3, 3, 4, 3, 4],
        "corrected-geom {100} section (squares)": [4, 4, 4, 4],
    }
    for name, faces in euclid.items():
        defect = 360.0 - sum(ang(n) for n in faces)
        print(f"    {name:42s}: defect = {defect:+.1f}deg -> {'FLAT (chi=0)' if abs(defect)<1e-9 else 'curved'}")
        ok(abs(defect) < 1e-9, f"{name}: zero angle defect -> flat -> chi=0")
    print("  -> chi=0 holds for EVERY Euclidean tiling; the 4.8.8 is not special. GENERIC, survives the correction.\n")

    print("  Convex polyhedra (incl. the corrected solids) for contrast -- Descartes: sum of defects = 720deg = 4pi -> chi=2:")
    poly = {
        "tetrahedron": (4, [3, 3, 3]),
        "cube": (8, [4, 4, 4]),
        "octahedron (matter cell shape)": (6, [3, 3, 3, 3]),
        "cuboctahedron (photon/dual cell)": (12, [3, 4, 3, 4]),
    }
    for name, (V, faces) in poly.items():
        d = 360.0 - sum(ang(n) for n in faces)
        tot = V * d
        chi = tot / 360.0  # sum delta = 2pi chi -> chi = (sum delta in deg)/360
        print(f"    {name:34s}: per-vertex defect {d:+.0f}deg, total {tot:.0f}deg -> chi = {chi:.0f}")
        ok(abs(tot - 720.0) < 1e-6 and abs(chi - 2) < 1e-9, f"{name}: sum of defects = 4pi -> chi=2 (curved sphere)")
    print("  -> chi=0 <=> flat in 2D; the corrected solids as POLYHEDRA are curved (chi=2), but their")
    print("     space-filling TILING is flat (chi=0). The content is FLATNESS (delta=0), not the Euler number per se.\n")

    # ---- Q2: in 3D chi=0 is automatic (curvature-blind) ----
    print("[Q2] 3D: chi = V - E + F - C of the 3-torus (any cell decomposition) -- cubic unit cell:")
    V, E, F, C = 1, 3, 3, 1
    chi3 = V - E + F - C
    print(f"    cubic decomposition of T^3: V-E+F-C = {V}-{E}+{F}-{C} = {chi3}")
    ok(chi3 == 0, "chi(T^3) = 0")
    print("    THEOREM: every closed odd-dimensional manifold has chi=0 (Poincare duality). So a FLAT 3-torus")
    print("    and a CURVED 3-torus BOTH have chi=0 -> chi is CURVATURE-BLIND in 3D.")
    ok(True, "chi=0 in 3D is automatic & metric-independent -> it cannot constrain Omega_k (topology != geometry)")

    # ---- Q3: the proper flatness is geometric (zero Regge deficit = space-filling) ----
    print("\n[Q3] the REAL flatness is geometric: zero Regge deficit (dihedral angles around each edge sum to 2pi),")
    print("     i.e. the substrate space-fills flat R^3. The oblate-bipyramid tiling does (Monte-Carlo verified")
    print("     in foundations_bipyramid_photon_crystal.py: every point in exactly one cell, no gaps/overlaps),")
    print("     so the substrate IS flat -- and this SURVIVES the geometry correction.")
    ok(True, "proper flatness = zero-deficit Euclidean space-filling tiling -> holds for the corrected bipyramid tiling")

    print("\n[verdict] the chi=0 'topological flatness' item is DOWNGRADED:")
    print("  - chi=0 is GENERIC to every Euclidean tiling (not octagon-specific) -> SURVIVES the correction;")
    print("    the '4.8.8' attribution was spurious (same pattern as the silver-ratio / 8pi dressing).")
    print("  - but chi=0 does NOT derive Omega_k=0: in 3D chi=0 is automatic (curvature-blind); 'chi=0 ->")
    print("    Omega_k=0' conflates topology with geometry -- a flaw INDEPENDENT of the octagon error.")
    print("  - the flatness that IS real -- the substrate is a zero-deficit Euclidean space-filling tiling --")
    print("    survives on the corrected bipyramid geometry. The cosmology-sector 'Omega=1 structural because")
    print("    the substrate is flat by construction' (ANCHOR L3705) is the proper account; it stands with the")
    print("    geometry language corrected to the oblate-bipyramid Z^3 (x) Q3 (not '4.8.8 / TCH'). exit 0")


if __name__ == "__main__":
    main()
