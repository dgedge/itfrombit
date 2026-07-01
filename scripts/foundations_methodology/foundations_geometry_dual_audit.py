#!/usr/bin/env python3
r"""foundations_geometry_dual_audit.py

DOWNSTREAM AUDIT of the 2026-06-25 geometry correction: does the corrected substrate (oblate-bipyramid
tiling, ANCHOR 1.2) contain the truncated-cube / octagon / 4.8.8 structure that a class of canon results
is built on (the silver-ratio node area 1.4, the 4.8.8 "Rosetta Stone" projections 1.5, the rho-meson
silver-ratio mass ratio 3+2sqrt2, the two-plaquette silver-ratio lock 7.3)?

Method: compute the DUAL of the bipyramid tiling directly. The dual cell of a primal vertex v = the convex
hull of the centres of the cells (bipyramids) meeting at v; a bipyramid centre = its cube-face centre.

RESULT:
  - dual cell of a matter CORNER (Z^3)        = CUBOCTAHEDRON   (12 vertices, 8 triangles + 6 SQUARES)
  - dual cell of a matter APEX  ((Z+1/2)^3)   = OCTAHEDRON      (6 vertices, 8 triangles)
  => dual of the bipyramid tiling = the RECTIFIED cubic honeycomb (cuboctahedra + octahedra).
  => NO octagonal faces anywhere -- in neither the primal (triangles+squares) nor the dual.

CONSEQUENCES:
  (GOOD)  the cuboctahedron of canon 7.2 IS the dual cell of a matter corner -- so 7.2's photon cell is
          DOUBLY grounded (L(Q3) of the 8 faces AND the dual cell). The C=-1 photon stands.
  (FLAG)  the truncated cube / octagon / 4.8.8 vertex figure is FOREIGN to the corrected substrate and its
          dual. Canon results that derive content FROM the octagon/4.8.8/silver-ratio are therefore on a
          geometry that is neither the matter primal nor its gauge dual -- their geometric grounding is IN
          QUESTION and needs re-derivation on the correct (octahedral/cuboctahedral) vertex figures:
            - 1.4 silver-ratio node area A_node = 1/(4 Lambda^2)  [octagon edge a=(sqrt2-1)/Lambda]
            - 1.5 the "z=3 degree-3 4.8.8 vertex" basis of the 8 pi G coefficient + the 1/3 transport
              penalty  (canon already flags this "vertex-figure-level, 3D-bulk open" -- now sharpened:
              the assumed 2D vertex figure is foreign to the substrate)
            - item 48 rho-meson silver-ratio mass ratio |E1|/|E2| = 3+2sqrt2
            - 7.3 two-plaquette geometric lock / silver-ratio suppression delta_S / beta_4
            - the 4.8.8 dispersion-isotropy / dimensional-reduction arguments
  (ROBUST, NOT octagon-based, unaffected): a_0 = hbar c/Lambda scale; the [8,4,4] 8-face alphabet; 3
          colours; z=6 matter coordination; the SC gauge web; f_pi = Lambda/sqrt(4 pi) (sphere projection,
          not octagon); Dashen pion splitting (semicircle pi/2); graviton; neutron colour capacity.

NB the silver-ratio NUMBER may yet re-derive on the correct geometry (sqrt2 is present: the bipyramid
equatorial diagonal = a sqrt2, the cuboctahedron square diagonal = edge x sqrt2) -- it is the ATTRIBUTION
to the 4.8.8 octagon that is broken, not necessarily the value. A focused reconciliation pass is the
follow-up.

Self-asserting; exit 0.
"""
import itertools
import numpy as np
from collections import Counter, defaultdict


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def hull_faces(pts):
    """return (n_vertices, {face_size: count}) for the convex hull of pts (coplanar simplices merged)."""
    from scipy.spatial import ConvexHull
    P = np.array(pts)
    h = ConvexHull(P)
    faces = defaultdict(set)
    for eq, simp in zip(h.equations, h.simplices):
        faces[tuple(np.round(eq, 4))] |= set(simp)
    return len(P), dict(Counter(sorted(len(s) for s in faces.values())))


def main():
    print("=== Dual of the oblate-bipyramid tiling (downstream geometry audit) ===\n")

    # all unit-square faces of the Z^3 cubic lattice, as (centre, set-of-4-corners). bipyramid centre = face centre.
    faces = []
    for nrm in range(3):
        for base in itertools.product(range(-2, 3), repeat=3):
            c = np.array(base, float) + 0.5
            c[nrm] = base[nrm]
            t = [a for a in range(3) if a != nrm]
            cs = []
            for du, dv in [(-.5, -.5), (-.5, .5), (.5, .5), (.5, -.5)]:
                p = c.copy(); p[t[0]] += du; p[t[1]] += dv; cs.append(p)
            faces.append((c, cs))

    def bip_centres_at(v):
        return [tuple(np.round(c, 6)) for c, cs in faces if any(np.allclose(v, x) for x in cs)]

    # dual cell of a matter CORNER (a cube vertex, Z^3)
    corner = np.array([0., 0., 0.])
    dc = sorted(set(bip_centres_at(corner)))
    nv, ff = hull_faces(dc)
    print(f"[corner Z^3]  bipyramids meeting = {len(dc)};  dual cell = {nv} vertices, faces {ff}")
    ok(nv == 12 and ff.get(3) == 8 and ff.get(4) == 6 and 8 not in ff,
       "dual cell of a matter corner = CUBOCTAHEDRON (12 v, 8 tri + 6 squares, NO octagons) = canon 7.2's photon cell")

    # dual cell of a matter APEX (a cube centre, (Z+1/2)^3): the 6 face-centres of its cube
    apex = np.array([0.5, 0.5, 0.5])
    da = sorted({tuple(np.round(c, 6)) for c, cs in faces
                 if np.max(np.abs(c - apex)) < 0.6 and np.min(np.abs(c - apex)) < 1e-9})
    # the 6 faces of the cube [0,1]^3 have centres at apex +/- 0.5 e_i
    da = [tuple(apex + 0.5 * s * np.eye(3)[i]) for i in range(3) for s in (-1, 1)]
    nv2, ff2 = hull_faces(da)
    print(f"[apex (Z+1/2)^3]  bipyramids meeting = {len(da)};  dual cell = {nv2} vertices, faces {ff2}")
    ok(nv2 == 6 and ff2.get(3) == 8 and 8 not in ff2,
       "dual cell of a matter apex = OCTAHEDRON (6 v, 8 tri, NO octagons)")

    print("\n[verdict] dual(bipyramid tiling) = RECTIFIED cubic honeycomb (cuboctahedra + octahedra):")
    print("  - GOOD: the cuboctahedron = the dual cell of a matter corner -> canon 7.2's photon cell is")
    print("    DOUBLY grounded (L(Q3) of the 8 faces AND the dual cell); the C=-1 photon stands;")
    print("  - FLAG: NO octagons appear in the primal OR the dual. The truncated-cube / octagon / 4.8.8")
    print("    structure is FOREIGN to the corrected substrate. Canon results derived FROM the octagon /")
    print("    4.8.8 vertex figure / silver-ratio (1.4 node area; 1.5 8piG z=3 vertex basis; item 48 rho")
    print("    3+2sqrt2; 7.3 two-plaquette silver lock; 4.8.8 dispersion-isotropy) are on a foreign geometry")
    print("    -> geometric grounding IN QUESTION, needs re-derivation on the octahedral/cuboctahedral");
    print("    vertex figures. (The silver number may survive; the octagon attribution does not.)")
    print("  - ROBUST (not octagon-based): a_0 scale; [8,4,4] alphabet; 3 colours; z=6; SC gauge web;")
    print("    f_pi=Lambda/sqrt(4pi) (sphere); Dashen (semicircle); graviton; neutron capacity. exit 0")


if __name__ == "__main__":
    main()
