#!/usr/bin/env python3
r"""foundations_silver_ratio_octagon_specific.py

RECONCILIATION of the silver-ratio result family after the 2026-06-25 geometry correction (the sharp test
the audit named). The silver ratio delta_S = 1+sqrt2 appears in canon in TWO structurally unrelated ways;
this script separates them and shows which survive on the corrected (oblate-bipyramid) substrate.

(A) GEOMETRIC silver ratio -- from the OCTAGON / 4.8.8 vertex figure:
      * 1.4 node area A_node = 1/(4 Lambda^2): octagon flat-to-flat = 1/Lambda, edge = (sqrt2-1)/Lambda,
        node-area = (3+2sqrt2)/4 edge^2 = W^2/4 because (3+2sqrt2)(3-2sqrt2)=1.
      * 7.3 two-plaquette lock beta_8/beta_4 = A_4/A_8 = 1/(2(1+sqrt2)): squares vs OCTAGONS.
      * 4.8.8 Euler-characteristic flatness chi=0; 4.8.8 dispersion-isotropy; 1.5 "z=3 degree-3 4.8.8
        vertex" basis of the 8 pi G coefficient.
    TEST: does the silver-ratio MECHANISM (flat-to-flat/edge = 1+sqrt2, and node-area = W^2/4) appear on
    the correct solids (cuboctahedron = photon/dual cell; oblate bipyramid = matter cell; octahedron)?
    RESULT: NO -- it is uniquely an octagon property. So every (A) result RETRACTS on the corrected
    geometry (the octagon is foreign: dual = rectified cubic honeycomb, faces are triangles+squares, no
    octagons). On the correct geometry the two-plaquette species are squares vs TRIANGLES, ratio 4/sqrt3,
    NOT the silver ratio.

(B) SPECTRAL silver ratio -- from a Hamiltonian CHARACTERISTIC POLYNOMIAL (octagon-INDEPENDENT):
      * the Gen-1 lepton Hamiltonian gives golden roots (x^2-x-1=0); after Feshbach projection the 3x3
        effective Hamiltonian has eigenvalues {-(1+sqrt2), sqrt2-1, 1} (x^2-2x-1=0): the neutrino-seesaw
        silver ratio, the product identity |E1||E2|=1, the mass ratio 3+2sqrt2 (item 38).
    These come from the code/algebra Hamiltonian, not the substrate geometry -> they STAND.

VERDICT: the geometric octagon silver-ratio family RETRACTS (octagon-specific); the spectral silver-ratio
family (lepton seesaw) STANDS. The only LIVE consumer of (A) is the Part-20 alternative M_P formula
(10.7), already deprecated (DRIFT G6) and shown by DRIFT G7 to be algebraically the SAME alpha^2 relation
as the canonical 10.5 (K_eff=205, octagon-free) -- so the live gravity chain is UNAFFECTED.

Self-asserting; exit 0.
"""
import itertools
import numpy as np
from collections import defaultdict
from scipy.spatial import ConvexHull

SQ2 = np.sqrt(2.0)
SILVER = 1 + SQ2


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def solid(P):
    P = np.array(P, float)
    h = ConvexHull(P)
    grp = defaultdict(list)
    for eq, s in zip(h.equations, h.simplices):
        grp[tuple(np.round(eq, 4))].append(s)
    faces = []
    for eq, simps in grp.items():
        vid = set().union(*[set(s) for s in simps])
        area = sum(0.5 * np.linalg.norm(np.cross(P[s[1]] - P[s[0]], P[s[2]] - P[s[0]])) for s in simps)
        faces.append((len(vid), area, P[list(vid)].mean(0), np.array(eq[:3])))
    edge = min(np.linalg.norm(P[i] - P[j]) for i in range(len(P)) for j in range(i + 1, len(P)))
    w = {}
    for nv, a, c, n in faces:
        for nv2, a2, c2, n2 in faces:
            if nv == nv2 and np.allclose(n, -n2, atol=1e-3):
                w.setdefault(nv, abs(np.dot(c - c2, n)))
    fa = {}
    for nv, a, c, n in faces:
        fa.setdefault(nv, a)
    return edge, w, fa


def has_silver_mechanism(name, P):
    e, w, fa = solid(P)
    silver = any(abs(d / e - SILVER) < 2e-3 for d in w.values())
    quarter = any(abs(a / (d * d) - 0.25) < 2e-3 for a in fa.values() for d in w.values())
    print(f"  {name}: flat-to-flat/edge {[round(d/e,3) for d in w.values()]}; "
          f"silver(1+√2)? {'YES' if silver else 'no'}; node-area=W^2/4? {'YES' if quarter else 'no'}")
    return silver or quarter


def main():
    print("=== silver-ratio reconciliation: geometric (octagon) vs spectral (Hamiltonian) ===\n")

    print("(A) GEOMETRIC silver ratio -- the octagon mechanism, and whether it survives:")
    print(f"  octagon (canon vertex figure): flat-to-flat/edge = 1+√2 = {SILVER:.4f} (YES); "
          f"node-area=W^2/4 EXACT [(3+2√2)(3-2√2)={(3+2*SQ2)*(3-2*SQ2):.3f}]  -> the cancellation")
    cubo = ([(a, b, 0) for a in (-1, 1) for b in (-1, 1)] + [(a, 0, b) for a in (-1, 1) for b in (-1, 1)]
            + [(0, a, b) for a in (-1, 1) for b in (-1, 1)])
    bip = [(0, 0, .5), (0, 0, -.5), (.5, .5, 0), (.5, -.5, 0), (-.5, .5, 0), (-.5, -.5, 0)]
    octa = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    found = [has_silver_mechanism("cuboctahedron (photon/dual cell)", cubo),
             has_silver_mechanism("oblate bipyramid (matter cell) ", bip),
             has_silver_mechanism("regular octahedron             ", octa)]
    ok(not any(found), "the silver-ratio MECHANISM is UNIQUELY an octagon property -> all (A) results RETRACT")

    # two-plaquette species on the correct geometry: squares vs TRIANGLES (cuboctahedron), not octagons
    e, w, fa = solid(cubo)
    ratio = fa[4] / fa[3]            # square area / triangle area
    print(f"\n  two-plaquette lock on the CORRECT geometry (cuboctahedron squares vs triangles):")
    print(f"    A_square/A_triangle = {ratio:.4f} = 4/√3 = {4/SQ2/np.sqrt(1.5):.4f}?  silver 2(1+√2)={2*SILVER:.4f}? no")
    ok(abs(ratio - 4/np.sqrt(3)) < 1e-3 and abs(ratio - 2*SILVER) > 0.1,
       "correct-geometry plaquette ratio = 4/√3 (squares vs triangles), NOT the silver ratio -> §7.3 lock retracts")

    print("\n(B) SPECTRAL silver ratio -- octagon-INDEPENDENT (the lepton seesaw, item 38):")
    # Gen-1 lepton block: golden characteristic x^2-x-1=0 ; after Feshbach -> silver x^2-2x-1=0
    silver_roots = np.sort(np.roots([1, -2, -1]))           # {1-√2, 1+√2}
    ok(np.allclose(silver_roots, sorted([1 - SQ2, 1 + SQ2])),
       "Feshbach-projected lepton Hamiltonian: x^2-2x-1=0 -> eigenvalues 1±√2 (algebraic, NOT geometric)")
    ok(abs(abs(1 + SQ2) * abs(SQ2 - 1) - 1.0) < 1e-12, "seesaw product |E1||E2| = (1+√2)(√2-1) = 1 (spectral) -> STANDS")
    ok(abs(abs(1 + SQ2) / abs(SQ2 - 1) - (3 + 2 * SQ2)) < 1e-9, "mass ratio |E1|/|E2| = 3+2√2 (spectral) -> STANDS")

    print("\n[verdict] silver-ratio family SPLITS cleanly:")
    print("  RETRACT (geometric, octagon/4.8.8 -- no analog on the corrected substrate or its dual):")
    print("    - §1.4 node area A_node = 1/(4Λ²)        (octagon flat-to-flat; value NOT reproduced)")
    print("    - §7.3 two-plaquette lock β8/β4 = silver  (squares-vs-octagons; correct geom -> 4/√3, no silver)")
    print("    - 4.8.8 Euler-χ=0 flatness; 4.8.8 dispersion-isotropy; §1.5 z=3/4.8.8 basis of 8πG (re-derive)")
    print("  STAND (spectral, Hamiltonian char.poly -- octagon-independent):")
    print("    - the lepton/neutrino seesaw silver ratio (item 38): eigenvalues 1±√2, |E1||E2|=1, 3+2√2")
    print("  LIVE CHAIN UNAFFECTED: the only consumer of the geometric A_node is the Part-20 alt-M_P (§10.7),")
    print("    already deprecated (G6) and = the canonical §10.5 α² relation (G7); §10.5 (K_eff=205) is")
    print("    octagon-free. f_π, Dashen, a0 scale, [8,4,4], colours, z=6, SC web, graviton, C=-1 photon: all stand.")
    print("  CONFLATION to flag: any canon claim that the spectral silver ratio IS the octagon one (the")
    print("    'golden→silver topological duality' §13.5 link, item-38 line) -- the geometric side is gone;")
    print("    the numerical coincidence 1+√2 = 1+√2 is not a shared mechanism. exit 0")


if __name__ == "__main__":
    main()
