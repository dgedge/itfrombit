#!/usr/bin/env python3
r"""foundations_tch_linegraph_literal.py

Extract the ACTUAL t{4,3,4} edge connectivity from the vertex coordinates and build the literal line
graph. Result (a correction): the literal L(TCH) is a 15-band, 8-REGULAR structure with Gamma-spectrum
{8, 2^5, -2^9} -- NOT canon 7.2's 4-regular 12-band cuboctahedron {4,2,2,2,0,0,0,-2^5}. The cuboctahedron
is the line graph of the gauge-cell CUBE GRAPH Q3 (the 8 gauge cells around a matter cell), a
per-matter-cell CLUSTER, which canon labels "L(TCH)" but which is a different object from the line graph
of the t{4,3,4} 1-skeleton.

GEOMETRY (item113): truncated cubes on Z^3 with octagonal faces at +-1/2 (truncation xi=(sqrt2-1)/2, all
edges length sqrt2-1); octahedra fill the (Z+1/2)^3 corner gaps. The t{4,3,4} vertex figure is a square
pyramid -> every vertex has degree 5 (verified). So L(TCH) is 8-regular (each edge shares an endpoint
with 4+4=8 others), with 6 vertices/primitive cell -> 15 edges/cell -> 15 bands.

  [1] truncated-cube vertices = perms(+-xi,+-1/2,+-1/2) (24); all edges length sqrt2-1 (36/cube).
  [2] in the assembled honeycomb every vertex has degree 5 (square-pyramid vertex figure), 6 vertices and
      15 edges per primitive cell -> L(TCH) is 8-regular, 15-band.
  [3] the literal 15-band Bloch H(k): Gamma-spectrum {8, 2^5, -2^9} (nine flat bands at -2). This is NOT
      the cuboctahedron {4,2,2,2,0,0,0,-2^5}.
  [4] CONSEQUENCE: canon 7.2's "cuboctahedral 12-band L(TCH)" is the gauge-link CLUSTER (L of the
      8-gauge-cell cube Q3 around one matter cell) -- a legitimate gauge model and the source of the
      12-eigenvalue cuboctahedron spectrum, but a per-cell cluster, NOT the literal line graph of the
      honeycomb. As a CRYSTAL the gauge links are 4-fold shared -> 3 bands/primitive cell (the SC macro
      photon, 7.3). So the chiral-edge question has TWO literal crystals -- the 15-band L(TCH) or the
      3-band macro photon -- not a 12-band cuboctahedron crystal; the 12-band cuboctahedron is the cluster.

Self-asserting; exit 0. Tier: literal geometry + 15-band/8-regular L(TCH) {8,2^5,-2^9} DERIVED from the
coordinates; the cuboctahedron is the per-matter-cell gauge-link cluster, not the literal L(TCH) -- a
geometric correction to the basis of canon 7.2's photon model.
"""
import itertools
from collections import defaultdict, Counter
import numpy as np

XI = (np.sqrt(2)-1)/2.0
LEN = np.sqrt(2)-1


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def tc_vertices():
    vs = []
    for sx in (1, -1):
        for sy in (1, -1):
            for sz in (1, -1):
                vs += [(sx*XI, sy*0.5, sz*0.5), (sx*0.5, sy*XI, sz*0.5), (sx*0.5, sy*0.5, sz*XI)]
    return np.array(vs)


def kf(p):
    return tuple(np.round(p, 6))


def main():
    print("=== Literal t{4,3,4} line graph from the vertex coordinates ===\n")
    TCV = tc_vertices()
    ok(len(TCV) == 24, "truncated cube has 24 vertices = perms(+-xi,+-1/2,+-1/2)")

    R = range(-1, 4)
    verts = {}
    for g in itertools.product(R, repeat=3):
        for v in TCV:
            p = np.array(g, float)+v
            verts.setdefault(kf(p), p)
    P = list(verts.values())
    idx = {kf(p): i for i, p in enumerate(P)}
    edges = set()
    for g in itertools.product(R, repeat=3):
        cvs = [np.array(g, float)+v for v in TCV]
        for i in range(24):
            for j in range(i+1, 24):
                if abs(np.linalg.norm(cvs[i]-cvs[j]) - LEN) < 1e-6:
                    a, b = idx[kf(cvs[i])], idx[kf(cvs[j])]
                    edges.add((min(a, b), max(a, b)))
    edges = sorted(edges)

    # [2] vertex degrees (interior) and per-cell counts
    deg = Counter()
    for a, b in edges:
        deg[a] += 1; deg[b] += 1
    gc = np.array([1, 1, 1], float)
    central_v = [idx[kf(gc+v)] for v in TCV]
    degs = Counter(deg[i] for i in central_v)
    print(f"[2] interior vertex degrees: {dict(degs)} (square-pyramid vertex figure -> degree 5)")
    ok(set(degs) == {5}, "every t{4,3,4} vertex has degree 5 -> L(TCH) is 8-regular")

    # [3] literal 15-band Bloch H(k)
    mid = {e: (np.array(P[e[0]])+np.array(P[e[1]]))/2 for e in edges}
    def cell(e): return tuple(np.floor(mid[e]+1e-6).astype(int))
    orb = {e: kf(np.array(mid[e])-np.array(cell(e))) for e in edges}
    orbits = sorted(set(orb.values())); oid = {o: i for i, o in enumerate(orbits)}
    nb = len(orbits)
    print(f"\n[3] literal L(TCH): {nb} edge-orbits/cell -> {nb} bands")
    ok(nb == 15, "L(TCH) has 15 bands (6 vertices/cell x degree 5 / 2)")
    vinc = defaultdict(list)
    for e in edges:
        vinc[e[0]].append(e); vinc[e[1]].append(e)
    central = [e for e in edges if cell(e) == (1, 1, 1)]
    hops = []
    for e in central:
        for v in e:
            for f in vinc[v]:
                if f != e:
                    hops.append((oid[orb[e]], oid[orb[f]], tuple(np.array(cell(f))-np.array(cell(e)))))

    def Hk(k):
        H = np.zeros((nb, nb), complex)
        for oi, oj, d in hops:
            H[oi, oj] += np.exp(1j*np.dot(k, d))
        return (H+H.conj().T)/2
    wG = np.sort(np.linalg.eigvalsh(Hk(np.zeros(3))))
    print(f"    Gamma spectrum: {np.round(wG,3)}")
    nflat = int(np.sum(np.abs(wG + 2) < 1e-6))
    ok(abs(wG[-1]-8) < 1e-6, "top band = 8 (8-regular line graph; NOT the cuboctahedron's 4)")
    ok(nflat == 9, "nine flat bands at -2 (line-graph flat-band signature)")
    cub = np.array(sorted([4, 2, 2, 2, 0, 0, 0, -2, -2, -2, -2, -2]))
    ok(not (len(wG) == 12 and np.allclose(wG, cub)), "literal Gamma {8,2^5,-2^9} != cuboctahedron {4,2,2,2,0,0,0,-2^5}")

    # [4] the cuboctahedron is the gauge-cell cube graph Q3 (per-matter-cell cluster), not L(TCH)
    print("\n[4] the cuboctahedron is L(Q3) of the 8 gauge cells around a matter cell (a per-cell CLUSTER):")
    cube_V = list(itertools.product((0, 1), repeat=3))
    cube_E = [(i, j) for i in range(8) for j in range(i+1, 8) if sum(a != b for a, b in zip(cube_V[i], cube_V[j])) == 1]
    Lq = np.zeros((12, 12))
    for a in range(12):
        for b in range(a+1, 12):
            if set(cube_E[a]) & set(cube_E[b]):
                Lq[a, b] = Lq[b, a] = 1
    sp = np.round(np.sort(np.linalg.eigvalsh(Lq)), 6)
    ok(np.allclose(sp, cub), f"L(Q3) (gauge-cell cube) = cuboctahedron {sp} = canon 7.2's 12-band CLUSTER")
    print("    => canon's 'cuboctahedral L(TCH)' is this per-matter-cell gauge-link cluster, not the literal")
    print("       line graph of the honeycomb; as a CRYSTAL the gauge links are 4-fold shared -> 3 bands")
    print("       (the SC macro photon, 7.3). The 12 eigenvalues are the cluster spectrum.")

    print("\n[verdict] LITERAL t{4,3,4} EDGE CONNECTIVITY EXTRACTED -- a geometric correction:")
    print("  - the literal L(TCH) is 15-band, 8-REGULAR, Gamma-spectrum {8, 2^5, -2^9} (degree-5 vertices,")
    print("    square-pyramid vertex figure) -- NOT a 12-band cuboctahedron;")
    print("  - the 12-band cuboctahedron is L(Q3) of the 8 gauge cells around a matter cell -- a per-cell")
    print("    CLUSTER (the source of the 12-eigenvalue spectrum), which canon mislabels 'L(TCH)';")
    print("  - as a CRYSTAL the gauge-link cuboctahedron is 4-fold-shared -> 3 bands (the SC macro photon);")
    print("  - so the chiral-edge question has two LITERAL crystals (15-band L(TCH) or 3-band macro photon),")
    print("    not a 12-band cuboctahedron crystal. The cuboctahedron + its Chern is a cluster idealization.")
    print("  The inter-cell assembly is thus geometrically pinned -- and it corrects the 12-band-crystal")
    print("  premise the prior turns carried. exit 0")


if __name__ == "__main__":
    main()
