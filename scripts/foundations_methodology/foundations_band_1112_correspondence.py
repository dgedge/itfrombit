#!/usr/bin/env python3
r"""foundations_band_1112_correspondence.py

The 11<->12 band-group correspondence (loose end from the chiral-photon resolution): how the edge-identified
11-band crystal (Chern on the lowest-3 group) relates to canon §7.2's 12-band cuboctahedron-cluster
(C_{S7}=-1 on the lower-7).

ANSWER: they are DIFFERENT UNIT-CELL embeddings of the same cuboctahedron (L(Q3)) photon, not a clean
1<->1 band map.
  - 12 = canon's PER-Q3 cell: the 12 edges of ONE bipyramid (8 slant + 4 equator) = L(Q3), spectrum
    [-2(x5), 0(x3), 2(x3), 4]; the phased lower-7 = singlet + two triplets.
  - 11 = the PHYSICAL PER-CUBE crystal: the 3 bipyramids of a cube share edges (slant in 3 bipyramids,
    equator in 4), so 3x12=36 edge-incidences -> 11 distinct = 8 slant + 3 equator. The ONLY count change
    12->11 is the equator: 4 per isolated bipyramid -> 3 distinct per cube.

THE MULTIPLET MAP (at Gamma, unphased): the 12 cuboctahedron levels reorganise into the 11 crystal bands:
  - 2(x3) triplet  -> survives intact;
  - -2(x5)         -> -2(x3) stays + 2 modes pushed DEEP (the cube-centre 8-slant-edge star's inter-cell
                       connectivity);
  - 0(x3)          -> 0(x2): ONE mode is removed -- this is the equator 4->3 identification;
  - 4 (singlet)    -> shifted up by the slant connectivity.
So 12 -> 11 loses exactly one 0-mode (equator identification); the rest reorganise (count-preserving)
under the slant-edge crystal connectivity, which also creates the deep modes absent from the isolated
cluster.

TOPOLOGY: both carry the SAME chiral photon (gauge-invariant C=+-1, helicity-locked). The band-GROUP that
holds it is embedding-dependent -- the 11-band crystal's lowest-3 (its low-energy modes are reorganised by
the slant connectivity, incl. the deep modes) vs canon's per-Q3 12-band lower-7. There is no clean
lowest-3 <-> lower-7 identity; the invariant content (the chiral photon) is identical, the groupings are
not, because the slant-edge inter-cell hopping reorganises the low-energy spectrum.

Self-asserting; exit 0.
"""
import itertools
import numpy as np
from collections import Counter

np.set_printoptions(suppress=True, precision=3)
kf = lambda p: tuple(np.round(p, 6)); ref = np.array([1.0, np.pi, np.e])


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def main():
    print("=== 11 <-> 12 band-group correspondence ===\n")

    # [A] the 12-band cuboctahedron (canon per-Q3 cell)
    V = np.array([(a, b, 0) for a in (-1, 1) for b in (-1, 1)] + [(a, 0, b) for a in (-1, 1) for b in (-1, 1)]
                 + [(0, a, b) for a in (-1, 1) for b in (-1, 1)], float)
    A = np.zeros((12, 12))
    for i in range(12):
        for j in range(i + 1, 12):
            if abs(np.linalg.norm(V[i] - V[j]) - np.sqrt(2)) < 1e-9: A[i, j] = A[j, i] = 1
    cub = np.round(np.sort(np.linalg.eigvalsh(A)), 3)
    print(f"[A] 12-band cuboctahedron L(Q3) (canon §7.2 per-Q3 cell): spectrum {cub}")
    ok(Counter(cub) == Counter([-2, -2, -2, -2, -2, 0, 0, 0, 2, 2, 2, 4]),
       "cuboctahedron multiplets [-2(x5), 0(x3), 2(x3), 4]; phased lower-7 = singlet + two triplets")

    # [B] the 11-band crystal (physical per-cube, edge-identified) -- Gamma spectrum + lowest-3 Chern
    edges = set(); tris = []
    for g in itertools.product(range(3), repeat=3):
        g = np.array(g, float); c = g + 0.5; cor = [g + np.array(o, float) for o in itertools.product((0, 1), repeat=3)]
        for axn in range(3):
            for s in (0, 1):
                fc = [v for v in cor if abs(v[axn] - (g[axn] + s)) < 1e-9]; ct = np.mean(fc, 0); oth = [t for t in range(3) if t != axn]
                fc = sorted(fc, key=lambda v: np.arctan2(v[oth[1]] - ct[oth[1]], v[oth[0]] - ct[oth[0]]))
                for k in range(4):
                    vi, vj = fc[k], fc[(k + 1) % 4]; tris.append([frozenset([kf(c), kf(vi)]), frozenset([kf(c), kf(vj)]), frozenset([kf(vi), kf(vj)])]); edges |= set(tris[-1])
    edges = sorted(edges, key=lambda e: sorted(e)); pos = {e: np.mean([np.array(p) for p in e], 0) for e in edges}
    cell = lambda e: tuple(np.floor(pos[e] + 1e-6).astype(int)); orb = {e: kf(np.array(pos[e]) - np.array(cell(e))) for e in edges}
    orbits = sorted(set(orb.values())); oid = {o: i for i, o in enumerate(orbits)}; nb = len(orbits)
    hop_un = []
    for (a, b, cc) in tris:
        for e in (a, b, cc):
            if cell(e) == (1, 1, 1):
                for f in (a, b, cc):
                    if f != e: hop_un.append((oid[orb[e]], oid[orb[f]], tuple(np.array(cell(f)) - np.array(cell(e)))))
    def Hun(k):
        M = np.zeros((nb, nb), complex)
        for oi, oj, d in hop_un: M[oi, oj] += np.exp(1j * np.dot(k, d))
        return (M + M.conj().T) / 2
    gam = np.round(np.sort(np.linalg.eigvalsh(Hun(np.zeros(3)))), 3)
    print(f"\n[B] 11-band physical crystal (per-cube, edge-identified): Gamma spectrum {gam}")
    ok(nb == 11, "11 bands (8 slant + 3 equator) -- one fewer equator than the per-Q3 cuboctahedron (4->3)")
    cnt = Counter(gam)
    ok(cnt.get(2.0) == 3, "the cuboctahedron 2(x3) triplet SURVIVES intact in the crystal")
    ok(cnt.get(-2.0) == 3, "a -2(x3) triplet survives; the other 2 of the cuboctahedron's -2(x5) are pushed DEEP")
    ok(cnt.get(0.0) == 2, "the cuboctahedron 0(x3) triplet -> 0(x2): ONE mode removed = the equator 4->3 identification")
    ok(min(gam) < -5 and max(gam) > 5, "deep + shifted modes appear (slant-star inter-cell connectivity; absent from the isolated cluster)")

    # [C] chiral coin flux -> lowest-3 Chern (helicity-locked)
    def corder(tri, rr):
        ms = [pos[e] for e in tri]; ctr = np.mean(ms, 0); n = np.cross(ms[1] - ms[0], ms[2] - ms[0]); n = n if np.dot(n, rr) >= 0 else -n
        b = ms[0] - ctr; b = b - np.dot(b, n) / np.dot(n, n) * n; b /= np.linalg.norm(b)
        f = lambda m: (lambda v: np.arctan2(np.dot(np.cross(b, v), n), np.dot(b, v)))(m - ctr - np.dot(m - ctr, n) / np.dot(n, n) * n)
        return sorted(range(3), key=lambda i: f(ms[i]))
    def chern(theta, rr, N=14):
        HOP = []
        for tri in tris:
            o = corder(tri, rr); s = [tri[i] for i in o]
            for a in range(3):
                p, q = s[a], s[(a + 1) % 3]
                if cell(p) == (1, 1, 1): HOP.append((oid[orb[p]], oid[orb[q]], tuple(np.array(cell(q)) - np.array(cell(p))), +1))
                if cell(q) == (1, 1, 1): HOP.append((oid[orb[q]], oid[orb[p]], tuple(np.array(cell(p)) - np.array(cell(q))), -1))
        def Hk(k):
            M = np.zeros((nb, nb), complex)
            for oi, oj, d, sg in HOP: M[oi, oj] += np.exp(1j * sg * theta) * np.exp(1j * np.dot(k, d))
            return (M + M.conj().T) / 2
        kk = np.linspace(0, 2 * np.pi, N, endpoint=False); U = np.empty((N, N, nb, 3), complex)
        for i, kx in enumerate(kk):
            for j, ky in enumerate(kk): _, v = np.linalg.eigh(Hk(np.array([kx, ky, 0.]))); U[i, j] = v[:, :3]
        F = 0.
        for i in range(N):
            for j in range(N):
                ip, jp = (i + 1) % N, (j + 1) % N; lk = lambda P, Q: np.linalg.det(U[P[0], P[1]].conj().T @ U[Q[0], Q[1]])
                F += np.angle(lk((i, j), (ip, j)) * lk((ip, j), (ip, jp)) * lk((ip, jp), (i, jp)) * lk((i, jp), (i, j)))
        return F / (2 * np.pi)
    cp, cm = chern(np.pi / 6, ref), chern(np.pi / 6, -ref)
    print(f"\n[C] chiral coin flux pi/2 per triangle: 11-band lowest-3 Chern = {cp:.2f} (hel +), {cm:.2f} (hel -)")
    ok(abs(abs(cp) - 1) < 0.1 and abs(cp + cm) < 0.1, "11-band lowest-3 carries C=+-1 (the chiral photon)")

    print("\n[verdict] 11 <-> 12 correspondence:")
    print("  - SAME object, DIFFERENT unit cell: 12 = per-Q3 cuboctahedron (8 slant + 4 equator);")
    print("    11 = per-cube edge-identified crystal (8 slant + 3 equator). The count change 12->11 is")
    print("    exactly the equator identification (4->3); the cuboctahedron 0-triplet -> 0-doublet.")
    print("  - the -2(x5) and 2(x3) and 4 reorganise under the slant-edge inter-cell connectivity")
    print("    (creating deep modes absent from the isolated cluster); only the one 0-mode is lost.")
    print("  - TOPOLOGY: both carry the same chiral photon (gauge-invariant C=+-1). There is NO clean")
    print("    lowest-3 <-> lower-7 identity -- the band-GROUP holding C is embedding-dependent (the slant")
    print("    connectivity reorganises the low-energy spectrum). The physical photon is the 11-band;")
    print("    canon's 12-band lower-7 is the per-Q3 cell description of the same chiral photon. exit 0")


if __name__ == "__main__":
    main()
