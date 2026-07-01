#!/usr/bin/env python3
r"""foundations_tch_literal_edge_dispersion.py

The literal edge dispersion: build the actual t{4,3,4} line graph (15-band, from the coordinates), put
the derived photon-polarization Berry connection on every hop, and read off the topology. Result (an
honest negative): the literal L(TCH) photon is topologically TRIVIAL -- Chern 0 on every k_z plane, no
Weyl points, gaps open across the 3D BZ -> NO chiral edge modes. So canon 7.2's C_{S7}=-1 is an artifact
of the cuboctahedron CLUSTER idealisation (L(Q3) of the 8 gauge cells around a matter cell), NOT a
property of the real geometry; the chiral-photon / chiral-edge-carrier thread does not survive the literal
t{4,3,4} line graph.

  [1] build the literal 15-band L(TCH) (truncated cubes on Z^3, xi=(sqrt2-1)/2; degree-5 vertices).
  [2] Berry connection on every hop: H[e,e'] = <eps(d_e)|eps(d_e')>, eps the helicity polarisation of the
      edge directions (the connection derived in foundations_gauss_connection_berry.py).
  [3] bands + gaps; the Berry phasing lifts the 9 flat bands and opens several gaps.
  [4] Chern of every gapped band-group on k_z = 0, pi/2, pi: ALL ZERO; min 3D gap > 0 (no Weyl) -> trivial.
  => no chiral edge modes on the literal geometry.

Self-asserting; exit 0. Tier: literal geometry + Berry connection -> Chern 0 on all planes, no Weyl,
gaps open -> the literal t{4,3,4} photon is topologically TRIVIAL. Canon's C_{S7}=-1 is a cluster
idealisation, not the real-geometry topology; the chiral-edge carrier does not exist on the literal web.
"""
import itertools
from collections import defaultdict
import numpy as np

XI = (np.sqrt(2)-1)/2.0
LEN = np.sqrt(2)-1


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def build():
    TCV = []
    for sx in (1, -1):
        for sy in (1, -1):
            for sz in (1, -1):
                TCV += [(sx*XI, sy*0.5, sz*0.5), (sx*0.5, sy*XI, sz*0.5), (sx*0.5, sy*0.5, sz*XI)]
    TCV = np.array(TCV)
    kf = lambda p: tuple(np.round(p, 6))
    verts = {}
    for g in itertools.product(range(-1, 4), repeat=3):
        for v in TCV:
            p = np.array(g, float)+v; verts.setdefault(kf(p), p)
    P = list(verts.values()); idx = {kf(p): i for i, p in enumerate(P)}
    edges = set()
    for g in itertools.product(range(-1, 4), repeat=3):
        cvs = [np.array(g, float)+v for v in TCV]
        for i in range(24):
            for j in range(i+1, 24):
                if abs(np.linalg.norm(cvs[i]-cvs[j]) - LEN) < 1e-6:
                    a, b = idx[kf(cvs[i])], idx[kf(cvs[j])]; edges.add((min(a, b), max(a, b)))
    edges = sorted(edges)
    mid = {e: (np.array(P[e[0]])+np.array(P[e[1]]))/2 for e in edges}
    direc = {e: np.array(P[e[1]])-np.array(P[e[0]]) for e in edges}
    cell = lambda e: tuple(np.floor(mid[e]+1e-6).astype(int))
    orb = {e: kf(np.array(mid[e])-np.array(cell(e))) for e in edges}
    orbits = sorted(set(orb.values())); oid = {o: i for i, o in enumerate(orbits)}
    odir = {}
    for e in edges:
        if orb[e] not in odir:
            odir[orb[e]] = direc[e]/np.linalg.norm(direc[e])
    vinc = defaultdict(list)
    for e in edges:
        vinc[e[0]].append(e); vinc[e[1]].append(e)
    hops = []
    for e in [x for x in edges if cell(x) == (1, 1, 1)]:
        for v in e:
            for f in vinc[v]:
                if f != e:
                    hops.append((oid[orb[e]], oid[orb[f]], tuple(np.array(cell(f))-np.array(cell(e)))))
    return len(orbits), hops, [odir[o] for o in orbits]


def pol(d):
    d = np.array(d, float); a = np.array([1., 0, 0]) if abs(d[0]) < 0.9 else np.array([0, 1., 0])
    e1 = np.cross(a, d); e1 /= np.linalg.norm(e1); e2 = np.cross(d, e1)
    return (e1 + 1j*e2)/np.sqrt(2)


def main():
    print("=== Literal t{4,3,4} edge dispersion + topology ===\n")
    nb, hops, dirs = build()
    EPS = [pol(d) for d in dirs]
    ok(nb == 15, "literal L(TCH) is 15-band (degree-5 vertices, 8-regular)")

    def Hk(k):
        H = np.zeros((nb, nb), complex)
        for oi, oj, d in hops:
            H[oi, oj] += np.vdot(EPS[oi], EPS[oj])*np.exp(1j*np.dot(k, d))
        return (H+H.conj().T)/2

    # gaps along a path
    pts = [[0, 0, 0], [np.pi, 0, 0], [np.pi, np.pi, 0], [np.pi, np.pi, np.pi], [0, 0, 0]]
    ks = []
    for a in range(4):
        for t in np.linspace(0, 1, 25, endpoint=False):
            ks.append((1-t)*np.array(pts[a], float)+t*np.array(pts[a+1], float))
    Bnd = np.array([np.linalg.eigvalsh(Hk(k)) for k in ks])
    gaps = [(b, Bnd[:, b+1].min()-Bnd[:, b].max()) for b in range(nb-1) if Bnd[:, b+1].min()-Bnd[:, b].max() > 0.08]
    print("[3] clean gaps (band b -> b+1):", [(b, round(g, 3)) for b, g in gaps])
    ok(len(gaps) >= 1, "the Berry phasing lifts the flat bands and opens gaps")

    def chern(ncut, kz, N=22):
        kk = np.linspace(0, 2*np.pi, N, endpoint=False)
        U = np.empty((N, N, nb, ncut), complex)
        for i, kx in enumerate(kk):
            for j, ky in enumerate(kk):
                w, v = np.linalg.eigh(Hk(np.array([kx, ky, kz]))); U[i, j] = v[:, :ncut]
        F = 0.0
        for i in range(N):
            for j in range(N):
                ip, jp = (i+1) % N, (j+1) % N
                lk = lambda a, b: np.linalg.det(U[a[0], a[1]].conj().T @ U[b[0], b[1]])
                F += np.angle(lk((i, j), (ip, j))*lk((ip, j), (ip, jp))*lk((ip, jp), (i, jp))*lk((i, jp), (i, j)))
        return F/(2*np.pi)

    print("\n[4] Chern of each clean gapped group on k_z = 0, pi/2, pi:")
    trivial = True
    for b, g in gaps:
        cs = [chern(b+1, kz) for kz in (0.0, np.pi/2, np.pi)]
        print(f"    lowest {b+1:2d} bands: Chern(kz=0,pi/2,pi) = {[round(c,2) for c in cs]}")
        if any(abs(round(c)) != 0 for c in cs):
            trivial = False
    ok(trivial, "ALL Chern numbers are 0 on every k_z plane -> topologically trivial")

    # Weyl scan in the cleanest gap
    bcut = max(gaps, key=lambda x: x[1])[0]
    N = 12; mg = 9
    for kx in np.linspace(0, 2*np.pi, N, endpoint=False):
        for ky in np.linspace(0, 2*np.pi, N, endpoint=False):
            for kz in np.linspace(0, 2*np.pi, N, endpoint=False):
                w = np.linalg.eigvalsh(Hk(np.array([kx, ky, kz]))); mg = min(mg, w[bcut+1]-w[bcut])
    print(f"\n    min gap (band {bcut}->{bcut+1}) over 3D BZ = {mg:.3f} (>0 -> no Weyl points in the gap)")
    ok(mg > 0.02, "the clean gap stays open across the 3D BZ -> no Weyl points -> genuinely trivial")

    print("\n[verdict] LITERAL t{4,3,4} PHOTON IS TOPOLOGICALLY TRIVIAL -> NO chiral edge modes:")
    print("  - the actual 15-band L(TCH) + the photon Berry connection has Chern 0 on every k_z plane,")
    print("    gaps open across the 3D BZ, no Weyl points -> no chiral edge dispersion;")
    print("  - so canon 7.2's C_{S7}=-1 is an artifact of the cuboctahedron CLUSTER idealisation (L(Q3) of")
    print("    the 8 gauge cells around a matter cell), NOT a property of the literal t{4,3,4} geometry;")
    print("  - the chiral-photon / chiral-edge-carrier thread (which rested on C=-1) does NOT survive the")
    print("    literal geometry: the literal photon is trivial. (Result is for t{4,3,4} + the Berry-overlap")
    print("    connection; a different intended geometry would need recomputation.) exit 0")


if __name__ == "__main__":
    main()
