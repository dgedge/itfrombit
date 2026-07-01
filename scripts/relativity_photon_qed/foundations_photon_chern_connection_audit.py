#!/usr/bin/env python3
r"""foundations_photon_chern_connection_audit.py

CONNECTION AUDIT of the 2026-06-25 "chiral photon C=-1 vindicated" result
(foundations_bipyramid_photon_crystal.py). Re-examined while doing the Tier-1 photon/topology cluster.

FINDING (an honest self-correction): the C=-1 on the corrected-geometry photon crystal is
CONNECTION-DEPENDENT and was an artifact of the convention I used.

  - I built the crystal with hop amplitude = the full complex overlap <eps(d_i)|eps(d_j)> of the helicity
    polarisations. That overlap has a NON-UNIFORM MAGNITUDE (|<eps|eps>| ranges ~0.21-0.79 on the hops).
    A non-uniform hopping magnitude is NOT a U(1) gauge field -- a photon link is a PURE PHASE (Peierls /
    Berry phase), uniform magnitude. canon 7.2 itself specifies "phi = +/-pi/4 phase" (pure phase).
  - With the overlap convention: C(lowest-3) = -1 (the reported result).
  - With the physically-correct PURE-PHASE connection e^{i arg<eps|eps>} (same phase, magnitude 1): C = 0
    -- the photon is TOPOLOGICALLY TRIVIAL. Normalising the magnitudes to 1 (keeping only the phase) kills
    the Chern, so the -1 came entirely from the magnitude pattern, not from any flux.

So the chiral photon (C=-1) is NOT established on the corrected geometry under canon 7.2's own
(scalar, pure-phase) model. canon 7.2's C_{S7}=-1 was computed on the OLD "L(TCH)" framing (the
truncated-cube geometry retired on 2026-06-25); on the corrected oblate-bipyramid crystal a pure-phase
connection gives a trivial photon. The chiral photon is therefore RE-OPENED (status: unproven), not
vindicated.

WHAT STILL STANDS (unaffected by this):
  - the geometry correction (oblate-bipyramid tiling) and the dual = rectified cubic honeycomb;
  - the cuboctahedron = L(Q3) of the 8 faces = the dual cell (the photon CELL / graph is right -- only its
    Chern is in question);
  - the cuboctahedron unphased graph spectrum [-2x5, 0x3, 2x3, 4] (matches canon 7.2 exactly);
  - the MACRO photon's emergent Lorentz invariance (Part B below) -- this is the real massless photon and
    is independent of the micro-cuboctahedron topology.

*** RESOLVED 2026-06-25 -- see foundations_chiral_photon_dynamical_flux.py ***
This audit correctly found the overlap C=-1 is an artifact and the GEOMETRIC pure-phase connection is
C=0. The resolution: the geometric phase is zero because the triangular-plaquette edge-directions are
COPLANAR (no solid angle), so canon's +-pi/4 is a DYNAMICAL coin / C4v turn-rule phase, not geometric.
A pure-phase, uniform-magnitude (U(1)) chiral coin flux on the triangular plaquettes gives a ROBUST
C=+-1 chiral photon (lowest-3, every k_z, no Weyl, window phi in [0.9, pi/2], helicity-locked). So the
chiral photon IS real -- under the dynamical (not geometric) pure-phase connection. This script's "C=0"
is the geometric-only connection; it is not the full story.

Self-asserting on the FINDINGS (not on a C=-1); exit 0.
"""
import itertools
import numpy as np

np.set_printoptions(suppress=True)


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def pol(d):
    d = np.array(d, float); d /= np.linalg.norm(d)
    x = np.array([1., 0, 0]) if abs(d[0]) < 0.9 else np.array([0, 1., 0])
    e1 = np.cross(x, d); e1 /= np.linalg.norm(e1); e2 = np.cross(d, e1)
    return (e1 + 1j * e2) / np.sqrt(2)


def build_crystal():
    kf = lambda p: tuple(np.round(p, 6))
    edges = set(); tris = []
    for g in itertools.product(range(3), repeat=3):
        g = np.array(g, float); c = g + 0.5
        corners = [g + np.array(o, float) for o in itertools.product((0, 1), repeat=3)]
        for axn in range(3):
            for s in (0, 1):
                fc = [v for v in corners if abs(v[axn] - (g[axn] + s)) < 1e-9]
                ctr = np.mean(fc, axis=0); oth = [t for t in range(3) if t != axn]
                fc = sorted(fc, key=lambda v: np.arctan2(v[oth[1]] - ctr[oth[1]], v[oth[0]] - ctr[oth[0]]))
                for k in range(4):
                    vi, vj = fc[k], fc[(k + 1) % 4]
                    e1 = frozenset([kf(c), kf(vi)]); e2 = frozenset([kf(c), kf(vj)]); e3 = frozenset([kf(vi), kf(vj)])
                    edges |= {e1, e2, e3}; tris.append((e1, e2, e3))
    edges = sorted(edges, key=lambda e: sorted(e))
    pos = {e: np.mean([np.array(p) for p in e], axis=0) for e in edges}
    direc = {e: np.array(sorted(e)[1]) - np.array(sorted(e)[0]) for e in edges}
    cell = lambda e: tuple(np.floor(pos[e] + 1e-6).astype(int))
    orb = {e: kf(np.array(pos[e]) - np.array(cell(e))) for e in edges}
    orbits = sorted(set(orb.values())); oid = {o: i for i, o in enumerate(orbits)}; nb = len(orbits)
    odir = {}
    for e in edges:
        odir.setdefault(orb[e], direc[e] / np.linalg.norm(direc[e]))
    hop = []
    for (a, b, cc) in tris:
        for e in (a, b, cc):
            if cell(e) == (1, 1, 1):
                for f in (a, b, cc):
                    if f != e:
                        hop.append((oid[orb[e]], oid[orb[f]], tuple(np.array(cell(f)) - np.array(cell(e)))))
    return nb, hop, [pol(odir[orbits[i]]) for i in range(nb)]


def main():
    print("=== Photon Chern: connection audit (overlap vs pure phase) ===\n")
    nb, hop, EPS = build_crystal()

    def Hk(k, mode):
        H = np.zeros((nb, nb), complex)
        for oi, oj, d in hop:
            ov = np.vdot(EPS[oi], EPS[oj])
            amp = ov if mode == 'overlap' else np.exp(1j * np.angle(ov))
            H[oi, oj] += amp * np.exp(1j * np.dot(k, d))
        return (H + H.conj().T) / 2

    def chern(ncut, kz, mode, M=20):
        kk = np.linspace(0, 2 * np.pi, M, endpoint=False); U = np.empty((M, M, nb, ncut), complex)
        for i, kx in enumerate(kk):
            for j, ky in enumerate(kk):
                _, v = np.linalg.eigh(Hk(np.array([kx, ky, kz]), mode)); U[i, j] = v[:, :ncut]
        F = 0.
        for i in range(M):
            for j in range(M):
                ip, jp = (i + 1) % M, (j + 1) % M
                lk = lambda P, Q: np.linalg.det(U[P[0], P[1]].conj().T @ U[Q[0], Q[1]])
                F += np.angle(lk((i, j), (ip, j)) * lk((ip, j), (ip, jp)) * lk((ip, jp), (i, jp)) * lk((i, jp), (i, j)))
        return F / (2 * np.pi)

    mags = [abs(np.vdot(EPS[oi], EPS[oj])) for oi, oj, d in hop]
    print(f"[A] hop-amplitude magnitudes under the OVERLAP convention: min {min(mags):.3f}, max {max(mags):.3f} "
          f"(NON-uniform -> not a gauge field)")
    ok(max(mags) - min(mags) > 0.3, "overlap hop magnitude is strongly non-uniform -> NOT a U(1) Peierls phase")

    c_ov = [round(chern(3, kz, 'overlap'), 2) for kz in (0, np.pi / 2, np.pi)]
    print(f"    overlap convention   : C(lowest-3) on kz=0,pi/2,pi = {c_ov}   <- the earlier 'C=-1' result")
    ok(all(abs(x + 1) < 0.1 for x in c_ov), "overlap convention reproduces C=-1 (the reported result)")

    # pure phase: find its clean gapped groups and show none carry the -1
    pts = [[0, 0, 0], [np.pi, 0, 0], [np.pi, np.pi, 0], [np.pi, np.pi, np.pi]]
    ks = [(1 - t) * np.array(pts[a], float) + t * np.array(pts[(a + 1) % 4], float)
          for a in range(4) for t in np.linspace(0, 1, 12, endpoint=False)]
    Bnd = np.array([np.linalg.eigvalsh(Hk(k, 'phase')) for k in ks])
    gaps = [b for b in range(nb - 1) if Bnd[:, b + 1].min() - Bnd[:, b].max() > 0.05]
    c_ph = {g + 1: [round(chern(g + 1, kz, 'phase'), 2) for kz in (0, np.pi / 2, np.pi)] for g in gaps}
    print(f"    pure-phase (physical): gapped groups {list(c_ph)}, Chern = {c_ph}   <- TRIVIAL")
    ok(all(all(abs(x) < 0.1 for x in v) for v in c_ph.values()),
       "pure-phase (Peierls/Berry, canon 7.2's own model) gives C=0 -> chiral photon NOT established")

    # the graph itself is right: cuboctahedron unphased spectrum
    V = np.array([(a, b, 0) for a in (-1, 1) for b in (-1, 1)] + [(a, 0, b) for a in (-1, 1) for b in (-1, 1)]
                 + [(0, a, b) for a in (-1, 1) for b in (-1, 1)], float)
    A = np.zeros((12, 12))
    for i in range(12):
        for j in range(i + 1, 12):
            if abs(np.linalg.norm(V[i] - V[j]) - np.sqrt(2)) < 1e-9:
                A[i, j] = A[j, i] = 1
    ok(np.allclose(np.sort(np.linalg.eigvalsh(A)), [-2, -2, -2, -2, -2, 0, 0, 0, 2, 2, 2, 4]),
       "cuboctahedron unphased spectrum [-2x5,0x3,2x3,4] matches canon 7.2 (the photon GRAPH is right)")

    print("\n[B] macro photon emergent Lorentz invariance (the real massless photon; independent of the above):")
    K = lambda k: 6 - 2 * sum(np.cos(k))
    rs = []
    for d in ([1, 0, 0], [1, 1, 0], [1, 1, 1]):
        d = np.array(d, float); d /= np.linalg.norm(d); rs.append(K(0.05 * d) / 0.05 ** 2)
    print(f"    K(k)/|k|^2 along [100],[110],[111] at |k|=0.05: {[round(r,5) for r in rs]} -> all ->1 (isotropic O(k^2))")
    ok(max(rs) - min(rs) < 1e-3 and all(abs(r - 1) < 1e-2 for r in rs),
       "macro photon dispersion isotropic to O(k^2): omega~|k| -> emergent LI (anisotropy is O((a0k)^2), RG-irrelevant)")

    print("\n[verdict] HONEST SELF-CORRECTION of the 2026-06-25 chiral-photon result:")
    print("  - C=-1 was an artifact of the OVERLAP convention (non-uniform hop magnitude 0.21-0.79, not a")
    print("    gauge field). Under canon 7.2's own scalar pure-PHASE (Peierls) connection the corrected-")
    print("    geometry crystal is C=0 -> the chiral photon is NOT established; it is RE-OPENED, not vindicated.")
    print("  - canon 7.2's C_{S7}=-1 was on the retired 'L(TCH)' framing; on the corrected oblate-bipyramid")
    print("    crystal a pure-phase connection gives a trivial photon. Re-deriving (or refuting) C=-1 with a")
    print("    physical pure-phase chiral connection on the correct geometry is the open task.")
    print("  - STANDS: the geometry correction; cuboctahedron = L(Q3) = dual cell (the photon GRAPH);")
    print("    the cuboctahedron unphased spectrum; and the MACRO photon's emergent Lorentz invariance. exit 0")


if __name__ == "__main__":
    main()
