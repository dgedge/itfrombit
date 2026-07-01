#!/usr/bin/env python3
r"""Z3 deconfinement on the bond-bipyramid TCH BULK x periodic Euclidean time.

Rung 3 of the finite-temperature programme, on the framework's OWN geometry
(the earlier torus rungs used a standard hypercubic lattice). The spatial
complex is the licensed bond-bipyramid bulk from
``record_grammar_tch_bond_bipyramid_bulk.build_bond_complex``; a periodic
Euclidean-time direction of extent Nt is added on top, giving the non-contractible
time cycle a Polyakov line needs.

Construction
------------
Spatial sites/edges/faces come from the bond complex.  Gauge links: a Z3 phase
on every spatial edge per time slice, plus a time-link on every spatial vertex
per slice (periodic in t).  Plaquettes: the spatial triangular faces (per slice)
and the time-plaquettes (spatial edge x time step).  Re(U_p)=cos is orientation-
even, so only the relative link signs within a loop matter.  Vectorised
checkerboard Metropolis (links sharing a plaquette are colour-separated).  The
order parameter is the volume-averaged Polyakov line
L = (1/Nv) sum_v w^{sum_t k_time(v,t)}.

Result
------
<|L|> melts from the finite-size floor 1/sqrt(Nv) (confined, Z3 unbroken) to ~1
(deconfined, Z3 broken), at a framework-specific beta_c that RISES with Nt:
    Nt=2 -> beta_c ~ 0.64,   Nt=4 -> beta_c ~ 0.73.
So the bond-bipyramid bulk supports a genuine finite-temperature
confinement/deconfinement transition with the correct Nt scaling, at its own
beta_c (close to but distinct from the hypercubic Z3 value ~0.67).

Boundary
--------
This is the Z3 centre on the framework bulk -- a consistency result (any 3+1D
gauge theory with a periodic time deconfines; the new content is that the
framework geometry does so, and its beta_c(Nt)).  It is NOT full SU(3) on the
bulk (the next rung), nor a continuum / finite-size-scaling study.  Companion to
``record_grammar_tch_z3_polyakov_torus.py`` (hypercubic Z3) and
``record_grammar_tch_su3_polyakov_torus.py`` (SU(3) + beta_c-vs-Nt, validated).
"""
from __future__ import annotations

import numpy as np

from record_grammar_tch_bond_bipyramid_bulk import build_bond_complex, bonds_box


def build(nx, ny, nz, Nt):
    C = build_bond_complex(f"slab{nx}{ny}{nz}", bonds_box(nx, ny, nz))
    verts = list(C.vertices)
    vidx = {v: i for i, v in enumerate(verts)}
    edges = list(C.edges)
    eidx = {e: i for i, e in enumerate(edges)}
    Nv, Ne = len(verts), len(edges)

    def slink(e, t):
        return e * Nt + (t % Nt)

    def tlink(v, t):
        return Ne * Nt + v * Nt + (t % Nt)

    Nlinks = Ne * Nt + Nv * Nt
    plaqs = []
    for (a, b, c) in C.faces:                          # spatial triangular faces (sorted a<b<c)
        eab, ebc, eac = (eidx[tuple(sorted((a, b)))],
                         eidx[tuple(sorted((b, c)))],
                         eidx[tuple(sorted((a, c)))])
        for t in range(Nt):
            plaqs.append(([slink(eab, t), slink(ebc, t), slink(eac, t), 0], [1, 1, -1, 0]))
    for ei, (u, w) in enumerate(edges):                # time-plaquettes
        iu, iw = vidx[u], vidx[w]
        for t in range(Nt):
            plaqs.append(([slink(ei, t), tlink(iw, t), slink(ei, t + 1), tlink(iu, t)], [1, 1, -1, -1]))

    plink = np.array([p[0] for p in plaqs], dtype=np.int64)
    psign = np.array([p[1] for p in plaqs], dtype=np.int64)
    Np = len(plaqs)

    inc = [[] for _ in range(Nlinks)]
    for p in range(Np):
        for s in range(4):
            if psign[p, s] != 0:
                inc[plink[p, s]].append((p, int(psign[p, s])))
    Lmax = max(len(x) for x in inc)
    lp = np.full((Nlinks, Lmax), -1, dtype=np.int64)
    ls = np.zeros((Nlinks, Lmax), dtype=np.int64)
    for l in range(Nlinks):
        for j, (p, s) in enumerate(inc[l]):
            lp[l, j], ls[l, j] = p, s

    adj = [set() for _ in range(Nlinks)]
    for p in range(Np):
        lk = [plink[p, s] for s in range(4) if psign[p, s] != 0]
        for i in range(len(lk)):
            for j in range(i + 1, len(lk)):
                adj[lk[i]].add(lk[j])
                adj[lk[j]].add(lk[i])
    color = -np.ones(Nlinks, dtype=np.int64)
    for l in range(Nlinks):
        used = {color[n] for n in adj[l] if color[n] >= 0}
        c = 0
        while c in used:
            c += 1
        color[l] = c
    groups = [np.where(color == c)[0] for c in range(int(color.max()) + 1)]
    tpoly = np.array([[tlink(v, t) for t in range(Nt)] for v in range(Nv)])
    return dict(Nlinks=Nlinks, Nv=Nv, Ne=Ne, Nt=Nt, Np=Np, plink=plink, psign=psign,
                lp=lp, ls=ls, groups=groups, tpoly=tpoly)


def plaq_exponent(k, L):
    return (L["psign"] * k[L["plink"]]).sum(1)


def mean_plaq(k, L):
    return float(np.cos(2 * np.pi * plaq_exponent(k, L) / 3.0).mean())


def sweep(k, beta, rng, L):
    for G in L["groups"]:
        E = plaq_exponent(k, L)
        lpG, lsG, kG = L["lp"][G], L["ls"][G], k[G]
        mask = (lpG >= 0)
        R = E[lpG] - lsG * kG[:, None]
        kp = (kG + rng.integers(1, 3, size=kG.shape)) % 3

        def cs(kk):
            return (mask * np.cos(2 * np.pi * (lsG * kk[:, None] + R) / 3.0)).sum(1)

        dS = -beta * (cs(kp) - cs(kG))
        k[G] = np.where(rng.random(kG.shape) < np.exp(-dS), kp, kG)
    return k


def polyakov_abs(k, L):
    return float(np.abs(np.exp(2j * np.pi * (k[L["tpoly"]].sum(1) % 3) / 3.0).mean()))


def scan(nx, ny, nz, Nt, betas, therm=200, meas=200, seed=7):
    L = build(nx, ny, nz, Nt)
    rng = np.random.default_rng(seed)
    k = rng.integers(0, 3, size=L["Nlinks"]).astype(np.int64)   # hot start, anneal up
    floor = 1.0 / np.sqrt(L["Nv"])
    print(f"  Nt={Nt}: Nv={L['Nv']} Ne={L['Ne']} links={L['Nlinks']} plaqs={L['Np']} "
          f"colours={len(L['groups'])} floor={floor:.3f}")
    rows = []
    for beta in betas:
        for _ in range(therm):
            k = sweep(k, beta, rng, L)
        m = np.array([polyakov_abs(sweep(k, beta, rng, L), L) for _ in range(meas)])
        rows.append((float(beta), float(m.mean())))
        print(f"    beta={beta:5.2f}  <|L|>={m.mean():.4f}  chi={L['Nv']*m.var():7.3f}")
    return rows, floor


def beta_c_halfcross(rows, floor):
    Ls = [r[1] for r in rows]
    thr = 0.5 * (np.mean(sorted(Ls)[-3:]) + floor)
    for i in range(1, len(rows)):
        if rows[i - 1][1] < thr <= rows[i][1]:
            (b0, l0), (b1, l1) = rows[i - 1], rows[i]
            return b0 + (b1 - b0) * (thr - l0) / (l1 - l0)
    return float("nan")


def assert_true(name, value):
    print(f"  {name:<60s} value={value}")
    if not value:
        raise AssertionError(name)


def main():
    print("Z3 deconfinement on the bond-bipyramid TCH bulk x periodic Euclidean time")
    print("=" * 100)
    rng = np.random.default_rng(0)
    L0 = build(3, 3, 2, 2)
    plaq0 = abs(mean_plaq(rng.integers(0, 3, size=L0["Nlinks"]).astype(np.int64), L0))
    print(f"  beta=0 random control: |<plaq>|={plaq0:.3f} (expect ~0)\n")

    print("  N_t=2 scan:")
    rows2, floor = scan(3, 3, 2, 2, np.round(np.arange(0.40, 1.31, 0.10), 2))
    print("  N_t=4 scan:")
    rows4, _ = scan(3, 3, 2, 4, np.round(np.arange(0.50, 1.01, 0.07), 2))

    bc2, bc4 = beta_c_halfcross(rows2, floor), beta_c_halfcross(rows4, floor)
    conf = np.mean([l for b, l in rows2 if b <= 0.5])
    deconf = np.mean([l for b, l in rows2 if b >= 1.1])
    print(f"\n  beta_c(N_t=2) = {bc2:.2f}   beta_c(N_t=4) = {bc4:.2f}   (hypercubic Z3 ~0.67)")

    print("\n[checks]")
    assert_true("beta=0 control: |<plaq>| < 0.10", plaq0 < 0.10)
    assert_true("confined <|L|> near floor (< 0.25)", conf < 0.25)
    assert_true("deconfined <|L|> O(1) (> 0.80)", deconf > 0.80)
    assert_true("beta_c(N_t=2) in (0.4, 0.9)", 0.4 < bc2 < 0.9)
    assert_true("beta_c(N_t=4) in (0.4, 0.95)", 0.4 < bc4 < 0.95)
    assert_true("beta_c rises with N_t on the bulk (bc4 > bc2 + 0.03)", bc4 > bc2 + 0.03)

    print(
        """
VERDICT:
  PASS.  The bond-bipyramid TCH bulk, given a periodic Euclidean-time direction,
  has a genuine Z3 deconfinement transition: the Polyakov order parameter melts
  from the finite-size floor (centre unbroken) to O(1) (centre broken), and its
  beta_c rises with N_t -- the correct finite-temperature scaling -- on the
  framework's own geometry, at a beta_c close to but distinct from the hypercubic
  value.  Rung 3 done.  This is the Z3 centre; full SU(3) on the bulk, and a
  continuum/finite-size-scaling study, remain the further rungs.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
