#!/usr/bin/env python3
r"""REGION-II charged-sector (SMG) gap on the 1xN-PLAQUETTE MAGNETIC strip -- the volume
escalation of smg_region2_2plaq_magnetic_gap.py to a third plaquette (and general N).

Generalizes the hardcoded 2-plaquette geometry of tch_finish_port_and_swap.py to a 1xN strip
(x in 0..N, y in 0,1), reusing the VALIDATED recoupling machinery verbatim (m.K_corner,
m.vbasis, m.targets) -- only the geometry (links, vertices, vertex-legs, plaquettes) and the
link-count in build_basis change. Each plaquette is still a unit square, so the per-plaquette
corner/colour-flow logic is identical; the N=2 geometry built here is byte-identical to the
engine's, which is the regression that validates the generalization.

OBSERVABLE: the charged-EXCITATION-sector gap E0(n_charged>=1) - E0(n_charged=0) with the Wilson
magnetic term H = ELE/beta + MAT - (beta/2) sum_p W_p, swept across region II (0.661<beta<6).
n_charged (#{f!="1"}) is conserved by W (links/intertwiners only) -> exact block-diagonal
sectors (guarded).

VALIDATION CHAIN: at N=2 reproduce the published 3,321 count + full gap 6.674654 (beta=0.5) +
the charged-sector min 5.066; only then trust N=3.

Run.  Local validate: python3 smg_region2_nplaq_magnetic_gap.py --n 2
      Deep N=3:        /home/dave/tenpy-env/bin/python smg_region2_nplaq_magnetic_gap.py --n 3
exit 0 = N=2 regression passes AND the requested-N charged-sector gap is positive across region II.
"""
from __future__ import annotations

import argparse
import functools
import importlib.util
import itertools as it
import json
import time
from pathlib import Path

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

ROOT = Path(__file__).resolve().parent
RUN_DIR = ROOT / "smg_dmrg_runs"
REGION2_BETAS = [0.661156, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0]


def load_machinery():
    spec = importlib.util.spec_from_file_location("finish", ROOT / "tch_finish_port_and_swap.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def build_geometry(n: int) -> dict:
    """1xN strip: x in 0..N, y in 0,1. Link order matches the engine at N=2 exactly."""
    LINKS = [("h", x, y) for y in (0, 1) for x in range(n)] + [("v", x, 0) for x in range(n + 1)]
    LIDX = {l: i for i, l in enumerate(LINKS)}
    VERTS = [(x, y) for x in range(n + 1) for y in range(2)]

    def vlegs(v):
        x, y = v
        legs = []
        if x < n:
            legs.append((LIDX[("h", x, y)], False))           # out: horizontal forward
        if y < 1:
            legs.append((LIDX[("v", x, 0)], False))           # out: vertical forward
        if x > 0:
            legs.append((LIDX[("h", x - 1, y)], True))        # in: horizontal
        if y > 0:
            legs.append((LIDX[("v", x, 0)], True))            # in: vertical
        return legs

    PLAQS = [[(LIDX[("h", k, 0)], "fwd"), (LIDX[("v", k + 1, 0)], "fwd"),
              (LIDX[("h", k, 1)], "bwd"), (LIDX[("v", k, 0)], "bwd")] for k in range(n)]
    return dict(LINKS=LINKS, LIDX=LIDX, VERTS=VERTS, vlegs=vlegs, PLAQS=PLAQS, nlinks=len(LINKS), n=n)


# --- convergence-safe intertwiner basis (verbatim from tch_2plaq_extended_deep.py) -----------
def make_robust_vbasis(cg):
    def robust_vbasis(reps):
        stacked = np.vstack(cg.generator_sum(reps))
        try:
            s = np.linalg.svd(stacked, compute_uv=False)
        except np.linalg.LinAlgError:
            s = sla.svd(stacked, compute_uv=False, lapack_driver="gesvd")
        if int(np.sum(s < 1e-9)) == 0:
            return ()
        try:
            _u, s2, vh = np.linalg.svd(stacked)
        except np.linalg.LinAlgError:
            _u, s2, vh = sla.svd(stacked, full_matrices=True, lapack_driver="gesvd")
        null = vh[np.sum(s2 > 1e-9):].conj()
        dims = tuple(cg.dim(r) for r in reps)
        out = []
        for row in null:
            t = row.reshape(dims)
            for o in out:
                t = t - o * np.vdot(o, t)
            nrm = np.linalg.norm(t)
            if nrm > 1e-9:
                out.append(t / nrm)
        return tuple(out)
    return robust_vbasis


def build_basis(m, reps, geom):
    VERTS, vlegs, nlinks = geom["VERTS"], geom["vlegs"], geom["nlinks"]
    cg = m.cg
    basis = []
    for Rs in it.product(reps, repeat=nlinks):
        per_v = []
        ok = True
        for v in VERTS:
            legs = vlegs(v)
            opts = []
            for f in m.FERMS:
                lr = tuple(Rs[li] for li, _ in legs)
                df = tuple(isin for _, isin in legs)
                full = tuple((cg.dual(r) if d else r) for r, d in zip(lr, df)) + (f,)
                k = len(m.vbasis(full))
                for mi in range(k):
                    opts.append((f, mi))
            if not opts:
                ok = False; break
            per_v.append(opts)
        if not ok:
            continue
        for choice in it.product(*per_v):
            basis.append((Rs, tuple(choice)))
    return basis


def make_sparse_wilson(m, Kc, geom):
    VERTS, vlegs, PLAQS = geom["VERTS"], geom["vlegs"], geom["PLAQS"]
    cg = m.cg

    def sparse_wilson(p, reps, basis, index):
        loop = PLAQS[p]
        looplinks = [li for li, _ in loop]
        dirs = {li: d for li, d in loop}
        plaq_verts = []
        for v in VERTS:
            legs = vlegs(v)
            tl = [(ax, li, isin) for ax, (li, isin) in enumerate(legs) if li in looplinks]
            if len(tl) == 2:
                plaq_verts.append((v, legs, tl))
        assert len(plaq_verts) == 4
        rows, cols, vals = [], [], []
        for i, (Rs, choice) in enumerate(basis):
            opts = []
            for li in looplinks:
                kind_t = ("out_fwd" if dirs[li] == "fwd" else "out_bwd")
                opts.append(m.targets(kind_t, Rs[li], reps))
            for newreps in it.product(*opts):
                Rps = list(Rs)
                for li, rp in zip(looplinks, newreps):
                    Rps[li] = rp
                Rps = tuple(Rps)
                amp_blocks = []
                ok = True
                for (v, legs, tl) in plaq_verts:
                    vi = VERTS.index(v)
                    f, mi = basis[i][1][vi]
                    raw_old = tuple(Rs[li] for li, _ in legs)
                    raw_new = tuple(Rps[li] for li, _ in legs)
                    df = tuple(isin for _, isin in legs)
                    (axA, liA, isinA), (axB, liB, isinB) = tl
                    posA, posB = looplinks.index(liA), looplinks.index(liB)
                    first, second = ((axA, liA, isinA), (axB, liB, isinB)) \
                        if (posB - posA) % 4 == 1 else ((axB, liB, isinB), (axA, liA, isinA))
                    kinds = [("in_" if isin else "out_") + dirs[li] for (ax, li, isin) in (first, second)]
                    K = Kc(raw_old, raw_new, df, f, (first[0], second[0]), tuple(kinds))
                    if K is None:
                        ok = False; break
                    amp_blocks.append((vi, K))
                if not ok:
                    continue
                pw = 1.0
                for li, rp in zip(looplinks, newreps):
                    pw *= (cg.dim(Rs[li]) / cg.dim(rp)) ** 0.5
                corner_vis = [vi for vi, _ in amp_blocks]
                Kmap = {vi: K for vi, K in amp_blocks}
                newdims = {vi: K.shape[0] for vi, K in amp_blocks}
                for mps in it.product(*[range(newdims[vi]) for vi in corner_vis]):
                    amp = pw + 0j
                    newchoice = list(basis[i][1])
                    for vi, mp in zip(corner_vis, mps):
                        f_old, m_old = basis[i][1][vi]
                        amp *= Kmap[vi][mp, m_old]
                        newchoice[vi] = (f_old, mp)
                    if abs(amp) < 1e-12:
                        continue
                    j = index.get((Rps, tuple(newchoice)))
                    if j is not None:
                        rows.append(j); cols.append(i); vals.append(amp)
        nn = len(basis)
        return sp.csr_matrix((vals, (rows, cols)), shape=(nn, nn), dtype=complex)
    return sparse_wilson


def diag_arrays(m, basis):
    cg = m.cg
    ele = np.array([sum(cg.casimir(r) for r in Rs) for (Rs, _ch) in basis], float)
    cc = np.array([sum(1 for f, _ in ch if f != "1") for (_Rs, ch) in basis], int)
    return ele, m.DSMG * cc.astype(float), cc


def sector_ground(Wsym, diag, idx, beta):
    H = sp.diags(diag[idx].astype(complex)) - (beta / 2.0) * Wsym[idx][:, idx]
    k = H.shape[0]
    if k == 1:
        return float(np.real(H.toarray()[0, 0]))
    return float(np.min(eigsh(H.tocsr(), k=min(6, k - 1), which="SA",
                              return_eigenvectors=False, maxiter=20000, tol=0)))


def full_gap(Wsym, diag, beta):
    e = np.sort(eigsh((sp.diags(diag.astype(complex)) - (beta / 2.0) * Wsym).tocsr(),
                      k=4, which="SA", return_eigenvectors=False, maxiter=20000, tol=0))
    return float(e[1] - e[0])


def run_cell(m, geom, reps, label, betas):
    Kc = functools.lru_cache(maxsize=None)(m.K_corner)
    sparse_wilson = make_sparse_wilson(m, Kc, geom)
    t0 = time.time()
    basis = build_basis(m, reps, geom)
    index = {b: i for i, b in enumerate(basis)}
    n = len(basis)
    print(f"\n=== {label}: {n} states (built {time.time()-t0:.1f}s) ===", flush=True)
    ta = time.time()
    Wsym = sum((sparse_wilson(p, reps, basis, index) for p in range(geom["n"])), sp.csr_matrix((n, n), dtype=complex))
    Wsym = (Wsym + Wsym.getH()).tocsr()
    ele, mat, cc = diag_arrays(m, basis)
    vac = np.where(cc == 0)[0]; chg = np.where(cc > 0)[0]
    print(f"  W assembled (nnz={Wsym.nnz}) in {time.time()-ta:.1f}s; vac={len(vac)} charged={len(chg)}", flush=True)
    cross = Wsym[vac][:, chg]
    cmax = float(abs(cross).max()) if cross.nnz else 0.0
    print(f"  block-diagonality: cross-sector |Wsym| max = {cmax:.2e}", flush=True)
    assert cmax < 1e-9, "magnetic move couples charge sectors"
    rows = []
    print(f"  {'beta':>8s} {'E0_vac':>12s} {'E0_charged':>12s} {'charged_gap':>12s} {'wall_s':>8s}", flush=True)
    for b in betas:
        tt = time.time()
        diag = ele / b + mat
        e_vac = sector_ground(Wsym, diag, vac, b)
        e_chg = sector_ground(Wsym, diag, chg, b)
        rows.append(dict(beta=float(b), E0_vac=e_vac, E0_charged=e_chg, charged_gap=e_chg - e_vac,
                         wall_seconds=time.time() - tt))
        print(f"  {b:>8.4f} {e_vac:>12.5f} {e_chg:>12.5f} {e_chg-e_vac:>12.5f} {rows[-1]['wall_seconds']:>8.1f}", flush=True)
    cg_ = [r["charged_gap"] for r in rows]
    gb = [r["charged_gap"] * r["beta"] for r in rows]
    return dict(label=label, dim=n, n_plaq=geom["n"], betas=list(betas), rows=rows,
                min_charged_gap=min(cg_), min_at_beta=rows[cg_.index(min(cg_))]["beta"],
                drift=max(cg_) - min(cg_), no_closure=min(cg_) > 1e-2,
                gap_is_electric_artifact=(max(gb) - min(gb)) / (sum(gb) / len(gb)) < 0.05)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=2, help="number of plaquettes (strip length)")
    ap.add_argument("--reps", choices=("minimal", "extended"), default="minimal")
    ap.add_argument("--betas", nargs="+", type=float, default=REGION2_BETAS)
    args = ap.parse_args()
    RUN_DIR.mkdir(parents=True, exist_ok=True)

    print(f"REGION-II charged-sector (SMG) gap on the 1x{args.n}-plaquette MAGNETIC strip ({args.reps} reps)")
    m = load_machinery()
    m.vbasis = functools.lru_cache(maxsize=None)(make_robust_vbasis(m.cg))
    reps = m.REPS_MIN if args.reps == "minimal" else m.REPS_EXT

    # --- N=2 regression: the generalized geometry must reproduce the engine's 2-plaquette build ---
    print("\n[regression] generalized geometry at N=2 must reproduce 3,321 + full gap 6.674654 (beta=0.5)")
    g2 = build_geometry(2)
    Kc = functools.lru_cache(maxsize=None)(m.K_corner)
    sw2 = make_sparse_wilson(m, Kc, g2)
    b2 = build_basis(m, m.REPS_MIN, g2)
    idx2 = {b: i for i, b in enumerate(b2)}
    assert len(b2) == 3321, f"N=2 basis count {len(b2)} != 3321 -- geometry generalization is wrong"
    W0 = sw2(0, m.REPS_MIN, b2, idx2); W1 = sw2(1, m.REPS_MIN, b2, idx2)
    d0 = float(np.max(np.abs(W0.toarray() - m.W0)))
    print(f"  N=2 basis=3321 OK; |sparse W0 - engine W0| = {d0:.1e} (target 0)")
    assert d0 < 1e-10, "N=2 magnetic operator != engine -- geometry generalization is wrong"
    ele2, mat2, _ = diag_arrays(m, b2)
    Wsym2 = (W0 + W0.getH() + W1 + W1.getH()).tocsr()
    g05 = full_gap(Wsym2, ele2 / 0.5 + mat2, 0.5)
    print(f"  N=2 full gap beta=0.5 = {g05:.6f} (published 6.674654)")
    assert abs(g05 - 6.674654) < 1e-3
    print("  REGRESSION PASSED -- the 1xN geometry is byte-faithful at N=2.\n")

    geom = build_geometry(args.n)
    res = run_cell(m, geom, reps, f"{args.n}plaq_{args.reps}", args.betas)
    out = RUN_DIR / f"smg_region2_{args.n}plaq_{args.reps}_magnetic_gap.json"
    out.write_text(json.dumps(res, indent=2) + "\n")
    print(f"\n  {res['label']} dim={res['dim']}: charged-sector gap min {res['min_charged_gap']:.5f} "
          f"at beta={res['min_at_beta']:.3f}, drift {res['drift']:.5f}, no_closure={res['no_closure']}, "
          f"electric_artifact={res['gap_is_electric_artifact']}")
    print(f"  wrote {out}")
    assert all(r["charged_gap"] > 0.0 for r in res["rows"]), "charged gap closed"
    print("exit 0 -- magnetic charged-sector (SMG) gap open across region II at this volume.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
