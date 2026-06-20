#!/usr/bin/env python3
r"""REGION-II charged-sector (SMG) gap on the 2-PLAQUETTE MAGNETIC cell -- the proper
VOLUME leg, WITH the Wilson magnetic term (the piece tch_nonreduced_spin_network_scaling
lacks, which made its charged gap collapse to the electric artifact C3/beta).

WHY this engine. smg_region2_beta_sweep.py established the region-II charged-sector gap on
the 336-state 1-plaquette magnetic cell (tch_nonreduced_magnetic_plaquette: H = ELE/beta +
MAT + t*HOP - (beta/2) W). The deep volume leg on the spin-network scaling engine measured
only ELE/beta (no magnetic W) -> gap = C3/beta -> 0 as beta->inf (a finite-cell artifact,
NOT the SMG mass). This script scales the CORRECT (magnetic) engine to the 2-plaquette cells
(3,321 minimal; 106,460 extended; tch_finish_port_and_swap + tch_2plaq_extended_deep) and runs
the SMG observable -- the charged-EXCITATION-sector gap E0(n_charged>=1) - E0(n_charged=0) --
across region II (0.661 < beta < 6).

The charged-vertex count n_charged = #{f != "1"} is conserved by W (which touches only link
reps + intertwiners, never the fermion reps) and by ELE/MAT (diagonal), so H is exactly
block-diagonal in n_charged -> the sector-restricted Lanczos is exact (guarded below).

DISCIPLINE: reuse the VALIDATED machinery verbatim (tch_finish_port_and_swap GATE A/B on
import; the sparse_wilson + robust SVD from tch_2plaq_extended_deep), regress the minimal
FULL gap against the published 6.674654 (beta=0.5) before trusting the charged-sector numbers.

Run.  Local:  python3 smg_region2_2plaq_magnetic_gap.py            # minimal 3,321 (regress + sweep)
      Deep :  /home/dave/tenpy-env/bin/python smg_region2_2plaq_magnetic_gap.py --extended  # + 106,460
exit 0 = regression passes AND the charged-sector gap is positive across region II.
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
    spec.loader.exec_module(m)            # runs GATE A (flux character) + GATE B (TCH spectrum)
    return m


# --- convergence-safe intertwiner basis (verbatim from tch_2plaq_extended_deep.py) -----------
def _robust_svd(a, compute_uv, cg):
    try:
        return np.linalg.svd(a, compute_uv=compute_uv)
    except np.linalg.LinAlgError:
        if compute_uv:
            return sla.svd(a, full_matrices=True, lapack_driver="gesvd")
        return sla.svd(a, compute_uv=False, lapack_driver="gesvd")


def make_robust_vbasis(cg):
    def robust_vbasis(reps):
        stacked = np.vstack(cg.generator_sum(reps))
        s = _robust_svd(stacked, False, cg)
        if int(np.sum(s < 1e-9)) == 0:
            return ()
        _u, s2, vh = _robust_svd(stacked, True, cg)
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


def make_sparse_wilson(m, Kc):
    """EXACT sparse copy of m.wilson_P over an arbitrary (basis, index)."""
    VERTS, vlegs, PLAQS = m.VERTS, m.vlegs, m.PLAQS

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
                    kinds = []
                    for (ax, li, isin) in (first, second):
                        d = dirs[li]
                        kinds.append(("in_" if isin else "out_") + d)
                    K = Kc(raw_old, raw_new, df, f, (first[0], second[0]), tuple(kinds))
                    if K is None:
                        ok = False; break
                    amp_blocks.append((vi, K))
                if not ok:
                    continue
                pw = 1.0
                for li, rp in zip(looplinks, newreps):
                    pw *= (m.cg.dim(Rs[li]) / m.cg.dim(rp)) ** 0.5
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
        n = len(basis)
        return sp.csr_matrix((vals, (rows, cols)), shape=(n, n), dtype=complex)
    return sparse_wilson


def diag_arrays(m, basis):
    cg = m.cg
    ele = np.array([sum(cg.casimir(r) for r in Rs) for (Rs, _ch) in basis], float)
    mat = np.array([m.DSMG * sum(1 for f, _ in ch if f != "1") for (_Rs, ch) in basis], float)
    cc = np.array([sum(1 for f, _ in ch if f != "1") for (_Rs, ch) in basis], int)
    return ele, mat, cc


def sector_ground(Wsym, diag, idx, beta):
    H = sp.diags(diag[idx].astype(complex)) - (beta / 2.0) * Wsym[idx][:, idx]
    k = H.shape[0]
    if k == 1:
        return float(np.real(H.toarray()[0, 0]))
    vals = eigsh(H.tocsr(), k=min(6, k - 1), which="SA", return_eigenvectors=False, maxiter=20000, tol=0)
    return float(np.min(vals))


def full_gap(Wsym, diag, beta):
    H = sp.diags(diag.astype(complex)) - (beta / 2.0) * Wsym
    e = np.sort(eigsh(H.tocsr(), k=4, which="SA", return_eigenvectors=False, maxiter=20000, tol=0))
    return float(e[1] - e[0])


def run_cell(m, sparse_wilson, reps, label, betas):
    basis = m.build_basis(reps)
    index = {b: i for i, b in enumerate(basis)}
    n = len(basis)
    print(f"\n=== {label}: {n} states ===", flush=True)
    t0 = time.time()
    W0 = sparse_wilson(0, reps, basis, index)
    W1 = sparse_wilson(1, reps, basis, index)
    Wsym = (W0 + W0.getH() + W1 + W1.getH()).tocsr()
    ele, mat, cc = diag_arrays(m, basis)
    vac = np.where(cc == 0)[0]
    chg = np.where(cc > 0)[0]
    print(f"  W assembled (nnz={Wsym.nnz}) in {time.time()-t0:.1f}s; vac={len(vac)} charged={len(chg)}", flush=True)
    # block-diagonality guard: W must preserve the charged-vertex count
    cross = Wsym[vac][:, chg]
    cross_max = float(abs(cross).max()) if cross.nnz else 0.0
    print(f"  block-diagonality: cross-sector |Wsym| max = {cross_max:.2e}", flush=True)
    assert cross_max < 1e-9, "magnetic move couples charge sectors -- sector restriction invalid"
    rows = []
    print(f"  {'beta':>8s} {'E0_vac':>12s} {'E0_charged':>12s} {'charged_gap':>12s} {'wall_s':>8s}", flush=True)
    for b in betas:
        tt = time.time()
        diag = ele / b + mat
        e_vac = sector_ground(Wsym, diag, vac, b)
        e_chg = sector_ground(Wsym, diag, chg, b)
        gap = e_chg - e_vac
        rows.append(dict(beta=float(b), E0_vac=e_vac, E0_charged=e_chg, charged_gap=gap,
                         wall_seconds=time.time() - tt))
        print(f"  {b:>8.4f} {e_vac:>12.5f} {e_chg:>12.5f} {gap:>12.5f} {rows[-1]['wall_seconds']:>8.1f}", flush=True)
    cg_ = [r["charged_gap"] for r in rows]
    min_cg = min(cg_); min_at = rows[cg_.index(min_cg)]["beta"]
    drift = max(cg_) - min(cg_)
    # is the gap electric (C3/beta, drift ~ full) or a real beta-independent SMG mass?
    gb = [r["charged_gap"] * r["beta"] for r in rows]
    electric_like = (max(gb) - min(gb)) / (sum(gb) / len(gb)) < 0.05
    return dict(label=label, dim=n, betas=list(betas), rows=rows, min_charged_gap=min_cg,
                min_at_beta=min_at, drift=drift, no_closure=min_cg > 1e-2,
                gap_is_electric_artifact=electric_like)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--extended", action="store_true", help="also run the 106,460 extended cell")
    ap.add_argument("--betas", nargs="+", type=float, default=REGION2_BETAS)
    args = ap.parse_args()
    RUN_DIR.mkdir(parents=True, exist_ok=True)

    print("REGION-II charged-sector (SMG) gap on the 2-PLAQUETTE MAGNETIC cell")
    m = load_machinery()
    Kc = functools.lru_cache(maxsize=None)(m.K_corner)
    m.vbasis = functools.lru_cache(maxsize=None)(make_robust_vbasis(m.cg))
    sparse_wilson = make_sparse_wilson(m, Kc)

    # --- regression: reproduce the published minimal FULL gap before trusting the sector gap ---
    print("\n[regression] minimal 3,321 FULL gap must reproduce the published 6.674654 (beta=0.5)")
    basis_min = m.build_basis(m.REPS_MIN)
    index_min = {b: i for i, b in enumerate(basis_min)}
    assert len(basis_min) == 3321
    W0 = sparse_wilson(0, m.REPS_MIN, basis_min, index_min)
    W1 = sparse_wilson(1, m.REPS_MIN, basis_min, index_min)
    d0 = float(np.max(np.abs(W0.toarray() - m.W0))); d1 = float(np.max(np.abs(W1.toarray() - m.W1)))
    print(f"  elementwise |sparse W - dense W| = {d0:.1e}, {d1:.1e} (target 0)")
    assert d0 < 1e-10 and d1 < 1e-10, "sparse assembler != validated dense machinery"
    Wsym = (W0 + W0.getH() + W1 + W1.getH()).tocsr()
    ele, mat, _cc = diag_arrays(m, basis_min)
    g05 = full_gap(Wsym, ele / 0.5 + mat, 0.5)
    print(f"  minimal full gap beta=0.5 = {g05:.6f} (published 6.674654)")
    assert abs(g05 - 6.674654) < 1e-3, "minimal full-gap regression failed"
    print("  REGRESSION PASSED.\n")

    results = [run_cell(m, sparse_wilson, m.REPS_MIN, "2plaq_minimal_3321", args.betas)]
    if args.extended:
        results.append(run_cell(m, sparse_wilson, m.REPS_EXT, "2plaq_extended_106460", args.betas))

    out = RUN_DIR / "smg_region2_2plaq_magnetic_gap.json"
    out.write_text(json.dumps(results, indent=2) + "\n")
    print("\n" + "=" * 70)
    for r in results:
        print(f"  {r['label']:24s} dim={r['dim']:>7d}: charged-sector gap min {r['min_charged_gap']:.5f} "
              f"at beta={r['min_at_beta']:.3f}, drift {r['drift']:.5f}, "
              f"no_closure={r['no_closure']}, electric_artifact={r['gap_is_electric_artifact']}")
    print(f"  wrote {out}")
    for r in results:
        assert all(row["charged_gap"] > 0.0 for row in r["rows"]), f"{r['label']}: charged gap closed"
    print("exit 0 -- magnetic 2-plaquette charged-sector (SMG) gap open across region II.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
