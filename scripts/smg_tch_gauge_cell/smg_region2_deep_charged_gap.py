#!/usr/bin/env python3
r"""REGION-II deep charged-sector (SMG) gap sweep -- the VOLUME leg of
smg_region2_beta_sweep.py, on the existing TCH gauge engine
(tch_nonreduced_spin_network_scaling), via SECTOR-RESTRICTED sparse eigensolves.

The SMG observable is the charged-excitation sector ground minus the vacuum sector ground,
where the sector is the CHARGED-SITE COUNT (model.charges != 0), matching the working 336-cell
sweep (smg_region2_beta_sweep.charged_counts). NB matter_diagonal is the matter ENERGY, which
is identically zero in this truncation -- it is NOT the charge; the conserved charge that
superselects the sectors is model.charges (Gauss law forces net triality 0, so the sector is
the count of charged excitations, not the net charge). H is block-diagonal in charged count
([N_charged,H]=0), so each sector is an exact sub-block. Only the electric diagonal depends on
beta (electric=(1/beta) sum C), so we BUILD ONCE and rescale per beta.

PERFORMANCE: the engine's matvec is pure-Python; for 1.47M states that makes eigsh
infeasible. So we assemble the transition list into a scipy sparse CSR ONCE (validated
against the engine's matvec, up to its own ~1e-8 non-Hermiticity, which our symmetrization
removes), take the exact charge-sector sub-blocks, and run a fast C-level eigsh per beta.
The model build (basis + transitions) is the dominant cost (~12 min for 242K, ~1.7 h for
1.47M), so we CACHE the assembled sparse H + diagonals + sector indices to disk; reruns
load the cache and do only the (fast) per-beta eigensolves.

Run. Local (small charged cell): python3 smg_region2_deep_charged_gap.py --case min_2x1
     Deep/volume:  /home/dave/tenpy-env/bin/python smg_region2_deep_charged_gap.py --case ext_2x1
exit 0 = sweep completes and the charged-sector gap is positive across region II.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

ROOT = Path(__file__).resolve().parent
MODEL_PATH = ROOT / "tch_nonreduced_spin_network_scaling.py"
RUN_DIR = ROOT / "smg_dmrg_runs"

CASES = {
    "smoke_2x1": dict(nx=2, ny=1, states_per_charge=3, energy_cutoff=2.0, rep_set="minimal",  hopping=0.2, n_lanczos=40, progress_every=10**9),
    "ext_1x1": dict(nx=1, ny=1, states_per_charge=3, energy_cutoff=4.0, rep_set="extended", hopping=0.2, n_lanczos=60, progress_every=10**9),
    "min_2x1": dict(nx=2, ny=1, states_per_charge=3, energy_cutoff=4.0, rep_set="minimal",  hopping=0.2, n_lanczos=60, progress_every=100000),
    "ext_2x1": dict(nx=2, ny=1, states_per_charge=3, energy_cutoff=4.0, rep_set="extended", hopping=0.2, n_lanczos=80, progress_every=100000),
}
REGION2_BETAS = [0.661156, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0]
SELFCHECK_TOL = 1e-6   # the engine matvec carries ~1e-8 non-Hermiticity; symmetrization removes it


def load_engine():
    spec = importlib.util.spec_from_file_location("tch", MODEL_PATH)
    m = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(m)
    return m


def assemble_offdiag(model) -> sp.csr_matrix:
    """Assemble the off-diagonal (hopping) part of H as a Hermitian sparse CSR."""
    rows, cols, data = [], [], []
    for col, trans in enumerate(model.transitions):
        for row, coeff in trans:
            rows.append(row); cols.append(col); data.append(coeff)
    H = sp.coo_matrix((np.asarray(data, dtype=complex), (rows, cols)),
                      shape=(model.dim, model.dim)).tocsr()
    return (H + H.getH()) * 0.5            # Hermitian part (engine matvec is symmetric to ~1e-8)


def build_or_load(case: str, eng):
    """Return (H_off csr, matter, casimir, vac_idx, ch_idx), caching the expensive build."""
    h_path = RUN_DIR / f"cache_{case}_H.npz"
    meta_path = RUN_DIR / f"cache_{case}_meta.npz"
    if h_path.exists() and meta_path.exists():
        H_off = sp.load_npz(h_path)
        meta = np.load(meta_path)
        print(f"  loaded validated cache: dim={H_off.shape[0]} nnz={H_off.nnz} ({h_path.name})", flush=True)
        return H_off, meta["matter"], meta["casimir"], meta["vac_idx"], meta["ch_idx"]

    template = dict(CASES[case]); template["beta"] = 1.0   # build at beta=1; rescale per beta below
    t0 = time.time()
    model = eng.SpinNetworkCell(**template)
    print(f"  built dim={model.dim} in {time.time()-t0:.1f}s; assembling sparse H ...", flush=True)
    ta = time.time()
    H_off = assemble_offdiag(model)
    matter = np.asarray(model.matter_diagonal, dtype=float)
    casimir = np.asarray(model.electric_diagonal, dtype=float)   # at beta=1: = sum Casimir
    print(f"  sparse H: nnz={H_off.nnz} in {time.time()-ta:.1f}s", flush=True)

    # correctness self-check: sparse H_off + diag(matter+casimir) == engine matvec (up to non-Herm noise)
    rng = np.random.default_rng(7)
    v = rng.standard_normal(model.dim) + 1j * rng.standard_normal(model.dim)
    ref = model.matvec(v.astype(complex))
    mine = H_off @ v + (matter + casimir) * v
    err = float(np.linalg.norm(mine - ref) / (np.linalg.norm(ref) + 1e-30))
    print(f"  self-check ||sparse - matvec|| / ||matvec|| = {err:.2e} (tol {SELFCHECK_TOL:.0e})", flush=True)
    assert err < SELFCHECK_TOL, "sparse assembly disagrees with engine matvec"

    # SMG sector split: the CHARGED-EXCITATION count (model.charges != 0), NOT matter_diagonal
    # (which is the matter ENERGY, identically zero here). Net triality is forced to 0 by Gauss
    # law, so the sector is the count of charged sites. Matches smg_region2_beta_sweep.
    charges = np.asarray(model.charges)
    qcount = np.fromiter(
        (sum(1 for s in config if charges[s] != 0) for (config, _r, _m) in model.basis),
        dtype=np.int32, count=model.dim,
    )
    vac_idx = np.where(qcount == 0)[0]
    ch_idx = np.where(qcount > 0)[0]
    print(f"  charge-sector split: vac={len(vac_idx)} charged={len(ch_idx)} "
          f"(max charged sites/state={int(qcount.max())})", flush=True)
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    sp.save_npz(h_path, H_off)
    np.savez(meta_path, matter=matter, casimir=casimir, vac_idx=vac_idx, ch_idx=ch_idx)
    print(f"  cached -> {h_path.name}, {meta_path.name}", flush=True)
    return H_off, matter, casimir, vac_idx, ch_idx


def sector_ground(H_off_sec: sp.csr_matrix, diag_sec: np.ndarray) -> float:
    H = H_off_sec + sp.diags(diag_sec.astype(complex))
    k = H.shape[0]
    if k == 1:
        return float(np.real(H[0, 0]))
    kk = min(4, k - 1)
    vals = eigsh(H, k=kk, which="SA", maxiter=20000, tol=1e-9, return_eigenvectors=False)
    return float(np.min(vals))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="min_2x1", choices=list(CASES))
    ap.add_argument("--betas", nargs="+", type=float, default=REGION2_BETAS)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    RUN_DIR.mkdir(parents=True, exist_ok=True)

    eng = load_engine()
    print(f"REGION-II deep charged-sector gap | case={args.case}", flush=True)
    H_off, matter, casimir, vac_idx, ch_idx = build_or_load(args.case, eng)
    dim = H_off.shape[0]

    print(f"  vacuum-sector states={len(vac_idx)}  charged-sector states={len(ch_idx)}", flush=True)
    assert len(vac_idx) and len(ch_idx), "need both charge sectors"
    # block-diagonality: the sector eigensolve neglects inter-sector hopping. The charged-site
    # count is broken only by 3(x)3bar->1 annihilation (a tiny matrix element; net triality is the
    # one exact charge, pinned to 0 by Gauss law). Require the coupling << gap scale so the block
    # restriction is accurate to second order (eigenvalue error ~ cross^2/gap).
    cross = H_off[vac_idx][:, ch_idx]
    cross_max = float(abs(cross).max()) if cross.nnz else 0.0
    print(f"  block-diagonality: cross-sector |H_off| max = {cross_max:.2e} "
          f"(restriction error ~ {cross_max**2:.1e}/gap)", flush=True)
    assert cross_max < 1e-3, "hopping strongly couples charge sectors -- sector restriction invalid"
    H_vac = H_off[vac_idx][:, vac_idx].tocsr()
    H_ch = H_off[ch_idx][:, ch_idx].tocsr()
    matter_vac, cas_vac = matter[vac_idx], casimir[vac_idx]
    matter_ch, cas_ch = matter[ch_idx], casimir[ch_idx]

    rows = []
    print(f"  {'beta':>8s} {'E0_vac':>13s} {'E0_charged':>13s} {'charged_gap':>12s} {'wall_s':>8s}", flush=True)
    for b in args.betas:
        tt = time.time()
        e_vac = sector_ground(H_vac, matter_vac + cas_vac / b)
        e_ch = sector_ground(H_ch, matter_ch + cas_ch / b)
        gap = e_ch - e_vac
        rows.append(dict(beta=float(b), E0_vac=e_vac, E0_charged=e_ch,
                         charged_gap=gap, wall_seconds=time.time() - tt))
        print(f"  {b:>8.4f} {e_vac:>13.5f} {e_ch:>13.5f} {gap:>12.5f} {rows[-1]['wall_seconds']:>8.1f}", flush=True)

    cg = [r["charged_gap"] for r in rows]
    min_cg = min(cg); min_at = rows[cg.index(min_cg)]["beta"]
    result = dict(case=args.case, dim=int(dim), betas=args.betas, rows=rows,
                  min_charged_gap=min_cg, min_at_beta=min_at, no_closure=min_cg > 1e-2)
    out = args.out or RUN_DIR / f"smg_region2_deep_charged_gap_{args.case}.json"
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(f"\n  region-II charged-sector gap: min {min_cg:.5f} at beta={min_at:.4f} "
          f"(dim={dim}); no_closure={min_cg > 1e-2}; wrote {out}", flush=True)
    assert all(r["charged_gap"] > 0.0 for r in rows), "charged-sector gap closed"
    print("exit 0 -- region-II charged-sector gap open across the window at this volume.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
