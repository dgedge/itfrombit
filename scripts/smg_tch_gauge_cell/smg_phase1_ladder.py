#!/usr/bin/env python3
r"""PHASE 1 -- volume ladder for the SMG electric-subtracted mirror gap.

REUSE-ONLY: every physics function is imported from the existing code --
  * smg_region2_nplaq_magnetic_gap.py : build_geometry, build_basis, make_sparse_wilson,
                                        diag_arrays, make_robust_vbasis, load_machinery
  * tch_2plaq_extended_hop.py         : sparse_hop, hop_tail, hop_head  (the VALIDATED
                                        gauge-covariant neutral fermion hop, |dHOP|<1e-10 vs oracle)
This driver only (a) dispatches dense-vs-sparse eigensolve by sector size (the runner's
sector_ground crashes for sectors with <=7 states -- ARPACK k>=N-1 -- so tiny N=1 sectors
use dense eigvalsh, exactly as the runner's own full_gap/regression does), and (b) wires the
generic geometry's link tail/head lists into the generic hop builder.

OBSERVABLES (per the scope doc and the prior code audit):
  Delta_raw(N)  = E0(n_charged>=1) - E0(n_charged=0)   at t=0   [the conserved-sector meson gap]
  Delta_mir(N)  = Delta_raw - C3/beta,  C3 = 4/3
  t>0 leg: with hopping ON, n_charged is NOT conserved (the hop sends 3<->1 at a vertex), so the
    clean charged-vs-vacuum sector gap is ill-defined. We therefore report the FULL spectral gap
    e1-e0 of H(t) on the SAME basis, and define the hopping observable consistently as the full
    gap at t=0 vs t>0. softening = (g_full(t=0) - g_full(t))/g_full(t=0). This is labelled
    'full_gap' (charge-mixed) and is the SAME quantity tch_2plaq_extended_hop.py reports.
    NOTE: g_full(t=0) is generally NOT Delta_raw (its first excited state is the cheapest
    excitation of ANY sector, not the charged sector), so the softening is reported on its own
    footing -- see the report.
"""
from __future__ import annotations
import argparse
import functools
import importlib.util
import io
import contextlib
import json
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

ROOT = Path("/home/dave/octahedrons/python_code")
RUN_DIR = ROOT / "smg_dmrg_runs"
C3 = 4.0 / 3.0


def load(name, fname, silence=True):
    spec = importlib.util.spec_from_file_location(name, ROOT / fname)
    mod = importlib.util.module_from_spec(spec)
    if silence:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            spec.loader.exec_module(mod)
    else:
        spec.loader.exec_module(mod)
    return mod


def load_hop_defs():
    """Exec ONLY the setup+definitions of tch_2plaq_extended_hop.py (lines before the STEP-1
    driver banner at line 144). The file has no `if __name__` guard, so a plain import would run
    the full 106,460 t-scan (~6 min) every time. Truncating before STEP 1 runs the file's OWN
    initialization (load finish, the 336 oracle, robust_vbasis, the IDENTICAL sparse_hop/hop_tail/
    hop_head definitions) and stops -- giving byte-faithful validated functions, no driver."""
    src = (ROOT / "tch_2plaq_extended_hop.py").read_text().splitlines()
    cut = next(i for i, ln in enumerate(src)
               if ln.startswith("# ===") and "STEP 1: validate hop" in ln)
    body = "\n".join(src[:cut])
    ns = {"__name__": "_hop_defs_only", "__file__": str(ROOT / "tch_2plaq_extended_hop.py")}
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        exec(compile(body, "tch_2plaq_extended_hop.py[defs]", "exec"), ns)
    return ns  # ns["sparse_hop"], ns["hop_tail"], ns["hop_head"]


def sector_ground_robust(H_sector):
    """Lowest eigenvalue; dense for tiny sectors, sparse Lanczos otherwise (matches runner intent)."""
    k = H_sector.shape[0]
    if k <= 8:
        return float(np.min(np.linalg.eigvalsh(H_sector.toarray())))
    return float(np.min(eigsh(H_sector.tocsr(), k=min(6, k - 1), which="SA",
                              return_eigenvectors=False, maxiter=40000, tol=0)))


def full_gap(H):
    k = H.shape[0]
    if k <= 8:
        e = np.sort(np.linalg.eigvalsh(H.toarray()))
    else:
        e = np.sort(eigsh(H.tocsr(), k=min(6, k - 1), which="SA",
                          return_eigenvectors=False, maxiter=40000, tol=0))
    return float(e[1] - e[0])


def link_tail_head(geom):
    """tail/head vertex of each link on the 1xN strip (matches tch_2plaq_extended_hop tail_head)."""
    LINKS = geom["LINKS"]
    tails, heads = [], []
    for (kind, x, y) in LINKS:
        if kind == "h":
            tails.append((x, y)); heads.append((x + 1, y))
        else:  # "v"
            tails.append((x, 0)); heads.append((x, 1))
    return tails, heads


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--reps", choices=("minimal", "extended"), required=True)
    ap.add_argument("--betas", nargs="+", type=float, default=[0.5, 1.0])
    ap.add_argument("--ts", nargs="+", type=float, default=[0.0, 1.0])
    ap.add_argument("--tag", default="")
    args = ap.parse_args()

    R = load("runner", "smg_region2_nplaq_magnetic_gap.py")
    m = R.load_machinery()
    m.vbasis = functools.lru_cache(maxsize=None)(R.make_robust_vbasis(m.cg))
    cg = m.cg
    reps = m.REPS_MIN if args.reps == "minimal" else m.REPS_EXT
    hopns = load_hop_defs() if any(t != 0.0 for t in args.ts) else None

    geom = R.build_geometry(args.n)
    t0 = time.time()
    basis = R.build_basis(m, reps, geom)
    index = {b: i for i, b in enumerate(basis)}
    n = len(basis)
    print(f"[{args.n}plaq_{args.reps}] basis={n} built {time.time()-t0:.1f}s", flush=True)

    # --- magnetic Wilson (reused) ---
    Kc = functools.lru_cache(maxsize=None)(m.K_corner)
    sparse_wilson = R.make_sparse_wilson(m, Kc, geom)
    ta = time.time()
    Wsym = sum((sparse_wilson(p, reps, basis, index) for p in range(geom["n"])),
               sp.csr_matrix((n, n), dtype=complex))
    Wsym = (Wsym + Wsym.getH()).tocsr()
    ele, mat, cc = R.diag_arrays(m, basis)
    vac = np.where(cc == 0)[0]; chg = np.where(cc > 0)[0]
    print(f"  Wsym nnz={Wsym.nnz} ({time.time()-ta:.1f}s); vac={len(vac)} charged={len(chg)}", flush=True)
    cross = Wsym[vac][:, chg]
    cmax = float(abs(cross).max()) if cross.nnz else 0.0
    assert cmax < 1e-9, f"magnetic move couples charge sectors: {cmax:.2e}"

    # --- neutral hop (reused), only if t>0 requested ---
    HOP = None
    if hopns is not None:
        th = time.time()
        tails, heads = link_tail_head(geom)
        HOP = hopns["sparse_hop"](reps, basis, index, geom["VERTS"], geom["vlegs"], tails, heads)
        herm = abs((HOP + HOP.getH()) - (HOP + HOP.getH()).getH()).max() if HOP.nnz else 0.0
        assert herm < 1e-9, f"HOP+HOPd not Hermitian: {herm:.2e}"
        print(f"  HOP nnz={HOP.nnz} ({time.time()-th:.1f}s) herm_resid={herm:.1e}", flush=True)

    rows = []
    for b in args.betas:
        diag = (ele / b + mat).astype(complex)
        Hmag = sp.diags(diag) - (b / 2.0) * Wsym
        # t=0 charged-sector gap (the primary Delta_raw)
        tt = time.time()
        e_vac = sector_ground_robust(Hmag[vac][:, vac])
        e_chg = sector_ground_robust(Hmag[chg][:, chg])
        delta_raw = e_chg - e_vac
        delta_mir = delta_raw - C3 / b
        row = dict(n_plaq=args.n, rep_set=args.reps, beta=float(b), dim=n,
                   E0_vac=e_vac, E0_charged=e_chg,
                   delta_raw=delta_raw, delta_mir=delta_mir, C3_over_beta=C3 / b,
                   pred_2Dsmg_plus_C3=2 * m.DSMG + C3 / b,
                   full_gap_by_t={}, wall_seconds_t0=time.time() - tt)
        # full spectral gap vs t (the hopping-softening leg)
        for t in args.ts:
            ts = time.time()
            Ht = Hmag if (t == 0.0 or HOP is None) else (Hmag + t * (HOP + HOP.getH())).tocsr()
            g = full_gap(Ht)
            row["full_gap_by_t"][f"{t}"] = dict(full_gap=g, wall_seconds=time.time() - ts)
            print(f"  beta={b:.3f} t={t:.2f}: full_gap={g:.6f}  "
                  f"[Delta_raw(t0)={delta_raw:.6f} Delta_mir={delta_mir:.6f}]", flush=True)
        # softening from full gap (t=0 baseline vs each t>0)
        g0 = row["full_gap_by_t"].get("0.0", {}).get("full_gap")
        if g0:
            for t in args.ts:
                if t != 0.0:
                    gt = row["full_gap_by_t"][f"{t}"]["full_gap"]
                    row["full_gap_by_t"][f"{t}"]["softening_vs_t0"] = (g0 - gt) / g0
        rows.append(row)

    out = {"meta": dict(n_plaq=args.n, rep_set=args.reps, dim=n, betas=args.betas, ts=args.ts,
                        DSMG=m.DSMG, C3=C3, tag=args.tag),
           "rows": rows}
    fn = RUN_DIR / f"phase1_{args.n}plaq_{args.reps}{('_' + args.tag) if args.tag else ''}.json"
    fn.write_text(json.dumps(out, indent=2) + "\n")
    print(f"  wrote {fn}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
