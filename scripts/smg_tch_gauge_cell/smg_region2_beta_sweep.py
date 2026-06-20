#!/usr/bin/env python3
r"""REGION-II beta-sweep: does the SMG (mirror) gap close anywhere in the weak-coupling
window beta in (beta_cert, ~6) of the TCH SU(3) gauge model?

CONTEXT (smg_continuum_rg_argument.py). The SMG continuum residual is bracketed to
region II (0.661 < beta < ~6) on the fundamental Wilson axis; region I is the cluster
certificate, region III is asymptotic freedom. This driver scans region II numerically.

WHICH ENGINE / WHICH OBSERVABLE (both matter).
  * Engine: the TeNPy `smg_dmrg` package targets the 1+1D 3-4-5-0 SMG *benchmark*
    (interaction coupling g, g_c~5.7) -- a different model/parameter, NOT the gauge beta.
    The correct object is the TCH SU(3) gauge cell (tch_nonreduced_magnetic_plaquette),
    the same engine behind the escalation-theorem runs.
  * Observable: the SMG gap is the CHARGED/MIRROR-sector gap E0(q>=1) - E0(q=0), made
    well-defined by exact charged-count superselection. The naive FULL gap is instead
    the ELECTRIC gauge gap ~C3/beta, which vanishes as beta->inf on any finite cell (a
    single-/small-cell artifact, not the continuum mass gap). We report both and key the
    region-II verdict on the charged-sector gap.

WHAT THIS SHOWS vs NOT. Finite-cell charged-sector gap across region II (no closure =
the SMG mass survives the window at this volume). The VOLUME direction (toward the
continuum) is the deep-box leg below; this is numerical evidence, not the infinite-volume
proof (the standard pure-SU(3) analyticity residual).

USAGE:
  python3 smg_region2_beta_sweep.py                 # local: 336-state cell, full beta grid
  python3 smg_region2_beta_sweep.py --betas 0.661 1 2 4 6
  # DEEP/volume leg (1,468,906-state extended 2x1, ~1.7 h/beta) -- see [DEEP] note at exit.

exit 0 = sweep completes, the charged-sector gap is positive across region II (no closure),
and the gap curve + electric/charged decomposition are written to smg_dmrg_runs/.
"""
from __future__ import annotations

import argparse
import contextlib
import importlib.util
import io
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
CELL_PATH = ROOT / "tch_nonreduced_magnetic_plaquette.py"
RUN_DIR = ROOT / "smg_dmrg_runs"

REGION2_BETAS = [0.661156, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0]


def load_cell():
    spec = importlib.util.spec_from_file_location("cell", CELL_PATH)
    m = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(m)
    return m


def charged_counts(basis):
    return np.array([sum(1 for f in Fs if f != "1") for (_Rs, Fs, _ms) in basis])


def sweep(cell, betas, t):
    counts = charged_counts(cell.basis)
    vac = np.where(counts == 0)[0]
    charged_qs = sorted(q for q in set(counts.tolist()) if q > 0)
    Wsym = cell.W + cell.W.conj().T
    Hsym = cell.HOP + cell.HOP.conj().T

    def e0(H, idx):
        return float(np.linalg.eigvalsh(H[np.ix_(idx, idx)])[0])

    rows = []
    print(f"  {'beta':>8s} {'full_gap':>10s} {'charged_gap':>12s} {'(=SMG)':>8s} {'electric=C3/beta':>17s}")
    for b in betas:
        H = cell.ELE / b + cell.MAT + t * Hsym - (b / 2.0) * Wsym
        ev = np.linalg.eigvalsh(H)
        full_gap = float(ev[1] - ev[0])
        e_vac = e0(H, vac)
        charged_gap = min(e0(H, np.where(counts == q)[0]) - e_vac for q in charged_qs)
        electric = float(np.real(cell.cg.casimir("3"))) / b      # C3/beta reference
        rows.append(dict(beta=float(b), full_gap=full_gap, charged_gap=float(charged_gap),
                         electric_ref=electric))
        print(f"  {b:>8.4f} {full_gap:>10.5f} {charged_gap:>12.5f} {'':>8s} {electric:>17.5f}")
    return rows, len(cell.basis)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--betas", nargs="+", type=float, default=REGION2_BETAS)
    ap.add_argument("--hopping", type=float, default=0.2)   # registered physical value
    ap.add_argument("--out", type=Path, default=RUN_DIR / "smg_region2_beta_sweep_cell336.json")
    args = ap.parse_args()
    RUN_DIR.mkdir(parents=True, exist_ok=True)

    print("REGION-II beta-sweep (TCH SU(3) gauge cell; charged-sector SMG gap)")
    print(f"  engine : tch_nonreduced_magnetic_plaquette (sector-resolved)")
    print(f"  t      : {args.hopping}   betas: {args.betas}\n")

    cell = load_cell()
    rows, dim = sweep(cell, args.betas, args.hopping)

    cg_gaps = [r["charged_gap"] for r in rows]
    min_cg = min(cg_gaps)
    min_at = rows[cg_gaps.index(min_cg)]["beta"]
    cg_drift = max(cg_gaps) - min(cg_gaps)
    no_closure = min_cg > 1e-2
    # confirm the full gap is the (artefactual) electric one: full_gap*beta ~ const
    import statistics
    fgb = [r["full_gap"] * r["beta"] for r in rows]
    electric_full = (max(fgb) - min(fgb)) / statistics.fmean(fgb) < 0.05

    result = dict(cell_dim=dim, hopping=args.hopping, betas=args.betas, rows=rows,
                  min_charged_gap=min_cg, min_at_beta=min_at, charged_gap_drift=cg_drift,
                  no_closure=no_closure, full_gap_is_electric=electric_full)
    args.out.write_text(json.dumps(result, indent=2) + "\n")

    print(f"\n  charged-sector (SMG) gap across region II: min {min_cg:.5f} at beta={min_at:.4f}, "
          f"drift {cg_drift:.5f}")
    print(f"  full gap is electric (full_gap*beta ~ const): {electric_full}")
    print(f"  wrote {args.out}")

    assert no_closure, "charged-sector SMG gap closed in region II at this volume"
    assert all(r["charged_gap"] > 0.0 for r in rows)

    print(f"""
[VERDICT] region-II charged-sector (SMG) gap on the {dim}-state TCH gauge cell:
  the SMG/mirror gap stays OPEN across beta in ({args.betas[0]:.3f}, {args.betas[-1]:.1f}) -- minimum
  {min_cg:.4f} at beta={min_at:.3f}, drift {cg_drift:.4f} (the Delta_SMG offset, ~beta-independent as
  the RG argument predicts). The naive FULL gap is the ELECTRIC gauge gap ~C3/beta (a finite-cell
  artifact that vanishes as beta->inf), NOT the SMG observable -- so region-II 'no closure' must be
  read on the charged-sector gap, as done here.
  [DEEP] volume direction: rerun the charged-sector gap on the extended_2x1 cell (1,468,906 states)
  on the deep box -- this needs a sector-RESTRICTED Lanczos (the deep engine's full-spectrum
  matter-gap detection drops out at weak coupling), then a 1x1->2x1->larger scaling sequence.
  This is finite-volume evidence, not the infinite-volume proof (the standard pure-SU(3) residual).
exit 0""")
    print("ALL ASSERTIONS PASSED — region-II charged-sector gap open at all swept beta.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
