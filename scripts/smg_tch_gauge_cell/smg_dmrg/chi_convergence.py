#!/usr/bin/env python3
"""
chi-convergence + free-edge gate for the 3-4-5-0 SMG DMRG.

Decides whether "edge A (light) is gapped" in the g-scan is a finite-chi artifact or real.

Part 1 (exact, fast): free (g=0) single-particle half-fill gap vs L. Gapless chiral edges
  => gap ~ 1/L -> 0. (Confirmed: gap ~ 3.0/L, model IS gapless.)
Part 2 (DMRG): edge-A fermion correlation length xi_A vs chi, at g=0 (control: DMRG should
  capture the now-proven-gapless free edge) and at g=6.5 (the interacting point). xi_A GROWING
  with chi => gapless captured (the scan's "gap" was finite-chi); xi_A SATURATING at small
  value => genuine gap OR the chiral-MPS wall (a net-chiral c_-!=0 edge cannot be an MPS GS).

Run (from repo root):
  ~/tenpy-env/bin/python python_code/smg_dmrg/chi_convergence.py
"""
import sys, os, time, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from smg_dmrg.model_3450 import (Model3450Spec, free_spectrum_3450,
                                 single_fermion_correlation_term, centered_left_position)
from smg_dmrg.tenpy_backend import run_ground_dmrg_3450

L = 16
CHIS = (256, 512, 1024)
GS = (0.0, 6.5)


def edgeA_xi(psi, spec, maxd=10):
    md = min(maxd, spec.edge_sites - 1)
    rs, mags = [], []
    for d in range(1, md + 1):
        left = centered_left_position(spec.edge_sites, d)
        vals = [abs(psi.expectation_value_term(single_fermion_correlation_term("A", fl, left, d), autoJW=True))
                for fl in ("3", "4", "5", "0")]
        rs.append(d); mags.append(float(np.mean(vals)))
    rs = np.array(rs, float); lm = np.log(np.maximum(mags, 1e-300))
    A = np.vstack([rs, np.ones(len(rs))]).T
    slope, _ = np.linalg.lstsq(A, lm, rcond=None)[0]
    return ((-1.0 / slope) if slope < 0 else float("inf")), mags


def trunc(info):
    try:
        ss = info.get("sweep_statistics", {})
        for k in ("max_trunc_err", "max_trunc", "trunc_err"):
            if k in ss and len(ss[k]):
                return float(ss[k][-1])
    except Exception:
        pass
    return None


print("PART 1 — exact free (g=0) half-fill gap vs L (gapless => ~1/L -> 0):", flush=True)
for LL in (8, 12, 16, 20, 24):
    sp = Model3450Spec(unit_cells=LL, t1=1.0, t2=0.5, g1=0.0, g2=0.0)
    ev, herm, gap = free_spectrum_3450(sp)
    print(f"  L={LL:>2}  gap={gap:.5f}  (gap*L={gap*LL:.2f})", flush=True)

print(f"\nPART 2 — DMRG edge-A xi vs chi, L={L}:", flush=True)
# Interruption-robust: each completed (g,chi) is saved immediately and skipped on resume,
# so a power cut / kill costs at most the one in-progress run, never the whole scan.
os.makedirs("smg_dmrg_runs", exist_ok=True)
RESULTS = "smg_dmrg_runs/chi_conv_results.json"
done = {}
if os.path.exists(RESULTS):
    done = {tuple(r["key"]): r for r in json.load(open(RESULTS))}
    print(f"  resuming: {len(done)} (g,chi) already saved -> will skip them", flush=True)


def save_results():
    json.dump(list(done.values()), open(RESULTS, "w"), indent=2)


print(f"  {'g':>4} {'chi':>5} {'chi_used':>8} {'energy':>12} {'trunc':>9} {'xi_A':>7}  time", flush=True)
for g in GS:
    for chi in CHIS:
        key = (g, chi)
        if key in done:
            r = done[key]
            print(f"  {g:>4} {chi:>5} {r['chi_used']:>8} {r['energy']:>12.4f} "
                  f"{'':>9} {r['xi_A']:>7.2f}  (cached)", flush=True)
            continue
        t0 = time.time()
        sp = Model3450Spec(unit_cells=L, t1=1.0, t2=0.5, g1=g, g2=g)
        try:
            res = run_ground_dmrg_3450(sp, chi_max=chi, max_sweeps=30, svd_min=1e-12,
                                       conserve="q3450", max_trunc_err=None, combine=True)
            psi = res["psi"]; E = res["energy"]; te = trunc(res["dmrg_info"]); cu = max(psi.chi)
            xi, mags = edgeA_xi(psi, sp)
            done[key] = {"key": [g, chi], "L": L, "chi_used": int(cu), "energy": float(E),
                         "trunc": (float(te) if te is not None else None), "xi_A": float(xi),
                         "mags": mags, "seconds": time.time() - t0}
            save_results()   # persist BEFORE printing, so it's on disk even if killed next
            print(f"  {g:>4} {chi:>5} {cu:>8} {E:>12.4f} "
                  f"{(te if te is not None else float('nan')):>9.1e} {xi:>7.2f}  ({time.time()-t0:.0f}s)", flush=True)
        except Exception as exc:
            print(f"  {g:>4} {chi:>5}  FAILED: {exc}", flush=True)

print("""
READ:
  g=0 control: if xi_A GROWS with chi, DMRG captures the (proven-gapless) free edge -> the scan's
    "edge A gapped" was a finite-chi artifact, and the scan needs re-running at higher chi.
    If xi_A SATURATES even at g=0, DMRG cannot represent the net-chiral gapless edge (c_-!=0
    chiral-MPS wall) -> correlator decay will NEVER show the light edge gapless; switch to the
    entanglement-spectrum / chiral-central-charge diagnostic.
  g=6.5: compare to g=0 to see whether the interaction changes the light edge.
""", flush=True)
