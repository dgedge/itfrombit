#!/usr/bin/env python3
"""
Converged SMG verdict driver for the 3-4-5-0 model: both halves in one resume-safe scan.

For each (g, chi) it runs DMRG (conserve='q3450') and measures, on BOTH edges:
  * fermion correlation length xi  (edge A light: want GROWING/large = gapless;
                                    edge B mirror at g>g_c: want SMALL/saturating = gapped)
  * the 8 Dirac+Majorana mass correlators -> condensate check (any large-r PLATEAU on edge B
    = bilinear condensate = SSB, which KILLS the SMG reading).

SMG pass (at converged chi) =
  edge A xi_A keeps growing (gapless)  AND  edge B xi_B saturates small (mirror gapped)
  AND  no edge-B mass plateau (no condensate);  with g=0 edge B gapless as the control.

Resume-safe: each (g,chi) saved immediately to smg_dmrg_runs/converged_verdict.json; reruns skip
done points. Ctrl-C / power loss costs at most the one in-progress run.

PROGRESS HEARTBEAT: this driver prints at every phase -- DMRG start, DMRG done (with timing), and
each measurement distance as it completes -- so a slow point is visibly *advancing*, never silent.
If the distance counter is ticking, it is working. Only a counter frozen on one line for many
minutes is a real stall. (DMRG completing prints TeNPy's "final DMRG state not in canonical
form ..." line and possibly "Maximum truncation error exceeded" -- both are normal, not errors.)

Run on the box (backgrounded + logged so it survives SSH drops):
  cd ~/octahedrons/python_code
  nohup ~/tenpy-env/bin/python -m smg_dmrg.converged_verdict > ~/converged.log 2>&1 &
  tail -f ~/converged.log
Timing rule of thumb (L=16): DMRG ~7 min @ chi=256, ~25-30 min @ chi=512; measurement a few min more.
"""
import sys, os, time, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from smg_dmrg.model_3450 import (Model3450Spec, single_fermion_correlation_term,
                                 mass_correlation_term, MASS_PAIRS, centered_left_position)
from smg_dmrg.tenpy_backend import run_ground_dmrg_3450
from smg_dmrg.analyze_correlations import classify, MASS_KEYS

import logging
# Per-sweep DMRG visibility: surface TeNPy's sweep checkpoints (energy, Delta E, trunc) so the long
# DMRG phase ticks sweep-by-sweep instead of going dark; suppress the verbose config-read spam.
# (DMRG converges in a few sweeps but runs all max_sweeps because the critical edge keeps S moving.)
logging.basicConfig(level=logging.WARNING, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
logging.getLogger("tenpy.algorithms.dmrg").setLevel(logging.INFO)

L = 16
CHIS = (256, 512, 1024, 2048)  # 256/512 cached; push to 1024/2048 to RESOLVE the mirror-gap mechanism
GS = (6.5, 0.0)              # SMG-relevant g=6.5 FIRST, then g=0 control (free; g_c~5.7)
MAXD = 12


def measure_edge(psi, spec, edge, t0, label):
    md = min(MAXD, spec.edge_sites - 1)
    rs = list(range(1, md + 1))
    ferm, mass = [], {k: [] for k in MASS_KEYS}
    for d in rs:
        left = centered_left_position(spec.edge_sites, d)
        fv = [abs(psi.expectation_value_term(single_fermion_correlation_term(edge, fl, left, d), autoJW=True))
              for fl in ("3", "4", "5", "0")]
        ferm.append(float(np.mean(fv)))
        for kind in ("dirac", "majorana"):
            for a, b in MASS_PAIRS:
                mass[f"{kind}_{a}{b}"].append(
                    abs(psi.expectation_value_term(mass_correlation_term(edge, kind, (a, b), left, d), autoJW=True)))
        print(f"      [{label}] edge {edge}: distance {d}/{md} done  (+{time.time()-t0:.0f}s)", flush=True)
    return rs, ferm, mass


def xi_exp(rs, mags):
    rs = np.asarray(rs, float); lm = np.log(np.maximum(mags, 1e-300))
    A = np.vstack([rs, np.ones(len(rs))]).T
    slope, _ = np.linalg.lstsq(A, lm, rcond=None)[0]
    return (-1.0 / slope) if slope < 0 else float("inf")


def trunc(info):
    try:
        ss = info.get("sweep_statistics", {})
        for k in ("max_trunc_err", "max_trunc", "trunc_err"):
            if k in ss and len(ss[k]):
                return float(ss[k][-1])
    except Exception:
        pass
    return None


os.makedirs("smg_dmrg_runs", exist_ok=True)
RESULTS = "smg_dmrg_runs/converged_verdict.json"
done = {}
if os.path.exists(RESULTS):
    done = {tuple(r["key"]): r for r in json.load(open(RESULTS))}
    print(f"resuming: {len(done)} (g,chi) already saved", flush=True)


def save():
    json.dump(list(done.values()), open(RESULTS, "w"), indent=2)


print(f"L={L}  CHIS={CHIS}  GS={GS}  (heartbeat ON: a ticking distance counter = working)", flush=True)
print(f"  {'g':>4} {'chi':>5} {'chi_use':>7} {'energy':>11} {'trunc':>8} {'xiA':>6} {'xiB':>6} {'B-condensate':>16}  time", flush=True)
for g in GS:
    for chi in CHIS:
        key = (g, chi)
        if key in done:
            r = done[key]
            print(f"  {g:>4} {chi:>5} {r['chi_used']:>7} {r['energy']:>11.3f} {'':>8} "
                  f"{r['xiA']:>6.2f} {r['xiB']:>6.2f} {str(r['B_condensate']):>16}  (cached)", flush=True)
            continue
        t0 = time.time()
        label = f"g={g} chi={chi}"
        sp = Model3450Spec(unit_cells=L, t1=1.0, t2=0.5, g1=g, g2=g)
        print(f"  [{label}] DMRG start ...", flush=True)
        try:
            res = run_ground_dmrg_3450(sp, chi_max=chi, max_sweeps=40, svd_min=1e-12,
                                       conserve="q3450", max_trunc_err=None, combine=True)
            psi = res["psi"]; E = res["energy"]; te = trunc(res["dmrg_info"]); cu = max(psi.chi)
            print(f"  [{label}] DMRG done  E={E:.3f}  chi_used={cu}  "
                  f"trunc={(te if te is not None else float('nan')):.1e}  (+{time.time()-t0:.0f}s); "
                  f"measuring 2 edges x {MAXD} distances ...", flush=True)
            rA, fA, mA = measure_edge(psi, sp, "A", t0, label)
            rB, fB, mB = measure_edge(psi, sp, "B", t0, label)
            xiA, xiB = xi_exp(rA, fA), xi_exp(rB, fB)
            # Surface the large-r tail of every edge-B mass channel (the condensate candidates).
            # The discriminator is chi-STABILITY: a real condensate's tail is flat across chi; an
            # under-convergence artifact swings by orders of magnitude (e.g. dirac_40: 1e-12 @256
            # -> 3e-3 @512). Watch the top channel hold steady (-> SSB) or collapse (-> SMG) vs chi.
            tails = {k: float(np.mean(mB[k][-3:])) for k in MASS_KEYS}
            top3 = sorted(tails.items(), key=lambda kv: -kv[1])[:3]
            print(f"  [{label}] edge-B mass large-r tail (top3, watch chi-stability): " +
                  ", ".join(f"{k}={v:.2e}" for k, v in top3), flush=True)
            # condensate: any edge-B mass correlator that the verified classifier flags as plateau
            plat = [k for k in MASS_KEYS if classify(rB, mB[k]).get("plateau")]
            done[key] = {"key": [g, chi], "L": L, "chi_used": int(cu), "energy": float(E),
                         "trunc": (float(te) if te is not None else None),
                         "xiA": float(xiA), "xiB": float(xiB),
                         "B_condensate": plat or "none",
                         "edgeA_ferm": fA, "edgeB_ferm": fB,
                         "edgeB_mass": {k: mB[k] for k in MASS_KEYS}, "seconds": time.time() - t0}
            save()
            print(f"  {g:>4} {chi:>5} {cu:>7} {E:>11.3f} "
                  f"{(te if te is not None else float('nan')):>8.1e} {xiA:>6.2f} {xiB:>6.2f} "
                  f"{str(plat or 'none'):>16}  ({time.time()-t0:.0f}s)  <-- SAVED", flush=True)
        except Exception as exc:
            print(f"  {g:>4} {chi:>5}  FAILED: {exc}", flush=True)

print("""
READ (at the highest converged chi):
  SMG-favorable  <=>  g=6.5: xiA grows/large (light edge gapless) AND xiB small/saturating
                      (mirror gapped) AND B-condensate = none (no SSB);  g=0: both edges gapless.
  SMG-fails      <=>  xiA saturates (light edge gapped) OR a persistent edge-B condensate (SSB).
  Convergence: chi_used < chi_max and energy/xi flattening between chi steps => converged.
""", flush=True)
