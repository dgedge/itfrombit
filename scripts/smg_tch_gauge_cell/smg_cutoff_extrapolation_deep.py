#!/usr/bin/env python3
r"""SMG/TCH continuum lift -- phase 1: cutoff(Lambda) + beta + O(1)-hopping sweep.

Extends the finite-cell mirror-gap evidence at the one-plaquette Z3 cell along the
three axes that the conditional continuum theorem needs a uniform lower bound over:
  * Lambda_cut  = states_per_charge in {2,3,4,5,6}   (local Fock truncation -> infinity)
  * beta        in {0.5,1,2,4}                        (pure Wilson axis, toward weak coupling)
  * hopping t   in {0, 0.5, 1.0}                      (O(1) -> a genuine gauge-dynamical test,
                                                       answering the script's own t-caveat)
The observable is the matter-resolved (mirror) gap: the gap to the first state whose
matter expectation exceeds the ground by >0.5 (a charged/mirror excitation). We then
extrapolate the gap to Lambda_cut -> infinity at each (beta,t) and report inf over the grid.

Reuses the committed one-plaquette Hamiltonian. numpy only; writes JSON.
"""
import importlib.util, json, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("plaq", HERE / "css_2d_mirror_fock_plaquette_lanczos.py")
plaq = importlib.util.module_from_spec(spec); spec.loader.exec_module(plaq)

SP   = [2, 3, 4, 5, 6]
BETA = [0.5, 1.0, 2.0, 4.0]
THOP = [0.0, 0.5, 1.0]

def mirror_gap(model, n_eigs=16, n_iter=110):
    evals, vecs = plaq.lanczos_lowest(model, n_iter=n_iter, n_eigs=n_eigs)
    g0 = model.matter_expectation(vecs[:, 0])
    full = float(evals[1] - evals[0])
    mir = float("inf")
    for i in range(1, vecs.shape[1]):
        if model.matter_expectation(vecs[:, i]) - g0 > 0.5:
            mir = float(evals[i] - evals[0]); break
    return full, mir

def main():
    results = []
    for sp in SP:
        for beta in BETA:
            for t in THOP:
                t0 = time.time()
                m = plaq.PlaquetteHamiltonian(beta=beta, hopping=t, states_per_charge=sp)
                full, mir = mirror_gap(m)
                dt = time.time() - t0
                row = dict(L=1, sp=sp, beta=beta, t=t, dim=int(m.dim),
                           full_gap=full, mirror_gap=(None if mir == float("inf") else mir), sec=dt)
                results.append(row)
                print(f"  sp={sp} beta={beta:<4} t={t:<4} dim={m.dim:<7} "
                      f"full={full:.4f} mirror={'--' if mir==float('inf') else f'{mir:.4f}'} ({dt:.1f}s)",
                      flush=True)
                json.dump(results, open(HERE / "smg_cutoff_sweep_results.json", "w"), indent=2)

    # cutoff extrapolation: mirror gap vs 1/sp at each (beta,t)
    print("\n  Lambda_cut -> infinity extrapolation (linear fit of mirror gap vs 1/sp):")
    extrap = {}
    for beta in BETA:
        for t in THOP:
            pts = [(1.0/r["sp"], r["mirror_gap"]) for r in results
                   if r["beta"] == beta and r["t"] == t and r["mirror_gap"] is not None]
            if len(pts) >= 3:
                x = np.array([p[0] for p in pts]); y = np.array([p[1] for p in pts])
                a, b = np.polyfit(x, y, 1)          # y = a*(1/sp) + b ; b = sp->inf intercept
                extrap[f"beta={beta},t={t}"] = float(b)
                print(f"    beta={beta:<4} t={t:<4}: gap(Lambda->inf) ~ {b:.4f}  (slope {a:+.3f})")
    if extrap:
        worst = min(extrap.values())
        print(f"\n  inf over grid of extrapolated mirror gap = {worst:.4f}")
        json.dump({"extrapolation": extrap, "inf_gap": worst},
                  open(HERE / "smg_cutoff_extrapolation.json", "w"), indent=2)
        print("  VERDICT:", "mirror gap stays OPEN as Lambda_cut->inf (positive uniform bound, Z3 proxy)"
              if worst > 0.2 else "gap closes -- continuum mirror pole NOT excluded")
    print("\nDONE.")

if __name__ == "__main__":
    main()
