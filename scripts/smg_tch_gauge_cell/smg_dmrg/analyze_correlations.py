#!/usr/bin/env python3
"""
Verdict tool for the 3-4-5-0 SMG DMRG correlation output (cli.py measure_3450_correlations).

Reads one or more correlation JSON files (one per edge/g) and, for each fermion and mass
correlator, returns: exponential-vs-power-law classification (gapped vs gapless), the fitted
xi / eta, and a large-r PLATEAU flag (the condensate / SSB signature). Then it combines edges
into the SMG verdict.

Usage:
  ~/tenpy-env/bin/python python_code/smg_dmrg/analyze_correlations.py FILE [FILE ...]
  ~/tenpy-env/bin/python python_code/smg_dmrg/analyze_correlations.py 'runs/*.json'

What the SMG phase requires (the figure in arXiv:2202.12355):
  * EDGE A (light/physical): fermion correlators POWER-LAW at all g  -> chiral mode stays gapless.
  * EDGE B (mirror): below g_c power-law (gapless); above g_c EXPONENTIAL (gapped) with NO plateau.
  * CRUCIAL: no mass-correlator PLATEAU on either side of g_c -> no bilinear condensation (= SMG,
    not SSB). A nonzero large-r plateau in a mass correlator = condensate = SSB, NOT SMG.

HONEST LIMITS (read these):
  * Needs >= ~4 distances to classify a decay; 2 points classify nothing.
  * A single system size cannot prove a plateau is a true condensate vs a finite-size tail.
    Run >=2 unit_cells values; a real condensate's plateau is L-independent, a finite-size
    artifact shrinks with L. This tool flags 'possible condensate -> confirm with L-scaling'.
  * It reports the data honestly; it does not decide SMG from one file.

stdlib + numpy.
"""
from __future__ import annotations
import sys, json, glob
import numpy as np

FLOOR = 1e-11          # below this a correlator is treated as decayed-to-zero (noise)
R2_OK = 0.85           # min R^2 to accept a fit form
MASS_KEYS = [f"{k}_{a}{b}" for k in ("dirac", "majorana") for a, b in
             (("3", "5"), ("3", "0"), ("4", "5"), ("4", "0"))]


def mag(pair):  # [re, im] -> |z|
    return float(np.hypot(pair[0], pair[1]))


def linfit(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    A = np.vstack([x, np.ones(len(x))]).T
    (m, c), *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = m * x + c
    ss_res = float(np.sum((y - yhat) ** 2)); ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-30 else (1.0 if ss_res < 1e-30 else 0.0)
    return m, c, r2


def classify(rs, mags):
    """Return dict: form, xi/eta, R2s, plateau flag, given arrays of distance and |C|."""
    rs = np.asarray(rs, float); mags = np.asarray(mags, float)
    keep = mags > FLOOR
    rk, mk = rs[keep], mags[keep]
    if len(mk) < 4:
        if len(mags) >= 3 and np.all(mags <= FLOOR):
            return {"form": "decayed~0 (gapped)", "detail": "all |C| below noise floor", "plateau": False}
        return {"form": "INSUFFICIENT", "detail": f"only {len(mk)} usable points (need >=4)", "plateau": None}
    me, ce, r2e = linfit(rk, np.log(mk))                 # exponential: ln|C| vs r
    mp, cp, r2p = linfit(np.log(rk), np.log(mk))         # power-law:  ln|C| vs ln r
    xi = (-1.0 / me) if me < 0 else float("inf")
    eta = -mp
    if max(r2e, r2p) < R2_OK:
        form, detail = "unclear/noisy", f"R2_exp={r2e:.2f} R2_pow={r2p:.2f}"
    elif r2e >= r2p:
        form, detail = "EXPONENTIAL (gapped)", f"xi={xi:.2f}  R2={r2e:.3f} (vs pow {r2p:.3f})"
    else:
        form, detail = "POWER-LAW (gapless)", f"eta={eta:.2f}  R2={r2p:.3f} (vs exp {r2e:.3f})"
    # plateau / condensate heuristic (single L): the decay STALLS at a nonzero value.
    # Measure the relative drop across the last third: ~0 => plateau; exp ~large; power ~moderate.
    # (Confirm any plateau with L-scaling: a true condensate's plateau is L-independent.)
    tail_n = max(3, len(mk) // 3)
    tail = mk[-tail_n:]
    tail_drop = (1.0 - tail[-1] / tail[0]) if tail[0] > FLOOR else 1.0
    head_decayed = mk[0] > 3.0 * tail[-1]          # not a uniformly-flat (trivial) dataset
    plateau = bool(tail[-1] > 100 * FLOOR and tail_drop < 0.15 and head_decayed)
    detail += f"  [tail drop {tail_drop*100:.0f}%/last {tail_n}]"
    return {"form": form, "detail": detail, "plateau": plateau,
            "tail": float(mk[-1]), "tail_drop": float(tail_drop), "r2e": r2e, "r2p": r2p}


def analyze_file(path):
    d = json.load(open(path))
    edge, g1, g2, L = d.get("edge"), d.get("g1"), d.get("g2"), d.get("unit_cells")
    rs = [e["distance"] for e in d["distances"]]
    print(f"\n=== {path}")
    print(f"    model={d.get('model')} edge={edge} unit_cells={L} g1={g1} g2={g2} "
          f"n_distances={len(rs)} (r={min(rs)}..{max(rs)})")
    out = {"edge": edge, "g1": g1, "g2": g2, "L": L, "fermion": {}, "mass": {}}
    for group, keys in (("fermion", ("3", "4", "5", "0")), ("mass", MASS_KEYS)):
        print(f"    [{group}]")
        for key in keys:
            mags = [mag(e[group][key]) for e in d["distances"]]
            res = classify(rs, mags)
            out[group][key] = res
            flag = "  <<< PLATEAU/possible-condensate" if res.get("plateau") else ""
            print(f"      {key:<12} {res['form']:<22} {res.get('detail','')}{flag}")
    return out


def verdict(results):
    """Combine per-file results into an SMG/SSB/gapless read where possible."""
    print("\n" + "=" * 70 + "\nVERDICT\n" + "=" * 70)
    any_plateau = False
    for r in results:
        for grp in ("fermion", "mass"):
            for k, v in r[grp].items():
                if v.get("plateau"):
                    any_plateau = True
                    print(f"  CONDENSATE FLAG: edge {r['edge']} g1={r['g1']} {grp}:{k} "
                          f"plateaus (tail={v.get('tail'):.2e}) -> possible SSB, NOT SMG. "
                          f"Confirm with L-scaling.")
    edgeB = [r for r in results if r["edge"] == "B"]
    edgeA = [r for r in results if r["edge"] == "A"]
    if edgeB:
        forms = [v["form"] for r in edgeB for v in r["mass"].values()]
        exp = sum("EXPONENTIAL" in f or "decayed" in f for f in forms)
        pw = sum("POWER-LAW" in f for f in forms)
        print(f"  edge B (mirror) mass correlators: {exp} gapped / {pw} gapless / "
              f"{len(forms)-exp-pw} unclear  (SMG phase wants all gapped, no plateau)")
    if edgeA:
        forms = [v["form"] for r in edgeA for v in r["fermion"].values()]
        pw = sum("POWER-LAW" in f for f in forms)
        print(f"  edge A (light) fermion correlators: {pw}/{len(forms)} power-law "
              f"(SMG wants light edge GAPLESS = power-law)")
    if not edgeA:
        print("  NOTE: no edge-A file supplied -> cannot confirm the light edge stays gapless,")
        print("        which is half the SMG claim. Measure edge A too.")
    if not any_plateau:
        print("  No condensate plateau detected in the supplied files (necessary, not sufficient,")
        print("  for SMG -- still need: full-distance fits, edge A gapless, and an L-scaling check).")


def main(argv):
    paths = []
    for a in argv:
        paths.extend(sorted(glob.glob(a)) or [a])
    if not paths:
        print(__doc__); raise SystemExit("supply one or more correlation JSON files")
    results = [analyze_file(p) for p in paths]
    verdict(results)


if __name__ == "__main__":
    main(sys.argv[1:])
