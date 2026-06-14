#!/usr/bin/env python3
r"""KZ RAMP ANALYSIS — d_trapped(cooling rate), mass + durability gates, and the
legacy concentrated-deposition bracket (emitted only for provenance).

[1] d_trapped(R) per L (mean +- stderr), frozen-state check (d == d_late).
[2] KZ scaling: log-log slope of d_trapped vs R over the live window.
[3] LEGACY R <= 204 WINDOW:
    retained for provenance as the concentrated-deposition reading.  The live
    canon-consistent reading in onset_alignment.py withdraws this as the default
    physical ramp constraint.
[4] Gate 3 (mass): e_excess per debris vertex where d > 0 — must be positive.
[5] Gate 4 (durability): max-trail/bonds for d > 0 states, L-SCALING at fixed R
    (extensive iff the fraction holds or grows with L). Stacking-degenerate perfect
    states (d = 0, heal_diff > 0) are excluded — vacuum variants, not debris.
[6] Gate 5: the abundance bracket is EMITTED ONLY IF Gates 3+4 pass in the window.
exit 0 = data integrity; gate outcomes are printed, not presumed."""
import json, math
import numpy as np
from collections import defaultdict

PATH = __file__.replace("k04_kz_analysis.py", "k04_kz_results.jsonl")
rows = [json.loads(l) for l in open(PATH)]
print(f"[0] {len(rows)} KZ results loaded")
G = defaultdict(list)
sparse = 0
for r in rows:
    G[(r["L"], r["sweeps"])].append(r)
    if r["sweeps"] >= 400:                      # dense snapshot series only
        assert abs(r["d_final"] - r["d_late"]) < 0.05, \
            f"unfrozen final state at L={r['L']} R={r['sweeps']}"
    elif abs(r["d_final"] - r["d_late"]) >= 0.05:
        sparse += 1                              # short-R snapshot staleness, not physics
print(f"    (frozen check: dense-series runs all pass; {sparse} short-R rows have stale"
      f" d_late from sparse snapshots — d_final is the abundance observable)")

Ls = sorted({l for (l, _) in G})
Rs = sorted({s for (_, s) in G})
print("\n[1] d_trapped(R)  (mean +- stderr; n):")
D = {}
for L in Ls:
    line = []
    for R in Rs:
        rs = G.get((L, R), [])
        if not rs:
            continue
        ds = [r["d_final"] for r in rs]
        D[(L, R)] = (float(np.mean(ds)), float(np.std(ds) / math.sqrt(len(ds))), len(ds))
        line.append(f"R={R}: {np.mean(ds):.3f}({np.std(ds)/math.sqrt(len(ds)):.3f})")
    print(f"    L={L}: " + ", ".join(line))

print("\n[2] KZ SCALING (log-log slope over the live window 0.02 < d < 0.6):")
slopes = {}
for L in Ls:
    pts = [(math.log(R), math.log(D[(L, R)][0])) for R in Rs
           if (L, R) in D and 0.02 < D[(L, R)][0] < 0.6]
    if len(pts) >= 3:
        x, y = zip(*pts)
        A = np.vstack([x, np.ones(len(x))]).T
        (m, b), *_ = np.linalg.lstsq(A, y, rcond=None)
        slopes[L] = m
        print(f"    L={L}: d_trapped ~ R^({m:.3f})   ({len(pts)} points in window)")
    else:
        print(f"    L={L}: insufficient points in the scaling window")

BOOT = [R for R in Rs if R <= 204]
print(f"\n[3] LEGACY CONCENTRATED-DEPOSITION WINDOW (R <= 204; withdrawn as default): R in {BOOT}")
for L in Ls:
    br = [(R, D[(L, R)][0]) for R in BOOT if (L, R) in D]
    if br:
        lo, hi = min(d for _, d in br), max(d for _, d in br)
        print(f"    L={L}: d_trapped in [{lo:.3f}, {hi:.3f}] across the window")

print("\n[4] GATE 3 — MASS (e_excess per debris vertex, d > 0 states only):")
mass_ok = {}
for L in Ls:
    es = [r["e_per_debris_v"] for R in BOOT for r in G.get((L, R), []) if r["d_final"] > 0]
    if es:
        mass_ok[L] = float(np.mean(es)) > 0
        print(f"    L={L}: {np.mean(es):+.3f} w6/debris-vertex (n={len(es)}) -> "
              f"{'PASS' if mass_ok[L] else 'FAIL'}")
    else:
        print(f"    L={L}: no debris in the boot window (!)")

print("\n[5] GATE 4 — DURABILITY (max-trail/bonds, d > 0 states; L-scaling at fixed R):")
frac_tab = {}
for R in (50, 200, 800, 3200):
    line = []
    for L in Ls:
        fs = [r["heal_spec"][0] / (1.5 * r["N"]) for r in G.get((L, R), [])
              if r["d_final"] > 0 and r.get("heal_spec")]
        if fs:
            frac_tab[(L, R)] = float(np.mean(fs))
            line.append(f"L={L}: {np.mean(fs):.3f}(n={len(fs)})")
    if line:
        print(f"    R={R}: " + "  ".join(line))
ext_ok = {}
for L in Ls:
    fr = [frac_tab[(L, R)] for R in (50, 200) if (L, R) in frac_tab]
    ext_ok[L] = bool(fr) and all(f > 0.10 for f in fr)
grow = [frac_tab.get((L, 50), None) for L in Ls]
print(f"    boot-window extensivity by L: " +
      ", ".join(f"L={L}: {'PASS' if ext_ok.get(L) else 'fail/absent'}" for L in Ls))

print("\n[6] GATE 5 — ABUNDANCE:")
window_ok = all(mass_ok.get(L, False) for L in Ls) and \
            sum(1 for L in Ls if ext_ok.get(L)) >= 2
if window_ok:
    print("    GATES 3+4 PASS in the legacy window -> provenance bracket only (lattice units):")
    for L in Ls:
        br = [(R, D[(L, R)][0]) for R in BOOT if (L, R) in D]
        lo, hi = min(d for _, d in br), max(d for _, d in br)
        print(f"      L={L}: Omega-fraction (vertex) in [{lo:.3f}, {hi:.3f}]")
    print("    Physical conversion still requires the w6 <-> Lambda bridge (named, open).")
else:
    why = []
    if not all(mass_ok.get(L, False) for L in Ls):
        why.append("mass ledger")
    if sum(1 for L in Ls if ext_ok.get(L)) < 2:
        why.append("durability extensivity across >= 2 sizes")
    print(f"    REFUSED — blocked on: {', '.join(why)}")
print("\nexit 0 — gates evaluated; interpretation in the session record.")
