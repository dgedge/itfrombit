#!/usr/bin/env python3
r"""MECHANISM-1 VERDICT: is the billable relic the winding sub-network?
Pre-registered primary: f_wind_B(R=496, L->inf) vs 9 alpha0 = 0.065693.
Requirements: (i) intensive in L (unlike f_line, which grew); (ii) the
extrapolant lands the gate. Secondary: R-trend at L=10; winding-trail counts.
exit 0 = analysis complete (verdict may be pass OR fail)."""
import json
import math
import statistics
import sys
from collections import defaultdict

ALPHA0 = 1 / 137.0
TARGET = 9 * ALPHA0
rows = [json.loads(l) for l in open(sys.argv[1] if len(sys.argv) > 1 else "d1_winding_results.jsonl")]
g = defaultdict(list)
for r in rows:
    g[(r["L"], r["R"])].append(r)
print("[1] PER-(L,R): winding observables (mean +/- SE over seeds):")
print(f"    {'L':>3s} {'R':>5s} {'n':>3s} {'f_line':>8s} {'f_wind_B':>9s} {'SE':>8s} {'n_wind':>7s}")
stats = {}
for (L, R) in sorted(g):
    v = g[(L, R)]
    fw = [x["f_wind_B"] for x in v]
    m = statistics.mean(fw)
    se = statistics.stdev(fw) / math.sqrt(len(fw)) if len(fw) > 1 else 0.0
    stats[(L, R)] = (m, se, len(v))
    print(f"    {L:>3d} {R:>5d} {len(v):>3d} {statistics.mean(x['f_line'] for x in v):>8.4f} "
          f"{m:>9.5f} {se:>8.5f} {statistics.mean(x['n_wind'] for x in v):>7.1f}")

print(f"\n[2] INTENSIVITY + GATE at R = 496 (target 9 alpha0 = {TARGET:.5f}):")
pts = [(L,) + stats[(L, 496)][:2] for L in (8, 10, 12, 16) if (L, 496) in stats]
for L, m, se in pts:
    print(f"    L = {L:>2d}: f_wind_B = {m:.5f} +/- {se:.5f}   ratio to gate = {m/TARGET:.3f}")
def wfit(xs, ys, ws):
    sw = sum(ws); sx = sum(w*x for w, x in zip(ws, xs)); sy = sum(w*y for w, y in zip(ws, ys))
    sxx = sum(w*x*x for w, x in zip(ws, xs)); sxy = sum(w*x*y for w, x, y in zip(ws, xs, ys))
    den = sw*sxx - sx*sx
    b = (sw*sxy - sx*sy) / den
    a = (sy - b*sx) / sw
    return a, b, math.sqrt(max(sxx/den, 0))
for nm, xf in (("1/L", lambda L: 1/L), ("1/L^2", lambda L: 1/L**2)):
    xs = [xf(L) for L, m, se in pts]; ys = [m for L, m, se in pts]
    ws = [1/se**2 if se > 0 else 1.0 for L, m, se in pts]
    a, b, ea = wfit(xs, ys, ws)
    pull = (a - TARGET)/ea if ea > 0 else float("inf")
    print(f"    {nm:>6s} extrapolation: f_wind_B(inf) = {a:.5f} +/- {ea:.5f}"
          f"  -> gate: {a/TARGET-1:+.1%} ({pull:+.1f} sigma)")

print(f"\n[3] R-TREND at L = 10 (KZ winding scaling):")
for R in (248, 496, 992):
    if (10, R) in stats:
        m, se, n = stats[(10, R)]
        print(f"    R = {R:>4d}: f_wind_B = {m:.5f} +/- {se:.5f}")
print("\n[4] verdicts are read from [2]: intensive + lands -> Mechanism 1 corroborated;")
print("    intensive + misses -> the winding density is the measured boundary condition")
print("    for whatever mechanism replaces it; non-intensive -> mechanism mis-posed.")
print("exit 0")
