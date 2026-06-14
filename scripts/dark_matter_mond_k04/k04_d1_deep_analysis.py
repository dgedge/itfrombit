#!/usr/bin/env python3
r"""D1' DEEP VERDICT — L-escalation of the protected-line deposition at the
canonical boot rate, against the zero-parameter gate
    f_line(R=496, L->inf) = 9 alpha0 e^{3/(2 phi)} = 0.743983.
Input: d1_deep_results.jsonl (the deep sweep).  Analysis: per-(L,R) means/SE;
1/L and 1/L^2 extrapolations at R=496; the gate verdict; R-trend at L=10.
exit 0 = analysis complete (verdict may be pass OR fail — pre-registered)."""
import json
import math
import statistics
import sys
from collections import defaultdict

ALPHA0 = 1 / 137.0
PHI = (math.sqrt(5) - 1) / 2
GATE = 9 * ALPHA0 * math.exp(3 / (2 * PHI))

rows = [json.loads(l) for l in open(sys.argv[1] if len(sys.argv) > 1 else "d1_deep_results.jsonl")]
g = defaultdict(list)
for r in rows:
    g[(r["L"], r["R"])].append(r["f_line"])
print(f"[1] PER-(L,R) protected-line density (n runs, mean +/- SE):")
print(f"    {'L':>3s} {'R':>5s} {'n':>3s} {'f_line':>9s} {'SE':>8s}")
stats = {}
for (L, R) in sorted(g):
    v = g[(L, R)]
    m = statistics.mean(v)
    se = statistics.stdev(v) / math.sqrt(len(v)) if len(v) > 1 else 0.0
    stats[(L, R)] = (m, se, len(v))
    print(f"    {L:>3d} {R:>5d} {len(v):>3d} {m:>9.5f} {se:>8.5f}")

print(f"\n[2] L-EXTRAPOLATION at the canonical rate R = 496:")
pts = [(L, stats[(L, 496)][0], stats[(L, 496)][1]) for L in (8, 10, 12, 16) if (L, 496) in stats]
# include L=6 from the laptop audit for context only (not in fits)
print(f"    laptop L=6 (context): 0.74815 +/- 0.05488")
for L, m, se in pts:
    print(f"    deep   L={L}: {m:.5f} +/- {se:.5f}")
def wfit(xs, ys, ws):
    sw = sum(ws); sx = sum(w*x for w, x in zip(ws, xs)); sy = sum(w*y for w, y in zip(ws, ys))
    sxx = sum(w*x*x for w, x in zip(ws, xs)); sxy = sum(w*x*y for w, x, y in zip(ws, xs, ys))
    den = sw*sxx - sx*sx
    b = (sw*sxy - sx*sy) / den
    a = (sy - b*sx) / sw
    var_a = sxx / den
    return a, b, math.sqrt(max(var_a, 0))
for nm, xf in (("1/L", lambda L: 1/L), ("1/L^2", lambda L: 1/L**2)):
    xs = [xf(L) for L, m, se in pts]
    ys = [m for L, m, se in pts]
    ws = [1/se**2 if se > 0 else 1.0 for L, m, se in pts]
    a, b, ea = wfit(xs, ys, ws)
    pull = (a - GATE) / ea if ea > 0 else float("inf")
    print(f"    {nm:>6s} extrapolation: f_line(inf) = {a:.5f} +/- {ea:.5f}"
          f"   gate {GATE:.5f}: {a/GATE-1:+.2%} ({pull:+.1f} sigma)")

print(f"\n[3] R-TREND at L = 10 (exclusivity check):")
for R in (350, 496, 700):
    if (10, R) in stats:
        m, se, n = stats[(10, R)]
        print(f"    R = {R}: f_line = {m:.5f} +/- {se:.5f}   x gamma_base = {m*math.exp(-3/(2*PHI)):.5f}"
              f"  (target 9a0 = {9*ALPHA0:.5f}, {m*math.exp(-3/(2*PHI))/(9*ALPHA0)-1:+.1%})")

print(f"\n[4] GATE: D1' requires f_line(496, inf) = {GATE:.6f} (zero free parameters).")
print("    Verdict per [2]; pre-registered: |pull| < 2 -> survives tighter gate;")
print("    |pull| > 3 with both extrapolants -> D1' refuted at deep precision.")
print("exit 0")
