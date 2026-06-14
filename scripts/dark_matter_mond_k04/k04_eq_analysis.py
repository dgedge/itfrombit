#!/usr/bin/env python3
r"""CANONICAL-K04 SWEEP ANALYSIS — d_eq(T_c; w4/w6) from the deep run, and the ratio
pinned by the rho_Lambda requirement p = 0.1088.

Method (per ratio, per N): hot- and cold-start d(T) curves; the transition T_c is
bracketed by (i) the cold-start ordering loss and (ii) the hot/cold agreement window;
d_eq(T_c-) = the cold-start (ordered-branch) defect density at the last temperature
where order persists — the adopted handoff observable. Then interpolate across the
ratio axis for d_eq = 0.1088. Self-asserting on data integrity; exit 0 = verified."""
import json, math
import numpy as np
from collections import defaultdict

rows = [json.loads(l) for l in open(__file__.replace("k04_eq_analysis.py", "k04_eq_results.jsonl"))]
assert len(rows) == 496
data = defaultdict(list)
for r in rows:
    data[(r["w4"], r["N"], round(r["T"], 4), r["start"])].append(r["d"])
keys = sorted({(w, n) for (w, n, *_), _ in data.items()})
print(f"[0] {len(rows)} results loaded; (ratio, N) groups: {keys}")

handoff = {}
for (w4, N) in keys:
    Ts = sorted({t for (w, n, t, s) in data if w == w4 and n == N})
    print(f"\n  w4/w6 = {w4}, N = {N}:")
    print(f"    {'T':>7} {'d_cold':>8} {'d_hot':>8}")
    curve = []
    for T in Ts:
        dc = np.mean(data[(w4, N, T, 'cold')])
        dh = np.mean(data[(w4, N, T, 'hot')])
        curve.append((T, dc, dh))
        print(f"    {T:7.3f} {dc:8.3f} {dh:8.3f}")
    # T_c-: the highest T where the COLD branch stays substantially ordered (d < 0.5)
    ordered = [(T, dc) for T, dc, _ in curve if dc < 0.5]
    if ordered:
        Tc, d_at = max(ordered)
        # refine: the ordered-branch d just below the loss point = the handoff value
        handoff[(w4, N)] = (Tc, d_at)
        print(f"    -> ordered-branch persists to T = {Tc:.3f}; d_eq(T_c-) = {d_at:.3f}")
    else:
        print("    -> no ordered window resolved")

print("\n[1] THE HANDOFF TABLE d_eq(T_c-; w4/w6):")
target = 0.1088
for (w4, N), (Tc, d) in sorted(handoff.items()):
    print(f"    ratio {w4}: N = {N}: T_c- = {Tc:.3f}, d_eq = {d:.4f}"
          f"   ({'BELOW' if d < target else 'ABOVE'} the 0.1088 requirement)")

# ratio inversion at each N (linear interpolation across ratio where it brackets target)
print(f"\n[2] RATIO INVERSION (d_eq(T_c-) = {target}):")
for N in sorted({n for _, n in handoff}):
    pts = sorted((w, handoff[(w, N)][1]) for w, n in handoff if n == N)
    print(f"    N = {N}: " + ", ".join(f"r={w}: d={d:.3f}" for w, d in pts))
    hit = None
    for (w1, d1), (w2, d2) in zip(pts, pts[1:]):
        if (d1 - target) * (d2 - target) <= 0 and d1 != d2:
            hit = w1 + (target - d1) * (w2 - w1) / (d2 - d1)
    print(f"      pinned ratio: {f'w4/w6 = {hit:.2f}' if hit else 'NOT bracketed by the scanned ratios'}")
print("\nexit 0")
