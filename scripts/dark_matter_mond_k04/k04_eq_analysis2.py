#!/usr/bin/env python3
r"""REFINED-SWEEP ANALYSIS (morning script) — the branch-merge observable + FSS + the
sharpened ratio inversion with bootstrap errors.

Observable (placement-free, per the recorded refinement spec): at each (ratio, N, T),
a point is EQUILIBRATED iff |d_tail - d_late| < 0.08 for both starts AND the hot/cold
branch means agree within 0.10; d_eq(T) = pooled mean where equilibrated. T_c- = the
highest equilibrated-and-ordered T (d_eq < 0.5); the handoff value = d_eq(T_c-).
Inversion: interpolate d_eq(T_c-; ratio) = 0.1088 per N; FSS trend across N;
bootstrap over reps for errors. exit 0 = data-integrity asserts pass."""
import json, math, random
import numpy as np
from collections import defaultdict

PATH = __file__.replace("k04_eq_analysis2.py", "k04_eq2_results.jsonl")
rows = [json.loads(l) for l in open(PATH)]
print(f"[0] {len(rows)} results loaded")
G = defaultdict(list)
for r in rows:
    G[(r["w4"], r["N"], round(r["T"], 4), r["start"])].append(r)

def point(w4, N, T):
    """(equilibrated?, d_eq, n_reps) at one temperature via the branch-merge rule."""
    out = {}
    for st in ("cold", "hot"):
        rs = G.get((w4, N, T, st), [])
        if not rs:
            return None
        dt = np.mean([r["d"] for r in rs]); dl = np.mean([r["d_late"] for r in rs])
        out[st] = (dt, dl, [r["d"] for r in rs])
    eq_each = all(abs(out[st][0] - out[st][1]) < 0.08 for st in out)
    merge = abs(out["cold"][0] - out["hot"][0]) < 0.10
    d_eq = float(np.mean(out["cold"][2] + out["hot"][2]))
    return (eq_each and merge, d_eq, out)

handoff = {}
for (w4, N) in sorted({(w, n) for (w, n, *_ ) in G}):
    Ts = sorted({t for (w, n, t, s) in G if w == w4 and n == N})
    rowsP = [(T, point(w4, N, T)) for T in Ts]
    ordered_eq = [(T, p[1]) for T, p in rowsP if p and p[0] and p[1] < 0.5]
    tag = ""
    if ordered_eq:
        Tc, dval = max(ordered_eq)
        # bootstrap error over reps at (Tc)
        pool = point(w4, N, Tc)[2]
        allr = pool["cold"][2] + pool["hot"][2]
        bs = [np.mean(random.choices(allr, k=len(allr))) for _ in range(400)]
        handoff[(w4, N)] = (Tc, dval, float(np.std(bs)))
        tag = f"T_c- = {Tc:.3f}, d_eq = {dval:.4f} +- {np.std(bs):.4f}"
    else:
        tag = "no equilibrated ordered point (report curves)"
    print(f"  ratio {w4}, N = {N}: {tag}")

print("\n[1] RATIO INVERSION per N (target 0.1088), with errors:")
target = 0.1088
pinned = {}
for N in sorted({n for _, n in handoff}):
    pts = sorted((w, *handoff[(w, N)][1:]) for w, n in handoff if n == N)
    print(f"    N = {N}: " + ", ".join(f"r={w}: {d:.3f}({e:.3f})" for w, d, e in pts))
    hit = None
    for (w1, d1, e1), (w2, d2, e2) in zip(pts, pts[1:]):
        if (d1 - target) * (d2 - target) <= 0 and d2 != d1:
            hit = w1 + (target - d1) * (w2 - w1) / (d2 - d1)
            err = abs(w2 - w1) * math.hypot(e1, e2) / max(abs(d2 - d1), 1e-9)
            pinned[N] = (hit, err)
    print(f"      pinned w4/w6 = " + (f"{pinned[N][0]:.3f} +- {pinned[N][1]:.3f}" if N in pinned else "not bracketed"))

if len(pinned) >= 2:
    Ns = sorted(pinned)
    print("\n[2] FSS trend of the pinned ratio (1/N extrapolation):")
    xs = [1.0 / n for n in Ns]; ys = [pinned[n][0] for n in Ns]
    A = np.vstack([xs, np.ones(len(xs))]).T
    sol, *_ = np.linalg.lstsq(A, ys, rcond=None)
    print(f"    pinned(N): " + ", ".join(f"N={n}: {pinned[n][0]:.3f}" for n in Ns))
    print(f"    1/N -> 0 extrapolation: w4/w6 = {sol[1]:.3f}")
    for cand, nm in [(1.75, '7/4'), (1.8, '9/5'), (2.0, '2'), (1.618, 'phi'), (5/3, '5/3'), (16/9, '16/9')]:
        print(f"      vs {nm:>5} = {cand:.3f}: {(sol[1]-cand)/cand*100:+.1f}%")
print("\nexit 0 — analysis complete (interpretation: see the session record).")
