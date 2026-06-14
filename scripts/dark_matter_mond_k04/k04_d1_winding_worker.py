#!/usr/bin/env python3
"""Winding-relic worker (Mechanism 1): canonical ramp -> trail census with
TOPOLOGICAL classification. A closed alternating trail is WINDING iff its
unwrapped displacement around the circuit is nonzero (it closes only through
the torus) — the content no contractible-healing process can remove.
PRE-REGISTERED primary observable: f_wind_B = occupied (B-side) bonds of
winding trails per site. Mechanism-1 gate: f_wind_B(496, L->inf) = 9 alpha0.
Usage: k04_d1_winding_worker.py L R seed -> one JSON line."""
import json
import math
import sys
import time
from pathlib import Path

L_ARG, R_ARG, SEED = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
src = (Path(__file__).parent / "k04_embedded_sweep.py").read_text()
prefix = src[: src.index('if START in ("hot", "ramp"):')]
sys.argv = ["k04", str(L_ARG), "1.7", "1.0", "0.5", "cold", "1", "0"]
env = {}
exec(compile(prefix, "k04_prefix", "exec"), env)
N, Lsz = env["N"], env["L"]
metropolis, census = env["metropolis"], env["census"]
wrap, add, AX = env["wrap"], env["add"], env["AX"]
B, c4, c6 = set(env["B"]), env["c4"], env["c6"]
env["random"].seed(100000 * L_ARG + 100 * R_ARG + SEED)

t0 = time.time()
B, c4, c6, _ = metropolis(B, c4, c6, math.inf, 20)
T_HI, T_LO = 6.0, 0.5
for s in range(R_ARG):
    Ts = T_HI * (T_LO / T_HI) ** (s / max(R_ARG - 1, 1))
    B, c4, c6, _ = metropolis(B, c4, c6, Ts, 1)
B, c4, c6, _ = metropolis(B, c4, c6, T_LO, max(20, R_ARG // 10))

# ---- trail census with winding classification (mirrors healing_spectrum) ----
best = None
for an in ((i, j, k) for i in (0, 1) for j in (0, 1) for k in (0, 1)):
    Tt = {(p, a) for (p, a) in ((s_, ax) for s_ in
          ((x, y, z) for x in range(Lsz) for y in range(Lsz) for z in range(Lsz))
          for ax in range(3)) if (p[a] - an[a]) % 2 == 0}
    D = B ^ Tt
    if best is None or len(D) < len(best):
        best = D
D = best
inc = {}
for bond in D:
    (p, a) = bond
    q = wrap(add(p, AX[a]))
    tag = "B" if bond in B else "T"
    inc.setdefault(p, []).append((bond, q, tag))
    inc.setdefault(q, []).append((bond, p, tag))
used, trails = set(), []
for v0 in list(inc):
    while True:
        start = next(((b, q, t) for (b, q, t) in inc.get(v0, []) if b not in used), None)
        if start is None:
            break
        ln = nB = 0
        disp = [0, 0, 0]
        cur, (b, q, t) = v0, start
        while True:
            used.add(b)
            ln += 1
            nB += (t == "B")
            (bp, ba) = b
            step = AX[ba] if cur == bp else (-AX[ba][0], -AX[ba][1], -AX[ba][2])
            disp = [disp[k] + step[k] for k in range(3)]
            cur, need = q, ("T" if t == "B" else "B")
            nxt = next(((bb, qq, tt) for (bb, qq, tt) in inc[cur]
                        if bb not in used and tt == need), None)
            if cur == v0 and ln % 2 == 0 and nxt is None:
                break
            b, q, t = nxt
        wind = any(d != 0 for d in disp)
        assert all(d % Lsz == 0 for d in disp)
        trails.append((ln, nB, wind))
assert sum(t[0] for t in trails) == len(D)
_, d = census(B)
f_line = sum(ln for ln, nB, w in trails if ln > 4) / N
f_wind_B = sum(nB for ln, nB, w in trails if w) / N
f_wind_tot = sum(ln for ln, nB, w in trails if w) / N
print(json.dumps(dict(L=L_ARG, R=R_ARG, seed=SEED, d_trap=round(d, 6),
                      f_line=round(f_line, 6), f_wind_B=round(f_wind_B, 6),
                      f_wind_tot=round(f_wind_tot, 6),
                      n_wind=sum(1 for *_, w in trails if w),
                      n_trails_gt4=sum(1 for ln, nB, w in trails if ln > 4),
                      secs=round(time.time() - t0, 1))))
