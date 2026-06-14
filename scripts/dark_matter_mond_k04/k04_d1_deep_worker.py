#!/usr/bin/env python3
"""D1' deep worker: canonical boot ramp -> protected-line census.
Usage: k04_d1_deep_worker.py L R seed   (one JSON line to stdout)
Protocol: geometric ramp T 6.0 -> 0.5 over R sweeps + max(20, R//10) hold.
Observables: d_trapped, |D|/N, f_4 (healable), f_line (protected, trails > 4).
D1' gate: f_line(R=496, L->inf) = 9*alpha0*exp(3/(2*phi)) = 0.743983."""
import json
import math
import sys
import time
from pathlib import Path

L_ARG, R_ARG, SEED = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
src = (Path(__file__).parent / "k04_embedded_sweep.py").read_text()
marker = 'if START in ("hot", "ramp"):'
prefix = src[: src.index(marker)]
sys.argv = ["k04", str(L_ARG), "1.7", "1.0", "0.5", "cold", "1", "0"]
env = {}
exec(compile(prefix, "k04_prefix", "exec"), env)
N = env["N"]
metropolis, census, healing = env["metropolis"], env["census"], env["healing_spectrum"]
B, c4, c6 = set(env["B"]), env["c4"], env["c6"]
env["random"].seed(100000 * L_ARG + 100 * R_ARG + SEED)

t0 = time.time()
B, c4, c6, _ = metropolis(B, c4, c6, math.inf, 20)
T_HI, T_LO = 6.0, 0.5
for s in range(R_ARG):
    Ts = T_HI * (T_LO / T_HI) ** (s / max(R_ARG - 1, 1))
    B, c4, c6, _ = metropolis(B, c4, c6, Ts, 1)
B, c4, c6, _ = metropolis(B, c4, c6, T_LO, max(20, R_ARG // 10))
_, d = census(B)
lengths, _, nD = healing(B)
f_line = sum(l for l in lengths if l > 4) / N
f_4 = sum(l for l in lengths if l == 4) / N
print(json.dumps(dict(L=L_ARG, R=R_ARG, seed=SEED, N=N, d_trap=round(d, 6),
                      fD=round(nD / N, 6), f4=round(f_4, 6), f_line=round(f_line, 6),
                      n_trails_gt4=sum(1 for l in lengths if l > 4),
                      max_trail=max(lengths) if lengths else 0,
                      secs=round(time.time() - t0, 1))))
