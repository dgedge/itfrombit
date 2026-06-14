#!/usr/bin/env python3
r"""CANONICAL-K04 EQUILIBRIUM SWEEP WORKER — one (N, w4, w6, T, start, rep) point.

Measures the equilibrium defect fraction d_eq under the ADOPTED canonical protocol
(ANCHOR 1.1 banner, 2026-06-11): Kempe-move cubic stratum, Metropolis at fixed T,
hot (random cubic) and cold (perfect Q3 tiling) starts, d time-averaged over the
final third of the run. Emits one JSON line.
Usage: k04_eq_sweep.py N w4 w6 T start sweeps rep
Exactness: incremental C4+C6 with drift checks (asserts)."""
import sys, json, random, math
import numpy as np

N = int(sys.argv[1]); W4 = float(sys.argv[2]); W6 = float(sys.argv[3])
T = float(sys.argv[4]); START = sys.argv[5]; SWEEPS = int(sys.argv[6]); REP = int(sys.argv[7])
seed = hash((N, W4, W6, round(T * 1000), START, REP)) % (2**31)
random.seed(seed); np.random.seed(seed % 65536)

def adj_from(E, n):
    A = [[] for _ in range(n)]
    for a, b in E:
        A[a].append(b); A[b].append(a)
    return A

def _paths(adj, start, goal, length):
    out, stack = [], [(start, (start,))]
    while stack:
        v, path = stack.pop()
        if len(path) == length:
            if goal in adj[v]:
                out.append(path)
            continue
        for w in adj[v]:
            if w not in path and w != goal:
                stack.append((w, path + (w,)))
    return out

def _canon(cycle):
    best = None
    for i in range(len(cycle)):
        r1 = cycle[i:] + cycle[:i]
        for rot in (r1, tuple(reversed(r1))):
            if best is None or rot < best:
                best = rot
    return best

def cycles_through(E, n, edges, k):
    adj = adj_from(E, n)
    cyc = set()
    for a, b in edges:
        if (min(a, b), max(a, b)) not in E:
            continue
        for path in _paths(adj, b, a, k - 1):
            cyc.add(_canon((a,) + path))
    return cyc

def count_all(E, n, k):
    return len(cycles_through(E, n, list(E), k))

CUBE_SPEC = sorted([3, 1, 1, 1, -1, -1, -1, -3])
def defect_fraction(E, n):
    adj = adj_from(E, n)
    seen, good = set(), 0
    for s in range(n):
        if s in seen:
            continue
        comp, stack = {s}, [s]
        while stack:
            v = stack.pop()
            for w in adj[v]:
                if w not in comp:
                    comp.add(w); stack.append(w)
        seen |= comp
        if len(comp) == 8 and all(len(adj[v]) == 3 for v in comp):
            idx = {v: i for i, v in enumerate(sorted(comp))}
            M = np.zeros((8, 8))
            for v in comp:
                for w in adj[v]:
                    M[idx[v], idx[w]] = 1.0
            if sorted(np.round(np.linalg.eigvalsh(M)).astype(int).tolist()) == CUBE_SPEC:
                good += 8
    return 1.0 - good / n

def random_cubic(n):
    while True:
        stubs = [v for v in range(n) for _ in range(3)]
        random.shuffle(stubs)
        E, ok = set(), True
        for i in range(0, len(stubs), 2):
            a, b = stubs[i], stubs[i + 1]
            if a == b or (min(a, b), max(a, b)) in E:
                ok = False; break
            E.add((min(a, b), max(a, b)))
        if ok:
            return E

def cold_start(n):
    E = set()
    for blk in range(n // 8):
        for u in range(8):
            for bit in (1, 2, 4):
                v = u ^ bit
                if u < v:
                    E.add((blk * 8 + u, blk * 8 + v))
    return E

E = cold_start(N) if START == "cold" else random_cubic(N)
c4, c6 = count_all(E, N, 4), count_all(E, N, 6)
acc, dsum, dn, dlate, dln = 0, 0.0, 0, 0.0, 0
tail = SWEEPS - SWEEPS // 3
late = SWEEPS - SWEEPS // 10
for sweep in range(SWEEPS):
    for _ in range(N):
        El = list(E)
        (a, b) = random.choice(El); (c, d) = random.choice(El)
        if len({a, b, c, d}) < 4:
            continue
        if random.random() < 0.5:
            n1, n2 = (min(a, c), max(a, c)), (min(b, d), max(b, d))
        else:
            n1, n2 = (min(a, d), max(a, d)), (min(b, c), max(b, c))
        if n1 in E or n2 in E:
            continue
        E2 = set(E)
        E2.discard((min(a, b), max(a, b))); E2.discard((min(c, d), max(c, d)))
        E2.add(n1); E2.add(n2)
        rm = [e for e in E if e not in E2]; ad = [e for e in E2 if e not in E]
        dc4 = len(cycles_through(E2, N, ad, 4)) - len(cycles_through(E, N, rm, 4))
        dc6 = len(cycles_through(E2, N, ad, 6)) - len(cycles_through(E, N, rm, 6))
        dE = -W4 * dc4 - W6 * dc6
        if dE <= 0 or random.random() < math.exp(-dE / T):
            E = E2; c4 += dc4; c6 += dc6; acc += 1
            if acc % 20000 == 0:
                assert c4 == count_all(E, N, 4) and c6 == count_all(E, N, 6)
    if sweep >= tail and sweep % 5 == 0:
        dd = defect_fraction(E, N)
        dsum += dd; dn += 1
        if sweep >= late:
            dlate += dd; dln += 1
assert c4 == count_all(E, N, 4) and c6 == count_all(E, N, 6)
print(json.dumps(dict(N=N, w4=W4, w6=W6, T=T, start=START, rep=REP,
                      d=dsum / max(dn, 1), d_late=dlate / max(dln, 1),
                      c4=c4, c6=c6, sweeps=SWEEPS)), flush=True)
