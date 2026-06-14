#!/usr/bin/env python3
r"""THE DECISIVE K04 ARREST-DENSITY COMPUTATION — paper-pinned protocol, scaled N.

PROTOCOL (pinned by crystallize/crystallisation.tex, found this session):
  moves     : parallel-swap (Kempe chain) double-edge swaps — cubic ensemble exact
              (the lambda term is a formal hard constraint, never dynamically active);
  T0        : adaptive — chosen so ~80% of initial uphill moves accept;
  cooling   : geometric, gamma = 0.995 per sweep (sweep = N proposals);
  reheat    : T *= 5 after R no-improvement sweeps (paper: 'periodic reheating');
  stop      : after MAXR reheats each followed by R no-improvement sweeps (paper:
              '2-3 reheating cycles suffice'), or hard sweep cap;
  weights   : ORDERING pinned (w4 > w6 > 0, 4-cycle dominant); values free — primary
              (w4, w6) = (2, 1), sensitivity (3, 1) at N = 24.
OBSERVABLE: d = per-vertex defective fraction at termination (vertex defective unless
in a perfect Q3 component — spectral check, stricter than the paper's girth test).
HARNESS GATE: the paper claims 94/100 perfect at N = 24 under this protocol.
Usage: k04_decisive.py [N trials w4 w6 seed tag]   (defaults: 24 10 2 1 11 local)
Self-asserting on exactness (incremental C4+C6 drift checks); results to JSONL."""
import sys, json, random, math
import numpy as np

args = sys.argv[1:]
N = int(args[0]) if len(args) > 0 else 24
TRIALS = int(args[1]) if len(args) > 1 else 10
W4 = float(args[2]) if len(args) > 2 else 2.0
W6 = float(args[3]) if len(args) > 3 else 1.0
SEED = int(args[4]) if len(args) > 4 else 11
TAG = args[5] if len(args) > 5 else "local"
random.seed(SEED); np.random.seed(SEED)
R_NOIMP, MAXR, CAP_SWEEPS = 1000, 3, 40000

def random_cubic(n):
    while True:
        stubs = [v for v in range(n) for _ in range(3)]
        random.shuffle(stubs)
        E = set(); ok = True
        for i in range(0, len(stubs), 2):
            a, b = stubs[i], stubs[i + 1]
            if a == b or (min(a, b), max(a, b)) in E:
                ok = False; break
            E.add((min(a, b), max(a, b)))
        if ok:
            return E

def adj_from(E, n):
    A = [[] for _ in range(n)]
    for a, b in E:
        A[a].append(b); A[b].append(a)
    return A

def _paths(adj, start, goal, length, blocked=None):
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
    """Canonical set of k-cycles through any of `edges` (k = 4 or 6)."""
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

def propose(E):
    El = list(E)
    (a, b) = random.choice(El)
    (c, d) = random.choice(El)
    if len({a, b, c, d}) < 4:
        return None
    if random.random() < 0.5:
        n1, n2 = (min(a, c), max(a, c)), (min(b, d), max(b, d))
    else:
        n1, n2 = (min(a, d), max(a, d)), (min(b, c), max(b, c))
    if n1 in E or n2 in E:
        return None
    E2 = set(E)
    E2.discard((min(a, b), max(a, b))); E2.discard((min(c, d), max(c, d)))
    E2.add(n1); E2.add(n2)
    return E2

def trial(seed):
    random.seed(seed)
    E = random_cubic(N)
    c4, c6 = count_all(E, N, 4), count_all(E, N, 6)
    cur = -W4 * c4 - W6 * c6
    # adaptive T0: ~80% uphill acceptance at start
    ups = []
    for _ in range(300):
        E2 = propose(E)
        if E2 is None:
            continue
        rm = [e for e in E if e not in E2]; ad = [e for e in E2 if e not in E]
        dc4 = len(cycles_through(E2, N, ad, 4)) - len(cycles_through(E, N, rm, 4))
        dc6 = len(cycles_through(E2, N, ad, 6)) - len(cycles_through(E, N, rm, 6))
        dE = -W4 * dc4 - W6 * dc6
        if dE > 0:
            ups.append(dE)
    T = -np.mean(ups) / math.log(0.8) if ups else 1.0
    T0 = T
    best, noimp, reheats, sweeps, acc = cur, 0, 0, 0, 0
    while sweeps < CAP_SWEEPS:
        for _ in range(N):
            E2 = propose(E)
            if E2 is None:
                continue
            rm = [e for e in E if e not in E2]; ad = [e for e in E2 if e not in E]
            dc4 = len(cycles_through(E2, N, ad, 4)) - len(cycles_through(E, N, rm, 4))
            dc6 = len(cycles_through(E2, N, ad, 6)) - len(cycles_through(E, N, rm, 6))
            dE = -W4 * dc4 - W6 * dc6
            if dE <= 0 or random.random() < math.exp(-dE / max(T, 1e-9)):
                E = E2; c4 += dc4; c6 += dc6; cur += dE; acc += 1
                if acc % 5000 == 0:
                    assert c4 == count_all(E, N, 4) and c6 == count_all(E, N, 6)
        sweeps += 1
        T *= 0.995
        if cur < best - 1e-9:
            best, noimp = cur, 0
        else:
            noimp += 1
        if noimp >= R_NOIMP:
            reheats += 1
            if reheats > MAXR:
                break
            T = max(T * 5.0, 0.6 * T0)
            noimp = 0
    assert c4 == count_all(E, N, 4) and c6 == count_all(E, N, 6)
    d = defect_fraction(E, N)
    return dict(N=N, seed=seed, d=d, perfect=(d == 0.0), E=cur, c4=c4, c6=c6,
                sweeps=sweeps, reheats=reheats, w4=W4, w6=W6, tag=TAG)

if __name__ == "__main__":
    out = []
    for t in range(TRIALS):
        r = trial(SEED * 1000 + t)
        out.append(r)
        print(json.dumps(r), flush=True)
    ds = [r["d"] for r in out]
    nperf = sum(r["perfect"] for r in out)
    print(f"# SUMMARY N={N} w=({W4},{W6}) trials={TRIALS}: perfect {nperf}/{TRIALS}, "
          f"mean d = {np.mean(ds):.4f} +- {np.std(ds)/math.sqrt(len(ds)):.4f}", flush=True)
