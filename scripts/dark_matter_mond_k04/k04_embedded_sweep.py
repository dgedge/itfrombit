#!/usr/bin/env python3
r"""K04 IN THE Z3-EMBEDDED ENSEMBLE — the substrate-faithful closing computation.

Replaces the configuration-model toy (k04_eq_sweep.py), whose ground state above
w4/w6 = 4/3 is the K33 crystal — proven NON-EMBEDDABLE in Z3 (bipartite => no K4;
max 2 common neighbours => no K33; debris_dark_matter_audit.py). Here:
  STATE : bond subset B of the periodic L^3 lattice with EVERY site degree exactly 3
          (each site keeps 3 of its 6 incident bonds).
  COLD  : exact cube tiling — bond (p, axis) occupied iff p[axis] is even (L even);
          components = disjoint perfect Q3 cells; C4 = 3N/4 (the proven maximum:
          each site has <= 3 axis-pairs, tiling completes all), C6 = 2N.
  HOT   : tiling randomised by T = infinity plaquette swaps.
  MOVES : plaquette (unit-square) swaps — a square with exactly one PARALLEL pair
          occupied and the other free swaps pairs; degree-preserving by construction.
  ENERGY: E = -w4 C4 - w6 C6 over CONTRACTIBLE cycles only — path search closes in
          UNWRAPPED coordinates, so torus-winding cycles never count, any L >= 4.
  OBS   : d = fraction of sites not in a perfect isolated Q3 cell (spectral check);
          d_late (last tenth) vs d (last third) for the branch-merge instrument;
          full debris census + excess energy at termination.
Usage: k04_embedded_sweep.py [L w4 w6 T start sweeps rep]   (one JSON line to stdout)
Self-asserting (exact incremental-count drift checks); exit 0 = verified."""
import sys, json, math, random
import numpy as np

args = sys.argv[1:]
L      = int(args[0])   if len(args) > 0 else 4
W4     = float(args[1]) if len(args) > 1 else 1.7
W6     = float(args[2]) if len(args) > 2 else 1.0
T      = float(args[3]) if len(args) > 3 else 2.0
START  = args[4]        if len(args) > 4 else "cold"
SWEEPS = int(args[5])   if len(args) > 5 else 2000
REP    = int(args[6])   if len(args) > 6 else 1
HOLD   = int(args[7])   if len(args) > 7 else -1     # ramp hold; -1 -> R/10 (aging knob)
assert L % 2 == 0 and START in ("cold", "hot", "ramp")
# ramp protocol (KZ): hot-init, geometric cool T_START -> T over SWEEPS sweeps
# (SWEEPS = R, the KZ control parameter; 1 sweep ~ 1 tick), then hold at T.
T_START_RAMP = 6.0
random.seed(hash((L, W4, T, START, REP, 20260611)) % (2**32))
N = L ** 3
AX = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]

def wrap(p):
    return (p[0] % L, p[1] % L, p[2] % L)
def add(p, s):
    return (p[0] + s[0], p[1] + s[1], p[2] + s[2])

# ---------------- state ----------------
def tiling():
    return {(p, a) for p in ((x, y, z) for x in range(L) for y in range(L) for z in range(L))
            for a in range(3) if p[a] % 2 == 0}
def degree_map(B):
    deg = {}
    for (p, a) in B:
        q = wrap(add(p, AX[a]))
        deg[p] = deg.get(p, 0) + 1
        deg[q] = deg.get(q, 0) + 1
    return deg

# ---------------- contractible cycle counting (unwrapped closure) ----------------
def nbrs(B, p):
    """occupied neighbours of wrapped site p, as (wrapped_q, unwrapped_step)."""
    out = []
    for a in range(3):
        if (p, a) in B:
            out.append((wrap(add(p, AX[a])), AX[a]))
        m = wrap(add(p, (-AX[a][0], -AX[a][1], -AX[a][2])))
        if (m, a) in B:
            out.append((m, (-AX[a][0], -AX[a][1], -AX[a][2])))
    return out

def cycles_through_bonds(B, bonds, k):
    """canonical set of contractible k-cycles through any of `bonds`."""
    cyc = set()
    for (p, a) in bonds:
        if (p, a) not in B:
            continue
        q = wrap(add(p, AX[a]))
        # paths of k-1 steps from q back to p, unwrapped displacement must cancel
        stack = [(q, AX[a], (p, q))]
        while stack:
            v, disp, path = stack.pop()
            for (w, s) in nbrs(B, v):
                nd = add(disp, s)
                if len(path) == k:
                    if w == p and nd == (0, 0, 0):       # k-th step closes, zero winding
                        best = None
                        for i in range(len(path)):
                            r1 = path[i:] + path[:i]
                            for rot in (r1, tuple(reversed(r1))):
                                if best is None or rot < best:
                                    best = rot
                        cyc.add(best)
                elif w not in path:
                    if abs(nd[0]) + abs(nd[1]) + abs(nd[2]) <= k - len(path):
                        stack.append((w, nd, path + (w,)))
    return cyc

def count_all(B, k):
    return len(cycles_through_bonds(B, list(B), k))

# ---------------- plaquette move ----------------
SQUARES = [(p, a, b) for p in ((x, y, z) for x in range(L) for y in range(L) for z in range(L))
           for a in range(3) for b in range(3) if a < b]
def square_bonds(p, a, b):
    pa = (wrap(add(p, AX[a])), b)
    pb = (wrap(add(p, AX[b])), a)
    return ((p, a), pb), ((p, b), pa)        # parallel pair along a, parallel pair along b

def propose(B):
    p, a, b = random.choice(SQUARES)
    pairA, pairB = square_bonds(p, a, b)
    inA = sum(1 for e in pairA if e in B)
    inB = sum(1 for e in pairB if e in B)
    if inA == 2 and inB == 0:
        return pairA, pairB
    if inB == 2 and inA == 0:
        return pairB, pairA
    return None

# ---------------- healing spectrum (readout #4: durability, distinct from mass) ----------------
def healing_spectrum(B):
    """alternating-trail lengths of B (+) nearest anchor tiling (see k04_ergodicity_audit)."""
    best = None
    for an in ((i, j, k) for i in (0, 1) for j in (0, 1) for k in (0, 1)):
        Tt = {(p, a) for (p, a) in ((s, ax) for s in
              ((x, y, z) for x in range(L) for y in range(L) for z in range(L))
              for ax in range(3)) if (p[a] - an[a]) % 2 == 0}
        D = B ^ Tt
        if best is None or len(D) < len(best[0]):
            best = (D, an)
    D, an = best
    inc = {}
    for bond in D:
        (p, a) = bond
        q = wrap(add(p, AX[a]))
        tag = "B" if bond in B else "T"
        inc.setdefault(p, []).append((bond, q, tag))
        inc.setdefault(q, []).append((bond, p, tag))
    for v, lst in inc.items():
        assert 2 * sum(1 for *_, t in lst if t == "B") == len(lst)
    used, lengths = set(), []
    for v0 in list(inc):
        while True:
            start = next(((b, q, t) for (b, q, t) in inc.get(v0, []) if b not in used), None)
            if start is None:
                break
            trail, v0_, (b, q, t) = 0, v0, start
            while True:
                used.add(b); trail += 1
                v, need = q, ("T" if t == "B" else "B")
                nxt = next(((bb, qq, tt) for (bb, qq, tt) in inc[v]
                            if bb not in used and tt == need), None)
                if v == v0_ and trail % 2 == 0 and nxt is None:
                    break
                b, q, t = nxt
            lengths.append(trail)
    assert sum(lengths) == len(D)
    return sorted(lengths, reverse=True), list(an), len(D)

# ---------------- percolation of the crystallised (register-bearing) network ----------------
def percolation(B):
    """vertex-level percolation of the good-cell material: adjacency = Z3 nearest
    neighbour between crystallised vertices (spatially adjacent cells are exactly
    what the gauge web can bridge — bond occupancy between cells not required).
    Returns (n_goodcells, largest-cluster fraction of good vertices, spans_torus);
    spanning = some vertex reached at two unwrapped offsets differing by >= L."""
    cells, _ = census(B, return_cells=True)
    good = set().union(*cells) if cells else set()
    if not good:
        return 0, 0.0, False
    seen, best, spans = set(), 0, False
    steps = [s for a in range(3) for s in (AX[a], (-AX[a][0], -AX[a][1], -AX[a][2]))]
    for v0 in good:
        if v0 in seen:
            continue
        offset = {v0: (0, 0, 0)}
        stack, comp = [v0], {v0}
        while stack:
            v = stack.pop()
            for s in steps:
                w = wrap(add(v, s))
                if w not in good:
                    continue
                off = add(offset[v], s)
                if w in comp:
                    d0 = offset[w]
                    if any(abs(off[k] - d0[k]) >= L for k in range(3)):
                        spans = True
                else:
                    comp.add(w); offset[w] = off; stack.append(w)
        seen |= comp
        best = max(best, len(comp))
    return len(cells), best / len(good), spans

# ---------------- census ----------------
CUBE_SPEC = sorted([3, 1, 1, 1, -1, -1, -1, -3])
def census(B, return_cells=False):
    adj = {}
    for (p, a) in B:
        q = wrap(add(p, AX[a]))
        adj.setdefault(p, []).append(q)
        adj.setdefault(q, []).append(p)
    seen, comps = set(), {}
    good, cells = 0, []
    for s in adj:
        if s in seen:
            continue
        comp, stack = {s}, [s]
        while stack:
            v = stack.pop()
            for w in adj[v]:
                if w not in comp:
                    comp.add(w); stack.append(w)
        seen |= comp
        kind = f"other{len(comp)}"
        if len(comp) == 8:
            idx = {v: i for i, v in enumerate(sorted(comp))}
            M = np.zeros((8, 8))
            for v in comp:
                for w in adj[v]:
                    M[idx[v], idx[w]] = 1.0      # set once per direction -> symmetric
            if sorted(np.round(np.linalg.eigvalsh(M)).astype(int).tolist()) == CUBE_SPEC:
                kind = "Q3"; good += 8
                cells.append(frozenset(comp))
        comps[kind] = comps.get(kind, 0) + 1
    if return_cells:
        return cells, 1.0 - good / N
    return comps, 1.0 - good / N

# ---------------- run ----------------
B = tiling()
deg = degree_map(B)
assert all(v == 3 for v in deg.values()) and len(deg) == N
c4, c6 = count_all(B, 4), count_all(B, 6)
assert c4 == 3 * N // 4 and c6 == 2 * N, f"tiling counts wrong: {c4}, {c6}"

def metropolis(B, c4, c6, T, sweeps, record=None):
    acc = 0
    for s in range(sweeps):
        for _ in range(N):
            mv = propose(B)
            if mv is None:
                continue
            rm, ad = mv
            dc4_rm = len(cycles_through_bonds(B, rm, 4))
            dc6_rm = len(cycles_through_bonds(B, rm, 6))
            B2 = set(B); B2.discard(rm[0]); B2.discard(rm[1]); B2.add(ad[0]); B2.add(ad[1])
            dc4_ad = len(cycles_through_bonds(B2, ad, 4))
            dc6_ad = len(cycles_through_bonds(B2, ad, 6))
            dc4, dc6 = dc4_ad - dc4_rm, dc6_ad - dc6_rm
            dE = -W4 * dc4 - W6 * dc6
            if T == math.inf or dE <= 0 or random.random() < math.exp(-dE / T):
                B = B2; c4 += dc4; c6 += dc6; acc += 1
                if acc % 20000 == 0:
                    assert c4 == count_all(B, 4) and c6 == count_all(B, 6)
        if record is not None and s % 50 == 0:
            record.append(census(B)[1])
    return B, c4, c6, acc

if START in ("hot", "ramp"):
    B, c4, c6, _ = metropolis(B, c4, c6, math.inf, max(1, 3000 // N + 1) * 20)
    deg = degree_map(B)
    assert all(v == 3 for v in deg.values())

series = []
if START == "ramp":
    acc = 0
    for s in range(SWEEPS):
        Ts = T_START_RAMP * (T / T_START_RAMP) ** (s / max(SWEEPS - 1, 1))
        B, c4, c6, a = metropolis(B, c4, c6, Ts, 1)
        acc += a
        if s % 50 == 0:
            series.append(census(B)[1])
    hold = HOLD if HOLD >= 0 else max(20, SWEEPS // 10)
    B, c4, c6, a = metropolis(B, c4, c6, T, hold, record=series)
    acc += a
else:
    B, c4, c6, acc = metropolis(B, c4, c6, T, SWEEPS, record=series)
assert c4 == count_all(B, 4) and c6 == count_all(B, 6)
deg = degree_map(B)
assert all(v == 3 for v in deg.values()) and len(deg) == N

comps, d_now = census(B)
n3 = max(1, len(series) // 3); n10 = max(1, len(series) // 10)
d_tail = float(np.mean(series[-n3:])) if series else d_now
d_late = float(np.mean(series[-n10:])) if series else d_now
E_conf = -(W4 * c4 + W6 * c6)
E_tile = -(W4 * 3 * N / 4 + W6 * 2 * N)
ndeb = round(d_now * N)
spec, anch, dsz = healing_spectrum(B)        # readout #4: durability vs mass
ncells, perc_frac, perc_spans = percolation(B)   # readout #5: visible-sector viability
print(json.dumps(dict(L=L, N=N, w4=W4, w6=W6, T=T, start=START, rep=REP,
                      d=d_tail, d_late=d_late, d_final=d_now, c4=c4, c6=c6,
                      census=comps, e_excess=E_conf - E_tile,
                      e_per_debris_v=(E_conf - E_tile) / ndeb if ndeb else 0.0,
                      heal_spec=spec[:30], heal_ntrails=len(spec),
                      heal_diff=dsz, heal_anchor=anch,
                      n_goodcells=ncells, perc_frac=perc_frac, perc_spans=perc_spans,
                      hold=(HOLD if HOLD >= 0 else max(20, SWEEPS // 10)) if START == "ramp" else 0,
                      sweeps=SWEEPS, acc=acc)), flush=True)
