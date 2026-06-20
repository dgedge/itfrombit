#!/usr/bin/env python3
"""
bh_fast_scrambling_expander_test.py

Black-hole sector (DRIFT B5): fast scrambling (Hayden-Preskill, t_* ~ log S)
requires the horizon-cell service graph to be an EXPANDER (spectral gap bounded
below, independent of N_H). The documented substrate horizon is a 2-sphere shell
of cells with LOCAL connectivity (V_Sch = direct sum over shell cells; KMS
service one-bit local). Does that local-surface graph fast-scramble?

Test: compute the Fiedler value (spectral gap) lambda_2 of the graph Laplacian
vs horizon size N_H for the substrate-faithful horizon graphs (2D torus shell;
icosphere = a genuine 2-sphere triangulation), against a random 3-regular
EXPANDER control at the same N and degree.

Diagnostic:
  expander  -> lambda_2 -> const  (gap survives; diameter ~ log N; t_* ~ log S)
  local 2D  -> lambda_2 ~ 1/N_H   (gap collapses; diameter ~ sqrt N; t_* ~ poly S)

A bounded-degree, sphere-embeddable (genus-0) graph CANNOT be an expander
(planar separator theorem => lambda_2 = O(1/N)), so the result is a theorem the
torus/icosphere only illustrate.

Self-asserting; numpy only (+ stdlib).
"""
import sys
from collections import deque
from itertools import combinations
import numpy as np

rng = np.random.default_rng(7)
_ok = True


def check(name, cond):
    global _ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    _ok = _ok and bool(cond)


def lambda2(N, edges):
    A = np.zeros((N, N))
    for u, v in edges:
        A[u, v] += 1.0
        A[v, u] += 1.0
    L = np.diag(A.sum(1)) - A
    ev = np.linalg.eigvalsh(L)
    return float(ev[1])                      # Fiedler value (2nd smallest)


def diameter(N, edges):
    adj = [[] for _ in range(N)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)

    def bfs(s):
        d = [-1] * N
        d[s] = 0
        q = deque([s])
        far = s
        while q:
            x = q.popleft()
            for y in adj[x]:
                if d[y] < 0:
                    d[y] = d[x] + 1
                    q.append(y)
                    far = y
        return far, max(d)
    a, _ = bfs(0)
    _, dd = bfs(a)                            # double-BFS approx diameter
    return dd


# ---- horizon models -------------------------------------------------------
def torus(L):
    def idx(i, j):
        return (i % L) * L + (j % L)
    E = []
    for i in range(L):
        for j in range(L):
            E += [(idx(i, j), idx(i + 1, j)), (idx(i, j), idx(i, j + 1))]
    return L * L, E


def icosphere(subdiv):
    t = (1 + 5 ** 0.5) / 2
    base = np.array([(-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0),
                     (0, -1, t), (0, 1, t), (0, -1, -t), (0, 1, -t),
                     (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)], float)
    base = base / np.linalg.norm(base[0])
    dmin = min(np.linalg.norm(base[i] - base[j]) for i, j in combinations(range(12), 2)) + 1e-6
    adj = [[j for j in range(12) if j != i and np.linalg.norm(base[i] - base[j]) < dmin]
           for i in range(12)]
    faces = [(i, j, k) for i, j, k in combinations(range(12), 3)
             if j in adj[i] and k in adj[j] and k in adj[i]]
    verts = [v for v in base]
    vmap = {tuple(np.round(v, 5)): i for i, v in enumerate(base)}

    def vid(p):
        key = tuple(np.round(p, 5))
        if key not in vmap:
            vmap[key] = len(verts)
            verts.append(p)
        return vmap[key]

    for _ in range(subdiv):
        nf = []
        for (a, b, c) in faces:
            pa, pb, pc = verts[a], verts[b], verts[c]
            m = []
            for p, q in ((pa, pb), (pb, pc), (pc, pa)):
                w = (p + q) / 2
                m.append(vid(w / np.linalg.norm(w)))
            iab, ibc, ica = m
            nf += [(a, iab, ica), (b, ibc, iab), (c, ica, ibc), (iab, ibc, ica)]
        faces = nf
    E = set()
    for (a, b, c) in faces:
        for u, v in ((a, b), (b, c), (c, a)):
            E.add((min(u, v), max(u, v)))
    return len(verts), list(E)


def random_regular(N, d=3):
    E = []
    for _ in range(d):
        p = rng.permutation(N)
        E += [(int(p[i]), int(p[i + 1])) for i in range(0, N, 2)]
    return N, E


# ===========================================================================
# 1. local-surface horizon (2D torus shell): lambda_2 ~ 1/N_H ?
# ===========================================================================
print("TORUS horizon shell (local, degree 4):")
Ns_t, l2_t = [], []
for L in (6, 8, 10, 14, 18, 22, 26):
    N, E = torus(L)
    g = lambda2(N, E)
    Ns_t.append(N)
    l2_t.append(g)
    print(f"  N_H={N:4d}  lambda_2={g:.5f}  lambda_2*N_H={g*N:.2f}  diam={diameter(N,E)}")
slope_t = np.polyfit(np.log(Ns_t), np.log(l2_t), 1)[0]
print(f"  -> log-log slope(lambda_2 vs N_H) = {slope_t:.3f}  (expect -1)")
check("torus horizon gap collapses as ~1/N_H (slope in [-1.2,-0.8])",
      -1.2 < slope_t < -0.8)

# ===========================================================================
# 2. genuine 2-sphere horizon (icosphere): same collapse
# ===========================================================================
print("\nICOSPHERE horizon (genus-0 2-sphere, ~degree 6):")
Ns_s, l2_s = [], []
for k in (1, 2, 3):
    N, E = icosphere(k)
    g = lambda2(N, E)
    Ns_s.append(N)
    l2_s.append(g)
    print(f"  N_H={N:4d}  lambda_2={g:.5f}  lambda_2*N_H={g*N:.2f}  diam={diameter(N,E)}")
slope_s = np.polyfit(np.log(Ns_s), np.log(l2_s), 1)[0]
print(f"  -> log-log slope = {slope_s:.3f}  (expect ~-1)")
check("icosphere (real horizon topology) gap collapses ~1/N_H (slope in [-1.3,-0.7])",
      -1.3 < slope_s < -0.7)

# ===========================================================================
# 3. EXPANDER control: random 3-regular -- gap survives
# ===========================================================================
print("\nRANDOM 3-REGULAR expander control (same degree class):")
l2_r = []
for N in (64, 256, 784):
    _, E = random_regular(N, 3)
    g = lambda2(N, E)
    l2_r.append(g)
    print(f"  N={N:4d}  lambda_2={g:.4f}  diam={diameter(N,E)}")
slope_r = np.polyfit(np.log([64, 256, 784]), np.log(l2_r), 1)[0]
print(f"  -> log-log slope = {slope_r:.3f}  (expect ~0, gap bounded below)")
check("random-regular gap STAYS O(1) (min lambda_2 > 0.1; slope > -0.2)",
      min(l2_r) > 0.1 and slope_r > -0.2)

# ===========================================================================
# 4. the diverging gap ratio + scrambling-time consequence
# ===========================================================================
# gap ratio expander/local: ~0.4 at N=64 (small tori are well-connected) but
# grows ~linearly in N as the local gap collapses (~39/N) and the expander gap
# stays ~0.17 -> ratio ~ N/229 -> infinity.
r_small = l2_r[0] / l2_t[1]                   # N=64
r_large = l2_r[-1] / l2_t[-1]                 # N~780
print(f"\nexpander/local gap ratio: {r_small:.2f} at N=64  ->  {r_large:.2f} at N~780  "
      f"(crosses 1 and diverges ~N; ~{0.17*1e6/39:.0f}x at N=1e6)")
check("gap ratio grows + crosses 1 (expander wins at large N -> divergence)",
      r_large > 1.0 and r_large > 3 * r_small)

print("\n--- VERDICT ---")
print("The substrate horizon-cell graph (local 2-sphere shell -- torus AND icosphere)")
print("has spectral gap lambda_2 ~ 1/N_H -> 0: it is NOT an expander. A random graph")
print("of the SAME degree keeps lambda_2 ~ O(1), so this is a property of the LOCAL")
print("2D-surface connectivity, not of the degree. (Rigorously: a bounded-degree")
print("genus-0 graph obeys the planar separator theorem, lambda_2 = O(1/N) -- a")
print("theorem the icosphere only illustrates.)")
print()
print("CONSEQUENCE: random-walk mixing time ~ 1/lambda_2 ~ N_H and graph diameter")
print("~ sqrt(N_H), so information scrambles in POWER-LAW time in the entropy")
print("(t_* ~ sqrt(S) ballistic to ~S diffusive), NOT the holographic t_* ~ log S.")
print()
print("FALSIFIABLE PREDICTION: with its documented LOCAL horizon couplings the")
print("framework's black holes are NOT fast scramblers -- a structural deviation")
print("from the Sekino-Susskind fast-scrambling conjecture. Fast scrambling would")
print("require NON-LOCAL horizon-cell connectivity (long-range links / all-to-all),")
print("which the substrate (direct-sum V_Sch, one-bit-local KMS service) does not")
print("supply. So DRIFT B5 closes NEGATIVE: no substrate-native fast scrambler.")
print("Scope: assumes scrambling propagates on the spatial horizon-cell adjacency")
print("(the documented coupling); a non-spatial global coupling, if later derived,")
print("would change this -- but none is in canon.")

print("\n" + ("ALL CHECKS PASSED" if _ok else "SOME CHECKS FAILED"))
sys.exit(0 if _ok else 1)
