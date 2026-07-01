#!/usr/bin/env python3
r"""foundations_8piG_dimensional_not_octagon.py

RE-DERIVATION of the §1.5 / Part-15 "8 pi G Einstein-tensor coefficient" after the 2026-06-25 geometry
correction retired the 4.8.8 octagon vertex figure.

The Part-15 derivation read: the 4.8.8 matter-gauge node has graph coordination z=3, a defect's delay eps
distributes as eps/3 (Ollivier-Ricci Dkappa = +/- eps/3), and "that fundamental topological 1/3 from the
degree-3 lattice vertex perfectly maps into the 1/3 volume factor of 3D space" (V_3 = 4 pi r^3 / 3), after
which Poisson reconstructs 8 pi.

The tell: the 1/3 in the 3-ball volume is 1/d (d = 3 SPATIAL DIMENSIONS), a dimensional fact -- it has
nothing to do with graph coordination. The Part-15 argument conflated d=3 (dimension) with z=3 (the 4.8.8
coordination) because they are numerically equal on that tiling. This script separates them.

PART A (rigorous): 8 pi is DIMENSIONAL.
  - the d-ball volume radial factor is 1/d (d=2 -> 1/2, d=3 -> 1/3, d=4 -> 1/4): tracks DIMENSION.
  - S_{d-1} (unit (d-1)-sphere area) = 4 pi for d=3 -> Newtonian Poisson grad^2 Phi = 4 pi G rho;
    the relativistic trace-reversal multiplies by 2 -> Einstein G_uv = 8 pi G T_uv. All dimensional; no
    lattice, no octagon.

PART B (computed): the Ollivier-Ricci shift is LINEAR in eps (the load-bearing proportionality, the only
  piece used downstream at 8.2/10.8), with an O(1), graph-dependent slope -- it is NEITHER 1/z NOR the
  dimensional 1/d. So the curvature does not supply the sphere's 1/3 at all; the "eps/3 = sphere-1/3"
  identification was unsupported on two counts.

VERDICT: the 8 pi G coefficient SURVIVES on the corrected (oblate-bipyramid) geometry -- it is a property
of 3+1D spacetime (still d=3+1 here), recoverable on ANY 3D substrate. The octagon z=3 was coincidental
dressing; removing it strengthens the result. The Ollivier-Ricci curvature LINEARITY survives; the
specific "1/3" transport coefficient was never load-bearing.

Self-asserting; exit 0.
"""
from math import gamma, pi
from itertools import product
from collections import deque
import numpy as np
from scipy.optimize import linprog


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def main():
    print("=== Re-deriving 8 pi G on the corrected geometry: dimensional, not octagon ===\n")

    # ---- PART A: 8 pi is dimensional ----
    print("[A] the '1/3' in V_3 = 4 pi r^3/3 is 1/d (dimension), not 1/z (coordination):")
    for d in (2, 3, 4):
        Sd1 = 2 * pi ** (d / 2) / gamma(d / 2)
        Vd = pi ** (d / 2) / gamma(d / 2 + 1)
        print(f"    d={d}: S_(d-1)={Sd1:.4f}, V_d/S_(d-1) = 1/d = {Vd/Sd1:.4f}")
    S2 = 2 * pi ** 1.5 / gamma(1.5)
    ok(abs((pi ** 1.5 / gamma(2.5)) / S2 - 1 / 3) < 1e-9, "3-ball radial factor = 1/d = 1/3 (DIMENSIONAL; d=2->1/2, d=4->1/4)")
    ok(abs(S2 - 4 * pi) < 1e-9, "S_2 = 4 pi (3D Gauss) -> Newtonian 4 pi G; Einstein 8 pi G = 2 x 4 pi (trace-reversal). No lattice.")

    # ---- PART B: Ollivier-Ricci shift is linear in eps; slope is O(1), not 1/z and not 1/d ----
    def bfs(adj, s):
        dist = {s: 0}; q = deque([s])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if v not in dist:
                    dist[v] = dist[u] + 1; q.append(v)
        return dist

    def kappa(adj, x, y, eps):
        nx, ny = list(adj[x]), list(adj[y])
        mx = {x: eps}
        for u in nx:
            mx[u] = mx.get(u, 0) + (1 - eps) / len(nx)
        my = {}
        for v in ny:
            my[v] = my.get(v, 0) + 1.0 / len(ny)
        A, B = list(mx), list(my)
        D = {a: bfs(adj, a) for a in set(A)}
        C = np.array([[D[a][b] for b in B] for a in A], float)
        nA, nB = len(A), len(B)
        Aeq, beq = [], []
        for i in range(nA):
            r = np.zeros(nA * nB); r[i * nB:(i + 1) * nB] = 1; Aeq.append(r); beq.append(mx[A[i]])
        for j in range(nB):
            r = np.zeros(nA * nB); r[j::nB] = 1; Aeq.append(r); beq.append(my[B[j]])
        res = linprog(C.ravel(), A_eq=np.array(Aeq), b_eq=np.array(beq), bounds=[(0, None)] * nA * nB)
        return 1.0 - res.fun

    def cube_graph():
        V = list(product((0, 1), repeat=3)); adj = {v: [] for v in V}
        for u in V:
            for i in range(3):
                w = list(u); w[i] ^= 1; adj[u].append(tuple(w))
        return adj

    def torus3(L):
        V = list(product(range(L), repeat=3)); adj = {v: [] for v in V}
        for u in V:
            for i in range(3):
                for s in (1, -1):
                    w = list(u); w[i] = (w[i] + s) % L; adj[u].append(tuple(w))
        return adj

    print("\n[B] Ollivier-Ricci curvature shift of a defect (laziness eps), via exact W1:")
    eps = np.array([0.02, 0.04, 0.06, 0.08, 0.10])
    for name, adj, z in [("octagon coordination z=3 (cube graph Q3)", cube_graph(), 3),
                         ("CORRECT matter lattice z=6 (3-torus)   ", torus3(4), 6)]:
        x = list(adj)[0]; y = adj[x][0]; k0 = kappa(adj, x, y, 0.0)
        dk = np.array([kappa(adj, x, y, e) - k0 for e in eps])
        slope = np.polyfit(eps, dk, 1)[0]
        resid = np.max(np.abs(dk - slope * eps))
        print(f"    {name}: slope dDkappa/deps = {slope:.3f}  (1/z={1/z:.3f}, 1/d=0.333); linear (max resid {resid:.1e})")
        ok(resid < 1e-3, f"Dkappa LINEAR in eps on {name.split('(')[0].strip()} (proportionality survives)")
        ok(abs(slope - 1 / z) > 0.2 and abs(slope - 1 / 3) > 0.2,
           "slope is O(1), NOT 1/z and NOT the dimensional 1/d -> the curvature does not supply the sphere's 1/3")

    print("\n[verdict] 8 pi G SURVIVES on the corrected geometry -- it is DIMENSIONAL:")
    print("  - the 1/3 in V_3=4pi r^3/3 is 1/d (d=3 spatial dims); 8pi = S_2(=4pi) x 2(trace-reversal).")
    print("    These hold on ANY 3+1D substrate, including the oblate-bipyramid tiling (still d=3+1).")
    print("  - the Part-15 'z=3 octagon -> 1/3 -> 8pi' conflated the dimensional 1/d with the 4.8.8")
    print("    coordination 1/z (equal only by coincidence, d=3=z=3 there). On the corrected geometry the")
    print("    matter coordination is z=6, so the coincidence is gone -- but 8pi does NOT depend on z, so")
    print("    it is UNAFFECTED, and now correctly attributed to dimension. The octagon was dressing.")
    print("  - the Ollivier-Ricci curvature LINEARITY (Dkappa proportional to eps -- the only piece used at")
    print("    8.2 Yukawa / 10.8 hierarchy) survives; the specific 'eps/3' coefficient is convention-")
    print("    dependent (standard OR gives slope ~1, not eps/3) and was never load-bearing. exit 0")


if __name__ == "__main__":
    main()
