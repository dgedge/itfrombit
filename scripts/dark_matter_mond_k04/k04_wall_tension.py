#!/usr/bin/env python3
r"""K04 DOMAIN-WALL TENSION sigma_wall — the energy scale for the defect-network
abundance (k04_defect_network_abundance.py), computed exactly.

Setup as in k04_kempe_locked_defect.py / k04_embedded_sweep.py: degree-3 bond
subsets of the periodic Z3 lattice; perfect cells = isolated 2x2x2 cubes; energy
E = -w4 C4 - w6 C6 over contractible cycles; the string tension is mu = w4 + 4 w6.

KEY STRUCTURE (surfaced by the construction, must be reported): the cube-crystal
ground state is DEGENERATE. Because the crystal is a set of ISOLATED Q3 cubes, the
energy is (N/8)(6 w4 + 16 w6) for ANY partition of Z3 into 2x2x2 cubes -- not just
the 8 axis-aligned tilings. So a domain wall between two cube-tilings is FREE
whenever the two tilings stay compatible across the interface (the cubes simply
don't straddle it). Only a FRUSTRATED interface -- one that forces cubes to be cut
-- costs energy. This script measures both, and that splits the abundance into a
free (entropic) part and a costed (energetic, durable) part.

  [1] FREE WALL  : flip the x-dimerisation phase inside a z-slab. The interface is
                  perpendicular to z and cube-aligned, so no cube straddles it.
                  Predict dE = 0 (a distinct ground state) -> degeneracy.
  [2] COSTED WALL: flip the y-dimerisation phase inside an x-region whose boundaries
                  fall MID-cube (odd x). Cut cubes line the interface -> dE > 0.
                  sigma_wall = dE / (2 * Ly * Lz). Homology [D]=(0,0,0): the wall is
                  an ENERGETIC, kinetically-frozen defect, not a topological one.
  [3] STRING     : mu = w4 + 4 w6 per step (the 1D costed defect), for comparison.

All dE via EXACT local cycle deltas through the toggled edges (unchanged-only cycles
cancel). exit 0 = free wall is free; costed wall tension measured; mu reproduced;
the degeneracy/abundance reframe stated.
"""
import math
import numpy as np

W4, W6 = 1.7, 1.0
AX = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]

def wrap(p, L): return (p[0] % L, p[1] % L, p[2] % L)
def add(p, s): return (p[0] + s[0], p[1] + s[1], p[2] + s[2])

def tiling(L):
    return {(p, a) for p in ((x, y, z) for x in range(L) for y in range(L) for z in range(L))
            for a in range(3) if p[a] % 2 == 0}

def degree_map(B, L):
    deg = {}
    for (p, a) in B:
        q = wrap(add(p, AX[a]), L)
        deg[p] = deg.get(p, 0) + 1; deg[q] = deg.get(q, 0) + 1
    return deg

def nbrs(B, p, L):
    out = []
    for a in range(3):
        if (p, a) in B:
            out.append((wrap(add(p, AX[a]), L), AX[a]))
        m = wrap(add(p, (-AX[a][0], -AX[a][1], -AX[a][2])), L)
        if (m, a) in B:
            out.append((m, (-AX[a][0], -AX[a][1], -AX[a][2])))
    return out

def cycles_through_bonds(B, bonds, k, L):
    cyc = set()
    for (p, a) in bonds:
        if (p, a) not in B:
            continue
        q = wrap(add(p, AX[a]), L)
        stack = [(q, AX[a], (p, q))]
        while stack:
            v, disp, path = stack.pop()
            for (w, s) in nbrs(B, v, L):
                nd = add(disp, s)
                if len(path) == k:
                    if w == p and nd == (0, 0, 0):
                        best = None
                        for i in range(len(path)):
                            r1 = path[i:] + path[:i]
                            for rot in (r1, tuple(reversed(r1))):
                                if best is None or rot < best:
                                    best = rot
                        cyc.add(best)
                elif w not in path and abs(nd[0]) + abs(nd[1]) + abs(nd[2]) <= k - len(path):
                    stack.append((w, nd, path + (w,)))
    return cyc

def winding(D):
    w = [0, 0, 0]
    for (p, a) in D:
        if p[a] == 0:
            w[a] ^= 1
    return tuple(w)

def delta_E(old, new, L):
    """exact dE = E(new) - E(old) via local cycle deltas through changed bonds."""
    added = [b for b in new if b not in old]
    removed = [b for b in old if b not in new]
    dC4 = len(cycles_through_bonds(new, added, 4, L)) - len(cycles_through_bonds(old, removed, 4, L))
    dC6 = len(cycles_through_bonds(new, added, 6, L)) - len(cycles_through_bonds(old, removed, 6, L))
    return -W4 * dC4 - W6 * dC6, dC4, dC6

def n_cubes(B, L):
    """number of perfect isolated Q3 cells (spectral test) -- defect diagnostic."""
    adj = {}
    for (p, a) in B:
        q = wrap(add(p, AX[a]), L)
        adj.setdefault(p, []).append(q); adj.setdefault(q, []).append(p)
    SPEC = sorted([3, 1, 1, 1, -1, -1, -1, -3])
    seen, good = set(), 0
    for s in adj:
        if s in seen:
            continue
        comp, st = {s}, [s]
        while st:
            v = st.pop()
            for w in adj[v]:
                if w not in comp:
                    comp.add(w); st.append(w)
        seen |= comp
        if len(comp) == 8:
            idx = {v: i for i, v in enumerate(sorted(comp))}
            M = np.zeros((8, 8))
            for v in comp:
                for w in adj[v]:
                    M[idx[v], idx[w]] = 1.0
            if sorted(np.round(np.linalg.eigvalsh(M)).astype(int).tolist()) == SPEC:
                good += 1
    return good

# ================= [1] FREE WALL: x-phase flip in a z-slab => dE = 0 (degeneracy) =================
print("[1] FREE WALL (x-dimerisation flipped in a z-slab; interface perp to z, cube-aligned):")
L = 8
CR = tiling(L)
slab = set(range(0, L // 2))                                    # z in [0, L/2)
free = set()
for (p, a) in [(p, a) for p in ((x, y, z) for x in range(L) for y in range(L) for z in range(L))
               for a in range(3)]:
    if a == 0:                                                 # x-edges: flip phase inside the z-slab
        present = (p[0] % 2 == 0) ^ (p[2] in slab)
    else:
        present = (p[a] % 2 == 0)
    if present:
        free.add((p, a))
assert all(v == 3 for v in degree_map(free, L).values())       # still a valid 3-factor
dEfree, _, _ = delta_E(CR, free, L)
print(f"    dE(free wall) = {dEfree:+.1f} w6   (cubes: crystal={n_cubes(CR,L)}, walled={n_cubes(free,L)})")
assert abs(dEfree) < 1e-9
print("    -> dE = 0: a DISTINCT ground state. The cube crystal is DEGENERATE (any 2x2x2")
print("       partition), so compatible domain walls cost NOTHING -- they are entropic, not")
print("       durable energetic defects.")

# ================= [2] COSTED WALL: y-phase flip in a mid-cube x-region => sigma_wall =================
print("\n[2] COSTED WALL (y-dimerisation flipped for x in {1,2,3,4}; interfaces at x=1,5 cut cubes):")
region = {1, 2, 3, 4}                                          # boundaries x=1 and x=5 are mid-cube (odd)
costed = set()
for (p, a) in [(p, a) for p in ((x, y, z) for x in range(L) for y in range(L) for z in range(L))
               for a in range(3)]:
    if a == 1:                                                 # y-edges: flip phase inside the x-region
        present = (p[1] % 2 == 0) ^ (p[0] in region)
    else:
        present = (p[a] % 2 == 0)
    if present:
        costed.add((p, a))
assert all(v == 3 for v in degree_map(costed, L).values())
wind = winding(set(costed) ^ CR)
dEc, dC4, dC6 = delta_E(CR, costed, L)
n_walls = 2                                                    # interfaces at x=1 and x=5
area = L * L                                                   # y-z plane area per wall
sigma = dEc / (n_walls * area)
print(f"    dE(costed wall) = {dEc:+.1f} w6 (dC4={dC4:+d}, dC6={dC6:+d}); cubes lost = {n_cubes(CR,L)-n_cubes(costed,L)}")
print(f"    homology [D] = {wind}  -> energetic defect, NOT topologically locked (heals if barrier crossed)")
print(f"    sigma_wall = dE / (2 walls x {L}x{L}) = {sigma:.4f} w6 per unit lattice area")
assert dEc > 0 and wind == (0, 0, 0)

# ================= [3] STRING tension for comparison =================
print("\n[3] STRING tension (1D costed defect), for the abundance bridge:")
Ls = 6
CRs = tiling(Ls)
line = {((x, 0, 0), 0) for x in range(Ls)}
slip = set(CRs) ^ line
mu, dc4s, dc6s = delta_E(CRs, slip, Ls)
mu /= Ls
assert all(v == 3 for v in degree_map(slip, Ls).values()) and winding(slip ^ CRs) == (1, 0, 0)
print(f"    mu_string = {mu:.2f} w6/step = w4 + 4 w6   (locked: [D]=(1,0,0))")
assert abs(mu - (W4 + 4 * W6)) < 1e-9

# ================= [4] verdict =================
sig_cell = sigma * 4                                           # per 2x2 cube face (area 4)
print(f"""
[4] VERDICT — sigma_wall, and the abundance split it forces:
  * sigma_wall = {sigma:.3f} w6 per unit lattice area = {sig_cell:.2f} w6 per 2x2 cube-face,
    for the FRUSTRATED (cube-cutting) wall. The string tension is mu = w4 + 4 w6 = {W4+4*W6:.2f}
    w6/step. Both are EXACT (local cycle deltas).
  * BUT the cube ground state is DEGENERATE ([1]): compatible walls are FREE (dE=0). So the
    KZ wall network counted in k04_defect_network_abundance.py SPLITS:
      - FREE walls (compatible tiling reconfigurations): configurational ENTROPY, not durable
        gravitating relics -- they reconfigure at no energy cost (no barrier to anneal);
      - COSTED walls (frustrated, sigma_wall>0) + winding STRINGS (mu): the durable ENERGETIC
        defects that gravitate and are kinetically frozen.
  * CORRECTION to the recast: rho_dark is driven by the COSTED-defect length/area density, NOT
    the full wall count. The energy scale is set (sigma_wall, mu); the surviving question is the
    costed-defect FRACTION of the KZ network (how often domains meet frustrated, ~ the geodesic-
    rule mismatch rate), times xi(R)^-1, times the w6<->Lambda bridge for absolute units.
  * Net: the defect-network abundance has a clean energy scale now (sigma_wall, mu both exact);
    the open pieces narrow to (a) the costed/frustrated fraction of KZ junctions, (b) xi(R), and
    (c) w6<->Lambda. The degeneracy finding is itself a result: 'frozen domain structure' is
    partly free reconfiguration (entropy) and partly costed frustration (the gravitating relic).
exit 0""")
print("ALL ASSERTIONS PASSED — free wall is free (degeneracy); costed sigma_wall measured; mu=w4+4w6.")
