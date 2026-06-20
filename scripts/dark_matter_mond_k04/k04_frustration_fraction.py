#!/usr/bin/env python3
r"""K04 FRUSTRATED FRACTION f_frust — what fraction of KZ domain junctions are COSTED.

Closes the last combinatoric in the defect-network abundance
(k04_defect_network_abundance.py): rho_dark ~ sigma_wall * f_frust * n_wall(xi(R)),
where k04_wall_tension.py showed only FRUSTRATED (cube-cutting) walls cost energy
(sigma_wall = 2w4+12w6/face) while COMPATIBLE walls are FREE (the cube ground state
is degenerate). f_frust is the costed fraction.

THE GEOMETRIC CRITERION (from the cube-boundary parity, verified below):
  the tiling phase phi = (phi_x,phi_y,phi_z) in Z2^3 shifts the cube lattice in each
  axis. The n-axis cube boundaries sit at v_n ODD for phi_n=0 and v_n EVEN for
  phi_n=1 -- DISJOINT position sets. So two domains share an n-cube-boundary at their
  interface IFF phi_n agrees. Therefore:

     a domain junction perpendicular to axis n is COSTED  <=>  phi_A,n != phi_B,n
                                                   (the NORMAL phase component differs);
     transverse-only differences (phi_a, a!=n) are FREE (transversely-shifted but
     still clean isolated cubes), confirmed by an explicit dE=0 build.

CONSEQUENCE (the number): over i.i.d. Z2^3 domain phases,
     P(junction is a wall, any phi difference)        = 1 - 1/8 = 7/8,
     P(junction perp n is costed, phi_n differs)      = 1/2,
  so the COSTED FRACTION OF WALLS is
     f_frust = (1/2)/(7/8) = 4/7 ~ 0.571,   (free fraction 3/7).
  The costed-junction area density scales as ~ (3/2)/xi (codim-1, KZ).

exit 0 = parity criterion verified (disjoint boundaries); transverse wall built free
         (dE=0); f_frust = 4/7 by enumeration AND Monte Carlo; costed density ~1/xi.
"""
import math
import random
import numpy as np

W4, W6 = 1.7, 1.0
AX = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]

def wrap(p, L): return (p[0] % L, p[1] % L, p[2] % L)
def add(p, s): return (p[0] + s[0], p[1] + s[1], p[2] + s[2])

# ---------------- [1] the cube-boundary parity criterion ----------------
print("[1] PARITY CRITERION: do two phases share an n-cube-boundary? (cube lattice positions)")
def n_boundaries(phi_n, L):
    """positions v_n (between cell v_n and v_n+1) where the n-bond is ABSENT = a cube
    boundary, for a 1D dimerisation of phase phi_n. n-bond present iff (v_n - phi_n) even."""
    return frozenset(v for v in range(L) if (v - phi_n) % 2 == 1)
L = 8
b0, b1 = n_boundaries(0, L), n_boundaries(1, L)
print(f"    phi_n=0 cube boundaries at v_n in {sorted(b0)}")
print(f"    phi_n=1 cube boundaries at v_n in {sorted(b1)}")
print(f"    shared (phi_n agree)  : {sorted(n_boundaries(0,L) & n_boundaries(0,L))} (identical)")
print(f"    shared (phi_n differ) : {sorted(b0 & b1)}  -> EMPTY => no cube-aligned interface")
assert b0 & b1 == frozenset() and n_boundaries(0, L) == n_boundaries(0, L)
print("    -> phi_n differ => NO shared cube boundary => the junction MUST cut cubes => COSTED.")
print("       phi_n agree   => shared boundaries => a cube-aligned (free) interface exists.")

# ---------------- [1b] anchor: a transverse-only wall is FREE (dE=0) ----------------
print("\n[1b] ANCHOR (energy build): a cube-aligned junction with only a TRANSVERSE phase")
print("     difference is FREE -- confirm dE=0 (the costed case sigma_wall>0 is k04_wall_tension.py):")
def tiling(L):
    return {(p, a) for p in ((x, y, z) for x in range(L) for y in range(L) for z in range(L))
            for a in range(3) if p[a] % 2 == 0}
def nbrs(B, p, L):
    out = []
    for a in range(3):
        if (p, a) in B:
            out.append((wrap(add(p, AX[a]), L), AX[a]))
        m = wrap(add(p, tuple(-x for x in AX[a])), L)
        if (m, a) in B:
            out.append((m, tuple(-x for x in AX[a])))
    return out
def cyc_through(B, bonds, k, L):
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
CR = tiling(L)
# transverse wall: flip x-phase (transverse) across a z-perpendicular cube-aligned slab
slab = set(range(0, L // 2))
trans = set()
for p in ((x, y, z) for x in range(L) for y in range(L) for z in range(L)):
    for a in range(3):
        present = ((p[0] % 2 == 0) ^ (p[2] in slab)) if a == 0 else (p[a] % 2 == 0)
        if present:
            trans.add((p, a))
added = [b for b in trans if b not in CR]; removed = [b for b in CR if b not in trans]
dC4 = len(cyc_through(trans, added, 4, L)) - len(cyc_through(CR, removed, 4, L))
dC6 = len(cyc_through(trans, added, 6, L)) - len(cyc_through(CR, removed, 6, L))
dE = -W4 * dC4 - W6 * dC6
print(f"    transverse cube-aligned wall: dE = {dE:+.1f} w6 -> FREE. (costed cube-cutting wall:")
print(f"    sigma_wall = 2w4+12w6/face = {2*W4+12*W6:.1f}, see k04_wall_tension.py.)")
assert abs(dE) < 1e-9

# ---------------- [2] f_frust by exhaustive enumeration over Z2^3 x Z2^3 ----------------
print("\n[2] f_frust by EXHAUSTIVE enumeration of all phase pairs (phi_A, phi_B) in (Z2^3)^2:")
phases = [(a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]
walls = costed = 0
for pa in phases:
    for pb in phases:
        for n in range(3):                          # a junction perpendicular to axis n
            is_wall = (pa != pb)
            is_costed = (pa[n] != pb[n])
            walls += is_wall
            costed += is_costed
print(f"    over {len(phases)**2 * 3} (pair x normal) junctions: walls = {walls}, costed = {costed}")
print(f"    P(wall) = {walls}/{len(phases)**2*3} = {walls/(len(phases)**2*3):.4f} (= 7/8 = {7/8:.4f})")
print(f"    P(costed) = {costed}/{len(phases)**2*3} = {costed/(len(phases)**2*3):.4f} (= 1/2)")
f_frust = costed / walls
print(f"    f_frust = costed/walls = {costed}/{walls} = {f_frust:.6f}  (= 4/7 = {4/7:.6f})")
assert abs(f_frust - 4 / 7) < 1e-9
assert abs(walls / (len(phases) ** 2 * 3) - 7 / 8) < 1e-9
assert abs(costed / (len(phases) ** 2 * 3) - 1 / 2) < 1e-9

# ---------------- [3] Monte-Carlo cross-check + costed density ~ 1/xi ----------------
print("\n[3] MONTE-CARLO cross-check + costed-junction density vs correlation length xi:")
rng = random.Random(20260614)
w = c = 0
for _ in range(300000):
    pa = tuple(rng.getrandbits(1) for _ in range(3))
    pb = tuple(rng.getrandbits(1) for _ in range(3))
    n = rng.randrange(3)
    if pa != pb:
        w += 1
        if pa[n] != pb[n]:
            c += 1
print(f"    MC f_frust = {c}/{w} = {c/w:.4f}  (target 4/7 = {4/7:.4f})")
assert abs(c / w - 4 / 7) < 0.01
# density on a cell lattice with i.i.d. xi-blocks
M = 12
def densities(xi, seed):
    r = random.Random(seed)
    nb = M // xi
    bp = {(i, j, k): tuple(r.getrandbits(1) for _ in range(3))
          for i in range(nb) for j in range(nb) for k in range(nb)}
    ph = {(x, y, z): bp[(x // xi, y // xi, z // xi)]
          for x in range(M) for y in range(M) for z in range(M)}
    nw = nc = 0
    for cpos in ph:
        for n in range(3):
            d = wrap(add(cpos, AX[n]), M)
            if ph[cpos] != ph[d]:
                nw += 1
                if ph[cpos][n] != ph[d][n]:
                    nc += 1
    V = M ** 3
    return nw / V, nc / V
print(f"    {'xi':>3s} {'n_wall':>8s} {'n_costed':>9s} {'costed/wall':>12s} {'xi*n_costed':>12s}")
xis, ncs = [], []
for xi in (1, 2, 3, 4, 6):
    nw = np.mean([densities(xi, s)[0] for s in range(6)])
    nc = np.mean([densities(xi, s)[1] for s in range(6)])
    xis.append(xi); ncs.append(nc)
    print(f"    {xi:>3d} {nw:>8.4f} {nc:>9.4f} {nc/nw:>12.4f} {xi*nc:>12.4f}")
slope, logC = np.polyfit(np.log(xis), np.log(ncs), 1)
print(f"    costed density fit: n_costed ~ {math.exp(logC):.3f} * xi^({slope:+.3f}) (KZ codim-1 => -1)")
assert -1.2 < slope < -0.8

# ---------------- [4] assemble the abundance ----------------
sig_face = 2 * W4 + 12 * W6
print(f"""
[4] VERDICT — f_frust = 4/7, and the assembled abundance law:
  * RULE (verified): a KZ domain junction perp axis n is COSTED iff the domains differ
    in the NORMAL phase component phi_n; transverse differences are FREE (dE=0, [1b]).
    Proof: n-cube-boundaries are disjoint in position for phi_n=0 vs 1 ([1]).
  * f_frust = 4/7 ~ 0.571 of all KZ domain walls are costed (energetic, durable); the
    other 3/7 are FREE compatible reconfigurations (configurational entropy). EXACT by
    enumeration, confirmed by Monte Carlo.
  * the COSTED defect (gravitating relic) density is n_costed ~ {math.exp(logC):.2f}/xi
    (slope {slope:+.2f}), so the abundance law assembles as
        rho_dark  ~  sigma_wall * n_costed(xi(R))  ~  ({sig_face:.1f} w6/face) x (4/7) x n_wall(xi(R))
                  ~  (1/xi(R))   in shape, with sigma_wall = 2w4+12w6/face fixed exactly.
  * The three loose ends are now TWO: only xi(R) (the boot cooling law, dark_sector
    open frontier #2) and the w6<->Lambda_QCD bridge (frontier #3) remain. The costed
    fraction f_frust=4/7 and the tension sigma_wall are CLOSED. The abundance is a clean
    rho ~ (4/7) sigma_wall / xi(R) -- a falsifiable scaling surface, no orphan-island count.
exit 0""")
print("ALL ASSERTIONS PASSED — parity criterion; transverse wall free; f_frust = 4/7 (enum+MC); n_costed ~ 1/xi.")
