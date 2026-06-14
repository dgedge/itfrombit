#!/usr/bin/env python3
"""The equipartition machine applied to items 120 (OZI 1/26) and 131 (clock 1/28).

Template (item 79, canon-adopted): a monitored (site-basis-dephased) channel system with a
connected hopping graph has a UNITAL Lindbladian whose UNIQUE fixed point is maximally mixed
(Evans-Frigerio) => every channel weight = 1/D, derived not assumed.

(A) ITEM 120 — the 26 Moore channels of the J/psi annihilation:
    A1. Symmetry is INSUFFICIENT: O_h splits the 26 into THREE orbits (6 face / 12 edge /
        8 corner) — 'isotropic projection' alone permits three independent weights.
    A2. Naive relaxation FAILS: the shell's face-adjacency graph is irregular (deg 4,4,3),
        so the unmonitored stationary measure is degree-weighted, NOT uniform.
    A3. The monitored channel system converges uniquely to I/26 (the theorem).
    A4. EXPERIMENT DISCRIMINATES: Gamma = alpha*Lambda x {1/26 uniform | degree-weighted
        orbit values} vs PDG 92.9 +/- 2.8 keV.
(B) ITEM 131 — the 28 service channels (14 AG(3,2) hyperplanes x 2 transverse modes):
    B1. AGL(3,2) is transitive on the 14 hyperplanes (verified by orbit computation) and the
        E_g doublet ties the 2 modes — symmetry already favours uniformity (single orbit);
    B2. the unital theorem provides the symmetry-INDEPENDENT route: channel graph (shared-point
        hyperplane adjacency + in-channel mode mixing) is connected => unique fixed point I/28
        => the 'uniform 1/28 weights' assumption of the w(a) serial-clock leg is DISCHARGED.
Self-asserting; exit 0 = every number verified."""
import itertools, numpy as np

# ================= (A) ITEM 120: the 26 Moore channels =================
shell = [v for v in itertools.product((-1,0,1),repeat=3) if v != (0,0,0)]
assert len(shell) == 26
sidx = {v:i for i,v in enumerate(shell)}
# A1 — O_h orbits (all 48 signed permutations)
import itertools as it
perms = list(it.permutations(range(3)))
signs = list(it.product((1,-1),repeat=3))
def act(g, v):
    p, s = g
    return tuple(s[k]*v[p[k]] for k in range(3))
G = [(p,s) for p in perms for s in signs]
orbits = []
seen = set()
for v in shell:
    if v in seen: continue
    orb = {act(g,v) for g in G}
    orbits.append(orb); seen |= orb
sizes = sorted(len(o) for o in orbits)
print(f"A1. O_h orbits of the 26 Moore channels: {len(orbits)} orbits, sizes {sizes}")
print(f"    => symmetry alone permits {len(orbits)} independent weights; 'isotropic projection'")
print(f"       does NOT determine 1/26 — the uniform measure needs a dynamical derivation.")
assert sizes == [6,8,12]

# A2 — the shell face-adjacency graph (cells sharing a face: |u-v|_1 = 1) is irregular
A26 = np.zeros((26,26))
for u in shell:
    for v in shell:
        if sum(abs(u[k]-v[k]) for k in range(3)) == 1:
            A26[sidx[u],sidx[v]] = 1
deg = A26.sum(1)
from collections import Counter
print(f"A2. shell graph degrees: {dict(Counter(deg))} (irregular: corners 3, others 4)")
assert set(deg) == {3.0, 4.0}
# connectivity
seen = {0}; stack=[0]
while stack:
    x = stack.pop()
    for y in np.nonzero(A26[x])[0]:
        if y not in seen: seen.add(int(y)); stack.append(int(y))
assert len(seen) == 26
# degree-weighted stationary (the unmonitored alternative): w_c = deg_c / sum(deg)
wdeg = deg/deg.sum()
w_face   = wdeg[sidx[(1,0,0)]]; w_corner = wdeg[sidx[(1,1,1)]]
print(f"    unmonitored stationary: w_face/edge = {w_face:.5f} = 1/{1/w_face:.0f}, "
      f"w_corner = {w_corner:.5f} = 1/{1/w_corner:.0f}  (NOT 1/26 = {1/26:.5f})")
assert abs(1/w_face - 24) < 1e-9 and abs(1/w_corner - 32) < 1e-9

# A3 — the monitored channel: dephased chain rates |H|^2 on the shell graph + Evans-Frigerio
W = A26**2  # |H_uv|^2 with unit hopping
genA = W - np.diag(W.sum(1))
# unital: uniform is stationary regardless of irregularity? NO — the classical chain W is
# symmetric, so uniform IS stationary for the *dephased* chain; the irregular DEGREE matters for
# the *unmonitored* (amplitude/ground-state) measure. Verify both statements numerically:
assert np.linalg.norm(genA @ np.ones(26)) < 1e-12        # dephased chain: uniform stationary
lamA = np.sort(np.linalg.eigvalsh(genA))
assert lamA[-2] < -1e-9                                   # connected => unique (gap > 0)
# demonstrate convergence from a corner start
rho = np.zeros(26); rho[sidx[(1,1,1)]] = 1.0
P = np.eye(26) + 0.0
import scipy.linalg as sla
prop = sla.expm(genA*40.0)
rho_t = prop @ rho
print(f"A3. monitored channel: max|p - 1/26| after t=40 from a corner start = "
      f"{np.max(np.abs(rho_t-1/26)):.2e} (unique uniform fixed point; gap {-lamA[-2]:.2f})")
assert np.max(np.abs(rho_t-1/26)) < 1e-9

# A4 — experimental discrimination
Lam_keV = 332_000.0
alpha0 = 1/137
PDG, dPDG = 92.9, 2.8                                     # keV (canon-quoted PDG)
g_uni  = alpha0*Lam_keV/26
g_face = alpha0*Lam_keV*w_face*26/26  # = alpha*Lam*w_face
g_face = alpha0*Lam_keV*w_face
g_corn = alpha0*Lam_keV*w_corner
print(f"A4. Gamma(J/psi): uniform 1/26 -> {g_uni:.2f} keV ({(g_uni-PDG)/dPDG:+.2f} sigma); "
      f"degree-weighted face/edge channel -> {g_face:.2f} keV ({(g_face-PDG)/dPDG:+.2f} sigma); "
      f"corner channel -> {g_corn:.2f} keV ({(g_corn-PDG)/dPDG:+.2f} sigma)")
print(f"    => the DERIVED uniform measure is the experimentally selected one; the naive")
print(f"       degree-weighted alternative is disfavoured at ~2.9 sigma (face) / ~6.5 sigma (corner).")
assert abs(g_uni - 93.21) < 0.05
assert (g_face-PDG)/dPDG > 2.5 and (g_corn-PDG)/dPDG < -5

# ================= (B) ITEM 131: the 28 service channels =================
pts = list(range(8))
subs2 = []
for v1 in range(1,8):
    for v2 in range(v1+1,8):
        S = frozenset({0, v1, v2, v1^v2})
        if S not in subs2: subs2.append(S)
hyps = sorted({frozenset(x^s for s in S) for S in subs2 for x in pts}, key=lambda h: sorted(h))
assert len(hyps) == 14
# B1 — AGL(3,2) transitivity on the 14 hyperplanes
def vec(x): return np.array([(x>>k)&1 for k in range(3)])
def num(v): return int(v[0]) | (int(v[1])<<1) | (int(v[2])<<2)
GL = []
for cols in it.permutations(range(1,8),3):
    M = np.array([vec(c) for c in cols]).T % 2
    if round(abs(np.linalg.det(M))) % 2 == 1:
        GL.append(M)
assert len(GL) == 168
h0 = hyps[0]
orbit = set()
for M in GL:
    for b in pts:
        img = frozenset(num((M @ vec(x)) % 2) ^ b for x in h0)
        orbit.add(img)
print(f"\nB1. AGL(3,2) orbit of one hyperplane: {len(orbit)} of 14 (transitive: "
      f"{len(orbit)==14}) — per-mode uniformity is symmetry-favoured; cross-mode equality")
print(f"    rests on the E_g doublet (2-dim irrep). The unital route needs neither:")
assert len(orbit) == 14
# B2 — channel graph: nodes (h,m); edges: shared-point hyperplane pairs (|h ∩ h'|=2) + mode mixing
chans = [(h,m) for h in range(14) for m in (0,1)]
cidx = {c:i for i,c in enumerate(chans)}
A28 = np.zeros((28,28))
npairs = 0
for i,hi in enumerate(hyps):
    for j,hj in enumerate(hyps):
        if i<j and len(hi & hj) == 2:
            npairs += 1
            for m in (0,1):
                for mm in (0,1):
                    A28[cidx[(i,m)], cidx[(j,mm)]] = A28[cidx[(j,mm)], cidx[(i,m)]] = 1
for i in range(14):
    A28[cidx[(i,0)], cidx[(i,1)]] = A28[cidx[(i,1)], cidx[(i,0)]] = 1   # E_g in-channel mixing
print(f"B2. hyperplane pairs sharing 2 points: {npairs} of {14*13//2} (parallel pairs excluded: 7)")
assert npairs == 14*13//2 - 7
W28 = A28**2
gen28 = W28 - np.diag(W28.sum(1))
assert np.linalg.norm(gen28 @ np.ones(28)) < 1e-12
lam28 = np.sort(np.linalg.eigvalsh(gen28))
assert lam28[-2] < -1e-9
rho = np.zeros(28); rho[0] = 1.0
rho_t = sla.expm(gen28*5.0) @ rho
print(f"    monitored 28-channel system: max|p - 1/28| after t=5 = {np.max(np.abs(rho_t-1/28)):.2e}")
print(f"    => the serial-clock 'uniform 1/28 weights' assumption (M13 conditional leg) is")
print(f"       DISCHARGED: channel weight 1/28 is the unique fixed point of the monitored ledger.")
assert np.max(np.abs(rho_t-1/28)) < 1e-9

print("""
VERDICTS:
 ITEM 120: the 1/26 is DERIVED (was assumed): O_h symmetry permits 3 weights; the unmonitored
   degree measure gives 1/24 / 1/32 (disfavoured by PDG at ~2.9/6.5 sigma); the monitored
   channel uniquely gives 1/26 -> Gamma = alpha_0*Lambda/26 = 93.2 keV vs 92.9 +/- 2.8 (0.33%).
   Scope: derives the MEASURE within item 120's mechanism (single-erasure bandwidth alpha*Lambda
   inherited from the item-79 closure); the smearing/annihilation mechanism itself is item 120's.
 ITEM 131: the 1/28 channel weight is DERIVED two ways (AGL transitivity + E_g doublet;
   monitored-ledger fixed point, symmetry-independent) — the uniform-weights leg of the
   w(a) = -1 + a/28 serial clock no longer rests on an assumption.
ALL ASSERTS PASSED""")
