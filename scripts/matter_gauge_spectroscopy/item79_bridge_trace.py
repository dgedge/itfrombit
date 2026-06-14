#!/usr/bin/env python3
"""Item 79 (i)-leg audit: does the WALK DYNAMICS on the 2-octagon gauge bridge actually
realise the equipartition that turns the count T(16)+1=137 into alpha = Tr_non-unit[W|K2]?

Method (vacpol-style: compute the actual object):
 1. Build the canonical 16-node interaction graph (two C8 octagons + 2 bridge edges):
    V=16, E=18, beta1=3, bridge-square interface n_b=4 — the Part-12 C2=-4!/7 ingredients.
 2. Counting audit: the three natural 2-particle state-counts -> only one yields 137.
 3. Equipartition audit: (a) classical stationary measure (degree-weighted, graph is
    irregular 3-vs-2); (b) quantum time-averaged (diagonal-ensemble) measure from the
    physical start (pair meeting at the bridge). Equipartition over the 136 confined
    states is what 1/137 branching NEEDS; measure its actual violation.
Self-asserting. exit 0 = every quoted number verified."""
import itertools, numpy as np

# ---- 1. the graph ----
N = 16
edges  = [(i,(i+1)%8) for i in range(8)]                     # octagon A: 0..7
edges += [(8+i,8+(i+1)%8) for i in range(8)]                 # octagon B: 8..15
edges += [(0,8),(1,9)]                                       # bridge (interface {0,1,8,9})
A = np.zeros((N,N))
for i,j in edges: A[i,j]=A[j,i]=1
deg = A.sum(1)
V,E = N, len(edges); beta1 = E-V+1
nb  = 4                                                       # bridge-square interface nodes
assert (V,E,beta1)==(16,18,3) and sorted(set(deg))==[2.0,3.0]
assert [deg[i] for i in (0,1,8,9)]==[3,3,3,3]
import math
assert math.factorial(nb)==24 and 2**beta1-1==7              # C2 = -24/7 ingredients
print(f"1. graph: V={V} E={E} beta1={beta1}; degrees: four 3s (interface), twelve 2s")
print(f"   Part-12 C2 ingredients verified: n_b!={math.factorial(nb)}, 2^beta1-1={2**beta1-1}")

# Laplacian spectral gap (the canon's cited justification)
L = np.diag(deg)-A
lam = np.sort(np.linalg.eigvalsh(L))
print(f"   Laplacian gap lambda_2 = {lam[1]:.4f} > 0 (connected/ergodic) — but gap>0 gives")
print(f"   ERGODICITY, not EQUIDISTRIBUTION: stationary measure is degree-weighted on an")
print(f"   irregular graph (uniformity needs regularity).")
assert lam[1] > 0.05

# ---- 2. counting audit ----
c_fermi  = N*(N-1)//2          # antisymmetric pairs (identical fermions)  = 120
c_tri    = N*(N+1)//2          # unordered WITH diagonal (T(16))           = 136
c_dist   = N*N                 # distinguishable pair (e+ e-)              = 256
print(f"\n2. counting audit (+1 emission channel):")
print(f"   antisym (Grassmann-consistent) {c_fermi} -> {c_fermi+1}")
print(f"   T(16) unordered-with-diagonal  {c_tri} -> {c_tri+1}   <- the ONLY one giving 137")
print(f"   distinguishable (e+e-)         {c_dist} -> {c_dist+1}")
assert (c_fermi+1, c_tri+1, c_dist+1) == (121, 137, 257)
print(f"   NOTE: the count branded 'Grassmann' uses the NON-fermionic convention (diagonal")
print(f"   pairs included); the fermionic count gives 121, not 137 — formula-freedom flag.")

# ---- 3a. classical equipartition ----
pi1 = deg/deg.sum()
ratio1 = pi1.max()/pi1.min()
print(f"\n3a. classical stationary measure: per-walker max/min = {ratio1:.2f}; pair-state")
print(f"    max/min = {ratio1**2:.4f} (= 9/4). Equipartition over pair states FAILS by 2.25x.")
assert abs(ratio1-1.5)<1e-12

# ---- 3b. quantum diagonal-ensemble measure on the 136 sym-pair states ----
pairs = [(i,j) for i in range(N) for j in range(i,N)]
idx   = {p:k for k,p in enumerate(pairs)}
assert len(pairs)==136
def sym(i,j):
    v = np.zeros(136)
    v[idx[(min(i,j),max(i,j))]] = 1.0
    return v
# two-walker Hamiltonian H2 = A(x)I + I(x)A restricted to the symmetric subspace:
# build the symmetric-subspace pair Hamiltonian cleanly from the full 256-dim lift:
H_full = np.kron(A,np.eye(N)) + np.kron(np.eye(N),A)
# symmetriser basis
Bas = np.zeros((256,136))
for (i,j),k in idx.items():
    v = np.zeros((N,N))
    if i==j: v[i,j]=1.0
    else:    v[i,j]=v[j,i]=1/np.sqrt(2)
    Bas[:,k] = v.reshape(256)
H2 = Bas.T @ H_full @ Bas
assert np.linalg.norm(H2-H2.T)<1e-12
w2, U2 = np.linalg.eigh(H2)
# physical start: the pair meeting ACROSS the gauge bridge: sym(0,8)
phi0 = Bas.T @ (np.kron(np.eye(N)[0],np.eye(N)[8]) + np.kron(np.eye(N)[8],np.eye(N)[0]))/np.sqrt(2)
phi0 /= np.linalg.norm(phi0)
# diagonal ensemble with degeneracy projectors: pbar[s] = sum_lambda |<s|P_lambda|phi0>|^2
tol = 1e-9
pbar = np.zeros(136); k0=0
while k0 < 136:
    k1 = k0
    while k1+1<136 and abs(w2[k1+1]-w2[k0])<tol: k1 += 1
    P = U2[:,k0:k1+1]
    proj = P @ (P.T @ phi0)
    pbar += proj**2
    k0 = k1+1
assert abs(pbar.sum()-1) < 1e-9
uni = 1/136
print(f"\n3b. quantum diagonal-ensemble occupation (start = pair at the bridge, sym(0,8)):")
print(f"    max p = {pbar.max():.5f} ({pbar.max()/uni:.1f}x uniform)   min p = {pbar.min():.2e} "
      f"({pbar.min()/uni:.3f}x uniform)")
print(f"    fraction of states below half-uniform: {(pbar<uni/2).mean()*100:.0f}%")
# the would-be emission weight if the portal is the bridge-interface pair class:
portal = [idx[(min(a,b),max(a,b))] for a in (0,1,8,9) for b in (0,1,8,9) if a<=b]
w_port = pbar[portal].sum()
print(f"    bridge-interface pair-class weight = {w_port:.4f} vs equipartition {len(portal)/136:.4f}"
      f" ({w_port/(len(portal)/136):.1f}x)")
assert pbar.max()/uni > 3.0          # strongly non-uniform
assert w_port/(len(portal)/136) > 2.0

print(f"""
VERDICT (item 79, (i)-leg):
  - The count T(16)+1=137 is clean as a COUNT, but its state-space convention
    (unordered-with-diagonal) is one of three natural choices and the only one giving 137;
    the fermionic ('Grassmann') convention gives 121.
  - The equipartition that converts the count into a dynamical branching 1/137 is NOT
    produced by the walk on the actual graph: classical stationary fails by 2.25x
    (irregular 3-vs-2 degrees), and the quantum diagonal ensemble from the physical
    bridge start is strongly non-uniform ({pbar.max()/uni:.0f}x peak, bridge class {w_port/(len(portal)/136):.1f}x over-weighted).
    The cited justification (positive spectral gap) gives ergodicity, not equidistribution.
  - Status mirror: K5 (dressed-alpha count-vs-integral). 'alpha = Tr_non-unit[W|K2]' is
    currently a COUNT IDENTIFICATION, not a trace-of-dynamics theorem.
  - Constructive promotion criterion (sharper than 'formal proof'): derive equipartition
    as the unique fixed point of the bridge-sector CPTP channel (the item-119 dissipator
    restricted to the bridge), i.e. show the OPEN dynamics, not the unitary walk,
    uniformises the 137-state sector. Unitary W cannot do it; the dissipator might.
ALL ASSERTS PASSED""")
