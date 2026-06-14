#!/usr/bin/env python3
"""
Scrutiny of the "SO(10) pseudospin capstone" claim: the inter-site 120-channel quartic
-U sum_a |psi_A^T A_a psi_B|^2 Fierz-collapses to H = -U(S^2-64) - B.S, B=(2t,0,U), a giant
spin-8 with gap sqrt(U^2+4t^2) claimed "topologically protected, survives kinetic hopping".

Reproduces the 1,617,804-state sector; ED gap = sqrt(U^2+4t^2) EXACTLY (the Fierz->SU(2)
algebra is correct). But the diagnostics show it is NOT chiral SMG, it is a trivial collective
(Lipkin / eta-pairing) model:
  * U=0 -> gap = 2t: the gap exists with ZERO interaction = the free bonding/antibonding
    kinetic gap. "Survives kinetics" is vacuous (at large t the gap IS the kinetic gap).
  * <N_A>=16 at t=0: ground state is a TRIVIAL polarised product state (A full, B empty),
    not a symmetric paired SMG vacuum.
  * interaction and kinetics live in the SAME collective SU(2) algebra (T_+-, N) -> components
    of one field B; they cannot compete. The collapse to one spin is the signature of an
    all-to-all collective model, not a local chiral lattice theory.
  * the interaction is a vector-like A<->B (Dirac) pairing that gaps everything -> opposite of
    SMG (which gaps a chiral set without a partner/condensate, sparing the physical mode).
numpy + scipy.
"""

import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh

weights = []
for i in range(32):
    b = [1 if (i & (1 << j)) else -1 for j in range(5)]
    if b[0]*b[1]*b[2]*b[3]*b[4] == 1: weights.append(b)
weights = np.array(weights)

popcounts = np.zeros(65536, dtype=np.int64); W_sums = np.zeros((65536, 5), dtype=np.int64)
for i in range(16):
    valid = (np.arange(65536) & (1 << i)) > 0
    popcounts[valid] += 1; W_sums[valid] += weights[i]

B_groups = {}
for sB in range(65536):
    B_groups.setdefault((int(popcounts[sB]), tuple(int(x) for x in W_sums[sB])), []).append(sB)
keys_list = []
for sA in range(65536):
    key = (16 - int(popcounts[sA]), tuple(-int(x) for x in W_sums[sA]))
    if key in B_groups:
        for sB in B_groups[key]: keys_list.append((sA << 16) | sB)
keys = np.array(keys_list, dtype=np.uint64); keys.sort()
dim = len(keys); print(f"sector dim = {dim}")

sA = (keys >> np.uint64(16)).astype(np.int64); sB = (keys & np.uint64(0xFFFF)).astype(np.int64)
N_A = popcounts[sA]
rows, cols, data = [], [], []
for i in range(16):
    mask = 1 << i
    valid = ((sA & mask) == 0) & ((sB & mask) > 0)
    idx = np.where(valid)[0]
    if len(idx) == 0: continue
    mb = mask - 1
    cA = popcounts[sA[idx] & mb]; cB = popcounts[sB[idx] & mb]
    sA_new = (sA[idx] | mask).astype(np.uint64); sB_new = (sB[idx] ^ mask).astype(np.uint64)
    keys_new = (sA_new << np.uint64(16)) | sB_new
    idx_new = np.searchsorted(keys, keys_new)
    signs = 1 - 2*((cA + cB + N_A[idx]) % 2)
    rows.append(idx_new); cols.append(idx); data.append(signs)
T_plus = sp.csr_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))), shape=(dim, dim))
T_minus = T_plus.T
NANB = sp.diags((N_A*(16-N_A)).astype(float), format='csr')

def run(U, t):
    H = -U*(T_plus @ T_minus - NANB) - t*(T_plus + T_minus)
    ev, vec = eigsh(H, k=2, which='SA')
    o = np.argsort(ev); ev = ev[o]; gs = vec[:, o[0]]
    return ev[0], ev[1]-ev[0], float(np.sum(np.abs(gs)**2 * N_A))

print(f"\n  {'U':>4} {'t':>5} {'E0':>10} {'ED gap':>8} {'sqrt(U^2+4t^2)':>14} {'<N_A>':>7}")
for U, t in [(1,0),(1,1),(1,4),(1,10),(0,1),(0,4)]:
    e0, gap, na = run(float(U), float(t))
    print(f"  {U:>4} {t:>5} {e0:>10.3f} {gap:>8.3f} {np.sqrt(U**2+4*t**2):>14.3f} {na:>7.2f}")
