#!/usr/bin/env python3
"""Item 79 question (a): is T(16)=136 derivable, or is the antisymmetric 120 right?

RESOLUTION CANDIDATE (tested here): K9's '"Grassmann" should give 120' objection neglected the
COIN. The framework's matter excitations carry a 2-dim internal doublet (the 3.5 spinor/chirality
coin). For TWO IDENTICAL fermionic defects each with a 2-dim coin, overall antisymmetry splits:

    antisym( 16 nodes x 2 coin ) = [Sym^2(16) (x) coin-SINGLET]  (+)  [Antisym^2(16) (x) coin-TRIPLET]
            C(32,2) = 496        =      136 x 1                  +        120 x 3

The coin-singlet (epsilon-contracted, scalar) channel — the one that couples to S-wave photon
emission (1.5's annihilation normalisation; C-parity) — has spatial dimension EXACTLY T(16)=136,
INCLUDING the 16 on-site states |i,up; i,down> (legal for fermions: same node, opposite coin —
no F2 double-occupancy). Grassmann statistics, done with the coin, FORCES the symmetric count.
Self-asserting; exit 0 = every number verified."""
import itertools, numpy as np

N, C = 16, 2
M = N*C                                     # one-particle modes (node, coin)
def mode(i,s): return 2*i+s

# two-fermion antisymmetric (wedge) basis
wedge = [(m1,m2) for m1 in range(M) for m2 in range(m1+1,M)]
W = len(wedge)
assert W == 496 == M*(M-1)//2

# embed wedge basis in the full 2-particle tensor space (dim M^2) to act with projectors
def wedge_vec(m1,m2):
    v = np.zeros(M*M)
    v[m1*M+m2] = 1/np.sqrt(2); v[m2*M+m1] = -1/np.sqrt(2)
    return v
B = np.column_stack([wedge_vec(*p) for p in wedge])      # (M^2, 496), orthonormal

# coin-singlet projector on the pair: |eps><eps| with eps = (|01>-|10>)/sqrt2 per node-pair,
# i.e. P_sing = I_space (x) |eps><eps| acting on (node1,coin1,node2,coin2)
P = np.zeros((M*M, M*M))
for i in range(N):
    for j in range(N):
        # basis vectors |i s1> |j s2>
        def t(s1,s2): return mode(i,s1)*M + mode(j,s2)
        # eps component on (coin1,coin2): (|01>-|10>)/sqrt2
        for (a1,a2),ca in (((0,1),1/np.sqrt(2)),((1,0),-1/np.sqrt(2))):
            for (b1,b2),cb in (((0,1),1/np.sqrt(2)),((1,0),-1/np.sqrt(2))):
                P[t(a1,a2), t(b1,b2)] += ca*cb
Psub = B.T @ P @ B                                       # singlet projector on the wedge space
assert np.linalg.norm(Psub@Psub - Psub) < 1e-10          # projector
rank_singlet = int(round(np.trace(Psub)))
rank_triplet = W - rank_singlet
print(f"two identical fermions, 16 nodes x 2-dim coin: antisym dim = {W}")
print(f"  coin-SINGLET sector: dim = {rank_singlet}   (T(16) = {N*(N+1)//2})")
print(f"  coin-TRIPLET sector: dim = {rank_triplet}   (3 x C(16,2) = {3*N*(N-1)//2})")
assert rank_singlet == 136 and rank_triplet == 360

# the 16 on-site states are in the singlet block (same node, opposite coin)
onsite = [wedge.index((mode(i,0),mode(i,1))) for i in range(N)]
ons_w = sum(Psub[k,k] for k in onsite)
print(f"  on-site pairs |i,up; i,down>: {len(onsite)} states, total singlet weight = {ons_w:.6f} (=16: all in singlet)")
assert abs(ons_w-16) < 1e-9

# spatial pair-walk (coin-blind) commutes with the singlet projector; singlet block == the
# symmetric-space H2 used in item79_unital_channel.py (same spectrum)
edges  = [(i,(i+1)%8) for i in range(8)] + [(8+i,8+(i+1)%8) for i in range(8)] + [(0,8),(1,9)]
A1 = np.zeros((N,N))
for i,j in edges: A1[i,j]=A1[j,i]=1
H1 = np.kron(A1, np.eye(C))                              # coin-blind one-particle walk on 32 modes
H2full = np.kron(H1,np.eye(M)) + np.kron(np.eye(M),H1)
Hw = B.T @ H2full @ B
assert np.linalg.norm(Hw@Psub - Psub@Hw) < 1e-9          # singlet sector invariant
# singlet-block spectrum vs the symmetric-space pair walk
evP, VP = np.linalg.eigh(Psub)
S = VP[:, evP > 0.5]                                     # orthonormal basis of the singlet block
assert S.shape[1] == 136
ev_sing = np.sort(np.linalg.eigvalsh(S.T @ Hw @ S))
pairs_s = [(i,j) for i in range(N) for j in range(i,N)]
Bas = np.zeros((N*N,len(pairs_s)))
for k,(i,j) in enumerate(pairs_s):
    v = np.zeros((N,N))
    if i==j: v[i,i]=1.0
    else: v[i,j]=v[j,i]=1/np.sqrt(2)
    Bas[:,k]=v.reshape(-1)
H2sym = Bas.T @ (np.kron(A1,np.eye(N))+np.kron(np.eye(N),A1)) @ Bas
ev_sym = np.sort(np.linalg.eigvalsh(H2sym))
assert np.max(np.abs(ev_sing-ev_sym)) < 1e-8
print(f"  singlet-block walk spectrum == symmetric-space H2 spectrum (max dev {np.max(np.abs(ev_sing-ev_sym)):.1e})")
print(f"  => the item79_unital_channel 'symmetric' computation WAS the coin-singlet sector.")

print(f"""
THE THREE READINGS, NOW WITH THEIR PREMISES EXPLICIT:
  reading                          space  alpha_0^-1  premise required
  coin-singlet (photon channel)     136      137      identical register defects + 2-dim coin
                                                      + photon couples to scalar/singlet (S-wave)
  spatial-antisym (coinless)        120      121      coinless point fermions (CONTRADICTS 3.5's
                                                      spinor structure on matter excitations)
  distinguishable species           256      257      e+/e- as distinct registers (CONTRADICTS
                                                      2.7/2.8 same-register CPT structure)
VERDICT (a): with the framework's own coin, Grassmann statistics FORCES the photon-channel pair
space to be Sym^2(16) = T(16) = 136 — the diagonal pairs are the on-site coin-singlets, and K9's
counting objection is RESOLVED in favour of 137 under named canon premises (3.5 coin; 2.7/2.8
same-register defects; 1.5 S-wave emission channel). The 120 exists — as the coin-triplet sector,
which does not couple to the scalar emission mode (ortho- vs para-positronium structure).
CAVEAT (kept): ANCHOR 5.4's Shell-1 capacity uses C(8,2)=28 (no diagonal) — a DIFFERENT pairing
convention than Shell-2's T(16); under this script's logic those are different objects (direction
pairs, no coin antisymmetrisation), but the per-shell convention difference deserves its own audit.
ALL ASSERTS PASSED""")
