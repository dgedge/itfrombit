#!/usr/bin/env python3
"""
Stress-test of a parallel claim that an SO(6)~=SU(4) quartic -|U| sum_a |psi^T A_a psi|^2
(the 6 antisymmetric vector-channel forms) "symmetrically gaps" a 2-cell strip with a gap
that "survives arbitrary kinetic hopping" -> SMG victory.

Reproduces the claim's spectrum AND adds the discriminating checks. Result: it is NOT SMG:
 (1) the gap exists ONLY in the attractive U<0 (condensate) branch; the SYMMETRIC U>0 branch
     (penalise pairing -> singlet/kernel sector) is GAPLESS (gap=0 at small t). So the gap is
     a pairing condensate <O_a> breaking SO(6)->SO(5) (SSB, "mass WITH symmetry breaking"),
     the NJL/BCS mechanism -- not SMG (mass WITHOUT symmetry breaking).
 (2) at large t the U=-1 gap == the U=0 FREE-fermion gap (88.1 ~ 100 at t=100): the surviving
     gap is trivial bonding-antibonding kinetics, not the interaction.
 (3) SO(6) vector 6 is antisymmetric in 4x4; the SO(10) vector 10 is SYMMETRIC (bilinear
     vanishes, cf. so10_smg_operator_check.py) -- SO(10)'s antisymmetric channel is the 120,
     not the vector. The analogy does not transfer.
Self-contained; numpy+scipy.
"""
import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh

def build_so6():
    X=np.array([[0,1],[1,0]]);Y=np.array([[0,-1j],[1j,0]]);Z=np.array([[1,0],[0,-1]]);I=np.eye(2)
    def kl(o):
        r=o[0]
        for x in o[1:]: r=np.kron(r,x)
        return r
    G=[]
    for k in range(3):
        G.append(kl([Z]*k+[X]+[I]*(2-k))); G.append(kl([Z]*k+[Y]+[I]*(2-k)))
    G7=kl([Z]*3); idx=np.where(np.diag(G7)>0)[0]; C=G[1]@G[3]@G[5]
    Ma=[]
    for i in range(6):
        M=(C@G[i])[np.ix_(idx,idx)]
        if np.linalg.norm(M+M.T)<1e-10: Ma.append(M)
    return Ma
def ac(s,m):
    if not(s&(1<<m)): return None,0
    return s&~(1<<m), (-1 if bin(s&((1<<m)-1)).count('1')%2 else 1)
def ad(s,m):
    if s&(1<<m): return None,0
    return s|(1<<m), (-1 if bin(s&((1<<m)-1)).count('1')%2 else 1)

Ma=build_so6(); print(f"# antisymmetric SO(6) matrices (the vector 6): {len(Ma)}")
dim=256
Hs=sp.lil_matrix((dim,dim),dtype=complex)
for M in Ma:
    for off in [0,4]:
        O=sp.lil_matrix((dim,dim),dtype=complex)
        for i in range(4):
            for j in range(4):
                if abs(M[i,j])>1e-10:
                    for s in range(dim):
                        s1,p1=ac(s,off+j)
                        if s1 is not None:
                            s2,p2=ac(s1,off+i)
                            if s2 is not None: O[s2,s]+=M[i,j]*p1*p2
        O=O.tocsr(); Hs+=O.T.conj()@O
Hh=sp.lil_matrix((dim,dim),dtype=complex)
for i in range(4):
    for s in range(dim):
        s1,p1=ac(s,4+i)
        if s1 is not None:
            s2,p2=ad(s1,i)
            if s2 is not None: Hh[s2,s]-=p1*p2
        s1,p1=ac(s,i)
        if s1 is not None:
            s2,p2=ad(s1,4+i)
            if s2 is not None: Hh[s2,s]-=p1*p2
Hs=Hs.tocsr(); Hh=Hh.tocsr()

def gaps(U,t,k=6):
    H=U*Hs+t*Hh
    ev=np.sort(eigsh(H+sp.eye(dim)*1e-6,k=k,which='SA',return_eigenvectors=False).real)
    return ev[0]-1e-6, ev

print("\n  t  | U=-1 (their condensate) E0/gap | U=0 (FREE) gap | U=+1 (SMG kernel) E0/gap")
for t in [0.0,1.0,4.0,16.0,32.0,100.0]:
    e_m,evm=gaps(-1.0,t); e_f,evf=gaps(0.0,t); e_p,evp=gaps(1.0,t)
    gm=evm[1]-evm[0]; gf=evf[1]-evf[0]; gp=evp[1]-evp[0]
    print(f" {t:6.1f} | E0={e_m:8.2f} gap={gm:7.3f} | free gap={gf:7.3f} | E0={e_p:7.3f} gap={gp:7.3f}")

# degeneracy at t=0, U=-1 (SSB signature = degenerate/near-degenerate ground space)
_,ev0=gaps(-1.0,0.0,k=12)
print(f"\n U=-1, t=0 lowest 12 eigenvalues: {np.round(ev0,3)}")
print(" (a unique symmetric ground state should be non-degenerate; an SSB-prone condensate")
print("  branch shows (near-)degeneracy = competing symmetry-broken directions.)")
