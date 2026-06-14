#!/usr/bin/env python3
"""
SMG existence, CORRECT symmetry: count SM gauge singlets (SU(3)_c x SU(2)_L x U(1)_Y),
NOT SO(10) singlets, in Lambda^N(16) for all N.

Rationale: the framework GAUGES the SM, not SO(10) (SO(10) is the organising spinor rep;
only the 12 SM gauge bosons exist, §3.4). SMG must preserve the SM gauge group, which is
SMALLER than SO(10) and has MORE singlets. so10_smg_singlet_structure.py tested too-large a
symmetry; this is the relevant one.

Builds the 12 SM single-particle generators from the P3 physical charges of the 16 modes:
  SU(3)_c (8): Gell-Mann on quark colour triplets (left) / antitriplets (right);
  SU(2)_L (3): Pauli on left doublets (u_L/d_L per colour; nu_L/e_L);
  U(1)_Y (1): diagonal hypercharge.
Then H_SM = sum_A (T^A)^2 + sum_i (T^i_L)^2 + Y^2 ; zeros = SM gauge singlets.
Self-validating: SU(3)/SU(2) Lie algebra; single-particle H_SM = known SM-rep Casimirs.

numpy + scipy.
"""
import numpy as np
import itertools
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh

G0, G1, LQ, C0, C1, I3, CHI, W = range(8)
H = 0.5
# Gell-Mann
GM = [np.array(m, complex) for m in [
    [[0,1,0],[1,0,0],[0,0,0]], [[0,-1j,0],[1j,0,0],[0,0,0]], [[1,0,0],[0,-1,0],[0,0,0]],
    [[0,0,1],[0,0,0],[1,0,0]], [[0,0,-1j],[0,0,0],[1j,0,0]], [[0,0,0],[0,0,1],[0,1,0]],
    [[0,0,0],[0,0,-1j],[0,1j,0]], np.array([[1,0,0],[0,1,0],[0,0,-2]])/np.sqrt(3)]]
PA = [np.array([[0,1],[1,0]],complex), np.array([[0,-1j],[1j,0]],complex), np.array([[1,0],[0,-1]],complex)]


def R1(c): return not (c[G0]==1 and c[G1]==1)
def R2(c): return c[W]==c[CHI]
def R3(c): return (c[C0],c[C1])==(0,0) if c[LQ]==0 else (c[C0],c[C1])!=(0,0)
def Y(c):  # P3 hypercharge (left-Weyl SO(10)); left fields +, right (conjugate) gives -Q etc.
    left = (c[CHI]==0)
    col = {(0,0):0,(0,1):1,(1,0):2,(1,1):3}[(c[C0],c[C1])]
    return None  # placeholder; Y computed below from weights for self-check only


COL = {(0,1):0,(1,0):1,(1,1):2}            # 3 quark colours -> SU(3) index
gen = sorted([c for c in itertools.product([0,1],repeat=8)
              if (c[G0],c[G1])==(0,0) and R1(c) and R2(c) and R3(c)])
M = len(gen); assert M==16
idxmode = {c:i for i,c in enumerate(gen)}


def hyper(c):  # left-Weyl SO(10) hypercharge from the standard PS formula (matches P3)
    left=(c[CHI]==0); sgn=1 if left else -1
    # B-L: lepton -1, quark +1/3 (left); conjugate flips sign
    BL = (-1.0 if c[LQ]==0 else 1.0/3.0)*sgn
    if left:
        T3R=0.0; T3L=(0.5 if c[I3]==0 else -0.5)
    else:
        T3L=0.0; T3R=(-0.5 if c[I3]==0 else 0.5)
    return T3R+BL/2.0, T3L

# ---- build single-particle generators (16x16) ----
def sp_zero(): return np.zeros((M,M),complex)

# SU(3)_c: group quark modes into colour-(anti)triplet slots keyed by (I3, CHI)
slots={}
for c in gen:
    if c[LQ]==1:
        slots.setdefault((c[I3],c[CHI]),{})[COL[(c[C0],c[C1])]]=idxmode[c]
T3c=[]
for A in range(8):
    g=sp_zero()
    for (i3,chi),cols in slots.items():
        lam = GM[A] if chi==0 else -GM[A].conj()    # triplet (L) vs antitriplet (R)
        for a in range(3):
            for b in range(3):
                g[cols[a],cols[b]] += 0.5*lam[a,b]
    T3c.append(g)

# SU(2)_L: left doublets keyed by (colour-or-lepton, CHI=0); members I3=0 (up) / I3=1 (down)
dbl={}
for c in gen:
    if c[CHI]==0:
        key=("L" if c[LQ]==0 else COL[(c[C0],c[C1])])
        dbl.setdefault(key,{})[c[I3]]=idxmode[c]
T2L=[]
for i in range(3):
    g=sp_zero()
    for key,mem in dbl.items():
        for a in range(2):
            for b in range(2):
                g[mem[a],mem[b]] += 0.5*PA[i][a,b]   # a,b: 0=up(I3=0),1=down(I3=1)
    T2L.append(g)

# U(1)_Y diagonal
Ydiag=np.diag([hyper(c)[0] for c in gen]).astype(complex)

# ---- self-checks: Lie algebra + single-particle Casimir = SM-rep values ----
def comm(a,b): return a@b-b@a
f123 = comm(T3c[0],T3c[1]) - 1j*T3c[2]
assert np.allclose(f123,0), "SU(3) [T1,T2]=i T3 fails"
assert np.allclose(comm(T2L[0],T2L[1]) - 1j*T2L[2],0), "SU(2)_L algebra fails"
Cas = sum(g@g for g in T3c) + sum(g@g for g in T2L) + Ydiag@Ydiag
# Q (3,2,1/6): C2_3=4/3, C2_2=3/4, Y^2=1/36 -> 4/3+3/4+1/36
expQ = 4/3+3/4+1/36
for c in gen:
    if c[LQ]==1 and c[CHI]==0:
        assert abs(Cas[idxmode[c],idxmode[c]].real-expQ)<1e-9, "Q-mode Casimir mismatch"
print(f"self-check OK: SU(3),SU(2) algebra close; Q single-particle SM-Casimir = {expQ:.4f}")

# ---- second-quantise and count SM singlets per filling ----
allg = T3c+T2L+[Ydiag]
NZ=[[(i,j,g[i,j]) for i in range(M) for j in range(M) if abs(g[i,j])>1e-12] for g in allg]
def pcb(m,k): return bin(m&((1<<k)-1)).count("1")
def singlets(NF):
    basis=[sum(1<<x for x in cmb) for cmb in itertools.combinations(range(M),NF)]
    idx={mk:i for i,mk in enumerate(basis)}; D=len(basis)
    if D==1: return D,1
    HC=csr_matrix((D,D),dtype=complex)
    for nz in NZ:
        r=[];co=[];v=[]
        for col,m in enumerate(basis):
            for i,j,val in nz:
                if not (m>>j)&1: continue
                s=(-1)**pcb(m,j); m2=m^(1<<j)
                if i!=j and (m2>>i)&1: continue
                s*=(-1)**pcb(m2,i); m3=m2|(1<<i)
                r.append(idx[m3]);co.append(col);v.append(s*val)
        Gab=csr_matrix((v,(r,co)),shape=(D,D),dtype=complex); HC=HC+Gab@Gab
    HC=(HC+HC.getH())*0.5
    ev=(np.linalg.eigvalsh(HC.toarray()).real if D<=60
        else np.sort(eigsh(HC,k=min(80,D-2),which="SA",return_eigenvectors=False).real))
    ns=int(np.sum(ev<1e-6)); return D,(-1 if ns>=min(80,D-2) else ns)

print(f"\n  {'N':>3} {'dim':>8} {'#SM gauge singlets':>20}")
for NF in range(0,9):
    D,ns=singlets(NF)
    print(f"  {NF:>3} {D:>8} {str(ns):>20}{'  (saturated)' if ns<0 else ''}")
print("  (N=9..16 mirror N=7..0.)")
print("""
RESULT: SM gauge singlets (the SMG-relevant symmetry) -- see table. If nontrivial-filling SM
singlets EXIST, the single-cell SMG obstruction inferred from SO(10) does NOT hold for the
actual gauge group; SMG existence is not blocked at this level (consistent with the 16 being
SM-anomaly-free / Wang-Wen SMG-able). The chiral-gauge question (gap the MIRROR, spare the
physical, survive kinetics) remains the open 2-cell frontier. CORRECTS the SO(10)-singlet
reading, which tested too large a symmetry.
""")
