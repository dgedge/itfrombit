#!/usr/bin/env python3
"""
L1 basis-robust confirmation for the ANCHOR item-96 caution.

tch_256_walk_mass_rank_audit.py tests the colour blocks of the documented
256-state walk mass operator M2=(U^dag U)^2 in the generation COMPUTATIONAL
(G0,G1) basis. The proposed 1/3 mechanism's colour blocks live in the generation
MASS-EIGENBASIS, so this script redoes the test there -- the fair decomposition.

Result (asserted): in the mass eigenbasis (eigvecs of the colour-traced generation
operator T), every per-generation conditional colour block is FULL RANK (three
comparable singular values) and far from proportional to I_3. No block is rank-one,
none is the identity -- so L1 (C_light rank-one, C_others=I_3) is NEGATIVE for the
documented operator, basis-independently. (T is also near-diagonal in the
computational basis, so the two bases nearly coincide -- the audit's basis is fine.)
"""
import numpy as np
G0,G1,LQ,C0,C1,I3,CHI,Wb=range(8); DELTA=2/9; N=256
def i2b(i): return [(i>>(7-b))&1 for b in range(8)]
def b2i(b):
    o=0
    for x in b: o=(o<<1)|x
    return o
def cnot_k(i,k):
    b=i2b(i); c=(2-k)%8; t=(5-k)%8
    if b[c]==1: b[t]^=1
    return b2i(b)
amp=[np.sqrt(1-DELTA)]+[np.sqrt(DELTA/7)*np.exp(1j*k*np.pi/4) for k in range(1,8)]
U=np.zeros((N,N),complex)
for s in range(N):
    for k,a in enumerate(amp): U[cnot_k(s,k),s]+=a
M2=(U.conj().T@U)@(U.conj().T@U)
GENS=[(0,0),(0,1),(1,0)]; COLS=[(0,1),(1,0),(1,1)]
def qbasis(iso):
    out=[]
    for g in GENS:
        for c in COLS:
            st=[0]*8; st[G0],st[G1]=g; st[LQ]=1; st[C0],st[C1]=c; st[I3]=iso
            out.append(b2i(st))
    return out
def analyze(name,iso):
    M9=M2[np.ix_(qbasis(iso),qbasis(iso))]; blk=M9.reshape(3,3,3,3)
    T=np.einsum('gcHc->gH',blk)/3
    offd=max(abs(T[a,b]) for a in range(3) for b in range(3) if a!=b)
    dmin=min(abs(T[a,a]) for a in range(3))
    w,V=np.linalg.eigh(T); order=np.argsort(w.real)
    print(f"\n{name}: T near-diagonal? off/diag={offd/dmin:.3f}; masses={np.round(w[order].real,4)}")
    for n,oi in enumerate(order):
        v=V[:,oi]; Cn=np.einsum('g,gcHd,H->cd',v.conj(),blk,v)
        sv=np.linalg.svd(Cn,compute_uv=False); rk=int(np.sum(sv>1e-9))
        tau=np.trace(Cn).real/3; resid=float(np.max(np.abs(Cn-tau*np.eye(3))))
        tag=["light","middle","heavy"][n]
        print(f"  {tag:<7} rank={rk} svals={np.round(sv,4)} tau={tau:.4f} ||C-tauI||={resid:.4f}")
        assert rk==3                       # full rank, NOT rank-one
        assert sv[-1]>0.5                  # smallest singular value far from 0
        assert resid>0.1                   # NOT proportional to I_3
    assert offd/dmin<0.06                  # bases nearly coincide (audit basis is fine)
analyze("up-left  (iso=0)",0)
analyze("down-left(iso=1)",1)
print("\nALL ASSERTS PASSED -> L1 negative for the documented operator, basis-independently.")
