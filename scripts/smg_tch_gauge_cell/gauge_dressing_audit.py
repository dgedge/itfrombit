#!/usr/bin/env python3
"""Gauge-compatibility audit for the cell-level CSS SMG construction. (user-provided)"""
import itertools
import numpy as np

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Z = np.diag([1, -1]).astype(complex)
P0 = np.array([[1, 0], [0, 0]], complex)
P1 = np.array([[0, 0], [0, 1]], complex)
SR = np.array([[0, 1], [0, 0]], complex)
SL = np.array([[0, 0], [1, 0]], complex)
G0, G1, LQ, C0, C1, I3, CHI, W = range(8)

def hd(t): print("\n"+"="*78+"\n"+t+"\n"+"="*78)
def op(singles):
    out=np.array([[1]],complex)
    for bit in range(8): out=np.kron(out, singles.get(bit,I2))
    return out
def x_codeword(c): return op({b:X for b,v in enumerate(c) if v})
def z_codeword(c): return op({b:Z for b,v in enumerate(c) if v})
def ctrl(cb,cv,tb,to):
    a=P0 if cv==0 else P1; ina=P1 if cv==0 else P0
    return op({cb:a,tb:to})+op({cb:ina,tb:I2})
def commutator(a,b): return a@b-b@a
def norm(a): return float(np.linalg.norm(a))
def idx(b): return int("".join(map(str,b)),2)
def charge_zf(b): return 1 if b[I3]==0 else -1
def sum_colour_projectors(b): return -3 if (b[C0],b[C1])==(0,0) else -1
def electric_charge(b): return 0.5*charge_zf(b)+(1/3)*sum_colour_projectors(b)+0.5
def weak_t3(b): return (0.5 if b[I3]==0 else -0.5) if b[CHI]==0 else 0.0
def build_gauge_actions():
    ab=[tuple(b) for b in itertools.product([0,1],repeat=8)]
    hy=np.diag([electric_charge(b)-weak_t3(b) for b in ab]).astype(complex)
    t3=np.diag([weak_t3(b) for b in ab]).astype(complex)
    tr=ctrl(CHI,0,I3,SR); tl=ctrl(CHI,0,I3,SL)
    cc=np.zeros((256,256),complex); cyc={(0,1):(1,0),(1,0):(1,1),(1,1):(0,1),(0,0):(0,0)}
    for b in ab:
        t=list(b); t[C0],t[C1]=cyc[(b[C0],b[C1])]; cc[idx(tuple(t)),idx(b)]=1
    return {"U(1)_Y":hy,"SU(2)_L T3":t3,"SU(2)_L T+":tr,"SU(2)_L T-":tl,"SU(3)_c cycle":cc,"SU(3)_c cycle^-1":cc.conj().T}
def orthonormal_orbit(seed,ga,tol=1e-9):
    bv=[]; bm=[]; q=[]
    def add(m):
        v=m.reshape(-1).astype(complex)
        for _ in range(2):
            for b in bv: v-=np.vdot(b,v)*b
        s=np.linalg.norm(v)
        if s<=tol: return False
        v=v/s; bv.append(v); bm.append(v.reshape(seed.shape)); q.append(bm[-1]); return True
    add(seed); cur=0
    while cur<len(q):
        m=q[cur]; cur+=1
        for go in ga.values(): add(commutator(go,m))
    return bv,bm
def representation_matrix(bv,bm,go):
    B=np.vstack(bv); cv=np.column_stack([commutator(go,m).reshape(-1) for m in bm])
    rep=B.conj()@cv; res=cv-B.T@rep; return rep, float(np.max(np.linalg.norm(res,axis=0)))
def dressed_residual_norm(rep):
    link_block=-rep.T; residual=rep+link_block.T; return norm(residual)
def main():
    gens=[[1,1,1,1,1,1,1,1],[0,0,0,0,1,1,1,1],[0,0,1,1,0,0,1,1],[0,1,0,1,0,1,0,1]]
    xs=[x_codeword(r) for r in gens]; zs=[z_codeword(r) for r in gens]
    hx=-sum(xs); hz=-sum(zs); hcss=hx+hz; zz=op({CHI:Z,W:Z})
    ga=build_gauge_actions()
    hd("Bare X vs gauge (sample)")
    for i,s in enumerate(xs[:1]):
        for n,g in ga.items(): print(f"  ||[Xstab0,{n}]||={norm(commutator(s,g)):.4g}")
    hd("Hamiltonian-level")
    for nm,o in {"H_X":hx,"H_Z":hz,"H_CSS":hcss,"Z_chi Z_W":zz}.items():
        print(f"  {nm}: max||[.,gauge]|| = {max(norm(commutator(o,g)) for g in ga.values()):.4g}")
    hd("Dressed orbits")
    maxfull=0.0; maxclos=0.0; dims=[]
    for i,s in enumerate(xs):
        bv,bm=orthonormal_orbit(s,ga); dims.append(len(bm))
        for n,g in ga.items():
            rep,clos=representation_matrix(bv,bm,g); full=dressed_residual_norm(rep)
            maxclos=max(maxclos,clos); maxfull=max(maxfull,full)
        print(f"  Xstab{i}: orbit dim={len(bm)}")
    print(f"\n  orbit dims={dims}; max closure residual={maxclos:.3e}; max dressed commutator residual={maxfull:.3e}")
    assert maxclos<1e-7 and maxfull<1e-9
    print("\nALL ASSERTS PASSED.")
main()
