#!/usr/bin/env python3
import itertools, numpy as np
G0,G1,LQ,C0,C1,I3,CHI,W=range(8)
GENERATORS=[[1,1,1,1,1,1,1,1],[0,0,0,0,1,1,1,1],[0,0,1,1,0,0,1,1],[0,1,0,1,0,1,0,1]]
ALL_BITS=[tuple(b) for b in itertools.product([0,1],repeat=8)]
OMEGA=np.exp(2j*np.pi/3)
CL={(0,0):"ell",(0,1):"r",(1,0):"g",(1,1):"b"}; ORD=["ell","r","g","b"]; CI={l:i for i,l in enumerate(ORD)}
def hd(t): print("\n"+"="*78+"\n"+t+"\n"+"="*78)
def idx(b): return int("".join(map(str,b)),2)
def clab(b): return CL[(b[C0],b[C1])]
def projs():
    P={}
    for l in ORD:
        m=np.zeros((4,4),complex); m[CI[l],CI[l]]=1; P[l]=m
    return P
def lift(local):
    M=np.zeros((256,256),complex)
    for b in ALL_BITS:
        s=CI[clab(b)]
        for tgt,co in enumerate(local[:,s]):
            if abs(co)<1e-14: continue
            tb=list(b); tl=ORD[tgt]
            tb[C0],tb[C1]={"ell":(0,0),"r":(0,1),"g":(1,0),"b":(1,1)}[tl]
            M[idx(tuple(tb)),idx(b)]+=co
    return M
def pflip(bit):
    m=np.zeros((4,4),complex)
    for c0,c1 in itertools.product([0,1],repeat=2):
        s=CL[(c0,c1)]; t=[c0,c1]; t[bit]^=1; m[CI[CL[tuple(t)]],CI[s]]=1
    return m
def ccycle():
    m=np.zeros((4,4),complex)
    for s,t in {"ell":"ell","r":"g","g":"b","b":"r"}.items(): m[CI[t],CI[s]]=1
    return m
def qshift():
    m=np.zeros((3,3),complex)
    for s in range(3): m[(s+1)%3,s]=1
    return m
def qops(): return np.diag([1,OMEGA,OMEGA**2]).astype(complex), qshift()
hd("A. F2 colour -> qutrit one-hot algebra")
P=projs(); pcol=P["r"]+P["g"]+P["b"]; Q=P["r"]+OMEGA*P["g"]+OMEGA**2*P["b"]
rec={"r":(pcol+Q+Q.conj().T)/3,"g":(pcol+OMEGA**2*Q+OMEGA*Q.conj().T)/3,"b":(pcol+OMEGA*Q+OMEGA**2*Q.conj().T)/3}
mrec=max(np.linalg.norm(rec[l]-P[l]) for l in "rgb"); print(f"  one-hot reconstruction residual={mrec:.1e}")
x0=pflip(0); x1=pflip(1); xb=x0@x1; lep=P["ell"]
leak={"X_C0":np.linalg.norm(pcol@x0@lep)+np.linalg.norm(lep@x0@pcol),"X_C1":np.linalg.norm(pcol@x1@lep)+np.linalg.norm(lep@x1@pcol),"X_C0X_C1":np.linalg.norm(pcol@xb@lep)+np.linalg.norm(lep@xb@pcol)}
print(f"  Pauli leakage: {dict((k,round(v,3)) for k,v in leak.items())}")
cyc=ccycle(); print(f"  cycle preserves coloured subspace={np.linalg.norm((np.eye(4)-pcol)@cyc@pcol):.1e}; coloured block==qutrit shift={np.linalg.norm(cyc[1:,1:]-qshift()):.1e}")
assert mrec<1e-12 and min(leak.values())>1.0
def xcw(cw):
    M=np.zeros((256,256),complex); r=tuple(cw)
    for b in ALL_BITS: M[idx(tuple(x^f for x,f in zip(b,r))),idx(b)]=1
    return M
def zcw(cw):
    z=np.diag([1,-1]).astype(complex); I=np.eye(2,dtype=complex); M=np.array([[1]],complex)
    for v in cw: M=np.kron(M,z if v else I)
    return M
def ccyc_u(): return lift(ccycle())
def adj_comp(op,U):
    pw=[np.linalg.matrix_power(U,n) for n in range(3)]; ip=[p.conj().T for p in pw]; comps=[]
    for ch in range(3):
        c=np.zeros_like(op)
        for n in range(3): c+=(OMEGA**(-ch*n))*pw[n]@op@ip[n]
        comps.append(c/3)
    return comps
def qdress(comps):
    clock,_=qops(); out=np.zeros((768,768),complex)
    for ch,c in enumerate(comps): out+=np.kron(c,np.linalg.matrix_power(clock,ch))
    return out
hd("C. Full dressed CSS (both halves) — computed local gap")
cyc_u=ccyc_u(); _,cs=qops(); tot=np.kron(cyc_u,cs); I768=np.eye(768,dtype=complex)
dz=[qdress(adj_comp(zcw(r),cyc_u)) for r in GENERATORS]
dx=[qdress(adj_comp(xcw(r),cyc_u)) for r in GENERATORS]
bz=[np.kron(zcw(r),np.eye(3,dtype=complex)) for r in GENERATORS]
mbzx=max(np.linalg.norm(z@x-x@z) for z in bz for x in dx)
alld=dz+dx; mcomm=max(np.linalg.norm(a@b-b@a) for a in alld for b in alld)
proj=I768.copy()
for s in alld: proj=proj@((I768+s)/2)
grank=round(float(np.real(np.trace(proj))))
gp=(I768+tot+tot@tot)/3; prank=round(float(np.real(np.trace(proj@gp))))
H=-sum(alld)-2.0*(tot+tot.conj().T)/2; ev=np.linalg.eigvalsh(H)
uniq=[]
for v in ev[:20]:
    if not uniq or abs(v-uniq[-1])>1e-8: uniq.append(float(v))
lgap=uniq[1]-uniq[0]; ldeg=int(np.sum(np.isclose(ev,ev[0],atol=1e-8)))
print(f"  ||[bare Z, dressed X]||={mbzx:.4g}  (dressing only X insufficient)")
print(f"  same-qutrit dressed comm={mcomm:.1e}; rank before Gauss={grank}; after orbit-Gauss={prank}")
print(f"  orbit-lock lambda=2: ground degeneracy={ldeg}, COMPUTED local gap={lgap:.4g}")
assert mbzx>1.0 and mcomm<1e-10 and grank==3 and prank==1 and ldeg==1 and lgap>1.9
def su3_reps(L): return [(p,q) for p in range(L+1) for q in range(L+1-p)]
def su3_cas(p,q): return (p*p+q*q+p*q+3*p+3*q)/3
def pauli_noflip(beta,level):
    reps=su3_reps(level); ix={r:i for i,r in enumerate(reps)}; d=len(reps); M=np.zeros((d,d),complex)
    for col,(p,q) in enumerate(reps):
        for t in [(p+1,q),(p-1,q+1),(p,q-1)]:
            if t in ix: M[ix[t],col]+=1
    W=(M+M.conj().T)/6; g2=1/beta
    def e0(src):
        H=np.diag([g2*su3_cas(*r) for r in reps]).astype(complex)-(beta+src)*W
        return float(np.linalg.eigvalsh(H)[0])
    return min(2.0, e0(0)-e0(4))
hd("D. Strong-coupling bookkeeping — IS the qutrit matter gap COMPUTED?")
print("  beta   Pauli no-flip (computed)   qutrit matter (?)")
for beta,level in [(0.5,18),(1,18),(2,20)]:
    pg=pauli_noflip(beta,level)
    qutrit_matter_gap=2.0     # <-- look: this is the line
    print(f"  {beta:<6g} {pg:<25.6g} {qutrit_matter_gap}")
print("  ^ the 'qutrit matter' value is the literal assignment qutrit_matter_gap=2.0 (HARDCODED).")
print("ALL ASSERTS PASSED (A real; C local gap computed; D qutrit column hardcoded).")
