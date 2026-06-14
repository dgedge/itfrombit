from collections import defaultdict
import itertools
import numpy as np
G0,G1,LQ,C0,C1,I3,CHI,W=range(8)
GENERATORS=[[1,1,1,1,1,1,1,1],[0,0,0,0,1,1,1,1],[0,0,1,1,0,0,1,1],[0,1,0,1,0,1,0,1]]
N_MODES=8; DIM=256; OMEGA=np.exp(2j*np.pi/3)
ALL_BITS=[tuple((s>>(N_MODES-1-b))&1 for b in range(N_MODES)) for s in range(DIM)]
_LC={}
def idx(bits): return sum(bits[b]<<(N_MODES-1-b) for b in range(N_MODES))
def jw_ann(bit):
    M=np.zeros((DIM,DIM),complex)
    for s,bits in enumerate(ALL_BITS):
        if bits[bit]==0: continue
        t=list(bits); t[bit]=0; M[idx(t),s]=(-1)**sum(bits[:bit])
    return M
def jw_maj(bit): return jw_ann(bit)+jw_ann(bit).conj().T
def parity(bit): return np.diag([1 if b[bit]==0 else -1 for b in ALL_BITS]).astype(complex)
def prod(ops):
    M=np.eye(DIM,dtype=complex)
    for o in ops: M=M@o
    return M
def mirror_idx(): return [i for i,b in enumerate(ALL_BITS) if b[CHI]!=b[W]]
def color_cycle_full():
    cyc={(0,0):(0,0),(0,1):(1,0),(1,0):(1,1),(1,1):(0,1)}; M=np.zeros((DIM,DIM),complex)
    for s,b in enumerate(ALL_BITS):
        t=list(b); t[C0],t[C1]=cyc[(b[C0],b[C1])]; M[idx(t),s]=1
    return M
def qops():
    clock=np.diag([1,OMEGA,OMEGA**2]).astype(complex); shift=np.zeros((3,3),complex)
    for s in range(3): shift[(s+1)%3,s]=1
    return clock,shift
def adj_comp(op,U):
    pw=[np.linalg.matrix_power(U,p) for p in range(3)]; comps=[]
    for ch in range(3):
        c=np.zeros_like(op)
        for p in range(3): c+=OMEGA**(-ch*p)*pw[p]@op@pw[p].conj().T
        comps.append(c/3)
    return comps
def local_block(spc):
    if spc in _LC: return _LC[spc]
    mir=mirror_idx(); cf=color_cycle_full(); color=cf[np.ix_(mir,mir)]; clock,shift=qops()
    tc=np.kron(color,shift); ld=len(mir)*3
    g=[jw_maj(b) for b in range(N_MODES)]; p=[parity(b) for b in range(N_MODES)]
    ds=[]
    for row in GENERATORS:
        sup=[b for b,v in enumerate(row) if v]
        for full in [prod([p[b] for b in sup]), prod([g[b] for b in sup])]:
            op=full[np.ix_(mir,mir)]; d=np.zeros((ld,ld),complex)
            for ch,comp in enumerate(adj_comp(op,color)): d+=np.kron(comp,np.linalg.matrix_power(clock,ch))
            ds.append(d)
    hl=-sum(ds); assert np.linalg.norm(hl@tc-tc@hl)<1e-10
    ev,evec=np.linalg.eigh(hl); locs=[]; cur=0
    while cur<len(ev):
        end=cur+1
        while end<len(ev) and abs(ev[end]-ev[cur])<1e-8: end+=1
        bv=evec[:,cur:end]; cb=bv.conj().T@tc@bv; cval,cvec=np.linalg.eig(cb)
        for pos,val in enumerate(cval):
            ch=int(np.round((np.angle(val)%(2*np.pi))/(2*np.pi/3)))%3
            vec=bv@cvec[:,pos]; vec/=np.linalg.norm(vec); locs.append((float(ev[cur]),ch,vec))
        cur=end
    sel=[]
    for ch in range(3):
        sec=sorted([s for s in locs if s[1]==ch],key=lambda x:x[0])[:spc]; sel.extend(sec)
    sel=sorted(sel,key=lambda x:(x[0],x[1])); basis=np.column_stack([s[2] for s in sel])
    energies=np.array([s[0] for s in sel]); charges=np.array([s[1] for s in sel],dtype=int)
    anns=[]
    for bit in range(6):
        cm=jw_ann(bit)[np.ix_(mir,mir)]; cl=np.kron(cm,np.eye(3,dtype=complex)); comps=[]
        for ch in range(3):
            c=np.zeros_like(cl); pw=[np.linalg.matrix_power(tc,p) for p in range(3)]
            for p in range(3): c+=OMEGA**(-ch*p)*pw[p]@cl@pw[p].conj().T
            c/=3; comps.append(basis.conj().T@c@basis)
        anns.append(comps)
    cc={ch:int(np.sum(charges==ch)) for ch in range(3)}; _LC[spc]=(energies,charges,anns,cc); return _LC[spc]
def link_fluxes(charges):
    f=[]; cum=0
    for ch in charges[:-1]: cum=(cum+ch)%3; f.append((-cum)%3)
    return tuple(f)
def chain_basis(charges,L):
    ld=len(charges); basis=[]
    for cfg in itertools.product(range(ld),repeat=L):
        if sum(int(charges[i]) for i in cfg)%3!=0: continue
        basis.append((cfg,link_fluxes(tuple(int(charges[i]) for i in cfg))))
    return basis,{it:p for p,it in enumerate(basis)}
def build_chain(L,beta,t,spc):
    energies,charges,anns,cc=local_block(spc); basis,index=chain_basis(charges,L); dim=len(basis); g2=1/beta
    H=np.zeros((dim,dim),complex); Hm=np.zeros((dim,dim),complex)
    for col,(cfg,fl) in enumerate(basis):
        me=sum(energies[s] for s in cfg); el=g2*sum(0.0 if f==0 else 1.0 for f in fl)
        Hm[col,col]+=me; H[col,col]+=me+el
    for col,(cfg,_) in enumerate(basis):
        cfg=list(cfg)
        for link in range(L-1):
            ls,rs=cfg[link],cfg[link+1]
            for mc in anns:
                for ch in range(3):
                    a=mc[ch]; lt=np.nonzero(np.abs(a[:,ls])>1e-12)[0]; rt=np.nonzero(np.abs(a[:,rs])>1e-12)[0]
                    for nl in lt:
                        al=np.conj(a[nl,ls])
                        for nr in rt:
                            ar=a[nr,rs]; nc=list(cfg); nc[link]=int(nl); nc[link+1]=int(nr)
                            if sum(int(charges[s]) for s in nc)%3!=0: continue
                            nf=link_fluxes(tuple(int(charges[s]) for s in nc)); row=index.get((tuple(nc),nf))
                            if row is not None: H[row,col]+=-t*al*ar
    H=(H+H.conj().T)/2; assert np.linalg.norm(H-H.conj().T)<1e-10
    return H,cc,{"matter":Hm}
def matter_gap(H,Hm):
    ev,evec=np.linalg.eigh(H); gm=float(np.real(evec[:,0].conj()@Hm@evec[:,0]))
    for i in range(1,len(ev)):
        if float(np.real(evec[:,i].conj()@Hm@evec[:,i]))-gm>0.5: return float(ev[i]-ev[0])
    return float('inf')
def full_gap(H): ev=np.linalg.eigvalsh(H); return float(ev[1]-ev[0])
