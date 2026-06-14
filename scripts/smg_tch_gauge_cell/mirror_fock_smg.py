import itertools, numpy as np
G0,G1,LQ,C0,C1,I3,CHI,W=range(8)
GENERATORS=[[1,1,1,1,1,1,1,1],[0,0,0,0,1,1,1,1],[0,0,1,1,0,0,1,1],[0,1,0,1,0,1,0,1]]
N=8; DIM=256
ALL_BITS=[tuple((s>>(N-1-b))&1 for b in range(N)) for s in range(DIM)]
def idx(bits): return sum(bits[b]<<(N-1-b) for b in range(N))
def jw_majorana(bit):
    M=np.zeros((DIM,DIM),complex)
    for src,bits in enumerate(ALL_BITS):
        t=list(bits); t[bit]^=1; sign=(-1)**sum(bits[:bit]); M[idx(t),src]=sign
    return M
def parity(bit): return np.diag([1 if b[bit]==0 else -1 for b in ALL_BITS]).astype(complex)
def prod(ops):
    M=np.eye(DIM,dtype=complex)
    for o in ops: M=M@o
    return M
def css():
    g=[jw_majorana(b) for b in range(N)]; p=[parity(b) for b in range(N)]
    xs=[]; zs=[]
    for row in GENERATORS:
        sup=[b for b,v in enumerate(row) if v]; xs.append(prod([g[b] for b in sup])); zs.append(prod([p[b] for b in sup]))
    return zs,xs
def mirror_idx(): return [i for i,b in enumerate(ALL_BITS) if b[CHI]!=b[W]]
def gap_of(H):
    ev=np.linalg.eigvalsh(H); 
    # gap to first level above ground
    lv=[]
    for v in ev:
        if not lv or abs(v-lv[-1])>1e-8: lv.append(float(v))
    deg=int(np.sum(np.isclose(ev,ev[0],1e-8)))
    return lv[1]-lv[0], lv[0], deg
zs,xs=css()
allst=zs+xs
mc=max(np.linalg.norm(a@b-b@a) for a in allst for b in allst); mh=max(np.linalg.norm(s-s.conj().T) for s in allst)
print(f"stabilizers: max commutator={mc:.1e}, max non-herm={mh:.1e}  (even fermion ops, commuting)")
H=-sum(allst)
g,e0,deg=gap_of(H); print(f"full 8-mode Fock CSS : dim=256 E0={e0:.4g} gap={g:.4g} deg={deg}")
mi=mirror_idx()
gm,e0m,degm=gap_of(H[np.ix_(mi,mi)]); print(f"mirror W!=chi (Z+X)  : dim={len(mi)} E0={e0m:.4g} gap={gm:.4g} deg={degm}")
xn=[np.linalg.norm(x[np.ix_(mi,mi)]) for x in xs]; print(f"  projected X norms: {np.round(xn,4)}")

print("\nDIAGNOSTIC — does the gap come from the SMG (X) interaction, or from Z alone?")
HZ=-sum(zs); HX=-sum(xs)
gz,e0z,degz=gap_of(HZ[np.ix_(mi,mi)]); print(f"  mirror, Z-only (NO SMG X): gap={gz:.4g} E0={e0z:.4g} ground degeneracy={degz}")
gx,e0x,degx=gap_of(HX[np.ix_(mi,mi)]); print(f"  mirror, X-only (SMG only): gap={gx:.4g} E0={e0x:.4g} ground degeneracy={degx}")
gzx,_,degzx=gap_of((HZ+HX)[np.ix_(mi,mi)]); print(f"  mirror, Z+X (full SMG)   : gap={gzx:.4g} ground degeneracy={degzx}")
print("  -> if Z-only is degenerate/gapless and Z+X is unique gap 2, the SMG X-interaction is")
print("     what gaps the mirror (genuine SMG), not the trivial diagonal stabilizer.")
