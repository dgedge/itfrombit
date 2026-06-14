import itertools, numpy as np
EDGES=[(0,1),(1,2),(3,2),(0,3)]
PERMS=list(itertools.permutations(range(3))); P_INDEX={p:i for i,p in enumerate(PERMS)}
IDENTITY=P_INDEX[(0,1,2)]; RESIDUAL=[IDENTITY,P_INDEX[(0,2,1)]]
NONIDENTITY=[i for i in range(6) if i!=IDENTITY]
def compose(l,r): p=PERMS[l]; q=PERMS[r]; return P_INDEX[tuple(p[q[i]] for i in range(3))]
def inverse(pi):
    p=PERMS[pi]; inv=[0,0,0]
    for s,t in enumerate(p): inv[t]=s
    return P_INDEX[tuple(inv)]
INV=[inverse(i) for i in range(6)]
def transform(cfg,vr): return tuple(compose(vr[t],compose(cfg[l],INV[vr[h]])) for l,(t,h) in enumerate(EDGES))
def orbit(seed): return {transform(seed,r) for r in itertools.product(RESIDUAL,repeat=4)}
def build_orbits():
    rem=set(itertools.product(range(6),repeat=4)); orbits=[]; c2o={}
    while rem:
        seed=next(iter(rem)); cur=orbit(seed); oi=len(orbits); orbits.append(sorted(cur))
        for c in cur: c2o[c]=oi
        rem-=cur
    return orbits,c2o
def trace_norm(pi): return sum(1 for i,v in enumerate(PERMS[pi]) if v==i)/3
def preserves_ref(pi): return 1.0 if PERMS[pi][0]==0 else 0.0
def plaq(cfg): return compose(cfg[0],compose(cfg[1],compose(INV[cfg[2]],INV[cfg[3]])))
def diag_E(cfg,beta,kappa): return -kappa*sum(preserves_ref(l) for l in cfg)-beta*trace_norm(plaq(cfg))
def shifted(cfg,l,g): o=list(cfg); o[l]=compose(g,o[l]); return tuple(o)
def H(beta,kappa=1.0,elec_scale=1.0):
    g2=1/beta; orbits,c2o=build_orbits(); d=len(orbits); h=np.zeros((d,d),complex)
    for si,configs in enumerate(orbits):
        osz=len(configs); h[si,si]+=diag_E(configs[0],beta,kappa)
        tc={}
        for cfg in configs:
            for l in range(4):
                for g in NONIDENTITY:
                    ti=c2o[shifted(cfg,l,g)]; tc[ti]=tc.get(ti,0)+1
        for ti,cnt in tc.items():
            h[ti,si]+=elec_scale*(-g2/len(NONIDENTITY))*cnt/np.sqrt(osz*len(orbits[ti]))
    return h,orbits
def gap(beta,kappa=1.0,elec_scale=1.0):
    h,_=H(beta,kappa,elec_scale); ev=np.linalg.eigvalsh(h); return float(ev[1]-ev[0])
orbits,_=build_orbits()
print(f"raw S3 link configs = {6**4}; residual Gauss orbits = {len(orbits)}")
print("\nReproduce gaps (kappa=1):")
for beta in [0.25,0.5,0.75,1.0]:
    print(f"  beta={beta}: gap={gap(beta):.5f}")
print("\nDIAGNOSTIC 1 — electric (gauge) scaling: gap vs 1/beta (g2)?")
for beta in [0.25,0.5,1.0]:
    print(f"  beta={beta}: gap={gap(beta):.4f}   g2=1/beta={1/beta:.4f}   gap*beta={gap(beta)*beta:.4f}")
print("\nDIAGNOSTIC 2 — turn off the electric term (elec_scale->0), beta=0.5:")
for es in [1.0,0.5,0.1,0.0]:
    print(f"  elec_scale={es}: gap={gap(0.5,elec_scale=es):.4f}")
print("\nDIAGNOSTIC 3 — content of H: terms = electric(S3 links) + Higgs(-kappa preserves-ref)")
print("  + plaquette(-beta Re Tr). NO fermion occupation, NO CSS X/Z stabilizers, NO mirror modes.")
print("  -> pure S3 gauge-Higgs. The gap is the S3 gauge/confinement gap, not a mirror gap.")
