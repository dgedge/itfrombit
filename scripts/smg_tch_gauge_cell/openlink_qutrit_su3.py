import itertools, numpy as np
OMEGA=np.exp(2j*np.pi/3)
EDGES=[(0,1),(1,2),(3,2),(0,3)]; ORIENTATION=[1,1,-1,-1]
def cos_z3(v): return float(np.cos(2*np.pi*(v%3)/3))
def split_config(c): return tuple(c[:4]),tuple(c[4:])
def covariant_links(matter,links):
    return tuple((matter[h]-matter[t]+links[l])%3 for l,(t,h) in enumerate(EDGES))
def plaquette_flux(b): return sum(s*v for s,v in zip(ORIENTATION,b))%3
def gauss_transform(c,sh):
    m,lk=split_config(c); nm=tuple((m[v]+sh[v])%3 for v in range(4))
    nl=[(v+sh[t]-sh[h])%3 for v,(t,h) in zip(lk,EDGES)]
    return nm+tuple(nl)
def canon(c): return min(gauss_transform(c,sh) for sh in itertools.product(range(3),repeat=4))
def b_basis():
    bs=list(itertools.product(range(3),repeat=4)); return bs,{s:i for i,s in enumerate(bs)}
def shift_b(s,l,a):
    o=list(s); o[l]=(o[l]+a)%3; return tuple(o)
def matter_shift_b(s,vx,a):
    o=list(s)
    for l,(t,h) in enumerate(EDGES):
        if h==vx: o[l]=(o[l]+a)%3
        if t==vx: o[l]=(o[l]-a)%3
    return tuple(o)
def hamiltonian(beta,kappa=1.0,eta=0.0,elec_scale=1.0):
    g2=1/beta; bs,ix=b_basis(); d=len(bs); h=np.zeros((d,d),complex)
    for col,s in enumerate(bs):
        h[col,col]+=-kappa*sum(cos_z3(v) for v in s) - beta*cos_z3(plaquette_flux(s))
        for l in range(4):
            for a in [1,-1]: h[ix[shift_b(s,l,a)],col]+=-elec_scale*g2/2
        for vx in range(4):
            for a in [1,-1]: h[ix[matter_shift_b(s,vx,a)],col]+=-eta/2
    return h
def gap(beta,kappa=1.0,eta=0.0,elec_scale=1.0):
    ev=np.linalg.eigvalsh(hamiltonian(beta,kappa,eta,elec_scale)); return float(ev[1]-ev[0])
# A. reproduce orbit audit
cfgs=list(itertools.product(range(3),repeat=8))
print("orbit audit: raw=%d, Gauss orbits=%d, covariant b=%d"%(
    len(cfgs), len({canon(c) for c in cfgs}), len({covariant_links(*split_config(c)) for c in cfgs})))
# B. reproduce the reported gaps
print("\nReproduce reported gaps:")
for beta,eta in [(0.5,0),(1.0,0),(0.5,0.5),(1.0,0.5)]:
    print(f"  beta={beta} eta={eta}: gap={gap(beta,eta=eta):.5f}")
# DIAGNOSTIC 1: is it the electric gap? compare to pure-electric Z3 link gap = 3*g2/2 = 3/(2beta)
print("\nDIAGNOSTIC 1 — is the gap the ELECTRIC gauge gap (= 3/(2beta))?")
for beta in [0.25,0.5,0.75,1.0]:
    print(f"  beta={beta}: reported gap(eta=0)={gap(beta,eta=0):.4f}   3/(2beta)={3/(2*beta):.4f}")
# DIAGNOSTIC 2: kill the electric term (elec_scale->0). If gap collapses, it WAS the electric gap.
print("\nDIAGNOSTIC 2 — gap with electric term scaled down (eta=0, beta=0.5):")
for es in [1.0,0.5,0.1,0.01]:
    print(f"  elec_scale={es:<5}: gap={gap(0.5,eta=0,elec_scale=es):.4f}")
print("\nDIAGNOSTIC 3 — what is in the Hamiltonian? terms = electric(links) + Higgs(cos b) +")
print("  magnetic(cos flux) + matter-orbit(eta). There are NO CSS X/Z stabilizers and NO")
print("  mirror fermions in this 81-dim space — 'matter' is the Z3 colour qutrit, projected into b_l.")
