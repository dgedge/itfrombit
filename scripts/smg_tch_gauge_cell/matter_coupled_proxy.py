#!/usr/bin/env python3
import itertools
import numpy as np
GROUND_X=(1,1,1,1); GLOBAL_FLIP_X=(-1,-1,-1,-1); BETA_VALUES=[0.5,1,2,5,10,20,50,100]
def hd(t): print("\n"+"="*78+"\n"+t+"\n"+"="*78)
def parity(xs):
    out=1
    for v in xs: out*=v
    return out
def x_sectors(closed_flux=None):
    s=list(itertools.product([1,-1],repeat=4))
    if closed_flux is not None: s=[xs for xs in s if parity(xs)==closed_flux]
    return s
def su3_reps(level): return [(p,q) for p in range(level+1) for q in range(level+1-p)]
def su3_casimir(p,q): return (p*p+q*q+p*q+3*p+3*q)/3
class U1Rotor:
    name="U(1) rotor"; center_quotient=True
    def __init__(self,n_max=120):
        ns=np.arange(-n_max,n_max+1); self.electric_diag=ns**2/2; d=len(ns)
        self.wilson=np.zeros((d,d),complex)
        for i in range(d-1): self.wilson[i,i+1]=0.5; self.wilson[i+1,i]=0.5
    def eigs(self,rho,g2):
        H=np.diag(g2*self.electric_diag).astype(complex)-rho*self.wilson
        return np.linalg.eigvalsh(H)[:2]
class SU2Character:
    name="SU(2) character"; center_quotient=True
    def __init__(self,n_max=240):
        self.electric_diag=np.array([(n/2)*(n/2+1) for n in range(n_max+1)])
        self.wilson=np.zeros((n_max+1,n_max+1),complex)
        for n in range(n_max): self.wilson[n,n+1]=0.5; self.wilson[n+1,n]=0.5
    def eigs(self,rho,g2):
        H=np.diag(g2*self.electric_diag).astype(complex)-rho*self.wilson
        return np.linalg.eigvalsh(H)[:2]
class SU3Character:
    name="SU(3) character shadow"; center_quotient=False
    def __init__(self,level):
        self.level=level; reps=su3_reps(level); idx={r:i for i,r in enumerate(reps)}; d=len(reps)
        M=np.zeros((d,d),complex)
        for col,(p,q) in enumerate(reps):
            for t in [(p+1,q),(p-1,q+1),(p,q-1)]:
                if t in idx: M[idx[t],col]+=1
        self.electric_diag=np.array([su3_casimir(*r) for r in reps]); self.wilson=(M+M.conj().T)/6; self.dim=d
    def eigs(self,rho,g2):
        H=np.diag(g2*self.electric_diag).astype(complex)-rho*self.wilson
        assert np.linalg.norm(H-H.conj().T)<1e-10
        return np.linalg.eigvalsh(H)[:2]
def gauge_spectra_by_source(model,beta,kappa):
    g2=1/beta
    return {src: model.eigs(beta+kappa*src,g2) for src in [-4,-2,0,2,4]}
def min_x_gap(spectra,cands):
    ge=float(spectra[4][0]); gaps=[]
    for xs in cands:
        if xs==GROUND_X: continue
        gaps.append(float(spectra[sum(xs)][0]-ge))
    return min(gaps)
def diagnostics(model,beta,kappa=1.0):
    spectra=gauge_spectra_by_source(model,beta,kappa); gs=spectra[4]
    gauge_gap=float(gs[1]-gs[0]); z=2.0
    raw=min_x_gap(spectra,x_sectors()); closed=min_x_gap(spectra,x_sectors(closed_flux=+1))
    qc=[xs for xs in x_sectors(closed_flux=+1) if xs not in {GROUND_X,GLOBAL_FLIP_X}]
    quo=min_x_gap(spectra,qc)
    return {"gauge_gap":gauge_gap,"raw_x_gap":raw,"closed_x_gap":closed,"quotient_x_gap":quo,
            "quotient_matter_gap":min(z,quo),"quotient_full_gap":min(gauge_gap,z,quo)}
def report_model(model,bv=BETA_VALUES):
    hd(model.name); print("  beta    gauge     raw-X    P_X=+X   quot-X    matter(q)  full(q)")
    mins={"gauge_gap":1e9,"raw_x_gap":1e9,"quotient_x_gap":1e9,"quotient_matter_gap":1e9}
    for beta in bv:
        r=diagnostics(model,beta)
        for k in mins: mins[k]=min(mins[k],r[k])
        print(f"  {beta:<7g} {r['gauge_gap']:<9.5g} {r['raw_x_gap']:<8.5g} {r['closed_x_gap']:<8.5g} {r['quotient_x_gap']:<9.5g} {r['quotient_matter_gap']:<10.5g} {r['quotient_full_gap']:<8.5g}")
    print(f"  minima: gauge={mins['gauge_gap']:.5g} raw-X={mins['raw_x_gap']:.5g} quot-X={mins['quotient_x_gap']:.5g} quot-matter={mins['quotient_matter_gap']:.5g}")
    return mins
def report_su3():
    hd("SU(3) character shadow"); print("  Z2 CSS signs vs Z3 SU(3) center: quotient only removes the Z2 global flip, not a real Z3 quotient.")
    print("  beta   lvl dim   gauge    raw-X   P_X=+X  noflip-X matter(noflip)")
    bl=[(0.5,16),(1,16),(2,18),(5,22),(10,26),(20,30),(30,34)]
    mins={"gauge_gap":1e9,"raw_x_gap":1e9,"quotient_x_gap":1e9,"quotient_matter_gap":1e9}
    for beta,level in bl:
        m=SU3Character(level=level); r=diagnostics(m,beta)
        for k in mins: mins[k]=min(mins[k],r[k])
        print(f"  {beta:<6g} {level:<3d} {m.dim:<5d} {r['gauge_gap']:<8.5g} {r['raw_x_gap']:<7.5g} {r['closed_x_gap']:<7.5g} {r['quotient_x_gap']:<8.5g} {r['quotient_matter_gap']:<8.5g}")
    print(f"  minima: gauge={mins['gauge_gap']:.5g} raw-X={mins['raw_x_gap']:.5g} noflip-X={mins['quotient_x_gap']:.5g} noflip-matter={mins['quotient_matter_gap']:.5g}")
    return mins
u1=report_model(U1Rotor()); su2=report_model(SU2Character()); su3=report_su3()
hd("asserts")
assert u1["raw_x_gap"]<1.0; assert su2["raw_x_gap"]<1.0
assert u1["quotient_matter_gap"]>1.9; assert su2["quotient_matter_gap"]>1.9
assert su3["quotient_matter_gap"]<1.0
print("U(1)/SU(2) quotient matter gap > 1.9 ; SU(3) quotient matter gap < 1.0  -> ALL ASSERTS PASSED")
