"""Extended Target A anomaly-compression scan.

Checks that role permutations are SM-isomorphic, and that R1 generation
choices and R2 weak-label choices are anomaly-neutral once the matter/colour
sector is fixed.
"""

from fractions import Fraction as F
from itertools import product, permutations, combinations

CW = list(product([0,1], repeat=8))

# A role map says which bit index plays each physical role.
# SM: I3=5, chi=6, colour=(3,4), LQ=2, gen=(0,1), W=7
def make(rm):
    I3,chi,LQ = rm['I3'], rm['chi'], rm['LQ']
    Ca,Cb = rm['C']; G0,G1 = rm['gen']; W = rm['W']
    def R1(c): return not (c[G0]==1 and c[G1]==1)
    def R2(c): return c[W]==c[chi]
    def R3(c):
        col=(c[Ca],c[Cb])
        return col==(0,0) if c[LQ]==0 else col!=(0,0)
    def R4(c): return not (c[LQ]==0 and c[I3]==0 and c[chi]==1)
    def Zf(c): return 1 if c[I3]==0 else -1
    def sZc(c): return -3 if (c[Ca],c[Cb])==(0,0) else -1
    def Q(c,b=F(1,3)): return F(1,2)*Zf(c) + b*sZc(c) + F(1,2)
    def T3(c): return F(1,2)*Zf(c) if c[chi]==0 else F(0)
    def Y(c,b=F(1,3)): return Q(c,b)-T3(c)
    def Ye(c,b=F(1,3)):
        y=Y(c,b); return y if c[chi]==0 else -y
    trip=lambda c:(c[Ca],c[Cb])!=(0,0); doub=lambda c:c[chi]==0
    return R1,R2,R3,R4,Q,Ye,trip,doub

def six(content,Ye,trip,doub):
    A1=sum(Ye(c) for c in content if trip(c))
    A2=F(1,2)*sum(Ye(c) for c in content if doub(c))
    A3=sum(Ye(c)**3 for c in content)
    A4=sum(Ye(c) for c in content)
    A5=sum(1 for c in content if trip(c) and not doub(c)==False and c) # placeholder
    return A1,A2,A3,A4
def afree(content,Ye,trip,doub):
    return six(content,Ye,trip,doub)==(0,0,0,0)

SMrm=dict(I3=5,chi=6,LQ=2,C=(3,4),gen=(0,1),W=7)
R1,R2,R3,R4,Q,Ye,trip,doub=make(SMrm)
SM=[c for c in CW if R1(c)and R2(c)and R3(c)and R4(c)]
print("=== sanity (SM roles) ===")
print("size",len(SM),"A1..A4",six(SM,Ye,trip,doub),"sumQ2",sum(Q(c)**2 for c in SM))

print("\n=== ROLE-PERMUTATION genericity (which bits = I3,chi,LQ,colour) ===")
bits=list(range(8)); tot=0; af=0; af_list=[]
for I3,chi,LQ in permutations(bits,3):
    rest=[b for b in bits if b not in (I3,chi,LQ)]
    for C in combinations(rest,2):
        rem=[b for b in rest if b not in C]
        gen=(rem[0],rem[1]); W=rem[2]            # spectators (anomaly-neutral)
        rm=dict(I3=I3,chi=chi,LQ=LQ,C=C,gen=gen,W=W)
        r1,r2,r3,r4,q,ye,tr,db=make(rm)
        cont=[c for c in CW if r1(c)and r2(c)and r3(c)and r4(c)]
        if not cont: continue
        tot+=1
        if afree(cont,ye,tr,db):
            af+=1; af_list.append(rm)
print("role assignments tested:",tot,"  anomaly-free:",af,f"({100*af/tot:.1f}%)")
# how many of the anomaly-free are SM-isomorphic (same charge multiset incl chirality)?
def sig(content,q,ye,db,tr):
    return tuple(sorted((str(q(c)),str(ye(c)),db(c),tr(c)) for c in content))
sm_sig=sig(SM,Q,Ye,doub,trip)
iso=0
for rm in af_list:
    r1,r2,r3,r4,q,ye,tr,db=make(rm)
    cont=[c for c in CW if r1(c)and r2(c)and r3(c)and r4(c)]
    if sig(cont,q,ye,db,tr)==sm_sig: iso+=1
print("  of those, SM-isomorphic (same (Q,Yeff,chir,colour) multiset):",iso,
      "| structurally-different anomaly-free:",af-iso)

print("\n=== R1 variants (all Boolean fns of the generation pair) ===")
n_af=0;n_tot=0
for f in product([0,1],repeat=4):   # f over (G0,G1) in {00,01,10,11}
    def r1f(c,f=f): return f[2*c[0]+c[1]]==1
    cont=[c for c in CW if r1f(c) and R2(c) and R3(c) and R4(c)]
    if not cont: continue
    n_tot+=1; 
    if afree(cont,Ye,trip,doub): n_af+=1
print(f"non-empty generation rules: {n_tot}, anomaly-free: {n_af}  -> R1/generations anomaly-NEUTRAL" if n_af==n_tot else f"{n_af}/{n_tot}")

print("\n=== R2 variants (relation of W to other bits) ===")
import itertools
rels={'W=chi(SM)':lambda c:c[7]==c[6],'W=~chi':lambda c:c[7]!=c[6],
      'W=I3':lambda c:c[7]==c[5],'W=0':lambda c:c[7]==0,'no R2':lambda c:True,
      'W=LQ':lambda c:c[7]==c[2]}
for name,r2f in rels.items():
    cont=[c for c in CW if R1(c) and r2f(c) and R3(c) and R4(c)]
    print(f"  {name:10s} size={len(cont):3d}  anomaly-free={afree(cont,Ye,trip,doub)}")
