from fractions import Fraction as F
from itertools import product
CW=list(product([0,1],repeat=8))
# bits: 0G0 1G1 2LQ 3C0 4C1 5I3 6chi 7W
def R1(c): return not(c[0]==1 and c[1]==1)
def R2(c): return c[7]==c[6]
def R3(c):
    col=(c[3],c[4]); return col==(0,0) if c[2]==0 else col!=(0,0)
def R4(c): return not(c[2]==0 and c[5]==0 and c[6]==1)
def valid(c): return R1(c)and R2(c)and R3(c)and R4(c)
def Zf(c): return 1 if c[5]==0 else -1
def sZc(c): return -3 if (c[3],c[4])==(0,0) else -1
def Q(c,Zp=1): return F(1,2)*Zf(c)+F(1,3)*sZc(c)+F(1,2)*Zp
def comp(c): return tuple(1-b for b in c)   # bitwise CPT (sec 2.7)

print("=== Q1: does bitwise-CPT give Q -> -Q ? (antiparticle = complement, Zp=-1) ===")
names={(0,0,0,0,0,0,0,0):'nu_L', (0,0,0,0,0,1,0,0):'e_L',
       (0,0,1,0,1,0,0,0):'u_L(r)', (0,0,1,0,1,1,0,0):'d_L(r)'}
for c,nm in names.items():
    cb=comp(c)
    print(f"  {nm:8s} Q={str(Q(c,1)):4s}  ->  complement Q(Zp=-1)={str(Q(cb,-1)):5s}  "
          f"want -Q={str(-Q(c,1)):4s}  match={Q(cb,-1)==-Q(c,1)}  (complement valid codeword? {valid(cb)})")

ok=bad=0
for c in CW:
    if not valid(c): continue
    if Q(comp(c),-1)==-Q(c,1): ok+=1
    else: bad+=1
print(f"\nover all {ok+bad} valid codewords: Q(c-bar,-1) == -Q(c) holds for {ok}, FAILS for {bad}")
print("reason: the one-hot colour term sumZc in {-3,-1} never maps to -sumZc under bit-complement")
print("        (lepton 00<->11 quark-colour; quark colour<->quark colour) -> NO 3<->3bar conjugation.")

print("\n=== Q2: is I3 the unique bit the weak coin (sec 3.1 zero-controlled CNOT) flips? ===")
# sec 3.1 coin: flip I3 when chi=0. Model: C flips bit5 (I3) iff bit6 (chi)==0.
def coin(c):
    c=list(c)
    if c[6]==0: c[5]^=1
    return tuple(c)
flipped=set()
for c in CW:
    d=coin(c)
    for i in range(8):
        if c[i]!=d[i]: flipped.add(i)
print(f"  bits the coin ever flips: {sorted(flipped)}  -> {'ONLY I3 (bit5)' if flipped=={5} else 'NOT only I3'}")
print("  => T3 target = I3 is structural (the CNOT acts on no other bit). Q2 sound (modulo definitional).")
