"""Target A anomaly-compression audit.

Reproduces the exact-rational checks used by the anomaly-compression paper:
the 45-codeword SM content, six anomaly sums, colour-weight uniqueness,
R4 chiral nontriviality, and the non-genericity scan over matter/colour
sector rules.
"""

from fractions import Fraction as F
from itertools import product

# Codeword = (G0,G1,LQ,C0,C1,I3,chi,W)  (per ANCHOR 2.1 octant addresses)
CW = list(product([0,1], repeat=8))
def R1(c): return not (c[0]==1 and c[1]==1)        # 3 generations
def R2(c): return c[7]==c[6]                        # W = chi
def R3(c):
    LQ,C0,C1 = c[2],c[3],c[4]
    return (C0,C1)==(0,0) if LQ==0 else (C0,C1)!=(0,0)
def R4(c): return not (c[2]==0 and c[5]==0 and c[6]==1)   # forbid nu_R

def Zf(c):  return 1 if c[5]==0 else -1             # Pauli-Z of I3 (0->+1)
def sumZc(c): return -3 if (c[3],c[4])==(0,0) else -1   # one-hot: singlet -3, triplet -1
def Q(c,a=F(1,2),b=F(1,3),Zp=1): return a*Zf(c) + b*sumZc(c) + F(1,2)*Zp
def T3(c):  return F(1,2)*Zf(c) if c[6]==0 else F(0)    # doublet (chi=0) only
def Y(c,a=F(1,2),b=F(1,3)): return Q(c,a,b) - T3(c)
def Yeff(c,a=F(1,2),b=F(1,3)):                      # all-left convention
    y = Y(c,a,b); return y if c[6]==0 else -y
def triplet(c): return (c[3],c[4])!=(0,0)          # colour-charged
def doublet(c): return c[6]==0

def six(content,a=F(1,2),b=F(1,3)):
    A1 = sum(Yeff(c,a,b)               for c in content if triplet(c))      # [SU3]^2 U1
    A2 = F(1,2)*sum(Yeff(c,a,b)        for c in content if doublet(c))      # [SU2]^2 U1
    A3 = sum(Yeff(c,a,b)**3            for c in content)                    # [U1]^3
    A4 = sum(Yeff(c,a,b)               for c in content)                    # grav-U1
    A5 = sum(1 for c in content if triplet(c) and c[6]==0) - \
         sum(1 for c in content if triplet(c) and c[6]==1)                  # [SU3]^3
    A6 = sum(1 for c in content if doublet(c)) // 2 % 2                     # Witten
    return (A1,A2,A3,A4,A5,A6)

SM = [c for c in CW if R1(c) and R2(c) and R3(c) and R4(c)]
print("=== SANITY GATE (SM content, sec 2.8 charges) ===")
print("content size :", len(SM), "(expect 45)")
print("six anomalies:", six(SM), "(expect all 0)")
print("sum Q  =", sum(Q(c) for c in SM), " sum Q^2 =", sum(Q(c)**2 for c in SM),
      " (expect 0, 16)")
print("sum QL^3 =", sum(Q(c)**3 for c in SM if c[6]==0),
      " sum QR^3 =", sum(Q(c)**3 for c in SM if c[6]==1), " (expect -2/9, -2/9)")

print("\n=== TEST 1: colour weight b scanned, is 1/3 forced? ===")
hits=[]
for num in range(-6,13):
    for den in range(1,13):
        b=F(num,den)
        if six(SM,b=b)[:4]==(0,0,0,0): hits.append(b)
hits=sorted(set(hits))
print("colour weights b in [-6,1] step rationals q<=12 giving A1..A4 = 0:", hits)
print("  -> uniquely b = 1/3" if hits==[F(1,3)] else "  -> NOT unique")

print("\n=== TEST 2: is the chiral content special? (R4 removes nu_R) ===")
V48=[c for c in CW if R1(c) and R2(c) and R3(c)]     # no R4 -> includes nu_R
print("R1&R2&R3 (48, vector-like, nu_R present):", six(V48), "<- cancels trivially (vector-like)")
print("R1&R2&R3&R4 (45, chiral, nu_R removed)  :", six(SM), "<- still cancels = the non-trivial fact")

print("\n=== TEST 3: genericity over R3-type sector rules ===")
# family: allow any nonempty set V of (LQ,C0,C1) triples; keep R1,R2; remove the
# 'up-type colour-singlet chi=1' state per generation (the R4 analogue) to stay chiral.
triples=list(product([0,1],repeat=3))
free=0; anomaly_free=0; chiral_af=0
import itertools
for r in range(1,len(triples)+1):
    for V in itertools.combinations(triples,r):
        Vs=set(V)
        cont=[c for c in CW if R1(c) and R2(c) and (c[2],c[3],c[4]) in Vs
              and not (c[2]==0 and c[5]==0 and c[6]==1)]
        if not cont: continue
        free+=1
        s=six(cont)
        if s==(0,0,0,0,0,0):
            anomaly_free+=1
            # chiral = not vector-like (some state lacks its chi-partner)
            partners=set((c[0],c[1],c[2],c[3],c[4],c[5],1-c[6]) for c in cont)
            if not set(cont)<=partners or any((c[0],c[1],c[2],c[3],c[4],c[5],1-c[6]) not in set(cont) for c in cont):
                chiral_af+=1
print("sector rules tested:", free)
print("fully anomaly-free  :", anomaly_free, f"({100*anomaly_free/free:.1f}%)")
print("  of which chiral    :", chiral_af)
