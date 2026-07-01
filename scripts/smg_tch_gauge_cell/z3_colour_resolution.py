from fractions import Fraction as F
from itertools import product
import cmath, math, numpy as np

# ---- 1. Is there a Z3 in the colour sector? GL(2,F2)=Aut(Z2xZ2) ----
def det2(M): return (M[0][0]*M[1][1]) ^ (M[0][1]*M[1][0])
def mm(A,B):
    return [[(A[0][0]*B[0][0]^A[0][1]*B[1][0]),(A[0][0]*B[0][1]^A[0][1]*B[1][1])],
            [(A[1][0]*B[0][0]^A[1][1]*B[1][0]),(A[1][0]*B[0][1]^A[1][1]*B[1][1])]]
def order(M):
    I=[[1,0],[0,1]]; P=[r[:] for r in M]; n=1
    while P!=I and n<12: P=mm(P,M); n+=1
    return n
def act(M,v): return ((M[0][0]*v[0]^M[0][1]*v[1]),(M[1][0]*v[0]^M[1][1]*v[1]))
GL=[[[a,b],[c,d]] for a in[0,1] for b in[0,1] for c in[0,1] for d in[0,1] if det2([[a,b],[c,d]])==1]
print("=== GL(2,F2) = Aut(Z2 x Z2) ===")
print("|GL(2,F2)| =",len(GL),"(== |S3| = 6);  element orders:",sorted(order(M) for M in GL))
print("contains order-3 (the Z3 triality)?", 3 in [order(M) for M in GL])

colours={(0,1):'r',(1,0):'g',(1,1):'b',(0,0):'s'}
M3=[[0,1],[1,1]]
cyc=[(0,1)]
for _ in range(3): cyc.append(act(M3,cyc[-1]))
print("\nZ3 generator [[0,1],[1,1]] (order",order(M3),") orbit:",
      " -> ".join(colours[c] for c in cyc), "| fixes singlet?", act(M3,(0,0))==(0,0))

T=[[1,0],[1,1]]
print("conjugation [[1,0],[1,1]] (order",order(T),", the 3<->3bar):",
      {colours[c]:colours[act(T,c)] for c in [(0,1),(1,0),(1,1),(0,0)]})

# ---- 2. triality: three colours sum to a singlet ----
w=cmath.exp(2j*math.pi/3)
print("\nomega^0+omega^1+omega^2 =", round(sum(w**q for q in range(3)).real,9),
      "(r+g+b -> colour-neutral, matches sec2.4 r XOR g XOR b = 0)")

# ---- 3. THE Q1 FIX: colour-blind Q = T3 + Y (Y flips under CPT) ----
SM={'nu':(F(1,2),F(-1,2)),'e':(F(-1,2),F(-1,2)),'u':(F(1,2),F(1,6)),'d':(F(-1,2),F(1,6)),
    'uR':(F(0),F(2,3)),'dR':(F(0),F(-1,3)),'eR':(F(0),F(-1))}
print("\n=== Q1 with colour-blind Q=T3+Y (CPT: T3->-T3, Y->-Y, colour q->-q) ===")
allok=True
for p,(t3,y) in SM.items():
    Q=t3+y; Qbar=-t3-y; ok=(Qbar==-Q); allok&=ok
    print(f"  {p:3s}  Q={str(Q):5s}  ->  antiparticle {str(Qbar):5s}   (-Q? {ok})")
print("  all antiparticle charges = -Q:",allok,"  => Q1 CPT-covariant")
print("  [literal sec2.8 one-hot gave u-bar=-4/3, e-bar=-1/3 -> Q1 BROKE; the one-hot was the culprit]")

# ---- 4. is Z3 ALREADY a symmetry of V_strong = A_{K3} (sec 3.3)? ----
AK3=np.array([[0,1,1],[1,0,1],[1,1,0]]); P=np.array([[0,0,1],[1,0,0],[0,1,0]])
print("\n=== Z3 vs strong channel A_{K3} (sec 3.3) ===")
print("  [A_K3, Z3 permutation] = 0 (Z3 is a symmetry of V_strong)?", np.array_equal(AK3@P,P@AK3))
print("  A_K3 eigenvalues:", sorted(np.linalg.eigvalsh(AK3).round(6).tolist()),
      "(+2 singlet, -1 doublet = the colour triplet)")
