#!/usr/bin/env python3
"""
Feasibility + faithfulness check for the 2-cell (32-mode) SO(10)/SM SMG test, BEFORE
committing the heavy ED. Two questions:
  (1) does naive SM-invariant 2-cell hopping give a chiral mode + mirror, or a TRIVIAL
      band insulator (in which case the interaction just competes with a trivial gap -- the
      SO(6) defect, no chiral mode to 'spare')?
  (2) how large is the SM-Cartan-neutral half-filling sector (the Lanczos target)?
"""
import numpy as np, itertools

G0,G1,LQ,C0,C1,I3,CHI,W = range(8)
def R1(c): return not(c[G0]==1 and c[G1]==1)
def R2(c): return c[W]==c[CHI]
def R3(c): return (c[C0],c[C1])==(0,0) if c[LQ]==0 else (c[C0],c[C1])!=(0,0)
gen=sorted([c for c in itertools.product([0,1],repeat=8)
            if (c[G0],c[G1])==(0,0) and R1(c) and R2(c) and R3(c)])
assert len(gen)==16

# integer SM weight per physical mode: (colour wx, wy ; 2*T3L ; 6*Y)
CW={ -1:(0,0) }  # placeholder
COLW={0:(1,0),1:(-1,1),2:(0,-1)}     # SU(3) fundamental weights (sum 0); antitriplet = negatives
CIDX={(0,1):0,(1,0):1,(1,1):2}
def weights(c):
    left=(c[CHI]==0); sgn=1 if left else -1
    if c[LQ]==0: cw=(0,0)
    else:
        w=COLW[CIDX[(c[C0],c[C1])]]; cw=w if left else (-w[0],-w[1])
    t3l = (1 if c[I3]==0 else -1) if left else 0        # 2*T3L (left doublets only)
    BL=(-1.0 if c[LQ]==0 else 1.0/3.0)*sgn
    T3R=0.0 if left else (-0.5 if c[I3]==0 else 0.5)
    Y=T3R+BL/2.0
    return (cw[0],cw[1],t3l,int(round(6*Y)))

phys=[weights(c) for c in gen]
mirror=[(-a,-b,-d,-e) for (a,b,d,e) in phys]   # mirror = conjugate 16-bar (negated charges)
modes=phys+mirror
print(f"32 modes built (16 physical + 16 mirror=conjugate).")

# (1) naive single-particle hopping spectrum: H = -t (|A,i><B,i| + h.c.), mode-diagonal
t=1.0
sp=np.zeros((32,32))
for i in range(16):
    sp[i,16+i]=sp[16+i,i]=-t
ev=np.sort(np.linalg.eigvalsh(sp))
print(f"\n(1) naive 2-cell single-particle spectrum (t=1): min={ev.min():.2f} max={ev.max():.2f}")
print(f"    distinct levels: {sorted(set(np.round(ev,3)))}  -> bonding/antibonding +-t, each x16")
print(f"    zero modes (|E|<1e-9): {int(np.sum(np.abs(ev)<1e-9))}")
print("    => half-filled = TRIVIAL band insulator (gap 2t), NO zero mode / NO chiral fermion.")
print("    A faithful chiral test needs a lattice-Dirac/Wilson kinetic term with a zero mode +")
print("    doubler -- a single naive hopping has none. Naive 2-cell repeats the SO(6) defect.")

# (2) SM-Cartan-neutral, half-filling (N=16) sector size via DP over modes
from collections import defaultdict
dp=defaultdict(int); dp[(0,0,0,0,0)]=1     # (N, cx,cy, 2T3L, 6Y)
for (a,b,d,e) in modes:
    nd=defaultdict(int)
    for (n,cx,cy,w,y),cnt in dp.items():
        nd[(n,cx,cy,w,y)]+=cnt                       # mode empty
        nd[(n+1,cx+a,cy+b,w+d,y+e)]+=cnt             # mode filled
    dp=nd
neutral_half=dp[(16,0,0,0,0)]
total_half=sum(v for (n,*r),v in dp.items() if n==16)
print(f"\n(2) half-filling (N=16) sector:  C(32,16) = {total_half:,}")
print(f"    SM-Cartan-neutral (cx=cy=2T3L=6Y=0) sub-sector: {neutral_half:,}")
print(f"    (Lanczos target; the true SM-SINGLET sector is smaller but >= this needs non-abelian projection.)")
