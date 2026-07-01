#!/usr/bin/env python3
"""
Willow Q4, attempt via the canonical S5.2 Boltzmann/frustration mechanism.
Question: does S5.2 frustration counting produce the quark Koide NUMERATORS
d=1 (down), d=2 (up) with a NAMED referent -- or not?

Canon used verbatim:
  S2.1 register: addr 000..111 -> G0,G1,LQ,C0,C1,I3,chi,W
  S5.2 F(c) = sum over Q3 cube edges (i,j) of (c_i XOR c_j)   [INTEGER 0..12]
        mass MAGNITUDE M = exp(phi F /2).  (a magnitude, NOT a phase)
  S5.8 Koide delta enters as RADIANS inside cos(delta + 2*pi*n/3).
Nothing fitted.
"""
from itertools import product
from fractions import Fraction as F

# address (b2,b1,b0) integer 0..7 -> semantic name (S2.1)
NAME = {0:"G0",1:"G1",2:"LQ",3:"C0",4:"C1",5:"I3",6:"chi",7:"W"}
idx  = {v:k for k,v in NAME.items()}
CUBE_EDGES = [(u,v) for u in range(8) for v in range(u+1,8)
              if bin(u^v).count("1")==1]   # 12 edges of Q3

def cw(**kw):
    a=[0]*8
    for k,v in kw.items(): a[idx[k]]=v
    return tuple(a)

def valid(a):
    G0,G1,LQ,C0,C1,I3,chi,Wb = (a[idx[n]] for n in
        ["G0","G1","LQ","C0","C1","I3","chi","W"])
    if G0==1 and G1==1: return False        # R1
    if Wb!=chi: return False                # R2
    if LQ==0 and (C0,C1)!=(0,0): return False   # R3
    if LQ==1 and (C0,C1)==(0,0): return False   # R3
    return True

def frust(a):
    return sum(1 for (u,v) in CUBE_EDGES if a[u]!=a[v])

codespace=[a for a in product((0,1),repeat=8) if valid(a)]

L="="*72
def Hd(s): print("\n"+L+"\n"+s+"\n"+L)

Hd("0. census")
print(f"  |codespace| = {len(codespace)}  (expect 48)")

# Gen-1 quark codewords: LQ=1, G0=G1=0. up: I3=0 ; down: I3=1.
g1q=[a for a in codespace if a[idx['LQ']]==1 and a[idx['G0']]==0 and a[idx['G1']]==0]
ups   =[a for a in g1q if a[idx['I3']]==0]
downs =[a for a in g1q if a[idx['I3']]==1]
def show(tag,lst):
    fs=sorted(frust(a) for a in lst)
    print(f"  {tag:14s} n={len(lst)}  F values={fs}  meanF={sum(fs)/len(fs):.3f}")
Hd("1. Frustration F on gen-1 quark codewords (colour x chirality)")
show("up   (I3=0)",ups)
show("down (I3=1)",downs)
# split by chirality chi
for tag,lst in (("up",ups),("down",downs)):
    for ch in (0,1):
        sub=[a for a in lst if a[idx['chi']]==ch]
        if sub: show(f"{tag} chi={ch}",sub)

Hd("2. THE CATEGORY CHECK (decisive)")
print("  S5.2 output is a MAGNITUDE exponent: M = exp(phi*F/2), F an integer.")
print("  S5.8 Koide delta is a PHASE in RADIANS: cos(delta + 2*pi*n/3).")
print("  A magnitude mechanism cannot emit a phase. So S5.2 frustration")
print("  CANNOT source the Koide delta for ANY sector -- including quarks.")
print("  (Consistent with S5.8 line 646: S5.2 supplies the quark MASS")
print("   HIERARCHY/magnitudes; the Koide circulant supplies PHASES.)")

Hd("3. Does any F-based integer accidentally equal the numerators (down,up)=(1,2)?")
# candidate named counts on a quark codeword:
def n_active_internal(a):  # count of logical-1 among internal EW bits {I3,chi,W} (S5.9 'active')
    return sum(a[idx[n]] for n in ["I3","chi","W"])
def n_colour_active(a):    # active colour bits
    return a[idx['C0']]+a[idx['C1']]
def dF_from_lepton(a):     # F(quark) - F(corresponding charged lepton baseline)
    return None
repU, repD = ups[0], downs[0]
print(f"  sample up   {repU}  F={frust(repU)} active_int={n_active_internal(repU)} colour={n_colour_active(repU)}")
print(f"  sample down {repD}  F={frust(repD)} active_int={n_active_internal(repD)} colour={n_colour_active(repD)}")
print("  Target numerators: down=1, up=2.")
print("  - active_int: depends on chi (gen/chirality), not a clean (1,2) split by up/down.")
print("  - colour_active: same multiset {1,2} for BOTH up and down (R,G have 1; B has 2).")
print("    => cannot distinguish up from down. DEAD as the (down,up)=(1,2) source.")
print("  - any F-difference is a MAGNITUDE exponent, wrong type for a phase (S2).")

Hd("4. VERDICT on Willow Q4 (attempt complete)")
print("""  S5.2 Boltzmann frustration does NOT supply the quark Koide numerators:
    (i) CATEGORY: it yields mass magnitudes exp(phi F/2), not phases; the
        Koide delta is a radian phase. Wrong type. (decisive)
    (ii) COUNTS: no F-derived integer cleanly splits up/down as (2,1) with a
        named referent; colour-activity gives the same {1,2} multiset to both.
  => The quark Koide PHASES delta_d=1/9, delta_u=2/27 have NO substrate
     mechanism at present:
       - nu_R Feshbach index  -> 0 for quarks (Colour Firewall; phi=LQ invariant)
       - S5.2 Boltzmann       -> magnitudes, not phases (this script)
     item-96's 'd/N_eff' for quarks BORROWS the lepton defect notation for a
     quantity neither quark mechanism produces. The equality delta_d=delta_l/2,
     delta_u=delta_l/N_c is an ASSERTED analogy, not a derivation.

  This SHARPENS (does not close) the open problem:
     The quark Koide phase requires a genuine PHASE source on the LQ=1 sector
     (candidates: a Berry/holonomy phase of the S5.2 dwell-time walk, or the
     CKM-side irreducible phase of S6.7 -- both produce radian phases). The
     magnitude mechanism is the wrong generator; the nu_R index is firewalled.
""")
