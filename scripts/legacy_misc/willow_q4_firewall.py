#!/usr/bin/env python3
"""
Willow Q4 examination + firewall basis check.

Built from canon (ANCHOR S2.1/2.2/2.3/2.5/2.8/2.9/2.14, item 51, item 96).
8-bit register, octant addresses 000..111 -> semantic names (S2.1):
   000 G0   001 G1   010 LQ   011 C0   100 C1   101 I3   110 chi   111 W

History note (2026-05-29): an earlier Part 1 of this file claimed to PROVE an
analytical firewall invariant phi(c)=LQ by finding the codespace splits into
LQ=const components. That was CIRCULAR -- it built its move set by explicitly
forbidding LQ-flips, then "discovered" LQ is constant. WITHDRAWN. Part 1 below
is replaced by the honest test: does the canonical O_h bridge rotation R_d
(S2.5) actually move the LQ face? It does -> phi=LQ is NOT a naive lab-frame
W-invariant, and the firewall rests on S2.14's exhaustive Feshbach enumeration
(H^n_ij = 0, n<=8), NOT on any one-line invariant from this script.
Nothing fitted.
"""
from fractions import Fraction as F
from itertools import product

G0,G1,LQ,C0,C1,I3,CHI,W = range(8)   # bit indices = octant addresses
NAME = {0:"G0",1:"G1",2:"LQ",3:"C0",4:"C1",5:"I3",6:"chi",7:"W"}

def valid(c):
    if c[G0]==1 and c[G1]==1: return False          # R1
    if c[W]!=c[CHI]: return False                    # R2
    if c[LQ]==0 and (c[C0],c[C1])!=(0,0): return False  # R3 lepton
    if c[LQ]==1 and (c[C0],c[C1])==(0,0): return False  # R3 quark
    return True

codespace=[c for c in product((0,1),repeat=8) if valid(c)]
def is_nuR(c): return c[LQ]==0 and c[I3]==0 and c[CHI]==1
def sector(c): return "lepton" if c[LQ]==0 else "quark"
def hd(a,b): return sum(x!=y for x,y in zip(a,b))

L="="*74
def H(s): print("\n"+L+"\n"+s+"\n"+L)

H("0. Codespace census (expect 48 = 12 lepton + 36 quark; 3 nu_R)")
print(f"  |codespace|          = {len(codespace)}")
print(f"  LQ=0 (lepton)        = {sum(1 for c in codespace if c[LQ]==0)}")
print(f"  LQ=1 (quark)         = {sum(1 for c in codespace if c[LQ]==1)}")
print(f"  nu_R pseudocodewords = {sum(1 for c in codespace if is_nuR(c))}")

# ---------------------------------------------------------------------------
H("1. HONEST firewall-basis check: does the canonical R_d move the LQ face?")
# Addresses are cube vertices (b2,b1,b0). O_h bridge rotations (S2.5) permute
# the 8 faces. Two canonical generators:
def trip(a): return ((a>>2)&1,(a>>1)&1,a&1)
def addr(b2,b1,b0): return (b2<<2)|(b1<<1)|b0
def C3(a):                       # 120deg about (1,1,1) body diagonal: cycle the 3 coord bits
    b2,b1,b0=trip(a); return addr(b1,b0,b2)
def C4z(a):                      # 90deg about z (the C4u used in T_y=R T_x R^-1): (x,y)->(y,-x)
    b2,b1,b0=trip(a); bx,by,bz=b0,b1,b2
    return addr(bz,1-bx,by)
print(f"  C3 about (111): LQ(010) -> {NAME[C3(LQ)]}   fixed? {C3(LQ)==LQ}")
print(f"  C4 about z    : LQ(010) -> {NAME[C4z(LQ)]}   fixed? {C4z(LQ)==LQ}")
print("  => The canonical bridge rotations MOVE the LQ face. Therefore the")
print("     lab-frame operator 'read face 010' does NOT commute with R_d, so")
print("     phi(c)=LQ is NOT a naive single-particle W-invariant.")
print("  Caveat: S2.5 proves [G0,W]=0 only for the rotation-FIXED face G0, on a")
print("     single tetrahedral sublattice, reflections open. LQ is the HARDER")
print("     (moved) case; a frame-tracked LQ invariant is neither proven nor")
print("     refuted here and would need the canonical 256-dim walk operator")
print("     (absent from the repo). Honest status: UNDETERMINED.")
print("  => Firewall basis = S2.14 exhaustive enumeration H^n_ij=0 (n<=8 = Q8")
print("     diameter), the framework's own canonical support. Not this script.")

# ---------------------------------------------------------------------------
H("2. Firewall CONSEQUENCE (independent of the invariant question)")
# Whatever the proof basis, the Feshbach single-flip channel count into nu_R is
# directly computable and is what the lepton defect d actually counts.
def feshbach_channels(c):
    return sum(1 for r in codespace if is_nuR(r) and hd(c,r)==1)
g1=lambda c: c[G0]==0 and c[G1]==0
eR  = next(c for c in codespace if c[LQ]==0 and c[I3]==0 and c[CHI]==0 and g1(c))
nuL = next(c for c in codespace if c[LQ]==0 and c[I3]==1 and c[CHI]==0 and g1(c)) \
      if any(c[LQ]==0 and c[I3]==1 and c[CHI]==0 and g1(c) for c in codespace) else None
ups   =[c for c in codespace if c[LQ]==1 and c[I3]==0 and g1(c)]
downs =[c for c in codespace if c[LQ]==1 and c[I3]==1 and g1(c)]
print(f"  lepton single-flip channels into nu_R (e_R rep): {feshbach_channels(eR)}  -> d/N = {F(feshbach_channels(eR),9)}")
print(f"  quark  single-flip channels into nu_R:")
print(f"        up   reps: {sorted(set(feshbach_channels(c) for c in ups))}")
print(f"        down reps: {sorted(set(feshbach_channels(c) for c in downs))}")
print("  => 0 for every quark. The lepton defect-count mechanism is identically")
print("     absent on LQ=1. (Consistency check on S2.14, not a new proof.)")
# nearest nu_R to a quark: what must change
q=ups[0] if ups else downs[0]
nearest=min((c for c in codespace if is_nuR(c)), key=lambda r: hd(q,r))
diff=[NAME[i] for i in range(8) if q[i]!=nearest[i]]
print(f"  sample quark {q} -> nearest nu_R {nearest}; bits to flip: {diff}")
print("  -> requires joint LQ-flip + colour reset = S4.4 GUT escape valve.")

# ---------------------------------------------------------------------------
H("3. WILLOW Q4: are the quark Koide numerators (down,up)=(1,2) derivable?")
print("  (a) nu_R Feshbach index : 0 for quarks (Part 2). Cannot give (1,2). [DEAD]")
print("  (b) S5.9 active-bit (I3): up=0, down=1 -> WRONG ORDER (need down<up). [DEAD]")
print("  (c) item-51 hops        : down=3, up=4 -> integers (3,4)!=(1,2);")
print("                            d=hops-2 unmotivated => integer-shopping. [DEAD]")
print("  See willow_q4_boltzmann.py for the decisive type argument: S5.2")
print("  frustration yields a mass MAGNITUDE exp(phi F/2), not the radian PHASE")
print("  delta. Wrong type for ANY sector's Koide phase.")

H("4. Q4 sub-points (2)-(4)")
print("  (2) delta_d = delta_l/2 : line 648 'structurally unresolved'; item 96")
print("      'single-path scaling' documents, does not derive.")
print("  (3) item 96 uses TWO knobs: up scales DENOM (9->27), down scales NUMER")
print("      (2->1). No single uniform (d,N) map gives all four. ")
print("  (4) chirality: item 51 is RH-specific (d_R,u_R); item 96 chirality-blind.")

H("5. NET")
print("""  Firewall: rests on S2.14 exhaustive enumeration (canonical). The
    analytical phi=LQ invariant is WITHDRAWN as circular; whether a
    frame-tracked LQ is W-conserved is UNDETERMINED here (needs the 256-dim
    operator). The Feshbach channel count is 0 for quarks regardless (Part 2).

  Q4 OPEN, sharply: lepton/nu numerators d=2,3 are named (Feshbach channels);
    quark numerators d=1,2 are produced by NO substrate mechanism (index=0,
    Boltzmann=magnitude-not-phase, step-count=wrong integers). item-96's quark
    d/N_eff borrows the lepton count. The quark Koide phase needs a type-correct
    radian-phase source on LQ=1 (Berry/holonomy of the S5.2 walk, or the S6.7
    CKM irreducible phase) -- recorded in ANCHOR S15 item 135 III(b) revision.""")
