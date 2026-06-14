#!/usr/bin/env python3
r"""The load-bearing baryogenesis step: are EXACTLY 3 of the 14 weight-4 codewords 'colour-singlet'?
Item 126 asserts '3 colour-symmetric weight-4 codewords = R1,R2,R3' out of 14, giving eta=(3/14)alpha^4.
The framework never defines 'colour-singlet' on the codewords, so we (a) build both the ideal 16-word
[8,4,4] and the framework's actual 48-word set, (b) try several natural definitions, (c) see whether 3
emerges robustly or is choice-dependent. Bit map (ANCHOR §2.1): 0=G0,1=G1,2=LQ,3=C0,4=C1,5=I3,6=chi,7=W;
colour = (C0,C1) = bits (3,4); 3 colours = (C0,C1) in {01,10,11}, (00)=colourless.
"""
import sys, itertools as it
import numpy as np
from collections import Counter

def wt(n): return bin(n).count("1")
def b(n,i): return (n>>i)&1
def colour(n): return (b(n,3), b(n,4))   # (C0,C1)

# ---------- the ideal [8,4,4] extended Hamming (standard generator), 14 weight-4 codewords ----------
G = np.array([[1,0,0,0,0,1,1,1],[0,1,0,0,1,0,1,1],[0,0,1,0,1,1,0,1],[0,0,0,1,1,1,1,0]])
cw16 = []
for bits in it.product([0,1],repeat=4):
    v = (np.array(bits)@G)%2
    cw16.append(sum(int(v[i])<<i for i in range(8)))     # pack to int (bit i = position i)
w4 = [n for n in cw16 if wt(n)==4]
print(f"[ideal [8,4,4]] {len(cw16)} codewords; weight-4 count = {len(w4)} (expect 14)")
print(f"  the 14 weight-4 codewords (binary, lo->hi bit) and their colour (C0,C1):")
for n in sorted(w4):
    print(f"    {format(n,'08b')[::-1]}   (C0,C1)={colour(n)}  {'colourless' if colour(n)==(0,0) else 'coloured'}")
col_dist = Counter(colour(n) for n in w4)
print(f"  colour distribution of the 14: {dict(col_dist)}")

# ---------- try several natural 'colour-singlet' definitions ----------
defs = {
  "colourless (C0=C1=0)":            lambda n: colour(n)==(0,0),
  "C0=C1 (swap-symmetric)":          lambda n: b(n,3)==b(n,4),
  "carries a definite colour (C0,C1)!=00": lambda n: colour(n)!=(0,0),
  "LQ=0 & colourless (R3 lepton-like)":    lambda n: b(n,2)==0 and colour(n)==(0,0),
  "W=chi (R2) among weight-4":       lambda n: b(n,7)==b(n,6),
}
print(f"\n[counts among the 14 weight-4 codewords, by definition]:")
hit3 = []
for name,f in defs.items():
    c = sum(1 for n in w4 if f(n))
    flag = "  <-- equals 3!" if c==3 else ""
    if c==3: hit3.append(name)
    print(f"    {name:<42} -> {c}{flag}")

# ---------- the framework's ACTUAL 48-word set: does it even have 14 weight-4? ----------
def is_cw48(n):
    if b(n,0) and b(n,1): return False                   # R1
    if b(n,7)!=b(n,6):    return False                   # R2
    cc=(b(n,3),b(n,4))                                   # R3
    if b(n,2)==0 and cc!=(0,0): return False
    if b(n,2)==1 and cc==(0,0): return False
    return True
P48=[n for n in range(256) if is_cw48(n)]
wdist48=Counter(wt(n) for n in P48)
print(f"\n[framework 48-set] weight distribution: {dict(sorted(wdist48.items()))}")
print(f"  weight-4 count in the 48-set = {wdist48[4]}  (the '14' is NOT a property of the actual codeword set)")

# ---------- can R1,R2,R3 (the rules) BE 3 weight-4 colour-singlet codewords? ----------
# R2 as a linear check W xor chi=0 -> its support is {chi,W}=bits 6,7 -> weight 2, not 4.
# R1 (G0.G1!=1) and R3 (LQ/colour conditional) are NON-linear -> not codewords/parity rows at all.
print(f"""
[R1,R2,R3 as 'weight-4 codewords']: R2 (W=chi) is the only linear rule; its parity support {{chi,W}} is
  weight-2, NOT weight-4. R1 (G0.G1!=1) and R3 (LQ<->colour conditional) are NON-LINEAR exclusions, not
  parity-check rows, so they are not codewords of ANY linear code. So 'R1,R2,R3 = 3 weight-4 codewords by
  inspection of the parity-check matrix' has no realisation in the actual rule set.""")

lq_of_colourless = [(format(n,'08b')[::-1], b(n,2)) for n in w4 if colour(n)==(0,0)]
print(f"""
=========================================================================================
VERDICT (does '3 colour-singlet of 14 weight-4' hold?) -- BALANCED, both halves:
  THE NUMBER REPRODUCES (credit): under the natural physical reading 'colour-singlet = colour-neutral
  (C0,C1)=(0,0)', EXACTLY 3 of the 14 weight-4 codewords are colourless -> 3/14 has a genuine code-count
  realisation, NOT arbitrary numerology. (Colour distribution {{(0,0):3,(0,1):4,(1,0):4,(1,1):3}}.)
  BUT the framework's stated DERIVATION fails on every specific:
   * Wrong mechanism: 'the 3 = R1,R2,R3 weight-4 codewords' has NO realisation -- R2 (W=chi) is weight-2,
     R1 (G0.G1!=1) and R3 are NON-linear (not codewords of any linear code). The number 3 is right; the
     stated reasoning is not.
   * Wrong physics: a colour-NEUTRAL codeword is LEPTON-like (no net colour), not a baryon (a baryon is the
     antisymmetric 3-COLOURED qqq singlet -- a combination of coloured codewords, not one colourless word).
     And the 3 colourless weight-4 words have LQ-bits {dict(lq_of_colourless) if False else [lq for _,lq in lq_of_colourless]} -- a mix, some with LQ=1
     (quark-like + colourless = R3-VIOLATING, not even in the 48-set). So they are not '3 baryons'.
   * Definition-dependent: colourless->3, but swap-symmetric(C0=C1)->6, coloured->11. The 3-giving reading
     was the favourable one.
   * Wrong code: the 14 is the IDEAL 16-word [8,4,4]; the framework's actual matter set is the 48-word
     NON-linear code, which has {wdist48[4]} weight-4 words, not 14.
  CONCLUSION: the count, attempted, gives 3 colourless of 14 weight-4 -- so 3/14 is a REAL code coincidence,
  not pure numerology -- but it is NOT the '3 colour-singlet baryons via R1,R2,R3' the framework claims: the
  stated mechanism is false, the physical identification (colourless=baryon) is backwards, and it borrows the
  ideal [8,4,4] not the 48-set. NET tier: the alpha^4 SCALE is derived; the 3/14 COEFFICIENT is a numerically-
  reproducible-but-mis-derived code coincidence -- 'consistent with eta via a real but mis-attributed count',
  not a clean parameter-free derivation. (Better than the gravity routes; short of an exhibited derivation.)
=========================================================================================""")
print("exit 0 -- count attempted: 3 colourless of 14 weight-4 DOES hold, but NOT via R1,R2,R3 and not on the 48-set.")
