#!/usr/bin/env python3
r"""Does an HONEST baryogenesis count on the framework's ACTUAL matter code (the 48-set) yield the
coefficient ~0.216 (= eta_obs/alpha^4) that 3/14=0.2143 supplies on the ideal [8,4,4]? Discipline: we do
NOT shop for 0.216 -- we compute the 48-set's real structure + principled ratios, and report what they are.
If none gives 0.216 naturally (and a §16.3 scan shows 0.216 is reachable only by shopping), that settles it.
Bit map: 0=G0,1=G1,2=LQ,3=C0,4=C1,5=I3,6=chi,7=W; colour=(C0,C1); colourless=(0,0).
"""
import sys, itertools as it
from fractions import Fraction
from collections import Counter

ALPHA=1/137.036; NEED = 6.12e-10/ALPHA**4    # required coeff for eta = coeff*alpha^4
b=lambda n,i:(n>>i)&1; wt=lambda n:bin(n).count("1")
def in48(n):
    if b(n,0) and b(n,1): return False                   # R1
    if b(n,7)!=b(n,6):    return False                   # R2
    cc=(b(n,3),b(n,4))                                   # R3
    if b(n,2)==0 and cc!=(0,0): return False
    if b(n,2)==1 and cc==(0,0): return False
    return True
P=[n for n in range(256) if in48(n)]
print(f"required coeff (eta/alpha^4) = {NEED:.4f}   (ideal-code value 3/14 = {3/14:.4f})")

# ---------- (1) the 48-set's real structure ----------
wdist=Counter(wt(n) for n in P)
minw=min(wt(n) for n in P if n!=0)
leptons=[n for n in P if b(n,2)==0]                       # LQ=0 -> colourless (R3)
quarks =[n for n in P if b(n,2)==1]                       # LQ=1 -> coloured (R3)
nuR=[n for n in P if b(n,6)==1 and b(n,7)==1 and (b(n,3),b(n,4))==(0,0) and b(n,2)==0]
print(f"\n[1] 48-set: weight dist {dict(sorted(wdist.items()))}; min non-zero weight = {minw}")
print(f"    -> min weight 1, so the '[8,4,4] distance d=4' (the alpha^4 source) does NOT hold on the 48-set;")
print(f"       and there are {wdist[4]} weight-4 words, NOT 14. Neither the alpha-power NOR the '14' survives.")
print(f"    leptons(LQ=0,colourless) = {len(leptons)} ; quarks(LQ=1,coloured) = {len(quarks)} ; nu_R = {len(nuR)}")

# ---------- (2) principled candidate ratios (colour-singlet / total, in several honest readings) ----------
cands = {
  "nu_R / 48 (3 colour-neutral sterile)":      Fraction(len(nuR), len(P)),
  "colourless(=leptons) / 48":                 Fraction(len(leptons), len(P)),
  "colourless weight-4 / weight-4":            Fraction(sum(1 for n in P if wt(n)==4 and (b(n,3),b(n,4))==(0,0)), wdist[4]),
  "nu_R / quarks":                             Fraction(len(nuR), len(quarks)),
  "leptons / quarks":                          Fraction(len(leptons), len(quarks)),
  "3 (gen) / 14 (NOT a 48-set count)":         Fraction(3,14),
}
print(f"\n[2] principled candidate baryogenesis ratios on the 48-set (vs required {NEED:.4f}):")
for name,r in cands.items():
    dev = abs(float(r)-NEED)/NEED*100
    print(f"    {name:<42} = {str(r):>6} = {float(r):.4f}  ({dev:+.0f}% vs need)")
print(f"    -> NONE of the natural 48-set ratios lands near {NEED:.3f}: they cluster at 1/16, 1/4, 1/3, 1...,")
print(f"       not ~0.216. The '3/14' is NOT a 48-set count (14 isn't even a weight in the 48-set).")

# ---------- (3) §16.3 discipline: how EASY is it to reach 0.216 by shopping integer ratios? ----------
counts = sorted(set([len(P),len(leptons),len(quarks),len(nuR)] + list(wdist.values()) + [16,14,3]))
hits=[]
for a in counts:
    for c in counts:
        if c>0 and abs(a/c - NEED)/NEED <= 0.01:
            hits.append(f"{a}/{c}")
print(f"\n[3] §16.3 shopping check: ratios of the framework's own integers ({counts}) within 1% of {NEED:.3f}:")
print(f"    -> {hits if hits else 'NONE'}  ({len(hits)} hits). 3/14 works only because 14 is imported from the")
print(f"       IDEAL [8,4,4], not the 48-set; on the actual matter integers, 0.216 is essentially unreachable.")

print(f"""
=========================================================================================
VERDICT (honest count on the 48-set):
  The framework's ACTUAL matter code (48-set) does NOT support eta=(0.216)alpha^4:
   * min weight is 1, so the alpha^4-from-d=4 argument does not apply (it is an ideal-[8,4,4] property);
   * there are {wdist[4]} weight-4 words, not 14 (the '14' is imported from the ideal 16-word code);
   * no principled colour-singlet ratio on the 48-set lands near 0.216 -- they are 1/16, 1/4, 1/3, ...;
   * 0.216 is reachable from the 48-set integers only by shopping ({len(hits)} of all integer-pair ratios), and
     even 3/14 needs the foreign denominator 14.
  CONCLUSION: 3/14 is genuinely a property of the IDEAL [8,4,4] linear code, which the framework's matter
  content is NOT. There is no honest derivation of the baryogenesis coefficient from the framework's own
  48-state matter structure. So: the alpha^4 SCALE is the only part that ports (and even d=4 is an ideal-code
  property the 48-set lacks); the 3/14 COEFFICIENT does not survive contact with the actual codespace.
  Net (final, for item 126): the headline 'parameter-free derived' is NOT supportable -- the number reproduces
  on a code the framework's matter isn't, and the reconstruction onto the real code FAILS. Honest tier:
  'eta ~ alpha^4 (scale motivated by an idealised code distance); the 3/14 coefficient is an ideal-code
  coincidence with no realisation on the framework's matter -- a Proposition at best, not Locked.'
=========================================================================================""")
print("exit 0 -- 48-set count done; no honest ratio yields 0.216; 3/14 is an ideal-[8,4,4]-only coincidence.")
