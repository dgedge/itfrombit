#!/usr/bin/env python3
r"""§16.3 + structural audit of ANCHOR §15 item 126: baryogenesis eta = (3/14)alpha^4 + (1/3)alpha^5.
Headline: "parameter-free" eta=6.145e-10 vs Planck 6.12+-0.04e-10, ~0.6 sigma. Derivation chain:
  alpha^4  <- [8,4,4] code distance d=4 (weight-4 logical fault, base per-bit leak rate alpha, item 79)
  /14      <- 14 weight-4 codewords of [8,4,4] (weight enumerator W(x)=1+14x^4+x^8)
  x3       <- "3 SU(3) colour-singlet configs (the R1,R2,R3 stabilizers) out of the 14"
  +(1/3)a^5<- weight-4 fault + 1-bit gauge fluctuation, isotropic /3 (added same day, improves the fit)
This script checks each piece that is cleanly checkable and flags the rest.
"""
import sys, itertools as it
from fractions import Fraction
import numpy as np

ALPHA = 1/137.036
ETA_OBS, ETA_ERR = 6.12e-10, 0.04e-10

# ---------- (1) verify the numbers ----------
a4 = ALPHA**4; a5 = ALPHA**5
eta_lo = (3/14)*a4
eta_full = (3/14)*a4 + (1/3)*a5
print(f"[1] alpha={ALPHA:.6f}  alpha^4={a4:.4e}  alpha^5={a5:.4e}")
print(f"    (3/14)alpha^4              = {eta_lo:.4e}  (resid {(eta_lo-ETA_OBS)/ETA_OBS*100:+.2f}%, "
      f"{abs(eta_lo-ETA_OBS)/ETA_ERR:.2f} sigma)")
print(f"    (3/14)alpha^4+(1/3)alpha^5 = {eta_full:.4e}  (resid {(eta_full-ETA_OBS)/ETA_OBS*100:+.2f}%, "
      f"{abs(eta_full-ETA_OBS)/ETA_ERR:.2f} sigma)")
assert abs(eta_lo-6.076e-10)<2e-13 and abs(eta_full-6.145e-10)<3e-13, "numbers must match the canon"

# ---------- (2) the [8,4,4] weight enumerator: is "14" real? ----------
G = np.array([[1,0,0,0,0,1,1,1],[0,1,0,0,1,0,1,1],[0,0,1,0,1,1,0,1],[0,0,0,1,1,1,1,0]])  # [8,4,4] gen
cws = []
for bits in it.product([0,1],repeat=4):
    cws.append(np.array(bits)@G % 2)
weights = sorted(int(c.sum()) for c in cws)
from collections import Counter
wc = Counter(weights)
print(f"\n[2] [8,4,4] extended Hamming: {len(cws)} codewords; weight enumerator = {dict(sorted(wc.items()))}")
print(f"    -> W(x) = 1 + {wc[4]}x^4 + 1x^8 : the '14 weight-4 codewords' IS a real code property. [OK]")
assert wc[0]==1 and wc[4]==14 and wc[8]==1, "[8,4,4] weight enumerator must be 1+14x^4+x^8"

# ---------- (3) BUT the framework's actual codeword set is the 48 (R1/R2/R3), NOT the 16-word [8,4,4] ----------
def is_cw48(n):  # the canonical 48-codeword rules (R1/R2/R3), from keff_construct
    b=lambda i:(n>>i)&1
    if b(0) and b(1): return False                       # R1
    if b(7)!=b(6):    return False                       # R2  (W=chi)
    cc=(b(3),b(4))                                        # R3
    if b(2)==0 and cc!=(0,0): return False
    if b(2)==1 and cc==(0,0): return False
    return True
P48=[n for n in range(256) if is_cw48(n)]
import math
is_pow2 = (len(P48) & (len(P48)-1))==0
print(f"\n[3] framework's physical codeword set: |P| = {len(P48)}  (a linear [8,k] code has 2^k codewords)")
print(f"    {len(P48)} is a power of two? {is_pow2}  ->  the 48-set is NOT a linear code; it is NOT the")
print(f"    16-word [8,4,4]. So the derivation mixes TWO different objects: the '14' is a property of the")
print(f"    *ideal* [8,4,4] (16 codewords), while 'R1,R2,R3 = 3 of the 14' refers to the *48-set* rules.")
print(f"    (48 = 3 x 16 suggests '3 generations of [8,4,4]', but then the weight-4 count per cell is still")
print(f"    14 of 16, and the 'colour-symmetric = 3' identification with R1/R2/R3 is asserted, not shown.)")
assert not is_pow2

# ---------- (4) §16.3: is the 3/14 coefficient special, or one of many simple rationals at alpha^4? ----------
need = ETA_OBS/a4                                        # required prefactor at alpha^4
tol = abs(eta_lo-ETA_OBS)/ETA_OBS                        # the ACTUAL leading-order match tolerance (0.7%)
comps=[]
for q in range(1,30):
    for p in range(1,q):
        r=Fraction(p,q)
        if abs(float(r)-need)/need <= tol:
            comps.append(r)
comps=sorted(set(comps), key=float)
print(f"\n[4] §16.3 competitor count at alpha^4: required prefactor eta_obs/alpha^4 = {need:.5f}")
print(f"    simple rationals p/q (q<30) within the leading-order tolerance {tol*100:.2f}%: {len(comps)}")
print(f"    -> {', '.join(str(c)+f'={float(c):.4f}' for c in comps[:8])}")
print(f"    SELF-CORRECTION: this is LOW ({len(comps)}), not high -- 3/14 is fairly UNIQUE at the 0.7% tolerance,")
print(f"    so (unlike the gravity alpha^2's 60 competitors) the coefficient is NOT numerology-exposed.")
print(f"    Caveat: q=14 sits just above the e1 simplicity ceiling (12), so it is borderline-simple, not tiny.")
# is alpha^4 itself forced? which power gives a clean O(1) prefactor?
print(f"    power check (eta_obs/alpha^n): a^3->{ETA_OBS/ALPHA**3:.4f}  a^4->{ETA_OBS/a4:.4f}  a^5->{ETA_OBS/a5:.2f}")
print(f"    -> only alpha^4 gives an O(1) clean prefactor (~0.22); so the MAGNITUDE itself selects n=4,")
print(f"       which COINCIDES with d=4 (a genuine motivation, not a reverse-solve -- credit the alpha^4).")

# ---------- (5) the alpha^5 term: post-hoc fit-improver ----------
print(f"\n[5] the (1/3)alpha^5 term: it shifts {eta_lo:.3e} -> {eta_full:.3e}, i.e. residual {abs(eta_lo-ETA_OBS)/ETA_OBS*100:.1f}% "
      f"({abs(eta_lo-ETA_OBS)/ETA_ERR:.2f}sigma) -> {abs(eta_full-ETA_OBS)/ETA_OBS*100:.1f}% ({abs(eta_full-ETA_OBS)/ETA_ERR:.2f}sigma).")
print(f"    It was 'added 2026-05-27' (the SAME day as the leading term) and IMPROVES the agreement -- the")
print(f"    hallmark of a post-hoc correction. Its 1/3 ('isotropic across 3 axes') is itself a fresh O(1)")
print(f"    coefficient choice; the leading (3/14)alpha^4 already matched at 0.7% (1 sigma).")

print(f"""
=========================================================================================
VERDICT (item 126 baryogenesis eta = (3/14)alpha^4 [+ (1/3)alpha^5]):
  GROUNDED: the alpha^4 power is the strongest part -- it is BOTH a genuine QEC code-distance argument
  (weight-d=4 logical fault, per-bit leak rate alpha) AND the unique power giving an O(1) prefactor.
  This is materially better-motivated than item 136's alpha^4 (which was a magnitude-selected Dirac exponent).
  NOT numerology (credit): unlike the gravity alpha^2, the 3/14 prefactor is fairly UNIQUE at the match
  tolerance ({len(comps)} competitor at 0.7%) -- it was NOT shopped from many options.
  STRUCTURAL ISSUES (where it actually weakens):
   * Code conflation: '14' is the weight-4 count of the IDEAL 16-word [8,4,4] LINEAR code, but the framework's
     actual codeword set is the 48-word (NON-linear, 48 != 2^k) R1/R2/R3 structure. The 48=3x16 reading
     ('3 generations of [8,4,4]') is plausible but unstated; either way '3 colour-symmetric weight-4 codewords
     = R1,R2,R3' is ASSERTED, never exhibited -- the load-bearing '3' is the one step not shown.
   * The (1/3)alpha^5 is a same-day post-hoc fit-improver (1.1 sigma -> 0.6 sigma) with its own fresh 1/3; the
     leading (3/14)alpha^4 already sat at 1.1 sigma, so the headline 0.6 sigma is dressed by the correction.
  NET: BETTER than the gravity routes -- the alpha^4 is genuinely motivated (d=4 AND clean-prefactor) and the
  3/14 is not numerology -- but it is NOT the clean 'parameter-free derived prediction' the headline implies:
  the load-bearing '3 of 14' is unshown and conflates a 16-word linear code with the 48-word non-linear set,
  and the sub-sigma match leans on a post-hoc alpha^5. Honest tier: the leading-order (3/14)alpha^4 = 6.08e-10
  at ~1.1 sigma is a code-combinatoric PROPOSITION (alpha^4 derived, 3/14 plausible-but-unexhibited); call it
  'consistent with, pending the explicit 3-colour-singlet weight-4 count on the framework's own codeword set.'
=========================================================================================""")
print("exit 0 -- baryogenesis numbers verified; [8,4,4] weight-enum real; 48!=16 conflation + §16.3 density flagged.")
