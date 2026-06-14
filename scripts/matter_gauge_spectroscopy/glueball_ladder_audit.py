#!/usr/bin/env python3
r"""Hadronic audit: the glueball ladder m_N = (2N-1)Lambda (ANCHOR §7.10/7.13/7.15, "Capstone").
The headline (L1534): "5 out of 5 LQCD glueball masses to 1-3% ... the ONLY continuous parameter
being the empirical Lambda_QCD." Audited on its own terms, against the canon's OWN quoted LQCD values
and its OWN status caveats. Three findings, each asserted:
  (1) SELF-CONTRADICTION: the 2 parity-odd channels carry a SECOND continuous fit parameter delta~0.155
      (L1520 "empirical fit parameter"; L1532 "rigorous derivation anchored as downstream target") --
      so "only Lambda" is false by the canon's own admission. Honest count: 3 at Lambda-only + 2 fitted.
  (2) INVERTED EPISTEMICS on the parity-even three: the boasted "strongest match" (1+- at 0.3%, L1393/
      L1431) is the ONE state whose N=5 is POST-HOC (L1446: N=5 rests on "topological intuition PLUS the
      0.2% LQCD match" -- circular). The two states whose N is FORCED by O_h geometry (0++ N=3 = 3 C4 axes;
      2++ N=4 = 4 C3 hexagons) are BOTH ~2.9% low at the canonical Lambda=332.
  (3) GENUINE CONTENT + SCALE TENSION: the (2N-1) ladder and the forced N-assignments ARE defensible; but
      the glueballs prefer Lambda~336.5 (all 3 to ~1.6%), NOT the nucleon-anchored 332 -- a ~3% inter-sector
      chiral-scale spread the single-Lambda framing absorbs silently.
LQCD values are quoted FROM CANON (auditing internal consistency); they themselves vary by reference
(Morningstar-Peardon vs Athenodorou-Teper) -- noted, not relied on for the structural points.
"""
import sys, math
from scipy.optimize import brentq

LAM = 332.0  # canonical chiral scale (ANCHOR §1.4), nucleon-anchored: M_N/(2 sqrt2) = 938.9/2.828 = 332.0
assert abs((938.272+939.565)/2/(2*math.sqrt(2)) - LAM) < 0.5, "Lambda=332 is the nucleon-anchored value"

# ---- the three parity-even states (ANCHOR L1437-1441, L1420-1421, L1393) ----
# (J^PC, N, N forced by O_h?, LQCD central, LQCD err, why-N)
PE = [
    ("0++", 3, True,  1710, 50, "3 orthogonal C4 face-axes (A_1g) -- forced by O_h"),
    ("2++", 4, True,  2390, 70, "4 C3 body-diagonal Petrie hexagons (T_2g) -- forced by O_h"),
    ("1+-", 5, False, 2980, 35, "N=5: 'minimal connected solenoid'; obstruction UNPROVEN (L1446)"),
]
print("(1) PARITY-EVEN three: m = (2N-1)Lambda at canonical Lambda=332, and the per-state effective Lambda:")
print(f"   {'JPC':>4} {'N':>2} {'forced?':>8} {'(2N-1)L':>9} {'LQCD':>10} {'dev':>7} {'L_eff=m/(2N-1)':>15}")
Leff = {}
for jpc, N, forced, m, err, why in PE:
    pred = (2*N-1)*LAM
    dev = (pred-m)/m*100
    Leff[jpc] = m/(2*N-1)
    print(f"   {jpc:>4} {N:>2} {str(forced):>8} {pred:>9.0f} {m:>6}±{err:<3} {dev:>+6.1f}% {Leff[jpc]:>14.1f}")
# the two FORCED states agree on L_eff; the POST-HOC one differs by ~3%
assert abs(Leff["0++"]-Leff["2++"]) < 1.0, "the two forced-N states agree on effective Lambda"
spread = (Leff["0++"]-Leff["1+-"])/Leff["1+-"]*100
print(f"   -> the two FORCED-N states agree: L_eff(0++)={Leff['0++']:.1f} ~ L_eff(2++)={Leff['2++']:.1f}.")
print(f"      the POST-HOC-N state wants L_eff(1+-)={Leff['1+-']:.1f} -- a {spread:+.1f}% disagreement.")
assert spread > 2.5, "post-hoc 1+- wants a >2.5% different Lambda than the forced pair"
print(f"   => canonical Lambda=332 NAILS the post-hoc-N 1+- (dev {((9*LAM-2980)/2980*100):+.1f}%) and leaves the")
print(f"      two FORCED states ~2.9% low. The canon boasts the 1+- as 'the strongest single match in the")
print(f"      entire glueball ladder' (L1393/L1431) -- i.e. the LEAST-constrained N gives the best number.")
print(f"      L1446 concedes N=5 rests on 'topological intuition PLUS the 0.2% LQCD match' -- the match is")
print(f"      cited as evidence for the N that produces the match. Circular.")

# ---- a single best-fit Lambda for all three: min-max fractional deviation ----
def dev_of(L, N, m): return (2*N-1)*L/m - 1
# binding extremes: 0++ (most negative) vs 1+- (most positive); equalise |dev|
f = lambda L: dev_of(L,3,1710) + dev_of(L,5,2980)   # 0++ low + 1+- high = 0 at the minimax
Lstar = brentq(f, 300, 360)
maxdev = max(abs(dev_of(Lstar,N,m)) for _,N,_,m,_,_ in PE)*100
print(f"\n   SINGLE best-fit Lambda for the 3 parity-even: Lambda* = {Lstar:.1f} MeV -> all three within {maxdev:.2f}%.")
print(f"   So 'glueballs to 1-3%' is achievable, but at Lambda*={Lstar:.1f}, NOT the nucleon-anchored 332.")
assert 335 < Lstar < 338 and maxdev < 2.0, "glueballs prefer Lambda~336.5, fitting all 3 to <2%"

# ---- inter-sector scale tension: each hadronic sector's preferred Lambda ----
print(f"\n(3) INTER-SECTOR chiral-scale spread (one Lambda is used everywhere; each sector's best-fit differs):")
L_nuc = (938.272+939.565)/2/(2*math.sqrt(2))      # nucleon M_N=2sqrt2 L
L_pi  = math.sqrt(4*math.pi)*92.28                 # f_pi = L/sqrt(4pi), PDG f_pi=92.28(10) MeV
L_glue = Lstar
for lbl, L in [("nucleon  M_N=2sqrt2 L", L_nuc), ("pion     L=sqrt(4pi) f_pi", L_pi),
               ("glueball best-fit (2N-1)L", L_glue)]:
    print(f"   {lbl:>28}: Lambda = {L:.1f} MeV")
span = (max(L_glue,L_nuc,L_pi)-min(L_glue,L_nuc,L_pi))/min(L_glue,L_nuc,L_pi)*100
print(f"   -> spread {min(L_pi,L_nuc,L_glue):.0f}-{max(L_pi,L_nuc,L_glue):.0f} MeV = {span:.1f}%. The single canonical 332 is a")
print(f"      nucleon-anchored COMPROMISE; per-sector it is individually off by up to ~{span:.0f}%. So the")
print(f"      hadronic masses match 'to a few %' BECAUSE they share one scale that is ~{span:.0f}% off per sector.")
assert span > 2.5, "inter-sector Lambda spread exceeds ~2.5%"

# ---- (2) the hidden SECOND continuous parameter delta in the 'Capstone' master formula ----
# master (L1527): m = (2N-1)L + N_shared*(1/3 - delta) L ; delta~0.155 fit to the 2 parity-odd channels.
# canon's beta values (L1520): 0.71, 0.54  <=>  N_shared in {4,3} times (1/3-delta):
print(f"\n(2) the 'Capstone' (L1527) for the 2 PARITY-ODD channels adds a SECOND continuous parameter delta:")
delta = 0.155
for Nsh, beta_claimed in [(4, 0.71), (3, 0.54)]:
    beta = Nsh*(1/3 - delta)
    print(f"   N_shared={Nsh}: beta = N_shared*(1/3 - delta) = {beta:.3f}  vs canon's quoted beta={beta_claimed} (L1520)")
    assert abs(beta - beta_claimed) < 0.02, "one universal delta~0.155 reproduces both quoted beta"
print(f"   -> beta(0.71),beta(0.54) ARE one universal delta=0.155 x geometric N_shared in {{3,4}} (good: 1 param,")
print(f"      not 2). BUT delta is STILL a continuous fit parameter: L1520 calls beta 'an empirical fit")
print(f"      parameter rather than a first-principles result'; L1532 anchors delta's 'rigorous derivation")
print(f"      as downstream target'. So L1534's 'the ONLY continuous parameter being Lambda' is FALSE by the")
print(f"      canon's OWN caveats: the 2 parity-odd of the '5 of 5' ride on a fitted delta.")

print(f"""
=========================================================================================
VERDICT (glueball ladder m_N=(2N-1)Lambda, "Capstone 5 of 5 to 1-3%, only-Lambda"):
  GENUINE CORE (credit where due): the (2N-1) ladder is real topological structure, and the N-assignments
  for the two FORCED states -- 0++ <-> N=3 (three C4 face-axes), 2++ <-> N=4 (four C3 Petrie hexagons) --
  ARE O_h geometry, not fits. The 7/5 tensor/scalar ratio is good (canon ALREADY demotes it to a §16.3
  class-1 consistency check, 10 competitors, L1274). This is structure, not numerology.
  BUT the "5 of 5 to 1-3% with only Lambda" headline (L1534) overstates three ways, each from canon's OWN text:
    (1) HIDDEN 2nd PARAMETER: the 2 parity-odd channels carry a fitted delta~0.155 (L1520 "empirical fit
        parameter"; L1532 derivation "downstream target"). Honest: 3 at Lambda-only + 2 with a fitted delta.
    (2) INVERTED EPISTEMICS: the boasted 0.3% "strongest match" (1+-, L1393/L1431) is the POST-HOC-N=5 state
        whose N is fixed BY that match (circular, L1446); the two FORCED-N states are BOTH ~2.9% low at 332.
    (3) SILENT SCALE SPREAD: the glueballs prefer Lambda~336.5 (all 3 to ~1.6%), the nucleon 332, the pion
        ~327 -- a ~3% inter-sector spread absorbed silently by the single canonical 332.
  NET: same headline-vs-footnote pattern as the cosmology sector, but MILDER -- here the core (ladder +
  forced O_h N) is genuine; the overstatement is in the precision claim ("1-3%, only-Lambda, strongest
  match"), not in the existence of structure. Honest calibration: "a topological glueball ladder whose
  forced ratio matches at ~3%, with a fitted delta for the parity-odd pair and a ~3% inter-sector scale
  spread" -- NOT "5 of 5 to 1-3% from one parameter."
=========================================================================================""")
print("exit 0 -- glueball ladder: genuine forced ladder + ratio; but hidden delta, post-hoc-N boast, and ~3% scale spread overstate the '5/5 only-Lambda' headline.")
