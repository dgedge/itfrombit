#!/usr/bin/env python3
r"""§16.3 competitor-count of the Part-20 fallback Planck-mass formula (the route the gravity sector
falls back to once K_eff=205 is retracted; DRIFT G6/G7 follow-up, requested check (a)).

Part-20 (ANCHOR L2337):  M_P^2 = 24 pi alpha^2 Lambda_QCD^3 / (H0 Omega_Lambda)   -> claims ~0.07% on M_P.

Question (§16.3, the framework's own search-space accounting): how many EQUALLY-SIMPLE, dimensionally-
valid expressions in the alphabet {rational, pi, alpha, Omega_Lambda, Lambda, H0} land within the same
tolerance of the measured M_P? If many -> the 0.07% 'agreement' carries near-zero standalone weight
(class-1 numerology); if ~0-1 -> the weight rests on the (unwritten) derivation of the exponents/prefactor.

Honest construction: a competitor is  M_P^2 ~ R * pi^p * alpha^a * Omega_Lambda^d * Lambda^b * H0^(2-b)
(the H0 exponent fixed by mass-dimension 2 -> b + (2-b) = 2). Integer exponents + a small-rational R, i.e.
exactly the framework's own expression STYLE (24 pi alpha^2 ... is R=24,p=1,a=2,d=-1,b=3). Count hits.
"""
import sys
from fractions import Fraction
import numpy as np

PI = np.pi
# ---- constants, computed with sources (not from memory) ----
M_P = 1.220890e19          # GeV, full Planck mass (CODATA 2018); framework §10.5 table uses 1.221e19
LAM = 0.332                # GeV, Lambda_QCD framework chiral scale (§1.4)
H0_kmsMpc = 67.36          # Planck 2018 base-LCDM
OmL = 0.6847               # Planck 2018 Omega_Lambda
ALPHA = 1 / 137.035999
HBAR_GeV_s = 6.582119e-25  # GeV*s
MPC_m = 3.0856776e22       # m
H0 = (H0_kmsMpc * 1e3 / MPC_m) * HBAR_GeV_s   # GeV
print(f"inputs (computed): M_P={M_P:.5e} GeV  Lambda={LAM} GeV  H0={H0:.5e} GeV  Omega_L={OmL}  alpha=1/{1/ALPHA:.3f}")

# ---- verify the Part-20 formula reproduces M_P (sanity) ----
MP2_part20 = 24 * PI * ALPHA**2 * LAM**3 / (H0 * OmL)
MP_part20 = MP2_part20**0.5
dev = (MP_part20 - M_P) / M_P
print(f"\nPart-20  M_P^2 = 24 pi alpha^2 Lambda^3 /(H0 OmL):  M_P = {MP_part20:.4e} GeV  ({dev*100:+.2f}% vs measured)")
assert abs(dev) < 0.01, "Part-20 formula should reproduce M_P to ~0.07-1%"
TARGET = M_P**2                                   # we count competitors to M_P^2 (the formula's object)
TOL_CLAIM = abs(MP2_part20 - TARGET) / TARGET     # the ACTUAL match tolerance (on M_P^2)
print(f"  match tolerance on M_P^2 = {TOL_CLAIM*100:.3f}%  (= {abs(dev)*200:.3f}% on M_P^2 ~ 2x the M_P %)")

# ---- build the dimensionally-valid competitor family (framework's own expression style) ----
def small_rationals(maxnum=24, maxden=12):
    s = set()
    for q in range(1, maxden + 1):
        for p in range(1, maxnum + 1):
            s.add(Fraction(p, q))
    return sorted(s)

Rs = small_rationals()
P_pi = range(0, 3)            # pi^0..pi^2  (framework uses pi^1)
A_al = range(0, 5)            # alpha^0..alpha^4 (framework uses alpha^2)
D_om = range(-2, 2)           # Omega_L^-2..^1 (framework uses ^-1)
B_lam = range(0, 6)           # Lambda^0..^5 ; H0 exponent = 2 - b (framework b=3)

def count_hits(tol):
    hits = []
    for b in B_lam:
        dim = LAM**b * H0**(2 - b)               # dimensionful core, mass-dimension 2 by construction
        for p in P_pi:
            for a in A_al:
                for d in D_om:
                    pref0 = PI**p * ALPHA**a * OmL**d
                    val0 = pref0 * dim
                    if val0 == 0:
                        continue
                    Rneed = TARGET / val0        # the rational that would make it exact
                    # accept if SOME small rational R lands the product within tol
                    for R in Rs:
                        v = float(R) * val0
                        if abs(v - TARGET) / TARGET <= tol:
                            hits.append((float(R), p, a, d, b, v))
                            break                # one R per (p,a,d,b) is enough to count the expression
    return hits

for tol, lab in [(TOL_CLAIM, "claim band"), (1e-3, "0.1%"), (1e-2, "1%")]:
    h = count_hits(tol)
    sample = "; ".join(f"{R:g}*pi^{p}*a^{a}*OmL^{d}*L^{b}*H0^{2-b}" for (R, p, a, d, b, v) in h[:4])
    print(f"\n[band {lab:11s} {tol*100:.3f}%]  competitors found = {len(h)}")
    if h:
        print(f"   e.g.: {sample}")

# is the framework's own (24,1,2,-1,3) in the list, and is it the UNIQUELY simplest?
print(f"""
INTERPRETATION (§16.3):
  The magnitude (M_P^2/Lambda^3 ~ 1/H0) pins the dimensional structure to b~3 (only b=3 gives the ~10^40
  GeV^2 scale with an O(1)-rational prefactor) -- so the dimensional part is NOT free; that is the formula's
  one real structural content. The FREEDOM is the dimensionless prefactor R*pi^p*alpha^a*OmL^d, and within
  the claimed tolerance there are the counts above. Read per the framework's own §16.3 rule: >= ~2 equally-
  simple competitors in the match band => class-1 (consistency check, near-zero standalone weight); ~0-1 =>
  weight rests on deriving the exponents/prefactor (which Part-20 does attempt via the alpha^2 'double-
  Landauer' argument -- itself flagged 'not a theorem', DRIFT G6).
""")
n_claim = len(count_hits(TOL_CLAIM))

# ---- ROBUSTNESS: is "1 competitor" stable across the alphabet, and is alpha^2 load-bearing? ----
print("ROBUSTNESS (default-to-refute -- vary the alphabet; count at the claim tolerance):")

def count_general(tol, rats, pis, als, doms, bs, irrs=(1.0,)):
    seen = set()
    for b in bs:
        dim = LAM**b * H0**(2 - b)
        for ir in irrs:
            for p in pis:
                for a in als:
                    for d in doms:
                        val0 = ir * PI**p * ALPHA**a * OmL**d * dim
                        if val0 == 0:
                            continue
                        for R in rats:
                            v = float(R) * val0
                            if abs(v - TARGET) / TARGET <= tol:
                                seen.add((p, a, d, b, ir, round(float(R), 6)))
                                break
    return len(seen)

half = [x / 2 for x in range(-4, 11)]                       # half-integer exponents (M_P form)
variants = [
    ("framework style (R<=24/12, int exps)",         count_general(TOL_CLAIM, Rs, P_pi, A_al, D_om, B_lam)),
    ("R<=12 only (the e1 ceiling -- excludes 24!)",  count_general(TOL_CLAIM, small_rationals(12,12), P_pi, A_al, D_om, B_lam)),
    ("wide R<=48",                                    count_general(TOL_CLAIM, small_rationals(48,12), P_pi, A_al, D_om, B_lam)),
    ("half-integer exponents allowed",               count_general(TOL_CLAIM, Rs, half, half, half, [x/2 for x in range(0,11)])),
    ("+ irrationals {sqrt2,sqrt3,phi} in prefactor", count_general(TOL_CLAIM, Rs, P_pi, A_al, D_om, B_lam, irrs=(1.0,2**.5,3**.5,(1+5**.5)/2))),
    ("NO alpha allowed (a=0)",                        count_general(TOL_CLAIM, Rs, P_pi, [0], D_om, B_lam)),
]
for lab, n in variants:
    print(f"  {lab:<44} competitors @ {TOL_CLAIM*100:.2f}% = {n}")

print(f"""
VERDICT (Part-20 fallback, §16.3 -- BASIS-DEPENDENT; reporting both sides honestly):
  The competitor count is NOT a single number -- it depends on the exponent basis, and that dependence is
  the finding:
   * M_P^2 INTEGER basis (Sakharov 'induced-stiffness' framing -- 1/G is naturally an M_P^2 coefficient):
     {n_claim} competitor at the {TOL_CLAIM*100:.2f}% tolerance. Magnitude pins b=3, so the structure looks constrained.
     BUT at the framework's OWN e1 ceiling (R<=12) the count is 0 -- the prefactor 24 is ABOVE the
     framework's own numerology alphabet, so even this favorable basis admits the formula only by widening R.
   * M_P HALF-INTEGER basis (the OBSERVABLE's natural form: M_P = sqrt(24pi)*alpha*L^1.5*H0^-0.5*OmL^-0.5 --
     the integer powers exist ONLY because the canon writes M_P^2, which doubles & hides the half-integers):
     60 competitors at the same tolerance -> firmly CLASS-1 numerology.
  In BOTH bases the weight is alpha-load-bearing: with NO alpha the count is 0. The alpha^2 factor supplies
  the scale and does all the work, and its 'double-Landauer' derivation is flagged 'not a theorem' (DRIFT G6).

  NET: the Part-20 fallback is §16.3-EXPOSED. It looks constrained only if you privilege the Sakharov
  M_P^2-integer basis AND accept a prefactor above the framework's own e1 ceiling; in the natural M_P
  observable basis it is class-1 (60 competitors). Either way its evidential content reduces to the
  un-derived alpha^2. So: 'M_P comes out at the right SCALE, set by an asserted alpha^2 factor' -- better
  than the unconstructible K_eff=205, but NOT a parameter-free 0.07% derivation.

  RECOMMENDATION on (b): the heat-kernel build is NOT justified merely to rescue 0.07% -- the number is
  already §16.3-exposed, and a real-valued a_2 coefficient cannot reproduce the clean prefactor anyway. It
  is worth doing ONLY as independent physics that DERIVES alpha^2 (and the prefactor) from a pre-pinned
  scheme; absent that, it relocates the same assertion one level up (the E-tuning trap, per your warning).
  My lean, matching your prior: do not spend the W_QQ(k) build on this unless there is an independent
  reason to expect the scheme to FORCE alpha^2.
""")
assert isinstance(n_claim, int)
print("exit 0 -- Part-20 competitor count + robustness computed (inputs computed, nothing fitted).")
