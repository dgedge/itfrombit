#!/usr/bin/env python3
r"""§16.3 + constructibility audit of ANCHOR §15 item 136 (the THIRD 'parameter-free Planck mass' route,
distinct from the 205 route and the Part-20 alpha^2 route; requested 2026-06-08).

Item 136 unified equation (ANCHOR L3102):
    M_Pl = M_0 * C * (L_H/a_0)^(1/4) * (1/alpha)^4 ,   M_0 ~ 2*Lambda_QCD ,  C = sqrt(2pi) .
Claims: ~3% match (1.26e19 vs 1.22e19), 'zero free parameters'; the ~1e8 residual in G 'matches' 1/alpha^4.
The two exponents are asserted independent: 4 = [8,4,4] code distance d (alpha^-4), 1/4 = conformal scaling
dim for 1D->4D coarse-graining. The audit tests whether that pair is DERIVED or magnitude-selected.

KEY structural fact (mass dimension): (L_H/a_0) and alpha are DIMENSIONLESS, so M_Pl = [one mass] x
[dimensionless]. The mass is forced to be ~Lambda (the only intrinsic scale); ALL predictive content sits
in the dimensionless (L_H/a_0)^p * alpha^-q. The magnitude is then ONE equation in (p,q) -> a LINE, not a point.
"""
import sys, math
from fractions import Fraction
import numpy as np

# ---- constants, computed with sources ----
M_P   = 1.220890e19            # GeV, full Planck mass (CODATA)
LAM   = 0.332                  # GeV, Lambda_QCD (framework chiral scale)
ALPHA = 1/137.035999
H0    = (67.36e3/3.0856776e22)*6.582119e-25     # GeV (Planck-2018 H0)
L_H   = 1.0/H0                 # GeV^-1, Hubble length c/H0
a0    = 1.0/LAM                # GeV^-1, natural lattice step (body-diagonal ~ 1/Lambda, sect 1.4)
RATIO = L_H/a0                 # = L_H*Lambda, the holographic ratio (item uses ~2.2e41)
M0    = 0.646                  # GeV, bare cluster mass ~ 2*Lambda (item 130)
C     = math.sqrt(2*math.pi)   # sqrt(2pi) CLT factor
print(f"inputs: M_P={M_P:.4e}  Lambda={LAM}  1/alpha={1/ALPHA:.3f}  L_H/a0={RATIO:.3e}  M0={M0} (~2Lambda={2*LAM})  C=sqrt(2pi)={C:.4f}")

# ---- T0: reproduce the item's number ----
M_136 = M0 * C * RATIO**0.25 * (1/ALPHA)**4
dev = (M_136 - M_P)/M_P
print(f"\n[T0] item-136 formula M0*sqrt(2pi)*(L_H/a0)^(1/4)*alpha^-4 = {M_136:.4e} GeV  ({dev*100:+.2f}% vs M_P)")
assert abs(dev) < 0.06, "should reproduce M_P to a few %"

# logs used throughout (base 10)
Lr = math.log10(RATIO)         # ~41.36
La = math.log10(1/ALPHA)       # ~2.137
def logpred(R, c, p, q, massGeV=LAM):
    return math.log10(R) + math.log10(c) + math.log10(massGeV) + p*Lr + q*La
target = math.log10(M_P)

# ---- T1: the (p,q) underdetermination -- magnitude is ONE line, item 117 sits at a DIFFERENT point ----
# required line (canonical prefactor R=2 i.e. M0=2Lambda, c=sqrt(2pi)):  p*Lr + q*La = target - log10(2*Lambda*C)
const = target - math.log10(2*LAM*C)
print(f"\n[T1] magnitude is ONE equation: {Lr:.3f}*p + {La:.3f}*q = {const:.3f}  (a LINE in (p,q), not a point).")
print(f"     required p for each integer alpha-power q (canonical prefactor 2*Lambda*sqrt(2pi)):")
for q in range(2, 7):
    p_req = (const - q*La)/Lr
    # nearest simple fraction
    best = min((Fraction(n, d) for d in range(1, 9) for n in range(1, d)), key=lambda fr: abs(float(fr)-p_req))
    err = abs(logpred(2, C, p_req, q) - target)
    # how far is the SIMPLE fraction from giving 3%?
    devp = (10**(logpred(2, C, float(best), q)) - M_P)/M_P
    flag = "  <-- ITEM 136 (1/4, 4)" if q==4 else ""
    print(f"       q={q}: p_required={p_req:.4f}  nearest simple={best} ({float(best):.4f}) -> M_P dev {devp*100:+6.1f}%{flag}")
print("  => item 117 (published, Zenodo) used (alpha^2, ratio^1); item 136 uses (alpha^4, ratio^1/4) -- the")
print("     SAME framework, TWO different exponent pairs, both called 'parameter-free'. The text ADMITS it:")
print("     'The alpha^4 replaces the alpha^2'. So the exponents are demonstrably ADJUSTABLE (internal >=2).")

# ---- T2: hypersensitivity -- the 'match' is carried by the EXACT rationals, not by a computed magnitude ----
dlnM_dp = math.log(RATIO)      # d ln M_P / dp  (natural log)
dp_for_3pct = 0.03/dlnM_dp
print(f"\n[T2] sensitivity d(ln M_P)/dp = ln(L_H/a0) = {dlnM_dp:.1f}.  To stay within 3%, p must be 1/4 +/- {dp_for_3pct:.4f}.")
print(f"     i.e. the exponent must be 0.2500 to ~4 sig figs. The prefactor M0*sqrt(2pi) is O(1): changing")
print(f"     M0 from 2*Lambda to Lambda (x2) shifts required p by only {math.log10(2)/Lr:.4f} -- cosmetic.")
print(f"  => the 3% 'agreement' does NOT test a computed magnitude (which is hypersensitive); it tests whether")
print(f"     the POSTULATED rationals (p=1/4, q=4) are exactly right. All the work is in guessing 1/4 and 4.")

# ---- T3: the 'double 4' (code distance d=4 AND spacetime D=4) is ONE root near 4, not two confirmations ----
# require D=d=n integer:  Lr/n + La*n = const  ->  La*n^2 - const*n + Lr = 0
A, B_, Cc = La, -const, Lr
disc = B_*B_ - 4*A*Cc
roots = [(-B_+s*math.sqrt(disc))/(2*A) for s in (+1,-1)] if disc>=0 else []
print(f"\n[T3] 'd=4 (code distance) AND 1/D=1/4 (spacetime) independently' -> set D=d=n: {La:.3f}n^2 - {const:.3f}n + {Lr:.3f}=0")
print(f"     roots n = {', '.join(f'{r:.2f}' for r in roots)}.  One root sits at n~4 -- so 'both are 4' is a SINGLE")
print(f"     tuned coincidence (the substrate ratio L_H/a0 and alpha put a root near the integer 4), presented")
print(f"     as two independent derivations. (There is also a second root ~{max(roots):.1f}.) The magnitude is one")
print(f"     constraint; two integer inputs hitting it to 3% is ~1-in-{int(1/ (2*dp_for_3pct/1)) if dp_for_3pct else 0}-ish selection, not double confirmation.")

# ---- T4: §16.3 competitor count -- (R, C, p, q) within the 3% claim tolerance ----
rats   = [Fraction(n,d) for d in (1,2,3) for n in range(1,7)]            # small prefactor rationals
irrs   = {"1":1.0, "sqrt2":2**.5, "2":2.0, "sqrt(pi)":math.pi**.5, "sqrt(2pi)":C, "pi":math.pi, "2pi":2*math.pi}
psimple= sorted({Fraction(n,d) for d in range(1,9) for n in range(1,d+1) if Fraction(n,d)<=1}, key=float)
qs     = range(2,7)
TOL    = abs(dev)                                                        # the item's ACTUAL match tolerance (~2-3%)
seen=set()
for R in rats:
    for cn,cv in irrs.items():
        for p in psimple:
            for q in qs:
                v = float(R)*cv*LAM*RATIO**float(p)*(1/ALPHA)**q
                if abs(v-M_P)/M_P <= TOL:
                    seen.add((R, cn, p, q))
n136 = len(seen)
print(f"\n[T4] §16.3 competitor count: simple expressions R*c*Lambda*(L_H/a0)^p*alpha^-q within the item's own")
print(f"     {TOL*100:.1f}% tolerance, over R(rational<=6/3), c in {{1,sqrt2,2,sqrt(pi),sqrt(2pi),pi,2pi}}, p simple<=1, q in 2..6:")
print(f"     competitors found = {n136}")
ex = sorted(seen, key=lambda t:(t[3], float(t[2])))[:6]
print("     e.g.: " + "; ".join(f"{R}*{cn}*L*(ratio)^{p}*a^-{q}" for (R,cn,p,q) in ex))

print(f"""
=========================================================================================
VERDICT (item 136, parameter-free Planck mass via alpha^-4 (L_H/a0)^1/4):
  (T0) the formula reproduces M_P to ~{abs(dev)*100:.1f}% (leaning on Hubble tension for the residual, per the item).
  (T1) mass dimension forces M_Pl = [~Lambda] x [dimensionless]; the magnitude is ONE equation in the two
       exponents (p on L_H/a0, q on 1/alpha) -> a LINE. The framework's OWN item 117 sits at a different point
       (alpha^2, ratio^1) and the text admits 'alpha^4 replaces alpha^2' -> exponents adjustable (internal >=2).
  (T2) hypersensitive: d ln M_P/dp = {dlnM_dp:.0f}, so 1/4 must hold to +/-{dp_for_3pct:.4f}. The 3% 'match' tests the
       GUESSED rationals (1/4, 4), not a computed magnitude; M0=2Lambda and sqrt(2pi) are O(1) cosmetic.
  (T3) 'code-distance 4 AND spacetime 4, independently' = ONE root of the single magnitude equation near n=4
       (second root ~{max(roots):.1f}); two integer inputs meeting one constraint is selection, not confirmation.
  (T4) {n136} equally-simple competitors land within the item's own {TOL*100:.1f}% tolerance -> firmly §16.3 class-1.
  Plus: 1/alpha^4 = {(1/ALPHA)**4:.2e} 'matching the ~1e8 residual' is a reverse-solve -- the residual is
  M_P/(everything else), and 'everything else' contains the tunable (L_H/a0)^1/4; only q=4 lands near 1e8, but
  that is one constrained choice riding on the tuned p=1/4.
  CONCLUSION: item 136 is the WEAKEST of the three gravity routes -- it is §16.3 class-1 numerology with two
  adjustable exponents (the framework's own history exhibits a second assignment), a 3% tolerance propped by
  Hubble tension, and 'derivations' (1/4 conformal dim; alpha^4 from d=4; sqrt(2pi) from CLT) that are
  post-hoc rationalisations of magnitude-selected values. Its own residual-polish list concedes the alpha^4
  rep-theory step is only 'plausible', not derived. It does NOT independently rescue the gravity headline.
  (Fairness: the [8,4,4] code distance IS genuinely 4, so alpha^-4 is a *motivated candidate*, not arbitrary;
  but motivated-candidate + free second exponent + O(1) prefactor + 3% + admitted reassignment = consistency
  check, not parameter-free prediction.)
=========================================================================================""")
assert isinstance(n136, int)
print("exit 0 -- item 136 audited (inputs computed, nothing fitted).")
