#!/usr/bin/env python3
r"""§16.3 + structural audit of the cosmological-constant row (ANCHOR §15 item 123 / 118 / §10.7).
Headline: rho_Lambda = (3/4) alpha Lambda^4 a_0 / (4 pi L_H) ~ 2.4e-47 GeV^4, "~4%" vs observed ~2.5e-47,
billed as resolving the 10^121 cosmological-constant problem. There is ALSO a Zel'dovich form
rho_Lambda = 9 alpha^2 Lambda^3 H_0 (L2346). We check: the numbers; the horizon-consuming Dirac structure;
the alpha^1-vs-alpha^2 internal inconsistency; the Friedmann relation to the Part-20 M_P; and §16.3.
"""
import sys, math
from fractions import Fraction

ALPHA=1/137.036; LAM=0.332
H0=(67.36e3/3.0856776e22)*6.582119e-25      # GeV
a0=1/LAM; L_H=1/H0                            # GeV^-1  (a0 = hbar c/Lambda = 1/Lambda)
OmL=0.6847
RHO_OBS=2.5e-47                               # GeV^4 (Omega_L * rho_crit; framework uses ~2.5e-47)
print(f"inputs: alpha=1/{1/ALPHA:.2f} Lambda={LAM} H0={H0:.3e} a0/L_H={a0/L_H:.3e}(=H0/Lambda={H0/LAM:.3e})")

# ---------- (1) the two formulas ----------
rho_123 = (3/4)*ALPHA*LAM**4*a0/(4*math.pi*L_H)     # item 123, alpha^1
rho_zel = 9*ALPHA**2*LAM**3*H0                       # Zel'dovich, alpha^2
print(f"\n[1] item 123  (3/4)alpha Lambda^4 a0/(4pi L_H) = {rho_123:.3e}  ({(rho_123-RHO_OBS)/RHO_OBS*100:+.1f}% vs obs)")
print(f"    Zel'dovich  9 alpha^2 Lambda^3 H0          = {rho_zel:.3e}  ({(rho_zel-RHO_OBS)/RHO_OBS*100:+.1f}% vs obs)")
print(f"    -> the SAME framework gives TWO rho_Lambda formulas with DIFFERENT alpha-powers (alpha^1 vs alpha^2).")

# ---------- (2) both are coeff * Lambda^3 * H0 -- HORIZON-CONSUMING (same as M_P) ----------
print(f"\n[2] both reduce to (coeff) x Lambda^3 x H0 (since a0/L_H = H0/Lambda):")
print(f"    item 123: coeff = 3/(16 pi) alpha = {3/(16*math.pi)*ALPHA:.4e} (alpha^1)")
print(f"    Zel'dovich: coeff = 9 alpha^2          = {9*ALPHA**2:.4e} (alpha^2)")
print(f"    rho_Lambda ~ Lambda_QCD^3 * H0 EATS THE HORIZON H0 (cosmological input) -- the Lambda<->horizon")
print(f"    Dirac coincidence, structurally identical to the M_P routes (DRIFT G7). NOT horizon-free.")

# ---------- (3) Friedmann: the alpha^2 (Zel'dovich) form IS the Part-20 M_P restated, not independent ----------
MP2_part20 = 24*math.pi*ALPHA**2*LAM**3/(H0*OmL)
rho_from_MP = 3*H0**2*OmL*MP2_part20/(8*math.pi)
print(f"\n[3] Friedmann rho_Lambda = 3 H0^2 OmL M_P^2/(8pi). Substituting the Part-20 M_P^2={MP2_part20:.3e}:")
print(f"    -> rho_Lambda = {rho_from_MP:.3e} = 9 alpha^2 Lambda^3 H0 (the Zel'dovich form) EXACTLY.")
print(f"    So the alpha^2 rho_Lambda is the Part-20 M_P expressed via Friedmann -- NOT an independent")
print(f"    prediction; it inherits the M_P routes' §16.3 exposure (asserted alpha^2). The alpha^1 item-123")
print(f"    form is a *different*, *worse*-fitting formula that CONFLICTS with it.")
assert abs(rho_from_MP-rho_zel)/rho_zel < 1e-9

# ---------- (4) which alpha-power does the magnitude prefer? ----------
need = RHO_OBS/(LAM**3*H0)                            # required coeff on Lambda^3 H0
print(f"\n[4] required coeff on Lambda^3 H0 = rho_obs/(Lambda^3 H0) = {need:.4e}. By alpha-power:")
for n in (0,1,2,3):
    c=need/ALPHA**n
    print(f"      alpha^{n}: coeff = {c:.4e}  ({'~'+str(round(c)) if 0.5<c<200 else 'tiny/huge'})")
print(f"    -> alpha^2 gives coeff ~ 9 (clean integer); alpha^1 gives ~0.065 (item 123's 3/16pi={3/(16*math.pi):.4f} is "
      f"{abs(3/(16*math.pi)-need/ALPHA)/(need/ALPHA)*100:.0f}% off). So the magnitude SELECTS alpha^2; the canonical")
print(f"    item-123 alpha^1 form is the INFERIOR fit and needs the '3/4' to be dressed toward agreement.")

# ---------- (5) §16.3 competitor count (coeff * alpha^n * Lambda^3 H0 within ~4%) ----------
tol=0.04
comps=[]
for n in (1,2,3):
    for q in range(1,25):
        for p in range(1,4*q):
            for extra,el in [(1,''),(math.pi,'pi'),(1/(4*math.pi),'/4pi'),(1/math.pi,'/pi')]:
                v=Fraction(p,q)*extra*ALPHA**n
                if v>0 and abs(float(v)-need)/need<=tol:
                    comps.append((f"{Fraction(p,q)}{el}*a^{n}", float(v), n))
comps=sorted(set((c[0],round(c[1],8),c[2]) for c in comps), key=lambda x:x[1])
ns=sorted(set(c[2] for c in comps))
print(f"\n[5] §16.3: simple coeff*alpha^n hitting rho_Lambda within {tol*100:.0f}%: {len(comps)} across alpha-powers {ns}")
print(f"    e.g.: {', '.join(c[0] for c in comps[:6])}")

print(f"""
=========================================================================================
VERDICT (cosmological-constant rho_Lambda = (3/4)alpha Lambda^4 a0/(4pi L_H)):
  * rho_Lambda ~ (coeff) Lambda_QCD^3 H0 is the SAME Lambda<->horizon DIRAC COINCIDENCE as the gravity M_P,
    HORIZON-CONSUMING (H0 is a cosmological INPUT) -- not an independent intrinsic prediction. In fact the
    alpha^2 (Zel'dovich) form IS the Part-20 M_P via Friedmann (verified exactly), so it is M_P restated,
    inheriting all of DRIFT G7's §16.3 exposure.
  * INTERNAL INCONSISTENCY: the framework asserts BOTH an alpha^1 form (item 123, (3/16pi)alpha) AND an
    alpha^2 form (Zel'dovich, 9 alpha^2) for the same rho_Lambda. They differ by ~a factor alpha and are
    numerically degenerate near alpha=1/137; the magnitude actually prefers alpha^2 (clean coeff ~9), making
    the canonical item-123 alpha^1 the inferior fit -- the '3/4 rule-class projection' is a coefficient
    fudge that dresses the worse alpha^1 form toward the ~4% match.
  * §16.3: {len(comps)} simple coeff*alpha^n competitors land within ~4% (across alpha^1/2/3) -- the alpha-power
    and coefficient are NOT pinned; the prefactor carries little standalone weight.
  * "Resolves the 10^121 problem": the 10^121 (or ~10^47 for QCD) suppression is supplied ENTIRELY by the
    horizon/lattice ratio a0/L_H = H0/Lambda ~ 1e-42 -- i.e. by inserting the cosmological input, the Dirac
    large-number, not by a derived mechanism.
  NET: the cosmological-constant 'prediction' is the Lambda_QCD<->horizon Dirac coincidence (horizon-input),
  the alpha^2 branch is the M_P restated (not independent), the alpha^1 branch conflicts with it and is the
  worse fit dressed by 3/4, and the alpha-power/coefficient are §16.3-open. Honest claim: "rho_Lambda comes
  out at the right SCALE via the Lambda<->horizon coincidence (horizon supplied), with the alpha-power and
  prefactor asserted" -- the same downgrade as the gravity headline, of which this is a Friedmann image.
=========================================================================================""")
print("exit 0 -- rho_Lambda audited: horizon-consuming Dirac coincidence = M_P-via-Friedmann; alpha^1-vs-alpha^2 inconsistency; §16.3-open.")
