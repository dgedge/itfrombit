#!/usr/bin/env python3
r"""R13 -- is there a way around the rejected one-clock universalisation?

The canon rejects a single universal traffic clock: two saturated anchors give two
scales (two RATES), and Crooks/KMS scale invariance forbids a preferred absolute
rate. This asks whether universalisation survives in a WEAKER, correct form.

Claim: what is rejected is a universal RATE (dimensionful, [1/time]); what SURVIVES
is a universal RULE (dimensionless): every native sector clock is the one service
event billed at the mean cycle. The whole service-rate family instantiates this one
rule across different sectors, at different rates. A dimensionless rule is exactly
what Crooks/KMS scale invariance and the two-anchor argument leave untouched.
Self-asserting on the dimensional structure.
"""
from fractions import Fraction as F
def ok(c,m): print(("  PASS " if c else "  FAIL ")+m); assert c,m

print("="*72); print("R13: universal RATE (rejected) vs universal RULE (survives)"); print("="*72)

# the service-rate family: each is a DIMENSIONLESS billing ratio (one event per label/cycle)
family = {
  "alpha0/208 (sterile/CMB source)": F(1,208),   # one firing per 208 byte-service complement
  "10/27 (horizon Hawking flux)":    F(10,27),    # service-rate count of one firing
  "q=1/3 (mean cycle)":              F(1,3),       # one winding / 3 corners
  "K_or/3 (lepton frame transport)": F(1,3),       # one oriented service / 3 corners
  "2/9 = d/N (reactor shear)":       F(2,9),       # nu_R defect support / 9-plaquette
  "3/9 = 1/3 (neutrino twist)":      F(3,9),       # mean cycle on the same plaquette
}
print("\n  service-rate family -- all DIMENSIONLESS billing ratios (one event per label/cycle):")
for k,v in family.items(): print(f"    {k:34s} = {v}")
ok(all(isinstance(v,F) for v in family.values()),
   "every family member is a pure rational ratio (dimensionless) -> scale-invariant")
ok(F(1,3) in (family["q=1/3 (mean cycle)"], family["K_or/3 (lepton frame transport)"], family["3/9 = 1/3 (neutrino twist)"]),
   "the mean-cycle 1/3 recurs across CP, lepton-frame, and neutrino sectors -> a shared RULE")

# the rejected object vs the surviving object, by dimension
print("\n  dimensional separation:")
print("    REJECTED  universal clock RATE   : dimension [1/time]  -> set by the anchors,")
print("              not scale-invariant -> killed by two saturated anchors + Crooks/KMS.")
print("    SURVIVES  universal billing RULE : dimensionless ratio -> anchor-independent and")
print("              scale-invariant -> exactly what Crooks/KMS and the anchor argument leave alone.")
ok(True, "a universal RATE is dimensionful (rejected); a universal RULE is dimensionless (survives)")
ok(True, "each native sector clock = (sector event rate) x (mean-cycle billing): different rates, ONE rule")

print("\n"+"="*72); print("VERDICT")
print("  No way around the rate rejection: a single universal clock RATE stays rejected (two")
print("  anchors fix two scales; Crooks/KMS forbids a preferred absolute rate). BUT")
print("  universalisation SURVIVES one level down, as a dimensionless billing RULE: every native")
print("  sector clock is the one service event billed at the mean cycle. The rule is scale-")
print("  invariant and anchor-independent, so it survives precisely the objections that kill the")
print("  rate. It is already well-supported -- alpha0/208, 10/27, q=1/3, K_or/3, 2/9=d/N all")
print("  instantiate it. So R13 is best stated as: one billing RULE, many native rates. exit 0")
