#!/usr/bin/env python3
"""Derive-or-retire the 5.7 bounds — the M5 decision gate.

5.7 states (underived): kappa in [ln3/2, ln2] ~ [0.549, 0.693], floor 'from diagonal-plane
parity checks', ceiling 'Shannon 2-bit detection evasion'; plus a dynamical selection kappa=phi
via the fixed point x = 1/(1+x). 5.2 fixes the units: M = 1/P with P a survival PROBABILITY,
P = p^F per frustrated edge => kappa = -ln(p), p a probability.

Reconstruction of the two bounds from their stated sources:
  GEOMETRIC: the canonical walk's hop amplitude is a = 1/sqrt(3) per axis (3.2: T_d = -i/sqrt3 ...;
    unitarity check: 3 x |1/sqrt3|^2 = 1). Stated floor ln3/2 = -ln|a| — an AMPLITUDE log.
    In 5.2's probability units the same physics gives p_geom = |a|^2 = 1/3 => kappa_geom = ln3.
  SHANNON: a binary parity check extracts <= 1 bit => survival probability >= 1/2
    => kappa <= ln2 — a PROBABILITY log.
Consistency test: a bound pair constrains one kappa only if both are in the same units.
Self-asserting; exit 0 = every number verified."""
import math

ln = math.log
phi  = (math.sqrt(5)-1)/2          # 0.618...
kII  = 1/(2*phi)                   # Group II exponent 0.809...

a_geom   = 1/math.sqrt(3)          # 3.2 hop amplitude (unitarity: 3*a^2 = 1)
assert abs(3*a_geom**2 - 1) < 1e-15
p_shan   = 0.5                     # binary-check survival floor (probability)

floor_amp,  ceil_amp  = -ln(a_geom), -ln(math.sqrt(p_shan))   # both as amplitude logs
floor_prob, ceil_prob = -ln(a_geom**2), -ln(p_shan)           # both as probability logs
floor_mixed, ceil_mixed = -ln(a_geom), -ln(p_shan)            # the PUBLISHED pair

print("the 5.7 bound pair under each unit convention (kappa = -ln p, M = 1/P per 5.2):")
print(f"  consistent PROBABILITY units: [ln3, ln2]   = [{floor_prob:.4f}, {ceil_prob:.4f}]  -> "
      f"{'EMPTY (floor > ceiling)' if floor_prob > ceil_prob else 'ok'}")
print(f"  consistent AMPLITUDE  units: [ln3/2, ln2/2] = [{floor_amp:.4f}, {ceil_amp:.4f}]  -> "
      f"{'EMPTY (floor > ceiling)' if floor_amp > ceil_amp else 'ok'}")
print(f"  MIXED (published)          : [ln3/2, ln2]   = [{floor_mixed:.4f}, {ceil_mixed:.4f}]  -> "
      f"{'non-empty' if floor_mixed < ceil_mixed else 'empty'}  <- amplitude floor + probability ceiling")
assert floor_prob > ceil_prob and floor_amp > ceil_amp and floor_mixed < ceil_mixed

print(f"""
=> The published interval [0.549, 0.693] exists ONLY by reading the floor in amplitude units
   and the ceiling in probability units. In EITHER consistent convention the interval is EMPTY
   (the geometric attenuation already exceeds the Shannon ceiling). The bounds are therefore not
   derivable as a consistent constraint on a single kappa — they compare incommensurable logs.""")

# the kappa = phi dynamical selection: category check
x = (math.sqrt(5)-1)/2             # fixed point of x = 1/(1+x)
assert abs(x - 1/(1+x)) < 1e-15
print(f"""the 5.7 dynamical selection: x = 1/(1+x) -> x = phi = {x:.4f}; 5.7 then SETS kappa = x.
  x is a probability/amplitude-like ratio (a survival fraction); kappa is a LOG-rate. Equating
  them is a category mismatch. The consistent reading kappa = -ln(x) = {-ln(x):.4f} violates even
  the published (mixed) floor {floor_mixed:.4f} — under no reading does the fixed-point argument
  select a kappa inside its own bounds.""")
assert -ln(x) < floor_mixed

# where the surviving candidates sit (for the record; no constraint applies once bounds retire)
print("candidate exponents for the record (no consistent bound constrains them):")
for name, k in (("Group II 1/(2phi)", kII), ("5.7 phi", phi)):
    print(f"  {name:18s} kappa = {k:.4f}")

print(f"""
DECISION (executing the pre-stated M5 audit rule, DRIFT M5 addendum 2026-06-10):
  'if the bounds are heuristic -> adopt Group II and retire the bounds.'
  The bounds FAIL derivation (unit-inconsistent; empty in any consistent convention), and the
  kappa=phi selection rests on a category mismatch + an uncited V_ub/V_cb scan. Therefore:
  - the 5.7 bounds and the kappa=phi selection are RETIRED (recorded, not deleted);
  - Group II  M(c) = exp(F/(2 phi)),  kappa = 1/(2 phi) = {kII:.5f}, unit-F multiplier
    exp(1/(2 phi)) = {math.exp(kII):.4f}, becomes the CANONICAL inscription — as the sole
    surviving candidate, supported by every numerical contact (9/4 at 0.19%, Yukawa 5.4%,
    cascade 1.002), each individually 16.3 class-1 (evidential status stays modest; the
    INSCRIPTION is canonical, not thereby 'confirmed').
ALL ASSERTS PASSED""")
