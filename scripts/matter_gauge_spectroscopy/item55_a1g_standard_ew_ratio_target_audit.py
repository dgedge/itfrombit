#!/usr/bin/env python3
"""
item55_a1g_standard_ew_ratio_target_audit.py

Section 229 / item 55 target:

    Can the matter-cell A_1g scalar mode plus standard electroweak mixing give
    a parameter-free Higgs/Z ratio?

This audits the sharpest tempting closure route.  At tree level in standard EW
normalisation,

    m_H^2 = 2 lambda_H v^2,
    m_Z^2 = (g^2 + g'^2) v^2 / 4,

so

    (m_H/m_Z)^2 = 8 lambda_H / g_Z^2,
    g_Z^2 = g^2 + g'^2,
    sin^2 theta_W = g'^2 / g_Z^2.

The A_1g construction can identify the scalar direction.  EW mixing identifies
the neutral gauge direction once g and g' are known.  The ratio still needs the
dimensionless scalar/gauge theorem lambda_H/g_Z^2, or equivalently
lambda_H/g^2 plus a same-scale value of sin^2 theta_W.

The attractive rational candidate is:

    sin^2 theta_W = 3/8       (GUT-normalised charge trace)
    lambda_H/g^2 = 3/8        (hypothetical A_1g scalar/gauge service theorem)

which would imply

    lambda_H/g_Z^2 = (3/8)(1 - 3/8) = 15/64,
    m_H/m_Z = sqrt(15/8).

The audit shows that this is a near-hit against pole masses, not a derivation:
it uses a high-scale charge-trace mixing value with low-energy pole masses, and
standard EW mixing alone never supplies lambda_H/g^2.
"""

from __future__ import annotations

import math
import sys
from fractions import Fraction


ok = True


def check(name: str, cond: bool) -> None:
    global ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    ok = ok and bool(cond)


def ratio_from_lambda_over_g2(lambda_over_g2: float, sin2: float) -> float:
    """Tree-level m_H/m_Z from lambda/g^2 and sin^2(theta_W)."""
    return math.sqrt(8.0 * lambda_over_g2 * (1.0 - sin2))


def required_lambda_over_g2(mh_over_mz: float, sin2: float) -> float:
    """Tree-level lambda/g^2 required by a measured m_H/m_Z at a chosen sin^2."""
    lambda_over_gz2 = mh_over_mz * mh_over_mz / 8.0
    return lambda_over_gz2 / (1.0 - sin2)


def pct(x: float) -> float:
    return 100.0 * x


# Same reference inputs already used in the EW audits in this repository.
m_h = 125.25
m_z = 91.1876
m_w = 80.3692
sin2_on_shell = 1.0 - (m_w / m_z) ** 2
sin2_msbar_mz = 0.23122
sin2_gut = Fraction(3, 8)

lambda_over_g2_candidate = Fraction(3, 8)
lambda_over_gz2_candidate = lambda_over_g2_candidate * (1 - sin2_gut)
ratio_candidate = math.sqrt(8.0 * float(lambda_over_gz2_candidate))

ratio_observed = m_h / m_z
lambda_over_gz2_observed = ratio_observed * ratio_observed / 8.0

print("=" * 88)
print("ITEM 55 / SECTION 229: A_1g + STANDARD EW MIXING RATIO AUDIT")
print("=" * 88)
print("\nStandard EW identity:")
print("  (m_H/m_Z)^2 = 8 lambda_H / g_Z^2")
print("  lambda_H/g_Z^2 = (lambda_H/g^2)(1 - sin^2 theta_W)")

print("\nReference pole-mass target used by the repo EW scripts:")
print(f"  m_H = {m_h:.4f} GeV")
print(f"  m_Z = {m_z:.4f} GeV")
print(f"  m_H/m_Z = {ratio_observed:.9f}")
print(f"  required lambda_H/g_Z^2 = {lambda_over_gz2_observed:.9f}")

print("\nTempting rational target:")
print(f"  sin^2 theta_W[GUT trace] = {sin2_gut} = {float(sin2_gut):.9f}")
print(f"  assumed A_1g service theorem lambda_H/g^2 = {lambda_over_g2_candidate}"
      f" = {float(lambda_over_g2_candidate):.9f}")
print(f"  implied lambda_H/g_Z^2 = {lambda_over_gz2_candidate}"
      f" = {float(lambda_over_gz2_candidate):.9f}")
print(f"  implied m_H/m_Z = sqrt(15/8) = {ratio_candidate:.9f}")

ratio_error = ratio_candidate / ratio_observed - 1.0
lambda_error = float(lambda_over_gz2_candidate) / lambda_over_gz2_observed - 1.0
print("\nNear-hit size:")
print(f"  ratio error = {pct(ratio_error):+.4f}%")
print(f"  lambda_H/g_Z^2 error = {pct(lambda_error):+.4f}%")

check("15/64 candidate is genuinely close at the sub-percent level",
      abs(ratio_error) < 0.01)
check("15/64 candidate is not an exact pole-mass closure",
      abs(ratio_error) > 0.002)

print("\nSame-scale EW-mixing audit:")
for label, sin2 in [
    ("GUT charge trace", float(sin2_gut)),
    ("on-shell at M_Z", sin2_on_shell),
    ("MSbar at M_Z", sin2_msbar_mz),
]:
    req = required_lambda_over_g2(ratio_observed, sin2)
    pred = ratio_from_lambda_over_g2(float(lambda_over_g2_candidate), sin2)
    print(f"  {label:16s}: sin^2={sin2:.9f}"
          f"  required lambda/g^2={req:.9f}"
          f"  if lambda/g^2=3/8 -> m_H/m_Z={pred:.9f}")

req_gut = required_lambda_over_g2(ratio_observed, float(sin2_gut))
req_on_shell = required_lambda_over_g2(ratio_observed, sin2_on_shell)
req_msbar = required_lambda_over_g2(ratio_observed, sin2_msbar_mz)
check("using GUT sin^2 makes the required lambda/g^2 look near 3/8",
      abs(req_gut / float(lambda_over_g2_candidate) - 1.0) < 0.01)
check("using low-energy on-shell sin^2 does not support lambda/g^2=3/8",
      abs(req_on_shell / float(lambda_over_g2_candidate) - 1.0) > 0.15)
check("using low-energy MSbar sin^2 does not support lambda/g^2=3/8",
      abs(req_msbar / float(lambda_over_g2_candidate) - 1.0) > 0.15)

print("\nInterpretation:")
print("  The A_1g mode supplies a scalar representation/eigenmode target.  Standard")
print("  EW mixing supplies the neutral-gauge rotation after g and g' are specified.")
print("  Neither operation derives lambda_H/g^2.  The 15/64 candidate becomes close")
print("  only when the high-scale 3/8 charge-trace mixing value is combined with")
print("  low-energy pole masses.  At the same low-energy schemes, the required")
print("  lambda_H/g^2 is about 0.304-0.307, not 3/8.")
print("\nVerdict:")
print("  Current framework status: not closed.  The real theorem target is now sharp:")
print("  derive a single EW service/Hessian or RG-matched boundary condition that fixes")
print("  lambda_H/g_Z^2 directly, or fixes lambda_H/g^2 at the same scale as the")
print("  chosen sin^2 theta_W.  Without that, A_1g + standard EW mixing gives a")
print("  useful structural map and a near rational target, not a parameter-free")
print("  prediction.")

if ok:
    print("\nALL CHECKS PASSED")
    sys.exit(0)

print("\nCHECKS FAILED")
sys.exit(1)
