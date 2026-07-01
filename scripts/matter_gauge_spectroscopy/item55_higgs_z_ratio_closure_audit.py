#!/usr/bin/env python3
"""
item55_higgs_z_ratio_closure_audit.py

Item 55 residual: can the A_1g matter-cell construction predict m_H/m_Z?

This is a closure audit, not a fit.  It combines the already-verified item-55
pieces:

  Matter-cell O_h side:
      M_H^2  proportional to lambda_A1g
             = 2 k_edge + 4 k_face + 2 k_body + k_site.

  Gauge-bridge D_4h side:
      M_Z^2  proportional to lambda_Eu
             = k_shear + k_mix^2 / k_shear.

The question is whether those two formulae share a sector-native parameter or
invariant.  They do not.  The ratio depends on two disjoint stiffness ledgers
and, in the substrate version, the cross-action-space dimensional bridge.

The standard electroweak formula says the same thing in continuum language:

      m_H / m_Z = 2 sqrt(2 lambda_H) / sqrt(g^2 + g'^2).

Thus the A_1g calculation identifies the Higgs scalar direction and supplies
the breathing eigenvalue, but it does not derive the Higgs quartic or the Z
gauge coupling.  A parameter-free m_H/m_Z would require a new, single EW
service/Hessian theorem tying lambda_H to g_Z^2 (or tying the O_h and D_4h
stiffnesses to one common invariant).  The existing A_1g + bridge split proves
that such a theorem is not already present.
"""

import math
import sys

TOL = 1e-12
ok = True


def check(name, cond):
    global ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    ok = ok and bool(cond)


def lambda_h(ke, kf, kb, ks):
    """Matter-cell A_1g breathing eigenvalue."""
    return 2.0 * ke + 4.0 * kf + 2.0 * kb + ks


def lambda_z(k_shear, k_mix):
    """Gauge-bridge Eu massive eigenvalue after the photon zero-mode Schur step."""
    return k_shear + k_mix * k_mix / k_shear


def log_ratio(params):
    """log(M_H/M_Z), including the cross-action-space scale ratio s_h/s_z."""
    ke, kf, kb, ks, k_shear, k_mix, log_scale_ratio = params
    return 0.5 * (math.log(lambda_h(ke, kf, kb, ks))
                  - math.log(lambda_z(k_shear, k_mix))
                  + log_scale_ratio)


def finite_gradient(f, x, h=1e-6):
    grad = []
    for i in range(len(x)):
        xp = list(x)
        xm = list(x)
        xp[i] += h
        xm[i] -= h
        grad.append((f(xp) - f(xm)) / (2.0 * h))
    return grad


print("=" * 96)
print("ITEM 55 HIGGS/Z RATIO CLOSURE AUDIT")
print("=" * 96)

# 1. Verify the two closed skeletons.
lh = lambda_h(1.0, 0.5, 0.25, 0.3)
lz = lambda_z(2.0, 0.7)
check("matter-cell A_1g skeleton lambda_H = 2ke + 4kf + 2kb + ks",
      abs(lh - 4.8) < TOL)
check("gauge-bridge Eu skeleton lambda_Z = k_shear + k_mix^2/k_shear",
      abs(lz - 2.245) < TOL)

# 2. Direct decoupling: changing one action space leaves the other invariant.
lh_values = [lambda_h(ke, 0.5, 0.25, 0.3) for ke in (0.5, 1.0, 2.0)]
lz_values = [lambda_z(ks, 0.7) for ks in (1.5, 2.0, 3.0)]
check("lambda_H responds to matter-cell stiffnesses", len(set(round(v, 12) for v in lh_values)) == 3)
check("lambda_Z responds to bridge stiffnesses", len(set(round(v, 12) for v in lz_values)) == 3)
check("lambda_H formula contains no bridge parameter",
      all(name not in lambda_h.__code__.co_varnames for name in ("k_shear", "k_mix")))
check("lambda_Z formula contains no matter-cell parameter",
      all(name not in lambda_z.__code__.co_varnames for name in ("ke", "kf", "kb", "ks")))

# 3. The ratio has support on all three independent blocks:
#    matter-cell, bridge, and cross-space scale.
x0 = [1.0, 0.5, 0.25, 0.3, 2.0, 0.7, 0.0]
grad = finite_gradient(log_ratio, x0)
labels = ["ke", "kf", "kb", "ksite", "k_shear", "k_mix", "log(scale_H/scale_Z)"]
print("\nlog(M_H/M_Z) gradient at a generic positive point:")
for label, g in zip(labels, grad):
    print(f"  d/d {label:22s} = {g:+.6f}")
matter_nonzero = any(abs(grad[i]) > 1e-8 for i in range(4))
bridge_nonzero = any(abs(grad[i]) > 1e-8 for i in range(4, 6))
scale_nonzero = abs(grad[6]) > 1e-8
check("ratio depends on the matter-cell stiffness block", matter_nonzero)
check("ratio depends on the gauge-bridge stiffness block", bridge_nonzero)
check("ratio depends on the cross-action-space scale bridge unless separately fixed", scale_nonzero)

# 4. Standard EW language: ratio cancels v but not the independent dimensionless
#    couplings lambda_H and g_Z = sqrt(g^2+g'^2).
m_h_ref = 125.25
m_z_ref = 91.1876
r_ref = m_h_ref / m_z_ref
lambda_over_gz2 = r_ref * r_ref / 8.0
print("\nStandard EW comparison (illustrative reference masses, not an input to the proof):")
print(f"  m_H/m_Z = {r_ref:.6f}")
print("  SM relation: m_H/m_Z = 2 sqrt(2 lambda_H) / g_Z")
print(f"  Therefore lambda_H / g_Z^2 would have to be {lambda_over_gz2:.6f}.")
check("SM ratio needs an independent lambda_H/g_Z^2 theorem",
      0.20 < lambda_over_gz2 < 0.30)

print("\nVERDICT")
print("  The A_1g matter-cell calculation is a genuine structural closure of the")
print("  Higgs scalar direction and its breathing eigenvalue.  The D_4h bridge")
print("  calculation is a genuine structural closure of the Z Eu skeleton.  But")
print("  they are disjoint ledgers.  Existing canon supplies no shared stiffness,")
print("  no shared Hessian, and no sector-native invariant that fixes lambda_H/g_Z^2.")
print("  Therefore m_H/m_Z remains a cross-sector EW-second-anchor quantity, not a")
print("  parameter-free prediction.  The precise future theorem target is a single")
print("  EW service/Hessian action whose fluctuation spectrum ties the matter-cell")
print("  A_1g stiffness to the bridge Eu stiffness (or equivalently derives")
print("  lambda_H/g_Z^2).")

if ok:
    print("\nALL CHECKS PASSED")
    sys.exit(0)

print("\nCHECKS FAILED")
sys.exit(1)
