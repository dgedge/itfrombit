#!/usr/bin/env python3
r"""ITEM 132 -- the cored-profile JEANS CONSISTENCY issue.

Two open sub-issues in dark_sector "Cored Profiles and Jeans Targets":
  (A) r_c selection: central-harmonic one-a0 rule gives r_c=r_M/3 (0.333),
      but the exact enclosed-field condition gives r_c/r_M=1-pi/4 (0.2146).
      Which (if either) is the dynamically consistent core radius?
  (B) Jeans support: a literal constant-tension Jeans model is retired.
      Is there a consistent (variable, physical) support, and does it close?

Method: compare BOTH candidate rules against the ACTUAL MOND phantom density
that a baryonic point mass produces (simple and standard interpolation), and
solve the anisotropic Jeans/static-stress equation for the cored profile.
Dimensionless throughout: u=r/r_M, y=g_N/a0=1/u^2, x=g/a0. Self-asserting.
"""
import math
import numpy as np

def check(ok, msg):
    print(("  PASS " if ok else "  FAIL ")+msg); assert ok, msg

print("="*72); print("CORED-PROFILE JEANS CONSISTENCY"); print("="*72)

# ---------- (A) the two r_c rules, exactly ----------
print("\n[A] r_c selection: inner-approx (1/3) vs exact enclosed-field (1-pi/4)")
# cored profile rho=A/(r^2+r_c^2): M(r)=4piA[r-r_c atan(r/r_c)], g_h=GM/r^2
# 4 pi G A = v_inf^2 ; condition g(r_c)=a0, with a0 = v_inf^2 / r_M.
# exact:  g_h(r_c) = (v_inf^2/r_c)(1-pi/4) = a0  -> r_c/r_M = 1-pi/4
rc_exact = 1.0 - math.pi/4.0
# inner quadratic: small-r g_inner=(4piGA/3 r_c^2) r = (v_inf^2/(3 r_c^2)) r
#   g_inner(r_c) = v_inf^2/(3 r_c) = a0 -> r_c/r_M = 1/3
rc_inner = 1.0/3.0
print(f"   exact enclosed-field  r_c/r_M = 1-pi/4 = {rc_exact:.4f}")
print(f"   inner-approx (canon)  r_c/r_M = 1/3    = {rc_inner:.4f}")
# the discrepancy is exactly the small-r approximation error of g_h at r=r_c:
#   g_inner(r_c)/g_h(r_c) = (1/3)/(1-pi/4)
ratio = rc_inner/rc_exact
check(abs(ratio - (1/3)/(1-math.pi/4)) < 1e-12,
      f"1/3 overshoots exact by {ratio:.3f}x = the quadratic-approx error at r=r_c")
# verify the inner expansion of the EXACT g_h reproduces the 1/3 coefficient:
# g_h(r) = (v_inf^2/r)[1-(r_c/r)atan(r/r_c)]; small x=r/r_c: 1-(1/x)atan(x) -> x^2/3
# so g_h -> v_inf^2 r/(3 r_c^2) = g_inner. Confirms 1/3 is the r<<r_c limit.
x=1e-3
inner_coeff = (1.0 - (1.0/x)*math.atan(x))/(x*x)   # -> 1/3 as x->0
check(abs(inner_coeff - 1/3) < 1e-5, "exact g_h has the 1/3 coefficient ONLY in the r<<r_c limit")

# ---------- (B) the ACTUAL MOND phantom density of a point mass ----------
print("\n[B] actual MOND phantom density: where does the real core sit?")
def x_simple(y):   # mu=x/(1+x): x^2-yx-y=0
    return 0.5*(y+np.sqrt(y*y+4*y))
def x_standard(y): # mu=x/sqrt(1+x^2): x^2/sqrt(1+x^2)=y -> solve
    # x^4 - y^2 x^2 - y^2 = 0 -> x^2=(y^2+sqrt(y^4+4y^2))/2
    return np.sqrt(0.5*(y*y+np.sqrt(y**4+4*y*y)))
u = np.logspace(-3, 3, 400001)
y = 1.0/(u*u)
for name, xf in (("simple", x_simple), ("standard", x_standard)):
    x = xf(y)
    # phantom rho_ph * r^2  proportional to  f(u)=d(u^2 x)/du  (asymptote -> 1)
    f = np.gradient(u*u*x, u)
    check(abs(f[-1]-1.0) < 0.02, f"[{name}] rho_ph r^2 -> A (isothermal tail) as u->inf")
    check(f[0] < 0.05, f"[{name}] phantom vanishes at center (baryon-dominated core)")
    # effective core: half-asymptote point f(u_c)=1/2  (matches pseudo-iso def)
    idx = np.argmin(np.abs(f-0.5))
    uc = u[idx]
    print(f"   mu={name:9s}: effective core u_c=r_c/r_M = {uc:.4f}  (f=1/2 half-asymptote)")
    # also peak of rho_ph (= peak of f/u^2)
    upk = u[np.argmax(f/(u*u))]
    print(f"               phantom-density peak at u = {upk:.4f}")

# ---------- (C) Jeans support: variable anisotropic tension ----------
print("\n[C] Jeans support: solve required tau(r) for the cored profile")
# anisotropic static stress p_r=-tau*rho, p_t=0:  tau' + [2 r_c^2/(r(r^2+r_c^2))] tau = g_h
# dimensionless x=r/r_c, g_h=(v_inf^2/r)[1-atan(x)/x]:
#   tau/v_inf^2 = (1+x^2)/x^2 * INT_0^x  u/(1+u^2) [1-atan(u)/u] du
xx = np.logspace(-3, 3, 200001)
integ_grid = (xx/(1+xx*xx))*(1.0-np.arctan(xx)/xx)
cum = np.concatenate([[0.0], np.cumsum(0.5*(integ_grid[1:]+integ_grid[:-1])*np.diff(xx))])
def tau_over_v2(xq):
    I = np.interp(xq, xx, cum); return (1+xq*xq)/(xq*xq)*I
for xq in (0.1, 1.0, 10.0, 100.0):
    print(f"   tau({xq:6.1f} r_c)/v_inf^2 = {tau_over_v2(xq):.4f}")
check(tau_over_v2(1e-2) < 1e-3, "tau -> 0 at center (regular, ~x^2/12): support is well-defined")
check(tau_over_v2(100.0) > tau_over_v2(10.0) > 0, "tau positive and growing (NOT constant): variable support")
g_log = tau_over_v2(100.0)-tau_over_v2(10.0)
check(abs(g_log-math.log(10.0)) < 0.2, "large-r tau grows as v_inf^2 ln(r): the log-tension support")

print("\n"+"="*72)
print("VERDICT")
print("  (A) 1/3 and 1-pi/4 are the SAME condition g(r_c)=a0 read two ways:")
print("      1-pi/4=0.215 uses the EXACT halo field (self-consistent global core);")
print("      1/3=0.333 is its r<<r_c quadratic limit, so it is a LOCAL central-cell")
print("      statement, NOT the global enclosed-field boundary. They must not be equated.")
print("  (B) the real MOND phantom core (half-asymptote) lands between them (~0.2-0.5 r_M),")
print("      i.e. neither rule is exact for a point mass; the pseudo-iso r_c is a regulator.")
print("  (C) Jeans support CLOSES with a variable anisotropic tension")
print("      tau(r)=v_inf^2[..]->v_inf^2 ln(r/r_c), regular (tau(0)=0), positive.")
print("      The retired 'constant-tension' ansatz was simply the wrong EOS. exit 0")
