#!/usr/bin/env python3
"""
Exact 3-parameter Koide fit (mu, R, delta) to the 3 charged-lepton masses.
Closed-form linear solve -- no RNG, no arccos branch ambiguity, no Monte-Carlo.

Supersedes the removed koide_delta_honest.py, whose "form-misfit spread" was
conceptually wrong: with R free, 3 parameters fit 3 masses EXACTLY, so there is
no misfit. The genuine tests are (a) R_fit vs sqrt2 [= Koide Q=2/3] and
(b) delta_fit vs 2/9, each carrying only the experimental error from m_tau.
"""
import math
M = {"e": 0.51099895000, "mu": 105.6583755, "tau": 1776.86}
dMtau = 0.12   # PDG 1-sigma on m_tau (dominant; e, mu negligible)

def fit(m):
    s = {k: math.sqrt(v) for k, v in m.items()}
    sbar = sum(s.values()) / 3.0               # = sqrt(mu) exactly (sum of cos = 0)
    x = {k: s[k] / sbar - 1 for k in m}        # x_n = R cos(delta + 2pi n/3)
    A = x["tau"]                               # R cos(delta)      (n=0)
    B = (x["mu"] - x["e"]) / math.sqrt(3)      # R sin(delta)      (n=1,2)
    return sbar**2, math.hypot(A, B), math.atan2(B, A), (x["e"] + x["mu"] + x["tau"])

mu, R, delta, sumx = fit(M)
print("Exact 3-param Koide fit to the 3 charged-lepton masses:")
print("  sqrt(mu)      = %.6f MeV^1/2   (sum x_n = %+.2e ~ 0 check)" % (math.sqrt(mu), sumx))
print("  R_fit         = %.7f   vs sqrt2 = %.7f" % (R, math.sqrt(2)))
print("  R_fit - sqrt2 = %+.2e        (this IS the Koide Q=2/3 test)" % (R - math.sqrt(2)))
print("  delta_fit     = %.7f rad" % delta)
print("  2/9           = %.7f rad" % (2 / 9))
print("  delta_fit-2/9 = %+.2e" % (delta - 2 / 9))

d_hi = fit({**M, "tau": M["tau"] + dMtau})[2]
d_lo = fit({**M, "tau": M["tau"] - dMtau})[2]
sig = abs(d_hi - d_lo) / 2
print("\n  sigma(delta) from m_tau +/- %.2f = %.2e" % (dMtau, sig))
print("  |delta_fit - 2/9| / sigma = %.2f" % (abs(delta - 2 / 9) / sig))
print("""
HONEST READING:
  delta_emp = %.6f +/- %.0e (experimental, m_tau-dominated). 2/9 = %.6f is
  %.2f sigma away: CONSISTENT with 2/9, NOT a proof of exactness. The strong,
  robust result is R_fit = sqrt2 to %.0e (Koide Q = 2/3); delta near 2/9 is the
  weaker claim riding on top. Branch (ii) predicts delta is transcendental, so
  sharper m_tau should reveal delta != 2/9 at the ~1e-5 level -- a real test.
""" % (delta, sig, 2 / 9, abs(delta - 2 / 9) / sig, abs(R - math.sqrt(2))))
