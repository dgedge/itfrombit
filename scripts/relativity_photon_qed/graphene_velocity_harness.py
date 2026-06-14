#!/usr/bin/env python3
r"""Graphene known-case harness for the Chadha-Nielsen velocity-renormalisation machinery (Item 115 loop).
Before trusting ANY c_4.8.8 sign, validate the self-energy + log-extraction logic on graphene, whose
1-loop velocity renormalisation has a textbook answer (Gonzalez-Guinea-Vozmediano 1994):
    v_F(mu) = v_F(0) [ 1 + (alpha_g/4) ln(Lambda/mu) ]   (velocity GROWS in the IR; coefficient = 1/4).

Model: massless Dirac H0 = v_F (qx sx + qy sy), v_F=1; instantaneous 2D Coulomb V(p) = 2 pi g / |p|.
Exchange (Fock) self-energy:  Sigma(k) = - INT d^2q/(2pi)^2 V(k-q) P_-(q),  P_-(q) = (I - qhat.sigma)/2.
The velocity shift is the coefficient of (khat.sigma)|k| in Sigma; here we take k=(k,0) and read the
sigma_x coefficient Sigma_x(k); delta_v(k) = Sigma_x(k)/k should equal (g/4) ln(Lambda/k).
PASS criterion declared in advance: the slope of Sigma_x(k)/k vs ln(Lambda/k) equals g/4 = 0.25 (g=1).
"""
import numpy as np

g = 1.0          # coupling (alpha_g = g/v_F, v_F=1)
Lam = 1.0        # UV cutoff
# polar grid for the loop momentum q; log-spaced in q to resolve the IR log, fine in theta
Nq, Nth = 4000, 2400
qs = np.exp(np.linspace(np.log(1e-5), np.log(Lam), Nq))
th = (np.arange(Nth)+0.5)/Nth*2*np.pi
dth = 2*np.pi/Nth
COS = np.cos(th)

def sigma_x(k):
    """sigma_x component of the exchange self-energy at external k=(k,0).
       Sigma_x = (g/(4 pi)) INT q dq dtheta cos(theta) / |k - q| ,  |k-q|=sqrt(k^2+q^2-2 k q cos theta)."""
    # integrate over q (trapz in log-spacing => dq = q d(lnq)) and theta
    dlnq = np.log(Lam/1e-5)/(Nq-1)
    tot = 0.0
    for i,q in enumerate(qs):
        denom = np.sqrt(k*k + q*q - 2*k*q*COS)
        ang = np.sum(COS/denom)*dth          # angular integral at this q
        tot += q*ang * q*dlnq                 # q dq = q*(q dlnq); integrand q dq dtheta cos/denom
    return g/(4*np.pi) * tot

ks = np.array([0.01, 0.02, 0.04, 0.08])
dv = np.array([sigma_x(k)/k for k in ks])     # delta_v(k) = Sigma_x(k)/k
logs = np.log(Lam/ks)
slope, intercept = np.polyfit(logs, dv, 1)
print("[graphene harness] static-Coulomb exchange self-energy, velocity renormalisation:")
print(f"  {'k':>7} {'ln(Lam/k)':>10} {'delta_v=Sig_x/k':>16}")
for k,l,d in zip(ks,logs,dv):
    print(f"  {k:>7.3f} {l:>10.3f} {d:>16.4f}")
print(f"\n  fit delta_v = slope*ln(Lam/k) + const :  slope = {slope:.4f}   (textbook g/4 = {g/4:.4f})")
err = abs(slope - g/4)/(g/4)
print(f"  deviation from g/4: {err*100:.1f}%")
ok = err < 0.08
print(f"  => harness {'PASSES' if ok else 'FAILS'}: the self-energy + log-extraction machinery reproduces the "
      f"known graphene coefficient 1/4 to {err*100:.0f}%.")

print(f"""
=========================================================================================
VERDICT (graphene harness):
  The discrete-quadrature exchange self-energy + linear-in-k velocity extraction + log-slope fit reproduces
  the textbook GGV coefficient g/4 = 0.25 (got {slope:.3f}). So the *architecture* (BZ-integration, projector
  onto the filled band, isolating the k-linear velocity shift, extracting the leading-log slope) is VALIDATED
  on a case with a known answer -- the sign (velocity grows toward the IR) and the magnitude both check out.
  CAVEATS carried forward to the 4.8.8 / K6 loop (NOT yet done):
   * Different propagator: graphene uses instantaneous 2D Coulomb V~1/q; the 4.8.8 photon is the RETARDED 3D
     SC propagator 1/(w^2 - q^2) with K(q)=6-2 sum cos. The harness validates the architecture, NOT the
     propagator handling -- the frequency integral d_omega and the 3D q-integral must be added and re-checked.
   * R_photon=1 is an ASSUMPTION (user's step 1), not validated here: the SC web is isotropic at O(q^2) but
     the discrete BZ-cutoff still gives a finite ratio vs the continuum at leading log; needs its own check.
   * The K6 vertex must be the SYMMETRY-CORRECT k.p (T_1u<->E_g, item102_k6_dispersion.py), NOT the user's
     explicit 6x6 (verified this turn to give the WRONG velocities 1.53/1.22/1.15 != sqrt(2/3)/.../...).
  So: machinery validated; three substrate-specific pieces (retarded 3D photon, R_photon, correct k.p vertex)
  remain before a defensible c_4.8.8 sign. No sign is claimed yet.
=========================================================================================""")
assert ok, "harness must reproduce g/4 to <8%"
print("exit 0 -- graphene velocity-renorm harness reproduces the textbook 1/4 coefficient; architecture validated.")
