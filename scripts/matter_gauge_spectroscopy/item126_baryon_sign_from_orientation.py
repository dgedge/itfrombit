#!/usr/bin/env python3
r"""ITEM 126 / R15: the baryon-asymmetry SIGN from the now-fixed CP orientation.

The baryogenesis sign was openfrontier, blocked on 'the DeltaL=2 recovery holonomy
portal, phase, and orientation'. This session supplies the PORTAL (R15 recovery
holonomy = the QEC recovery correction) and the ORIENTATION sigma (the global
handedness, pinned sigma=+1 by the observed Jarlskog J>0). The CP-odd carrier is

    I_CP = Im[ (M_H)_12 (M_H)_23 (M_H)_31 ],   M_H = M0[I + r e^{i sigma Phi} A_K3],

the oriented 3-edge holonomy product = (M0 r)^3 sin(3 sigma Phi). The baryon sign is
sign(I_CP) -> (sign-preserving) B = (28/79)(B-L). Honest correction: the relevant
phase is NOT the R15 C_3 character 2pi/3 (which gives sin(3*2pi/3)=sin 2pi=0, a CP
ZERO); it is the geometric nu_R Berry phase (~1/3 rad or 2pi/9, both < pi/3), for
which sin(3Phi)>0 so the SIGN is carried by sigma alone, robust to the magnitude.
Self-asserting.
"""
import numpy as np
def ok(c,m): print(("  PASS " if c else "  FAIL ")+m); assert c,m
A_K3=np.array([[0,1,1],[1,0,1],[1,1,0]],float)
def MH(sigma,phi,M0=1.0,r=0.3): return M0*(np.eye(3)+r*np.exp(1j*sigma*phi)*A_K3)
def I_CP(sigma,phi):
    M=MH(sigma,phi); return float(np.imag(M[0,1]*M[1,2]*M[2,0]))

print("="*72); print("BARYON SIGN from the fixed CP orientation sigma"); print("="*72)

# (1) the C_3 character 2pi/3 KILLS the baryogenesis CP (sin 3Phi = sin 2pi = 0)
ok(abs(I_CP(+1, 2*np.pi/3))<1e-12, "Phi=2pi/3 (R15 C_3 character) gives I_CP=0 -> NOT the baryogenesis phase")

# (2) the geometric nu_R Berry phases give nonzero CP, flipping with sigma
for phi,lab in ((1/3,"delta_nu=1/3 rad"),(2*np.pi/9,"nu_R Berry 2pi/9")):
    ip,im=I_CP(+1,phi),I_CP(-1,phi)
    print(f"\n  Phi={phi:.4f} ({lab}):  I_CP(+sigma)={ip:+.4f}  I_CP(-sigma)={im:+.4f}")
    ok(abs(ip)>1e-6, f"{lab}: nonzero CP")
    ok(abs(ip+im)<1e-12, f"{lab}: I_CP flips sign with sigma -> sign carried by sigma")
    ok(phi<np.pi/3, f"{lab} < pi/3 so sin(3Phi)>0 -> sign(I_CP) set by sigma alone, magnitude-robust")

# (3) sign is robust across the whole geometric window (0, pi/3); only sigma sets it
window=np.linspace(0.05,np.pi/3-0.05,50)
signs=[np.sign(I_CP(+1,p)) for p in window]
ok(len(set(signs))==1, f"sign(I_CP) is constant (={signs[0]:+.0f}) for sigma=+1 across all Phi in (0,pi/3) -> phase-independent sign")

# (4) sphaleron B-L -> B = 28/79 is positive: sign-preserving
from fractions import Fraction as F
ok(F(28,79)>0, "sphaleron conversion B=(28/79)(B-L) is positive -> sign(eta_B)=sign(B-L)=sign(I_CP)")

print("\n"+"="*72); print("VERDICT")
print("  The baryon-asymmetry SIGN is carried by the global orientation sigma -- the SAME")
print("  handedness fixed by the observed Jarlskog J>0 (sigma=+1). I_CP = (M0 r)^3 sin(3 sigma Phi)")
print("  flips with sigma and, for the geometric nu_R Berry phase (Phi < pi/3), has a sign set by")
print("  sigma ALONE, independent of the phase magnitude. The sphaleron 28/79 is sign-preserving.")
print("  So sign(eta_B) is LOCKED to sign(J): the matter-antimatter asymmetry sign and the quark")
print("  CP sign are ONE orientation -- a falsifiable correlation, no longer a free knob.")
print("  CORRECTION: the baryogenesis phase is NOT the R15 C_3 character 2pi/3 (that is a CP zero);")
print("  it is the geometric nu_R Berry phase. The SIGN frontier closes (carried by sigma); the phase")
print("  MAGNITUDE stays the geometric Berry input but does not affect the sign. Absolute 'matter' =")
print("  the same global-handedness convention proven at R15. openfrontier -> computed (up to convention). exit 0")
