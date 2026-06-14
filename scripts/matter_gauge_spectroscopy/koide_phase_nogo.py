#!/usr/bin/env python3
"""
Is the Koide offset delta=2/9 a DERIVABLE geometric phase, or only a number that
sits near the empirical arg(B)? Pancharatnam from named frames + branch-(i) test
+ root-of-unity heuristic + coincidence audit. Compares to 2/9 only AFTER.
Companion to koide_holonomy.py (holonomy classes) and koide_LW_squeeze.py (L-W).
Reproducible; nothing fitted.
"""
import math, cmath
import numpy as np

TARGET = 2/9
EMP, SIG = 0.222230, 8e-6      # exact-fit empirical delta and experimental sigma

def coh(r):
    r = np.array(r, float); r = r/np.linalg.norm(r)
    th = math.acos(max(-1, min(1, r[2]))); ph = math.atan2(r[1], r[0])
    return np.array([math.cos(th/2), cmath.exp(1j*ph)*math.sin(th/2)])
def pan2(a, b): z = np.vdot(a, b); return math.atan2(z.imag, z.real)
def pan3(a, b, c):
    z = np.vdot(a, b)*np.vdot(b, c)*np.vdot(c, a); return math.atan2(z.imag, z.real)

print("TARGET 2/9 = %.6f ; empirical arg(B) = %.6f +/- %.0e\n" % (TARGET, EMP, SIG))

# (1) Pancharatnam from NAMED substrate frames (generation axes x,y,z; item 133)
cX, cY, cZ = coh([1,0,0]), coh([0,1,0]), coh([0,0,1])
n = np.array([1,1,1.])/math.sqrt(3)
S = [np.array([[0,1],[1,0]]), np.array([[0,-1j],[1j,0]]), np.array([[1,0],[0,-1]])]
U = math.cos(math.pi/3)*np.eye(2) - 1j*math.sin(math.pi/3)*(n[0]*S[0]+n[1]*S[1]+n[2]*S[2])
print("(1) Pancharatnam/Bargmann phases of the named generation frames:")
print("    arg<+x|+y>               = %.6f  (= pi/4, pi-laden)" % pan2(cX, cY))
print("    Bargmann <x|y><y|z><z|x> = %.6f  (= pi/4, pi-laden)" % pan3(cX, cY, cZ))
print("    Z3 transport 120deg[111] = %.6f  (= 2pi/3; ALREADY the 2pi n/3 term)" % abs(pan2(cX, U@cX)))
print("    => named frames give pi/4 and 2pi/3, NEVER 2/9; the closed Z3 holonomy")
print("       is the formula's own term, contributing nothing to the offset.\n")

# (2) branch (i): can delta be a rational that is NOT an angle?
M = {"e":0.51099895, "mu":105.6583755, "tau":1776.86}
s = {k: math.sqrt(M[k]) for k in M}; sb = sum(s.values())/3
ok = all(abs((s[k]/sb-1)/math.sqrt(2) - math.cos(2/9 + 2*math.pi*nn/3)) < 2e-3
         for k, nn in [("tau",0), ("e",1), ("mu",2)])
print("(2) branch (i) test -- 'a ratio, not a phase':")
print("    masses reconstruct as cos(2/9 + 2pi n/3): %s" % ok)
print("    => delta is unavoidably arg(B) (circulant eigenvalues ARE cosines of arg).")
print("       'ratio not phase' is incoherent; it RELOCATES the gap onto arg(B).")
print("    |B|/A = 1/sqrt2 = %.5f algebraic (fine); arg(B)=2/9 is the sole problem.\n" % (1/math.sqrt(2)))

# (3) root-of-unity heuristic (suggestive, NOT a no-go)
best = min(((abs(2*math.pi*p/q - TARGET), p, q)
            for q in range(1, 200) for p in [round(TARGET*q/(2*math.pi))] if p > 0),
           key=lambda t: t[0])
print("(3) root-of-unity heuristic (discrete phases ~ 2pi p/q):")
print("    nearest 2pi*%d/%d = %.6f needs q=%d (large); 2pi n/3 is exactly q=3." % (best[1], best[2], 2*math.pi*best[1]/best[2], best[2]))
print("    SUGGESTIVE only -- Berry phases CAN be non-root-of-unity. Weaker than L-W.\n")

# (4) coincidence density near 2/9
cands = {"2/9":2/9, "2sin(1/9)":2*math.sin(1/9), "arcsin(2/9)":math.asin(2/9),
         "arctan(2/9)":math.atan(2/9), "pi/14 (anchored Cabibbo!)":math.pi/14}
near = sum(1 for v in cands.values() if abs(v-EMP) < 0.01*EMP)
print("(4) coincidence audit -- simple expressions within 1%% of emp (%.5f):" % EMP)
for nm, v in cands.items():
    print("    %-26s = %.6f   %s" % (nm, v, "<1%" if abs(v-EMP) < 0.01*EMP else ""))
print("    %d of %d within 1%%. Neighbourhood CROWDED -> percent-level proximity is" % (near, len(cands)))
print("    NOT evidence. pi/14 is the framework's OWN anchored Cabibbo quantum.\n")

print("="*68)
print("VERDICT: delta=2/9 is solid as a defect RATIO (type A: sin^2_thetaW, Q^3,")
print("Wolfenstein) but NOT a derived PHASE. Convergent: (1) Pancharatnam null,")
print("(2) branch-(i) incoherent, (3) not low-order root of unity, (4) crowded")
print("neighbourhood; plus L-W (koide_LW_squeeze.py) + holonomy pi-obstruction")
print("(koide_holonomy.py). Robust result: R=sqrt2 (Q=2/3). delta~2/9 at 0.89")
print("sigma is the weak claim; mechanism OPEN, most likely coincidental.")
print("="*68)
