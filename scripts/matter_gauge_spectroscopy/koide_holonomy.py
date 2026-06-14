#!/usr/bin/env python3
"""
Koide phase delta: holonomy attempt + pi-obstruction theorem.
Tests whether delta_l = 2/9 rad (empirical, 0.02%) can be a geometric phase.
Reproducible. Nothing fitted.
"""
import math, cmath
import numpy as np

def bargmann(states):
    z = 1 + 0j
    n = len(states)
    for i in range(n):
        z *= np.vdot(states[i], states[(i + 1) % n])
    return math.atan2(z.imag, z.real)

def coh(theta, phi):
    return np.array([math.cos(theta / 2), cmath.exp(1j * phi) * math.sin(theta / 2)])

print("TARGET delta_l = 2/9 =", round(2/9, 6), "rad (empirical from PDG, 0.02%)\n")

# (A) continuous Berry/Bargmann phase: spin-1/2 coherent states on x,y,z
X = coh(math.pi/2, 0); Y = coh(math.pi/2, math.pi/2); Z = coh(0, 0)
gA = bargmann([X, Y, Z])
print("(A) closed Bargmann phase, axes x/y/z: %.6f rad = %.4f*pi  (=pi/4: solid-angle/2)"
      % (gA, gA/math.pi))

# (B) Gauss-Bonnet: continuous closed geometric phase always carries pi
print("\n(B) Gauss-Bonnet: closed phase = -Omega/2, Omega = (angle sum) - pi.")
print("    pi-free rational delta would need angle-sum = pi + 2*delta -> carries pi.")
print("    => continuous closed holonomy is pi-LADEN by construction.")

# (C) cyclic single-valued loop (AB-style): 2pi-quantized
print("\n(C) cyclic loop, 3-hop generation ring, single-valuedness => total = 2pi*(d/N):")
for d, N, name in [(2,9,'lepton'),(3,9,'nu'),(1,9,'down'),(2,27,'up')]:
    dN = d/N; ab = 2*math.pi*dN/3
    print("    %-7s target/hop d/N=%.5f  closed/hop=%.5f  ratio=%.5f"
          % (name, dN, ab, ab/dN))
print("    uniform ratio 2pi/3 = %.5f  => every sector off by 2pi (the AB failure)." % (2*math.pi/3))

# (D) is 2/9 a pi-multiple at all?
print("\n(D) nearest pi-multiple to 2/9 (scan p/q, q<50):")
best = None
for q in range(1, 50):
    for p in range(1, 4*q):
        for base, bn in [(2*math.pi,'2pi'), (math.pi,'pi')]:
            v = base*p/q
            if abs(v - 2/9) < 2e-3 and (best is None or abs(v-2/9) < best[0]):
                best = (abs(v-2/9), "%s*%d/%d=%.5f" % (bn, p, q, v))
print("    ", best if best else "none within 0.2% for q<50")

print("""
VERDICT (corrected -- earlier draft OVERCLAIMED a universal obstruction):
  A general closed Berry/Bargmann phase = -Omega/2 and Omega (solid angle) can
  be ANYTHING in [0,4pi); a loop enclosing Omega=4/9 sr gives delta=2/9 pi-free.
  So 'all closed holonomies are pi-obstructed' is FALSE. WITHDRAWN.
  What IS true and robust:
    * Cube-symmetric / natural loops give pi-RATIONAL phases: the octant geodesic
      triangle (axes x,y,z) gives exactly pi/4 (B,A above), NOT 2/9.
    * AB single-valued quantization gives 2pi*(d/N), off by 2pi/3 uniformly (C).
    * To get delta=2/9 from a closed loop you must enclose Omega=4/9 steradian --
      a solid angle with NO cube-symmetric origin (those are pi-rational).
  So the natural geometric candidates are excluded; a bespoke Omega=4/9 loop is
  not excluded but has no substrate motivation (and see the L-W squeeze).

The decisive, holonomy-independent result is the Lindemann-Weierstrass fork:
see koide_LW_squeeze.py. delta=2/9-exact <=> transcendental circulant amplitude,
so 'exact rational 2/9' and 'derived geometric phase from algebraic substrate
amplitudes' cannot both hold. That fork does not depend on any holonomy claim.
""")
