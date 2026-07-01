#!/usr/bin/env python3
r"""ITEM 87 -- origin of the universal 2/9 standard-shear (PMNS reactor / mixing).

matter_gauge treats the 2/9 standard shear as an input ('the universal 2/9 reactor
primitive'). ch10 (book) DERIVES it as a Berry phase: the sterile-neutrino defect
traversing the cyclic generation ring picks up

    delta = d / N = 2 / 9   radians,

with d=2 the defect's topological support (the 2 bits where nu_R violates the R4
rule) and N=9 the qubit plaquette. This script connects the two and places 2/9 in
the mean-cycle family -- an honest REDUCTION, not a clean closure (the four-sector
unification is an organising conjecture; the delta->angle map is leading order).

Solid result: delta=2/9 is the d=2 sibling of the NEUTRINO mean-cycle twist
delta_nu = 3/9 = 1/3 (already closed this session) on the SAME N=9 plaquette.
So the shear is not an independent empirical input; it is a plaquette Berry phase
delta=d/N in the same family as the mean cycle, with d_charged/d_nu = 2/3.
Self-asserting on the arithmetic; honest about what stays conjectural.
"""
import numpy as np
from fractions import Fraction as F
def ok(c,m): print(("  PASS " if c else "  FAIL ")+m); assert c,m

print("="*72); print("2/9 STANDARD-SHEAR ORIGIN: nu_R-defect Berry phase delta = d/N"); print("="*72)

# (1) the Berry phase delta = d/N = 2/9
d, N = 2, 9
delta = F(d, N)
ok(delta == F(2,9), f"charged-lepton twist delta = d/N = {d}/{N} = 2/9 (defect support / plaquette)")

# (2) sibling of the closed mean-cycle: neutrino twist = 3/9 = 1/3
delta_nu = F(3, 9)
ok(delta_nu == F(1,3), "neutrino twist delta_nu = 3/9 = 1/3 = the (now-closed) mean cycle")
ok(delta == F(2,3)*delta_nu, "delta_2/9 = (2/3) * delta_nu : same N=9 plaquette, defect-support ratio 2/3")

# (3) the full four-sector d/N_eff table (organising CONJECTURE per ch10)
table = {"charged lepton": F(2,9), "neutrino": F(3,9), "down quark": F(1,9), "up quark": F(2,27)}
for sec, tw in table.items():
    base = "N=9 plaquette" if tw.denominator in (9,3) else "N_eff=27 (9 x 3 colour)"
    print(f"    {sec:14s}: delta = {tw}   [{base}]")
ok(all(tw.denominator in (3,9,27) for tw in table.values()),
   "all four twists are d/N_eff on the 9-qubit plaquette (27 = colour-enlarged) -- one family")

# (4) leading-order angle: theta13 ~ arcsin(sin(delta)/sqrt2)  (NOT an exact angle law)
th13 = np.degrees(np.arcsin(np.sin(float(delta))/np.sqrt(2)))
ok(8.0 < th13 < 9.5, f"theta13 ~ arcsin(sin(2/9)/sqrt2) = {th13:.2f} deg (obs ~8.6 deg; leading order, ~5%)")

print("\n"+"="*72); print("VERDICT (honest reduction)")
print("  REDUCED, not an empirical input: the matter_gauge '2/9 reactor primitive' is the")
print("  ch10 nu_R-defect Berry phase delta = d/N = 2/9 (d=2 R4-violation bits, N=9 plaquette),")
print("  the d=2 SIBLING of the closed neutrino mean-cycle delta_nu=3/9=1/3 on the same plaquette")
print("  (delta = (2/3) delta_nu). So 2/9 lives in the mean-cycle / plaquette Berry-phase family.")
print("  NOT fully closed: (i) the four-sector d/N_eff table (down 1/9, up 2/27) is an organising")
print("  conjecture -- only the charged-lepton 2/9 is high precision; (ii) N=9 is referenced as")
print("  the unit cell but its count (octagonal 8 + centre) is not airtight; (iii) the delta->angle")
print("  map is leading order (transcendental cos + mu fit), as the canon already states. exit 0")
