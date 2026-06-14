#!/usr/bin/env python3
"""
Confirm the CORRECTED §5.4 / Part 12 Dyson-Schwinger equation and pin down the
exact canonical form, before any edit to ANCHOR/DRIFT.

The published 2-loop term  (24/7)*(1/(2 pi alpha^-1))^2  is QUADRATIC in
alpha^-1 and yields 137.0360037 (the 1-loop value). The headline numbers come
from the term being LINEAR in alpha^-1. The natural loop expansion makes this
unambiguous: writing delta = alpha^-1 - 137,

    delta = (N1/2pi) alpha - (|C2|/(2pi)^2) alpha^2 + (|C3|/(2pi)^3) alpha^3 - ...

i.e. the k-loop term is order alpha^k in delta -> the 2-loop term must be
O(alpha^2), which forces ONE power of alpha (= 1/alpha^-1) in the RHS term:

    alpha^-1 (alpha^-1 - 137) = N1/(2pi)
                                - (|C2|/(2pi)^2) (1/alpha^-1)
                                + (|C3|/(2pi)^3) (1/alpha^-1)^2 - ...
    with N1=31, |C2|=24/7, |C3|=24/49 (C_k = 24/7^{k-1}, geometric ratio 1/7).

This script verifies that THIS form reproduces ALL THREE quoted numbers:
  1-loop 137.036004 ; 2-loop 137.035999077 ; 3-loop 137.035999078 .
It also checks whether the '1/7 geometric series' claim is about the
COEFFICIENTS C_k (true) vs the loop-to-loop SHIFTS (carries an extra alpha/2pi).
Self-asserting. exit 0 == corrected form reproduces all three quoted values.
"""
import mpmath as mp
mp.mp.dps = 40
twopi = 2*mp.pi
N1 = mp.mpf(31)

def solve(nloops):
    """Solve corrected DS truncated at nloops. term_k (k>=2) = (-1)^(k-1) (24/7^(k-1))/((2pi)^k (a)^(k-1))."""
    a = mp.mpf('137.036')
    for _ in range(300):
        rhs = N1/twopi
        for k in range(2, nloops+1):
            Ck = mp.mpf(24)/mp.mpf(7)**(k-1)
            rhs += (-1)**(k-1) * Ck/(twopi**k * a**(k-1))
        an = 137 + rhs/a
        if abs(an-a) < mp.mpf(10)**(-35): a = an; break
        a = an
    return a

x1 = solve(1); x2 = solve(2); x3 = solve(3)
print("CORRECTED (linear-power) Dyson-Schwinger equation, truncations:")
print(f"  1-loop : {mp.nstr(x1,13)}   quoted 137.036004")
print(f"  2-loop : {mp.nstr(x2,13)}   quoted 137.035999077")
print(f"  3-loop : {mp.nstr(x3,13)}   quoted 137.035999078")

# Loop-to-loop shifts and the role of the 1/7 ratio.
s12 = x1 - x2; s23 = x2 - x3
print(f"\n  shift 1->2 = {mp.nstr(s12,4)} ;  shift 2->3 = {mp.nstr(s23,4)}")
print(f"  term-shift ratio |s12/s23| = {mp.nstr(abs(s12/s23),5)}")
alpha = 1/mp.mpf('137.036')
print(f"  predicted term-shift ratio = 7 / (alpha/2pi) = {mp.nstr(7/(alpha/twopi),5)}")
print( "  => the loop-to-loop TERM ratio is 7*(2pi/alpha) ~ 6000, NOT 7. The '1/7")
print( "     geometric series' is a statement about the COEFFICIENTS C_k=24/7^(k-1)")
print( "     (correct); each higher loop also carries an extra alpha/2pi, so the")
print( "     three quoted numbers ARE mutually consistent. [self-correction: an")
print( "     earlier audit compared term-ratio to 7 and wrongly called it inconsistent.]")

# Asserts: corrected form reproduces all three quoted values.
assert abs(x1 - mp.mpf('137.036004'))    < mp.mpf('5e-7'), "1-loop mismatch"
assert abs(x2 - mp.mpf('137.035999077')) < mp.mpf('5e-10'), "2-loop mismatch"
assert abs(x3 - mp.mpf('137.035999078')) < mp.mpf('5e-10'), "3-loop mismatch"
print("\nexit 0 == corrected linear-power equation reproduces all three quoted numbers")
print("           (137.036004 / 137.035999077 / 137.035999078) to the quoted digits.")
