#!/usr/bin/env python3
"""
a_e = (g-2)/2 reachability test for §15 item 138 / §5.4.

Question: can the framework's §5.4 (Part 12) rational graph-invariant machinery
produce the QED anomalous moment a_e = (g-2)/2 = alpha/2pi + c2 (alpha/pi)^2 + ...?

Two independent obstructions, both demonstrated here:
  (A) WRONG GREEN'S FUNCTION — §5.4 computes the photon self-energy Pi(q^2)
      (running of alpha); a_e is the electron vertex form factor F2(0).
  (B) RATIONAL vs TRANSCENDENTAL — §5.4 coefficients are in Q (N1=31, C2=-24/7,
      C_k=(-1)^{k-1} 24 (1/7)^{k-1}); the QED c2 is transcendental (pi^2, zeta(3)).

Result: a_e is NOT reachable. Corrects the earlier "comes from §5.4" framing.
Needs only mpmath. Nothing fitted.
"""
import mpmath as mp
mp.mp.dps = 40

pi = mp.pi; ln2 = mp.log(2); z3 = mp.zeta(3)
alpha = 1 / mp.mpf('137.035999177')
x = alpha / pi
a_exp = mp.mpf('0.00115965218059')   # measured electron a_e (Fan et al 2023 / PDG)

# QED mass-independent coefficients (electron a_e), known closed forms:
c1 = mp.mpf(1) / 2                                              # Schwinger 1947, RATIONAL
c2 = mp.mpf(197)/144 + pi**2/12 - (pi**2/2)*ln2 + mp.mpf(3)/4*z3  # Sommerfield/Petermann 1957
c3 = mp.mpf('1.181241456587')                                  # Laporta-Remiddi 1996 (transcendental)

a1 = c1*x
a2 = a1 + c2*x**2
a3 = a2 + c3*x**3

print("="*66)
print("a_e QED tower vs experiment  (a_e^exp = %s)" % mp.nstr(a_exp, 12))
print("="*66)
print("  c1 = 1/2                = %s   [RATIONAL]" % mp.nstr(c1, 6))
print("  c2 = 197/144+pi^2/12-(pi^2/2)ln2+(3/4)z3 = %s   [TRANSCENDENTAL: pi^2, zeta(3)]" % mp.nstr(c2, 10))
print("  c3 = %s   [TRANSCENDENTAL]" % mp.nstr(c3, 10))
print()
print("  1-loop (Schwinger only): a_e = %s   dev = %s %%" % (mp.nstr(a1,12), mp.nstr((a1-a_exp)/a_exp*100,3)))
print("  + 2-loop:                a_e = %s   dev = %s %%" % (mp.nstr(a2,12), mp.nstr((a2-a_exp)/a_exp*100,3)))
print("  + 3-loop:                a_e = %s   dev = %s %%" % (mp.nstr(a3,12), mp.nstr((a3-a_exp)/a_exp*100,3)))
print("  => Schwinger's RATIONAL 1/2 alone is +0.15%% high (a FAILING precision result);")
print("     the 10-digit agreement is supplied entirely by the TRANSCENDENTAL tower.")

print()
print("="*66)
print("OBSTRUCTION (A): wrong Green's function")
print("="*66)
print("  §5.4 Dyson-Schwinger  -> photon self-energy Pi(q^2)  -> running of alpha.")
print("  a_e                   -> electron vertex form factor F2(0)  -> different diagram.")
print("  Deriving alpha does NOT derive a_e.")

print()
print("="*66)
print("OBSTRUCTION (B): rational graph invariants cannot give transcendentals")
print("="*66)
# §5.4 series sums to a rational; show its limit is in Q and != c2.
S = sum(mp.mpf((-1)**(k-1)) * 24 * mp.mpf(1)/7**(k-1) for k in range(1, 40))  # 24/(1+1/7)=21
print("  §5.4 geometric series  sum_k (-1)^{k-1} 24 (1/7)^{k-1} = %s = 24/(8/7) = 21  [RATIONAL]" % mp.nstr(S, 10))
print("  QED 2-loop coefficient c2 = %s  [TRANSCENDENTAL]" % mp.nstr(c2, 10))
print("  A geometric series in 1/7 with integer numerators is rational; c2 is not.")
print("  => §5.4's machinery can match c1=1/2 at most; the transcendental tower is")
print("     unreachable. Same obstruction class as the L-W Koide no-go (item 135 III(b))")
print("     and the §16.3 search-space null.")

print()
print("VERDICT: a_e (the genuine QED precision test) is NOT reachable from §5.4 —")
print("wrong diagram AND rational-vs-transcendental wall. The framework gets tree-level")
print("g=2 (Dirac baseline, item 138) but not the anomalous moment.")
