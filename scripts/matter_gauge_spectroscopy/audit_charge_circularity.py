#!/usr/bin/env python3
"""
CIRCULARITY AUDIT of ANCHOR §2.8 electric-charge formula + the §14 "anomaly
cancellation (sum Q=0, sum Q^2=16)" claim. EXPLORATORY (uncommitted), 2026-05-30.

Worry (from the outside-auditor pass): are Q = (1/2)Z_f + (1/3) sum Z_ci + (1/2)Z_p,
weights {1/2,1/3,1/2}, REVERSE-ENGINEERED to make the SM charges (and sum Q=0,
sum Q^2=16) come out -- or forced by independent upstream structure?

Method: separate INPUT (chosen) from DERIVED (forced) by counterfactual freedom-count.
Self-asserting where exact. Needs only fractions + numpy.

TWO self-corrections folded in (caught during the audit):
 (a) sum Q^2 needs the 15-Weyl-per-generation count (each Weyl once, nu_R excluded):
     16/3 per generation x3 = 16. An 8-Dirac count gives 8/3 x3 = 8 (wrong convention).
 (b) for sum Q, the right-handed states must enter as LEFT-HANDED CONJUGATES with
     CONJUGATE charges (u_R^c has Q=-2/3, not +2/3). Q^2 is insensitive to this, but
     sum Q is not -- so the construction below uses proper conjugate charges.
"""
from fractions import Fraction as F
import numpy as np
def H(t): print("\n"+"="*76+"\n"+t+"\n"+"="*76)

# canonical §2.8 assignments: (Z_f, sum_Zc, Z_p) -> Q_SM
rows={'u':((+1,-1,+1),F(2,3)),'d':((-1,-1,+1),F(-1,3)),
      'e':((-1,-3,+1),F(-1)),'nu':((+1,-3,+1),F(0)),'eplus':((+1,+3,-1),F(1))}
def Q(zf,zc,zp,a,b,c): return a*zf+b*zc+c*zp

H("(1) Are the weights {1/2,1/3,1/2} FORCED by the charge assignments, or fitted?")
M=np.array([[1,-1,1],[-1,-1,1],[-1,-3,1]],float); y=np.array([2/3,-1/3,-1.0])
abc=np.linalg.solve(M,y)
print(f"  solve (a,b,c) from u,d,e rows: a={abc[0]:.4f} b={abc[1]:.4f} c={abc[2]:.4f}  (= 1/2,1/3,1/2)")
ok=all(abs(float(Q(*rows[p][0],F(1,2),F(1,3),F(1,2)))-float(rows[p][1]))<1e-12 for p in rows)
print(f"  all 5 rows reproduced by {{1/2,1/3,1/2}}: {ok}")
print("  READING: 3 assignments fix 3 coefficients exactly (square solve); nu and e+ are then")
print("  genuine over-determined checks. Weights are FORCED by the Z-assignments. Question")
print("  moves upstream: are the Z-ASSIGNMENTS themselves free?")
assert ok and np.allclose(abc,[0.5,1/3,0.5])

H("(2) Are the Z-assignments independently motivated, or chosen to hit Q?")
print("  Z_f = sgn(I3): I3 is register bit c5 (§2.1), fixed upstream (weak isospin). INDEP of Q. CLEAN.")
print("  Z_p = CPT tag: §2.7 bitwise-NOT involution, global Z2 of the code. INDEP of Q. CLEAN.")
print("  one-hot colour Z_ci (lepton |00|-> -3, quark -> -1): the LOAD-BEARING choice; it is")
print("  what makes integer-lepton vs third-quark charge work. Canon §2.8 ALREADY flags")
print("  '(the F2-closure -> one-hot dual mapping ... is open)'. Disclosed, not hidden.")

H("(3) Does sum Q = 0 / sum Q^2 = 16 survive the proper 15-Weyl/generation count?")
def q(zf,zc): return F(1,2)*zf+F(1,3)*zc+F(1,2)          # matter (Z_p=+1)
def qbar(Qval): return -Qval                              # LH conjugate of a RH state
gen=[]
# Q_L doublet, 3 colours: u_L=+2/3, d_L=-1/3  (6 states)
for _ in range(3): gen += [q(+1,-1), q(-1,-1)]
# RH singlets entered as LH conjugates, 3 colours: u_R^c=-2/3, d_R^c=+1/3  (6 states)
for _ in range(3): gen += [qbar(q(+1,-1)), qbar(q(-1,-1))]
# L doublet: nu_L=0, e_L=-1  (2 states)
gen += [q(+1,-3), q(-1,-3)]
# e_R^c = +1  (1 state)
gen += [qbar(q(-1,-3))]
sumQ=sum(gen); sumQ2=sum(x*x for x in gen)
print(f"  15-Weyl/generation (RH as LH conjugates): count={len(gen)}")
print(f"    sum Q   = {sumQ}   (charge neutrality of a generation)")
print(f"    sum Q^2 = {sumQ2} = {float(sumQ2):.4f} = 16/3")
print(f"  x3 generations: sum Q = {3*sumQ}, sum Q^2 = {3*sumQ2} (= {float(3*sumQ2):.0f}) over {3*len(gen)} codewords")
print("  => BOTH hold under the 15-Weyl/gen convention. sum Q=0 is robust (charge neutrality);")
print("     sum Q^2=16 is correct but CONVENTION-LADEN -- should be quoted WITH '15 Weyl/gen,")
print("     nu_R excluded', not as a bare number. (8-Dirac/gen would give 8; the convention matters.)")
assert sumQ==0 and 3*sumQ2==16 and len(gen)==15

H("VERDICT -- §2.8 + anomaly circularity audit: GREEN-AMBER")
for ln in [
 "NOT the damaging circularity the assessment feared; foundations substantially clean,",
 "with ONE disclosed soft spot and ONE labelling tighten.",
 "",
 " * weights {1/2,1/3,1/2}: FORCED by the assignments (exact solve), not independently fitted.",
 " * Z_f, Z_p: genuine UPSTREAM inputs (§2.1 register, §2.7 CPT) -- not reverse-engineered. CLEAN.",
 "   => the 'charges reverse-engineered' worry is DISCONFIRMED for Z_f, Z_p, and the weights.",
 " * colour one-hot (the -3/-1 split giving 1/3 quark charge): the one real load-bearing",
 "   assumption -- and canon ALREADY flags its F2-closure->one-hot derivation OPEN. So the",
 "   quark 1/3 charge is honestly CONDITIONAL on that disclosed-open step.",
 " * sum Q = 0 per generation: FORCED by colour-triplet + lepton-doublet structure. ROBUST.",
 " * sum Q^2 = 16: CORRECT, but convention-laden (15 Weyl/gen, nu_R excluded). Quote with",
 "   the convention. (A first-draft mis-count -- 8-Dirac/gen -> 8, and particle-charges for",
 "   conjugates -- was self-caught; the proper 15-Weyl LH-conjugate count gives the canon 16.)",
 "",
 "RECOMMENDED canon action (minor, honest-labelling, NOT a foundational collapse):",
 " (a) note the quark 1/3 charge is conditional on the open colour-one-hot step (§2.8);",
 " (b) state the sum Q^2=16 counting convention at the §14 anomaly row.",
 "The audit DISCONFIRMS the reverse-engineering worry for the weights and 2 of 3 Z-inputs,",
 "and CONFIRMS one disclosed-open assumption (colour normalisation) that canon already owns.",
]: print(ln)
print("\nALL ASSERTS PASSED.")
