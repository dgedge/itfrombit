#!/usr/bin/env python3
"""
COUNTERFACTUAL for ANCHOR item 138: is a_mu (or the differential a_mu - a_e)
reachable, given item 138 already retracted a_e? Tests the user's reading (b):
the per-generation Hamming-weight differential might give a_mu - a_e even if neither
absolute is reachable.

Canon-consistent framing (item 138 sub-notes, lines ~3867/3869/3871):
  - item 138 already established a_e is NOT reachable: (A) WRONG GREEN'S FUNCTION
    (sec5.4 computes the self-energy Pi(q^2)/running-alpha, not the vertex F2(0) that
    IS the anomalous moment), and a 'rational-vs-transcendental' point;
  - but item 138 ITSELF softened the transcendental point: the zeta(3) proof-of-route
    (zeta_lattice.py) shows a lattice iterated integral DOES reach periods -- the wall
    is 'shortcut-specific, not substrate-fundamental'. So transcendentality per se is
    NOT the binding obstruction; the BINDING one is (A) no F2(0) vertex route exists
    (the period-route is machinery-class only; it does not supply the vertex graph).
This script asks what the DIFFERENTIAL adds, and whether it escapes (A). Needs mpmath.
"""
import mpmath as mp
mp.mp.dps = 40
pi=mp.pi; ln2=mp.log(2); z3=mp.zeta(3)
alpha=1/mp.mpf('137.035999177'); x=alpha/pi
ae_exp=mp.mpf('0.00115965218059')                 # electron a_e (Fan 2023)
amu_exp=mp.mpf('116592057e-11')                   # muon a_mu (FNAL 2023)
mmu_me=mp.mpf('206.7682830')                      # mass ratio
phi=(mp.sqrt(5)-1)/2                              # framework reciprocal golden ratio 0.618

def hd(t): print("\n"+"="*70+"\n"+t+"\n"+"="*70)

hd("(0) the target difference")
diff=amu_exp-ae_exp
print(f"  a_mu^exp - a_e^exp = {mp.nstr(diff,8)}  (~6.3e-6); precision bar wants 1e-9")

hd("(1) the UNIVERSAL QED tower is identical for e and mu -> it CANCELS in the difference")
c1=mp.mpf(1)/2
c2=mp.mpf(197)/144+pi**2/12-(pi**2/2)*ln2+mp.mpf(3)/4*z3
print(f"  c1=1/2, c2={mp.nstr(c2,8)}, ... are MASS-INDEPENDENT (same number for e and mu).")
print(f"  => a_mu - a_e gets ZERO from the universal tower. Reading (b)'s cancellation is REAL:")
print(f"     a_e's universal-tower point does NOT, by itself, hit the difference.")

hd("(2) but the RESIDUAL a_mu - a_e is mass-dependent QED + NON-PERTURBATIVE hadronic")
lead = x**2 * ( mp.mpf(1)/3*mp.log(mmu_me) - mp.mpf(25)/36 )   # e-loop vac-pol in the mu vertex
print(f"  leading mass-dep 2-loop QED (e-loop VAP): (a/pi)^2[(1/3)ln(mmu/me) - 25/36]")
print(f"     = {mp.nstr(lead,8)}  -> {mp.nstr(lead/diff*100,4)}% of the whole difference")
print(f"     it carries ln(mmu/me)={mp.nstr(mp.log(mmu_me),8)} and dilogs at higher order")
print(f"     -> a PERIOD/iterated-integral object (item 138 zeta(3) class: reachable IN PRINCIPLE")
print(f"        only via the integral that yields F2(0) -- which is exactly obstruction (A), unbuilt).")
hvp_TI=mp.mpf('6.845e-8'); hvp_BMW=mp.mpf('7.075e-8')
print(f"  hadronic VP (NON-PERTURBATIVE; the 4.2sigma lives here -- a NEW obstruction vs a_e):")
print(f"     TI dispersive = {mp.nstr(hvp_TI,4)} ; BMW lattice = {mp.nstr(hvp_BMW,4)} ; split ~{mp.nstr(hvp_BMW-hvp_TI,3)}")
print(f"     this is beyond even the period-reaching machinery class: it is a non-perturbative")
print(f"     QCD spectral integral, the very quantity lattice/dispersive QCD itself splits on.")

hd("(3) ingredient (i): the Hamming-weight differential is a COUNT (rational label), not an integral")
r = mp.log(mmu_me)/phi
print(f"  ln(mmu/me)/phi = {mp.nstr(r,8)}  -> not a clean integer/half-integer weight difference;")
print(f"  and the physical mass ratio carries the sec5.4/5.5 macroscopic dressing, not bare phi*F.")
print(f"  A Hamming-weight tier is a defect COUNT; the difference is a period-integral + a")
print(f"  non-perturbative integral. Not 'transcendentality is forbidden' (item 138 softened that)")
print(f"  -- it is the WRONG KIND of object: a count is not the F2(0) vertex integral.")

hd("VERDICT  (canon-consistent with item 138)")
print("""  (A) WRONG DIAGRAM -- BINDING and unchanged: a_mu, like a_e, is the vertex form factor
      F2(0); sec5.4 computes the self-energy Pi(q^2). Item 138 rejected the 'items 93+79
      supply the vertex' framing, and its zeta(3) proof-of-route is machinery-class only --
      it does NOT supply the F2(0) vertex graph. No vertex route exists for e OR mu.
  (B) THE DIFFERENTIAL DOES NOT ESCAPE:
      - reading (b) is RIGHT that the universal tower cancels in a_mu - a_e;
      - but the residual is 93% a mass-dependent QED period (needs the unbuilt F2(0) integral,
        i.e. obstruction A) + a NON-PERTURBATIVE hadronic VP (the 4.2sigma piece) that is a
        NEW obstruction not present for a_e and beyond the period-reaching machinery class;
      - ingredient (i)'s Hamming-weight differential is a rational COUNT, the wrong kind of
        object for a period-integral + non-perturbative integral.
  => a_mu (absolute AND differential) is NOT reachable at the 1e-9 precision bar -- SAME
     STATUS as a_e (item 138). Binding obstruction: no F2(0) vertex route (A); for the
     difference, additionally the non-perturbative HVP. Reading (a) is essentially correct;
     reading (b)'s universal-tower cancellation is real but insufficient. The transcendental
     'wall' is NOT invoked as fundamental (item 138 softened it); the no-go is wrong-diagram
     + non-perturbative, not impossibility-of-periods.""")
print("\nNOT reachable at the precision bar. (Counterfactual done BEFORE committing time, per the brief.)")
