#!/usr/bin/env python3
r"""Locked-tier sweep, hadronic/particle side (NOT horizon-consuming, unlike the cosmology claims).
Three claims, verified on their own terms -- fairly, since these may actually hold up:
  (A) anomaly cancellation: sum_45 Q = 0, sum_45 Q^2 = 16
  (B) nucleon mass M_N = 2 sqrt2 Lambda
  (C) Koide tau-mass: m_tau predicted from m_e, m_mu via Q = 2/3
"""
import sys, math
from fractions import Fraction

# ============================ (A) anomaly cancellation ============================
# SM Weyl fermions per generation (15, no nu_R), electric charges Q:
gen = ([Fraction(2,3)]*3 + [Fraction(-1,3)]*3 + [Fraction(-2,3)]*3 + [Fraction(1,3)]*3
       + [Fraction(0)] + [Fraction(-1)] + [Fraction(1)])     # u^3, d^3, uc^3, dc^3, nu, e, ec
assert len(gen)==15
states45 = gen*3
sumQ  = sum(states45)
sumQ2 = sum(q*q for q in states45)
print("(A) anomaly cancellation over the 45 SM Weyl fermions:")
print(f"    sum Q   = {sumQ}   (claim 0)")
print(f"    sum Q^2 = {sumQ2}  (claim 16; = 16/3 per generation x 3)")
assert sumQ==0 and sumQ2==16
print(f"    VERIFIED. But this is the STANDARD SM/SO(10) anomaly-freedom (the 16-spinor minus nu_R): a")
print(f"    correct CONSISTENCY CHECK that the framework's fermion content matches the SM -- NOT a novel")
print(f"    parameter-free prediction (the SM is anomaly-free by construction; sum Q^2=16/3 per gen is fixed")
print(f"    by the standard charges). Solid + real, but 'reproduces the SM', not 'predicts something new'.")

# ============================ (B) nucleon mass M_N = 2 sqrt2 Lambda ============================
LAM=332.0  # MeV, framework chiral scale
# C_8 cycle-graph eigenvalues lambda_k = 2 cos(2 pi k/8); k=1,7 are the first harmonic
l1 = 2*math.cos(2*math.pi*1/8); l7 = 2*math.cos(2*math.pi*7/8)
MN = (l1+l7)*LAM
MN_obs = (938.272+939.565)/2
print(f"\n(B) nucleon mass: C_8 eigenvalues lambda_1=lambda_7=2cos(pi/4)={l1:.4f}")
print(f"    M_N = (lambda_1+lambda_7) Lambda = 2 sqrt2 Lambda = {MN:.2f} MeV  vs obs {MN_obs:.2f}  "
      f"({(MN-MN_obs)/MN_obs*100:+.2f}%)")
print(f"    The sqrt2 IS genuine graph theory (C_8 first harmonic). BUT two inputs:")
print(f"     * the factor '2' = LINEAR superposition of the 2 modes -- an ANSATZ the canon ITSELF admits is")
print(f"       'not uniquely forced' (L1779/item 47): orthogonal-quadrature would give a different factor")
print(f"       (mesons use orthogonal-quadrature -> sqrt2*phi; baryons use linear -> 2sqrt2). Different")
print(f"       ansatz per particle, chosen to fit each.")
print(f"     * the scale Lambda=332 (chiral input). M_N/Lambda = {MN_obs/LAM:.4f} -> needs Lambda~{MN_obs/(2*math.sqrt(2)):.1f} MeV for exact 2sqrt2.")
# §16.3: how unique is 2sqrt2 as M_N/Lambda?
need=MN_obs/LAM
comps=[]
for label,val in [("2sqrt2",2*math.sqrt(2)),("sqrt8",math.sqrt(8)),("e",math.e),("17/6",17/6),
                  ("sqrt2*2",2*math.sqrt(2)),("3-1/6",3-1/6),("2+5/6",2+5/6),("sqrt(8.0)",2.828)]:
    if abs(val-need)/need<0.02: comps.append((label,val))
print(f"    §16.3: simple values within 2% of M_N/Lambda={need:.3f}: {[c[0] for c in comps]} (2sqrt2 is fairly clean)")
print(f"    NET: real graph-theory sqrt2 + admitted-ansatz factor 2 + chiral scale -> SEMI-derived. The 0.02%")
print(f"    match is genuine but the factor 2 rests on the per-particle ansatz choice (item 47, open).")

# ============================ (C) Koide tau-mass ============================
me, mmu = 0.51099895, 105.6583755           # MeV (PDG, essentially exact)
def koideQ(ms): return sum(ms)/(sum(math.sqrt(m) for m in ms))**2
# predict m_tau from m_e, m_mu at Q=2/3: (se+smu+st)^2 = (3/2)(se^2+smu^2+st^2); solve for st (larger root)
se,smu=math.sqrt(me),math.sqrt(mmu)
# (3/2)(se^2+smu^2+st^2) - (se+smu+st)^2 = 0 -> 0.5 st^2 -2(se+smu) st + [0.5(se^2+smu^2)-2 se smu]=0
A=0.5; B=-2*(se+smu); C=0.5*(se**2+smu**2)-2*se*smu
st=(-B+math.sqrt(B*B-4*A*C))/(2*A)
mtau_pred=st**2
print(f"\n(C) Koide charged-lepton: Q = (sum m)/(sum sqrt m)^2")
for lbl,mt in [("PDG-old 1776.86",1776.86),("Belle II 2023 1777.09",1777.09),("world-avg ~1776.97",1776.97)]:
    print(f"    Q[{lbl}] = {koideQ([me,mmu,mt]):.6f}  (2/3 = {2/3:.6f})")
print(f"    m_tau PREDICTED from m_e,m_mu at Q=2/3: {mtau_pred:.3f} MeV")
print(f"      vs Belle II 1777.09+-0.14: {(mtau_pred-1777.09)/0.14:+.1f} sigma ; vs world-avg ~1776.97: "
      f"{(mtau_pred-1776.97)/0.13:+.1f} sigma")
# is Q=2/3 automatic for the (1+sqrt2 cos theta)^2 form? test across random theta
auto=True; tested=0
for i in range(629):                      # theta in [0,2pi)
    th=i*0.01
    roots=[1+math.sqrt(2)*math.cos(th+2*math.pi*n/3) for n in range(3)]
    if all(r>0 for r in roots):           # physical regime (sqrt m = 1+sqrt2 cos, not |.|)
        if abs(koideQ([r**2 for r in roots])-2/3)>1e-9: auto=False
        tested+=1
print(f"    Q=2/3 is an ALGEBRAIC IDENTITY of (1+sqrt2 cos)^2 in the all-positive-root regime (the leptons':")
print(f"    fitted theta~0.222 gives roots 2.38/0.58/0.034, all >0): holds for all {tested} valid theta = {auto}")
assert auto
print(f"    So Q=2/3 is NOT a prediction -- it is automatic for ANY (1+sqrt2 cos)^2 parametrisation. The")
print(f"    framework's content is (i) deriving R=sqrt2 (the circulant amplitude) and (ii) m_tau from m_e,m_mu.")
print(f"    The m_tau prediction IS Koide's 1981 formula -- a real, forward-testable, currently-matching")
print(f"    relation (0 sigma at world avg). Genuinely the strongest of the three; but it is Koide's relation")
print(f"    grounded via R=sqrt2 (whose own derivation is the load-bearing piece), not a new mass formula.")

print(f"""
=========================================================================================
VERDICT (Locked-tier hadronic sweep) -- these HOLD UP markedly better than the cosmology claims:
  (A) anomaly sum Q=0, sum Q^2=16: VERIFIED EXACT, but it is the standard SM/SO(10) anomaly-freedom
      reproduced -- a correct consistency check, not a novel prediction (the canon's own §14 note says so).
  (B) M_N = 2 sqrt2 Lambda: 0.02% match, the sqrt2 IS real graph theory (C_8 first harmonic), NOT a Dirac
      coincidence and NOT horizon-consuming -- but the factor 2 is an admitted, not-forced per-particle
      ansatz (item 47) and the scale is the chiral input. SEMI-derived.
  (C) Koide tau = 1776.97: a real, forward-testable, 0-sigma-matching relation -- but it is Koide's 1981
      formula (Q=2/3 is functional-form-automatic), grounded via R=sqrt2; the value-add is R=sqrt2's
      derivation, not the mass relation itself.
  NET: UNLIKE the cosmology sector (Dirac coincidences + ideal-code coincidences), the hadronic Locked tier
  contains GENUINE content -- an exact SM consistency check, a real graph-theory mass with one ansatz, and a
  real matching mass relation. The honest calibration is 'verified/semi-derived/real-but-standard', NOT the
  headline-vs-footnote over-claim pattern. This is where the framework's substance actually lives.
=========================================================================================""")
print("exit 0 -- (A) exact SM consistency check; (B) semi-derived (real sqrt2 + ansatz 2); (C) Koide's relation via R=sqrt2.")
