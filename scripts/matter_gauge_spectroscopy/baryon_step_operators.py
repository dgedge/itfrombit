#!/usr/bin/env python3
"""
Can the three decuplet 'NOT closed' items be closed? EXPLORATORY (uncommitted).

Method: the De Rujula-Georgi-Glashow (DGG) constituent quark model is the minimal
operator decomposition that QCD actually uses for baryon masses:
    M = sum_i m_i  +  a * sum_{i<j} (S_i . S_j)/(m_i m_j)
The second term is the chromomagnetic HYPERFINE. We fit (m_u, m_s, a) to THREE
masses (N, Delta, Omega) -- no curvature freedom left -- and then ask the model what
kappa and the strange step actually ARE as operators. This decides items (1) and (2)
without any framework assumption; (3) we just pin numerically.

S_i.S_j sums (3 quarks): sum_{i<j} S_i.S_j = [S(S+1) - 9/4]/2.
  octet  S=1/2 -> -3/4 ;  decuplet S=3/2 -> +3/4 (so each aligned pair contributes +1/4).
"""
import mpmath as mp
mp.mp.dps = 25
def H(s): print("\n"+"="*72+"\n"+s+"\n"+"="*72)

MN, MD, MO = mp.mpf('939.0'), mp.mpf('1232.0'), mp.mpf('1672.45')
PDG_DEC = {'Delta':1232.0,'Sigma*':1384.6,'Xi*':1533.4,'Omega':1672.45}
MPI = mp.mpf('139.57'); L = mp.mpf('332.0')

# ---- fit (m_u, a) from N, Delta ; then m_s from Omega ----
# M_D - M_N = (3/2) a/u^2 ; M_D + M_N = 6u
u  = (MN + MD)/6
Au = (MD - MN)*2/3            # = a/u^2
a  = Au * u**2
# Omega: 3s + (3/4)a/s^2 = M_O  -> bisection
f  = lambda s: 3*s + mp.mpf('0.75')*a/s**2 - MO
lo, hi = mp.mpf('400'), mp.mpf('700')
for _ in range(200):
    mid=(lo+hi)/2
    if f(lo)*f(mid)<=0: hi=mid
    else: lo=mid
s = (lo+hi)/2
kappa = (MD - MN)/2
print(f"DGG fit to (N,Delta,Omega):  m_u={float(u):.1f}  m_s={float(s):.1f}  "
      f"m_s-m_u={float(s-u):.1f}  a/u^2={float(Au):.1f}  kappa=(MD-MN)/2={float(kappa):.1f}")

# pair-sum helper in units of 1/u^2 : n_s strange, (3-n_s) light
def Hdec(ns):
    r=u/s
    ll = (3-ns)*(2-ns)//2 if ns<=1 else (1 if ns==1 else 0)  # C(3-ns,2)
    ll = max((3-ns)*(3-ns-1)//2,0); ls=(3-ns)*ns; ss=ns*(ns-1)//2
    return (a/4)*(ll/u**2 + ls/(u*s) + ss/s**2)
def Mdec(ns): return (3-ns)*u + ns*s + Hdec(ns)

H("(1) Is P_s = kappa a structural identity, or a coincidence of two operators?")
print("   In the DGG decomposition the two numbers are DIFFERENT operators:")
print(f"     kappa     = (3/4) a/u^2 = pure N-Delta CHROMOMAGNETIC hyperfine   = {float(kappa):.1f} MeV")
dstep=(MO-Mdec(0))/3
print(f"     dec step  = (m_s-m_u) + (a/4)(1/s^2 - 1/u^2)")
print(f"               = constituent MASS-DIFF {float(s-u):.1f} + hyperfine corr {float((a/4)*(1/s**2-1/u**2)):.1f} = {float(s-u+(a/4)*(1/s**2-1/u**2)):.1f} MeV")
print(f"   kappa is 100%% spin-spin; the strange step is ~90%% mass-difference. They coincide")
print(f"   only because (m_s-m_u) {float(s-u):.0f} - hyperfine {abs(float((a/4)*(1/s**2-1/u**2))):.0f} ~ pure-hyperfine {float(kappa):.0f}: a NUMERICAL accident.")
print(f"   FALSIFICATION via the OCTET (universal-kappa would step by kappa there too):")
print(f"     octet N->Lambda step = {1115.7-939.0:.1f} (empirical) / 166 (item-129 4/8)  vs kappa {float(kappa):.1f}")
print(f"     => octet step != kappa by ~20-30 MeV. 'one tension kappa' is FALSIFIED.")
print("   VERDICT (1): CLOSED -- as REFUTED. Not a structural identity; a numerical near-")
print("   coincidence of distinct operators, non-universal (octet breaks it). The empirical")
print("   postdiction 'decuplet step = kappa to 0.2%' survives, but only as empirical.")

H("(2) The ~14 MeV curvature: which way does leading hyperfine actually bend it?")
steps_pdg=[PDG_DEC['Sigma*']-PDG_DEC['Delta'],PDG_DEC['Xi*']-PDG_DEC['Sigma*'],PDG_DEC['Omega']-PDG_DEC['Xi*']]
steps_dgg=[float(Mdec(1)-Mdec(0)),float(Mdec(2)-Mdec(1)),float(Mdec(3)-Mdec(2))]
print(f"   PDG steps      : {[f'{x:.1f}' for x in steps_pdg]}  (DECREASING; curvature {steps_pdg[0]-steps_pdg[2]:+.0f})")
print(f"   DGG hyperfine  : {[f'{x:.1f}' for x in steps_dgg]}  (INCREASING; curvature {steps_dgg[0]-steps_dgg[2]:+.0f})")
print(f"   DGG predicts intermediate Sigma*/Xi* LOW: "
      f"Sig* {float(Mdec(1)):.0f} (PDG {PDG_DEC['Sigma*']:.0f}), Xi* {float(Mdec(2)):.0f} (PDG {PDG_DEC['Xi*']:.0f})")
print("   => leading chromomagnetic hyperfine bends the ladder the WRONG WAY. The observed")
print("   curvature is a HIGHER-ORDER effect (config mixing / relativistic) that even standard")
print("   QCD's leading correction misses. So the framework's EQUAL-SPACING (linear kappa) is")
print("   the better LEADING description, not a defect.")
print("   VERDICT (2): NOT closed, but REFRAMED -- not a framework-specific failure; the 14 MeV")
print("   curvature is second-order and open for ALL leading models. Closing it needs a genuine")
print("   spin-spin operator (= item 138), with sign/magnitude beyond leading DGG.")

H("(3) Why kappa ~ 146.5 ? M_pi vs Lambda_QCD/2 vs neither")
print(f"   kappa = {float(kappa):.1f}")
print(f"     M_pi          = {float(MPI):.1f}   (kappa is +{float((kappa-MPI)/MPI*100):.1f}% above)")
print(f"     Lambda_QCD/2  = {float(L/2):.1f}   (kappa is {float((kappa-L/2)/(L/2)*100):.1f}% below)")
print(f"   kappa sits BETWEEN M_pi and Lambda/2, ~5% from the pion. The framework's 'kappa = M_pi")
print(f"   (the D=1 meson)' is 5% loose; 'kappa = Lambda/2' is 12% loose. Neither is clean.")
print("   VERDICT (3): NOT closed -- genuinely recurses to item-130 / §9.9 GMOR (why the J^P")
print("   tension sits near M_pi). No new strange mass; the open scale is kappa/M_pi itself.")

H("SUMMARY")
print("""  (1) CLOSED (refuted): P_s=kappa is a coincidence of two distinct operators
      (pure hyperfine vs mass-difference), non-universal -- octet step != kappa.
      Keep the empirical postdiction; drop the 'single tension' interpretation.
  (2) NOT closed, REFRAMED: leading hyperfine bends the curvature the WRONG way, so
      equal-spacing is the better leading description; the 14 MeV residual is 2nd-order
      and open for everyone (needs item-138 spin operator). Not a framework-specific gap.
  (3) NOT closed: recurses to GMOR/item-130 (kappa ~ M_pi, 5% loose). No strange mass.
  NET: 1 of 3 resolvable now (as a refutation), 1 reframed in the framework's favour,
  1 genuinely open and correctly located upstream.""")
